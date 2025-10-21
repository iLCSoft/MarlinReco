#ifndef RUNKUTODESOLVER_H
#define RUNKUTODESOLVER_H 1

#include "Colours.h"
#include <cstdlib>
#include <math.h>

// Include Marlin
#include <streamlog/streamlog.h>

namespace sistrip {

#define MAXSTEPS 1000
#define FCE2(x, y) (_pTemplate->*_fce)(x, y)

//! Double precision first-order ODE (IVP, i.e. initial-value problem) solver,
//! that utilizes 4th order Runge-Kutta algorithm with or without adaptive
//! control of its stepsize or so-called 5th(4th) order Runge-Kutta-Fehlberg
//! algorithm with adaptive control of its stepsize (The template class represents
//! class whose method corresponds to a right-hand side FCE2(x,y) of the ODE,
//! defined as dy/dx = FCE2(x,y)!). The template class doesn't have *.cpp file!
//!
//! @see RKCSolve
//! @see RKFSolve
//! @see _eps
//!
//! @author Z. Drasal, Charles University Prague
//!

template <class T>
class RunKutODESolver {

public:
  //! Constructor - sets ODE right-hand side and required absolut accumulated precision
  RunKutODESolver(double (T::*fce)(double, double), T* pTemplate, float eps)
      : _fce(fce), _pTemplate(pTemplate), _eps(eps) {
    ;
  }

  //! Destructor
  ~RunKutODESolver() { ; }

  //! This method returns a solution y(x) to the so-called ODE initial value problem
  //!(IVP), using 4th order Runge-Kutta method. Where IVP is defined as follows:
  //!
  //!   &nbsp;&nbsp; y' = dy/dx = f(x,y) with initial condition f(x0) = y0
  //!
  //!  Calculation (w(x,h) corresponds to the best approximation of y(x)):
  //!
  //!   &nbsp;&nbsp; w_0   = y0
  //!
  //!   &nbsp;&nbsp; for i = 0,1,2 ... (i-th step) <br>
  //!   &nbsp;&nbsp; w_i+1 = w_i + h*Phi(x_i,w_i;h;f) <br>
  //!   &nbsp;&nbsp; x_i+1 = x_i + h, x_0 = x0
  //!
  //!   &nbsp;&nbsp; Euler method gives Phi(x,y;h;f) = f(x,y)
  //!
  //! This algorithm doesn't implement adaptive stepsize control, each stepsize is
  //! calculated as relative precision * (x-xo). (Epsilon ~ 0.001, i.e 1000 iterations,
  //! gives approximately good results.)
  //!
  //! For more details see J.Stoer, R.Bulirsch: Introduction to Numerical Analysis,
  //! chap. 7, 3rd ed., and <a href=http://www.nr.com/oldverswitcher.html>www.nr.com</a>
  //!
  double RKSolve(double x0, double y0, double x);

  //! This method returns a solution y(x) to the so-called ODE initial value problem
  //!(IVP), using 4th order Runge-Kutta method with adaptive stepsize control. The
  //! stepsize control algorithm is performed as follows: first, an approximative
  //! solution w(x0+h;h) using stepsize h is calculated and then is compared to the
  //! solution w(x0+h;h/2), obtained in the same way, but using two consequent
  //! stepsizes h/2. If the difference is not similar to the required epsilon,
  //! then the ratio of a new proposed stepsize and the old stepsize must be much
  //! bigger than 2 (i.e. take for example bigger or equal than 3), else it's
  //! approximately equal to 2. Hence, in both cases the stepsize is set as 2 times
  //! stepsize H, that corresponds to required precision, but in the first case, one
  //! must perform the algorithm once again, while in the latter the approximative
  //! solution w_i+1 and a stepsize for new iteration are found. The epsilon
  //! corresponds to precision in each step!!! The IVP problem is defined as:
  //!
  //!   &nbsp;&nbsp; y' = dy/dx = f(x,y) with initial condition f(x0) = y0
  //!
  //!  Calculation (w(x,h) corresponds to the best approximation of y(x)):
  //!
  //!   &nbsp;&nbsp; w_0   = y0
  //!
  //!   &nbsp;&nbsp; for i = 0,1,2 ... (i-th step) <br>
  //!   &nbsp;&nbsp; w_i+1 = w_i + h*Phi(x_i,w_i;h;f) <br>
  //!   &nbsp;&nbsp; x_i+1 = x_i + h, x_0 = x0
  //!
  //!   &nbsp;&nbsp; Euler method gives Phi(x,y;h;f) = f(x,y)
  //!
  //!  Stepsize calculation (H denotes required stepsize, h current stepsize):
  //!
  //!   &nbsp;&nbsp; h/H = (16/15*(|w_i+1(x_i+h;h) - w_i+1(x_i+h;h/2)|/eps))^0.20 <br>
  //!   &nbsp;&nbsp; if h/H >> 2 -> h = 2*H -> required precision not achieved -> stepsize recalculation <br>
  //!   &nbsp;&nbsp; if h/H ~  2 -> h = 2*H -> required precision achieved -> use h for another iteration
  //!
  //! For more details see J.Stoer, R.Bulirsch: Introduction to Numerical Analysis,
  //! chap. 7, 3rd ed., and <a href=http://www.nr.com/oldverswitcher.html>www.nr.com</a>
  //!
  double RKCSolve(double x0, double y0, double x, double hTry);

  //! This method returns a solution y(x) to the so-called ODE initial value problem
  //!(IVP), using 5th(4th) order Runge-Kutta-Fehlberg method with adaptive stepsize
  //! control. The stepsize control algorithm is performed as follows: first two
  //! approximations to the exact solution w(x_i+1,h) (O(h6)) and v(x_i+1,h) (O(h5))
  //! are calculated and then the difference proportional to h*C(x) is expressed
  //! and used for calculation of a new stepsize, i.e. used as an estimate of the
  //! error. Due to a fact that the required precision _eps represents an accumulated
  //! error, permitted in maximum, C(x) is also dependent on h and thus expressed
  //! in more complicated way (Two coefficients: 0.25 and 0.20 are used), see
  //! Stepsize calculation. The IVP problem is defined as:
  //!
  //!   &nbsp;&nbsp; y' = dy/dx = f(x,y) with initial condition f(x0) = y0
  //!
  //!  Calculation (w(x,h), resp. v(x,h), correspond to the best approximation of y(x)):
  //!
  //!   &nbsp;&nbsp; w_0   = y0 <br>
  //!   &nbsp;&nbsp; v_0   = y0
  //!
  //!   &nbsp;&nbsp; for i = 0,1,2 ... (i-th step) <br>
  //!   &nbsp;&nbsp; w_i+1 = w_i + h*PhiI(x_i,w_i;h;f) <br>
  //!   &nbsp;&nbsp; v_i+1 = v_i + h*PhiII(x_i,w_i;h;f) <br>
  //!   &nbsp;&nbsp; x_i+1 = x_i + h, x_0 = x0
  //!
  //!  Stepsize calculation (H denotes required stepsize, h current stepsize):
  //!
  //!   &nbsp;&nbsp; fraction = |eps*h/(x - x0)/(h*PhiI - h*PhiII)| <br>
  //!   &nbsp;&nbsp; if fraction <= 1. -> reduce stepsize   h = 0.9*h*fraction^0.25 <br>
  //!   &nbsp;&nbsp; if fraction >  1. -> increase stepsize h = 0.9*h*fraction^0.2
  //!
  //! For more details see J.Stoer, R.Bulirsch: Introduction to Numerical Analysis,
  //! chap. 7, 3rd ed., and <a href=http://www.nr.com/oldverswitcher.html>www.nr.com</a>
  //!
  double RKFSolve(double x0, double y0, double x, double hTry);

protected:
  //! This method based on 4th order Runge-Kutta algorithm returns an approximation
  //! w(x_i+1,h) to the exact solution y(x_i+1) of the so-called initial value
  //! problem (IVP) at position x_i+1 = x_i + h and using given "integration" stepsize
  //! h, where the IVP is defined as y' = f(x,y), y(x0) = y0 and for equidistant
  //! stepsizes h the coordinate x_i can be expressed as x_i = x0 + ih.
  //!
  //!  General numerical algorithm for approximate values' calculation:
  //!
  //!   &nbsp;&nbsp; w_0   = y0
  //!   &nbsp;&nbsp; w_i+1 = w_i + h*Phi(x_i,w_i;h;f)
  //!   &nbsp;&nbsp; x_i+1 = x_i + h, x_0 = x0
  //!
  //!   &nbsp;&nbsp; Euler method gives Phi(x,y;h;f) = f(x,y)
  //!
  //!  Calculation using Runge-Kutta algorithm:
  //!
  //!   &nbsp;&nbsp; k1 = f(x_i, w_i)
  //!   &nbsp;&nbsp; k2 = f(x_i + h/2, w_i + h/2*k1)
  //!   &nbsp;&nbsp; k3 = f(x_i + h/2, w_i + h/2*k2)
  //!   &nbsp;&nbsp; k4 = f(x_i + h, w_i + h*k3)
  //!
  //!   &nbsp;&nbsp; Phi(x_i,w_i;h;f) = 1/6*[k1 + 2*k2 + 2*k3 + k4]
  //!
  //! Parameters: position x_i, approximation w_i, stepsize h and evaluated righ-hand
  //! side of the ODE dydx_i (It may be called more than once, thus given as param!)
  //!
  //! For more details see J.Stoer, R.Bulirsch: Introduction to Numerical Analysis,
  //! chap. 7, 3rd ed., and <a href=http://www.nr.com/oldverswitcher.html>www.nr.com</a>
  //!
  double RK4Method(double x_i, double w_i, double h, double dydx_i);

  //! This method based on 5th(4th) order Runge-Kutta-Fehlberg algorithm returns
  //! two approximations: w(x_i+1,h) (O(h6)) and v(x_i+1,h) (O(h5)), that approximate
  //! the exact solution y(x_i+1) of given initial value problem (IVP) at position
  //! x_i+1 = x_i + h and using an integration stepsize h. These approximative solutions
  //! are effectively used, when controlling required precision. The difference is
  //! proportional to h*C(x) and hence can be used to estimate an error and find
  //! a more suitable stepsize H. Due to a fact, that different schemes can be used to
  //! calculate w or v, we've chosen the one of Cash and Karp. The IVP is defined as
  //! y' = f(x,y), y(x0) = y0.
  //!
  //!  Rung-Kutta-Fehlberg algorithm for approximate values' calculation:
  //!
  //!   &nbsp;&nbsp; w_0   = y0 <br>
  //!   &nbsp;&nbsp; v_0   = y0 <br>
  //!   &nbsp;&nbsp; w_i+1 = w_i + h*PhiI(x_i,w_i;h;f) <br>
  //!   &nbsp;&nbsp; v_i+1 = v_i + h*PhiII(x_i,w_i;h;f) <br>
  //!   &nbsp;&nbsp; x_i+1 = x_i + h, x_0 = x0
  //!
  //!   &nbsp;&nbsp; PhiI(x_i,w_i;h;f)  = Sum(k=0, 5){c_kI  * f_k(x_i, w_i; h)} <br>
  //!   &nbsp;&nbsp; PhiII(x_i,w_i;h;f) = Sum(k=0, 5){c_kII * f_k(x_i, w_i; h)} <br>
  //!   &nbsp;&nbsp; f_k(x_i,w_i;h) = f(x_i + a_k*h, w_i + h * Sum(l=0, k-1){b_kl * f_l)
  //!
  //!   &nbsp;&nbsp; a_k, b_kl, c_k ...... Cash-Karp parameters
  //!
  //!   <table style="text-align: center">
  //!   <tr><td> k </td><td> a_k  </td><td>    b_k0    </td><td>   b_k1  </td><td>   b_k2    </td><td>     b_k3
  //!   </td><td>   b_k4   </td><td>   c_kI   </td><td>    c_kII    </td></tr> <tr><td> 0 </td><td>  0   </td><td>
  //!   </td><td>         </td><td>           </td><td>              </td><td>          </td><td>  37/378  </td><td>
  //!   2825/27648 </td></tr> <tr><td> 1 </td><td> 1/5  </td><td>    1/5     </td><td>         </td><td> </td><td>
  //!   </td><td>          </td><td>    0     </td><td>      0      </td></tr> <tr><td> 2 </td><td> 3/10 </td><td> 3/40
  //!   </td><td>   9/40  </td><td>           </td><td>              </td><td>          </td><td> 250/621  </td><td>
  //!   18575/48384 </td></tr> <tr><td> 3 </td><td> 3/5  </td><td>    3/10    </td><td>  -9/10  </td><td>   6/5
  //!   </td><td>              </td><td>          </td><td> 125/594  </td><td> 13525/55296 </td></tr> <tr><td> 4
  //!   </td><td>  1   </td><td>  -11/54    </td><td>   5/2   </td><td> -70/27    </td><td>    35/27     </td><td>
  //!   </td><td>    0     </td><td>   277/14336 </td></tr> <tr><td> 5 </td><td> 7/8  </td><td> 1631/55296 </td><td>
  //!   175/512 </td><td> 575/13824 </td><td> 44275/110592 </td><td> 253/4096 </td><td> 512/1771 </td><td>     1/4
  //!   </td></tr>
  //!   </table>
  //!
  //! Parameters: position x_i, approximation w_i, stepsize h, evaluated righ-hand
  //! side of the ODE dydx_i and estimated error between 5th and 4th order (i.e.
  //! the difference between w_i+1 and v_i+1, defined as wDiff).
  //!
  //! For more details see J.Stoer, R.Bulirsch: Introduction to Numerical Analysis,
  //! chap. 7, 3rd ed., and <a href=http://www.nr.com/oldverswitcher.html>www.nr.com</a>
  //!
  double RKFMethod(double x_i, double w_i, double h, double dydx_i, double* wDiffPointer);

  double (T::*_fce)(double, double); //!< Relative pointer to an integrated function

  T* _pTemplate; //!< Pointer to the template class

  float _eps; //!< Absolute precision (to hold "global" accumulation of errors constant)

}; // Class

//
// 4th order Runge-Kutta solver without adaptive stepsize control
//
template <class T>
double RunKutODESolver<T>::RKSolve(double x0, double y0, double x) {
  // Approximation to the required y(x_i) for given x_i
  double x_i = x0;
  double w_i = y0;

  // Set stepsize and number of iterations
  int nIter = (int)(1 / _eps);
  double h = (x - x0) / nIter;

  // Perform n steps
  for (int i = 0; i <= nIter; i++) { // nSteps; i++) {

    // Calculate new approximation w_i+1 (denoted again as w_i)
    w_i = RK4Method(x_i, w_i, h, FCE2(x_i, w_i));

    // Verify that step is not beyond machine precision
    if ((double)(x_i + h) == x_i) {
      //       std::cout << "Error - RunKutSolver::RKSolve: Stepsize is too small!!!"
      streamlog_out(ERROR) << "RunKutSolver::RKSolve: Stepsize is too small!!!" << std::endl;
      exit(1);
    } else
      x_i = x_i + h;
  }

  // std::cout << "Number of iterations: " << nIter+1 << std::endl;
  streamlog_out(MESSAGE0) << "   ODESolver - number of iterations: " << nIter + 1 << std::endl;

  // Result
  return w_i;
}

//
// 4th order Runge-Kutta solver with adaptive stepsize control
//
template <class T>
double RunKutODESolver<T>::RKCSolve(double x0, double y0, double x, double hTry) {
  // Approximation to the required y(x_i), resp. y(x_i+1), for given x_i, resp. x_i+1, and using stepsize hTry
  double x_i = x0;
  double w_i = y0;
  double w_i1 = 0.;

  // Approximation to the required y(x_i+1) for given x_i+1 and using 2 * half-stepsize hTry/2
  double w2_i1 = 0.;

  // Stepsize (h = actual stepsize, H = new proposed stepsize)
  double h = hTry;
  double h2 = h / 2.;
  double H = 0.;

  // Number of iterations
  int nIter = 0;

  // Allow MAXSTEPS steps in maximum
  for (int i = 0; i < MAXSTEPS; i++) {

    // Evaluate ODE right-hand side
    double dydx_i = FCE2(x_i, w_i);

    // Actual position + stepsize cannot be bigger than the final position x
    if ((x_i + h - x) * (x_i + h - x0) > 0.)
      h = x - x_i;

    // Adapt stepsize and calculate new approximation w_i+1
    bool adapt = true;
    while (adapt) {

      nIter++;

      // Calculate w_i+1(x_i + h; h) (denoted as w_i1) using stepsize h
      w_i1 = RK4Method(x_i, w_i, h, dydx_i);

      // Calculate w_i+1(x_i + h; h/2) (denoted as w2_i1) using 2 * stepsize h/2
      w2_i1 = RK4Method(x_i, w_i, h2, dydx_i);
      w2_i1 = RK4Method(x_i + h2, w2_i1, h2, FCE2(x_i + h2, w2_i1));

      // Verify if w_i+1(x_i + h; h) - w_i+1(x_i + h; h/2) is within required
      // precision
      double fraction = fabs(_eps / (w_i1 - w2_i1));
      H = h * pow(15. / 16. * fraction, 0.2);

      // If not, try again
      if ((h / H) > 2.5) {
        // Reduce stepsize, no more than 10 times
        h = (((2 * H) > (0.1 * h)) ? (2 * H) : (0.1 * h));
        h2 = h / 2.;

        adapt = true;
      } else {
        // New step
        hTry = 2 * H;

        adapt = false;
      }
    } // While

    // Verify that the step was not beyond machine precision
    if ((double)(x_i + h) == x_i) {
      //       std::cout << "Error - RunKutSolver::RKSolve: Stepsize is too small!!!"
      streamlog_out(ERROR) << "RunKutSolver::RKCSolve: Stepsize is too small!!!" << std::endl;
      exit(1);
    }
    // Calculate x_i+1 position
    x_i += h;

    // Set new w_i+1 position
    w_i = w2_i1;

    // Final position x achived?
    if ((x_i - x) * (x - x0) >= 0.) {

      //       std::cout << "Number of iterations: " << nIter+1 << std::endl;
      streamlog_out(MESSAGE0) << "   ODESolver - number of iterations: " << nIter + 1 << std::endl;
      return w2_i1;
    }

    // If not, set new stepsize
    h = hTry;
    h2 = hTry / 2.;

  } // For

  // Too many steps peformed
  // std::cout << "Error - RunKutSolver::RKSolve: Stepsize is too small!!!"
  streamlog_out(ERROR) << "RunKutSolver::RKCSolve: Too many steps performed!!!" << std::endl;
  exit(1);
}

//
// 5th(4th) order Runge-Kutta-Fehlberg solver with adaptive stepsize control
//
template <class T>
double RunKutODESolver<T>::RKFSolve(double x0, double y0, double x, double hTry) {
  // Approximation to the required y(x_i), resp. y(x_i+1), for given x_i, resp. x_i+1
  double x_i = x0;
  double w_i = y0;
  double w_i1 = 0;

  // Estimated precision of current step
  double wDiff = 0.;

  // Stepsize (h = actual stepsize, H = new stepsize)
  double h = hTry;
  double H = 0.;

  // Number of iterations
  int nIter = 0;

  // Allow MAXSTEPS steps in maximum
  for (int i = 0; i < MAXSTEPS; i++) {

    // Evaluate ODE right-hand side
    double dydx_i = FCE2(x_i, w_i);

    // Actual position + stepsize cannot be bigger than the final position x
    if ((x_i + h - x) * (x_i + h - x0) > 0.)
      h = x - x_i;

    // Adapt stepsize and calculate new approximation w_i+1
    bool adapt = true;
    while (adapt) {

      nIter++;

      // Calculate w_i+1(x_i + h; h) (denoted as w_i) using stepsize h and get
      // estimated precision of obtained result
      w_i1 = RKFMethod(x_i, w_i, h, dydx_i, &wDiff);

      // Verify if wDiff is within required precision
      double fraction = fabs((_eps * h / (x - x0)) / wDiff);

      // If not, try again
      if (fraction <= 1.) {

        // Reduce stepsize
        H = 0.9 * h * pow(fraction, 0.25);

        // No more than 10 times
        h = ((H > (0.1 * h)) ? H : (0.1 * h));

        adapt = true;
      }
      // Step suceeded
      else {

        // Suggest new step
        H = 0.9 * h * pow(fraction, 0.2);
        hTry = H;

        adapt = false;
      }
    } // While

    // Verify that the step was not beyond machine precision
    if ((double)(x_i + h) == x_i) {
      //       std::cout << "Error - RunKutSolver::RKSolve: Stepsize is too small!!!"
      streamlog_out(ERROR) << "RunKutSolver::RKFSolve: Stepsize is too small!!!" << std::endl;
      exit(1);
    }
    // Calculate x_i+1 position
    x_i += h;

    // Set new w_i+1 positio
    w_i = w_i1;

    // Final position x achived?
    if ((x_i - x) * (x - x0) >= 0.) {

      //       std::cout << "Number of iterations: " << nIter+1 << std::endl;
      streamlog_out(MESSAGE0) << "   ODESolver - number of iterations: " << nIter + 1 << std::endl;
      return w_i1;
    }

    // If not, set new stepsize
    h = hTry;

  } // For

  // Too many steps peformed
  // std::cout << "Error - RunKutSolver::RKSolve: Stepsize is too small!!!"
  streamlog_out(ERROR) << "RunKutSolver::RKFSolve: Too many steps performed!!!" << std::endl;
  exit(1);
}

//
// 4th order Runge-Kutta method - calculates i+1-th step, when solving 1st order ODE
//
template <class T>
double RunKutODESolver<T>::RK4Method(double x_i, double w_i, double h, double dydx_i) {
  double h2 = h / 2.;
  double k1, k2, k3, k4;
  double Phi;

  // First step
  k1 = dydx_i;

  // Second step
  k2 = FCE2(x_i + h2, w_i + h2 * k1);

  // Third step
  k3 = FCE2(x_i + h2, w_i + h2 * k2);

  // Fourth step
  k4 = FCE2(x_i + h, w_i + h * k3);

  // Phi(x,y;h)
  Phi = 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4);

  // Result w(x_i+1,h) = w_i+1
  return (w_i + h * Phi);
}

//
// 5th(4th) order Runge-Kutta-Fehlberg method - calculates i+1-th step, when solving 1st order ODE
//
template <class T>
double RunKutODESolver<T>::RKFMethod(double x_i, double w_i, double h, double dydx_i, double* wDiffPointer) {
  static float a1 = 0.2, a2 = 0.3, a3 = 0.6, a4 = 1., a5 = 0.875;
  static float b10 = 0.2, b20 = 0.075, b30 = 0.3, b40 = -11. / 54., b50 = 1631. / 55296., b21 = 0.225, b31 = -0.9,
               b41 = 2.5, b51 = 175. / 512., b32 = 1.2, b42 = -70. / 27., b52 = 575. / 13824., b43 = 35. / 27.,
               b53 = 44275. / 110592., b54 = 253. / 4096.;

  static float c0I = 37. / 378., c2I = 250. / 621., c3I = 125. / 594., c4I = 0., c5I = 512. / 1771.,
               c0II = 2825. / 27648., c2II = 18575. / 48384., c3II = 13525. / 55296., c4II = 277. / 14336., c5II = 0.25,
               difc0 = c0I - c0II, difc2 = c2I - c2II, difc3 = c3I - c3II, difc4 = c4I - c4II, difc5 = c5I - c5II;

  double f0, f1, f2, f3, f4, f5;
  double PhiI;

  // First step
  f0 = dydx_i;

  // Second step
  f1 = FCE2(x_i + a1 * h, w_i + h * (b10 * f0));

  // Third step
  f2 = FCE2(x_i + a2 * h, w_i + h * (b20 * f0 + b21 * f1));

  // Fourth step
  f3 = FCE2(x_i + a3 * h, w_i + h * (b30 * f0 + b31 * f1 + b32 * f2));

  // Fifth step
  f4 = FCE2(x_i + a4 * h, w_i + h * (b40 * f0 + b41 * f1 + b42 * f2 + b43 * f3));

  // Sixth step
  f5 = FCE2(x_i + a5 * h, w_i + h * (b50 * f0 + b51 * f1 + b52 * f2 + b53 * f3 + b54 * f4));

  // Phi(x,y;h) (c1I = 0 and c4I = 0)
  PhiI = c0I * f0 + c2I * f2 + c3I * f3 + c5I * f5;

  // Estimated precision wDiff = w(x_i+1,h) - v(x_i+1,h)
  *wDiffPointer = h * (difc0 * f0 + difc2 * f2 + difc3 * f3 + difc4 * f4 + difc5 * f5);

  // Result w(x_i+1,h) = w_i+1
  return (w_i + h * PhiI);
}

} // namespace sistrip

#endif // RUNKUTODESOLVER_H
