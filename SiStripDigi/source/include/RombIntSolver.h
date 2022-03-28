#ifndef ROMBINTSOLVER_H
#define ROMBINTSOLVER_H 1

#include <math.h>
#include "Colours.h"

// Include Marlin
#include <streamlog/streamlog.h>

/**
\addtogroup SiStripDigi SiStripDigi
@{
*/

namespace sistrip {

#define EPS_R       1E-4
#define NSTEPSMAX_R 20
#define NSTEPSMIN_R 5
#define FCE(x)      (_pTemplate->*_fce)(x)

//! Double precision Romberg integration solver (Template class represents
//! class whose method is to be integrate!). The template class doesn't have
//! *.cpp file!
//!
//! @see Integrate
//!
//! @author Z. Drasal, Charles University in Prague
//!

template <class T> class RombIntSolver {

 public:

//!Constructor - sets function to be integrated
   RombIntSolver(double (T::* fce)(double), T* pTemplate) : _eps(EPS_R),
                 _nStepsMin(NSTEPSMIN_R), _nStepsMax(NSTEPSMAX_R),
                 _fce(fce), _pTemplate(pTemplate) {;}

//!Constructor - sets function to be integrated and integration precision
   RombIntSolver(double (T::* fce)(double), T* pTemplate, float eps) : _eps(eps),
                 _nStepsMin(NSTEPSMIN_R), _nStepsMax(NSTEPSMAX_R),
                 _fce(fce), _pTemplate(pTemplate) {;}

//!Destructor
   ~RombIntSolver() {;}

//!This method returns a definite integral of function fce calculated
//!in interval (a,b) using so-called Romberg algorithm. The algorithm
//!exploites a very general idea of Richardson's deffered approach to the
//!limit, where the function is integrated using the trapezium rule first
//!(with h subsequently set to (a-b), (a-b)/2, (a-b)/4) and then the result
//!is extrapolated to the limit. For n=0 (h=(a-b)) it corresponds to m=0;
//!for n=1 (h=(a-b)/2) it corresponds to m=0 and m=1; ..., where m represents
//!the m-th step to the limit (Richardson extrapolation).<br>
//!
//! Calculation (R(0,0); R(1,0)->R(1,1); R(2,0->R(2,1)->R(2,2)):<br>
//!
//!  &nbsp;&nbsp; R(0,0) = 1/2*(b-a)*(f(a)+f(b)<br>
//!  &nbsp;&nbsp; R(n,0) = 1/2*R(n-1,0) + h*Sum(i=1,2^(n-1)){f(a+(2*i-1)*h)}<br>
//!
//!  &nbsp;&nbsp; R(n,m) ... m-th iteration - achieved precision O(h^(2^(m+1)))<br>
//!  &nbsp;&nbsp; R(n,m) = R(n,m-1) + (R(n,m-1)-R(n-1,m-1))/(4^m-1)<br>
//!
//! Total relative precision (influenced by maximum number of steps permitted):<br>
//!
//!  &nbsp;&nbsp; epsilon = |R(n,n) - R(n,n-1)|/R(n,n) <br>
//!
//!For more details see J.Stoer, R.Bulirsch: Introduction to Numerical Analysis,
//!chap. 2 and 3, 3rd ed., and <a href=http://www.nr.com/oldverswitcher.html>www.nr.com</a>
//!
   double Integrate(double endPointA, double endPointB);

 protected:

//!This method returns a definite integral of function fce in interval (a,b)
//!and is based on classical extended trapezium rule.
//!
//! Calculation (n-th iteration):
//!
//!  &nbsp;&nbsp; Int(a,b){f(x)*dx} = h/2*[f(a)+f(b)] + h*Sum(i=1,n-1){f(a+i*h)} - O(h^2)
//!
//!  &nbsp;&nbsp; I(n) = 1/2*I(n-1) + h*Sum(i=1,2^(n-1)){f(a+(2*i-1)*h)}<br>
//!
//!Where n=0 corresponds to h=(a-b), n=1 to h=(a-b)/2 and the error of the
//!aproximation can be derived from Euler-MacLaurin summation formula and is
//!following:<br>
//!
//!  &nbsp;&nbsp; -B2*h^2/2!*(f'(b)-f'(a)) - B4*h^4/4!*(f'''(b)-f'''(a)) - ...
//!
//!For more details see J.Stoer, R.Bulirsch: Introduction to Numerical Analysis,
//!chap. 2 and 3, 3rd ed., and <a href=http://www.nr.com/oldverswitcher.html>www.nr.com</a>
//!
   double TrapezRuleInt(int n);

 private:

   double _a;//!<Left integration endpoint
   double _b;//!<Right integration endpoint

   float _eps;//!<Total integration precision (Caution: Don't overestimate this value!)

   int _nStepsMin;//!<Minimum number of integration steps permitted (Avoid spurious early convergence!)
   int _nStepsMax;//!<Maximum number of integration steps permitted

   double (T::* _fce)(double);//!<Relative pointer to an integrated function

   T * _pTemplate;//!< Pointer to the template class

}; // Class

//
// Romberg method calculating a definite integral in interval (a,b)
//
template <class T>
double RombIntSolver<T>::Integrate(double endPointA, double endPointB)
{
// Set endpoints
   _a = endPointA;
   _b = endPointB;

// Degree of refinement n, m, m4 = 4^m
   int i, n, m, m4;

// Integration sesult for given n, m
   double R[NSTEPSMAX_R][NSTEPSMAX_R];

// Calculate integral using Trapezium method
   for (n=0; n<_nStepsMax; n++) {
      R[n][0] = TrapezRuleInt(n);

//    std::cout << std::endl;
//    std::cout << "Step n=" << n << "R[" << n << "][0]" << R[n][0];

   // Use Richardson approach to increase precision
      for (m=1; m<=n; m++) {
         for (m4=1, i=0; i<(2*m); i++) m4<<=1;
         R[n][m] = R[n][m-1] + (R[n][m-1] - R[n-1][m-1])/(m4 - 1);

//       std::cout << std::endl;
//       std::cout << "Step n,m=" << n << "," << m << "R[" << n << "]"
//                 << "[" << m << "]" <<R[n][m];
      }

   // Find out precision
      if ((n>=_nStepsMin)&&((fabs(R[n][n] - R[n][n-1])/R[n][n])<_eps)) return R[n][n];
//    std::cout << std::endl;
   }

   // Too many steps to obtain required precision
   //std::cout << "Warning - Too many steps in Romberg integr., required "
   //          << "precision wasn't achieved" << std::endl;
   streamlog_out(WARNING) << "Warning - Too many steps in Romberg integr., required "
                          << "precision wasn't achieved" << std::endl;

   return R[n][n];
}

//
// Trapezium rule - calculates a definite integral in interval (a,b) - n-th iteration
//
template <class T>
double RombIntSolver<T>::TrapezRuleInt(int n)
{
   static double trapSum = 0;
   double x, h, h2, subSum;
   int i, iN;

// 1st iteration (n = 0)
   if (n == 0) {
      trapSum = 0.5 * (_b - _a) * (FCE(_a)+FCE(_b));

      return trapSum;
   }
// n-th iteration
   else {

   // Calculate iN = 2^n
      for (iN=1, i=0; i<n; i++) iN<<=1;

   // Calculate h
      h = (_b - _a)/iN;
      h2= 2*h;

   // Calculate subSum (Add to trapSum values that hasn't been calculated yet)
      x = _a + h;
      for (subSum=0, i=0; i<(iN/2); i++, x+=h2) subSum += FCE(x);

   // Calculate final trapSum
      trapSum = 0.5*trapSum + subSum*h;

      return trapSum;
   }
}

} // Namespace

/** @} */

#endif // ROMBINTSOLVER_H
