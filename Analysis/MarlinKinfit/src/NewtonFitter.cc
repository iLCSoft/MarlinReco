////////////////////////////////////////////////////////////////
// Class NewtonFitter
//
// Author: Benno List
// Last update: $Date: 2008-02-25 08:23:56 $
//          by: $Author: gaede $
// 
// Description: kinematic fit using Newton method
//               
////////////////////////////////////////////////////////////////

#undef NDEBUG

#include "NewtonFitter.h" 

#include<iostream>
#include<cmath>
#include<cassert>

#include "BaseFitObject.h"
#include "BaseHardConstraint.h"
#include "BaseSoftConstraint.h"
#include "ftypes.h"
#include "cernlib.h"

using std::cout;
using std::cerr;
using std::endl;

// constructor
NewtonFitter::NewtonFitter() 
: npar (0), ncon (0), nsoft (0) 
{}

// destructor
NewtonFitter::~NewtonFitter() {}

double NewtonFitter::fit() {

  // cout statements
  int debug = 0;
  int nitdebug = -1;

  // order parameters etc
  initialize();
  
  // initialize eta, etasv, y 
  int idim = npar+ncon;
  FDouble *x  = new FDouble [idim];
  FDouble *dx = new FDouble [idim];
  FDouble *y  = new FDouble [idim];
  FDouble *M  = new FDouble [idim*idim];
  FInteger *ir = new FInteger[idim];
  double *perr = new double[idim];
  
  for (int i = 0; i < idim; ++i) x[i] = y[i] = 0;
  
  // Get starting values
  for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
    BaseFitObject *fo = *i;
    assert (fo);
    for (int ilocal = 0; ilocal < fo->getNPar(); ++ilocal) {
      int iglobal = fo->getGlobalParNum (ilocal);
      assert (iglobal >= 0 && iglobal < npar);
      x[iglobal] = fo->getParam (ilocal);
    }
  }
  

//   // Set initial lambda values all to 1
//   for (unsigned int icon = 0; icon < constraints.size(); ++icon) {
//     BaseHardConstraint *c = constraints[icon];
//     assert (c);
//     int kglobal = c->getGlobalNum();
//     if (debug) cout << "Constraint " << kglobal << ", idim=" << idim << endl;
//     assert (kglobal >= 0 && kglobal < idim);
//     x[kglobal] = 1;
//   }

  // For step size limitation: get errors of all parameters
  for (int i = 0; i<idim; ++i) perr[i] = 0;
  for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
    BaseFitObject *fo = *i;
    assert (fo);
    for (int ilocal = 0; ilocal < fo->getNPar(); ++ilocal) {
      int iglobal = fo->getGlobalParNum (ilocal);
      assert (iglobal >= 0 && iglobal < npar);
      perr[iglobal] = fo->getError (ilocal);
    }
  }
  

   
  
  // LET THE GAMES BEGIN
  
  bool converged = 0;
  ierr = 0;
  
  double chi2new = calcChi2();
  nit = 0;
  if (debug>1) {
    cout << "Fit objects:\n";
    for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
      BaseFitObject *fo = *i;
      assert (fo);
      cout << fo->getName() << ": " << *fo << ", chi2=" << fo->getChi2() << endl;
    }
    cout << "constraints:\n";
    for (ConstraintIterator i = constraints.begin(); i != constraints.end(); ++i) {
      BaseHardConstraint *c = *i;
      assert (c);
      cout << i-constraints.begin() << " " << c->getName() << ": " << c->getValue() << "+-" << c->getError() << endl;
    }
    cout << "soft constraints:\n";
    for (SoftConstraintIterator i = softconstraints.begin(); i != softconstraints.end(); ++i) {
      BaseSoftConstraint *c = *i;
      assert (c);
      cout << i-softconstraints.begin() << " " << c->getName() << ": " << c->getValue() << "+-" << c->getError() 
           << ", chi2: " << c->getChi2() << endl;
    }
  }
  
  do {
  
    double chi2old = chi2new;
    
    if (debug>1 && (nit==0 || nit>nitdebug)) cout << "===================\nStarting iteration " << nit << endl;
    if (debug>2 && (nit==0 || nit>nitdebug)) cout << "Fit objects:\n";
    for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
      BaseFitObject *fo = *i;
      assert (fo);
      if (debug>2 && (nit==0 || nit>nitdebug))  cout << fo->getName() << ": " << *fo << ", chi2=" << fo->getChi2() << endl;
    }
    if (debug>2 && (nit==0 || nit>nitdebug)) cout << "constraints:\n";
    for (ConstraintIterator i = constraints.begin(); i != constraints.end(); ++i) {
      BaseHardConstraint *c = *i;
      assert (c);
      if (debug>2 && (nit==0 || nit>nitdebug)) cout << i-constraints.begin() << " " << c->getName() << ": " << c->getValue() << "+-" << c->getError() << endl;
    }
    if (debug>2 && (nit==0 || nit>nitdebug)) cout << "soft constraints:\n";
    for (SoftConstraintIterator i = softconstraints.begin(); i != softconstraints.end(); ++i) {
      BaseSoftConstraint *c = *i;
      assert (c);
      if (debug>2 && (nit==0 || nit>nitdebug)) cout << i-softconstraints.begin() << " " << c->getName() << ": " << c->getValue() << "+-" << c->getError() 
           << ", chi2: " << c->getChi2() << endl;
    }
 
    // Set M to 0
    for (int i = 0; i < idim*idim; ++i) M[i] = 0;
    // Reset y
    for (int i = 0; i < idim; ++i) y[i] = 0;
    // Compose M: 
    // First, all terms d^2 chi^2/dx1 dx2
    for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
      BaseFitObject *fo = *i;
      assert (fo);
      fo->addToGlobalChi2DerMatrix (M, idim);
    }
    
//     cout << "After adding covariances from fit ojects:\n";
//     printMy (M, y, idim);
    // Second, all terms d^2 chi^2/dlambda dx, 
    // i.e. the first derivatives of the contraints,
    // plus the second derivatives times the lambda values
    for (unsigned int k = 0; k < constraints.size(); ++k) {
      BaseHardConstraint *c = constraints[k];
      assert (c);
      int kglobal = c->getGlobalNum();
      assert (kglobal >= 0 && kglobal < idim);
      c->add1stDerivativesToMatrix (M, idim);
      c->add2ndDerivativesToMatrix (M, idim, x[kglobal]);
    }
//     cout << "After adding derivatives of constraints::\n";
//     printMy (M, y, idim);
    
    // Now, calculate the result vector y with the values of the derivatives
    // d chi^2/d x
    // First, for the parameters
    for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
      BaseFitObject *fo = *i;
      assert (fo);
      fo->addToGlobalChi2DerVector (y, idim);
    }
//     cout << "After adding fo derivatives to y::\n";
//     printMy (M, y, idim);
    
    // Now add lambda*derivatives of constraints,
    // And finally, the derivatives w.r.t. to the constraints, i.e. the constraints themselves
    for (unsigned int k = 0; k < constraints.size(); ++k) {
      BaseHardConstraint *c = constraints[k];
      assert (c);
      int kglobal = c->getGlobalNum();
      assert (kglobal >= 0 && kglobal < idim);
      c->addToGlobalChi2DerVector (y, idim, x[kglobal]);
      y[kglobal] = c->getValue();
    }
  
    // Finally, treat the soft constraints

    for (SoftConstraintIterator i = softconstraints.begin(); i != softconstraints.end(); ++i) {
      BaseSoftConstraint *bsc = *i;
      assert (bsc);
      bsc->add2ndDerivativesToMatrix (M, idim);
      bsc->addToGlobalChi2DerVector (y, idim);
    }
  
    // from x_(n+1) = x_n - y/y' = x_n - M^(-1)*y we have M*(x_n-x_(n+1)) = y, 
    // which we solve for dx = x_n-x_(n+1) and hence x_(n+1) = x_n-dx
  
    for (int i = 0; i < idim; ++i) dx[i] = y[i];
    
//     cout << "Complete system:\n";
//     printMy(M, y, idim);
    
    FInteger ifail = 0;
    // dseqn  (idim, M, idim, ifail, 1, dx);
//     double Morig[idim*idim];
//     for (int i = 0; i < idim*idim; ++i) Morig[i]=M[i];
    deqn  (idim, M, idim, ir, ifail, 1, dx);
    if (ifail != 0) cout << "NewtonFitter::fit: ifail=" << ifail << endl;
    
//     cout << "Result of dseqn: " << ifail << endl;
//     cout << "dx = (";
//    for (int i = 0; i < idim; ++i) cout << dx[i] << (i+1<idim?", ":")\n");
    
//     double ycalc[idim];
//     for (int i = 0; i < idim; ++i) {
//       ycalc[i] = 0;
//       for (int j = 0; j < idim; ++j) {
//         ycalc[i] += Morig[idim*i+j]*dx[j];
//       }
//       cout << "i=" << i << ", y=" << y[i] << " - ycalc=" << ycalc[i] << " = diff: " << y[i]-ycalc[i] << endl;
//     }
    
    
  
    // Update vector x to new values
//      cout << "y = (";
//      for (int i = 0; i < idim; ++i) cout << y[i] << (i+1<idim?", ":")\n");
    // Update vector x to new values
//      cout << "old x = (";
//      for (int i = 0; i < idim; ++i) cout << x[i] << (i+1<idim?", ":")\n");

    // Step size limitation: Calculate max number of sigmas
    double maxstep = 4.0;
    double maxsigma = 0;
    for (int i = 0; i < idim; ++i) {
      if(perr[i]>0 && std::abs(dx[i])>maxsigma*perr[i]) 
        maxsigma=std::abs(dx[i])/perr[i];
        if (debug && nit>nitdebug && perr[i]>0 && std::abs(dx[i]) > maxstep*perr[i]) 
          cout << "step size for parameter " << i << ": " << dx[i]
               << " / " << perr[i] << " = " << dx[i]/perr[i] << endl;
    }
    double scale = 1.0;    
    if (maxsigma > maxstep) {
      if (debug && nit>nitdebug) cout << "Step size limitation: maxsigma=" << maxsigma << endl;
      scale = maxstep/maxsigma;
    }    
    
    for (int i = 0; i < idim; ++i) x[i] -= scale*dx[i];
    
//      cout << "new x = (";
//      for (int i = 0; i < idim; ++i) cout << x[i] << (i+1<idim?", ":")\n");
  
    // Update values in Fitobjects
    bool significant = false;
    for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
      BaseFitObject *fo = *i;
      assert (fo);
      bool s = fo->updateParams (x, idim);
      significant |=  s;
      if (debug && nit>nitdebug && s) 
        cout << "Significant update for FO " << i-fitobjects.begin() << " (" 
             << fo->getName() << ")\n";
    }
  
//   *-- Convergence criteria 

    chi2new = calcChi2();
    if (debug && nit>nitdebug) cout << "old chi2: " << chi2old << ", new chi2: " << chi2new << ", diff=" << chi2old-chi2new << endl;
    converged = !significant;
    bool constraintsok = true;
    if (true || converged || debug>1) {
      for (ConstraintIterator ci = constraints.begin(); ci != constraints.end(); ++ci) {
        BaseHardConstraint *c = *ci;
        assert (c);
        constraintsok &= (std::abs (c->getValue()) <= 1E-2*c->getError());
        if (debug && nit>nitdebug && std::abs (c->getValue()) >= 1E-2*c->getError())
          cout << "Constraint " << ci-constraints.begin() << " " << c->getName()
               << ": not fulfilled: " << c->getValue() << "+-" << c->getError() << endl;
          if (debug > 2 && nit>nitdebug) cout << "Constraint " << ci-constraints.begin() << ": " << c->getValue()
                              << "+-" << c->getError() << endl;
      }
    }
    converged = !significant && constraintsok;
    if (debug && nit > nitdebug) {
      cout << "nit=" << nit << ", significant: " << significant << ", constraintsok=" <<  constraintsok << endl;
    }
    if (debug && nit > nitdebug) {
      cout << "Fit probability: " << prob(chi2new,ncon+nsoft-nunm) << endl;
    }
    ++nit;
    if (nit > 200) ierr = 1;
  
  } while (!(converged || ierr));
  
// *-- End of iterations - calculate errors.

//    ====> HERE GOES ERROR CALCULATION <======

  if (debug>1) {
    cout << "========= END =========\n";
    cout << "Fit objects:\n";
    for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
      BaseFitObject *fo = *i;
      assert (fo);
      cout << fo->getName() << ": " << *fo << ", chi2=" << fo->getChi2() << endl;
    }
    cout << "constraints:\n";
    for (ConstraintIterator i = constraints.begin(); i != constraints.end(); ++i) {
      BaseHardConstraint *c = *i;
      assert (c);
      cout << i-constraints.begin() << " " << c->getName() << ": " << c->getValue() << "+-" << c->getError() << endl;
    }
    cout << "=============================================\n";
  }

  delete[] x;
  delete[] dx;
  delete[] y;
  delete[] M;
  delete[] perr;
  delete[] ir;
  

// *-- Turn chisq into probability.
  FReal chi = FReal(chi2new);

  return fitprob = prob(chi ,ncon+nsoft-nunm);
    
}

bool NewtonFitter::initialize() {
//  bool debug = true;

  // tell fitobjects the global ordering of their parameters:
  int iglobal = 0;
  // measured parameters first!
  for (unsigned int ifitobj = 0; ifitobj < fitobjects.size(); ++ifitobj) {
    for (int ilocal = 0; ilocal < fitobjects[ifitobj]->getNPar(); ++ilocal) {
      fitobjects[ifitobj]->setGlobalParNum (ilocal, iglobal);
      ++iglobal;
    }
  }
  npar = iglobal;
  nunm = 0;
  // now unmeasured parameters!
  for (unsigned int ifitobj = 0; ifitobj < fitobjects.size(); ++ifitobj) {
    for (int ilocal = 0; ilocal < fitobjects[ifitobj]->getNPar(); ++ilocal) {
      if (!fitobjects[ifitobj]->isParamMeasured(ilocal)) ++nunm;
    }
  }
  
  // set number of constraints
  ncon = constraints.size();
  // Tell the constraints their numbers
  for (unsigned int icon = 0; icon < constraints.size(); ++icon) {
    BaseHardConstraint *c = constraints[icon];
    assert (c);
    c->setGlobalNum (npar+icon);
//    if (debug) cout << "Constraint " << icon << " -> global " << c->getGlobalNum() << endl;
  }
  
  nsoft = softconstraints.size();
  
  if (nunm > ncon+nsoft) {
    cerr << "NewtonFitter::initialize: nunm=" << nunm << " > ncon+nsoft=" 
         << ncon << "+" << nsoft << endl;
  }
  
  return true;

}
  
double NewtonFitter::calcChi2() {
  chi2 = 0;
  for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
    BaseFitObject *fo = *i;
    assert (fo);
    chi2 += fo->getChi2();
  }
  for (SoftConstraintIterator i = softconstraints.begin(); i != softconstraints.end(); ++i) {
    BaseSoftConstraint *bsc = *i;
    assert (bsc);
    chi2 += bsc->getChi2();
  }
  return chi2;
}

void NewtonFitter::printMy (double M[], double y[], int idim) {
  for (int i = 0; i < idim; ++i) {
    cout << i << "  [ " << M[idim*i+0];
      for (int j = 1; j<idim; ++j) cout << ", " << M[idim*i+j];
      cout << "]  [" << y[i] << "]\n";
  }
}

int NewtonFitter::getError() const {return ierr;}
double NewtonFitter::getProbability() const {return fitprob;}
double NewtonFitter::getChi2() const {return chi2;}
int NewtonFitter::getDoF() const {return ncon+nsoft-nunm;}
int NewtonFitter::getIterations() const {return nit;}
