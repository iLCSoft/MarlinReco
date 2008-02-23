////////////////////////////////////////////////////////////////
// Class OPALFitter
//
// Author: Jenny Boehme
// Last update: $Date: 2008-02-23 11:18:39 $
//          by: $Author: listj $
// 
// Description: kinematic fit a la WWFGO
// DISCLAIMER: the only object oriented stuff in here is the
//             interface to fit objects (jets, neutrinos,....)
//             and constraints (px, py, M_W, ....)
//             which replaces the calls to WWKCNS.
//             The OPALFitter::fit() method is almost an
//             'F to C' translation of WWFGO. It is NOT
//             considered to be good C++, but was done
//             on purpose as a first implementation.
//             An OO version might follow later!     
//
// Changed 6.12.04 B. List: Added additional term in the definition
//             of S to make S nonsingular also in the case when
//             some constraint depends only on unmeasured quantities.
//               
////////////////////////////////////////////////////////////////

#include<iostream>
#include<cmath>
#include<cassert>

#include "OPALFitter.h" 

#include "BaseFitObject.h"
#include "BaseHardConstraint.h"
#include "ftypes.h"
#include "cernlib.h"

using std::cout;
using std::cerr;
using std::endl;
using std::abs;

static int debug = 0;

// constructor
OPALFitter::OPALFitter() : npar(0), nmea(0), nunm(0), ncon(0), ierr (0), nit (0)  {};

// destructor
OPALFitter::~OPALFitter() {std::cout << "destroying OPALFitter" << std::endl;};

// do it (~ transcription of WWFGO as of ww113)
double OPALFitter::fit() {

  // cout statements
  int inverr = 0;

  // order parameters etc
  initialize();
  
  // initialize eta, etasv, y 
  FDouble eta[NPARMAX];
  FDouble etasv[NPARMAX];
  FDouble y[NPARMAX];           
  
  for (unsigned int i = 0; i < fitobjects.size(); ++i) {
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ++ilocal) {
      int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
      if (iglobal >= 0) {
        assert (iglobal >= 0 && iglobal < NPARMAX);
        y[iglobal] = eta[iglobal] = fitobjects[i]->getParam(ilocal);
      }
      //if (debug) cout << "eta[" << iglobal << "] = " << eta[iglobal] 
      //                << " for jet " << i << " and ilocal = " << ilocal << endl;
    }
  }

  assert (ncon <= NCONMAX);
  /// initialize dfeta ( = d F / d eta)
  FDouble dfeta[NCONMAX][NPARMAX];   //  npar x nconstraint
  for (int k=0; k < ncon; k++) {
    for (int j=0; j < npar; j++) dfeta[k][j]=0;
    constraints[k]->getDerivatives(NPARMAX, dfeta[k]);
    if (debug>1) for (int j=0; j < npar; j++) 
      if (dfeta[k][j]!= 0) cout << "1: dfeta[" << k << "][" << j << "] = " << dfeta[k][j] << endl;
  }
  

  // R & S & W1
  FDouble  R[NCONMAX];                      // list of constraints
  FDouble  S[NCONMAX][NCONMAX];             // nconstraint  x  nconstraint
  FDouble  W1[NUNMMAX][NUNMMAX];            // nunmeasured  x  nunmeasured
  FDouble  V[NPARMAX][NPARMAX];             // error matrix of parameters
  FDouble  VINV[NPARMAX][NPARMAX];          // inverse error matrix of parameters
  FDouble  dxi[NUNMMAX];                    // shift of unmeasured parameters
  FDouble  alam[NCONMAX];                   // Lagrange multipliers

  // chi2's, step size, # iterations
  double chinew=0, chit=0, chik=0;
  double alph = 1.;
  nit = 0;
  // convergence criteria (as in WWINIT)
  // int nitmax = 20;
  int nitmax = 200;
  double chik0 = 100.;
  double chit0 = 100.;
  double dchikc = 1.0E-3;
  double dchitc = 1.0E-4;
  double dchikt = 1.0E-2;
  double dchik  = 1.05;
  double chimxw = 10000.; // 1000.;
  double almin = 0.05;
 

  // repeat with or with out smaller steps size 
  bool repeat = true;
  bool scut = false;
  bool calcerr = true;
 
  // start of iterations
  while (repeat) {
    bool updatesuccess = true;
  
//*-- If necessary, retry smaller step, same direction
    if (scut) {
      for (int j=0; j < npar; j++) eta[j] = etasv[j];
      updatesuccess = updateFitObjects (eta);
      if (!updatesuccess) {
        std::cerr << "OPALFitter::fit: old parameters are garbage!" << std::endl;
        return -1;
      }
      for (int k = 0; k < ncon; ++k)  {
        for (int j=0; j < npar; j++) dfeta[k][j]=0;
        constraints[k]->getDerivatives(NPARMAX, dfeta[k]);
        if (debug>1) for (int j=0; j < npar; j++) 
          if (dfeta[k][j]!= 0) cout << "2: dfeta[" << k << "][" << j << "] = " << dfeta[k][j] << endl;
      }
    } 
    else {
      for (int j=0; j < npar; j++) etasv[j] = eta[j];
      chik0 = chik;
      chit0 = chit;
    }
    
    // Get covariance matrix
    for (int i = 0; i < npar; ++i) 
      for (int j = 0; j < npar; ++j) V[i][j] = 0;
    for (unsigned int ifitobj = 0; ifitobj<fitobjects.size(); ++ifitobj) {
      fitobjects[ifitobj]->addToGlobCov (&V[0][0], NPARMAX);
    }  
    if (debug>1) 
      for (int  i = 0; i < npar; ++i) 
        for (int j = 0; j < npar; ++j)
          if (V[i][j] != 0)
            cout << "V[" << i << "][" << j << "]=" << V[i][j] << endl;;
    for (int i = 0; i < nmea; ++i) 
      for (int j = 0; j < nmea; ++j) VINV[i][j] = V[i][j];
      
    if (debug>1) 
      for (int i = 0; i < npar; ++i) 
        for (int j = 0; j < npar; ++j) 
          if (V[i][j] != 0) cout << "V[" << i << "][" << j << "] = " << V[i][j] << endl;
    
    // invert covariance matrix (needed for chi2 calculation later)
    dsinv(nmea, &VINV[0][0], NPARMAX, inverr);
    if (inverr != 0) {
      if (debug) cerr << "VINV: dsinv error " << inverr << endl;
      ierr = 6;
      calcerr = false;
      if (debug) 
        for (int i = 0; i < nmea; ++i) {
          if (V[i][i] == 0) cout << "  zero element! i = " << i << "\n";
          for (int j = i; j < nmea; ++j) {
            if (V[i][j] != 0) 
              cout << "...V[" << i << "][" << j << "] = " << V[i][j]
                   << ", r = " <<  V[i][j]/std::sqrt(V[i][i]*V[j][j])
                   << endl;
            if (std::abs(V[i][j]-V[j][i])>1E-6*std::abs(V[i][j]+V[j][i]))
              cout << "asymmetric: V[i][j]-V[j][i] = " << V[i][j]-V[j][i]
                   << ", V[i][j]+V[j][i]= " << V[i][j]+V[j][i] << endl;
          } 
        }
        
      // If dsinv fails, try dinv...  
      FInteger ir[NPARMAX];  
      dinv(nmea, &VINV[0][0], NPARMAX, ir, inverr);
      if (inverr != 0) {
        std::cerr << "VINV: dsinv and dinv failed!\n";
        break;
      }
      else {
        if (debug) std::cerr << "... but dinv succeeded! => continue.\n";
        ierr = 0;
        calcerr = true;
      }
    }
    
// *-- Evaluate r and S.
    for (int k = 0; k < ncon; ++k) {
      R[k] = constraints[k]->getValue();
      if (debug>1) cout << "F[" << k << "] = " << R[k] << endl;
      for (int j = 0; j < nmea; ++j) {
        R[k] += dfeta[k][j]*(y[j]-eta[j]);
      }
      if (debug>1) cout << "R[" << k << "] = " << R[k] << endl;
      for (int l = 0; l < ncon; ++l) {
        S[k][l] = 0;
        for (int i = 0; i < nmea; ++i) {
          for (int j = 0; j < nmea; ++j) {
            S[k][l] += dfeta[k][i] * V[i][j] * dfeta[l][j]; 
          }
        }
        // New invention by B. List, 6.12.04:
        // add F_xi * F_xi^T to S, to make method work when some
        // constraints do not depend on any measured parameter
        for (int j = 0; j < nunm; ++j) {
          S[k][l] += dfeta[k][nmea+j] * dfeta[l][nmea+j]; 
        }
      }
    }
    
    
    
    if (debug>1) 
      for (int k = 0; k < ncon; ++k) 
        for (int l = 0; l < ncon; ++l) 
          if (S[k][l] != 0)
            cout << "S[" << k << "][" << l << "]=" << S[k][l] << endl;;
    
// *-- Invert S, testing for singularity first.
// skip singularity test, since it only fixes some machine dependency
// of DSINV/RSINV (-> cernlib)

// S is positive definite, and we can use dsinv
   dsinv(ncon, &S[0][0], NCONMAX, inverr);
   if (inverr != 0) {
     cerr << "S: dsinv error " << inverr << endl;
     ierr = 7;
     calcerr = false;
     break;
   }

// *-- Calculate new unmeasured quantities.
    if (nunm > 0) {
      for (int i = 0; i < nunm; ++i) {
        for (int j = 0; j < nunm; ++j) {
          W1[i][j] = 0;
          for (int k = 0; k < ncon; ++k) {
            for (int l = 0; l < ncon; ++l) {
              W1[i][j] += dfeta[k][nmea+i] * S[k][l] * dfeta[l][nmea+j];
            }
          }
          if (debug>1) cout << "W1[" << i << "][" << j << "] = " << W1[i][j] << endl;
        }
      }
      if (debug > 1) {
        // Check symmetry of W1
        for (int i = 0; i < nunm; ++i) {
          for (int j = 0; j < nunm; ++j) {
            if (abs(W1[i][j]-W1[j][i]) > 1E-3*abs(W1[i][j]+W1[j][i]))
              cout << "W1[" << i << "][" << j << "] = " << W1[i][j] 
                   << "   W1[" << j << "][" << i << "] = " << W1[j][i] 
                   << "   => diff=" << abs(W1[i][j]-W1[j][i])
                   << "   => tol=" << 1E-3*abs(W1[i][j]+W1[j][i])
                   << endl;
          }
       }
      }
      
      
      // invert W1
      // Note added 23.12.04: W1 is symmetric and positive definite,
      dsinv(nunm, &W1[0][0], NUNMMAX, inverr);
      if (inverr != 0) {
        cerr << "W1: dsinv error " << inverr << endl;
        ierr = 8;
        calcerr = false;
        break;
      }
      // calculate shift of unmeasured parameters
      for (int i = 0; i < nunm; ++i) {
        dxi[i] = 0;
        for (int j = 0; j < nunm; ++j) {
          for (int k = 0; k < ncon; ++k) {
            for (int l = 0; l < ncon; ++l) {
              dxi[i] -= alph * W1[i][j] * dfeta[k][nmea+j]* S[k][l] * R[l];
            }
          }
        }
        if (debug>1) cout << "dxi[" << i << "] = " << dxi[i] << endl;
      }
    }
    
// *-- And now update unmeasured parameter.
    for (int i = 0; i < nunm; ++i) {
      eta[nmea+i] += dxi[i];
    }
    
// *-- Calculate new Lagrange multipliers.
    for (int k = 0; k < ncon; ++k) {
      alam[k] = 0.;
      for (int l = 0; l < ncon; ++l) {
        alam[k] += S[k][l] * R[l];
        if (debug>2) cout << "alam[" << k << "] = " << alam[k] << endl;
        for (int j = 0; j < nunm; ++j) {
          alam[k] += S[k][l] * dfeta[l][nmea+j] * dxi[j];
          if (debug>2) cout << "change of alam[" << k << 
                             "] for j,l = " << j << " , " << l << 
                             " : " << S[k][l] * dfeta[l][nmea+j] * dxi[j] << endl;
        }
        if (debug>2) cout << "alam[" << k << "] = " << alam[k] << endl;
      }
      if (debug>1) cout << "alam[" << k << "] = " << alam[k] << endl;
    }
    
// *-- Calculate new fitted parameters.
    for (int i = 0; i < nmea; ++i) {
      eta[i] = y[i];
      for (int j = 0; j < nmea; ++j) {
        for (int k = 0; k < ncon; ++k) {
          eta[i] -= V[i][j] * dfeta[k][j] * alam[k];
        }
      }
      if (debug>1) cout << "updated eta[" << i << "] = " << eta[i] << endl;
    }
    
// *-- Calculate constraints and their derivatives.
    // since the constraints ask the fitobjects for their parameters, 
    // we need to update the fitobjects first!
    // COULD BE DONE: update also ERRORS! (now only in the very end!)
    updatesuccess = updateFitObjects (eta);
         
    if (debug) {
      cout << "After adjustment of all parameters:\n";
      for (int k = 0; k < ncon; ++k) {
        cout << "Value of constraint " << k << " = " << constraints[k]->getValue()
             << endl;
      }
    }
    for (int k=0; k < ncon; k++) {
      for (int j=0; j < npar; j++) dfeta[k][j]=0;
      constraints[k]->getDerivatives(NPARMAX, dfeta[k]);
      if (debug>1) 
        for (int j=0; j < npar; j++) 
          if (dfeta[k][j]!= 0) cout << "dfeta[" << k << "][" << j << "] = " << dfeta[k][j] << endl;;
    }

// *-- Calculate new chisq.
    chit = 0;
    for (int i = 0; i < nmea; ++i) 
      for (int j = 0; j < nmea; ++j) {
        double dchit = (y[i]-eta[i]) * VINV[i][j] * (y[j]-eta[j]);
        chit  +=  dchit;
        if (debug>1 && dchit != 0)
          cout << "chit for i,j = " << i << " , " << j << " = " 
               << (y[i]-eta[i]) * VINV[i][j] * (y[j]-eta[j]) << endl;
      }
    chik = 0;
    for (int k = 0; k < ncon; ++k) chik += std::abs(2*alam[k]*constraints[k]->getValue());
    
    chinew = chit + chik;
    
//*-- Calculate change in chisq, and check constraints are satisfied.
    nit++;
    
    bool sconv = (std::abs(chik-chik0) < dchikc) 
              && (std::abs(chit-chit0) < dchitc*chit) 
              && (chik < dchikt*chit);
    // Second convergence criterium:
    // If all constraints are fulfilled to better than 1E-8,
    // and all parameters have changed by less than 1E-8,
    // assume convergence
    // This criterium assumes that all constraints and all parameters
    // are "of order 1", i.e. their natural values are around 1 to 100,
    // as for GeV or radians
    double eps = 1E-6;
    bool sconv2 = true;
    for (int k = 0; sconv2 && (k < ncon); ++k) 
      sconv2 &= (std::abs(constraints[k]->getValue()) < eps);
    if (sconv2 && debug) 
      cout << "All constraints fulfilled to better than " << eps << endl;
       
    for (int j = 0; sconv2 && (j < npar); ++j) 
      sconv2 &= (std::abs(eta[j] - etasv[j]) < eps);
    if (sconv2 && debug) 
      cout << "All parameters stable to better than " << eps << endl;
    sconv |= sconv2;
             
    bool sbad  = (chik > dchik*chik0) 
              && (chik > dchikt*chit)
              && (chik > chik0 + 1.E-10);
              
    scut = false;
           
    if (nit > nitmax) {
// *-- Out of iterations
      repeat = false;
      ierr = 1;
    }  
    else if (sconv && updatesuccess) {
// *-- Converged
      repeat = false;
      ierr = 0;
    }  
    else if ( nit > 2 && chinew > chimxw  && updatesuccess) {
// *-- Chi2  crazy ?
      repeat = false;
      calcerr = false;
      ierr = 2;
    }  
    else if (sbad && nit > 1 || !updatesuccess) {
// *-- ChiK increased, try smaller step
      if ( alph == almin ) {
        repeat = true;   // false;
        calcerr = false;
        ierr = 3;
      }  
      else {
        alph  =  std::max (almin, 0.5 * alph);
        scut  =  true;
        repeat = true;
        ierr = 4;
      }
    }    
    else {
// *-- Keep going..
      alph  =  std::min (alph+0.1, 1. );
      repeat = true;
      ierr = 5;
    }
    
    if (debug) cout << "======== NIT = " << nit << ",  CHI2 = " << chinew 
                                     << ",  ierr = " << ierr << ", alph=" << alph << endl;
                                     
    if (debug) 
      for (unsigned int i = 0; i < fitobjects.size(); ++i) 
        cout << "fitobject " << i << ": " << *fitobjects[i] << endl;                                 

  }   // end of while (repeat)
  
// *-- End of iterations - calculate errors.
  FDouble VETA[NPARMAX][NPARMAX];
  for (int i = 0; i < npar; ++i) 
    for (int j = 0; j < npar; ++j) VETA[i][j] = 0;
 
  if (calcerr) {
  
// *-- Evaluate S and invert.
    for (int k = 0; k < ncon; ++k) {
      for (int l = 0; l < ncon; ++l) {
        S[k][l] = 0;
        for (int i = 0; i < npar; ++i) {
          for (int j = 0; j < npar; ++j) {
            S[k][l] += dfeta[k][i] * V[i][j] * dfeta[l][j]; 
          }
        }
        // New invention by B. List, 6.12.04:
        // add F_xi * F_xi^T to S, to make method work when some
        // constraints do not depend on any measured parameter
        for (int j = 0; j < nunm; ++j) {
          S[k][l] += dfeta[k][nmea+j] * dfeta[l][nmea+j]; 
        }
      }
    }
   
// *-- Invert S, testing for singularity first.
// skip singularity test, since it only fixes some machine dependency
// of DSINV/RSINV (-> cernlib)

// S is positive definite, and we can use dsinv
    dsinv(ncon, &S[0][0], NCONMAX, inverr);
    if (inverr != 0) cerr << "S(2): dsinv error " << inverr << endl;

// *-- Calculate G.  
//  (same as W1, but for measured parameters) 
    FDouble G[NPARMAX][NPARMAX];
    for (int i = 0; i < nmea; ++i) {
      for (int j = 0; j < nmea; ++j) {
        G[i][j] = 0;
        for (int k = 0; k < ncon; ++k) {
          for (int l = 0; l < ncon; ++l) {
            G[i][j] += dfeta[k][i] * S[k][l] * dfeta[l][j];
          }
        }
      }
    }
    
    if (nunm > 0) {

// *-- Calculate H.
      FDouble H[NPARMAX][NUNMMAX];
      for (int i = 0; i < nmea; ++i) {
        for (int j = 0; j < nunm; ++j) {
          H[i][j] = 0;
          for (int k = 0; k < ncon; ++k) {
            for (int l = 0; l < ncon; ++l) {
              H[i][j] += dfeta[k][i] * S[k][l] * dfeta[l][nmea+j];
            }
          }
        }
      }
      
// *-- Calculate U**-1 and invert.
//   (same as W1)
      FDouble U[NUNMMAX][NUNMMAX];
      for (int i = 0; i < nunm; ++i) {
        for (int j = 0; j < nunm; ++j) {
          U[i][j] = 0;
          for (int k = 0; k < ncon; ++k) {
            for (int l = 0; l < ncon; ++l) {
              U[i][j] += dfeta[k][nmea+i] * S[k][l] * dfeta[l][nmea+j];
            }
          }
        }
      }
      dsinv(nunm, &U[0][0], NUNMMAX, inverr);
      if (inverr != 0) {
        cerr << "U: dsinv error " << inverr << endl;
        return -1;
      }
     
// *-- U is now error matrix of unmeasured parameters
      for (int i = 0; i < nunm; ++i) {
        for (int j = 0; j < nunm; ++j) {
          VETA[nmea+i][nmea+j] = U[i][j];
        }
      }

// *-- Covariance matrix between measured and unmeasured parameters.
      for (int i = 0; i < nmea; ++i) {
        for (int j = 0; j < nunm; ++j) {
          VETA[i][nmea+j] = 0.;
          for (int ii = 0; ii < nmea; ++ii) {
            for (int jj = 0; jj < nunm; ++jj) {
              VETA[i][nmea+j] -= V[i][ii] * H[ii][jj] * U[jj][j];
            }
          }
        }
      }
      
// *-- Fill in symmetric part:
      for (int i = 0; i < nmea; ++i) {
        for (int j = 0; j < nunm; ++j) {
          VETA[nmea+j][i] = VETA[i][nmea+j];
        }
      }
      
// *-- Calculate G-HUH.
      for (int i = 0; i < nmea; ++i) {
        for (int j = 0; j < nmea; ++j) {
          for (int ii = 0; ii < nunm; ++ii) {
            for (int jj = 0; jj < nunm; ++jj) {
              G[i][j] -= H[i][ii] * U[ii][jj] * H[j][jj];
            }
          }
        }
      }
      
    }  // endif nunm > 0

// *-- Calculate I-GV.
    FDouble GV[NPARMAX][NPARMAX];
    for (int i = 0; i < nmea; ++i) {
      for (int j = 0; j < nmea; ++j) {
        if (i == j) GV[i][j] = 1.;
        else GV[i][j] = 0.;
        for (int ii = 0; ii < nmea; ++ii) {
          GV[i][j] -= G[i][ii] * V[ii][j];
        }
      }
    }

// *-- And finally error matrix on fitted parameters.
    for (int i = 0; i < nmea; ++i) {
      for (int j = 0; j < nmea; ++j) {
        VETA[i][j] = 0.;
        for (int ii = 0; ii < nmea; ++ii) {
          VETA[i][j] += V[i][ii] * GV[ii][j];
        }
      }
    }

    // update errors in fitobjects
    // (consider only diagonal elements of VETA for the moment...)
    for (unsigned int ifitobj = 0; ifitobj < fitobjects.size(); ++ifitobj) {
      for (int ilocal = 0; ilocal < fitobjects[ifitobj]->getNPar(); ++ilocal) {
        int iglobal = fitobjects[ifitobj]->getGlobalParNum (ilocal); 
        for (int jlocal = ilocal; jlocal < fitobjects[ifitobj]->getNPar(); ++jlocal) {
          int jglobal = fitobjects[ifitobj]->getGlobalParNum (jlocal); 
          if (iglobal >= 0 && jglobal >= 0) 
            fitobjects[ifitobj]->setCov(ilocal, jlocal, VETA[iglobal][jglobal]); 
        }
      }
    }
    
  } // endif calcerr == true

// *-- Turn chisq into probability.
  FReal chi = FReal(chinew);
  fitprob = (ncon-nunm > 0) ? prob(chi,ncon-nunm) : 0.5;
  chi2 = chinew;

  return fitprob;
    
};

bool OPALFitter::initialize() {
  // tell fitobjects the global ordering of their parameters:
  int iglobal = 0;
  // measured parameters first!
  for (unsigned int ifitobj = 0; ifitobj < fitobjects.size(); ++ifitobj) {
    for (int ilocal = 0; ilocal < fitobjects[ifitobj]->getNPar(); ++ilocal) {
      if (fitobjects[ifitobj]->isParamMeasured(ilocal) &&
          !fitobjects[ifitobj]->isParamFixed(ilocal)) {
        fitobjects[ifitobj]->setGlobalParNum (ilocal, iglobal);
        if (debug) 
          cout << "Object " << fitobjects[ifitobj]->getName()
               << " Parameter " << fitobjects[ifitobj]->getParamName(ilocal)
               << " is measured, global number " << iglobal << endl;
        ++iglobal;
      }
    }
  }
  nmea = iglobal;
  // now  unmeasured parameters!
  for (unsigned int ifitobj = 0; ifitobj < fitobjects.size(); ++ifitobj) {
    for (int ilocal = 0; ilocal < fitobjects[ifitobj]->getNPar(); ++ilocal) {
      if (!fitobjects[ifitobj]->isParamMeasured(ilocal) &&
          !fitobjects[ifitobj]->isParamFixed(ilocal)) {
        fitobjects[ifitobj]->setGlobalParNum (ilocal, iglobal);
        if (debug) 
          cout << "Object " << fitobjects[ifitobj]->getName()
               << " Parameter " << fitobjects[ifitobj]->getParamName(ilocal)
               << " is unmeasured, global number " << iglobal << endl;
        ++iglobal;
      }
    }
  }
  npar = iglobal;
  assert (npar <= NPARMAX);
  nunm = npar - nmea;
  assert (nunm <= NUNMMAX);
  
  // set number of constraints
  ncon = constraints.size();
  assert (ncon <= NCONMAX);
  
  return true;

};
  
bool OPALFitter::updateFitObjects (double eta[]) {
  bool debug = false;
  bool result = true;
  for (unsigned int ifitobj = 0; ifitobj < fitobjects.size(); ++ifitobj) {
    for (int ilocal = 0; ilocal < fitobjects[ifitobj]->getNPar(); ++ilocal) {
      fitobjects[ifitobj]->updateParams (eta, NPARMAX);
//       int iglobal = fitobjects[ifitobj]->getGlobalParNum (ilocal); 
//       if (!fitobjects[ifitobj]->isParamFixed (ilocal) && iglobal >= 0) {
//         if (debug) cout << "Parameter " << iglobal 
//                         << " (" << fitobjects[ifitobj]->getName()
//                         << ": " << fitobjects[ifitobj]->getParamName(ilocal)
//                         << ") set to " << eta[iglobal];
//         result &= fitobjects[ifitobj]->setParam(ilocal, eta[iglobal]); 
//         eta[iglobal] = fitobjects[ifitobj]->getParam(ilocal);
//         if (debug) cout << " => " << eta[iglobal] << endl;
//       }
    }
  }
  return result;
};

int OPALFitter::getError() const {return ierr;}
double OPALFitter::getProbability() const {return fitprob;}
double OPALFitter::getChi2() const {return chi2;}
int OPALFitter::getDoF() const {return ncon-nunm;}
int OPALFitter::getIterations() const {return nit;}
