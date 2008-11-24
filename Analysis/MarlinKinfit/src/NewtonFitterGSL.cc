/*! \file 
 *  \brief Implements class NewtonFitterGSL
 *
 * Author: Benno List
 * $Date: 2008-11-24 11:01:01 $
 * $Author: beckmann $
 *
 * \b Changelog:
 * - 26.9.08 BL: Correct parameter counting (discared fixed parameters)
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.1  2008/10/16 08:13:44  blist
 * - New versions of OPALfitter and Newtonfitter using GSL
 * -
 * - Revision 1.4  2008/09/26 09:58:10  boehmej
 * - removed ~100 semicolons after } at end of function implementation :)
 * -
 * - Revision 1.3  2008/09/26 07:30:54  blist
 * - Added documentation,
 * - removed bug in Newtonfitter (could not handle fixed parameters correctly)
 *
 */ 

#undef NDEBUG

#include "NewtonFitterGSL.h" 

#include<iostream>
#include<cmath>
#include<cassert>

#include "BaseFitObject.h"
#include "BaseHardConstraint.h"
#include "BaseSoftConstraint.h"
#include "ftypes.h"
#include "cernlib.h"

#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using std::cout;
using std::cerr;
using std::endl;

static int debug = 0;
static int nitdebug = 100;
static int nitcalc = 0;
static int nitsvd = 0;

// constructor
NewtonFitterGSL::NewtonFitterGSL() 
: npar (0), ncon (0), nsoft (0),
  x(0), dx(0), y(0), perr(0), v1 (0), v2(0), Meval (0), 
  M(0), M1(0), Mevec (0), permM(0), ws (0), wsdim (0), idim (0)
{}

// destructor
NewtonFitterGSL::~NewtonFitterGSL() {

  if (x) gsl_vector_free (x);               x=0;
  if (dx) gsl_vector_free (dx);             dx=0;
  if (y) gsl_vector_free (y);               y=0;
  if (perr) gsl_vector_free (perr);         perr=0;
  if (v1) gsl_vector_free (v1);             v1=0;
  if (v2) gsl_vector_free (v2);             v2=0;
  if (Meval) gsl_vector_free (Meval);       Meval=0;
  if (M) gsl_matrix_free (M);               M=0;
  if (M1) gsl_matrix_free (M1);             M1=0;
  if (Mevec) gsl_matrix_free (Mevec);       Mevec=0;
  if (permM) gsl_permutation_free (permM);  permM=0;
  if (ws) gsl_eigen_symmv_free (ws);        ws=0; wsdim=0;
}



double NewtonFitterGSL::fit() {

  // order parameters etc
  initialize();
  
  // initialize eta, etasv, y   
  assert (x && x->size == idim);
  assert (dx && dx->size == idim);
  assert (y && y->size == idim);
  assert (perr && perr->size == idim);
  assert (v1 && v1->size == idim);
  assert (v2 && v2->size == idim);
  assert (Meval && Meval->size == idim);
  assert (M && M->size1 == idim && M->size1 == idim);
  assert (M1 && M1->size1 == idim && M1->size1 == idim);
  assert (Mevec && Mevec->size1 == idim && Mevec->size1 == idim);
  assert (permM && permM->size == idim);
  
  gsl_vector_set_zero (x);
  gsl_vector_set_zero (y);
  gsl_vector_set_all (perr, 1);
  
  // Get starting values
  for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
    BaseFitObject *fo = *i;
    assert (fo);
    for (int ilocal = 0; ilocal < fo->getNPar(); ++ilocal) {
      if (!fo->isParamFixed(ilocal)) {
        int iglobal = fo->getGlobalParNum (ilocal);
        assert (iglobal >= 0 && iglobal < npar);
        gsl_vector_set (x, iglobal, fo->getParam (ilocal));
        double e = std::abs(fo->getError (ilocal));
        gsl_vector_set (perr, iglobal, e ? e : 1);
      }
    }
  }
  for (ConstraintIterator i = constraints.begin(); i != constraints.end(); ++i) {
    BaseHardConstraint *c = *i;
    assert (c);
    int iglobal = c->getGlobalNum ();
    assert (iglobal >= 0 && iglobal < (int)idim);
    double e =  c->getError();
    gsl_vector_set (perr, iglobal, e ? e : 1);
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
      for (int ilocal = 0; ilocal < fo->getNPar(); ++ilocal) {
        if (!fo->isParamFixed(ilocal)) {
          int iglobal = fo->getGlobalParNum (ilocal);
          cout << "  par " << fo->getParamName(ilocal) << ": local: " << ilocal << ": global: " << iglobal
               << " value=" << fo->getParam (ilocal) << " +- " << fo->getError(ilocal);
          if (fo->isParamMeasured(ilocal)) 
            cout << " measured: " << fo->getMParam (ilocal);
          cout << endl;
        }
        else {
          cout << "  par " << fo->getParamName(ilocal) << ": local: " << ilocal << " -- fixed -- " 
               << " value=" << fo->getParam (ilocal) << " +- " << fo->getError(ilocal)
               << endl;
        
        }
      }
    }
    cout << "constraints:\n";
    for (ConstraintIterator i = constraints.begin(); i != constraints.end(); ++i) {
      BaseHardConstraint *c = *i;
      assert (c);
      cout << i-constraints.begin() << " " << c->getName() << ": " << c->getValue() << "+-" << c->getError() << endl;
      int kglobal = c->getGlobalNum();
      cout << "  global number: " << kglobal << endl;
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
    
    if (debug>1 && (nit==0 || nit<nitdebug)) cout << "===================\nStarting iteration " << nit << endl;
    if (debug>2 && (nit==0 || nit<nitdebug)) {
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
 
    // Set M to 0
    gsl_matrix_set_zero (M);
    // Reset y
    gsl_vector_set_zero (y);
    
    // Compose M: 
    // First, all terms d^2 chi^2/dx1 dx2
    for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
      BaseFitObject *fo = *i;
      assert (fo);
      fo->addToGlobalChi2DerMatrix (M->block->data, M->tda);
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
      assert (kglobal >= 0 && kglobal < (int)idim);
      c->add1stDerivativesToMatrix (M->block->data, M->tda);
      c->add2ndDerivativesToMatrix (M->block->data, M->tda, gsl_vector_get (x, kglobal));
    }
//     cout << "After adding derivatives of constraints::\n";
//     printMy (M, y, idim);
    
    // Now, calculate the result vector y with the values of the derivatives
    // d chi^2/d x
    // First, for the parameters
    for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
      BaseFitObject *fo = *i;
      assert (fo);
      fo->addToGlobalChi2DerVector (y->block->data, y->size);
    }
//     cout << "After adding fo derivatives to y::\n";
//     printMy (M, y, idim);
    
    // Now add lambda*derivatives of constraints,
    // And finally, the derivatives w.r.t. to the constraints, i.e. the constraints themselves
    for (unsigned int k = 0; k < constraints.size(); ++k) {
      BaseHardConstraint *c = constraints[k];
      assert (c);
      int kglobal = c->getGlobalNum();
      assert (kglobal >= 0 && kglobal < (int)idim);
      c->addToGlobalChi2DerVector (y->block->data, y->size, gsl_vector_get (x, kglobal));
      gsl_vector_set (y, kglobal, c->getValue());
    }
  
    // Finally, treat the soft constraints

    for (SoftConstraintIterator i = softconstraints.begin(); i != softconstraints.end(); ++i) {
      BaseSoftConstraint *bsc = *i;
      assert (bsc);
      bsc->add2ndDerivativesToMatrix (M->block->data, M->tda);
      bsc->addToGlobalChi2DerVector (y->block->data, y->size);
    }
    
    if (debug>3 && (nit==0 || nit<nitdebug)) {
      cout << "After setting up equations: \n";
      debug_print (M, "M");
      debug_print (y, "y");
      debug_print (x, "x");
    }
  
    int ifail = calcDx ();
    if (ifail != 0) {
      cout << "NewtonFitterGSL::fit: calcDx: ifail=" << ifail << endl;
      ierr = 2;
      break;
    }
    
    gsl_vector_sub (x, dx);
    
    if (debug>3 && (nit==0 || nit<nitdebug)) {
      cout << "After solving equations: \n";
      debug_print (dx, "dx");
      debug_print (x, "x");
    }
    
//      cout << "new x = (";
//      for (int i = 0; i < idim; ++i) cout << x[i] << (i+1<idim?", ":")\n");
  
    // Update values in Fitobjects
    bool significant = false;
    for (FitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
      BaseFitObject *fo = *i;
      assert (fo);
      bool s = fo->updateParams (x->block->data, x->size);
      significant |=  s;
      if (debug && nit<nitdebug && s) 
        cout << "Significant update for FO " << i-fitobjects.begin() << " (" 
             << fo->getName() << ")\n";
    }
  
//   *-- Convergence criteria 

    chi2new = calcChi2();
    if (debug && nit<nitdebug) cout << "old chi2: " << chi2old << ", new chi2: " << chi2new << ", diff=" << chi2old-chi2new << endl;
    converged = !significant;
    bool constraintsok = true;
    if (true || converged || debug>1) {
      for (ConstraintIterator ci = constraints.begin(); ci != constraints.end(); ++ci) {
        BaseHardConstraint *c = *ci;
        assert (c);
        constraintsok &= (std::abs (c->getValue()) <= 1E-2*c->getError());
        if (debug && nit<nitdebug && std::abs (c->getValue()) >= 1E-2*c->getError())
          cout << "Constraint " << ci-constraints.begin() << " " << c->getName()
               << ": not fulfilled: " << c->getValue() << "+-" << c->getError() << endl;
          if (debug > 2 && nit<nitdebug) cout << "Constraint " << ci-constraints.begin() << ": " << c->getValue()
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

// *-- Turn chisq into probability.
  FReal chi = FReal(chi2new);

  return fitprob = prob(chi ,ncon+nsoft-nunm);
    
}

bool NewtonFitterGSL::initialize() {
//  bool debug = true;

  // tell fitobjects the global ordering of their parameters:
  npar = 0;
  nunm = 0;
  // 
  for (unsigned int ifitobj = 0; ifitobj < fitobjects.size(); ++ifitobj) {
    for (int ilocal = 0; ilocal < fitobjects[ifitobj]->getNPar(); ++ilocal) {
      if (!fitobjects[ifitobj]->isParamFixed(ilocal)) {
        fitobjects[ifitobj]->setGlobalParNum (ilocal, npar);
        ++npar;
        if (!fitobjects[ifitobj]->isParamMeasured(ilocal)) ++nunm;
      }
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
    cerr << "NewtonFitterGSL::initialize: nunm=" << nunm << " > ncon+nsoft=" 
         << ncon << "+" << nsoft << endl;
  }
  
  idim = npar+ncon;
  
  ini_gsl_vector (x, idim);
  ini_gsl_vector (dx, idim);
  ini_gsl_vector (y, idim);
  ini_gsl_vector (perr, idim);
  ini_gsl_vector (v1, idim);
  ini_gsl_vector (v2, idim);
  ini_gsl_vector (Meval, idim);
  
  ini_gsl_matrix (M, idim, idim);
  ini_gsl_matrix (M1, idim, idim);
  ini_gsl_matrix (Mevec, idim, idim);
  
  ini_gsl_permutation (permM, idim);
  
  if (ws && wsdim != idim) {
    gsl_eigen_symmv_free (ws); 
    ws = 0;
  }
  if (ws == 0) ws = gsl_eigen_symmv_alloc (idim); 
  wsdim = idim;
 
  return true;

}
  
double NewtonFitterGSL::calcChi2() {
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

void NewtonFitterGSL::printMy (double M[], double y[], int idim) {
  for (int i = 0; i < idim; ++i) {
    cout << i << "  [ " << M[idim*i+0];
      for (int j = 1; j<idim; ++j) cout << ", " << M[idim*i+j];
      cout << "]  [" << y[i] << "]\n";
  }
}

int NewtonFitterGSL::getError() const {return ierr;}
double NewtonFitterGSL::getProbability() const {return fitprob;}
double NewtonFitterGSL::getChi2() const {return chi2;}
int NewtonFitterGSL::getDoF() const {return ncon+nsoft-nunm;}
int NewtonFitterGSL::getIterations() const {return nit;}

int NewtonFitterGSL::calcDx () const {
    if (debug>1)cout << "entering calcDx" << endl;
  
    nitcalc++;
    // from x_(n+1) = x_n - y/y' = x_n - M^(-1)*y we have M*(x_n-x_(n+1)) = y, 
    // which we solve for dx = x_n-x_(n+1) and hence x_(n+1) = x_n-dx
  
    gsl_matrix_memcpy (M1, M);
    
//     cout << "Complete system:\n";
//     printMy(M, y, idim);
    
    FInteger ifail = 0;
    
    int signum;
    int result = gsl_linalg_LU_decomp (M1, permM, &signum);
    if (debug>1)cout << "calcDx: gsl_linalg_LU_decomp result=" << result << endl;
    // SOlve M1*dx = y
    ifail = gsl_linalg_LU_solve (M1, permM, y, dx);
    if (debug>1)cout << "calcDx: gsl_linalg_LU_solve result=" << ifail << endl;
    
    if (ifail != 0) {
      cerr << "NewtonFitterGSL::calcDx: ifail from gsl_linalg_LU_solve=" << ifail << endl;
      return calcDxSVD ();
    }
//     cout << "Result of dseqn: " << ifail << endl;
//     cout << "dx = (";
//    for (int i = 0; i < idim; ++i) cout << gsl_vector_get (dx, i) << (i+1<idim?", ":")\n");
    
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
    for (unsigned int i = 0; i < idim; ++i) {
      if(gsl_vector_get (perr, i)>0 && std::abs(gsl_vector_get (dx, i))>maxsigma*gsl_vector_get (perr, i)) 
        maxsigma=std::abs(gsl_vector_get (dx, i))/gsl_vector_get (perr, i);
        if (debug && nitcalc<nitdebug && gsl_vector_get (perr, i)>0 && std::abs(gsl_vector_get (dx, i)) > maxstep*gsl_vector_get (perr, i)) 
          cout << "step size for parameter " << i << ": " << gsl_vector_get (dx, i)
               << " / " << gsl_vector_get (perr, i) << " = " << gsl_vector_get (dx, i)/gsl_vector_get (perr, i) << endl;
    }
    double scale = 1.0;    
    if (maxsigma > maxstep) {
      if (debug && nitcalc<nitdebug) cout << "Step size limitation: maxsigma=" << maxsigma << endl;
      scale = maxstep/maxsigma;
    }   
    
    if (scale < 0.5) {
      if (debug > 1) cout << "NewtonFitterGSL::calcDx: reverting to calcDxSVD\n";
      return calcDxSVD ();
    }
    
    gsl_blas_dscal  (scale, dx);
        
    return 0;
}

int NewtonFitterGSL::calcDxSVD () const {
    if (debug>1)cout << "entering calcDxSVD" << endl;

    nitsvd++;
    // from x_(n+1) = x_n - y/y' = x_n - M^(-1)*y we have M*(x_n-x_(n+1)) = y, 
    // which we solve for dx = x_n-x_(n+1) and hence x_(n+1) = x_n-dx
  
    for (unsigned int i = 0; i < idim; ++i) assert (gsl_vector_get (perr, i) > 0);
    
    // Rescale columns and rows by perr
    for (unsigned int i = 0; i < idim; ++i) 
      for (unsigned int j = 0; j < idim; ++j) 
        gsl_matrix_set (M1, i, j,  
          gsl_vector_get (perr, i)*gsl_vector_get (perr, j)*gsl_matrix_get (M, i, j));
    
//     cout << "Complete system:\n";
//     printMy(M, y, idim);
     // Get eigenvalues and eigenvectors of M1
     int ierr=0;
     ierr = gsl_eigen_symmv (M1, Meval, Mevec, ws); 
     if (ierr != 0) {
       cerr << "NewtonFitterGSL::calcDxSVD: ierr=" << ierr << "from gsl_eigen_symmv!\n";
     }
     // Sort the eigenvalues and eigenvectors in descending order in magnitude
     ierr = gsl_eigen_symmv_sort (Meval, Mevec, GSL_EIGEN_SORT_ABS_DESC);
     if (ierr != 0) {
       cerr << "NewtonFitterGSL::calcDxSVD: ierr=" << ierr << "from gsl_eigen_symmv_sort!\n";
     }
     
     
     // The eigenvectors are stored in the columns of Mevec;
     // the eigenvectors are orthonormal, therefore Mevec^T = Mevec^-1
     // Therefore, M1 = Mevec * diag(Meval)* Mevec^T, and
     // M1^-1 = (Mevec^T)^-1 * diag(Meval)^-1 * Mevec^-1 
     //       =  Mevec * diag(Meval)^-1 * Mevec^T
     // So, the solution of M1*dx = y is given by
     // dx = M1^-1 * y = Mevec * diag(Meval)^-1 * Mevec^-1 *y
     //    = Mevec * v2
     // For the pseudoinverse, the last elements of Meveal^-1 are set
     // to 0, therefore the last elements of v2 will be 0, 
     // therefore we can restrict the calculation of Mevec * v2
     // to the first ndim rows.
     // So, we calculate v2 only once, with only the inverse of zero eigenvalues
     // set to 0, and then calculate Mevec * v2 for fewer and fewer rows
   
//    cout << "calcDxSVD: info = " << info << endl;
//    cout << " s = ";
//    for (int i = 0; i < idim; ++i) cout << s[i] << ", ";
//    cout << endl;
   
   // Now M1 = U * s * V^T
   // We want to solve M1*dx = y, hence dx = V * s^-1 * U^T * y
   // Calculate UTy first; we need only the first ndim entries
   unsigned int ndim = 0;
   while (ndim<idim && gsl_vector_get (Meval, ndim) != 0) ++ndim;
   
   if (ndim < idim) {
     cout << "calcDxSVD: idim = " << idim << " > ndim = " << ndim << endl;
     cout << " Meval = ";
     for (unsigned int i = 0; i < idim; ++i) cout << gsl_vector_get (Meval, ndim) << ", ";
     cout << endl;
   }
   
   // Calculate v1 = y*perr (component wise)
   gsl_vector_memcpy (v1, y);
   gsl_vector_mul (v1, perr);
   // Calculate v2 = 1*Mevec^T*v1 + 0*v2
    gsl_blas_dgemv (CblasTrans, 1, Mevec, v1, 0, v2);
    
   // Divide by nonzero eigenvalues
   for (unsigned int i = 0; i<idim; ++i) {
     if (double e = gsl_vector_get (Meval, i)) gsl_vector_set (v2, i, gsl_vector_get (v2, i)/e);
     else gsl_vector_set (v2, i, 0);
   }
      
   double maxstep = 4.0;
   double maxsigma = 0;
   double scale = 1.0;    
   do { 
   
     gsl_vector_view v2part = gsl_vector_subvector (v2, 0, ndim);
     gsl_matrix_view Mevecpart = gsl_matrix_submatrix (Mevec, 0, 0, idim, ndim);
     
     // Calculate dx = 1*Mevecpart^T*v2 + 0*dx
     gsl_blas_dgemv (CblasNoTrans, 1, &Mevecpart.matrix, &v2part.vector, 0, dx);
     
      // Step size limitation: Calculate max number of sigmas
      maxsigma = 0;
      for (unsigned int i = 0; i < idim; ++i) {
        if(std::abs(gsl_vector_get (dx, i))>maxsigma) 
          maxsigma=std::abs(gsl_vector_get (dx, i));
          if (debug && nitsvd<nitdebug && std::abs(gsl_vector_get (dx, i)) > maxstep) 
            cout << "step size for parameter " << i << ": " << gsl_vector_get (dx, i)*gsl_vector_get (perr, i)
                 << " / " << gsl_vector_get (perr, i) << " = " << gsl_vector_get (dx, i) << endl;
      }
      scale = 1.0;    
      if (maxsigma > maxstep) {
        if (debug && nitsvd<nitdebug) cout << "Step size limitation: maxsigma=" << maxsigma << endl;
        scale = maxstep/maxsigma;
      }  
      --ndim;
      
      if (debug>1 && (scale < 0.5 || ndim < idim-1)) {
        cout << "ndim=" << ndim << ", scale=" << scale << endl;
      }
      
    } 
    while (ndim > 0 && scale < 0.5);
    
    // dx = dx*perr (component wise)
    gsl_vector_mul (dx, perr);
     // dx = dx*scale 
    gsl_blas_dscal  (scale, dx);
   
    return 0;
}

void NewtonFitterGSL::ini_gsl_permutation (gsl_permutation *&p, unsigned int size) {
  if (p) {
    if (p->size != size) {
      gsl_permutation_free (p);
      if (size > 0) p = gsl_permutation_alloc (size);
      else p=0;
    }
  }
  else 
    if (size > 0) p = gsl_permutation_alloc (size);
}

void NewtonFitterGSL::ini_gsl_vector (gsl_vector *&v, unsigned int size) {
  
  if (v) {
    if (v->size != size) {
      gsl_vector_free (v);
      if (size>0) v = gsl_vector_alloc (size);
      else v=0;
    }
  }
  else 
    if (size > 0) v = gsl_vector_alloc (size);
}

void NewtonFitterGSL::ini_gsl_matrix (gsl_matrix *&m, unsigned int size1, unsigned int size2) {
  if (m) {
    if (m->size1 != size1 || m->size2 != size2) {
      gsl_matrix_free (m);
      if (size1*size2 > 0) m = gsl_matrix_alloc (size1, size2);
      else m=0;
    }
  }
  else 
    if (size1*size2 > 0) m = gsl_matrix_alloc (size1, size2);
}

void NewtonFitterGSL::debug_print (gsl_matrix *m, const char *name) {
  for (unsigned int  i = 0; i < m->size1; ++i) 
    for (unsigned int j = 0; j < m->size2; ++j)
      if (gsl_matrix_get (m, i, j) != 0)
        cout << name << "[" << i << "][" << j << "]=" << gsl_matrix_get (m, i, j) << endl;
}

void NewtonFitterGSL::debug_print (gsl_vector *v, const char *name) {
  for (unsigned int  i = 0; i < v->size; ++i) 
      if (gsl_vector_get (v, i) != 0)
        cout << name << "[" << i << "]=" << gsl_vector_get (v, i) << endl;
}

int NewtonFitterGSL::getNcon() const {return ncon;}
int NewtonFitterGSL::getNsoft() const {return nsoft;}
int NewtonFitterGSL::getNunm() const {return nunm;}
int NewtonFitterGSL::getNpar() const {return npar;}
