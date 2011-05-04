/*! \file 
 *  \brief Implements class ParticleFitObject
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: ParticleFitObject.cc,v $
 * - Revision 1.7  2009/02/17 12:46:35  blist
 * - Improved version of NewtonFitterGSL, JetFitObject changed
 * -
 * - Revision 1.6  2009/02/11 15:33:49  mbeckman
 * - Bug fixes: mass initialization in ParticleFitObject, parameter handling in PhotonFitObjectPxyg
 * -
 * - Revision 1.5  2008/11/23 17:53:41  mbeckman
 * - Fixed minor bug in ParticleFitObject.cc
 * -
 * - Revision 1.4  2008/10/17 13:17:17  blist
 * - Avoid variable-size arrays
 * -
 * - Revision 1.3  2008/10/16 08:13:44  blist
 * - New versions of OPALfitter and Newtonfitter using GSL
 * -
 * - Revision 1.2  2008/09/26 09:58:11  boehmej
 * - removed ~100 semicolons after } at end of function implementation :)
 * -
 * - Revision 1.1  2008/02/12 10:19:09  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.7  2008/02/04 17:30:54  blist
 * - NewtonFitter works now!
 * -
 * - Revision 1.6  2008/01/30 21:48:03  blist
 * - Newton Fitter still doesnt work :-(
 * -
 * - Revision 1.5  2008/01/30 09:14:54  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.4  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.3  2007/09/13 13:33:06  blist
 * - Print methods return os
 * -
 * - Revision 1.2  2007/09/13 08:09:51  blist
 * - Updated 2nd derivatives for px,py,pz,E constraints, improved header documentation
 * -
 *
 */ 
 
#include "ParticleFitObject.h"
#include "cernlib.h"

#include <iostream>
#include <cassert>
#include <cmath>
using std::isfinite;
using std::cout; 
using std::endl;


ParticleFitObject::ParticleFitObject()
: mass (0)
{
  for (int ilocal = 0; ilocal < NPAR; ++ilocal) globalParNum[ilocal] = -1;
  for (int ilocal = 0; ilocal < NPAR; ++ilocal) fixed[ilocal] = false;
  for (int ilocal = 0; ilocal < NPAR; ++ilocal) 
    for (int jlocal = 0; jlocal < NPAR; ++jlocal) cov[ilocal][jlocal] = 0; 
}

ParticleFitObject::~ParticleFitObject()
{}

bool ParticleFitObject::setParam (int ilocal, double par_, 
                                    bool measured_, bool fixed_) {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (measured[ilocal] != measured_ || fixed[ilocal] != fixed_) invalidateCache();
  measured[ilocal] = measured_;
  fixed[ilocal] = fixed_;
  return setParam (ilocal, par_);
}  

bool ParticleFitObject::setParam (int ilocal, double par_ ) {
  if (!isfinite(par_)) return true;
  assert (ilocal >= 0 && ilocal < NPAR);
  if (par[ilocal] == par_) return false;
  invalidateCache();
  bool result = (par_-par[ilocal])*(par_-par[ilocal]) > eps2*cov[ilocal][ilocal]; 
  par[ilocal] = par_;
  return result;
}
  
bool ParticleFitObject::setMParam (int ilocal, double mpar_ ) {
  if (!isfinite(mpar_)) return false;
  assert (ilocal >= 0 && ilocal < NPAR);
  if (mpar[ilocal] == mpar_) return true;
  invalidateCache();
  mpar[ilocal] = mpar_;
  return true;
}

bool ParticleFitObject::setError (int ilocal, double err_) {
  if (!isfinite(err_)) return false;
  assert (ilocal >= 0 && ilocal < NPAR);
  invalidateCache();
  covinvvalid = false;
  cov[ilocal][ilocal] = err_*err_;
  return true;
}

bool ParticleFitObject::setCov (int ilocal, int jlocal, double cov_) {
  if (!isfinite(cov_)) return false;
  assert (ilocal >= 0 && ilocal < NPAR);
  assert (jlocal >= 0 && jlocal < NPAR);
  invalidateCache();
  covinvvalid = false;
  cov[ilocal][jlocal] = cov[jlocal][ilocal] = cov_;
  return true;
}
bool ParticleFitObject::setMass (double mass_) {
  if (!isfinite(mass_)) return false;
  if (mass == mass_) return true;
  invalidateCache();
  mass = std::abs(mass_);
  return true;
}
double ParticleFitObject::getMass () const {
  return mass;
}
bool ParticleFitObject::fixParam (int ilocal, bool fix) {
  assert (ilocal >= 0 && ilocal < NPAR);
  return fixed [ilocal] = fix;
}

bool ParticleFitObject::setGlobalParNum (int ilocal, int iglobal) {
  if (ilocal < 0 || ilocal >= NPAR) return false;
  globalParNum[ilocal] = iglobal;
  return true;
}
int  ParticleFitObject::getGlobalParNum(int ilocal) const {
  if (ilocal < 0 || ilocal >= NPAR) return -1;
  return globalParNum[ilocal];
}

double ParticleFitObject::getParam (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  return par[ilocal];
}
double ParticleFitObject::getMParam (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  return mpar[ilocal];
}

double ParticleFitObject::getError (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  return std::sqrt(cov[ilocal][ilocal]);
}
double ParticleFitObject::getCov (int ilocal, int jlocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  assert (jlocal >= 0 && jlocal < NPAR);
  return cov[ilocal][jlocal];
}
bool ParticleFitObject::isParamMeasured (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  return measured[ilocal];
}

bool ParticleFitObject::isParamFixed (int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  return fixed[ilocal];
}
    
std::ostream&  ParticleFitObject::print4Vector(std::ostream& os) const {
  os << "[" << getE() << ", " << getPx() << ", " 
      << getPy() << ", "  << getPz() << "]";
  return os;
}

std::ostream&  ParticleFitObject::print (std::ostream& os) const {
  printParams(os);
  os << " => ";
  print4Vector(os);
  return os;
}
bool ParticleFitObject::calculateCovInv() const {
  int n = getNPar();
  int idim = 0;
  for (int i = 0; i < n; ++i) {
    if (isParamMeasured (i)) {
      idim = i;
      for (int j = 0; j < n; ++j) {
        covinv[i][j] = isParamMeasured (j) ? cov[i][j] : 0;
      }
    }
    else {
      for (int j = 0; j < n; ++j) {
        covinv[i][j] = static_cast<double>(i == j);
      }
    }
  }
  int ierr = (idim == 0) ? 0 : dsinv(idim+1, &covinv[0][0], NPAR);
  if (ierr != 0) {
    //std::cerr << "ParticleFitObject::calculateCovInv: Error "
    //          << ierr << " from dsinv! Object " << getName() << std::endl;
    // printCov (std::cerr);         
  }
  return covinvvalid = (ierr == 0);
}


double ParticleFitObject::getChi2() const {
  if (!covinvvalid) calculateCovInv();
  if (!covinvvalid) return -1;
  double chi2 = 0;
  static double resid[NPAR];
  static bool chi2contr[NPAR];
  for (int i = 0; i < getNPar(); ++i) {
    resid[i] = par[i]-mpar[i];
    if (chi2contr[i] = isParamMeasured(i) && !isParamFixed(i)) {
      chi2 += resid[i]*covinv[i][i]*resid[i];
      for (int j = 0; j < i; ++j) {
        if (chi2contr[j]) chi2 += 2*(resid[i])*covinv[i][j]*(resid[j]);
      }
    }
  }
  return chi2;
}


double ParticleFitObject::getDChi2DParam(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (isParamFixed(ilocal) || !isParamMeasured(ilocal)) return 0;
  if (!covinvvalid) calculateCovInv();
  if (!covinvvalid) return 0;
  double result = 0;
  for (int jlocal = 0; jlocal < getNPar(); jlocal++) 
    if (!isParamFixed(jlocal) && isParamMeasured(jlocal))
      result += covinv[ilocal][jlocal]*(par[jlocal]-mpar[jlocal]);
  return 2*result;
}

void ParticleFitObject::addToGlobalChi2DerVector (double *y, int idim) const {
  assert (getNPar() <= NPAR);
  if (!covinvvalid) calculateCovInv();
  assert (covinvvalid);
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    if (!isParamFixed(ilocal) && isParamMeasured(ilocal)) {
      int iglobal = getGlobalParNum (ilocal);
      assert (iglobal>= 0 && iglobal < idim);
      for (int jlocal = 0; jlocal < getNPar(); jlocal++) {
        if (!isParamFixed(jlocal) && isParamMeasured(jlocal)) {
//           cout << "ParticleFitObject::addToGlobalChi2DerVector: fo=" << getName()
//                << ", ilocal=" << ilocal << ", jlocal=" << jlocal
//                << ", covinv=" << covinv[ilocal][jlocal]
//                << ", resid = " << par[jlocal]-mpar[jlocal]
//                << ", term=" << 2*covinv[ilocal][jlocal]*(par[jlocal]-mpar[jlocal])
//                << ", par=" << par[jlocal]
//                << ", mpar=" << mpar[jlocal]
//                << endl;
          y[iglobal] += 2*covinv[ilocal][jlocal]*(par[jlocal]-mpar[jlocal]);
        }
      }
    }
  }
}

void ParticleFitObject::addToGlobalChi2DerVectorNum (double *y, int idim, double eps)  {
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    int iglobal = getGlobalParNum(ilocal);
    y[iglobal] += num1stDerivative (ilocal, eps);
  }
}
    
double ParticleFitObject::getD2Chi2DParam2(int ilocal, int jlocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  assert (jlocal >= 0 && jlocal < NPAR);
  if (isParamFixed(ilocal) || !isParamMeasured(ilocal) || 
      isParamFixed(jlocal) || !isParamMeasured(jlocal))
    return 0;
  if (!covinvvalid) calculateCovInv();
  if (!covinvvalid) return 0;
  return 2*covinv[ilocal][jlocal];
}
      
void ParticleFitObject::addToGlobalChi2DerMatrix (double *M, int idim) const {
  if (!covinvvalid) calculateCovInv();
  if (!covinvvalid) return;
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    if (!isParamFixed(ilocal) && isParamMeasured(ilocal)) {
      int iglobal = getGlobalParNum (ilocal);
      assert (iglobal >= 0 && iglobal < idim);
      int ioffs = idim*iglobal;
      for (int jlocal = 0; jlocal < getNPar(); ++jlocal) {
        if (!isParamFixed(jlocal) && isParamMeasured(jlocal)) {
          int jglobal = getGlobalParNum (jlocal);
          assert (jglobal >= 0 && jglobal < idim);
          M[ioffs+jglobal] += 2*covinv[ilocal][jlocal];
        }
      }
    }
  }
}
      
void ParticleFitObject::addToGlobalChi2DerMatrixNum (double *M, int idim, double eps) {
  for (int ilocal1 = 0; ilocal1 < getNPar(); ++ilocal1) {
    int iglobal1 = getGlobalParNum (ilocal1);
    for (int ilocal2 = ilocal1; ilocal2 < getNPar(); ++ilocal2) {
      int iglobal2 = getGlobalParNum (ilocal2);
      M[idim*iglobal1 + iglobal2]+= num2ndDerivative (ilocal1, eps, ilocal2, eps);
    }
  }
}

void ParticleFitObject::addToGlobCov(double *globCov, int idim) const {
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    if (!isParamFixed(ilocal) && isParamMeasured(ilocal)) {
      int iglobal = getGlobalParNum (ilocal);
      assert (iglobal >= 0 && iglobal < idim);
      int ioffs = idim*iglobal;
      for (int jlocal = 0; jlocal < getNPar(); ++jlocal) {
        if (!isParamFixed(jlocal) && isParamMeasured(jlocal)) {
          int jglobal = getGlobalParNum (jlocal);
          assert (jglobal >= 0 && jglobal < idim);
          globCov[ioffs+jglobal] += cov[ilocal][jlocal];
        }
      }
    }
  }
}

void ParticleFitObject::getDerivatives (double der[], int idim) const {
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    assert (ilocal < idim);
    der [4*ilocal]   = getDE (ilocal);
    der [4*ilocal+1] = getDPx (ilocal);
    der [4*ilocal+2] = getDPy (ilocal);
    der [4*ilocal+3] = getDPz (ilocal);
  }
}


void ParticleFitObject::test1stDerivatives () {
  cout << "ParticleFitObject::test1stDerivatives, object " << getName() << "\n";
  double ycalc[100],ynum[100];
  for (int i = 0; i < 100; ++i) ycalc[i]=ynum[i]=0;
  addToGlobalChi2DerVector (ycalc, 100);
  double eps = 0.00001;
  addToGlobalChi2DerVectorNum (ynum, 100, eps);
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    int iglobal = getGlobalParNum(ilocal);
    double calc = ycalc[iglobal];
    double num = ynum[iglobal];
    cout << "fo: " << getName() << " par " << ilocal << "/" 
         << iglobal << " ("<< getParamName(ilocal)
         << ") calc: " << calc << " - num: " << num << " = " << calc-num
         << endl;
  }
}

void ParticleFitObject::test2ndDerivatives () {
  cout << "ParticleFitObject::test2ndDerivatives, object " << getName() << "\n";
  const int idim=100;
  double *Mnum = new double[idim*idim];
  double *Mcalc = new double[idim*idim];
  for (int i = 0; i < idim*idim; ++i) Mnum[i]=Mcalc[i]=0;
  addToGlobalChi2DerMatrix (Mcalc, idim);
  double eps = 0.0001;
  cout << "eps=" << eps << endl;
  addToGlobalChi2DerMatrixNum (Mnum, idim, eps);
  for (int ilocal1 = 0; ilocal1 < getNPar(); ++ilocal1) {
    int iglobal1 = getGlobalParNum (ilocal1);
    for (int ilocal2 = ilocal1; ilocal2 < getNPar(); ++ilocal2) {
      int iglobal2 = getGlobalParNum (ilocal2);
      double calc = Mcalc[idim*iglobal1 + iglobal2];
      double num = Mnum[idim*iglobal1 + iglobal2];
      cout << "fo: " << getName() << " par " << ilocal1 << "/" 
           << iglobal1 << " ("<< getParamName(ilocal1)
           << "), par " << ilocal2 << "/" 
           << iglobal2 << " ("<< getParamName(ilocal2)
           << ") calc: " << calc << " - num: " << num << " = " << calc-num
           << endl;
    }
  }
  delete[] Mnum;
  delete[] Mcalc;
}

double ParticleFitObject::num1stDerivative (int ilocal, double eps) {
    double save = getParam (ilocal);
    setParam (ilocal, save+eps);
    double v1 = getChi2();
    setParam (ilocal, save-eps);
    double v2 = getChi2();
    double result = (v1-v2)/(2*eps);
    setParam (ilocal, save);
    return result;
}

double ParticleFitObject::num2ndDerivative (int ilocal1, double eps1,
                                            int ilocal2, double eps2) {
  double result;

  if (ilocal1 == ilocal2) {
    double save = getParam (ilocal1);
    double v0 = getChi2();
    setParam (ilocal1, save+eps1);
    double v1 = getChi2();
    setParam (ilocal1, save-eps1);
    double v2 = getChi2();
    result = (v1+v2-2*v0)/(eps1*eps1);
    setParam (ilocal1, save);
  }
  else {
    double save1 = getParam (ilocal1);
    double save2 = getParam (ilocal2);
    setParam (ilocal1, save1+eps1);
    setParam (ilocal2, save2+eps2);
    double v11 = getChi2();
    setParam (ilocal2, save2-eps2);
    double v12 = getChi2();
    setParam (ilocal1, save1-eps1);
    double v22 = getChi2();
    setParam (ilocal2, save2+eps2);
    double v21 = getChi2();
    result = (v11+v22-v12-v21)/(4*eps1*eps2);
    setParam (ilocal1, save1);
    setParam (ilocal2, save2);
  }
  return result;
}
