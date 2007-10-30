/*! \file 
 *  \brief Implements class ParticleFitObject
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
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


ParticleFitObject::ParticleFitObject() {
  for (int ilocal = 0; ilocal < NPAR; ++ilocal) globalParNum[ilocal] = -1;
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
  if (!isfinite(par_)) return false;
  assert (ilocal >= 0 && ilocal < NPAR);
  if (par[ilocal] == par_) return true;
  invalidateCache();
  par[ilocal] = par_;
  return true;
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
  int ierr = (idim == 0) ? 0 : dsinv(idim, &covinv[0][0], NPAR);
  if (ierr != 0) {
    std::cerr << "ParticleFitObject::calculateCovInv: Error "
              << ierr << " from dsinv! Object " << getName() << std::endl;
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
      result += covinv[ilocal][jlocal]*(par[jlocal]-measured[jlocal]);
  return 2*result;
}
    
double ParticleFitObject::getD2Chi2DParam2(int ilocal, int jlocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  assert (jlocal >= 0 && jlocal < NPAR);
  if (isParamFixed(ilocal) || !isParamMeasured(ilocal) && 
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

