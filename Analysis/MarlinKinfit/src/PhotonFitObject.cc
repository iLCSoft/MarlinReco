/*! \file 
 *  \brief Implements class PhotonFitObject
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.6  2009/02/23 12:04:05  mbeckman
 * - - PhotonFitObject:     bug fix (1/0), removed dispensable variables
 * - - PhotonFitObjectPxyg: bug fixes (1/0, order of computing variables), modified parametrization
 * - - JetFitObject:        added start parameter check (inf, nan)
 * -
 * - Revision 1.5  2009/02/18 11:56:22  mbeckman
 * - PhotonFitObject*.cc: documentation, debug output
 * - NewtonFitterGSL.cc:  bug fix (Lagrange multipliers not initialized), debug output
 * - JetFitObject.cc:     bug fix: division by 0, if energy <= mass
 * -
 *
 */ 

#include "PhotonFitObject.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "marlin/Processor.h"

using std::sqrt;
using std::cout; 
using std::endl;

// constructor
PhotonFitObject::PhotonFitObject(double px, double py, double pz,  
                           double Dpz) {
  initCov();                         
  setParam (0, px, true, true);
  setParam (1, py, true, true);
  setParam (2, pz, true);
  setMParam (0, px);
  setMParam (1, py);
  setMParam (2, pz);
  setError (2, Dpz);
  setMass (0.);
  invalidateCache();
}

// destructor
PhotonFitObject::~PhotonFitObject() {}

const char *PhotonFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "P_x";
    case 1: return "P_y";
    case 2: return "P_z";
  }
  return "undefined";
}

// needed for constructor!
bool PhotonFitObject::setParam (int ilocal, double par_, 
                                    bool measured_, bool fixed_) {
  assert (ilocal >= 0 && ilocal < 3);
  if (measured[ilocal] != measured_ || fixed[ilocal] != fixed_) invalidateCache();
  measured[ilocal] = measured_;
  fixed[ilocal] = fixed_;
// this doesn't work inn constructor, since par[i], mass etc are not initialized yet!!!!  
//  return setParam (ilocal, par_);
// old version of bool ParticleFitObject::setParam (int ilocal, double par_ )
//  if (!isfinite(par_)) return false;  doesn't exist anymore?
//  assert (ilocal >= 0 && ilocal < NPAR);   done before
//  if (par[ilocal] == par_) return true;    
//  invalidateCache();                     done in constructor anyhow
  par[ilocal] = par_;
  return true;
} 

bool PhotonFitObject::setParam (int i, double par_ ) {
  invalidateCache();
//   if (i==0) {
//      std::cout << "setParam: par_ = " << par_ << endl;
//      std::cout << "setParam: par[0] = " << par[0] << endl;
//   }   
  bool result = (par_-par[i])*(par_-par[i]) > eps2*cov[i][i]; 
  switch (i) {
    // p_x
    case 0: par[0] = par_;
            //std::cout << "setParam: par[0] = " << par[0] << endl;
            break;
    // p_y
    case 1: par[1] = par_;
            break;          
    // p_z
    case 2: par[2] = par_;
            break;          
    default: std::cerr << "PhotonFitObject::setParam: Illegal i=" << i << std::endl;
  }
  return result;
}
 
bool PhotonFitObject::updateParams (double p[], int idim) {
  invalidateCache();
  
  int i2 = getGlobalParNum(2);
// std::cout << "updateParams: i2 = " << i2 << "\n";
// std::cout << "updateParams: idim = " << idim << "\n";
  assert (i2 >= 0 && i2 < idim);
  
  double p2 = p[i2];
// std::cout << "updateParams: p2 = " << p[i2] << "   par[2] = " << par[2] << "\n";
  
  bool result = ((p2-par[2])*(p2-par[2]) > eps2*cov[2][2]);

  par[2] = p2;
  p[i2] = par[2];         
  return result;
}  

// these depend on actual parametrisation!
double PhotonFitObject::getPx() const {return par[0];}

double PhotonFitObject::getPy() const {return par[1];}

double PhotonFitObject::getPz() const {return par[2];}

double PhotonFitObject::getE() const {
  if (!cachevalid) updateCache();
  return p;
}
double PhotonFitObject::getP() const {
  if (!cachevalid) updateCache();
  return p;
}
double PhotonFitObject::getP2() const {
  if (!cachevalid) updateCache();
  return p2;
}
double PhotonFitObject::getPt() const {
  if (!cachevalid) updateCache();
  return std::sqrt(pt2);
}
double PhotonFitObject::getPt2() const {
  if (!cachevalid) updateCache();
  return pt2;
}

double PhotonFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
//   if (!cachevalid) updateCache();
  switch (ilocal) {
//     case 0: return dpx0;
//     case 1: return dpx1;
//     case 2: return dpx2;
    case 0: return 1.;
    case 1: return 0.;
    case 2: return 0.;
  }
  return 0; 
}

double PhotonFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
//   if (!cachevalid) updateCache();
  switch (ilocal) {
//     case 0: return dpy0;
//     case 1: return dpy1;
//     case 2: return dpy2;
    case 0: return 0.;
    case 1: return 1.;
    case 2: return 0.;
  }
  return 0; 
}

double PhotonFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
//   if (!cachevalid) updateCache();
  switch (ilocal) {
//     case 0: return dpz0;
//     case 1: return dpz1;
//     case 2: return dpz2;
    case 0: return 0.;
    case 1: return 0.;
    case 2: return 1.;
  }
  return 0; 
}

double PhotonFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dE0;
    case 1: return dE1;
    case 2: return dE2;
  }
  return 0; 
}
 
void   PhotonFitObject::addToDerivatives (double der[], int idim, 
                                       double efact, double pxfact, 
                                       double pyfact, double pzfact) const {
  int i0 = globalParNum[0];
  int i1 = globalParNum[1];
  int i2 = globalParNum[2];
  assert (i0 >= 0 && i0 < idim);
  assert (i1 >= 0 && i1 < idim);
  assert (i2 >= 0 && i2 < idim);
  
  if (!cachevalid) updateCache();
  // for numerical accuracy, add up derivatives first,
  // then add them to global vector
//   double der0 = efact*dE0;
//   double der1 = efact*dE1;
//   double der2 = efact*dE2;
  
//   if (pxfact != 0) {
//     der0 += pxfact*dpx0;
//     der1 += pxfact*dpx1;
//     der2 += pxfact*dpx2;
//   }
//   if (pyfact != 0) {
//     der0 += pyfact*dpy0;
//     der1 += pyfact*dpy1;
//     der2 += pyfact*dpy2;
//   }
//   if (pzfact != 0) {
//     der0 += pzfact*dpz0;
//     der1 += pzfact*dpz1;
//     der2 += pzfact*dpz2;
//   }

//   der[i0] += der0;
//   der[i1] += der1;
//   der[i2] += der2;
  der[i0] += (efact*dE0+pxfact);
  der[i1] += (efact*dE1+pyfact);
  der[i2] += (efact*dE2+pzfact);
}
    
void   PhotonFitObject::addTo2ndDerivatives (double der2[], int idim, 
                                          double efact, double pxfact, 
                                          double pyfact, double pzfact) const {
  int i2 = globalParNum[2];
  assert (i2 >= 0 && i2 < idim);
  
  if (!cachevalid) updateCache();
  // if p==0, 2nd derivates are zero (catch up 1/0)
  if (p==0 || pt2==0) return;
  
  double der22 = pt2/p/p/p*efact;   // NOT using /p2/p to avoid numerical problems
  
  der2[idim*i2+i2] += der22;
}
    
void   PhotonFitObject::addTo2ndDerivatives (double M[], int idim,  double lambda, double der[]) const {
  addTo2ndDerivatives (M, idim, lambda*der[0], lambda*der[1], lambda*der[2], lambda*der[3]);
}

void   PhotonFitObject::addTo1stDerivatives (double M[], int idim, double der[], int kglobal) const {
  assert (kglobal >= 0 && kglobal < idim);
  int i2 = globalParNum[2];
  assert (i2 >= 0 && i2 < idim);
  
  if (!cachevalid) updateCache();
  
//   double d2 = der[0]*dE2 + der[1]*dpx2 + der[2]*dpy2 + der[3]*dpz2;
  double d2 = der[0]*dE2 + der[3];
  
  M[idim*kglobal + i2] += d2; 
  M[idim*i2 + kglobal] += d2; 
}

         
void PhotonFitObject::initCov() {
  for (int i = 0; i < NPAR; ++i) {
    for (int j = 0; j < NPAR; ++j) {
      cov[i][j] = static_cast<double>(i == j);
    }
  }    
}

void PhotonFitObject::invalidateCache() const {
  cachevalid = false;
}

void PhotonFitObject::addToGlobalChi2DerVector (double *y, int idim, 
                                             double lambda, double der[]) const {
  int i2 = globalParNum[2];
  assert (i2 >= 0 && i2 < idim);
  
  if (!cachevalid) updateCache();
  
//   double d2 = der[0]*dE2 + der[1]*dpx2 + der[2]*dpy2 + der[3]*dpz2;

//   y[i2] += lambda*d2;
  y[i2] += lambda*(der[0]*dE2 + der[3]);
}

void PhotonFitObject::updateCache() const {
  // std::cout << "PhotonFitObject::updateCache" << std::endl;
  double px = par[0];
  double py = par[1];
  double pz = par[2];

  pt2 = px*px+py*py;
  p2  = pt2+pz*pz;
  p   = std::sqrt(p2);
  
//   dpx0 = 1.;
//   dpx1 = 0.;
//   dpx2 = 0.;
//   dpy0 = 0.;
//   dpy1 = 1.;
//   dpy2 = 0.;
//   dpz0 = 0.;
//   dpz1 = 0.;
//   dpz2 = 1.;
  // if p==0, derivates are zero (catch up division by zero)
  if(p){
    dE0  = px/p;
    dE1  = py/p;
    dE2  = pz/p;
  }
  else{
    dE0  = 0.;
    dE1  = 0.;
    dE2  = 0.;
  }

  cachevalid = true;
}

double PhotonFitObject::getError2 (double der[]) const {
  if (!cachevalid) updateCache();
  double cov4[4][4]; // covariance of E, px, py, px
//   cov4[0][0] =                                       // E, E
//         dE0*(dE0*cov[0][0]+2*dE1*cov[0][1]+2*dE2*cov[0][2])+dE1*(dE1*cov[1][1]+2*dE2*cov[1][2])+dE2*dE2*cov[2][2];
//   cov4[0][1] = cov4[1][0] =                          // E, px 
//         dE0*dpx0*cov[0][0]+(dE0*dpx1+dE1*dpx0)*cov[0][1]+(dE0*dpx2+dE2*dpx0)*cov[0][2]
//         +dE1*dpx1*cov[1][1]+(dE1*dpx2+dE2*dpx1)*cov[1][2]+dE2*dpx2*cov[2][2];
//   cov4[0][2] = cov4[2][0] =                          // E, py
//         dE0*dpy0*cov[0][0]+(dE0*dpy1+dE1*dpy0)*cov[0][1]+(dE0*dpy2+dE2*dpy0)*cov[0][2]
//         +dE1*dpy1*cov[1][1]+(dE1*dpy2+dE2*dpy1)*cov[1][2]+dE2*dpy2*cov[2][2];
//   cov4[0][3] = cov4[3][0] =                          // E, pz
//         dE0*dpz0*cov[0][0]+(dE0*dpz1+dE1*dpz0)*cov[0][1]+(dE0*dpz2+dE2*dpz0)*cov[0][2]
//         +dE1*dpz1*cov[1][1]+(dE1*dpz2+dE2*dpz1)*cov[1][2]+dE2*dpz2*cov[2][2];
//   cov4[1][1] =                                       // px, px
//         dpx0*(dpx0*cov[0][0]+2*dpx1*cov[0][1]+2*dpx2*cov[0][2])+dpx1*(dpx1*cov[1][1]+2*dpx2*cov[1][2])+dpx2*dpx2*cov[2][2];
//   cov4[1][2] = cov4[2][1] =                          // px, py
//         dpx0*dpy0*cov[0][0]+(dpx0*dpy1+dpx1*dpy0)*cov[0][1]+(dpx0*dpy2+dpx2*dpy0)*cov[0][2]
//         +dpx1*dpy1*cov[1][1]+(dpx1*dpy2+dpx2*dpy1)*cov[1][2]+dpx2*dpy2*cov[2][2];
//   cov4[1][3] = cov4[3][1] =                          // px, pz
//         dpx0*dpz0*cov[0][0]+(dpx0*dpz1+dpx1*dpz0)*cov[0][1]+(dpx0*dpz2+dpx2*dpz0)*cov[0][2]
//         +dpx1*dpz1*cov[1][1]+(dpx1*dpz2+dpx2*dpz1)*cov[1][2]+dpx2*dpz2*cov[2][2];
//   cov4[2][2] =                                       // py, py
//         dpy0*(dpy0*cov[0][0]+2*dpy1*cov[0][1]+2*dpy2*cov[0][2])+dpy1*(dpy1*cov[1][1]+2*dpy2*cov[1][2])+dpy2*dpy2*cov[2][2];
//   cov4[2][3] = cov4[3][2] =                          // py, pz
//         dpy0*dpz0*cov[0][0]+(dpy0*dpz1+dpy1*dpz0)*cov[0][1]+(dpy0*dpz2+dpy2*dpz0)*cov[0][2]
//         +dpy1*dpz1*cov[1][1]+(dpy1*dpz2+dpy2*dpz1)*cov[1][2]+dpy2*dpz2*cov[2][2];
//   cov4[3][3] =                                       // pz, pz
//         dpz0*(dpz0*cov[0][0]+2*dpz1*cov[0][1]+2*dpz2*cov[0][2])+dpz1*(dpz1*cov[1][1]+2*dpz2*cov[1][2])+dpz2*dpz2*cov[2][2];
  cov4[0][0] =                                       // E, E
        dE0*(dE0*cov[0][0]+2*dE1*cov[0][1]+2*dE2*cov[0][2])+dE1*(dE1*cov[1][1]+2*dE2*cov[1][2])+dE2*dE2*cov[2][2];
  cov4[0][1] = cov4[1][0] =                          // E, px 
        dE0*cov[0][0]+dE1*cov[0][1]+dE2*cov[0][2];
  cov4[0][2] = cov4[2][0] =                          // E, py
        dE0*cov[0][1]+dE1*cov[1][1]+dE2*cov[1][2];
  cov4[0][3] = cov4[3][0] =                          // E, pz
        dE0*cov[0][2]+dE1*cov[1][2]+dE2*cov[2][2];
  cov4[1][1] =               cov[0][0];              // px, px
  cov4[1][2] = cov4[2][1] =  cov[0][1];              // px, py
  cov4[1][3] = cov4[3][1] =  cov[0][2];              // px, pz
  cov4[2][2] =               cov[1][1];              // py, py
  cov4[2][3] = cov4[3][2] =  cov[1][2];              // py, pz
  cov4[3][3] =               cov[2][2];              // pz, pz
  
  return der[0]*(der[0]*cov4[0][0] + 2*der[1]*cov4[0][1] + 2*der[2]*cov4[0][2] + 2*der[3]*cov4[0][3])
                             + der[1]*(der[1]*cov4[1][1] + 2*der[2]*cov4[1][2] + 2*der[3]*cov4[1][3])
                                                   + der[2]*(der[2]*cov4[2][2] + 2*der[3]*cov4[2][3])
                                                                          + der[3]*der[3]*cov4[3][3];
}
