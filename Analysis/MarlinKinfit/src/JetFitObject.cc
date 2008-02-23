/*! \file 
 *  \brief Implements class JetFitObject
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.10  2008/02/07 08:15:25  blist
 * - error calculation of constraints fixed
 * -
 * - Revision 1.9  2008/02/04 17:30:53  blist
 * - NewtonFitter works now!
 * -
 * - Revision 1.8  2008/01/30 21:48:02  blist
 * - Newton Fitter still doesnt work :-(
 * -
 * - Revision 1.7  2008/01/30 09:14:53  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.6  2008/01/29 17:22:59  blist
 * - new addTo2ndDerivatives and setParam
 * -
 * - Revision 1.5  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.4  2007/09/13 08:09:50  blist
 * - Updated 2nd derivatives for px,py,pz,E constraints, improved header documentation
 * -
 *
 */ 

#include "JetFitObject.h"
#include <cmath>
#include <cassert>
#include <iostream>

using std::sqrt;
using std::sin;
using std::cos;
using std::cout; 
using std::endl;

// constructor
JetFitObject::JetFitObject(double E, double theta, double phi,  
                           double DE, double Dtheta, double Dphi, 
                           double m) {
  initCov();                         
  setParam (0, E, true);
  setParam (1, theta, true);
  setParam (2, phi, true);
  setMParam (0, E);
  setMParam (1, theta);
  setMParam (2, phi);
  setError (0, DE);
  setError (1, Dtheta);
  setError (2, Dphi);
  setMass (m);
  invalidateCache();
//   std::cout << "JetFitObject::JetFitObject: E = " << E << std::endl;
//   std::cout << "JetFitObject::JetFitObject: getParam(0) = " << getParam(0) << std::endl;
//   std::cout << "JetFitObject::JetFitObject: " << *this << std::endl;
//   std::cout << "mpar= " << mpar[0] << ", " << mpar[1] << ", " << mpar[2] << std::endl;
}

// destructor
JetFitObject::~JetFitObject() {}

const char *JetFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "E";
    case 1: return "theta";
    case 2: return "phi";
  }
  return "undefined";
}

// needed for constructor!
bool JetFitObject::setParam (int ilocal, double par_, 
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
};  

bool JetFitObject::setParam (int i, double par_ ) {
  invalidateCache();
//   if (i==0) {
//      std::cout << "setParam: par_ = " << par_ << endl;
//      std::cout << "setParam: par[0] = " << par[0] << endl;
//   }   
  bool result = (par_-par[i])*(par_-par[i]) > eps2*cov[i][i]; 
  switch (i) {
    // Energy: positive, greater than mass
    case 0: par[0] = (par_ >= mass) ? par_ : mass;
            //std::cout << "setParam: par[0] = " << par[0] << endl;
            break;
    // theta: between 0 and pi
    case 1: par[1] = (par_ >= 0 && par_ < M_PI) ? 
                      par_ : std::acos (std::cos (par_));
            break;          
    // phi: any value
    case 2: par[2] = par_;
            break;          
    default: std::cerr << "JetFitObject::setParam: Illegal i=" << i << std::endl;
  }
  return result;
}; 
 
bool JetFitObject::updateParams (double p[], int idim) {
  invalidateCache();
  
  int iE  = getGlobalParNum(0);
  int ith = getGlobalParNum(1);
  int iph = getGlobalParNum(2);
  assert (iE  >= 0 && iE  < idim);
  assert (ith >= 0 && ith < idim);
  assert (iph >= 0 && iph < idim);
  
  double e  = p[iE];
  double th = p[ith];
  double ph = p[iph];
  
  if (e<0) {
    // cout << "JetFitObject::updateParams: mirrored E!\n";
    e  = -e;
    th = M_PI-th;
    ph = M_PI+ph;
  }
  if (th<0 || th>M_PI) {
    // cout << "JetFitObject::updateParams: mirrored theta!\n";
    th = M_PI-th;
    ph = M_PI+ph;
  }
  
//   if ((e -par[0])*(e -par[0]) > eps2*cov[0][0])
//     cout << getName() << ": significant change of E: " << par[0] << "->" << e 
//          << ", diff=" << e -par[0] << ", error: " << std::sqrt(cov[0][0]) << endl;
//   if ((th-par[1])*(th-par[1]) > eps2*cov[1][1])
//     cout << getName() << ": significant change of theta: " << par[1] << "->" << th 
//          << ", diff=" << th -par[1]<< ", error: " << std::sqrt(cov[1][1]) << endl;
//   if ((ph-par[2])*(ph-par[2]) > eps2*cov[2][2])
//     cout << getName() << ": significant change of phi: " << par[2] << "->" << ph 
//          << ", diff=" << ph -par[2]<< ", error: " << std::sqrt(cov[2][2]) << endl;
  
  bool result = ((e -par[0])*(e -par[0]) > eps2*cov[0][0]) ||
                ((th-par[1])*(th-par[1]) > eps2*cov[1][1]) ||
                ((ph-par[2])*(ph-par[2]) > eps2*cov[2][2]);
  
  par[0] = (e >= mass) ? e : mass;
  par[1] = (th >= 0 && th < M_PI) ? 
            th : std::acos (std::cos (th));
  double dphi = ph-mpar[2];
  if (std::abs(dphi) > M_PI) dphi = atan2 (sin(dphi), cos (dphi));
  par[2] = mpar[2]+dphi;
  p[iE]  = par[0];         
  p[ith] = par[1];         
  p[iph] = par[2];         
  return result;
}  

// these depend on actual parametrisation!
double JetFitObject::getPx() const {
  if (!cachevalid) updateCache();
  return px;
}
double JetFitObject::getPy() const {
  if (!cachevalid) updateCache();
  return py;
}
double JetFitObject::getPz() const {
  if (!cachevalid) updateCache();
  return pz;
}
double JetFitObject::getE() const {return par[0];}

double JetFitObject::getP() const {
  if (!cachevalid) updateCache();
  return p;
}
double JetFitObject::getP2() const {
  if (!cachevalid) updateCache();
  return p2;
}
double JetFitObject::getPt() const {
  if (!cachevalid) updateCache();
  return pt;
}
double JetFitObject::getPt2() const {
  if (!cachevalid) updateCache();
  return pt*pt;
}

double JetFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpxdE;
    case 1: return dpxdtheta;
    case 2: return -py;
  }
  return 0; 
}

double JetFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpydE;
    case 1: return dpydtheta;
    case 2: return px;
  }
  return 0; 
}

double JetFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpzdE;
    case 1: return -pt;
    case 2: return 0;
  }
  return 0; 
}

double JetFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  switch (ilocal) {
    case 0: return 1;
    case 1: return 0;
    case 2: return 0;
  }
  return 0; 
}
 
void   JetFitObject::addToDerivatives (double der[], int idim, 
                                       double efact, double pxfact, 
                                       double pyfact, double pzfact) const {
  int i_E     = globalParNum[0];
  int i_theta = globalParNum[1];
  int i_phi   = globalParNum[2];
  assert (i_E     >= 0 && i_E     < idim);
  assert (i_theta >= 0 && i_theta < idim);
  assert (i_phi   >= 0 && i_phi   < idim);
  
  if (!cachevalid) updateCache();
  // for numerical accuracy, add up derivatives first,
  // then add them to global vector
  double der_E = efact;
  double der_theta = 0;
  double der_phi = 0;
  
  if (pxfact != 0) {
    der_E     += pxfact*dpxdE;
    der_theta += pxfact*dpydtheta;
    der_phi   -= pxfact*py;
  }
  if (pyfact != 0) {
    der_E     += pyfact*dpydE;
    der_theta += pyfact*dpydtheta;
    der_phi   += pyfact*px;
  }
  if (pzfact != 0) {
    der_E     += pzfact*dpdE*ctheta;
    der_theta -= pzfact*pt;
  }
  
  der[i_E]     += der_E;
  der[i_theta] += der_theta;
  der[i_phi]   += der_phi;
}
    
void   JetFitObject::addTo2ndDerivatives (double der2[], int idim, 
                                          double efact, double pxfact, 
                                          double pyfact, double pzfact) const {
  int i_E  = globalParNum[0];
  int i_th = globalParNum[1];
  int i_ph= globalParNum[2];
  assert (i_E  >= 0 && i_E  < idim);
  assert (i_th >= 0 && i_th < idim);
  assert (i_ph >= 0 && i_ph < idim);
  
  if (!cachevalid) updateCache();
  // for numerical accuracy, add up derivatives first,
  // then add them to global vector
  double der_EE   = 0;
  double der_Eth  = 0;
  double der_Eph  = 0;
  double der_thth = 0;
  double der_thph = 0;
  double der_phph = 0;
  
  double d2pdE2 = (mass != 0) ? -mass*mass/(p*p*p) : 0;
  double d2ptdE2 = d2pdE2*stheta;
  
  if (pxfact != 0) {
    der_EE   += pxfact*d2ptdE2*cphi;
    der_Eth  += pxfact*dpzdE*cphi;
    der_Eph  -= pxfact*dpydE;
    der_thth -= pxfact*px;
    der_thph -= pxfact*dpydtheta;
    der_phph -= pxfact*px;
  }
  if (pyfact != 0) {
    der_EE   += pyfact*d2ptdE2*sphi;
    der_Eth  += pyfact*dpzdE*sphi;
    der_Eph  += pyfact*dpxdE;
    der_thth -= pyfact*py;
    der_thph += pyfact*dpxdtheta;
    der_phph -= pyfact*py;
  }
  if (pzfact != 0) {
    der_EE   += pzfact*d2pdE2*ctheta;
    der_Eth  -= pzfact*dptdE;
    der_thth -= pzfact*pz;
  }
  
  der2[idim*i_E+i_E]   += der_EE;
  der2[idim*i_E+i_th]  += der_Eth;
  der2[idim*i_E+i_ph]  += der_Eph;
  der2[idim*i_th+i_E]  += der_Eth;
  der2[idim*i_th+i_th] += der_thth;
  der2[idim*i_th+i_ph] += der_thph;
  der2[idim*i_ph+i_E]  += der_Eph;
  der2[idim*i_ph+i_th] += der_thph;
  der2[idim*i_ph+i_ph] += der_phph;
}
    
void   JetFitObject::addTo2ndDerivatives (double M[], int idim,  double lambda, double der[]) const {
  addTo2ndDerivatives (M, idim, lambda*der[0], lambda*der[1], lambda*der[2], lambda*der[3]);
}

void   JetFitObject::addTo1stDerivatives (double M[], int idim, double der[], int kglobal) const {
  assert (kglobal >= 0 && kglobal < idim);
  int i_E  = globalParNum[0];
  int i_th = globalParNum[1];
  int i_ph= globalParNum[2];
  assert (i_E  >= 0 && i_E  < idim);
  assert (i_th >= 0 && i_th < idim);
  assert (i_ph >= 0 && i_ph < idim);
  
  if (!cachevalid) updateCache();
  
  double dE  = der[0] + der[1]*dpxdE     + der[2]*dpydE     + der[3]*dpzdE;
  double dth =          der[1]*dpxdtheta + der[2]*dpydtheta - der[3]*pt;
  double dph =        - der[1]*py        + der[2]*px;
  
  M[idim*kglobal + i_E]  += dE;  
  M[idim*kglobal + i_th] += dth; 
  M[idim*kglobal + i_ph] += dph; 
  M[idim*i_E  + kglobal] += dE;  
  M[idim*i_th + kglobal] += dth; 
  M[idim*i_ph + kglobal] += dph; 
}

         
void JetFitObject::initCov() {
  for (int i = 0; i < NPAR; ++i) {
    for (int j = 0; j < NPAR; ++j) {
      cov[i][j] = static_cast<double>(i == j);
    }
  }    
}

void JetFitObject::invalidateCache() const {
  cachevalid = false;
}

void JetFitObject::addToGlobalChi2DerVector (double *y, int idim, 
                                             double lambda, double der[]) const {
  int i_E  = globalParNum[0];
  int i_th = globalParNum[1];
  int i_ph = globalParNum[2];
  assert (i_E  >= 0 && i_E  < idim);
  assert (i_th >= 0 && i_th < idim);
  assert (i_ph >= 0 && i_ph < idim);
  
  if (!cachevalid) updateCache();
  
  y[i_E]  += lambda*(der[0] + der[1]*dpxdE     + der[2]*dpydE     + der[3]*dpzdE);
  y[i_th] += lambda*(         der[1]*dpxdtheta + der[2]*dpydtheta - der[3]*pt);
  y[i_ph] += lambda*(       - der[1]*py        + der[2]*px);
}

void JetFitObject::updateCache() const {
  // std::cout << "JetFitObject::updateCache" << std::endl;
  double e     = par[0];
  double theta = par[1];
  double phi   = par[2];

  ctheta = cos(theta);
  stheta = sin(theta);
  cphi   = cos(phi);
  sphi   = sin(phi);

  p2 = std::abs(e*e-mass*mass);
  p = std::sqrt(p2);
  pt = p*stheta;

  px = pt*cphi;
  py = pt*sphi;
  pz = p*ctheta;
  dpdE = e/p;
  dptdE = dpdE*stheta;
  dpxdE = dptdE*cphi;
  dpydE = dptdE*sphi;
  dpzdE = dpdE*ctheta;
  dpxdtheta = pz*cphi;
  dpydtheta = pz*sphi;
  // dpzdtheta = -p*stheta = -pt
  // dpxdphi  = -pt*sphi  = -py
  // dpydphi  =  pt*cphi  = px
//   d2pdE2 
//   d2ptsE2

  cachevalid = true;
}

double JetFitObject::getError2 (double der[]) const {
  if (!cachevalid) updateCache();
  double cov4[4][4]; // covariance of E, px, py, px
  cov4[0][0] =                      cov[0][0];       // E, E
  cov4[0][1] = cov4[1][0] =                          // E, px 
        dpxdE*cov[0][0] + 2*dpxdtheta*cov[0][1] - 2*py*cov[0][2]; 
  cov4[0][2] = cov4[2][0] =                          // E, py
        dpydE*cov[0][0] + 2*dpydtheta*cov[0][1] + 2*px*cov[0][2]; 
  cov4[0][3] = cov4[3][0] =                          // E, pz
        dpzdE*cov[0][0] - 2*pt*cov[0][1]; 
  cov4[1][1] =                                       // px, px
       dpxdE*(dpxdE*cov[0][0] + 2*dpxdtheta*cov[0][1] - 2*py*cov[0][2])
     + dpxdtheta*(dpxdtheta*cov[1][1] - 2*py*cov[1][2]) + py*py*cov[2][2]; 
  cov4[1][2] = cov4[2][1] =                         // px, py
       dpxdE*(dpydE*cov[0][0] + 2*dpydtheta*cov[0][1] + 2*px*cov[0][2])
     + dpxdtheta*(dpydtheta*cov[1][1] + 2*px*cov[1][2]) - py*px*cov[2][2]; 
  cov4[1][3] = cov4[3][1] =                         // px, pz
       dpxdE*(dpzdE*cov[0][0] - 2*pt*cov[0][1])
     - dpxdtheta*pt*cov[1][1]; 
  cov4[2][2] =                                       // py, py
       dpydE*(dpydE*cov[0][0] + 2*dpydtheta*cov[0][1] + 2*px*cov[0][2])
     + dpydtheta*(dpydtheta*cov[1][1] + 2*px*cov[1][2]) + px*px*cov[2][2]; 
  cov4[2][3] = cov4[3][2] =                          // py, pz
       dpydE*(dpzdE*cov[0][0] - 2*pt*cov[0][1])
     - dpydtheta*pt*cov[1][1]; 
  cov4[3][3] =                                      // pz, pz
       dpzdE*(dpzdE*cov[0][0] - 2*pt*cov[0][1]) + pt*pt*cov[1][1]; 
  return der[0]*(der[0]*cov4[0][0] + 2*der[1]*cov4[0][1] + 2*der[2]*cov4[0][2] + 2*der[3]*cov4[0][3])
                             + der[1]*(der[1]*cov4[1][1] + 2*der[2]*cov4[1][2] + 2*der[3]*cov4[1][3])
                                                   + der[2]*(der[2]*cov4[2][2] + 2*der[3]*cov4[2][3])
                                                                          + der[3]*der[3]*cov4[3][3];
}
