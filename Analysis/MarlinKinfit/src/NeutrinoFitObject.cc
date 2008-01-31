////////////////////////////////////////////////////////////////
// Class NeutrinoFitObject
//
// Author: Jenny Boehme
// Last update: $Date: 2008-01-31 13:01:55 $
//          by: $Author: listj $
// 
// Description: class for neutrinos with (E, theta, phi) in kinematic fits
//              parametrisation, derivatives etc are the same as for jets,
//              but here E, theta and phi are declared "unmeasured"  
//               
////////////////////////////////////////////////////////////////

#include"NeutrinoFitObject.h"
#include <cmath>
#include <cassert>

using std::sqrt;
using std::sin;
using std::cos;

// constructor
NeutrinoFitObject::NeutrinoFitObject(double E, double theta, double phi,  
                           double DE, double Dtheta, double Dphi, 
                           double m) {
  initCov();                         
  setParam (0, E, false);
  setParam (1, theta, false);
  setParam (2, phi, false);
  setMParam (0, E);
  setMParam (1, theta);
  setMParam (2, phi);
  setError (0, DE);
  setError (1, Dtheta);
  setError (2, Dphi);
  setMass (m);
  invalidateCache();
}

// destructor
NeutrinoFitObject::~NeutrinoFitObject() {}

const char *NeutrinoFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "E";
    case 1: return "theta";
    case 2: return "phi";
  }
  return "undefined";
}

// these depend on actual parametrisation!
double NeutrinoFitObject::getPx() const {
  if (!cachevalid) updateCache();
  return px;
}
double NeutrinoFitObject::getPy() const {
  if (!cachevalid) updateCache();
  return py;
}
double NeutrinoFitObject::getPz() const {
  if (!cachevalid) updateCache();
  return pz;
}
double NeutrinoFitObject::getE() const {return par[0];}

double NeutrinoFitObject::getP() const {
  if (!cachevalid) updateCache();
  return p;
}
double NeutrinoFitObject::getP2() const {
  if (!cachevalid) updateCache();
  return p2;
}
double NeutrinoFitObject::getPt() const {
  if (!cachevalid) updateCache();
  return pt;
}
double NeutrinoFitObject::getPt2() const {
  if (!cachevalid) updateCache();
  return pt*pt;
}

double NeutrinoFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpxdE;
    case 1: return dpxdtheta;
    case 2: return -py;
  }
  return 0; 
}

double NeutrinoFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpydE;
    case 1: return dpydtheta;
    case 2: return px;
  }
  return 0; 
}

double NeutrinoFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpzdE;
    case 1: return -pt;
    case 2: return 0;
  }
  return 0; 
}

double NeutrinoFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return 1;
    case 1: return 0;
    case 2: return 0;
  }
  return 0; 
}
    
void   NeutrinoFitObject::addToDerivatives (double der[],
                                       int idim,
                                       double pxfact, 
                                       double pyfact, 
                                       double pzfact, 
                                       double efact) const {
  int i_E     = globalParNum[0];
  int i_theta = globalParNum[1];
  int i_phi   = globalParNum[2];
  assert (i_E     >= 0 && i_E     < idim);
  assert (i_theta >= 0 && i_theta < idim);
  assert (i_phi   >= 0 && i_phi   < idim);
  
  if (!cachevalid) updateCache();
  // for numerical accuracy, add up derivatives first,
  // then add them to global vector
  double der_E = 0;
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
  der_E     += efact;
  
  der[i_E]     += der_E;
  der[i_theta] += der_theta;
  der[i_phi]   += der_phi;
}
    
void   NeutrinoFitObject::addTo2ndDerivatives (double der2[],
                                          int idim,
                                          double pxfact, 
                                          double pyfact, 
                                          double pzfact, 
                                          double efact) const {
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


          
void NeutrinoFitObject::addToGlobalDerMatrix (int idim, double c, double *M) const {
  // add second derivatives to global matrix
  if (c == 0) return;
  assert (0);
  
  
}

void NeutrinoFitObject::initCov() {
  for (int i = 0; i < NPAR; ++i) {
    for (int j = 0; j < NPAR; ++j) {
      cov[i][j] = static_cast<double>(i == j);
    }
  }    
}

void NeutrinoFitObject::invalidateCache() const {
  cachevalid = false;
}

void NeutrinoFitObject::updateCache() const {
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
//   d2pdE2 
//   d2ptsE2

  cachevalid = true;
}
