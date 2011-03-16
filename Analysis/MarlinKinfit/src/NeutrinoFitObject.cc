////////////////////////////////////////////////////////////////
// Class NeutrinoFitObject
//
// Author: Jenny Boehme
// Last update: $Date: 2011/03/03 15:03:03 $
//          by: $Author: blist $
// 
// Description: class for neutrinos with (E, theta, phi) in kinematic fits
//               
////////////////////////////////////////////////////////////////

#include "NeutrinoFitObject.h"
#include <cmath>
#include <cassert>
#include <algorithm>

using std::sqrt;
using std::sin;
using std::cos;
using std::cout; 
using std::endl;

// constructor
NeutrinoFitObject::NeutrinoFitObject(double E, double theta, double phi, 
                                     double DE, double Dtheta, double Dphi) {
  setMass (0);
  setParam (0, E, false);
  setParam (1, theta, false);
  setParam (2, phi, false);
  setError (0, DE);
  setError (1, Dtheta);
  setError (2, Dphi);
  invalidateCache();
}

// destructor
NeutrinoFitObject::~NeutrinoFitObject() {}

NeutrinoFitObject *NeutrinoFitObject::copy() const {
  return new NeutrinoFitObject (*this);
}
    
NeutrinoFitObject& NeutrinoFitObject::assign (const BaseFitObject& source) {
  if (const NeutrinoFitObject *psource = dynamic_cast<const NeutrinoFitObject *>(&source)) {
    if (psource != this) *this = *psource;
  }
  else {
    assert (0);
  }
  return *this;
}

const char *NeutrinoFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "E";
    case 1: return "theta";
    case 2: return "phi";
  }
  return "undefined";
}

// needed for constructor!
bool NeutrinoFitObject::setParam (int ilocal, double par_, 
                                    bool measured_, bool fixed_) {
  assert (ilocal >= 0 && ilocal < 3);
  if (measured[ilocal] != measured_ || fixed[ilocal] != fixed_) invalidateCache();
  measured[ilocal] = measured_;
  fixed[ilocal] = fixed_;
  return setParam (ilocal, par_);
}  

bool NeutrinoFitObject::setParam (int i, double par_ ) {
  invalidateCache();
  bool result = (par_-par[i])*(par_-par[i]) > eps2*cov[i][i]; 
  switch (i) {
    // Energy: positive, greater than mass
    case 0: par[0] = (par_ >= 0) ? par_ : 0;
            break;
    // theta: between 0 and pi
    case 1: par[1] = (par_ >= 0 && par_ < M_PI) ? 
                      par_ : std::acos (std::cos (par_));
            break;          
    // phi: any value
    case 2: par[2] = par_;
            break;          
    default: std::cerr << "NeutrinoFitObject::setParam: Illegal i=" << i << std::endl;
  }
  return result;
}  
 
bool NeutrinoFitObject::updateParams (double p[], int idim) {

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
    // cout << "NeutrinoFitObject::updateParams: mirrored E!\n";
    e  = -e;
    th = M_PI-th;
    ph = M_PI+ph;
  }
  if (th<0 || th>M_PI) {
    // cout << "NeutrinoFitObject::updateParams: mirrored theta!\n";
    th = M_PI-th;
    ph = M_PI+ph;
  }
  
  bool result = (e -par[0])*(e -par[0]) > eps2*cov[0][0] ||
                (th-par[1])*(th-par[1]) > eps2*cov[1][1] ||
                (ph-par[2])*(ph-par[2]) > eps2*cov[2][2];
  par[0] = e;
  par[1] = (th >= 0 && th < M_PI) ? 
            th : std::acos (std::cos (th));
  if (std::abs(ph) > M_PI) ph = atan2 (sin(ph), cos (ph));          
  par[2] = ph; 
  p[iE]  = par[0];         
  p[ith] = par[1];         
  p[iph] = par[2];         
  return result;
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
double NeutrinoFitObject::getE() const {
  return par[0];
}

double NeutrinoFitObject::getP() const {
  return par[0];
}

double NeutrinoFitObject::getP2() const {
  return par[0]*par[0];
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
    case 0: return ctheta;
    case 1: return -pt;
    case 2: return 0;
  }
  return 0; 
}

double NeutrinoFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  switch (ilocal) {
    case 0: return 1;
    case 1: return 0;
    case 2: return 0;
  }
  return 0; 
}

// 
// double NeutrinoFitObject::getDPx(int ilocal) const {
//   double result;
//   if (ilocal == 0) {                                 // d px / dE
//     result = par[0]/getP() * cos(par[2]) * sin(par[1]);
//   } 
//   else if (ilocal == 1) {                                 // d px / d theta
//     result = getP() * cos(par[2]) * cos(par[1]);
//   }  
//   else if (ilocal == 2) {                                 // d px / d phi
//     result =  -getP() * sin(par[2]) * sin(par[1]);
//   }  
//   return result; 
// }
// 
// double NeutrinoFitObject::getDPy(int ilocal) const {
//   double result;
//   if (ilocal == 0) {                                 // d py / dE
//     result = par[0]/getP() * sin(par[2]) * sin(par[1]);
//   } 
//   else if (ilocal == 1) {                                 // d py / d theta
//     result =  getP() * sin(par[2]) * cos(par[1]);
//   }  
//   else if (ilocal == 2) {                                 // d py / d phi
//     result =  getP() * cos(par[2]) * sin(par[1]);
//   }  
//   return result; 
// }
// 
// double NeutrinoFitObject::getDPz(int ilocal) const {
//   double result;
//   if (ilocal == 0) {                                 // d pz / dE
//     result = par[0]/getP() * cos(par[1]);
//   } 
//   else if (ilocal == 1) {                                 // d pz / d theta
//     result = -getP() * sin(par[1]); 
//   }  
//   else if (ilocal == 2) {                                 // d pz / d phi
//     result = 0; 
//   }  
//   return result; 
// }
// double NeutrinoFitObject::getD2Px(int ilocal1, int ilocal2) const {
//   double result = 0;
//   if (ilocal1 > ilocal2) std::swap (ilocal1, ilocal2);
//   assert (ilocal1 >= 0 && ilocal1 <= ilocal2);
//   assert (ilocal2 >= 0 && ilocal2 <= NPAR);
//   if (ilocal1 == 0) { 
//     if (ilocal2 == 0) {                            
//       result = -mass/std::pow(par[0]*par[0]-mass*mass,1.5) * cos(par[2]) * sin(par[1]);
//     }
//     else if (ilocal2 == 1) { 
//       result = par[0]/getP() * cos(par[2]) * cos(par[1]);
//     }
//     else if (ilocal2 == 2) { 
//       result = -par[0]/getP() * sin(par[2]) * sin(par[1]);
//     }
//   } 
//   else if (ilocal1 == 1) { 
//     if (ilocal2 == 1) {                             
//     result = -getP() * cos(par[2]) * sin(par[1]);
//     }
//     else if (ilocal2 == 2) {
//       result = -getP() * sin(par[2]) * cos(par[1]);
//     }
//   }  
//   else if (ilocal1 == 2) {  
//     if (ilocal2 == 2) {                                
//       result =  -getP() * cos(par[2]) * sin(par[1]);
//     }
//   }  
//   return result;
// }
// double NeutrinoFitObject::getD2Py(int ilocal1, int ilocal2) const {
//   double result = 0;
//   if (ilocal1 > ilocal2) std::swap (ilocal1, ilocal2);
//   assert (ilocal1 >= 0 && ilocal1 <= ilocal2);
//   assert (ilocal2 >= 0 && ilocal2 <= NPAR);
//   if (ilocal1 == 0) { 
//     if (ilocal2 == 0) {                            
//       result = -mass/std::pow(par[0]*par[0]-mass*mass, 1.5) * sin(par[2]) * sin(par[1]);
//     }
//     else if (ilocal2 == 1) { 
//       result = par[0]/getP() * sin(par[2]) * cos(par[1]);
//     }
//     else if (ilocal2 == 2) { 
//       result = par[0]/getP() * cos(par[2]) * sin(par[1]);
//     }
//   } 
//   else if (ilocal1 == 1) { 
//     if (ilocal2 == 1) {                             
//       result = -getP() * sin(par[2]) * sin(par[1]);
//     }
//     else if (ilocal2 == 2) {
//       result =  getP() * cos(par[2]) * cos(par[1]);
//     }
//   }  
//   else if (ilocal1 == 2) {  
//     if (ilocal2 == 2) {                                
//       result = -getP() * sin(par[2]) * sin(par[1]);
//     }
//   }  
//   return result;
// }
// double NeutrinoFitObject::getD2Pz(int ilocal1, int ilocal2) const {
//   double result = 0;
//   if (ilocal1 > ilocal2) std::swap (ilocal1, ilocal2);
//   assert (ilocal1 >= 0 && ilocal1 <= ilocal2);
//   assert (ilocal2 >= 0 && ilocal2 <= NPAR);
//   if (ilocal1 == 0) { 
//     if (ilocal2 == 0) {                            
//       result = -mass/std::pow(std::abs(par[0]*par[0]-mass*mass), 1.5) * cos(par[1]);
//     }
//     else if (ilocal2 == 1) { 
//       result = -par[0]/getP() * sin(par[1]);
//     }
//   } 
//   else if (ilocal1 == 1) { 
//     if (ilocal2 == 1) {                             
//       result = -getP() * cos(par[1]); 
//     }
//   }  
//   return result;
// }
// double NeutrinoFitObject::getD2E (int ilocal1, int ilocal2) const {
//   double result = 0;
//   if (ilocal1 == 0 && ilocal2 == 0) result = 1;
//   return result;
// }
// 
// double NeutrinoFitObject::getDE(int ilocal) const {
//   double result;
//   if (ilocal == 0) {                                 // d E / dE
//     result = 1.; 
//   }  
//   else  {                                 // d E / d theta, phi
//     result = 0;
//   } 
//   return result; 
// }
          
void   NeutrinoFitObject::addToDerivatives (double der[], int idim, 
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
    der_E     += pzfact*ctheta;
    der_theta -= pzfact*pt;
  }
  
  der[i_E]     += der_E;
  der[i_theta] += der_theta;
  der[i_phi]   += der_phi;
}
    
void   NeutrinoFitObject::addTo2ndDerivatives (double der2[], int idim, 
                                          double efact, double pxfact, 
                                          double pyfact, double pzfact) const {
  int i_E  = globalParNum[0];
  int i_th = globalParNum[1];
  int i_ph = globalParNum[2];
  assert (i_E  >= 0 && i_E  < idim);
  assert (i_th >= 0 && i_th < idim);
  assert (i_ph >= 0 && i_ph < idim);
  
  if (!cachevalid) updateCache();
  // for numerical accuracy, add up derivatives first,
  // then add them to global vector
  double der_Eth  = 0;
  double der_Eph  = 0;
  double der_thth = 0;
  double der_thph = 0;
  double der_phph = 0;
  
  if (pxfact != 0) {
    der_Eth  += pxfact*ctheta*cphi;
    der_Eph  -= pxfact*dpydE;
    der_thth -= pxfact*px;
    der_thph -= pxfact*dpydtheta;
    der_phph -= pxfact*px;
  }
  if (pyfact != 0) {
    der_Eth  += pyfact*ctheta*sphi;
    der_Eph  += pyfact*dpxdE;
    der_thth -= pyfact*py;
    der_thph += pyfact*dpxdtheta;
    der_phph -= pyfact*py;
  }
  if (pzfact != 0) {
    der_Eth  -= pzfact*stheta;
    der_thth -= pzfact*pz;
  }
  
  der2[idim*i_E+i_th]  += der_Eth;
  der2[idim*i_E+i_ph]  += der_Eph;
  der2[idim*i_th+i_E]  += der_Eth;
  der2[idim*i_th+i_th] += der_thth;
  der2[idim*i_th+i_ph] += der_thph;
  der2[idim*i_ph+i_E]  += der_Eph;
  der2[idim*i_ph+i_th] += der_thph;
  der2[idim*i_ph+i_ph] += der_phph;
}
    
void   NeutrinoFitObject::addTo2ndDerivatives (double M[], int idim,  double lambda, double der[]) const {
  if (lambda == 0) return;
  double pxfact = lambda*der[1];
  double pyfact = lambda*der[2];
  double pzfact = lambda*der[3];
  
  int i_E  = globalParNum[0];
  int i_th = globalParNum[1];
  int i_ph = globalParNum[2];
  assert (i_E  >= 0 && i_E  < idim);
  assert (i_th >= 0 && i_th < idim);
  assert (i_ph >= 0 && i_ph < idim);
  
  if (!cachevalid) updateCache();
  // for numerical accuracy, add up derivatives first,
  // then add them to global vector
  double der_Eth  = 0;
  double der_Eph  = 0;
  double der_thth = 0;
  double der_thph = 0;
  double der_phph = 0;
  
  if (pxfact != 0) {
    der_Eth  += pxfact*ctheta*cphi;
    der_Eph  -= pxfact*dpydE;
    der_thth -= pxfact*px;
    der_thph -= pxfact*dpydtheta;
    der_phph -= pxfact*px;
  }
  if (pyfact != 0) {
    der_Eth  += pyfact*ctheta*sphi;
    der_Eph  += pyfact*dpxdE;
    der_thth -= pyfact*py;
    der_thph += pyfact*dpxdtheta;
    der_phph -= pyfact*py;
  }
  if (pzfact != 0) {
    der_Eth  -= pzfact*stheta;
    der_thth -= pzfact*pz;
  }
  
  M[idim*i_E+i_th]  += der_Eth;
  M[idim*i_E+i_ph]  += der_Eph;
  M[idim*i_th+i_E]  += der_Eth;
  M[idim*i_th+i_th] += der_thth;
  M[idim*i_th+i_ph] += der_thph;
  M[idim*i_ph+i_E]  += der_Eph;
  M[idim*i_ph+i_th] += der_thph;
  M[idim*i_ph+i_ph] += der_phph;
}
void   NeutrinoFitObject::addTo1stDerivatives (double M[], int idim, double der[], int kglobal) const {
  assert (kglobal >= 0 && kglobal < idim);
  int i_E  = globalParNum[0];
  int i_th = globalParNum[1];
  int i_ph= globalParNum[2];
  assert (i_E  >= 0 && i_E  < idim);
  assert (i_th >= 0 && i_th < idim);
  assert (i_ph >= 0 && i_ph < idim);
  
  if (!cachevalid) updateCache();
  
  double dE  = der[0] + der[1]*dpxdE     + der[2]*dpydE     + der[3]*ctheta;
  double dth =          der[1]*dpxdtheta + der[2]*dpydtheta - der[3]*pt;
  double dph =        - der[1]*py        + der[2]*px;
  
  M[idim*kglobal + i_E]  += dE;  
  M[idim*kglobal + i_th] += dth; 
  M[idim*kglobal + i_ph] += dph; 
  M[idim*i_E  + kglobal] += dE;  
  M[idim*i_th + kglobal] += dth; 
  M[idim*i_ph + kglobal] += dph; 
}

void NeutrinoFitObject::addToGlobalChi2DerVector (double *y, int idim, 
                                                  double lambda, double der[]) const {
  int i_E  = globalParNum[0];
  int i_th = globalParNum[1];
  int i_ph= globalParNum[2];
  assert (i_E  >= 0 && i_E  < idim);
  assert (i_th >= 0 && i_th < idim);
  assert (i_ph >= 0 && i_ph < idim);
  
  if (!cachevalid) updateCache();
  
  y[i_E]  += lambda*(der[0] + der[1]*dpxdE     + der[2]*dpydE     + der[3]*ctheta);
  y[i_th] += lambda*(         der[1]*dpxdtheta + der[2]*dpydtheta - der[3]*pt);
  y[i_ph] += lambda*(       - der[1]*py        + der[2]*px);
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

  pt = e*stheta;

  px = pt*cphi;
  py = pt*sphi;
  pz = e*ctheta;
  dpxdE = stheta*cphi;
  dpydE = stheta*sphi;
  dpxdtheta = pz*cphi;
  dpydtheta = pz*sphi;
  // dpzdtheta = -p*stheta = -pt
  // dpxdphi  = -pt*sphi  = -py
  // dpydphi  =  pt*cphi  = px

  cachevalid = true;
}

double NeutrinoFitObject::getError2 (double der[]) const {
  if (!cachevalid) updateCache();
  double cov4[4][4]; // covariance of E, px, py, px
  cov4[0][0] =                      cov[0][0];       // E, E
  cov4[0][1] = cov4[1][0] =                          // E, px 
        dpxdE*cov[0][0] + 2*dpxdtheta*cov[0][1] - 2*py*cov[0][2]; 
  cov4[0][2] = cov4[2][0] =                          // E, py
        dpydE*cov[0][0] + 2*dpydtheta*cov[0][1] + 2*px*cov[0][2]; 
  cov4[0][3] = cov4[3][0] =                          // E, pz
        ctheta*cov[0][0] - 2*pt*cov[0][1]; 
  cov4[1][1] =                                       // px, px
       dpxdE*(dpxdE*cov[0][0] + 2*dpxdtheta*cov[0][1] - 2*py*cov[0][2])
     + dpxdtheta*(dpxdtheta*cov[1][1] - 2*py*cov[1][2]) + py*py*cov[2][2]; 
  cov4[1][2] = cov4[2][1] =                         // px, py
       dpxdE*(dpydE*cov[0][0] + 2*dpydtheta*cov[0][1] + 2*px*cov[0][2])
     + dpxdtheta*(dpydtheta*cov[1][1] + 2*px*cov[1][2]) - py*px*cov[2][2]; 
  cov4[1][3] = cov4[3][1] =                         // px, pz
       dpxdE*(ctheta*cov[0][0] - 2*pt*cov[0][1])
     - dpxdtheta*pt*cov[1][1]; 
  cov4[2][2] =                                       // py, py
       dpydE*(dpydE*cov[0][0] + 2*dpydtheta*cov[0][1] + 2*px*cov[0][2])
     + dpydtheta*(dpydtheta*cov[1][1] + 2*px*cov[1][2]) + px*px*cov[2][2]; 
  cov4[2][3] = cov4[3][2] =                          // py, pz
       dpydE*(ctheta*cov[0][0] - 2*pt*cov[0][1])
     - dpydtheta*pt*cov[1][1]; 
  cov4[3][3] =                                      // pz, pz
       ctheta*(ctheta*cov[0][0] - 2*pt*cov[0][1]) + pt*pt*cov[1][1]; 
  return der[0]*(der[0]*cov4[0][0] + 2*der[1]*cov4[0][1] + 2*der[2]*cov4[0][2] + 2*der[3]*cov4[0][3])
                             + der[1]*(der[1]*cov4[1][1] + 2*der[2]*cov4[1][2] + 2*der[3]*cov4[1][3])
                                                   + der[2]*(der[2]*cov4[2][2] + 2*der[3]*cov4[2][3])
                                                                          + der[3]*der[3]*cov4[3][3];
}
