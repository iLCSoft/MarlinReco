////////////////////////////////////////////////////////////////
// Class NeutrinoFitObject
//
// Author: Jenny Boehme
// Last update: $Date: 2007-10-30 15:51:14 $
//          by: $Author: gaede $
// 
// Description: class for jets with (E, theta, phi) in kinematic fits
//               
////////////////////////////////////////////////////////////////

#include "NeutrinoFitObject.h"
#include <cmath>
#include <cassert>
#include <algorithm>

using std::sqrt;
using std::sin;
using std::cos;

// constructor
NeutrinoFitObject::NeutrinoFitObject(double E, double theta, double phi,  
                           double DE, double Dtheta, double Dphi, 
                           double m) {
  setMass (m);
  setParam (0, E, false);
  setParam (1, theta, false);
  setParam (2, phi, false);
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

//setters
bool NeutrinoFitObject::setParam (int i, double par_, bool measured_ ) {
  if (i < 0 || i > 2) {
    std::cerr << "NeutrinoFitObject::setParam: Illegal i=" << i << std::endl;
    return false;
  }
  measured[i] = measured_;
  return setParam (i, par_);
}  

bool NeutrinoFitObject::setParam (int i, double par_ ) {
  switch (i) {
    // Energy: positive, greater than mass
    case 0: par[0] = (par_ >= mass) ? par_ : mass;
            return (par_ >= mass);
    // theta: between 0 and pi
    case 1: par[1] = (par_ >= 0 && par_ < M_PI) ? 
                      par_ : std::acos (std::cos (par_));
            return true;          
    // phi: between -pi and pi
    case 2: par[2] = (std::abs(par_) <= M_PI) ? 
                     par_ : std::atan2 (std::sin (par_), std::cos (par_));
            return true;          
    default: std::cerr << "NeutrinoFitObject::setParam: Illegal i=" << i << std::endl;
  }
  return false;
} 

// these depend on actual parametrisation!
double NeutrinoFitObject::getPx() const {
  return getP() * cos(par[2]) * sin(par[1]);
}
double NeutrinoFitObject::getPy() const {
  return getP() * sin(par[2]) * sin(par[1]);
}
double NeutrinoFitObject::getPz() const {
  return getP() * cos(par[1]);
}
double NeutrinoFitObject::getE() const {return std::abs(par[0]);}


double NeutrinoFitObject::getDPx(int ilocal) const {
  double result = 0;
  if (ilocal == 0) {                                 // d px / dE
    result = par[0]/getP() * cos(par[2]) * sin(par[1]);
  } 
  else if (ilocal == 1) {                                 // d px / d theta
    result = getP() * cos(par[2]) * cos(par[1]);
  }  
  else if (ilocal == 2) {                                 // d px / d phi
    result =  -getP() * sin(par[2]) * sin(par[1]);
  }  
  return result; 
}

double NeutrinoFitObject::getDPy(int ilocal) const {
  double result = 0;
  if (ilocal == 0) {                                 // d py / dE
    result = par[0]/getP() * sin(par[2]) * sin(par[1]);
  } 
  else if (ilocal == 1) {                                 // d py / d theta
    result =  getP() * sin(par[2]) * cos(par[1]);
  }  
  else if (ilocal == 2) {                                 // d py / d phi
    result =  getP() * cos(par[2]) * sin(par[1]);
  }  
  return result; 
}

double NeutrinoFitObject::getDPz(int ilocal) const {
  double result = 0;
  if (ilocal == 0) {                                 // d pz / dE
    result = par[0]/getP() * cos(par[1]);
  } 
  else if (ilocal == 1) {                                 // d pz / d theta
    result = -getP() * sin(par[1]); 
  }  
  else if (ilocal == 2) {                                 // d pz / d phi
    result = 0; 
  }  
  return result; 
}
double NeutrinoFitObject::getD2Px(int ilocal1, int ilocal2) const {
  double result = 0;
  if (ilocal1 > ilocal2) std::swap (ilocal1, ilocal2);
  assert (ilocal1 >= 0 && ilocal1 <= ilocal2);
  assert (ilocal2 >= 0 && ilocal2 <= NPAR);
  if (ilocal1 == 0) { 
    if (ilocal2 == 0) {                            
      result = -mass/std::pow(par[0]*par[0]-mass*mass,1.5) * cos(par[2]) * sin(par[1]);
    }
    else if (ilocal2 == 1) { 
      result = par[0]/getP() * cos(par[2]) * cos(par[1]);
    }
    else if (ilocal2 == 2) { 
      result = -par[0]/getP() * sin(par[2]) * sin(par[1]);
    }
  } 
  else if (ilocal1 == 1) { 
    if (ilocal2 == 1) {                             
    result = -getP() * cos(par[2]) * sin(par[1]);
    }
    else if (ilocal2 == 2) {
      result = -getP() * sin(par[2]) * cos(par[1]);
    }
  }  
  else if (ilocal1 == 2) {  
    if (ilocal2 == 2) {                                
      result =  -getP() * cos(par[2]) * sin(par[1]);
    }
  }  
  return result;
}
double NeutrinoFitObject::getD2Py(int ilocal1, int ilocal2) const {
  double result = 0;
  if (ilocal1 > ilocal2) std::swap (ilocal1, ilocal2);
  assert (ilocal1 >= 0 && ilocal1 <= ilocal2);
  assert (ilocal2 >= 0 && ilocal2 <= NPAR);
  if (ilocal1 == 0) { 
    if (ilocal2 == 0) {                            
      result = -mass/std::pow(par[0]*par[0]-mass*mass, 1.5) * sin(par[2]) * sin(par[1]);
    }
    else if (ilocal2 == 1) { 
      result = par[0]/getP() * sin(par[2]) * cos(par[1]);
    }
    else if (ilocal2 == 2) { 
      result = par[0]/getP() * cos(par[2]) * sin(par[1]);
    }
  } 
  else if (ilocal1 == 1) { 
    if (ilocal2 == 1) {                             
      result = -getP() * sin(par[2]) * sin(par[1]);
    }
    else if (ilocal2 == 2) {
      result =  getP() * cos(par[2]) * cos(par[1]);
    }
  }  
  else if (ilocal1 == 2) {  
    if (ilocal2 == 2) {                                
      result = -getP() * sin(par[2]) * sin(par[1]);
    }
  }  
  return result;
}
double NeutrinoFitObject::getD2Pz(int ilocal1, int ilocal2) const {
  double result = 0;
  if (ilocal1 > ilocal2) std::swap (ilocal1, ilocal2);
  assert (ilocal1 >= 0 && ilocal1 <= ilocal2);
  assert (ilocal2 >= 0 && ilocal2 <= NPAR);
  if (ilocal1 == 0) { 
    if (ilocal2 == 0) {                            
      result = -mass/std::pow(std::abs(par[0]*par[0]-mass*mass), 1.5) * cos(par[1]);
    }
    else if (ilocal2 == 1) { 
      result = -par[0]/getP() * sin(par[1]);
    }
  } 
  else if (ilocal1 == 1) { 
    if (ilocal2 == 1) {                             
      result = -getP() * cos(par[1]); 
    }
  }  
  return result;
}
double NeutrinoFitObject::getD2E (int ilocal1, int ilocal2) const {
  double result = 0;
  if (ilocal1 == 0 && ilocal2 == 0) result = 1;
  return result;
}

double NeutrinoFitObject::getDE(int ilocal) const {
  double result;
  if (ilocal == 0) {                                 // d E / dE
    result = 1.; 
  }  
  else  {                                 // d E / d theta, phi
    result = 0;
  } 
  return result; 
}
          
void NeutrinoFitObject::addToGlobalDerMatrix (int idim, double c, double *M) const {
  assert (0);

}
