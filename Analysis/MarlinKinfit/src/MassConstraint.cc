////////////////////////////////////////////////////////////////
// Class MassConstraint
//
// Author: Jenny Boehme
// Last update: $Date: 2007-10-30 15:51:14 $
//          by: $Author: gaede $
// 
// Description: M(p1,p2)=M_W constraint
//               
////////////////////////////////////////////////////////////////

#include "MassConstraint.h"
#include "ParticleFitObject.h"

#include<iostream>
#include<cmath>
#include<cassert>

using std::cerr;
using std::cout;
using std::endl;

// constructor
MassConstraint::MassConstraint (double mass_) : mass(mass_) {}

// destructor
MassConstraint::~MassConstraint () {}

// calulate current value of constraint function
double MassConstraint::getValue() const {
  double totE[2] = {0,0};
  double totpx[2] = {0,0}; 
  double totpy[2] = {0,0}; 
  double totpz[2] = {0,0}; 
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    int index = (flags[i] == 1) ? 0 : 1; // default is 1, but 2 may indicate fitobjects for a second W -> equal mass constraint!
    totE[index] += fitobjects[i]->getE(); 
    totpx[index] += fitobjects[i]->getPx(); 
    totpy[index] += fitobjects[i]->getPy(); 
    totpz[index] += fitobjects[i]->getPz(); 
  }
  double result = -mass;
  result += std::sqrt(std::abs(totE[0]*totE[0]-totpx[0]*totpx[0]-totpy[0]*totpy[0]-totpz[0]*totpz[0]));
  result -= std::sqrt(std::abs(totE[1]*totE[1]-totpx[1]*totpx[1]-totpy[1]*totpy[1]-totpz[1]*totpz[1]));
  return result;
}

// calculate vector/array of derivatives of this contraint 
// w.r.t. to ALL parameters of all fitobjects
// here: d M /d par(j) 
//          = d M /d p(i) * d p(i) /d par(j)
//          =  +-1/M * p(i) * d p(i) /d par(j)
void MassConstraint::getDerivatives(int idim, double der[]) const {
  double totE[2] = {0,0};
  double totpx[2] = {0,0}; 
  double totpy[2] = {0,0}; 
  double totpz[2] = {0,0}; 
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    int index = (flags[i]==1) ? 0 : 1; // default is 1, but 2 may indicate fitobjects for a second W -> equal mass constraint!
    totE[index]  += fitobjects[i]->getE(); 
    totpx[index] += fitobjects[i]->getPx(); 
    totpy[index] += fitobjects[i]->getPy(); 
    totpz[index] += fitobjects[i]->getPz(); 
  }
  double m2[2]; 
  double m_inv[2] = {0,0}; 
  for (int index = 0; index < 2; ++index) {
    m2[index] = totE[index]*totE[index] - totpx[index]*totpx[index]
                - totpy[index]*totpy[index] - totpz[index]*totpz[index];
    if (m2[index] < 0) cerr << "MassConstraint::getDerivatives: m2<0!" << endl;
    if (m2[index] == 0) cerr << "MassConstraint::getDerivatives: m2==0!" << endl;
    else m_inv[index] = 1/std::sqrt (std::abs(m2[index]));
  }
  
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    int index = (flags[i]==1) ? 0 : 1; // default is 1, but 2 may indicate fitobjects for a second W -> equal mass constraint!
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
      if (!fitobjects[i]->isParamFixed(ilocal)) {
        int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
        assert (iglobal >= 0 && iglobal < idim);
        der[iglobal] =   totE[index]  * fitobjects[i]->getDE (ilocal)
                       - totpx[index] * fitobjects[i]->getDPx (ilocal)
                       - totpy[index] * fitobjects[i]->getDPy (ilocal)
                       - totpz[index] * fitobjects[i]->getDPz (ilocal);
        der[iglobal] *= m_inv[index];
        if (index == 1) der[iglobal] *= -1.;
      }
    }
  }
}
  
double MassConstraint::getMass (int flag) {
  double totE = 0;
  double totpx = 0; 
  double totpy = 0; 
  double totpz = 0; 
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    if (flags[i] == flag) {
      totE += fitobjects[i]->getE(); 
      totpx += fitobjects[i]->getPx(); 
      totpy += fitobjects[i]->getPy(); 
      totpz += fitobjects[i]->getPz(); 
    }
  }
  return std::sqrt(std::abs(totE*totE-totpx*totpx-totpy*totpy-totpz*totpz));
}

void MassConstraint::add1stDerivativesToMatrix(int idim, double *M) const {
  assert (false);
}

void MassConstraint::add2ndDerivativesToMatrix(int idim, double *M, double lambda) const {
  assert (false);
}


