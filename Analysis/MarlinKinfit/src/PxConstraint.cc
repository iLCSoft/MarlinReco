////////////////////////////////////////////////////////////////
// Class PxConstraint
//
// Author: Jenny Boehme
// Last update: $Date: 2007-10-30 15:51:14 $
//          by: $Author: gaede $
// 
// Description: sum (px) = 0 constraint
//               
////////////////////////////////////////////////////////////////

#include "PxConstraint.h"
#include "ParticleFitObject.h"

#include<iostream>
#include<cassert>

using std::cout;
using std::endl;

// constructor
PxConstraint::PxConstraint () {}

// destructor
PxConstraint::~PxConstraint () {}

// calulate current value of constraint function
double PxConstraint::getValue() const {
  double totpx = 0;
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    totpx += fitobjects[i]->getPx(); 
  }
  return totpx;
}

// calculate vector/array of derivatives of this contraint 
// w.r.t. to ALL parameters of all fitobjects
// here: d sum(px) /d par(i,j) 
//                      = d sum(px) /d px(i) * d px(i) /d par(i, j)
//                                      =  1 * d px(i) /d par(i, j)
void PxConstraint::getDerivatives(int idim, double der[]) const {
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
      if (!fitobjects[i]->isParamFixed(ilocal)) {
        int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
        assert (iglobal >= 0 && iglobal < idim);
        der[iglobal] = fitobjects[i]->getDPx(ilocal);
//       cout << "Px: der[" << iglobal << "] = " << der[iglobal] 
//            << " for jet " << i << " and ilocal = " << ilocal << endl;
      }
    }
  }
}
  
void PxConstraint::add1stDerivativesToMatrix(int idim, double *M) const {
  assert (M);
  int kglobal = getGlobalNum();
  assert (kglobal >= 0 && kglobal < idim);
  
  for (ConstFitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
    const ParticleFitObject *fo = *i;
    assert (fo);
    for (int ilocal = 0; ilocal < fo->getNPar(); ++ilocal) {
      if (!fo->isParamFixed(ilocal)) {
        int iglobal = fo->getGlobalParNum (ilocal);
        assert (iglobal >= 0 && iglobal < idim);
        double d = fo->getDPx(ilocal);
        M[idim*iglobal+kglobal] += d;
        M[idim*kglobal+iglobal] += d;
      }
    }
  }
}
  
void PxConstraint::add2ndDerivativesToMatrix(int idim, double *M, double lambda) const {
  assert (M);
  // Nothing to be done: All second derivatives are 0!
}

