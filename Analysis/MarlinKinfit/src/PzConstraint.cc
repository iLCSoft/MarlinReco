/*! \file 
 *  \brief Implements class PzConstraint
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.3  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 * - Revision 1.2  2007/09/13 08:09:51  blist
 * - Updated 2nd derivatives for px,py,pz,E constraints, improved header documentation
 * -
 *
 */ 

#include "PzConstraint.h"
#include "ParticleFitObject.h"

#include<iostream>
#include<cassert>

using std::cout;
using std::endl;

// constructor
PzConstraint::PzConstraint () {}

// destructor
PzConstraint::~PzConstraint () {}

// calulate current value of constraint function
double PzConstraint::getValue() const {
  double totpz = 0;
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    totpz += fitobjects[i]->getPz(); 
  }
  return totpz;
}

// calculate vector/array of derivatives of this contraint 
// w.r.t. to ALL parameters of all fitobjects
// here: d sum(px) /d par(i,j) 
//                      = d sum(px) /d px(i) * d px(i) /d par(i, j)
//                                      =  1 * d px(i) /d par(i, j)
void PzConstraint::getDerivatives(int idim, double der[]) const {
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
      if (!fitobjects[i]->isParamFixed(ilocal)) {
        int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
        assert (iglobal >= 0 && iglobal < idim);
        der[iglobal] = fitobjects[i]->getDPz(ilocal);
//       cout << "Pz: der[" << iglobal << "] = " << der[iglobal] 
//            << " for jet " << i << " and ilocal = " << ilocal << endl;
      }
    }
  }
}
  
void PzConstraint::add1stDerivativesToMatrix(int idim, double *M) const {
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
        double d = fo->getDPz(ilocal);
        M[idim*iglobal+kglobal] += d;
        M[idim*kglobal+iglobal] += d;
      }
    }
  }
}
  
void PzConstraint::add2ndDerivativesToMatrix(int idim, double *M, double lambda) const {
  assert (M);
  for (ConstFitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
    const ParticleFitObject *fo = *i;
    assert (fo);
    fo->addTo2ndDerivatives(M, idim, 0, 0, lambda, 0);
  }
}

