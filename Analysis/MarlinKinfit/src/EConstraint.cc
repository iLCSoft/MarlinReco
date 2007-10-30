/*! \file 
 *  \brief Implements class EConstraint
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.2  2007/09/13 08:09:50  blist
 * - Updated 2nd derivatives for px,py,pz,E constraints, improved header documentation
 * -
 *
 */ 

#include "EConstraint.h"
#include "ParticleFitObject.h"

#include<iostream>
#include<cassert>

using std::cout;
using std::endl;

// constructor
EConstraint::EConstraint (double ecm_) : ecm(ecm_) {}

// destructor
EConstraint::~EConstraint () {}

// calulate current value of constraint function
double EConstraint::getValue() const {
  double tote = -getEcm();
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    tote += fitobjects[i]->getE(); 
  }
  return tote;
}

// calculate vector/array of derivatives of this contraint 
// w.r.t. to ALL parameters of all fitobjects
// here: d sum(E) /d par(i,j) 
//                      = d sum(E) /d E(i) * d E(i) /d par(i, j)
//                                      =  1 * d E(i) /d par(i, j)
void EConstraint::getDerivatives(int idim, double der[]) const {
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
      if (!fitobjects[i]->isParamFixed(ilocal)) {
        int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
        assert (iglobal >= 0 && iglobal < idim);
        der[iglobal] = fitobjects[i]->getDE(ilocal);
//       cout << "E: der[" << iglobal << "] = " << der[iglobal] 
//            << " for jet " << i << " and ilocal = " << ilocal << endl;
      }
    }
  }
}
  
void EConstraint::add1stDerivativesToMatrix(int idim, double *M) const {
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
        double d = fo->getDE(ilocal);
        M[idim*iglobal+kglobal] += d;
        M[idim*kglobal+iglobal] += d;
      }
    }
  }
}
  
void EConstraint::add2ndDerivativesToMatrix(int idim, double *M, double lambda) const {
  assert (M);
  for (ConstFitObjectIterator i = fitobjects.begin(); i != fitobjects.end(); ++i) {
    const ParticleFitObject *fo = *i;
    assert (fo);
    fo->addTo2ndDerivatives(M, idim, 0, 0, 0, lambda);
  }
}

