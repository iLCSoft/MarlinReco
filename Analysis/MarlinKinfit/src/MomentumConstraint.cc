/*! \file 
 *  \brief Implements class MomentumConstraint
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.1  2008/02/18 09:59:35  blist
 * - MomentumConstraint and SoftGaussMomentumCOnstraint added; PConstraint is obsolete
 * -
 * -
 *
 */ 

#include "MomentumConstraint.h"
#include "ParticleFitObject.h"

#include<iostream>
#include<cassert>

using std::cout;
using std::endl;

MomentumConstraint::MomentumConstraint (double efact_, double pxfact_, double pyfact_, 
                                        double pzfact_, double value_) 
: efact (efact_),
  pxfact (pxfact_),
  pyfact (pyfact_),
  pzfact (pzfact_),
  value (value_),
  cachevalid(false)
{}

// destructor
MomentumConstraint::~MomentumConstraint () {
  //std::cout << "destroying MomentumConstraint" << std::endl;
}

// calculate current value of constraint function
double MomentumConstraint::getValue() const {
  double totpx = 0;
  double totpy = 0;
  double totpz = 0;
  double totE = 0;
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    if (pxfact != 0) totpx += fitobjects[i]->getPx(); 
    if (pyfact != 0) totpy += fitobjects[i]->getPy(); 
    if (pzfact != 0) totpz += fitobjects[i]->getPz(); 
    if (efact  != 0) totE  += fitobjects[i]->getE(); 
  }
  return pxfact*totpx + pyfact*totpy + pzfact*totpz + efact*totE - value;
}

// calculate vector/array of derivatives of this contraint 
// w.r.t. to ALL parameters of all fitobjects
// here: d sum(px) /d par(i,j) 
//                      = d sum(px) /d px(i) * d px(i) /d par(i, j)
//                                      =  1 * d px(i) /d par(i, j)
void MomentumConstraint::getDerivatives(int idim, double der[]) const {
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
      if (!fitobjects[i]->isParamFixed(ilocal)) {
        int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
        assert (iglobal >= 0 && iglobal < idim);
        double d = 0;
        if (pxfact != 0) d += pxfact*fitobjects[i]->getDPx (ilocal);
        if (pyfact != 0) d += pyfact*fitobjects[i]->getDPy (ilocal);
        if (pzfact != 0) d += pzfact*fitobjects[i]->getDPz (ilocal);
        if (efact  != 0) d +=  efact*fitobjects[i]->getDE  (ilocal);
        der[iglobal] = d;
      }
    }
  }
}
  
    
void MomentumConstraint::addToGlobalDerMatrix (double lambda, int idim, double *M) const {

  assert (0);
  // Add lambda*d^2 g / d x_i dx_j to global matrix
  
  if (lambda == 0) return;
  
  // d^2 g / (dx_i dx_j) = 
  //   = sum_k,l d^2 g/(dpx_k dpx_l) * dpx_k/dx_i dpx_l/dx_j
  //     + sum_k dg/dpx_k * d^2 px_k/(dx_i dx_j)
  //   = sum_k,l      1              * dpx_k/dx_i dpx_l/dx_j
  
  // assume here that different 4-vectors always depend on 
  // different parameters!
  
  if (!cachevalid) updateCache();
  
  int *globalParNum = new int[nparams];
  double *der = new double[nparams];
  
//   ipar = 0;
//   for (int i = 0; i < fitobjects.size(); i++) {
//     for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
//       int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
//       if (iglobal >= 0) {
//         assert (ipar < nparams);
//         globalParNum[ipar] = iglobal;
//         der[ipar] = fitobjects[i]->getDPx (ilocal);
//         ipar++;
//       }
//   }
  
  for (int ipar = 0; ipar < nparams; ipar++) {
    int iglobal = globalParNum[ipar];
    double der_i = der[ipar];
    for (int jpar = ipar; jpar < nparams; jpar++) {
      int jglobal = globalParNum[ipar];
      double der_j = der[jpar];
      double l_der_ij = lambda*der_i*der_j;
      M[idim*iglobal+jglobal] += l_der_ij;
      if (ipar != jpar) M[idim*jglobal+iglobal] += l_der_ij;
    }
  }       
}

void MomentumConstraint::invalidateCache() const {
  cachevalid = false;
}

void MomentumConstraint::updateCache() const {
  nparams = 0;
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
      int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
      if (!fitobjects[i]->isParamFixed(ilocal)) {
        assert (iglobal >= 0);
        nparams++;
      }
    }
  }
  cachevalid = true;
}
  
bool MomentumConstraint::secondDerivatives (int i, int j, double *derivatives) const {
  return false;
}  
  
bool MomentumConstraint::firstDerivatives (int i, double *derivatives) const {
  derivatives[0] = efact;
  derivatives[1] = pxfact;
  derivatives[2] = pyfact;
  derivatives[3] = pzfact;
  return true;
}
