////////////////////////////////////////////////////////////////
// Class PyConstraint
//
// Author: Jenny Boehme
// Last update: $Date: 2007-10-30 15:51:14 $
//          by: $Author: gaede $
// 
// Description: sum (px) = 0 constraint
//               
////////////////////////////////////////////////////////////////

#include "PyConstraint.h"
#include "ParticleFitObject.h"

#include<iostream>
#include<cassert>

using std::cout;
using std::endl;

// constructor
PyConstraint::PyConstraint () {}

// destructor
PyConstraint::~PyConstraint () {}

// calulate current value of constraint function
double PyConstraint::getValue() const {
  double totpy = 0;
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    totpy += fitobjects[i]->getPy(); 
  }
  return totpy;
}

// calculate vector/array of derivatives of this contraint 
// w.r.t. to ALL parameters of all fitobjects
// here: d sum(py) /d par(i,j) 
//                      = d sum(py) /d py(i) * d py(i) /d par(i, j)
//                                      =  1 * d py(i) /d par(i, j)
void PyConstraint::getDerivatives(int idim, double der[]) const {
  for (unsigned int i = 0; i < fitobjects.size(); i++) {
    for (int ilocal = 0; ilocal < fitobjects[i]->getNPar(); ilocal++) {
      int iglobal = fitobjects[i]->getGlobalParNum (ilocal);
      if (!fitobjects[i]->isParamFixed(ilocal)) {
        assert (iglobal >= 0 && iglobal < idim);
        der[iglobal] = fitobjects[i]->getDPy(ilocal);
//       cout << "Py: der[" << iglobal << "] = " << der[iglobal] 
//            << " for jet " << i << " and ilocal = " << ilocal << endl;
      }
    }
  }
}
  

  

