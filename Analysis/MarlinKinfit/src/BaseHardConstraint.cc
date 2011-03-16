/*! \file 
 *  \brief Implements class BaseHardConstraint
 *
 * \b Changelog:
 * - 15.11.2010 First version
 *
 *
 * \b CVS Log messages:
 * - $Log: BaseHardConstraint.cc,v $
 * - Revision 1.1  2011/03/03 15:03:02  blist
 * - Latest version, with NewFitterGSL
 * -
 *
 */ 
 
#include "BaseHardConstraint.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include <cmath>

BaseHardConstraint::~BaseHardConstraint()
{}

double BaseHardConstraint::dirDer (double *p, double *w, int idim, double mu) {
  double *pw, *pp;
  for (pw = w; pw < w+idim; *(pw++) = 0);
  addToGlobalChi2DerVector (w, idim, mu);
  double result = 0;
  for (pw = w, pp = p; pw < w+idim; result += *(pp++) * *(pw++));
  return mu*result;
}

double BaseHardConstraint::dirDerAbs (double *p, double *w, int idim, double mu) {
  double val = getValue();
  if (val == 0) return mu*std::fabs(dirDer (p, w, idim, 1));
  else if (val > 0) return mu*dirDer (p, w, idim, 1);
  else return -mu*dirDer (p, w, idim, 1);
}
