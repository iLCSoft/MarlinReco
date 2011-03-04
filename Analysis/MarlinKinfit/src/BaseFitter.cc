/*! \file 
 *  \brief Implements class BaseFitter
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: BaseFitter.cc,v $
 * - Revision 1.3  2011/03/03 15:03:02  blist
 * - Latest version, with NewFitterGSL
 * -
 * - Revision 1.2  2009/09/01 09:48:13  blist
 * - Added tracer mechanism, added access to fit covariance matrix
 * -
 * - Revision 1.1  2008/02/13 12:37:38  blist
 * - new file BaseFitter.cc
 * -
 * -
 *
 */ 
 
#include "BaseFitter.h"
#include "BaseSoftConstraint.h"
#include "BaseHardConstraint.h"
#include <cassert>

BaseFitter::BaseFitter()  
: covDim (0), cov(0), covValid (false)
#ifndef FIT_TRACEOFF    
  , tracer (0)
#endif FIT_TRACEOFF    
{}

BaseFitter::~BaseFitter()  
{
  delete[] cov;
  cov = 0;
}

void BaseFitter::addFitObject (BaseFitObject* fitobject_)  
{ 
  covValid = false;
  fitobjects.push_back(fitobject_);
}

void BaseFitter::addFitObject (BaseFitObject& fitobject_)  
{
  covValid = false;
  fitobjects.push_back(&fitobject_);
}

void BaseFitter::addConstraint (BaseConstraint* constraint_)  
{
  covValid = false;
  if (BaseHardConstraint *hc = dynamic_cast<BaseHardConstraint *>(constraint_))
    constraints.push_back(hc);
  else if (BaseSoftConstraint *sc = dynamic_cast<BaseSoftConstraint *>(constraint_))
    softconstraints.push_back(sc);
  else {
    // illegal constraint
    assert (0);
  }
}

void BaseFitter::addConstraint (BaseConstraint& constraint_)  
{
  covValid = false;
  if (BaseHardConstraint *hc = dynamic_cast<BaseHardConstraint *>(&constraint_)) 
    constraints.push_back(hc);
  else if (BaseSoftConstraint *sc = dynamic_cast<BaseSoftConstraint *>(&constraint_)) 
    softconstraints.push_back(sc);
}

void BaseFitter::addHardConstraint (BaseHardConstraint* constraint_)  
{
  covValid = false;
  constraints.push_back(constraint_);
}

void BaseFitter::addHardConstraint (BaseHardConstraint& constraint_) {
  covValid = false;
  constraints.push_back(&constraint_);
}

void BaseFitter::addSoftConstraint (BaseSoftConstraint* constraint_)  
{
  covValid = false;
  softconstraints.push_back(constraint_);
}

void BaseFitter::addSoftConstraint (BaseSoftConstraint& constraint_)  
{
  covValid = false;
  softconstraints.push_back(&constraint_);
}

std::vector<BaseFitObject*>* BaseFitter::getFitObjects() 
{
  return &fitobjects;
}

std::vector<BaseHardConstraint*>* BaseFitter::getConstraints()  
{
  return &constraints;
}

std::vector<BaseSoftConstraint*>* BaseFitter::getSoftConstraints()  
{
  return &softconstraints;
}

void BaseFitter::reset() 
{
  fitobjects.resize(0);
  constraints.resize(0);
  softconstraints.resize(0);
  covValid = false;
}  
    
BaseTracer *BaseFitter::getTracer() { 
  return tracer; 
}
const BaseTracer *BaseFitter::getTracer() const { 
  return tracer; 
}
void BaseFitter::setTracer(BaseTracer *newTracer) {
  tracer = newTracer; 
}
void BaseFitter::setTracer(BaseTracer& newTracer) {
  tracer = &newTracer; 
}

const double *BaseFitter::getGlobalCovarianceMatrix (int& idim) const {
  if (covValid && cov) {
    idim = covDim;
    return cov;
  }
  idim = 0;
  return 0;
}                                          

double *BaseFitter::getGlobalCovarianceMatrix (int& idim) {
  if (covValid && cov) {
    idim = covDim;
    return cov;
  }
  idim = 0;
  return 0;
}                                          
