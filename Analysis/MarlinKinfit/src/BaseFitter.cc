/*! \file 
 *  \brief Implements class BaseFitter
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
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
{}

BaseFitter::~BaseFitter()  
{}

void BaseFitter::addFitObject (BaseFitObject* fitobject_)  
{ 
  fitobjects.push_back(fitobject_);
}

void BaseFitter::addFitObject (BaseFitObject& fitobject_)  
{
  fitobjects.push_back(&fitobject_);
}

void BaseFitter::addConstraint (BaseConstraint* constraint_)  
{
  if (BaseHardConstraint *hc = dynamic_cast<BaseHardConstraint *>(constraint_))
    constraints.push_back(hc);
  else if (BaseSoftConstraint *sc = dynamic_cast<BaseSoftConstraint *>(constraint_))
    softconstraints.push_back(sc);
}

void BaseFitter::addConstraint (BaseConstraint& constraint_)  
{
  if (BaseHardConstraint *hc = dynamic_cast<BaseHardConstraint *>(&constraint_)) 
    constraints.push_back(hc);
  else if (BaseSoftConstraint *sc = dynamic_cast<BaseSoftConstraint *>(&constraint_)) 
    softconstraints.push_back(sc);
}

void BaseFitter::addHardConstraint (BaseHardConstraint* constraint_)  
{
  constraints.push_back(constraint_);
}

void BaseFitter::addHardConstraint (BaseHardConstraint& constraint_) {
  constraints.push_back(&constraint_);
}

void BaseFitter::addSoftConstraint (BaseSoftConstraint* constraint_)  
{
  softconstraints.push_back(constraint_);
}

void BaseFitter::addSoftConstraint (BaseSoftConstraint& constraint_)  
{
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
}  
