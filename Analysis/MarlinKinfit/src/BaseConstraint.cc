/*! \file 
 *  \brief Implements class BaseConstraint
 *
 * \b Changelog:
 *
 *
 * \b CVS Log messages:
 *
 */ 
 
#include "BaseConstraint.h"
#include <cassert>
#include <cstring>
#include <iostream>

BaseConstraint::BaseConstraint()
:name (0)
{
  setName("???");
}

BaseConstraint::BaseConstraint (const BaseConstraint& rhs)
: name(0)
{
  if (rhs.name) setName(rhs.name);
  else setName ("???");
}
BaseConstraint& BaseConstraint::operator= (const BaseConstraint& rhs) {
  if (this != &rhs) {
    if (rhs.name) setName(rhs.name);
    else setName ("???");
  }
  return *this;
}
    

BaseConstraint::~BaseConstraint()
{
  //std::cout << "destroying BaseConstraint with name" << name << std::endl;
  delete[] name;
}

const char*BaseConstraint::getName() const {
  return name;
}

void  BaseConstraint::setName (const char * name_) {
  if (name_ == 0) return;
  size_t l = strlen(name_);
  if (name) delete[] name;
  name = new char[l+1];
  strcpy (name, name_);
}

double BaseConstraint::getError() const {
  assert (false); 
  return 0;
}
