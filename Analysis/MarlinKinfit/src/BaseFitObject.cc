/*! \file 
 *  \brief Implements class BaseFitObject
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.1  2008/02/12 10:19:08  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.8  2008/02/04 17:30:53  blist
 * - NewtonFitter works now!
 * -
 * - Revision 1.7  2008/01/29 17:17:33  blist
 * - implemented setname
 * -
 * - Revision 1.6  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.5  2007/09/13 13:33:06  blist
 * - Print methods return os
 * -
 *
 */ 
 
#include "BaseFitObject.h"
#include <cassert>
#include <cstring>
#include <iostream>

BaseFitObject::BaseFitObject(): name(0) {
  setName ("???");
}

BaseFitObject::BaseFitObject (const BaseFitObject& rhs)
: name(0)
{
  //std::cout << "copying BaseFitObject with name" << rhs.name << std::endl;
  if (rhs.name) setName(rhs.name);
  else setName ("???");
}
BaseFitObject& BaseFitObject::operator= (const BaseFitObject& rhs) {
  if (this != &rhs) {
    if (rhs.name) setName(rhs.name);
    else setName ("???");
  }
  return *this;
}
    
BaseFitObject::~BaseFitObject() {
  //std::cout << "destroying BaseFitObject with name" << name << std::endl;
  delete[] name;
}

const double BaseFitObject::eps2 = 0.00001;

void  BaseFitObject::setName (const char * name_) {
  if (name_ == 0) return;
  size_t l = strlen(name_);
  if (name) delete[] name;
  name = new char[l+1];
  strcpy (name, name_);
}

int BaseFitObject::getNMeasured() const {
  int nmeasured = 0;
  for (int i = 0; i < getNPar(); ++i) if (isParamMeasured(i) && !isParamFixed(i)) ++nmeasured;
  return nmeasured;
}
int BaseFitObject::getNUnmeasured() const {
  int nunmeasrd = 0;
  for (int i = 0; i < getNPar(); ++i) if (!isParamMeasured(i) && !isParamFixed(i)) ++nunmeasrd;
  return nunmeasrd;
}
int BaseFitObject::getNFree() const {
  int nfree = 0;
  for (int i = 0; i < getNPar(); ++i) if (!isParamFixed(i)) ++nfree;
  return nfree;
}
int BaseFitObject::getNFixed() const {
  int nfixed = 0;
  for (int i = 0; i < getNPar(); ++i) if (isParamFixed(i)) ++nfixed;
  return nfixed;
}
    
std::ostream& BaseFitObject::printParams(std::ostream& os) const {
  os << "(";
  for (int i = 0; i < getNPar(); ++i) {
    if (i>0) os << ", ";
    os << getParam(i);
    if (isParamFixed (i))  os << " fix";
    else if (getError(i)>0) os << " \261 " << getError(i);
  }
  os << ")";  
  return os;
}
bool BaseFitObject::updateParams (double p[], int idim) {
  bool result = false;
  invalidateCache();
  for (int ilocal = 0; ilocal < getNPar(); ++ilocal) {
    int iglobal = getGlobalParNum (ilocal);
    assert (iglobal >= 0 && iglobal < idim);
    result = result || setParam (ilocal, p[iglobal]);
    // if illegal value: read back legal value
    p[iglobal] = getParam (ilocal);
  }
  return result;
}  
                                
