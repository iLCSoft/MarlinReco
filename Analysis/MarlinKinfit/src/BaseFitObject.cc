/*! \file 
 *  \brief Implements class BaseFitObject
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
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
