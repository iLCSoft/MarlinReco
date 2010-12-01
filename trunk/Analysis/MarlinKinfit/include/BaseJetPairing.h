/*! \file 
 *  \brief Declares class BaseJetPairing
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: BaseJetPairing.h,v $
 * - Revision 1.1  2008/02/12 10:19:05  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.2  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 *
 */ 
 
#ifndef __BASEJETPAIRING_H
#define __BASEJETPAIRING_H

class JetFitObject;

//  Class BaseConstraint:
/// Abstract base class for jet pairings
/**
 *
 * Author: Jenny List, Benno List
 * Last update: $Date: 2008/02/12 10:19:05 $
 *          by: $Author: blist $
 *
 */
class BaseJetPairing {
  public:
    virtual ~BaseJetPairing() {};
    
    // setters
    virtual void reset() {iperm = 0;};
    
    // getters
    virtual int getNPerm() const = 0;
    
    // does the job
    virtual int nextPermutation (JetFitObject *permObjects[]) = 0;
    
  protected:
    int iperm; 

};
    
#endif // __BASEJETPAIRING_H

