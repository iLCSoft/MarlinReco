/*! \file 
 *  \brief Declares class FourJetPairing
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.3  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 *
 */ 

#ifndef __FOURJETPAIRING_H
#define __FOURJETPAIRING_H

#include <iostream>
#include "BaseJetPairing.h"
#include "JetFitObject.h"

//  Class FourJetPairing:
/// Class to handle permutations of 2b jets and 4 light jets
/**
 *
 * Author: Jenny List, Anca Siebel
 * Last update: $Date: 2007-10-30 15:51:14 $
 *          by: $Author: gaede $
 *
 */

class FourJetPairing : public BaseJetPairing {
  public:
    /// constructor
    FourJetPairing (JetFitObject *jets_[]);
    
    /// Virtual destructor
    virtual ~FourJetPairing() {};    
        
    /// Number of permutaions
    virtual int getNPerm() const {return NPERM;};
    
    /// does the job
    virtual int nextPermutation (JetFitObject *permObjects[]);
    
  protected:
    enum {NPERM = 3};
    enum {NJETS = 4};
    JetFitObject *jets[NJETS]; 
    int permutations [NPERM][NJETS];

};
    
#endif // __FOURJETPAIRING_H

