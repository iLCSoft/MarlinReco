/*! \file 
 *  \brief Implements class TopEventILCSoft
 *
 * Author: Benno List
 * Last update: $Date: 2011/05/03 13:18:48 $
 *          by: $Author: blist $
 * $Date: 2011/05/03 13:18:48 $
 * $Author: blist $
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: TopEventILCSoft.h,v $
 * - Revision 1.3  2011/05/03 13:18:48  blist
 * - New test stuff
 * -
 * -
 *
 */ 

////////////////////////////////////////////////////////////////
// Class TopEventILCSoft
//
// Author: Benno List, Jenny Boehme
// Last update: $Date: 2011/05/03 13:18:48 $
//          by: $Author: blist $
// 
// Description: class to generate and fit top pair events at ILC
//               
////////////////////////////////////////////////////////////////
#ifndef __TOPEVENTILCSOFT_H
#define __TOPEVENTILCSOFT_H

#include "BaseEvent.h"
#include "JetFitObject.h"
#include "PConstraint.h"
#include "SoftGaussMassConstraint.h"
#include "SoftBWMassConstraint.h"


// Class TopEventILCSoft
/// Class to generate and fit top pair events at ILC, using soft constraints
/**
 *
 * Author: Benno List
 * Last update: $Date: 2011/05/03 13:18:48 $
 *          by: $Author: blist $
 *
 * \b Changelog:
 * - 15.11.2010 First version
 *
 */
class TopEventILCSoft : public BaseEvent {
  public: 
    TopEventILCSoft();
    virtual ~TopEventILCSoft();
    virtual void genEvent();
    virtual int fitEvent (BaseFitter& fitter);

    double bwrandom (double r, double e0, double gamma, double emin, double emax) const;
    
    BaseConstraint& getPxConstraint() {return pxc;};
    BaseConstraint& getPyConstraint() {return pyc;};
    BaseConstraint& getW1Constraint() {return w1;};
    BaseConstraint& getW2Constraint() {return w2;};
    BaseConstraint& getTopConstraint() {return w;};
    
    double getW1Mass()  {return w1.getMass();};
    double getW2Mass()  {return w2.getMass();};
    double getTopMass(int flag)  {return w.getMass(flag);};
    
    ParticleFitObject* getTrueFitObject (int i) {return bfo[i];};
    ParticleFitObject* getSmearedFitObject (int i) {return bfosmear[i];};
    FourVector* getTrueFourVector (int i) {return fv[i];};
    
    bool leptonic;
    
  protected:
    enum {NFV = 11, NBFO = 6};
    FourVector *fv[NFV];
    FourVector *fvsmear[NFV];
    FourVector *fvfinal[NFV];
    ParticleFitObject *bfo[NBFO];
    ParticleFitObject *bfosmear[NBFO];
    
    PConstraint pxc;
    PConstraint pyc;
    PConstraint pzc;
    PConstraint ec;
//     SoftBWMassConstraint w1;
//     SoftBWMassConstraint w2;
//     SoftBWMassConstraint w;
    SoftGaussMassConstraint w1;
    SoftGaussMassConstraint w2;
    SoftGaussMassConstraint w;
    
    

};


#endif // __TOPEVENTILCSOFT_H
