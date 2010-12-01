/*! \file 
 *  \brief Declares class BaseEvent
 *
 * \b Changelog:
 * - 7.6.04 JB: First doxygen docu
 *
 */ 

#ifndef __BASEEVENT_H
#define __BASEEVENT_H

#include "FourVector.h"
#include "BaseFitObject.h"
#include "BaseFitter.h"

// class BaseEvent
/// Abstract base class for different kinds of events
/**
 * This class defines the minimal functionality of a class describing
 * an event hypothesis that is ment to be tested by a kinematic fit.
 * It is an optional class - the event hypothesis can also be defined 
 * directly in the main program.
 * 
 *
 * Author: Jenny Böhme, Benno List
 * $Date: 2008/02/12 10:19:05 $
 * $Author: blist $
 *
 */

class BaseEvent  {
  public: 
    virtual ~BaseEvent() {};
    /// provides four-momenta (i.e. read values from ntuple, run toy MC, ...)
    virtual void genEvent() = 0;
    /// do it!
    virtual int fitEvent (BaseFitter& fitter) = 0;
  
};


#endif // __BASEEVENT_H
