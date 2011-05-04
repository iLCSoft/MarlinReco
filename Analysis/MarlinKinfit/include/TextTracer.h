/*! \file 
 *  \brief Declares class TextTracer
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: TextTracer.h,v $
 * - Revision 1.2  2011/05/03 13:18:48  blist
 * - New test stuff
 * -
 * - Revision 1.1  2009/09/01 09:48:13  blist
 * - Added tracer mechanism, added access to fit covariance matrix
 * -
 * - Revision 1.1  2009/05/22 12:12:07  blist
 * -
 *
 */ 

#ifndef __TEXTTRACER_H
#define __TEXTTRACER_H

#include <iostream>
#include "BaseTracer.h"

class BaseFitter;

//  Class TextTracer:
/// Class to produce text output during kinematic fits
/**
 *
 * Author: Benno List
 * Last update: $Date: 2011/05/03 13:18:48 $
 *          by: $Author: blist $
 *
 */
 
class BaseFitter; 
 
class TextTracer: public BaseTracer {
  public:
    TextTracer(std::ostream& os_);
    virtual ~TextTracer();
    
    /// Called at the start of a new fit (during initialization)
    virtual void initialize (BaseFitter& fitter);
    /// Called at the end of each step
    virtual void step (BaseFitter& fitter);
    /// Called at intermediate points during a step
    virtual void substep (BaseFitter& fitter,
                          int flag
                          );
    /// Called at the end of a fit
    virtual void finish (BaseFitter& fitter);
    
    void printFitObjects (BaseFitter& fitter);
    void printConstraints (BaseFitter& fitter);
    void printTraceValues (BaseFitter& fitter);
    void printSums (BaseFitter& fitter);
    
  protected:
    std::ostream& os;
    
    int istep;
    int isubstep;
    double chi2fo;
    double chi2sc;
    double sumhc;
    double sumhcscal;
};

#endif // __TEXTTRACER_H
