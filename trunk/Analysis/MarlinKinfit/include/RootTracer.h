/*! \file 
 *  \brief Declares class RootTracer
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: RootTracer.h,v $
 * - Revision 1.1  2010/05/25 13:23:56  boehmej
 * - fixed RootTracer
 * -
 * - Revision 1.2  2009/09/01 08:34:24  blist
 * - Kinfit update
 * -
 * - Revision 1.1  2009/05/22 12:12:07  blist
 * - New tracers
 * -
 * -
 *
 */ 

#ifndef __ROOTTRACER_H
#define __ROOTTRACER_H

#include <iostream>
#include "BaseTracer.h"

class BaseFitter;
class TFile;

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//  Class RootTracer:
/// Class to produce text output during kinematic fits
/**
 *
 * Author: Benno List
 * Last update: $Date: 2010/05/25 13:23:56 $
 *          by: $Author: boehmej $
 *
 */
 
class BaseFitter; 
 
class RootTracer: public BaseTracer {
  public:
    RootTracer(const char* filename="trace.root", const char *option="RECREATE");
    virtual ~RootTracer();
    
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
    
  protected:
    void SetBranchAddresses();
    void CreateBranches();
    void CreateEventBranches(BaseFitter& fitter);
    void FillParameterValues(BaseFitter& fitter);
  
    TFile *file;
    TTree *tree; 
    TTree *eventtree; 
    
    int istep;
    int isubstep;
    
    Int_t eventnumber;
    Int_t stepnumber;
    Int_t substepnumber;
    Double_t chi2;
    
    enum {NPARMAX = 100};
    Double_t parvalue[NPARMAX];
    
};

#endif // __ROOTTRACER_H
