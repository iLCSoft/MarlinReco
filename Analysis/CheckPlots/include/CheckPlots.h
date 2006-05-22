#ifndef CheckPlots_h
#define CheckPlots_h 1

#include <iostream>
#include <string>
#include <vector>

#include "marlin/Processor.h"
#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/CalorimeterHit.h>


#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#endif


using namespace lcio ;
using namespace marlin ;


/**
 *    This processor provides check plots. The plots distributed in the following different categories:
 *    <br>
 *    1. MC particle related plots <br>
 *    2. Plots related to the simulated hits in tracking and calorimeter devices <br>
 *    3. Plots related to the hits in tracking and calorimeter devices <br>
 *    4. more to come <br>
 *    <br>
 *    Steering parameters: <br>
 *    FillMC: toggles the check plots for the MC particles with generator status == 1 ( 0 or 1 ) <br>
 *    FillMCSim : toggles the check plots for the MC particles with generator status != 1 ( 0 or 1 ) <br>
 *    FillSimCalo : produces check plots for the simulated calorimeter hits <br>
 *    SimECut : cut on the simulated energy in the cell. Only energies above this cut are taken into account. <br>
 *    FillCalo : produces check plots for the calorimeter hits <br>
 *    ECut : cut on the energy in the cell. Only energies above this cut are taken into account. <br>
 *
 *
 *    @author O. Wendt (DESY)
 *    @version $Id: CheckPlots.h,v 1.1 2006-05-22 13:16:44 owendt Exp $
 *
 */
class CheckPlots : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CheckPlots ; }
    
  CheckPlots() ;
  
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  
  
 private:

  int _nRun;
  int _nEvt;

  int _fillMC;
  int _fillMCSim;

  int _fillSimCalo;
  float _simECut;
  int _fillCalo;
  float _ECut;

  void fillMCCheckPlots(LCEvent * evt);

  //  void fillSimTrackerCheckPlots(LCEvent * evt);
  void fillSimCaloCheckPlots(LCEvent * evt);

  //  void fillTrackerCheckPlots(LCEvent * evt);
  void fillCaloCheckPlots(LCEvent * evt);

} ;

#endif



