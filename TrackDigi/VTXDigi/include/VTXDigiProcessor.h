/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef VTXDigiProcessor_h
#define VTXDigiProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

using namespace lcio ;
using namespace marlin ;


/**  Produces VTX TrackerHit collection from SimTrackerHit collection. 
 * @author A. Raspereza, MPI (Munich)
 */
class VTXDigiProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new VTXDigiProcessor ; }
  
  
  VTXDigiProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection names.
   */
  std::string _colNameVTX ;
  std::string _colNameSIT ;
  std::string _outColNameVTX;
  std::string _outColNameSIT;

  int _nRun ;
  int _nEvt ;
  int _debug;
  int _removeDRays;
  float _pointResoRPhi,_pointResoRPhi_VTX,_pointResoRPhi_SIT;
	float _pointResoZ,_pointResoZ_VTX,_pointResoZ_SIT;
  float _momCut;

  gsl_rng * r ;


} ;

#endif



