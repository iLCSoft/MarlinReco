/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef FTDDigiProcessor_h
#define FTDDigiProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

using namespace lcio ;
using namespace marlin ;


/**  Produces FTD TrackerHit collection from SimTrackerHit collection. At present simple gaussian smearing is applied
 * @author S. Aplin, DESY
 */
class FTDDigiProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new FTDDigiProcessor ; }
  
  
  FTDDigiProcessor() ;
  
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

  /** Input collection name.
   */
  std::string _colName ;
  std::string _outColName ;

  int _nRun ;
  int _nEvt ;
  float _pointReso;
  float _momCut;
  
  // gsl random number generator
  gsl_rng * r ;

} ;

#endif



