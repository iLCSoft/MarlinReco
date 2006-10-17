/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef VTXDigiProcessor_h
#define VTXDigiProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;


/**  Produces VTX TrackerHit collection from SimTrackerHit collection. At present no smearing is applied
 * @author S. Aplin, DESY
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

  /** Input collection name.
   */
  std::string _colName ;
  std::string _colNameSIT ;
  std::string _outColName ;
  std::string _outColNameSIT ;

  int _nRun ;
  int _nEvt ;
} ;

#endif



