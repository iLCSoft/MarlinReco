/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef TPCDigiProcessor_h
#define TPCDigiProcessor_h 1

#include <marlin/Processor.h>
#include <lcio.h>
#include <string>
#include <gsl/gsl_rng.h>


using namespace lcio ;
using namespace marlin ;



/**  Produces TPC TrackerHit collection from SimTrackerHit collection, smeared in RPhi and Z. 
 * Double hits are identified but are currently not added to the collection. This may be change 
 * at a later date when criteria for their seperation is defined. The resolutions are defined in 
 * the GEAR stearing file.  
 * @author S. Aplin, DESY
 */
class TPCDigiProcessor : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new TPCDigiProcessor ; }
  
  
  TPCDigiProcessor() ;
  
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
  std::string _TPCTrackerHitsCol ;

  int _nRun ;
  int _nEvt ;

  // gsl random number generator
  gsl_rng * _random ;

} ;

#endif



