#ifndef FORMTRUETRACKSAR_H
#define FORMTRUETRACKSAR_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;


/** Example processor for marlin. If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 */
class TrackCheater : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TrackCheater ; }
  
  
  TrackCheater() ;
  
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

  int _nRun ;
  int _nEvt ;

  std::string _trueTracksCollection;
  std::vector<std::string> _trackerHitCollections;
  std::string _relCollection;
 
  float _bField;
  float _eCut;
 
} ;

#endif



