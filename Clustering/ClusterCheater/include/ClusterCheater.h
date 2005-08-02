#ifndef FORMTRUECLUSTERSAR_H
#define FORMTRUECLUSTERSAR_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include "HelixClass.h"
#include <EVENT/CalorimeterHit.h>
#include <EVENT/MCParticle.h>

using namespace lcio ;
using namespace marlin ;


/** Example processor for marlin. If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 */
class ClusterCheater : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ClusterCheater ; }
  
  
  ClusterCheater() ;
  
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
  int _ifBrahms;

  float _proximityCut;
  float _bField;
  int _minimal_hits;

  std::string _trueClustCollection;
  std::vector<std::string> _caloCollections;
  std::string _relCollection;

  float DistanceToChargeParticle(HelixClass * helix, CalorimeterHit * hit);
  float DistanceToNeutralParticle(MCParticle * par, CalorimeterHit * hit);
  HelixClass* AssignHelixToMCP( MCParticle * par);

} ;

#endif



