#ifndef Use_TrueJet_h
#define Use_TrueJet_h 1

//#include "TrueJet.h"
#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/LCRelation.h>
#include <UTIL/LCRelationNavigator.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/ParticleIDImpl.h>
#include <string>
#include "TrueJet_Parser.h"

using namespace lcio ;
using namespace marlin ;




class Use_TrueJet : public Processor , public TrueJet_Parser {

 public:
  
  virtual Processor*  newProcessor() { return new Use_TrueJet ; }
  
  
  Use_TrueJet();
  
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
  
  
  virtual void check( LCEvent * evt) ; 
  

  /** Called after data processing for clean up.
   */
  virtual void end() ;

  std::string get_recoMCTruthLink(){ return _recoMCTruthLink  ; } ;

 protected:

  /** Input collection name.
   */

  int _nRun ;
  int _nEvt ;

private:


  std::string  _MCParticleColllectionName ;
  std::string _recoParticleCollectionName ;
  std::string _recoMCTruthLink ;

} ;

#endif



