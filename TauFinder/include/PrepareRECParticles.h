#ifndef PrepareRECParticles_h
#define PrepareRECParticles_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include "TTree.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>

using namespace lcio ;
using namespace marlin ;


/**  PreProcessor for TauFinder to provide necessary information
 *  and create universal input for TauFinder, so that the TauFinder
 *  can run on various input information
 *
 * @author A. Muennich, CERN
 *
 */

class PrepareRECParticles : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new PrepareRECParticles ; }
  
  
  PrepareRECParticles() ;
  
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
  std::string _colNameMC,_colNameTrack,_colNamePFO,_colNamePFOMCTruth,;
  std::string _colNameMCTruth, _colNameTrackTruth,_colNamePFOTruth;
  std::string _outcolMC, _outcolTracks ,_outcolPFO;
  double _bField;
  int _nRun ;
  int _nEvt ;

  
} ;

#endif



