#ifndef EvaluateTauFinder_h
#define EvaluateTauFinder_h 1

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
#include "UTIL/LCRelationNavigator.h"

using namespace lcio ;
using namespace marlin ;

struct MyParticle;

/**  Evaluation processor for TauFinder
 * 
 * @author A. Muennich, CERN
 */

class EvaluateTauFinder : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new EvaluateTauFinder ; }
  
  
  EvaluateTauFinder() ;
  
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
  virtual void LoopDaughters(MCParticle *particle,double &Evis,double &ptvis,double &pvis);
  virtual void LoopDaughtersRelation(MCParticle *particle,LCRelationNavigator* relationNavigatorTau ,
				     LCRelationNavigator* relationNavigatorMC ,bool &ralToTau);
  
 protected:

  /** Input collection name.
   */
  std::string _colNameMC, _colNameRECO,_colNameTrack;
  std::string _colNameMCTruth,  _incol,_colNamePFORecLink;
  std::string _colNameMCRecLink,_colNameTracksRecLink, _colNameTauRecLink;
  int _nRun ;
  int _nEvt ;
  
  double _ntot_rec;
  double _ntot_mc;
  double _ntau_correct;
  double _dEsum;
  double _dEsumsq;
  int _ndE;

  float _bField;
  TFile *rootfile;
  TNtuple *tautuple;
  TNtuple *mcmisstuple;
  TNtuple *taumatchtuple;
  TNtuple *tauexacttuple;
  TNtuple *evtuple;
} ;

#endif



