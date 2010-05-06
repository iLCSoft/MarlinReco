#ifndef TauFinder_h
#define TauFinder_h 1

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


/** TauFinder processor for marlin.
 * 
 * @author A. Muennich, CERN
 * 
 */

class TauFinder : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TauFinder ; }
  
  
  TauFinder() ;
  
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
  std::string _colNameMC, _colNameRECO, _incol;
  std::string _colNameMCTruth, _colNameTauRecLink,  _outcol;

  int _nRun ;
  int _nEvt ;

  float _bField;
  float _ptcut,_ptseed;
  float _coneAngle,_isoAngle,_isoE;
  float _D0seed,_minv;

  int _fail_minv,_fail_minv_neg,_fail_Qtr,_fail_isoE;
   
  bool FindTau(std::vector<ReconstructedParticle*> &Qvec,std::vector<ReconstructedParticle*> &Nvec,
	       std::vector<std::vector<ReconstructedParticle*> > &tauvec);
  
} ;

#endif



