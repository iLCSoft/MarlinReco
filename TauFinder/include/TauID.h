#ifndef TauID_h
#define TauID_h 1

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
#include <IMPL/ReconstructedParticleImpl.h>

using namespace lcio ;
using namespace marlin ;


/** TauID processor for marlin.
 * 
 * @author A. Muennich, CERN
 * 
 */

class TauID : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TauID ; }
  
  
  TauID() ;
  
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
  std::string  _incol, _outcol;

  int _nRun ;
  int _nEvt ;

  float _bField;
  float _ptcut,_ptseed;
  float _coneAngle,_isoAngle,_isoE;
  float _D0seedmin, _D0seedmax,_minv;

  int _fail_minv,_fail_minv_neg,_fail_Qtr,_fail_isoE, _fail_Q;
   
  bool FindTau(std::vector<ReconstructedParticle*> &Qvec,std::vector<ReconstructedParticle*> &Nvec,
	       std::vector<ReconstructedParticleImpl* > &tauvec);
  
} ;

#endif



