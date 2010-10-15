#ifndef EvaluateTauID_h
#define EvaluateTauID_h 1

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

/**  Evaluation processor for TauID
 * 
 * @author A. Muennich, CERN
 */

class EvaluateTauID : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new EvaluateTauID ; }
  
  
  EvaluateTauID() ;
  
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
  virtual void LoopDaughters(MCParticle *particle,double &Evis,double &ptvis,int &nQ, int &nN, double &Qw, double &pQ);
 
  
 protected:

  /** Input collection name.
   */
  std::string _colNameMC, _colNameRECO;
  std::string _colNameMCTruth,  _incol,_colNamePFORecLink;
  std::string _colNameMCRecLink, relcol;
  int _nRun ;
  int _nEvt ;
  
  float _bField;
  double _kappa;
  TFile *rootfile;
  
  TNtuple *mctuple;
  TNtuple *taumatchtuple;
  TNtuple *evtuple;
  TNtuple *faketuple;
  TNtuple *Qtuple;
  
} ;

#endif



