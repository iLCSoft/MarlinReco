#ifndef TrackProcessor_h
#define TrackProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#ifdef MARLIN_USE_ROOT
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#endif

using namespace lcio ;
using namespace marlin ;

/** Plots the error on the inverse of the momentum for tracks reconstructed in the TPC.
* @author S. Aplin, DESY
*/

class TrackProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TrackProcessor ; }
  
  
  TrackProcessor() ;
  
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
  
 private:

  std::string _incol,_rootfilename;
  std::string _roothistname;
  int _rootnbins;
  float _roothistlow,_roothisthigh;  

#ifdef MARLIN_USE_ROOT
  TH1F *_massHisto;
  TFile *_histoFile;
#endif
 
 protected:

  /** Input collection name.
   */
  std::string _tpcTrackColName ;
  std::string _mcTrackRelColName ;
  int _nRun ;
  int _nEvt ;
} ;

#endif



