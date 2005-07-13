#ifndef TrackingPerformanceProcessor_h
#define TrackingPerformanceProcessor_h 1

#include <marlin/Processor.h>
#include <lcio.h>
#include <string>

#ifdef MARLIN_USE_ROOT
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#endif

using namespace lcio ;
using namespace marlin ;

class TrackingPerformanceProcessor : public Processor {

 public:

  virtual Processor* newProcessor() { return new TrackingPerformanceProcessor ; }

  TrackingPerformanceProcessor() ;

  /** Called at eh begin of the job before anything is read.
   *  Used to initialize the processor, e.g. book histograms.
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

#ifdef MARLIN_USE_ROOT
  TH1I *_efficiencyHisto;
  TH1I *_efficiencyHisto1;
  TH1I *_efficiencyHisto2;
  TH1I *_efficiencyHisto3;
  TH1I *_efficiencyHisto4;
  TH1F *_invpHisto1;
  TH1F *_invpHisto2;
  TH1F *_invpHisto3;
  TH1F *_invpHisto4;
  TH1I *_ghostHisto;

  TFile *_histoFile;

#endif
  
 protected:

   /** Input collection name.
   */
  std::string _colName ;

  int _nRun ;
  int _nEvt ;
} ;

#endif
