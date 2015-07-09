#ifndef IsolatedLeptonTaggingProcessor_h
#define IsolatedLeptonTaggingProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include "TMVA/Reader.h"

/**  processor for isolated lepton tagging.
 * @version $Id: IsolatedLeptonTaggingProcessor.h,v 1.0 2013-10-31 12:57:39 Junping Exp $ 
 */

class IsolatedLeptonTaggingProcessor : public marlin::Processor {
  
 public:
  
  virtual marlin::Processor*  newProcessor() { return new IsolatedLeptonTaggingProcessor ; }
  
  
  IsolatedLeptonTaggingProcessor() ;
  
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
  std::string _colPFOs ;
  std::string _colNewPFOs ;
  std::string _colLeptons ;
  std::string _isolated_electron_weights ;
  std::string _isolated_muon_weights ;

  bool _is_one_isolep ;
  float _minEOverPForElectron ;
  float _maxEOverPForElectron ;
  float _minEecalOverTotEForElectron ;
  float _minPForElectron ;
  float _maxD0SigForElectron ;
  float _maxZ0SigForElectron ;  
  float _maxEOverPForMuon ;
  float _minEyokeForMuon ;
  float _minPForMuon ;  
  float _maxD0SigForMuon ;
  float _maxZ0SigForMuon ;  

  float _cosConeSmall ;
  float _cosConeLarge ;  
  
  std::vector<TMVA::Reader*> _readers;
  Float_t _coneec, _coneen, _momentum, _coslarcon, _energyratio;
  Float_t _ratioecal, _ratiototcal, _nsigd0, _nsigz0, _yokeenergy, _totalcalenergy;

  float _mvaCutForElectron, _mvaCutForMuon;

} ;

#endif



