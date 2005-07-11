#ifndef MyProcessor_h
#define MyProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;



/** Example processor for marlin. If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * @author P.K , DESY
 */
class MokkaCaloDigi : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new MokkaCaloDigi ; }
  
  
  MokkaCaloDigi() ;
  
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
  std::vector<std::string>  _ecalCollections  ;
  std::vector<std::string>  _hcalCollections  ;
  std:: string _newCollNameHCAL;
  std:: string _newCollNameECAL;
  int _nRun ;
  int _nEvt ;
  int cell_size;

  float _thresholdEcal;
  float _thresholdHcal;

  int _digitalEcal;
  int _digitalHcal;


  std::vector<float> _calibrCoeffEcal;
  std::vector<float> _calibrCoeffHcal;

  std::vector<int> _ecalLayers;
  std::vector<int> _hcalLayers;

} ;

#endif



