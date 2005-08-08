#ifndef DIGITIZERAR_H
#define DIGITIZERAR_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;


/** === SimpleCaloDigi Processor === <br>
 *  Simple calorimeter digitizer Processor. <br>
 *  Takes SimCalorimeterHit Collections and <br>
 *  produces CalorimeterHit Collections. <br>
 *  Simulated energy depositions in active <br>
 *  layers of calorimeters are <br>
 *  converted into physical energy. This is done <br>
 *  taking into account sampling fractions of <br>
 *  ECAL and HCAL. <br>
 *  User has to specify ECAL and HCAL SimCalorimeterHit <br>
 *  collections with processor parameters <br>
 *  HCALCollections and ECALCollections. <br>
 *  The names of the output CalorimeterHit Collections <br>
 *  are specified with processor parameters <br>
 *  ECALOutputCollection and HCALOutputCollection. <br>
 *  Conversion factors for ECAL and HCAL <br>
 *  are specified via processor parameters  <br>
 *  CalibrECAL and CalibrHCAL. <br>
 *  It should be noted that ECAL and HCAL may consist <br>
 *  of several sections with different sampling fractions. <br>
 *  To handle this situation, calibration coefficients for <br>
 *  ECAL and HCAL are passed as arrays of floats with each element <br>
 *  in this array corresponding to certain section with <br>
 *  a given sampling fraction. <br>
 *  List of layer numbers terminating each section are given through <br>
 *  processor parameters ECALLayers and HCALLayers <br>
 *  There is an option to perform digitization of <br> 
 *  both ECAL and HCAL in a digital mode. <br>
 *  Digital digitization is activated by  <br>
 *  setting processor parameters <br>
 *  IfDigitalEcal / IfDigitalHcal to 1. <br>
 *  In this case CalibrECAL / CalibrHCAL will  <br>
 *  convert the number of hits into physical energy. <br>
 *  Thresholds on hit energies in ECAL and HCAL <br>
 *  are set with processor parameters <br>
 *  ECALThreshold and HCALThreshold.  <br>
 *  Relations between CalorimeterHits and SimCalorimeterHits <br>
 *  are held in the corresponding relation collection. <br>
 *  The name of this relation collection is specified <br>
 *  via processor parameter RelationOutputCollection. <br> 
 *  @author A. Raspereza (DESY) <br>
 *  @version $Id: SimpleCaloDigi.h,v 1.3 2005-08-07 16:22:39 gaede Exp $ <br>
 */
class SimpleCaloDigi : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SimpleCaloDigi ; }
  
  
  SimpleCaloDigi() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  virtual void end() ;
  
  
 protected:

  int _nRun ;
  int _nEvt ;
  
  std::vector<std::string> _ecalCollections;
  std::vector<std::string> _hcalCollections;

  std::string _outputEcalCollection;
  std::string _outputHcalCollection;
  std::string _outputRelCollection;


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



