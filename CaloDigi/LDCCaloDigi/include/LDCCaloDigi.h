#ifndef DIGITIZER_LDCCALODIGI_H
#define DIGITIZER_LDCCALODIGI_H 1

#include "marlin/Processor.h"
#include <IMPL/CalorimeterHitImpl.h>
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;

const int MAX_LAYERS = 200;
const int MAX_STAVES =  16;

/** === LDCCaloDigi Processor === <br>
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
 *  Digital mode is activated by  <br>
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
 *  <h4>Input collections and prerequisites</h4>
 *  SimCalorimeterHit collections <br>
 *  <h4>Output</h4>
 *  CalorimeterHit collections for ECal and HCal. <br>
 *  Collection of relations <br>
 *  between CalorimeterHits and SimCalorimeterHits. <br> 
 *  For ECal Calorimeter hits the variable type is set to 0, <br>
 *  whereas for HCal Calorimeter hits the type is set to 1 <br>
 *  @author A. Raspereza (DESY) <br>
 *  @author M. Thomson (DESY) <br>
 *  @version $Id$ <br>
 */
class LDCCaloDigi : public Processor {
  
 public:
  
  LDCCaloDigi(const LDCCaloDigi&) = delete;
  LDCCaloDigi& operator=(const LDCCaloDigi&) = delete;

  virtual Processor*  newProcessor() { return new LDCCaloDigi ; }
  
  
  LDCCaloDigi() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
   
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;

  virtual void fillECALGaps() ;
  
  
 protected:

  int _nRun{};
  int _nEvt{};
  
  std::vector<std::string> _ecalCollections{};
  std::vector<std::string> _hcalCollections{};

  std::string _outputEcalCollection{};
  std::string _outputHcalCollection{};
  std::string _outputRelCollection{};

  float _thresholdEcal{};
  float _thresholdHcal{};

  int _digitalEcal{};
  int _digitalHcal{};

  std::vector<float> _calibrCoeffEcal{};
  std::vector<float> _calibrCoeffHcal{};

  std::vector<int> _ecalLayers{};
  std::vector<int> _hcalLayers{};

  int _ecalGapCorrection{};
  float _ecalGapCorrectionFactor{};
  float _ecalModuleGapCorrectionFactor{};
  float _ecalEndcapCorrectionFactor{};

  std::vector<CalorimeterHitImpl*> _calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS];
  std::vector<int> _calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS];

  float _zOfEcalEndcap{};
  float _barrelPixelSizeT[MAX_LAYERS]{};
  float _barrelPixelSizeZ[MAX_LAYERS]{};
  float _endcapPixelSizeX[MAX_LAYERS]{};
  float _endcapPixelSizeY[MAX_LAYERS]{};
  float _barrelStaveDir[MAX_STAVES][2]{};

} ;

#endif
