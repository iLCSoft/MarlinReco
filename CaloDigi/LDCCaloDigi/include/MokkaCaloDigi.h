#ifndef MokkaCaloDigi1_h
#define MokkaCaloDigi1_h 1

#include "marlin/Processor.h"
#include "EVENT/SimCalorimeterHit.h"
#include "IMPL/CalorimeterHitImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;

/**
\addtogroup CaloDigi CaloDigi
@{

\addtogroup MokkaCaloDigi MokkaCaloDigi
@{
 Calorimeter digitizer Processor for LCIO files produced by Mokka.
=== MokkaCaloDigi Processor === <br>
 *  Calorimeter digitizer Processor for LCIO <br>
 *  files produced by Mokka. <br>
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
 *  NewECALCollName and NewHCALCollName. <br>
 *  Processor performs merging of neighboring virtual cells <br>
 *  in calorimeters into a larger cells. <br>
 *  Processor is meant to perform digitization of <br>
 *  calorimeter hits for an arbitrary detector geometry <br>
 *  parameter NewHCALCellSize and should be multiple of 10 <br>
 *  as virtual cell size used in Mokka is 10x10 mm2. <br>
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
 *  @authors A. Raspereza and P. Krstonosic (DESY) <br>
 *  @version $Id$ <br>
 */

struct MyHit {
  CalorimeterHitImpl * hit{};
  std::vector<SimCalorimeterHit*> simHits{};
};


class MokkaCaloDigi : public Processor {
  
 public:
  
  MokkaCaloDigi(const MokkaCaloDigi&) = delete;
  MokkaCaloDigi& operator=(const MokkaCaloDigi&) = delete;

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


  MyHit * ProcessHitInBarrel(SimCalorimeterHit * hit);
  MyHit * ProcessHitInEndcap(SimCalorimeterHit * hit);

 protected:

  std::vector<std::string>  _ecalCollections{};
  std::vector<std::string>  _hcalCollections{};
  std::string _newCollNameHCAL{};
  std::string _newCollNameECAL{};
  std::string _relationCollName{};
  int _nRun{};
  int _nEvt{};

  float _thresholdEcal{};
  float _thresholdHcal{};

  int _digitalEcal{};
  int _digitalHcal{};

  std::vector<int> _ecalLayers{};
  std::vector<int> _hcalLayers{};
  std::vector<float> _calibrCoeffEcal{};
  std::vector<float> _calibrCoeffHcal{};
  float * _endBarrelChamberLength{};
  float * _barrelLateralWidth{};
  float * _barrelOffsetMaxX{};
  float * _endBarrelOffsetMaxZ{};
  float _regularBarrelOffsetMaxZ{};
  float _lateralPlateThickness{};
  float _modulesGap{};
  float _innerHcalRadius{};
  int _numberOfHcalLayers{};
  int _nStaves{};
  int _nModules{};
  int _cellScaleX{};
  int _cellScaleZ{};
  float _newCellSizeX{};
  float _newCellSizeZ{};
  float _hcalLayerThickness{};
  float _hcalAbsorberThickness{};
  float _hcalSensitiveThickness{};
  float _virtualCellSizeX{};
  float _virtualCellSizeZ{};
  float _regularBarrelModuleLength{};
  float _regularBarrelChamberLength{};
  float _deltaPhi{};
  std::vector< std::vector<MyHit*> > _calorimeterHitVec{};
  LCCollectionVec * _relationCollection{};
  int _startIEndcap{};
  int _startJEndcap{};
  float _startXEndcap{};
  float _startZEndcap{};



} ;

/** @} @}*/

#endif
