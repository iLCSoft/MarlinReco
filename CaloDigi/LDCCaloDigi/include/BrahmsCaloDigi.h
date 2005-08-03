#ifndef DIGITIZERBRAHMS_H
#define DIGITIZERBRAHMS_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "DigiHitExtended.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;

/** === BrahmsCaloDigi Processor === <br>
 *  Calorimeter digitizer Processor for LCIO <br>
 *  files produced by Brahms. <br>
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
 *  Processor performs merging of neighboring virtual cells <br>
 *  in calorimeters into a larger cells. <br>
 *  Transverse cell size (in mm) is specified with processor <br>
 *  parameter TileSize and should be multiple of 10 <br>
 *  as virtual cell size used in Brahms is 10x10 mm2. <br>
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
 *  @version $ld: $ <br>
 */

class BrahmsCaloDigi : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new BrahmsCaloDigi ; }
  
  
  BrahmsCaloDigi() ;
  
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

  int _nRun ;
  int _nEvt ;
  
  std::vector<std::string> _ecalCollections;
  std::vector<std::string> _hcalCollections;

  std::string _outputEcalCollection;
  std::string _outputHcalCollection;
  std::string _outputRelCollection;

  float _thresholdEcal;
  float _thresholdHcal;

  float _tWall;
  float _zBarrelHcal;
  float _zEndcapHcal;
  float _r1BarrelHcal;
  float _r2BarrelHcal;
  float _halfEndcapHole;
  float _layerThickHcal;
  float _zInLayer;

  float _y0BarrelModule;
  float _y1BarrelModule;
  float _x1BarrelModule;
  float _x2BarrelModule;

  int _digitalHcal;
  int _digitalEcal;


  int _nLayerBarrel1;
  int _nLayerBarrel;

  int _tileSize;

  float _const_pi;
  float _const_pi4;
  float _const_pi8;
  float _const_twopi;

  std::vector<float> _calibrCoeffEcal;
  std::vector<float> _calibrCoeffHcal;

  std::vector<int> _ecalLayers;
  std::vector<int> _hcalLayers;

  std::vector<DigiHitExtended*> _HitVector;

  void getCell(float xhit, float yhit, float zhit, int & Module, int & Stave, int & SubModule, int & icell, int & jcell, int & layer);

  void getCellEndcap(float xhit, float yhit, float zhit, int & Stave, int & icell, int & jcell, int & layer);


  void getCellBarrel(float xhit, float yhit, float zhit, int & Stave, int & SubModule, int & icell, int & jcell, int & layer);

  void getCoordinates(int Module, int Stave, int SubModule, int icell, int jcell, int layer, float & xnew, float & ynew, float & znew);

  void getCoordinatesEndcap(int Module, int Stave, int icell, int jcell, int layer, float & xnew, float & ynew, float & znew);

  void getCoordinatesBarrel(int Module, int Stave, int SubModule, int icell, int jcell, int layer, float & xnew, float & ynew, float & znew);


  int EncodeCellID(int Module, int Stave, int SubModule, int icell, int jcell, int layer);

  void CleanUp();


};

#endif



