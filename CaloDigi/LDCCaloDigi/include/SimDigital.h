#ifndef SimDigital_HHH
#define SimDigital_HHH

#include "marlin/Processor.h"
#include "lcio.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <utility>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <marlin/Global.h>
#include <gear/GearParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>
#include <TROOT.h>
#include <TFile.h>
#include <TF2.h>
#include <TH1.h>
#include <TRandom3.h>
#include "TH1F.h"

#include "CalorimeterHitType.h" //in MarlinUtil

class TH1F;
class TF1;
class TTree;

using namespace lcio ;
using namespace marlin ;

const int MAX_LAYERS = 200;
const int MAX_STAVES =  16;


/** Digitization for the SDHcal - based on NewLDCCaloDigi. 
 *
 *  @author  G.Grenier, INPL
 *  @author  R.Han, INPL
 *  @version $Id$
 */

class SimDigital : public Processor {
 public:
  virtual Processor*  newProcessor() { return new SimDigital;}
  SimDigital();

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
  virtual void fillECALGaps() ;
 
  
 
 private:
  std::vector<std::string> _ecalCollections;
  std::string _outputEcalCollection0;
  std::string _outputEcalCollection1;
  std::string _outputEcalCollection2;
  std::vector<std::string> _outputEcalCollections;
  std::vector<float> _calibrCoeffEcal;
  std::vector<int> _ecalLayers;
  int _digitalEcal;
  float _thresholdEcal;
  int _ecalGapCorrection;
  float _ecalGapCorrectionFactor;
  float _ecalModuleGapCorrectionFactor;
  float _ecalEndcapCorrectionFactor;
  std::vector<CalorimeterHitImpl*> _calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS];
  std::vector<int> _calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS];

  float _zOfEcalEndcap;
  float _barrelPixelSizeT[MAX_LAYERS];
  float _barrelPixelSizeZ[MAX_LAYERS];
  float _endcapPixelSizeX[MAX_LAYERS];
  float _endcapPixelSizeY[MAX_LAYERS];
  float _barrelStaveDir[MAX_STAVES][2];
  //std::vector<int> _hcalLayers;
  //std::vector<std::string> _calorimeterHitCollections;
  std::vector<std::string> _hcalCollections;
  std::vector<std::string> _outputHcalCollections;
  std::map<std::string, int> _counters;
  std::vector<float> _thresholdHcal;
  std::vector<float>_calibrCoeffHcal;
  std::string _outputRelCollection;
  bool _printSimDigital;
  TF1 * _QPolya;
  //TF2 *c2d;  //NOT USED FOR THE MOMENT : WAITING LCIO v1.6 ? CHECK WITH RAN
  
  LCCollectionVec* _relcol;
  void processECAL(LCEvent* evt, LCFlagImpl& flag);
  void processHCAL(LCEvent* evt, LCFlagImpl& flag);

  //intermediate storage class
  struct hitMemory
  {
    hitMemory() : ahit(0),relatedHits(), rawHit(-1) {} 
    CalorimeterHitImpl *ahit;
    std::set<int> relatedHits;
    int rawHit;
  };
  //helper class to manage multiplicity
  struct multiplicityChargeSplitter
  {
    float cellSize,edgeSize;
    typedef std::pair<int,int> LeftRight_LowerUpper_Coordinates;
    enum DIRECTION {Central=0,Lower=-1,Upper=1,Left=-1,Right=1};
    std::map<LeftRight_LowerUpper_Coordinates,float> chargeMap;
    multiplicityChargeSplitter(float cell_Size, float edge_Size);
    void addCharge(float charge, float posLeftRight, float posLowerUpper);
  };
  TRandom3 _Ranm;
  float _edgedistance;
  bool _doThresholds;
  //float _RPC_PadSeparation; //distance between cells

  //predicate class to remove potential hit below threshold
  class ThresholdIsBelow
  {
    float value;
  public:
    ThresholdIsBelow(float f) : value(f) {;}
    bool operator()(std::pair<int,hitMemory> f) { return f.second.ahit->getEnergy()<value;}
  };

  CHT::Layout _currentHCALCollectionCaloLayout;
  const gear::LayerLayout& getGearLayout();

  LCCollectionVec * processHCALCollection(LCCollection * col ,LCFlagImpl& flag);
  void createPotentialOutputHits(std::map<int,hitMemory>& myHitMap, LCCollection *col, CellIDEncoder<CalorimeterHitImpl>& encodid);
  void removeHitsBelowThreshold(std::map<int,hitMemory>& myHitMap, float threshold);
  void applyThresholds(std::map<int,hitMemory>& myHitMap);

};
  
#endif
