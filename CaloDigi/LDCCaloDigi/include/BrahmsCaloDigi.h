#ifndef DIGITIZERBRAHMS_H
#define DIGITIZERBRAHMS_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "DigiHitExtended.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;

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



