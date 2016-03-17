#ifndef DIGITIZER_REALISTICCALORECO_H
#define DIGITIZER_REALISTICCALORECO_H 1

#include "marlin/Processor.h"
#include "lcio.h"

#include <EVENT/CalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>

#include <UTIL/CellIDDecoder.h>

#include <IMPL/LCFlagImpl.h>

#include "CalorimeterHitType.h"

#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;

const int MAX_LAYERS  = 35;
const int MAX_STAVES  = 12;
const int MAX_MODULES = 10;
const int MAX_WAFERS  = 100;

/** === RealisticCaloReco Processor === <br>
    realistic reconstruction of calorimeter hits
    e.g. apply sampling fraction correction
    virtual class, technology indenpendent
    D.Jeans 02/2016.
*/


class RealisticCaloReco : virtual public Processor {
  
 public:
  RealisticCaloReco() ;
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;

 protected:

  float getLayerCalib( int ilayer );
  virtual float reconstructEnergy(const CalorimeterHit* hit)=0;  // to be overloaded, technology-specific

  // geometry info, eg for gap filling
  // can be overlaoded if wanted
  virtual void getGeometryInformation() {} // called for every collection
  virtual void resetGaps() {} // for every event
  virtual void prepareForGaps(const CalorimeterHit* hit) {} // for every hit
  virtual void fillGaps() {}  // to be overlaoded

  // parameters
  std::vector <std::string> _inputHitCollections;
  std::vector <std::string> _outputHitCollections;
  std::vector <std::string> _outputRelCollections;
  std::vector <float> _calibrCoeff;
  std::vector <int> _calLayers;

  int   _gapCorrection;
  std::string _cellIDLayerString;
  std::string _cellIDModuleString;
  std::string _cellIDStaveString;
  std::string _cellIDWaferString;
  std::string _cellIDTowerString;
  std::string _cellIDIndexIString;
  std::string _cellIDIndexJString;

  // internal variables
  LCFlagImpl _flag;
  CellIDDecoder<CalorimeterHit> * _idDecoder;

  CHT::Layout   _cht_layout;
  CHT::CaloID   _cht_caloid;
  CHT::CaloType _cht_calotype;

  std::vector <CalorimeterHitImpl*> _calHitsByStaveModuleLayer[MAX_STAVES][MAX_MODULES][MAX_LAYERS];

  int _countWarnings;

} ;

#endif



