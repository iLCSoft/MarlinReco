#ifndef BRUTEFORCEECAL_GAPFILLER_H
#define BRUTEFORCEECAL_GAPFILLER_H 1

#include "marlin/Processor.h"
#include <IMPL/LCFlagImpl.h>
#include <EVENT/CalorimeterHit.h>

#include <IMPL/LCCollectionVec.h>

#include "lcio.h"
#include <string>
#include <vector>

#include "DDRec/DetectorData.h"


using namespace lcio ;
using namespace marlin ;

class BruteForceEcalGapFiller : public Processor {

 public:

  virtual Processor*  newProcessor() { return new BruteForceEcalGapFiller ; }


  BruteForceEcalGapFiller( ) ;
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) {}
  virtual void processEvent( LCEvent * evt ) ;
  virtual void check( LCEvent * evt ) {}
  virtual void end() ;

 protected:

  enum {MAXMODULE=10, MAXSTAVE=15, MAXLAYER=50};
  std::vector < CalorimeterHit* > hitsByLayerModuleStave[MAXLAYER][MAXSTAVE][MAXMODULE];

  std::string _inputHitCollection;
  std::string _outputHitCollection;

  std::string _cellIDLayerString;
  std::string _cellIDModuleString;
  std::string _cellIDStaveString;

  LCFlagImpl _flag;

  enum {ECALENDCAP, ECALBARREL};

  int _currentLayout;

  float _interModuleDist;

  dd4hep::rec::LayeredCalorimeterData* _caloGeomData;

  void getGeometryData(const int ihitType);

  void addIntraModuleGapHits( LCCollectionVec* newcol );
  void addInterModuleGapHits( LCCollectionVec* newcol );
 
} ;

#endif
