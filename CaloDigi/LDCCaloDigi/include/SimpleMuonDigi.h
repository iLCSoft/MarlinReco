#ifndef SimpleMuonDigi_H
#define SimpleMuonDigi_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

#include "CalorimeterHitType.h"


using namespace lcio ;
using namespace marlin ;


/** === SimpleMuonDigi Processor === <br>
 *  Simple calorimeter digitizer for the muon detectors.
 *  Converts SimCalorimeterHit collections to one 
 *  CalorimeterHit collection applying a threshold and an calibration constant...
 * 
 *  @version $Id: SimpleMuonDigi.h,v 1.1 2008-01-11 15:10:39 gaede Exp $
 */
class SimpleMuonDigi : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SimpleMuonDigi ; }
  
  
  SimpleMuonDigi() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  virtual void end() ;

  bool useLayer(CHT::Layout caloLayout, unsigned int layer) ;
  
 protected:

  int _nRun ;
  int _nEvt ;

  IntVec _layersToKeepBarrelVec, _layersToKeepEndcapVec;
  std::vector<bool>  _useLayersBarrelVec, _useLayersEndcapVec;

  std::vector<std::string> _muonCollections;

  std::string _outputMuonCollection;
  std::string _outputRelCollection;


  float _thresholdMuon;
  float _calibrCoeffMuon;
  float _maxHitEnergyMuon;


} ;

#endif



