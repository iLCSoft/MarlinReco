#ifndef DIGITIZERLHCAL_H
#define DIGITIZERLHCAL_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;


/** === SimpleLHCalDigi Processor === <br>
 *  Simple calorimeter digitizer for the LCal Processor. <br>
 *  Converts SimCalorimeterHit collections to one 
 *  CalorimeterHit collection applying a threshold and an calibration constant...
 */

class SimpleLHCalDigi : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SimpleLHCalDigi ; }
  
  
  SimpleLHCalDigi() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  virtual void end() ;
  
  
 protected:

  int _nRun ;
  int _nEvt ;
  
  std::vector<std::string> _lhcalCollections;

  std::string _outputLHCalCollection;
  std::string _outputRelCollection;

  std::string _cellIDLayerString ;

  std::string  _defaultEncoding;
  std::string  _caloLayout;
  std::string  _caloID;
  //std::string  _caloType;
  
  // jl: replace by string as soon as MarlinUtil is updated
  int _caloType;
  
  float _thresholdLHCal;
  float _calibrCoeffLHCal;
  bool _fixLCalHits; 

  float xing_angle ;
  float zMin       ;
  float dZ         ;
  float rMin       ;
  float cellDimR   ;
  float cellDimPhi ;
  float WThickness ;
 
} ;

#endif



