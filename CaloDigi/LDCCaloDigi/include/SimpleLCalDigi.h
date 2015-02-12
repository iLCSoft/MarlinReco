
#ifndef DIGITIZERLCAL_H
#define DIGITIZERLCAL_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;


/** === SimpleLCalDigi Processor === <br>
 *  Simple calorimeter digitizer for the LCal Processor. <br>
 *  Simple calorimeter digitizer for the muon detectors.
 *  Converts SimCalorimeterHit collections to one 
 *  CalorimeterHit collection applying a threshold and an calibration constant...
 * 
 *  @version $Id$
 */

class SimpleLCalDigi : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SimpleLCalDigi ; }
  
  
  SimpleLCalDigi() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  virtual void end() ;
  
  
 protected:

  int _nRun ;
  int _nEvt ;
  
  std::vector<std::string> _lcalCollections;

  std::string _outputLCalCollection;
  std::string _outputRelCollection;

  std::string _cellIDLayerString ;

  float _thresholdLCal;
  float _calibrCoeffLCal;


} ;

#endif



