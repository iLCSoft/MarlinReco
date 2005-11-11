#ifndef FORMTRUECLUSTERSKP_H
#define FORMTRUECLUSTERSKP_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include "HelixClass.h"
#include <EVENT/CalorimeterHit.h>
#include <EVENT/MCParticle.h>

using namespace lcio ;
using namespace marlin ;


/** Cluster Cheater <br>
 *  This processor constructs true clusters.<br>
 *  All the hits are collected. <br>
 *  Uses gear to get inner radius and z of ecal. <br>
 *    @author P. Krstonosic (DESY)<br>
 *    @version $ld: $<br>
 */
class ClusterCheater : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ClusterCheater ; }
  
  
  ClusterCheater() ;

  /** Initialization
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


  std::string _trueClustCollection;
  std::vector<std::string> _caloCollections;
  std::string _relCollection;
  std::string _trueClustToMCP;


} ;

#endif



