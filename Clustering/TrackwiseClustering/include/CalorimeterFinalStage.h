#ifndef CalorimeterFinalStage_h
#define CalorimeterFinalStage_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;



/** Computes cluster properties from clustered hits.
 * Reads a collection of type Cluster and creates a new Cluster
 * collection that holds the final clusters created from the internediate 
 * input clusters. Could for example be used if another clustering algorithm 
 * than Trackwise Clustering had been used to find cluster the hits.
 * 
 *  @author A. Raspereza, DESY 
 *  @version $Id: CalorimeterFinalStage.h,v 1.1 2005-08-07 15:57:07 gaede Exp $ 
 */
class CalorimeterFinalStage : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CalorimeterFinalStage ; }
  
  
  CalorimeterFinalStage() ;
  
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

  std::string _clusterInput;
  std::string _clusterOutput;
  int _nhit_minimal;
  int _nRun ;
  int _nEvt ;
} ;

#endif



