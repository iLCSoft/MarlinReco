#ifndef NNClusterProcessor_h
#define NNClusterProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;



/** Simple nearest neighbour clustering ....
 */
class NNClusterProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new NNClusterProcessor ; }
  
  
  NNClusterProcessor() ;
  
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

  /** Input collection name.
   */
  StringVec _colNames ;

  std::string _outputColName ;

  float _distCut ;
  float _eCut ;

  int _nThetaPhi ;

  int _nRun ;
  int _nEvt ;

//   NNClusterer* _clusterer ;

} ;

#endif



