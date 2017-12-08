#ifndef NNClusterProcessor_h
#define NNClusterProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;


/** Example processor that does a simple nearest neighbour (NN) clustering on one or more CalorimeterHit 
 *  collections. It uses a simple euclididan distance cut.
 * 
 * @param HitCollections    - Name of the input collection(s) (CalorimeterHit)
 * @param OutputCollection  - Name of the output collection (Cluster)
 * @param DistanceCut       - Cut for distance between hits in mm
 * @param EnergyCut         - Cut for hit energy in GeV
 *
 *  @author F.Gaede (DESY)
 *  @version $Id$
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
  StringVec _colNames{};

  std::string _outputColName{};

  float _distCut{};
  float _eCut{};

  int _nThetaPhi{};

  int _nRun{};
  int _nEvt{};

//   NNClusterer* _clusterer ;

} ;

#endif



