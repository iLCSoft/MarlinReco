#ifndef TOFEstimators_h
#define TOFEstimators_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

using namespace lcio ;
using namespace marlin ;


class TH2F ;

/** Compute various estimators for the time of flight from the CalorimeterHits in Clusters.
 *  
 * 
 * @author F.Gaede, DESY, April 2018
 */

class TOFEstimators : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TOFEstimators ; }
  
  TOFEstimators(const TOFEstimators&) = delete;
  TOFEstimators& operator=(const TOFEstimators&) = delete;

  TOFEstimators() ;
  
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

  /** Input collection name with ReconstructedParticles (PFOs).
   */
  std::string _colNamePFO{};
  int   _maxLayerNum{} ;
  float _resolution{} ;
  
  std::vector<std::string> _TOFNames {} ; 

  int _nRun{};
  int _nEvt{};

  gsl_rng* _rng = nullptr ;

  std::vector<TH2F*> _h{};

} ;

#endif



