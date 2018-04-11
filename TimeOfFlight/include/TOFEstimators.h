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

/** Compute estimators for the time of flight from the CalorimeterHits in Clusters.
 *  The estimators are stored in a PID object with the name of the processor.
 *  The following parameters are stored:
 *   TOFFirstHit         -  first hit ( closest to calo entry point)
 *   TOFClosestHits      -  closest hit in every layer (<lMax)
 *   TOFClosestHitsError -  error of above 
 *   TOFFlightLength     -  trajectory length to reference point
 *   TOFLastTrkHit       -  (unsmeared) time of last tracker hit (in SET)
 *   TOFLastTrkHitFlightLength -  trajectory length to last trk hit
 *
 *  processor parameters:
 *  MaxLayerNumber  - restrict the hits to the first MaxLayerNumbers
 *  ReconstructedParticleCollection - input collection
 *  TimeResolution  - assumed single hit resolution in ps
 *
 *  Only Ecal hits are considered.
 *  For charged particles the TrackState at the calorimeter is used as reference point and for the 
 *  extrapolation into the calorimeter. For neutral particles the position of the first hit (the one closest to the IP) 
 *  is used and a straight flight path from the IP is taken for the extrapolation.
 *  The hit time is corrected for the flight time from the reference assuming speed of light.
 *  
 *  See ../scripts/tofestimators.xml for example steering file.
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



