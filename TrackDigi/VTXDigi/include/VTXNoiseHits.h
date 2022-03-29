#ifndef VTXNoiseHits_h
#define VTXNoiseHits_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>


using namespace lcio ;
using namespace marlin ;


/**
\addtogroup TrackDigi TrackDigi
@{
\addtogroup VTXNoiseHits VTXNoiseHits
@{
Adds random noise hits to collection of TrackerHits of the vertex detector.
======= VTXNoiseHits ========== <br>
 * Adds random noise hits to collection of TrackerHits of the vertex detector.
 * The number of noise hits are given by the parameter HitDensityPerLayer (hits/cm^2).
 * The noise hits are created with a uniform distribution over the ladder
 * surface.
 * 
 * @param HitDensityPerLayer  hit densities (hits/cm^2) per VXD layer
 * @param PointResolutionRPhi_VTX Point resolution in r-phi for the vertex detector (in mm) <br>
 * (default value 0.0027) <br>
 * @param PointResolutionZ_VTX Point resolution in z for the vertex detector (in mm) <br>
 * (default value 0.0027) <br>
 * 
 * <br>
 * @version $Id$
 * @author F.Gaede, DESY
 */
class VTXNoiseHits : public Processor {
  
 public:
  VTXNoiseHits(const VTXNoiseHits&) = delete;
  VTXNoiseHits& operator=(const VTXNoiseHits&) = delete; 
 
  virtual Processor*  newProcessor() { return new VTXNoiseHits ; }
  
  
  VTXNoiseHits() ;
  
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

  std::string _colNameVTX{};
  FloatVec _densities{};
  float _pointResoRPhiVTX{};
  float _pointResoZVTX{};

  int _nRun{};
  int _nEvt{};

  gsl_rng * r{};


} ;
/** @} @} */
#endif
