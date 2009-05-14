#ifndef VTXNoiseClusters_h
#define VTXNoiseClusters_h 1

#include "marlin/Processor.h"
#include "marlin/EventModifier.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>


using namespace lcio ;
using namespace marlin ;

class VXDGeometry ;


/** ======= VTXNoiseClusters ========== <br>
 * Adds random noise hits to collection of SimTrackerHits of the vertex detector.
 * The number of noise hits are given by the parameter HitDensityPerLayer (hits/cm^2).
 * The noise hits are created with a uniform distribution over the ladder
 * surface. An object of type VXDClusterParameters is added to every hit to describe the 
 * extension of the cluster on the ladder surface. 
 * 
 * 
 * @param HitDensityPerLayer  hit densities (hits/cm^2) per VXD layer
 * @param PointResolutionRPhi_VTX Point resolution in r-phi for the vertex detector (in mm) <br>
 * (default value 0.0027) <br>
 * @param PointResolutionZ_VTX Point resolution in z for the vertex detector (in mm) <br>
 * (default value 0.0027) <br>
 * 
 * <br>
 * @version $Id: VTXNoiseClusters.h,v 1.1 2009-05-14 07:24:19 gaede Exp $
 * @author F.Gaede, DESY
 */

class VTXNoiseClusters : public  Processor, public EventModifier{
  
 public:
  
  virtual Processor*  newProcessor() { return new VTXNoiseClusters ; }
  
  
  VTXNoiseClusters() ;
  
  virtual const std::string & name() const { return Processor::name() ; }

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  void modifyEvent( LCEvent * evt ) ; 
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;


 protected:

  std::string _colNameVTX ;
  FloatVec _densities ;
//   float _pointResoRPhiVTX ;
//   float _pointResoZVTX ;
  int   _ranSeed  ;

  int _nRun ;
  int _nEvt ;

  gsl_rng* _rng ;

  VXDGeometry* _vxdGeo ;


} ;

#endif



