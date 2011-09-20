#ifndef VTXNoiseClusters_h
#define VTXNoiseClusters_h 1


//#ifdef USE_ROOT 
// we need some root histograms with cluster size distribution for this processor

#include "marlin/Processor.h"
#include "marlin/EventModifier.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

#include "TH2F.h"
#include "TFile.h"

#ifdef MARLIN_USE_AIDA
#include <AIDA/AIDA.h>
typedef std::vector< AIDA::IHistogram1D* > Hist1DVec ;
typedef std::vector< AIDA::IHistogram2D* > Hist2DVec ;
#endif

using namespace lcio ;
using namespace marlin ;

class VXDGeometry ;


/** ======= VTXNoiseClusters ========== <br>
 * Adds random noise hits to collection of SimTrackerHits of the vertex detector.
 * The number of noise hits are given by the parameter HitDensityPerLayer (hits/cm^2).
 * The noise hits are created with a uniform distribution over the ladder
 * surface. An object of type VXDClusterParameters is added to every hit to describe the 
 * extension of the cluster on the ladder surface. 
 * The distribution of the cluster sizes is read from ROOT histograms - one per layer.
 * 
 * @param RootHistograms      root file name and histogram names (one per layer)
 * @param HitDensityPerLayer  hit densities (hits/cm^2) per VXD layer
 * @param VTXCollectionName   collection of VXD SimTrackerhits
 * @param RandomSeed          random number seed
 * 
 * <br>
 * @version $Id$
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
  StringVec _rootNames ;  
//   float _pointResoRPhiVTX ;
//   float _pointResoZVTX ;
  int   _ranSeed  ;

  int _nRun ;
  int _nEvt ;

  gsl_rng* _rng ;
  VXDGeometry* _vxdGeo ;

  std::vector<TH2F*> _hist ;
  TFile* _hfile ;

#ifdef MARLIN_USE_AIDA
   Hist1DVec _hist1DVec ;
   Hist2DVec _hist2DVec ;
#endif

} ;

#endif
//#endif // USE_ROOT



