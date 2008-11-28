/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef ETDDigiProcessor_h
#define ETDDigiProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

using namespace lcio ;
using namespace marlin ;


/** ======= ETDigiProcessor ========== <br>
 * Produces ETD TrackerHit collection from SimTrackerHit collection. <br> 
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits in (x,y) plane according to the specified point resolution. 
 * Each created TrackerHit is assigned the type via method TrackerHitImpl::setType(int type).
 * The ETD TrackerHit type is encoded in the following way : <br>
 * type = 200 + layer_index_ftd  for vertex hits (layer_index_ftd = 1...7) <br>
 * To access this type use method TrackerHit::getType() <br> 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collections of SimTrackerHits in ETD <br>
 * <h4>Output</h4>
 * Processor produces collection of digitized TrackerHits in ETD <br>
 * @param CollectionName The name of input collection of ETD SimTrackerHits <br>
 * (default name ftd01_VXD) <br>
 * @param OutputCollectionName The name of output collection of TrackerHits <br>
 * (default name ETDTrackerHits) <br> 
 * @param PointResolution Point resolution in (x,y) for the ETD planar detectors (in mm) <br>
 * (default value 0.01) <br>
 * @param RemoveDrays When this flag is set to 1 hits produced by delta-electrons are removed
 * from output collections <br>
 * (default value 0) <br>
 * @param MomentumCutForDRays The upper cut on delta-electron momentum (in MeV) <br>
 * (default value 10) <br>
 * <br>
 * 
 * [F.Gaede: this is a one to one copy of the original FTDDigiProcessor]
 * F.Gaede: 2008-11-28
 *    added parameter ActiveLayers: only hits from these layers will be digitized
 *    -> used to mimic the stereo layers strip detectors (eg. use only hits from layers -2,2)
 * 
 * @author A. Raspereza, MPI (Munich)
 */


class ETDDigiProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ETDDigiProcessor ; }
  
  
  ETDDigiProcessor() ;
  
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
  std::string _colName ;
  std::string _outColName ;

  std::vector<int> _activeLayers ;

  int _nRun ;
  int _nEvt ;
  float _pointReso;
  float _momCut;

  int _removeDrays;
  
  // gsl random number generator
  gsl_rng * r ;

} ;

#endif



