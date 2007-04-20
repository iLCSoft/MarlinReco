/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef FTDDigiProcessor_h
#define FTDDigiProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

using namespace lcio ;
using namespace marlin ;


/** ======= FTDigiProcessor ========== <br>
 * Produces FTD TrackerHit collection from SimTrackerHit collection. <br> 
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits in (x,y) plane according to the specified point resolution. 
 * Each created TrackerHit is assigned the type via method TrackerHitImpl::setType(int type).
 * The FTD TrackerHit type is encoded in the following way : <br>
 * type = 200 + layer_index_ftd  for vertex hits (layer_index_ftd = 1...7) <br>
 * To access this type use method TrackerHit::getType() <br> 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collections of SimTrackerHits in FTD <br>
 * <h4>Output</h4>
 * Processor produces collection of digitized TrackerHits in FTD <br>
 * @param CollectionName The name of input collection of FTD SimTrackerHits <br>
 * (default name ftd01_VXD) <br>
 * @param OutputCollectionName The name of output collection of TrackerHits <br>
 * (default name FTDTrackerHits) <br> 
 * @param PointResolution Point resolution in (x,y) for the FTD planar detectors (in mm) <br>
 * (default value 0.01) <br>
 * @param RemoveDrays When this flag is set to 1 hits produced by delta-electrons are removed
 * from output collections <br>
 * (default value 0) <br>
 * @param MomentumCutForDRays The upper cut on delta-electron momentum (in MeV) <br>
 * (default value 10) <br>
 * <br>
 * @author A. Raspereza, MPI (Munich)
 */


class FTDDigiProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new FTDDigiProcessor ; }
  
  
  FTDDigiProcessor() ;
  
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

  int _nRun ;
  int _nEvt ;
  float _pointReso;
  float _momCut;

  int _removeDrays;
  
  // gsl random number generator
  gsl_rng * r ;

} ;

#endif



