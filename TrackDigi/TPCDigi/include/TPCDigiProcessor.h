/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef TPCDigiProcessor_h
#define TPCDigiProcessor_h 1

#include <marlin/Processor.h>
#include <lcio.h>
#include <string>
#include <gsl/gsl_rng.h>


using namespace lcio ;
using namespace marlin ;



/** ====== TPCDigiProcessor ====== <br>
 * Caution: This digitiser presently does not process space-point like SimTrackerHits which have been flagged with CellIDs set to the negetive row number. This must be implemented in future. 
 *Produces TPC TrackerHit collection from SimTrackerHit collection, smeared in r-phi and z. 
 * Double hits are identified but are currently not added to the collection. This may be change 
 * at a later date when criteria for their seperation is defined. The resolutions are defined in 
 * the GEAR stearing file. Resolution in r-phi is calculated according to the formular <br>
 * sigma(r-phi) = sqrt(const**2 + diffusion*z_drift) <br>
 * where 'const' stands for the constant term and 'diffusion' stands for the diffusion term,
 * 'z_drift' is the drift length <br>
 * At the moment resolution in z assumed to be independent of drift length. <br>
 * The type of TPC TrackerHit is set to 500 via method TrackerHitImpl::setType(int type) <br>
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collections of SimTrackerHits in TPC <br>
 * <h4>Output</h4>
 * Processor produces collection of digitized TrackerHits in TPC <br>
 * @param CollectionName The name of input SimTrackerHit collection <br>
 * (default name STpc01_TPC)
 * @param TPCTrackerHitsCol The name of output collection of TrackerHits <br>
 * (default name TPCTrackerHits) <br>
 * <br>
 * @authors S. Aplin, DESY and A.Raspereza, MPI
 *
 * Changed 7/9/07 so that the const and diffusion resolution terms are taken as processor parameters rather than the gear file.
 * The parameters _pixZ and pixRP were also changed from gear parameters to processor parameters
 * clare.lynch@bristol.ac.uk
 *
 */
class TPCDigiProcessor : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new TPCDigiProcessor ; }
  
  
  TPCDigiProcessor() ;
  
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
  std::string _TPCTrackerHitsCol ;

  int _nRun ;
  int _nEvt ;

  // gsl random number generator
  gsl_rng * _random ;



  float _pointResoRPhi;
  float _pointResoZ;
  float _diffRPhi;
  float _pixZ;
  float _pixRP;
} ;

#endif



