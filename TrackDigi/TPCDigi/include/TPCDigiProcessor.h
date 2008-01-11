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
 * @param RejectCellID0 Whether or not to reject SimTrackerHits with Cell ID 0. Mokka drivers
 * TPC00-TPC03 encoded the pad row number in the cell ID, which should always be non-zero anyway.
 * Drivers TPC04 and TPC05 do not simulate pad rows and thus have the cell ID set to zero for all hits.
 * You will need to set RejectCellID0 to 0 in order to use this processor with these drivers, but note
 * that the implications for track reconstruction are not strictly defined. Mokka driver TPC06 uses
 * a mixed approach with one hit per pad row having non-zero cell ID, extra hits having 0 cell ID.
 * Typically, unless you use TPC04 or TPC05, you should not touch this parameter. <br>
 * (default value 1)
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

  int _rejectCellID0;

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



