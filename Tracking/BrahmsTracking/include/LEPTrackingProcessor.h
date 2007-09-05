#ifndef LEPTrackingProcessor_h
#define LEPTrackingProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;



/** Produces Track collection from TPC TrackerHit collections using LEP tracking algorithms. 
* The Geometric information via GEAR. The f77 tracking code previously relied on ZEBRA banks, 
* these have been replace by C++ structures implemented in tkhitbank.h etc. There must be 
* instanitated and deleted at the begining and the end of the processor respectively.
*
* Currently only the TPC patrec and fitting is implemented, though plans to incorporate VTX
* hits into Tracks exist.
*
* The name of the TPC Track collection is specified in the steering file, togeather with the name 
* of the hit collections needed as input.
*
* For the track collection the following applys:
*
* The reference point is the Point of Closest Approach.
*
* All parameters are defined there.
*
* Phi is defined between -PI and +PI
*
* Z0 is the z coordinate of the PCA in the R-Phi plane.
*
* For the relations, the weighs are calculated as the percentage of hits that a given MC particle contributes to the reconstucted track's hit collection
*
* At present the covariance matrix is not filled.
*
*
* The LEP algorithms are taken from Brahms but does not contain the DELPHI ambiguity resolver. 
* The track finding is based on  
*  1) SUBROUTINE CIRCLE  (N.CHERNOV, G. OSOSKOV )  
*     REFERENCE:  COMPUTER PHYSICS COMMUNICATIONS VOL 33,P329   
*  2) 3-DIMENSIONAL ITERACTION  (MARTIN POPPE)     
*     REFERENCE:  ALEPH NOTE 87-102
*
* The final track fitting is done using a Kalman filter. At present the both the magnetic field 
* and the material desciption of the TPC is hard coded into the Fortran code and C++ code, 
* and will be improved when GEAR becomes more evolved.
* <h4>Input collections and prerequisites</h4> 
* Processor requires collection of digitized TPC tracker hits. 
* If such a collections with the user specified names do not exist 
* processor takes no action. Processor still attempts to assign 
* VTX and SIT tracker hits to the found TPC tracks and produce 
* combined Si-TPC tracks. Hence, optionally an user can provide 
* the names of the VTX and SIT tracker hit collections. A more 
* efficient algorithm of combining
* information from TPC and silicon detectors <br>
* <h4>Output</h4>
* Processor produces collection of TPC tracks<br>
* @param TPCTrackerHitCollectionName Name of the TPC TrackerHit collection <br>
* (default value is TPCTrackerHits) <br>
* @param VTXTrackerHitCollectionName Name of the VTX TrackerHit collection <br>
* (default value is VTXTrackerHits) <br>
* @param SITTrackerHitCollectionName Name of the SIT TrackerHit collection <br>
* (default value is SITTrackerHits) <br>
* @param TPCTrackCollectionName Name of the output TPC Track collection <br>
* (default value is TPCTracks) <br>
* @param TrackCollectionName Name of the combined Si-TPC Track collection <br>
* (default value is Tracks) <br>
* @param MCTPCTrackRelCollectionName Name of the TPC Track MC Relation collection <br>
* (default value is TPCTracksMCP) <br>
* @param MCTrackRelCollectionName Name of the Track MC Relation collection <br>
* (default value is TracksMCP) <br>
* @authors S. Aplin, DESY, A. Raspereza MPI-Munich    
*/
class LEPTrackingProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new LEPTrackingProcessor ; }
  
  
  LEPTrackingProcessor() ;
  
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

  /** Input collections name.
   */
  std::string _colNameTPC ;
  std::string _colNameVTX ;
  std::string _colNameSIT ;
  std::string _colNameTPCTracks ;
  std::string _colNameTracks ;
  std::string _colNameMCTPCTracksRel ;
  std::string _colNameMCTracksRel ;

  int _nRun ;
  int _nEvt ;
} ;

#endif



