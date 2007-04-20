#ifndef FULLLDCTRACKING_H
#define FULLLDCTRACKING_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include "ClusterExtended.h"
#include "TrackExtended.h"
#include "TrackerHitExtended.h"
#include "HelixClass.h"
#include "GroupTracks.h"
#include "MarlinTrackFit.h"

using namespace lcio ;
using namespace marlin ;

/** === FullLDCTracking Processor === <br>
 * Processor performing track finding procedure in 
 * the entire LDC detector by linking track segments
 * found by the SiliconTracking module in the silicon detectors
 * and by the LEPTracking module in TPC. 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collection of digitized vertex, sit, ftd & tpc tracker hits 
 * and also the collections of tracks found in the silicon detectors
 * and in TPC.
 * <h4>Output</h4>
 * Processor produces an LCIO collection of the Tracks. Each track is characterised by 
 * five parameters : Omega (signed curvuture), Tan(lambda) where
 * lambda is the dip angle, Phi (azimuthal angle @ point of closest approach), D0 (signed impact parameter),
 * Z0 (displacement along z axis at the point of closest approach to IP). 
 * Covariance matrix for these parameters is also provided.
 * Only lower left corner of the covariance matrix is stored. The sequence of the covariance matrix elements 
 * assigned to track is the following: <br>
 * (D0,D0) <br>
 * (Phi,D0), (Phi,Phi) <br>
 * (Omega,D0), (Omega,Phi), (Omega,Omega) <br>
 * (Z0,D0), (Z0,Phi), (Z0,Omega), (Z0,Z0) <br>
 * (TanL,D0), (TanL,Phi), (TanL,Omega), (TanL,Z0), (TanL,TanL) <br>
 * The number of hits in the different subdetectors associated
 * with each track can be accessed via method Track::getSubdetectorHitNumbers().
 * This method returns vector of integers : <br>
 * number of VXD hits in track is the first element in this vector  
 * (Track::getSubdetectorHitNumbers()[0]) <br>
 * number of FTD hits in track is the second element in this vector  
 * (Track::getSubdetectorHitNumbers()[1]) <br>
 * number of SIT hits in track is the third element in this vector  
 * (Track::getSubdetectorHitNumbers()[2]) <br>
 * number of TPC hits in track is the fourth element in this vector  
 * (Track::getSubdetectorHitNumbers()[3]) <br>
 * Output track collection has by default a name "LDCTracks". 
 * In addition collection of relations of the tracks to MCParticles is stored if flag CreateMap is set to 1. 
 * Collection of relations has by default a name "LDCTracksMCP" 
 * @param VXDHitCollectionName name of input VXD TrackerHit collection <br>
 * (default parameter value : "VTXTrackerHits") <br>
 * @param FTDHitCollectionName name of input FTD TrackerHit collection <br>
 * (default parameter value : "FTDTrackerHits") <br>
 * @param SITHitCollectionName name of input SIT TrackerHit collection <br>
 * (default parameter value : "SITTrackerHits") <br>
 * @param TPCHitCollectionName name of input TPC TrackerHit collection <br>
 * (default parameter value : "TPCTrackerHits") <br>
 * @param TPCTracks collection name of TPC tracks <br>
 * (default parameter value : "TPCTracks") <br>
 * @param SiTracks collection name of Si tracks <br>
 * (default parameter value : "SiTracks") <br>
 * @param Chi2FitCut cut on the Chi2/ndf of the track fit <br>
 * (default parameter value : 100.0) <br>
 * @param Chi2PrefitCut cut on the prefit Chi2 of track candidate, 
 * prefit is done with the simple helix hypothesis <br>
 * (default parameter value : 1e+10) <br>
 * @param CreateMap flag to create relations between Tracks and MCParticles, 
 * if set to 1, relations collection is created and stored in event <br>
 * (default parameter value : 1)
 * @param AngleCutForMerging  cut on opening angle between 
 * particle momentum reconstructed with TPC and momentum reconstructed
 * with the Silicon detectors. If the opening angle is smaller that this cut
 * the track segment in Silicon trackers and in TPC are tested for their
 * compatibility <br>
 * (default parameter value : 0.2) <br>
 * @param RefitTPCTracks flag to refit TPC tracks,
 * if set to 1 TPC tracks are refitted
 * (default parameter value : 0)
 * @param RefitSiTracks flag to refit Si Tracks,
 * if set to 1 Si tracks are refitted
 * (default parameter value : 0)
 * @param Debug flag to allow for printout of debug information,
 * if set to 1 debugging printout is done
 * (default parameter value : 0)
 * @param OptPrefit option for track candidate prefitting <br>
 * if OptPrefit=0 - FORTRAN code tfithl is invoked, <br> 
 * if OptPrefit=1 - prefit of track candidate is done with the fitting method of ClusterShapes class <br>
 * (default parameter value : 1) <br>
 * @param UseExtraPoint This flag is used to steer DELPHI fitting code. If set to 1, additional 
 * artificial mesurement point at PCA is introduced with relatively large errors. <br>
 * (default parameter value : 0)
 * @param ForceMerging This flag steers merging of Si and TPC track segments. If ForceMerging=1
 * Si and TPC track segments are forced to be merged if the opening angle between Si track 
 * momentum and TPC track momentum
 * is less than AngleCutForForcedMerging (see below) and difference in tracks 
 * parameters Omega is less than OmegaCutForForcedMerging (see below) <br>
 * (default parameter value : 1)
 * @param AngleCutForForcedMerging cut on opening angle between Si track momentum and
 * TPC track momentum. Used to steer forced merging of Si and TPC track segments. <br>
 * (default parameter value : 0.07)
 * @param OmegaCutForForcedMerging cut on the difference between Si and TPC tracks parameter
 * Omega. Used to steer forced merging of Si and TPC track segments. <br>
 * (default parameter value : 0.1)
 * <br>
 * @author A. Raspereza (MPI Munich)<br>
 */

class FullLDCTracking : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new FullLDCTracking ; }  
  FullLDCTracking() ;  
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  void prepareVectors( LCEvent * evt );
  void CleanUp();
  void MergeTPCandSiTracks();
  TrackExtended * CombineTracks(TrackExtended * tpcTrk, TrackExtended * siTrk);
  void Sorting(TrackExtendedVec & trackVec);
  void SelectCombinedTracks();
  void AddNotCombinedTracks();
  void AddTrackColToEvt(LCEvent * evt, TrackExtendedVec & trkVec, 
			std::string TrkColName, std::string RelColName);

 protected:

  int _nRun ;
  int _nEvt ;

  std::string _TPCTrackCollection;
  std::string _SiTrackCollection;
  std::string _VTXTrackerHitCollection;
  std::string _SITTrackerHitCollection;
  std::string _FTDTrackerHitCollection;
  std::string _TPCTrackerHitCollection;

  TrackExtendedVec _allSiTracks;
  TrackExtendedVec _allTPCTracks;
  TrackExtendedVec _allCombinedTracks;
  TrackExtendedVec _trkImplVec;
  TrackerHitExtendedVec _allTPCHits;
  TrackerHitExtendedVec _allVTXHits;
  TrackerHitExtendedVec _allFTDHits;
  TrackerHitExtendedVec _allSITHits;

  float _resolutionRPhi_TPC,_resolutionZ_TPC;
  float _resolutionRPhi_VTX,_resolutionZ_VTX;
  float _resolutionRPhi_FTD,_resolutionZ_FTD;
  float _resolutionRPhi_SIT,_resolutionZ_SIT;
  float PI, PIOVER2, TWOPI;
  float _deltaPhiForTracks,_deltaQForTracks;
  float _angleCutForMerging;
  float _bField;
  float _chi2PrefitCut;
  float _chi2FitCut;
  int _debug;
  int _createMap;
  int _useExtraPoint,_optFit;
  int _forceMerging;
  float _angleCutForForcedMerging;
  float _omegaCutForForcedMerging;
  MarlinTrackFit _trackFit;

  int _refitTPCTracks;
  int _refitSiTracks;
  int _storeRefittedTPCTracks;
  int _storeRefittedSiTracks;

} ;

#endif



