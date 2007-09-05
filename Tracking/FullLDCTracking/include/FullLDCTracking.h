#ifndef FULLLDCTRACKING_H
#define FULLLDCTRACKING_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include "ClusterExtended.h"
#include "TrackExtended.h"
#include "TrackerHitExtended.h"
#include "TrackHitPair.h"
#include "HelixClass.h"
#include "GroupTracks.h"
#include "../../BrahmsTracking/include/MarlinTrackFit.h"

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
 * number of VTX hits used in the track fit is the first element in this vector  
 * (Track::getSubdetectorHitNumbers()[0]) <br>
 * number of FTD hits used in the track fit is the second element in this vector  
 * (Track::getSubdetectorHitNumbers()[1]) <br>
 * number of SIT hits used in the track fit is the third element in this vector  
 * (Track::getSubdetectorHitNumbers()[2]) <br>
 * number of TPC hits used in the track fit is the forth element in this vector  
 * (Track::getSubdetectorHitNumbers()[3]) <br>
 * total number of VTX hits in track is the fifth element in this vector 
 * (Track::getSubdetectorHitNumbers()[4]) <br>
 * total number of FTD hits in track is the sixth element in this vector
 * (Track::getSubdetectorHitNumbers()[5]) <br>
 * total number of SIT hits in track is the seventh element in this vector
 * (Track::getSubdetectorHitNumbers()[6]) <br>
 * total number of TPC hits in track is the eighth element in this vector
 * (Track::getSubdetectorHitNumbers()[7]) <br>
 * Output track collection has by default a name "LDCTracks". 
 * In addition collection of relations of the tracks to MCParticles is stored if flag CreateMap is set to 1. 
 * Collection of relations has by default a name "LDCTracksMCP" 
 * @param VTXHitCollection name of input VTX TrackerHit collection <br>
 * (default parameter value : "VTXTrackerHits") <br>
 * @param FTDHitCollection name of input FTD TrackerHit collection <br>
 * (default parameter value : "FTDTrackerHits") <br>
 * @param SITHitCollection name of input SIT TrackerHit collection <br>
 * (default parameter value : "SITTrackerHits") <br>
 * @param TPCHitCollection name of input TPC TrackerHit collection <br>
 * (default parameter value : "TPCTrackerHits") <br>
 * @param TPCTracks collection name of TPC tracks <br>
 * (default parameter value : "TPCTracks") <br>
 * @param TPCTracksMCPRelColl Name of input TPC track to MC particle relation collection <br>
 * (default parameter value : "TPCTracksMCP") <br>
 * @param SiTracks collection name of Si tracks <br>
 * (default parameter value : "SiTracks") <br>
 * @param SiTracksMCPRelColl Name of input Si track to MC particle relation collection <br>
 * (default parameter value : "SiTracksMCP") <br> 
 * @param LDCTrackCollection name of the output LDC track collection <br>
 * (default parameter value : "LDCTracks") <br>
 * @param LDCTrackMCPRelCollection name of the output LDC track to MC Particle relation collection <br>
 * (default parameter value : "LDCTracksMCP") <br>
 * @param ReffitedTPCTrackCollection name of the output collection of the refitted TPC tracks <br>
 * (default parameter value : "RefittedTPCTracks") <br>
 * @param RefittedTPCTrackMCPRelCollection name of the output refitted TPC track to MC Particle relation collection <br>
 * (default parameter value : "RefittedTPCTracksMCP") <br>
 * @param ReffitedSiTrackCollection name of the output collection of the refitted Si tracks <br>
 * (default parameter value : "RefittedSiTracks") <br>
 * @param RefittedSiTrackMCPRelCollection name of the output refitted Si track to MC Particle relation collection <br>
 * (default parameter value : "RefittedSiTracksMCP") <br>
 * @param Chi2FitCut cut on the Chi2/Ndf of the track fit <br>
 * (default parameter value : 100.0) <br>
 * @param Chi2PrefitCut cut on the prefit Chi2 of the track candidate, 
 * prefit is done with the simple helix hypothesis <br>
 * (default parameter value : 1e+5) <br>
 * @param CreateMap flag to create relations between Tracks and MCParticles, 
 * if set to 1, relations collection is created and stored in an event <br>
 * (default parameter value : 1) <br>
 * @param AngleCutForMerging  cut on opening angle between 
 * particle momentum reconstructed with TPC and momentum reconstructed
 * with the Silicon detectors. If the opening angle is smaller that this cut
 * the track segment in Silicon trackers and in TPC are tested for their
 * compatibility <br>
 * (default parameter value : 0.10) <br>
 * @param OmegaCutForMerging  cut on the relative difference in the track Omega
 * parameter reconstructed with TPC and with Si detectors. If the relative difference is smaller
 * than this cut, the track segments in TPC and Si are tested for their compatibility <br>
 * (default parameter value : 0.25) <br>
 * @param D0CutForMerging Upper cutoff on the difference in D0 [mm] to allow for merging 
 * of the Si and TPC segments <br>
 * (default parameter value : 500) <br>
 * @param Z0CutForMerging Upper cutoff on the difference in Z0 [mm] to allow for merging
 * of the Si and TPC segments <br>
 * (default parameter value : 1000) <br>
 * @param RefitTPCTracks flag to refit TPC tracks,
 * if set to 1 TPC tracks are refitted <br>
 * (default parameter value : 1) <br>
 * @param RefitSiTracks flag to refit Si Tracks,
 * if set to 1 Si tracks are refitted <br>
 * (default parameter value : 0) <br>
 * @param StoreRefittedTPCTracks flag to store refitted TPC tracks in additional 
 * LCIOTrack collection named RefittedTPCTracks. Corresponding track MCParticle relations
 * are stored in additional collection named RefittedTPCTracksMCP.
 * If set to 1 refitted TPC tracks are stored in the separate collection <br>
 * (default parameter value : 0) <br>
 * @param StoreRefittedSiTracks flag to store refitted Si tracks in additional 
 * LCIOTrack collection named RefittedSiTracks. Corresponding track MCParticle relations
 * are stored in additional collection named RefittedSiTracksMCP.
 * If set to 1 refitted Si tracks are stored in the separate collection <br>
 * (default parameter value : 0) <br>
 * @param Debug flag to allow for printout of debug information,
 * if set to 1 debugging printout is activated
 * (default parameter value : 1) <br>
 * @param OptFit option for track candidate prefitting <br>
 * if OptFit=0 - FORTRAN code tfithl is invoked, <br> 
 * if OptFit=1 - prefit of track candidate is done with the fitting method of ClusterShapes class <br>
 * if OptFit=2 - initial track parameters d0, z0, phi0 and tan(lambda) are taken from the Si 
 * and parameter Omega for the TPC track segments <br>
 * if OptFit=3 - sophisticated iterative prefit, improving cov matrix estimate <br>
 * if OptFit=4 - track parameters are determined from the two separate fits of the Si
 * and TPC track segments (recommended) <br>
 * (default parameter value : 4) <br>
 * @param OptFitTPC option for TPC track refit. Options are the same as for the combined LDC track prefit <br>
 * (default parameter value : 2) <br>
 * @param OptFitSi option for Si track refit. Options are the same as for the combined LDC track prefit <br>
 * (default parameter value : 2) <br>
 * @param UseExtraPoint This flag is used to steer DELPHI fitting code. If set to 1, additional 
 * artificial mesurement point at PCA is introduced with relatively large errors [OBSOLETE] <br>
 * (default parameter value : 0) <br>
 * @param ForceSiTPCMerging This flag steers merging of Si and TPC track segments. If ForceMerging=1
 * Si and TPC track segments are forced to be merged if the opening angle between Si track 
 * momentum and TPC track momentum
 * is less than AngleCutForForcedMerging (see below) and difference in tracks 
 * parameters Omega is less than OmegaCutForForcedMerging (see below) <br>
 * (default parameter value : 0)
 * @param AngleCutForForcedMerging cut on opening angle between Si track momentum and
 * TPC track momentum. Used to steer forced merging of Si and TPC track segments. <br>
 * (default parameter value : 0.05)
 * @param OmegaCutForForcedMerging cut on the difference between Si and TPC tracks parameter
 * Omega. Used to steer forced merging of Si and TPC track segments. Relative 
 * errors are compared. <br>
 * (default parameter value : 0.15) <br>
 * @param D0CutForForcedMerging Upper cutoff on the difference in D0 to allow for forced
 * merging of the Si and TPC segments <br>
 * (default parameter value : 50) <br>
 * @param Z0CutForForcedMerging Upper cutoff on the difference in Z0 to allow for forced
 * merging of the Si and TPC segments <br>
 * (default parameter value : 200) <br>
 * @param ForceTPCSegmentsMerging If this flag is set to 1, the code attempts to 
 * merge TPC segments from the low pt splitted loopers <br>
 * (default parameter value : 1) <br>
 * @param D0CutToMergeTPCSegments cut on the difference in the track parameter
 * d0 [mm] to allow for merging TPC segments <br>
 * (default parameter value : 100) <br>
 * @param Z0CutToMergeTPCSegments cut on the difference in the track parameter
 * z0 [mm] to allow for merging TPC segments <br>
 * (default parameter value : 5000) <br> 
 * @param DeltaPCutMergeTPCSegments cut on the magnitude [GeV/c] of the vectorial difference
 * of the momentum vectors, associated with TPC segments, to allow for merging 
 * of these segments. <br>
 * (default parameter value : 0.1) <br>
 * @param AssignTPCHits If this flag is set to 1, the code attempts to assign left-over 
 * TPC hits to the accepted track candidates. No track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param TPCHitToTrackDistance Cut on the distance between left-over TPC hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 25.0) <br>
 * @param CutOnTPCHits minimal number of TPC hits, used in the track fit, which is 
 * required for tracks which have no hits from the Si detectors <br>
 * (default parameter value : 20) <br> 
 * @param CutOnTrackD0 cut on the d0 parameter of the track. If the d0 parameter is greater that 
 * this cut, track is rejected <br>
 * (default parameter value : 500) <br>
 * @param CutOnTrackZ0 cut on the z0 parameter of the track. If the z0 parameter is greater that 
 * this cut, track is rejected <br>
 * (default parameter value : 500) <br>
 * @param ForbidOverlapInZTPC If this flag is set to 1 then merging of the TPC semiloops is 
 * forbiden for segment overlapping in z <br>
 * (default parameter value : 0) <br>
 * @param ForbidOverlapInZComb If this flag is set to 1 then merging of left-over TPC semiloop and
 * combined Si-TPC track is their segments overlap in z <br>
 * (default parameter value : 0) <br>
 * @param aParameterForIPError parameter a defining minimal IP resolution
 * according to the formular sigma[IP] = a + b/[P*sin^3/2{Q}]^s where P is the particle momentum
 * and Q is the polar angle. If the resolution on IP, calculated with DELPHI fitting routine,
 * is smaller than minimal allowed resolution, the corresponding track parameter 
 * covariance matrix element is set to the minimal allowed value <br>
 * (default is 0.002 [mm]) <br>
 * @param bParameterForIPError parameter b in the parametrisation of the minimal IP resolution
 * sigma[IP] = a + b/[P*sin^3/2{Q}]^s <br>
 * (default is 0.0076 [mm]) <br>
 * @param sParameterForIPError parameter s in the parametrisation of the minimal IP resolution
 * sigma[IP] = a + b/[P*sin^3/2{Q}]^s <br>
 * (default is 0.75 [mm]) <br>
 * @param StoreHitsInFit if set to 1 only hits used in the track fit are stored in the 
 * corresponding associated vector of TrackerHits <br>
 * (default is 0) <br>
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
 protected:

  void prepareVectors( LCEvent * evt );
  void CleanUp();
  void MergeTPCandSiTracks();
  TrackExtended * CombineTracks(TrackExtended * tpcTrk, TrackExtended * siTrk);
  void Sorting(TrackExtendedVec & trackVec);
  void SelectCombinedTracks();
  void AddNotCombinedTracks();
  void AddNotAssignedHits();
  void AddTrackColToEvt(LCEvent * evt, TrackExtendedVec & trkVec, 
			std::string TrkColName, std::string RelColName);
  float CompareTrk(TrackExtended * first, TrackExtended * second, 
		     float d0Cut, float z0Cut, int iopt);
  
  float CompareTrkII(TrackExtended * first, TrackExtended * second, 
		     float d0Cut, float z0Cut, int iopt, float &Angle);

  void SortingTrackHitPairs(TrackHitPairVec & trackHitPairVec);

  void AssignSiHitsToTracks(TrackerHitExtendedVec hitVec,
			    float dcut);

  void AssignTPCHitsToTracks(TrackerHitExtendedVec hitVec,
			     float dcut);

  void PrintOutMerging(TrackExtended * firstTrackExt, TrackExtended * secondTrackExt, 
		       int iopt);

  int _nRun ;
  int _nEvt ;

  std::string _TPCTrackCollection;
  std::string _SiTrackCollection;
  std::string _TPCTrackMCPCollName;
  std::string _SiTrackMCPCollName;
  std::string _VTXTrackerHitCollection;
  std::string _SITTrackerHitCollection;
  std::string _FTDTrackerHitCollection;
  std::string _TPCTrackerHitCollection;
  std::string _LDCTrackCollection;
  std::string _LDCTrackMCPCollection;
  std::string _RefittedTPCTrackCollection;
  std::string _RefittedTPCTrackMCPCollection;
  std::string _RefittedSiTrackCollection;
  std::string _RefittedSiTrackMCPCollection;

  TrackExtendedVec _allSiTracks;
  TrackExtendedVec _allTPCTracks;
  TrackExtendedVec _allCombinedTracks;
  TrackExtendedVec _allNonCombinedTPCTracks;
  TrackExtendedVec _allNonCombinedSiTracks;
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

  float _bField;
  float _chi2PrefitCut;
  float _chi2FitCut;

  int _debug;
  int _createMap;
  int _useExtraPoint,_optFit;

  float _dPCutForMerging;
  float _d0CutForMerging;
  float _z0CutForMerging;
  float _dOmegaForMerging;
  float _angleForMerging;


  int _forceMerging;
  float _dPCutForForcedMerging;
  float _d0CutForForcedMerging;
  float _z0CutForForcedMerging;
  float _dOmegaForForcedMerging;
  float _angleForForcedMerging;


  int _mergeTPCSegments;
  float _dPCutToMergeTPC;
  float _d0CutToMergeTPC;
  float _z0CutToMergeTPC;

  MarlinTrackFit _trackFit;

  int _refitTPCTracks;
  int _refitSiTracks;
  int _storeRefittedTPCTracks;
  int _storeRefittedSiTracks;
  int _storeHitsInFit;

  int _cutOnTPCHits;

  float _aParIpReso,_bParIpReso,_sParIpReso;

  int _assignVTXHits,_assignFTDHits,_assignSITHits,_assignTPCHits;

  float _distCutForVTXHits,_distCutForFTDHits,_distCutForSITHits,_distCutForTPCHits;

  int _optFitTPC,_optFitSi;

  float _d0TrkCut,_z0TrkCut;

  int _forbidOverlapInZTPC,_forbidOverlapInZComb;

  LCEvent * _evt;

} ;

#endif



