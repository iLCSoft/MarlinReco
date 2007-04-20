#ifndef SILICONTRACKING_H
#define SILICONTRACKING_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <IMPL/TrackImpl.h>
#include "ClusterExtended.h"
#include "TrackExtended.h"
#include "TrackerHitExtended.h"
#include "HelixClass.h"
#include "MarlinTrackFit.h"

using namespace lcio ;
using namespace marlin ;


/** === Silicon Tracking Processor === <br>
 * Processor performing stand-alone pattern recognition
 * in the vertex detector (VXD), forward tracking disks and SIT. <br>
 * The procedure consists of three steps : <br> 
 * 1) Tracking in VXD and SIT ; <br>
 * 2) Tracking in FTD ; <br>
 * 3) Merging compatible track segments reconstructed in VXD and FTD <br>
 * STEP 1 : TRACKING IN VXD and SIT <br>
 * Algorithm starts with finding of hit triplets satisfying helix hypothesis <br> 
 * in three different layers. Two layers of SIT are effectively considered as outermost <br>
 * layers of the vertex detector. To accelerate procedure, the 4-pi solid angle
 * is divided in NDivisionsInTheta and NDivisionsInPhi sectors in cosQ and Phi, 
 * respectively. Triplets are looked for in 2x2 window of adjacent sectors. 
 * Once triplet is found, attempt is made to associate additional hits to 
 * track. Combinatin of hits is accepted for further analysis if the Chi2 
 * of the fit is less than certain predefined threshold. All accepted 
 * combinations are sorted in ascending order of their Chi2. First track candidate 
 * in the sorted array is automatically accepted. The hits belonging to this track are 
 * marked as used, and track candidates sharing these hits are discarded.
 * The procedure proceeds with increasing index of track candidate in the sorted 
 * array until all track candidate have been output or discarded. <br>
 * STEP 2 : TRACKING IN FTD <br>
 * In the next step tracking in FTD is performed. The strategy of tracking in the FTD 
 * is the same as used for tracking in the VXD+SIT. <br>
 * STEP 3 : MERGING TRACK SEGMENTS FOUND IN FTD AND VXD+SIT <br>
 * In the last step, track segments reconstructed in the FTD and VXD+SIT, belonging to the
 * same track  are identified and merged into one track. All possible 
 * pairings are tested for their compatibility.
 * The number of pairings considered is Ntrk_VXD_SIT*Ntrk_FTD, where Ntrk_VXD_SIT is the number of 
 * track segments reconstructed in the first step in VXD+SIT (segments containing solely VXD and SIT hits) and  
 * Ntrk_FTD is the number of track segments reconstructed in the second step 
 * (segments containing solely FTD hits).
 * Pair of segments is accepted for further examination if the angle between track segments and 
 * than certain specified threshold.
 * Pairing satisfying this condition is subjected for 
 * addtitional test. The fit is performed on unified array of hits belonging to both segments. 
 * If the chi2 of the fit does not exceed predefined cut value two segments are unified into 
 * one track. 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collection of digitized vertex, sit and ftd tracker hits. <br>
 * If such a collections with the user specified names do not exist 
 * processor takes no action. <br>
 * <h4>Output</h4>
 * Processor produces an LCIO collection of the Tracks. Each track is characterised by 
 * five parameters : Omega (signed curvuture), Tan(lambda) where
 * lambda is the dip angle, Phi (azimuthal angle @ point of closest approach), D0 (signed impact parameter),
 * Z0 (displacement along z axis at the point of closest approach to IP). Covariance matrix for these parameters is also provided.
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
 * Output track collection has a name "SiTracks". <br>
 * In addition collection of relations of the tracks to MCParticles is stored 
 * if flag CreateMap is set to 1. <br>
 * Collection of relations has a name "SiTracksMCP" <br>
 * @param VXDHitCollectionName name of input VXD TrackerHit collection <br>
 * (default parameter value : "VTXTrackerHits") <br>
 * @param FTDHitCollectionName name of input FTD TrackerHit collection <br>
 * (default parameter value : "FTDTrackerHits") <br>
 * @param SITHitCollectionName name of input SIT TrackerHit collection <br>
 * (default parameter value : "SITTrackerHits") <br>
 * @param LayerCombinations combinations of layers used to search for hit triplets in VXD+SIT <br>
 * (default parameters : 6 4 3  6 4 2  5 4 3  5 4 2  4 3 2  4 3 1  4 2 1  3 2 1) <br> 
 * Note that in the VXD+SIT system the first and the second layers of SIT have indicies nLayerVTX and nLayerVTX+1. 
 * Combination given above means that triplets are looked first in layers 6 4 3, and then 
 * in 6 4 2;  5 4 3;  5 4 2;  4 3 2; 4 3 1; 4 2 1 and finally 3 2 1. NOTE THAT LAYER INDEXING STARTS FROM 0.
 * LAYER 0 is the innermost layer  <br>
 * @param LayerCombinationsFTD combinations of layers used to search for hit triplets in FTD <br>
 * (default parameters 5 4 3 5 4 2 5 4 1 5 3 2 5 3 1 5 2 1 4 3 2 4 3 1 4 3 0 4 2 0 4 1 0 3 2 1 3 2 0 2 1 0). 
 * NOTE THAT TRACKS IN FTD ARE SEARCHED ONLY IN ONE HEMISPHERE. TRACK IS NOT 
 * ALLOWED TO HAVE HITS BOTH IN BACKWARD AND FORWARD PARTS OF FTD SIMULTANEOUSLY. 
 * @param NDivisionsInPhi Number of divisions in Phi for tracking in VXD+SIT <br>
 * (default value is 60) <br>
 * @param NDivisionsInTheta Number of divisions in cosQ for tracking in VXD+SIT <br>
 * (default value is 60) <br>
 * @param NDivisionsInPhiFTD Number of divisions in Phi for tracking in FTD <br>
 * (default value is 3) <br>
 * @param Chi2WRphiTriplet weight on chi2 in R-Phi plane for track with 3 hits <br>
 * (default value is 1) <br>
 * @param Chi2WZTriplet weight on chi2 in S-Z plane for track with 3 hits <br>
 * (default value is 0.5) <br>
 * @param Chi2WRphiQuartet weight on chi2 in R-Phi plane to accept track with 4 hits <br>
 * (default value is 1) <br>
 * @param Chi2WZQuartet weight on chi2 in S-Z plane for track with 4 hits <br>
 * (default value is 0.5) <br>
 * @param Chi2WRphiSeptet weight on chi2 in R-Phi plane for track with 5 and more hits <br>
 * (default value is 1) <br>
 * @param Chi2WZSeptet Cut on chi2 in S-Z plane for track with 5 and more hits <br>
 * (default value is 0.5) <br>
 * @param Chi2FitCut Cut on chi2/ndf to accept track candidate
 * (default value is 50.)
 * @param ResolutionRPhiVTX Point resolution in R-Phi plane in mm for VXD <br> 
 * (default value is 0.004)
 * @param ResolutionZVTX Point resolution along Z coordinate in mm for VXD <br> 
 * (default value is 0.004)
 * @param ResolutionRPhiFTD Point resolution in R-Phi plane in mm for FTD <br> 
 * (default value is 0.01)
 * @param ResolutionZFTD Point resolution along Z coordinate in mm for FTD <br> 
 * (default value is 0.1)
 * @param ResolutionRPhiSIT Point resolution in R-Phi plane in mm for SIT <br> 
 * (default value is 0.01) 
 * @param ResolutionZSIT Point resolution along Z coordinate in mm for SIT <br> 
 * (default value is 0.01)  
 * @param AngleCutForMerging cut on the angle between two track segments.  
 * If the angle is greater than this cut, segments are not allowed to be merged. <br>
 * (default value is 0.1) <br>
 * @param MinDistCutAttach cut on the distance (in mm) from hit to the helix. This parameter is used
 * to decide whether hit can be attached to the track. If the distance is less than 
 * cut value. The track is refitted with a given hit being added to the list of hits already 
 * assigned for the track. Additional hit is assigned if chi2 of the new fit has good chi2. <br>
 * (default value is 1 ) <br>
 * @param CutOnZ0 cut on Z0 parameter of track (in mm). If abs(Z0) is greater than the cut value, track is 
 * discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 100) <br>
 * @param CutOnD0 cut on D0 parameter of track (in mm). If abs(D0) is greater than the cut value, track is 
 * discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 100) <br>
 * @param CutOnPt cut on Pt (GeV/c). If Pt is less than this cut, track is discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 0.1) <br>
 * @param MinimalHits minimal number of hits in track required <br>
 * (default value is 3) <br>
 * @param FastAttachment if this flag is set to 1, less accurate but 
 * fast procedure to merge additional hits to tracks is used <br> 
 * if set to 0, a more accurate, but slower procedure is invoked <br>
 * (default value is 1) <br>
 * @param OptPrefit Option for prefit of the track with simple helix model. If 
 * set to 0, helix fit based on FORTRAN code tfithl is used, when set to 1 ClusterShapes class
 * is used to fit track with the simple helix model <br>
 * (default value is 0) <br>
 * @param SimpleHelixFit Flag to enable fast procedure of track fitting based on simple helix fit.
 * This procedure is used when performing pattern recognition.<br>
 * (default value is 1)
 * @param UseSIT When this flag is set to 1, SIT is included in pattern recognition. When this flag is set
 * to 0, SIT is excluded from the procedure of pattern recognition <br>
 * (default value is 1) <br>
 * @param FinalRefit If set to 1, final track candidates are refitted using DELPHI fitting code, which 
 * accounts for effects of multiple scattering and energy loss <br>
 * (default value is 1) <br>
 * @param CreateMap When this flag is set to 1 collection of relations between tracks and MCParticles is 
 * created <br>
 * (default value is 1) <br>
 * @param UseExtraPoint This flag is used to steer DELPHI fitting code. If set to 1, additional 
 * artificial mesurement point at PCA is introduced with relatively large errors. This helps
 * to improve resolution on D0 and Z0 for fitted track. <br>
 * (default value 0)
 * <br>
 * @author A. Raspereza (MPI Munich)<br>
 */
class SiliconTracking : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SiliconTracking ; }
  
  
  SiliconTracking() ;
  
  /**  
   * Initialization
   */
  virtual void init() ;
  
  /** Run header processor.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Event processor.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  int _nRun ;
  int _nEvt ;

  int _nDivisionsInPhi;
  int _nDivisionsInTheta;
  int _nLayers;

  int _nLayersFTD;
  int _nLayersVTX;
  int _nLayersSIT;
  int _nPhiFTD;

  std::string _VTXHitCollection;
  std::string _FTDHitCollection;
  std::string _SITHitCollection;
  std::string _relCollection;
  
  std::vector<TrackerHitExtendedVec> _sectors;
  std::vector<TrackerHitExtendedVec> _sectorsFTD;
  
  TrackExtendedVec _tracks5Hits;
  TrackExtendedVec _tracks4Hits;
  TrackExtendedVec _tracks3Hits;

  int InitialiseVTX(LCEvent * evt);
  int InitialiseFTD(LCEvent * evt);
  void ProcessOneSector(int iSectorPhi, int iSectorTheta);
  void CleanUp();
  TrackExtended * TestTriplet(TrackerHitExtended * outerHit, 
			      TrackerHitExtended * middleHit,
			      TrackerHitExtended * innerHit,
			      HelixClass & helix);
  int BuildTrack(TrackerHitExtended * outerHit, 
		 TrackerHitExtended * middleHit,
		 TrackerHitExtended * innerHit,
		 HelixClass & helix, 
		 int innerlayer,
		 int iPhiLow, int iPhiUp,
		 int iTheta, int iThetaUp,
		 TrackExtended * trackAR);

  void Sorting( TrackExtendedVec & trackVec);
  void CreateTrack(TrackExtended * trackAR );
  void AttachRemainingVTXHitsSlow();
  void AttachRemainingFTDHitsSlow();
  void AttachRemainingVTXHitsFast();
  void AttachRemainingFTDHitsFast();
  void TrackingInFTD();
  int BuildTrackFTD(TrackExtended * trackAR, int * nLR, int iS);
  int AttachHitToTrack(TrackExtended * trackAR, TrackerHitExtended * hit);
  void FinalRefit();

  float _bField;
  float _chi2WRPhiTriplet;
  float _chi2WRPhiQuartet;
  float _chi2WRPhiSeptet;
  float _chi2WZTriplet;
  float _chi2WZQuartet;
  float _chi2WZSeptet;
  float _minDistCutAttach;
  int _minimalLayerToAttach;

  double PI,TWOPI,PIOVER2;
  double _dPhi;
  double _dTheta;
  double _dPhiFTD;

  std::vector<float> _zLayerFTD;
  std::vector<int> _Combinations;
  std::vector<int> _CombinationsFTD;

  float _resolutionRPhiVTX;
  float _resolutionZVTX;
  
  float _resolutionRPhiFTD;
  float _resolutionZFTD;

  float _resolutionRPhiSIT;
  float _resolutionZSIT;

  float _phiCutForMerging;
  float _tanlambdaCutForMerging;
  float _angleCutForMerging;

  int _print;

  float _distRPhi;
  float _distZ;
  float _chi2FitCut;

  float _chi2PrefitCut;

  TrackExtendedVec _trackImplVec;

  float _cutOnD0, _cutOnZ0, _cutOnOmega, _cutOnPt;

  int _minimalHits;
  int _attachFast;

  int _nTotalVTXHits,_nTotalFTDHits,_nTotalSITHits;
  int _optFit,_simpleHelixFit;
  int _useSIT;
  int _finalRefit;
  int _createMap;
  int _useExtraPoint;
  MarlinTrackFit _trackFit;


} ;

#endif



