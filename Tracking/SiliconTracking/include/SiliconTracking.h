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

using namespace lcio ;
using namespace marlin ;


/** === Silicon Tracking Processor === <br>
 * Processor performing stand-alone pattern recognition
 * in the vertex detector (VXD) and forward tracking disks. <br>
 * The procedure consists of three steps : <br> 
 * 1) Tracking in VXD ; <br>
 * 2) Tracking in FTD ; <br>
 * 3) Merging compatible track segments reconstructed in VXD and FTD <br>
 * STEP 1 : TRACKING IN VXD <br>
 * Algorithm starts with finding of hit triplets satisfying helix hypothesis 
 * in three different layers. To accelerate procedure, the 4-pi solid angle
 * is divided in NDivisionsInTheta and NDivisionsInPhi sectors in cosQ and Phi, 
 * respectively. Triplets are looked for in 2x2 window of adjacent sectors. 
 * Once triplet is found, attempt is made to associate additional hits to 
 * track. Combinatin of hits is accepted for further analysis if the Chi2 
 * of the helix fit is less than certain predefined threshold. All accepted 
 * combinations are sorted in ascending order of their Chi2. First track candidate 
 * in the sorted array is automatically accepted. The hits belonging to this track are 
 * marked as used, and track candidates sharing these hits are discarded.
 * The procedure proceeds with increasing index of track candidate in the sorted 
 * array until all track candidate have been output or discarded. <br>
 * STEP 2 : TRACKING IN FTD <br>
 * In the next step tracking in FTD is performed. The strategy of tracking in the FTD 
 * is the same as used for tracking in the VXD.
 * STEP 3 : MERGING TRACK SEGMENTS FOUND IN FTD AND VXD <br>
 * In the last step, track segments reconstructed in the FTD and VXD, belonging to the
 * same track  are identified and merged into one track. All possible 
 * pairings are tested for their compatibility.
 * The number of pairings considered is Ntrk_VXD*Ntrk_FTD, where Ntrk_VXD is the number of 
 * track segments reconstructed in the first step in VXD (segments containing solely VXD hits) and  
 * Ntrk_FTD is the number of track segments reconstructed in the second step 
 * (segments containing solely FTD hits).
 * Pair of segments is accepted for further examination if the difference in tan(dip angle) 
 * and azimuthal angle
 * between segments is less than certain . Pairing satisfying this condition is subjected for 
 * addtitional test. The fit is performed on unified array of hits belonging to both segments. 
 * If the chi2 of the fit does not exceed predefined cut value two segments are unified into 
 * one track.     
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collection of digitized vertex tracker hits. <br>
 * If such a collection with the user specified name does not exist <br>
 * processor takes no action.
 * <h4>Output</h4>
 * Processor produces an LCIO collection of the Tracks. Each track is characterised by 
 * five parameters : Omega (signed curvuture), Phi (azimuthal angle), D0 (signed impact parameter),
 * Z0 (displacement along z axis at the point of closest approach to IP) and Tan(lambda) where
 * lambda is the dip angle. Covariance matrix for these parameters is also provided.
 * Output track collection has a name "VertexTracks". <br>
 * @param VXDHitCollectionName name of input VXD TrackerHit collection <br>
 * (default parameter value : "VTXTrackerHits") <br>
 * @param FTDHitCollectionName name of input FTD TrackerHit collection <br>
 * (default parameter value : "FTDTrackerHits") <br>
 * @param NumberOfLayers number of layers in the vertex detector <br>
 * (default parameter value : 5 ) <br>
 * @param NumberOfFTDLayers number of layers in the forward tracking detector <br>
 * (default parameter value : 7 ) <br>
 * @param LayerCombinations combinations of layers used to search for hit triplets in VXD <br>
 * (default parameters : 4 3 2 4 3 1 4 2 1 3 2 1; this means that triplets are looked for 
 * first in layers 4 3 2, then in 4 3 1; 4 2 1 and finally 3 2 1. NOTE THAT LAYER INDEXING STARTS FROM 0.
 * LAYER 0 is the innermost layer ) <br>
 * @param LayerCombinationsFTD combinations of layers used to search for hit triplets in FTD <br>
 * (default parameters 5 4 3 5 4 2 5 4 1 5 3 2 5 3 1 5 2 1 4 3 2 4 3 1 4 3 0 4 2 0 4 1 0 3 2 1 3 2 0 2 1 0). 
 * NOTE THAT TRACKS IN FTD ARE SEARCHED ONLY IN ONE HEMISPHERE. TRACK IS NOT 
 * ALLOWED TO HAVE HITS BOTH IN BACKWARD AND FORWARD PARTS OF FTD SIMULTANEOUSLY. 
 * @param NDivisionsInPhi Number of divisions in Phi for tracking in VXD <br>
 * (default value is 60) <br>
 * @param NDivisionsInTheta Number of divisions in cosQ for tracking in VXD <br>
 * (default value is 60) <br>
 * @param NDivisionsInPhiFTD Number of divisions in Phi for tracking in FTD <br>
 * (default value is 40) <br>
 * @param Chi2WRphiTriplet Cut on chi2 in R-Phi plane to accept track with 3 hits <br>
 * (default value is 10) <br>
 * @param Chi2WZTriplet Cut on chi2 in S-Z plane to accept track with 3 hits <br>
 * (default value is 20) <br>
 * @param Chi2WRphiQuartet Cut on chi2 in R-Phi plane to accept track with 4 hits <br>
 * (default value is 10) <br>
 * @param Chi2WZQuartet Cut on chi2 in S-Z plane to accept track with 4 hits <br>
 * (default value is 20) <br>
 * @param Chi2WRphiSeptet Cut on chi2 in R-Phi plane to accept track with 5 and more hits <br>
 * (default value is 10) <br>
 * @param Chi2WZSeptet Cut on chi2 in S-Z plane to accept track with 5 and more hits <br>
 * (default value is 20) <br>
 * @param ResolutionRPhi Point resolution in R-Phi plane in mm for VXD <br> 
 * (default value is 0.005)
 * @param ResolutionZ Point resolution along Z coordinate in mm for VXD <br> 
 * (default value is 0.005)
 * @param ResolutionRPhiFTD Point resolution in R-Phi plane in mm for FTD <br> 
 * (default value is 0.1)
 * @param TanLambdaCutForMerging cut on the difference in tan(lambda) for two track segments.  
 * If the difference is greater than this cut, segments are not allowed to be merged. <br>
 * (default value is 0.05) <br>
 * @param PhiCutForMerging cut on the difference in azimuthal angle for two track segments.  
 * If the difference is greater than this cut, segments are not allowed to be merged. <br>
 * (default value is 0.05) <br>
 * @param MinDistCutAttach cut on the distance (in mm) from hit to the helix. This parameter is used
 * to decide whether hit can be attached to the track. If the distance is less than 
 * cut value. The track is refitted with a given hit being added to the list of hits already 
 * assigned for the track. Additional hit is assigned if chi2 of the new fit has good chi2. <br>
 * (default value is 1 ) <br>
 * @param CutOnZ0 cut on Z0 parameter of track (in mm). If abs(Z0) is greater than the cut value, track is 
 * discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 50) <br>
 * @param CutOnD0 cut on D0 parameter of track (in mm). If abs(D0) is greater than the cut value, track is 
 * discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 50) <br>
 * @param CutOnPt cut on Pt (GeV/c). If Pt is less than this cut, track is discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 0.1) <br>
 * @param MinimalHits minimal number of hits in track required <br>
 * (default value is 4) <br>
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
  std::string _SiTracks;
  std::string _SiTracksMCP;
  
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
  int _simpleHelixFit;

  int _useSIT;

  int _finalRefit;
  
  int _createMap;

} ;

#endif



