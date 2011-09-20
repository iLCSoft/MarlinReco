#ifndef FORMTRUETRACKSAR_H
#define FORMTRUETRACKSAR_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include "../../BrahmsTracking/include/MarlinTrackFit.h"
#include <EVENT/TrackerHit.h>

using namespace lcio ;
using namespace marlin ;

/** === TrackCheater Processor === <br>
 *  Constructs true tracks. True track is considered to comprise 
 *  all the TrackerHits attributable to the same particle. 
 *  An user has to provide the names of TrackerHit collections 
 *  via processor parameter TrackerHitCollections. 
 *  An output collection of tracks is specified 
 *  with processor parameter TrueTrackCollection. Furthermore,
 *  relations between tracks and MCParticles are created and stored for each event. 
 *  There is an option to perform a fit of TrackerHits assigned to 
 *  track. This is done if processor parameter FitTrueTrack is set 
 *  to 1. In this case parameters extracted from the fit are used 
 *  to define track parameters, namely Omega (signed curvature), tan(lambda),
 *  where lamda is the track dip angle, Phi (azimuthal angle of track momentum at 
 *  the point of closest approach to IP), D0 (signed impact parameter) and Z0 
 *  (displacement along z axis at the point of closest approach to IP).
 *  Covariance matrix for these parameters is also provided.
 *  Only lower left corner of the covariance matrix is stored. The sequence of the covariance matrix elements 
 *  assigned to track is the following: <br>
 *  (Omega,Omega) <br>
 *  (Omega,TanLambda), (TanLambda,TanLambda) <br>
 *  (Omega,Phi), (TanLamda,Phi), (Phi,Phi) <br>
 *  (Omega,D0), (TanLambda,D0), (Phi,D0), (D0,D0) <br>
 *  (Omega,Z0), (TanLambda,Z0), (Phi,Z0), (D0,Z0), (Z0,Z0) <br>
 *  If FitTrueTrack is set to 0 then true Monte Carlo information, namely 
 *  MCParticle momentum and vertex, is used to define 
 *  track parameters. The number of hits in the different subdetectors associated
 *  with each track can be accessed via method Track::getSubdetectorHitNumbers().
 *  This method returns vector of integers : <br>
 *  number of VXD hits used in the track fit is the first element in this vector  
 *  (Track::getSubdetectorHitNumbers()[0]) <br>
 *  number of FTD hits used in the track fit is the second element in this vector  
 *  (Track::getSubdetectorHitNumbers()[1]) <br>
 *  number of SIT hits used in the track fit is the third element in this vector  
 *  (Track::getSubdetectorHitNumbers()[2]) <br>
 *  number of TPC hits used in the track fit is the forth element in this vector  
 *  (Track::getSubdetectorHitNumbers()[3]) <br>
 *  total number of VXD hits in track is the fifth element in this vector 
 *  (Track::getSubdetectorHitNumbers()[4]) <br>
 *  total number of FTD hits in track is the sixth element in this vector
 *  (Track::getSubdetectorHitNumbers()[5]) <br>
 *  total number of SIT hits in track is the seventh element in this vector
 *  (Track::getSubdetectorHitNumbers()[6]) <br>
 *  total number of TPC hits in track is the eighth element in this vector
 *  (Track::getSubdetectorHitNumbers()[7]) <br>
 *  <br>
 *  <h4>Input collections and prerequisites</h4> 
 *  Processor requires collections of digitized tracker hits 
 *  in tracker subdetectors 
 *  namely VXD, FTD,  SIT and TPC. <br>
 *  If such collections with the user specified names do not exist, 
 *  processor takes no action. <br>
 *  <h4>Output</h4>
 *  Processor produces collection of true MC tracks and collection of relations between
 *  true MC tracks and MCParticles. <br>
 *  @param TrackerHitCollections Vector of strings containing collection names of TrackerHits <br>
 *  (default names VTXTrackerHits FTDTrackerHits SITTrackerHits TPCTrackerHits) <br>
 *  @param TrueTrackCollection The name of the output collection of true MC tracks <br>
 *  (default name TrueTracks) <br>
 *  @param MCTrueTrackRelCollectionName The name of the output collection 
 *  of relations between true MC tracks and MCParticles <br>
 *  (default name TrueTracksMCP) <br>
 *  @param HitToHelixDist Cut on distance between hit and the true MC track 
 *  helix (in mm). True MC track helix is 
 *  defined by MCParticle momentum and vertex. If the distance between 
 *  hit and helix is greater than this cut, hit is 
 *  not assigned to track <br>
 *  (default value is 100) <br>
 *  @param HitToHelixInFit Cut on distance between hit and true MC 
 *  track helix (in mm) to include hit into fitting procedure. 
 *  True MC track helix is 
 *  defined by MCParticle momentum and vertex. If the distance between 
 *  hit and helix is greater than this cut, hit is 
 *  not used in the track fit<br>
 *  (default value is 50) <br>
 *  @param ECut Lower cut on the energy of MCParticle (in GeV) <br>
 *  (default value is 0.1) <br>
 *  @param FitTrueTrack When this flag is set to 1 tracks are fitted and covariance matrix for each
 *  track is calculated. If this flag is set to 0 then track parameters are calculated using 
 *  true MC information at the generator level. No covariance matricies are provided in this case <br>
 *  (default value is 1) <br>
 *  @param Chi2Cut Cut on the chi2/ndf for the track fit. Tracks failing this cut are dropped
 *  from output collection. This occurs only if the flag FitTrueTrack is set to 1 <br>
 *  (default value is 100) <br>
 *  @param CutOnD0 Cut on the reconstructed D0 parameter [mm] of the track candidate. 
 *  If D0 is greater than this cut, track candidate is rejected <br>
 *  (default value 500) <br>
 *  @param CutOnZ0 Cut on the reconstructed Z0 parameter [mm] of the track candidate. 
 *  If Z0 is greater than this cut, track candidate is rejected <br>
 *  (default value 500) <br>
 *  @param CutOnTPCHits minimal number of TPC hits required for tracks which have no 
 *  hits from the Si detectors <br>
 *  (default value 40) <br> 
 *  @param MinimalHits Minimal required number of hits in track. If number of hits
 *  assigned to true track is less than MinimalHits, track is dropped from output collection <br>
 *  (default value is 3) <br>
 *  @param UseExtraPoint This flag is used to steer DELPHI fitting code. If set 
 *  to 1, an additional 
 *  artificial mesurement point at PCA is introduced with relatively 
 *  large errors (OBSOLETE : used in early versions of the code to improve d0 and z0
 *  resolutions, for the current version of the LDC Tracking code this parameter is 
 *  recommended to be set to 0). <br> 
 *  (default value is 0) <br>
 *  @param OptFit Option for the track fit. <br>
 *  0 - simple helix model is used for track prefit, prefit is done with tfithl 
 *  routine of the DELPHI code, <br>
 *  1 - simple helix model is used for track prefit, prefit is done with the 
 *  ClusterShapes class, <br>
 *  2 - true Monte Carlo parameters of the charged particle at the generator
 *  level are used as initial approximation for the track fit, <br>
 *  3 - a sophisticated multiiterative procedure is employed to determine
 *  initial track parameters passed to the Kalman filter, <br>
 *  4 - track parameters are determined in an optimal way by combining
 *  separate fits of Si and TPC track segments. Parameters d0, z0, phi0 and
 *  tan(lambda) are determinde from the fit of the Si segment (if it contains 
 *  3 or more hits), while parameter omega is determined from the fit of 
 *  the entire set of hits contributing to a given track (recommended option) <br>
 *  (default value is 3) <br>
 *  @param StoreHitsInFit if set to 1 only hits used in fit are stored
 *  in the LCIO object Track,
 *  if StoreHitsInFit = 0 all hits belonging to a given track are stored
 *  in the object Track <br>
 *  (defalut value is 0) <br>
 *  @param aParameterForIPError parameter a, defining minimal allowed IP resolution
 *  according to the formular sigma[IP] = a + b/[P*sin^3/2{Q}]^s where P is the particle momentum
 *  and Q is the polar angle <br>
 *  (default value is 0.002 [mm]) <br>
 *  @param bParameterForIPError parameter b in the parametrisation of the minimal IP resolution
 *  sigma[IP] = a + b/[P*sin^3/2{Q}]^s <br>
 *  (default value is 0.0076 [mm]) <br>
 *  @param sParameterForIPError parameter s in the parametrisation of the minimal IP resolution
 *  sigma[IP] = a + b/[P*sin^3/2{Q}]^s <br>
 *  (default value is 0.75 [mm]) <br>
 *  @param Debug defines level of details for the printout information (meant for debuging purposes).
 *  If set to 0, printout is suppressed <br>
 *  (default value is 1)
 *  <br>
 *  @author A. Raspereza (MPI Munich)
 *  @version $Id$<br>
 */
class TrackCheater : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TrackCheater ; }
  
  
  TrackCheater() ;
  
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

  void SortTrackerHitsByRadius(TrackerHitVec & trackerHitVec);

  int _nRun ;
  int _nEvt ;

  int _fitTrueTrack;
  int _minimal_hits;
  float _cutOnChi2;

  int _cutOnTPCHits;
  float _cutOnD0,_cutOnZ0;

  std::string _trueTracksCollection;
  std::vector<std::string> _trackerHitCollections;
  std::string _colNameMCTrueTracksRel;
  std::vector<float> _deviations;
  
  float _bField;
  float _eCut;
  float _hitToHelixCut;
  float _hitToHelixInFit;

  float PI;
  int _useOnlyOneLoop;
  int _useExtraPoint;
  int _optFit;
  int _storeHitsInFit;
  
  float _aParIpReso,_bParIpReso,_sParIpReso;

  int _debug;

  MarlinTrackFit _trackFit;

 
} ;

#endif



