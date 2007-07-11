#ifndef TrackBasedPFlow_h
#define TrackBasedPFlow_h 1

#include <iostream>
#include <iomanip>

#include <vector>
#include <utility>
#include <algorithm>


#include <marlin/Processor.h>
#include <marlin/Exceptions.h>
#include <marlin/Global.h>

#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <IMPL/ClusterImpl.h>
#include <EVENT/Cluster.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <UTIL/LCRelationNavigator.h>
#include <LCRTRelations.h>

#include <MarlinUtil.h>

#include <gear/GEAR.h>
#include <gear/TPCParameters.h>

#include "constants.h"
#include "LCGeometryTypes.h"
#include "Trajectory.h"
#include "SimpleHelix.h"
#include "ClusterShapes.h"
#include "CalorimeterHitWithAttributes.h"
#include "ClusterImplWithAttributes.h"
#include "NNClusters.h"
#include "TrackwiseClusters.h"






// MarlinCED is only used for debugging
#include <MarlinCED.h>
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/ICloud2D.h>
#endif








using namespace lcio;
using namespace marlin;
using namespace constants;



// integer runtime extension which flags CalorimeterHits in EMShowers
struct isPartOfEMShowerCandidate : LCIntExtension<isPartOfEMShowerCandidate> {};




/**
 *    This processor is an initial version of a Track-Based Particle Flow approach. More docu will come soon.
 *    <br>
 *
 *    @author O. Wendt (DESY)
 *    @version
 *
 */
class TrackBasedPFlow : public Processor {
  
 public:

  virtual Processor* newProcessor() { return new TrackBasedPFlow; }

  TrackBasedPFlow();
  
  virtual void init();
  virtual void processRunHeader( LCRunHeader* run );
  virtual void processEvent( LCEvent* evt ); 
  virtual void check( LCEvent* evt ); 
  virtual void end();

  
 private:

  
 protected:

  // 'global constants
  float _cutOnAbsMomentumForFitOnOuterMostTrackerHits;
  

  float _bField;

  std::string _colNameTracks;
  std::string _colNameECAL;
  std::string _colNameHCAL;
  std::string _colNameEMShowerCandidates;
  std::string _colNameRelationTrackToMCP;
  std::string _colNameRelationCaloHitToSimCaloHit;
  std::string _reconstructedClusterCollectionName;
  std::string _reconstructedParticleCollectionName;

  std::vector<float> _calibrCoeffECAL;
  std::vector<float> _calibrCoeffHCAL;

  double _cutOnPt;
  double _absD0Cut;
  double _absZ0Cut;
  int _minNTPCHits;
  int _minNNonTPCHits;
  int _nOfTrackerHitsUsedForExtrapolation;
  int _nOfTrackerHitsOutsideCylindricalCut;
  double _rMinCutHelixExtrapolation;
  double _zMinCutHelixExtrapolation;
  double _rMinCylindricalCut;
  double _zMinCylindricalCut;
  double _openingAngleConeTube;
  double _maximalConeTubeLength;
  double _maximalPathLengthForMIPLikeStub;
  std::vector<float> _maximalRadiusOfInnerTubeForMIPLikeStub;
  std::vector<float> _minimalRadiusOfOuterTubeForMIPLikeStub;
  std::vector<float> _maximalDistanceToHelixToAssignCluster;
  std::vector<float> _maximalDistanceOfElectronShowerPosition;
  double _fractionEM;
  double _outputConditionLimit;
  std::vector<float> _mipCoeffEcal;
  std::vector<float> _mipCoeffHcal;
  //  int _doComparisonWithMC;
  int _drawOnCED;
  int _debugLevel;


  // parameters needed for Trackwise Clustering
  std::vector<float> _distanceTrackBack;
  std::vector<float> _stepTrackBack;
  std::vector<float> _resolutionParameter;
  std::vector<float> _resolutionParameterForNeutrals;
  std::vector<float> _distanceMergeForward;
  float _distanceToTrackSeed;
  float _distanceToDefineDirection;
  float _resolutionToMerge;
  int _nhit_merge_forward;
  int _use_tracks;
  int _nhit_minimal;
  int _nhit_neutral_minimal;
  int _typeOfGenericDistance;
  int _typeOfGenericDistanceForNeutrals;

  int _doMerging;
  int _doMergingForward;
  int _displayClusters;
  int _NDefineSP;
  int _nScanToMergeForward;

  float _zofendcap;
  float _rofbarrel;
  float _phiofbarrel;
  int _nsymmetry;
  float _thetaofendcap;
  float _weightForReso;
  float _weightForDist;
  float _weightForResoForNeutrals;
  float _weightForDistForNeutrals;


  // debug
  MCParticleHelper* _mcParticleHelper;

  int _nOfRealMIPStubs;
  int _nOfFoundMIPStubs;
  int _nOfChargedObjectsExtrapolatedIntoCalorimeter;

  std::vector<Track*> _tracksExtrapolatedIntoCalorimeter;

  std::vector<Track*> _tracksNotExtrapolatedIntoCalorimeter;
  std::vector<Track*> _tracksNotFulFillingPtCut;
  std::vector<Track*> _tracksWithTooFewTrackerHits;
  std::vector<Track*> _tracksNotFulfillingCombinedCylinderShellCuts;

  std::vector<Track*> _tracksNotExtrapolatedWithEnoughSiliconAndTPCHits;
  std::vector<Track*> _tracksNotExtrapolatedComingFromIP;
  std::vector<Track*> _tracksDiscarded;

  std::vector<Track*> _tracksWhichWouldReachCalorimeter;
  std::vector<Track*> _tracksNotReachingTheCalorimeter;

  std::vector<Cluster*> _emShowerCandidatesRecognisedAsCharged;



  #ifdef MARLIN_USE_AIDA
  AIDA::IHistogram1D* _hPtRealMIPStubs;
  AIDA::IHistogram1D* _hPtFoundMIPStubs;
  AIDA::ICloud1D* _cPtRealMIPStubs;
  AIDA::ICloud1D* _cPtFoundMIPStubs;
  AIDA::ICloud2D* _cPtNRealMIPStubs;
  AIDA::ICloud2D* _cPtNFoundMIPStubs;
  AIDA::IHistogram1D* _hCosThRealMIPStubs;
  AIDA::IHistogram1D* _hCosThFoundMIPStubs;
  AIDA::ICloud1D* _cCosThRealMIPStubs;
  AIDA::ICloud1D* _cCosThFoundMIPStubs;
  AIDA::ICloud2D* _cCosThNRealMIPStubs;
  AIDA::ICloud2D* _cCosThNFoundMIPStubs;
  AIDA::ICloud1D* _cNExtrapolatedObjects;

  AIDA::ICloud1D* _cSECalo;
  AIDA::ICloud1D* _cSEReco;
  AIDA::ICloud1D* _cInvMass;
  AIDA::ICloud1D* _cSEMC;
  AIDA::ICloud1D* _cDSERecoMC;
  AIDA::ICloud1D* _cDSERecoCalo;
  AIDA::ICloud1D* _cDSECaloMC;
  AIDA::ICloud1D* _cAlphaRecoMC;
  AIDA::ICloud1D* _cAlphaRecoCalo;
  AIDA::ICloud1D* _cAlphaCaloMC;

  AIDA::ICloud1D* _cNTracksWithEnergyAssignedByMC;
  AIDA::ICloud1D* _cNTracksNotPassingCylinderCuts;
  AIDA::ICloud1D* _cPTracksNotPassingCylinderCuts;
  AIDA::ICloud1D* _cSumPTracksNotPassingCylinderCuts;
  AIDA::ICloud2D* _cSumPTracksNotPassingCylinderCutsvsDSERecoMC;

  AIDA::ICloud1D* _cNNeutralHitsAssignedToCharged;
  AIDA::ICloud1D* _cENeutralHitsAssignedToCharged;
  AIDA::ICloud1D* _cSumNNeutralHitsAssignedToCharged;
  AIDA::ICloud1D* _cSumENeutralHitsAssignedToCharged;
  AIDA::ICloud2D* _cSumENeutralHitsAssignedToChargedvsDSERecoMC;

  AIDA::ICloud1D* _cNChargedHitsAssignedToNeutral;
  AIDA::ICloud1D* _cEChargedHitsAssignedToNeutral;
  AIDA::ICloud1D* _cSumNChargedHitsAssignedToNeutrals;
  AIDA::ICloud1D* _cSumEChargedHitsAssignedToNeutrals;
  AIDA::ICloud2D* _cSumEChargedHitsAssignedToNeutralsvsDSERecoMC;

  AIDA::ICloud2D* _cDEdcEwavsDSERecoMC;
  AIDA::ICloud2D* _cDEdcEwaMinusDSERecoMCvsDSERecoMC;

  AIDA::IHistogram1D* _hChi2;
  AIDA::IHistogram2D* _hChi2DEP;
  #endif


  int _nRun;
  int _nEvt;
 

  const TrackerHitVec getOuterTrackerHits(const Track* track, unsigned int n=12);
  const std::vector<CalorimeterHitWithAttributes*> getRelatedCalorimeterHits(const LCEvent* evt, const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix, 
									     int trackNumber/*only for debugging*/);
  
  void getRelatedClusters(const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, const TrackerHitVec outermostTrackerHits,
			  Trajectory* fittedHelix, ClusterImplWithAttributes* mipStub, std::vector<ClusterImplWithAttributes*>& clusters);

  void getRelatedClusterPerfectly(const LCEvent* evt, Track* track, const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix,
				  std::vector<ClusterImplWithAttributes*>& clusters);

  void getMIPStub(ClusterImplWithAttributes* clusterWithAttributes, const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, 
		  const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix);

// je: for gcc4 compatibility
// const std::vector<CalorimeterHitWithAttributes*> TrackBasedPFlow::removeMIPStub(std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes,
//										  ClusterImplWithAttributes* mipStub);
// 
  const std::vector<CalorimeterHitWithAttributes*> removeMIPStub(std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes,
										  ClusterImplWithAttributes* mipStub);
  
  std::vector<ClusterImplWithAttributes*> doTrackwiseClustering(const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix,
								const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, 
								const double* startPoint, const float pathLengthOnHelixOfStartPoint, const float distanceToHelixOfStartPoint,
								const double* startDirection);

  std::vector<ClusterImplWithAttributes*> doConeClustering(const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix,
							   const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, 
							   const double* startPoint, const double* startDirection);

  std::vector<ClusterImplWithAttributes*> doNNClustering(const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix,
							 const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, 
							 const double* startPoint, const double* startDirection);

  std::vector<ClusterImplWithAttributes*> doTrackBasedClustering(const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix,
								 const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, 
								 const double* startPoint, const double* startDirection);


  std::vector<ClusterImplWithAttributes*> doTrackwiseClusteringForNeutrals(const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes,
									   const double* startPoint, const double* startDirection);

  ReconstructedParticleImpl* assignClustersToTrack(Track* track, ClusterImplWithAttributes* mipStub, std::vector<ClusterImplWithAttributes*> clusters,
						   const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix);

  ReconstructedParticleImpl* assignClustersToTrackPerfectly(Track* track, std::vector<ClusterImplWithAttributes*> clusters);
  
  ReconstructedParticleImpl* assignNeutralClusterToReconstructedParticle(ClusterImplWithAttributes* cluster);

  void assignClusterProperties(ClusterImpl* cluster); // basic properties of the cluster such as energy, position, angles etc.
  void assignClusterAttributes(const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix, std::vector<ClusterImplWithAttributes*> resultingClusterImplsWithAttributes,
			       const double* startPoint); // additional cluster attributes, such as start point, start direction etc.

  void doPID(ReconstructedParticleImpl* recoParticle ,bool isMuon);

  bool isRealMIPStub(const LCEvent* evt, Track* track, const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix);

  void getRelatedCalorimeterHitsPerfectly(const LCEvent* evt, Track* track, ClusterImpl* clusterRealEnergy, ClusterImpl* clusterPerfectEnergy);
  void drawRelatedCalorimeterHits(const LCEvent* evt,Track* track);
  void drawEMShowerCandidates(const LCEvent* evt);

  std::vector<ClusterImplWithAttributes*> assignClusters(ClusterImplWithAttributes* cluster, std::vector<ClusterImplWithAttributes*> clusters, LCVector3D referencePosition,
							 const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix);

  std::vector<ClusterImplWithAttributes*> assignAdditionalClusters(ClusterImplWithAttributes* mipStub, Trajectory* fittedHelix, 
								   std::vector<ClusterImplWithAttributes*> clustersAlreadyAssigned,
								   std::vector<ClusterImplWithAttributes*> clustersToCheck, Track* track);

  
  std::vector<CalorimeterHit*> getNeutralHitsAssignedToChargedParticle(LCEvent* evt, ReconstructedParticle* recoParticle, int& n, double& energy, double hitEnergyFraction=0.5);
  std::vector<CalorimeterHit*> getChargedHitsAssignedToNeutralParticle(LCEvent* evt, ReconstructedParticle* recoParticle, int& n, double& energy, double hitEnergyFraction=0.5);
  double getNotAssignedCalorimeterEnergy(LCEvent* evt, LCCollectionVec* reconstructedParticles);

  int getTypeOfPositionOfCluster(ClusterImpl* cluster);
  int getTypeOfPositionOfCluster(Cluster* cluster);

  double getDistanceToHelix(double* point, Trajectory* helix);
  double getDistanceToHelix(std::vector<double> point, Trajectory* helix);
  double getDistanceToHelix(LCVector3D point, Trajectory* helix);

  double* getProjectedPointOnHelix(double* point, Trajectory* helix);
  std::vector<double> getProjectedPointOnHelix(std::vector<double> point, Trajectory* helix);
  LCVector3D getProjectedPointOnHelix(LCVector3D point, Trajectory* helix);

  double getPathLengthOnHelix(double* point, Trajectory* helix);
  double getPathLengthOnHelix(std::vector<double> point, Trajectory* helix);
  double getPathLengthOnHelix(LCVector3D point, Trajectory* helix);
  
  double getPathLengthOnHelix(double* point1, double* point2, Trajectory* helix);
  double getPathLengthOnHelix(std::vector<double> point1, std::vector<double> point2, Trajectory* helix);
  double getPathLengthOnHelix(LCVector3D point1, LCVector3D point2, Trajectory* helix);

} ;

#endif
