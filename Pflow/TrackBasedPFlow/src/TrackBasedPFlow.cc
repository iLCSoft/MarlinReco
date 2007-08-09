#include "TrackBasedPFlow.h"

using namespace lcio;
using namespace marlin;
using namespace constants;


TrackBasedPFlow aTrackBasedPFlow ;



TrackBasedPFlow::TrackBasedPFlow() : Processor("TrackBasedPFlow")
{
  // modify processor description
  _description = "simple track-based particle flow processor" ;

  // FIXME: remove hard coded number, if absP > 10 GeV do not perform a fir on outermost tracker hits but take track parameter as parameters for trajectory
  _cutOnAbsMomentumForFitOnOuterMostTrackerHits = 10.0;
  

  registerProcessorParameter( "colNameTracks" ,
			      "name of the Track collection" ,
			      _colNameTracks ,
			      std::string("Tracks") ) ;

  registerProcessorParameter( "colNameECAL" , 
			      "ECAL Collection Name"  ,
			      _colNameECAL,
			      std::string("ECAL") );

  registerProcessorParameter( "colNameHCAL" , 
			      "HCAL Collection Name"  ,
			      _colNameHCAL,
			      std::string("HCAL") );

  registerProcessorParameter( "colNameEMShowerCandidates" , 
			      "collection name of candidates of EM showers in the ECAL"  ,
			      _colNameEMShowerCandidates,
			      std::string("EMShowerCandidates") );

  registerProcessorParameter( "colNameRelationTrackToMCP" , 
			      "name of the LC Relation collection between Tracks and MC particles"  ,
			      _colNameRelationTrackToMCP,
			      std::string("TrackToMCP") );

  registerProcessorParameter( "colNameRelationCaloHitToSimCaloHit" , 
			      "name of the LC Relation collection between Calorimeterhits and SimCalorimeterhits"  ,
			      _colNameRelationCaloHitToSimCaloHit,
			      std::string("RelationCaloHit") );

  registerProcessorParameter("reconstructedClusterCollectionName",
			     "name of the collection of clusters assigned to the reconstructed particles",
			     _reconstructedClusterCollectionName,
			     std::string("ClustersFromTrackBasedPFlow"));

  registerProcessorParameter("reconstructedParticleCollectionName",
			     "name of the collection of reconstructed particles",
			     _reconstructedParticleCollectionName,
			     std::string("RecoParticlesFromTrackBasedPFlow"));

  std::vector<float> calibrECAL;
  calibrECAL.push_back(33.0235);
  calibrECAL.push_back(93.5682);
  registerProcessorParameter("calibrCoeffECAL" , 
			     "Calibration coefficients for ECAL" ,
			     _calibrCoeffECAL,
			     calibrECAL);
  
  std::vector<float> calibrHCAL;
  calibrHCAL.push_back(21.19626);
  registerProcessorParameter("calibrCoeffHCAL" , 
			     "Calibration coefficients for HCAL" ,
			     _calibrCoeffHCAL,
			     calibrHCAL);
  
  registerProcessorParameter( "cutOnPt",
			      "cut on pt, use only tracks with pt larger than cut (in GeV)",
			      _cutOnPt,
			      (double)0.1 ) ; // 100 MeV

  registerProcessorParameter( "absD0Cut",
			      "cut on d0 of the tracks, use only tracks with d0 smaller than cut (in mm)",
			      _absD0Cut,
			      (double)5.0 ) ; // 5 mm

  registerProcessorParameter( "absZ0Cut",
			      "cut on z0 of the tracks, use only tracks with z0 smaller than cut (in mm))",
			      _absZ0Cut,
			      (double)5.0 ) ; // 5 mm

  registerProcessorParameter( "minNTPCHits",
			      "cut on minimal number of TPC hits, use only tracks with more TPC tracker hits than cut for the tracks which were not extrapolated",
			      _minNTPCHits,
			      (int)5 ) ;

  registerProcessorParameter( "minNNonTPCHits",
			      "cut on minimal number of non TPC hits, use only tracks with more non TPC tracker hits than cut for the tracks which were not extrapolated",
			      _minNNonTPCHits,
			      (int)1 ) ;

  registerProcessorParameter( "nOfTrackerHitsUsedForExtrapolation",
			      "number of outermost Tracker hits used for the extrapolation of the track into the calorimeter",
			      _nOfTrackerHitsUsedForExtrapolation,
			      (int)90);

  registerProcessorParameter( "rMinCutHelixExtrapolation",
			      "cut on the minimal radius of the tracker hits to be taken into account in the helix extrapolation (in cylindrical coordinates)",
			      _rMinCutHelixExtrapolation,
			      (double)600.0) ;
  
  registerProcessorParameter( "zMinCutHelixExtrapolation",
			      "cut on the minimal z coordiante of the tracker hits to be taken into account in the helix extrapolation (in cylindrical coordinates)",
			      _zMinCutHelixExtrapolation,
			      (double)900.0) ;

  registerProcessorParameter( "nOfTrackerHitsOutsideCylindricalCut",
			      "number of Tracker hits (as a subset of nOfTrackerHitsUsedForExtrapolation) which are located outside cylindrical cut",
			      _nOfTrackerHitsOutsideCylindricalCut,
			      (int)25);

  registerProcessorParameter( "rMinCylindricalCut",
			      "cut on the minimal radius of the 'nOfTrackerHitsOutsideCylindricalCut' tracker hits",
			      _rMinCylindricalCut,
			      (double)800.0) ;
  
  registerProcessorParameter( "zMinCylindricalCut",
			      "cut on the minimal z coordiante of the 'nOfTrackerHitsOutsideCylindricalCut' tracker hits",
			      _zMinCylindricalCut,
			      (double)1300.0) ;
  
  registerProcessorParameter( "openingAngleConeTube",
			      "opening angle (in degree) of the cone-like tube around the track extrapolation (only calorimeter hits insider this tube are taken into account)",
			      _openingAngleConeTube,
			      (double)120.0) ;

  registerProcessorParameter( "maximalConeTubeLength",
			      "maximal length of the cone-like tube around the track extrapolation where calorimeter hits around are taken into account",
			      _maximalConeTubeLength,
			      (double)4000.0) ; // could also be implemented as a function of B field and particle momentum and type

  registerProcessorParameter( "maximalPathLengthForMIPLikeStub",
			      "maximal path length on the extrapolated helix for calorimeter hits to be taken into account as contributors to the MIP like stub",
			      _maximalPathLengthForMIPLikeStub,
			      (double)4000.0) ; // could also be implemented as a function of B field and particle momentum and type

  std::vector<float> maximalRadiusOfInnerTubeForMIPLikeStubInit;
  maximalRadiusOfInnerTubeForMIPLikeStubInit.push_back(15.0);
  maximalRadiusOfInnerTubeForMIPLikeStubInit.push_back(17.0);
  maximalRadiusOfInnerTubeForMIPLikeStubInit.push_back(40.0);

  registerProcessorParameter( "maximalRadiusOfInnerTubeForMIPLikeStub",
			      "maximal radius of an inner cylindrical tube around the extralpolated helix for Calorimeter hits to be assigned as MIP like stub (in 3 zones)",
			      _maximalRadiusOfInnerTubeForMIPLikeStub,
			      maximalRadiusOfInnerTubeForMIPLikeStubInit);// could also be implemented as a function of B field and particle momentum and type??

  std::vector<float> minimalRadiusOfOuterTubeForMIPLikeStubInit;
  minimalRadiusOfOuterTubeForMIPLikeStubInit.push_back(40.0);
  minimalRadiusOfOuterTubeForMIPLikeStubInit.push_back(45.0);
  minimalRadiusOfOuterTubeForMIPLikeStubInit.push_back(110.0);

  registerProcessorParameter( "minimalRadiusOfOuterTubeForMIPLikeStub",
			      "minimal radius of an outer cylindrical tube around the extralpolated helix for Calorimeter hits to be excluded as MIP like stub (in 3 zones)",
			      _minimalRadiusOfOuterTubeForMIPLikeStub,
			      minimalRadiusOfOuterTubeForMIPLikeStubInit);// could also be implemented as a function of B field and particle momentum and type??

  std::vector<float> maximalDistanceToHelixToAssignClusterInit;
  maximalDistanceToHelixToAssignClusterInit.push_back(200.0);
  maximalDistanceToHelixToAssignClusterInit.push_back(50.0);
  maximalDistanceToHelixToAssignClusterInit.push_back(60.0);
  maximalDistanceToHelixToAssignClusterInit.push_back(80.0);

  registerProcessorParameter( "maximalDistanceToHelixToAssignCluster",
			      "maximal distance to the extrapolated helix of the Calorimeter hit with the smallest 3-dim distance to MIP stub (in 4 zones)",
			      _maximalDistanceToHelixToAssignCluster,
			       maximalDistanceToHelixToAssignClusterInit);// could also be implemented as a function of B field and particle momentum and type, 
                                                                          // and especially the detector geometry and material


  std::vector<float> maximalDistanceOfElectronShowerPositionInit;
  maximalDistanceOfElectronShowerPositionInit.push_back(20.0);
  maximalDistanceOfElectronShowerPositionInit.push_back(25.0);


  registerProcessorParameter( "maximalDistanceOfElectronShowerPosition" ,
			      "maximal distance of electron shower position to extrapolated track",
			      _maximalDistanceOfElectronShowerPosition,
 			       maximalDistanceOfElectronShowerPositionInit);


  registerProcessorParameter( "fractionEM" ,
			      "fraction of EM Energy",
			      _fractionEM,
			      (double)0.95); // used for a simple PID
  
  
  registerProcessorParameter( "outputConditionLimit",
			      "sets processor output condition 'to true' if fabs(a) >= 'outputConditionLimit' (where a is a/sqrt(E)), otherwise 'false'",
			      _outputConditionLimit,
			      double(1.0) );
  

  std::vector<float> mipCoeffEcal;
  mipCoeffEcal.push_back(0.007);
  mipCoeffEcal.push_back(0.022);

  registerProcessorParameter("MIPCoeffEcal", 
			     "Coefficients for the MIP calibration in the ECAL in GeV/MIP",
			     _mipCoeffEcal,
			     mipCoeffEcal);

  

  std::vector<float> mipCoeffHcal;
  mipCoeffHcal.push_back(0.03);

  registerProcessorParameter("MIPCoeffHcal",
			     "Coefficients for the MIP calibration in the HCAL in GeV/MIP",
			     _mipCoeffHcal,
			     mipCoeffHcal);

  /*
  // not used at the moment
  registerProcessorParameter( "DoComparisonWithMC",
			      "toggles wether a comparison with MC tree is done (DoComparisonWithMC = 1) or not (DoComparisonWithMC = 0).",
			      _doComparisonWithMC,
			      int(1) );
  */

  registerProcessorParameter( "DrawOnCED",
			      "draw objects on CED",
			      _drawOnCED,
			      int(0) );


  registerProcessorParameter( "DebugLevel",
			      "limits the amount of information written to std out (0 - none, 9 - maximal information)",
			      _debugLevel,
			      int(0) );






  // parameters needed for Trackwise Clustering
  // Not used
  registerProcessorParameter( "DistanceForDirection", 
			      "Distance to Define Direction", 
			      _distanceToDefineDirection,
			      (float)1.0);

  // Not used
  registerProcessorParameter( "DistanceToTrackSeed", 
			      "Distance to Track Seed", 
			      _distanceToTrackSeed,
			      (float)25.0);    
  
  // RCutMax
  std::vector<float>  distanceTrackBack;
  distanceTrackBack.push_back(50.0);
  distanceTrackBack.push_back(100.0);
  
  registerProcessorParameter( "DistanceTrackBack" , 
			      "Distance to Track Back "  ,
			      _distanceTrackBack,
			      distanceTrackBack); 
  
  // RCut
  std::vector<float>  stepTrackBack;
  stepTrackBack.push_back(40.0);
  stepTrackBack.push_back(160.0);
  
  registerProcessorParameter( "StepTrackBack" , 
			      "Step to Track Back "  ,
			      _stepTrackBack,
			      stepTrackBack); 
  
  // SCut (merging parameter cut)
  std::vector<float>  resolutionParameter;
  resolutionParameter.push_back(20.0);
  resolutionParameter.push_back(60.0);
  
  registerProcessorParameter( "ResolutionParameter" , 
			      "Resolution Parameter "  ,
			      _resolutionParameter,
			      resolutionParameter); 

  // SCut (merging parameter cut for neutrals)
  std::vector<float>  resolutionParameterForNeutrals;
  resolutionParameterForNeutrals.push_back(30.0);
  resolutionParameterForNeutrals.push_back(100.0);
  
  registerProcessorParameter( "ResolutionParameterForNeutrals" , 
			      "Resolution Parameter for neutral particles"  ,
			      _resolutionParameterForNeutrals,
			      resolutionParameterForNeutrals); 

  std::vector<float>  distanceMergeForward;
  distanceMergeForward.push_back(50.0);
  distanceMergeForward.push_back(100.0);
  
  registerProcessorParameter( "DistanceMergeForward" , 
			      "Distance To Merge Forward" ,
			      _distanceMergeForward,
			      distanceMergeForward); 
  
  
  registerProcessorParameter( "NToDefineSP",
			      "N hits to define SP " , 
			      _NDefineSP, 
			      3);
  
  registerProcessorParameter( "NScanToMergeForward",
			      "N hits scan to merge forward " , 
			      _nScanToMergeForward, 
			      40);
  
  registerProcessorParameter( "TypeOfGenericDistance" , 
			      "Type of Generic Distance "  ,
			      _typeOfGenericDistance,
			      2);

  registerProcessorParameter( "TypeOfGenericDistanceForNeutrals" , 
			      "Type of Generic Distance for neutral particles"  ,
			      _typeOfGenericDistanceForNeutrals,
			      1);

  registerProcessorParameter( "MinimalHitsInCluster" ,
			      "Minimal allowed hits in cluster" , 
			      _nhit_minimal, 
			      3);

  registerProcessorParameter( "MinimalNeutralHitsInCluster" ,
			      "Minimal neutral hits allowed in cluster" , 
			      _nhit_neutral_minimal, 
			      3);
  
  registerProcessorParameter( "MaximalHitsToMerge" ,
			      "Maximal Hits To Merge" , 
			      _nhit_merge_forward, 
			      50);
  
  registerProcessorParameter( "UseTracking" ,
			      "Use tracks to seed clusters" , 
			      _use_tracks, 
			      0);
  
  registerProcessorParameter( "DoMergingLowMultiplicity" , 
			      "merging low multiplicity clusters?", 
			      _doMerging,
			       1);
  
  registerProcessorParameter( "DoMergingForward" , 
			      "merging clusters forward-wise?", 
			      _doMergingForward,
			      0);

  registerProcessorParameter( "DisplayClusterInfo",
			      "Display Info on Clusters",
			      _displayClusters,
			      0);
  
  registerProcessorParameter( "ResolutionToMerge",
			      "Resolution To Merge Halo Hits",
			      _resolutionToMerge,
			      (float)400.);

  
  registerProcessorParameter( "WeightForResolution",
			      "Weight For Resolution",
			      _weightForReso,
			      (float)4.0);
  
  registerProcessorParameter( "WeightForDistance",
			      "Weight For Distance",
			      _weightForDist,
			      (float)5.0);

  
  registerProcessorParameter( "WeightForResolutionForNeutrals",
			      "Weight For Resolution for Neutrals",
			      _weightForResoForNeutrals,
			      (float)1.0);
  
  registerProcessorParameter( "WeightForDistanceForNeutrals",
			      "Weight For Distance for Neutrals",
			      _weightForDistForNeutrals,
			      (float)1.0);


}



void TrackBasedPFlow::init()
{
  
  // usually a good idea to 
  printParameters();
 
  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;

  _bField = gearTPC.getDoubleVal("BField");


  // only the absolut value is used
  _zMinCutHelixExtrapolation = fabs(_zMinCutHelixExtrapolation);
  _zMinCylindricalCut = fabs(_zMinCylindricalCut);

  if (_nOfTrackerHitsOutsideCylindricalCut > _nOfTrackerHitsUsedForExtrapolation) _nOfTrackerHitsOutsideCylindricalCut = _nOfTrackerHitsUsedForExtrapolation;
  if (_rMinCylindricalCut < _rMinCutHelixExtrapolation) _rMinCylindricalCut = _rMinCutHelixExtrapolation;
  if (_zMinCylindricalCut < _zMinCutHelixExtrapolation) _zMinCylindricalCut = _zMinCutHelixExtrapolation;
  

  // convert into radian
  _openingAngleConeTube = (constants::twopi/360.0) * _openingAngleConeTube;



  // debug
  if (_drawOnCED)  MarlinCED::init(this);


  _mcParticleHelper = new MCParticleHelper();

  _nOfRealMIPStubs = 0;
  _nOfFoundMIPStubs = 0;
  _nOfChargedObjectsExtrapolatedIntoCalorimeter = 0;

  #ifdef MARLIN_USE_AIDA
  _hPtRealMIPStubs = AIDAProcessor::histogramFactory(this)->createHistogram1D( "PtRealMIPStubs","PtRealMIPStubs",10,0.0,20.0 );
  _hPtFoundMIPStubs = AIDAProcessor::histogramFactory(this)->createHistogram1D( "PtFoundMIPStubs","PtFoundMIPStubs",10,0.0,20.0 );
  _cPtRealMIPStubs = AIDAProcessor::histogramFactory(this)->createCloud1D( "cPtRealMIPStubs","cPtRealMIPStubs", -1 );
  _cPtFoundMIPStubs = AIDAProcessor::histogramFactory(this)->createCloud1D( "cPtFoundMIPStubs","cPtFoundMIPStubs", -1 );
  _cPtNRealMIPStubs = AIDAProcessor::histogramFactory(this)->createCloud2D( "PtNRealMIPStubs","PtNRealMIPStubs", -1 );
  _cPtNFoundMIPStubs = AIDAProcessor::histogramFactory(this)->createCloud2D( "PtNFoundMIPStubs","PtNFoundMIPStubs", -1 );
  _hCosThRealMIPStubs = AIDAProcessor::histogramFactory(this)->createHistogram1D( "CosThRealMIPStubs","CosThRealMIPStubs",10,-1.0,1.0 );
  _hCosThFoundMIPStubs = AIDAProcessor::histogramFactory(this)->createHistogram1D( "CosThFoundMIPStubs","CosThFoundMIPStubs",10,-1.0,1.0 );
  _cCosThRealMIPStubs = AIDAProcessor::histogramFactory(this)->createCloud1D( "cCosThRealMIPStubs","cCosThRealMIPStubs", -1 );
  _cCosThFoundMIPStubs = AIDAProcessor::histogramFactory(this)->createCloud1D( "cCosThFoundMIPStubs","cCosThFoundMIPStubs", -1 );
  _cCosThNRealMIPStubs = AIDAProcessor::histogramFactory(this)->createCloud2D( "CosThNRealMIPStubs","CosThNRealMIPStubs", -1 );
  _cCosThNFoundMIPStubs = AIDAProcessor::histogramFactory(this)->createCloud2D( "CosThNFoundMIPStubs","CosThNFoundMIPStubs", -1 );
  _cNExtrapolatedObjects = AIDAProcessor::histogramFactory(this)->createCloud1D( "NExtrapolatedObjects","NExtrapolatedObjects", -1 );

  _cSECalo = AIDAProcessor::histogramFactory(this)->createCloud1D( "SECalo", "SECalo", -1 ); 
  _cSEReco = AIDAProcessor::histogramFactory(this)->createCloud1D( "SEReco", "SEReco", -1 );
  _cInvMass = AIDAProcessor::histogramFactory(this)->createCloud1D( "InvMass", "InvMass", -1 );
  _cSEMC = AIDAProcessor::histogramFactory(this)->createCloud1D( "SEMC", "SEMC", -1 );
  _cDSERecoMC = AIDAProcessor::histogramFactory(this)->createCloud1D( "DSERecoMC", "DSERecoMC", -1 );
  _cDSERecoCalo = AIDAProcessor::histogramFactory(this)->createCloud1D( "DSERecoCalo", "DSERecoCalo", -1 );
  _cDSECaloMC = AIDAProcessor::histogramFactory(this)->createCloud1D( "DSECaloMC", "DSECaloMC", -1 );
  _cAlphaRecoMC = AIDAProcessor::histogramFactory(this)->createCloud1D( "AlphaRecoMC", "AlphaRecoMC", -1 );
  _cAlphaRecoCalo = AIDAProcessor::histogramFactory(this)->createCloud1D( "AlphaRecoCalo", "AlphaRecoCalo", -1 );
  _cAlphaCaloMC = AIDAProcessor::histogramFactory(this)->createCloud1D( "AlphaCaloMC", "AlphaCaloMC", -1 );

  _cNTracksWithEnergyAssignedByMC = AIDAProcessor::histogramFactory(this)->createCloud1D( "NTracksWithEnergyAssignedByMC", "NTracksWithEnergyAssignedByMC", -1 );
  _cNTracksNotPassingCylinderCuts = AIDAProcessor::histogramFactory(this)->createCloud1D( "NTracksNotPassingCylinderCuts", "NTracksNotPassingCylinderCuts", -1 );
  _cPTracksNotPassingCylinderCuts = AIDAProcessor::histogramFactory(this)->createCloud1D( "PTracksNotPassingCylinderCuts", "PTracksNotPassingCylinderCuts", -1 );;
  _cSumPTracksNotPassingCylinderCuts = AIDAProcessor::histogramFactory(this)->createCloud1D( "SumPTracksNotPassingCylinderCuts", "SumPTracksNotPassingCylinderCuts", -1 );;
  _cSumPTracksNotPassingCylinderCutsvsDSERecoMC = AIDAProcessor::histogramFactory(this)->createCloud2D( "SumPTracksNotPassingCylinderCutsvsDSERecoMC", 
													"SumPTracksNotPassingCylinderCutsvsDSERecoMC", -1 );;
  _hChi2 = AIDAProcessor::histogramFactory(this)->createHistogram1D( "hChi2","hChi2",50,0.0,5.0 );
  _hChi2DEP = AIDAProcessor::histogramFactory(this)->createHistogram2D( "hChi2DEP","hChi2DEP",25,0.0,5.0,100,-50.0,50.0 );
  
  _cNNeutralHitsAssignedToCharged = AIDAProcessor::histogramFactory(this)->createCloud1D( "NNeutralHitsAssignedToCharged", "NNeutralHitsAssignedToCharged", -1 );
  _cENeutralHitsAssignedToCharged = AIDAProcessor::histogramFactory(this)->createCloud1D( "ENeutralHitsAssignedToCharged", "ENeutralHitsAssignedToCharged", -1 );
  _cSumNNeutralHitsAssignedToCharged = AIDAProcessor::histogramFactory(this)->createCloud1D( "SumNNeutralHitsAssignedToCharged", "SumNNeutralHitsAssignedToCharged", -1 );
  _cSumENeutralHitsAssignedToCharged = AIDAProcessor::histogramFactory(this)->createCloud1D( "SumENeutralHitsAssignedToCharged", "SumENeutralHitsAssignedToCharged", -1 );
  _cSumENeutralHitsAssignedToChargedvsDSERecoMC = AIDAProcessor::histogramFactory(this)->createCloud2D( "SumENeutralHitsAssignedToChargedvsDSERecoMC",
													"SumENeutralHitsAssignedToChargedvsDSERecoMC", -1 );

  _cNChargedHitsAssignedToNeutral = AIDAProcessor::histogramFactory(this)->createCloud1D( "NChargedHitsAssignedToNeutral", "NChargedHitsAssignedToNeutral", -1 );
  _cEChargedHitsAssignedToNeutral = AIDAProcessor::histogramFactory(this)->createCloud1D( "EChargedHitsAssignedToNeutral", "EChargedHitsAssignedToNeutral", -1 );
  _cSumNChargedHitsAssignedToNeutrals = AIDAProcessor::histogramFactory(this)->createCloud1D( "SumNChargedHitsAssignedToNeutrals", "SumNChargedHitsAssignedToNeutrals", -1 );
  _cSumEChargedHitsAssignedToNeutrals = AIDAProcessor::histogramFactory(this)->createCloud1D( "SumEChargedHitsAssignedToNeutrals", "SumEChargedHitsAssignedToNeutrals", -1 );
  _cSumEChargedHitsAssignedToNeutralsvsDSERecoMC = AIDAProcessor::histogramFactory(this)->createCloud2D( "SumEChargedHitsAssignedToNeutralsvsDSERecoMC",
													 "SumEChargedHitsAssignedToNeutralsvsDSERecoMC", -1 );

  _cDEdcEwavsDSERecoMC = AIDAProcessor::histogramFactory(this)->createCloud2D( "DEdcEwavsDSERecoMC", "DEdcEwavsDSERecoMC", -1 );
  _cDEdcEwaMinusDSERecoMCvsDSERecoMC = AIDAProcessor::histogramFactory(this)->createCloud2D( "DEdcEwaMinusDSERecoMCvsDSERecoMC", "DEdcEwaMinusDSERecoMCvsDSERecoMC", -1 );


  #endif




  _nRun = -1 ;
  _nEvt = 0 ;


}



void TrackBasedPFlow::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
}



void TrackBasedPFlow::processEvent( LCEvent * evt )
{
  static bool firstEvent = true ;
  
  if(firstEvent==true) {
    std::cout << "TrackBasedPFlow called for first event" << std::endl;
  }

  // extrapolate tracks to calo
  // - get outermost hits per track
  // - fit new helix
 
  // find calohits close to track extrapolation (find MIP stub)
  // - cone-like tube cut around the helix with start point in the last point of the track, take only hits inside the cone 
  //   (cone angle as a funtion of particle momentum?)
  // - calc distance of such a hit to track extrapolation (put in plot?)
  // - if distance is smaller than limit1, assign hit to track
  // - else search furtheron for hits
  // - if no hits found release distance criteria? or at the end assign no hits to track
  // - if a few hits are found (MIP stub) define point where nuclear interaction takes place
  // - start (i) trackwise clustering, (ii) nn-clustering, ... with this point as seed (and perhaps with direction as well)
  

  // divide up energy of showers which are assigned to two or more tracks
  // - first approach: divide up linearly
  // - perform PID (as good as possible at this point)
  // - second approach: devide up by their predicted energy
  // - perform PID (as good as possible at this point)
  // - if ( ( ECalo - n*SigmaECaloIntrinsic ) < (ECalo - (p^2 + m^2)) / ( p^2 + m^2 ) < ( ECalo + n*SigmaECaloIntrinsic ) ) keep keep track and calo entries and
  //   buid RecoParticle and 'remove' hits (mark hits not to be used again)
  // - else 'remove' tracks and use the intrinsic calorimeter resolution

  // on remaining hits perform (i) trackwise clustering, (ii) nn-clustering, ...
  // - perform PID (as good as possible at this point)
  // - 
  // - buid neutral RecoParticles
  // - divide up energy of showers which are assigned to two or more clusters:
  // - first approach: divide up linearly
  // - perform PID (as good as possible at this point)





  // debug
  int NTracks = 0;
  int NTracksWithEnergyAssignedByMC = 0;
  double SEReco  = 0.0;
  double SpxReco = 0.0;
  double SpyReco = 0.0;
  double SpzReco = 0.0;


  _tracksExtrapolatedIntoCalorimeter.clear();

  _tracksNotExtrapolatedIntoCalorimeter.clear();
  _tracksNotFulFillingPtCut.clear();
  _tracksWithTooFewTrackerHits.clear();
  _tracksNotFulfillingCombinedCylinderShellCuts.clear();

  _tracksNotExtrapolatedWithEnoughSiliconAndTPCHits.clear();
  _tracksNotExtrapolatedComingFromIP.clear();
  _tracksDiscarded.clear();

  _tracksWhichWouldReachCalorimeter.clear();
  _tracksNotReachingTheCalorimeter.clear();

  _emShowerCandidatesRecognisedAsCharged.clear();


  int NTracksNotPassingCylinderCuts = 0;
  double SumPTracksNotPassingCylinderCuts = 0.0;

  int nOfChargedObjectsExtrapolatedIntoCalorimeterPerEvent = 0;
  
  std::vector<CalorimeterHit*> neutralCaloHitsAssignedToCharged;
  std::vector<CalorimeterHit*> chargedCaloHitsAssignedToNeutral;

  
  LCCollectionVec* reconstructedParticles = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec* reconstructedClusters = new LCCollectionVec(LCIO::CLUSTER);

  // set flag to store more information in the output file
  LCFlagImpl flag;
  flag.setBit(LCIO::CLBIT_HITS);
  reconstructedClusters->setFlag(flag.getFlag());

  try {
    
    std::vector< std::string >::const_iterator iter;
    const std::vector< std::string >* ColNames = evt->getCollectionNames();

    for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
      
      LCCollection* col = evt->getCollection( *iter ) ;
      
      if ( (col->getTypeName() == LCIO::TRACK) && (*iter == _colNameTracks) ) {

	// FIXME: Here is a problem, if there are more than one Track Collections
	NTracks = col->getNumberOfElements();
	
	// debug
	if ( _debugLevel > 1 ) std::cout << "event " << _nEvt << "  " << "NTracks = " << NTracks << std::endl << std::endl << std::endl;
	
	
	TrackerHitVec outermostTrackerHits;
	
	for(int j=0; j<NTracks; ++j){
	  
	  Track* track = dynamic_cast<Track*>(col->getElementAt(j));

	  const double absP = MarlinUtil::getAbsMomentum(track,_bField);
	  const double absD0 = fabs(track->getD0());
	  const double absZ0 = fabs(track->getZ0());



	  // debug
	  
	  const double* p = MarlinUtil::getMomentum(track,_bField);
	  double cosTh = p[2]/absP;
	  double pt = sqrt( p[0]*p[0] + p[1]*p[1] );


	  /*
	  std::cout << std::endl
		    << "New Track: " << std::endl
		    << "track momentum: |p| = " << absP << "  " << "cosTh = " << cosTh << "  " 
		    << "momentum cut: pt = " << pt << " < " << << _cutOnPt << std::endl;
	  */


	  // debug

	  if ( _debugLevel > 4 )  MarlinUtil::printTrack(track);	        	 	    
	  if (_drawOnCED) {
	    MarlinCED::newEvent(this,0);
	    // MarlinCED::drawMCParticleTree(evt,"MCParticle",0.05,4.0,15.5,50.0,1626.0,2500.0);
	    MarlinCED::drawTrack(track,3,1,0xff0000,3);
	  }
	  


	  bool ptCut = (pt >= _cutOnPt);
	  bool minNTPCHitsReached = false;
	  bool minNNonTPCHitsReached = false;
	  bool minNHitCut = hasTrackSufficientNumberOfHits(track,minNTPCHitsReached,minNNonTPCHitsReached);
	  bool absD0Cut = (absD0 < _absD0Cut);
	  bool absZ0Cut = (absZ0 < _absZ0Cut);
	  
	  if ( ptCut && minNHitCut && absD0Cut &&  absZ0Cut ) {
	    
	    // 1. get outermost hits of the track	  
	    outermostTrackerHits = getOuterTrackerHits(track,_nOfTrackerHitsUsedForExtrapolation);


	    if (!outermostTrackerHits.empty()) {


	      // 2. proceed only if all of the _nOfTrackerHitsUsedForExtrapolation hits of the track fullfill the cut (rMinCutHelixExtrapolation,zMinCutHelixExtrapolation)   

	      // debug
	      if ( _debugLevel > 7 ) {
		std::cout << "Outermost " << _nOfTrackerHitsUsedForExtrapolation << " Trackerhits. " << "Cut at: " 
			  << "rMinCutHelixExtrapolation = " << _rMinCutHelixExtrapolation << "  " << "zMinCutHelixExtrapolation = " << _zMinCutHelixExtrapolation 
			  << std::endl;
	      }
	     
	      bool allOuterTrackerHitsFullfillRPhiCutForHelixExtrapolation = true;
	      for(TrackerHitVec::const_iterator k = outermostTrackerHits.begin(); k != outermostTrackerHits.end(); ++k){
		
		float x = (*k)->getPosition()[0];
		float y = (*k)->getPosition()[1];
		float z = (*k)->getPosition()[2];
		
		float r = sqrt( (x*x) + (y*y) );


		// debug
		if ( _debugLevel > 7 ) {
		  std::cout << k - outermostTrackerHits.begin() << "  " << "x = " << x << "  " << "y = " << y << "  " << "r = " << r << "  " << "z = " << z << "  " 
			    << "allOuterTHitsFullfillRPhiCut = " << allOuterTrackerHitsFullfillRPhiCutForHelixExtrapolation << std::endl;
		}


		if ( ( r < _rMinCutHelixExtrapolation ) && ( fabs(z) < _zMinCutHelixExtrapolation ) ) {
		  
		  allOuterTrackerHitsFullfillRPhiCutForHelixExtrapolation = false;
		  break;

		}
		
	      }
	      

	      // 3. proceed only if _nOfTrackerHitsOutsideCylindricalCut of the _nOfTrackerHitsUsedForExtrapolation hits fullfill the cut 
	      //    (rMinCylindricalCut,zMinCylindricalCut). Based on the ordered collection 'outermostTrackerHits'

	      // debug
	      if ( _debugLevel > 7 ) {
		std::cout << "Outermost " << _nOfTrackerHitsOutsideCylindricalCut << " Trackerhits. " << "Cut at: " 
			  << "rMinCylindricalCut = " <<  _rMinCylindricalCut<< "  " << "zMinCylindricalCut = " <<  _zMinCylindricalCut
			  << std::endl;
	      }


	      bool outerTrackerHitsFullfillCylindricalCut = true;
	      int index = 0;
	      for(TrackerHitVec::const_iterator k = outermostTrackerHits.begin(); k != outermostTrackerHits.end(); ++k){
		       		
		float x = (*k)->getPosition()[0];
		float y = (*k)->getPosition()[1];
		float z = (*k)->getPosition()[2];
		
		float r = sqrt( (x*x) + (y*y) );
		
		// debug
		if ( _debugLevel > 7 ) {
		  std::cout << k - outermostTrackerHits.begin() << "  " << "x = " << x << "  " << "y = " << y << "  " << "r = " << r << "  " << "z = " << z << "  " 
			    << "outerTrackerHitsFullfillCylindricalCut = " << outerTrackerHitsFullfillCylindricalCut << std::endl;
		}


	
		if ( ( r < _rMinCylindricalCut ) && ( fabs(z) < _zMinCylindricalCut ) ) {
		  
		  outerTrackerHitsFullfillCylindricalCut = false;
		  break;
		  
		}
		
		if (index > _nOfTrackerHitsOutsideCylindricalCut) break;

		++index;

	      }
	      
	      
	      if ( allOuterTrackerHitsFullfillRPhiCutForHelixExtrapolation && outerTrackerHitsFullfillCylindricalCut ) {


		// debug
		
		if (_drawOnCED) {
		  int index = 0;
		  for (TrackerHitVec::const_iterator k = outermostTrackerHits.begin(); k != outermostTrackerHits.end(); ++k) {
		    
		    float x = (*k)->getPosition()[0];
		    float y = (*k)->getPosition()[1];
		    float z = (*k)->getPosition()[2];
		    
		    int color = 0;
		    if (index < _nOfTrackerHitsOutsideCylindricalCut) color = 0x9100ff;
		    else color = 0xffffff;
		    
		    ced_hit ( x,y,z, 0 | 4 << CED_LAYER_SHIFT, 3, color );		
		    ++index;

		  }
	      
		}
		
		// prepare arrays for helix fit on the outermost tracker hits
		float* xCoordiantesOfOutermostHits = new float[_nOfTrackerHitsUsedForExtrapolation];
		float* yCoordiantesOfOutermostHits = new float[_nOfTrackerHitsUsedForExtrapolation];
		float* zCoordiantesOfOutermostHits = new float[_nOfTrackerHitsUsedForExtrapolation];
		float* aCoordiantesOfOutermostHits = new float[_nOfTrackerHitsUsedForExtrapolation];




		
		for(int k=0; k<_nOfTrackerHitsUsedForExtrapolation; ++k){
		
		  // transform coordinates, set origin at outermost tracker hit
		  xCoordiantesOfOutermostHits[k] = (float)( (outermostTrackerHits.at(k)->getPosition()[0]) - (outermostTrackerHits.at(0)->getPosition()[0]) );
		  yCoordiantesOfOutermostHits[k] = (float)( (outermostTrackerHits.at(k)->getPosition()[1]) - (outermostTrackerHits.at(0)->getPosition()[1]) );
		  zCoordiantesOfOutermostHits[k] = (float)( (outermostTrackerHits.at(k)->getPosition()[2]) - (outermostTrackerHits.at(0)->getPosition()[2]) );
		  aCoordiantesOfOutermostHits[k] = 0.0;
		  
		  // debug
		  /*
		  std::cout << "x = " << xCoordiantesOfOutermostHits[k] << "  " << "y = " << yCoordiantesOfOutermostHits[k] << "  " << "z = " << zCoordiantesOfOutermostHits[k] 
		            << "  " << "a = " << aCoordiantesOfOutermostHits[k] << std::endl;
		  */
		
		}
		
	      
		// 4. fit helix on outermost hits, if |p| < _cutOnAbsMomentumForFitOnOuterMostTrackerHits

		// canonical parametrisation, parameter vector: (z0,Phi0,omega,d0,tanL)				
		double z0    = 0.0;
		double Phi0  = 0.0;
		double omega = 0.0;
		double d0    = 0.0;
		double tanL  = 0.0;
		
		double chi2  = 0.0;

		LCVector3D referencePoint;



		if ( absP < _cutOnAbsMomentumForFitOnOuterMostTrackerHits ) { 

		  ClusterShapes* shape = new ClusterShapes(_nOfTrackerHitsUsedForExtrapolation,aCoordiantesOfOutermostHits,
							   xCoordiantesOfOutermostHits,yCoordiantesOfOutermostHits,zCoordiantesOfOutermostHits);
	      
		  double par[5];
		  double dpar[5];
		  double distmax;
		  int status = 0;

		  int direction = 1;

		  // FIXME: replace hard coded number, needed to avoid instable fits
		  if ( fabs(track->getOmega()) < 0.0001 ) direction = 1; 
		  else direction = (int)( (track->getOmega())/fabs(track->getOmega()));
		  
		  shape->FitHelix(1000,status,3,par,dpar,chi2,distmax,direction); // canonical parametrisation

		  z0    = par[0];
		  Phi0  = par[1];
		  omega = par[2];
		  d0    = par[3];
		  tanL  = par[4];

		  referencePoint.setX(outermostTrackerHits.at(0)->getPosition()[0]);
		  referencePoint.setY(outermostTrackerHits.at(0)->getPosition()[1]);
		  referencePoint.setZ(outermostTrackerHits.at(0)->getPosition()[2]);

		  delete shape;
		  shape = 0;
		
		  // debug
		  if ( _debugLevel > 1 ) std::cout << "Result of the helix fit on outermost tracker hits: chi2 = " << chi2 << std::endl;
		  
		}
		else {
		  
		  z0    = track->getZ0();
		  Phi0  = track->getPhi();
		  omega = track->getOmega();
		  d0    = track->getD0();
		  tanL  = track->getTanLambda();

		  referencePoint.setX(0.0);
		  referencePoint.setY(0.0);
		  referencePoint.setZ(0.0);

		  // debug
		  if ( _debugLevel > 1 ) std::cout << "Track momentum |p| = " << absP << " larger than cut. Track parameter are used for track extrapolation" << std::endl;

		}

	       

		// Trajectory interface
	
		
		SimpleHelix* fittedHelix = new SimpleHelix(d0,Phi0,omega,z0,tanL,referencePoint);

		LCVector3D firstPointInHelixFit(xCoordiantesOfOutermostHits[_nOfTrackerHitsUsedForExtrapolation-1] 
						+ outermostTrackerHits.at(0)->getPosition()[0],
						yCoordiantesOfOutermostHits[_nOfTrackerHitsUsedForExtrapolation-1]
						+ outermostTrackerHits.at(0)->getPosition()[1],
						zCoordiantesOfOutermostHits[_nOfTrackerHitsUsedForExtrapolation-1]
						+ outermostTrackerHits.at(0)->getPosition()[2] );

		LCVector3D lastPointInHelixFit(xCoordiantesOfOutermostHits[0]
					       + outermostTrackerHits.at(0)->getPosition()[0],
					       yCoordiantesOfOutermostHits[0]
					       + outermostTrackerHits.at(0)->getPosition()[1],
					       zCoordiantesOfOutermostHits[0]
					       + outermostTrackerHits.at(0)->getPosition()[2] );
		  
		double sOfFirstPointInHelixFit = getPathLengthOnHelix(firstPointInHelixFit,fittedHelix);
		double sOfLastPointInHelixFit  = getPathLengthOnHelix(lastPointInHelixFit,fittedHelix);

		
		
		fittedHelix->setStart(sOfFirstPointInHelixFit);
		fittedHelix->setEnd(sOfLastPointInHelixFit+_maximalConeTubeLength);
		
		

		// debug
				
		if (_drawOnCED) {
		  
		  // MarlinCED::drawTrajectory(fittedHelix,0 | 3 << CED_LAYER_SHIFT,1,0xff75fa);


		  // self made draw helix (for debugging of MarlinCED::drawHelix() or MarlinCED::drawTrajectory() method)		 		  

		  if ( _debugLevel > 1 ) std::cout << "sOfFirstPointInHelixFit = " << sOfFirstPointInHelixFit << "  "
						   << "sOfLastPointInHelixFit = " << sOfLastPointInHelixFit << "  "
						   << "maxPathLengthInConeTube = " << sOfLastPointInHelixFit+_maximalConeTubeLength << std::endl;	     
		  
		  LCVector3D firstPointInHelixFitProjected = getProjectedPointOnHelix(firstPointInHelixFit,fittedHelix);
		  LCVector3D lastPointInHelixFitProjected  = getProjectedPointOnHelix(lastPointInHelixFit,fittedHelix);

		  
		  ced_hit ( firstPointInHelixFit.x(),firstPointInHelixFit.y(),firstPointInHelixFit.z(), 0 | 1 << CED_LAYER_SHIFT, 4, 0x60ff4b );
		  ced_hit ( firstPointInHelixFitProjected.x(),firstPointInHelixFitProjected.y(),firstPointInHelixFitProjected.z(), 
			    0 | 2 << CED_LAYER_SHIFT, 4, 0xf7ff86 );
		  ced_hit ( lastPointInHelixFit.x(),lastPointInHelixFit.y(),lastPointInHelixFit.z(), 0 | 1 << CED_LAYER_SHIFT, 4, 0x60ff4b );
		  ced_hit ( lastPointInHelixFitProjected.x(),lastPointInHelixFitProjected.y(),lastPointInHelixFitProjected.z(), 
			    0 | 2 << CED_LAYER_SHIFT, 4, 0xf7ff86 );

		  /*
		  std::cout << "first point in helix fit: " << std::endl
			    << "x = " << "(" << firstPointInHelixFit.x() << "," <<  firstPointInHelixFit.y() << "," <<  firstPointInHelixFit.z() << ")" << "  "
			    << "xproj = " << "(" << firstPointInHelixFitProjected.x() << "," <<  firstPointInHelixFitProjected.y() << "," <<  firstPointInHelixFit.z() 
			    << ")" << std::endl << std::endl; 		  
		  std::cout << "last point in helix fit: " << std::endl
			    << "x = " << "(" << lastPointInHelixFit.x() << "," <<  lastPointInHelixFit.y() << "," <<  lastPointInHelixFit.z() << ")" << "  "
			    << "xproj = " << "(" << lastPointInHelixFitProjected.x() << "," <<  lastPointInHelixFitProjected.y() << "," <<  lastPointInHelixFit.z() 
			    << ")" << std::endl << std::endl;
		  */
   		  
		  
		  double step = 10.0; // mm
		  int nOfSteps = (int)std::ceil(_maximalConeTubeLength/step);
		  
		  for(int k=0; k<nOfSteps; ++k) {		    
		    
		    LCVector3D pointAtS   = fittedHelix->getPosition(k*step+sOfFirstPointInHelixFit);
		    LCVector3D pointAtSp1 = fittedHelix->getPosition((k+1)*step+sOfFirstPointInHelixFit);
		    
		    ced_line( pointAtS.x(),pointAtS.y(),pointAtS.z(), pointAtSp1.x(),pointAtSp1.y(),pointAtSp1.z(), 0 | 6 << CED_LAYER_SHIFT, 1,0xffffff );

		  }		  
		  
		  


		  // draw hits and projections
		  for(int k=_nOfTrackerHitsUsedForExtrapolation-1; k>=0; --k) {
		    
		    LCVector3D trackerHit(xCoordiantesOfOutermostHits[k] + outermostTrackerHits.at(0)->getPosition()[0],
					  yCoordiantesOfOutermostHits[k] + outermostTrackerHits.at(0)->getPosition()[1],
					  zCoordiantesOfOutermostHits[k] + outermostTrackerHits.at(0)->getPosition()[2]);
		    double sOfTrackerHit = getPathLengthOnHelix(trackerHit,fittedHelix);
		    LCVector3D trackerHitProjected = getProjectedPointOnHelix(trackerHit,fittedHelix);
		    double sOfTrackerHitProjected = getPathLengthOnHelix(trackerHitProjected,fittedHelix);
		    LCVector3D trackerHitRecalc = fittedHelix->getPosition(sOfTrackerHit);
		    LCVector3D trackerHitProjectedRecalc = fittedHelix->getPosition(sOfTrackerHitProjected);

		    // debug
		    /*
		    std::cout << "Hit = " << "(" << trackerHit.x() << "," <<  trackerHit.y() << "," <<  trackerHit.z() << ")" << "  " << "sOfTrackerHit = " << sOfTrackerHit 
			      << "  " << "Hit(sOfTrackerHit) = " << "(" << trackerHitRecalc.x() << "," <<  trackerHitRecalc.y() << "," <<  trackerHitRecalc.z() << ")" 
			      << "  " << std::endl
			      << "Hitproj = " << "(" << trackerHitProjected.x() << "," <<  trackerHitProjected.y() << "," <<  trackerHitProjected.z() << ")" 
			      << "  " << "sOfTrackerHitProjected = " << sOfTrackerHitProjected 
			      << "  " << "Hit(sOfTrackerHitProjected) = " << "(" << trackerHitProjectedRecalc.x() << "," <<  trackerHitProjectedRecalc.y() << "," 
			      << trackerHitProjectedRecalc.z() << ")"
			      << std::endl << std::endl;
		    */

		    ced_hit ( trackerHit.x(),trackerHit.y(),trackerHit.z(), 0 | 1 << CED_LAYER_SHIFT, 2, 0xff0004 );
		    ced_hit ( trackerHitProjected.x(),trackerHitProjected.y(),trackerHitProjected.z(), 0 | 1 << CED_LAYER_SHIFT, 2, 0x17ff06 );
		    ced_line( trackerHit.x(),trackerHit.y(),trackerHit.z(),
			      trackerHitProjected.x(),trackerHitProjected.y(),trackerHitProjected.z(), 0 | 1 << CED_LAYER_SHIFT, 1, 0xffffff );
      
		  }


		  
		  // manually draw helix only within the _nOfTrackerHitsUsedForExtrapolation Tracker Hits
		  int nOfStepsHelixExtrapolatedTrackerHits = 4*_nOfTrackerHitsUsedForExtrapolation; // hard coded four times more steps than _nOfTrackerHitsUsedForExtrapolation
		  double ds = (sOfLastPointInHelixFit - sOfFirstPointInHelixFit)/nOfStepsHelixExtrapolatedTrackerHits;
		  
		  for(int k=0; k<nOfStepsHelixExtrapolatedTrackerHits; ++k){		    
		    
		    LCVector3D pointAtS   = fittedHelix->getPosition(k*ds+sOfFirstPointInHelixFit);
		    LCVector3D pointAtSp1 = fittedHelix->getPosition((k+1)*ds+sOfFirstPointInHelixFit);

		    ced_line( pointAtS.x(),pointAtS.y(),pointAtS.z(), pointAtSp1.x(),pointAtSp1.y(),pointAtSp1.z(), 0 | 5 << CED_LAYER_SHIFT, 1, 0xffca67 );

		  }
		  
		    
				  
		}
		
		// debug
		
		if ( _debugLevel > 1 ) std::cout << "Fit (d0,z0,phi0,omega,tanL) = " << "(" << d0 << "," << z0 << "," << Phi0 << "," << omega << "," << tanL << ")" 
						 << std::endl;

		
		SimpleHelix* helixTrajectory = dynamic_cast<SimpleHelix*>(fittedHelix);
		if ( (_debugLevel > 1) && (helixTrajectory != 0)) helixTrajectory->printProperties(); // check if dynamic_cast has been successful
		
		
		
		
		  
		// also used to cut off GEANT4 bugs
		drawRelatedCalorimeterHits(evt,track);
		// debug
		// FIXME: to take out GEANT4 bugs in this way takes far to much time, therefore it is done only for _drawOnCED; find an other way
		if (_drawOnCED){ 
		  drawEMShowerCandidates(evt);		  
		}




		// 5. collect Calorimeter hits in a cone-like tube around the helix extrapolation		
		
		const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes = getRelatedCalorimeterHits(evt,outermostTrackerHits,fittedHelix,
															   j/*only for debugging*/);
		


		// 6. collect MIP like stub, if exists, and perform trackwise clustering => result are clusters related to this track
		
		ClusterImplWithAttributes* mipStub = new ClusterImplWithAttributes();
		std::vector<ClusterImplWithAttributes*> clusters;
		
		getRelatedClusters(calorimeterHitsWithAttributes,outermostTrackerHits,fittedHelix,mipStub,clusters);
		//getRelatedClusterPerfectly(evt,track,outermostTrackerHits,fittedHelix,clusters);




		// 7. assign clusters to reconstructed particle (always assign the MIP stub, if it exists, and assign other clusters by proximity criteria)
		
		ReconstructedParticleImpl* recoParticle = assignClustersToTrack(track,mipStub,clusters,outermostTrackerHits,fittedHelix);
		//ReconstructedParticleImpl* recoParticle = assignClustersToTrackPerfectly(track,clusters);
		



		// FIXME: very simple approach to find mistakes in the clustering and im the assignment of clusters
		// simply redo the full procedure with smaller resolution parameters of the trackwise clustering
		// this can also been done iteratively!!!
	      	      		
		double energyAssignedToRecoParticle = 0.0;
		const ClusterVec ClustersAssignedToRecoParticle = recoParticle->getClusters();
		int nOfClustersAssignedToRecoParticle = ClustersAssignedToRecoParticle.size();
 
		for (int iOfClustersAssignedToRecoParticle = 0; iOfClustersAssignedToRecoParticle<nOfClustersAssignedToRecoParticle; ++iOfClustersAssignedToRecoParticle) {
		  
		  Cluster* cl = recoParticle->getClusters()[iOfClustersAssignedToRecoParticle];
		   energyAssignedToRecoParticle+= cl->getEnergy();
		  
		}


		// get track energy for comparison
		double energyOfTrack = sqrt( absP*absP + 0.13957*0.13957); // FIXME: at the moment the pion hypothesis is used

		if ( fabs(energyAssignedToRecoParticle-energyOfTrack) > (0.6*sqrt(energyOfTrack)) ) {
		  
		  if ( _debugLevel > 5 ) {
		  
		    std::cout << "RECLUSTERING => => =>" << std::endl
			      << "absP of track: " << absP << "  " << "energy of track: " << energyOfTrack << "  " 
			      << "energyAssignedToRecoParticle: " << energyAssignedToRecoParticle << std::endl
			      << "fabs(energyAssignedToRecoParticle-energyOfTrack) = " << fabs(energyAssignedToRecoParticle-energyOfTrack) << "  "
			      << "0.6*sqrt(energyOfTrack) = " << 0.6*sqrt(energyOfTrack) << std::endl
			      << "Full Condition: fabs(energyAssignedToRecoParticle-energyOfTrack) > (0.6*sqrt(energyOfTrack)) : " 
			      << ( fabs(energyAssignedToRecoParticle-energyOfTrack) > (0.6*sqrt(energyOfTrack)) ) << std::endl << std::endl;
		    
		  }
		  
		  delete mipStub;
		  mipStub = 0;

		  mipStub = new ClusterImplWithAttributes();
		  
		  clusters.clear();

		  delete recoParticle;
		  recoParticle = 0;
		  

		  
		  float weightForResoStore = _weightForReso;
		  _weightForReso = 5.0;

		  for(std::vector<float>::iterator k = _resolutionParameter.begin(); k != _resolutionParameter.end(); ++k) {

		    // FIXME: hard coded number
		    (*k) -= 10.0;
		    
		     if ( _debugLevel > 5 ) {
		       
		       std::cout << "resolutionParameter: " << (*k) << std::endl;
		       
		     }

		  }

		  if (_drawOnCED) {
		    MarlinCED::newEvent(this,0);
		    // MarlinCED::drawMCParticleTree(evt,"MCParticle",0.05,4.0,15.5,50.0,1626.0,2500.0);
		    MarlinCED::drawTrack(track,3,1,0xff0000,3);
		  }


		  getRelatedClusters(calorimeterHitsWithAttributes,outermostTrackerHits,fittedHelix,mipStub,clusters);
		  //getRelatedClusterPerfectly(evt,track,outermostTrackerHits,fittedHelix,clusters);

		  recoParticle = assignClustersToTrack(track,mipStub,clusters,outermostTrackerHits,fittedHelix);
		  //ReconstructedParticleImpl* recoParticle = assignClustersToTrackPerfectly(track,clusters);


		  _weightForReso = weightForResoStore;

		  for(std::vector<float>::iterator k = _resolutionParameter.begin(); k != _resolutionParameter.end(); ++k) {
		    
		    // FIXME: hard coded number
		    (*k) += 10.0;

		  }

		}



		// fill cluster collection	      
		for(ClusterVec::const_iterator k = recoParticle->getClusters().begin(); k != recoParticle->getClusters().end(); ++k) {

		  reconstructedClusters->addElement(*k);

		}



	      	      	  


		// in between: calculate the MIP-Stub finding efficiency
			      
		bool realMIPStub = isRealMIPStub(evt,track,outermostTrackerHits,fittedHelix);

		// hard-coded 5, see isRealMIPStub()
		bool foundMIPStub = ( (mipStub->getClusterImpl()->getCalorimeterHits().size()) >= 5 ) && ( (mipStub->getEnergy()) > 0.0 );
		
		
		if (realMIPStub) ++_nOfRealMIPStubs;
		
		if (foundMIPStub) ++_nOfFoundMIPStubs;

		  
		// debug
		/*		
		  std::cout << nOfChargedObjectsExtrapolatedIntoCalorimeterPerEvent << "  " << "REAL MIP STUB: " << realMIPStub << "  " << "FOUND MIP STUB: " 
		  << foundMIPStub << "  " << "pt = " << pt << "  " << "cosTh = " << cosTh << std::endl;
		*/
		  
		++nOfChargedObjectsExtrapolatedIntoCalorimeterPerEvent;
		++_nOfChargedObjectsExtrapolatedIntoCalorimeter;
		  
		
		
                #ifdef MARLIN_USE_AIDA
		if (realMIPStub) {
		  
		  _hPtRealMIPStubs->fill(pt);
		  _cPtRealMIPStubs->fill(pt);
		  _cPtNRealMIPStubs->fill(pt,1);
		  
		  _hCosThRealMIPStubs->fill(cosTh);
		  _cCosThRealMIPStubs->fill(cosTh);
		  _cCosThNRealMIPStubs->fill(cosTh,1);
		    
		}
		else {
		    
		  _cPtNRealMIPStubs->fill(pt,0);
		  _cCosThNRealMIPStubs->fill(cosTh,0);
		  
		}
		if (foundMIPStub) {
		    
		  _hPtFoundMIPStubs->fill(pt);
		  _cPtFoundMIPStubs->fill(pt);
		  _cPtNFoundMIPStubs->fill(pt,1);
		  
		  _hCosThFoundMIPStubs->fill(cosTh);
		  _cCosThFoundMIPStubs->fill(cosTh);
		  _cCosThNFoundMIPStubs->fill(cosTh,1);
		  
		}
		else {
		  
		  _cPtNFoundMIPStubs->fill(pt,0);
		  _cCosThNFoundMIPStubs->fill(cosTh,0);
		  
		}
                #endif

		

	       
		// 8. perform PID and add to reco particle collection

		doPID(recoParticle,mipStub->isMuon());

		
		// debug
		if (_drawOnCED) {

		  MarlinCED::drawSimCalorimeterHits(evt,0,1,0xff0e12,1);
		  MarlinCED::drawCalorimeterHits(evt,2,4,0x22ff5d,1);
	    
		}

		if ( _debugLevel > 1 ) MarlinUtil::printRecoParticle(recoParticle,_bField); // to check the CED output


		// debug, draw clusters assigned to track
		if (_drawOnCED) {
		  
		  ClusterVec reconstructedClusters = recoParticle->getClusters();
		  for(ClusterVec::const_iterator i = reconstructedClusters.begin(); i != reconstructedClusters.end(); ++i) MarlinCED::drawCluster((*i),2,10,0xff0004,9); 
        
		}
  

		

                #ifdef MARLIN_USE_AIDA

		double recoParticleEnergyInCluster = 0.0;
		const ClusterVec Clusters = recoParticle->getClusters();
		int nOfClusters = Clusters.size();
		
		for (int iOfClusters = 0; iOfClusters<nOfClusters; ++iOfClusters) {
		  
		  Cluster* cluster = recoParticle->getClusters()[iOfClusters];
		   recoParticleEnergyInCluster += cluster->getEnergy();
		  
		}


		double m2RecoParticle = recoParticleEnergyInCluster*recoParticleEnergyInCluster - ( (recoParticle->getMomentum()[0])*(recoParticle->getMomentum()[0]) + 
												    (recoParticle->getMomentum()[1])*(recoParticle->getMomentum()[1]) +
												    (recoParticle->getMomentum()[2])*(recoParticle->getMomentum()[2]) );
		// debug
		// std::cout << chi2 << "  " << m2RecoParticle << std::endl;

		_hChi2->fill(chi2);
		_hChi2DEP->fill(chi2,m2RecoParticle);	       
                #endif
		


		// check plots
		int nOfNeutralHits = 0;
		double eOfNeutralHits = 0.0;
		std::vector<CalorimeterHit*> neutralCaloHits = getNeutralHitsAssignedToChargedParticle(evt,recoParticle,nOfNeutralHits,eOfNeutralHits,0.5);
		
                #ifdef MARLIN_USE_AIDA
		_cNNeutralHitsAssignedToCharged->fill(nOfNeutralHits);
		_cENeutralHitsAssignedToCharged->fill(eOfNeutralHits);
                #endif

		// debug
		// std::cout << "NEUTRAL HITS: n = " << nOfNeutralHits << "  " << "E = " << eOfNeutralHits << std::endl;
		
		for (std::vector<CalorimeterHit*>::const_iterator i = neutralCaloHits.begin(); i != neutralCaloHits.end(); ++i) {
		  
		  std::vector<CalorimeterHit*>::const_iterator position = find(neutralCaloHitsAssignedToCharged.begin(),neutralCaloHitsAssignedToCharged.end(),(*i));
		  
		  if ( position == neutralCaloHitsAssignedToCharged.end() ) neutralCaloHitsAssignedToCharged.push_back(*i);
		  
		}
		neutralCaloHits.clear();
		
		if ( _drawOnCED ) {
		  
		  for (std::vector<CalorimeterHit*>::const_iterator i = neutralCaloHitsAssignedToCharged.begin(); i != neutralCaloHitsAssignedToCharged.end(); ++i) {
		    
		    ced_hit ( (*i)->getPosition()[0], (*i)->getPosition()[1], (*i)->getPosition()[2], 0 | 0 << CED_LAYER_SHIFT, 4, 0xff12e3 );
		    
		  }
		  
		}




		// 9. add recoParticle to collection		
		reconstructedParticles->addElement(recoParticle);
		
		// track is added to the _tracksExtrapolatedIntoCalorimeter collection assuming that the assignment of energy was successful
		_tracksExtrapolatedIntoCalorimeter.push_back(track);

		
		delete fittedHelix;
		fittedHelix = 0;	       		
		
		delete[] xCoordiantesOfOutermostHits;
		xCoordiantesOfOutermostHits = 0;
		delete[] yCoordiantesOfOutermostHits;
		yCoordiantesOfOutermostHits = 0;
		delete[] zCoordiantesOfOutermostHits;
		zCoordiantesOfOutermostHits = 0;
		delete[] aCoordiantesOfOutermostHits;
		aCoordiantesOfOutermostHits = 0;     
		
	      }
	      else { // track does not fulfil 'combined cylinder shell cuts'
		
		if ( _debugLevel > 5 ) {
		  
		std::cout << "Track has NOT been extrapolated ==> Track is added to _tracksNotExtrapolated collection" << std::endl 
			  << "Reason: Track does not fulfil 'combined cylinder shell cuts'" << std::endl
			  << "allOuterTrackerHitsFullfillRPhiCutForHelixExtrapolation = " << allOuterTrackerHitsFullfillRPhiCutForHelixExtrapolation << "  "
			  << "outerTrackerHitsFullfillCylindricalCut = " << outerTrackerHitsFullfillCylindricalCut << std::endl
			  << "Full Condition: (allOuterTrackerHitsFullfillRPhiCutForHelixExtrapolation && outerTrackerHitsFullfillCylindricalCut) = "
			  << (allOuterTrackerHitsFullfillRPhiCutForHelixExtrapolation && outerTrackerHitsFullfillCylindricalCut) << std::endl << std::endl;
		
		}
		

		_tracksNotExtrapolatedIntoCalorimeter.push_back(track);
		_tracksNotFulfillingCombinedCylinderShellCuts.push_back(track);


		++NTracksNotPassingCylinderCuts;

                #ifdef MARLIN_USE_AIDA
		_cPTracksNotPassingCylinderCuts->fill(absP);
                #endif
		SumPTracksNotPassingCylinderCuts += absP;
	      
	      }
	      
	    }
	    else { // too few tracker hits

	      if ( _debugLevel > 5 ) {
		
		std::cout << "Track has NOT been extrapolated ==> Track is added to _tracksNotExtrapolated collection" << std::endl 
			  << "Reason: too few tracker hits (" << track->getTrackerHits().size() << "). " << track->getTrackerHits().size() 
			  << " < " << _nOfTrackerHitsUsedForExtrapolation << std::endl << std::endl;
		
	      }

	      _tracksNotExtrapolatedIntoCalorimeter.push_back(track);
	      _tracksWithTooFewTrackerHits.push_back(track);
	      
	    }
	    
	  }
	  else { // cut on pt



	  bool ptCut = (pt >= _cutOnPt);
	  bool minNTPCHitsReached = false;
	  bool minNNonTPCHitsReached = false;
	  bool minNHitCut = hasTrackSufficientNumberOfHits(track,minNTPCHitsReached,minNNonTPCHitsReached);
	  bool absD0Cut = (absD0 < _absD0Cut);
	  bool absZ0Cut = (absZ0 < _absZ0Cut);
	  


	    if ( _debugLevel > 5 ) {
 
	      std::cout << "Track has NOT been extrapolated ==> Track is added to _tracksNotExtrapolated collection" << std::endl 
			<< "Reason: pt = " << pt << " < " << _cutOnPt << "  " << "ptCut = " << ptCut << std::endl
			<< "        minNTPCHitsReached = " << minNTPCHitsReached << "  " 
			<< "minNNonTPCHitsReached = " << minNNonTPCHitsReached << "  " << "minNHitCut = " << minNHitCut << std::endl
			<< "        absD0 = " << absD0 << " < " << _absD0Cut << "  " << "absZ0 = " << absZ0 << " < " << _absZ0Cut
			<< std::endl
			<< "Full Condition: ( ptCut && minNHitCut && absD0Cut &&  absZ0Cut ) = " 
			<< ( ptCut && minNHitCut && absD0Cut &&  absZ0Cut ) << std::endl << std::endl;

	    }

	    _tracksNotExtrapolatedIntoCalorimeter.push_back(track);
	    if (ptCut)_tracksNotFulFillingPtCut.push_back(track);

	  }


	  // debug
	  if ( _drawOnCED ) MarlinCED::draw(this,true);


	} // loop over tracks
	
      }
      
    }
    


    
      
    double pSumOfTracksNotExtrapolatedIntoCalorimeter = 0.0;
      
    for ( std::vector<Track*>::const_iterator i = _tracksNotExtrapolatedIntoCalorimeter.begin(); i != _tracksNotExtrapolatedIntoCalorimeter.end(); ++i ) {
      
      pSumOfTracksNotExtrapolatedIntoCalorimeter += MarlinUtil::getAbsMomentum((*i),_bField);
      
    }
      

      
    // debug
    if ( _debugLevel > 5 ) {
      std::cout << "Failed to extrapolate "<< _tracksNotExtrapolatedIntoCalorimeter.size() << " tracks with a momentum sum of " 
		<< pSumOfTracksNotExtrapolatedIntoCalorimeter << " into the calorimeter" << std::endl << std::endl
		<< "Loop over these tracks:" << std::endl << std::endl;
      
    }
    
    

    // have a look at the charged particles, which do not reach the calorimeter 
    // FIXME: A realistic kink, V0, pair production and 'hard delta ray' finding is needed

    double pSumOfTracksWithReassignedEnergyByMC = 0.0;
    double ESumOfTracksWithReassignedEnergyByMC = 0.0;
      
    for ( std::vector<Track*>::const_iterator i = _tracksNotExtrapolatedIntoCalorimeter.begin(); i != _tracksNotExtrapolatedIntoCalorimeter.end(); ++i ) {

      bool minNTPCHitsReached = false;
      bool minNNonTPCHitsReached = false;;

      bool trackHasSufficientNumberOfHits = hasTrackSufficientNumberOfHits(*i,minNTPCHitsReached,minNNonTPCHitsReached);
    
      // discard track if it has to few TPC or Silicon Tracker Hits
      if ( trackHasSufficientNumberOfHits ) {
      
	const double absD0ofNotExtrapolatedTrack = fabs((*i)->getD0());
	const double absZ0ofNotExtrapolatedTrack = fabs((*i)->getZ0());

	if ( (absD0ofNotExtrapolatedTrack < _absD0Cut) && (absZ0ofNotExtrapolatedTrack < _absZ0Cut) ) {

	  // debug
	  if ( _debugLevel > 5 ) std::cout << "Track will be discarded, particle will be measured as neutral." << std::endl;
	  _tracksDiscarded.push_back(*i);




	  /*
	  // took out the assignment of energy to tracks by MC, particles will be measured as neutrals
	  
	  // debug
	  if ( _debugLevel > 5 ) std::cout << "Energy of the following track will be assigned by MC information: " << std::endl;


	  // FIXME: there seems to be a bug in the filling of _tracksWhichWouldReachCalorimeter, use method 'getRelatedCalorimeterHitsPerfectly()' instead
	  // But this takes more computing time
	  // std::vector<Track*>::const_iterator index = find(_tracksWhichWouldReachCalorimeter.begin(),_tracksWhichWouldReachCalorimeter.end(),(*i));

	  ClusterImpl* clusterRealEnergy = new ClusterImpl();
	  ClusterImpl* clusterPerfectEnergy = new ClusterImpl();
	
	  // also used to cut off GEANT4 bugs
	  // FIXME: to take out GEANT4 bugs in this way takes far to much time, therefore it is done only for _drawOnCED; find an other way
	  getRelatedCalorimeterHitsPerfectly(evt,(*i),clusterRealEnergy,clusterPerfectEnergy);


	  ReconstructedParticleImpl* recoParticle = new ReconstructedParticleImpl();
	    
	  recoParticle->addTrack(*i);
	  recoParticle->addCluster(clusterRealEnergy);
	  // do not delete clusters which are assigned to an reconstructed particle
	  
	  doPID(recoParticle,false);
	  
	  reconstructedClusters->addElement(clusterRealEnergy);
	  
	  reconstructedParticles->addElement(recoParticle);
	  

	  ++NTracksWithEnergyAssignedByMC;

	  pSumOfTracksWithReassignedEnergyByMC += MarlinUtil::getAbsMomentum((*i),_bField);
	  
	 
	  int nCaloHitsAssigned = clusterRealEnergy->getCalorimeterHits().size();
	  for (int j = 0; j < nCaloHitsAssigned; ++j) {

	    CalorimeterHit* caloHit = clusterRealEnergy->getCalorimeterHits().at(j);
	  
	    ESumOfTracksWithReassignedEnergyByMC += caloHit->getEnergy();
	    
	  }
	  
	  // FIXME: perhaps put these clusters into a separate lcio collection
	  // delete clusters which are NOT assigned to any reconstructed particle
      
	  //  delete clusterRealEnergy;
	  //  clusterRealEnergy = 0;
  
	  delete clusterPerfectEnergy;
	  clusterPerfectEnergy = 0;
	  
	  */
	
	  
	  /*
	    bool isTrackReachingCalorimeter = false;
	    if ( clusterRealEnergy->getCalorimeterHits().size() > 0 ) isTrackReachingCalorimeter = true;
	    
	    // take only tracks into account which would not reach the Calorimeter and which fullfill the d0, z0, number of TPC hits and number of non TPC hits cuts 
	    
	    if ( (!isTrackReachingCalorimeter) && (absD0 < _absD0Cut) && (absZ0 < _absZ0Cut) && 
	    minNTPCHitsReached && minNNonTPCHitsReached ) {
	    
	    _tracksNotReachingTheCalorimeter.push_back(*i);	
	    
	    if ( _debugLevel > 5 ) {
	
	    std::cout << "Track with no energy deposition in Calorimeter found. Added to the collection of reconstructed particles..." << std::endl;
	    MarlinUtil::printTrack(*i);
	    
	    }
	    
	    if (_drawOnCED) {
	    
	    MarlinCED::newEvent(this,0);
	    MarlinCED::drawMCParticleTree(evt,"MCParticle",0.05,4.0,15.5,50.0,1626.0,2500.0);
	    MarlinCED::drawTrack((*i),3,1,0xff0000,7);
	    drawRelatedCalorimeterHits(evt,(*i));
	    MarlinCED::draw(this,true);
	    
	    }
	    
	    
	
	    ReconstructedParticleImpl* recoParticle = new ReconstructedParticleImpl();
	    
	    recoParticle->addTrack(*i);
	    doPID(recoParticle,false);
	    
	    reconstructedParticles->addElement(recoParticle);
		
	    
	    }
	  */


	
	  
	  _tracksNotExtrapolatedComingFromIP.push_back(*i);
	  
	}
	else { // track (not extrapolated) is not comming from IP
	  
	  if ( _debugLevel > 5 ) {
	    
	    std::cout << "Track (not extrapolated) is not coming from IP ==> Track will be discarded" << std::endl 
		      << "Reason: absD0 = " << absD0ofNotExtrapolatedTrack << " < " << _absD0Cut << "  " << "absZ0 = " << absZ0ofNotExtrapolatedTrack << " < " 
		      << _absZ0Cut << std::endl
		      << "Full Condition: (absD0 < absD0Cut) && (absZ0 < absZ0Cut) = "
		      << ( (absD0ofNotExtrapolatedTrack < _absD0Cut) && (absZ0ofNotExtrapolatedTrack < _absZ0Cut) ) << std::endl << std::endl;
	    
	  }
	  
	  _tracksDiscarded.push_back(*i);

	}
	
	
	_tracksNotExtrapolatedWithEnoughSiliconAndTPCHits.push_back(*i);

      }
      else { // too few TPC or Silicon Tracker Hits of this not extrapolated track 

	if ( _debugLevel > 5 ) {

	  std::cout << "Track (not extrapolated) has too few TPC or Silicon Tracker hits ==> Track will be discarded" << std::endl 
		    << "Reason: minNTPCHitsReached = " << minNTPCHitsReached << "  " << "minNSiliconHitsReached = " << minNNonTPCHitsReached << std::endl
		    << "Full Condition: (minNTPCHitsReached && minNNonTPCHitsReached) = " << ( minNTPCHitsReached && minNNonTPCHitsReached ) << std::endl << std::endl;
	  
	}

	_tracksDiscarded.push_back(*i);
	
      }
            
      // debug  
      if (_drawOnCED) {
	MarlinCED::newEvent(this,0);   
	MarlinCED::drawTrack((*i),3,1,0xff0000,3);
	MarlinCED::draw(this,true);
      }
	        
    }



    
    // debug
    if ( _debugLevel > 5 ) {
      std::cout << "Recovered "<< NTracksWithEnergyAssignedByMC << " tracks where the energy is assigned by MC. Total momentum sum: " 
		<< pSumOfTracksWithReassignedEnergyByMC << "  " << "Total energy sum: " << ESumOfTracksWithReassignedEnergyByMC << std::endl << std::endl;
      
    }
    
    

















    /*
    
    // debug
    if ( _debugLevel > 5 ) {
    
      std::cout << "Track Summary:" << std::endl 
		<< "N: " << NTracks << "  " << "N failed to extrapolate: " << _tracksNotExtrapolatedIntoCalorimeter.size() << "  "
		<< "N would reach Calo: " << _tracksWhichWouldReachCalorimeter.size() << "  " 
		<< "N not reaching Calo: " << _tracksNotReachingTheCalorimeter.size() << std::endl;
      
      
      std::cout << "TracksWhichWouldReachCalorimeter:" << std::endl;
      for ( std::vector<Track*>::const_iterator i = _tracksWhichWouldReachCalorimeter.begin(); i != _tracksWhichWouldReachCalorimeter.end(); ++i ) {
	MarlinUtil::printTrack(*i,_bField);
      }
      
      std::cout << "TracksNotReachingTheCalorimeter:" << std::endl;
      for ( std::vector<Track*>::const_iterator i = _tracksNotReachingTheCalorimeter.begin(); i != _tracksNotReachingTheCalorimeter.end(); ++i ) {
	MarlinUtil::printTrack(*i,_bField);
      }
      
    }
    
    */









    // carry on with the remaining neutral particles
    
    std::vector<CalorimeterHitWithAttributes*> remainingCalorimeterHits;
    
    for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
      
      LCCollection* col = evt->getCollection( *iter ) ;
      
      if ( (col->getTypeName() == LCIO::CALORIMETERHIT) && ( (*iter == _colNameECAL) || (*iter == _colNameHCAL) ) ) {
	       

	// 10. remove assigned hits to charged particles from copy of calorimeter hit collection
	
	int NCalorimeterHits = col->getNumberOfElements();

	for(int j=0; j<NCalorimeterHits; ++j){
	  
	  CalorimeterHit* calorimeterHit = dynamic_cast<CalorimeterHit*>(col->getElementAt(j));
       
	  bool found = false;
	  
	  for (LCCollectionVec::const_iterator k = reconstructedParticles->begin(); k != reconstructedParticles->end(); ++k) {

	    ReconstructedParticle* recoParticle = dynamic_cast<ReconstructedParticle*>( (*k) );
	  
	    if ( recoParticle->getTracks().size() != 0 ) {
	      
	      for (ClusterVec::const_iterator l = recoParticle->getClusters().begin(); l != recoParticle->getClusters().end(); ++l) {
		
		for (CalorimeterHitVec::const_iterator m = (*l)->getCalorimeterHits().begin(); m != (*l)->getCalorimeterHits().end(); ++m) {
		  
		  if ( calorimeterHit == (*m) ) { 
		    
		    found = true;
		    break;
		    
		  }
		  
		}
		
		if (found) break;
		
	      }

	    }

	    if (found) break;

	  }

	  if (!found) {

	    CalorimeterHitWithAttributes* calorimeterHitWithAttributes = new CalorimeterHitWithAttributes(calorimeterHit,0.0,0.0);
	    remainingCalorimeterHits.push_back(calorimeterHitWithAttributes);

	  }

	}
	
      }

    }



    
    // 11. perform a trackwise clustering on the remaining hits 
    double startPoint[3] = {0.0,0.0,0.0}; // not used at the moment
    double startDirection[3] = {0.0,0.0,0.0}; // not used at the moment

    std::vector<ClusterImplWithAttributes*> remainingClusters = doTrackwiseClusteringForNeutrals(remainingCalorimeterHits,startPoint,startDirection);
       
    for (std::vector<ClusterImplWithAttributes*>::const_iterator i = remainingClusters.begin(); i != remainingClusters.end(); ++i) {
      
      
      // 12. build PFlow objects

      ReconstructedParticleImpl* recoParticle = assignNeutralClusterToReconstructedParticle(*i);

        
      // 13. perfom PID on them and add them to reco particle collection
      
      doPID(recoParticle,false);



      // check plots
      int nOfChargedHits = 0;
      double eOfChargedHits = 0.0;
      std::vector<CalorimeterHit*> chargedCaloHits = getChargedHitsAssignedToNeutralParticle(evt,recoParticle,nOfChargedHits,eOfChargedHits,0.5);
      
      #ifdef MARLIN_USE_AIDA
      _cNChargedHitsAssignedToNeutral->fill(nOfChargedHits);
      _cEChargedHitsAssignedToNeutral->fill(eOfChargedHits);
      #endif

      // debug
      // std::cout << "CHARGED HITS: n = " << nOfChargedHits << "  " << "E = " << eOfChargedHits << std::endl;
      
      for (std::vector<CalorimeterHit*>::const_iterator i = chargedCaloHits.begin(); i != chargedCaloHits.end(); ++i) {
	
	std::vector<CalorimeterHit*>::const_iterator position = find(chargedCaloHitsAssignedToNeutral.begin(),chargedCaloHitsAssignedToNeutral.end(),(*i));
	
	if ( position == chargedCaloHitsAssignedToNeutral.end() ) chargedCaloHitsAssignedToNeutral.push_back(*i);
	
      }
      chargedCaloHits.clear();
      
      if (_drawOnCED) {
	
	for (std::vector<CalorimeterHit*>::const_iterator i = chargedCaloHitsAssignedToNeutral.begin(); i != chargedCaloHitsAssignedToNeutral.end(); ++i) {
	  
	  ced_hit ( (*i)->getPosition()[0], (*i)->getPosition()[1], (*i)->getPosition()[2], 0 | 0 << CED_LAYER_SHIFT, 4, 0xe655ff );
	  
	}
	
      }
      
      
      
      // debug      
      /*
      if ( _debugLevel > 5 )MarlinUtil::printRecoParticle(recoParticle,_bField);
      
      if (_drawOnCED) { 

	MarlinCED::newEvent(this,0);
	//MarlinCED::drawMCParticleTree(evt,"MCParticle",0.05,4.0,15.5,50.0,1626.0,2500.0);       
	MarlinCED::drawRecoParticle(recoParticle,0,1,0xfff641,9 );
	MarlinCED::draw(this,true);

      }
      */
      



      // fill cluster collection     
      for(ClusterVec::const_iterator j = recoParticle->getClusters().begin(); j != recoParticle->getClusters().end(); ++j) {

	reconstructedClusters->addElement(*j);

      }
      

      // 14. add recoParticle to collection
		
      reconstructedParticles->addElement(recoParticle);
      
      
      
    }

    

    
    // debug 
    // calculate check numbers to fill histograms
    
    int numberOfReconstructedElectrons = 0;
    int numberOfReconstructedChargedHadrons = 0;
    int numberOfReconstructedPhotons = 0;
    int numberOfReconstructedNeutralHadrons = 0;
    int numberOfReconstructedMuons = 0;

    int numberOfChargedReconstructedParticles = 0;
    int numberOfNeutralReconstructedParticles = 0;
    int numberOfReconstructedNoTypeParticles = 0;  
    
    double SumEnergyOfReconstructedElectrons = 0.0;
    double SumEnergyOfReconstructedChargedHadrons = 0.0;
    double SumEnergyOfReconstructedPhotons = 0.0;
    double SumEnergyOfReconstructedNeutralHadrons = 0.0;
    double SumEnergyOfReconstructedMuons = 0.0;

    double SumEnergyOfChargedReconstructedParticles = 0.0;
    double SumEnergyOfNeutralReconstructedParticles = 0.0;
    double SumEnergyOfReconstructedNoTypeParticles = 0.0;  
    


    /*if ( _debugLevel > 1 )*/ std::cout << "collection of reconstructed particles (Size: " << reconstructedParticles->getNumberOfElements() << "):" << std::endl;
    for (LCCollectionVec::const_iterator i = reconstructedParticles->begin(); i != reconstructedParticles->end(); ++i) {

      ReconstructedParticle* recoParticle = dynamic_cast<ReconstructedParticle*>( (*i) );
      SEReco  += recoParticle->getEnergy();
      SpxReco += recoParticle->getMomentum()[0];
      SpyReco += recoParticle->getMomentum()[1];
      SpzReco += recoParticle->getMomentum()[2];

      // debug
      if ( _debugLevel > 5 ) MarlinUtil::printRecoParticle(recoParticle,_bField);

      
      switch ( recoParticle->getType() ) {
	
      case 1 : 
	++numberOfReconstructedElectrons;
	SumEnergyOfReconstructedElectrons += recoParticle->getEnergy();
	
	++numberOfChargedReconstructedParticles;
	SumEnergyOfChargedReconstructedParticles += recoParticle->getEnergy();
	break;
	
      case 2 : 
	++numberOfReconstructedChargedHadrons;
	SumEnergyOfReconstructedChargedHadrons += recoParticle->getEnergy();
	
	++numberOfChargedReconstructedParticles;
	SumEnergyOfChargedReconstructedParticles += recoParticle->getEnergy();
	break;
	
      case 3 : 
	++numberOfReconstructedPhotons;
	SumEnergyOfReconstructedPhotons += recoParticle->getEnergy();
	
	++numberOfNeutralReconstructedParticles;
	SumEnergyOfNeutralReconstructedParticles += recoParticle->getEnergy();
	break;

      case 4 : 
	++numberOfReconstructedNeutralHadrons;
	SumEnergyOfReconstructedNeutralHadrons += recoParticle->getEnergy();
	
	++numberOfNeutralReconstructedParticles;
	SumEnergyOfNeutralReconstructedParticles += recoParticle->getEnergy();
	break;

      case 5 : 
	++numberOfReconstructedMuons;
	SumEnergyOfReconstructedMuons += recoParticle->getEnergy();
	
	++numberOfChargedReconstructedParticles;
	SumEnergyOfChargedReconstructedParticles += recoParticle->getEnergy();
	break;

      default :	
	++numberOfReconstructedNoTypeParticles;
	SumEnergyOfReconstructedNoTypeParticles += recoParticle->getEnergy();

      }
          
    }
    
  
    // debug
    // summary of reconstructed particles
    //    if ( _debugLevel > 1 ) {
      
      std::cout << "N of reconstructed electrons " << numberOfReconstructedElectrons << "  " << "with energy " << SumEnergyOfReconstructedElectrons << std::endl
		<< "N of reconstructed charged hadrons " << numberOfReconstructedChargedHadrons << "  " << "with energy " << SumEnergyOfReconstructedChargedHadrons << std::endl
		<< "N of reconstructed photons " << numberOfReconstructedPhotons << "  " << "with energy " << SumEnergyOfReconstructedPhotons << std::endl
		<< "N of reconstructed neutral hadrons " << numberOfReconstructedNeutralHadrons << "  " << "with energy " << SumEnergyOfReconstructedNeutralHadrons << std::endl
    		<< "N of reconstructed muons " << numberOfReconstructedMuons << "  " << "with energy " << SumEnergyOfReconstructedMuons << std::endl
		<< std::endl
		<< "N of charged reconstructed particles " << numberOfChargedReconstructedParticles << "  " << "with energy " 
		<< SumEnergyOfChargedReconstructedParticles << std::endl
		<< "N of neutral reconstructed particles " << numberOfNeutralReconstructedParticles << "  " << "with energy " 
		<< SumEnergyOfNeutralReconstructedParticles << std::endl
		<< "N of reconstructed particles w/o type " <<  numberOfReconstructedNoTypeParticles << "  " << "with energy " 
		<< SumEnergyOfReconstructedNoTypeParticles << std::endl
		<< std::endl;
      
      //    }
    

    evt->addCollection(reconstructedClusters,_reconstructedClusterCollectionName.c_str());
    evt->addCollection(reconstructedParticles,_reconstructedParticleCollectionName.c_str());



  }
  
  catch(DataNotAvailableException &e){std::cout << "no valid collection in event " << _nEvt << std::endl; };



  // debug  
  /*
  std::cout << "# REAL MIP STUB: " << _nOfRealMIPStubs << "  " << "# FOUND MIP STUB: " << _nOfFoundMIPStubs << "  " << std::endl
	    << "# CHARGED OBJECTS EXTRAPOLATED INTO CALORIMETER in this event: " << nOfChargedObjectsExtrapolatedIntoCalorimeterPerEvent << "  " 
	    << "# CHARGED OBJECTS EXTRAPOLATED INTO CALORIMETER: " << _nOfChargedObjectsExtrapolatedIntoCalorimeter << std::endl;	
  */



  // fill check plots
  #ifdef MARLIN_USE_AIDA
  double SECalo = MarlinUtil::getEnergyDepositedInFullCalorimeter(evt);
  _cSECalo->fill(SECalo);
  _cSEReco->fill(SEReco);

  double Minv = sqrt( pow(SEReco,2) - pow(SpxReco,2) - pow(SpyReco,2) - pow(SpzReco,2) );
  _cInvMass->fill(Minv);

  double accumulatedEnergies[21] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  MarlinUtil::getMC_Balance(evt,accumulatedEnergies);
  double SEMC = accumulatedEnergies[20];
  _cSEMC->fill(SEMC);
  
  _cDSERecoMC->fill(SEReco - SEMC);
  _cDSERecoCalo->fill(SEReco - SECalo);
  _cDSECaloMC->fill(SECalo - SEMC);
  
  double alphaRecoMC = ((SEReco - SEMC)/sqrt(SEMC))*100.0;  
  double alphaRecoCalo = ((SEReco - SECalo)/sqrt(SECalo))*100.0;
  double alphaCaloMC = ((SECalo - SEMC)/sqrt(SEMC))*100.0;

  _cAlphaRecoMC->fill(alphaRecoMC);
  _cAlphaRecoCalo->fill(alphaRecoCalo);
  _cAlphaCaloMC->fill(alphaCaloMC);

  _cNTracksWithEnergyAssignedByMC->fill(NTracksWithEnergyAssignedByMC);
  _cNTracksNotPassingCylinderCuts->fill(NTracksNotPassingCylinderCuts);
  _cSumPTracksNotPassingCylinderCuts->fill(SumPTracksNotPassingCylinderCuts);
  _cSumPTracksNotPassingCylinderCutsvsDSERecoMC->fill(SumPTracksNotPassingCylinderCuts,SEReco - SEMC);

  _cNExtrapolatedObjects->fill(nOfChargedObjectsExtrapolatedIntoCalorimeterPerEvent);


  double sumENeutralHitsAssignedToCharged = 0.0;
  for (std::vector<CalorimeterHit*>::const_iterator i = neutralCaloHitsAssignedToCharged.begin(); i != neutralCaloHitsAssignedToCharged.end(); ++i) {

    sumENeutralHitsAssignedToCharged += (*i)->getEnergy();

  }

  _cSumNNeutralHitsAssignedToCharged->fill(neutralCaloHitsAssignedToCharged.size());
  _cSumENeutralHitsAssignedToCharged->fill(sumENeutralHitsAssignedToCharged);

  _cSumENeutralHitsAssignedToChargedvsDSERecoMC->fill(sumENeutralHitsAssignedToCharged,SEReco - SEMC);



  double sumEChargedHitsAssignedToNeutrals = 0.0;
  for (std::vector<CalorimeterHit*>::const_iterator i = chargedCaloHitsAssignedToNeutral.begin(); i != chargedCaloHitsAssignedToNeutral.end(); ++i) {

    sumEChargedHitsAssignedToNeutrals += (*i)->getEnergy();

  }

  _cSumNChargedHitsAssignedToNeutrals->fill(chargedCaloHitsAssignedToNeutral.size());
  _cSumEChargedHitsAssignedToNeutrals->fill(sumEChargedHitsAssignedToNeutrals);

  _cSumEChargedHitsAssignedToNeutralsvsDSERecoMC->fill(sumEChargedHitsAssignedToNeutrals,SEReco - SEMC);

  
  _cDEdcEwavsDSERecoMC->fill(sumEChargedHitsAssignedToNeutrals-sumENeutralHitsAssignedToCharged,SEReco - SEMC);
  _cDEdcEwaMinusDSERecoMCvsDSERecoMC->fill(((sumEChargedHitsAssignedToNeutrals-sumENeutralHitsAssignedToCharged)-(SEReco - SEMC)),SEReco - SEMC);


  if ( fabs(alphaRecoMC) >= _outputConditionLimit*100.0 ) setReturnValue(true);
  else setReturnValue(false);

  // debug
  /*
  std::cout << "fabs(alphaRecoMC) = " << fabs(alphaRecoMC) << "  " << "_outputConditionLimit*100.0 = " << _outputConditionLimit*100.0 << "  " 
	    <<  (fabs(alphaRecoMC) >= _outputConditionLimit*100.0)  << std::endl << std::endl;
  */


  /*
  _debugLevel = 7;
  _drawOnCED = 1;
  */

  double notAssignedCaloEnergy = getNotAssignedCalorimeterEnergy(evt,reconstructedParticles);

  // debug
  std::cout << "SpxReco = " << SpxReco << "  " << "SpyReco = " << SpyReco << "  " << "SpzReco = " << SpzReco << std::endl
	    << "SECalo = " << SECalo << "  " << "SEReco = " << SEReco << " ( Minv = " << Minv << " )  " << "SEMC = " << SEMC << std::endl
	    << "SEdc   = " << sumEChargedHitsAssignedToNeutrals << "  " << "SEwa = " << sumENeutralHitsAssignedToCharged << "  " << "SEnotAssigned = "
	    << notAssignedCaloEnergy << std::endl
	    << "SEReco-SEMC   = " << SEReco-SEMC << "  " << "aRecoMC   = " << alphaRecoMC << "%" << std::endl
	    << "SEReco-SECalo = " << SEReco-SECalo << "  " << "aRecoCalo = " << alphaRecoCalo << "%" << std::endl
	    << "SECalo-SEMC   = " << SECalo-SEMC << "  " << "aCaloMC   = " << alphaCaloMC << "%" << std::endl
	    << "DSEdcSEwa     = " << sumEChargedHitsAssignedToNeutrals-sumENeutralHitsAssignedToCharged << "  " 
	    << "(DSEdcSEwa)/(SEReco-SEMC) = " << (sumEChargedHitsAssignedToNeutrals-sumENeutralHitsAssignedToCharged)/(SEReco-SEMC) << std::endl;
  #endif


  if (_drawOnCED) MarlinCED::draw(this,true);

  /*
  _debugLevel = -7;
  _drawOnCED = 0;
  */

  _nEvt ++;
  firstEvent = false ;

}



void TrackBasedPFlow::check( LCEvent * evt )
{
  // nothing to check here
}



void TrackBasedPFlow::end()
{

  delete _mcParticleHelper;
  _mcParticleHelper = 0;

}



const TrackerHitVec TrackBasedPFlow::getOuterTrackerHits(const Track* track, unsigned int n) {


  unsigned int nHits = track->getTrackerHits().size();

  // returning an empty vector if less tracker hits than outer hits demanded
  if ( nHits < n ) {

    const TrackerHitVec emptyTrackerHitVec;
    return emptyTrackerHitVec;

  }
  else {

    LCVector3D refPoint(track->getReferencePoint()[0],track->getReferencePoint()[1],track->getReferencePoint()[2]);

    Trajectory* helix = new SimpleHelix(track->getD0(),track->getPhi(),track->getOmega(),track->getZ0(),track->getTanLambda(),refPoint);

    //double distanceRefPoint = getDistanceToHelix(refPoint,helix); // not used at the moment
    LCVector3D refPointProjected = getProjectedPointOnHelix(refPoint,helix);
    


    // debug
    /*
    std::cout << std::setprecision(6)
	      << "D0 = " <<  track->getD0() << "  " << "Z0 = " << track->getZ0() << "  " << "phi = " << track->getPhi() << "  " << "omega = " << track->getOmega() << "  " 
	      << "tanlambda = " << track->getTanLambda() << "  " << "refPoint = " << "(" << refPoint.x() << "," << refPoint.y() << "," << refPoint.z() << ")" << "  "
	      << "BField = " << _bField << std::endl;
    */


    // debug
    /*
    std::cout << std::setprecision(6)
	      << "nHits = " << nHits << "  "
	      << "refPoint = " << "(" << refPoint.x() << "," << refPoint.y() << "," << refPoint.z() << ")" << "  "
	      << "refPointProjected = " << "(" << refPointProjected.x() << "," << refPointProjected.y() << "," << refPointProjected.z() << ")" << "  "
	      << "d = " << distanceRefPoint << std::endl;
    */
    


    // debug
    /*   
    double d0    = (double)(track->getD0());
    double z0    = (double)(track->getZ0());
    double Phi0  = (double)(track->getPhi());
    double omega = (double)(track->getOmega());
    double tanL  = (double)(track->getTanLambda());

    double X0 = (1/omega -  d0)*sin(Phi0);
    double Y0 = (-1.0)*(1/omega -  d0)*cos(Phi0);
    double R0 = 1/fabs(omega);
    double bz = (-1.0)*omega/tanL;
    double phi0 = (z0*omega)/tanL + Phi0 + (omega*(acos(-1.0)))/(2.0*fabs(omega));

    std::cout << "High Precision:" << "  " << "DBL_EPSILON = " << DBL_EPSILON << "  " << "DBL_DIG = " << DBL_DIG << std::endl;
    std::cout << std::setprecision(DBL_DIG) 
	      << "d0    = " << d0 << std::endl
	      << "z0    = " << z0 << std::endl
	      << "Phi0  = " << Phi0 << std::endl
	      << "omega = " << omega << std::endl 
	      << "tanL  = " << tanL << std::endl
	      << std::endl  << std::endl
	      << "X0    = " << X0 << std::endl
	      << "Y0    = " << Y0 << std::endl
	      << "R0    = " << R0 << std::endl
	      << "bz    = " << bz << std::endl
	      << "phi0  = " << phi0 << std::endl
  	      << std::endl  << std::endl;    
    */    

    // debug
    // std::cout << std::setprecision(6) << "Track (X0,Y0,R0,bz,phi0) = " << "(" << X0 << "," << Y0 << "," << R0 << "," << bz << "," << phi0 << ")" << std::endl;  



    

    TrackerHitVec outermostTrackerHits;
    std::vector< std::pair<TrackerHit*,float> > trackerHitsWithDistance;
    
    for(unsigned int i=0; i<nHits; ++i){
      
      TrackerHit* trackerHit = track->getTrackerHits().at(i);

      LCVector3D trackerHitPoint(trackerHit->getPosition()[0],trackerHit->getPosition()[1],trackerHit->getPosition()[2]);
      LCVector3D trackerHitPointProjected = getProjectedPointOnHelix(trackerHitPoint,helix);

      // double distanceTrackerHitPoint = getDistanceToHelix(trackerHitPoint,helix); // not used at the moment
      double s = getPathLengthOnHelix(refPointProjected,trackerHitPointProjected,helix);


      
      
      // debug 

      if (i==0) {

	// double rmin = track->getRadiusOfInnermostHit();

	if (_drawOnCED) {
	  MarlinCED::drawTrajectory(helix,0 | ( 5 << CED_LAYER_SHIFT ),1,0x51fcff);
	
	  // or self made draw helix (see above)

	}
	/*
	std::cout << "trackerHitPoint = " << "(" << trackerHitPoint.x() << "," << trackerHitPoint.y() << "," << trackerHitPoint.z() << ")" << "  " 
		  << "trackerHitPointProjected = " << "(" << trackerHitPointProjected.x() << "," << trackerHitPointProjected.y() << "," << trackerHitPointProjected.z() << ")" 
		  << "  " << "d = " << distanceTrackerHitPoint << "  " << "s = " << s << std::endl;
	*/
      }      
      if (_drawOnCED) {
		      
	ced_hit ( trackerHitPoint.x(),trackerHitPoint.y(),trackerHitPoint.z(), 0 | 2 << CED_LAYER_SHIFT, 2, 0xff0004 );
	ced_hit ( trackerHitPointProjected.x(),trackerHitPointProjected.y(),trackerHitPointProjected.z(), 0 | 2 << CED_LAYER_SHIFT, 2, 0x17ff06 );
	ced_line( trackerHitPoint.x(),trackerHitPoint.y(),trackerHitPoint.z(),
		  trackerHitPointProjected.x(),trackerHitPointProjected.y(),trackerHitPointProjected.z(),0 | 2 << CED_LAYER_SHIFT, 1, 0xffffff );
	
      }
      
      
      std::pair<TrackerHit*,float> p(trackerHit,s);
      
      if ( i == 0) {
	
	trackerHitsWithDistance.push_back(p);
	
      }
      else {
	
	bool trackerHitInserted = false;

	for(std::vector< std::pair<TrackerHit*,float> >::iterator j = trackerHitsWithDistance.begin(); j != trackerHitsWithDistance.end(); ++j){
	  
	  if ( s > (*j).second) {
	    
	    trackerHitsWithDistance.insert(j,p);
	    trackerHitInserted = true;
	    break;
	    
	  }
	  
	}
	
	if (!trackerHitInserted) trackerHitsWithDistance.push_back(p);

      }
      
    }
    
    trackerHitsWithDistance.resize(n); // take only the largest n trackerHitsWithDistance
    
    // fill the output vector
    for(std::vector< std::pair<TrackerHit*,float> >::iterator j = trackerHitsWithDistance.begin(); j != trackerHitsWithDistance.end(); ++j){


      // debug
      /*      
      std::cout << "n outermost hits = " << trackerHitsWithDistance.size() << "  " << "hit = " << (*j).first << "  " 
		<< "projposref = " << "( " << refPointProjected.x() << ", " << refPointProjected.y() << ", " << refPointProjected.z() << " )" << "  "
		<< "pos = " << "( " << (*j).first->getPosition()[0] << ", " << (*j).first->getPosition()[1] << ", " << (*j).first->getPosition()[2] << " )" << "  "
		<< "path length = " << (*j).second << std::endl;
      */


      outermostTrackerHits.push_back( (*j).first );
      
    }
    
    delete helix;
    helix = 0;
    
    return outermostTrackerHits;
    
  }
  
}




const std::vector<CalorimeterHitWithAttributes*> TrackBasedPFlow::getRelatedCalorimeterHits(const LCEvent* evt, const TrackerHitVec outermostTrackerHits, 
											    Trajectory* fittedHelix,
											    int trackNumber/*only for debugging*/) {

  
  std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes;


  try {

    LCVector3D referencePosition(outermostTrackerHits.at(0)->getPosition()[0],outermostTrackerHits.at(0)->getPosition()[1],outermostTrackerHits.at(0)->getPosition()[2]);
    LCVector3D referencePositionProjected = getProjectedPointOnHelix(referencePosition,fittedHelix);

    double sOfReferencePositionProjected = getPathLengthOnHelix(referencePositionProjected,fittedHelix);

    SimpleHelix* helixTrajectory = dynamic_cast<SimpleHelix*>(fittedHelix);
    helixTrajectory->setStart(sOfReferencePositionProjected);
    helixTrajectory->setEnd(sOfReferencePositionProjected+_maximalConeTubeLength);


    //debug
    if ( _debugLevel > 7 ) {
      std::cout << "event " << _nEvt << "  " << "referencePosition = " << "( " << referencePosition.x() << ", " << referencePosition.y() << ", " 
		<< referencePosition.z() << " )" << std::endl;
    }


    /*
    #ifdef MARLIN_USE_AIDA
    static AIDA::ICloud2D* cDistanceVSPathlength;
    std::ostringstream ostr;
    ostr << "_" << _nEvt << "_" << trackNumber;
    std::string cloudName = "Track " + ostr.str();
    
    std::cout<< cloudName <<  std::endl;
    
    cDistanceVSPathlength = AIDAProcessor::histogramFactory(this)->createCloud2D( cloudName, "d = d(s) for this particular track", -1 );
    #endif
    */


    float energyInConelikeTube = 0.0;

    double meanEnergyECALSampling1 = 0.0;
    double meanEnergyECALSampling2 = 0.0;
    double meanEnergyHCALSampling1 = 0.0;

    int nHitsECALSampling1 = 0;
    int nHitsECALSampling2 = 0;
    int nHitsHCALSampling1 = 0;
    

    std::vector< std::string >::const_iterator iter;
    const std::vector< std::string >* ColNames = evt->getCollectionNames();
    
    for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
	
      LCCollection* col = evt->getCollection( *iter ) ;

      if ( (col->getTypeName() == LCIO::CALORIMETERHIT) && ( (*iter == _colNameECAL) || (*iter == _colNameHCAL) ) ) {

	int NCalorimeterHits = col->getNumberOfElements();

	
	// debug  
	if ( _debugLevel > 7 ) std::cout << "collection: " << *iter << "  " << "NCalorimeterHits in this collection = " << NCalorimeterHits << std::endl;


	// debug
	float Emax = 0.0; // get maximal deposited energy in whole Calorimeter for the rainbow scale
	if ( _drawOnCED ) {

	  for(int j=0; j<NCalorimeterHits; ++j) {
	  
	    CalorimeterHit* calorimeterHit = dynamic_cast<CalorimeterHit*>(col->getElementAt(j));
	    if ( (calorimeterHit->getEnergy()) > Emax ) Emax = calorimeterHit->getEnergy();
	  
	  }
	
	}	


	for(int j=0; j<NCalorimeterHits; ++j) {
	
	  CalorimeterHit* calorimeterHit = dynamic_cast<CalorimeterHit*>(col->getElementAt(j));

	
	  LCVector3D calorimeterHitPosition(calorimeterHit->getPosition()[0],calorimeterHit->getPosition()[1],calorimeterHit->getPosition()[2]);
	  LCVector3D calorimeterHitPositionProjected = getProjectedPointOnHelix(calorimeterHitPosition,fittedHelix);

	  double pathLengthOnHelix = getPathLengthOnHelix(referencePositionProjected,calorimeterHitPositionProjected,fittedHelix);
	  double distance = getDistanceToHelix(calorimeterHitPosition,fittedHelix);



	  // debug
	  if ( _debugLevel > 7 ) {
	    std::cout << "hit = " << calorimeterHit << "  " << "pos = " << "(" << calorimeterHitPosition.x() << "," << calorimeterHitPosition.y() << "," 
		      << calorimeterHitPosition.z() << ")" << "  " << "pos proj = " << "(" << calorimeterHitPositionProjected.x() << "," 
		      << calorimeterHitPositionProjected.y() << "," << calorimeterHitPositionProjected.z() << ")" << "  " << "ref pos on helix = " << "(" 
		      << referencePosition.x() << "," << referencePosition.y() << "," << referencePosition.z() << ")" << "  " 
		      << std::endl 
		      << "                   " << "s = " << pathLengthOnHelix << "  " << "d = " << distance << "  " << "alphaTube = " << _openingAngleConeTube << "  "
		      << "tan(alpha/2)*pathLengthOnHelix = " << tan(_openingAngleConeTube/2) * pathLengthOnHelix << "  " 
		      << "E = " << calorimeterHit->getEnergy() << std::endl;
	  }



	  // cut on cone-like tube
	  if ( ( distance <= fabs(tan(_openingAngleConeTube/2) * pathLengthOnHelix) ) && 
	       ( pathLengthOnHelix <= _maximalConeTubeLength ) && 
	       ( pathLengthOnHelix >= 0.0 ) ) {	    


  
	    // have a look if the CalorimeterHit is part of a EM ShowerCandidate
	    if ( calorimeterHit->ext<isPartOfEMShowerCandidate>() == 1 ) {
	      
	      // if the distance of the center of the EMShower Candidate to the trajectory extrapolation is small (electron candidate), continue
	      LCCollection* collectionOfEMShowerCandidates = evt->getCollection(_colNameEMShowerCandidates.c_str()) ;
	      
	      if( collectionOfEMShowerCandidates != 0 ) {
		
		// find cluster which the CalorimeterHit is assigned to
		Cluster* clusterOfCalorimeterHit = 0;
		bool clusterOfCalorimeterHitFound = false;
		
		for(int k=0; k<collectionOfEMShowerCandidates->getNumberOfElements(); ++k) {
		  
		  Cluster* cluster = dynamic_cast<Cluster*>(collectionOfEMShowerCandidates->getElementAt(k));
		  
		  CalorimeterHitVec::const_iterator position = find(cluster->getCalorimeterHits().begin(),cluster->getCalorimeterHits().end(),calorimeterHit);
		  
		  if ( position != cluster->getCalorimeterHits().end() ) {
		    
		    // each CalorimeterHit can belong to one cluster only
		    clusterOfCalorimeterHit = cluster;
		    clusterOfCalorimeterHitFound =true;
		    break;
		    
		  }
		  
		}
		
		LCVector3D clusterPosition(clusterOfCalorimeterHit->getPosition()[0],clusterOfCalorimeterHit->getPosition()[1],clusterOfCalorimeterHit->getPosition()[2]);
		
		double distanceOfClusterPositionToTrajectory = getDistanceToHelix(clusterPosition,fittedHelix);
		
		int typeOfClusterPosition = getTypeOfPositionOfCluster(clusterOfCalorimeterHit);
		double maximalDistanceToCompareClusterPositionWith = _maximalDistanceOfElectronShowerPosition.at(typeOfClusterPosition);
		
		// continue loop, i.e. do not take into account this CalorimeterHit in the following if distanceOfClusterPositionToTrajectory is larger than cut
		if ( distanceOfClusterPositionToTrajectory > maximalDistanceToCompareClusterPositionWith ) continue;
		
	      }
	      
	    }



	    energyInConelikeTube += calorimeterHit->getEnergy();

	    
	    // debug 
	    float E = calorimeterHit->getEnergy();	    
	    if ( _debugLevel > 7 ) {
	      std::cout << "HIT PASSES CUT" << std::endl
			<< "hit = " << calorimeterHit << "  " << "pos = " << "(" << calorimeterHitPosition.x() << "," << calorimeterHitPosition.y() << "," 
			<< calorimeterHitPosition.z() << ")" << "  " << "pos proj = " << "(" << calorimeterHitPositionProjected.x() << "," 
			<< calorimeterHitPositionProjected.y() << "," << calorimeterHitPositionProjected.z() << ")" << "  " << "ref pos on helix = " << "(" 
			<< referencePosition.x() << "," << referencePosition.y() << "," << referencePosition.z() << ")" << "  " 
			<< std::endl 
			<< "                   " << "s = " << pathLengthOnHelix << "  " << "d = " << distance << "  " 
			<< "tan(alpha/2)*pathLengthOnHelix = " << tan(_openingAngleConeTube/2) * pathLengthOnHelix << "  " 
			<< "E = " << E << std::endl << std::endl;
	    }    

  
	    if ( _drawOnCED ) {
	      
	      ced_hit ( calorimeterHitPosition.x(),calorimeterHitPosition.y(),calorimeterHitPosition.z(), 0 | 10 << CED_LAYER_SHIFT, 4, 0xff0004 );
	      ced_hit ( calorimeterHitPositionProjected.x(),calorimeterHitPositionProjected.y(),calorimeterHitPositionProjected.z(), 0 | 10 << CED_LAYER_SHIFT, 4,
			0x17ff06 );
	      ced_line( calorimeterHitPosition.x(),calorimeterHitPosition.y(),calorimeterHitPosition.z(),
			calorimeterHitPositionProjected.x(),calorimeterHitPositionProjected.y(),calorimeterHitPositionProjected.z(),0 | 10 << CED_LAYER_SHIFT, 1,
			0xffffff );

	      
	      ced_hit ( referencePosition.x(),referencePosition.y(),referencePosition.z(), 0 | 4 << CED_LAYER_SHIFT, 6, 0x08f51f );
	      unsigned int color = MarlinDrawUtil::getColorAmplitude(E,Emax,"rainbow",0.4);

	      // draw all calorimeter hits in cone-like tube in rainbow color
	      ced_hit ( calorimeterHitPosition.x(),calorimeterHitPosition.y(),calorimeterHitPosition.z(), 0 | 7 << CED_LAYER_SHIFT, 2, color );

	      // draw all calorimeter hits in cone-like tube in one color
	      // ced_hit ( calorimeterHitPosition.x(),calorimeterHitPosition.y(),calorimeterHitPosition.z(), 2 | 8 << CED_LAYER_SHIFT, 8, 0xffff00 );

	      
	      double hitEnergyInMIPs = 0.0;
		
	      switch ( calorimeterHit->getType() ) {
		
	      case 0 : 
		hitEnergyInMIPs = (calorimeterHit->getEnergy())/_mipCoeffEcal.at(0);
		meanEnergyECALSampling1 += (calorimeterHit->getEnergy());
		++nHitsECALSampling1;
		break;
	      case 1 : 
		hitEnergyInMIPs = (calorimeterHit->getEnergy())/_mipCoeffEcal.at(1);
		meanEnergyECALSampling2 += (calorimeterHit->getEnergy());
		++nHitsECALSampling2;
		break;
	      case 2 : 
		hitEnergyInMIPs = (calorimeterHit->getEnergy())/_mipCoeffHcal.at(0);
		meanEnergyHCALSampling1 += (calorimeterHit->getEnergy());
		++nHitsHCALSampling1;
		break;
	      
	      }

	      
	      // std::cout << calorimeterHit->getEnergy() << "  " << hitEnergyInMIPs << "  " << calorimeterHit->getType() << std::endl;
	      
    
	      if (hitEnergyInMIPs < 0.5) { 

		ced_hit( calorimeterHitPosition.x(),calorimeterHitPosition.y(),calorimeterHitPosition.z(), 0 | 2 << CED_LAYER_SHIFT, 2, 0xff696c );

	      }
	      else if ( (hitEnergyInMIPs >= 0.5) && (hitEnergyInMIPs < 1.7) ) {

		ced_hit( calorimeterHitPosition.x(),calorimeterHitPosition.y(),calorimeterHitPosition.z(), 0 | 2 << CED_LAYER_SHIFT, 2, 0xff0000 );

	      }
	      else if ( (hitEnergyInMIPs >= 1.7) && (hitEnergyInMIPs < 3.5) ) {

		ced_hit( calorimeterHitPosition.x(),calorimeterHitPosition.y(),calorimeterHitPosition.z(), 0 | 3 << CED_LAYER_SHIFT, 2, 0x0dff00 );

	      }
	      else {

		ced_hit( calorimeterHitPosition.x(),calorimeterHitPosition.y(),calorimeterHitPosition.z(), 0 | 4 << CED_LAYER_SHIFT, 2, 0x2f6dff );
		
	      }

	    }

	    // debug
	    /*
            #ifdef MARLIN_USE_AIDA
	    cDistanceVSPathlength->fill(pathLengthOnHelix,distance);
	    #endif
	    */
	    // end debug




	    CalorimeterHitWithAttributes* calorimeterHitWithAttributes = new CalorimeterHitWithAttributes(calorimeterHit,distance,pathLengthOnHelix);

	    bool calorimeterHitInserted = false;
	      
	    for(std::vector<CalorimeterHitWithAttributes*>::iterator l = calorimeterHitsWithAttributes.begin(); l !=  calorimeterHitsWithAttributes.end(); ++l){
		
	      if ( pathLengthOnHelix > (*l)->getPathLengthOnHelix() ) {
		
		calorimeterHitsWithAttributes.insert(l,calorimeterHitWithAttributes);
		calorimeterHitInserted = true;
		break;
		  
	      }
		
	    }
	      
	    if (!calorimeterHitInserted) calorimeterHitsWithAttributes.push_back(calorimeterHitWithAttributes);
	    	    
	  }
	  
	}
	
	/*
	std::cout << "meanEnergyECALSampling1 = " << meanEnergyECALSampling1/nHitsECALSampling1 << "  " 
		  << "meanEnergyECALSampling2 = " << meanEnergyECALSampling2/nHitsECALSampling2 << "  " 
		  << "meanEnergyHCALSampling1 = " << meanEnergyHCALSampling1/nHitsHCALSampling1 << std::endl;
	*/


      }
      
    }

    // hits with the smallest path length on helix first    
    reverse(calorimeterHitsWithAttributes.begin(), calorimeterHitsWithAttributes.end());


    // debug 
    if ( _debugLevel > 7 ) {
      std::cout << "collected Calorimeter Hits sorted by their path length on the helix: " << std::endl;
      for (std::vector< CalorimeterHitWithAttributes* >::const_iterator i = calorimeterHitsWithAttributes.begin(); i != calorimeterHitsWithAttributes.end(); ++i) {

	std::cout << "Calorimeter Hit: " << (*i)->getCalorimeterHit() << "  " << "( " << (*i)->getCalorimeterHit()->getPosition()[0] << "," 
		  << (*i)->getCalorimeterHit()->getPosition()[1] << "," << (*i)->getCalorimeterHit()->getPosition()[2] << " )" << "  " 
		  << "l on helix = " << (*i)->getPathLengthOnHelix() << "  " << "d to helix = " << (*i)->getDistanceToHelix() << std::endl;
	
      }
    }

    // debug 
    if ( _debugLevel > 1 ) std::cout << "energy in cone-like tube = " << energyInConelikeTube << "  " << "n of hits in conelike tube = " 
				     << calorimeterHitsWithAttributes.size() << std::endl;
    


    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // needs to be done, but somewhwere else
    // delete calorimeterHitsWithAttributes;

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

  }
  
  catch(DataNotAvailableException &e){std::cout << "no valid collection in event " << _nEvt << std::endl; };

  return calorimeterHitsWithAttributes;

}




void TrackBasedPFlow::getRelatedClusters(const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, const TrackerHitVec outermostTrackerHits,
					 Trajectory* fittedHelix, ClusterImplWithAttributes* mipStub, std::vector<ClusterImplWithAttributes*>& clusters) {



  // set path start and end to the MIP tube limits
  LCVector3D referencePosition(outermostTrackerHits.at(0)->getPosition()[0],outermostTrackerHits.at(0)->getPosition()[1],outermostTrackerHits.at(0)->getPosition()[2]);
  LCVector3D referencePositionProjected = getProjectedPointOnHelix(referencePosition,fittedHelix);
  
  double sOfReferencePositionProjected = getPathLengthOnHelix(referencePositionProjected,fittedHelix);

  double relevantPathLength = _maximalPathLengthForMIPLikeStub;
  if ( _maximalPathLengthForMIPLikeStub > _maximalConeTubeLength ) relevantPathLength = _maximalConeTubeLength;

  SimpleHelix* helixTrajectory = dynamic_cast<SimpleHelix*>(fittedHelix);
  helixTrajectory->setStart(sOfReferencePositionProjected);
  helixTrajectory->setEnd(sOfReferencePositionProjected+relevantPathLength);


  // find mip stub
  getMIPStub(mipStub,calorimeterHitsWithAttributes,outermostTrackerHits,fittedHelix);


  // reset path start and end to the cone-like tube limits
  helixTrajectory->setStart(sOfReferencePositionProjected); // not necessary, just for safety
  helixTrajectory->setEnd(sOfReferencePositionProjected+_maximalConeTubeLength);




  // debug
  if ( _debugLevel > 1 ) {
    float mipStubEnergy = mipStub->getEnergy();
    std::cout << "MIP stub " << mipStub->getClusterImpl() << " with energy E = " << mipStubEnergy << "  " << "n of hits = " 
	      << mipStub->getClusterImpl()->getCalorimeterHits().size() << " and positions: " << std::endl
	      << "xStart = " << "( " << mipStub->getPositionStartHit()[0] << ", " << mipStub->getPositionStartHit()[1] << ", " << mipStub->getPositionStartHit()[2] << " )" 
	      << " " << "type = " << mipStub->getTypeStartHit() << "  "
	      << "xPos = "<< "( " << mipStub->getPosition()[0] << ", " << mipStub->getPosition()[1] << ", " << mipStub->getPosition()[2] << " )" << "  "
	      << "xEnd = " << "( " << mipStub->getPositionEndHit()[0] << ", " << mipStub->getPositionEndHit()[1] << ", " << mipStub->getPositionEndHit()[2] << " )" 
	      << " " << "type = " << mipStub->getTypeEndHit() 
	      << std::endl << std::endl;
  
  /*
  for (CalorimeterHitVec::const_iterator i = mipStub->getClusterImpl()->getCalorimeterHits().begin(); i != mipStub->getClusterImpl()->getCalorimeterHits().end(); ++i ) {
  
    std::cout << "Calorimeter Hit: " << (*i) << "  " << std::endl;
    
  }
  */
  }

  if (_drawOnCED) {
    MarlinCED::drawClusterImpl(mipStub->getClusterImpl(),0,1,0x08fff7,0); 
    ced_hit ( mipStub->getPositionStartHit()[0], mipStub->getPositionStartHit()[1], mipStub->getPositionStartHit()[2], 0 | 0 << CED_LAYER_SHIFT, 6, 0x08fff7 );
    ced_hit ( mipStub->getPosition()[0], mipStub->getPosition()[1], mipStub->getPosition()[2], 0 | 0 << CED_LAYER_SHIFT, 4, 0x08fff7 );
    ced_hit ( mipStub->getPosition()[0], mipStub->getPosition()[1], mipStub->getPosition()[2], 2 | 0 << CED_LAYER_SHIFT, 12, 0x08fff7 );
    ced_hit ( mipStub->getPositionEndHit()[0], mipStub->getPositionEndHit()[1], mipStub->getPositionEndHit()[2], 0 | 0 << CED_LAYER_SHIFT, 6, 0x08fff7 );

  }



  
  if (!mipStub->isMuon()) {

    // perform clustering on the remaining hits

    
    // 1. method: cone clustering
    /*
    if ( mipStub->isMIPStub() ) {

      std::vector<CalorimeterHitWithAttributes*> remainingCalorimeterHitsWithAttributes = removeMIPStub(calorimeterHitsWithAttributes,mipStub);

      clusters = doConeClustering(outermostTrackerHits,fittedHelix,remainingCalorimeterHitsWithAttributes, mipStub->getPositionEndHit(),
				  mipStub->getDirectionEndHit());

    }
    else {

      clusters = doConeClustering(outermostTrackerHits,fittedHelix,calorimeterHitsWithAttributes,mipStub->getPositionEndHit(),
				  mipStub->getDirectionEndHit());

    }
    */

        
    // 2. method: apply trackwise clustering

    LCVector3D startPosition(mipStub->getPositionEndHit()[0],mipStub->getPositionEndHit()[1],mipStub->getPositionEndHit()[2]);
    LCVector3D startPositionProjected = getProjectedPointOnHelix(startPosition,fittedHelix);

    double pathLengthOnHelixOfStartPositionProjected  = getPathLengthOnHelix(referencePositionProjected,startPositionProjected,fittedHelix);
    double distanceToHelixOfStartPosition = getDistanceToHelix(startPosition,fittedHelix);
            
    if ( mipStub->isMIPStub() ) {

      std::vector<CalorimeterHitWithAttributes*> remainingCalorimeterHitsWithAttributes = removeMIPStub(calorimeterHitsWithAttributes,mipStub);
      	
      clusters = doTrackwiseClustering(outermostTrackerHits,fittedHelix,remainingCalorimeterHitsWithAttributes, mipStub->getPositionEndHit(),
				       pathLengthOnHelixOfStartPositionProjected, distanceToHelixOfStartPosition, mipStub->getDirectionEndHit());

    }
    else {
      
      clusters = doTrackwiseClustering(outermostTrackerHits,fittedHelix,calorimeterHitsWithAttributes,mipStub->getPositionEndHit(),
				       pathLengthOnHelixOfStartPositionProjected, distanceToHelixOfStartPosition, mipStub->getDirectionEndHit());
      
    }
    
    
    
    // 3. method: perform NN clustering
    /*
    if ( mipStub->isMIPStub() ) {
      
      std::vector<CalorimeterHitWithAttributes*> remainingCalorimeterHitsWithAttributes = removeMIPStub(calorimeterHitsWithAttributes,mipStub);
      
      clusters = doNNClustering(outermostTrackerHits,fittedHelix,remainingCalorimeterHitsWithAttributes, mipStub->getPositionEndHit(),
				mipStub->getDirectionEndHit());
      
    }
    else {
      
      clusters = doNNClustering(outermostTrackerHits,fittedHelix,calorimeterHitsWithAttributes,mipStub->getPositionEndHit(),
				mipStub->getDirectionEndHit());
      
    }
    */


    // 4. method: perform track based clustering
    /*
    if ( mipStub->isMIPStub() ) {
      
      std::vector<CalorimeterHitWithAttributes*> remainingCalorimeterHitsWithAttributes = removeMIPStub(calorimeterHitsWithAttributes,mipStub);
      
      clusters = doTrackBasedClustering(outermostTrackerHits,fittedHelix,remainingCalorimeterHitsWithAttributes, mipStub->getPositionEndHit(),
					mipStub->getDirectionEndHit());
      
    }
    else {
      
      clusters = doTrackBasedClustering(outermostTrackerHits,fittedHelix,calorimeterHitsWithAttributes,mipStub->getPositionEndHit(),
					mipStub->getDirectionEndHit());
      
    }
    */


 
    // debug
    
    float clustersEnergy = 0.0;
    /*
    std::cout << std::endl << "all remaining clusters, start position = ( " << mipStub->getPositionEndHit()[0] << "," << mipStub->getPositionEndHit()[1] << "," 
	      << mipStub->getPositionEndHit()[2] << " ), start direction = ( " << mipStub->getDirectionEndHit()[0] << "," << mipStub->getDirectionEndHit()[1] << ","
	      << mipStub->getDirectionEndHit()[2] << " ):" << std::endl;
    */

    int iCluster = 0; 
    for (std::vector<ClusterImplWithAttributes*>::const_iterator i = clusters.begin(); i != clusters.end(); ++i) {
      clustersEnergy += (*i)->getEnergy();
      if (_drawOnCED) {
	int color = 0x8f57ff*(iCluster+32);

	if (iCluster < 8) { // draw all remaining clusters on layer 19 (hard coded) if number of clusters is bigger than 8
	  MarlinCED::drawClusterImpl((*i)->getClusterImpl(),0,1,color,11+iCluster);
	  ced_hit( (*i)->getPositionStartHit()[0], (*i)->getPositionStartHit()[1], (*i)->getPositionStartHit()[2], 0 | (11+iCluster) << CED_LAYER_SHIFT, 6, color );
	  ced_hit( (*i)->getPosition()[0], (*i)->getPosition()[1], (*i)->getPosition()[2], 0 | (11+iCluster) << CED_LAYER_SHIFT, 4, color );
	  ced_hit( (*i)->getPosition()[0], (*i)->getPosition()[1], (*i)->getPosition()[2], 2 | (11+iCluster) << CED_LAYER_SHIFT, 12, color );
	  ced_hit( (*i)->getPositionEndHit()[0], (*i)->getPositionEndHit()[1], (*i)->getPositionEndHit()[2], 0 | (11+iCluster) << CED_LAYER_SHIFT, 6, color );
	}
	else {
	  MarlinCED::drawClusterImpl((*i)->getClusterImpl(),0,1,color,19);
	  ced_hit( (*i)->getPositionStartHit()[0], (*i)->getPositionStartHit()[1], (*i)->getPositionStartHit()[2], 0 | 19 << CED_LAYER_SHIFT, 6, color );
	  ced_hit( (*i)->getPosition()[0], (*i)->getPosition()[1], (*i)->getPosition()[2], 0 | 19 << CED_LAYER_SHIFT, 4, color );
	  ced_hit( (*i)->getPosition()[0], (*i)->getPosition()[1], (*i)->getPosition()[2], 2 | 19 << CED_LAYER_SHIFT, 12, color );
	  ced_hit( (*i)->getPositionEndHit()[0], (*i)->getPositionEndHit()[1], (*i)->getPositionEndHit()[2], 0 | 19 << CED_LAYER_SHIFT, 6, color );
	}

      }


      if ( _debugLevel > 2 ) {
      
	float distanceToEndPointOfMIP = sqrt( pow( ( (*i)->getPositionStartHit()[0] - mipStub->getPositionEndHit()[0] ) ,2) + 
					      pow( ( (*i)->getPositionStartHit()[1] - mipStub->getPositionEndHit()[1] ) ,2) + 
					      pow( ( (*i)->getPositionStartHit()[2] - mipStub->getPositionEndHit()[2] ) ,2) );
	
	std::cout << "cluster " << iCluster << " with energy E = " << (*i)->getEnergy() << "  " << "n of hits = " 
		  << (*i)->getClusterImpl()->getCalorimeterHits().size() << " and positions: " << std::endl 
		  << "xStart = " << "( " << (*i)->getPositionStartHit()[0] << ", " << (*i)->getPositionStartHit()[1] << ", " << (*i)->getPositionStartHit()[2] << " )" 
		  << " " << "type = " << (*i)->getTypeStartHit() << "  "
		  << "xPos = "<< "( " << (*i)->getPosition()[0] << ", " << (*i)->getPosition()[1] << ", " << (*i)->getPosition()[2] << " )" << "  "
		  << "xEnd = " << "( " << (*i)->getPositionEndHit()[0] << ", " << (*i)->getPositionEndHit()[1] << ", " << (*i)->getPositionEndHit()[2] << " )" 
		  << " " << "type = " << (*i)->getTypeEndHit() << "  "
		  << std::endl
		  << "3d distance to endpoint of MIP = " << distanceToEndPointOfMIP 	
		  << std::endl << std::endl;

      }
	
      ++iCluster;

    }
    
    // std::cout << "sum of energy in all clusters: " << mipStubEnergy+clustersEnergy << std::endl;




    
    
  }
  
}


								      
void TrackBasedPFlow::getRelatedClusterPerfectly(const LCEvent* evt, Track* track, const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix,
						 std::vector<ClusterImplWithAttributes*>& clusters) {

  std::vector<ClusterImpl*> resultingClusterImpls;

  ClusterImpl* clusterRealEnergy = new ClusterImpl();
  ClusterImpl* clusterPerfectEnergy = new ClusterImpl();

  getRelatedCalorimeterHitsPerfectly(evt,track,clusterRealEnergy,clusterPerfectEnergy);


  // assignClusterProperties(clusterRealEnergy);
  assignClusterProperties(clusterPerfectEnergy);

  //  if ( clusterRealEnergy->getCalorimeterHits().size() > 0 ) resultingClusterImpls.push_back(clusterRealEnergy);
  if ( clusterPerfectEnergy->getCalorimeterHits().size() > 0 ) resultingClusterImpls.push_back(clusterPerfectEnergy);



  // FIXME build function/class for all stuff below  
  // assign properties of the attributes of the cluster

  for (std::vector<ClusterImpl*>::const_iterator i = resultingClusterImpls.begin(); i != resultingClusterImpls.end(); ++i) {
    
    ClusterImplWithAttributes* clusterImplWithAttributes = new ClusterImplWithAttributes();
    
    clusterImplWithAttributes->setClusterImpl(*i);

    clusters.push_back(clusterImplWithAttributes);
    
  }



  // fill attributes of clusters
  double startPoint[3] = {0.0,0.0,0.0};
  assignClusterAttributes(outermostTrackerHits,fittedHelix,clusters,startPoint);

  return;
  
}



//  _________________________________________________________________________________________________________________________________________________________________________


// build class for MIP like stub etc.

void TrackBasedPFlow::getMIPStub(ClusterImplWithAttributes* clusterWithAttributes, const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, 
				 const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix) {


  double energy = 0.0;

  ClusterImpl* cluster = new ClusterImpl();
  clusterWithAttributes->setClusterImpl(cluster);

  
  std::vector<CalorimeterHitWithAttributes*>::const_iterator i = calorimeterHitsWithAttributes.begin();

  // no calorimeter hits in conlike tube or the the closest by calorimeter hit is already beyond 'maximalPathLengthForMIPLikeStub' -> take outermost tracker projection hit 
  // as a reference point
  if ( ( calorimeterHitsWithAttributes.size() == 0 ) || ( (*i)->getPathLengthOnHelix() > _maximalPathLengthForMIPLikeStub ) ) {
    
    LCVector3D outermostTrackerHit(outermostTrackerHits.at(0)->getPosition()[0],outermostTrackerHits.at(0)->getPosition()[1],outermostTrackerHits.at(0)->getPosition()[2]);
    LCVector3D outermostTrackerHitProjected = getProjectedPointOnHelix(outermostTrackerHit,fittedHelix);
    double s = getPathLengthOnHelix(outermostTrackerHitProjected,fittedHelix);
    LCVector3D tangent = fittedHelix->getDirection(s);

    clusterWithAttributes->setPositionEndHitLCVec(outermostTrackerHitProjected);
    clusterWithAttributes->setPositionStartHitLCVec(outermostTrackerHitProjected);
    
    clusterWithAttributes->setDirectionEndHitLCVec(tangent);
    clusterWithAttributes->setDirectionStartHitLCVec(tangent);
    
    clusterWithAttributes->setPositionLCVec(outermostTrackerHitProjected);

    clusterWithAttributes->setTypeEndHit(-1);
    clusterWithAttributes->setTypeStartHit(-1);   

    clusterWithAttributes->setIsMIPStub(false);
    clusterWithAttributes->setIsMuon(false);

    clusterWithAttributes->setEnergy(energy);



    //debug
    
    if (_drawOnCED){
      ced_hit ( outermostTrackerHitProjected.x(),outermostTrackerHitProjected.y(),outermostTrackerHitProjected.z(), 0 | 10 << CED_LAYER_SHIFT, 3, 0xff26c9);        
    }

    if ( _debugLevel > 5 ) {
      std::cout << "1. case: ( calorimeterHitsWithAttributes.size() == 0 ) || ( (*i)->getPathLengthOnHelix() > _maximalPathLengthForMIPLikeStub ) :" << std::endl
		<< "pos outermost tracker hit = " << "( " << outermostTrackerHit.x() << "," << outermostTrackerHit.y()  << "," << outermostTrackerHit.z() << " )" << "  " 
		<< "pos hit projected = " << "( " << outermostTrackerHitProjected.x() << "," << outermostTrackerHitProjected.y()  << "," 
		<< outermostTrackerHitProjected.z() << " )" << std::endl 
		<< "type start hit = " << clusterWithAttributes->getTypeStartHit() << "  " << "type end hit = " << clusterWithAttributes->getTypeEndHit()
		<< std::endl << std::endl;
    }


    return;
    
  }
  

  CalorimeterHitWithAttributes* storeLastHit = (*i);
  CalorimeterHitWithAttributes* storeFirstHit = (*i); 
  bool firstHitFound = false;
  int numberOfHitsInVetoTube = 0;
  double maximalRadiusOfInnerTubeForMIPLikeStubCompare = 0.0;
  double minimalRadiusOfOuterTubeForMIPLikeStubCompare = 0.0;


  // debug
  int index = 0;
  
  for (i = calorimeterHitsWithAttributes.begin(); i!=calorimeterHitsWithAttributes.end(); ++i) {


    /*
    double hitEnergyInMIPs = 0.0;
		
    switch ( (*i)->getCalorimeterHit()->getType() ) {
		
    case 0 : 
      hitEnergyInMIPs = ((*i)->getCalorimeterHit()->getEnergy())/_mipCoeffEcal.at(0);
      break;
    case 1 : 
      hitEnergyInMIPs = ((*i)->getCalorimeterHit()->getEnergy())/_mipCoeffEcal.at(1);
      break;
    case 2 : 
      hitEnergyInMIPs = ((*i)->getCalorimeterHit()->getEnergy())/_mipCoeffHcal.at(0);
      break;
      
    }
    */



    //debug: draws the sequence of hits collected for the MIP stub
    /*
    if (_drawOnCED){
      ced_hit ( (*i)->getCalorimeterHit()->getPosition()[0],(*i)->getCalorimeterHit()->getPosition()[1],(*i)->getCalorimeterHit()->getPosition()[2], 
		0 | 10 << CED_LAYER_SHIFT, 3, 0xff26c9);
    }
    */
    /*
    std::cout << index << "  " << "hit = " << (*i)->getCalorimeterHit() << "  " << "( " << (*i)->getCalorimeterHit()->getPosition()[0] << "," 
	      << (*i)->getCalorimeterHit()->getPosition()[1] << "," << (*i)->getCalorimeterHit()->getPosition()[2] << " )" << "  "
	      << "l on helix = " << (*i)->getPathLengthOnHelix() << "  " << "d to helix = " 
	      << (*i)->getDistanceToHelix() << std::endl;
    */
    ++index;
    


    // set up proper cylindrical tube in the 3 calorimeter zones to collect the MIP stub hits
    switch ( (*i)->getCalorimeterHit()->getType() ) {

      case 0 : maximalRadiusOfInnerTubeForMIPLikeStubCompare = _maximalRadiusOfInnerTubeForMIPLikeStub.at(0);
               minimalRadiusOfOuterTubeForMIPLikeStubCompare = _minimalRadiusOfOuterTubeForMIPLikeStub.at(0);
	       break;
      case 1 : maximalRadiusOfInnerTubeForMIPLikeStubCompare = _maximalRadiusOfInnerTubeForMIPLikeStub.at(1);
	       minimalRadiusOfOuterTubeForMIPLikeStubCompare = _minimalRadiusOfOuterTubeForMIPLikeStub.at(1);
	       break;
      case 2 : maximalRadiusOfInnerTubeForMIPLikeStubCompare = _maximalRadiusOfInnerTubeForMIPLikeStub.at(2);
	       minimalRadiusOfOuterTubeForMIPLikeStubCompare = _minimalRadiusOfOuterTubeForMIPLikeStub.at(2);
	       break;

    }


    // this calorimeter hits belongs to the MIP stub -> collect it 
    if ( (*i)->getDistanceToHelix() <= maximalRadiusOfInnerTubeForMIPLikeStubCompare ) {

      // debug
      // std::cout << "LAST HIT TAKEN" << std::endl << std::endl;

      
      // set up 'storeFirstHit' if the first hit in MIP stub is found
      if ( !firstHitFound ) {
	
	storeFirstHit = (*i); 
	firstHitFound = true;
	
      }
      
      clusterWithAttributes->addHit((*i)->getCalorimeterHit(),1.0);
      
      energy += (*i)->getCalorimeterHit()->getEnergy();
      storeLastHit = (*i);      
      
    }

    else {
      
      // don't look at hits outside 'veto-tube'
      if ( (*i)->getDistanceToHelix() > minimalRadiusOfOuterTubeForMIPLikeStubCompare )  {

	// debug
	// std::cout << "last hit outside VETO cylinder tube" << std::endl << std::endl;
	
	continue;
	
      }



      /*
      // already the closest hit to track is located outside the inner cylinder and inside the outer cylinder. -> take the this hit projection as start and end point
      if (i == calorimeterHitsWithAttributes.begin()) {

	LCVector3D caloHitPosition((*i)->getCalorimeterHit()->getPosition()[0],(*i)->getCalorimeterHit()->getPosition()[1],(*i)->getCalorimeterHit()->getPosition()[2]);
	LCVector3D caloHitPositionProjected = getProjectedPointOnHelix(caloHitPosition,fittedHelix);
	double s = getPathLengthOnHelix(caloHitPositionProjected,fittedHelix);
	LCVector3D tangent = fittedHelix->getDirection(s);

	clusterWithAttributes->setPositionEndHitLCVec(caloHitPositionProjected);
	clusterWithAttributes->setPositionStartHitLCVec(caloHitPositionProjected);
    
	clusterWithAttributes->setDirectionEndHitLCVec(tangent);
	clusterWithAttributes->setDirectionStartHitLCVec(tangent);
	
	clusterWithAttributes->setPositionLCVec(caloHitPositionProjected);

	clusterWithAttributes->setTypeEndHit((*i)->getCalorimeterHit()->getType());
	clusterWithAttributes->setTypeStartHit((*i)->getCalorimeterHit()->getType());   

	clusterWithAttributes->setIsMIPStub(false);
	clusterWithAttributes->setIsMuon(false);

	clusterWithAttributes->setEnergy(0.0);





	//debug
	if (_drawOnCED){
	  ced_hit ( caloHitPositionProjected.x(),caloHitPositionProjected.y(),caloHitPositionProjected.z(), 0 | 10 << CED_LAYER_SHIFT, 3, 0xff26c9);        
	}

	if ( _debugLevel > 5 ) {
	  std::cout << "2. case: ( (*i)->getDistanceToHelix() > maximalRadiusOfInnerTubeForMIPLikeStubCompare ) && (i == calorimeterHitsWithAttributes.begin()):" << std::endl
		    << "pos calorimeter hit = " << "( " << caloHitPosition.x() << "," << caloHitPosition.y()  << "," << caloHitPosition.z() << " )" << "  " 
		    << "pos hit projected = " << "( " << caloHitPositionProjected.x() << "," << caloHitPositionProjected.y()  << "," 
		    << caloHitPositionProjected.z() << " )" 
		    << std::endl << "type start hit = " << clusterWithAttributes->getTypeStartHit() << "  " << "type end hit = " << clusterWithAttributes->getTypeEndHit()
		    << std::endl << std::endl;
	}


	return;

      }

      */

      // a calorimeter hit located outside the inner cylinder and inside the outer cylinder was found. -> take the information of the calorimeter hit before, which
      // was located in the MIP stub and take this hit as end point and the first hit as a start point 
      // else {
     
      
      ++numberOfHitsInVetoTube;

      // FIXME: remove hard coded number for number of hits allowed in veto tube
      if ( numberOfHitsInVetoTube > 0 ) { // if more than zero hit is found in the veto tube break MIP stub searching
	// break if no calorimeter hit has fullfilled the MIP criteria
	if ( (clusterWithAttributes->getClusterImpl()->getCalorimeterHits().size()) == 0 ) break;

	LCVector3D caloHitPositionFirst(storeFirstHit->getCalorimeterHit()->getPosition()[0],
					storeFirstHit->getCalorimeterHit()->getPosition()[1],
					storeFirstHit->getCalorimeterHit()->getPosition()[2]);
	LCVector3D caloHitPositionProjected = getProjectedPointOnHelix(caloHitPositionFirst,fittedHelix);
	double sFirst = getPathLengthOnHelix(caloHitPositionProjected,fittedHelix);
	LCVector3D tangentFirst = fittedHelix->getDirection(sFirst);

	clusterWithAttributes->setPositionStartHitLCVec(caloHitPositionFirst);
	clusterWithAttributes->setDirectionStartHitLCVec(tangentFirst);
	clusterWithAttributes->setTypeStartHit(storeFirstHit->getCalorimeterHit()->getType());


	LCVector3D caloHitPositionLast(storeLastHit->getCalorimeterHit()->getPosition()[0],
				       storeLastHit->getCalorimeterHit()->getPosition()[1],
				       storeLastHit->getCalorimeterHit()->getPosition()[2]);
	LCVector3D caloHitPositionLastProjected = getProjectedPointOnHelix(caloHitPositionLast,fittedHelix);
	double sLast = getPathLengthOnHelix(caloHitPositionLastProjected,fittedHelix);
	LCVector3D tangentLast = fittedHelix->getDirection(sLast);

	clusterWithAttributes->setPositionEndHitLCVec(caloHitPositionLast);
	clusterWithAttributes->setDirectionEndHitLCVec(tangentLast);
	clusterWithAttributes->setTypeEndHit(storeLastHit->getCalorimeterHit()->getType());

       
	unsigned int n = clusterWithAttributes->getClusterImpl()->getCalorimeterHits().size();

	float* x = new float[n];
	float* y = new float[n];
	float* z = new float[n];
	float* a = new float[n];

	for (unsigned int j = 0; j < n; ++j) {

	  x[j] = clusterWithAttributes->getClusterImpl()->getCalorimeterHits().at(j)->getPosition()[0];
	  y[j] = clusterWithAttributes->getClusterImpl()->getCalorimeterHits().at(j)->getPosition()[1];
	  z[j] = clusterWithAttributes->getClusterImpl()->getCalorimeterHits().at(j)->getPosition()[2];
	  a[j] = clusterWithAttributes->getClusterImpl()->getCalorimeterHits().at(j)->getEnergy();

	}

      
	ClusterShapes* clusterShape = new ClusterShapes(n,a,x,y,z);


	LCVector3D position(clusterShape->getCentreOfGravity()[0],clusterShape->getCentreOfGravity()[1],clusterShape->getCentreOfGravity()[2]);
	clusterWithAttributes->setPositionLCVec(position);

	clusterWithAttributes->setIsMIPStub(true);
	clusterWithAttributes->setIsMuon(false);

	clusterWithAttributes->setEnergy(energy);

	delete clusterShape;
	clusterShape = 0;

	delete[] x;
	x = 0;
	delete[] y;
	y = 0;
	delete[] z;
	z = 0;
	delete[] a;
	a = 0;



	//debug
	if (_drawOnCED) {
	  ced_hit ( caloHitPositionFirst.x(),caloHitPositionFirst.y(),caloHitPositionFirst.z(), 0 | 10 << CED_LAYER_SHIFT, 3, 0xff26c9);
	  ced_hit ( position.x(),position.y(),position.z(), 0 | 10 << CED_LAYER_SHIFT, 3, 0xff26c9);
	  ced_hit ( caloHitPositionLast.x(),caloHitPositionLast.y(),caloHitPositionLast.z(), 0 | 10 << CED_LAYER_SHIFT, 3, 0xff26c9);
	}

	if ( _debugLevel > 5 ) {
	  std::cout << "3. case: ( (*i)->getDistanceToHelix() > maximalRadiusOfInnerTubeForMIPLikeStubCompare ) && (i != calorimeterHitsWithAttributes.begin()):" << std::endl
		    << "pos first hit = " << "( " << caloHitPositionFirst.x() << "," << caloHitPositionFirst.y()  << "," << caloHitPositionFirst.z() << " )" << "  " 
		    << "pos = " << "( " << position.x() << "," << position.y()  << "," << position.z() << " )" << "  " 
		    << "pos last hit = " << "( " << caloHitPositionLast.x() << "," << caloHitPositionLast.y()  << "," << caloHitPositionLast.z() << " )"
		    << std::endl << "type start hit = " << clusterWithAttributes->getTypeStartHit() << "  " << "type end hit = " << clusterWithAttributes->getTypeEndHit()
		    << std::endl << std::endl;
	}


	return;
       
      }

    }
    
  }

  // no calorimeter hits was found fullfilling the MIP stub conditions -> take hit with closest path length inside relaxed MIP stub criteria
  if ( clusterWithAttributes->getClusterImpl()->getCalorimeterHits().size() == 0 ) {

    // initialise with the closest by hit in the cone like tube
    CalorimeterHitWithAttributes* calorimeterHitWithAttributes = calorimeterHitsWithAttributes.at(0);

    for (i = calorimeterHitsWithAttributes.begin(); i!=calorimeterHitsWithAttributes.end(); ++i) {

      if ( (*i)->getDistanceToHelix() <= 2.0*minimalRadiusOfOuterTubeForMIPLikeStubCompare ) { // FIXME: hard coded factor

	calorimeterHitWithAttributes = (*i);
	break;
	
      }

    }

    LCVector3D caloHitPosition(calorimeterHitWithAttributes->getCalorimeterHit()->getPosition()[0],
			       calorimeterHitWithAttributes->getCalorimeterHit()->getPosition()[1],
			       calorimeterHitWithAttributes->getCalorimeterHit()->getPosition()[2]);
    LCVector3D caloHitPositionProjected = getProjectedPointOnHelix(caloHitPosition,fittedHelix);
    double s = getPathLengthOnHelix(caloHitPositionProjected,fittedHelix);
    LCVector3D tangent = fittedHelix->getDirection(s);
    
    clusterWithAttributes->setPositionEndHitLCVec(caloHitPositionProjected);
    clusterWithAttributes->setPositionStartHitLCVec(caloHitPositionProjected);
    
    clusterWithAttributes->setDirectionEndHitLCVec(tangent);
    clusterWithAttributes->setDirectionStartHitLCVec(tangent);
    
    clusterWithAttributes->setPositionLCVec(caloHitPositionProjected);
    
    clusterWithAttributes->setTypeEndHit(calorimeterHitWithAttributes->getCalorimeterHit()->getType());
    clusterWithAttributes->setTypeStartHit(calorimeterHitWithAttributes->getCalorimeterHit()->getType());   
    
    clusterWithAttributes->setIsMIPStub(false);
    clusterWithAttributes->setIsMuon(false);
    
    clusterWithAttributes->setEnergy(0.0);
    
    
    //debug
    if (_drawOnCED){ 
      ced_hit ( caloHitPositionProjected.x(),caloHitPositionProjected.y(),caloHitPositionProjected.z(), 0 | 10 << CED_LAYER_SHIFT, 3, 0xff26c9);
    }

    if ( _debugLevel > 5 ) {
      std::cout << "4. case: ( clusterWithAttributes->getClusterImpl()->getCalorimeterHits().size() == 0 ):" << std::endl
		<< "pos calorimeter hit = " << "( " << caloHitPosition.x() << "," << caloHitPosition.y()  << "," << caloHitPosition.z() << " )" << "  " 
		<< "pos hit projected = " << "( " << caloHitPositionProjected.x() << "," << caloHitPositionProjected.y()  << "," 
		<< caloHitPositionProjected.z() << " )" 
		<< std::endl << "type start hit = " << clusterWithAttributes->getTypeStartHit() << "  " << "type end hit = " << clusterWithAttributes->getTypeEndHit()
		<< std::endl << std::endl;
    }


    return;
     
  }


  // whole procedure was passed. -> all hits show a MIP like pattern -> MUON

  LCVector3D caloHitPositionFirst(storeFirstHit->getCalorimeterHit()->getPosition()[0],
				  storeFirstHit->getCalorimeterHit()->getPosition()[1],
				  storeFirstHit->getCalorimeterHit()->getPosition()[2]);
  LCVector3D caloHitPositionProjected = getProjectedPointOnHelix(caloHitPositionFirst,fittedHelix);
  double sFirst = getPathLengthOnHelix(caloHitPositionProjected,fittedHelix);
  LCVector3D tangentFirst = fittedHelix->getDirection(sFirst);
  
  clusterWithAttributes->setPositionStartHitLCVec(caloHitPositionFirst);
  clusterWithAttributes->setDirectionStartHitLCVec(tangentFirst);
  clusterWithAttributes->setTypeStartHit(storeFirstHit->getCalorimeterHit()->getType());
  
  
  LCVector3D caloHitPositionLast(storeLastHit->getCalorimeterHit()->getPosition()[0],
				 storeLastHit->getCalorimeterHit()->getPosition()[1],
				 storeLastHit->getCalorimeterHit()->getPosition()[2]);
  LCVector3D caloHitPositionLastProjected = getProjectedPointOnHelix(caloHitPositionLast,fittedHelix);
  double sLast = getPathLengthOnHelix(caloHitPositionLastProjected,fittedHelix);
  LCVector3D tangentLast = fittedHelix->getDirection(sLast);
  
  clusterWithAttributes->setPositionEndHitLCVec(caloHitPositionLast);
  clusterWithAttributes->setDirectionEndHitLCVec(tangentLast);
  clusterWithAttributes->setTypeEndHit(storeLastHit->getCalorimeterHit()->getType());
  
  
  unsigned int n = clusterWithAttributes->getClusterImpl()->getCalorimeterHits().size();
  
  float* x = new float[n];
  float* y = new float[n];
  float* z = new float[n];
  float* a = new float[n];
  
  for (unsigned int j = 0; j < n; ++j) {
    
    x[j] = clusterWithAttributes->getClusterImpl()->getCalorimeterHits().at(j)->getPosition()[0];
    y[j] = clusterWithAttributes->getClusterImpl()->getCalorimeterHits().at(j)->getPosition()[1];
    z[j] = clusterWithAttributes->getClusterImpl()->getCalorimeterHits().at(j)->getPosition()[2];
    a[j] = clusterWithAttributes->getClusterImpl()->getCalorimeterHits().at(j)->getEnergy();
    
  }

      
  ClusterShapes* clusterShape = new ClusterShapes(n,a,x,y,z);
  
  
  LCVector3D position(clusterShape->getCentreOfGravity()[0],clusterShape->getCentreOfGravity()[1],clusterShape->getCentreOfGravity()[2]);
  
  clusterWithAttributes->setPositionLCVec(position);
  
  clusterWithAttributes->setIsMIPStub(true);
  clusterWithAttributes->setIsMuon(false);
  
  clusterWithAttributes->setEnergy(energy);
  
  delete clusterShape;
  clusterShape = 0;
  
  delete[] x;
  x = 0;
  delete[] y;
  y = 0;
  delete[] z;
  z = 0;
  delete[] a;
  a = 0;


  
  //debug
  if (_drawOnCED){
    ced_hit ( caloHitPositionFirst.x(),caloHitPositionFirst.y(),caloHitPositionFirst.z(), 0 | 10 << CED_LAYER_SHIFT, 3, 0xff26c9);
    ced_hit ( position.x(),position.y(),position.z(), 0 | 10 << CED_LAYER_SHIFT, 3, 0xff26c9);
    ced_hit ( caloHitPositionLast.x(),caloHitPositionLast.y(),caloHitPositionLast.z(), 0 | 10 << CED_LAYER_SHIFT, 3, 0xff26c9);
  }

  if ( _debugLevel > 5 ) {
  
    std::cout << "5. case:  whole procedure was passed. -> all hits show a MIP like pattern -> MUON:" << std::endl
	      << "pos first hit = " << "( " << caloHitPositionFirst.x() << "," << caloHitPositionFirst.y()  << "," << caloHitPositionFirst.z() << " )" << "  " 
	      << "pos = " << "( " << position.x() << "," << position.y()  << "," << position.z() << " )" << "  " 
	      << "pos last hit = " << "( " << caloHitPositionLast.x() << "," << caloHitPositionLast.y()  << "," << caloHitPositionLast.z() << " )"
	      << std::endl << "type start hit = " << clusterWithAttributes->getTypeStartHit() << "  " << "type end hit = " << clusterWithAttributes->getTypeEndHit()
	      << std::endl << std::endl;
    
  }

}


//  _________________________________________________________________________________________________________________________________________________________________________





const std::vector<CalorimeterHitWithAttributes*> TrackBasedPFlow::removeMIPStub(std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes,
										ClusterImplWithAttributes* mipStub) {


  std::vector<CalorimeterHitWithAttributes*> remainingCalorimeterHitWithAttributes;

  for (std::vector<CalorimeterHitWithAttributes*>::const_iterator i = calorimeterHitsWithAttributes.begin(); i != calorimeterHitsWithAttributes.end(); ++i) {
  
    CalorimeterHitVec::const_iterator position = find(mipStub->getClusterImpl()->getCalorimeterHits().begin(),mipStub->getClusterImpl()->getCalorimeterHits().end(),
						      (*i)->getCalorimeterHit());

    if (position == mipStub->getClusterImpl()->getCalorimeterHits().end()) {

      remainingCalorimeterHitWithAttributes.push_back(*i);

       //debug
      if (_drawOnCED){

	ced_hit ( (*i)->getCalorimeterHit()->getPosition()[0],
		  (*i)->getCalorimeterHit()->getPosition()[1],
		  (*i)->getCalorimeterHit()->getPosition()[2], 2 | 5 << CED_LAYER_SHIFT, 8, 0x3bff89);


      }

    }
    
  }

  return remainingCalorimeterHitWithAttributes;
  
}



//  _________________________________________________________________________________________________________________________________________________________________________


std::vector<ClusterImplWithAttributes*> TrackBasedPFlow::doTrackwiseClustering(const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix,
									       const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes,
									       const double* startPoint, const float pathLengthOnHelixOfStartPoint,
									       const float distanceToHelixOfStartPoint, const double* startDirection) {


  std::vector<ClusterImpl*> resultingClusterImpls;
  std::vector<ClusterImplWithAttributes*> resultingClusterImplsWithAttributes;


  TrackwiseClustersParameters trackwiseClustersParameters;

  trackwiseClustersParameters.distanceTrackBack         = _distanceTrackBack;
  trackwiseClustersParameters.stepTrackBack             = _stepTrackBack;
  trackwiseClustersParameters.resolutionParameter       = _resolutionParameter;
  trackwiseClustersParameters.distanceMergeForward      = _distanceMergeForward;
  trackwiseClustersParameters.distanceToTrackSeed       = _distanceToTrackSeed;
  trackwiseClustersParameters.distanceToDefineDirection = _distanceToDefineDirection;
  trackwiseClustersParameters.resolutionToMerge         = _resolutionToMerge;
  trackwiseClustersParameters.nhit_merge_forward        = _nhit_merge_forward; 
  trackwiseClustersParameters.nhit_minimal              = _nhit_minimal;
  trackwiseClustersParameters.typeOfGenericDistance     = _typeOfGenericDistance;
  trackwiseClustersParameters.doMerging                 = _doMerging;
  trackwiseClustersParameters.doMergingForward          = _doMergingForward;
  trackwiseClustersParameters.displayClusters           = _displayClusters;
  trackwiseClustersParameters.NDefineSP                 = _NDefineSP;
  trackwiseClustersParameters.nScanToMergeForward       = _nScanToMergeForward;

  const TrackwiseClustersParameters* pTrackwiseClustersParameters = &trackwiseClustersParameters;



  TrackwiseClustersGeometryParameters trackwiseClustersGeometryParameters;
 
  trackwiseClustersGeometryParameters.zofendcap         = _zofendcap;
  trackwiseClustersGeometryParameters.rofbarrel         = _rofbarrel;
  trackwiseClustersGeometryParameters.phiofbarrel       = _phiofbarrel;
  trackwiseClustersGeometryParameters.nsymmetry         = _nsymmetry;
  trackwiseClustersGeometryParameters.thetaofendcap     = _thetaofendcap;

  trackwiseClustersGeometryParameters.weightForReso     = _weightForReso;
  trackwiseClustersGeometryParameters.weightForDist     = _weightForDist;
  trackwiseClustersGeometryParameters.bField            = _bField;

  const TrackwiseClustersGeometryParameters* pTrackwiseClustersGeometryParameters = &trackwiseClustersGeometryParameters;


  // FIXME: need (Clustering which is able to take doubles)
  
  float startP[3];
  float startDir[3];
  for (int i = 0; i < 3; ++i) {

    startP[i] = (float)(startPoint[i]);
    startDir[i] = (float)(startDirection[i]);

    // debug
    // std::cout << "TrackwiseClustering for charged particles : " << startP[i] << "  " << startDir[i] << std::endl;

  }

  
  TrackwiseClusters* clusters = new TrackwiseClusters(calorimeterHitsWithAttributes,startP,pathLengthOnHelixOfStartPoint,distanceToHelixOfStartPoint,
						      startDir,pTrackwiseClustersParameters,pTrackwiseClustersGeometryParameters);


  resultingClusterImpls = clusters->doClustering();



  // FIXME build function/class for all stuff below  
  // assign properties of the attributes of the cluster

  for (std::vector<ClusterImpl*>::const_iterator i = resultingClusterImpls.begin(); i != resultingClusterImpls.end(); ++i) {
                                             
    ClusterImplWithAttributes* clusterImplWithAttributes = new ClusterImplWithAttributes();

    clusterImplWithAttributes->setClusterImpl(*i);

    resultingClusterImplsWithAttributes.push_back(clusterImplWithAttributes);
    
  }

  delete clusters;
  clusters = 0;


  // fill attributes of clusters
  assignClusterAttributes(outermostTrackerHits,fittedHelix,resultingClusterImplsWithAttributes,startPoint);

  return resultingClusterImplsWithAttributes;

}


std::vector<ClusterImplWithAttributes*> TrackBasedPFlow::doConeClustering(const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix,
									  const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, 
									  const double* startPoint, const double* startDirection) {


  std::vector<ClusterImpl*> resultingClusterImpls;
  std::vector<ClusterImplWithAttributes*> resultingClusterImplsWithAttributes;

  float openingAngleConeTubeRed = (constants::twopi/360.0) * 60.0; // initial: 100 degrees
  float maxPathLengthInTubeRed = 1000.0; // initial: 1000 mm
  float openingAngleConeTubeGreen = (constants::twopi/360.0) * 40.0; // initial: 80 degrees
  float maxPathLengthInTubeGreen = 1000.0; // initial: 1000 mm
  float openingAngleConeTubeBlue = (constants::twopi/360.0) * 20.0; // initial: 40 degrees
  float maxPathLengthInTubeBlue = 1000.0; // initial: 1000 mm


  LCVector3D referencePosition(outermostTrackerHits.at(0)->getPosition()[0],outermostTrackerHits.at(0)->getPosition()[1],outermostTrackerHits.at(0)->getPosition()[2]);
  LCVector3D referencePositionProjected = getProjectedPointOnHelix(referencePosition,fittedHelix);

  LCVector3D startPosition(startPoint[0],startPoint[1],startPoint[2]);
  LCVector3D startPositionProjected = getProjectedPointOnHelix(startPosition,fittedHelix);
  
  double pathLengthOnHelixOfStartPoint = getPathLengthOnHelix(referencePositionProjected,startPositionProjected,fittedHelix);



  // debug

  if (_drawOnCED) ced_hit ( startPosition.x(), startPosition.y(), startPosition.z(), 0 | 11 << CED_LAYER_SHIFT, 6, 0xff8ed2 );
  //  double distanceOfStartPoint = getDistanceToHelix(startPosition,fittedHelix);
  //  std::cout << "pathLengthOnHelixOfStartPoint = " << pathLengthOnHelixOfStartPoint << "  " << "distanceOfStartPointToHelix = " << distanceOfStartPoint << std::endl;



  ClusterImpl* cluster = new ClusterImpl(); // only one cluster is assigned to each track

  // assign blue hits to cluster
  for (std::vector<CalorimeterHitWithAttributes*>::const_iterator i = calorimeterHitsWithAttributes.begin(); i!= calorimeterHitsWithAttributes.end(); ++i) {
    
    float distanceOfHit = (*i)->getDistanceToHelix();
    float pathLengthOfHitToStartPoint = ( (*i)->getPathLengthOnHelix() ) - pathLengthOnHelixOfStartPoint;
    
    // std::cout << "pathLengthOfHitToStartPoint = " << pathLengthOfHitToStartPoint << "  " << "distanceOfHit = " << distanceOfHit << "  " <<std::endl;




    double hitEnergyInMIPs = 0.0;
		
    switch ( (*i)->getCalorimeterHit()->getType() ) {
      
    case 0 : 
      hitEnergyInMIPs = ((*i)->getCalorimeterHit()->getEnergy())/_mipCoeffEcal.at(0);
      break;
    case 1 : 
      hitEnergyInMIPs = ((*i)->getCalorimeterHit()->getEnergy())/_mipCoeffEcal.at(1);
      break;
    case 2 : 
      hitEnergyInMIPs = ((*i)->getCalorimeterHit()->getEnergy())/_mipCoeffHcal.at(0);
      break;
      
    }
    
   
    if (hitEnergyInMIPs < 0.5) { 

      // do nothing?

    }
  
    else if ( (hitEnergyInMIPs >= 0.5) && (hitEnergyInMIPs < 1.7) ) {
              
      if ( ( distanceOfHit <= tan(openingAngleConeTubeRed/2) * pathLengthOfHitToStartPoint ) && ( pathLengthOfHitToStartPoint > 0.0 ) && 
	   ( pathLengthOfHitToStartPoint <= maxPathLengthInTubeRed ) ) {	
	
	// std::cout << "FOUND RED HIT" << std::endl;
	
	CalorimeterHit* caloHit = (*i)->getCalorimeterHit();      
	cluster->addHit(caloHit,(float)1.0);       
	 
	// debug
	if (_drawOnCED) {

	  ced_hit ( caloHit->getPosition()[0], caloHit->getPosition()[1], caloHit->getPosition()[2], 1 | 12 << CED_LAYER_SHIFT, 6, 0xff0000 );

	}

      }

    }
 
    else if ( (hitEnergyInMIPs >= 1.7) && (hitEnergyInMIPs < 3.5) ) {
      
      if ( ( distanceOfHit <= tan(openingAngleConeTubeGreen/2) * pathLengthOfHitToStartPoint ) && ( pathLengthOfHitToStartPoint > 0.0 ) && 
	   ( pathLengthOfHitToStartPoint <= maxPathLengthInTubeGreen ) ) {	
	
	// std::cout << "FOUND GREEN HIT" << std::endl;
	
	CalorimeterHit* caloHit = (*i)->getCalorimeterHit();      
	cluster->addHit(caloHit,(float)1.0);       
	 
	// debug
	if (_drawOnCED) {

	  ced_hit ( caloHit->getPosition()[0], caloHit->getPosition()[1], caloHit->getPosition()[2], 1 | 13 << CED_LAYER_SHIFT, 6, 0x0dff00 );

	}

      }
      
    }

    else {
      
      if ( ( distanceOfHit <= tan(openingAngleConeTubeBlue/2) * pathLengthOfHitToStartPoint ) && ( pathLengthOfHitToStartPoint > 0.0 ) && 
	   ( pathLengthOfHitToStartPoint <= maxPathLengthInTubeBlue ) ) {	
	
	// std::cout << "FOUND BLUE HIT" << std::endl;
	
	CalorimeterHit* caloHit = (*i)->getCalorimeterHit();      
	cluster->addHit(caloHit,(float)1.0);       
	 
	// debug
	if (_drawOnCED) {

	  ced_hit ( caloHit->getPosition()[0], caloHit->getPosition()[1], caloHit->getPosition()[2], 1 | 14 << CED_LAYER_SHIFT, 6, 0x2f6dff );

	}

      }

    }
    
  }
  
  assignClusterProperties(cluster);
  
  if ( cluster->getCalorimeterHits().size() > 0 ) resultingClusterImpls.push_back(cluster);
  

  // FIXME build function/class for all stuff below  
  // assign properties of the attributes of the cluster

  for (std::vector<ClusterImpl*>::const_iterator i = resultingClusterImpls.begin(); i != resultingClusterImpls.end(); ++i) {
    
    ClusterImplWithAttributes* clusterImplWithAttributes = new ClusterImplWithAttributes();
    
    clusterImplWithAttributes->setClusterImpl(*i);

    resultingClusterImplsWithAttributes.push_back(clusterImplWithAttributes);
    
  }


  // fill attributes of clusters
  assignClusterAttributes(outermostTrackerHits,fittedHelix,resultingClusterImplsWithAttributes,startPoint);

  return resultingClusterImplsWithAttributes;

}


std::vector<ClusterImplWithAttributes*> TrackBasedPFlow::doNNClustering(const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix,
									const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, 
									const double* startPoint, const double* startDirection) {


  std::vector<ClusterImpl*> resultingClusterImpls;
  std::vector<ClusterImplWithAttributes*> resultingClusterImplsWithAttributes;



  // FIXME: delete this somewhere????
  LCCollectionVec* CaloHits = new LCCollectionVec(LCIO::CALORIMETERHIT);


  for (std::vector<CalorimeterHitWithAttributes*>::const_iterator i = calorimeterHitsWithAttributes.begin(); i!= calorimeterHitsWithAttributes.end(); ++i) {
    
    CaloHits->addElement( (*i)->getCalorimeterHit());

  }




  // taken from NNClusterProcesssor.cc

  LCCollectionVec* lcioClusters = new LCCollectionVec( LCIO::CLUSTER )  ;
  
  GenericHitVec<CalorimeterHit> h ;

  GenericClusterVec<CalorimeterHit> cl ;
  
  EnergyCut<CalorimeterHit> eCut( 0.0 ) ; // parameter to set (default 0.0)
  
  ZIndex<CalorimeterHit,100> zIndex( -4300. , 4300. ) ; 

  NNDistance< CalorimeterHit, float> dist( _resolutionParameter.at(0) )  ; // parameter to set (default 40.0)

  LCIOCluster<CalorimeterHit> converter ;
  

  // create a vector of generic hits from the collection applying an energy cut
  // int nHit = CaloHits->getNumberOfElements();

  addToGenericHitVec( h , CaloHits , eCut , zIndex ) ;

  // cluster the hits with a nearest neighbour condition
  cluster( h.begin() , h.end() , std::back_inserter( cl )  , &dist ) ;


  /*
  std::cout << "  passing " << h.size() << " of " << nHit  
	    << "  hits to clustering (E_cut: " << eCut << ") " 
	    << "  found  " << cl.size() << " clusters " << std::endl ;
  */

  // create lcio::Clusters from the clustered GenericHits
  std::transform( cl.begin(), cl.end(), std::back_inserter( *lcioClusters ) , converter ) ;

  // END (taken from NNClusterProcesssor.cc)




  // FIXME: not very elegant below, change this

  for (int i = 0; i < lcioClusters->getNumberOfElements(); ++i) {
  
    ClusterImpl* cluster = new ClusterImpl();
    Cluster* clusterFromNNClustering = dynamic_cast<Cluster*>(lcioClusters->getElementAt(i));
    
    for (unsigned int j = 0; j < clusterFromNNClustering->getCalorimeterHits().size(); ++j) {
      
      CalorimeterHit* caloHit = clusterFromNNClustering->getCalorimeterHits().at(j);
      cluster->addHit(caloHit,float(1.0));

    }

    assignClusterProperties(cluster);
    resultingClusterImpls.push_back(cluster);

  }

  delete lcioClusters;
  lcioClusters = 0;


  
  for (std::vector<ClusterImpl*>::const_iterator i = resultingClusterImpls.begin(); i != resultingClusterImpls.end(); ++i) {
                                             
    ClusterImplWithAttributes* clusterImplWithAttributes = new ClusterImplWithAttributes();

    clusterImplWithAttributes->setClusterImpl(*i);

    resultingClusterImplsWithAttributes.push_back(clusterImplWithAttributes);
    
  }



  // fill attributes of clusters
  assignClusterAttributes(outermostTrackerHits,fittedHelix,resultingClusterImplsWithAttributes,startPoint);

  return resultingClusterImplsWithAttributes;


}




std::vector<ClusterImplWithAttributes*> TrackBasedPFlow::doTrackBasedClustering(const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix,
										const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, 
										const double* startPoint, const double* startDirection) {


  std::vector<ClusterImplWithAttributes*> returnValue;
  
  return returnValue;

}




std::vector<ClusterImplWithAttributes*> TrackBasedPFlow::doTrackwiseClusteringForNeutrals(const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, 
											  const double* startPoint, const double* startDirection) {
  

  std::vector<ClusterImpl*> resultingClusterImpls;
  std::vector<ClusterImplWithAttributes*> resultingClusterImplsWithAttributes;


  TrackwiseClustersParameters trackwiseClustersParameters;

  trackwiseClustersParameters.distanceTrackBack         = _distanceTrackBack;
  trackwiseClustersParameters.stepTrackBack             = _stepTrackBack;
  trackwiseClustersParameters.resolutionParameter       = _resolutionParameterForNeutrals;
  trackwiseClustersParameters.distanceMergeForward      = _distanceMergeForward;
  trackwiseClustersParameters.distanceToTrackSeed       = _distanceToTrackSeed;
  trackwiseClustersParameters.distanceToDefineDirection = _distanceToDefineDirection;
  trackwiseClustersParameters.resolutionToMerge         = _resolutionToMerge;
  trackwiseClustersParameters.nhit_merge_forward        = _nhit_merge_forward; 
  trackwiseClustersParameters.nhit_minimal              = _nhit_neutral_minimal;
  trackwiseClustersParameters.typeOfGenericDistance     = _typeOfGenericDistanceForNeutrals;
  trackwiseClustersParameters.doMerging                 = _doMerging;
  trackwiseClustersParameters.doMergingForward          = _doMergingForward;
  trackwiseClustersParameters.displayClusters           = _displayClusters;
  trackwiseClustersParameters.NDefineSP                 = _NDefineSP;
  trackwiseClustersParameters.nScanToMergeForward       = _nScanToMergeForward;

  const TrackwiseClustersParameters* pTrackwiseClustersParameters = &trackwiseClustersParameters;



  TrackwiseClustersGeometryParameters trackwiseClustersGeometryParameters;
 
  trackwiseClustersGeometryParameters.zofendcap         = _zofendcap;
  trackwiseClustersGeometryParameters.rofbarrel         = _rofbarrel;
  trackwiseClustersGeometryParameters.phiofbarrel       = _phiofbarrel;
  trackwiseClustersGeometryParameters.nsymmetry         = _nsymmetry;
  trackwiseClustersGeometryParameters.thetaofendcap     = _thetaofendcap;

  trackwiseClustersGeometryParameters.weightForReso     = _weightForResoForNeutrals;
  trackwiseClustersGeometryParameters.weightForDist     = _weightForDistForNeutrals;
  trackwiseClustersGeometryParameters.bField            = _bField;

  const TrackwiseClustersGeometryParameters* pTrackwiseClustersGeometryParameters = &trackwiseClustersGeometryParameters;


  // FIXME: need (Clustering which is able to take doubles)
  
  float startP[3];
  float startDir[3];
  for (int i = 0; i < 3; ++i) {
    
    startP[i] = (float)(startPoint[i]);
    startDir[i] = (float)(startDirection[i]);

    // debug
    // std::cout << "TrackwiseClustering for neutrals : " << startP[i] << "  " << startDir[i] << std::endl;

  }
  

  
  TrackwiseClusters* clusters = new TrackwiseClusters(calorimeterHitsWithAttributes,startP,0.0,0.0,startDir,
						      pTrackwiseClustersParameters,pTrackwiseClustersGeometryParameters);


  resultingClusterImpls = clusters->doClustering();



  // FIXME build function/class for all stuff below  
  // assign properties of the attributes of the cluster

  for (std::vector<ClusterImpl*>::const_iterator i = resultingClusterImpls.begin(); i != resultingClusterImpls.end(); ++i) {

    ClusterImplWithAttributes* clusterImplWithAttributes = new ClusterImplWithAttributes();

    clusterImplWithAttributes->setClusterImpl(*i);
    
    resultingClusterImplsWithAttributes.push_back(clusterImplWithAttributes);
    
  }

  delete clusters;
  clusters = 0;

  return resultingClusterImplsWithAttributes;


}




ReconstructedParticleImpl* TrackBasedPFlow::assignClustersToTrack(Track* track, ClusterImplWithAttributes* mipStub, std::vector<ClusterImplWithAttributes*> clusters, 
								  const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix) {


  ReconstructedParticleImpl* recoParticle = new ReconstructedParticleImpl();
  
  recoParticle->addTrack(track);

  
  if ( mipStub->isMuon() ) {

    recoParticle->addCluster(mipStub->getClusterImpl());
    return recoParticle;
    
  }

  
  if ( mipStub->isMIPStub() ) {

    recoParticle->addCluster(mipStub->getClusterImpl());

  }
  
  std::vector<ClusterImplWithAttributes*> relatedClusters = assignClusters(mipStub,clusters,mipStub->getPositionEndHitLCVec(),outermostTrackerHits,fittedHelix);



  // find clusters which are not assigned to track
  std::vector<ClusterImplWithAttributes*> notAssignedClusters;
  
  for(std::vector<ClusterImplWithAttributes*>::const_iterator i = clusters.begin(); i != clusters.end(); ++i) {
     
    std::vector<ClusterImplWithAttributes*>::const_iterator position = find(relatedClusters.begin(),relatedClusters.end(),(*i));   
    if ( position == relatedClusters.end() ) notAssignedClusters.push_back(*i);
    
  }



  // find additional Clusters
  std::vector<ClusterImplWithAttributes*> additionalRelatedClusters = assignAdditionalClusters(mipStub,fittedHelix,relatedClusters,notAssignedClusters,track);
  


  // remove MIP stub again, because it already has been added to the RecoParticle
  relatedClusters.erase(relatedClusters.begin());


  // add related Clusters to RecoParticle and to Cluster collection
  for(std::vector<ClusterImplWithAttributes*>::const_iterator i = relatedClusters.begin(); i != relatedClusters.end(); ++i) {
  
    recoParticle->addCluster((*i)->getClusterImpl());

  }
  

  // add additionally found clusters
  for(std::vector<ClusterImplWithAttributes*>::const_iterator i = additionalRelatedClusters.begin(); i != additionalRelatedClusters.end(); ++i) {
  
    recoParticle->addCluster((*i)->getClusterImpl());

  }

  return recoParticle;

}


ReconstructedParticleImpl* TrackBasedPFlow::assignClustersToTrackPerfectly(Track* track, std::vector<ClusterImplWithAttributes*> clusters) {


  ReconstructedParticleImpl* recoParticle = new ReconstructedParticleImpl();


  // perfect PFlow for charged particles -> assign all clusters and the track to the reco particle

  recoParticle->addTrack(track);

  for(std::vector<ClusterImplWithAttributes*>::const_iterator i = clusters.begin(); i != clusters.end(); ++i) {

    recoParticle->addCluster((*i)->getClusterImpl());

  }

  return recoParticle;

}



ReconstructedParticleImpl* TrackBasedPFlow::assignNeutralClusterToReconstructedParticle(ClusterImplWithAttributes* cluster) {


  ReconstructedParticleImpl* recoParticle = new ReconstructedParticleImpl();
  
  recoParticle->addCluster(cluster->getClusterImpl());
  
  return recoParticle;
  
}


void TrackBasedPFlow::assignClusterProperties(ClusterImpl* cluster) {

  int nCaloHitsAssigned = cluster->getCalorimeterHits().size();

  if ( nCaloHitsAssigned > 0 ) {

    float totalEnergy = 0.0;
    float totalEnergyECAL = 0.0;
    float totalEnergyHCAL = 0.0;
    
    float* xHits = new float[nCaloHitsAssigned];
    float* yHits = new float[nCaloHitsAssigned];
    float* zHits = new float[nCaloHitsAssigned];
    float* aHits = new float[nCaloHitsAssigned];
  
    for (int i = 0; i < nCaloHitsAssigned; ++i) {

      CalorimeterHit* caloHit = cluster->getCalorimeterHits().at(i);
      
      xHits[i] = caloHit->getPosition()[0];
      yHits[i] = caloHit->getPosition()[1];
      zHits[i] = caloHit->getPosition()[2];
      
      // float energyContributionOfCalorimeterHit = cluster->getHitContributions().at(i); // in case of using this the energy contributions 
                                                                                          // of the CaloHits need to be set correctely
      float energyOfIntrinsicCalorimeterHit = caloHit->getEnergy();

      
      // debug, compare hit contribution with 'intrinsic' energy of the corresponding hit
      /*      
      std::cout << "Energy contribution of the Calorimeter hit : " << energyContributionOfCalorimeterHit << "  " 
		<< "'intrinsic' energy of the Calorimeter hit : " << energyOfIntrinsicCalorimeterHit << std::endl;
      */

      aHits[i] = energyOfIntrinsicCalorimeterHit;
      totalEnergy += energyOfIntrinsicCalorimeterHit;

      if ( (caloHit->getType() == 0) || (caloHit->getType() == 1) ) totalEnergyECAL += energyOfIntrinsicCalorimeterHit;
      else totalEnergyHCAL += energyOfIntrinsicCalorimeterHit;
      
    }
    
    cluster->setEnergy(totalEnergy);
    // in a cluster there only exist two types: type == 1 for the whole ECAL and type == 2 for the HCAL
    cluster->subdetectorEnergies().resize(2);
    cluster->subdetectorEnergies()[0] = totalEnergyECAL;
    cluster->subdetectorEnergies()[1] = totalEnergyHCAL;
    
    ClusterShapes* shape = new ClusterShapes(nCaloHitsAssigned,aHits,xHits,yHits,zHits);	    
    
    cluster->setPosition(shape->getCentreOfGravity());
    float phiCluster = atan2(shape->getEigenVecInertia()[1],shape->getEigenVecInertia()[0]);
    float thetaCluster = acos(shape->getEigenVecInertia()[2]);
    cluster->setIPhi(phiCluster);
    cluster->setITheta(thetaCluster);	    
    
    delete shape;
    
    delete[] xHits;
    xHits = 0;
    delete[] yHits;
    yHits = 0;
    delete[] zHits;
    zHits = 0;
    delete[] aHits;
    aHits = 0;
    
  }

}


void TrackBasedPFlow::assignClusterAttributes(const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix,
					      std::vector<ClusterImplWithAttributes*> resultingClusterImplsWithAttributes, 
					      const double* startPoint) {


  // find the hit in each cluster with the smallest path length on the helix and store its distance to the helix

  LCVector3D referencePosition(outermostTrackerHits.at(0)->getPosition()[0],outermostTrackerHits.at(0)->getPosition()[1],outermostTrackerHits.at(0)->getPosition()[2]);
  LCVector3D referencePositionProjected = getProjectedPointOnHelix(referencePosition,fittedHelix);

      
  for (std::vector<ClusterImplWithAttributes*>::const_iterator i = resultingClusterImplsWithAttributes.begin(); i!= resultingClusterImplsWithAttributes.end(); ++i) {

    double smallestPathLength = DBL_MAX;
    double largestPathLength = (-1.0)*DBL_MAX;
    double smallestDistanceToStartPoint = DBL_MAX;
    double largestDistanceToStartPoint = (-1.0)*DBL_MAX;
    
    CalorimeterHit* closeByHitPathLength;
    CalorimeterHit* farAwayHitPathLength;
    CalorimeterHit* closeByHitDistanceToStartPoint;
    CalorimeterHit* farAwayHitDistanceToStartPoint;   


    for(CalorimeterHitVec::const_iterator j = ((*i)->getClusterImpl()->getCalorimeterHits().begin()); j != ((*i)->getClusterImpl()->getCalorimeterHits().end()); ++j ) {
	

      LCVector3D calorimeterHitPosition((*j)->getPosition()[0],(*j)->getPosition()[1],(*j)->getPosition()[2]);
      LCVector3D calorimeterHitPositionProjected = getProjectedPointOnHelix(calorimeterHitPosition,fittedHelix);      

      double pathLength = getPathLengthOnHelix(referencePositionProjected,calorimeterHitPositionProjected,fittedHelix);

      double distanceToStartPoint = sqrt( pow( (startPoint[0] - (*j)->getPosition()[0]), 2) + pow( (startPoint[1] - (*j)->getPosition()[1]), 2) + 
					  pow( (startPoint[2] - (*j)->getPosition()[2]), 2) );
	

      if ( pathLength < smallestPathLength ) {

	smallestPathLength = pathLength;
	closeByHitPathLength = (*j);

      }

      if ( pathLength > largestPathLength ) {

	largestPathLength = pathLength;
	farAwayHitPathLength = (*j);

      }
      
      if ( distanceToStartPoint < smallestDistanceToStartPoint ) {

	smallestDistanceToStartPoint = distanceToStartPoint;
	closeByHitDistanceToStartPoint = (*j);

      }

      if ( distanceToStartPoint > largestDistanceToStartPoint ) {

	largestDistanceToStartPoint = distanceToStartPoint;
	farAwayHitDistanceToStartPoint = (*j);

      }


    }

    /*
    LCVector3D CaloHitPositionFirst(closeByHitDistanceToStartPoint->getPosition()[0],
				    closeByHitDistanceToStartPoint->getPosition()[1],
				    closeByHitDistanceToStartPoint->getPosition()[2]);
    */
    LCVector3D CaloHitPositionFirst(closeByHitPathLength->getPosition()[0],
				    closeByHitPathLength->getPosition()[1],
				    closeByHitPathLength->getPosition()[2]);


    LCVector3D CaloHitPositionFirstProjected = getProjectedPointOnHelix(CaloHitPositionFirst,fittedHelix);
    double sFirst = getPathLengthOnHelix(CaloHitPositionFirstProjected,fittedHelix);
    LCVector3D tangentFirst = fittedHelix->getDirection(sFirst);

    /*
    LCVector3D closeByHitDistanceToStartPointPosition(closeByHitDistanceToStartPoint->getPosition()[0],
						      closeByHitDistanceToStartPoint->getPosition()[1],
						      closeByHitDistanceToStartPoint->getPosition()[2]);
    */
    LCVector3D closeByHitDistanceToStartPointPosition(closeByHitPathLength->getPosition()[0],
						      closeByHitPathLength->getPosition()[1],
						      closeByHitPathLength->getPosition()[2]);

    (*i)->setPositionStartHitLCVec(closeByHitDistanceToStartPointPosition);
    (*i)->setDirectionStartHitLCVec(tangentFirst);
    (*i)->setTypeStartHit(closeByHitDistanceToStartPoint->getType());

    /*    
    LCVector3D CaloHitPositionSecond(farAwayHitDistanceToStartPoint->getPosition()[0],
				     farAwayHitDistanceToStartPoint->getPosition()[1],
				     farAwayHitDistanceToStartPoint->getPosition()[2]);
    */
    LCVector3D CaloHitPositionSecond(farAwayHitPathLength->getPosition()[0],
				     farAwayHitPathLength->getPosition()[1],
				     farAwayHitPathLength->getPosition()[2]);



    LCVector3D CaloHitPositionSecondProjected = getProjectedPointOnHelix(CaloHitPositionSecond,fittedHelix);
    double sLast = getPathLengthOnHelix(CaloHitPositionSecondProjected,fittedHelix);
    LCVector3D tangentLast = fittedHelix->getDirection(sLast);

    /*
    LCVector3D farAwayHitDistanceToStartPointPosition(farAwayHitDistanceToStartPoint->getPosition()[0],
						      farAwayHitDistanceToStartPoint->getPosition()[1],
						      farAwayHitDistanceToStartPoint->getPosition()[2]);
    */

    LCVector3D farAwayHitDistanceToStartPointPosition(farAwayHitPathLength->getPosition()[0],
						      farAwayHitPathLength->getPosition()[1],
						      farAwayHitPathLength->getPosition()[2]);

    (*i)->setPositionEndHitLCVec(farAwayHitDistanceToStartPointPosition);
    (*i)->setDirectionEndHitLCVec(tangentLast);
    (*i)->setTypeEndHit(farAwayHitDistanceToStartPoint->getType());

    (*i)->setIsMIPStub(false);
    (*i)->setIsMuon(false);

    
    // debug 
    /*    
      std::cout << "cluster : " << (*i) << std::endl << "smallestPathLength = " << smallestPathLength << "  " << "largestPathLength = " << largestPathLength << std::endl;
    
    if (_drawOnCED) {	
      ced_hit(closeByHitDistanceToStartPoint->getPosition()[0], closeByHitDistanceToStartPoint->getPosition()[1], closeByHitDistanceToStartPoint->getPosition()[2], 
      0 | 9 << CED_LAYER_SHIFT, 6, 0x08fff7 );
      ced_hit(farAwayHitDistanceToStartPoint->getPosition()[0], farAwayHitDistanceToStartPoint->getPosition()[1], farAwayHitDistanceToStartPoint->getPosition()[2], 
      0 | 9 << CED_LAYER_SHIFT, 6, 0x08fff7 );
    }
    */
    
  }

}

void TrackBasedPFlow::doPID(ReconstructedParticleImpl* recoParticle, bool isMuon) {


  if ( recoParticle->getTracks().size() > 1 ) {

    std::cout << "More than one track per reconstructed particle -> PID not possible at the moment" << std::endl;
    return;
    
  }

  // Muons
  if ( isMuon ) {

    Track* track = recoParticle->getTracks().at(0);
    
    float charge = ( track->getOmega() ) / std::fabs( track->getOmega() );
    float mass = 0.105658369;
      
    float momentum[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) momentum[i] = (float) ( (MarlinUtil::getMomentum(track,_bField))[i]);
      
    float energy = 0.0;
    for (int i = 0; i < 3; ++i) energy += momentum[i]*momentum[i];
    energy = sqrt( energy + mass*mass);
      

    recoParticle->setType(5);
    recoParticle->setCharge(charge);
    recoParticle->setMass(mass);
      
    recoParticle->setMomentum(momentum);
    recoParticle->setEnergy(energy);

  }


  // Charged Particles without energy deposition in calorimeter -> pion hypothesis
  if ( (recoParticle->getTracks().size() == 1) && (recoParticle->getClusters().size() == 0) ) {

    Track* track = recoParticle->getTracks().at(0);
    
    float charge = ( track->getOmega() ) / std::fabs( track->getOmega() );
    float mass = 0.13957;
      
    float momentum[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) momentum[i] = (float) ( (MarlinUtil::getMomentum(track,_bField))[i]);
      
    float energy = 0.0;
    for (int i = 0; i < 3; ++i) energy += momentum[i]*momentum[i];
    energy = sqrt( energy + mass*mass);
      

    recoParticle->setType(2);
    recoParticle->setCharge(charge);
    recoParticle->setMass(mass);
      
    recoParticle->setMomentum(momentum);
    recoParticle->setEnergy(energy);

  }


  // charged particles
  if ( recoParticle->getTracks().size() == 1 ) {

    Track* track = recoParticle->getTracks().at(0);
    
    float charge = ( track->getOmega() ) / std::fabs( track->getOmega() );

    float ecalEnergy = 0.0;
    float hcalEnergy = 0.0;
    float totalEnergy = 0.0;
    
    for (ClusterVec::const_iterator i = recoParticle->getClusters().begin(); i != recoParticle->getClusters().end(); ++i) {

      for (CalorimeterHitVec::const_iterator j = (*i)->getCalorimeterHits().begin(); j != (*i)->getCalorimeterHits().end(); ++j) {

	if ( ( (*j)->getType() == 0 ) || ( (*j)->getType() == 1 ) ) ecalEnergy += (*j)->getEnergy();
	
	else hcalEnergy += (*j)->getEnergy();
	
      }
      
    }
    
    totalEnergy = ecalEnergy + hcalEnergy;

    // if >95% of total energy is in the ecal -> electron 
    float fraction = ecalEnergy/fmax(totalEnergy,1.0e-6);
    if ( fraction > _fractionEM) {

      float mass = 5.11E-4;
      
      float momentum[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < 3; ++i) momentum[i] = (float) ( (MarlinUtil::getMomentum(track,_bField))[i]);
      
      float energy = 0.0;
      for (int i = 0; i < 3; ++i) energy += momentum[i]*momentum[i];
      energy = sqrt( energy + mass*mass);
      

      recoParticle->setType(1);
      recoParticle->setCharge(charge);
      recoParticle->setMass(mass);
      
      recoParticle->setMomentum(momentum);
      recoParticle->setEnergy(energy);
      
    }
    // else -> pion
    else {

      float mass = 0.13957;

      float momentum[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < 3; ++i) momentum[i] = (float) ( (MarlinUtil::getMomentum(track,_bField))[i]);

      float energy = 0.0;
      for (int i = 0; i < 3; ++i) energy += momentum[i]*momentum[i];
      energy = sqrt( energy + mass*mass);
      

      recoParticle->setType(2);
      recoParticle->setCharge(charge);
      recoParticle->setMass(mass);
      
      recoParticle->setMomentum(momentum);
      recoParticle->setEnergy(energy);
      
    }
        
  }


  // neutral particles
  else {

    float ecalEnergy = 0.0;
    float hcalEnergy = 0.0;
    float totalEnergy = 0.0;

    float position[3] = {0.0, 0.0, 0.0};
    

    for (ClusterVec::const_iterator i = recoParticle->getClusters().begin(); i != recoParticle->getClusters().end(); ++i) {


      // position as vector sum of positions
      for (int j = 0; j < 3; ++j) position[j] += (*i)->getPosition()[j];
      
      
      for (CalorimeterHitVec::const_iterator j = (*i)->getCalorimeterHits().begin(); j != (*i)->getCalorimeterHits().end(); ++j) {

	if ( ( (*j)->getType() == 0 ) || ( (*j)->getType() == 1 ) ) ecalEnergy += (*j)->getEnergy();
	
	else hcalEnergy += (*j)->getEnergy();
	
      }
      
    }
    
    totalEnergy = ecalEnergy + hcalEnergy;
    
    // if >95% of total energy is in the ecal -> gamma 
    float fraction = ecalEnergy/fmax(totalEnergy,1.0e-6);
    if ( fraction > _fractionEM) {
      
      float mass = 0.0;
      float energy = totalEnergy;
      
      float absPostition = 0.0;
      for (int i = 0; i < 3; ++i) absPostition += position[i]*position[i];
      absPostition = sqrt(absPostition);
      
      float momentum[3] = {0.0, 0.0, 0.0};
      if (totalEnergy > mass) {

	for (int i = 0; i < 3; ++i) momentum[i] = ( sqrt( pow(totalEnergy,2) - pow(mass,2) )/absPostition ) * position[i];

      }
      else {

	for (int i = 0; i < 3; ++i) momentum[i] = ( sqrt( pow(totalEnergy,2) )/absPostition ) * position[i];

      }

      recoParticle->setType(3);
      recoParticle->setCharge(0.0);
      recoParticle->setMass(mass);
      
      recoParticle->setMomentum(momentum);
      recoParticle->setEnergy(energy);
      
    }
    // else -> K0sl
    else {
  
      float mass = 0.49765;
      float energy = totalEnergy;
      
      float absPostition = 0.0;
      for (int i = 0; i < 3; ++i) absPostition += position[i]*position[i];
      absPostition = sqrt(absPostition);
      
      float momentum[3] = {0.0, 0.0, 0.0};
      if (totalEnergy > mass) {

	for (int i = 0; i < 3; ++i) momentum[i] = ( sqrt( pow(totalEnergy,2) - pow(mass,2) )/absPostition ) * position[i];

      }
      else {

	for (int i = 0; i < 3; ++i) momentum[i] = ( sqrt( pow(totalEnergy,2) )/absPostition ) * position[i];

      }      

      recoParticle->setType(4);
      recoParticle->setCharge(0.0);
      recoParticle->setMass(mass);
      
      recoParticle->setMomentum(momentum);
      recoParticle->setEnergy(energy);
            
    }
    
  }

}


bool TrackBasedPFlow::isRealMIPStub(const LCEvent* evt, Track* track, const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix) {


  std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes;

  try {
    
    const LCCollection* LCRcol = evt->getCollection(_colNameRelationTrackToMCP);
    LCRelationNavigator* nav = new LCRelationNavigator(LCRcol);

    const LCObjectVec& relMCParticlesToTrack = nav->getRelatedToObjects(track); 
    if ( relMCParticlesToTrack.size() > 1 ) std::cout << "Error: More than one MCParticle related to track." << std::endl;

    // FIXME: here is only the 0th contribution taken into account
    MCParticle* mcpOfTrack = dynamic_cast<MCParticle*>(relMCParticlesToTrack.at(0));
    

    std::vector< std::string >::const_iterator iter;
    const std::vector< std::string >* ColNames = evt->getCollectionNames();
    
    for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
	
      LCCollection* col = evt->getCollection( *iter ) ;

      if ( (col->getTypeName()) == LCIO::SIMCALORIMETERHIT ) {

	int nHits = col->getNumberOfElements();

	for(int k = 0; k < nHits; ++k) {
	  
	  SimCalorimeterHit* simCaloHit = dynamic_cast<SimCalorimeterHit*>(col->getElementAt(k));

	  for ( int l = 0; l < simCaloHit->getNMCContributions(); ++l ) {
	    
	    MCParticle* mcpOfCalo = simCaloHit->getParticleCont(l);


	    if ( mcpOfCalo == mcpOfTrack ) {
 	    
	      CalorimeterHitImpl* caloHitImpl = new CalorimeterHitImpl();

	      caloHitImpl->setEnergy(simCaloHit->getEnergy());

	      float pos[3];
	      for ( int m = 0; m < 3; ++m ) pos[m] = simCaloHit->getPosition()[m];	      
	      caloHitImpl->setPosition(pos);
	      
	      CalorimeterHit* caloHit = static_cast<CalorimeterHit*>(caloHitImpl);


	      LCVector3D referencePosition(outermostTrackerHits.at(0)->getPosition()[0],
					   outermostTrackerHits.at(0)->getPosition()[1],
					   outermostTrackerHits.at(0)->getPosition()[2]);
	      LCVector3D referencePositionProjected = getProjectedPointOnHelix(referencePosition,fittedHelix);

	      LCVector3D calorimeterHitPosition(caloHit->getPosition()[0],caloHit->getPosition()[1],caloHit->getPosition()[2]);
	      LCVector3D calorimeterHitPositionProjected = getProjectedPointOnHelix(calorimeterHitPosition,fittedHelix);
	      
	      double pathLengthOnHelix = getPathLengthOnHelix(referencePositionProjected,calorimeterHitPositionProjected,fittedHelix);
	      double distance = getDistanceToHelix(calorimeterHitPosition,fittedHelix);


	      CalorimeterHitWithAttributes* calorimeterHitWithAttributes = new CalorimeterHitWithAttributes(caloHit,distance,pathLengthOnHelix);
	       
	      bool calorimeterHitInserted = false;
	    
	      for(std::vector<CalorimeterHitWithAttributes*>::iterator m = calorimeterHitsWithAttributes.begin(); m !=  calorimeterHitsWithAttributes.end(); ++m){
		
		if ( pathLengthOnHelix > (*m)->getPathLengthOnHelix() ) {
		
		  calorimeterHitsWithAttributes.insert(m,calorimeterHitWithAttributes);
		  calorimeterHitInserted = true;
		  break;
		  
		}
		
	      }
	      
	      if (!calorimeterHitInserted) calorimeterHitsWithAttributes.push_back(calorimeterHitWithAttributes);
	      
	    
	    
	    }
	  
	  }

	}

      }
      
    }

    delete nav;
    nav = 0;
    
  }  
  catch(DataNotAvailableException &e){
    std::cout << "Collection " << _colNameRelationTrackToMCP << " not available in event." << std::endl;
  };


  // hits with the smallest path length on helix first    
  reverse(calorimeterHitsWithAttributes.begin(), calorimeterHitsWithAttributes.end());
  

  bool returnValue = false;
  bool caloHitBeforeWithinInnerTube = true;
  int nOfSuccessiveCaloHitsInTube = 0;
  int limitNOfSuccessiveCaloHitsInTube = 5; // in the first implementation hard coded
 
  float rOfInnerTube = _maximalRadiusOfInnerTubeForMIPLikeStub.at(0);
  float rOfOuterTube = _minimalRadiusOfOuterTubeForMIPLikeStub.at(0);


  for (std::vector<CalorimeterHitWithAttributes*>::const_iterator i = calorimeterHitsWithAttributes.begin(); i!=calorimeterHitsWithAttributes.end(); ++i) {

    if ( ( (*i)->getDistanceToHelix() > rOfOuterTube ) && caloHitBeforeWithinInnerTube ) {

      //debug
      if (_drawOnCED){
	ced_hit ( (*i)->getCalorimeterHit()->getPosition()[0],(*i)->getCalorimeterHit()->getPosition()[1],(*i)->getCalorimeterHit()->getPosition()[2], 
		  0 | 1 << CED_LAYER_SHIFT, 3, 0xa4ff2d);
      }

      continue;

    }
        
    //debug
    if (_drawOnCED){
      ced_hit ( (*i)->getCalorimeterHit()->getPosition()[0],(*i)->getCalorimeterHit()->getPosition()[1],(*i)->getCalorimeterHit()->getPosition()[2], 
		0 | 1 << CED_LAYER_SHIFT, 3, 0xff1842);
    }


    if ( ( (*i)->getDistanceToHelix() <= rOfInnerTube ) && caloHitBeforeWithinInnerTube ) ++nOfSuccessiveCaloHitsInTube;
    else caloHitBeforeWithinInnerTube = false;

    if ( !caloHitBeforeWithinInnerTube ) break;

    if ( nOfSuccessiveCaloHitsInTube >= limitNOfSuccessiveCaloHitsInTube ) {

      returnValue = true;
      break;
      
    }

  }



  for (std::vector<CalorimeterHitWithAttributes*>::iterator i = calorimeterHitsWithAttributes.begin(); i!=calorimeterHitsWithAttributes.end(); ++i) {

    CalorimeterHit* caloHit = (*i)->getCalorimeterHit();
    
    delete caloHit;
    caloHit = 0;
    
    delete (*i);
    (*i) = 0;
   
  }

  return returnValue;

}



// FIXME: put this as a function in MarlinUtil (see processor 'CutOnGEANT4Bug' as well)
void TrackBasedPFlow::getRelatedCalorimeterHitsPerfectly(const LCEvent* evt, Track* track, ClusterImpl* clusterRealEnergy, ClusterImpl* clusterPerfectEnergy) {

  
  // calibration of the calorimeter
  std::vector<float> calibration;
  for ( unsigned int i = 0; i < _calibrCoeffECAL.size(); ++i ) calibration.push_back(_calibrCoeffECAL.at(i));
  for ( unsigned int i = 0; i < _calibrCoeffHCAL.size(); ++i ) calibration.push_back(_calibrCoeffHCAL.at(i));
  

  try {

    // set up relations for MC particle to track
    LCCollection* LCRcolTracks = evt->getCollection(_colNameRelationTrackToMCP);

    LCRelationNavigator* navTracks = new LCRelationNavigator(LCRcolTracks);
    const LCObjectVec& relMCParticlesToTrack = navTracks->getRelatedToObjects(track); 

    if ( relMCParticlesToTrack.size() > 1 ) std::cout << "Warning: More than one MCParticle related to track." << std::endl;
    
    // set up relations for CalorimeterHit to SimCalorimeterHit
    LCCollection* LCRcolCalorimeter = evt->getCollection(_colNameRelationCaloHitToSimCaloHit);    
    LCRelationNavigator* navCalorimeter = new LCRelationNavigator(LCRcolCalorimeter);

    MCParticle* mcpOfTrack = 0;

    // container for the CalorimeterHits which are related to the track with hit and sub-hit energy accuracy, both run in parallel, i.e. they have the same size
    std::vector< std::pair<CalorimeterHit*,float> > collectedCalorimeterHitsWithEnergies;
    std::vector< std::pair<CalorimeterHit*,float> > collectedSubCalorimeterHitsWithEnergies;

    unsigned int index = 0;
    bool alreadyCollected = false;	

    float ESumCalorimeterHits = 0.0; // accumulated energy of calorimeter hits where the track contributes
    float ESumSubCalorimeterHits = 0.0; // accumulated energy of calorimeter hit energies, but only the part which originates from the track,i.e. sub-hit accuracy


    for(unsigned int i = 0; i < relMCParticlesToTrack.size(); ++i) {

      mcpOfTrack = dynamic_cast<MCParticle*>(relMCParticlesToTrack.at(i)); 

      MCParticleVec allMCPsOfTrack = MarlinUtil::getAllMCDaughters(mcpOfTrack);

      for(unsigned int j = 0; j < allMCPsOfTrack.size(); ++j) {
      
	MCParticle* mcp = allMCPsOfTrack.at(j);
	
	std::vector< std::string >::const_iterator iter;
	const std::vector< std::string >* ColNames = evt->getCollectionNames();
    
	for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
	
	  LCCollection* col = evt->getCollection( *iter ) ;
	
	  if ( (col->getTypeName()) == LCIO::SIMCALORIMETERHIT ) {
	  
	    int nHits = col->getNumberOfElements();	 

	    for(int k = 0; k < nHits; ++k) {
	    
	      SimCalorimeterHit* simCaloHit = dynamic_cast<SimCalorimeterHit*>(col->getElementAt(k));
	    
	      // debug
	      // std::cout << "simCaloHit->getNMCContributions() = " << simCaloHit->getNMCContributions() << std::endl;

	      for ( int l = 0; l < simCaloHit->getNMCContributions(); ++l ) {

		MCParticle* mcpOfCalo = simCaloHit->getParticleCont(l);
	      
		if ( mcpOfCalo == mcp ) {
		  
		  // float ESimHit = simCaloHit->getEnergy(); // only for debugging
		  float ESimHitContribution = simCaloHit->getEnergyCont(l);

		  const LCObjectVec& relCaloHitsToSimCaloHit = navCalorimeter->getRelatedFromObjects(simCaloHit);

		  // there should only be one CalorimeterHit related to one SimCalorimeterHit since the CalorimeterHit consists (is related to) of several SimCalorimeterHit, 
		  // but the SimCalorimeterHit is only related to one CalorimeterHit (by the ganging in the Calorimeter digitize processor)		
		  if (relCaloHitsToSimCaloHit.size() > 1 ) {
		  
		    std::cout << "Warning: More than one (" << relCaloHitsToSimCaloHit.size() << ") CalorimeterHit related to one SimCalorimeterHit. " << std::endl;
		  
		  }


		  // debug
		  // std::cout << "relCaloHitsToSimCaloHit.size() = " << relCaloHitsToSimCaloHit.size() << std::endl;


		  for ( unsigned int m = 0; m < relCaloHitsToSimCaloHit.size(); ++m ) {
		  
		    CalorimeterHit* caloHit = dynamic_cast<CalorimeterHit*>(relCaloHitsToSimCaloHit.at(m)); 	        		    

		    int type   = caloHit->getType();
		    float EHit = caloHit->getEnergy();
		    // float EHitCalculated   = ESimHit*calibration.at(type); // only for debugging
		    float EHitContribution = ESimHitContribution*calibration.at(type);
		  

		    // search if caloHit has already been assigned to track	  
		    for (unsigned int n = 0; n < collectedCalorimeterHitsWithEnergies.size(); ++n) {
		      
		      if (collectedCalorimeterHitsWithEnergies.at(n).first == caloHit) {
		      
			index = n;
			alreadyCollected = true;
			break;

		      }
		      else {

			index = 0;
			alreadyCollected = false;

		      }

		    
		    }
		  
		    
		    if ( !alreadyCollected ) { 		    
		  		   			    	    
		      std::pair<CalorimeterHit*,float> calorimeterHitWithEnergy(caloHit,EHit);
		      collectedCalorimeterHitsWithEnergies.push_back(calorimeterHitWithEnergy);
		      
		      std::pair<CalorimeterHit*,float> calorimeterSubHitWithEnergy(caloHit,EHitContribution);
		      collectedSubCalorimeterHitsWithEnergies.push_back(calorimeterSubHitWithEnergy);
		      
		    }
		    // calorimeter hit has already been assigned
		    else {
		    
		      collectedSubCalorimeterHitsWithEnergies.at(index).second += EHitContribution;
		    
		    }
		  

		    // debug
		    /*
		    std::cout << "Related hit to track " << track << "  " << "MCP of track = " << MarlinUtil::getMCName(mcp->getPDG()) 
			      << " ( " << mcp->getPDG() << " )" << std::endl
			      << "SimCaloHit " << simCaloHit << "  "
			      << "Pos SimCaloHit = " << "( " << simCaloHit->getPosition()[0] << "," << simCaloHit->getPosition()[1] << "," << simCaloHit->getPosition()[2] 
			      << " )" << "  " << "CaloHit " << caloHit << "  "
			      << "Pos CaloHit = " << "( " << caloHit->getPosition()[0] << "," << caloHit->getPosition()[1] << "," << caloHit->getPosition()[2] << " )" 
			      <<  std::endl
			      << "E CaloHit (full hit) = " << EHit << "  " << "CaloHit type = " << type << "  " << "E SimCaloHit = " << ESimHit << "  " << "EHitCalculated = " 
			      << EHitCalculated << "  " << "EHitContributed = " << EHitContribution << std::endl << std::endl;		    
		    */
		  
		  
		  }
	      
		}
		
	      }

	    }
	      
	  }
	  
	}
	
      }
      
    }
    

    


    // fill cluster and determine ESumCalorimeterHits and ESumSubCalorimeterHits
    
    //debug
    /*
    std::cout << "collectedCalorimeterHitsWithEnergies.size() = " << collectedCalorimeterHitsWithEnergies.size() << "  " 
	      << "collectedSubCalorimeterHitsWithEnergies.size() = " << collectedSubCalorimeterHitsWithEnergies.size() << std::endl;
    */

    for(unsigned int i = 0; i < collectedCalorimeterHitsWithEnergies.size(); ++i) {
      
      CalorimeterHit* caloHit = collectedCalorimeterHitsWithEnergies.at(i).first;
      float EHit = collectedCalorimeterHitsWithEnergies.at(i).second;
      
      CalorimeterHit* caloSubHit = collectedSubCalorimeterHitsWithEnergies.at(i).first;
      float ESubHit = collectedSubCalorimeterHitsWithEnergies.at(i).second;

      
      ESumCalorimeterHits += EHit;
      ESumSubCalorimeterHits += ESubHit;

      clusterRealEnergy->addHit(caloHit,EHit);
      clusterPerfectEnergy->addHit(caloSubHit,ESubHit);

      // debug
      // std::cout << "caloHit = " << caloHit << "  " << "EHit = " << EHit << "  " << "caloSubHit = " << caloSubHit << "  " << "ESubHit = " << ESubHit << std::endl;

    }


    // debug
    if ( _debugLevel > 5 ) { 
     
      std::cout << "MCP = " << MarlinUtil::getMCName(mcpOfTrack->getPDG()) << " ( " << mcpOfTrack->getPDG() << " )" << "  " << mcpOfTrack << "  " 
		<< "Generator Status: " << mcpOfTrack->getGeneratorStatus() << std::endl
		<< "E related to track (|p| = " << MarlinUtil::getAbsMomentum(track,_bField) << ") (hit accuracy) = " << ESumCalorimeterHits << "  " 
		<< "E related to track perfectly (sub-hit accuracy) = " << ESumSubCalorimeterHits << std::endl;	      
      
    }


    /*
    if ( abs(mcpOfTrack->getPDG()) == 13 ) _drawOnCED = 1;
    else _drawOnCED = 0;
    */




    // cut on GEANT4 bug
    // FIXME: should be done in a stand alone processor or somewhere else
    double eMCP = mcpOfTrack->getEnergy();
	    
    if ( ( eMCP > 5.0 ) && ( ESumSubCalorimeterHits > (1.75 * eMCP)) ) {
	      
      std::cout << std::endl
		<< "--------------------------------------------------------------------------------------------------------------------------------------" 
		<< std::endl << std::endl
		<< " ==> EVENT WITH GEANT4 BUG FOUND <=="
		<< std::endl << std::endl
		<< "     EVENT WILL BE DISCARDED ..."
		<< std::endl << std::endl
		<< "--------------------------------------------------------------------------------------------------------------------------------------" 
		<< std::endl << std::endl;
	      
      throw SkipEventException(this);

    }

    delete navCalorimeter;
    navCalorimeter = 0;   

    delete navTracks;
    navTracks = 0;
		
  }
  catch(DataNotAvailableException &e){
    std::cout << "Collection " << _colNameRelationTrackToMCP << " or " << _colNameRelationCaloHitToSimCaloHit  << " not available in event." << std::endl;
  };

}



void TrackBasedPFlow::drawRelatedCalorimeterHits(const LCEvent* evt,Track* track) {
  
  ClusterImpl* clusterRealEnergy = new ClusterImpl();
  ClusterImpl* clusterPerfectEnergy = new ClusterImpl();
	
 // also used to cut off GEANT4 bugs
  // FIXME: to take out GEANT4 bugs in this way takes far to much time, therefore it is done only for _drawOnCED; find an other way
  getRelatedCalorimeterHitsPerfectly(evt,track,clusterRealEnergy,clusterPerfectEnergy);

  if (_drawOnCED){ 

    for (CalorimeterHitVec::const_iterator iter = clusterRealEnergy->getCalorimeterHits().begin(); iter !=  clusterRealEnergy->getCalorimeterHits().end(); ++iter) {
      
      ced_hit( (*iter)->getPosition()[0], (*iter)->getPosition()[1], (*iter)->getPosition()[2], 0 | 8 << CED_LAYER_SHIFT, 4, 0xf7ff12 );
      
    }

  }

  delete clusterRealEnergy;
  clusterRealEnergy = 0;
  
  delete clusterPerfectEnergy;
  clusterPerfectEnergy = 0;
  



}


void TrackBasedPFlow::drawEMShowerCandidates(const LCEvent* evt) {


  try {
    
    std::vector< std::string >::const_iterator iter;
    const std::vector< std::string >* ColNames = evt->getCollectionNames();

    for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
      
      LCCollection* col = evt->getCollection( *iter ) ;
      
      if ( (col->getTypeName() == LCIO::CLUSTER) && (*iter == _colNameEMShowerCandidates) ) {

	int NCores = col->getNumberOfElements();
	
	for(int i=0; i<NCores; ++i){
	  
	  int color = 0x8f57ff*(i+32);

	  Cluster* cluster = dynamic_cast<Cluster*>(col->getElementAt(i));
	  
	  MarlinCED::drawCluster(cluster,2,5,color,6); 

	}

      }

    }
    
  }  
  catch(DataNotAvailableException &e){std::cout << "no valid collection of EM cores in event " << _nEvt << std::endl; };

}


std::vector<ClusterImplWithAttributes*> TrackBasedPFlow::assignClusters(ClusterImplWithAttributes* cluster, std::vector<ClusterImplWithAttributes*> clusters,
									LCVector3D referencePosition, const TrackerHitVec outermostTrackerHits, Trajectory* fittedHelix) {


  std::vector<ClusterImplWithAttributes*> result;
  result.push_back(cluster);

  // debug
  if ( _debugLevel > 7 ) { 

    std::cout << "cluster to start with: " << cluster << "  " << "( E = " << cluster->getEnergy() << " n = " << cluster->getClusterImpl()->getCalorimeterHits().size() 
	      << " )" << std::endl;
    std::cout << "list of clusters: " << std::endl;
    for(std::vector<ClusterImplWithAttributes*>::const_iterator i = clusters.begin(); i != clusters.end(); ++i) {

      std::cout << (*i) << "  " << "( E = " << (*i)->getEnergy() << " n = " << (*i)->getClusterImpl()->getCalorimeterHits().size() << " )" << std::endl;

    }

    std::cout << "resulting cluster list: " << std::endl;
    for(std::vector<ClusterImplWithAttributes*>::const_iterator i = result.begin(); i != result.end(); ++i) {

      std::cout << (*i) << "  " << "( E = " << (*i)->getEnergy() << " n = " << (*i)->getClusterImpl()->getCalorimeterHits().size() << " )" << std::endl;

    }

  }


  std::vector<ClusterImplWithAttributes*>::iterator iter = find(clusters.begin(),clusters.end(),cluster); 

  if ( iter != clusters.end() ) {
    
    // debug
    if ( _debugLevel > 7 ) std::cout << "start cluster in list of clusters to compare with => ERASE IT" << std::endl;
    
    clusters.erase(iter);
    
  }

  // debug
  if ( _debugLevel > 7 ) { 

    std::cout << "list of clusters used for comparission: " << std::endl;
    for(std::vector<ClusterImplWithAttributes*>::const_iterator i = clusters.begin(); i != clusters.end(); ++i) {

      std::cout << (*i) << "  " << "( E = " << (*i)->getEnergy() << " n = " << (*i)->getClusterImpl()->getCalorimeterHits().size() << " )" << std::endl;

    }

  }


  double maximalDistanceToCompareStartHitWith = 0.0;
  double maximalDistanceToCompareEndHitWith = 0.0;
  double maximalDistanceToCompareCoGWith = 0.0;
  int typeOfCoG = 0;

  for(std::vector<ClusterImplWithAttributes*>::const_iterator i = clusters.begin(); i != clusters.end(); ++i) {    
        
    if ( (cluster->getTypeEndHit()) == ((*i)->getTypeStartHit()) ) {
    

      switch (cluster->getTypeStartHit()) {

        case -1 : maximalDistanceToCompareStartHitWith = _maximalDistanceToHelixToAssignCluster.at(0); break;
        case  0 : maximalDistanceToCompareStartHitWith = _maximalDistanceToHelixToAssignCluster.at(1); break;
        case  1 : maximalDistanceToCompareStartHitWith = _maximalDistanceToHelixToAssignCluster.at(2); break;
        case  2 : maximalDistanceToCompareStartHitWith = _maximalDistanceToHelixToAssignCluster.at(3); break;
        default : maximalDistanceToCompareStartHitWith = _maximalDistanceToHelixToAssignCluster.at(0); break;
	  
      }

      switch (cluster->getTypeEndHit()) {

        case -1 : maximalDistanceToCompareEndHitWith = _maximalDistanceToHelixToAssignCluster.at(0); break;
        case  0 : maximalDistanceToCompareEndHitWith = _maximalDistanceToHelixToAssignCluster.at(1); break;
        case  1 : maximalDistanceToCompareEndHitWith = _maximalDistanceToHelixToAssignCluster.at(2); break;
        case  2 : maximalDistanceToCompareEndHitWith = _maximalDistanceToHelixToAssignCluster.at(3); break;
        default : maximalDistanceToCompareEndHitWith = _maximalDistanceToHelixToAssignCluster.at(0); break;
	  
      }
      
    } 

    else if ( (cluster->getTypeEndHit()) < ((*i)->getTypeStartHit()) ) {

      maximalDistanceToCompareStartHitWith = _maximalDistanceToHelixToAssignCluster.at( (cluster->getTypeStartHit()) + 1);
      maximalDistanceToCompareEndHitWith = _maximalDistanceToHelixToAssignCluster.at( (cluster->getTypeEndHit()) + 1);

    }
    else {
      
      maximalDistanceToCompareStartHitWith = _maximalDistanceToHelixToAssignCluster.at( (cluster->getTypeStartHit()) + 1);
      maximalDistanceToCompareEndHitWith = _maximalDistanceToHelixToAssignCluster.at( (cluster->getTypeEndHit()) + 1);

    }

    
    typeOfCoG = getTypeOfPositionOfCluster(cluster->getClusterImpl());

    // std::cout << "typeOfCoG = " << typeOfCoG << std::endl;

    maximalDistanceToCompareCoGWith = _maximalDistanceToHelixToAssignCluster.at(typeOfCoG+1);
    



    
    double distanceEndToStartHit = sqrt( pow( ( (*i)->getPositionStartHit()[0] - cluster->getPositionEndHit()[0] ) ,2) + 
					 pow( ( (*i)->getPositionStartHit()[1] - cluster->getPositionEndHit()[1] ) ,2) + 
					 pow( ( (*i)->getPositionStartHit()[2] - cluster->getPositionEndHit()[2] ) ,2) );
    
    double distanceStartToEndHit = sqrt( pow( ( (*i)->getPositionEndHit()[0] - cluster->getPositionStartHit()[0] ) ,2) + 
					 pow( ( (*i)->getPositionEndHit()[1] - cluster->getPositionStartHit()[1] ) ,2) + 
					 pow( ( (*i)->getPositionEndHit()[2] - cluster->getPositionStartHit()[2] ) ,2) );



    double distanceCoGToCoG =  sqrt( pow( ( (*i)->getPosition()[0] - cluster->getPosition()[0] ) ,2) + 
				     pow( ( (*i)->getPosition()[1] - cluster->getPosition()[1] ) ,2) + 
				     pow( ( (*i)->getPosition()[2] - cluster->getPosition()[2] ) ,2) );



    double distanceEndToCoG = sqrt( pow( ( (*i)->getPosition()[0] - cluster->getPositionEndHit()[0] ) ,2) + 
				    pow( ( (*i)->getPosition()[1] - cluster->getPositionEndHit()[1] ) ,2) + 
				    pow( ( (*i)->getPosition()[2] - cluster->getPositionEndHit()[2] ) ,2) );

    double distanceStartToCoG = sqrt( pow( ( (*i)->getPosition()[0] - cluster->getPositionStartHit()[0] ) ,2) + 
				      pow( ( (*i)->getPosition()[1] - cluster->getPositionStartHit()[1] ) ,2) + 
				      pow( ( (*i)->getPosition()[2] - cluster->getPositionStartHit()[2] ) ,2) );



    double distanceCoGToStartHit =  sqrt( pow( ( (*i)->getPositionStartHit()[0] - cluster->getPosition()[0] ) ,2) + 
					  pow( ( (*i)->getPositionStartHit()[1] - cluster->getPosition()[1] ) ,2) + 
					  pow( ( (*i)->getPositionStartHit()[2] - cluster->getPosition()[2] ) ,2) );

    double distanceCoGToEndHit =  sqrt( pow( ( (*i)->getPositionEndHit()[0] - cluster->getPosition()[0] ) ,2) + 
					pow( ( (*i)->getPositionEndHit()[1] - cluster->getPosition()[1] ) ,2) + 
					pow( ( (*i)->getPositionEndHit()[2] - cluster->getPosition()[2] ) ,2) );



    /*
    LCVector3D clusterEndPoint(cluster->getPositionEndHit()[0],cluster->getPositionEndHit()[1],cluster->getPositionEndHit()[2]);
    LCVector3D clusterEndPointProjected = getProjectedPointOnHelix(clusterEndPoint,fittedHelix);
    double pathLengthEndPointProjected  = getPathLengthOnHelix(clusterEndPointProjected,fittedHelix) - pathLengthOfReferencePositionProjected;
    */

    // get path length of reference point. This needs as the origin of the 'path-length coordinate system' the last tracker hit as s=0
    LCVector3D zeroPathLength(outermostTrackerHits.at(0)->getPosition()[0],outermostTrackerHits.at(0)->getPosition()[1],outermostTrackerHits.at(0)->getPosition()[2]);
    LCVector3D zeroPathLengthProjected = getProjectedPointOnHelix(zeroPathLength,fittedHelix);    
    double pathLengthOfZeroPathLengthProjected = getPathLengthOnHelix(zeroPathLengthProjected,fittedHelix);

    LCVector3D referencePositionProjected = getProjectedPointOnHelix(referencePosition,fittedHelix);
    double pathLengthOfReferencePositionProjected = getPathLengthOnHelix(referencePositionProjected,fittedHelix) - pathLengthOfZeroPathLengthProjected;

    LCVector3D positionCoG((*i)->getPosition()[0],(*i)->getPosition()[1],(*i)->getPosition()[2] );
    LCVector3D positionCoGProjected = getProjectedPointOnHelix(positionCoG,fittedHelix);

    double pathLengthRefPositionToCoG = getPathLengthOnHelix(referencePositionProjected,positionCoGProjected,fittedHelix);

    double distanceExtrapolatedTrajectoryToCoG = getDistanceToHelix(positionCoG,fittedHelix);


    // debug
    if ( _debugLevel > 7 ) { 

      std::cout << "cluster: " << cluster << "  " << "( E = " << cluster->getEnergy() << " n = " << cluster->getClusterImpl()->getCalorimeterHits().size() << " )" << "  " 
		<< "cluster to compare: " << (*i) << "  " << "( E = " << (*i)->getEnergy() << " n = " << (*i)->getClusterImpl()->getCalorimeterHits().size() << " )" << "  " 
		<< "conditions to assign cluster:" << std::endl
		<< "distanceEndToStartHit < maximalDistanceToCompareEndHitWith (" << distanceEndToStartHit << " < " << maximalDistanceToCompareEndHitWith << "): " 
		<< (distanceEndToStartHit < maximalDistanceToCompareEndHitWith) << "  " << std::endl
		<< "distanceStartToEndHit < maximalDistanceToCompareStartHitWith (" << distanceStartToEndHit << " < " << maximalDistanceToCompareStartHitWith << "): " 
		<< (distanceStartToEndHit < maximalDistanceToCompareStartHitWith) << "  " << std::endl 
		<< "distanceCoGToCoG < maximalDistanceToCompareCoGWith (" << distanceCoGToCoG << " < " << maximalDistanceToCompareCoGWith << "):  " 
		<< (distanceCoGToCoG < maximalDistanceToCompareCoGWith) << "  " << std::endl
		<< "distanceEndToCoG < maximalDistanceToCompareEndHitWith (" << distanceEndToCoG << " < " << maximalDistanceToCompareEndHitWith << "):  " 
		<< (distanceEndToCoG < maximalDistanceToCompareEndHitWith) << "  " << std::endl
		<< "distanceStartToCoG < maximalDistanceToCompareStartHitWith (" << distanceStartToCoG << " < " << maximalDistanceToCompareStartHitWith << "): " 
		<< (distanceStartToCoG < maximalDistanceToCompareStartHitWith) << "  " << std::endl	
		<< "distanceCoGToStartHit < maximalDistanceToCompareCoGWith (" << distanceCoGToStartHit << " < " << maximalDistanceToCompareCoGWith << "): "
		<< (distanceCoGToStartHit < maximalDistanceToCompareCoGWith) << "  " << std::endl
		<< "distanceCoGToEndHit < maximalDistanceToCompareCoGWith (" << distanceCoGToEndHit << " < " << maximalDistanceToCompareCoGWith << "): " 
		<< (distanceCoGToEndHit < maximalDistanceToCompareCoGWith) << "  " << std::endl
		<< "pathLengthRefPositionToCoG (" << pathLengthRefPositionToCoG << ") < ( _maximalConeTubeLength (" << _maximalConeTubeLength 
		<< ") - pathLengthRefPosition (" << pathLengthOfReferencePositionProjected << ") ) (" 
		<< _maximalConeTubeLength-pathLengthOfReferencePositionProjected << "): " 
		<< (pathLengthRefPositionToCoG < (_maximalConeTubeLength-pathLengthOfReferencePositionProjected)) << "  " << std::endl
		<< "pathLengthRefPositionToCoG (" << pathLengthRefPositionToCoG << ") > -0.1 * ( _maximalConeTubeLength (" << _maximalConeTubeLength 
		<< ") - pathLengthRefPosition (" << pathLengthOfReferencePositionProjected << ") ) (" 
		<< -0.1*(_maximalConeTubeLength-pathLengthOfReferencePositionProjected)<< "): " 
		<< (pathLengthRefPositionToCoG > -0.1*(_maximalConeTubeLength-pathLengthOfReferencePositionProjected))
		<< "  " << std::endl 
		<< "distanceExtrapolatedTrajectoryToCoG (" << distanceExtrapolatedTrajectoryToCoG <<") < maximalDistanceToCompareCoGWith (" 
		<< maximalDistanceToCompareCoGWith << "): " << (distanceExtrapolatedTrajectoryToCoG < maximalDistanceToCompareCoGWith ) << std::endl;

    }


    
    bool conditionToAssign = 
      ( distanceEndToStartHit < maximalDistanceToCompareEndHitWith ) ||  ( distanceStartToEndHit < maximalDistanceToCompareStartHitWith ) || 
      ( distanceCoGToCoG < maximalDistanceToCompareCoGWith ) ||
      ( distanceEndToCoG < maximalDistanceToCompareEndHitWith ) ||  ( distanceStartToCoG < maximalDistanceToCompareStartHitWith ) ||
      ( distanceCoGToStartHit < maximalDistanceToCompareCoGWith ) ||  ( distanceCoGToEndHit < maximalDistanceToCompareCoGWith ) ||
      ( ( pathLengthRefPositionToCoG < (_maximalConeTubeLength-pathLengthOfReferencePositionProjected) ) 
	&& ( pathLengthRefPositionToCoG > -0.1 * (_maximalConeTubeLength-pathLengthOfReferencePositionProjected)  ) 
	&& ( distanceExtrapolatedTrajectoryToCoG < maximalDistanceToCompareCoGWith ) );

    

    if (conditionToAssign) {


      std::vector<ClusterImplWithAttributes*>::const_iterator isInResult = find(result.begin(),result.end(),(*i));      
      
      if (isInResult == result.end() ) {
      
	std::vector<ClusterImplWithAttributes*> relatedClusters = assignClusters( (*i),clusters,referencePosition,outermostTrackerHits,fittedHelix );

	for(std::vector<ClusterImplWithAttributes*>::const_iterator j = relatedClusters.begin(); j != relatedClusters.end(); ++j) {
	  
	  std::vector<ClusterImplWithAttributes*>::const_iterator relatedClusterInResult = find(result.begin(),result.end(),(*j));
	  
	  if ( relatedClusterInResult == result.end() ) result.push_back( (*j) );

	}

      }

    }
    
  }
  

  // debug
  if ( _debugLevel > 7 ) { 

    std::cout << "result : " << std::endl;
    for(std::vector<ClusterImplWithAttributes*>::const_iterator i = result.begin(); i != result.end(); ++i) std::cout << (*i) << "  " << (*i)->getEnergy() << std::endl;

  }  


  return result;
  
}




std::vector<ClusterImplWithAttributes*> TrackBasedPFlow::assignAdditionalClusters(ClusterImplWithAttributes* mipStub, Trajectory* fittedHelix,
										  std::vector<ClusterImplWithAttributes*> clustersAlreadyAssigned,
										  std::vector<ClusterImplWithAttributes*> clustersToCheck, 
										  Track* track) {
  


  std::vector<ClusterImplWithAttributes*> clustersToAssignAdditionally;



  // FIXME: this is only a simple approach to get the energy resolution of the assigned energy, we need a much more elaborated procedure here

  
  // get path lenght of end point of MIP stub as a reference
  
  LCVector3D mipStubPosition(mipStub->getPosition()[0],mipStub->getPosition()[1],mipStub->getPosition()[2]);
  LCVector3D mipStubPositionProjected = getProjectedPointOnHelix(mipStubPosition,fittedHelix);
  double pathLengthMIPStubPositionProjected = getPathLengthOnHelix(mipStubPositionProjected,fittedHelix);
  


  // get track energy for comparison
  double absTrackMomentum = MarlinUtil::getAbsMomentum(track,_bField);
  double energyOfTrack = sqrt( absTrackMomentum*absTrackMomentum + 0.13957*0.13957); // FIXME: at the moment the pion hypothesis is used
  
  // calculate already assigned energy
  double energyOfAssignedClusters = 0.0;
  for(std::vector<ClusterImplWithAttributes*>::const_iterator i = clustersAlreadyAssigned.begin(); i != clustersAlreadyAssigned.end(); ++i) {    
    
    for(CalorimeterHitVec::const_iterator j = (*i)->getClusterImpl()->getCalorimeterHits().begin(); j != (*i)->getClusterImpl()->getCalorimeterHits().end(); ++j) {    

      energyOfAssignedClusters += (*j)->getEnergy();
     
    }

  }
  // use loose cuts if the MIP stub is the only cluster which has been assigned
  bool mipStubIsTheOnlyCluster = false;
  bool mipStubHasBeenTheOnlyCluster = false;

  if ( (clustersAlreadyAssigned.size() == 1) && ( clustersAlreadyAssigned.at(0)->isMIPStub()) ) mipStubIsTheOnlyCluster = true;



  double maximalDistanceToCompareWith = 0.0;
  
  for(std::vector<ClusterImplWithAttributes*>::const_iterator i = clustersToCheck.begin(); i != clustersToCheck.end(); ++i) { 
    
    double smallestDistanceBetweenClusters = DBL_MAX;
    double energyOfClusterToCheck = 0.0;

    LCVector3D startPositionOfClusterToCheck((*i)->getPosition()[0],(*i)->getPosition()[1],(*i)->getPosition()[2]);
    LCVector3D startPositionOfClusterToCheckProjected = getProjectedPointOnHelix(startPositionOfClusterToCheck,fittedHelix);
    double pathLengthOfStartPositionProjected = getPathLengthOnHelix(startPositionOfClusterToCheckProjected,fittedHelix);

  
    // debug
    if ( _debugLevel > 5 ) { 
      std::cout << "pathLengthMIPStubPositionProjected: " << pathLengthMIPStubPositionProjected << "  " 
		<< "pathLengthOfStartPositionProjected: " << pathLengthOfStartPositionProjected << "  " 
		<< ( pathLengthOfStartPositionProjected < pathLengthMIPStubPositionProjected );
      
      if ( pathLengthOfStartPositionProjected < pathLengthMIPStubPositionProjected ) std::cout << "  " << "s is too small => cluster is skipped";
      
      std::cout << std::endl;
    }


    // have a look at the next cluster if path length of start position of the cluster to check is smaller than path length of position of MIP stub    
    if ( pathLengthOfStartPositionProjected < pathLengthMIPStubPositionProjected ) continue;



    for(CalorimeterHitVec::const_iterator j = (*i)->getClusterImpl()->getCalorimeterHits().begin(); j != (*i)->getClusterImpl()->getCalorimeterHits().end(); ++j) { 

      energyOfClusterToCheck += (*j)->getEnergy();
      
    }

    ClusterImpl* clusterToCheck = (*i)->getClusterImpl();

    
    // assign additional clusters if already at least one cluster has been assigned or if no cluster has been assigned search with loose cuts (see 'else' case)
    if ( (energyOfAssignedClusters > 0) && (!mipStubIsTheOnlyCluster) ) {
      
      for(CalorimeterHitVec::const_iterator j = clusterToCheck->getCalorimeterHits().begin(); j != clusterToCheck->getCalorimeterHits().end(); ++j) {    
	
	float positionCaloHitToCheck[3];    
	for (int k = 0; k < 3; ++k) positionCaloHitToCheck[k] = (*j)->getPosition()[k];

	std::vector<ClusterImplWithAttributes*> clustersAssignedToTrack;

	if ( !mipStubHasBeenTheOnlyCluster ) clustersAssignedToTrack = clustersAlreadyAssigned;
	else clustersAssignedToTrack = clustersToAssignAdditionally;
	  
	for(std::vector<ClusterImplWithAttributes*>::const_iterator k = clustersAssignedToTrack.begin(); k != clustersAssignedToTrack.end(); ++k) { 
	  
	  ClusterImpl* assignedCluster = (*k)->getClusterImpl();
	  
	  for(CalorimeterHitVec::const_iterator l = assignedCluster->getCalorimeterHits().begin(); l != assignedCluster->getCalorimeterHits().end(); ++l) {    
	    
	    float positionCaloHitAlreadyAssigned[3];
	    
	    for (int m = 0; m < 3; ++m) positionCaloHitAlreadyAssigned[m] = (*l)->getPosition()[m];
	    
	      
	    double distance = sqrt( pow( ( positionCaloHitAlreadyAssigned[0] - positionCaloHitToCheck[0] ) ,2) + 
				    pow( ( positionCaloHitAlreadyAssigned[1] - positionCaloHitToCheck[1] ) ,2) + 
				    pow( ( positionCaloHitAlreadyAssigned[2] - positionCaloHitToCheck[2] ) ,2) );
	    
	    
	    if ( distance < smallestDistanceBetweenClusters ) { 
	      
	      /*
	      // check if hit of 'cluster to check' is isolated
	      double smallestDistanceToNearestHit = DBL_MAX;
	      double maximalDistanceToCompareNearestHitWith = 0.0;
	      
	      for(CalorimeterHitVec::const_iterator m = clusterToCheck->getCalorimeterHits().begin(); m != clusterToCheck->getCalorimeterHits().end(); ++m) {
	      
	      float positionCaloHitNearestBy[3];    
	      for (int q = 0; q < 3; ++q) positionCaloHitNearestBy[q] = (*m)->getPosition()[q];
	      
	      double distanceToHit = sqrt( pow( ( positionCaloHitNearestBy[0] - positionCaloHitToCheck[0] ) ,2) + 
	      pow( ( positionCaloHitNearestBy[1] - positionCaloHitToCheck[1] ) ,2) + 
	      pow( ( positionCaloHitNearestBy[2] - positionCaloHitToCheck[2] ) ,2) );
	      
	      
	      if ( ( distanceToHit < smallestDistanceToNearestHit ) && ( m != j ) ) {
	      
	      smallestDistanceToNearestHit = distanceToHit;
	      
	      int typeOfHitToCheck = (*j)->getType();
	      
	      switch ( typeOfHitToCheck ) {
		    
	      case -1 : maximalDistanceToCompareNearestHitWith = _maximalDistanceToHelixToAssignCluster.at(0); break;
	      case  0 : maximalDistanceToCompareNearestHitWith = _maximalDistanceToHelixToAssignCluster.at(1); break;
	      case  1 : maximalDistanceToCompareNearestHitWith = _maximalDistanceToHelixToAssignCluster.at(2); break;
	      case  2 : maximalDistanceToCompareNearestHitWith = _maximalDistanceToHelixToAssignCluster.at(3); break;
	      default : maximalDistanceToCompareNearestHitWith = _maximalDistanceToHelixToAssignCluster.at(0); break;
	      
	      } 		  		
	      
	      }
	      
	      }
	      */
	      
	      if ( true /*smallestDistanceToNearestHit < (maximalDistanceToCompareNearestHitWith/1.5)*/ ) {
		
		smallestDistanceBetweenClusters = distance;
		
		int typeOfHitAlreadyAssigned = (*l)->getType();
		
		switch ( typeOfHitAlreadyAssigned ) {
		  
		case -1 : maximalDistanceToCompareWith = _maximalDistanceToHelixToAssignCluster.at(0); break;
		case  0 : maximalDistanceToCompareWith = _maximalDistanceToHelixToAssignCluster.at(1); break;
		case  1 : maximalDistanceToCompareWith = _maximalDistanceToHelixToAssignCluster.at(2); break;
		case  2 : maximalDistanceToCompareWith = _maximalDistanceToHelixToAssignCluster.at(3); break;
		default : maximalDistanceToCompareWith = _maximalDistanceToHelixToAssignCluster.at(0); break;
		  
		}
		
	      }
	      
	    }
	    
	  }
	  
	}
	               
      }
    
    }
    else {

      for(CalorimeterHitVec::const_iterator j = clusterToCheck->getCalorimeterHits().begin(); j != clusterToCheck->getCalorimeterHits().end(); ++j) {    
	
	float positionCaloHitToCheck[3];    
	for (int k = 0; k < 3; ++k) positionCaloHitToCheck[k] = (*j)->getPosition()[k];

	double distance = sqrt( pow( ( mipStub->getPositionEndHit()[0] - positionCaloHitToCheck[0] ) ,2) + 
				pow( ( mipStub->getPositionEndHit()[1] - positionCaloHitToCheck[1] ) ,2) + 
				pow( ( mipStub->getPositionEndHit()[2] - positionCaloHitToCheck[2] ) ,2) );

	if ( distance < smallestDistanceBetweenClusters ) smallestDistanceBetweenClusters = distance; 


      }

      maximalDistanceToCompareWith = _maximalDistanceToHelixToAssignCluster.at(0); // largest distance for assignment

    }
    



    //debug
    if ( _debugLevel > 5 ) { 

      if ( (energyOfAssignedClusters > 0)  && (!mipStubIsTheOnlyCluster) ) {

	std::cout << "Try to assign additional Clusters: " << std::endl
		  << "energyOfClusterToCheck = " << energyOfClusterToCheck << "  " << "0.5*energyOfAssignedClusters = " << 0.5*energyOfAssignedClusters << "  "
		  << (  energyOfClusterToCheck < 0.5*energyOfAssignedClusters ) << std::endl
		  << "energyOfAssignedClusters = " << energyOfAssignedClusters << "  " << "energyOfClusterToCheck = " << energyOfClusterToCheck << "  " 
		  << "energyOfTrack = " << energyOfTrack << std::endl 
		  << "(energyOfAssignedClusters+energyOfClusterToCheck-energyOfTrack) = " << (energyOfAssignedClusters+energyOfClusterToCheck-energyOfTrack) << "  " 
		  << "2.0*0.3*sqrt(energyOfTrack) = " << 2.0*0.3*sqrt(energyOfTrack) << "  " 
		  << ( (energyOfAssignedClusters+energyOfClusterToCheck-energyOfTrack) < 2.0*0.3*sqrt(energyOfTrack) ) << "  " << std::endl
		  << "smallestDistanceBetweenClusters = " << smallestDistanceBetweenClusters << "  " << "1.0*maximalDistanceToCompareWith = "
		  << 1.0*maximalDistanceToCompareWith << "  " << (smallestDistanceBetweenClusters < 1.0*maximalDistanceToCompareWith) << std::endl		
		  << "Full Condition: " << ( ( energyOfClusterToCheck < 0.5*energyOfAssignedClusters ) && 
					     ( (energyOfAssignedClusters+energyOfClusterToCheck-energyOfTrack) < 2.0*0.3*sqrt(energyOfTrack) ) && 
					     ( smallestDistanceBetweenClusters < 1.0*maximalDistanceToCompareWith ) ) << std::endl << std::endl;

      }
      else {

	std::cout << "Only a MIP stub or even no cluster has been assigned to track within the first step. Trying to find additional clusters around ... " << std::endl
		  << "Try to assign additional Clusters with loose cuts: " << std::endl
		  << "energyOfClusterToCheck = " << energyOfClusterToCheck << "  " << "energyOfTrack = " << energyOfTrack << std::endl
		  << "(energyOfClusterToCheck-energyOfTrack) = " << (energyOfClusterToCheck-energyOfTrack) << "  " 
		  << "3.0*0.3*sqrt(energyOfTrack) = " << 3.0*0.3*sqrt(energyOfTrack) << "  "
		  << ( (energyOfClusterToCheck-energyOfTrack) < 3.0*0.3*sqrt(energyOfTrack) ) << "  " << std::endl
		  << "smallestDistanceBetweenClusters = " << smallestDistanceBetweenClusters << "  " << "1.5*maximalDistanceToCompareWith = "
		  << 1.5*maximalDistanceToCompareWith << "  " << (smallestDistanceBetweenClusters < 1.5*maximalDistanceToCompareWith) << std::endl	
		  << "Full Condition: " << ( ( (energyOfClusterToCheck-energyOfTrack) < 3.0*0.3*sqrt(energyOfTrack) ) && 
					     ( smallestDistanceBetweenClusters < 1.5*maximalDistanceToCompareWith ) ) << std::endl << std::endl;
      }
      
    }
    

    if ( (energyOfAssignedClusters > 0)  && (!mipStubIsTheOnlyCluster) ) { 

      // FIXME: This is only a first, simple ansatz with hard coded factors
      if ( ( energyOfClusterToCheck < 0.5*energyOfAssignedClusters ) && 
	   ( (energyOfAssignedClusters+energyOfClusterToCheck-energyOfTrack) < 2.0*0.3*sqrt(energyOfTrack) ) && 
	   (smallestDistanceBetweenClusters < 1.0*maximalDistanceToCompareWith) ) {
      
	clustersToAssignAdditionally.push_back(*i);     
	energyOfAssignedClusters += (*i)->getEnergy(); // FIXME: this is not the best way, it would be better to look at all permutations

      }
    
    }
    else {

      if ( ( (energyOfClusterToCheck-energyOfTrack) < 3.0*0.3*sqrt(energyOfTrack) ) && 
	   ( smallestDistanceBetweenClusters < 1.5*maximalDistanceToCompareWith ) ) {

	clustersToAssignAdditionally.push_back(*i);     
	energyOfAssignedClusters += (*i)->getEnergy(); // FIXME: this is not the best way, it would be better to look at all permutations
	mipStubIsTheOnlyCluster = false; // FIXME: this is not the best way, it would be better to look at all permutations
	mipStubHasBeenTheOnlyCluster = true;
      }

    }

  }
    

  return clustersToAssignAdditionally;

}


  
std::vector<CalorimeterHit*> TrackBasedPFlow::getNeutralHitsAssignedToChargedParticle(LCEvent* evt, ReconstructedParticle* recoParticle, int& n, double& energy,
										      double hitEnergyFraction) {



  std::vector<CalorimeterHit*> collectedCalotimeterHits;

  if (recoParticle->getTracks().size() == 0 ) {

    // debug
    std::cout << "ERROR: Charged reconstructed particle without track found in method 'getNeutralHitsAssignedToChargedParticle()'. Cannot calculate neutral energy." 
	      << std::endl << std::endl;

    return collectedCalotimeterHits;

  }
  else if (recoParticle->getTracks().size() > 1 ) {

    // debug
    std::cout << "WARNING: More than one track assigned to reconstructed particle found in method 'getNeutralHitsAssignedToChargedParticle()'." << std::endl 
	      << "Continue with the first track only" 
	      << std::endl << std::endl;

  }
  

  Track* track = recoParticle->getTracks().at(0);

  // set up relations for MC particle to track
  LCCollection* LCRcolTracks = evt->getCollection(_colNameRelationTrackToMCP);  
  LCRelationNavigator* navTracks = new LCRelationNavigator(LCRcolTracks);
  const LCObjectVec& relMCParticlesToTrack = navTracks->getRelatedToObjects(track); 
  
  if ( relMCParticlesToTrack.size() > 1 ) std::cout << "Warning: More than one MCParticle related to track." << std::endl;
    

  // set up relations for CalorimeterHit to SimCalorimeterHit
  LCCollection* LCRcolCalorimeter = evt->getCollection(_colNameRelationCaloHitToSimCaloHit);  
  LCRelationNavigator* navCalorimeter = new LCRelationNavigator(LCRcolCalorimeter);


  // calibration of the calorimeter
  std::vector<float> calibration;
  for ( unsigned int i = 0; i < _calibrCoeffECAL.size(); ++i ) calibration.push_back(_calibrCoeffECAL.at(i));
  for ( unsigned int i = 0; i < _calibrCoeffHCAL.size(); ++i ) calibration.push_back(_calibrCoeffHCAL.at(i));


  for(ClusterVec::const_iterator i = recoParticle->getClusters().begin(); i != recoParticle->getClusters().end(); ++i) {
    
    for(CalorimeterHitVec::const_iterator j = (*i)->getCalorimeterHits().begin(); j != (*i)->getCalorimeterHits().end(); ++j) {
      
      double energyCaloHit = (*j)->getEnergy();
      double relatedNeutralEnergyPerCaloHit = 0.0;
   
      const LCObjectVec& relCaloHitsToSimCaloHit = navCalorimeter->getRelatedToObjects(*j);
      
      for ( unsigned int k = 0; k < relCaloHitsToSimCaloHit.size(); ++k ) {
		  
	SimCalorimeterHit* simCaloHit = dynamic_cast<SimCalorimeterHit*>(relCaloHitsToSimCaloHit.at(k));
	
	for ( int l = 0; l < simCaloHit->getNMCContributions(); ++l ) {
	  
	  MCParticle* mcpOfSimCaloHit = simCaloHit->getParticleCont(l);
	  double eSimHitContribution = simCaloHit->getEnergyCont(l);

	  MCParticleVec parentsOfMCofSimCaloHit = MarlinUtil::getAllMCParents(mcpOfSimCaloHit);

	  // debug
	  if ( _debugLevel > 5 ) {

	    if ( parentsOfMCofSimCaloHit.size() > 1 ) {

	      std::cout << "Warning: More than one (" << parentsOfMCofSimCaloHit.size() << ") MC parent found for MCParticle " << mcpOfSimCaloHit 
			<< " (" << mcpOfSimCaloHit->getPDG() << ") " << " in SimCalorimeterHit " << simCaloHit 
			<< " in method 'getNeutralHitsAssignedToChargedParticle()'" << std::endl;
	      
		}

	  }
	  
	  for ( MCParticleVec::const_iterator m = parentsOfMCofSimCaloHit.begin(); m != parentsOfMCofSimCaloHit.end(); ++m ) {
	 	    
	    MCParticle* mcpOfTrack = 0;
	    
	    for(unsigned int n = 0; n < relMCParticlesToTrack.size(); ++n) {
	    
	      mcpOfTrack = dynamic_cast<MCParticle*>(relMCParticlesToTrack.at(n)); 

	      // debug
	      /*
	      std::cout << "mcpOfSimCaloHit: " << mcpOfSimCaloHit << "  " << mcpOfSimCaloHit->getVertex()[0] << "|" << mcpOfSimCaloHit->getVertex()[1] << "|" 
			<< mcpOfSimCaloHit->getVertex()[2] << std::endl;
	      */

	      
	      bool isMCOfTrackDaughterOfParentOfMCOfSimCaloHit = MarlinUtil::isDaughterOf( mcpOfTrack, (*m) );
	    
	      if ( (mcpOfTrack != (*m)) && !isMCOfTrackDaughterOfParentOfMCOfSimCaloHit ) { // SimCaloHit found which has no contribution of mcpOfTrack
		
		// take only neutral particles causing this SimCaloHit into account
		std::string charge = _mcParticleHelper->getMCCharge( (*m)->getPDG() );


		// debug
		/*
		std::cout << "CaloHit: " << (*j) << "  " << "E = " << (*j)->getEnergy() << "  " 
			  << "SimCaloHit: " << simCaloHit << "  " << "E = " << eSimHitContribution*calibration.at((*j)->getType()) << std::endl 
			  << "MCContrib: " << l+1 << "/" << simCaloHit->getNMCContributions() << "  " << mcpOfSimCaloHit << "  " << mcpOfSimCaloHit->getPDG() << "  " 
			  << mcpOfSimCaloHit->getEnergy() << "  " << "GeneratorStatus: " << mcpOfSimCaloHit->getGeneratorStatus() << "  " 
			  << "CreatedInSimulation: " << mcpOfSimCaloHit->isCreatedInSimulation() << "  " << "DecayedInTracker: " << mcpOfSimCaloHit->isDecayedInTracker() 
			  << std::endl
			  << "MCParents: " << ( m-parentsOfMCofSimCaloHit.begin()+1 ) << "/" << parentsOfMCofSimCaloHit.size() << "  " << (*m) << "  " << (*m)->getPDG() 
			  << "  " << (*m)->getEnergy() << "  " << charge << "  " << "GeneratorStatus: " << (*m)->getGeneratorStatus() << "  " 
			  << "CreatedInSimulation: " << (*m)->isCreatedInSimulation() << "  "
			  << "DecayedInTracker: " << (*m)->isDecayedInTracker() << std::endl << std::endl;
		*/

	 
		if (charge == "0") { // neutral MC contribution found
	    
		  int type   = (*j)->getType();	     
		  double eHitContribution = eSimHitContribution*calibration.at(type);
		  
		  relatedNeutralEnergyPerCaloHit += eHitContribution;
	      
		}
		
	      }
	      
	    }
	    
	  }	  
	  
	}
	
      }

      if ( relatedNeutralEnergyPerCaloHit > (hitEnergyFraction*energyCaloHit) ) {
      
	collectedCalotimeterHits.push_back(*j);
	++n;
	energy += energyCaloHit;

      }

    

    }
    
  }

  delete navCalorimeter;
  navCalorimeter = 0;

  delete navTracks;
  navTracks = 0;
  
  return collectedCalotimeterHits;
  
}


std::vector<CalorimeterHit*> TrackBasedPFlow::getChargedHitsAssignedToNeutralParticle(LCEvent* evt, ReconstructedParticle* recoParticle, int& n, double& energy,
										      double hitEnergyFraction) {
  

  std::vector<CalorimeterHit*> collectedCalotimeterHits;

  if (recoParticle->getTracks().size() > 0 ) {

    // debug
    std::cout << "ERROR: Neutral reconstructed particle with track found in method 'getChargedHitsAssignedToNeutralParticle()'. Cannot calculate charged energy." 
	      << std::endl << std::endl;

    return collectedCalotimeterHits;

  }

  // set up relations for MC particle to track
  LCCollection* LCRcolTracks = evt->getCollection(_colNameRelationTrackToMCP);
  LCRelationNavigator* navTracks = new LCRelationNavigator(LCRcolTracks);

  // set up relations for CalorimeterHit to SimCalorimeterHit
  LCCollection* LCRcolCalorimeter = evt->getCollection(_colNameRelationCaloHitToSimCaloHit);
  LCRelationNavigator* navCalorimeter = new LCRelationNavigator(LCRcolCalorimeter);


  // calibration of the calorimeter
  std::vector<float> calibration;
  for ( unsigned int i = 0; i < _calibrCoeffECAL.size(); ++i ) calibration.push_back(_calibrCoeffECAL.at(i));
  for ( unsigned int i = 0; i < _calibrCoeffHCAL.size(); ++i ) calibration.push_back(_calibrCoeffHCAL.at(i));


  for(ClusterVec::const_iterator i = recoParticle->getClusters().begin(); i != recoParticle->getClusters().end(); ++i) {

    for(CalorimeterHitVec::const_iterator j = (*i)->getCalorimeterHits().begin(); j != (*i)->getCalorimeterHits().end(); ++j) {
      
      double energyCaloHit = (*j)->getEnergy();
      double relatedChargedEnergyPerCaloHit = 0.0;

      const LCObjectVec& relCaloHitsToSimCaloHit = navCalorimeter->getRelatedToObjects(*j);
            
      for ( unsigned int k = 0; k < relCaloHitsToSimCaloHit.size(); ++k ) {
		  
	SimCalorimeterHit* simCaloHit = dynamic_cast<SimCalorimeterHit*>(relCaloHitsToSimCaloHit.at(k));
		
	for ( int l = 0; l < simCaloHit->getNMCContributions(); ++l ) {
	  
	  MCParticle* mcpOfSimCaloHit = simCaloHit->getParticleCont(l);
	  double eSimHitContribution = simCaloHit->getEnergyCont(l);

	  // search for SimCalorimeterHit related to not extrapolated tracks
	  bool isMCOfSimCalorimeterHitDaughtherOfNotExtrapolatedTrack = false;

	   
	  for ( std::vector<Track*>::const_iterator m = _tracksNotExtrapolatedIntoCalorimeter.begin(); m != _tracksNotExtrapolatedIntoCalorimeter.end(); ++m) {
	    
	    const LCObjectVec& relMCParticlesToTrack = navTracks->getRelatedToObjects(*m); 
	      
	    if ( relMCParticlesToTrack.size() > 1 ) std::cout << "Warning: More than one MCParticle related to track." << std::endl;
		
	    MCParticle* mcpOfTrack = 0;
	    
	    for(unsigned int n = 0; n < relMCParticlesToTrack.size(); ++n) {
	    
	      mcpOfTrack = dynamic_cast<MCParticle*>(relMCParticlesToTrack.at(n)); 
		
	      isMCOfSimCalorimeterHitDaughtherOfNotExtrapolatedTrack = isMCOfSimCalorimeterHitDaughtherOfNotExtrapolatedTrack || 
		MarlinUtil::isDaughterOf(mcpOfSimCaloHit,mcpOfTrack) || (mcpOfSimCaloHit == mcpOfTrack);
	      
	      if (isMCOfSimCalorimeterHitDaughtherOfNotExtrapolatedTrack) break;
	      
	    }
		
	    if (isMCOfSimCalorimeterHitDaughtherOfNotExtrapolatedTrack) {

	      std::vector<Track*>::const_iterator position = find(_tracksWhichWouldReachCalorimeter.begin(),_tracksWhichWouldReachCalorimeter.end(),(*m));
	       
	      if ( position == _tracksWhichWouldReachCalorimeter.end() ) _tracksWhichWouldReachCalorimeter.push_back(*m);

	      break;

	    }
	    
	  }

	  if ( !isMCOfSimCalorimeterHitDaughtherOfNotExtrapolatedTrack ) {


	    MCParticleVec parentsOfMCofSimCaloHit = MarlinUtil::getAllMCParents(mcpOfSimCaloHit);


	    // debug
	    if ( _debugLevel > 5 ) {
	    
	      if ( parentsOfMCofSimCaloHit.size() > 1 ) {
	      
		std::cout << "Warning: More than one (" << parentsOfMCofSimCaloHit.size() << ") MC parent found for MCParticle " << mcpOfSimCaloHit 
			  << " (" << mcpOfSimCaloHit->getPDG() << ") " << " in SimCalorimeterHit " << simCaloHit 
			  << " in method 'getChargedHitsAssignedToNeutralParticle()'" << std::endl;	      
	      
	      }
	    
	    }

	  
	    for ( MCParticleVec::const_iterator m = parentsOfMCofSimCaloHit.begin(); m != parentsOfMCofSimCaloHit.end(); ++m ) {
	    
	      // take only charged particles causing this SimCaloHit into account
	      std::string charge = _mcParticleHelper->getMCCharge( (*m)->getPDG() );
   

	      // debug
	      /*
	      std::cout << "CaloHit: " << (*j) << "  " << "E = " << (*j)->getEnergy() << "  " 
			<< "SimCaloHit: " << simCaloHit << "  " << "E = " << eSimHitContribution*calibration.at((*j)->getType()) << std::endl 
			<< "MCContrib: " << l+1 << "/" << simCaloHit->getNMCContributions() << "  " << mcpOfSimCaloHit << "  " << mcpOfSimCaloHit->getPDG() << "  " 
			<< mcpOfSimCaloHit->getEnergy() << "  " << "GeneratorStatus: " << mcpOfSimCaloHit->getGeneratorStatus() << "  " 
			<< "CreatedInSimulation: " << mcpOfSimCaloHit->isCreatedInSimulation() << "  " << "DecayedInTracker: " << mcpOfSimCaloHit->isDecayedInTracker() 
			<< std::endl
			<< "MCParents: " << ( m-parentsOfMCofSimCaloHit.begin()+1 ) << "/" << parentsOfMCofSimCaloHit.size() << "  " << (*m) << "  " << (*m)->getPDG() 
			<< "  " << (*m)->getEnergy() << "  " << charge << "  " << "GeneratorStatus: " << (*m)->getGeneratorStatus() << "  " 
			<< "CreatedInSimulation: " << (*m)->isCreatedInSimulation() << "  "
			<< "DecayedInTracker: " << (*m)->isDecayedInTracker() << std::endl << std::endl;
	      */
	    

	  

	      // charged MC contribution
	      if ( charge != "0" ) { 	   	   
		
		int type   = (*j)->getType();	     
		double eHitContribution = eSimHitContribution*calibration.at(type);
		
		relatedChargedEnergyPerCaloHit += eHitContribution;
	      
	      }
	      
	    }
	  
	  }

	}
	
      }

      if ( relatedChargedEnergyPerCaloHit > (hitEnergyFraction*energyCaloHit) ) {
	
	collectedCalotimeterHits.push_back(*j);
	++n;
	energy += energyCaloHit;

      }
     
    }

  }

  delete navCalorimeter;
  navCalorimeter = 0;
	    
  delete navTracks;
  navTracks = 0;

  return collectedCalotimeterHits;

}


double TrackBasedPFlow::getNotAssignedCalorimeterEnergy(LCEvent* evt, LCCollectionVec* reconstructedParticles) {

  double notAssignedCalorimeterEnergy = 0.0;
  double fullCalorimeterEnergy = 0.0;

  std::vector<CalorimeterHit*> allCalorimeterHits;


  try {
    
    std::vector< std::string >::const_iterator iter;
    const std::vector< std::string >* ColNames = evt->getCollectionNames();
    
    for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
      
      LCCollection* col = evt->getCollection( *iter ) ;
      
      if ( (col->getTypeName() == LCIO::CALORIMETERHIT) && ( (*iter == _colNameECAL) || (*iter == _colNameHCAL) ) ) {
	
	int NCalorimeterHits = col->getNumberOfElements();

	for(int j=0; j<NCalorimeterHits; ++j) {
	
	  CalorimeterHit* calorimeterHit = dynamic_cast<CalorimeterHit*>(col->getElementAt(j));

	  fullCalorimeterEnergy += calorimeterHit->getEnergy();
	  allCalorimeterHits.push_back(calorimeterHit);

	}

      }

    }

  }
  catch(DataNotAvailableException &e){std::cout << "no valid calorimeter hit collection in event " << _nEvt << std::endl; };



  // debug
  if (_drawOnCED) {
    MarlinCED::init(this);
    MarlinCED::newEvent(this,0);
  }


      
  for (LCCollectionVec::const_iterator i = reconstructedParticles->begin(); i != reconstructedParticles->end(); ++i) {

    ReconstructedParticle* recoParticle = dynamic_cast<ReconstructedParticle*>( (*i) );


    // debug
    if (_drawOnCED) { 

      unsigned int color = 0;
      
      // draw charged reconstructed particles
      if ( recoParticle->getCharge() != 0.0 ) {
	
	switch ( recoParticle->getType() ) {
	  
	case 1 : // e+/-
	  color = MarlinDrawUtil::getColor(11);
	  break;
	case 2 : // pi+/-
	  color = MarlinDrawUtil::getColor(211);
	  break;
	case 5 : // mu+/-
	  color = MarlinDrawUtil::getColor(13);
	  break;
	default : // assume it is a pi+/-
	  color = MarlinDrawUtil::getColor(211);
	  break;
	}
	
	MarlinCED::drawRecoParticle(recoParticle,0,2,color,8);
	
      }
      else {

	switch ( recoParticle->getType() ) {
		
	case 3 : // gamma
	  color = MarlinDrawUtil::getColor(22);
	  break;
	case 4 : // K0L, that means neutral hadron
	  color = MarlinDrawUtil::getColor(130);
	  break;
	default : // assume it is a K0L, that means a neutral hadron
	  color = MarlinDrawUtil::getColor(130);
	  break;
	}
	
	MarlinCED::drawRecoParticle(recoParticle,0,2,color,9);

      }
    }
    

    for (ClusterVec::const_iterator j = recoParticle->getClusters().begin(); j != recoParticle->getClusters().end(); ++j) {
	
      for (CalorimeterHitVec::const_iterator k = (*j)->getCalorimeterHits().begin(); k != (*j)->getCalorimeterHits().end(); ++k) {
	
	std::vector<CalorimeterHit*>::iterator position = find(allCalorimeterHits.begin(),allCalorimeterHits.end(),(*k));
       
	  if ( position !=  allCalorimeterHits.end() ) { 

	    allCalorimeterHits.erase(position);
	    	    
	  }
	  
      }

    }

  }


  // debug
  if (_drawOnCED) MarlinCED::drawMCParticleTree(evt,"MCParticle",0.05,4.0,15.5,50.0,1626.0,2500.0);
  

  for (std::vector<CalorimeterHit*>::const_iterator i = allCalorimeterHits.begin(); i != allCalorimeterHits.end(); ++i) {

    notAssignedCalorimeterEnergy += (*i)->getEnergy();

    if (_drawOnCED) {

      double hitEnergyInMIPs = 0.0;
    
      switch ( (*i)->getType() ) {
		
      case 0 : 
	hitEnergyInMIPs = ((*i)->getEnergy())/_mipCoeffEcal.at(0);
	break;
      case 1 : 
	hitEnergyInMIPs = ((*i)->getEnergy())/_mipCoeffEcal.at(1);
	break;
      case 2 : 
	hitEnergyInMIPs = ((*i)->getEnergy())/_mipCoeffHcal.at(0);
	break;
	
      }
    
      // debug
      // std::cout << (*i)->getEnergy() << "  " << hitEnergyInMIPs << "  " << (*i)->getType() << std::endl;
	      
      
      if (hitEnergyInMIPs < 0.5) { 
	
	ced_hit( (*i)->getPosition()[0],(*i)->getPosition()[1],(*i)->getPosition()[2], 0 | 4 << CED_LAYER_SHIFT, 2, 0xff696c );
      
      }
      else if ( (hitEnergyInMIPs >= 0.5) && (hitEnergyInMIPs < 1.7) ) {
	
	ced_hit( (*i)->getPosition()[0],(*i)->getPosition()[1],(*i)->getPosition()[2], 0 | 5 << CED_LAYER_SHIFT, 2, 0xff0000 );
	
      }
      else if ( (hitEnergyInMIPs >= 1.7) && (hitEnergyInMIPs < 3.5) ) {
	
	ced_hit( (*i)->getPosition()[0],(*i)->getPosition()[1],(*i)->getPosition()[2], 0 | 6 << CED_LAYER_SHIFT, 2, 0x0dff00 );
	
      }
      else {
	
	ced_hit( (*i)->getPosition()[0],(*i)->getPosition()[1],(*i)->getPosition()[2], 0 | 7 << CED_LAYER_SHIFT, 2, 0x2f6dff );
      
      }
      
      
    }
    
  }
  
  return notAssignedCalorimeterEnergy;
    
}


int TrackBasedPFlow::getTypeOfPositionOfCluster(ClusterImpl* cluster) {

  int type = 0;
  double smallestDistance = DBL_MAX;
  
  for(CalorimeterHitVec::const_iterator i = cluster->getCalorimeterHits().begin(); i != cluster->getCalorimeterHits().end(); ++i) {
    
    double distanceHitToPosition = sqrt( pow( ( cluster->getPosition()[0] - (*i)->getPosition()[0] ) ,2) + 
					 pow( ( cluster->getPosition()[1] - (*i)->getPosition()[1] ) ,2) + 
					 pow( ( cluster->getPosition()[2] - (*i)->getPosition()[2] ) ,2) );
    
    if ( distanceHitToPosition < smallestDistance ) {

      smallestDistance = distanceHitToPosition;
      type = (*i)->getType();

    }
    
  }

  return type;

}


int TrackBasedPFlow::getTypeOfPositionOfCluster(Cluster* cluster) {

  int type = 0;
  double smallestDistance = DBL_MAX;
  
  for(CalorimeterHitVec::const_iterator i = cluster->getCalorimeterHits().begin(); i != cluster->getCalorimeterHits().end(); ++i) {
    
    double distanceHitToPosition = sqrt( pow( ( cluster->getPosition()[0] - (*i)->getPosition()[0] ) ,2) + 
					 pow( ( cluster->getPosition()[1] - (*i)->getPosition()[1] ) ,2) + 
					 pow( ( cluster->getPosition()[2] - (*i)->getPosition()[2] ) ,2) );
    
    if ( distanceHitToPosition < smallestDistance ) {

      smallestDistance = distanceHitToPosition;
      type = (*i)->getType();

    }
    
  }

  return type;

}


bool TrackBasedPFlow::hasTrackSufficientNumberOfHits(Track* track, bool& minNTPCHitsReached, bool& minNNonTPCHitsReached) {
  
  // check if track fulfilling cuts on number of silicon and TPC hits
  int nOfTPCHits = 0;
  int nOfNonTPCHits = 0;
  bool trackHasSufficientNumberOfHits = false;
      
  for ( TrackerHitVec::const_iterator i = track->getTrackerHits().begin(); i != track->getTrackerHits().end(); ++i ) {
	
    if ( ((*i)->getType()) == 500 ) { // FIXME: hard coded number for type describing the TPC hits
      
      ++nOfTPCHits;
	  
      if ( nOfTPCHits >= _minNTPCHits) minNTPCHitsReached = true;	  	  

    }
	
    if ( ((*i)->getType()) != 500 ) { // FIXME: hard coded number for type describing the non TPC hits
	  
      ++nOfNonTPCHits;
	  
      if ( nOfNonTPCHits >= _minNNonTPCHits ) minNNonTPCHitsReached = true;
	  
    }
		
  }
      
  trackHasSufficientNumberOfHits = minNTPCHitsReached && minNNonTPCHitsReached;

  return trackHasSufficientNumberOfHits;

}




// FIXME: the following methods should be placed somewhere in a helix util class


double TrackBasedPFlow::getDistanceToHelix(double* point, Trajectory* helix) {

  LCVector3D pointVector(point[0],point[1],point[2]);  
  return getDistanceToHelix(pointVector,helix);

}

double TrackBasedPFlow::getDistanceToHelix(std::vector<double> point, Trajectory* helix) {

  LCVector3D pointVector(point.at(0),point.at(1),point.at(2));  
  return getDistanceToHelix(pointVector,helix);
  
}

double TrackBasedPFlow::getDistanceToHelix(LCVector3D point, Trajectory* helix) {

  // real code

  double sOfPoint = helix->getPathAt(point);
  LCVector3D pointProjected = helix->getPosition(sOfPoint);
  double distance = LCVector3D( pointProjected - point ).mag();

  return distance;

}

double* TrackBasedPFlow::getProjectedPointOnHelix(double* point, Trajectory* helix) {
  
  LCVector3D pointVector(point[0],point[1],point[2]);
  LCVector3D pointVectorProjected = getProjectedPointOnHelix(pointVector,helix);

  // user needs to care about delete
  double* pointReturn = new double[3];
  pointReturn[0] = pointVectorProjected.x();
  pointReturn[1] = pointVectorProjected.y();
  pointReturn[2] = pointVectorProjected.z();

  return pointReturn;

}

std::vector<double> TrackBasedPFlow::getProjectedPointOnHelix(std::vector<double> point, Trajectory* helix) {

  LCVector3D pointVector(point.at(0),point.at(1),point.at(2));
  LCVector3D pointVectorProjected = getProjectedPointOnHelix(pointVector,helix);

  std::vector<double> pointReturn;
  pointReturn.push_back(pointVectorProjected.x());
  pointReturn.push_back(pointVectorProjected.y());
  pointReturn.push_back(pointVectorProjected.z());

  return pointReturn;

}

LCVector3D TrackBasedPFlow::getProjectedPointOnHelix(LCVector3D point, Trajectory* helix) {

  double s = helix->getPathAt(point);  
  return helix->getPosition(s);

}

double TrackBasedPFlow::getPathLengthOnHelix(double* point, Trajectory* helix) {
  
  LCVector3D pointVector(point[0],point[1],point[2]);
  return getPathLengthOnHelix(pointVector,helix);
  
}

double TrackBasedPFlow::getPathLengthOnHelix(std::vector<double> point, Trajectory* helix) {
 
  LCVector3D pointVector(point.at(0),point.at(1),point.at(2));
  return getPathLengthOnHelix(pointVector,helix);

}

double TrackBasedPFlow::getPathLengthOnHelix(LCVector3D point, Trajectory* helix) {

  double s = helix->getPathAt(point);
  return s;

}

double TrackBasedPFlow::getPathLengthOnHelix(double* point1, double* point2, Trajectory* helix) {

  LCVector3D pointVector1(point1[0],point1[1],point1[2]);
  LCVector3D pointVector2(point2[0],point2[1],point2[2]);

  return getPathLengthOnHelix(pointVector1,pointVector2,helix);

}

double TrackBasedPFlow::getPathLengthOnHelix(std::vector<double> point1, std::vector<double> point2, Trajectory* helix) {

  LCVector3D pointVector1(point1.at(0),point1.at(1),point1.at(2));
  LCVector3D pointVector2(point2.at(0),point2.at(1),point2.at(2));

  return getPathLengthOnHelix(pointVector1,pointVector2,helix);

}

double TrackBasedPFlow::getPathLengthOnHelix(LCVector3D point1, LCVector3D point2, Trajectory* helix) {

  double s1 = helix->getPathAt(point1);
  double s2 = helix->getPathAt(point2);
  double ds = s2 - s1;

  return ds;
  
}
