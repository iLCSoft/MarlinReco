#include "FullLDCTracking.h"
#include <iostream>
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <iostream>
#include <math.h>
#include <map>
#include <marlin/Global.h>
#include "ClusterShapes.h"
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/VXDParameters.h>
#include <gear/GearParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/BField.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>

using namespace lcio ;
using namespace marlin ;

FullLDCTracking aFullLDCTracking ;

FullLDCTracking::FullLDCTracking() : Processor("FullLDCTracking") {  
  _description = "Performs full tracking in LDC detector" ;  


  // Input tracker hit collections

  registerInputCollection(LCIO::TRACKERHIT,
			  "FTDHitCollection",
			  "FTD Hit Collection Name",
			  _FTDTrackerHitCollection,
			  std::string("FTDTrackerHits"));  

  registerInputCollection(LCIO::TRACKERHIT,
			  "VTXHitCollection",
			  "VTX Hit Collection Name",
			  _VTXTrackerHitCollection,
			  std::string("VTXTrackerHits"));  

  registerInputCollection(LCIO::TRACKERHIT,
			  "SITHitCollection",
			  "SIT Hit Collection Name",
			  _SITTrackerHitCollection,
			  std::string("SITTrackerHits"));
  
  registerInputCollection(LCIO::TRACKERHIT,
			  "SETHitCollection",
			  "SET Hit Collection Name",
			  _SETTrackerHitCollection,
			  std::string("SETTrackerHits"));
  
  registerInputCollection(LCIO::TRACKERHIT,
			  "ETDHitCollection",
			  "ETD Hit Collection Name",
			  _ETDTrackerHitCollection,
			  std::string("ETDTrackerHits"));
  
  registerInputCollection(LCIO::TRACKERHIT,
			  "TPCHitCollection",
			  "TPC Hit Collection Name",
			  _TPCTrackerHitCollection,
			  std::string("TPCTrackerHits"));
  

  // Input track collections

  registerInputCollection(LCIO::TRACK,
			  "TPCTracks",
			  "TPC Track Collection",
			  _TPCTrackCollection,
			  std::string("TPCTracks"));

  registerInputCollection(LCIO::TRACK,
			  "SiTracks",
			  "Si Track Collection",
			  _SiTrackCollection,
			  std::string("SiTracks"));

  // Input relation collections

  registerInputCollection(LCIO::LCRELATION,
			  "TPCTracksMCPRelColl",
			  "TPC Track to MCP Relation Collection Name",
			  _TPCTrackMCPCollName,
			  std::string("TPCTracksMCP"));

  registerInputCollection(LCIO::LCRELATION,
			  "SiTracksMCPRelColl",
			  "Si Track to Collection",
			  _SiTrackMCPCollName,
			  std::string("SiTracksMCP"));

  // Output track collection
  registerOutputCollection(LCIO::TRACK,
			   "LDCTrackCollection",
			   "LDC track collection name",
			   _LDCTrackCollection,
			   std::string("LDCTracks"));


  // Output relation collection
  registerOutputCollection(LCIO::LCRELATION,
			   "LDCTrackMCPRelCollection",
			   "Collection name for the LDC track to MCParticle relations",
			   _LDCTrackMCPCollection,
			   std::string("LDCTracksMCP"));

  // Refitted track collection
  registerOutputCollection(LCIO::TRACK,
			   "ReffitedTPCTrackCollection",
			   "Refitted TPC track collection name",
			   _RefittedTPCTrackCollection,
			   std::string("RefittedTPCTracks"));

  registerOutputCollection(LCIO::LCRELATION,
			   "RefittedTPCTrackMCPRelCollection",
			   "Collection name for the refitted TPC track to MCParticle relations",
			   _RefittedTPCTrackMCPCollection,
			   std::string("RefittedTPCTracksMCP"));

  registerOutputCollection(LCIO::TRACK,
			   "RefittedSiTrackCollection",
			   "Refitted Si track collection name",
			   _RefittedSiTrackCollection,
			   std::string("RefittedSiTracks"));

  registerOutputCollection(LCIO::LCRELATION,
			   "RefittedSiTrackMCPRelCollection",
			   "Collection name for the refittedSi track to MCParticle relations",
			   _RefittedSiTrackMCPCollection,
			   std::string("RefittedSiTracksMCP"));

  
  // steering parameters

  registerProcessorParameter("D0CutForMerging",
			     "Cut on D0 difference for merging of Si and TPC segments",
			     _d0CutForMerging,
			     float(500.0));

  registerProcessorParameter("Z0CutForMerging",
			     "Cut on Z0 difference for merging of Si and TPC segments",
			     _z0CutForMerging,
			     float(1000.0));

  registerProcessorParameter("OmegaCutForMerging",
			     "Cut on Omega difference for merging Si and TPC segments",
			     _dOmegaForMerging,
			     float(0.25));

  registerProcessorParameter("AngleCutForMerging",
			     "Cut on Opening Angle for merging Si and TPC segments",
			     _angleForMerging,
			     float(0.10));

  registerProcessorParameter("Chi2FitCut",
			     "Cut on fit Chi2",
			     _chi2FitCut,
			     float(100.0));

  registerProcessorParameter("Chi2PrefitCut",
			     "Cut on fit Chi2",
			     _chi2PrefitCut,
			     float(1.0e+5));

  registerProcessorParameter("CreateMap",
			     "Create Track to MCP Relations",
			     _createMap,
			     int(1));

  registerProcessorParameter("Debug",
			     "Activate debugging?",
			     _debug,
			     int(0));

  registerProcessorParameter("UseExtraPoint",
			     "Use Extra Point in Fit",
			     _useExtraPoint,
			     int(0));

  registerProcessorParameter("OptFit",
			     "Option for the LDC Track fit",
			     _optFit,
			     int(4));

  registerProcessorParameter("ForceSiTPCMerging",
			     "Force merging of Si and TPC segments?",
			     _forceMerging,
			     int(0));

  registerProcessorParameter("D0CutForForcedMerging",
			     "Cut on D0 difference for forced merging of Si and TPC segments",
			     _d0CutForForcedMerging,
			     float(50.));

  registerProcessorParameter("Z0CutForForcedMerging",
			     "Cut on Z0 difference for forced merging of Si and TPC segments",
			     _z0CutForForcedMerging,
			     float(200.));

  registerProcessorParameter("OmegaCutForForcedMerging",
			     "Cut on Omega difference for forced merging of Si and TPC segments",
			     _dOmegaForForcedMerging,
			     float(0.15));

  registerProcessorParameter("AngleCutForForcedMerging",
			     "Cut on Opening Angle for forced merging of Si and TPC segments",
			     _angleForForcedMerging,
			     float(0.05));

  registerProcessorParameter("RefitTPCTracks",
			     "Refit TPC Tracks ?",
			     _refitTPCTracks,
			     int(1));

  registerProcessorParameter("RefitSiTracks",
			     "Refit Si Tracks ?",
			     _refitSiTracks,
			     int(0));

  registerProcessorParameter("StoreRefittedTPCTracks",
			     "Store Refitted TPC Tracks ?",
			     _storeRefittedTPCTracks,
			     int(0));

  registerProcessorParameter("StoreRefittedSiTracks",
			     "Store Refitted Si Tracks ?",
			     _storeRefittedSiTracks,
			     int(0));

  registerProcessorParameter("ForceTPCSegmentsMerging",
			     "Force merging of TPC Segments?",
			     _mergeTPCSegments,
			     int(1));

  registerProcessorParameter("D0CutToMergeTPCSegments",
			     "Cut on D0 difference for merging TPC segments",
			     _d0CutToMergeTPC,
			     float(100.));

  registerProcessorParameter("Z0CutToMergeTPCSegments",
			     "Cut on Z0 difference for merging TPC segments",
			     _z0CutToMergeTPC,
			     float(5000.0));

  registerProcessorParameter("DeltaPCutToMergeTPCSegments",
			     "Cut on dP/P difference for merging TPC segments",
			     _dPCutToMergeTPC,
			     float(0.1));

  registerProcessorParameter("PtCutToMergeTPCSegments",
			     "Cut on Pt of tracks for merging TPC segments",
			     _PtCutToMergeTPC,
			     float(1.2));



  registerProcessorParameter("cosThetaCutHighPtMerge",
			     "Cut on cos theta between the two momentum vectors when considering merger of high Pt tracks",
			     _cosThetaCutHighPtMerge,
			     float(0.99));

  registerProcessorParameter("cosThetaCutSoftHighPtMerge",
			     "cut on cos theta between the two momentum vectors when considering merger of high Pt tracks for softer dp/p cut",
			     _cosThetaCutSoftHighPtMerge,
			     float(0.998));

  registerProcessorParameter("momDiffCutHighPtMerge",
			     "cut on dp/p when considering merger of high Pt tracks",
			     _momDiffCutHighPtMerge,
			     float(0.01));

  registerProcessorParameter("momDiffCutSoftHighPtMerge",
			     "softer cut on dp/p when considering merger of high Pt tracks",
			     _momDiffCutSoftHighPtMerge,
			     float(0.25));

  registerProcessorParameter("hitDistanceCutHighPtMerge",
			     "cut on 3D distance between hit and helix extrapolation when considering merger of high Pt tracks",
			     _hitDistanceCutHighPtMerge,
			     float(25.0));

  registerProcessorParameter("maxHitDistanceCutHighPtMerge",
			     "cut for max 3D distance between any hit and helix extrapolation when considering merger of high Pt tracks",
			     _maxHitDistanceCutHighPtMerge,
			     float(50.0));

  registerProcessorParameter("maxFractionOfOutliersCutHighPtMerge",
			     "cut on maximum fraction of outliers when considering merger of high Pt tracks",
			     _maxFractionOfOutliersCutHighPtMerge,
			     float(0.95));





  registerProcessorParameter("CutOnTPCHits",
			     "Cut on the number of the TPC hits for tracks with no Si hits",
			     _cutOnTPCHits,
			     int(35));
  registerProcessorParameter("CutOnSiHits",
			     "Cut on the number of the Si hits for tracks with no TPC hits",
			     _cutOnSiHits,
			     int(4));

  registerProcessorParameter("aParameterForIPError",
			     "Parameter a to define minimal IP error",
			     _aParIpReso,
			     float(0.002));
 
  registerProcessorParameter("bParameterForIPError",
			     "Parameter b to define minimal IP error",
			     _bParIpReso,
			     float(0.0076));

  registerProcessorParameter("sParameterForIPError",
			     "Parameter s to define minimal IP error",
			     _sParIpReso,
			     float(0.75));

  registerProcessorParameter("AssignVTXHits",
			     "Assign left over VTX hits",
			     _assignVTXHits,
			     int(1));

  registerProcessorParameter("AssignFTDHits",
			     "Assign left over FTD hits",
			     _assignFTDHits,
			     int(1));

  registerProcessorParameter("AssignSITHits",
			     "Assign left over SIT hits",
			     _assignSITHits,
			     int(1));

  registerProcessorParameter("AssignTPCHits",
			     "Assign left over TPC hits",
			     _assignTPCHits,
			     int(1));

  registerProcessorParameter("AssignSETHits",
			     "Assign SET Hits",
			     _assignSETHits,
			     int(1));


  registerProcessorParameter("AssignETDHits",
			     "Assign ETD Hits",
			     _assignETDHits,
			     int(1));

  registerProcessorParameter("NHitsExtrapolation",
			     "number of hits for outward extrapolation",
			     _nHitsExtrapolation,
			     int(35));

  registerProcessorParameter("StoreHitsInFit",
			     "Store only hits used in fit?",
			     _storeHitsInFit,
			     int(0));


  registerProcessorParameter("VTXHitToTrackDistance",
			     "Cut on distance between track and VTX hits",
			     _distCutForVTXHits,
			     float(1.5));


  registerProcessorParameter("FTDHitToTrackDistance",
			     "Cut on distance between track and FTD hits",
			     _distCutForFTDHits,
			     float(2.0));


  registerProcessorParameter("SITHitToTrackDistance",
			     "Cut on distance between track and SIT hits",
			     _distCutForSITHits,
			     float(2.0));

  registerProcessorParameter("SETHitToTrackDistance",
			     "Cut on distance between track and SET hits",
			     _distCutForSETHits,
			     float(2.0));


  registerProcessorParameter("ETDHitToTrackDistance",
			     "Cut on distance between track and ETD hits",
			     _distCutForETDHits,
			     float(10.0));


  registerProcessorParameter("TPCHitToTrackDistance",
			     "Cut on distance between track and TPC hits",
			     _distCutForTPCHits,
			     float(15.0));


  registerProcessorParameter("OptFitTPC",
			     "Option for TPC tracks refitting",
			     _optFitTPC,
			     int(2));

  registerProcessorParameter("OptFitSi",
			     "Option for Si tracks refitting",
			     _optFitSi,
			     int(2));


  registerProcessorParameter("CutOnTrackD0",
			     "Cut on the track parameter D0",
			     _d0TrkCut,
			     float(500.));


  registerProcessorParameter("CutOnTrackZ0",
			     "Cut on the track parameter Z0",
			     _z0TrkCut,
			     float(500.));


  registerProcessorParameter("ForbidOverlapInZTPC",
			     "Forbid overlap in Z for the merged TPC segments",
			     _forbidOverlapInZTPC,
			     int(0));

  registerProcessorParameter("ForbidOverlapInZComb",
			     "Forbid overlap in Z for combining TPC segments with tracks having Si hits",
			     _forbidOverlapInZComb,
			     int(0));




}



void FullLDCTracking::init() { 

  printParameters();  
  _nRun = -1 ;
  _nEvt = 0 ;
  PI = acos(-1.);
  PIOVER2 = 0.5*PI;
  TWOPI = 2*PI;

}

void FullLDCTracking::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
  _nEvt = 0;
  std::cout << std::endl;
  std::cout << "FullLDCTracking ---> new run : run number = " << _nRun << std::endl;

} 

void FullLDCTracking::processEvent( LCEvent * evt ) { 

  _evt = evt;

  //  _bField = Global::parameters->getFloatVal("BField");
  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;

  float parIpReso[3];
  parIpReso[0] = _aParIpReso;
  parIpReso[1] = _bParIpReso;
  parIpReso[2] = _sParIpReso;
  _trackFit.setParametersForIPErrors(parIpReso);

  if (_debug>=1) {
      std::cout << std::endl;
      std::cout << "FullLDCTracking -> run = " << _nRun 
		<< "  event = " << _nEvt << std::endl;
      std::cout << std::endl;
  }

  prepareVectors( evt );
  //  std::cout << "prepareVectors done..." << std::endl;
  if (_storeRefittedTPCTracks>0 && _refitTPCTracks>0) {
    AddTrackColToEvt(evt,_allTPCTracks,
		     _RefittedTPCTrackCollection,_RefittedTPCTrackMCPCollection);
    // std::cout << "Collection of refitted TPC tracks is added to event..." << std::endl;
  }
  if (_storeRefittedSiTracks>0 && _refitSiTracks>0) {
     AddTrackColToEvt(evt,_allSiTracks,
		      _RefittedSiTrackCollection,_RefittedSiTrackMCPCollection);
     // std::cout << "Collection of refitted Si tracks is added to event..." << std::endl; 
  }
  if (_debug>=2)std::cout << "************************************Merge TPC/Si" << std::endl;
  MergeTPCandSiTracks();
  if (_debug>=2)std::cout << "************************************Merging done..." << std::endl;
  MergeTPCandSiTracksII();
  if (_debug>=2)std::cout << "************************************Merging II done..." << std::endl;
  Sorting(_allCombinedTracks);
  if (_debug>=2)std::cout << "************************************Sorting done..." << std::endl;
  SelectCombinedTracks();
  if (_debug>=2)std::cout << "************************************Selection of combined tracks done..." << std::endl;
  AddNotCombinedTracks( );
  if (_debug>=2)std::cout << "************************************Not combined tracks added..." << std::endl;
  //CheckTracks( );

  AddNotAssignedHits();
  if (_debug>=2)std::cout << "***********************************Not assigned hits added..." << std::endl;
  AddTrackColToEvt(evt,_trkImplVec,
		   _LDCTrackCollection,_LDCTrackMCPCollection);
  //  std::cout << "Collections added to event..." << std::endl;
  CleanUp();
  //  std::cout << "Cleanup is done..." << std::endl;
  _nEvt++;
  //  getchar();
  if (_debug>=2) {
    std::cout << std::endl;
    std::cout << std::endl;
  }
}

void FullLDCTracking::AddTrackColToEvt(LCEvent * evt, TrackExtendedVec & trkVec, 
				       std::string TrkColName, std::string RelColName) {
  
  LCCollectionVec * colTRK = new LCCollectionVec(LCIO::TRACK);
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  colTRK->setFlag( trkFlag.getFlag()  ) ;  

  if (_debug >= 2)
    std::cout << "Collection " << TrkColName
	      << " is being added to event " << std::endl;

  LCCollectionVec * colRel = NULL;
  
  if (_createMap) {
    colRel = new LCCollectionVec(LCIO::LCRELATION);
    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    colRel->setFlag( lcFlag.getFlag()  ) ;
  }
  
  int nTrkCand = int(trkVec.size());
  
  int nTotTracks = 0;
  float eTot = 0.0;
  float pxTot = 0.0;
  float pyTot = 0.0;
  float pzTot = 0.0;

  for (int iTRK=0;iTRK<nTrkCand;++iTRK) {
    TrackExtended * trkCand = trkVec[iTRK];
    TrackImpl * track = new TrackImpl();
    track->setOmega(trkCand->getOmega());
    track->setTanLambda(trkCand->getTanLambda());
    track->setPhi(trkCand->getPhi());
    track->setZ0(trkCand->getZ0());
    track->setD0(trkCand->getD0());
    track->setChi2(trkCand->getChi2());
    track->setNdf(trkCand->getNDF());
    track->setCovMatrix(trkCand->getCovMatrix());
    TrackerHitExtendedVec hitVec = trkCand->getTrackerHitExtendedVec();
    int nHits = int(hitVec.size());
    std::vector<MCParticle*> mcPointers ;
    std::vector<int> mcHits ;
    mcPointers.clear();
    mcHits.clear();
    int nHitsVTX = 0;
    int nHitsFTD = 0;
    int nHitsSIT = 0;
    int nHitsTPC = 0;
    int nHitsSET = 0;
    int nHitsETD = 0;
    int nHitsVTXInFit = 0;
    int nHitsFTDInFit = 0;
    int nHitsSITInFit = 0;
    int nHitsTPCInFit = 0;
    int nHitsSETInFit = 0;
    int nHitsETDInFit = 0;
    float r2Min = 1.0e+30;
    for (int iH=0;iH<nHits;++iH) {
      TrackerHitExtended * hitExt = hitVec[iH];
      TrackerHit * hit = hitExt->getTrackerHit();
      bool isUsedInFit = hitExt->getUsedInFit();
      if (_storeHitsInFit==0)
	  track->addHit(hit);
      int det = hitExt->getDet();
      if (det == 3)
	nHitsVTX++;
      if (det == 2)
	nHitsFTD++;
      if (det == 7)
	nHitsSIT++;
      if (det == 1)
	nHitsTPC++;
      if (det == 8)
	nHitsSET++;
      if (det == 9)
	nHitsETD++;
      float hitX = float(hit->getPosition()[0]);
      float hitY = float(hit->getPosition()[1]);
      float hitR2 = hitX*hitX+hitY*hitY;
      if (hitR2<r2Min)
	  r2Min = hitR2;
      if (isUsedInFit) {
	  if (_storeHitsInFit!=0)
	    track->addHit(hit);
	  if (det == 3)
	    nHitsVTXInFit++;
	  if (det == 2)
	    nHitsFTDInFit++;
	  if (det == 7)
	    nHitsSITInFit++;
	  if (det == 1)
	    nHitsTPCInFit++;
	  if (det == 8)
	    nHitsSETInFit++;
	  if (det == 9)
	    nHitsETDInFit++;
      }


      if (_createMap > 0) {
	int nSH = int(hit->getRawHits().size());
	for (int ish=0;ish<nSH;++ish) {
	  SimTrackerHit * simHit = dynamic_cast<SimTrackerHit*>(hit->getRawHits()[ish]);
	  MCParticle * mcp = simHit->getMCParticle();
	  bool found = false;
	  int nMCP = int(mcPointers.size());
	  for (int iMCP=0;iMCP<nMCP;++iMCP) {
	    if (mcp == mcPointers[iMCP]) {
	      found = true;
	      mcHits[iMCP]++;
	      break;
	    }
	  }
	  if (!found) {
	    mcPointers.push_back(mcp);
	    mcHits.push_back(1);
	  }
	}
      }
    }
    GroupTracks * group = trkCand->getGroupTracks();
    if (group != NULL) {
      TrackExtendedVec trkVecGrp = group->getTrackExtendedVec();
      int nGrTRK = int(trkVecGrp.size());
      for (int iGr=0;iGr<nGrTRK;++iGr) {
	TrackExtended * subTrack = trkVecGrp[iGr];
	track->addTrack(subTrack->getTrack());
      }
    }

    float d0TrkCand = trkCand->getD0();
    float z0TrkCand = trkCand->getZ0();
    float phi0TrkCand = trkCand->getPhi();

    float RefPoint[3];

    RefPoint[0] = -d0TrkCand*sin(phi0TrkCand);
    RefPoint[1] =  d0TrkCand*cos(phi0TrkCand);
    RefPoint[2] =  z0TrkCand;

    track->setReferencePoint(RefPoint);
    track->setIsReferencePointPCA( true );
    track->setRadiusOfInnermostHit(sqrt(r2Min));

    track->subdetectorHitNumbers().resize(12);
    track->subdetectorHitNumbers()[0] = nHitsVTXInFit;
    track->subdetectorHitNumbers()[1] = nHitsFTDInFit;
    track->subdetectorHitNumbers()[2] = nHitsSITInFit;
    track->subdetectorHitNumbers()[3] = nHitsTPCInFit;
    track->subdetectorHitNumbers()[4] = nHitsSETInFit;
    track->subdetectorHitNumbers()[5] = nHitsETDInFit;    
    track->subdetectorHitNumbers()[6] = nHitsVTX;
    track->subdetectorHitNumbers()[7] = nHitsFTD;
    track->subdetectorHitNumbers()[8] = nHitsSIT;
    track->subdetectorHitNumbers()[9] = nHitsTPC;
    track->subdetectorHitNumbers()[10]= nHitsSET;
    track->subdetectorHitNumbers()[11]= nHitsETD;
    
    
    int nHitsSiInFit = nHitsVTXInFit+nHitsFTDInFit+nHitsSITInFit;
    bool rejectTrack = (nHitsTPCInFit<_cutOnTPCHits) && (nHitsSiInFit<=0);
    rejectTrack = rejectTrack || ( (nHitsTPCInFit<=0) && (nHitsSiInFit<_cutOnSiHits) );
    rejectTrack = rejectTrack || ( fabs(d0TrkCand) > _d0TrkCut ) || ( fabs(z0TrkCand) > _z0TrkCut );

    if ( rejectTrack ) {
      delete track;
    }
    else {
      float omega = trkCand->getOmega();
      float tanLambda = trkCand->getTanLambda();
      float phi0 = trkCand->getPhi();
      float d0 = trkCand->getD0();
      float z0 = trkCand->getZ0();
      HelixClass helix;
      helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
      float trkPx = helix.getMomentum()[0];
      float trkPy = helix.getMomentum()[1];
      float trkPz = helix.getMomentum()[2];
      float trkP = sqrt(trkPx*trkPx+trkPy*trkPy+trkPz*trkPz);
      eTot += trkP;
      pxTot += trkPx;
      pyTot += trkPy;
      pzTot += trkPz;	
      nTotTracks++;
      colTRK->addElement(track);
      if (_createMap > 0) {
	int nRel = int(mcPointers.size());
	for (int k=0;k<nRel;++k) {
	  LCRelationImpl* tpclcRel = new LCRelationImpl;
	  MCParticle * mcp = mcPointers[k];
	  tpclcRel->setFrom (track);
	  tpclcRel->setTo (mcp);
	  float weight = (float)(mcHits[k])/(float)(track->getTrackerHits().size());
	  tpclcRel->setWeight(weight);
	  colRel->addElement(tpclcRel);
	}
      }
    }
  }
  if (_debug >= 1) {
    std::cout << std::endl;
    std::cout << "Number of accepted " << TrkColName << " = " 
	      << nTotTracks << std::endl;
    std::cout << "Total 4-momentum of " << TrkColName << " : E = " << eTot
	      << " Px = " << pxTot
	      << " Py = " << pyTot
	      << " Pz = " << pzTot << std::endl;
    std::cout << std::endl;
  }  

  evt->addCollection(colTRK,TrkColName.c_str());
  if (_createMap)
    evt->addCollection(colRel,RelColName.c_str());

}


void FullLDCTracking::prepareVectors(LCEvent * event ) {


  _allTPCHits.clear();
  _allVTXHits.clear();
  _allFTDHits.clear();
  _allSITHits.clear();
  _allSETHits.clear();
  _allETDHits.clear();
  _allTPCTracks.clear();
  _allSiTracks.clear();
  _allCombinedTracks.clear();
  _allNonCombinedTPCTracks.clear();
  _allNonCombinedSiTracks.clear();
  _trkImplVec.clear();
  _candidateCombinedTracks.clear();

  //  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  //  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;


  std::map <TrackerHit*,TrackerHitExtended*> mapTrackerHits;

  // Reading TPC hits
  try {
    LCCollection * col = event->getCollection(_TPCTrackerHitCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
//       double tpcRPhiResConst = gearTPC.getDoubleVal("tpcRPhiResConst");
//       double tpcRPhiResDiff  = gearTPC.getDoubleVal("tpcRPhiResDiff");
//       double aReso = tpcRPhiResConst*tpcRPhiResConst;
//       double driftLenght = gearTPC.getMaxDriftLength() - fabs(hit->getPosition()[2]);
//       if (driftLenght <0) { 
//         driftLenght = 0.10;
//       }
//       double bReso = tpcRPhiResDiff*tpcRPhiResDiff;
//       double tpcRPhiRes = sqrt(aReso + bReso*driftLenght);
//       double tpcZRes = gearTPC.getDoubleVal("tpcZRes");
//       hitExt->setResolutionRPhi(float(tpcRPhiRes));
//       hitExt->setResolutionZ(float(tpcZRes));      

      // Covariance Matrix in LCIO is defined in XYZ convert to R-Phi-Z
      // For no error in r

      //      double phi = atan2(hit->getPosition()[1],hit->getPosition()[0]); 
      double tpcRPhiRes = sqrt(hit->getCovMatrix()[0] + hit->getCovMatrix()[2]);
      double tpcZRes = sqrt(hit->getCovMatrix()[5]);

 
      // f77 tracking code works in cm
//       tpcRPhiRes = 0.1 * tpcRPhiRes;
//       tpcZRes = 0.1 * tpcZRes;

      hitExt->setResolutionRPhi(float(tpcRPhiRes));
      hitExt->setResolutionZ(float(tpcZRes));

      hitExt->setType(int(3));
      hitExt->setDet(int(1));
      _allTPCHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    if (_debug >= 2)
      std::cout << _TPCTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  };

  // Reading FTD Hits
  try {
    LCCollection * col = event->getCollection(_FTDTrackerHitCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      hitExt->setResolutionRPhi(float(sqrt(hit->getCovMatrix()[0])));
      hitExt->setResolutionZ(0.1);
      //      hitExt->setResolutionRPhi(_resolutionRPhi_FTD);
      //      hitExt->setResolutionZ(_resolutionZ_FTD);
      hitExt->setType(int(2));
      hitExt->setDet(int(2));
      _allFTDHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    if (_debug >= 2)
      std::cout << _FTDTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading ETD Hits
  try {
    LCCollection * col = event->getCollection(_ETDTrackerHitCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      hitExt->setResolutionRPhi(float(sqrt(hit->getCovMatrix()[0])));
      hitExt->setResolutionZ(0.1);
      //      hitExt->setResolutionRPhi(_resolutionRPhi_FTD);
      //      hitExt->setResolutionZ(_resolutionZ_FTD);
      hitExt->setType(int(2));
      hitExt->setDet(int(9));
      _allETDHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    if (_debug >= 2)
      std::cout << _ETDTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading SIT Hits
  try {
    LCCollection * col = event->getCollection(_SITTrackerHitCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
       hitExt->setResolutionRPhi(float(sqrt(hit->getCovMatrix()[2])));
       hitExt->setResolutionZ(float(sqrt(hit->getCovMatrix()[5])));
//      hitExt->setResolutionRPhi(_resolutionRPhi_SIT);
//      hitExt->setResolutionZ(_resolutionZ_SIT);
      hitExt->setType(int(3));
      hitExt->setDet(int(7));
      _allSITHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    if (_debug >=2 )
      std::cout << _SITTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading SET Hits
//    printf("SET hits ----->\n");
//    printf(" id    r_hit   phi_hit    z_hit  e(r-phi)  e(z)  \n");
  //     "  0   1807.64   -1.97     26.45   0.010   0.010
  //     "  1   1812.64   -1.97     26.51   0.010   0.010
  try {
    LCCollection * col = event->getCollection(_SETTrackerHitCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      hitExt->setResolutionRPhi(float(sqrt(hit->getCovMatrix()[2])));
      hitExt->setResolutionZ(float(sqrt(hit->getCovMatrix()[5])));
//       float x_hit = float(hit->getPosition()[0]);
//       float y_hit = float(hit->getPosition()[1]);
//       float z_hit = float(hit->getPosition()[2]);
//       float phi_hit = atan2( y_hit, x_hit );
//       float r_hit = sqrt( x_hit*x_hit + y_hit*y_hit );
//       printf("%3i  %8.2f  %6.3f  %8.2f  %6.3f  %6.3f   ", 
//  	     ielem, r_hit, phi_hit, z_hit, 
//  	     hitExt->getResolutionRPhi(), 
//  	     hitExt->getResolutionZ() );
//       std::cout << hit << std::endl;
      hitExt->setType(int(3));
      hitExt->setDet(int(8));
      _allSETHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    if (_debug >=2 )
      std::cout << _SETTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading VTX Hits
  try {
    LCCollection * col = event->getCollection(_VTXTrackerHitCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      hitExt->setResolutionRPhi(float(sqrt(hit->getCovMatrix()[2])));
      hitExt->setResolutionZ(float(sqrt(hit->getCovMatrix()[5])));
      //      hitExt->setResolutionRPhi(_resolutionRPhi_VTX);
      //      hitExt->setResolutionZ(_resolutionZ_VTX);
      hitExt->setType(int(3));
      hitExt->setDet(int(3));
      _allVTXHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    if (_debug >=2 )
      std::cout << _VTXTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }


  // Reading TPC Tracks
  try {
    LCCollection * col = event->getCollection(_TPCTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    if (_debug >= 2) {
      std::cout << std::endl;
      std::cout << "Number of TPC Tracks = " << nelem << std::endl;
      std::cout << " Trk       p          D0         Z0       Px       Py       Pz    ntpc ndf Chi2/ndf" << std::endl;
      //           "  0  1.111   0.059      0.022    -0.54     0.61    -0.45    0.185
    }
    for (int iTrk=0; iTrk<nelem; ++iTrk) {
      Track * tpcTrack = dynamic_cast<Track*>(col->getElementAt(iTrk));
      TrackExtended * trackExt = new TrackExtended( tpcTrack );
           TrackerHitVec hitVec = tpcTrack->getTrackerHits();
      int nHits = int(hitVec.size());
      trackExt->setOmega(tpcTrack->getOmega());
      trackExt->setTanLambda(tpcTrack->getTanLambda());
      trackExt->setPhi(tpcTrack->getPhi());
      trackExt->setD0(tpcTrack->getD0());
      trackExt->setZ0(tpcTrack->getZ0());
      float cov[15];
      float param[5];
      float reso[4];
      param[0] = tpcTrack->getOmega();
      param[1] = tpcTrack->getTanLambda();
      param[2] = tpcTrack->getPhi();
      param[3] = tpcTrack->getD0();
      param[4] = tpcTrack->getZ0();
      // Setting crude resolutions for TPC Track parameters 
      reso[0] = 1e-3;  // dPt/Pt = 1e-3*Pt
      reso[1] = 0.001; // angular resolution = 1mrad
      reso[2] = 1.0;   // D0 reso = 1mm 
      reso[3] = 1.0;   // Z0 reso = 1mm	
      _trackFit.CrudeErrorEstimates(_bField,reso,param,cov);
      trackExt->setCovMatrix(cov);
      trackExt->setNDF(tpcTrack->getNdf());
      trackExt->setChi2(tpcTrack->getChi2());            
      char strg[200];
      HelixClass helixTPC;
      for (int iHit=0;iHit<nHits;++iHit) {
	TrackerHit * hit = hitVec[iHit];
	TrackerHitExtended * hitExt = mapTrackerHits[hit];
	hitExt->setTrackExtended( trackExt );
	trackExt->addTrackerHitExtended( hitExt );	
      }      
      if (_refitTPCTracks > 0) { // refit TPC tracks
	float * x_h = new float[nHits];
	float * y_h = new float[nHits];
	float * z_h = new float[nHits];
	int * idet_h = new int[nHits];
	int * ityp_h = new int[nHits];
	int * lhits = new int[nHits];
	float * rR_h = new float[nHits];
	float * rZ_h = new float[nHits];
	for (int iHit=0;iHit<nHits;++iHit) {
	  TrackerHit * hit = hitVec[iHit];
	  TrackerHitExtended * hitExt = mapTrackerHits[hit];
	  x_h[iHit] = float(hit->getPosition()[0]);
	  y_h[iHit] = float(hit->getPosition()[1]);
	  z_h[iHit] = float(hit->getPosition()[2]);
	  idet_h[iHit] = hitExt->getDet();
	  ityp_h[iHit] = hitExt->getType();
	  rR_h[iHit] = hitExt->getResolutionRPhi();
	  rZ_h[iHit] = hitExt->getResolutionZ();
	}
	int NPT = nHits;
	float chi2_D;
	int ndf_D;
	float chi2rphi,chi2z;
	float par[5] = {0.,0.,0.,0.,0.};
	for (int iPar=0;iPar<5;++iPar)
	  par[iPar] = param[iPar];
	float epar[15];
	float refPoint[3];	
	helixTPC.Initialize_Canonical(par[2],par[3],par[4],par[0],par[1],_bField);
	float momX = helixTPC.getMomentum()[0];
	float momY = helixTPC.getMomentum()[1];
	float momZ = helixTPC.getMomentum()[2];
	int ierr = _trackFit.DoFitting(_useExtraPoint,_optFitTPC,NPT,_bField,idet_h,ityp_h,_chi2PrefitCut,
				     x_h,y_h,z_h,rR_h,rZ_h,
				     par,epar,refPoint,chi2_D,ndf_D,chi2rphi,chi2z,lhits);
	helixTPC.Initialize_Canonical(par[2],par[3],par[4],par[0],par[1],_bField);
	float momXNew = helixTPC.getMomentum()[0];
	float momYNew = helixTPC.getMomentum()[1];
	float momZNew = helixTPC.getMomentum()[2];
	float dPx = momX - momXNew;
	float dPy = momY - momYNew;
	float dPz = momZ - momZNew;
	float dMom = sqrt(dPx*dPx+dPy*dPy+dPz*dPz);
	if (ierr>=0 && dMom < 0.05) { // redefine track parameters if fit is stable	  
	  trackExt->setCovMatrix(epar);
	  trackExt->setOmega(par[0]);	  
	  trackExt->setTanLambda(par[1]);
	  trackExt->setPhi(par[2]);
	  trackExt->setD0(par[3]);
	  trackExt->setZ0(par[4]);
	  trackExt->setChi2(chi2_D);
	  trackExt->setNDF(ndf_D);
	}
	delete[] idet_h;
	delete[] ityp_h;
	delete[] lhits;
	delete[] x_h;
	delete[] y_h;
	delete[] z_h;
	delete[] rZ_h;
	delete[] rR_h;
      }  
      if (_debug>=2) {
	float d0TPC = trackExt->getD0();
	float z0TPC = trackExt->getZ0();
	float omegaTPC = trackExt->getOmega();
	float phi0TPC = trackExt->getPhi();
	float tanLTPC = trackExt->getTanLambda();
	float Chi2TPC = trackExt->getChi2()/float(trackExt->getNDF());
	const int ndfTPC = trackExt->getNDF();
	helixTPC.Initialize_Canonical(phi0TPC,d0TPC,z0TPC,omegaTPC,tanLTPC,_bField);
	float pxTPC = helixTPC.getMomentum()[0];
	float pyTPC = helixTPC.getMomentum()[1];
	float pzTPC = helixTPC.getMomentum()[2];
	const float ptot = sqrt(pxTPC*pxTPC+pyTPC*pyTPC+pzTPC*pzTPC);
	sprintf(strg,"%3i  %9.3f  %9.3f  %9.3f  %7.2f  %7.2f  %7.2f %4i %4i %8.3f",iTrk,
		ptot, d0TPC,z0TPC,pxTPC,pyTPC,pzTPC,nHits,ndfTPC,Chi2TPC);
	std::cout << strg << std::endl;
      }
      _allTPCTracks.push_back( trackExt );                
    }      
  }
  catch ( DataNotAvailableException &e) {
    if (_debug >= 2 )
      std::cout << _TPCTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading Si Tracks
  try {
    LCCollection * col = event->getCollection(_SiTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    if (_debug >= 2) {
      std::cout << std::endl;
      std::cout << "Number of Si Tracks = " << nelem << std::endl;
      std::cout << " Trk       p          D0         Z0       Px       Py       Pz   hitsSi ndf Chi2/ndf" << std::endl;
    }
    for (int iTrk=0; iTrk<nelem; ++iTrk) {
      Track * siTrack = dynamic_cast<Track*>(col->getElementAt(iTrk));
      TrackExtended * trackExt = new TrackExtended( siTrack );
      TrackerHitVec hitVec = siTrack->getTrackerHits();
      int nHits = int(hitVec.size());
      trackExt->setOmega(siTrack->getOmega());
      trackExt->setTanLambda(siTrack->getTanLambda());
      trackExt->setPhi(siTrack->getPhi());
      trackExt->setD0(siTrack->getD0());
      trackExt->setZ0(siTrack->getZ0());
      float cov[15];
      float param[5];
      float reso[4];
      param[0] = siTrack->getOmega();
      param[1] = siTrack->getTanLambda();
      param[2] = siTrack->getPhi();
      param[3] = siTrack->getD0();
      param[4] = siTrack->getZ0();      
      if (_refitSiTracks < 0) { // Setting crude resolutions for the Si Track parameters
	reso[0] = 1e-2;  // dPt/Pt = 1e-2*Pt
	reso[1] = 0.001; // angular resolution = 1mrad
	reso[2] = 0.1;   // D0 reso = 1mm 
	reso[3] = 0.1;   // Z0 reso = 1mm
	_trackFit.CrudeErrorEstimates(_bField,reso,param,cov);
      }
      if (_refitSiTracks == 0) { // using cov matrix obtained from Si track fits
	const FloatVec Cov = siTrack->getCovMatrix();
	int NC = int(Cov.size());
	for (int ic=0;ic<NC;ic++) {
	  cov[ic] =  Cov[ic];
	}	
      }
      trackExt->setCovMatrix(cov);
      trackExt->setNDF(siTrack->getNdf());
      trackExt->setChi2(siTrack->getChi2());      
      char strg[200];
      HelixClass helixSi;
      for (int iHit=0;iHit<nHits;++iHit) {
	TrackerHit * hit = hitVec[iHit];
	TrackerHitExtended * hitExt = mapTrackerHits[hit];
	hitExt->setTrackExtended( trackExt );
	trackExt->addTrackerHitExtended( hitExt );	
      }
      if (_refitSiTracks > 0) { // refit Si tracks
	float * x_h = new float[nHits];
	float * y_h = new float[nHits];
	float * z_h = new float[nHits];
	int * idet_h = new int[nHits];
	int * ityp_h = new int[nHits];
	int * lhits = new int[nHits];
	float * rR_h = new float[nHits];
	float * rZ_h = new float[nHits];
	for (int iHit=0;iHit<nHits;++iHit) {
	  TrackerHit * hit = hitVec[iHit];
	  TrackerHitExtended * hitExt = mapTrackerHits[hit];
	  x_h[iHit] = float(hit->getPosition()[0]);
	  y_h[iHit] = float(hit->getPosition()[1]);
	  z_h[iHit] = float(hit->getPosition()[2]);
	  idet_h[iHit] = hitExt->getDet();
	  ityp_h[iHit] = hitExt->getType();
	  rR_h[iHit] = hitExt->getResolutionRPhi();
	  rZ_h[iHit] = hitExt->getResolutionZ();
	}
	int NPT = nHits;
	float chi2_D;
	int ndf_D;
	float chi2rphi,chi2z;
	float par[5];
	float epar[15];
	float refPoint[3];
	for (int iPar=0;iPar<5;++iPar)
	  par[iPar] = param[iPar];
	int ierr = _trackFit.DoFitting(_useExtraPoint,_optFitSi,NPT,_bField,idet_h,ityp_h,_chi2PrefitCut,
				       x_h,y_h,z_h,rR_h,rZ_h,
				       par,epar,refPoint,chi2_D,ndf_D,chi2rphi,chi2z,lhits);
	if (ierr>=0) {	  
	  trackExt->setCovMatrix(epar);
	  trackExt->setOmega(par[0]);	  
	  trackExt->setTanLambda(par[1]);
	  trackExt->setPhi(par[2]);
	  trackExt->setD0(par[3]);
	  trackExt->setZ0(par[4]);
	  trackExt->setChi2(chi2_D);
	  trackExt->setNDF(ndf_D);
	}
	delete[] idet_h;
	delete[] ityp_h;
	delete[] x_h;
	delete[] y_h;
	delete[] z_h;
	delete[] rZ_h;
	delete[] rR_h;
	delete[] lhits;
      }
      if (_debug>=2) {
	float d0Si = trackExt->getD0();
	float z0Si = trackExt->getZ0();
	float omegaSi = trackExt->getOmega();
	float phi0Si = trackExt->getPhi();
	float tanLSi = trackExt->getTanLambda();
	helixSi.Initialize_Canonical(phi0Si,d0Si,z0Si,omegaSi,tanLSi,_bField);
	float pxSi = helixSi.getMomentum()[0];
	float pySi = helixSi.getMomentum()[1];
	float pzSi = helixSi.getMomentum()[2];
	const float pTot = sqrt(pxSi*pxSi+pySi*pySi+pzSi*pzSi);
	const int ndfSi = trackExt->getNDF();
	float Chi2Si = trackExt->getChi2()/float(trackExt->getNDF());
	sprintf(strg,"%3i  %9.3f  %9.3f  %9.3f  %7.2f  %7.2f  %7.2f %4i %4i %8.3f",iTrk,
		pTot, d0Si,z0Si,pxSi,pySi,pzSi,nHits, ndfSi, Chi2Si);
	std::cout << strg << std::endl;
      }
      if(nHits>0){
	_allSiTracks.push_back( trackExt );
      }else{
	delete trackExt;
      }
    }
    if (_debug>=2)
      std::cout << std::endl;
  }
  catch ( DataNotAvailableException &e) {
    if (_debug >= 2)
      std::cout << _SiTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }



}

void FullLDCTracking::CleanUp(){

  int nNonCombTpc = int(_allNonCombinedTPCTracks.size());  
  for (int i=0;i<nNonCombTpc;++i) {
      TrackExtended * trkExt = _allNonCombinedTPCTracks[i];
      GroupTracks * group = trkExt->getGroupTracks();
      delete group;
  }
  _allNonCombinedTPCTracks.clear();

  int nNonCombSi = int(_allNonCombinedSiTracks.size());  
  for (int i=0;i<nNonCombSi;++i) {
      TrackExtended * trkExt = _allNonCombinedSiTracks[i];
      GroupTracks * group = trkExt->getGroupTracks();
      delete group;
  }
  _allNonCombinedSiTracks.clear();
  
  int nSITHits = int(_allSITHits.size());
  for (int i=0;i<nSITHits;++i) {
    TrackerHitExtended * hitExt = _allSITHits[i];
    delete hitExt;
  }
  _allSITHits.clear();

  int nSETHits = int(_allSETHits.size());
  for (int i=0;i<nSETHits;++i) {
    TrackerHitExtended * hitExt = _allSETHits[i];
    delete hitExt;
  }
  _allSETHits.clear();

  int nTPCHits = int(_allTPCHits.size());
  for (int i=0;i<nTPCHits;++i) {
    TrackerHitExtended * hitExt = _allTPCHits[i];
    delete hitExt;
  }
  _allTPCHits.clear();

  int nFTDHits = int(_allFTDHits.size());
  for (int i=0;i<nFTDHits;++i) {
    TrackerHitExtended * hitExt = _allFTDHits[i];
    delete hitExt;
  }
  _allFTDHits.clear();

  int nETDHits = int(_allETDHits.size());
  for (int i=0;i<nETDHits;++i) {
    TrackerHitExtended * hitExt = _allETDHits[i];
    delete hitExt;
  }
  _allETDHits.clear();

  int nVTXHits = int(_allVTXHits.size());
  for (int i=0;i<nVTXHits;++i) {
    TrackerHitExtended * hitExt = _allVTXHits[i];
    delete hitExt;
  }
  _allVTXHits.clear();

  int nSiTrk = int(_allSiTracks.size());
  for (int i=0;i<nSiTrk;++i) {
    TrackExtended * trkExt = _allSiTracks[i];
    delete trkExt;
  }
  _allSiTracks.clear();

  int nTPCTrk = int(_allTPCTracks.size());
  for (int i=0;i<nTPCTrk;++i) {
    TrackExtended * trkExt = _allTPCTracks[i];
    delete trkExt;
  }
  _allTPCTracks.clear();

  int nCombTrk = int(_allCombinedTracks.size());
  for (int i=0;i<nCombTrk;++i) {
    TrackExtended * trkExt = _allCombinedTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    delete group;
    delete trkExt;    
  }
  _allCombinedTracks.clear();

//   int nImplTrk = int(_trkImplVec.size());
//   for (int i=0;i<nImplTrk;++i) {
//     TrackExtended * trkImpl = _trkImplVec[i];
//     delete trkImpl;
//   }
  _trkImplVec.clear();

  //AS: Dont delete the individual entries, some of them are cleared elsewhere, I think
  _candidateCombinedTracks.clear();
}

void FullLDCTracking::MergeTPCandSiTracks() {

  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks  = int(_allSiTracks.size());

  for (int iTPC=0;iTPC<nTPCTracks;++iTPC) {
    TrackExtended * tpcTrackExt = _allTPCTracks[iTPC];
    for (int iSi=0;iSi<nSiTracks;++iSi) {
      TrackExtended * siTrackExt = _allSiTracks[iSi];
      int iComp = 0;
      float angle = 0;
      float dOmega = CompareTrkII(siTrackExt,tpcTrackExt,_d0CutForMerging,_z0CutForMerging,iComp,angle);
      if ( (dOmega<_dOmegaForMerging) && (angle<_angleForMerging) && !VetoMerge(siTrackExt,tpcTrackExt)) {
	TrackExtended *combinedTrack = CombineTracks(tpcTrackExt,siTrackExt);       
	if (combinedTrack != NULL) {
	  _allCombinedTracks.push_back( combinedTrack );
	  _candidateCombinedTracks.insert(tpcTrackExt);
	  _candidateCombinedTracks.insert(siTrackExt);
	  if (_debug == 3 ) {
	    int iopt = 1;
	    PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
	  }
	}else{
	  if (_debug == 3 ) {
	    int iopt = 6;
	    PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
	  }
	}
      }
      else {
	if (_debug == 3 ) {
	  int iopt = 6;
	  PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
	}
      }
    }
  }


}



void FullLDCTracking::MergeTPCandSiTracksII() {

  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks  = int(_allSiTracks.size());

  for (int iTPC=0;iTPC<nTPCTracks;++iTPC) {
    TrackExtended * tpcTrackExt = _allTPCTracks[iTPC];
    if(_candidateCombinedTracks.find(tpcTrackExt) != _candidateCombinedTracks.end() )continue;

    for (int iSi=0;iSi<nSiTracks;++iSi) {

      TrackExtended * siTrackExt = _allSiTracks[iSi];
      if(_candidateCombinedTracks.find(siTrackExt)!= _candidateCombinedTracks.end() )continue;
      int iComp = 0;
      float angleSignificance = 0;
      float significance = CompareTrkIII(siTrackExt,tpcTrackExt,_d0CutForMerging,_z0CutForMerging,iComp,angleSignificance);
      if ( (significance<10) && (angleSignificance<5) && !VetoMerge(tpcTrackExt,siTrackExt) ) {
	TrackExtended * combinedTrack = CombineTracks(tpcTrackExt,siTrackExt);

	if (combinedTrack != NULL) {

	  _allCombinedTracks.push_back( combinedTrack );
	  if (_debug == 3 ) {
	    int iopt = 1;
	    PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
	  }
	}else{
	  if (_debug == 3 ) {
	    int iopt = 6;
	    PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
	  }
	}
      }
      else {
	if (_debug == 3 ) {
	  int iopt = 6;
	  PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
	}
      }
    }
  }
}




TrackExtended * FullLDCTracking::CombineTracks(TrackExtended * tpcTrack, TrackExtended * siTrack) {

  TrackExtended * OutputTrack = NULL;

  TrackerHitExtendedVec siHitVec = siTrack->getTrackerHitExtendedVec();
  TrackerHitExtendedVec tpcHitVec = tpcTrack->getTrackerHitExtendedVec();

  int nSiHits = int(siHitVec.size());
  int nTPCHits = int(tpcHitVec.size());

  int nHits = nTPCHits + nSiHits;

  float * xh = new float[nHits];
  float * yh = new float[nHits];
  float * zh = new float[nHits];
  float * rphi_reso = new float[nHits];
  float * z_reso = new float[nHits];
  int * ityp_h = new int[nHits];
  int * idet_h = new int[nHits];
  int * lhits = new int[nHits];
  float par[5];
  float epar[15];
  for (int ih=0;ih<nSiHits;++ih) {
    TrackerHit * trkHit = siHitVec[ih]->getTrackerHit();
    float rR = siHitVec[ih]->getResolutionRPhi();
    float rZ = siHitVec[ih]->getResolutionZ();
    xh[ih] = float(trkHit->getPosition()[0]);
    yh[ih] = float(trkHit->getPosition()[1]);
    zh[ih] = float(trkHit->getPosition()[2]);
    rphi_reso[ih] = rR;
    z_reso[ih] = rZ;
    ityp_h[ih] =  siHitVec[ih]->getType();
    idet_h[ih] =  siHitVec[ih]->getDet();
  }      
  for (int ih=0;ih<nTPCHits;++ih) {
    TrackerHit * trkHit = tpcHitVec[ih]->getTrackerHit();
    float rR = tpcHitVec[ih]->getResolutionRPhi();
    float rZ = tpcHitVec[ih]->getResolutionZ();
    xh[ih+nSiHits] = float(trkHit->getPosition()[0]);
    yh[ih+nSiHits] = float(trkHit->getPosition()[1]);
    zh[ih+nSiHits] = float(trkHit->getPosition()[2]);
    rphi_reso[ih+nSiHits] = rR;
    z_reso[ih+nSiHits] = rZ;
    ityp_h[ih+nSiHits] =  tpcHitVec[ih]->getType();
    idet_h[ih+nSiHits] =  tpcHitVec[ih]->getDet();
  }      

  int NPT = nHits;
  float chi2_D;
  int ndf_D;
  float chi2rphi,chi2z;
  par[0] = tpcTrack->getOmega();
  par[1] = siTrack->getTanLambda();
  par[2] = siTrack->getPhi();
  par[3] = siTrack->getD0();
  par[4] = siTrack->getZ0();
  float refPoint[3];
  int ierr = _trackFit.DoFitting(_useExtraPoint,_optFit,NPT,_bField,
				 idet_h,ityp_h,_chi2PrefitCut,
				 xh,yh,zh,rphi_reso,z_reso,
				 par,epar,refPoint,chi2_D,ndf_D,
				 chi2rphi,chi2z,lhits);

  
  float chiQ = chi2_D/float(ndf_D);
  if (ierr == 0 && chiQ > 0.001) {
    float omega = par[0];
    float tanlambda = par[1];
    float phi0 = par[2];
    float d0 = par[3];
    float z0 = par[4];    
    float chi2Fit = chi2_D/float(ndf_D);
    if (chi2Fit > _chi2FitCut) {
      delete[] xh;
      delete[] yh;
      delete[] zh;
      delete[] idet_h;
      delete[] ityp_h;
      delete[] rphi_reso;
      delete[] z_reso;      
      delete[] lhits;
      return OutputTrack;
    }    
    OutputTrack = new TrackExtended();
    GroupTracks * group = new GroupTracks();
    group->addTrackExtended(siTrack);
    group->addTrackExtended(tpcTrack);
    OutputTrack->setGroupTracks(group);
    OutputTrack->setOmega(omega);
    OutputTrack->setTanLambda(tanlambda);
    OutputTrack->setPhi(phi0);
    OutputTrack->setZ0(z0);
    OutputTrack->setD0(d0);
    OutputTrack->setChi2(chi2_D);
    OutputTrack->setNDF(ndf_D);
    OutputTrack->setCovMatrix(epar);
    if (_optFit == 4) {
      // Preserving positiveness of cov. matrix
      float eD0 = siTrack->getCovMatrix()[0];
      float eZ0 = siTrack->getCovMatrix()[9];
      float sD0 = sqrt(eD0/epar[0]);
      float sZ0 = sqrt(eZ0/epar[9]);
      epar[0]  = eD0;
      epar[9]  = eZ0;
      epar[1]  = sD0*epar[1];
      epar[3]  = sD0*epar[3];
      epar[10] = sD0*epar[10];
      epar[7]  = sZ0*epar[7];
      epar[8]  = sZ0*epar[8];
      epar[13] = sZ0*epar[13];
      epar[6]  = sD0*sZ0*epar[6];
      OutputTrack->setZ0(siTrack->getZ0());
      OutputTrack->setD0(siTrack->getD0());
      OutputTrack->setCovMatrix(epar);
    }
    int iLH = 0;
    for (int i=0;i<nSiHits;++i) {
      TrackerHitExtended * hitExt = siHitVec[i];
      OutputTrack->addTrackerHitExtended(hitExt);
      if (lhits[iLH]>0)
	  hitExt->setUsedInFit(true);
      else 
	  hitExt->setUsedInFit(false);
      iLH++;
    }
    for (int i=0;i<nTPCHits;++i) {
      TrackerHitExtended * hitExt = tpcHitVec[i];
      OutputTrack->addTrackerHitExtended(hitExt);
      if (lhits[iLH]>0)
	  hitExt->setUsedInFit(true);
      else 
	  hitExt->setUsedInFit(false);
      iLH++;
    }
  }
  delete[] xh;
  delete[] yh;
  delete[] zh;
  delete[] idet_h;
  delete[] ityp_h;
  delete[] rphi_reso;
  delete[] z_reso;
  delete[] lhits;
  
  return OutputTrack;
}


TrackExtended * FullLDCTracking::TrialCombineTracks(TrackExtended * tpcTrack, TrackExtended * siTrack) {

  TrackExtended * OutputTrack = NULL;

  TrackerHitExtendedVec siHitVec = siTrack->getTrackerHitExtendedVec();
  TrackerHitExtendedVec tpcHitVec = tpcTrack->getTrackerHitExtendedVec();

  int nSiHits = int(siHitVec.size());
  int nTPCHits = int(tpcHitVec.size());

  int nHits = nTPCHits + nSiHits;

  float * xh = new float[nHits];
  float * yh = new float[nHits];
  float * zh = new float[nHits];
  float * rphi_reso = new float[nHits];
  float * z_reso = new float[nHits];
  int * ityp_h = new int[nHits];
  int * idet_h = new int[nHits];
  int * lhits = new int[nHits];
  float par[5];
  float epar[15];
  for (int ih=0;ih<nSiHits;++ih) {
    TrackerHit * trkHit = siHitVec[ih]->getTrackerHit();
    float rR = siHitVec[ih]->getResolutionRPhi();
    float rZ = siHitVec[ih]->getResolutionZ();
    xh[ih] = float(trkHit->getPosition()[0]);
    yh[ih] = float(trkHit->getPosition()[1]);
    zh[ih] = float(trkHit->getPosition()[2]);
    rphi_reso[ih] = rR;
    z_reso[ih] = rZ;
    ityp_h[ih] =  siHitVec[ih]->getType();
    idet_h[ih] =  siHitVec[ih]->getDet();
  }      
  for (int ih=0;ih<nTPCHits;++ih) {
    TrackerHit * trkHit = tpcHitVec[ih]->getTrackerHit();
    float rR = tpcHitVec[ih]->getResolutionRPhi();
    float rZ = tpcHitVec[ih]->getResolutionZ();
    xh[ih+nSiHits] = float(trkHit->getPosition()[0]);
    yh[ih+nSiHits] = float(trkHit->getPosition()[1]);
    zh[ih+nSiHits] = float(trkHit->getPosition()[2]);
    rphi_reso[ih+nSiHits] = rR;
    z_reso[ih+nSiHits] = rZ;
    ityp_h[ih+nSiHits] =  tpcHitVec[ih]->getType();
    idet_h[ih+nSiHits] =  tpcHitVec[ih]->getDet();
  }      

  int NPT = nHits;
  float chi2_D;
  int ndf_D;
  float chi2rphi,chi2z;
  par[0] = tpcTrack->getOmega();
  par[1] = siTrack->getTanLambda();
  par[2] = siTrack->getPhi();
  par[3] = siTrack->getD0();
  par[4] = siTrack->getZ0();
  float refPoint[3];
  int ierr = _trackFit.DoFitting(_useExtraPoint,_optFit,NPT,_bField,
				 idet_h,ityp_h,_chi2PrefitCut,
				 xh,yh,zh,rphi_reso,z_reso,
				 par,epar,refPoint,chi2_D,ndf_D,
				 chi2rphi,chi2z,lhits);

  
  float chiQ = chi2_D/float(ndf_D);
  if (ierr == 0 && chiQ > 0.001) {
    float omega = par[0];
    float tanlambda = par[1];
    float phi0 = par[2];
    float d0 = par[3];
    float z0 = par[4];    
    float chi2Fit = chi2_D/float(ndf_D);
    if (chi2Fit > _chi2FitCut) {
      delete[] xh;
      delete[] yh;
      delete[] zh;
      delete[] idet_h;
      delete[] ityp_h;
      delete[] rphi_reso;
      delete[] z_reso;      
      delete[] lhits;
      return OutputTrack;
    }    
    OutputTrack = new TrackExtended();
    GroupTracks * group = new GroupTracks();
    group->addTrackExtended(siTrack);
    group->addTrackExtended(tpcTrack);
    OutputTrack->setGroupTracks(group);
    OutputTrack->setOmega(omega);
    OutputTrack->setTanLambda(tanlambda);
    OutputTrack->setPhi(phi0);
    OutputTrack->setZ0(z0);
    OutputTrack->setD0(d0);
    OutputTrack->setChi2(chi2_D);
    OutputTrack->setNDF(ndf_D);
    OutputTrack->setCovMatrix(epar);
    if (_optFit == 4) {
      // Preserving positiveness of cov. matrix
      float eD0 = siTrack->getCovMatrix()[0];
      float eZ0 = siTrack->getCovMatrix()[9];
      float sD0 = sqrt(eD0/epar[0]);
      float sZ0 = sqrt(eZ0/epar[9]);
      epar[0]  = eD0;
      epar[9]  = eZ0;
      epar[1]  = sD0*epar[1];
      epar[3]  = sD0*epar[3];
      epar[10] = sD0*epar[10];
      epar[7]  = sZ0*epar[7];
      epar[8]  = sZ0*epar[8];
      epar[13] = sZ0*epar[13];
      epar[6]  = sD0*sZ0*epar[6];
      OutputTrack->setZ0(siTrack->getZ0());
      OutputTrack->setD0(siTrack->getD0());
      OutputTrack->setCovMatrix(epar);
    }
  }

  delete[] xh;
  delete[] yh;
  delete[] zh;
  delete[] idet_h;
  delete[] ityp_h;
  delete[] rphi_reso;
  delete[] z_reso;
  delete[] lhits;
  
  return OutputTrack;
}


void FullLDCTracking::SortingTrackHitPairs(TrackHitPairVec & trackHitPairVec) {

  int sizeOfVector = int(trackHitPairVec.size());
  TrackHitPair *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++)
      {
	one = trackHitPairVec[j];
	two = trackHitPairVec[j+1];
	float oneQ = one->getDistance();
	float twoQ = two->getDistance();
	if( oneQ > twoQ )
	  {
	    Temp = trackHitPairVec[j];
	    trackHitPairVec[j] = trackHitPairVec[j+1];
	    trackHitPairVec[j+1] = Temp;
	  }
      }  
    

}

void FullLDCTracking::Sorting(TrackExtendedVec & trackVec) {

  int sizeOfVector = int(trackVec.size());
  TrackExtended *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++)
      {
	one = trackVec[j];
	two = trackVec[j+1];
	float oneQ = one->getChi2()/float(one->getNDF());
	float twoQ = two->getChi2()/float(two->getNDF());
	if( oneQ > twoQ )
	  {
	    Temp = trackVec[j];
	    trackVec[j] = trackVec[j+1];
	    trackVec[j+1] = Temp;
	  }
      }  
}

void FullLDCTracking::SelectCombinedTracks() {

  int nCombTrk = int(_allCombinedTracks.size());

  for (int i=0; i<nCombTrk;++i) {
    TrackExtended * trkExt = _allCombinedTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    TrackExtendedVec tracks = group->getTrackExtendedVec();
    int nTracks = int(tracks.size());
    if (nTracks == 2) {
      TrackExtended * firstTrack = tracks[0];
      TrackExtended * secondTrack = tracks[1];
      if ((firstTrack->getGroupTracks() == NULL) &&
	  (secondTrack->getGroupTracks() == NULL) ) {
	firstTrack->setGroupTracks(group);
	secondTrack->setGroupTracks(group);	
	TrackerHitExtendedVec firstVec = firstTrack->getTrackerHitExtendedVec();
	TrackerHitExtendedVec secondVec = secondTrack->getTrackerHitExtendedVec();
	int nFirst = int(firstVec.size());
	int nSecond = int(secondVec.size());
	float edges[2];
	edges[0] = 1.0e+20;
	edges[1] = -1.0e+20;
	for (int iF=0;iF<nFirst;++iF) {
	  TrackerHitExtended * trkHitExt = firstVec[iF];
	  TrackerHit * trkHit = trkHitExt->getTrackerHit();
	  float zpos = float(trkHit->getPosition()[2]);
	  if (zpos>edges[1])
	    edges[1] = zpos;
	  if (zpos<edges[0])
	    edges[0] = zpos;	  
	}
	for (int iS=0;iS<nSecond;++iS) {
	  TrackerHitExtended * trkHitExt = secondVec[iS];
	  TrackerHit * trkHit = trkHitExt->getTrackerHit();
	  float zpos = float(trkHit->getPosition()[2]);
	  if (zpos>edges[1])
	    edges[1] = zpos;
	  if (zpos<edges[0])
	    edges[0] = zpos;
	}
	group->setEdges(edges);
	_trkImplVec.push_back(trkExt);	
	
	if (_debug == 3) {
	  int iopt = 1;
	  PrintOutMerging(secondTrack,firstTrack,iopt);
	}	
      }
    }else{
      if(nTracks>2 && _debug>1)std::cout << " MORE THAN TWO TRACKS " << std::endl;
    }
  }


}

void FullLDCTracking::AddNotCombinedTracks() {  

  int nTPCTrk = int(_allTPCTracks.size());
  int nSiTrk = int(_allSiTracks.size());

  // we need some buffer vector
  TrackExtendedVec allMergedTracks;
  allMergedTracks.clear();

  // forcing merging of Si and TPC track segments
  if (_forceMerging==1) { 
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExtTPC = _allTPCTracks[i];
      GroupTracks * groupTPC = trkExtTPC->getGroupTracks();
      if (groupTPC == NULL) {
	float diffMin = 1.0e+20;  
	TrackExtended * siTrkToAttach = NULL;
	for (int j=0;j<nSiTrk;++j) {
	  TrackExtended * trkExtSi = _allSiTracks[j];
	  GroupTracks * groupSi = trkExtSi->getGroupTracks();
	  if (groupSi == NULL) {
	    int iComp = 0;
	    //	    float deltaP = CompareTrk(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp);
	    float angle(0.);
	    float angleSignificance(0.);

	    float dOmega = CompareTrkII(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp,angle);
	    float significance = CompareTrkIII(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp,angleSignificance);
	    //	    if (deltaP < _dPCutForForcedMerging) {
	    if ( ((dOmega<_dOmegaForForcedMerging) && (angle<_angleForForcedMerging)) ||
		 ((significance<5)                 && (angleSignificance<5))
		 ) {
	      float chi2O = dOmega/_dOmegaForForcedMerging;
	      float chi2A = angle/_angleForForcedMerging;
	      float deltaP = chi2O*chi2O + chi2A*chi2A; 
	      if (deltaP<diffMin) {
		diffMin = deltaP;
		siTrkToAttach = trkExtSi;
	      }
	    }else{
	      if (_debug==3) {
		int  iopt = 7;
		std::cout << significance << " " << angleSignificance << std::endl;
		PrintOutMerging(trkExtTPC,trkExtSi,iopt);
	      }
	    }
	  }
	}

	if (siTrkToAttach!=NULL) {
	    TrackExtended * trkExtSi = siTrkToAttach; 
	    TrackExtended * OutputTrack = new TrackExtended();
	    GroupTracks * group = new GroupTracks();
	    group->addTrackExtended(trkExtSi);
	    group->addTrackExtended(trkExtTPC);
	    OutputTrack->setGroupTracks(group);
//	    trkExtSi->setGroupTracks(group);
//	    trkExtTPC->setGroupTracks(group);	      
	    OutputTrack->setOmega(trkExtTPC->getOmega());
	    OutputTrack->setTanLambda(trkExtSi->getTanLambda());
	    OutputTrack->setPhi(trkExtSi->getPhi());
	    OutputTrack->setZ0(trkExtSi->getZ0());
	    OutputTrack->setD0(trkExtSi->getD0());
	    float covMatTPC[15];
	    float covMatSi[15];
	    float covMat[15];
	    for (int iCov=0;iCov<15;++iCov) {
		covMatTPC[iCov] = trkExtTPC->getCovMatrix()[iCov];
		covMatSi[iCov] = trkExtSi->getCovMatrix()[iCov];		    
		covMat[iCov] = covMatSi[iCov];
	    }
	    float scaling = sqrt(covMatTPC[5]/covMatSi[5]);
	    covMat[5] = covMatTPC[5];
	    covMat[3] = scaling*covMatSi[3];
	    covMat[4] = scaling*covMatSi[4];
	    covMat[8] = scaling*covMatSi[8];
	    covMat[12] = scaling*covMatSi[12];	      
	    OutputTrack->setCovMatrix(covMat);
	    TrackerHitExtendedVec tpcHitVec = trkExtTPC->getTrackerHitExtendedVec();
	    TrackerHitExtendedVec siHitVec = trkExtSi->getTrackerHitExtendedVec();	      
	    int nTPCHits = int( tpcHitVec.size());
	    int nSiHits = int( siHitVec.size());	      
	    float edges[2];
	    edges[0] = 1.0e+20;
	    edges[1] = -1.0e+20;
	    for (int iH=0;iH<nSiHits;++iH) {
		TrackerHitExtended * hitExt = siHitVec[iH];
		OutputTrack->addTrackerHitExtended(hitExt);
		hitExt->setUsedInFit(true);
		TrackerHit * hit = hitExt->getTrackerHit();
		float zpos = float(hit->getPosition()[2]);
		if (zpos<edges[0])
		  edges[0] = zpos;
		if (zpos>edges[1])
		  edges[1] = zpos;
	    }	  
	    for (int iH=0;iH<nTPCHits;++iH) {
		TrackerHitExtended * hitExt = tpcHitVec[iH];
		OutputTrack->addTrackerHitExtended(hitExt);
		hitExt->setUsedInFit(true); 
		TrackerHit * hit = hitExt->getTrackerHit();
		float zpos = float(hit->getPosition()[2]);
		if (zpos<edges[0])
		  edges[0] = zpos;
		if (zpos>edges[1])
		  edges[1] = zpos;
	    }
	    group->setEdges(edges);
	    OutputTrack->setChi2(diffMin); // will be replaced if necessary
	    OutputTrack->setNDF(int(1));   // will be replaced if necessary
	    _allCombinedTracks.push_back( OutputTrack );
	    allMergedTracks.push_back( OutputTrack );
	}	    
      }
    }


    int nMerged = int(allMergedTracks.size());
    if (nMerged>0) {
	Sorting(allMergedTracks);	
	for (int iM=0;iM<nMerged;++iM) {
	    TrackExtended * mergedTrack = allMergedTracks[iM];
	    GroupTracks * grpTrk = mergedTrack->getGroupTracks();
	    TrackExtendedVec trkVec = grpTrk->getTrackExtendedVec();
	    TrackExtended * trkTPC = NULL;
	    TrackExtended * trkSi = NULL;
	    int nT = int(trkVec.size());
	    if (nT==2) {
		trkTPC = trkVec[0];
		trkSi = trkVec[1];
		GroupTracks * groupTPC = trkTPC->getGroupTracks();
		GroupTracks * groupSi  = trkSi->getGroupTracks();
		if (groupTPC == NULL && groupSi == NULL) {
		    trkTPC->setGroupTracks( grpTrk );
		    trkSi->setGroupTracks( grpTrk );
		    TrackerHitExtendedVec hitVec = mergedTrack->getTrackerHitExtendedVec();
		    int nhits = int(hitVec.size());
		    int totNdf = 2*nhits - 5;
		    float totChi2 = trkTPC->getChi2() + trkSi->getChi2();
		    mergedTrack->setNDF( totNdf );
		    mergedTrack->setChi2( totChi2 );
		    if (_debug == 3) {
		      int iopt = 2;
		      PrintOutMerging(trkTPC,trkSi,iopt);
		    }
		    _trkImplVec.push_back( mergedTrack );
		}
	    }
	}
    }
  }



  // clear buffer vector
  allMergedTracks.clear();
  
  // merging splitted TPC segments
  if (_mergeTPCSegments) {
    std::vector<GroupTracks*> TPCSegments;
    TPCSegments.clear();
    int nNonAssignedTPCSeg = 0;
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExt = _allTPCTracks[i];
      GroupTracks * group = trkExt->getGroupTracks();
      if (group == NULL) {
	TrackerHitExtendedVec currentVec = trkExt->getTrackerHitExtendedVec();
	int nCur = int(currentVec.size());
	float zmin = 1e+20;
	float zmax = -1e+20;
	for (int iCur=0;iCur<nCur;++iCur) {
	  TrackerHitExtended * curTrkHitExt = currentVec[iCur];
	  TrackerHit * curTrkHit = curTrkHitExt->getTrackerHit();
	  float zpos = float(curTrkHit->getPosition()[2]);
	  if (zpos < zmin)
	    zmin = zpos;
	  if (zpos > zmax)
	    zmax = zpos;
	}
	nNonAssignedTPCSeg++;
	int nGroups = int(TPCSegments.size());
	float dPtMin = 1.0e+10;
	GroupTracks * groupToAttach = NULL;
	TrackExtended * trkToAttach = NULL;
	for (int iG=0;iG<nGroups;++iG) {
	  GroupTracks * segments = TPCSegments[iG];
	  TrackExtendedVec segVec = segments->getTrackExtendedVec();
	  int nTrk = int(segVec.size());
	  bool consider = true;
	  if (_forbidOverlapInZTPC==1) { // if overlap in Z of the two segments is forbidden
	    for (int iTrk=0;iTrk<nTrk;++iTrk) {
	      TrackExtended * trkInGroup = segVec[iTrk];
	      TrackerHitExtendedVec hitInGroupVec = trkInGroup->getTrackerHitExtendedVec();
	      int nHitsInGrp = int(hitInGroupVec.size());
	      for (int iHitInGrp=0;iHitInGrp<nHitsInGrp;iHitInGrp++) {
		TrackerHitExtended * xTrkExt = hitInGroupVec[iHitInGrp];
		TrackerHit * xTrk = xTrkExt->getTrackerHit();
		float xZ = float(xTrk->getPosition()[2]);
		if (xZ>zmin&&xZ<zmax) {
		  consider = false;
		  break;
		}
	      }
	      if (!consider)
		break;
	    }
	  }
	  if (consider) {
	    for (int iTrk=0;iTrk<nTrk;++iTrk) {
	      TrackExtended * trkInGroup = segVec[iTrk];
	      int iComp = 1;
	      float dPt = CompareTrk(trkExt,trkInGroup,_d0CutToMergeTPC,_z0CutToMergeTPC,iComp);
	      if (dPt < dPtMin && !VetoMerge(trkExt,trkInGroup)) {
		dPtMin = dPt;
		groupToAttach = segments;
		trkToAttach = trkInGroup;
		if (_debug==3) {
		  int iopt = 6;
		  PrintOutMerging(trkExt,trkInGroup,iopt);
		}
	      }
	      else {
		if (_debug==3) {
		  int iopt = 9;
		  PrintOutMerging(trkExt,trkInGroup,iopt);
		}
	      }
	    }
	  }
	  else {
	    if (_debug == 3) {
	      int iopt = 9;
	      for (int iTrk=0;iTrk<nTrk;++iTrk) {
		TrackExtended * trkInGroup = segVec[iTrk];
		int iComp = 1;
		float dPt = CompareTrk(trkExt,trkInGroup,_d0CutToMergeTPC,_z0CutToMergeTPC,iComp);		
		if (dPt >= dPtMin) {	      
		  PrintOutMerging(trkExt,trkInGroup,iopt);
		}
	      }
	    }
	  }
	}
	if (dPtMin < _dPCutToMergeTPC && groupToAttach != NULL) {

	   

	  groupToAttach->addTrackExtended(trkExt);
	  trkExt->setGroupTracks(groupToAttach);
	  float zminGroup = groupToAttach->getEdges()[0]; 
	  float zmaxGroup = groupToAttach->getEdges()[1];
	  float edges[2];
	  edges[0] = zmin;
	  if (zminGroup<zmin)
	    edges[0] = zminGroup;
	  edges[1] = zmax;
	  if (zmaxGroup>zmax)
	    edges[1] = zmaxGroup;
	  groupToAttach->setEdges(edges);
	  if (_debug==3) {
	    int iopt = 3;
	    PrintOutMerging(trkExt,trkToAttach,iopt);
	  }
	}
	else {
	  GroupTracks * newSegment = new GroupTracks(trkExt);
	  trkExt->setGroupTracks(newSegment);
	  TPCSegments.push_back(newSegment);
	  float edges[2];
	  edges[0] = zmin;
	  edges[1] = zmax;
	  newSegment->setEdges(edges);
	}
      }
    }
  
    // combining splitted TPC segments with the reconstructed tracks having Si hits
    int nCombTrk = int(_trkImplVec.size());
    int nSegments = int(TPCSegments.size());
    //    std::cout << "Combined tracks = " << nCombTrk << std::endl;
    //    std::cout << "nSegments = " << nSegments << std::endl;
    for (int iS=0;iS<nSegments;++iS) {
	GroupTracks * segments = TPCSegments[iS];
	TrackExtendedVec segVec = segments->getTrackExtendedVec();
	float zminTPCSeg = segments->getEdges()[0];
	float zmaxTPCSeg = segments->getEdges()[1];
	int nTrk = int(segVec.size());
	TrackExtended * CombTrkToAttach = NULL;
	TrackExtended * keyTrack = NULL;
	float deltaPtMin = _dPCutToMergeTPC;
	for (int iCTrk=0;iCTrk<nCombTrk;++iCTrk) {
	  TrackExtended * combTrk = _trkImplVec[iCTrk];
	  GroupTracks * groupComb = combTrk->getGroupTracks();
	  bool consider = true;
	  if (_forbidOverlapInZComb==1) { // if overlap in Z of the two segments is forbidden
	    float zminComb = groupComb->getEdges()[0];
	    float zmaxComb = groupComb->getEdges()[1];
	    consider = (zminTPCSeg>zmaxComb) || (zmaxTPCSeg<zminComb);
	  }
	  if (consider) {
	    for (int iTrk=0;iTrk<nTrk;++iTrk) {
	      TrackExtended * trk = segVec[iTrk];
	      int iopt = 0;
	      float dPt = CompareTrk(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt);
	      float angleSignificance(0.);
	      float significance = CompareTrkIII(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt,angleSignificance);
	      if ( (dPt<deltaPtMin || significance <5 ) ) {
		if(VetoMerge(trk,combTrk)==false){
		  CombTrkToAttach = combTrk;
		  keyTrack = trk;
		  deltaPtMin = dPt;
		}
	      }
	      else {
		if (_debug==3) {
		  GroupTracks * groupCur = combTrk->getGroupTracks();
		  TrackExtended * dummySi = groupCur->getTrackExtendedVec()[0];
		  int iopt = 8;
		  PrintOutMerging(trk,dummySi,iopt);
		}
	      }
	    }
	  }
	  else {
	    if (_debug==3) {
	      for (int iTrk=0;iTrk<nTrk;++iTrk) {
		TrackExtended * trk = segVec[iTrk];
		int iopt = 0;
		float dPt = CompareTrk(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt);
		if (dPt>deltaPtMin) {
		  GroupTracks * groupCur = combTrk->getGroupTracks();
		  TrackExtended * dummySi = groupCur->getTrackExtendedVec()[0];
		  int iopt = 8;
		  PrintOutMerging(trk,dummySi,iopt);
		}
	      }
	    }
	  }
	}
	
	if (CombTrkToAttach != NULL) { // attach TPC segment to existing Comb Track
	  GroupTracks * groupToAttach = CombTrkToAttach->getGroupTracks();	  
	  TrackExtended * SiCombTrk = groupToAttach->getTrackExtendedVec()[0];
	  TrackExtended * TpcCombTrk = groupToAttach->getTrackExtendedVec()[1];
	  if (_debug==3) {
	    int iopt = 4;
	    PrintOutMerging(keyTrack,SiCombTrk,iopt);
	    iopt = 5;
	    PrintOutMerging(keyTrack,TpcCombTrk,iopt);	  
	  }
	  for (int iTrk=0;iTrk<nTrk;iTrk++) {
	    TrackExtended * trk = segVec[iTrk];
	    groupToAttach->addTrackExtended( trk );
	    trk->setGroupTracks( groupToAttach );
	    TrackerHitExtendedVec hitVec = trk->getTrackerHitExtendedVec();
	    int nHitS = int(hitVec.size());	       
	    for (int iHS=0;iHS<nHitS;++iHS) {
	      TrackerHitExtended * hitExt = hitVec[iHS];
	      hitExt->setUsedInFit(false);
	      CombTrkToAttach->addTrackerHitExtended( hitExt );
	    }
	  }
	}
	else {
	  if (nTrk==1) { //
	    GroupTracks * newGrp = new GroupTracks();
	    segVec[0]->setGroupTracks(newGrp);
	    newGrp->addTrackExtended(segVec[0]);
	    TrackerHitExtendedVec TpcHitVec = segVec[0]->getTrackerHitExtendedVec();
	    int nTpcH = int(TpcHitVec.size());
	    for (int iTpcH=0;iTpcH<nTpcH;++iTpcH) {
	      TpcHitVec[iTpcH]->setUsedInFit( true );
	    }
	    _trkImplVec.push_back(segVec[0]);
	    _allNonCombinedTPCTracks.push_back(segVec[0]);
	  }
	  else {
	    float zMin = 1.0e+20;
	    TrackExtended * chosenTrack = NULL;
	    for (int iTrk=0;iTrk<nTrk;++iTrk) {
	      TrackExtended * trk = segVec[iTrk];
	      Track * track = trk->getTrack();
	      TrackerHitVec hitVec = track->getTrackerHits();
	      int nHits = int(hitVec.size());
	      for (int iH=0;iH<nHits;++iH) {
		float zPosi = fabs(hitVec[iH]->getPosition()[2]);
		if (zPosi<zMin) {
		  chosenTrack = trk;
		  zMin = zPosi;
		  break;
		}
	      }
	    }
	    if (chosenTrack!=NULL) {
	      GroupTracks * newGroup = new GroupTracks();
	      for (int iTrk=0;iTrk<nTrk;++iTrk) {
		TrackExtended * trk = segVec[iTrk];
		trk->setGroupTracks( newGroup );
		newGroup->addTrackExtended( trk );			
		TrackerHitExtendedVec hitVecS = trk->getTrackerHitExtendedVec();
		int nHitS = int(hitVecS.size());			
		for (int iH=0;iH<nHitS;++iH) {
		  TrackerHitExtended * trkHitExt = hitVecS[iH];
		  if (trk!=chosenTrack) {
		    trkHitExt->setUsedInFit( false );
		    chosenTrack->addTrackerHitExtended( trkHitExt );				
		  }
		  else {
		    trkHitExt->setUsedInFit( true );
		  }
		}
	      }
	      _allNonCombinedTPCTracks.push_back(chosenTrack);
	      _trkImplVec.push_back(chosenTrack);
	    }
	  }
	}
    }
    for (int iS=0;iS<nSegments;++iS) {
      GroupTracks * segments = TPCSegments[iS];
      delete segments;
    }
    TPCSegments.clear();
  }
  else { // adding all TPC segments to the list of tracks (track splitting is allowed)
      for (int i=0;i<nTPCTrk;++i) {
	  TrackExtended * trkExt = _allTPCTracks[i];
	  Track * track = trkExt->getTrack();
	  GroupTracks * group = trkExt->getGroupTracks();
	  if (group == NULL) {
	      TrackerHitVec hitVec = track->getTrackerHits();       
	      _trkImplVec.push_back(trkExt);
	      _allNonCombinedTPCTracks.push_back( trkExt );
	      GroupTracks * newGrp = new GroupTracks();
	      newGrp->addTrackExtended( trkExt );
	      trkExt->setGroupTracks( newGrp );
	  }
      }    
  }
  
  for (int i=0;i<nSiTrk;++i) { // adding left-over Si segments to the list of tracks
      TrackExtended * trkExt = _allSiTracks[i];
      GroupTracks * group = trkExt->getGroupTracks();
      if (group == NULL) {
	  TrackerHitExtendedVec hitVec = trkExt->getTrackerHitExtendedVec();
	  int nHSi = int(hitVec.size());
	  for (int iHSi=0;iHSi<nHSi;++iHSi) {
	    hitVec[iHSi]->setUsedInFit(true);
	  }
	  _trkImplVec.push_back(trkExt);
	  GroupTracks * newGrp = new GroupTracks();
	  newGrp->addTrackExtended( trkExt );
	  trkExt->setGroupTracks( newGrp );	  
	  _allNonCombinedSiTracks.push_back( trkExt );
      }
  }
  
}

void FullLDCTracking::CheckTracks() {  

  for(unsigned int i = 0; i< _trkImplVec.size();i++){
    TrackExtended *first = _trkImplVec[i];
    if(first==NULL)continue;
    float d0First = first->getD0();
    float z0First = first->getZ0();
    float omegaFirst = first->getOmega();
    float tanLFirst = first->getTanLambda();
    float phiFirst = first->getPhi();
    HelixClass helixFirst;
    helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
    float momFirst[3];
    momFirst[0]= helixFirst.getMomentum()[0];
    momFirst[1]= helixFirst.getMomentum()[1];
    momFirst[2]= helixFirst.getMomentum()[2];
    float pFirst    = sqrt(momFirst[0]*momFirst[0]+momFirst[1]*momFirst[1]+momFirst[2]*momFirst[2]);
    if(isnan(pFirst))continue;
    TrackerHitExtendedVec firstHitVec  = first->getTrackerHitExtendedVec();
    if(firstHitVec.size()<1)continue;
    
    for(unsigned int j = i+1; j<_trkImplVec.size();j++){
      TrackExtended *second = _trkImplVec[j];
      if(second==NULL)continue;
      float d0Second = second->getD0();
      float z0Second = second->getZ0();
      float omegaSecond = second->getOmega();
      float tanLSecond = second->getTanLambda();
      float phiSecond = second->getPhi();
      HelixClass helixSecond;
      helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
      float momSecond[3];
      momSecond[0] = helixSecond.getMomentum()[0];
      momSecond[1] = helixSecond.getMomentum()[1];
      momSecond[2] = helixSecond.getMomentum()[2];
      float pSecond    = sqrt(momSecond[0]*momSecond[0]+momSecond[1]*momSecond[1]+momSecond[2]*momSecond[2]);
      if(isnan(pSecond))continue;
      TrackerHitExtendedVec secondHitVec  = second->getTrackerHitExtendedVec();
      if(secondHitVec.size()<1)continue;
      if(firstHitVec.size()+secondHitVec.size()<10)continue;


      float pdot = (momFirst[0]*momSecond[0]+momFirst[1]*momSecond[1]+momFirst[2]*momSecond[2])/pFirst/pSecond;
      if(pdot<0.999)continue;
      // const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
      // const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
      // const float deltaP = fabs(pFirst-pSecond);
      // const float sigmaDeltaP = sqrt(pFirst*sigmaPOverPFirst*pFirst*sigmaPOverPFirst+pSecond*sigmaPOverPSecond*pSecond*sigmaPOverPSecond);
      //      const float significance = deltaP/sigmaDeltaP;
      
      TrackExtended * combinedTrack = TrialCombineTracks(first,second);
      if(combinedTrack != NULL){
	const int minHits = std::min(firstHitVec.size(),secondHitVec.size());
	const int maxHits = std::max(firstHitVec.size(),secondHitVec.size());
       
	if( combinedTrack->getNDF() <= 2*maxHits+minHits-5){
	  delete combinedTrack->getGroupTracks();
	  delete combinedTrack;
	  continue;
	}

	float d0    = combinedTrack->getD0();
	float z0    = combinedTrack->getZ0();
	float omega = combinedTrack->getOmega();
	float tanL  = combinedTrack->getTanLambda();
	float phi   = combinedTrack->getPhi();

	HelixClass helix;
	helix.Initialize_Canonical(phi,d0,z0,omega,tanL,_bField);
	float mom[3];
	mom[0]  = helix.getMomentum()[0];
	mom[1]  = helix.getMomentum()[1];
	mom[2]  = helix.getMomentum()[2];
	// float p = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
	// float chi2Sig =  (combinedTrack->getChi2() - combinedTrack->getNDF());
	// chi2Sig = chi2Sig/sqrt(combinedTrack->getNDF()*2);



	int nTpcFirst(0);
	int nUsedFirst(0);
	for(unsigned int ihit = 0;ihit<firstHitVec.size();ihit++){
	  if(firstHitVec[ihit]->getDet()==1)nTpcFirst++;
	  if(firstHitVec[ihit]->getUsedInFit()==true)nUsedFirst++;
	}


	int nTpcSecond(0);
	int nUsedSecond(0);
	for(unsigned int ihit = 0;ihit<secondHitVec.size();ihit++){
	  if(secondHitVec[ihit]->getDet()==1)nTpcSecond++;
	  if(secondHitVec[ihit]->getUsedInFit()==true)nUsedSecond++;
	}
	delete combinedTrack->getGroupTracks();
	delete combinedTrack;
      }
    }
  }
  
  
}





float FullLDCTracking::CompareTrkII(TrackExtended * first, TrackExtended * second, 
				    float d0Cut, float z0Cut,int iopt,float & Angle) {


  float result = 1.0e+20;
  Angle  = 1.0e+20; 
  float omegaFirst = first->getOmega();
  float omegaSecond = second->getOmega();
  float deltaOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
  if(deltaOmega> 2*_dOmegaForMerging)return result;


  float d0First = first->getD0();
  float z0First = first->getZ0();
  float d0Second = second->getD0();
  float z0Second = second->getZ0();

  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);


  if(!isCloseInIP)return result;

  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  float qFirst = PIOVER2 - atan(tanLFirst);
  float qSecond = PIOVER2 - atan(tanLSecond);

  Angle = (cos(phiFirst)*cos(phiSecond)+sin(phiFirst)*sin(phiSecond))*
             sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
  Angle = acos(Angle);

  result = deltaOmega;
  
  return result;

}




float FullLDCTracking::CompareTrkIII(TrackExtended * first, TrackExtended * second, 
				    float d0Cut, float z0Cut,int iopt, float & AngleSignificance) {


  float result = 1.0e+20;

  float d0First = first->getD0();
  float z0First = first->getZ0();
  float omegaFirst = first->getOmega();
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  float qFirst = PIOVER2 - atan(tanLFirst);


  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  float omegaSecond = second->getOmega();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  float qSecond = PIOVER2 - atan(tanLSecond);


  //MB 2010 03
  float d0ErrFirst = sqrt(first->getCovMatrix()[0]);
  float z0ErrFirst = sqrt(first->getCovMatrix()[9]);
  //  float omegaErrFirst = sqrt(first->getCovMatrix()[5]);
  float phiErrFirst = sqrt(first->getCovMatrix()[2]);
  float qErrFirst = sqrt(cos(qFirst)*cos(qFirst)*first->getCovMatrix()[14]);
  //MB END
  //MB 2010 03
  float d0ErrSecond = sqrt(second->getCovMatrix()[0]);
  float z0ErrSecond = sqrt(second->getCovMatrix()[9]);
  //  float omegaErrSecond = sqrt(second->getCovMatrix()[5]);
  float phiErrSecond = sqrt(second->getCovMatrix()[2]);
  float qErrSecond = sqrt(cos(qSecond)*cos(qSecond)*second->getCovMatrix()[14]);
  //MB END


  //  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  //isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);

  //MB 2010 03
  bool isCloseInIP = (fabs(d0First-d0Second)/sqrt(d0ErrFirst*d0ErrFirst+d0ErrSecond*d0ErrSecond)<d0Cut);
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)/sqrt(z0ErrFirst*z0ErrFirst+z0ErrSecond*z0ErrSecond)<z0Cut);

  if (!isCloseInIP)return result;

  float Angle = (cos(phiFirst)*cos(phiSecond)+sin(phiFirst)*sin(phiSecond))*
    sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);


  HelixClass helixFirst;
  helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
  HelixClass helixSecond;
  helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);

  float pFirst[3];
  float pSecond[3];
  float momFirst = 0;
  float momSecond = 0;
  
  for (int iC=0;iC<3;++iC) {
    pFirst[iC] = helixFirst.getMomentum()[iC];
    pSecond[iC] = helixSecond.getMomentum()[iC];
    momFirst += pFirst[iC]* pFirst[iC];
    momSecond += pSecond[iC]*pSecond[iC];
  }
  momFirst = sqrt(momFirst);
  momSecond = sqrt(momSecond);

  
  float pdot = (pFirst[0]*pSecond[0]+pFirst[1]*pSecond[1]+pFirst[2]*pSecond[2])/momFirst/momSecond;


  const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
  const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
  const float deltaP = fabs(momFirst-momSecond);
  const float sigmaPFirst = momFirst*sigmaPOverPFirst;
  const float sigmaPSecond = momSecond*sigmaPOverPSecond;
  const float sigmaDeltaP = sqrt(sigmaPFirst*sigmaPFirst+sigmaPSecond*sigmaPSecond);
  const float significance = deltaP/sigmaDeltaP;

  //MB 2010 03
    float errorAngle =sin(phiFirst)*sin(phiFirst)*phiErrFirst*phiErrFirst*cos(phiSecond)*cos(phiSecond)+
      sin(phiSecond)*sin(phiSecond)*phiErrSecond*phiErrSecond*cos(phiFirst)*cos(phiFirst)+
      sin(qFirst)*sin(qFirst)*qErrFirst*qErrFirst*cos(qSecond)*cos(qSecond)+
      sin(qSecond)*sin(qSecond)*qErrSecond*qErrSecond*cos(qFirst)*cos(qFirst)+
      cos(phiFirst)*cos(phiFirst)*phiErrFirst*phiErrFirst*(sin(phiSecond)*sin(qFirst)*sin(qSecond))*(sin(phiSecond)*sin(qFirst)*sin(qSecond))+
      cos(phiSecond)*cos(phiSecond)*phiErrSecond*phiErrSecond*(sin(phiFirst)*sin(qFirst)*sin(qSecond))*(sin(phiFirst)*sin(qFirst)*sin(qSecond))+
      cos(qFirst)*cos(qFirst)*qErrFirst*qErrFirst*(sin(phiFirst)*sin(phiSecond)*sin(qSecond))*(sin(phiFirst)*sin(phiSecond)*sin(qSecond))+
      cos(qSecond)*cos(qSecond)*qErrSecond*qErrSecond*(sin(phiFirst)*sin(phiSecond)*sin(qFirst))*(sin(phiFirst)*sin(phiSecond)*sin(qFirst));
    
    if(Angle<1.){
      errorAngle = sqrt(1./(1.-Angle*Angle)*errorAngle);
    }else{
      errorAngle = sqrt(errorAngle);
    }
    
    if(errorAngle<1.e-6)errorAngle=1.e-6;

    AngleSignificance = fabs(acos(Angle)/errorAngle);
    if(pdot<0.999)return result;
        
    return significance;
   
}


float FullLDCTracking::CompareTrk(TrackExtended * first, TrackExtended * second, 
				  float d0Cut, float z0Cut,int iopt) {
  
  float result = 1.0e+20;
  
  float d0First = first->getD0();
  float z0First = first->getZ0();
  float omegaFirst = first->getOmega();
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  
  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  float omegaSecond = second->getOmega();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  
  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  
  if (iopt>0) isCloseInIP = isCloseInIP || (fabs(d0First+d0Second)<d0Cut);
  
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);
  
  
  HelixClass helixFirst;
  helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
  HelixClass helixSecond;
  helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
  
  float pFirst[3];
  float pSecond[3];
  float dPminus[3];
  float dPplus[3];
  float momFirst = 0;
  float momSecond = 0;
  float momMinus = 0;
  float momPlus = 0;

  if ( isCloseInIP ) {

    for (int iC=0;iC<3;++iC) {
      pFirst[iC] = helixFirst.getMomentum()[iC];
      pSecond[iC] = helixSecond.getMomentum()[iC];
      momFirst += pFirst[iC]* pFirst[iC];
      momSecond += pSecond[iC]*pSecond[iC];
      dPminus[iC] = pFirst[iC] - pSecond[iC];
      dPplus[iC] = pFirst[iC] + pSecond[iC];
      momMinus += dPminus[iC]*dPminus[iC];
      momPlus += dPplus[iC]*dPplus[iC];
    }
    momFirst = sqrt(momFirst);
    momSecond = sqrt(momSecond);
    
    float ptFirst = sqrt(pFirst[0]*pFirst[0]+pFirst[1]*pFirst[1]);
    float ptSecond = sqrt(pSecond[0]*pSecond[0]+pSecond[1]*pSecond[1]);


    if ( (ptFirst<_PtCutToMergeTPC) && (ptSecond<_PtCutToMergeTPC) ) {
      
      momMinus = sqrt(momMinus);
      momPlus = sqrt(momPlus);
      float nom = momMinus;
      if (momPlus<nom && iopt>0)
	nom = momPlus;
      float den = momFirst;
      if (momSecond<momFirst)
	den = momSecond;
      
      result = nom/den;     

    }

    else {

      
    // check for cases where PatRec splits non-looping TPC tracks 
    // look for two tracks where total tpc hits are not more than total number
    // of pad rows and that the hits on one track are close to the helix of the
    // other track

      float dpOverP = 2.0*fabs(momFirst-momSecond)/(momFirst+momSecond);
      const float pdot = (pFirst[0]*pSecond[0]+pFirst[1]*pSecond[1]+pFirst[2]*pSecond[2])/momFirst/momSecond;
      const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
      const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
      const float deltaP = fabs(momFirst-momSecond);
      const float sigmaPFirst = momFirst*sigmaPOverPFirst;
      const float sigmaPSecond = momSecond*sigmaPOverPSecond;
      const float sigmaDeltaP = sqrt(sigmaPFirst*sigmaPFirst+sigmaPSecond*sigmaPSecond);
      const float significance = deltaP/sigmaDeltaP;


    //compare angle between the two vectors (cos theta) and their momentum
    //    if( ( pdot>0.99 && dpOverP<0.01 ) || ( pdot>0.998 && dpOverP<0.25 ) ){
      if( (pdot>_cosThetaCutHighPtMerge && dpOverP<_momDiffCutHighPtMerge) 
	  || 
	  (pdot>_cosThetaCutSoftHighPtMerge && dpOverP<_momDiffCutSoftHighPtMerge) 
	  || (pdot > 0.9999 && significance <10) 
        ){

	
	int nTrkGrpFirst = 0;
	int nTrkGrpSecond = 0;
	TrackerHitVec hitvecFirst;
	TrackerHitVec hitvecSecond;
	GroupTracks * groupFirst = first->getGroupTracks();
	GroupTracks * groupSecond = second->getGroupTracks();
	
	if(groupFirst!=NULL){
	
	  TrackExtendedVec tracksInGroupFirst = groupFirst->getTrackExtendedVec();
	  nTrkGrpFirst = int(tracksInGroupFirst.size());
	
	  for (int iTrkGrp=0;iTrkGrp<nTrkGrpFirst;++iTrkGrp) {
	
	    TrackExtended * trkGrp = tracksInGroupFirst[iTrkGrp];
	    TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();
	  
	    for(unsigned int i =0; i<hitVec.size(); ++i){
	      hitvecFirst.push_back(hitVec[i]->getTrackerHit());	  
	    }
	  }
	}
      
	if(groupSecond!=NULL){

	  TrackExtendedVec tracksInGroupSecond = groupSecond->getTrackExtendedVec();
	  nTrkGrpSecond = int(tracksInGroupSecond.size());
	
	  for (int iTrkGrp=0;iTrkGrp<nTrkGrpSecond;++iTrkGrp) {
	    TrackExtended * trkGrp = tracksInGroupSecond[iTrkGrp];
	    TrackerHitExtendedVec hitVec = 
	      trkGrp->getTrackerHitExtendedVec();
	    
	    for(unsigned int i=0;i<hitVec.size();++i){
	      hitvecSecond.push_back(hitVec[i]->getTrackerHit());
	    }
	  }
	}
	
	const gear::TPCParameters& pTPC = Global::GEAR->getTPCParameters();
	const gear::PadRowLayout2D& pTPCpads = pTPC.getPadLayout();
	
	std::vector<double>tpcExtent = pTPCpads.getPlaneExtent();
      
	float tpcInnerR = tpcExtent[0]; 
	int tpcMaxRow = pTPCpads.getNRows();
      
	// for non-looping tracks 
	int nhitsFirst  = (int)hitvecFirst.size();
	int nhitsSecond = (int)hitvecSecond.size();
	int ntpcFirst   = 0;
	int ntpcSecond  = 0;
	float hitxyz[3];
	float dist[3];
	float maxdistFirst=0;
	float maxdistSecond=0;
	int ncloseFirst = 0;
	int ncloseSecond = 0;
	float zminFirst = 99999;
	float zminSecond = 99999;
	float zmaxFirst = -99999;
	float zmaxSecond = -99999;
	
	for(int ih =0;ih<nhitsFirst;++ih){
	
	  float x = (float) hitvecFirst[ih]->getPosition()[0];
	  float y = (float) hitvecFirst[ih]->getPosition()[1];
	  float z = (float) hitvecFirst[ih]->getPosition()[2];
	  
	  if(fabs(z)<zminFirst) zminFirst=fabs(z);
	  if(fabs(z)>zmaxFirst) zmaxFirst=fabs(z);
	  
	  float r = sqrt(x*x+y*y);
	  
	  if(r>tpcInnerR) ntpcFirst++;

	  hitxyz[0]=x;
	  hitxyz[1]=y;
	  hitxyz[2]=z;
	  helixSecond.getDistanceToPoint(hitxyz, dist);
	  
	  // compare 3D distance between hit and extrapolation
	  if(dist[2]>maxdistFirst) maxdistFirst=dist[2];
	  if(dist[2]<_hitDistanceCutHighPtMerge) ncloseFirst++;
	}
	
	for(int ih =0;ih<nhitsSecond;++ih){
	  
	  float x = (float)hitvecSecond[ih]->getPosition()[0];
	  float y = (float)hitvecSecond[ih]->getPosition()[1];
	  float z = (float)hitvecSecond[ih]->getPosition()[2];
	  
	  if(fabs(z)<zminSecond) zminSecond=fabs(z);
	  if(fabs(z)>zmaxSecond) zmaxSecond=fabs(z);
	  
	  float r = sqrt(x*x+y*y);
	  
	  if(r>tpcInnerR) ntpcSecond++;

	  hitxyz[0]=x;
	  hitxyz[1]=y;
	  hitxyz[2]=z;
	  helixFirst.getDistanceToPoint(hitxyz, dist);
	  
	  // compare 3D distance between hit and extrapolation
	  if(dist[2]>maxdistSecond) maxdistSecond=dist[2];
	  if(dist[2]<_hitDistanceCutHighPtMerge) ncloseSecond++;
	}
	
	float fcloseFirst  = (float)ncloseFirst/(float)nhitsFirst;
	float fcloseSecond = (float)ncloseSecond/(float)nhitsSecond;



	bool split = false;
	//std::cout << "Momenta = " << momFirst << " " << momSecond << std::endl;
 	//std::cout << "MaxDist = " << maxdistSecond << " " << maxdistFirst << " " << _maxHitDistanceCutHighPtMerge << std::endl;
 	//std::cout << "close   = " << fcloseSecond << " " << fcloseFirst << " " << _maxFractionOfOutliersCutHighPtMerge << std::endl;
 	//std::cout << "ntpc    = " << ntpcFirst << " " << ntpcSecond << " " << tpcMaxRow+10 << std::endl;
 	//std::cout << "pdot    = " << pdot << " significance " << significance << std::endl;
	
	TrackExtended * combinedTrack = TrialCombineTracks(first,second);
	
	if(combinedTrack != NULL){
	  //std::cout << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << first->getNDF()+second->getNDF()+5 << std::endl;
	  if(combinedTrack->getNDF()+10>first->getNDF()+second->getNDF()+5){
	    split = true;
	    dpOverP = 0;
	    //std::cout << " Forcing MERGE " << std::endl;
	  }
	  delete combinedTrack->getGroupTracks();
	  delete combinedTrack;
	}else{
	  //std::cout << "Could not combine track " << std::endl;
	  if(significance<5 && fcloseFirst>_maxFractionOfOutliersCutHighPtMerge){
	    split = true;
	    dpOverP = 0;
	    //	    int overlap = SegmentRadialOverlap(first,second);
	    //std::cout << " Forcing MERGE " << overlap << std::endl;
	  }
	}

	// criteria for split track
	// old criterion
	if( maxdistSecond < _maxHitDistanceCutHighPtMerge && maxdistFirst < _maxHitDistanceCutHighPtMerge 
	    && 
	    (fcloseSecond > _maxFractionOfOutliersCutHighPtMerge || fcloseFirst > _maxFractionOfOutliersCutHighPtMerge) 
	    && 
	    ntpcFirst+ntpcSecond < tpcMaxRow+10.)
	  split = true;
	
	if(split){
	  result = dpOverP;
	}
	
      }
    }
  }

  return result;
  
}

void FullLDCTracking::AddNotAssignedHits() {

    const gear::VXDParameters& pVXDDet = Global::GEAR->getVXDParameters();
    const gear::VXDLayerLayout& pVXDLayerLayout = pVXDDet.getVXDLayerLayout();
    int nLayersVTX = int(pVXDLayerLayout.getNLayers());

    const gear::GearParameters& pFTDDet = Global::GEAR->getGearParameters("FTD");  
    int nLayersFTD = int(pFTDDet.getDoubleVals("FTDZCoordinate").size());

    const gear::GearParameters& pSITDet = Global::GEAR->getGearParameters("SIT");  
    int nLayersSIT = int(pSITDet.getDoubleVals("SITLayerRadius").size());

    // Creating helix extrapolations
    if (_assignSETHits>0||_assignETDHits>0)
      CreateExtrapolations();

    if (_assignSETHits>0) { // Assignment of SET Hits

      const gear::GearParameters& pSETDet = Global::GEAR->getGearParameters("SET");  
      int nLayersSET = int(pSETDet.getDoubleVals("SETLayerRadius").size());

      int nSETHits = _allSETHits.size();
      std::vector<TrackerHitExtendedVec> SETHits;
      SETHits.resize(nLayersSET);

      for (int iSET=0;iSET<nSETHits;++iSET) {
	TrackerHitExtended * trkHit = _allSETHits[iSET];
	TrackerHit * hit = trkHit->getTrackerHit();
	int layer = hit->getType() - 401;
	if (layer>=0&&layer<nLayersSET) 
	  SETHits[layer].push_back(trkHit);	
      }
      for (int iL=0; iL<nLayersSET; ++iL) { // loop over SET layers
	TrackerHitExtendedVec hitVec = SETHits[iL];
	int refit = 1;
	AssignOuterHitsToTracks(hitVec,_distCutForSETHits,refit);
      }
    }

    if (_assignETDHits>0) { // Assignment of ETD Hits
      
      const gear::GearParameters& pETDDet = Global::GEAR->getGearParameters("ETD");  
      int nLayersETD = int(pETDDet.getDoubleVals("ETDLayerZ").size());

      int nETDHits = _allETDHits.size();
      std::vector<TrackerHitExtendedVec> ETDHits;
      ETDHits.resize(nLayersETD);

      for (int iETD=0;iETD<nETDHits;++iETD) {
	TrackerHitExtended * trkHit = _allETDHits[iETD];
	TrackerHit * hit = trkHit->getTrackerHit();
	int layer = hit->getType() - 201;
	if (layer>=0 && layer < nLayersETD) 
	  ETDHits[layer].push_back(trkHit);
      }
      for (int iL=0; iL<nLayersETD; ++iL) { // loop over ETD layers
	TrackerHitExtendedVec hitVec = ETDHits[iL];
	int refit = 0;
	AssignOuterHitsToTracks( hitVec, _distCutForETDHits, refit );
      }

    }
    
    // Cleaning up helix extrapolations
    if (_assignSETHits>0||_assignETDHits>0)
      CleanUpExtrapolations();

    if (_assignSITHits>0) { // Treatment of left-over SIT hits 
	std::vector<TrackerHitExtendedVec> nonAssignedSITHits;    
	nonAssignedSITHits.resize(nLayersSIT);
	int nSITHits = int(_allSITHits.size());
	for (int iH=0;iH<nSITHits;++iH) {
	    TrackerHitExtended * trkHitExt = _allSITHits[iH];
	    TrackExtended * trkExt = trkHitExt->getTrackExtended();
	    if (trkExt == NULL) {
		TrackerHit * trkHit = trkHitExt->getTrackerHit();
		int layer = trkHit->getType() - 401;
		if (layer >=0 && layer < nLayersSIT) {
		    nonAssignedSITHits[layer].push_back(trkHitExt);
		}
	    }
	}	
	for (int iL=nLayersSIT-1;iL>=0;--iL) { // loop over layers in Si
	    TrackerHitExtendedVec hitVec = nonAssignedSITHits[iL];
	    AssignSiHitsToTracks(hitVec,
				 _distCutForSITHits);
	}
    }
    
    if (_assignFTDHits>0) { // Treatment of left-over FTD hits
	std::vector<TrackerHitExtendedVec> nonAssignedFTDHits;
	nonAssignedFTDHits.resize(nLayersFTD);
	int nFTDHits = int(_allFTDHits.size());
	for (int iH=0;iH<nFTDHits;++iH) {
	    TrackerHitExtended * trkHitExt = _allFTDHits[iH];
	    TrackExtended * trkExt = trkHitExt->getTrackExtended();
	    if (trkExt == NULL) {
		TrackerHit * trkHit = trkHitExt->getTrackerHit();
		int layer = trkHit->getType() - 201;
		if (layer >=0 && layer < nLayersFTD)
		    nonAssignedFTDHits[layer].push_back(trkHitExt);
	    }
	}
	for (int iL=nLayersFTD-1;iL>=0;--iL) {
	    TrackerHitExtendedVec hitVec = nonAssignedFTDHits[iL];
	    AssignSiHitsToTracks(hitVec,
				 _distCutForFTDHits);	
	}
    }
    

    if (_assignVTXHits>0) { // Treatment of left-over VTX hits
	std::vector<TrackerHitExtendedVec> nonAssignedVTXHits;
	nonAssignedVTXHits.resize(nLayersVTX);
	int nVTXHits = int(_allVTXHits.size());
	for (int iH=0;iH<nVTXHits;++iH) {
	    TrackerHitExtended * trkHitExt = _allVTXHits[iH];
	    TrackExtended * trkExt = trkHitExt->getTrackExtended();
	    if (trkExt == NULL) {
		TrackerHit * trkHit = trkHitExt->getTrackerHit();
		int layer = trkHit->getType() - 101;
		if (layer >=0 && layer < nLayersVTX)
		    nonAssignedVTXHits[layer].push_back(trkHitExt);
	    }
	}
	for (int iL=nLayersVTX-1;iL>=0;--iL) {
	    TrackerHitExtendedVec hitVec = nonAssignedVTXHits[iL];
	    AssignSiHitsToTracks(hitVec,
				 _distCutForVTXHits);	
	}
    }
    
    if(_debug>1)std::cout << "Assign TPC hits *********************************" << std::endl;

    if (_assignTPCHits) {// Treatment of left-over TPC hits
	TrackerHitExtendedVec nonAssignedTPCHits;
	int nTPCHits = int(_allTPCHits.size());
	for (int iH=0;iH<nTPCHits;++iH) {
	    TrackerHitExtended * trkHitExt = _allTPCHits[iH];
	    TrackExtended * trkExt = trkHitExt->getTrackExtended();
	    if (trkExt == NULL) {
		nonAssignedTPCHits.push_back(trkHitExt);
	    }
	}
	AssignTPCHitsToTracks(nonAssignedTPCHits,
			      _distCutForTPCHits);
    }


}


void FullLDCTracking::CreateExtrapolations() {

  _trackExtrapolatedHelix.clear();

  int nTrk = int(_trkImplVec.size());
  for (int iTrk=0; iTrk<nTrk; ++iTrk) {
    TrackExtended * trk = _trkImplVec[iTrk];
    HelixClass * helix = GetExtrapolationHelix( trk );
    _trackExtrapolatedHelix[trk] = helix;
  }

}

void FullLDCTracking::CleanUpExtrapolations() {

  int nTrk = int(_trkImplVec.size());
  for (int iTrk=0; iTrk<nTrk; ++iTrk) {
    TrackExtended * trk = _trkImplVec[iTrk];
    HelixClass * helix =  _trackExtrapolatedHelix[trk];
    delete helix;
  }  
  _trackExtrapolatedHelix.clear();

}


void FullLDCTracking::AssignOuterHitsToTracks(TrackerHitExtendedVec hitVec, float dcut, int refit) {
  
    int nHits = int(hitVec.size());
    int nTrk = int(_trkImplVec.size());

    std::map <TrackExtended*,bool> flagTrack;
    std::map <TrackerHitExtended*,bool> flagHit;
    TrackHitPairVec pairs;
    flagTrack.clear();
    flagHit.clear();
    pairs.clear();

    for (int iH=0;iH<nHits;++iH) {
      float pos[3];
      TrackerHitExtended * trkHitExt = hitVec[iH];
      TrackerHit * hit = trkHitExt->getTrackerHit();
      for (int ip=0;ip<3;++ip)
	pos[ip] = float(hit->getPosition()[ip]);
      //      float r_hit = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
      //      float phi_hit = atan2(pos[1],pos[0]);
      //      float z_hit = pos[2];
      //      std::cout << r_hit << " " << phi_hit << " "
      // 		<< z_hit << " " 
      // 		<< trkHitExt->getResolutionRPhi() << " "
      // 		<< trkHitExt->getResolutionZ() << std::endl;

      for (int iT=0;iT<nTrk;++iT) {
	TrackExtended * trkExt = _trkImplVec[iT];
	float tanLambda = trkExt->getTanLambda();	    
	float product = pos[2]*tanLambda;
	if (product>0) {
	  HelixClass * helix = _trackExtrapolatedHelix[trkExt];
	  float distance = helix->getDistanceToPoint(pos,dcut);
	  //	  std::cout << "Dist = " << dist[2] << std::endl;
	  if (distance<dcut) {
	    TrackHitPair * trkHitPair = 
	      new TrackHitPair(trkExt,trkHitExt,distance);
	    pairs.push_back(trkHitPair);
	    flagTrack[trkExt] = true;
	    flagHit[trkHitExt] = true;
	  }
	}
      }
    }

    int nPairs = int(pairs.size());
    if (nPairs>0) {
      SortingTrackHitPairs(pairs);
      for (int iP=0;iP<nPairs;++iP) {
	TrackHitPair * trkHitPair = pairs[iP];
	TrackExtended * trkExt = trkHitPair->getTrackExtended();
	TrackerHitExtended * trkHitExt = 
	  trkHitPair->getTrackerHitExtended();	   
	if (flagTrack[trkExt] && flagHit[trkHitExt]) {
	  if (refit==0) {
	    trkExt->addTrackerHitExtended( trkHitExt );
	    trkHitExt->setUsedInFit( false );
	    trkHitExt->setTrackExtended( trkExt );
	  }
	  else {
	    TrackerHitExtendedVec hitsInTrack = 
	      trkExt->getTrackerHitExtendedVec();
	    int nTotH = int(hitsInTrack.size());
	    int nHitsInFit = 0;
	    for (int iTH=0;iTH<nTotH;++iTH) {
	      TrackerHitExtended * hitInTrack = hitsInTrack[iTH];
	      if (hitInTrack->getUsedInFit())
		nHitsInFit++;
	    }
	    float * x_h = new float[nHitsInFit+1];
	    float * y_h = new float[nHitsInFit+1];
	    float * z_h = new float[nHitsInFit+1];
	    int * idet_h = new int[nHitsInFit+1];
	    int * ityp_h = new int[nHitsInFit+1];
	    int * lhits = new int[nHitsInFit+1];
	    float * rR_h = new float[nHitsInFit+1];
	    float * rZ_h = new float[nHitsInFit+1];
	    int iHitInFit = 0;
	    for (int iHit=0;iHit<nTotH;++iHit) {
	      TrackerHitExtended * hitInTrack = hitsInTrack[iHit];
	      if (hitInTrack->getUsedInFit()) {
		TrackerHit * hit = hitInTrack->getTrackerHit();
		x_h[iHitInFit] = float(hit->getPosition()[0]);
		y_h[iHitInFit] = float(hit->getPosition()[1]);
		z_h[iHitInFit] = float(hit->getPosition()[2]);
		idet_h[iHitInFit] = hitInTrack->getDet();
		ityp_h[iHitInFit] = hitInTrack->getType();
		rR_h[iHitInFit] = hitInTrack->getResolutionRPhi();
		rZ_h[iHitInFit] = hitInTrack->getResolutionZ();
		iHitInFit++;
	      }
	    }
	    TrackerHit * remainHit = trkHitExt->getTrackerHit();
	    x_h[iHitInFit] = float(remainHit->getPosition()[0]);
	    y_h[iHitInFit] = float(remainHit->getPosition()[1]);
	    z_h[iHitInFit] = float(remainHit->getPosition()[2]);
	    idet_h[iHitInFit] = trkHitExt->getDet();
	    ityp_h[iHitInFit] = trkHitExt->getType();
	    rR_h[iHitInFit] = trkHitExt->getResolutionRPhi();
	    rZ_h[iHitInFit] = trkHitExt->getResolutionZ();	    
	    iHitInFit++;

	    int NPT = iHitInFit;
	    float chi2_D;
	    int ndf_D;
	    float chi2rphi,chi2z;
	    float par[5];
	    float epar[15];
	    float refPoint[3];

 	    int ierr = _trackFit.DoFitting(_useExtraPoint,_optFit,NPT,
 					   _bField,idet_h,ityp_h,
 					   _chi2PrefitCut,
 					   x_h,y_h,z_h,rR_h,rZ_h,
 					   par,epar,refPoint,chi2_D,ndf_D,
 					   chi2rphi,chi2z,lhits);	       

	    float chi2ndf = chi2_D/float(ndf_D);
	    if (ierr==0 && chi2ndf<_chi2FitCut) {
	      trkExt->addTrackerHitExtended( trkHitExt );
	      trkExt->setCovMatrix(epar);
	      trkExt->setOmega(par[0]);	  
	      trkExt->setTanLambda(par[1]);
	      trkExt->setPhi(par[2]);
	      trkExt->setD0(par[3]);
	      trkExt->setZ0(par[4]);
	      trkExt->setChi2(chi2_D);
	      trkExt->setNDF(ndf_D);
	      trkHitExt->setTrackExtended( trkExt );
	      trkHitExt->setUsedInFit( true );
	      flagTrack[trkExt] = false;
	      flagHit[trkHitExt] = false;
	      if (_optFit == 4) {
		GroupTracks * groupTracks = trkExt->getGroupTracks();
		if (groupTracks!=NULL) {
		  TrackExtendedVec group = groupTracks->getTrackExtendedVec();
		  int nTracks = int(group.size());
		  if (nTracks>0) {
		    TrackExtended * siTrack = group[0];
		    // preserving positiveness of matrix
		    float eD0 = siTrack->getCovMatrix()[0];
		    float eZ0 = siTrack->getCovMatrix()[9];
		    float sD0 = sqrt(eD0/epar[0]);
		    float sZ0 = sqrt(eZ0/epar[9]);
		    epar[0]  = eD0;
		    epar[9]  = eZ0;
		    epar[1]  = sD0*epar[1];
		    epar[3]  = sD0*epar[3];
		    epar[10] = sD0*epar[10];
		    epar[7]  = sZ0*epar[7];
		    epar[8]  = sZ0*epar[8];
		    epar[13] = sZ0*epar[13];
		    epar[6]  = sD0*sZ0*epar[6];
		    trkExt->setZ0(siTrack->getZ0());
		    trkExt->setD0(siTrack->getD0());
		    trkExt->setCovMatrix(epar);
		  }
		}
	      }
 	    }
	    delete[] x_h;
	    delete[] y_h;
	    delete[] z_h;
	    delete[] rR_h;
	    delete[] rZ_h;
	    delete[] idet_h;
	    delete[] ityp_h;
	    delete[] lhits;
	  }
	}
      }
	
      for (int iP=0;iP<nPairs;++iP) {
	TrackHitPair * trkHitPair = pairs[iP];
	delete trkHitPair;
      }
      
      pairs.clear();
      
    }
  

}

HelixClass * FullLDCTracking::GetExtrapolationHelix( TrackExtended * track) {
  
  TrackerHitExtendedVec hitVec = track->getTrackerHitExtendedVec();
  int nHits = int(hitVec.size());

  HelixClass * helix = new HelixClass();

  int nHitsFit = nHits;
  if (nHits>_nHitsExtrapolation)
    nHitsFit = _nHitsExtrapolation;

  float * ampl = new float[nHitsFit];
  float * xhit = new float[nHitsFit];
  float * yhit = new float[nHitsFit];
  float * zhit = new float[nHitsFit];
    

  if (nHits<=_nHitsExtrapolation) {
    for (int iH=0;iH<nHitsFit;++iH) {
      TrackerHit * trkHit = hitVec[iH]->getTrackerHit();
      ampl[iH] = 1.0;
      xhit[iH] = (trkHit->getPosition()[0]);
      yhit[iH] = (trkHit->getPosition()[1]);
      zhit[iH] = (trkHit->getPosition()[2]);
    }    

  }
  else {
    int * index = new int[nHits];
    float * zcoor = new float[nHits];        
    for (int iH=0;iH<nHits;++iH) {
      TrackerHit * trkHit = hitVec[iH]->getTrackerHit();
      zcoor[iH] = fabs(float(trkHit->getPosition()[2]));
    }
    int order = 1;
    GeneralSorting(index,zcoor,order,nHits);
    for (int iH=0;iH<nHitsFit;++iH) {
      int idx = index[iH];
      TrackerHit * trkHit = hitVec[idx]->getTrackerHit();
      ampl[iH] = 1.0;
      xhit[iH] = (trkHit->getPosition()[0]);
      yhit[iH] = (trkHit->getPosition()[1]);
      zhit[iH] = (trkHit->getPosition()[2]);
    }
    delete[] index;
    delete[] zcoor;
  }
 
  ClusterShapes * shapes = new ClusterShapes(nHitsFit,ampl,xhit,yhit,zhit);

  float parSh[5];
  float dparSh[5];
  float distmax = 0;
  float chi2Sh = 0;
  // do fitting
  shapes->FitHelix(500, 0, 1, parSh, dparSh, chi2Sh, distmax);  

  float x0Sh = parSh[0];
  float y0Sh = parSh[1];
  float r0Sh = parSh[2];
  float bzSh = parSh[3];
  float phi0Sh = parSh[4];
  float signPz = 1;
  float zBegin = zhit[0];
  float zEnd   = zhit[0];
  float zMax   = fabs(zBegin);
  float zMin   = fabs(zEnd);
  for (int iH=1;iH<nHitsFit;++iH) {
    float zCurrent = fabs(zhit[iH]);
    if (zCurrent>zMax) {
      zMax = zCurrent;
      zEnd = zhit[iH];
    }
    if (zCurrent<zMin) {
      zMin = zCurrent;
      zBegin = zhit[iH];
    }      
  }

  delete shapes;

  if (zEnd<zBegin)
    signPz = -1;
  helix->Initialize_BZ(x0Sh, y0Sh, r0Sh, 
		       bzSh, phi0Sh, _bField,signPz,
		       zBegin);
  
  delete[] ampl;
  delete[] xhit;
  delete[] yhit;
  delete[] zhit;
  
  return helix;

}


void FullLDCTracking::AssignTPCHitsToTracks(TrackerHitExtendedVec hitVec,
					    float dcut) {

    int nHits = int(hitVec.size());
    int nTrk = int(_trkImplVec.size());

    for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
      TrackExtended * foundTrack = _trkImplVec[iT];
      GroupTracks * group = foundTrack->getGroupTracks();
      TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
      int nTrkGrp = int(tracksInGroup.size());
      for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
	TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
	TrackerHitExtendedVec hitVecGrp = trkGrp->getTrackerHitExtendedVec();
	int nHits = int(hitVecGrp.size());
	float zMin = 1.0e+20;
	float zMax = -1.0e+20;
	float startPoint[3] = {0.,0.,0.};
	float endPoint[3]   = {0.,0.,0.};
	for (int iH=0;iH<nHits;++iH) {
	  TrackerHitExtended * trkHitExt = hitVecGrp[iH];
	  float pos[3] = {0.,0.,0.};
	  for (int iC=0;iC<3;++iC) 
	    pos[iC] = float(trkHitExt->getTrackerHit()->getPosition()[iC]);	  
	  if (pos[2]>zMax) {
	    zMax = pos[2];
	    for (int iC=0;iC<3;++iC)
	      endPoint[iC] = pos[iC];	      
	  }
	  if (pos[2]<zMin) {
	    zMin = pos[2];
	    for (int iC=0;iC<3;++iC)
	      startPoint[iC] = pos[iC];
	  }
	}
	trkGrp->setStart(startPoint);
	trkGrp->setEnd(endPoint);
      }
    }
    

    // replace previous version with faster loop ordering

    std::vector<float> minDistances(nHits, dcut);
    std::vector<TrackExtended*> tracksToAttach(nHits);
    std::vector< std::vector<float> > HitPositions(nHits);
    std::vector<int> HitSign(nHits);//Positive or Negative side
    for (int iH=0;iH<nHits;++iH) { // loop over leftover TPC hits
      tracksToAttach[iH]=NULL;
      //Get all TrackerHit positions, so we only have to get them once
      TrackerHit* temphit = hitVec[iH]->getTrackerHit();
      const double *temppos = temphit->getPosition();
      HitPositions[iH].push_back(float(temppos[0]));
      HitPositions[iH].push_back(float(temppos[1]));
      HitPositions[iH].push_back(float(temppos[2]));
      HitSign[iH]=signbit(temppos[2]);
    }    

    if(_debug>0)std::cout << " Starting loop " << nTrk << " tracks   and  " << nHits << " hits" << std::endl;

    for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
      TrackExtended * foundTrack = _trkImplVec[iT];
      int tanlambdaSign = signbit(foundTrack->getTanLambda());//we only care about positive or negative
      GroupTracks * group = foundTrack->getGroupTracks();
      TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
      int nTrkGrp = int(tracksInGroup.size());

      for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
	TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
	float tanLambda = trkGrp->getTanLambda();
	float omega = trkGrp->getOmega();
	float d0 = trkGrp->getD0();
	float z0 = trkGrp->getZ0();
	float phi0 = trkGrp->getPhi();
	float startPointZ = trkGrp->getStart()[2];
	float endPointZ   = trkGrp->getEnd()[2];
	HelixClass helix;
	helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
	float OnePFivehalfPeriodZ = 1.5*fabs(acos(-1.)*tanLambda/omega);

	for (int iH=0;iH<nHits;++iH) { // loop over leftover TPC hits

	  //check if the hit and the track or on the same side
	  //xor return 1, if hits are different
	  if ( tanlambdaSign^HitSign[iH] ) continue;

	  float DeltaStart = fabs(HitPositions[iH][2]-startPointZ);
	  float DeltaEnd = fabs(HitPositions[iH][2]-endPointZ);
	  bool consider = DeltaStart <= OnePFivehalfPeriodZ;
	  consider = consider || (DeltaEnd <= OnePFivehalfPeriodZ);
	  consider = consider || ( (HitPositions[iH][2]>=startPointZ) && (HitPositions[iH][2]<=endPointZ) );
	  
	  if(consider){
	    float distance = helix.getDistanceToPoint(HitPositions[iH], minDistances[iH]);
	    if (distance < minDistances[iH]) {
	      minDistances[iH] = distance;
	      tracksToAttach[iH] = foundTrack;
	    }
	  }
	} // loop over leftover TPC hits
      } //groups in tracks
    } // loop over all tracks

    for (int iH=0;iH<nHits;++iH) {
      TrackerHitExtended * trkHitExt = hitVec[iH];
      if (tracksToAttach[iH]!=NULL) {
	tracksToAttach[iH]->addTrackerHitExtended(trkHitExt);
	trkHitExt->setTrackExtended( tracksToAttach[iH] );
	trkHitExt->setUsedInFit( false );
      }
    }
	
    if(_debug>0)std::cout << " Fast loop done " << std::endl;
    

//     for (int iH=0;iH<nHits;iH++) { // loop over leftover TPC hits
// 	TrackerHitExtended * hitExt = hitVec[iH];
// 	float pos[3];
// 	for (int ip=0;ip<3;++ip) 
// 	    pos[ip] = float(hitExt->getTrackerHit()->getPosition()[ip]);
// 	float minDist = 1.0e+20;
// 	TrackExtended * trackToAttach = NULL;
// 	for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
// 	    TrackExtended * foundTrack = _trkImplVec[iT];
// 	    float tanLambdaFound = foundTrack->getTanLambda();
// 	    float product = tanLambdaFound*pos[2];
// 	    if (product>0) {
// 	      GroupTracks * group = foundTrack->getGroupTracks();
// 	      TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
// 	      int nTrkGrp = int(tracksInGroup.size());
// 	      for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
// 		TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
// 		float tanLambda = trkGrp->getTanLambda();
// 		float omega = trkGrp->getOmega();
// 		float d0 = trkGrp->getD0();
// 		float z0 = trkGrp->getZ0();
// 		float phi0 = trkGrp->getPhi();
// 		float dist[3];
// 		float startPoint[3];
// 		float endPoint[3];
// 		for (int iC=0;iC<3;++iC) {
// 		  startPoint[iC] = trkGrp->getStart()[iC];
// 		  endPoint[iC] = trkGrp->getEnd()[iC];
// 		}
// 		HelixClass helix;
// 		helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
// 		float halfPeriodZ = fabs(acos(-1.)*tanLambda/omega);
// 		helix.getDistanceToPoint(pos,dist);
// 		float DeltaStart = fabs(pos[2]-startPoint[2]);
// 		float DeltaEnd = fabs(pos[2]-endPoint[2]);
// 		bool consider = DeltaStart <= 1.5*halfPeriodZ;
// 		consider = consider || (DeltaEnd <= 1.5*halfPeriodZ);
// 		consider = consider || ( (pos[2]>=startPoint[2]) && (pos[2]<=endPoint[2]) );
// // 		float ZMin = DeltaStart;
// // 		if (DeltaEnd<ZMin)
// // 		  ZMin = DeltaEnd;
// 		if (dist[2]<dcut && consider && dist[2]<minDist) {
// 		  minDist = dist[2];
// 		  trackToAttach = foundTrack;
// 		}
// 	      }
// 	    }
// 	}

	//	if (trackToAttach!=NULL) {
	//   trackToAttach->addTrackerHitExtended(hitExt);
	//  hitExt->setTrackExtended( trackToAttach );
	//  hitExt->setUsedInFit( false );
	//  if(trackToAttach!=tracksToAttach[iH])std::cout << " Check Failed" << trackToAttach << "  " << tracksToAttach[iH] << std::endl;
	//
	//}
	//else {
	///	  std::cout << iH << " hit is not assigned : distance to closest track = " << minDist << std::endl;
	///}
    //}
    //std::cout << " Slow loop done " << std::endl;

}

void FullLDCTracking::AssignSiHitsToTracks(TrackerHitExtendedVec hitVec,
					   float dcut) {

    int nHits = int(hitVec.size());
    int nTrk = int(_allNonCombinedTPCTracks.size());

    std::map <TrackExtended*,bool> flagTrack;
    std::map <TrackerHitExtended*,bool> flagHit;
    TrackHitPairVec pairs;
    flagTrack.clear();
    flagHit.clear();
    pairs.clear();

    for (int iH=0;iH<nHits;++iH) {
	float pos[3];
	TrackerHitExtended * trkHitExt = hitVec[iH];
	TrackerHit * hit = trkHitExt->getTrackerHit();
	for (int ip=0;ip<3;++ip)
	    pos[ip] = float(hit->getPosition()[ip]);
	for (int iT=0;iT<nTrk;++iT) {
	    TrackExtended * trkExt = _allNonCombinedTPCTracks[iT];
	    float tanLambda = trkExt->getTanLambda();	    
	    float product = pos[2]*tanLambda;
	    if (product>0) {
		float d0 = trkExt->getD0();
		float z0 = trkExt->getZ0();
		float tanLambda = trkExt->getTanLambda();
		float phi0 = trkExt->getPhi();
		float omega = trkExt->getOmega();
		HelixClass helix;
		helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
		float distance = helix.getDistanceToPoint(pos,dcut);
		if (distance<dcut) {
		    TrackHitPair * trkHitPair = 
			new TrackHitPair(trkExt,trkHitExt,distance);
		    pairs.push_back(trkHitPair);
		    flagTrack[trkExt] = true;
		    flagHit[trkHitExt] = true;
		}
	    }
	}
    }

    int nPairs = int(pairs.size());
    if (nPairs>0) {
	SortingTrackHitPairs(pairs);
	for (int iP=0;iP<nPairs;++iP) {
	    TrackHitPair * trkHitPair = pairs[iP];
	    TrackExtended * trkExt = trkHitPair->getTrackExtended();
	    TrackerHitExtended * trkHitExt = 
		trkHitPair->getTrackerHitExtended();
	    if (flagTrack[trkExt] && flagHit[trkHitExt]) {		
		TrackerHitExtendedVec hitsInTrack = 
		    trkExt->getTrackerHitExtendedVec();
		int nTotH = int(hitsInTrack.size());
		int nHitsInFit = 0;
		for (int iTH=0;iTH<nTotH;++iTH) {
		    TrackerHitExtended * hitInTrack = hitsInTrack[iTH];
		    if (hitInTrack->getUsedInFit())
			nHitsInFit++;
		}
		float * x_h = new float[nHitsInFit+1];
		float * y_h = new float[nHitsInFit+1];
		float * z_h = new float[nHitsInFit+1];
		int * idet_h = new int[nHitsInFit+1];
		int * ityp_h = new int[nHitsInFit+1];
		int * lhits = new int[nHitsInFit+1];
		float * rR_h = new float[nHitsInFit+1];
		float * rZ_h = new float[nHitsInFit+1];
		int iHitInFit = 0;
		for (int iHit=0;iHit<nTotH;++iHit) {
		    TrackerHitExtended * hitInTrack = hitsInTrack[iHit];
		    if (hitInTrack->getUsedInFit()) {
			TrackerHit * hit = hitInTrack->getTrackerHit();
			x_h[iHitInFit] = float(hit->getPosition()[0]);
			y_h[iHitInFit] = float(hit->getPosition()[1]);
			z_h[iHitInFit] = float(hit->getPosition()[2]);
			idet_h[iHitInFit] = hitInTrack->getDet();
			ityp_h[iHitInFit] = hitInTrack->getType();
			rR_h[iHitInFit] = hitInTrack->getResolutionRPhi();
			rZ_h[iHitInFit] = hitInTrack->getResolutionZ();
			iHitInFit++;
		    }
		}
		TrackerHit * remainHit = trkHitExt->getTrackerHit();
		x_h[iHitInFit] = float(remainHit->getPosition()[0]);
		y_h[iHitInFit] = float(remainHit->getPosition()[1]);
		z_h[iHitInFit] = float(remainHit->getPosition()[2]);
		idet_h[iHitInFit] = trkHitExt->getDet();
		ityp_h[iHitInFit] = trkHitExt->getType();
		rR_h[iHitInFit] = trkHitExt->getResolutionRPhi();
		rZ_h[iHitInFit] = trkHitExt->getResolutionZ();		
		iHitInFit++;
		int NPT = iHitInFit;
		float chi2_D;
		int ndf_D;
		float chi2rphi,chi2z;
		float par[5];
		float epar[15];
		float refPoint[3];
		int ierr = _trackFit.DoFitting(_useExtraPoint,_optFit,NPT,
					       _bField,idet_h,ityp_h,
					       _chi2PrefitCut,
					       x_h,y_h,z_h,rR_h,rZ_h,
					       par,epar,refPoint,chi2_D,ndf_D,
					       chi2rphi,chi2z,lhits);	       
		float chi2ndf = chi2_D/float(ndf_D);
		if (ierr==0 && chi2ndf<_chi2FitCut) {
		  trkExt->addTrackerHitExtended( trkHitExt );
		  trkExt->setCovMatrix(epar);
		  trkExt->setOmega(par[0]);	  
		  trkExt->setTanLambda(par[1]);
		  trkExt->setPhi(par[2]);
		  trkExt->setD0(par[3]);
		  trkExt->setZ0(par[4]);
		  trkExt->setChi2(chi2_D);
		  trkExt->setNDF(ndf_D);
		  if (_optFit == 4) {
		    GroupTracks * group = trkExt->getGroupTracks();
		    if (group!=NULL) {
		      TrackExtendedVec trks = group->getTrackExtendedVec();
		      int nTrks = int(trks.size());
		      if (nTrks > 0) {
			TrackExtended * siTrack = trks[0];
			// preserving positiveness of cov. matrix
			float eD0 = siTrack->getCovMatrix()[0];
			float eZ0 = siTrack->getCovMatrix()[9];
			float sD0 = sqrt(eD0/epar[0]);
			float sZ0 = sqrt(eZ0/epar[9]);
			epar[0]  = eD0;
			epar[9]  = eZ0;
			epar[1]  = sD0*epar[1];
			epar[3]  = sD0*epar[3];
			epar[10] = sD0*epar[10];
			epar[7]  = sZ0*epar[7];
			epar[8]  = sZ0*epar[8];
			epar[13] = sZ0*epar[13];
			epar[6]  = sD0*sZ0*epar[6];
			trkExt->setZ0(siTrack->getZ0());
			trkExt->setD0(siTrack->getD0());
			trkExt->setCovMatrix(epar);			
		      }
		    }
		  }
		  trkHitExt->setTrackExtended( trkExt );
		  trkHitExt->setUsedInFit( true );
		  flagTrack[trkExt] = false;
		  flagHit[trkHitExt] = false;
		}
		delete[] x_h;
		delete[] y_h;
		delete[] z_h;
		delete[] rR_h;
		delete[] rZ_h;
		delete[] idet_h;
		delete[] ityp_h;
		delete[] lhits;
	    }
	}
	
	for (int iP=0;iP<nPairs;++iP) {
	    TrackHitPair * trkHitPair = pairs[iP];
	    delete trkHitPair;
	}

	pairs.clear();

    }
}

void FullLDCTracking::PrintOutMerging(TrackExtended * firstTrackExt, TrackExtended * secondTrackExt, int iopt) {
  // iopt = 1 false Si and TPC merging
  // iopt = 2 false Si and TPC forced merging
  // iopt = 3 false TPC segments merging
  // iopt = 4 false Comb Si and TPC merging
  // iopt = 5 false Comb TPC and TPC merging
  // iopt = 6 unmerged TPC and Si segments (soft merging)
  // iopt = 7 unmerged TPC and Si segments (forced merging)
  // iopt = 8 unmerged Comb and TPC
  // iopt = 9 unmerged TPC segments


  try {

    Track * firstTrack = firstTrackExt->getTrack();
    Track * secondTrack = secondTrackExt->getTrack();
    
    std::string firstColName = _TPCTrackMCPCollName;
    std::string secondColName = _TPCTrackMCPCollName;

    if (iopt==1) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==2) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==3) {
      secondColName = _TPCTrackMCPCollName;
    }
    else if (iopt==4) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==5) {
      secondColName = _TPCTrackMCPCollName;
    }
    else if (iopt==6) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==7) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==8) {
      secondColName = _SiTrackMCPCollName;
    }    
    else {
      secondColName = _TPCTrackMCPCollName;
    }

    
    LCCollection * firstCol = _evt->getCollection(firstColName.c_str());
    LCCollection * secondCol = _evt->getCollection(secondColName.c_str());

    LCRelationNavigator firstNav(firstCol);
    LCRelationNavigator secondNav(secondCol);
    LCObjectVec firstVec = firstNav.getRelatedToObjects(firstTrack);
    FloatVec firstWeights = firstNav.getRelatedToWeights(firstTrack);
    LCObject * firstMCP = NULL;
    float firstWght = 0;
    int nObj = firstVec.size();
    for (int iObj=0;iObj<nObj;++iObj) {
      if (firstWeights[iObj]>firstWght) {
	firstWght = firstWeights[iObj];
	firstMCP = firstVec[iObj];
      }
    }
    
    LCObjectVec secondVec = secondNav.getRelatedToObjects(secondTrack);
    FloatVec secondWeights = secondNav.getRelatedToWeights(secondTrack);
    LCObject * secondMCP = NULL;
    float secondWght = 0;
    nObj = secondVec.size();
    for (int iObj=0;iObj<nObj;++iObj) {
      if (secondWeights[iObj]>secondWght) {
	secondWght = secondWeights[iObj];
	secondMCP = secondVec[iObj];
      }
    }

  
    float d0First = firstTrackExt->getD0();
    float z0First = firstTrackExt->getZ0();
    float omegaFirst = firstTrackExt->getOmega();
    float tanLFirst = firstTrackExt->getTanLambda();
    float phi0First = firstTrackExt->getPhi();
    
    float d0Second = secondTrackExt->getD0();
    float z0Second = secondTrackExt->getZ0();
    float omegaSecond = secondTrackExt->getOmega();
    float tanLSecond = secondTrackExt->getTanLambda();
    float phi0Second = secondTrackExt->getPhi();	    

    HelixClass firstHelix;
    firstHelix.Initialize_Canonical(phi0First,d0First,z0First,omegaFirst,tanLFirst,_bField);
    float pxFirst = firstHelix.getMomentum()[0];
    float pyFirst = firstHelix.getMomentum()[1];
    float pzFirst = firstHelix.getMomentum()[2];	    
    
    HelixClass secondHelix;
    secondHelix.Initialize_Canonical(phi0Second,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
    float pxSecond = secondHelix.getMomentum()[0];
    float pySecond = secondHelix.getMomentum()[1];
    float pzSecond = secondHelix.getMomentum()[2];	    
    
    float dPx = pxFirst - pxSecond;
    float dPy = pyFirst - pySecond;
    float dPz = pzFirst - pzSecond;
    
    float dPplus = sqrt(dPx*dPx+dPy*dPy+dPz*dPz);
    
    dPx = pxFirst + pxSecond;
    dPy = pyFirst + pySecond;
    dPz = pzFirst + pzSecond;
    
    float dPminus = sqrt(dPx*dPx+dPy*dPy+dPz*dPz);
    
    float pFirst = sqrt(pxFirst*pxFirst+pyFirst*pyFirst+pzFirst*pzFirst);
    float pSecond = sqrt(pxSecond*pxSecond+pySecond*pySecond+pzSecond*pzSecond);

    if(pFirst < 1.5 && pSecond < 1.5 )return; 

    const float sigmaPOverPFirst  = sqrt(firstTrackExt->getCovMatrix()[5])/fabs(omegaFirst);
    const float sigmaPOverPSecond = sqrt(secondTrackExt->getCovMatrix()[5])/fabs(omegaSecond);
    //    const float deltaP = fabs(pFirst-pSecond);
    const float sigmaPFirst = pFirst*sigmaPOverPFirst;
    const float sigmaPSecond = pSecond*sigmaPOverPSecond;
    //const float sigmaDeltaP = sqrt(sigmaPFirst*sigmaPFirst+sigmaPSecond*sigmaPSecond);
    //    const float significance = deltaP/sigmaDeltaP;


    float den = pFirst;
    if (pSecond<pFirst)
      den = pSecond;
    
    dPplus  = dPplus/den;
    dPminus = dPminus/den; 
    
    if (firstMCP!=secondMCP && iopt < 6) {
      
      if (iopt==1) {
	std::cout << "Erroneous combining Si and TPC segments (iopt=1) --->" << std::endl;
      }
      else if (iopt==2) {
	std::cout << "Erroneous merging of Si and TPC segments (iopt=2) --->" << std::endl;
      }
      else if (iopt==3) {
	std::cout << "Erroneous merging of TPC segments (iopt=3) ---> " << std::endl; 
      }
      else if (iopt==4) {
	std::cout << "Erroneous merging of combSi segment with uncombTPC segment (iopt=4) ---> " << std::endl;
      }
      else {
	std::cout << "Erroneous merging of combTPC segment with uncombTPC segment (iopt=5) --->" << std::endl;
      }
      
      printf("%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
	     pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      printf("  %5.3f\n",firstWght);
      
      printf("%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
	     pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      printf("  %5.3f\n",secondWght);
      
      std::cout << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;


      TrackExtended * combinedTrack = TrialCombineTracks(firstTrackExt,secondTrackExt);
      
      if(combinedTrack != NULL){
	std::cout << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
	delete combinedTrack->getGroupTracks();
	delete combinedTrack;
      }else{
	std::cout << "Could not combine track " << std::endl;
      }
      std::cout << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;

      std::cout << std::endl;
    
    }
    else if (firstMCP==secondMCP && ( (iopt==8) || (iopt==9) ) ) {
      
      float deltaOmega = _dOmegaForMerging;
      float deltaAngle = _angleForMerging;
      if (iopt==8) {
	std::cout << "Unmerged TPC and Comb segments (iopt=8) --->" << std::endl;
      }
      else {
	std::cout << "Unmerged TPC segments (iopt=9) --->" << std::endl;
      }
      
      float qFirst = PIOVER2 - atan(tanLFirst);
      float qSecond = PIOVER2 - atan(tanLSecond);
      float dOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
      float angle = (cos(phi0First)*cos(phi0Second)+sin(phi0First)*sin(phi0Second))*
	sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
      angle = acos(angle);
      
      
      std::cout << " dOmegaCut = " << deltaOmega
		<< " AngleCut = " << deltaAngle
		<< " dOmega = " << dOmega
		<< " angle = " << angle << std::endl;
      

      printf("%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
	     pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      printf("  %5.3f\n",firstWght);
      
      printf("%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
	     pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      printf("  %5.3f\n",secondWght);

      std::cout << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;

      TrackExtended * combinedTrack = TrialCombineTracks(firstTrackExt,secondTrackExt);
      if(combinedTrack != NULL){
	std::cout << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
	delete combinedTrack->getGroupTracks();
	delete combinedTrack;
      }else{
	std::cout << "Could not combine track " << std::endl;
      }
      std::cout << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      
      std::cout << std::endl;
      
    }    
    else if (firstMCP==secondMCP && ( (iopt == 6) || (iopt == 7) ) ) {
      
      float deltaOmega = _dOmegaForMerging;
      float deltaAngle = _angleForMerging;
  
      if (iopt ==6) {
	std::cout << "Unmerged TPC and Si segments (iopt=6) --->" << std::endl;
      }
      else {
	std::cout << "Unmerged TPC and Si segments (iopt=7) --->" << std::endl;
	deltaOmega = _dOmegaForForcedMerging;
	deltaAngle = _angleForForcedMerging;
      }
      
      float qFirst = PIOVER2 - atan(tanLFirst);
      float qSecond = PIOVER2 - atan(tanLSecond);
      
      float dOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
      float angle = (cos(phi0First)*cos(phi0Second)+sin(phi0First)*sin(phi0Second))*
	sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
      angle = acos(angle);
      
      
      std::cout << " dOmegaCut = " << deltaOmega
		<< " AngleCut = " << deltaAngle
		<< " dOmega = " << dOmega
		<< " angle = " << angle << std::endl;
      printf("%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
	     pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      printf("  %5.3f\n",firstWght);
      
      printf("%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
	     pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      printf("  %5.3f\n",secondWght);
      
      std::cout << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;      
      TrackExtended * combinedTrack = TrialCombineTracks(firstTrackExt,secondTrackExt);
      
      if(combinedTrack != NULL){
	std::cout << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
	delete combinedTrack->getGroupTracks();
	delete combinedTrack;
      }else{
	std::cout << "Could not combine track " << std::endl;
      }
      std::cout << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      std::cout << std::endl;      
      
    }else if (firstMCP==secondMCP && iopt < 6) {
      return;
      if (iopt==1) {
	std::cout << "Correct combining Si and TPC segments (iopt=1) --->" << std::endl;
      }
      else if (iopt==2) {
	std::cout << "Correct merging of Si and TPC segments (iopt=2) --->" << std::endl;
      }
      else if (iopt==3) {
	std::cout << "Correct merging of TPC segments (iopt=3) ---> " << std::endl; 
      }
      else if (iopt==4) {
	std::cout << "Correct merging of combSi segment with uncombTPC segment (iopt=4) ---> " << std::endl;
      }
      else {
	std::cout << "Correct merging of combTPC segment with uncombTPC segment (iopt=5) --->" << std::endl;
      }
      
      printf("%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
	     pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      printf("  %5.3f\n",firstWght);
      
      printf("%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
	     pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      printf("  %5.3f\n",secondWght);
      
      std::cout << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;
      TrackExtended * combinedTrack = TrialCombineTracks(firstTrackExt,secondTrackExt);
      
      if(combinedTrack != NULL){
	std::cout << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
      }else{
	std::cout << "Could not combine track " << std::endl;
      }
      delete combinedTrack->getGroupTracks();
      delete combinedTrack;
      std::cout << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      
      std::cout << std::endl;
    
    }
  }

  catch(DataNotAvailableException &e){};

}    

void FullLDCTracking::GeneralSorting(int * index, float * val, int direct, int nVal) {
  /**
     Sorting of index vector in ascending (0) /descending (!=0) order of val
  */

  float valOne, valTwo, valTemp;
  int   indTemp;
  for (int i=0; i<nVal; ++i) {
    index[i] = i;
  }

  for (int i = 0 ; i < nVal-1; i++) {
    for (int j = 0; j < nVal-i-1; j++) {      
      valOne = val[j];
      valTwo = val[j+1];
      bool order = valOne > valTwo;
      if (direct>0) 
	order = valOne <= valTwo;
      if( order )
	{
	  valTemp = val[j];
	  val[j] = val[j+1];
	  val[j+1] = valTemp;
	  indTemp = index[j];
	  index[j] = index[j+1];
	  index[j+1] = indTemp;
	}
    }  
  }


}

int FullLDCTracking::SegmentRadialOverlap(TrackExtended* first, TrackExtended* second){

  const gear::TPCParameters& pTPC = Global::GEAR->getTPCParameters();
  const gear::PadRowLayout2D& pTPCpads = pTPC.getPadLayout();
  static const std::vector<double>tpcExtent = pTPCpads.getPlaneExtent();
  static const float tpcInnerR = tpcExtent[0]; 
  static const float tpcOuterR = tpcExtent[1];
  //static const float tpcZmax   = pTPC.getMaxDriftLength();
  static const float tpcMaxRow = pTPCpads.getNRows();


  float padHeight = (tpcOuterR-tpcInnerR)/(float)tpcMaxRow;


  int nTrkGrpFirst = 0;
  int nTrkGrpSecond = 0;
  TrackerHitVec hitvecFirst;
  TrackerHitVec hitvecSecond;
  GroupTracks * groupFirst = first->getGroupTracks();
  GroupTracks * groupSecond = second->getGroupTracks();
	
  if(groupFirst!=NULL){
	
    TrackExtendedVec tracksInGroupFirst = groupFirst->getTrackExtendedVec();
    nTrkGrpFirst = int(tracksInGroupFirst.size());
	
    for (int iTrkGrp=0;iTrkGrp<nTrkGrpFirst;++iTrkGrp) {
	
      TrackExtended * trkGrp = tracksInGroupFirst[iTrkGrp];
      TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();
	  
      for(unsigned int i =0; i<hitVec.size(); ++i){
	hitvecFirst.push_back(hitVec[i]->getTrackerHit());	  
      }
    }
  }
  if(groupSecond!=NULL){
    
    TrackExtendedVec tracksInGroupSecond = groupSecond->getTrackExtendedVec();
    nTrkGrpSecond = int(tracksInGroupSecond.size());
    
    for (int iTrkGrp=0;iTrkGrp<nTrkGrpSecond;++iTrkGrp) {
      TrackExtended * trkGrp = tracksInGroupSecond[iTrkGrp];
      TrackerHitExtendedVec hitVec = 
	trkGrp->getTrackerHitExtendedVec();
      
      for(unsigned int i=0;i<hitVec.size();++i){
	hitvecSecond.push_back(hitVec[i]->getTrackerHit());
      }
    }
  }
	

  int nhitsFirst = (int)hitvecFirst.size();
  int nhitsSecond = (int)hitvecSecond.size();
  int count = 0;
  for(int i =0;i<nhitsFirst;++i){
    float xi = (float)hitvecFirst[i]->getPosition()[0];
    float yi = (float)hitvecFirst[i]->getPosition()[1];
    float ri = sqrt(xi*xi+yi*yi);
    if(ri < tpcInnerR || ri > tpcOuterR)continue;
    for(int j =0;j<nhitsSecond;++j){
      float xj = (float)hitvecSecond[j]->getPosition()[0];
      float yj = (float)hitvecSecond[j]->getPosition()[1];
      float rj = sqrt(xj*xj+yj*yj);
      if(fabs(ri-rj)<padHeight/2.0)count++;
    }
  }  
  return count;
}

bool FullLDCTracking::VetoMerge(TrackExtended* firstTrackExt, TrackExtended* secondTrackExt){
  

  const float d0First = firstTrackExt->getD0();
  const float z0First = firstTrackExt->getZ0();
  const float omegaFirst = firstTrackExt->getOmega();
  const float tanLFirst = firstTrackExt->getTanLambda();
  const float phi0First = firstTrackExt->getPhi();
  
  const float d0Second = secondTrackExt->getD0();
  const float z0Second = secondTrackExt->getZ0();
  const float omegaSecond = secondTrackExt->getOmega();
  const float tanLSecond = secondTrackExt->getTanLambda();
  const float phi0Second = secondTrackExt->getPhi();	    
  
  HelixClass firstHelix;
  firstHelix.Initialize_Canonical(phi0First,d0First,z0First,omegaFirst,tanLFirst,_bField);
  const float pxFirst = firstHelix.getMomentum()[0];
  const float pyFirst = firstHelix.getMomentum()[1];
  const float pzFirst = firstHelix.getMomentum()[2];	    
    
  HelixClass secondHelix;
  secondHelix.Initialize_Canonical(phi0Second,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
  const float pxSecond = secondHelix.getMomentum()[0];
  const float pySecond = secondHelix.getMomentum()[1];
  const float pzSecond = secondHelix.getMomentum()[2];	    
  const float pFirst = sqrt(pxFirst*pxFirst+pyFirst*pyFirst+pzFirst*pzFirst);
  const float pSecond = sqrt(pxSecond*pxSecond+pySecond*pySecond+pzSecond*pzSecond);
  if(pFirst<2.5 || pSecond<2.5)return false;

  TrackExtended * combinedTrack = TrialCombineTracks(firstTrackExt,secondTrackExt);
  bool veto = false;
  if(combinedTrack!=NULL){
    if(combinedTrack->getNDF()+15<firstTrackExt->getNDF()+secondTrackExt->getNDF()+5)veto=true;
    delete combinedTrack->getGroupTracks();
    delete combinedTrack;
  }else{
    veto = true;
  }
  if(SegmentRadialOverlap(firstTrackExt,secondTrackExt)>10)veto=true;
  return veto;

}


void FullLDCTracking::check(LCEvent * evt) { }
void FullLDCTracking::end() {}


