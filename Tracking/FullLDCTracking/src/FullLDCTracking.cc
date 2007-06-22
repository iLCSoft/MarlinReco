#include "FullLDCTracking.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <iostream>
#include <math.h>
#include <map>
#include <marlin/Global.h>
#include "ClusterShapes.h"
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>

using namespace lcio ;
using namespace marlin ;

// extern "C" {
//   struct {
//     float consb;
//   } fkfild_;
// }

FullLDCTracking aFullLDCTracking ;

FullLDCTracking::FullLDCTracking() : Processor("FullLDCTracking") {  
  _description = "Performs full tracking in LDC detector" ;  

  registerInputCollection(LCIO::TRACKERHIT,
			  "FTDHitCollection",
			  "FTD Hit Collection Name",
			  _FTDTrackerHitCollection,
			  std::string("FTDTrackerHits"));  

  registerInputCollection(LCIO::TRACKERHIT,
			  "VXDHitCollection",
			  "VXD Hit Collection Name",
			  _VTXTrackerHitCollection,
			  std::string("VTXTrackerHits"));  

  registerInputCollection(LCIO::TRACKERHIT,
			  "SITHitCollection",
			  "SIT Hit Collection Name",
			  _SITTrackerHitCollection,
			  std::string("SITTrackerHits"));
  
  registerInputCollection(LCIO::TRACKERHIT,
			  "TPCHitCollection",
			  "TPC Hit Collection Name",
			  _TPCTrackerHitCollection,
			  std::string("TPCTrackerHits"));
  
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

//   registerProcessorParameter("ResolutionRPhi_VTX",
// 			     "R-Phi Resolution for VTX",
// 			     _resolutionRPhi_VTX,
// 			     float(0.004));

//   registerProcessorParameter("ResolutionZ_VTX",
// 			     "Z Resolution for VTX",
// 			     _resolutionZ_VTX,
// 			     float(0.004));

//   registerProcessorParameter("ResolutionRPhi_FTD",
// 			     "R-Phi Resolution for FTD",
// 			     _resolutionRPhi_FTD,
// 			     float(0.01));

//   registerProcessorParameter("ResolutionZ_FTD",
// 			     "Z Resolution for FTD",
// 			     _resolutionZ_FTD,
// 			     float(0.10));

//   registerProcessorParameter("ResolutionRPhi_SIT",
// 			     "R-Phi Resolution for SIT",
// 			     _resolutionRPhi_SIT,
// 			     float(0.01));

//   registerProcessorParameter("ResolutionZ_SIT",
// 			     "Z Resolution for SIT",
// 			     _resolutionZ_SIT,
// 			     float(0.01));

  registerProcessorParameter("AngleCutForMerging",
			     "Angle Cut For Merging",
			     _angleCutForMerging,
			     float(0.20));

  registerProcessorParameter("OmegaCutForMerging",
			     "Omega Cut For Merging",
			     _omegaCutForMerging,
			     float(0.25));

  registerProcessorParameter("Chi2FitCut",
			     "Cut on fit Chi2",
			     _chi2FitCut,
			     float(100.0));

  registerProcessorParameter("Chi2PrefitCut",
			     "Cut on fit Chi2",
			     _chi2PrefitCut,
			     float(1.0e+4));

  registerProcessorParameter("CreateMap",
			     "Create Track to MCP Relations",
			     _createMap,
			     int(1));

  registerProcessorParameter("Debug",
			     "Debug?",
			     _debug,
			     int(0));

  registerProcessorParameter("UseExtraPoint",
			     "Use Extra Point in Fit",
			     _useExtraPoint,
			     int(0));

  registerProcessorParameter("OptPrefit",
			     "Option of prefit ?",
			     _optFit,
			     int(4));

  registerProcessorParameter("ForceSiTPCMerging",
			     "Force merging of Si and TPC segments?",
			     _forceMerging,
			     int(1));

  registerProcessorParameter("AngleCutForForcedMerging",
			     "Angle Cut For Forced Merging",
			     _angleCutForForcedMerging,
			     float(0.05));

  registerProcessorParameter("OmegaCutForForcedMerging",
			     "Omega Cut For Forced Merging",
			     _omegaCutForForcedMerging,
			     float(0.25));

  registerProcessorParameter("RefitTPCTracks",
			     "Refit TPC Tracks ?",
			     _refitTPCTracks,
			     int(0));

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

  registerProcessorParameter("DeltaD0ToMergeTPCSegments",
			     "Cut on D0 difference for merging TPC segments",
			     _deltaD0ToMerge,
			     float(200));

  registerProcessorParameter("DeltaPToMergeTPCSegments",
			     "Cut on dP/P difference for merging TPC segments",
			     _deltaPtToMerge,
			     float(0.05));

  registerProcessorParameter("CutOnTPCHits",
			     "Cut on the number of TPC hits",
			     _cutOnTPCHits,
			     int(30));
  

}



void FullLDCTracking::init() { 

  printParameters();  
  _nRun = -1 ;
  _nEvt = 0 ;
  _bField = Global::parameters->getFloatVal("BField");
  PI = acos(-1.);
  PIOVER2 = 0.5*PI;
  TWOPI = 2*PI;

}

void FullLDCTracking::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
  std::cout << "New Run : run number = " << _nRun << std::endl;
  _bField = Global::parameters->getFloatVal("BField");
//   fkfild_.consb = 2.997924e-3*_bField;

} 

void FullLDCTracking::processEvent( LCEvent * evt ) { 

  std::cout << "Full LDC Tracking -> run " << _nRun 
	    << " event " << _nEvt << std::endl;
  prepareVectors( evt );
  //  std::cout << "prepareVectors done..." << std::endl;
  if (_storeRefittedTPCTracks>0 && _refitTPCTracks>0) {
    AddTrackColToEvt(evt,_allTPCTracks,
		     std::string("RefittedTPCTracks"),std::string("RefittedTPCTracksMCP"));
    // std::cout << "Collection of refitted TPC tracks is added to event..." << std::endl;
  }
  if (_storeRefittedSiTracks>0 && _refitSiTracks>0) {
     AddTrackColToEvt(evt,_allSiTracks,
		      std::string("RefittedSiTracks"),std::string("RefittedSiTracksMCP"));
     // std::cout << "Collection of refitted Si tracks is added to event..." << std::endl; 
  }
  MergeTPCandSiTracks();
  //  std::cout << "Merging done..." << std::endl;
  Sorting(_allCombinedTracks);
  //  std::cout << "Sorting done..." << std::endl;
  SelectCombinedTracks();
  //  std::cout << "Selection of combined tracks done..." << std::endl;
  AddNotCombinedTracks();
  //  std::cout << "Not combined tracks added..." << std::endl;
  AddTrackColToEvt(evt,_trkImplVec,
		   std::string("LDCTracks"),std::string("LDCTracksMCP"));
  //  std::cout << "Collections added to event..." << std::endl;

  CleanUp();

  //  std::cout << "Cleanup is done..." << std::endl;
  
  _nEvt++;
  //  getchar();

}

void FullLDCTracking::AddTrackColToEvt(LCEvent * evt, TrackExtendedVec & trkVec, 
				       std::string TrkColName, std::string RelColName) {
  
  LCCollectionVec * colTRK = new LCCollectionVec(LCIO::TRACK);

  std::cout << "Collection " << TrkColName
	    << " is added to event " << std::endl;

  LCCollectionVec * colRel = NULL;
  if (_createMap > 0)
    colRel = new LCCollectionVec(LCIO::LCRELATION);

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
    for (int iH=0;iH<nHits;++iH) {
      TrackerHitExtended * hitExt = hitVec[iH];
      TrackerHit * hit = hitExt->getTrackerHit();
      track->addHit(hit);
      int det = hit->getType()/100;
      if (det == 1)
	nHitsVTX++;
      if (det == 2)
	nHitsFTD++;
      if (det == 4)
	nHitsSIT++;
      if (det == 5)
	nHitsTPC++;
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
      TrackExtendedVec trkVec = group->getTrackExtendedVec();
      int nGrTRK = int(trkVec.size());
      for (int iGr=0;iGr<nGrTRK;++iGr) {
	TrackExtended * subTrack = trkVec[iGr];
	track->addTrack(subTrack->getTrack());
      }
    }

    track->subdetectorHitNumbers().resize(4);
    track->subdetectorHitNumbers()[0] = nHitsVTX;
    track->subdetectorHitNumbers()[1] = nHitsFTD;
    track->subdetectorHitNumbers()[2] = nHitsSIT;
    track->subdetectorHitNumbers()[3] = nHitsTPC;

    if ((nHitsTPC<_cutOnTPCHits) && ((nHitsVTX+nHitsFTD+nHitsSIT)<=0) ) {
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
  std::cout << "Number of " << TrkColName << " = " 
	    << nTotTracks << std::endl;
  std::cout << "Total 4-momentum of " << TrkColName << " : E = " << eTot
	    << " Px = " << pxTot
	    << " Py = " << pyTot
	    << " Pz = " << pzTot << std::endl;
  
  evt->addCollection(colTRK,TrkColName.c_str());
  if (_createMap)
    evt->addCollection(colRel,RelColName.c_str());

}


void FullLDCTracking::prepareVectors(LCEvent * event ) {


  _allTPCHits.clear();
  _allVTXHits.clear();
  _allFTDHits.clear();
  _allSITHits.clear();
  _allTPCTracks.clear();
  _allSiTracks.clear();
  _allCombinedTracks.clear();
  _trkImplVec.clear();

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
      hitExt->setResolutionRPhi(float(sqrt(hit->getCovMatrix()[2])));
      hitExt->setResolutionZ(float(sqrt(hit->getCovMatrix()[5])));
      hitExt->setType(int(3));
      hitExt->setDet(int(1));
      _allTPCHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
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
    std::cout << _FTDTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
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
    std::cout << _SITTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
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
    std::cout << _VTXTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }


  // Reading TPC Tracks
  try {
    LCCollection * col = event->getCollection(_TPCTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    if (_debug == 1) { 
      std::cout << "Number of TPC Tracks = " << nelem << std::endl;
      std::cout << "Trk    Omega      TanLambda    Phi         D0          Z0" << std::endl;
      //           "  0  0.00012154     0.000   1.555527     -0.006      0.002"
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
      reso[0] = 1e-3; // dPt/Pt = 1e-3*Pt
      reso[1] = 0.001; // angular resolution = 1mrad
      reso[2] = 1.0; // D0 reso = 1mm 
      reso[3] = 5.0; // Z0 reso = 5mm
      param[0] = tpcTrack->getOmega();
      param[1] = tpcTrack->getTanLambda();
      param[2] = tpcTrack->getPhi();
      param[3] = tpcTrack->getD0();
      param[4] = tpcTrack->getZ0();
      _trackFit.CrudeErrorEstimates(_bField,reso,param,cov);
      const FloatVec Cov = tpcTrack->getCovMatrix();
      int NC = int(Cov.size());
      for (int ic=0;ic<NC;ic++) {
	cov[ic] =  Cov[ic];
      }
      trackExt->setCovMatrix(cov);
      trackExt->setNDF(tpcTrack->getNdf());
      trackExt->setChi2(tpcTrack->getChi2());            
      char strg[200];
      sprintf(strg,"%3i %11.8f  %8.3f  %9.6f  %9.3f  %9.3f",iTrk,
	      tpcTrack->getOmega(),
	      tpcTrack->getTanLambda(),
	      tpcTrack->getPhi(),
	      tpcTrack->getD0(),
	      tpcTrack->getZ0());
      if (_debug==1) std::cout << strg << std::endl;
      for (int iHit=0;iHit<nHits;++iHit) {
	TrackerHit * hit = hitVec[iHit];
	TrackerHitExtended * hitExt = mapTrackerHits[hit];
	hitExt->setTrackExtended( trackExt );
	trackExt->addTrackerHitExtended( hitExt );	
      }      
      if (_refitTPCTracks > 0) {
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
	int ierr = _trackFit.DoFitting(_useExtraPoint,_optFit,NPT,_bField,idet_h,ityp_h,_chi2PrefitCut,
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
	delete[] lhits;
	delete[] x_h;
	delete[] y_h;
	delete[] z_h;
	delete[] rZ_h;
	delete[] rR_h;
      }      
      _allTPCTracks.push_back( trackExt );      
    }    
  }
  catch ( DataNotAvailableException &e) {
    std::cout << _TPCTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading Si Tracks
  try {
    LCCollection * col = event->getCollection(_SiTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    if (_debug == 1) { 
      std::cout << "Number of Si Tracks = " << nelem << std::endl;
      std::cout << "Trk    Omega      TanLambda    Phi         D0          Z0" << std::endl;
      //           "  0  0.00012154     0.000   1.555527     -0.006      0.002"
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
      reso[0] = 1e-2; // dPt/Pt = 1e-2*Pt
      reso[1] = 0.001; // angular resolution = 1mrad
      reso[2] = 1.0; // D0 reso = 1mm 
      reso[3] = 5.0; // Z0 reso = 5mm
      param[0] = siTrack->getOmega();
      param[1] = siTrack->getTanLambda();
      param[2] = siTrack->getPhi();
      param[3] = siTrack->getD0();
      param[4] = siTrack->getZ0();
      _trackFit.CrudeErrorEstimates(_bField,reso,param,cov);
      const FloatVec Cov = siTrack->getCovMatrix();
      int NC = int(Cov.size());
      for (int ic=0;ic<NC;ic++) {
	cov[ic] =  Cov[ic];
      }
      trackExt->setCovMatrix(cov);
      trackExt->setNDF(siTrack->getNdf());
      trackExt->setChi2(siTrack->getChi2());      
      char strg[200];
      sprintf(strg,"%3i %11.8f  %8.3f  %9.6f  %9.3f  %9.3f",iTrk,
	      siTrack->getOmega(),
	      siTrack->getTanLambda(),
	      siTrack->getPhi(),
	      siTrack->getD0(),
	      siTrack->getZ0());
      if (_debug==1) std::cout << strg << std::endl;
      for (int iHit=0;iHit<nHits;++iHit) {
	TrackerHit * hit = hitVec[iHit];
	TrackerHitExtended * hitExt = mapTrackerHits[hit];
	hitExt->setTrackExtended( trackExt );
	trackExt->addTrackerHitExtended( hitExt );	
      }
      if (_refitSiTracks > 0) {
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
	int ierr = _trackFit.DoFitting(_useExtraPoint,_optFit,NPT,_bField,idet_h,ityp_h,_chi2PrefitCut,
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
      _allSiTracks.push_back( trackExt );      
    }    
  }
  catch ( DataNotAvailableException &e) {
    std::cout << _SiTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }



};

void FullLDCTracking::CleanUp(){

  int nSITHits = int(_allSITHits.size());
  for (int i=0;i<nSITHits;++i) {
    TrackerHitExtended * hitExt = _allSITHits[i];
    delete hitExt;
  }
  _allSITHits.clear();

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

};

void FullLDCTracking::MergeTPCandSiTracks() {

  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks  = int(_allSiTracks.size());

  for (int iTPC=0;iTPC<nTPCTracks;++iTPC) {
    TrackExtended * tpcTrackExt = _allTPCTracks[iTPC];
    Track * tpcTrack = tpcTrackExt->getTrack();
    float phiTPC = tpcTrack->getPhi();
    float tanLambdaTPC = tpcTrack->getTanLambda();
    float qTPC = PIOVER2 - atan(tanLambdaTPC);
    float oTPC = tpcTrack->getOmega();
    for (int iSi=0;iSi<nSiTracks;++iSi) {
      TrackExtended * siTrackExt = _allSiTracks[iSi];
      Track * siTrack = siTrackExt->getTrack();
      float phiSi = siTrack->getPhi();
      float tanLambdaSi = siTrack->getTanLambda();
      float qSi = PIOVER2 - atan(tanLambdaSi);
      float oSi = siTrack->getOmega();
      float angle = (cos(phiTPC)*cos(phiSi)+sin(phiTPC)*sin(phiSi))*sin(qTPC)*sin(qSi)+cos(qTPC)*cos(qSi);
      angle = acos(angle);
      float deltaOmega = fabs((oSi-oTPC)/oTPC);
      if (angle < _angleCutForMerging && deltaOmega < _omegaCutForMerging) {
	TrackExtended * combinedTrack = CombineTracks(tpcTrackExt,siTrackExt);
	if (combinedTrack != NULL) {
	  _allCombinedTracks.push_back( combinedTrack );
	}
      }
    }
  }


};

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
  par[0] = siTrack->getOmega();
  par[1] = siTrack->getTanLambda();
  par[2] = siTrack->getPhi();
  par[3] = siTrack->getD0();
  par[4] = siTrack->getZ0();
  float refPoint[3];
  int ierr = _trackFit.DoFitting(_useExtraPoint,_optFit,NPT,_bField,idet_h,ityp_h,_chi2PrefitCut,
				 xh,yh,zh,rphi_reso,z_reso,
				 par,epar,refPoint,chi2_D,ndf_D,chi2rphi,chi2z,lhits);
  
  float chiQ = chi2_D/float(ndf_D);
  if (ierr >= 0 && chiQ > 0.001) {
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
    group->addTrackExtended(tpcTrack);
    group->addTrackExtended(siTrack);
    OutputTrack->setGroupTracks(group);
    OutputTrack->setOmega(omega);
    OutputTrack->setTanLambda(tanlambda);
    OutputTrack->setPhi(phi0);
    OutputTrack->setZ0(z0);
    OutputTrack->setD0(d0);
    OutputTrack->setChi2(chi2_D);
    OutputTrack->setNDF(ndf_D);
    OutputTrack->setCovMatrix(epar);
    for (int i=0;i<nTPCHits;++i) {
      TrackerHitExtended * hitExt = tpcHitVec[i];
      OutputTrack->addTrackerHitExtended(hitExt);
    }
    for (int i=0;i<nSiHits;++i) {
      TrackerHitExtended * hitExt = siHitVec[i];
      OutputTrack->addTrackerHitExtended(hitExt);
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
	_trkImplVec.push_back(trkExt);	
      }
    }
  }


}

void FullLDCTracking::AddNotCombinedTracks() {  

  int nTPCTrk = int(_allTPCTracks.size());
  int nSiTrk = int(_allSiTracks.size());

  if (_forceMerging==1) {
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExtTPC = _allTPCTracks[i];
      GroupTracks * groupTPC = trkExtTPC->getGroupTracks();
      if (groupTPC == NULL) {
	for (int j=0;j<nSiTrk;++j) {
	  TrackExtended * trkExtSi = _allSiTracks[j];
	  GroupTracks * groupSi = trkExtSi->getGroupTracks();
	  if (groupSi == NULL) {
	    float omegaTPC = trkExtTPC->getOmega();
	    float omegaSi = trkExtSi->getOmega();
	    float dOmega = fabs((omegaTPC-omegaSi)/omegaTPC);
	    float phiSi = trkExtSi->getPhi();
	    float tanLambdaSi = trkExtSi->getTanLambda();
	    float qSi = PIOVER2 - atan(tanLambdaSi);
	    float phiTPC = trkExtTPC->getPhi();
	    float tanLambdaTPC = trkExtTPC->getTanLambda();
	    float qTPC = PIOVER2 - atan(tanLambdaTPC);
	    float angle = (cos(phiTPC)*cos(phiSi)+sin(phiTPC)*sin(phiSi))*sin(qTPC)*sin(qSi)+cos(qTPC)*cos(qSi);
	    angle = acos(angle);
	    if (dOmega<_omegaCutForForcedMerging) {
	      if (angle<_angleCutForForcedMerging) {
		TrackExtended * OutputTrack = new TrackExtended();
		GroupTracks * group = new GroupTracks();
		group->addTrackExtended(trkExtTPC);
		group->addTrackExtended(trkExtSi);
		OutputTrack->setGroupTracks(group);
		trkExtSi->setGroupTracks(group);
		trkExtTPC->setGroupTracks(group);	      
		OutputTrack->setOmega(trkExtSi->getOmega());
		OutputTrack->setTanLambda(trkExtSi->getTanLambda());
		OutputTrack->setPhi(trkExtSi->getPhi());
		OutputTrack->setZ0(trkExtSi->getZ0());
		OutputTrack->setD0(trkExtSi->getD0());
		OutputTrack->setChi2(trkExtTPC->getChi2());
		OutputTrack->setNDF(trkExtTPC->getNDF());
		OutputTrack->setCovMatrix(trkExtSi->getCovMatrix());
		TrackerHitExtendedVec tpcHitVec = trkExtTPC->getTrackerHitExtendedVec();
		TrackerHitExtendedVec siHitVec = trkExtSi->getTrackerHitExtendedVec();	      
		int nTPCHits = int( tpcHitVec.size());
		int nSiHits = int( siHitVec.size());	      
		for (int iH=0;iH<nTPCHits;++iH) {
		  TrackerHitExtended * hitExt = tpcHitVec[iH];
		  OutputTrack->addTrackerHitExtended(hitExt);
		}
		for (int iH=0;iH<nSiHits;++iH) {
		  TrackerHitExtended * hitExt = siHitVec[iH];
		  OutputTrack->addTrackerHitExtended(hitExt);
		}	  
		_trkImplVec.push_back(OutputTrack);
		break;
	      }	    
	    }
	  }
	}
      }
    }  
  }

  if (_mergeTPCSegments) {
    std::vector<GroupTracks*> TPCSegments;
    TPCSegments.clear();
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExt = _allTPCTracks[i];
      GroupTracks * group = trkExt->getGroupTracks();
      if (group == NULL) {
	int nGroups = int(TPCSegments.size());
	float dPtMin = 1.0e+10;
	GroupTracks * groupToAttach = NULL;
	for (int iG=0;iG<nGroups;++iG) {
	  GroupTracks * segments = TPCSegments[iG];
	  TrackExtendedVec segVec = segments->getTrackExtendedVec();
	  int nTrk = int(segVec.size());
	  for (int iTrk=0;iTrk<nTrk;++iTrk) {
	    TrackExtended * trkInGroup = segVec[iTrk];
	    float dPt = CompareTrk(trkExt,trkInGroup);
	    if (dPt < dPtMin) {
	      dPtMin = dPt;
	      groupToAttach = segments;
	      break;
	    }
	  }
	}
	if (dPtMin < _deltaPtToMerge && groupToAttach != NULL) {
	  groupToAttach->addTrackExtended(trkExt);
	  trkExt->setGroupTracks(groupToAttach);
	}
	else {
	  GroupTracks * newSegment = new GroupTracks(trkExt);
	  trkExt->setGroupTracks(newSegment);
	  TPCSegments.push_back(newSegment);
	}
      }
    }
    int nCombTrk = int(_trkImplVec.size());
    int nSegments = int(TPCSegments.size());
    for (int iS=0;iS<nSegments;++iS) {
      GroupTracks * segments = TPCSegments[iS];
      TrackExtendedVec segVec = segments->getTrackExtendedVec();
      int nTrk = int(segVec.size());
      int found = 0;
      for (int iCTrk=0;iCTrk<nCombTrk;++iCTrk) {
	TrackExtended * combTrk = _trkImplVec[iCTrk];
	for (int iTrk=0;iTrk<nTrk;++iTrk) {
	  TrackExtended * trk = segVec[iTrk];
	  float dPt = CompareTrk(trk,combTrk);
	  if (dPt<_deltaPtToMerge) {
	    found = 1;
	    break;
	  }
	}
	if (found == 1)
	  break;
      }
      if (found == 0) {
	if (nTrk==1) {
	  _trkImplVec.push_back(segVec[0]);
	  segVec[0]->setGroupTracks(NULL);
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
	    _trkImplVec.push_back(chosenTrack);
	    for (int iTrk=0;iTrk<nTrk;++iTrk) {
		TrackExtended * trk = segVec[iTrk];
		if (trk!=chosenTrack) {
		    TrackerHitExtendedVec hitVec = trk->getTrackerHitExtendedVec();
		    int nHits = int(hitVec.size());
		    for (int iH=0;iH<nHits;++iH) {
			TrackerHitExtended * trkHitExt = hitVec[iH];
			chosenTrack->addTrackerHitExtended(trkHitExt);
		    }
		}
	    }
	    chosenTrack->setGroupTracks(NULL);
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
  else {
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExt = _allTPCTracks[i];
      Track * track = trkExt->getTrack();
      GroupTracks * group = trkExt->getGroupTracks();
      if (group == NULL) {
	TrackerHitVec hitVec = track->getTrackerHits();       
	_trkImplVec.push_back(trkExt);
      }
    }    
  }

  for (int i=0;i<nSiTrk;++i) {
    TrackExtended * trkExt = _allSiTracks[i];
    Track * track = trkExt->getTrack();
    GroupTracks * group = trkExt->getGroupTracks();
    if (group == NULL) {
      TrackerHitVec hitVec = track->getTrackerHits();
      _trkImplVec.push_back(trkExt);
    }
  }

}


float FullLDCTracking::CompareTrk(TrackExtended * first, TrackExtended * second) {

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

  if (fabs(d0First-d0Second)<_deltaD0ToMerge || fabs(d0First+d0Second)<_deltaD0ToMerge) {

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
    momMinus = sqrt(momMinus);
    momPlus = sqrt(momPlus);
    float nom = momMinus;
    if (momPlus<nom)
      nom = momPlus;
    float den = momFirst;
    if (momSecond<momFirst)
      den = momSecond;

    result = nom/den;
//     if (momMinus<0.1 || momPlus <0.1) {
//       std::cout << pFirst[0] << " "
// 		<< pFirst[1] << " "
// 		<< pFirst[2] << std::endl;
//       std::cout << pSecond[0] << " "
// 		<< pSecond[1] << " "
// 		<< pSecond[2] << std::endl;
//       std::cout << result << std::endl;	
//     }

  }

  return result;

}

void FullLDCTracking::check(LCEvent * evt) { };
void FullLDCTracking::end() {};
    
