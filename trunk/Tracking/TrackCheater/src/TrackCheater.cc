#include "TrackCheater.h"
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include "HelixClass.h"
#include "ClusterShapes.h"
#include <iostream>
#include <map>
#include <vector>
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/BField.h>

extern "C" {
  void tksetr_();
}

extern "C" {
  void tkinit_();
}



using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

typedef std::map <const MCParticle*,TrackImpl*> map_MCP_Track;
typedef std::map <const MCParticle*,HelixClass*> map_MCP_Helix;

TrackCheater aTrackCheater ;


TrackCheater::TrackCheater() : Processor("TrackCheater") {

  _description = "Creates true tracks..." ;
  
  registerOutputCollection( LCIO::TRACK,
			    "TrueTrackCollection",
			    "Collection of True Clusters",
			    _trueTracksCollection ,
			    std::string("TrueTracks"));
  
  std::vector<std::string> trackerHitCollections;  
  trackerHitCollections.push_back(std::string("VTXTrackerHits"));
  trackerHitCollections.push_back(std::string("FTDTrackerHits"));
  trackerHitCollections.push_back(std::string("SITTrackerHits"));
  trackerHitCollections.push_back(std::string("TPCTrackerHits"));
  
  registerInputCollections( LCIO::TRACKERHIT,
 			    "TrackerHitCollections",
			    "Tracker Hit Collection Names",
			    _trackerHitCollections ,
			    trackerHitCollections);
  
  registerOutputCollection( LCIO::LCRELATION, 
			    "MCTrueTrackRelCollectionName" , 
			    "Name of the TrueTrack MC Relation collection"  ,
			    _colNameMCTrueTracksRel ,
			    std::string("TrueTracksMCP") );
  
  registerProcessorParameter("ECut",
			     "Energy Cut",
			     _eCut ,
			     (float)0.1);

  registerProcessorParameter("HitToHelixDist",
			     "Cut on distance from hit to helix",
			     _hitToHelixCut,
			     (float)500.0);

  registerProcessorParameter("HitToHelixInFit",
			     "Cut on distance from hit to helix in fitting",
			     _hitToHelixInFit,
			     (float)20.0);

  registerProcessorParameter("FitTrueTrack",
			     "Flag to Fit True Track",
			     _fitTrueTrack,
			     (int)1);

  registerProcessorParameter("Chi2Cut",
			     "Cut On Fit Chi2",
			     _cutOnChi2,
			     (float)100.);    

  registerProcessorParameter("MinimalHits",
			     "Minimal Hits in Track Cluster",
			     _minimal_hits,
			     (int)3);

  registerProcessorParameter("UseExtraPoint",
			     "Use Extra Point in Fit?",
			     _useExtraPoint,
			     int(0));

  registerProcessorParameter("OptFit",
			     "Track Fit Option",
			     _optFit,
			     int(4));

  registerProcessorParameter("CutOnD0",
			     "Cut on d0 to accept track",
			     _cutOnD0,
			     float(500.));

  registerProcessorParameter("CutOnZ0",
			     "Cut on Z0 to accept track",
			     _cutOnZ0,
			     float(500.));

  registerProcessorParameter("CutOnTPCHits",
			     "Cut on TPC hits for tracks with no Si hits",
			     _cutOnTPCHits,
			     int(35));

  registerProcessorParameter("StoreHitsInFit",
			     "Store only hits used in fit?",
			     _storeHitsInFit,
			     int(0));

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

  registerProcessorParameter("Debug",
			     "Level of the printout info for the debuging purposes",
			     _debug,
			     int(1));
  

}

void TrackCheater::init() {
    _nRun = -1;
    _nEvt = 0;
    PI = acos(-1.0);
}


void TrackCheater::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 

void TrackCheater::processEvent( LCEvent * evt ) { 

  //    _bField = Global::parameters->getFloatVal("BField");
    _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;

    float parIpReso[3];
    parIpReso[0] = _aParIpReso;
    parIpReso[1] = _bParIpReso;
    parIpReso[2] = _sParIpReso;
    _trackFit.setParametersForIPErrors(parIpReso);
  
    map_MCP_Track _mcp_track;
    map_MCP_Helix _mcp_helix;

    int nonAssignedHits = 0;
    int assignedToNeutrals = 0;

    for (unsigned int i(0) ; i < _trackerHitCollections.size(); ++i) {    
	try {
	    LCCollection * col = 
	      evt->getCollection(_trackerHitCollections[i].c_str());
	    int nelem = col->getNumberOfElements();
	    for (int ielem(0); ielem < nelem; ++ielem) {
		TrackerHit * trkHit = 
		  dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
		LCObjectVec objVec =  trkHit->getRawHits();
		int nInVec = objVec.size();
		if (nInVec > 0 ) {
		    SimTrackerHit * simTrkHit = 
		      dynamic_cast<SimTrackerHit*>(objVec[0]);
		    const MCParticle * mcp = simTrkHit->getMCParticle();       
		    if ( mcp != NULL) {
		      bool _cuts = mcp->getEnergy() > _eCut ;
		      _cuts = _cuts && fabs(mcp->getCharge()) > 0.2;
		      if (fabs(mcp->getCharge())<0.2)
			assignedToNeutrals++;
		      if ( _cuts ) {
			if (_mcp_track[mcp] == NULL) {
			  TrackImpl * track = new TrackImpl();
			  float mom[3];
			  float ver[3];
			  HelixClass * helix = new HelixClass();
			  for (int icomp=0; icomp<3; ++icomp) {
			    mom[icomp]=(float)mcp->getMomentum()[icomp];
			    ver[icomp]=(float)mcp->getVertex()[icomp];
			  }
			  float charge = mcp->getCharge(); 
			  helix->Initialize_VP(ver,mom,charge,_bField);
			  _mcp_track[mcp] = track;
			  _mcp_helix[mcp] = helix;
			  float xPoint[3];
			  float Dist[3];
			  xPoint[0] = (float)trkHit->getPosition()[0];
			  xPoint[1] = (float)trkHit->getPosition()[1];
			  xPoint[2] = (float)trkHit->getPosition()[2];
			  helix->getDistanceToPoint(xPoint,Dist);
			  if (Dist[2] < _hitToHelixCut)
			    track->addHit( trkHit );
			}
			else {
			  TrackImpl * track  = _mcp_track[mcp];
			  HelixClass * helix = _mcp_helix[mcp];
			  float xPoint[3];
			  float Dist[3];
			  xPoint[0] = (float)trkHit->getPosition()[0];
			  xPoint[1] = (float)trkHit->getPosition()[1];
			  xPoint[2] = (float)trkHit->getPosition()[2];
			  helix->getDistanceToPoint(xPoint,Dist);
			  if (Dist[2]<_hitToHelixCut)
			    track->addHit( trkHit );
			}
		      }
		    }
		    else {
		      nonAssignedHits++;
		    }
		}
		
	    }
	}	
	catch(DataNotAvailableException &e){ 
	}
    }

    // debug    
    //     std::cout << std::endl;
    //     std::cout << std::endl;
    //     std::cout << "TrackCheater : " <<  nonAssignedHits 
    // 	      << " hits non pointing to MCParticle ; " 
    // 	      << assignedToNeutrals 
    // 	      << " hits assigned to neutrals " << std::endl;
    
    int nTotTracks = 0;
    float eTot = 0.;
    float pxTot = 0.;
    float pyTot = 0.;
    float pzTot = 0.;

    LCCollectionVec * trkcol = new LCCollectionVec(LCIO::TRACK);
    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit( LCIO::TRBIT_HITS ) ;
    trkcol->setFlag( trkFlag.getFlag()  ) ;


    LCCollectionVec * relationcol = new LCCollectionVec(LCIO::LCRELATION);
    map_MCP_Track::iterator pos;
    int itk = 0;
    int nlost = 0;
    int nincl = 0;
    int trkEnum = 0;
    if (_debug>=2) 
	std::cout << "Trk  # hits  MCParticle       Px       Py       Pz         D0         Z0    Si  TPC  chi2   ier Fit " << std::endl;
        //           "  0   191     0xd6ccf48     -1.861   -1.584   -1.071       0.00       0.00   0   0   0.00    0   0
    
    for (pos = _mcp_track.begin(); pos != _mcp_track.end(); ++pos) {
      const MCParticle * mcp = pos->first;
      TrackImpl * track = pos->second;
      TrackImpl * newTrack = new TrackImpl();
      // initialize variables that might not be filled in all cases
      newTrack->setNdf(0);
      newTrack->setChi2(0);
      newTrack->setdEdx(0);
      newTrack->setdEdxError(0);
      float zeromatrix[15]={0};
      newTrack->setCovMatrix(zeromatrix);
      HelixClass * helix = _mcp_helix[mcp];
      TrackerHitVec hitvec = track->getTrackerHits();
      int nhits = (int)hitvec.size();
      if (_debug>=2) {
	printf("%3i  %4i",trkEnum,nhits);
	std::cout << "     " << mcp;
      }
      trkEnum++;

      // debug	
      // 	std::cout << "Track " << itk << " Q = " << mcp->getCharge() 
      // 		  << " ; Px, Py, Pz = " 
      // 		  << "  " << mcp->getMomentum()[0]
      // 		  << "  " << mcp->getMomentum()[1]
      // 		  << "  " << mcp->getMomentum()[2] 
      // 		  << "  nhits = " << nhits << std::endl;

      bool storeTrack = true;
      float Pos[3];
      float Mom[3];
      Pos[0] = (float)mcp->getVertex()[0];
      Pos[1] = (float)mcp->getVertex()[1];
      Pos[2] = (float)mcp->getVertex()[2];
      Mom[0] = mcp->getMomentum()[0];
      Mom[1] = mcp->getMomentum()[1];
      Mom[2] = mcp->getMomentum()[2];

      int nHitsVTX = 0;
      int nHitsFTD = 0;
      int nHitsSIT = 0;
      int nHitsTPC = 0;
      int nHitsVTXFit = 0;
      int nHitsFTDFit = 0;
      int nHitsSITFit = 0;
      int nHitsTPCFit = 0;
      int nHitsSiFit = 0;

      int nHitsInFit=0;
      float chi2ndfTrk = 0;
      int status = 0;



      if (_fitTrueTrack==0) { // store True MC Information
	float d0Cur = helix->getD0();
	float z0Cur = helix->getZ0();
	if (d0Cur>_cutOnD0 || z0Cur>_cutOnZ0) {
	  storeTrack = false;
	}
	else {
	  if (nhits > _minimal_hits) {
	    nincl += nhits ;	  
	    newTrack->setReferencePoint(Pos);
	    newTrack->setD0(helix->getD0());
	    newTrack->setPhi(helix->getPhi0());
	    newTrack->setZ0(helix->getZ0());
	    newTrack->setOmega(helix->getOmega());
	    newTrack->setTanLambda(helix->getTanLambda());
	    int * lh = new int[nhits];
	    for (int ihit=0;ihit<nhits;++ihit) {
	      lh[ihit] = 1;
	      TrackerHit * trkHit = hitvec[ihit];
	      int det = trkHit->getType()/100;
	      if (det <= 4) {
		for (int lhit=0;lhit<ihit;++lhit) {
		  TrackerHit * trkHitS = hitvec[lhit];
		  if ((trkHitS->getType()==trkHit->getType()) && (lh[lhit]==1)) {
		    float xP[3];
		    float xPS[3];
		    for (int iC=0;iC<3;++iC) {
		      xP[iC] = float(trkHit->getPosition()[iC]);
		      xPS[iC] = float(trkHitS->getPosition()[iC]);
		    }
		    float Point[3];
		    float PointS[3];
		    if (det == 2) {
		      float time = helix->getPointInZ(xP[2],Pos,Point);
		      time = helix->getPointInZ(xPS[2],Pos,PointS);
		    }
		    else {
		      float RAD = sqrt(xP[0]*xP[0]+xP[1]*xP[1]);
		      float RADS = sqrt(xPS[0]*xPS[0]+xPS[1]*xPS[1]);
		      float time = helix->getPointOnCircle(RAD,Pos,Point);
		      time = helix->getPointOnCircle(RADS,Pos,PointS);
		    }
		    float DIST = 0;
		    float DISTS = 0;
		    for (int iC=0;iC<3;++iC) {
		      DIST += (Point[iC]-xP[iC])*(Point[iC]-xP[iC]);
		      DISTS += (PointS[iC]-xPS[iC])*(PointS[iC]-xPS[iC]);
		    }
		    if (DIST < DISTS) {
		      lh[lhit] = 0;
		    }
		    else {
		      lh[ihit] = 0;
		    }
		    break;
		  }
		}
	      }
	    }
	    double MinRSqHit = 99999;
	    for (int ihit=0;ihit<nhits;++ihit) {
	      if (lh[ihit] == 1) {
		TrackerHit * trkHit = hitvec[ihit];
		double RSqHit=trkHit->getPosition()[0]*trkHit->getPosition()[0]
		  +trkHit->getPosition()[1]*trkHit->getPosition()[1];
		if (RSqHit<MinRSqHit) MinRSqHit=RSqHit;
		newTrack->addHit( trkHit );
		int det = trkHit->getType()/100;
		if (det == 1) // VTX
		  nHitsVTX++;
		if (det == 2) // FTD
		  nHitsFTD++;
		if (det == 4) // SIT
		  nHitsSIT++;
		if (det == 5) // TPC
		  nHitsTPC++;
	      }
	    }
	    newTrack->setRadiusOfInnermostHit(sqrt(MinRSqHit));
	    newTrack->setReferencePoint(Pos);
	    newTrack->subdetectorHitNumbers().resize(8);
	    newTrack->subdetectorHitNumbers()[0] = nHitsVTX;
	    newTrack->subdetectorHitNumbers()[1] = nHitsFTD;
	    newTrack->subdetectorHitNumbers()[2] = nHitsSIT;
	    newTrack->subdetectorHitNumbers()[3] = nHitsTPC;
	    newTrack->subdetectorHitNumbers()[4] = nHitsVTX;
	    newTrack->subdetectorHitNumbers()[5] = nHitsFTD;
	    newTrack->subdetectorHitNumbers()[6] = nHitsSIT;
	    newTrack->subdetectorHitNumbers()[7] = nHitsTPC;
	    float PCA = sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]+Pos[2]*Pos[2]);
	    bool isRefPCA = false;
	    if (PCA<0.1)
	      isRefPCA = true;
	    newTrack->setIsReferencePointPCA(isRefPCA);
	    storeTrack = true;
	    int ntot = nHitsVTX + nHitsFTD + nHitsSIT + nHitsTPC;
	    if (ntot < _minimal_hits)
	      storeTrack = false;
	    delete[] lh; 
	  }
	  else {
	    storeTrack = false;
	  }
	}
      }
      else { // Track is supposed to be fitted
	float PeriodZ = fabs(2*PI*helix->getRadius()*helix->getTanLambda());
	TrackerHitVec hitsInFit;
	hitsInFit.clear();
	storeTrack = false;
	TrackerHitVec hitvec = track->getTrackerHits();	  
	float z1 = helix->getZ0() + PeriodZ*float(round((Pos[2]-helix->getZ0())/PeriodZ));
	float z2 = z1 + PI*helix->getRadius()*helix->getTanLambda();
	float zMin = z1;
	float zMax = z2;
	if (z1>z2) { 
	  zMin = z2;
	  zMax = z1;
	}

	int * lh = new int[nhits];

	// choose hits which are within the first track semi-loop
	// and close to the MC prediction for track tracjectory
	for (int ihit = 0; ihit < nhits; ++ihit) {
	  float xPoint[3];
	  TrackerHit * trkHit = hitvec[ihit];
	  xPoint[0] = (float)trkHit->getPosition()[0];
	  xPoint[1] = (float)trkHit->getPosition()[1];
	  xPoint[2] = (float)trkHit->getPosition()[2];
	  float distToH[3];
	  helix->getDistanceToPoint(xPoint,distToH);
	  if ((xPoint[2]>zMin+1.) && (xPoint[2]<zMax-1.) && distToH[2] < _hitToHelixInFit) {
	    lh[ihit] = 1;
	  }
	  else {
	    lh[ihit] = 0;
	  }
	}	  

	// check if there double hits in one Si layer
	for (int ihit=0;ihit<nhits;++ihit) {
	  TrackerHit * trkHit = hitvec[ihit];
	  int det = trkHit->getType()/100;
	  if (det <= 4) {
	    for (int lhit=0;lhit<ihit;++lhit) {
	      TrackerHit * trkHitS = hitvec[lhit];
	      if ((trkHitS->getType() == trkHit->getType()) && (lh[lhit] == 1) && lhit!=ihit) {
		float xP[3];
		float xPS[3];
		for (int iC=0;iC<3;++iC) {
		  xP[iC] = float(trkHit->getPosition()[iC]);
		  xPS[iC] = float(trkHitS->getPosition()[iC]);
		}
		float Point[3];
		float PointS[3];
		if (det == 2) {
		  float time = helix->getPointInZ(xP[2],Pos,Point);
		  time = helix->getPointInZ(xPS[2],Pos,PointS);
		}
		else {
		  float RAD = sqrt(xP[0]*xP[0]+xP[1]*xP[1]);
		  float RADS = sqrt(xPS[0]*xPS[0]+xPS[1]*xPS[1]);
		  float time = helix->getPointOnCircle(RAD,Pos,Point);
		  time = helix->getPointOnCircle(RADS,Pos,PointS);
		}
		float DIST = 0;
		float DISTS = 0;
		for (int iC=0;iC<3;++iC) {
		  DIST += (Point[iC]-xP[iC])*(Point[iC]-xP[iC]);
		  DISTS += (PointS[iC]-xPS[iC])*(PointS[iC]-xPS[iC]);
		}
		if (DIST < DISTS) {
		  lh[lhit] = 0;
		}
		else {
		  lh[ihit] = 0;
		}
	      }
	    }
	  }
	}

	for (int ihit=0;ihit<nhits;++ihit) {
	  if (lh[ihit]==1) {
	    TrackerHit * trkHit = hitvec[ihit];
	    nHitsInFit++;
	    hitsInFit.push_back(trkHit);	    
	  }
	}
	


	delete[] lh;

	if (nHitsInFit>=3) {
	  storeTrack = true;
	  float * xh = new float[nHitsInFit];
	  float * yh = new float[nHitsInFit];
	  float * zh = new float[nHitsInFit];
	  int   * idet = new int[nHitsInFit];
	  int   * ityp = new int[nHitsInFit];
	  float * rR = new float[nHitsInFit];
	  float * rZ = new float[nHitsInFit];
	  int * lhits = new int[nHitsInFit];
	  float parSi[5];
	  float eparSi[15];
	  int nSi = 0;
	  int foundSiSeg = 0;
	  if (_optFit >= 4) { // fitting first Si segment if found
	      for (int ihit = 0; ihit < nHitsInFit; ++ihit) {
		  float xPoint[3];
		  TrackerHit * trkHit = hitsInFit[ihit];
		  int det = trkHit->getType()/100;
		  if (det<5) {
		      xPoint[0] = (float)trkHit->getPosition()[0];
		      xPoint[1] = (float)trkHit->getPosition()[1];
		      xPoint[2] = (float)trkHit->getPosition()[2];	    
		      xh[nSi]=xPoint[0];
		      yh[nSi]=xPoint[1];
		      zh[nSi]=xPoint[2];
		      idet[nSi] = det;
		      if (det==2) { // FTD Planar Detector
			  ityp[nSi] = int(2);
			  rR[nSi] = sqrt(trkHit->getCovMatrix()[0]);
			  rZ[nSi] = 0.1;
		      }
		      else { // Cyllindrical Detector
			  ityp[nSi] = int(3);		  
			  rR[nSi] = sqrt(trkHit->getCovMatrix()[2]);
			  rZ[nSi] = sqrt(trkHit->getCovMatrix()[5]);		  
		      }
		      nSi++;
		  }
	      }
	      if (nSi>2) {
		  float chi2PrefitCutSi = 1.0e+10;
		  float chi2_DSi;
		  int ndf_DSi;
		  float chi2rphiSi,chi2zSi;
		  float refPointSi[3];
		  if (_trackFit.DoFitting(_useExtraPoint,_optFit,nSi,_bField,idet,ityp,chi2PrefitCutSi,
					  xh,yh,zh,rR,rZ,parSi,eparSi,refPointSi,chi2_DSi,ndf_DSi,
					  chi2rphiSi,chi2zSi,lhits)>=0)	      
		  foundSiSeg = 1;
	      }
	  }

	  for (int ihit = 0; ihit < nHitsInFit; ++ihit) {
	    float xPoint[3];
	    TrackerHit * trkHit = hitsInFit[ihit];
	    xPoint[0] = (float)trkHit->getPosition()[0];
	    xPoint[1] = (float)trkHit->getPosition()[1];
	    xPoint[2] = (float)trkHit->getPosition()[2];	    
	    xh[ihit]=xPoint[0];
	    yh[ihit]=xPoint[1];
	    zh[ihit]=xPoint[2];
	    int det = trkHit->getType()/100;
	    idet[ihit] = det;
	    if (det==2) { // FTD Planar Detector
	      ityp[ihit] = int(2);
	      rR[ihit] = sqrt(trkHit->getCovMatrix()[0]);
	      rZ[ihit] = 0.1;
	    }
	    else { // Cyllindrical Detector
	      ityp[ihit] = int(3);		  
	      rR[ihit] = sqrt(trkHit->getCovMatrix()[2]);
	      rZ[ihit] = sqrt(trkHit->getCovMatrix()[5]);		  
	    }
	  }
	  float chi2PrefitCut = 1.0e+10;
	  float par[5];
	  float epar[15];
	  float chi2_D;
	  int ndf_D;
	  float chi2rphi,chi2z;
	  int NPT = nHitsInFit;
	  float refPoint[3];
	  status = _trackFit.DoFitting(_useExtraPoint,_optFit,NPT,_bField,idet,ityp,chi2PrefitCut,
					   xh,yh,zh,rR,rZ,par,epar,refPoint,chi2_D,ndf_D,chi2rphi,chi2z,lhits);

	  chi2ndfTrk = fabs(chi2_D/float(ndf_D));

	  // hits in fit
	  for (int ihit=0;ihit<nHitsInFit;++ihit) {
	    if (lhits[ihit] > 0) {
	      TrackerHit * trkHit = hitsInFit[ihit];
	      int det = trkHit->getType()/100;
	      if (det == 1) // VTX
		nHitsVTXFit++;
	      if (det == 2) // FTD
		nHitsFTDFit++;
	      if (det == 4) // SIT
		nHitsSITFit++;
	      if (det == 5) // TPC
		nHitsTPCFit++;
	      if (_storeHitsInFit != 0)
		newTrack->addHit( trkHit );
	    }
	  }

	  bool rejectTrk = status < 0; // status of fit
	  rejectTrk = rejectTrk || (chi2ndfTrk > _cutOnChi2); // chi2 cut
	  rejectTrk = rejectTrk || (fabs(par[3]) > _cutOnD0); // cut on D0
	  rejectTrk = rejectTrk || (fabs(par[4]) > _cutOnZ0); // cut on Z0
	  
	  nHitsSiFit = nHitsVTXFit + nHitsFTDFit + nHitsSITFit;
	  
	  rejectTrk = rejectTrk || (nHitsSiFit==0 && (nHitsTPCFit<_cutOnTPCHits));

	  if (rejectTrk) {
	    storeTrack = false;
	  }	  
	  else {
	    if (foundSiSeg == 1) {
	      par[1] = parSi[1];
	      par[2] = parSi[2];
	      par[3] = parSi[3];
	      par[4] = parSi[4];
	      float scaling = sqrt(epar[5]/eparSi[5]);
	      float e2Omega = epar[5];
	      for (int iP=0;iP<15;++iP)
		epar[iP] = eparSi[iP];
	      epar[5] = e2Omega;
	      epar[3] = scaling*eparSi[3];
	      epar[4] = scaling*eparSi[4];
	      epar[8] = scaling*eparSi[8];
	      epar[12] = scaling*eparSi[12];	      
	    }
	    newTrack->setD0(par[3]);
	    newTrack->setZ0(par[4]);
	    newTrack->setPhi(par[2]);
	    newTrack->setTanLambda(par[1]);
	    newTrack->setOmega(par[0]);	  
	    newTrack->setCovMatrix(epar);
	    newTrack->setChi2(chi2_D);
	    newTrack->setNdf(ndf_D);
	    newTrack->setReferencePoint(refPoint);
	    newTrack->setIsReferencePointPCA(true);
	    int nHitsVTX = 0;
	    int nHitsFTD = 0;
	    int nHitsSIT = 0;
	    int nHitsTPC = 0;
	    double MinRSqHit = 99999;
	    for (int ihit=0;ihit<nhits;++ihit) {
	      TrackerHit * trkHit = hitvec[ihit];
	      double RSqHit=trkHit->getPosition()[0]*trkHit->getPosition()[0]
		+trkHit->getPosition()[1]*trkHit->getPosition()[1];
	      if (RSqHit<MinRSqHit) MinRSqHit=RSqHit;
	      int det = trkHit->getType()/100;
	      if (det == 1) // VTX
		nHitsVTX++;
	      if (det == 2) // FTD
		nHitsFTD++;
	      if (det == 4) // SIT
		nHitsSIT++;
	      if (det == 5) // TPC
		nHitsTPC++;	    
	      if (_storeHitsInFit == 0)
		newTrack->addHit( trkHit );	      
	    }
//************************************************************
//        Debug printout --->
// 	  std::cout << " VTX hits = " << nHitsVTX
// 		    << " FTD hits = " << nHitsFTD
// 		    << " SIT hits = " << nHitsSIT
// 		    << " TPC hits = " << nHitsTPC
// 		    << " N outliers = " << noutl << std::endl;
//*************************************************************
	    newTrack->setRadiusOfInnermostHit(sqrt(MinRSqHit));
	    newTrack->subdetectorHitNumbers().resize(8);
	    newTrack->subdetectorHitNumbers()[0] = nHitsVTXFit;
	    newTrack->subdetectorHitNumbers()[1] = nHitsFTDFit;
	    newTrack->subdetectorHitNumbers()[2] = nHitsSITFit;
	    newTrack->subdetectorHitNumbers()[3] = nHitsTPCFit;	  
	    newTrack->subdetectorHitNumbers()[4] = nHitsVTX;
	    newTrack->subdetectorHitNumbers()[5] = nHitsFTD;
	    newTrack->subdetectorHitNumbers()[6] = nHitsSIT;
	    newTrack->subdetectorHitNumbers()[7] = nHitsTPC;	  
	  }
	  delete[] xh;
	  delete[] yh;
	  delete[] zh;
	  delete[] idet;
	  delete[] ityp;
	  delete[] rR;
	  delete[] rZ;
	  delete[] lhits;
	}
	else {
	  storeTrack = false;
	}
      }
      if (storeTrack) {
	MCParticle * particle=0;
	try {
	  LCCollection * colpart = evt->getCollection("MCParticle");
	  int nelem = colpart->getNumberOfElements();
	  for (int ielem(0); ielem < nelem; ++ielem) {
	    particle = 
	      dynamic_cast<MCParticle*>(colpart->getElementAt(ielem));
	    if (particle == mcp) 
	      break;
	  }
	}
	catch(DataNotAvailableException &e){}
	float omegaT = newTrack->getOmega();
	float tanLambdaT = newTrack->getTanLambda();
	float phi0T = newTrack->getPhi();
	float d0T = newTrack->getD0();
	float z0T = newTrack->getZ0();
	HelixClass helixFinal;
	helixFinal.Initialize_Canonical(phi0T,d0T,z0T,omegaT,tanLambdaT,_bField);
	float trkPx = helixFinal.getMomentum()[0];
	float trkPy = helixFinal.getMomentum()[1];
	float trkPz = helixFinal.getMomentum()[2];
	float trkP = sqrt(trkPx*trkPx+trkPy*trkPy+trkPz*trkPz);
	eTot += trkP;
	pxTot += trkPx;
	pyTot += trkPy;
	pzTot += trkPz;	
	nTotTracks++;	
	trkcol->addElement( newTrack );
	LCRelationImpl * rel = new LCRelationImpl(newTrack,particle,(float)1.0);
	relationcol->addElement( rel );
	float d0Lost = helixFinal.getD0();
	float z0Lost = helixFinal.getZ0();
	if (_debug>=2) {
	    printf("   %8.3f %8.3f %8.3f %10.2f %10.2f %3i %3i %6.2f %4i %3i",
		   trkPx,trkPy,trkPz,d0Lost,z0Lost,nHitsSiFit,nHitsTPCFit,chi2ndfTrk,status,nHitsInFit);
	    std::cout << std::endl;
	}
      }
      else { // track is lost
	float d0Lost = helix->getD0();
	float z0Lost = helix->getZ0();
	if (_debug>=2) {
	    printf("   %8.3f %8.3f %8.3f %10.2f %10.2f %3i %3i %6.2f %4i %3i",
		   Mom[0],Mom[1],Mom[2],d0Lost,z0Lost,nHitsSiFit,nHitsTPCFit,chi2ndfTrk,status,nHitsInFit);
	    std::cout << " *** " << std::endl;
	}
	nlost += nhits;
	delete newTrack;
      }
      delete helix;
      delete track;
      itk++;
    
    // debug
    //std::cout << std::endl;
    }

    // debug
    //std::cout << std::endl;

    if (_debug>=1) {
	std::cout << "TrackCheater -> run " << _nRun
		  << " event " << _nEvt << std::endl; 
	std::cout << "Number of cheated tracks = " 
		  << nTotTracks << std::endl;
	std::cout << "Total 4-momentum of cheated tracks  : E = " << eTot
		  << " Px = " << pxTot
		  << " Py = " << pyTot
		  << " Pz = " << pzTot << std::endl;
	std::cout << std::endl;
    }
    evt->addCollection(trkcol,_trueTracksCollection);
    evt->addCollection(relationcol,_colNameMCTrueTracksRel);

  _nEvt++;

}

void TrackCheater::SortTrackerHitsByRadius(TrackerHitVec & trackerHitVec) {

  int sizeOfVector = int(trackerHitVec.size());
  TrackerHit *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++)
      {
	one = trackerHitVec[j];
	two = trackerHitVec[j+1];
	
	float xOne = float(one->getPosition()[0]);
	float xTwo = float(two->getPosition()[0]);
	float yOne = float(one->getPosition()[1]);
	float yTwo = float(two->getPosition()[1]);
	float rOne = xOne*xOne+yOne*yOne;
	float rTwo = xTwo*xTwo+yTwo*yTwo;
	if( rOne > rTwo )
	  {
	    Temp = trackerHitVec[j];
	    trackerHitVec[j] = trackerHitVec[j+1];
	    trackerHitVec[j+1] = Temp;
	  }
      }  

}


void TrackCheater::check( LCEvent * evt ) { }
  
void TrackCheater::end(){ } 
