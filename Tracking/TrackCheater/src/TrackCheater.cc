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
			     (float)10.0);

  registerProcessorParameter("FitTrueTrack",
			     "Flag to Fit True Track",
			     _fitTrueTrack,
			     (int)1);


  registerProcessorParameter("Chi2Cut",
			     "Cut On Fit Chi2",
			     _chi2Cut,
			     (float)100.);    

  registerProcessorParameter("MinimalHits",
			     "Minimal Hits in Track Cluster",
			     _minimal_hits,
			     (int)3);

  registerProcessorParameter("UseExtraPoint",
			     "Use Extra Point in Fit?",
			     _useExtraPoint,
			     int(0));

  registerProcessorParameter("OptPrefitFit",
			     "Prefit Option",
			     _optFit,
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
  _bField = Global::parameters->getFloatVal("BField");
} 

void TrackCheater::processEvent( LCEvent * evt ) { 


  // set flag to store more information in the output file
  LCFlagImpl flag;
  flag.setBit(LCIO::TRBIT_HITS);


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
			  if (fabs(Dist[2])<_hitToHelixCut)
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
    

    LCCollectionVec * trkcol = new LCCollectionVec(LCIO::TRACK);
    trkcol->setFlag(flag.getFlag());

    
    LCCollectionVec * relationcol = new LCCollectionVec(LCIO::LCRELATION);
    map_MCP_Track::iterator pos;
    int itk = 0;
    int nlost = 0;
    int nincl = 0;
    for (pos = _mcp_track.begin(); pos != _mcp_track.end(); ++pos) {
      const MCParticle * mcp = pos->first;
      TrackImpl * track = pos->second;
      TrackImpl * newTrack = new TrackImpl();
      HelixClass * helix = _mcp_helix[mcp];
      TrackerHitVec hitvec = track->getTrackerHits();
      int nhits = (int)hitvec.size();

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
      if (_fitTrueTrack==0) { // store True MC Information
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
	  int nHitsVTX = 0;
	  int nHitsFTD = 0;
	  int nHitsSIT = 0;
	  int nHitsTPC = 0;
	  for (int ihit=0;ihit<nhits;++ihit) {
	    if (lh[ihit] == 1) {
	      TrackerHit * trkHit = hitvec[ihit]; 
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
	  newTrack->subdetectorHitNumbers().resize(4);
	  newTrack->subdetectorHitNumbers()[0] = nHitsVTX;
	  newTrack->subdetectorHitNumbers()[1] = nHitsFTD;
	  newTrack->subdetectorHitNumbers()[2] = nHitsSIT;
	  newTrack->subdetectorHitNumbers()[3] = nHitsTPC;	  
	  newTrack->setIsReferencePointPCA(false);
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
      else { // Track is supposed to be fitted
	int nHitsInFit=0;
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
	for (int ihit=0;ihit<nhits;++ihit) {
	  lh[ihit] = 1;
	  TrackerHit * trkHit = hitvec[ihit];
	  int det = trkHit->getType()/100;
	  if (det <= 4) {
	    for (int lhit=0;lhit<ihit;++lhit) {
	      TrackerHit * trkHitS = hitvec[lhit];
	      if ((trkHitS->getType() ==trkHit->getType()) && (lh[lhit] == 1)) {
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

	for (int ihit = 0; ihit < nhits; ++ihit) {
	  if (lh[ihit] == 1) {
	    float xPoint[3];
	    TrackerHit * trkHit = hitvec[ihit];
	    xPoint[0] = (float)trkHit->getPosition()[0];
	    xPoint[1] = (float)trkHit->getPosition()[1];
	    xPoint[2] = (float)trkHit->getPosition()[2];
	    if ((xPoint[2]>zMin+1) && (xPoint[2]<zMax-1)) {
	      nHitsInFit++;
	      hitsInFit.push_back(trkHit);
	    }
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
	  _trackFit.DoFitting(_useExtraPoint,_optFit,NPT,_bField,idet,ityp,chi2PrefitCut,
			      xh,yh,zh,rR,rZ,par,epar,refPoint,chi2_D,ndf_D,chi2rphi,chi2z,lhits);
	  newTrack->setD0(par[3]);
	  newTrack->setZ0(par[4]);
	  newTrack->setPhi(par[2]);
	  newTrack->setTanLambda(par[1]);
	  newTrack->setOmega(par[0]);	  
	  newTrack->setCovMatrix(epar);
	  newTrack->setChi2(chi2_D);
	  newTrack->setNdf(ndf_D);
	  newTrack->setReferencePoint(refPoint);
	  newTrack->setIsReferencePointPCA(false);
	  if (_useExtraPoint == 1)
	    newTrack->setIsReferencePointPCA(true);
	  int nHitsVTX = 0;
	  int nHitsFTD = 0;
	  int nHitsSIT = 0;
	  int nHitsTPC = 0;
	  for (int ihit=0;ihit<nHitsInFit;++ihit) {
	    if (lhits[ihit] > 0) {
	      TrackerHit * trkHit = hitsInFit[ihit];
	      int det = trkHit->getType()/100;
	      if (det == 1) // VTX
		nHitsVTX++;
	      if (det == 2) // FTD
		nHitsFTD++;
	      if (det == 4) // SIT
		nHitsSIT++;
	      if (det == 5) // TPC
		nHitsTPC++;	    
	      newTrack->addHit( trkHit );
	    }
	  }
	  newTrack->subdetectorHitNumbers().resize(4);
	  newTrack->subdetectorHitNumbers()[0] = nHitsVTX;
	  newTrack->subdetectorHitNumbers()[1] = nHitsFTD;
	  newTrack->subdetectorHitNumbers()[2] = nHitsSIT;
	  newTrack->subdetectorHitNumbers()[3] = nHitsTPC;	  
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
	MCParticle * particle;
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
	trkcol->addElement( newTrack );
	LCRelationImpl * rel = new LCRelationImpl(newTrack,particle,(float)1.0);
	relationcol->addElement( rel );	  
      }
      else { // track is lost
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
