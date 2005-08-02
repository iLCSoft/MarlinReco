#include "TrackCheater.h"
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include "HelixClass.h"
#include "ClusterShapes.h"
#include <iostream>
#include <map>
#include <vector>

using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

typedef std::map <const MCParticle*,TrackImpl*> map_MCP_Track;
typedef std::map <const MCParticle*,HelixClass*> map_MCP_Helix;

TrackCheater aTrackCheater ;


TrackCheater::TrackCheater() : Processor("TrackCheater") {

  _description = "Creates true tracks..." ;
  
  registerProcessorParameter("TrueTrackCollection",
			     "Collection of True Clusters",
			     _trueTracksCollection ,
			     std::string("TrueTracks"));

  std::vector<std::string> trackerHitCollections;

  trackerHitCollections.push_back(std::string("TPCTrackerHits"));

  registerProcessorParameter("TrackerHitCollections",
			     "Tracker Hit Collection Names",
			     _trackerHitCollections ,
			     trackerHitCollections);

  registerProcessorParameter("BField",
			     "Magnetic Field",
			     _bField ,
			     (float)4.0);

  registerProcessorParameter("ECut",
			     "Energy Cut",
			     _eCut ,
			     (float)0.4);

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
			     (float)30.);    

  registerProcessorParameter("MinimalHits",
			     "Minimal Hits in Track Cluster",
			     _minimal_hits,
			     (int)10);
  

}

void TrackCheater::init() {
    _nRun = -1;
    _nEvt = 0;
}


void TrackCheater::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 

void TrackCheater::processEvent( LCEvent * evt ) { 

    map_MCP_Track _mcp_track;
    map_MCP_Helix _mcp_helix;

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
			  if (charge == 0.) {
			    std::cout << "Track from particle with 0 charge "
				      << std::endl;
			    
			  }
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
			  TrackImpl * track = _mcp_track[mcp];
			  HelixClass * helix = _mcp_helix[mcp];
			  float xPoint[3];
			  float Dist[3];
			  xPoint[0] = (float)trkHit->getPosition()[0];
			  xPoint[1] = (float)trkHit->getPosition()[1];
			  xPoint[2] = (float)trkHit->getPosition()[2];
			  helix->getDistanceToPoint(xPoint,Dist);
			  if (Dist[2] < _hitToHelixCut)
			    track->addHit( trkHit );
			}
		      }
		    }
		    else {
		      std::cout 
			<< "FORMTRUETRACKS : " 
			<< " Warning -> no " 
			<< "pointer to MCParticle " << std::endl;
		    }
		}
		
	    }
	}	
	catch(DataNotAvailableException &e){ 
	}
    }

    LCCollectionVec * trkcol = new LCCollectionVec(LCIO::TRACK);
    LCCollectionVec * relationcol = new LCCollectionVec(LCIO::LCRELATION);
    map_MCP_Track::iterator pos;
    int itk = 0;
    int nlost = 0;
    int nincl = 0;
    for (pos = _mcp_track.begin(); pos != _mcp_track.end(); ++pos) {
	const MCParticle * mcp = pos->first;
	TrackImpl * track = pos->second;
	HelixClass * helix = _mcp_helix[mcp];
	TrackerHitVec hitvec = track->getTrackerHits();
	int nhits = (int)hitvec.size();
	if (nhits > _minimal_hits) {
	  nincl += nhits ;	  
	  float Pos[3];
	  float Mom[3];
	  Pos[0] = (float)mcp->getVertex()[0];
	  Pos[1] = (float)mcp->getVertex()[1];
	  Pos[2] = (float)mcp->getVertex()[2];
	  Mom[0] = mcp->getMomentum()[0];
	  Mom[1] = mcp->getMomentum()[1];
	  Mom[2] = mcp->getMomentum()[2];
	  float charge = mcp->getCharge();
	  if (fabs(charge) < 0.5) 
	    std::cout << "Warning : Track points to neutral particle !!!! "
		      << std::endl;
	  track->setReferencePoint(Pos);
	  track->setD0(helix->getD0());
	  track->setPhi(helix->getPhi0());
	  track->setZ0(helix->getZ0());
	  track->setOmega(helix->getOmega());
	  track->setTanLambda(helix->getTanLambda());
	  TrackerHitVec hitvec = track->getTrackerHits();
	  if (_fitTrueTrack > 0 && nhits > 6) {
	    TrackerHitVec hitvec = track->getTrackerHits();
	    int nhits = (int)hitvec.size();
	    float * xh = new float[nhits];
	    float * yh = new float[nhits];
	    float * zh = new float[nhits];
	    float * ah = new float[nhits];
	    for (int ihit = 0; ihit < nhits; ++ihit) {
	      xh[ihit]=(float)hitvec[ihit]->getPosition()[0];
	      yh[ihit]=(float)hitvec[ihit]->getPosition()[1];
	      zh[ihit]=(float)hitvec[ihit]->getPosition()[2];
	      ah[ihit]=0.0;	    
	    }
	    ClusterShapes * shapes = new ClusterShapes(nhits,ah,xh,yh,zh);
	    float par[5];
	    float dpar[5];
	    float chi2;
	    float distmax;
	    shapes->FitHelix(500, 0, 1, par, dpar, chi2, distmax);
	    if (chi2 < _chi2Cut ) {
	      std::cout << "Successfull fit : " << itk << std::endl; 
	      HelixClass * fittedHelix = new HelixClass();
	      float x0 = par[0];
	      float y0 = par[1];
	      float r0 = par[2];
	      float bz = par[3];
	      float phi0 = par[4];
	      fittedHelix->Initialize_BZ(x0,y0,r0,bz,phi0,
					 _bField,Mom[2],
					 Pos[2]);
	      track->setD0(fittedHelix->getD0());
	      track->setZ0(fittedHelix->getZ0());
	      track->setPhi(fittedHelix->getPhi0());
	      track->setTanLambda(fittedHelix->getTanLambda());
	      track->setOmega(fittedHelix->getOmega());	  
	      float RefPoint[3];
	      RefPoint[0] = fittedHelix->getReferencePoint()[0];
	      RefPoint[1] = fittedHelix->getReferencePoint()[1];
	      RefPoint[2] = fittedHelix->getReferencePoint()[2];	    
	      track->setReferencePoint(RefPoint);
	      delete fittedHelix;
	    }
	    delete shapes;
	    delete[] xh;
	    delete[] yh;
	    delete[] zh;
	    delete[] ah;
	  }
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
	  trkcol->addElement( track );
	  LCRelationImpl * rel = new LCRelationImpl(track,particle,(float)1.0);
	  relationcol->addElement( rel );	  
	}
	else {
	  nlost += nhits;
	  std::cout << "Track is lost : # of hits = " << nhits << std::endl;
	  std::cout << "Px, Py, Pz = " 
		    << " " << mcp->getMomentum()[0]
		    << " " << mcp->getMomentum()[1]
		    << " " << mcp->getMomentum()[2]
		    << std::endl;
	}
	delete helix;
	itk++;
    }


    evt->addCollection(trkcol,_trueTracksCollection.c_str());
    evt->addCollection(relationcol,"TrueTrackToMCP");

    std::cout << "Lost hits = " << nlost << " Accepted hits = " << nincl << std::endl;

  _nEvt++;

}


void TrackCheater::check( LCEvent * evt ) { }
  
void TrackCheater::end(){ } 
