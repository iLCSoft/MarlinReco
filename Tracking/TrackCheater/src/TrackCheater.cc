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
#include <iostream>
#include <map>
#include <vector>

using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

typedef std::map <const MCParticle*,TrackImpl*> map_MCP_Track;

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
			     float(4.0));

  registerProcessorParameter("ECut",
			     "Energy Cut",
			     _eCut ,
			     float(0.4));

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

    
    LCCollectionVec * trkcol = new LCCollectionVec(LCIO::TRACK);

    int itrack(0);

    for (unsigned int i(0) ; i < _trackerHitCollections.size(); ++i) {    
	try {
	    LCCollection * col = evt->getCollection(_trackerHitCollections[i].c_str());
	    int nelem = col->getNumberOfElements();

	    for (int ielem(0); ielem < nelem; ++ielem) {
		TrackerHit * trkHit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));

		LCObjectVec objVec =  trkHit->getRawHits();
		int nInVec = objVec.size();
		if (nInVec > 0 ) {
		    SimTrackerHit * simTrkHit = dynamic_cast<SimTrackerHit*>(objVec[0]);
		    const MCParticle * mcp = simTrkHit->getMCParticle();
		    
		    if ( mcp != NULL) {
			bool _cuts = mcp->getEnergy() > _eCut ;
			if ( _cuts ) {
			    if (_mcp_track[mcp] == NULL) {
				itrack++;
				TrackImpl * track = new TrackImpl();
				trkcol->addElement( track );
				_mcp_track[mcp] = track;
				track->addHit( trkHit );
			    }
			    else {
				TrackImpl * track = _mcp_track[mcp];
				track->addHit(trkHit);
			    }
			}
		    }
		    else {
			std::cout << "FORMTRUETRACKS : Warning -> no pointer to MCParticle " << std::endl;
		    }
		}

	    }
	}	
	catch(DataNotAvailableException &e){ 
	}
    }


    LCCollectionVec * relationcol = new LCCollectionVec(LCIO::LCRELATION);
    map_MCP_Track::iterator pos;
    for (pos = _mcp_track.begin(); pos != _mcp_track.end(); ++pos) {
	const MCParticle * mcp = pos->first;
	TrackImpl * track = pos->second;
	float Pos[3];
	float Mom[3];
	Pos[0] = (float)mcp->getVertex()[0];
	Pos[1] = (float)mcp->getVertex()[1];
	Pos[2] = (float)mcp->getVertex()[2];
	Mom[0] = mcp->getMomentum()[0];
	Mom[1] = mcp->getMomentum()[1];
	Mom[2] = mcp->getMomentum()[2];
	float charge = mcp->getCharge(); 
	track->setReferencePoint(Pos);
	HelixClass * helix = new HelixClass();
	helix->Initialize_VP(Pos,Mom,charge,_bField);
	track->setD0(helix->getD0());
	track->setPhi(helix->getPhi0());
	track->setZ0(helix->getZ0());
	track->setOmega(helix->getOmega());
	track->setTanLambda(helix->getTanLambda());
	MCParticle * particle;
	try {
	    LCCollection * colpart = evt->getCollection("MCParticle");
	    int nelem = colpart->getNumberOfElements();
	    for (int ielem(0); ielem < nelem; ++ielem) {
		particle = dynamic_cast<MCParticle*>(colpart->getElementAt(ielem));
		if (particle == mcp) 
		    break;
	    }
	}
	catch(DataNotAvailableException &e){}
	LCRelationImpl * rel = new LCRelationImpl(track,particle,(float)1.0);
	relationcol->addElement( rel );
	delete helix;
    }


    evt->addCollection(trkcol,_trueTracksCollection.c_str());
    evt->addCollection(relationcol,"TrueTrackToMCP");

  _nEvt++;

}


void TrackCheater::check( LCEvent * evt ) { }
  
void TrackCheater::end(){ } 
