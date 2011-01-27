#include "CLICTrackSelector.h"
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
#include <gear/CalorimeterParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/VXDParameters.h>
#include <gear/GearParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/BField.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>

using namespace lcio ;
using namespace marlin ;

CLICTrackSelector aCLICTrackSelector ;

CLICTrackSelector::CLICTrackSelector() : Processor("CLICTrackSelector") {  
  _description = "Performs full tracking in LDC detector" ;  


  // Input track collections

  registerInputCollection(LCIO::TRACK,
			  "InputTracks",
			  "Input Track Collection",
			  _inputTrackCollection,
			  std::string("LDCTracks"));

  // Input relation collections

  registerInputCollection(LCIO::LCRELATION,
			  "InputTracksMCPRelColl",
			  "Input Track to MCP Relation Collection Name",
			  _inputTrackMCPCollName,
			  std::string("TPCTracksMCP"));

  
  // Output track collection
  registerOutputCollection(LCIO::TRACK,
			   "SelectedTrackCollection",
			   "Selected track collection name",
			   _selectedTrackCollection,
			   std::string("SelectedLDCTracks"));


  // Output relation collection
  registerOutputCollection(LCIO::LCRELATION,
			   "SelectedTracksMCPRelCollection",
			   "Collection name for the selected track to MCParticle relations",
			   _selectedTrackMCPCollection,
			   std::string("SelectedLDCTracksMCP"));
  
  registerProcessorParameter("CreateMap",
			     "Create Track to MCP Relations",
			     _createMap,
			     int(1));

  registerProcessorParameter("Debug",
			     "Activate debugging?",
			     _debug,
			     int(0));

  registerProcessorParameter("CutOnTPCHits",
			     "Cut on the number of the TPC hits for tracks with no Si hits",
			     _cutOnTPCHits,
			     int(30));

  registerProcessorParameter("CutOnSiHits",
			     "Cut on the number of the Si hits for tracks with no TPC hits",
			     _cutOnSiHits,
			     int(4));
}



void CLICTrackSelector::init() { 

  printParameters();  
  _nRun = -1 ;
  _nEvt = 0 ;
  PI = acos(-1.);
  PIOVER2 = 0.5*PI;
  TWOPI = 2*PI;

}

void CLICTrackSelector::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
  _nEvt = 0;
  std::cout << std::endl;
  std::cout << "CLICTrackSelector ---> new run : run number = " << _nRun << std::endl;

} 

void CLICTrackSelector::processEvent( LCEvent * evt ) { 

  _evt = evt;

  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;

  if (_debug>=1) {
      std::cout << std::endl;
      std::cout << "CLICTrackSelector -> run = " << _nRun 
		<< "  event = " << _nEvt << std::endl;
      std::cout << std::endl;
  }

  LCCollectionVec * colTRK = new LCCollectionVec(LCIO::TRACK);
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  colTRK->setFlag( trkFlag.getFlag()  ) ;  

  LCCollectionVec * colRel = NULL;
  
  if (_createMap) {
    colRel = new LCCollectionVec(LCIO::LCRELATION);
    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    colRel->setFlag( lcFlag.getFlag()  ) ;
  }
  

  // Reading LDC Tracks
  try {
    LCCollection * col = evt->getCollection(_inputTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    if (_debug >= 2) {
      std::cout << std::endl;
      std::cout << "Number of TPC Tracks = " << nelem << std::endl;
      std::cout << " Trk    p     D0         Z0       Px       Py       Pz     Chi2/ndf" << std::endl;
      //           "  0  1.111   0.059      0.022    -0.54     0.61    -0.45    0.185
    }
    int nDropped(0);
    for (int iTrk=0; iTrk<nelem; ++iTrk) {
      bool passTrackSelection = true;
      Track * track = dynamic_cast<Track*>(col->getElementAt(iTrk));
      TrackerHitVec hitVec = track->getTrackerHits();
      int nHits = int(hitVec.size());
      int nHitsVTX = track->getSubdetectorHitNumbers()[6];
      int nHitsFTD = track->getSubdetectorHitNumbers()[7];
      int nHitsSIT = track->getSubdetectorHitNumbers()[8];
      int nHitsTPC = track->getSubdetectorHitNumbers()[9];
      int nHitsSET = track->getSubdetectorHitNumbers()[10];
      int nHitsETD = track->getSubdetectorHitNumbers()[11];
      float r2Min = 1E10; 
      float d0    = track->getD0();
      float z0    = track->getZ0();
      float omega = track->getOmega();
      float tanL  = track->getTanLambda();
      float phi0  = track->getPhi();
      
      HelixClass helix;
      helix.Initialize_Canonical(phi0,d0,z0,omega,tanL,_bField);
      float px = helix.getMomentum()[0];
      float py = helix.getMomentum()[1];
      float pz = helix.getMomentum()[2];
      float pT = sqrt(px*px+py*py);
      float pTot = sqrt(pT*pT+pz*pz);

      float time = this->TimeAtEcal(track);
     
      for (int iH=0;iH<nHits;++iH) {
	TrackerHit * hit = hitVec[iH];
	float hitX = float(hit->getPosition()[0]);
	float hitY = float(hit->getPosition()[1]);
	float hitR2 = hitX*hitX+hitY*hitY;
	if (hitR2<r2Min)r2Min = hitR2;
      }
      
      // general hit count cut
      if( nHitsTPC<_cutOnTPCHits && nHitsVTX+nHitsSIT+nHitsFTD < _cutOnSiHits)passTrackSelection = false;  

      // specific for TPC tracks
      if( nHitsTPC>=_cutOnTPCHits && nHitsVTX+nHitsSIT+nHitsFTD < _cutOnSiHits){	
	int nPass = 0;
	if(nHitsSET>0)nPass++;
	if(nHitsETD>0)nPass++;
	if(nHitsSIT>0)nPass+=2;
	if(nHitsFTD>4)nPass+=2;
	if(nHitsVTX>4)nPass+=2;
	if(r2Min > 500*500)nPass++;
	if(nPass == 0)passTrackSelection = false;
      }

      if(pT<0.100)passTrackSelection = false;
      // require arrival time at ECAL to be consistent with calorimeter hit times
      if(fabs(time)>20)passTrackSelection = false;

      if(!passTrackSelection){
	nDropped++;
	if (_debug >= 1) std::cout << " Dropped " << pTot << " pT = " << pT << " rmin = " << sqrt(r2Min) << " nTPC = " << nHitsTPC << " nVTX = " << nHitsVTX << " nFTD = " << nHitsFTD << " nSIT = " << nHitsSIT << " nSET " << nHitsSET << " nETD " << nHitsETD << " time = " << time  << std::endl;   
	
      }

      if(passTrackSelection){

	if (_debug >= 1) std::cout << " Accepted " << pTot << " pT = " << pT << " rmin = " << sqrt(r2Min) << " nTPC = " << nHitsTPC << " nVTX = " << nHitsVTX << " nFTD = " << nHitsFTD << " nSIT = " << nHitsSIT << " nSET " << nHitsSET << " nETD " << nHitsETD << " time = " << time <<  std::endl;   

	TrackImpl * selectedTrack = new TrackImpl();
	selectedTrack->setOmega(track->getOmega());
	selectedTrack->setTanLambda(track->getTanLambda());
	selectedTrack->setPhi(track->getPhi());
	selectedTrack->setZ0(track->getZ0());
	selectedTrack->setD0(track->getD0());
	selectedTrack->setChi2(track->getChi2());
	selectedTrack->setNdf(track->getNdf());
	selectedTrack->setCovMatrix(track->getCovMatrix());
	int nHits = int(hitVec.size());	
	std::vector<MCParticle*> mcPointers ;
	std::vector<int> mcHits ;
	mcPointers.clear();
	mcHits.clear();
	for (int iH=0;iH<nHits;++iH) {
	  TrackerHit * hit = hitVec[iH];
	  selectedTrack->addHit(hit);
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
      
    

	float RefPoint[3];
	RefPoint[0] = track->getReferencePoint()[0];
	RefPoint[1] = track->getReferencePoint()[1];
	RefPoint[2] = track->getReferencePoint()[2];
	
	selectedTrack->setReferencePoint(RefPoint);
	selectedTrack->setIsReferencePointPCA( track->isReferencePointPCA());
	selectedTrack->setRadiusOfInnermostHit(track->getRadiusOfInnermostHit());
	selectedTrack->subdetectorHitNumbers().resize(12);
	for(unsigned int i=0; i<12; i++)selectedTrack->subdetectorHitNumbers()[i] = track->getSubdetectorHitNumbers()[i];
    
	if (_createMap > 0) {
	  int nRel = int(mcPointers.size());
	  for (int k=0;k<nRel;++k) {
	    LCRelationImpl* tpclcRel = new LCRelationImpl;
	    MCParticle * mcp = mcPointers[k];
	    tpclcRel->setFrom (selectedTrack);
	    tpclcRel->setTo (mcp);
	    float weight = (float)(mcHits[k])/(float)(selectedTrack->getTrackerHits().size());
	    tpclcRel->setWeight(weight);
	    colRel->addElement(tpclcRel);
	  }
	}	
	colTRK->addElement(selectedTrack);
      }

    }
  }
  catch( DataNotAvailableException &e ) {
    std::cout << _inputTrackCollection.c_str() << " collection is unavailable" << std::endl;
  };


  evt->addCollection(colTRK,_selectedTrackCollection.c_str());
  if (_createMap) evt->addCollection(colRel,_selectedTrackMCPCollection.c_str());

  CleanUp();

  _nEvt++;

}



float CLICTrackSelector::TimeAtEcal(Track* pTrack){


  const gear::CalorimeterParameters& pEcalBarrel = Global::GEAR->getEcalBarrelParameters();
  const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
  // determine pseudo-geometry of detector (determined by ECAL barrel)
  // symmetry 0 =cylinder, 1=prorotype, >1 = polygon
  const int symmetry = pEcalBarrel.getSymmetryOrder();
  const float zOfEndCap = (float)pEcalEndcap.getExtent()[2];
  const float phi0 = (float)pEcalBarrel.getPhi0();
  const float rBarrel = (float)pEcalBarrel.getExtent()[0];

  
  HelixClass helix;
  helix.Initialize_Canonical(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(), pTrack->getTanLambda(), _bField);
  
  const EVENT::TrackerHitVec &trackerHitvec(pTrack->getTrackerHits());
  float zMin(std::numeric_limits<float>::max()), zMax(-std::numeric_limits<float>::max());
  
  for (int iz = 0, nTrackHits = trackerHitvec.size(); iz < nTrackHits - 1; ++iz)
    {
      const float hitZ(trackerHitvec[iz]->getPosition()[2]);

      if (hitZ > zMax)
	zMax = hitZ;
      
      if (hitZ < zMin)
	zMin = hitZ;
    }
  
  const int signPz(fabs(zMin) < fabs(zMax) ? 1 : -1);
  
  float referencePoint[3];
  referencePoint[0] = helix.getReferencePoint()[0];
  referencePoint[1] = helix.getReferencePoint()[1];
  referencePoint[2] = helix.getReferencePoint()[2];

  
  // First project to endcap
  float minTime(std::numeric_limits<float>::max());
  bool isProjectedToEndCap(true);
  
  float bestECalProjection[3];
  minTime =  helix.getPointInZ(static_cast<float>(signPz) * zOfEndCap, referencePoint, bestECalProjection);

  // Then project to barrel surface(s)
  static const float pi(acos(-1.));
  float barrelProjection[3];
  
  // n-sided Polygon
  float twopi_n = 2. * pi / (static_cast<float>(symmetry));
  
  for (int i = 0; i < symmetry; ++i)
    {
      float time(std::numeric_limits<float>::max());
      const float phi(twopi_n * static_cast<float>(i) + phi0);
      
      time = helix.getPointInXY(rBarrel * cos(phi), rBarrel * sin(phi),
			   cos(phi + 0.5 * pi), sin(phi + 0.5 * pi), referencePoint, barrelProjection);
      
      if ((time < minTime))
	{
	    minTime = time;
	    isProjectedToEndCap = false;
	    bestECalProjection[0] = barrelProjection[0];
	    bestECalProjection[1] = barrelProjection[1];
	    bestECalProjection[2] = barrelProjection[2];

	}
    }
  
  float tof    = sqrt( bestECalProjection[0]*bestECalProjection[0]+
                       bestECalProjection[1]*bestECalProjection[1]+
                       bestECalProjection[2]*bestECalProjection[2])/300;

  float px = helix.getMomentum()[0];
  float py = helix.getMomentum()[1];
  float pz = helix.getMomentum()[2];
  float E = sqrt(px*px+py*py+pz*pz+0.139*0.139);
  minTime = minTime/300*E-tof;

  return minTime;
  
}


void CLICTrackSelector::CleanUp(){
 

}

void CLICTrackSelector::check(LCEvent * evt) { }

void CLICTrackSelector::end() {

}


