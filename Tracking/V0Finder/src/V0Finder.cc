#include "V0Finder.h"
#include "marlin/Global.h"
#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/VertexImpl.h"
#include <math.h>
#include <gear/GEAR.h>
#include <gear/BField.h>
#include "HelixClass.h"

using namespace lcio ;
using namespace marlin ;


V0Finder aV0Finder ;


V0Finder::V0Finder() : Processor("V0Finder") {

  _description = "V0 Finder Processor " ;
  
  registerInputCollection(LCIO::TRACK,
			  "TrackCollection",
			  "Name of input collection of reconstructed particles",
			  _trackColName,
			  std::string("LDCTracks"));

  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "RecoParticleCollection",
			  "Name of output collection of reconstructed particles",
			  _recoPartColName,
			  std::string("V0RecoParticles"));

  registerOutputCollection(LCIO::VERTEX,
			  "VertexCollection",
			  "Name of output collection of neutral vertices",
			  _vertexColName,
			  std::string("V0Vertices"));

  registerProcessorParameter("CutOnRadius",
			     "Cut on V0 radius",
			     _rVertCut,
			     float(50.));

  registerProcessorParameter("CutOnTrkDistance",
			     "Cut on two track distance",
			     _dVertCut,
			     float(1.0));



}

void V0Finder::init() {

  MASSProton = 0.93827203;
  MASSPion   = 0.13957018;
  MASSLambda0 = 1.115683;
  MASSK0S     = 0.497648;
  MASSGamma   = 0;

  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
  _nRun = -1;
  _nEvt = 0;

}


void V0Finder::processRunHeader( LCRunHeader* run) { 

  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
  _nRun++ ;
  _nEvt = 0;

} 

void V0Finder::processEvent( LCEvent * evt ) { 

  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;

  try {
  
    LCCollection * col = evt->getCollection( _trackColName.c_str() );

    int nelem = col->getNumberOfElements();

    TrackPairVec  trkPairs;
    trkPairs.clear();

    std::map<Track*,int> trackUsed;

    for (int i=0;i<nelem;++i) {
      Track * trk = dynamic_cast<Track*>(col->getElementAt(i));
      trackUsed[trk] = 0;
    }
    
    for (int i=0;i<nelem-1;++i) {
      Track * firstTrack = dynamic_cast<Track*>(col->getElementAt(i));
      HelixClass * firstHelix = new HelixClass();
      float d01 = firstTrack->getD0();
      float z01 = firstTrack->getZ0();
      float phi1 = firstTrack->getPhi();
      float tanLambda1 = firstTrack->getTanLambda();
      float omega1 = firstTrack->getOmega();
      firstHelix->Initialize_Canonical(phi1,d01,z01,omega1,tanLambda1,_bField);
      float charge1 = firstHelix->getCharge();
      for (int j=i+1;j<nelem;++j) {
	Track * secondTrack = dynamic_cast<Track*>(col->getElementAt(j));
	HelixClass * secondHelix = new HelixClass();
	float d02 = secondTrack->getD0();
	float z02 = secondTrack->getZ0();
	float phi2 = secondTrack->getPhi();
	float tanLambda2 = secondTrack->getTanLambda();
	float omega2 = secondTrack->getOmega();
	secondHelix->Initialize_Canonical(phi2,d02,z02,omega2,tanLambda2,_bField);
	float charge2 = secondHelix->getCharge();
	float prodCharge = charge1*charge2;
	if (prodCharge<0) { // two tracks with opposite charges
	  
	  float px1 = firstHelix->getMomentum()[0];
	  float py1 = firstHelix->getMomentum()[1];
	  float pz1 = firstHelix->getMomentum()[2];
	  float pp1 = sqrt(px1*px1+py1*py1+pz1*pz1);
	  
	  float px2 = secondHelix->getMomentum()[0];
	  float py2 = secondHelix->getMomentum()[1];
	  float pz2 = secondHelix->getMomentum()[2];
	  float pp2 = sqrt(px2*px2+py2*py2+pz2*pz2);
	  
	  float distV0;
	  float momentum[3];
	  float vertex[3];
	  
	  if (pp1>pp2) {
	    distV0 = firstHelix->getDistanceToHelix(secondHelix, vertex, momentum);
	  }
	  else {
	    distV0 = secondHelix->getDistanceToHelix(firstHelix, vertex, momentum);
	  }
	  
	  float radV0 = sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]);
	  
	  if (distV0 < _dVertCut && radV0 > _rVertCut) { // cut on vertex radius and track misdistance
	    
	    TrackPair * trkPair = new TrackPair();
	    trkPair->setFirstTrack( firstTrack );
	    trkPair->setSecondTrack( secondTrack );
	    trkPair->setDistance( distV0 );
	    trkPair->setVertex( vertex );
	    trkPair->setMomentum( momentum );
	    
	    // Testing K0 hypothesis
	    float energy1 = sqrt(pp1*pp1+MASSPion*MASSPion);
	    float energy2 = sqrt(pp2*pp2+MASSPion*MASSPion);
	    float energyV0 = energy1 + energy2;	  
	    float massK02 = energyV0*energyV0-momentum[0]*momentum[0]-momentum[1]*momentum[1]-momentum[2]*momentum[2];
	    
	    // Testing L0 hypothesis
	    if (charge1<0) {
	      energy1 = sqrt(pp1*pp1+MASSPion*MASSPion);
	      energy2 = sqrt(pp2*pp2+MASSProton*MASSProton);
	    }
	    else {
	      energy1 = sqrt(pp1*pp1+MASSProton*MASSProton);
	      energy2 = sqrt(pp2*pp2+MASSPion*MASSPion);
	    }
	    energyV0 = energy1 + energy2;	  
	    float massL02 = energyV0*energyV0-momentum[0]*momentum[0]-momentum[1]*momentum[1]-momentum[2]*momentum[2];
	  
	    // Testing photon hypothesis	  
	    energyV0 = pp1 + pp2;
	    float massGamma2 = energyV0*energyV0-momentum[0]*momentum[0]-momentum[1]*momentum[1]-momentum[2]*momentum[2];

	    float deltaK0 = fabs(massK02 - MASSK0S*MASSK0S);
	    float deltaL0 = fabs(massL02 - MASSLambda0*MASSLambda0);
	    float deltaGm = fabs(massGamma2 - MASSGamma*MASSGamma);
	    
	    int code = 22;

	    if (deltaGm<deltaL0&&deltaGm<deltaK0) {
	      code = 22;
	    }
	    else if (deltaK0<deltaL0) {
	      code = 310;
	    }
	    else {
	      code = 3122;
	    }
	    trkPair->setCode( code );
	    
	    trkPairs.push_back( trkPair );
	    
	  } 
	}
	delete secondHelix;
      }
      delete firstHelix;
    }
    
    // Sorting of all vertices in ascending order of the track misdistance
    
    int nTrkPairs = int(trkPairs.size());
    
    if (nTrkPairs>0) { // V0s are present in event
      
      Sorting( trkPairs );
      
      // Declaration of the output collections
      LCCollectionVec * colRecoPart = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
      LCCollectionVec * colVertex   = new LCCollectionVec(LCIO::VERTEX);
      
      for (int iTrkP=0;iTrkP<nTrkPairs;++iTrkP) {
	TrackPair * pair = trkPairs[iTrkP];
	Track * firstTrack = pair->getFirstTrack();
	Track * secondTrack = pair->getSecondTrack();
	if (trackUsed[firstTrack]==0&&trackUsed[secondTrack]==0) {
	  
	  ReconstructedParticleImpl * part = new ReconstructedParticleImpl();
	  VertexImpl * vtx = new VertexImpl();
	  
	  float vertex[3];
	  float momentum[3];
	  int code = pair->getCode();
	  for (int iC=0;iC<3;++iC) {
	    vertex[iC] = pair->getVertex()[iC];
	    momentum[iC] = pair->getMomentum()[iC];
	  }
	  
	  float distance = pair->getDistance();	
	  vtx->setPosition( vertex );
	  vtx->addParameter( distance );
	  
	  part->setMomentum( momentum );
	  part->setType( code );
	  
	  
	  float mass = 0;
	  if ( code == 22)
	    mass = 0;
	  else if ( code == 310 )
	    mass = MASSK0S;
	  else 
	    mass = MASSLambda0;
	  
	  part->setMass( mass );	
	  vtx->setAssociatedParticle( part );
	  part->setStartVertex( vtx );
	  
	  colRecoPart->addElement( part );
	  colVertex->addElement( vtx );
	  
	  trackUsed[firstTrack] = 1;
	  trackUsed[secondTrack] = 1;
	}
      }
      
      evt->addCollection( colRecoPart,_recoPartColName.c_str() );
      evt->addCollection( colVertex, _vertexColName.c_str() );
      
    }
    
    // Clean up memory
    for (int iTrkP=0;iTrkP<nTrkPairs;++iTrkP) {
      TrackPair * trkPair = trkPairs[iTrkP];
      delete trkPair;
    }
    trkPairs.clear();
    
  }
  catch(DataNotAvailableException &e) {}

  _nEvt++;

}


void V0Finder::check( LCEvent * evt ) { }
  
void V0Finder::end(){ } 

void V0Finder::Sorting( TrackPairVec & trkPairVec ) {

  int sizeOfVector = int(trkPairVec.size());
  TrackPair *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++)
      {
	one = trkPairVec[j];
	two = trkPairVec[j+1];
	if( one->getDistance() > two->getDistance() )
	  {
	    Temp = trkPairVec[j];
	    trkPairVec[j] = trkPairVec[j+1];
	    trkPairVec[j+1] = Temp;
	  }
      }  

}
