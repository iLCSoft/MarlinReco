#include "V0Finder.h"
#include "marlin/Global.h"
#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/VertexImpl.h"
#include "UTIL/Operators.h"
#include <math.h>
#include <gear/GEAR.h>
#include <gear/BField.h>
#include <gearimpl/Vector3D.h>

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
  
//   std::vector<float> rVertCut;
//   rVertCut.push_back(14.);
//   rVertCut.push_back(60.);
//   rVertCut.push_back(320.);
//   rVertCut.push_back(1600.);

  registerProcessorParameter("CutOnRadius",
			     "Cuts on V0 radius",
			     _rVertCut,
			     float(5.0));

//   std::vector<float> dVertCut;
//   dVertCut.push_back(0.2);
//   dVertCut.push_back(1.0);
//   dVertCut.push_back(1.5);

  
  registerProcessorParameter("CutOnTrkDistance",
			     "Cut on two track distance",
			     _dVertCut,
			     float(1.5));

  registerProcessorParameter("MinimumTrackHitRatio",
			     "Minimum ratio of inner track hit radius to reconstructed vertex radius",
			     _minTrackHitRatio,
			     float(0.7));

  registerProcessorParameter("MassRangeGamma",
			     "Maximal deviation in mass for photon candidate",
			     _deltaMassGamma,
			     float(0.01));

  registerProcessorParameter("MassRangeK0S",
			     "Maximal deviation in mass for K0S candidate",
			     _deltaMassK0S,
			     float(0.01));

  registerProcessorParameter("MassRangeL0",
			     "Maximal deviation in mass for Lamda0 candidate",
			     _deltaMassL0,
			     float(0.008));

  registerProcessorParameter("RxyCutGamma",
			     "Minimum radius in xy plane for photon candidate",
			     _rxyCutGamma,
			     float(10.0));
  registerProcessorParameter("RxyCutK0S",
			     "Minimum radius in xy plane for K0S candidate",
			     _rxyCutK0S,
			     float(30.0));
  registerProcessorParameter("RxyCutGamma",
			     "Minimum radius in xy plane for Lambda0 candidate",
			     _rxyCutLambda,
			     float(50.0));



}

void V0Finder::init() {

  MASSProton  = 0.93827203;
  MASSPion    = 0.13957018;
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
      float d01 = firstTrack->getD0();
      float z01 = firstTrack->getZ0();
      float phi1 = firstTrack->getPhi();
      float tanLambda1 = firstTrack->getTanLambda();
      float omega1 = firstTrack->getOmega();
      HelixClass firstHelix;
      firstHelix.Initialize_Canonical(phi1,d01,z01,omega1,tanLambda1,_bField);
      float charge1 = firstHelix.getCharge();

      float r1 = firstTrack->getRadiusOfInnermostHit();

      for (int j=i+1;j<nelem;++j) {
	Track * secondTrack = dynamic_cast<Track*>(col->getElementAt(j));
	float r2 = secondTrack->getRadiusOfInnermostHit();

	float d02 = secondTrack->getD0();
	float z02 = secondTrack->getZ0();
	float phi2 = secondTrack->getPhi();
	float tanLambda2 = secondTrack->getTanLambda();
	float omega2 = secondTrack->getOmega();
	HelixClass secondHelix;
	secondHelix.Initialize_Canonical(phi2,d02,z02,omega2,tanLambda2,_bField);
	float charge2 = secondHelix.getCharge();
	float prodCharge = charge1*charge2;
	if (prodCharge<0) { // two tracks with opposite charges
	  
	  float px1 = firstHelix.getMomentum()[0];
	  float py1 = firstHelix.getMomentum()[1];
	  float pz1 = firstHelix.getMomentum()[2];
	  float pp1 = sqrt(px1*px1+py1*py1+pz1*pz1);
	  
	  float px2 = secondHelix.getMomentum()[0];
	  float py2 = secondHelix.getMomentum()[1];
	  float pz2 = secondHelix.getMomentum()[2];
	  float pp2 = sqrt(px2*px2+py2*py2+pz2*pz2);
	  
	  float distV0;
	  float momentum[3];
	  float vertex[3];
	  
	  if (pp1>pp2) {
	    distV0 = firstHelix.getDistanceToHelix(&secondHelix, vertex, momentum);
	  }
	  else {
	    distV0 = secondHelix.getDistanceToHelix(&firstHelix, vertex, momentum);
	  }
	  
	  float radV0 = sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]);
	  
// 	  int nRV = int(_rVertCut.size())-1;	  
// 	  float dCut = _dVertCut[0];
// 	  for (int iRV=0;iRV<nRV;++iRV) {
// 	    if (radV0>=_rVertCut[iRV]&&radV0<_rVertCut[iRV+1]) {
// 	      dCut =  _dVertCut[iRV];
// 	      break;
// 	    }
// 	  }

	  // streamlog_out( DEBUG4 ) << " **** found vertex for tracks : " << gear::Vector3D( (const float*) vertex ) 
	  // 		      << " t1 " << lcshort( firstTrack ) << "\n"  
	  // 		      << " t2 " << lcshort( secondTrack )  << std::endl ;



	  // check to ensure there are no hits on tracks at radii significantly smaller than reconstructed vertex
	  // TO DO: should be done more precisely using helices
	  if(r1/radV0<_minTrackHitRatio)continue;
	  if(r2/radV0<_minTrackHitRatio)continue;
	 

	  //	  if (distV0 < _dVertCut && radV0 > _rVertCut ) { // cut on vertex radius and track misdistance
	  if (radV0 > _rVertCut  ) { 

	    streamlog_out( DEBUG4 ) << " ***************** found vertex for tracks : " << gear::Vector3D( (const float*) vertex ) 
				    << " t1 " << lcshort( firstTrack ) << "\n"  
				    << " t2 " << lcshort( secondTrack )  
				    << " distV0 " << distV0 
				    << std::endl ;

	    if( distV0 < _dVertCut ) { // cut on vertex radius and track misdistance
	    

	    // streamlog_out( DEBUG4 ) << " **** found vertex for tracks : " << gear::Vector3D( (const float*) vertex ) 
	    // 			    << " t1 " << lcshort( firstTrack ) << "\n"  
	    // 			    << " t2 " << lcshort( secondTrack )  << std::endl ;

	    streamlog_out( DEBUG ) << "  ***** testing various hypotheses " << std::endl ;

	    // Testing K0 hypothesis
	    float energy1 = sqrt(pp1*pp1+MASSPion*MASSPion);
	    float energy2 = sqrt(pp2*pp2+MASSPion*MASSPion);
	    float energyV0 = energy1 + energy2;	  
	    float massK0 = sqrt(energyV0*energyV0-momentum[0]*momentum[0]-momentum[1]*momentum[1]-momentum[2]*momentum[2]);
	    
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
	    float massL0 = sqrt(energyV0*energyV0-momentum[0]*momentum[0]-momentum[1]*momentum[1]-momentum[2]*momentum[2]);
	  
           // Testing L0bar hypothesis                                             
            if (charge1>0) {
              energy1 = sqrt(pp1*pp1+MASSPion*MASSPion);
              energy2 = sqrt(pp2*pp2+MASSProton*MASSProton);
            }
            else {
              energy1 = sqrt(pp1*pp1+MASSProton*MASSProton);
              energy2 = sqrt(pp2*pp2+MASSPion*MASSPion);
            }
            energyV0 = energy1 + energy2;
            float massL0bar = sqrt(energyV0*energyV0-momentum[0]*momentum[0]-moment\
um[1]*momentum[1]-momentum[2]*momentum[2]);



	    // Testing photon hypothesis	  
	    energyV0 = pp1 + pp2;
	    float massGamma = sqrt(energyV0*energyV0-momentum[0]*momentum[0]-momentum[1]*momentum[1]-momentum[2]*momentum[2]);

	    float deltaK0 = fabs(massK0 - MASSK0S);
	    float deltaL0 = fabs(massL0 - MASSLambda0);
	    float deltaGm = fabs(massGamma - MASSGamma);
	    float deltaL0bar = fabs(massL0bar - MASSLambda0);
	    if(radV0<_rxyCutGamma )deltaGm    = 100000.;
	    if(radV0<_rxyCutK0S   )deltaK0    = 100000.;
	    if(radV0<_rxyCutLambda)deltaL0    = 100000.;
	    if(radV0<_rxyCutLambda)deltaL0bar = 100000.;
	    
	    int code = 22;
	    bool massCondition = false;

           if (deltaGm<deltaL0&&deltaGm<deltaK0&&deltaGm<deltaL0bar) {
              code = 22;
              massCondition = deltaGm < _deltaMassGamma;
            }
            else if (deltaK0<deltaL0 && deltaK0<deltaL0bar) {
              code = 310;
              massCondition = deltaK0 < _deltaMassK0S;
            }
            else{
              if (deltaL0<deltaL0bar ) {
                code = 3122;
                massCondition = deltaL0 < _deltaMassL0;
              }else{
                code = -3122;
                massCondition = deltaL0bar < _deltaMassL0;
              }
            }

	   streamlog_out( DEBUG ) << "  ***** mass condition :  " <<  massCondition 
				  << "  code : " << code  << std::endl ;

	    if (massCondition) {
	      bool ok = true;
	      if(r1/radV0<_minTrackHitRatio|| r2/radV0<_minTrackHitRatio){
		r1 = this->Rmin(firstTrack);
		r2 = this->Rmin(secondTrack);
		if(r1/radV0<_minTrackHitRatio || r2/radV0<_minTrackHitRatio)ok = false;
		//std::cout << " V0X: " << ok << " r = " << radV0 << " r1 = " << r1 << " r2 = " << r2 << std::endl;
	      }
	      if(!ok)continue;
	      TrackPair * trkPair = new TrackPair();
	      trkPair->setFirstTrack( firstTrack );
	      trkPair->setSecondTrack( secondTrack );
	      trkPair->setDistance( distV0 );
	      trkPair->setVertex( vertex );
	      trkPair->setMomentum( momentum );	    
	      trkPair->setCode( code );
	      trkPairs.push_back( trkPair );
	      
	    }
	    else {
// 	      std::cout << "Rejected vertex : V = (" 
// 			<< vertex[0] << ","
// 			<< vertex[1] << ","
// 			<< vertex[2] << ")" << std::endl;
	    }
	    
// 	    std::cout << "Code = " << code << std::endl;
// 	    std::cout << "Vertex = " << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;
// 	    std::cout << "Momentum = " << momentum[0] << " " << momentum[1] << " " << momentum[2] << std::endl;

	    
	  } 
	}
      }
    }

    }//DEBUG ------

//     std::cout << std::endl;

    // Sorting of all vertices in ascending order of the track misdistance
    
    int nTrkPairs = int(trkPairs.size());
    
    if (nTrkPairs>0) { // V0s are present in event

      //      std::cout << "Number of track pairs = " << nTrkPairs << std::endl;
      
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
	  
// 	  std::cout << "Code = " << code << "  Distance = " << distance << std::endl;
// 	  std::cout << "Vertex = (" 
// 		    << vertex[0] << "," 
// 		    << vertex[1] << ","
// 		    << vertex[2] << ")" << std::endl;

// 	  std::cout << "Momentum = ("
// 		    << momentum[0] << ","
// 		    << momentum[1] << ","
// 		    << momentum[2] << ")" << std::endl;
// 	  std::cout << firstTrack << " " << secondTrack << std::endl;

	  
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
	  part->addTrack( firstTrack );
	  part->addTrack( secondTrack );

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

    //    getchar();
    
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

float V0Finder::Rmin( Track* track ) {

   // find track extrema

  float rmin = 1000000.;
  TrackerHitVec hitvec = track->getTrackerHits();
  int nhits = (int)hitvec.size();
  float zmax =-99999.;
  float zmin =99999.;
  for(int ih =0;ih<nhits;++ih){
    float z = (float)hitvec[ih]->getPosition()[2];
    if(z<zmin)zmin=z;
    if(z>zmax)zmax=z;
  }
  float tanLambda = track->getTanLambda();
  //  std::cout << " V0 Check : " << tanLambda << " z = " << zmin << " - " << zmax << std::endl; 
  float zzz = zmin;
  if(tanLambda<0)zzz=zmax;

  float zstart = 0;
  if(fabs(zmin)<fabs(zmax))zstart = zmin;
  if(fabs(zmax)<fabs(zmin))zstart = zmax;
  //std::cout << " V0 Check " << zstart << " - " << zzz << std::endl;
  for(int ih =0;ih<nhits;++ih){
    float z = (float)hitvec[ih]->getPosition()[2];
    if(fabs(z-zstart)<250){
      float x = (float)hitvec[ih]->getPosition()[0];
      float y = (float)hitvec[ih]->getPosition()[1];
      float r2 = x*x+y*y;
      float  r = sqrt(r2);
      if(r<rmin)rmin = r;
    }
    
  }
  
  return rmin;

}
