#include "Wolf.h"
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Cluster.h>
#include <EVENT/Track.h>
#include <EVENT/MCParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <iostream>
#include "ClusterShapes.h"
#include <ced_cli.h>
#include "HelixClass.h"
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include <math.h>

using namespace lcio ;
using namespace marlin ;


Wolf aWolf ;


Wolf::Wolf() : Processor("Wolf") {

  _description = "Particle Reconstruction" ;
  
  registerProcessorParameter("TrackCollection",
			     "Track Collection Name",
			     _trackCollection,
			     std::string("TPC_Tracks"));

  registerProcessorParameter("ClusterCollection",
			     "Cluster Collection Name",
			     _clusterCollection,
			     std::string("ClustersAR"));

  registerProcessorParameter("ParticleCollection",
			     "Particle Collection Name",
			     _particleCollection,
			     std::string("RecoParticles"));

  registerProcessorParameter( "ZOfEndcap" ,
			      "Z coordinate of Endcap" ,
			      _zofendcap,
			      (float)2820.);
    
  registerProcessorParameter( "ROfBarrel" , 
			      "Radius of Barrel" , 
			      _rofbarrel,
			      (float)1700.);


  registerProcessorParameter( "NSymmetry" , 
			      "N Fold Symmetry",
			      _nSymmetry,
			      (int)8);

  registerProcessorParameter( "GlobalPhi" , 
			      "Global Phi angle" ,
			      _phiofbarrel,
			      (float)0.);

  registerProcessorParameter( "DistanceTrackToCluster" ,
			      "Distance from Track Seed to Cluster",
			      _distTrackToCluster,
			      (float)50.);

  registerProcessorParameter( "FractionEM" ,
			      "Fraction of EM Energy",
			      _fractionEM,
			      (float)0.95);

  registerProcessorParameter( "RPhiCut", 
			      "Cut on D0 for tracks",
			      _rPhiCut,
			      (float)50.0); 

  registerProcessorParameter( "ZCut", 
			      "Cut on Z0 for tracks",
			      _zCut,
			      (float)50.0); 




}

void Wolf::init() {

  _bField = 4.0;
  _nRun = -1;
  _nEvt = 0;

}


void Wolf::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 

void Wolf::processEvent( LCEvent * evt ) { 


  initialiseEvent( evt );  
  ClusterTrackMatching();
  createPartCollection( evt );
  CleanUp();

  _nEvt++;


}


void Wolf::check( LCEvent * evt ) { }
  
void Wolf::end(){ } 

void Wolf::initialiseEvent( LCEvent * evt ) {

  _clusterVec.clear();
  _trackVec.clear();


  // Reading Track collection
  try {
    LCCollection * col = evt->getCollection(_trackCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0; ielem < nelem; ++ielem) {
      Track * track = dynamic_cast<Track*>(col->getElementAt(ielem));
      TrackExtended * trackAR = new TrackExtended( track );
      _trackVec.push_back( trackAR );
    }    
    //    std::cout << "Tracks : " << nelem  << std::endl;
  }
  catch( DataNotAvailableException &e){}

  // Reading Cluster collection
  try {
    LCCollection * col = evt->getCollection(_clusterCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0; ielem < nelem; ++ielem) {
      Cluster * cluster = dynamic_cast<Cluster*>(col->getElementAt(ielem));
      ClusterExtended * clusterAR = new ClusterExtended( cluster );
      _clusterVec.push_back( clusterAR );
    }
    //    std::cout << "Clusters : " << nelem  << std::endl;
  }
  catch( DataNotAvailableException &e){}


}

void Wolf::createPartCollection( LCEvent * evt) {

    LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);


    int nClusters = (int)_clusterVec.size();
    int nTracks = (int)_trackVec.size();

    //    std::cout << "Tracks : " << nTracks  << std::endl;
    //    std::cout << "Clusters : " << nClusters  << std::endl;



    for (int i=0; i < nTracks; ++i) {
      TrackExtended * trackAR = _trackVec[i];
      Track * track = trackAR->getTrack();
      float d0 = track->getD0();
      float z0 = track->getZ0();
      float omega = track->getOmega();
      float phi0  = track->getPhi();
      float tanlambda = track->getTanLambda();

      if (fabs(d0) < _rPhiCut && fabs(z0) < _zCut) {      
	HelixClass * helix = new HelixClass();
	helix->Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _bField);
	ReconstructedParticleImpl * recPart
	  = new ReconstructedParticleImpl();	
	float Mom[3];
	float mass = 0.140;
	float energy = 0.0;
	int IDPart = 2;
	ClusterExtended * clusterAR = trackAR->getCluster();
	if (clusterAR != NULL) {
	  Cluster * cluster = clusterAR->getCluster();
	  if (cluster != NULL) {
	    recPart->addCluster( cluster );
	    if (cluster->getSubdetectorEnergies().size() == 2) {
	      float EcalEnergy = cluster->getSubdetectorEnergies()[0];
	      float totEnergy = cluster->getEnergy();
	      float fraction = EcalEnergy/fmax(totEnergy,1.0e-6);
	      if (fraction > _fractionEM) {
		IDPart = 1;
		mass = 0.0;
	      }
	    }
	  }
	}
	for (int j=0; j < 3; ++j) {
	  Mom[j] = helix->getMomentum()[j];
	  energy += Mom[j]*Mom[j];
	}	
	energy = sqrt(energy+mass*mass);
	float charge = omega/fabs(omega); 	
	recPart->setMomentum( Mom );
	recPart->setEnergy( energy );
	recPart->setMass( mass );
	recPart->setCharge( charge );
	recPart->setType( IDPart );
	recPart->addTrack( track );
	recparcol->addElement( recPart );
	delete helix;
      }

    }


    for (int i=0; i < nClusters; ++i) {
      ClusterExtended * clusterAR = _clusterVec[i];
      TrackExtendedVec trkvec = clusterAR->getTrackExtendedVec();
      int nTrk = (int)trkvec.size();
      if (nTrk == 0) {
	Cluster * cluster =  clusterAR->getCluster();
	ReconstructedParticleImpl * recPart 
	  = new ReconstructedParticleImpl();
	float Mom[3];
	float mass = 0.498;
	float totGravity = 0.0;
	float totene = cluster->getEnergy();
	for (int ii(0); ii < 3; ++ii) {
	  float xgrav = cluster->getPosition()[ii];
	  Mom[ii] = totene*xgrav;
	  totGravity += xgrav*xgrav;
	}
	totGravity = sqrt(totGravity);
	Mom[0] = Mom[0]/totGravity;
	Mom[1] = Mom[1]/totGravity;
	Mom[2] = Mom[2]/totGravity;
	int IDPart = 4;
	if (cluster->getSubdetectorEnergies().size() == 2) {
	  float EcalEnergy = cluster->getSubdetectorEnergies()[0];
	  float totEnergy = cluster->getEnergy();
	  float fraction = EcalEnergy/fmax(totEnergy,1.0e-6);
	  if (fraction > _fractionEM) { 
	    IDPart = 3;
	    mass = 0.0;
	  }
	}
	float energy = sqrt(Mom[0]*Mom[0]+
			    Mom[1]*Mom[1]+
			    Mom[2]*Mom[2]+
			    mass*mass);
	recPart->setMomentum( Mom );
	recPart->setEnergy( energy );
	recPart->setMass( mass );
	recPart->setCharge(0.);
	recPart->setType(IDPart);
	recPart->addCluster( cluster );
	recparcol->addElement( recPart );	      
	
      }
    }

    evt->addCollection( recparcol , _particleCollection.c_str() );


}

void  Wolf::defineIntersection( TrackExtended * track) {


  HelixClass * helix = new HelixClass();

  TrackerHitVec hitvec = track->getTrack()->getTrackerHits();
  int nhits = (int)hitvec.size();

  float * xh = new float[nhits];
  float * yh = new float[nhits];
  float * zh = new float[nhits];
  float * ah = new float[nhits];

  float zmin = 1.0e+10;
  float zmax = -1.0e+10;
  for (int ihit=0; ihit < nhits; ++ihit) {
    xh[ihit] = (float)hitvec[ihit]->getPosition()[0];
    yh[ihit] = (float)hitvec[ihit]->getPosition()[1];
    zh[ihit] = (float)hitvec[ihit]->getPosition()[2];
    ah[ihit] = 0.0;
    if (zh[ihit] < zmin) 
      zmin = zh[ihit];
    if (zh[ihit] > zmax)
      zmax = zh[ihit];
  }
  float zBegin, zEnd;
  if (fabs(zmin)<fabs(zmax)) {
    zBegin = zmin;
    zEnd   = zmax;
  }
  else {
    zBegin = zmax;
    zEnd   = zmin;
  }
  float signPz = zEnd - zBegin;		  
  ClusterShapes * shapes = new ClusterShapes(nhits,ah,xh,yh,zh);
  float par[5];
  float dpar[5];
  float chi2;
  float distmax;
  shapes->FitHelix(500, 0, 1, par, dpar, chi2, distmax);
  float x0 = par[0];
  float y0 = par[1];
  float r0 = par[2];
  float bz = par[3];
  float phi0 = par[4]; 
  float ref[3];
  float seed[3];
  helix->Initialize_BZ(x0, y0, r0, 
		       bz, phi0, _bField,signPz,
		       zBegin);  

  // * OLD IMPLEMENTATION
  //float phi0  = track->getTrack()->getPhi();
  //float d0    = track->getTrack()->getD0();
  //float omega = track->getTrack()->getOmega();
  //float z0    = track->getTrack()->getZ0();
  //float tanlambda = track->getTrack()->getTanLambda();
  //  for (int j=0; j < 3; ++j)
  //   ref[j]=track->getTrack()->getReferencePoint()[j];
  // helix->Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _bField);  
  //  ref[0] = d0*sin(phi0);
  //  ref[1] = -d0*cos(phi0);
  //  ref[2] = z0;


  // * NEW IMPLEMENTATION
  
  float D0   = helix->getD0();
  float PHI0 = helix->getPhi0();
  float Z0   = helix->getZ0();
  ref[0] = D0*sin(PHI0);
  ref[1] = -D0*cos(PHI0);
  ref[2] = Z0;

  float zLine;

  if (signPz > 0) {
    zLine = _zofendcap;    
  }
  else {
    zLine = -_zofendcap;
  }

  float time_min = helix->getPointInZ(zLine, ref, seed);
  float _const_pi = acos(-1.0);
  float _const_2pi = 2.0*_const_pi;
  float _const_2pi_n = _const_2pi/((float)_nSymmetry);

  if (_nSymmetry > 0) {
    for (int i = 0; i < _nSymmetry; ++i) {
      float phi = _const_2pi_n*((float)i) + _phiofbarrel;
      float xx = _rofbarrel*cos(phi);
      float yy = _rofbarrel*sin(phi);
      float ax = cos(phi + 0.5*_const_pi);
      float ay = sin(phi + 0.5*_const_pi);
      float point[3];
      float tt = helix->getPointInXY(xx,yy,ax,ay,ref,point);
      if (tt < time_min) {
	time_min = tt;
	seed[0] = point[0];
	seed[1] = point[1];
	seed[2] = point[2];      
      }
    }
  }

  delete shapes;
  delete[] xh;
  delete[] yh;
  delete[] zh;
  delete[] ah;
  delete helix;
  track->setSeedPosition(seed);
}


float Wolf::DistanceBetweenPoints(float * x1, float * x2 ) {

  float Dist = 0.0;

  for (int i=0; i < 3; ++i) {
    Dist = Dist + (x1[i]-x2[i])*(x1[i]-x2[i]);
  }

  Dist = sqrt(Dist);

  return Dist;

}

void Wolf::ClusterTrackMatching() {

  int nTracks = (int)_trackVec.size();
  int nClusters = (int)_clusterVec.size();

  for (int i(0); i < nTracks; ++i) {
    TrackExtended * track = _trackVec[i];
    defineIntersection( track );    
    float seed[3];
    seed[0] = track->getSeedPosition()[0];
    seed[1] = track->getSeedPosition()[1];
    seed[2] = track->getSeedPosition()[2];
    float minDist = 10000.;
    ClusterExtended * clusterToAttach = NULL;
    for (int j(0); j < nClusters; ++j) {
      float point[3];
      ClusterExtended * clusterAR = _clusterVec[j];
      CalorimeterHitVec hitvec = clusterAR->getCluster()->getCalorimeterHits();
      int nHits = (int)hitvec.size();
      for (int k(0); k < nHits; ++k) {
	CalorimeterHit * caloHit = hitvec[k];
	point[0] = caloHit->getPosition()[0];
	point[1] = caloHit->getPosition()[1];
	point[2] = caloHit->getPosition()[2];
	float dist = DistanceBetweenPoints(seed,point);	
	if (dist < minDist) {
	  minDist = dist;
	  clusterToAttach = clusterAR;	
	}
      }

    }
    if (minDist < _distTrackToCluster && clusterToAttach != NULL) {
      track->setCluster(clusterToAttach);
      clusterToAttach->addTrackExtended( track );
    }

  }

}

void Wolf::CleanUp() {

  int nTracks = (int)_trackVec.size();
  int nClusters = (int)_clusterVec.size();

  for (int iTrk = 0; iTrk < nTracks; ++iTrk) {
    TrackExtended * track = _trackVec[iTrk];
    delete track;
  }
  _trackVec.clear();

  for (int iClust = 0; iClust < nClusters; ++iClust) {
    ClusterExtended * cluster = _clusterVec[iClust];
    delete cluster;
  }
  _clusterVec.clear();




}
