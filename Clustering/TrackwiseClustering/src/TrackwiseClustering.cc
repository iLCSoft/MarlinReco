#include "TrackwiseClustering.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/Track.h"
#include "IMPL/ClusterImpl.h"
#include "IMPL/LCCollectionVec.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include "ClusterShapes.h"
// GEAR include files
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/CalorimeterParameters.h>
#include "random.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;


TrackwiseClustering aTrackwiseClustering;

TrackwiseClustering::TrackwiseClustering() : Processor("TrackwiseClustering") {


  _description = "Performs clustering in a track-wise manner..." ;


     // Not used  
     registerProcessorParameter( "DistanceForDirection", 
				"Distance to Define Direction", 
				_distanceToDefineDirection,
				(float)5.);
     // Not used
     registerProcessorParameter( "DistanceToTrackSeed", 
				"Distance to Track Seed", 
				_distanceToTrackSeed,
				(float)25.);    


     // RCutMax
    std::vector<float>  distanceTrackBack;
    distanceTrackBack.push_back(100.);
    distanceTrackBack.push_back(500.);
    
    registerProcessorParameter( "DistanceTrackBack" , 
				"Distance to Track Back "  ,
				_distanceTrackBack,
				 distanceTrackBack); 


    // RCut
    std::vector<float>  stepTrackBack;
    stepTrackBack.push_back(10.0);
    stepTrackBack.push_back(100.0);
    
    registerProcessorParameter( "StepTrackBack" , 
				"Step to Track Back "  ,
				_stepTrackBack,
				 stepTrackBack); 

    // SCut (merging parameter cut)
    std::vector<float>  resolutionParameter;
    resolutionParameter.push_back(20.0);
    resolutionParameter.push_back(80.0);
    
    registerProcessorParameter( "ResolutionParameter" , 
				"Resolution Parameter "  ,
				_resolutionParameter,
				 resolutionParameter); 


    std::vector<float>  distanceMergeForward;
    distanceMergeForward.push_back(50.0);
    distanceMergeForward.push_back(100.0);
    
    registerProcessorParameter( "DistanceMergeForward" , 
				"Distance To Merge Forward" ,
				_distanceMergeForward,
				 distanceMergeForward); 


    registerProcessorParameter( "NToDefineSP",
				"N hits to define SP " , 
				_NDefineSP, 
				5);

    registerProcessorParameter( "NScanToMergeForward",
				"N hits scan to merge forward " , 
				_nScanToMergeForward, 
				10);

    registerProcessorParameter( "TypeOfGenericDistance" , 
				"Type of Generic Distance "  ,
				_typeOfGenericDistance,
				 0);


    std::vector<std::string> EcalCollections;
    EcalCollections.push_back(std::string("ECAL"));
    
    registerInputCollections( LCIO::CALORIMETERHIT,
			     "EcalCollections" , 
			     "Ecal Collection Names "  ,
			     _ecalCollections,
			     EcalCollections);

    std::vector<std::string> HcalCollections;
    HcalCollections.push_back(std::string("HCAL"));
    
    registerInputCollections( LCIO::CALORIMETERHIT,
			     "HcalCollections" , 
			     "Hcal Collection Names "  ,
			     _hcalCollections,
			     HcalCollections);
    
    std::vector<std::string> TrackCollections;
    TrackCollections.push_back(std::string("Track"));
    
    registerInputCollections( LCIO::TRACK,
			      "TrackCollections" , 
			      "Track Collection Names "  ,
			      _trackCollections,
			      TrackCollections);
    
    
    registerOutputCollection( LCIO::CLUSTER,
			      "ClusterCollection" , 
			      "Cluster Collection Name "  ,
			      _clusterCollection,
			      std::string("ClustersAR"));

    
   registerProcessorParameter( "MinimalHitsInCluster" ,
			       "Minimal allowed hits in cluster" , 
				_nhit_minimal, 
				10);
    
   registerProcessorParameter( "MaximalHitsToMerge" ,
			       "Maximal Hits To Merge" , 
				_nhit_merge_forward, 
				50);
    
   registerProcessorParameter( "UseTracking" ,
			       "Use tracks to seed clusters" , 
				_use_tracks, 
				0);

   registerProcessorParameter( "DoMergingLowMultiplicity" , 
			       "merging low multiplicity clusters?", 
			       _doMerging,
			       1);

   registerProcessorParameter( "DoMergingForward" , 
			       "merging clusters forward-wise?", 
			       _doMergingForward,
			       1);

   registerProcessorParameter( "DisplayClusterInfo",
			       "Display Info on Clusters",
			       _displayClusters,
			       0);

   registerProcessorParameter( "ResolutionToMerge",
			       "Resolution To Merge Halo Hits",
			       _resolutionToMerge,
			       (float)400.);


   registerProcessorParameter( "WeightForResolution",
			       "Weight For Resolution",
			       _weightForReso,
			       (float)1.0);

   registerProcessorParameter( "WeightForDistance",
			       "Weight For Distance",
			       _weightForDist,
			       (float)1.0);


   registerProcessorParameter( "BField",
			       "Magnetic Field (in TESLA)",
			       _bField,
			       float(4.0));

}

void TrackwiseClustering::init() {

    _nRun = -1;

}

void TrackwiseClustering::processRunHeader( LCRunHeader* run ) {

    _nRun++;
    _nEvent = 0;
    _allHits.clear();
    _allTracks.clear();
    _allSuperClusters.clear();
    _allClusters.clear();

}

void TrackwiseClustering::processEvent( LCEvent * evt ) {

    initialiseEvent( evt );
    GlobalSorting();
    GlobalClustering();
    //    propertiesForAll();

    if (_doMergingForward == 1) {    
      //      MergeTrackSegments();
      mergeForward();
    }

    if (_doMerging == 1) {
      mergeLowMultiplicity();
    }

    //    propertiesForAll();

    if (_displayClusters == 1)
      DisplayClusters(_allClusters);
    CreateClusterCollection(evt,_allClusters);
    CleanUp();
    _nEvent++;

}

void TrackwiseClustering::check( LCEvent * evt ) {}

void TrackwiseClustering::end() {}

void TrackwiseClustering::initialiseEvent( LCEvent * evt ) {

    int Total(0);

    int NOfHits(0);

  const gear::CalorimeterParameters& pEcalBarrel = Global::GEAR->getEcalBarrelParameters();

  const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();

  _rofbarrel = (float)pEcalBarrel.getExtent()[0];
  _phiofbarrel = (float)pEcalBarrel.getPhi0();
  _nsymmetry = pEcalBarrel.getSymmetryOrder();
  _zofendcap = (float)pEcalEndcap.getExtent()[2];
  _const_pi    = acos(-1.);
  _const_2pi = 2.0*_const_pi;
  _const_pi_n  = _const_pi/float(_nsymmetry);
  _const_2pi_n = 2.0*_const_pi/float(_nsymmetry);
  _thetaofendcap = (float)atan((double)(_rofbarrel/_zofendcap));


    
    _xmin_in_distance = 1.0e+10;
    _xmax_in_distance = -1.0e+10;


//  Reading Ecal hits
    for (unsigned int i(0); i < _ecalCollections.size(); ++i) {
	try {
	    LCCollection * col = evt->getCollection( _ecalCollections.at(i).c_str() );
	    if (col != 0) {
		int nelem = col->getNumberOfElements()  ;
		Total += nelem;
		for (int j(0); j < nelem; ++j) {
		    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(j));
		    CaloHitExtended *calohit = new CaloHitExtended(hit,0);
		    float dist[2];
		    CalculateGenericDistance(calohit, dist);		    
		    if (_typeOfGenericDistance == 0) 
		      calohit->setGenericDistance(dist[0]);
		    else 
		      calohit->setGenericDistance(dist[1]);
		    calohit->setDistanceToCalo(dist[1]);
		    float distance = calohit->getGenericDistance();
		    if (distance < _xmin_in_distance) 
			_xmin_in_distance = distance;
		    if (distance > _xmax_in_distance) 
			_xmax_in_distance = distance;		
		    _allHits.push_back(calohit);
		    NOfHits++;
		}
	    }
	}
	catch(DataNotAvailableException &e) {};

    }

// Reading Hcal hits
    for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
	try{
	    LCCollection * col = evt->getCollection( _hcalCollections.at(i).c_str() );
	    if (col != 0) {
		int nelem = col->getNumberOfElements()  ;
		Total+= nelem;
		for (int j(0); j < nelem; ++j) {
		    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(j));
		    CaloHitExtended * calohit = new CaloHitExtended(hit,1);
		    float  dist[2] ;
		    CalculateGenericDistance(calohit,dist);
		    if (_typeOfGenericDistance == 0)
		      calohit->setGenericDistance(dist[0]);
		    else 
		      calohit->setGenericDistance(dist[1]);
		    calohit->setDistanceToCalo(dist[1]);
		    float distance = calohit->getGenericDistance();
		    if (distance < _xmin_in_distance) 
			_xmin_in_distance = distance;
		    if (distance > _xmax_in_distance) 
			_xmax_in_distance = distance;
		    _allHits.push_back(calohit);
		    NOfHits++;
		}
	    }
	}
	catch(DataNotAvailableException &e) {};

    }


// Reading Track collection
    if (_use_track != 0) {
	for (unsigned int i(0); i < _trackCollections.size(); ++i) {
	    try {
		LCCollection * col = evt->getCollection( _trackCollections.at(i).c_str() ); 
		if (col != 0) {
		    int nelem = col->getNumberOfElements()  ;
		    for (int j(0); j < nelem; j++) {
			Track * track = dynamic_cast<Track*>(col->getElementAt(j));
			TrackExtended * trackAR = new TrackExtended(track);
			_allTracks.push_back(trackAR);
		    }
		}
	    }
	    catch(DataNotAvailableException &e) {};
	} 
    }

}

float TrackwiseClustering::findResolutionParameter(CaloHitExtended * fromHit, CaloHitExtended * toHit) {

    float xdistvec[3];
    float dirvec[3];
    float xdist(0.);
    float product(0.);
    float dir(0.);
    for (int i(0); i < 3; i++) {
	xdistvec[i] = toHit->getCalorimeterHit()->getPosition()[i] - fromHit->getCalorimeterHit()->getPosition()[i];
	xdist += xdistvec[i]*xdistvec[i];
	dirvec[i] = fromHit->getDirVec()[i];
	dir += dirvec[i]*dirvec[i];
	product += xdistvec[i]*dirvec[i]; 
    }

    xdist = sqrt(xdist);
    dir = sqrt(dir);
    product=product/fmax(1.0e-6,(xdist*dir));

    if (product > 1.) {
      product = 0.999999;
    }
    if (product < -1.) {
      product = -0.999999;
    }

    float angle = acos(product);
    
    return xdist*angle;
 
}

void TrackwiseClustering::CalculateGenericDistance(CaloHitExtended *calohit, float * dist) {
    float xDistance =0.0;
    float rDistance =0.0;

    for (int i(0); i < 3; ++i) {
	float x = calohit->getCalorimeterHit()->getPosition()[i];
	rDistance += x*x; 	
    }
    rDistance = sqrt(rDistance);
    
    float x = calohit->getCalorimeterHit()->getPosition()[0];
    float y = calohit->getCalorimeterHit()->getPosition()[1];
    float z = calohit->getCalorimeterHit()->getPosition()[2];
    float phi = atan2(y,x) - _phiofbarrel + _const_pi_n;
    int nZone = (int)(phi/_const_2pi_n);
    if (phi < 0.)
      phi = phi + _const_2pi;
    phi = phi - nZone * _const_2pi_n - _const_pi_n;
    float radius = sqrt(x*x + y*y);
    float rdist = radius * cos(phi) - _rofbarrel;
    float zdist = fabs(z) - _zofendcap;
    if (rdist > 0 && zdist < 0) {
      xDistance = rdist;
    }
    else if (rdist < 0 && zdist > 0 ) {
      xDistance = zdist;
    }
    else {
      float theta = (float)atan((float)(rdist/zdist));
      if (theta > _thetaofendcap) {
	xDistance = rdist;
      }
      else {
	xDistance = zdist;
      }

    }
    xDistance = xDistance + 1.0e-10*rDistance ;
    dist[0] = rDistance;
    dist[1] = xDistance;
}



float TrackwiseClustering::DistanceBetweenPoints(float *x1, float *x2) {
    float xDistance(0.);
    for (int i(0); i < 3; i++) {
	float xx = x1[i] - x2[i];
	xDistance += xx*xx ;
    }
    xDistance = sqrt(xDistance);
    return xDistance;
}

void TrackwiseClustering::BubbleSort(CaloHitExtendedVec & input)
{
    unsigned int sizeOfVector = input.size();
    CaloHitExtended *one,*two,*Temp;

    for (unsigned int i = 0 ; i < sizeOfVector-1; i++)
	for (unsigned int j = 0; j < sizeOfVector-i-1; j++)
	{
	    one = input.at(j);
	    two = input.at(j+1);
	    if(one->getGenericDistance() > two->getGenericDistance())
	    {
		Temp = input[j];
		input[j] = input[j+1];
		input[j+1] = Temp;
	    }
	}
}

void TrackwiseClustering::CleanUp() {

    for (unsigned int i(0); i < _allHits.size(); ++i) {
	CaloHitExtended * calohit = _allHits.at(i);
	delete calohit;
    }
    
    _allHits.clear();

    for (unsigned int i(0); i < _allTracks.size(); ++i) {
	TrackExtended * track = _allTracks.at(i);
	delete track;
    }
    
    _allTracks.clear();


    for (unsigned int i(0); i < _allSuperClusters.size(); ++i) {
	ClusterExtended * cluster = _allSuperClusters.at(i);
	delete cluster;
    }

    _allSuperClusters.clear();

    for (unsigned int i(0); i < _allClusters.size(); ++i) {
	ClusterExtended * cluster = _allClusters.at(i);
	delete cluster;
    }

    _allClusters.clear();

}

void TrackwiseClustering::propertiesForAll() {
  int nclusters = (int)_allClusters.size();
  for (int i=0;i<nclusters;++i) {
    ClusterExtended * Cl = _allClusters[i];
    calculateProperties(Cl);
  }

}

void TrackwiseClustering::calculateProperties(ClusterExtended * Cl) {

  CaloHitExtendedVec calohitvec = Cl->getCaloHitExtendedVec();
  int nhcl = (int)calohitvec.size();
  if (nhcl > 0) {
    float * xhit = new float[nhcl];
    float * yhit = new float[nhcl];
    float * zhit = new float[nhcl];
    float * ahit = new float[nhcl];
    float * exhit = new float[nhcl];
    float * eyhit = new float[nhcl];
    float * ezhit = new float[nhcl];    
    float totene = 0.0;
    float totecal = 0.0;
    float tothcal = 0.0;
    RandomNumberGenerator random;
    float zmin = 1.0e+20;
    float zmax = -1.0e+20;
    int jhit = 0;
    for (int ihit(0); ihit < nhcl; ++ihit) {
      CalorimeterHit * calhit = 
	calohitvec[ihit]->getCalorimeterHit();
      if (calohitvec[ihit]->getDistanceToNearestHit() < 100.) {
	xhit[jhit] = calhit->getPosition()[0] + random.EqualDistribution(1.0)[0];
	yhit[jhit] = calhit->getPosition()[1] + random.EqualDistribution(1.0)[0];
	zhit[jhit] = calhit->getPosition()[2];
	ahit[jhit] = calhit->getEnergy();
	exhit[jhit] = 4.0;
	eyhit[jhit] = 4.0;
	ezhit[jhit] = 4.0;
	totene += ahit[jhit];
	if (calohitvec[jhit]->getType() == 0) {
	  totecal += ahit[jhit];
	}
	else {
	  tothcal += ahit[jhit];
	}	
	if (zhit[jhit]<zmin )
	  zmin = zhit[jhit];
	if (zhit[jhit]>zmax)
	  zmax = zhit[jhit];	
	jhit++;
      }	      
    }
    float zBeg;
    float signPz = 1.0;
    if (fabs(zmin)>fabs(zmax)) {
      signPz = -1.0;
      zBeg = zmax;
    }
    else {
      zBeg = zmin;
    }


    ClusterShapes * shapes 
      = new ClusterShapes(jhit,ahit,xhit,yhit,zhit);	    
    shapes->setErrors(exhit,eyhit,ezhit);
    float par[5];
    float dpar[5];
    float chi2 = 1.0e+10;
    float distmax = 1.0e+20;
    float x0 = 1;
    float y0 = 1;
    float r0 = 1;
    float bz = 1;
    float phi0 = 1;
    if (jhit > 3) {
      shapes->FitHelix(500, 0, 1, par, dpar, chi2, distmax);
      x0 = par[0];
      y0 = par[1];
      r0 = par[2];
      bz = par[3];
      phi0 = par[4];
    }
    HelixClass helix;
    helix.Initialize_BZ(x0, y0, r0, 
			bz, phi0, _bField, signPz,
			zBeg);
    
    Cl->setHelix(helix);
    Cl->setHelixChi2R(chi2);
    Cl->setHelixChi2Z(distmax);

    float axis[3];
    float pos[3];
    float axisMod = 0.0;
    float low[3];
    float up[3];
    for (int i=0; i<3; ++i) {
      pos[i]  = shapes->getCentreOfGravity()[i];
      axis[i] = shapes->getEigenVecInertia()[i];
      axisMod += axis[i]*axis[i];
      low[i] = Cl->getLowEdge()[i];
      up[i]  = Cl->getUpEdge()[i];
    }
    axisMod = sqrt(axisMod);
    Cl->setAxis(axis);
    Cl->setPosition(pos);	    
    Cl->setEccentricity(shapes->getElipsoid_eccentricity());

//    Cl->setLowEdge(low);
//    Cl->setUpEdge(up);
    if (nhcl > 40  ) {      
      //debug
      /*
	std::cout << nhcl << " " << shapes->getElipsoid_eccentricity() << std::endl;
	std::cout << "Chi2 : " << chi2 << " DistMax = " << distmax << std::endl;
	std::cout << "Pos : " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	std::cout << "Low : " << low[0] << " " << low[1] << " " << low[2] << std::endl;
	std::cout << "Up : " << up[0] << " " << up[1] << " " << up[2] << std::endl;     
	std::cout << "Par : " << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << std::endl;
      */

       float xx[3];
       int nn = nhcl / 2;
       xx[0] = xhit[nn];
       xx[1] = yhit[nn];
       xx[2] = zhit[nn];
       float dd[3];
       float time = helix.getDistanceToPoint(xx,dd);
       float dist1 = dd[2];
    }

    delete shapes;
    delete[] exhit;
    delete[] eyhit;
    delete[] ezhit;
    delete[] xhit;
    delete[] yhit;
    delete[] zhit;
    delete[] ahit;
    
  }

}

void TrackwiseClustering::CreateClusterCollection(LCEvent * evt, ClusterExtendedVec clusterVec) {

    int nclust = (int)clusterVec.size();

    LCCollectionVec * clucol = new LCCollectionVec(LCIO::CLUSTER);

    for (int iclust(0); iclust < nclust; ++iclust) {
	ClusterExtended * Cl = clusterVec[iclust];
	CaloHitExtendedVec calohitvec = Cl->getCaloHitExtendedVec();
	int nhcl = (int)calohitvec.size();
	if (nhcl > _nhit_minimal) {
	    ClusterImpl * cluster = new ClusterImpl();
	    float * xhit = new float[nhcl];
	    float * yhit = new float[nhcl];
	    float * zhit = new float[nhcl];
	    float * ahit = new float[nhcl];
	    float totene = 0.0;
	    float totecal = 0.0;
	    float tothcal = 0.0;
	    for (int ihit(0); ihit < nhcl; ++ihit) {
		CalorimeterHit * calhit = 
		  calohitvec[ihit]->getCalorimeterHit();
		cluster->addHit(calhit,(float)1.0);
		xhit[ihit] = calhit->getPosition()[0];
		yhit[ihit] = calhit->getPosition()[1];
		zhit[ihit] = calhit->getPosition()[2];
		ahit[ihit] = calhit->getEnergy();
		totene += ahit[ihit];
		if (calohitvec[ihit]->getType() == 0) {
		  totecal += ahit[ihit];
		}
		else {
		  tothcal += ahit[ihit];
		}		
	    }
	    
	    ClusterShapes * shape 
	      = new ClusterShapes(nhcl,ahit,xhit,yhit,zhit);	    
	    cluster->setEnergy(totene);
	    cluster->subdetectorEnergies().resize(2);
	    cluster->subdetectorEnergies()[0] = totecal;
	    cluster->subdetectorEnergies()[1] = tothcal;	    
	    cluster->setPosition(shape->getCentreOfGravity());
	    float PhiCluster = atan2(shape->getEigenVecInertia()[1],shape->getEigenVecInertia()[0]);
	    float ThetaCluster = acos(shape->getEigenVecInertia()[2]);
	    cluster->setIPhi(PhiCluster);
	    cluster->setITheta(ThetaCluster);	    
	    clucol->addElement(cluster);

	    delete shape;
	    delete[] xhit;
	    delete[] yhit;
	    delete[] zhit;
	    delete[] ahit;

	}

    }
    evt->addCollection(clucol,_clusterCollection.c_str());

}


void TrackwiseClustering::DisplayClusters(ClusterExtendedVec clusterVec) {

    cout << " " << endl;
    cout << "Debugging --> Event : " << _nEvent << endl;
    int nhits  = _allHits.size();
    cout << "Total number of Hits : " << nhits << endl;
    int nclust = clusterVec.size();
    cout << "Number of clusters : " << nclust << endl;
    int ntot(0);
    int ninbig(0);

    

    float totenergy = 0.0;
    float energyinbig = 0.0;
    for (int iclust(0); iclust < nclust; ++iclust) {
	ClusterExtended * Cl = clusterVec[iclust];
	CaloHitExtendedVec calohitvec = Cl->getCaloHitExtendedVec();
	int nhcl = calohitvec.size();
	ntot += nhcl;
	float ene=0.0;
	for (int i=0; i<nhcl; ++i) {
	  ene += calohitvec[i]->getCalorimeterHit()->getEnergy();
	}
	totenergy += ene;
	if (nhcl > _nhit_minimal) {
	    cout << "Cluster  " << iclust << " Number of hits = " << nhcl << endl;
	    ninbig += nhcl;
	    energyinbig += ene;
	}

    }

    float fraction = energyinbig/totenergy;
    cout << endl;
    cout << "Hits in big clusters     : " << ninbig << endl;
    cout << "Total energy : " << totenergy << std::endl;
    cout << "Fraction in big clusters : " << fraction << std::endl;
    cout << "Sumcheck: number of hits : " <<  ntot << endl;



}



void TrackwiseClustering::GlobalSorting() {

    int _NGLAYERS = 1000;

    float inverse = ((float)_NGLAYERS)/(_xmax_in_distance - _xmin_in_distance);

    _allSuperClusters.clear();

    for (int i=0; i < _NGLAYERS; ++i) {
	ClusterExtended * cluster = new ClusterExtended();
	_allSuperClusters.push_back(cluster);
    }

    for (unsigned int i = 0; i < _allHits.size(); ++i) {
	CaloHitExtended * calohit = _allHits[i];
	float distToIP = calohit->getGenericDistance();
	int index = (int)((distToIP - _xmin_in_distance) * inverse) ;
	if (index >= _NGLAYERS)
	    index = _NGLAYERS - 1;
	if (index < 0) 
	    index = 0;

	_allSuperClusters[index]->addCaloHitExtended(calohit);
    }

    _allHits.clear();

    int counter = 0;

    for (int i(0); i < _NGLAYERS; ++i) {
	ClusterExtended * cluster = _allSuperClusters[i];
	CaloHitExtendedVec calohitvec = cluster->getCaloHitExtendedVec();
	int nhits = calohitvec.size();
	if (nhits > 0) {
	    BubbleSort(calohitvec);
	    for (int ihit(0); ihit < nhits; ++ihit) {
		CaloHitExtended * calohit = calohitvec[ihit];
		calohit->setIndex(counter);
		_allHits.push_back(calohit);
		counter++;
	    }
	}

    }

}


void TrackwiseClustering::GlobalClustering() {

    _allClusters.clear();

    for (unsigned int ihitTo(0); ihitTo < _allHits.size(); ++ihitTo) {

	CaloHitExtended * CaloHitTo = _allHits[ihitTo];
	int ihitFrom = ihitTo - 1;
	int ifound = 0;
	int idTo = CaloHitTo->getType();

	float r_step = _stepTrackBack[idTo];
	float r_min  = r_step;
	float r_dist = _distanceTrackBack[idTo]; 

	float YResMin = 1.0e+10;
	float YResCut = _resolutionParameter[idTo];
	float YDistMin = 1.0e+10;

	while (ihitFrom >=0) {

	    CaloHitExtended * CaloHitFrom = _allHits[ihitFrom];
	    float dist_in_generic = CaloHitTo->getGenericDistance()-CaloHitFrom->getGenericDistance();
	    if (dist_in_generic > r_min) {
		if (ifound ==1)
		    break;
		r_min +=  r_step; 
	    }

	    if (dist_in_generic > r_dist)
		break;
	    
	    float pos1[3];
	    float pos2[3];
	    for (int iposi=0; iposi<3; ++iposi) {
	      pos1[iposi] = 
		(float)CaloHitTo->getCalorimeterHit()->getPosition()[iposi];
	      pos2[iposi] = 
		(float)CaloHitFrom->getCalorimeterHit()->getPosition()[iposi];
	    }

	    float XDist = DistanceBetweenPoints(pos1,pos2);
	    float YRes = findResolutionParameter(CaloHitFrom, CaloHitTo);
	    if (YRes < 0.)
		std::cout << "Resolution parameter < 0" << std::endl; 


	    float YDist = 1 + _weightForReso*YRes + _weightForDist*XDist;

	    bool proxCriterion = YDist < YDistMin;

	    if (proxCriterion) {
	      YResMin = YRes;
	      YDistMin = YDist;
	      CaloHitTo->setCaloHitFrom(CaloHitFrom);
	      CaloHitTo->setYresFrom(YRes);
	    }

	    if (proxCriterion && YRes < YResCut) {
	      YResMin = YRes;
	      YDistMin = YDist;
	      CaloHitTo->setCaloHitFrom(CaloHitFrom);
	      CaloHitTo->setYresFrom(YRes);		
	    }


	    if (YRes < YResCut)
		ifound = 1;
	    
	    ihitFrom --;		
	}


	if (ifound == 1) { // Attach to already existing cluster
	    CaloHitExtended * calohit_AttachTo = CaloHitTo->getCaloHitFrom();
	    ClusterExtended * cluster = calohit_AttachTo->getClusterExtended();
	    CaloHitExtendedVec calohitvec = cluster->getCaloHitExtendedVec();
	    CaloHitTo->setClusterExtended(cluster);
	    cluster->addCaloHitExtended(CaloHitTo);
	    float distanceToHit = 0.0;
	    for (int ii=0;ii<3;++ii) {
	      float xx =  CaloHitTo->getCalorimeterHit()->getPosition()[ii]
		-  calohit_AttachTo->getCalorimeterHit()->getPosition()[ii];
	      distanceToHit += xx*xx;
	    }
	    distanceToHit = sqrt(distanceToHit);
	    CaloHitTo->setDistanceToNearestHit(distanceToHit);
	    float xDir[3];
	    bool redefineSP ;
	    float dif_in_dist = CaloHitTo->getGenericDistance() - calohit_AttachTo->getGenericDistance();	    
	    if (_typeOfGenericDistance == 0) {

		redefineSP = (int)calohitvec.size() < _NDefineSP;

	    }

	    else {
		
		redefineSP = dif_in_dist < _distanceToDefineDirection;

	    }

	    if (redefineSP ) {
		float xx = 0.;
		float yy = 0.;
		float zz = 0.;
		float ee = 0.;
		for (unsigned int i(0); i < calohitvec.size(); ++i) {
		    CaloHitExtended * chit = calohitvec[i];
		    float ene = chit->getCalorimeterHit()->getEnergy();
		    xx += chit->getCalorimeterHit()->getPosition()[0]*ene;
		    yy += chit->getCalorimeterHit()->getPosition()[1]*ene;
		    zz += chit->getCalorimeterHit()->getPosition()[2]*ene;
		    ee += ene;
		}		
		float xSP[3];
		xSP[0] = xx/ee;
		xSP[1] = yy/ee;
		xSP[2] = zz/ee;
		cluster->setStartingPoint(xSP);
	    }

	    float dist_to_SP(0.);
	    if (_typeOfGenericDistance == 0) {
	      //dist_to_SP = CaloHitTo->getDistanceToCalo();
		for (int i(0); i < 3; ++i) {
		    float xx = CaloHitTo->getCalorimeterHit()->getPosition()[i]-cluster->getStartingPoint()[i];
		    dist_to_SP += xx*xx;
		}
		dist_to_SP = sqrt(dist_to_SP);
	    }
	    else {
	      dist_to_SP = dif_in_dist;
	    }
		

	    if (dist_to_SP < _distanceToDefineDirection ) {
		for (int i(0); i < 3; ++i)
		    xDir[i] =  CaloHitTo->getCalorimeterHit()->getPosition()[i];		
	    }
	    else {
		for (int i(0); i < 3; ++i) 
		    xDir[i] = CaloHitTo->getCalorimeterHit()->getPosition()[i] - cluster->getStartingPoint()[i];
	    }

	      

	    CaloHitTo->setDirVec(xDir);
	}
	else { // Create new cluster
	    ClusterExtended * cluster = new ClusterExtended(CaloHitTo);
	    CaloHitTo->setClusterExtended(cluster);
	    CaloHitTo->setDistanceToNearestHit(0.0);
	    float xDir[3];
	    for (int i(0); i < 3; ++i)
		xDir[i] = CaloHitTo->getCalorimeterHit()->getPosition()[i];
	    CaloHitTo->setDirVec(xDir);	    
	    _allClusters.push_back( cluster );
	}


    }


}

void TrackwiseClustering::mergeForward() {

  int nClusters = (int)_allClusters.size();
  int nTotHits = (int)_allHits.size();

  int iCluster = 0;

  while (iCluster < nClusters) {
    ClusterExtended * clusterAR = _allClusters[iCluster];
    CaloHitExtendedVec hitvec = clusterAR->getCaloHitExtendedVec();
    int nHits = (int)hitvec.size();
    int iforw(0);
    int iback(0);
    CaloHitExtended * calohitAttachTo ;
    if (nHits > _nhit_minimal && nHits < _nhit_merge_forward) {
      //      std::cout << "attempt to merge forward" << std::endl;
      //      for (int i=1; i<nHits; ++i) {
      //      	float vec[3];
      //	for (int j=0;j<3;++j) 
      //	  vec[j] = hitvec[i]->getCalorimeterHit()->getPosition()[j]-
      //	    hitvec[0]->getCalorimeterHit()->getPosition()[j];
      //	hitvec[i]->setDirVec(vec);
      //      }
      int LowerBound = min(_nScanToMergeForward,nHits);
      for (int iCounterHit=0; iCounterHit < LowerBound; ++iCounterHit) {
	CaloHitExtended * calohit = hitvec[nHits-iCounterHit-1];
	int index = calohit->getIndex() + 1;
	int type = calohit->getType();
	float distance = 0.0;
	//	std::cout << " " << index << " " << nTotHits << std::endl; 
	while (distance < _distanceMergeForward[type] && index < nTotHits) {
	  CaloHitExtended * calohitTo = _allHits[index];	
	  distance = calohitTo->getGenericDistance() - calohit->getGenericDistance();
	  ClusterExtended * cluster_dummy = calohitTo->getClusterExtended();
	  CaloHitExtendedVec dummy = cluster_dummy->getCaloHitExtendedVec();
	  float yres = findResolutionParameter(calohit, calohitTo);
	  bool considerHit = yres < 2.0*_resolutionParameter[type]; 
	  //      int ndummy = (int)dummy.size();
	  //	  considerHit = considerHit && (ndummy > _nhit_minimal);
	  considerHit = considerHit && (cluster_dummy != clusterAR);
	  if (considerHit) {
	    calohitAttachTo = calohitTo;
	    calohit->setCaloHitTo(calohitTo);
	    calohit->setYresTo(yres);
	    iforw = 1;
	  }
	  if (iforw == 1)
	    break;
	  index++;
	}
	if (iforw == 1)
	  break;
      }
      if (iforw == 0) {
	for (int iCounterHit=0; iCounterHit < LowerBound; ++iCounterHit) {
	  CaloHitExtended * calohit = hitvec[nHits-iCounterHit-1];
	  int index = calohit->getIndex() - 1;
	  int type = calohit->getType();
	  float distance = 0.0;
	  //	std::cout << " " << index << " " << nTotHits << std::endl; 
	  while (distance < _stepTrackBack[type] && index >= 0) {
	    CaloHitExtended * calohitTo = _allHits[index];	
	    distance = - calohitTo->getGenericDistance() + calohit->getGenericDistance();
	    ClusterExtended * cluster_dummy = calohitTo->getClusterExtended();
	    CaloHitExtendedVec dummy = cluster_dummy->getCaloHitExtendedVec();
	    float yres = findResolutionParameter(calohitTo, calohit);
	    bool considerHit = yres < 2.0*_resolutionParameter[type]; 
	    //      int ndummy = (int)dummy.size();
	    //	  considerHit = considerHit && (ndummy > _nhit_minimal);
	    considerHit = considerHit && (cluster_dummy != clusterAR);
	    if (considerHit) {
	      calohitAttachTo = calohitTo;
	      calohit->setCaloHitTo(calohitTo);
	      calohit->setYresTo(yres);
	      iback = 1;
	    }
	    if (iback == 1)
	      break;
	    index--;
	  }
	  if (iback == 1)
	    break;
	}
      }
      if (iforw == 1) {
	//	std::cout << "Merging forward " << std::endl;
	ClusterExtended * clusterTo = calohitAttachTo->getClusterExtended();
	CaloHitExtendedVec hitvecAttached = clusterTo->getCaloHitExtendedVec();
	int nHitsAttached = (int)hitvecAttached.size();
	for (int jHit = 0; jHit < nHitsAttached; ++jHit) {
	  CaloHitExtended * hitToAttach = hitvecAttached[jHit];
	  clusterAR->addCaloHitExtended( hitToAttach );
	  hitToAttach->setClusterExtended( clusterAR );
	}
      clusterTo->Clear();
      }
      if (iback == 1) {
	//     std::cout << "Merging backward " << std::endl;
	ClusterExtended * clusterTo = calohitAttachTo->getClusterExtended();
	CaloHitExtendedVec hitvecAttached = clusterTo->getCaloHitExtendedVec();
	int nHitsAttached = (int)hitvecAttached.size();
	clusterAR->Clear();
	for (int jHit = 0; jHit < nHitsAttached; ++jHit) {
	  CaloHitExtended * hitToAttach = hitvecAttached[jHit];
	  clusterAR->addCaloHitExtended( hitToAttach );
	  hitToAttach->setClusterExtended( clusterAR );
	}
	for (int jHit =0; jHit < nHits; ++jHit) {
	  CaloHitExtended * hitToAttach = hitvec[jHit];
	  clusterAR->addCaloHitExtended( hitToAttach );
	  hitToAttach->setClusterExtended( clusterAR );
	}
      clusterTo->Clear();
      }
      if (iforw == 0 && iback == 0) {
	iCluster++;
      }
    }
    else {
      iCluster++;
    }
  }



}
void TrackwiseClustering::mergeLowMultiplicity() {

  int nClusters = (int)_allClusters.size();
  int nTotHits = (int)_allHits.size();
  int iCluster = 0;

  while (iCluster < nClusters) {
    ClusterExtended * clusterAR = _allClusters[iCluster];
    CaloHitExtendedVec hitvec = clusterAR->getCaloHitExtendedVec();
    int nHits = (int)hitvec.size(); 
    
    if (nHits > 0 && nHits < _nhit_minimal) {
      //      std::cout << "attempt to merge forward" << std::endl;
      //      for (int i=1; i<nHits; ++i) {
      //	float vec[3];
      //      	for (int j=0;j<3;++j) 
      //      	  vec[j] = hitvec[i]->getCalorimeterHit()->getPosition()[j]-
      //          hitvec[0]->getCalorimeterHit()->getPosition()[j];
      //      	hitvec[i]->setDirVec(vec);
      //      }
      //  std::cout << "attempt to merge low multiplicity cluster " << iCluster << std::endl;
      CaloHitExtended * first = hitvec[0];
      CaloHitExtended * last = hitvec[nHits - 1];
      int index = last->getIndex() + 1;
      float yres_min = 1.0e+20;
      CaloHitExtended * hitToAttach;
      if (first->getCaloHitFrom() != NULL) {
	yres_min = first->getYresFrom();
	hitToAttach = first->getCaloHitFrom();
	ClusterExtended * cluster_dummy =
	  hitToAttach->getClusterExtended();
	if (cluster_dummy == clusterAR) {
	  yres_min = 1.0e+20;
	}
      }
      else {
	yres_min = 1.0e+20;
      } 
      float distance = 0.0;
      int type = last->getType();
      int ifound = 0;
      while (distance < _distanceMergeForward[type] && index < nTotHits) {
	CaloHitExtended * hitTo = _allHits[index];
	distance = hitTo->getGenericDistance() - last->getGenericDistance();
	float yres =  findResolutionParameter(last, hitTo);
	ClusterExtended * cluster_dummy = 
	  hitTo->getClusterExtended();
	if (yres < yres_min && clusterAR != cluster_dummy ){
	  hitToAttach = hitTo;
	  yres_min = yres;
	  ifound = 1;
	}
	++index;
      }
      
      if (yres_min < _resolutionToMerge) {
	//	std::cout << "merging is done : cluster " << iCluster 
	//                << " is merged " << std::endl;
	if (ifound == 0) {
	  ClusterExtended * clusterTo = hitToAttach->getClusterExtended();	  
	  for (int jHit = 0; jHit < nHits; ++jHit){
	    CaloHitExtended * hh = hitvec[jHit];
	    clusterTo->addCaloHitExtended(hh);
	    hh->setClusterExtended(clusterTo);
	  }
	  clusterAR->Clear();
	}
	else {
	  ClusterExtended * clusterTo = hitToAttach->getClusterExtended();
	  CaloHitExtendedVec dummyVec = clusterTo->getCaloHitExtendedVec();
	  int nDummy = (int)dummyVec.size();
	  for (int jHit = 0; jHit < nDummy; ++jHit) {
	    CaloHitExtended * hh = dummyVec[jHit];
	    clusterAR->addCaloHitExtended( hh );
	    hh->setClusterExtended( clusterAR );
	  }
	  clusterTo->Clear();
	  if (nDummy == 0) // Just protection
	    iCluster++;
	}
      }
      else {	
	iCluster++;	
      }	
    }
    else {
      iCluster++;
    }
  }


}

void TrackwiseClustering::MergeTrackSegments() {

  int nClusters = int(_allClusters.size());
  int iCluster = 0;

  while (iCluster < nClusters) {
    ClusterExtended * Cluster = _allClusters[iCluster];
    float chi2ph = Cluster->getHelixChi2R();
    float chi2z  = Cluster->getHelixChi2Z();
    CaloHitExtendedVec hitVec = Cluster->getCaloHitExtendedVec();
    int nHits = int(hitVec.size());
    if (chi2ph < 3.0 && chi2z < 10.0 && nHits > 4) {
      int ifound = 0;
      int iforw = 0;
      float xEnd[3];
      float xBeg[3];
      ClusterExtended * ClusterAttach; 
      for (int j=0; j<3; ++j) {
	xBeg[j] = Cluster->getLowEdge()[j];
	xEnd[j] = Cluster->getUpEdge()[j];
      }
      for (int i=0; i<nClusters; ++i) {
	ClusterExtended * ClusterTo = _allClusters[i];
	CaloHitExtendedVec hitVecTo = ClusterTo->getCaloHitExtendedVec();
	int nHitsTo = int(hitVecTo.size());
	if (Cluster != ClusterTo && nHitsTo > 0) {
	  //	  std::cout << "Condition 1 " << std::endl;
	  float posCluster[3];
	  float posClusterTo[3];
	  for (int j=0;j<3;++j) {
	    posCluster[j] = Cluster->getPosition()[j];
	    posClusterTo[j] = ClusterTo->getPosition()[j];
	  }
	  float dist = DistanceBetweenPoints(posCluster,posClusterTo);
	  if (dist < 500.) {
	    //	    std::cout << "Condition 2 " << std::endl;
	    for (int ihit=0;ihit<nHitsTo;++ihit) {
	      float hpos[3];
	      for (int j=0;j<3;++j) 
		hpos[j] = hitVecTo[ihit]->getCalorimeterHit()->getPosition()[j];
	      
	      float distEnd = DistanceBetweenPoints(xEnd,hpos);
	      float distBeg = DistanceBetweenPoints(xBeg,hpos);
	      float distMin = distEnd ;
	      iforw = 1;
	      if (distBeg < distMin) {
		iforw = 0;
		distMin = distBeg;
	      }
	      float dd[3];
	      HelixClass helix = Cluster->getHelix();
	      float time = helix.getDistanceToPoint(hpos,dd);
	      float dist1 = dd[2];
	      float phi0 = helix.getPhi0();
	      float d0 = helix.getD0();
	      float z0 = helix.getZ0();
	      float tanlambda = helix.getTanLambda();
	      float omega = helix.getOmega();
	      float B = 4.0;
	      helix.Initialize_Canonical(phi0, -d0, z0, omega, 
			      tanlambda, B);
	      for (int j=0;j<3;++j) 
		hpos[j] = hitVecTo[ihit]->getCalorimeterHit()->getPosition()[j];

	      time = helix.getDistanceToPoint(hpos,dd);
	      float dist2 = dd[2];
	      float dist = fmin(dist1,dist2);
	      //	      std::cout << "Check distances" << dist1 << " " << dist2 << " " << dist << std::endl;
	      if (dist < 50.0 && distMin) {
		ClusterAttach = ClusterTo;
		ifound = 1;
		break;
	      }
	      
	      

	    }
	  }
	}
	if (ifound == 1)
	  break;
      }
      if (ifound == 1) {
	//	std::cout << "Merging Track Segment " << std::endl;
	CaloHitExtendedVec hitVecTo = ClusterAttach->getCaloHitExtendedVec();
	int nHitsTo = int(hitVecTo.size());
	if (iforw == 1)
	  for (int jHit = 0; jHit < nHitsTo; ++jHit) {
	    CaloHitExtended * hh = hitVecTo[jHit];
	    Cluster->addCaloHitExtended( hh );
	    hh->setClusterExtended( Cluster );
	  }
	else {
	  Cluster->Clear();
	  for (int jHit = 0; jHit < nHitsTo; ++jHit) {
	    CaloHitExtended * hh = hitVecTo[jHit];
	    Cluster->addCaloHitExtended( hh );
	    hh->setClusterExtended( Cluster );
	  }
	  for (int jHit = 0; jHit < nHits; ++jHit) {
	    CaloHitExtended * hh = hitVec[jHit];
	    Cluster->addCaloHitExtended( hh );
	    hh->setClusterExtended( Cluster );
	  }	  
	}
	ClusterAttach->Clear();
	if (nHitsTo == 0) // Just protection
	  iCluster++;
	calculateProperties(Cluster);
	//debug
	/*
	std::cout << "After merging " << Cluster->getHelixChi2R() << " " 
		  << Cluster->getHelixChi2Z() << std::endl;
	std::cout << "              " << ClusterAttach->getHelixChi2R() << " " 
		  << ClusterAttach->getHelixChi2Z() << std::endl;
	*/
      }
      else {
	iCluster++;
      }
    }
    else {
      iCluster++;
    }
  }

}
