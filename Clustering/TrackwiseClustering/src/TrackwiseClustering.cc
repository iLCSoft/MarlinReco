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

using namespace lcio ;
using namespace marlin ;
using namespace std;


TrackwiseClustering aTrackwiseClustering;

TrackwiseClustering::TrackwiseClustering() : Processor("TrackwiseClustering") {

    registerProcessorParameter( "DistanceForDirection", 
				"Distance to Define Direction", 
				_distanceToDefineDirection,
				(float)25.);

     registerProcessorParameter( "DistanceToTrackSeed", 
				"Distance to Track Seed", 
				_distanceToTrackSeed,
				(float)25.);    


    std::vector<float>  distanceTrackBack;
    distanceTrackBack.push_back(100.);
    distanceTrackBack.push_back(500.);
    
    registerProcessorParameter( "DistanceTrackBack" , 
				"Distance to Track Back "  ,
				_distanceTrackBack,
				 distanceTrackBack); 

    std::vector<float>  stepTrackBack;
    stepTrackBack.push_back(10.0);
    stepTrackBack.push_back(100.0);
    
    registerProcessorParameter( "StepTrackBack" , 
				"Step to Track Back "  ,
				_stepTrackBack,
				 stepTrackBack); 

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
    
    registerProcessorParameter( "EcalCollections" , 
				"Ecal Collection Names "  ,
				_ecalCollections,
				 EcalCollections);

    std::vector<std::string> HcalCollections;
    HcalCollections.push_back(std::string("HCAL"));
    
    registerProcessorParameter( "HcalCollections" , 
				"Hcal Collection Names "  ,
				_hcalCollections,
				 HcalCollections);

    std::vector<std::string> TrackCollections;
    TrackCollections.push_back(std::string("Track"));
    
    registerProcessorParameter( "TrackCollections" , 
				"Track Collection Names "  ,
				_trackCollections,
				 TrackCollections);

    
    registerProcessorParameter( "ClusterCollection" , 
				"Cluster Collection Name "  ,
				_clusterCollection,
				 std::string("ClustersAR"));

    
    registerProcessorParameter( "ZOfEndcap" ,
                                "Z coordinate of Endcap" ,
				_zofendcap,
				(float)2820.);
    
    registerProcessorParameter( "ROfBarrel" , 
                                "Radius of Barrel" , 
				_rofbarrel,
				(float)1700.);
				
    registerProcessorParameter( "GlobalPhi" , 
				"Global Phi angle" ,
				_phiofbarrel,
				(float)0.);

    registerProcessorParameter( "NFoldSymmetry" ,
				"Global N-fold Symmetry" , 
				_nsymmetry, 
				8);
    
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

   registerProcessorParameter( "DoMerging" , 
			       "Do Merging", 
			       _doMerging,
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



}

void TrackwiseClustering::init() {

    _const_pi    = acos(-1.);
    _const_2pi = 2.0*_const_pi;
    _const_pi_n  = _const_pi/float(_nsymmetry);
    _const_2pi_n = 2.0*_const_pi/float(_nsymmetry);

    _thetaofendcap = (float)atan((double)(_rofbarrel/_zofendcap));

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
    if (_doMerging == 1) {
      mergeForward();
      mergeLowMultiplicity();
    }
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
		    float dist = CalculateGenericDistance(calohit, _typeOfGenericDistance);
		    calohit->setGenericDistance(dist);
		    if (dist < _xmin_in_distance) 
			_xmin_in_distance = dist;
		    if (dist > _xmax_in_distance) 
			_xmax_in_distance = dist;		
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
		    float dist = CalculateGenericDistance(calohit, _typeOfGenericDistance);
		    calohit->setGenericDistance(dist);
		    if (dist < _xmin_in_distance) 
			_xmin_in_distance = dist;
		    if (dist > _xmax_in_distance) 
			_xmax_in_distance = dist;
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
      product = 0.999*fabs(product)/product;
    }

    float angle = acos(product);
    
    return xdist*angle;
 
}

float TrackwiseClustering::CalculateGenericDistance(CaloHitExtended *calohit, int itype) {
    float xDistance =0.0;
    float rDistance =0.0;
    
    for (int i(0); i < 3; ++i) {
	float x = calohit->getCalorimeterHit()->getPosition()[i];
	rDistance += x*x; 	
    }
    rDistance = sqrt(rDistance);
    
    if (itype == 0) {
	xDistance = rDistance;
    }
    else {
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
    }
    return xDistance;
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
	    float xDir[3];
	    bool redefineSP ;
	    float dif_in_dist = CaloHitTo->getGenericDistance() - calohit_AttachTo->getGenericDistance();	    
	    if (_typeOfGenericDistance == 0) {

		redefineSP = (int)calohitvec.size() < _NDefineSP;

	    }

	    else {
		
		redefineSP = dif_in_dist < _distanceToDefineDirection;

	    }

	    if (redefineSP) {
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

  for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
    ClusterExtended * clusterAR = _allClusters[iCluster];
    CaloHitExtendedVec hitvec = clusterAR->getCaloHitExtendedVec();
    int nHits = (int)hitvec.size();
    int ifound(0);
    CaloHitExtended * calohitAttachTo ;
    if (nHits > _nhit_minimal && nHits < _nhit_merge_forward) {
      int LowerBound = min(_nScanToMergeForward,nHits);
      for (int iCounterHit=0; iCounterHit < LowerBound; ++iCounterHit) {
	CaloHitExtended * calohit = hitvec[nHits-iCounterHit-1];
	int index = calohit->getIndex() + 1;
	int type = calohit->getType();
	float distance = 0.0;
	while (distance < _distanceMergeForward[type] && index < nTotHits) {
	  CaloHitExtended * calohitTo = _allHits[index];	
	  distance = calohitTo->getGenericDistance() - calohit->getGenericDistance();
	  ClusterExtended * cluster_dummy = calohitTo->getClusterExtended();
	  CaloHitExtendedVec dummy = cluster_dummy->getCaloHitExtendedVec();
	  int ndummy = (int)dummy.size();
	  float yres = findResolutionParameter(calohit, calohitTo);
	  bool considerHit = yres < _resolutionParameter[type]; 
	  considerHit = considerHit && (ndummy > _nhit_minimal);
	  considerHit = considerHit && (cluster_dummy != clusterAR);
	  if (considerHit) {
	    calohitAttachTo = calohitTo;
	    calohit->setCaloHitTo(calohitTo);
	    calohit->setYresTo(yres);
	    ifound = 1;
	  }
	  if (ifound == 1)
	    break;
	  ++index;
	}
	if (ifound == 1)
	  break;
      }
      if (ifound == 1) {
	std::cout << "Merging forward " << std::endl;
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
      CaloHitExtended * first = hitvec[0];
      CaloHitExtended * last = hitvec[nHits - 1];
      int index = last->getIndex() + 1;
      float yres_min;
      CaloHitExtended * hitToAttach;
      if (first->getCaloHitFrom() != NULL) {
	yres_min = first->getYresFrom();
	hitToAttach = first->getCaloHitFrom();
	ClusterExtended * cluster_dummy =
	  hitToAttach->getClusterExtended();
	if (cluster_dummy == clusterAR) {
	  yres_min = 1.0e+20;
	  std::cout << "already existing cluster" << std::endl;
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
	if (yres < yres_min && clusterAR != cluster_dummy) {
	  hitToAttach = hitTo;
	  yres_min = yres;
	  ifound = 1;
	}
	++index;
      }

      if (yres_min < _resolutionToMerge) {
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
	  if (nDummy == 0)
	    ++iCluster;
	}
      }
      else {
	
	++iCluster;
	
      }	
    }
    else {
      ++iCluster;
    }

  }


}
