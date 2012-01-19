#include "EMShowerFinder.h"

using namespace lcio ;
using namespace marlin ;

EMShowerFinder aEMShowerFinder ;



EMShowerFinder::EMShowerFinder() : Processor("EMShowerFinder") {

  // modify processor description
  _description = "a photon finder processor based on the KIT and KITutil 'classes and fuctions'" ;

  registerProcessorParameter( "colNameECAL" , 
			      "ECAL Collection Name"  ,
			      _colNameECAL,
			      std::string("ECAL") );

  registerProcessorParameter( "collectionNameOfEMShowerCandidates",
			      "Name of the collection of EM shower candidates",
			      _collectionNameOfEMShowerCandidates,
			      std::string("EMShowerCandidates"));

  registerProcessorParameter( "Cleaning",
			      "To do the cleaning on hits or not ",
			      _ToClean,
			      std::string("YES"));
  registerProcessorParameter( "TopologicalCut",
			      "At which number of neighbors to put the threshold, condition is < so you need to put N+1 ",
			       _CleanCut,
			      (int)5);

  registerProcessorParameter( "NumberOfLevels",
			      "Number of levels for central loop ",
			       _N,
			      (int)10);
  vector<float> miipstep;
  miipstep.push_back(0.1);
  miipstep.push_back(1.5);
  miipstep.push_back(2.5);
  miipstep.push_back(4.0);
  miipstep.push_back(6.0);
  miipstep.push_back(9.0);
  miipstep.push_back(16.0);
  miipstep.push_back(26.0);
  miipstep.push_back(41.0);
  miipstep.push_back(65.0);

  registerProcessorParameter( "Levels",
			     "Levels for central loop in MIP ",
			      _miipstep,
			       miipstep);

  registerProcessorParameter( "MinHit0",
			     "Minimal Number of hits for ground level cluster ",
			      _MinHit0,
			      (int)4);

  registerProcessorParameter( "MinHitSplit",
			      "Minimal Number of hits for i-th level cluster ",
			      _MinHitSplit,
			      (int)2);

  registerProcessorParameter( "Rcut",
			      "Fluctuation suprresion cut",
			      _Rcut,
			      (double)0.4);
  registerProcessorParameter( "Distcut",
			      "Square of distance cut for merging ",
			      _Distcut,
			      (double)35.0);
  registerProcessorParameter( "Coscut",
			      "Cosine of the angle for merging ",
			      _Coscut,
			      (double)0.95);

  registerProcessorParameter( "energyDeviationCut",
			      "cut on energy deviation of em shower candidates from estimated energy ( abs( (Ecluster - Eestimated)/Eestimated ) < cut )",
			      _energyDeviationCut,
			      (double)0.25);

  registerProcessorParameter( "probabilityDensityCut",
			      "cut on the probability density to assign hits to shower cores",
			      _probabilityDensityCut,
			      (double)1E-3);

  registerProcessorParameter( "DebugLevel",
			      "limits the amount of information written to std out (0 - none, 9 - maximal information)",
			      _debugLevel,
			      int(0) );

  registerProcessorParameter( "DrawOnCED",
			      "draw objects on CED",
			      _drawOnCED,
			      int(0) );


}


void EMShowerFinder::init() {
  
  // usually a good idea to 
  printParameters();

  // FIXME: hard coded cell id's for old Mokka (e.g. Mokka v5.4) versions)
  CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");


  // debug
  if (_drawOnCED)  MarlinCED::init(this);

  _nRun = 0 ;
  _nEvt = 0 ;


}


void EMShowerFinder::processRunHeader( LCRunHeader* run) {

  ++_nRun;

}


void EMShowerFinder::processEvent( LCEvent * evt ) {


  // debug
  if (_drawOnCED) {

    MarlinCED::newEvent(this,0);
    MarlinCED::drawMCParticleTree(evt,"MCParticle",0.05,4.0,15.5,50.0,1626.0,2500.0);

  }


  try {
    
    LCCollection* colt = evt->getCollection(_colNameECAL.c_str()) ;
    CellIDDecoder<CalorimeterHit>CDECAL =CellIDDecoder<CalorimeterHit>(colt);

    LCCollectionVec* clscol = new LCCollectionVec(LCIO::CLUSTER);

    if( colt!=0 ) {

      unsigned int nelem=colt->getNumberOfElements(); 
      // MAIN CONTAINER OF SHITS 
      vector<Superhit2*> calo[10]; 
      
      // creating all superhits 
      CreateAllShits2(colt,CDECAL,calo); 
      
      // precalculation
      TotalPrecalc2(calo,CDECAL,nelem,_CleanCut);  
      
      // setting the parameters of the alghorithm
      vector <PROTSEED2> prs2;
      CoreCut2 Ccut;
      Ccut.Rcut=_Rcut;
      Ccut.Distcut=_Distcut;
      Ccut.Coscut=_Coscut;
      Ccut.MinHit0=(unsigned int) _MinHit0;
      Ccut.MinHitSplit=(unsigned int) _MinHitSplit;
      const unsigned int N=_N;
      vector<Tmpcl2*> bbb[MAXARRAYSIZE];
	 
      // finding cores.
      if( _ToClean=="YES" || _ToClean=="yes")
	{
	  // cout << " da koristim cat " << endl;
	  FindCores2(&(calo[4]), bbb , &prs2,_N,_miipstep,Ccut);     
	}else{
	FindCores2(&(calo[0]), bbb , &prs2,_N,_miipstep,Ccut);   
      }


      // container to store information to assign hits to photon candidates
      std::vector<ECALHitWithAttributes> ECALHitsWithAttributes;

      std::vector<CoreCalib2> coreCalibrationLDC00;
      CreateCalibrationLDC00(&coreCalibrationLDC00);

      
      for(unsigned int i=0;i<prs2.size();i++)
	{
	  if(prs2[i].active==true)
	    {

	      double centerPosition[3];
	      centerPosition[0]=(double)prs2[i].cl->getCenter()[0];
	      centerPosition[1]=(double)prs2[i].cl->getCenter()[1];
	      centerPosition[2]=(double)prs2[i].cl->getCenter()[2];

	      double directionOfCluster[3];
	      prs2[i].cl->findInertia();
	      directionOfCluster[0] = prs2[i].cl->direction[0]; // already normalised
	      directionOfCluster[1] = prs2[i].cl->direction[1]; // already normalised
	      directionOfCluster[2] = prs2[i].cl->direction[2]; // already normalised

	      
	      double coreEnergy = 0.0;
	      double startPosition[3];
	      double startPositionMinProjection[3];
	      double startPositionMaxProjection[3];
	      double minProjectionAlongMainPrincipalAxis = DBL_MAX;
	      double maxProjectionAlongMainPrincipalAxis = 0.0;


	      for(unsigned int j=0;j<prs2[i].cl->hits.size();j++) {

		// sum up core energy
		coreEnergy += prs2[i].cl->hits[j]->chit->getEnergy();

		// debug 
		if (_drawOnCED) ced_hit( prs2[i].cl->hits[j]->chit->getPosition()[0], prs2[i].cl->hits[j]->chit->getPosition()[1], 
					 prs2[i].cl->hits[j]->chit->getPosition()[2], 0 | 4 << CED_LAYER_SHIFT, 2, 0xff47e6);
		


		double vCenterToHit[3];		
		for(unsigned int k=0;k<3;k++) vCenterToHit[k] = centerPosition[k] - prs2[i].cl->hits[j]->chit->getPosition()[k];
		
		double projectionAlongMainPrincipalAxis = 0.0;
		for(unsigned int k=0;k<3;k++) projectionAlongMainPrincipalAxis += directionOfCluster[k]*vCenterToHit[k];

		if ( projectionAlongMainPrincipalAxis < minProjectionAlongMainPrincipalAxis ) {
		  
		  minProjectionAlongMainPrincipalAxis = projectionAlongMainPrincipalAxis;
		  for(unsigned int k=0;k<3;k++) startPositionMinProjection[k] = centerPosition[k] - 2.0*minProjectionAlongMainPrincipalAxis*directionOfCluster[k];
		  //for(unsigned int k=0;k<3;k++) startPositionMinProjection[k] = prs2[i].cl->hits[j]->chit->getPosition()[k];

		}
		
		if ( projectionAlongMainPrincipalAxis > maxProjectionAlongMainPrincipalAxis ) {
		
		  maxProjectionAlongMainPrincipalAxis = projectionAlongMainPrincipalAxis;
		  for(unsigned int k=0;k<3;k++) startPositionMaxProjection[k] = centerPosition[k] - 2.0*maxProjectionAlongMainPrincipalAxis*directionOfCluster[k];
		  //for(unsigned int k=0;k<3;k++) startPositionMaxProjection[k] = prs2[i].cl->hits[j]->chit->getPosition()[k];

		}


		double rStartPositionMinProjection = sqrt( pow(startPositionMinProjection[0],2) + pow(startPositionMinProjection[1],2) + pow(startPositionMinProjection[2],2));
		double rStartPositionMaxProjection = sqrt( pow(startPositionMaxProjection[0],2) + pow(startPositionMaxProjection[1],2) + pow(startPositionMaxProjection[2],2));

		if ( rStartPositionMinProjection < rStartPositionMaxProjection ) {
		  
		  startPosition[0] = startPositionMinProjection[0];
		  startPosition[1] = startPositionMinProjection[1];
		  startPosition[2] = startPositionMinProjection[2];
		  
		}
		else {

		  startPosition[0] = startPositionMaxProjection[0];
		  startPosition[1] = startPositionMaxProjection[1];
		  startPosition[2] = startPositionMaxProjection[2];
		  
		}



	      }


	      double photonE = giveMeEEstimate2(prs2[i].level,prs2[i].cl->getEnergy(),coreCalibrationLDC00);





	      // debug
	      if ( _debugLevel > 5 ) {
		std::cout << "EM Shower Core found by the KIT package: " << std::endl
			  << "CoreEnergy: " << coreEnergy << "  " << " estimated 'em' particle energy: " <<  photonE << std::endl 
			  << "centerPos: " << "(" << centerPosition[0] << ", " << centerPosition[1] << ", " << centerPosition[2] << ")" 
			  << "  "	<< "startPos: " << "(" << startPosition[0] << ", " << startPosition[1] << ", " << startPosition[2] << ")" 
			  << std::endl << std::endl;
	      }

	      // debug 
	      if (_drawOnCED) {
		ced_hit(centerPosition[0],centerPosition[1],centerPosition[2], 2 | 4 << CED_LAYER_SHIFT, 6, 0xffffff);
		ced_hit(startPosition[0],startPosition[1],startPosition[2], 2 | 4 << CED_LAYER_SHIFT, 6, 0xff0004);
	      }


	      
	      
	      
	      if ( photonE > 0.0 ) {

	      
		Photon2* photonFinder = new Photon2(photonE,centerPosition,startPosition);

		// loop over all ECAL hits
		for(unsigned int j = 0; j < nelem; ++j) {

		  double probabilityAndDistance[2];
		  CalorimeterHit* ECALHit = dynamic_cast<CalorimeterHit*>(colt->getElementAt(j));
		  
		  photonFinder->Prob(ECALHit,_probabilityDensityCut,probabilityAndDistance);
		  
		  bool isECALHitAlreadyAssigned = false;
		  int indexOfECALHitAlreadyAssigned = 0;
		  for(unsigned int k = 0; k < ECALHitsWithAttributes.size(); ++k) {
		    
		    if ( ECALHitsWithAttributes.at(k).ECALHit == ECALHit) {
		      
		      isECALHitAlreadyAssigned = true;
		      indexOfECALHitAlreadyAssigned = k;
		      break;
		      
		    }
		    
		  }


		
		  if ( !isECALHitAlreadyAssigned ) {
		    
		    if ( probabilityAndDistance[0] > 0.0 ) {
		      
		      ECALHitWithAttributes hitWithAttributes;
		      
		      hitWithAttributes.ECALHit = ECALHit;
		      hitWithAttributes.relatedCores.push_back(&(prs2[i]));
		      hitWithAttributes.probabilitiesForThisECALHit.push_back(probabilityAndDistance[0]);
		      hitWithAttributes.distancesToCoresForThisECALHit.push_back(probabilityAndDistance[1]);
		      hitWithAttributes.estimatedEnergyPerCore.push_back(photonE);
		      
		      ECALHitsWithAttributes.push_back(hitWithAttributes);
		      
		    }
		    
		  }
		  else {
		    
		    if ( probabilityAndDistance[0] > 0.0 ) {
		    
		      ECALHitsWithAttributes.at(indexOfECALHitAlreadyAssigned).relatedCores.push_back(&(prs2[i]));
		      ECALHitsWithAttributes.at(indexOfECALHitAlreadyAssigned).probabilitiesForThisECALHit.push_back(probabilityAndDistance[0]);
		      ECALHitsWithAttributes.at(indexOfECALHitAlreadyAssigned).distancesToCoresForThisECALHit.push_back(probabilityAndDistance[1]);
		      ECALHitsWithAttributes.at(indexOfECALHitAlreadyAssigned).estimatedEnergyPerCore.push_back(photonE);
		      
		    }
		    
		  }
		  
		}
		
		delete photonFinder;
		photonFinder = 0;
		
	      }
	    }
	}
      
	

      
      // debug
      
      if ( _debugLevel > 5 ) {
	std::cout << "CalorimeterHits in ECAL and their probability to belong to EM showers: " << std::endl;
	for(unsigned int i = 0; i  <  ECALHitsWithAttributes.size(); ++i) {	
	  std::cout << "ECAL Hit: " << ECALHitsWithAttributes.at(i).ECALHit << " (" << ECALHitsWithAttributes.at(i).ECALHit->getPosition()[0] << ", " 
		    << ECALHitsWithAttributes.at(i).ECALHit->getPosition()[1] << ", " << ECALHitsWithAttributes.at(i).ECALHit->getPosition()[2] << ")" << "  " 
		    <<"E = " << ECALHitsWithAttributes.at(i).ECALHit->getEnergy() << "  " << std::endl;
	  std::cout << "Cores: " << "(";
	  for(unsigned int j = 0; j  <  ECALHitsWithAttributes.at(i).relatedCores.size(); ++j) { 
	    std::cout << ECALHitsWithAttributes.at(i).relatedCores.at(j);
	    if (j  <  ECALHitsWithAttributes.at(i).relatedCores.size() - 1) std::cout << ",";
	  }
	  std::cout << ")" << "  " << std::endl;
	  std::cout << "Probs: " << "(";
	  for(unsigned int j = 0; j  <  ECALHitsWithAttributes.at(i).probabilitiesForThisECALHit.size(); ++j) { 
	    std::cout << ECALHitsWithAttributes.at(i).probabilitiesForThisECALHit.at(j);
	    if (j  <  ECALHitsWithAttributes.at(i).probabilitiesForThisECALHit.size() - 1) std::cout << ",";
	  }
	  std::cout << ")" << "  " << std::endl;
	  std::cout << "Dists: " << "(";
	  for(unsigned int j = 0; j  <  ECALHitsWithAttributes.at(i).distancesToCoresForThisECALHit.size(); ++j) { 
	    std::cout << ECALHitsWithAttributes.at(i).distancesToCoresForThisECALHit.at(j);
	    if (j  <  ECALHitsWithAttributes.at(i).distancesToCoresForThisECALHit.size() - 1) std::cout << ",";
	  }
	  std::cout << ")" << "  " << std::endl;
	  std::cout << "Estimated EM Energy: " << "(";
	  for(unsigned int j = 0; j  <  ECALHitsWithAttributes.at(i).estimatedEnergyPerCore.size(); ++j) { 
	    std::cout << ECALHitsWithAttributes.at(i).estimatedEnergyPerCore.at(j);
	    if (j  <  ECALHitsWithAttributes.at(i).estimatedEnergyPerCore.size() - 1) std::cout << ",";
	  }
	  std::cout << ")" << std::endl << std::endl;
	}
      }
      
	  
      
      
      
      std::map<PROTSEED2*,ClusterImpl*> mapCoreToCluster;
      std::map<PROTSEED2*,double> mapCoreToEstimatedEnergy;
      
      for(unsigned int i = 0; i < ECALHitsWithAttributes.size(); ++i) {
		
	PROTSEED2* coreToAddHitTo;
	double estimatedEnergyOfCoreToAddHitTo = 0.0;
	
	if (ECALHitsWithAttributes.at(i).relatedCores.size() == 1) {
	  
	  coreToAddHitTo = ECALHitsWithAttributes.at(i).relatedCores.at(0);	
	  estimatedEnergyOfCoreToAddHitTo = ECALHitsWithAttributes.at(i).estimatedEnergyPerCore.at(0);
	  
	}
	else {
	  
	  unsigned int indexOfMaxProbability = 0;
	  double maxProbability = 0.0;
	  
	  for(unsigned int j = 0; j  <  ECALHitsWithAttributes.at(i).relatedCores.size(); ++j){
	    
	    if ( ECALHitsWithAttributes.at(i).probabilitiesForThisECALHit.at(j) > maxProbability ) {
	      
	      indexOfMaxProbability = j;
	      maxProbability = ECALHitsWithAttributes.at(i).probabilitiesForThisECALHit.at(j);
	      
	    }
	    
	  }
	  
	  coreToAddHitTo = ECALHitsWithAttributes.at(i).relatedCores.at(indexOfMaxProbability);
	  estimatedEnergyOfCoreToAddHitTo = ECALHitsWithAttributes.at(i).estimatedEnergyPerCore.at(indexOfMaxProbability);
	  
	}

	    
	std::map<PROTSEED2*,ClusterImpl*>::iterator position = mapCoreToCluster.find(coreToAddHitTo);
	    
	if ( position ==  mapCoreToCluster.end() ) { 
	      
	  mapCoreToCluster[coreToAddHitTo] = new ClusterImpl();
	  mapCoreToEstimatedEnergy[coreToAddHitTo] = estimatedEnergyOfCoreToAddHitTo;
	  position = mapCoreToCluster.find(coreToAddHitTo);
	      
	}	
	    
	(*position).second->addHit(ECALHitsWithAttributes.at(i).ECALHit,(float)1.0);
	    
      }
	  

      int iCluster = 0;

      for(std::map<PROTSEED2*,ClusterImpl*>::iterator i = mapCoreToCluster.begin(); i != mapCoreToCluster.end(); ++i) {
       	    
	ClusterImpl* cluster = (*i).second;
	    
	double estimatedEnergyOfCluster = mapCoreToEstimatedEnergy[(*i).first];
	    
	// set energy and position of cluster
	unsigned int n = cluster->getCalorimeterHits().size();
	    
	float* x = new float[n];
	float* y = new float[n];
	float* z = new float[n];
	float* a = new float[n];
	float energy = 0.0;
	    
	for (unsigned int j = 0; j < n; ++j) {
	      
	  x[j] = cluster->getCalorimeterHits().at(j)->getPosition()[0];
	  y[j] = cluster->getCalorimeterHits().at(j)->getPosition()[1];
	  z[j] = cluster->getCalorimeterHits().at(j)->getPosition()[2];
	  a[j] = cluster->getCalorimeterHits().at(j)->getEnergy();
	      
	  energy += cluster->getCalorimeterHits().at(j)->getEnergy();
	      
	}
	    
	    
	// cut of EM candidtate clusters which differ more than _energyDeviationCut from enery estimate
	double energyDeviation = fabs( (energy-estimatedEnergyOfCluster)/estimatedEnergyOfCluster );
	
	
	// simple cut on the starting layer of the shower
	// FIXME: this is not working properly for the edges of the ECAL and HCAL => need pseudo layer !!!
	const unsigned int cutOnLayer = 10; // FIXME: remove hard-coded cut	
	const unsigned int contOnNumberOfHitsInLayerCut = 4; // FIXME: remove hard-coded cut

	unsigned int nHitsInLayerCut = 0;
	bool hitWithinCutRangeFound = false;

	for (unsigned int j = 0; j < cluster->getCalorimeterHits().size(); ++j) {

	  unsigned int layer = CDECAL(cluster->getCalorimeterHits().at(j))["K-1"];
	  
	  if ( layer <= cutOnLayer ) ++nHitsInLayerCut;

	  if ( nHitsInLayerCut >= contOnNumberOfHitsInLayerCut ) {
	    hitWithinCutRangeFound = true;
	    break;
	    
	  }

	}
	


	   
	// debug
	if ( _debugLevel > 5 ) {
	  std::cout << std::endl << "Cluster Size:" << cluster->getCalorimeterHits().size() << std::endl
		    << "energy of cluster: " << energy << "  " << "estimated energy: " << estimatedEnergyOfCluster << "  " << "deviation: " << energyDeviation << std::endl
		    << "Calorimeter Hist in this Cluster:" << std::endl;
	  for (unsigned int j = 0; j < cluster->getCalorimeterHits().size(); ++j) {
	    
	    std::cout << "Position: " << "(" << cluster->getCalorimeterHits().at(j)->getPosition()[0] << "," << cluster->getCalorimeterHits().at(j)->getPosition()[1] << "," 
		      << cluster->getCalorimeterHits().at(j)->getPosition()[2] << ")" << "  " << "E = " << cluster->getCalorimeterHits().at(j)->getEnergy() << "  "
		      << "layer: " << CDECAL(cluster->getCalorimeterHits().at(j))["K-1"] << std::endl;
	    
	    
	  }
	}
 
	    
	if ( (energyDeviation < _energyDeviationCut) &&  hitWithinCutRangeFound ) {


	      
	  // flag hits 
	  for (unsigned int j = 0; j < cluster->getCalorimeterHits().size(); ++j) cluster->getCalorimeterHits().at(j)->ext<isPartOfEMShowerCandidate>() = 1;
		      
	  ClusterShapes* clusterShape = new ClusterShapes(n,a,x,y,z);
	      
	  float position[3];
	  position[0]=clusterShape->getCentreOfGravity()[0];
	  position[1]=clusterShape->getCentreOfGravity()[1];
	  position[2]=clusterShape->getCentreOfGravity()[2];
	      
	      
	  cluster->setPosition(position); 
	  cluster->setEnergy(energy);
	  clscol->addElement(cluster);
	      
	      
	  // debug
	  if ( _debugLevel > 5 ) {
	    std::cout << "EM Shower Candidate FOUND : " << "cluster " << cluster << " with energy E = " << cluster->getEnergy() << "  " << "n of hits = " 
		      << cluster->getCalorimeterHits().size() << std::endl;
	  }
	      
	      
	  if (_drawOnCED) {
	    
	    int color = 0x8f57ff * ( iCluster + 32);	   
	    MarlinCED::drawClusterImpl(cluster,2,2,color,5); 
	    ced_send_event();
	    ++iCluster;
	    getchar();

	  }
	      
	      
	      
	  delete clusterShape;
	  clusterShape = 0;
	      
	}
	    
	    
	    
	delete[] x;
	x = 0;
	delete[] y;
	y = 0;
	delete[] z;
	z = 0;
	delete[] a;
	a = 0;
	    
      }
      
      mapCoreToCluster.clear();
      ECALHitsWithAttributes.clear();



      
      // for strong memory and nice dreams ..
      for(unsigned int i=0;i<N;i++)
	{
	  if( bbb[i].size()!=0)
	    for(unsigned int im=0;im<bbb[i].size();im++)
	      {
		delete bbb[i][im];
	      }
	}
      
      for(unsigned int im=0;im<2;im++)
	{        
	  if(calo[im].size()!=0)
	    for( unsigned int iij=0;iij<calo[im].size();iij++)
	      delete (calo[im])[iij];
	}
      
    }      
    
    evt->addCollection(clscol,_collectionNameOfEMShowerCandidates.c_str());
   
  }catch(DataNotAvailableException &e) {std::cout << "no valid ECAL collection in event " << _nEvt << std::endl; }
  

  // debug
  if (_drawOnCED) MarlinCED::draw(this,true);


  ++_nEvt;

}


void EMShowerFinder::check( LCEvent * evt ) {
 
}


void EMShowerFinder::end() {

}
