#include "CalorimeterFinalStage.h"
#include "marlin/Global.h"
#include <math.h>

#include "IMPL/LCCollectionVec.h"
#include "IMPL/ClusterImpl.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCRelationImpl.h"
#include "EVENT/Cluster.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCCollection.h"
#include "ClusterShapes.h"

using namespace lcio;
using namespace marlin;
using namespace std;


CalorimeterFinalStage aCalorimeterFinalStage;


CalorimeterFinalStage::CalorimeterFinalStage() : Processor("CalorimeterFinalStage") {
  
  // Processor description
  _description = "CalorimeterFinalStage ";
  

  // Register steering parameters: name, description, class-variable, default value

  registerInputCollection(LCIO::CLUSTER,
			  "ClusterInputCollection", 
			  "Cluster Input Collection Name",
			  _clusterInput,
			  std::string("CalorimeterStage3Clusters"));
  
 

  registerOutputCollection(LCIO::CLUSTER,
			   "ClusterOutputCollection", 
			   "Cluster Output Collection Name",
			   _clusterOutput,
			   std::string("ClustersMagic"));
  
  
  registerProcessorParameter("MinimalHits",
			     "Minimal Hits in Cluster",
			     _nhit_minimal,
			     (int)0);

}


void CalorimeterFinalStage::init() { 

  // Print processor parameters
  printParameters();

  // Set number of runs and events to 0
  _nRun = 0;
  _nEvt = 0;
  
}

void CalorimeterFinalStage::processRunHeader( LCRunHeader* run) { 

  _nRun++;
} 

void CalorimeterFinalStage::processEvent( LCEvent * evt ) { 


  try{
    LCCollection * inputCollection = evt->getCollection(_clusterInput.c_str());
    int nclust = inputCollection->getNumberOfElements();
    LCCollectionVec * clucol = new LCCollectionVec(LCIO::CLUSTER);    
    for (int iclust(0); iclust < nclust; ++iclust) {
	Cluster * clusterInput = 
	  dynamic_cast<Cluster*>(inputCollection->getElementAt(iclust));
	CalorimeterHitVec calohitvec = clusterInput->getCalorimeterHits();	  
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
		CalorimeterHit * calhit = calohitvec[ihit];
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
    evt->addCollection(clucol,_clusterOutput.c_str());
  }
  catch(DataNotAvailableException &e){}
  

  _nEvt ++;
}

void CalorimeterFinalStage::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void CalorimeterFinalStage::end(){ 
  
  std::cout << "CalorimeterFinalStage::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    << std::endl;

}

