#include "ClusterCheater.h"
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include <iostream>
#include "ClusterShapes.h"
#include <map>
#include <vector>
#include <math.h>

using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

typedef std::map <MCParticle*,ClusterImpl*> map_MCP_Clust;
typedef std::map <MCParticle*,HelixClass*> map_MCP_Helix;

ClusterCheater aClusterCheater ;


ClusterCheater::ClusterCheater() : Processor("ClusterCheater") {

  _description = "Creates true clusters..." ;

  registerProcessorParameter("TrueClusterCollection",
			     "Collection of True Clusters",
			     _trueClustCollection ,
			     std::string("TrueClusters"));

  std::vector<std::string> caloCollections;

  caloCollections.push_back(std::string("ECAL"));
  caloCollections.push_back(std::string("HCAL"));
  

  registerProcessorParameter("CaloCollections",
			     "Calorimeter Collection Names",
			     _caloCollections ,
			     caloCollections);

  registerProcessorParameter("RelCollection",
			     "SimCaloHit to CaloHit Relations Collection Name",
			     _relCollection ,
			     std::string("RelationCaloHit"));

  registerProcessorParameter("TraceWholeShower",
			     "Trace Whole Shower Tree",
			     _ifBrahms,
			     (int)1);

  registerProcessorParameter("ProximityCut",
			     "Hit To Cluster Proximity Cut",
			     _proximityCut,
			     (float)1200.);

  registerProcessorParameter("MinimalHits",
			     "Minimal Hits in Cluster",
			     _minimal_hits,
			     (int)10);

  registerProcessorParameter("MagneticField",
			     "Magnetic Field",
			     _bField,
			     (float)4.0);
    

}

void ClusterCheater::init() {
    _nRun = -1;
    _nEvt = 0;
}


void ClusterCheater::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 

void ClusterCheater::processEvent( LCEvent * evt ) { 
    const LCCollection * relcol ;

    map_MCP_Clust _mcp_clust;
    map_MCP_Helix _mcp_helix;
    
    try {
	relcol = evt->getCollection(_relCollection.c_str());	
    }
    catch(DataNotAvailableException &e){
	std::cout << "Cluster Cheater : " 
		  << "No relation exist " 
		  << "between Sim and Calo hit" << std::endl; 	
	std::cout << "Cluster Cheater : Quitting program" << std::endl; 
	exit(1);
    }

    LCRelationNavigator navigate(relcol); 

    for (unsigned int i(0) ; i < _caloCollections.size(); ++i) {    
	try {
	    LCCollection * col = evt->getCollection(_caloCollections[i].c_str());
	    int nelem = col->getNumberOfElements();
	    for (int i(0); i < nelem; ++i) {
		CalorimeterHit * hit = 
		  dynamic_cast<CalorimeterHit*>(col->getElementAt(i));
		LCObjectVec objectVec = navigate.getRelatedToObjects(hit);
		if (objectVec.size() > 0) {
		    SimCalorimeterHit * shit = 
		      dynamic_cast<SimCalorimeterHit*>(objectVec[0]);
		    if (shit->getNMCParticles() > 0) {
			MCParticle * par = shit->getParticleCont(0);
			if (par != NULL ) {
			  if (_ifBrahms == 1) {			    
			    bool swi = par->isCreatedInSimulation();
			    if (swi) {
			    }
			    else {
			      if (_mcp_clust[par] == NULL) {
				ClusterImpl * cluster = new ClusterImpl();
				_mcp_clust[par] = cluster;
				float charge = par->getCharge();
				float Dist;
				if (fabs(charge)>0.5) {
				  HelixClass * helix = AssignHelixToMCP(par);
				  _mcp_helix[par] = helix;
				  Dist = DistanceToChargeParticle(helix,hit);
				}
				else {
				  Dist = DistanceToNeutralParticle(par,hit);
				}
				if (Dist < _proximityCut)
				  cluster->addHit(hit,(float)1.0);
			      }
			      else {
				ClusterImpl * cluster = _mcp_clust[par];
				float charge = par->getCharge();
				float Dist;
				if (fabs(charge)>0.5) {
				  HelixClass * helix = _mcp_helix[par];
				  Dist = DistanceToChargeParticle(helix,hit);
				}
				else {
				  Dist = DistanceToNeutralParticle(par,hit);
				}
				if (Dist < _proximityCut)
				  cluster->addHit(hit,(float)1.0);
			      }
			    }
			  }
			  else {
			    bool loop = 1;
			    while (loop) {
			      int nparents = par->getNumberOfParents();
			      bool _isDecayedInTracker = 0;
			      MCParticle * parent ;
			      if (nparents != 0) {
				parent = par->getParent(0);
				_isDecayedInTracker = 
				  parent->isDecayedInTracker();
			      }
			      if (nparents == 0 || _isDecayedInTracker ) {
				if (_mcp_clust[par] == NULL) {
				  ClusterImpl * cluster = new ClusterImpl();
				  _mcp_clust[par] = cluster;
				  float charge = par->getCharge();
				  float Dist;
				  if (fabs(charge)>0.5) {
				    HelixClass * helix = AssignHelixToMCP(par);
				    _mcp_helix[par] = helix;
				    Dist = DistanceToChargeParticle(helix,hit);
				  }
				  else {
				    Dist = DistanceToNeutralParticle(par,hit);
				  }
				  if (Dist < _proximityCut)
				    cluster->addHit(hit,(float)1.0);
				}
				else {
				  ClusterImpl * cluster = _mcp_clust[par];
				  float charge = par->getCharge();
				  float Dist;
				  if (fabs(charge)>0.5) {
				    HelixClass * helix = _mcp_helix[par];
				    Dist = DistanceToChargeParticle(helix,hit);
				  }
				  else {
				    Dist = DistanceToNeutralParticle(par,hit);
				  }
				  if (Dist < _proximityCut)
				    cluster->addHit(hit,(float)1.0);
				}
				loop = 0;
				break;
			      }
			    par = parent;
			    } 
			  }
			}
			else {
			  std::cout << "Hit is lost !!!!! " << std::endl;
			}
		    }
		    else {
		      std::cout << "Hit is lost !!!!! " << std::endl;
		    }
		}
		else {
		  std::cout << "Hit is lost !!!!! " << std::endl;
		}
	    }	
	}
	catch(DataNotAvailableException &e){ 
	}
    }

    LCCollectionVec * clscol = new LCCollectionVec(LCIO::CLUSTER);
    LCCollectionVec * relationcol = new LCCollectionVec(LCIO::LCRELATION);

    map_MCP_Clust::iterator pos;
    for (pos = _mcp_clust.begin(); pos != _mcp_clust.end(); ++pos) {
      MCParticle * mcp = pos->first;
      ClusterImpl * cluster = pos->second;
      if (fabs(mcp->getCharge()) > 0.5) {
	HelixClass * helix = _mcp_helix[mcp];
	delete helix;
      }
      CalorimeterHitVec calohitvec = cluster->getCalorimeterHits();
      int nhcl = (int)calohitvec.size();
      if (nhcl > _minimal_hits) {
	float * xhit = new float[nhcl];
	float * yhit = new float[nhcl];
	float * zhit = new float[nhcl];
	float * ahit = new float[nhcl];
	float totene = 0.0;
	float totecal = 0.0;
	float tothcal = 0.0;
	for (int ihit(0); ihit < nhcl; ++ihit) {
	  CalorimeterHit * calhit = calohitvec[ihit];
	  xhit[ihit] = calhit->getPosition()[0];
	  yhit[ihit] = calhit->getPosition()[1];
	  zhit[ihit] = calhit->getPosition()[2];
	  ahit[ihit] = calhit->getEnergy();
	  totene += ahit[ihit];
	  if (calhit->getType() == 0) {
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
	float PhiCluster = 
	  atan2(shape->getEigenVecInertia()[1],
		shape->getEigenVecInertia()[0]);
	float ThetaCluster = acos(shape->getEigenVecInertia()[2]);
	cluster->setIPhi(PhiCluster);
	cluster->setITheta(ThetaCluster);	    
	clscol->addElement(cluster);
	delete shape;
	delete[] xhit;
	delete[] yhit;
	delete[] zhit;
	delete[] ahit;	      
	LCRelationImpl * rel = new LCRelationImpl(cluster,mcp,(float)1.0);
	relationcol->addElement( rel );
      }
    }
    evt->addCollection(clscol,"TrueClusters");
    evt->addCollection(relationcol,"TrueClusterToMCP");
    _nEvt++;
    
}


void ClusterCheater::check( LCEvent * evt ) { }
  
void ClusterCheater::end(){ } 

float ClusterCheater::DistanceToChargeParticle(HelixClass * helix, 
					       CalorimeterHit * hit) {
  float Distance[3];
  float pos[3];
  for (int i=0; i<3; ++i)
    pos[i]=(float)hit->getPosition()[i];    
  float Time = helix->getDistanceToPoint(pos,Distance);

  float Dist = 1.0e+20;
  if (Time > 0.0)
    Dist = Distance[2];
  return Dist;
}


float ClusterCheater::DistanceToNeutralParticle(MCParticle * par, 
						CalorimeterHit * hit) {
  float distance[3];
  float momentum[3];
  float absDistance = 0.0;
  float absMomentum = 0.0;
  float product = 0.0;
  for (int i=0; i<3; ++i) {
    distance[i] = (float)hit->getPosition()[i]-(float)par->getVertex()[i];
    momentum[i] = (float)par->getMomentum()[i];
    product += distance[i]*momentum[i];
    absDistance += distance[i]*distance[i];
    absMomentum += momentum[i]*momentum[i];
  }
  absDistance = sqrt(absDistance);
  absMomentum = sqrt(absMomentum);
  float cosAngle = product/fmax(1.0e-30,absDistance*absMomentum);
  if (cosAngle > 1.0) 
    cosAngle = 0.999999;  
  if (cosAngle < -1.0)
    cosAngle = -0.999999;
  float angle = acos(cosAngle);
  return absDistance*angle;

}


HelixClass*  ClusterCheater::AssignHelixToMCP( MCParticle * par) {
  HelixClass * helix = new HelixClass();
  float vertex[3];
  float momentum[3];
  for (int i=0; i<3; ++i) {
    vertex[i] = (float)par->getVertex()[i];
    momentum[i] = (float)par->getMomentum()[i];
  }
  float charge = par->getCharge();
  helix->Initialize_VP(vertex,momentum,charge,_bField);
  return helix;

}
