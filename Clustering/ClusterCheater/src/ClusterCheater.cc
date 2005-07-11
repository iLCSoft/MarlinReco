#include "ClusterCheater.h"
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include <iostream>
#include <map>
#include <vector>

using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

typedef std::map <MCParticle*,ClusterImpl*> map_MCP_Clust;

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
  

  registerProcessorParameter("ECALCollection",
			     "ECAL Collection Name",
			     _caloCollections ,
			     caloCollections);

  registerProcessorParameter("RelCollection",
			     "Relation Collection Name",
			     _relCollection ,
			     std::string("RelationCaloHit"));

  registerProcessorParameter("IfBrahms",
			     "If Brahms",
			     _ifBrahms,
			     (int)1);

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
    
    try {
	relcol = evt->getCollection(_relCollection.c_str());	
    }
    catch(DataNotAvailableException &e){
	std::cout << "No relation exist between Sim and Calo hit" << std::endl; 	std::cout << "Quitting program" << std::endl; 
	exit(1);
    }

    LCRelationNavigator navigate(relcol); 

    LCCollectionVec * clscol = new LCCollectionVec(LCIO::CLUSTER);

    for (unsigned int i(0) ; i < _caloCollections.size(); ++i) {    
	try {
	    LCCollection * col = evt->getCollection(_caloCollections[i].c_str());
	    int nelem = col->getNumberOfElements();
	    for (int i(0); i < nelem; ++i) {
		CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(i));
		LCObjectVec objectVec = navigate.getRelatedToObjects(hit);
		if (objectVec.size() > 0) {
		    SimCalorimeterHit * shit = dynamic_cast<SimCalorimeterHit*>(objectVec[0]);		    
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
				clscol->addElement(cluster);
				_mcp_clust[par] = cluster;
				cluster->addHit(hit,(float)1.0);
			      }
			      else {
				ClusterImpl * cluster = _mcp_clust[par];
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
				_isDecayedInTracker = parent->isDecayedInTracker();
			      }
			      if (nparents == 0 || _isDecayedInTracker ) {
				if (_mcp_clust[par] == NULL) {
				  ClusterImpl * cluster = new ClusterImpl();
				  clscol->addElement(cluster);
				  _mcp_clust[par] = cluster;
				  cluster->addHit(hit,(float)1.0);
				}
				else {
				  ClusterImpl * cluster = _mcp_clust[par];
				  cluster->addHit(hit,(float)1.0);
				}
				loop = 0;
				break;
			      }
			    par = parent;
			    } 
			  }
			}
		    }
		}
	    }	
	}
	catch(DataNotAvailableException &e){ 
	}
    }
    evt->addCollection(clscol,"TrueClusters");

    LCCollectionVec * relationcol = new LCCollectionVec(LCIO::LCRELATION);

    map_MCP_Clust::iterator pos;
    for (pos = _mcp_clust.begin(); pos != _mcp_clust.end(); ++pos) {
	MCParticle * mcp = pos->first;
	ClusterImpl * cluster = pos->second;
	LCRelationImpl * rel = new LCRelationImpl(cluster,mcp,(float)1.0);
	relationcol->addElement( rel );
    }
    evt->addCollection(relationcol,"TrueClusterToMCP");

  _nEvt++;

}


void ClusterCheater::check( LCEvent * evt ) { }
  
void ClusterCheater::end(){ } 
