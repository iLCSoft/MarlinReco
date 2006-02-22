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
// STUFF needed for GEAR
#include <marlin/Global.h>
#include "gear/GEAR.h"
//#include "gear/TPCParameters.h"
//#include "gear/PadRowLayout2D.h"
#include "gear/CalorimeterParameters.h"

#define MASK_K (unsigned int) 0x3F000000
#define SHIFT_K 24
#define SHIFT_M 0
#define MASK_M (unsigned int) 0x00000007
using namespace std ;
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
 registerProcessorParameter("TrueClusterToMCPCollection",
			     "Relation Collection Cluster to MCP",
			     _trueClustToMCP,
			     std::string("TrueClusterToMCP"));
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
 

  float totalUnrecoverableOverlapEcal1max=0.0;
  float totalUnrecoverableOverlapEcal2max=0.0;
  float totalUnrecoverableOverlapHcalmax=0.0;
  float totalUnrecoverableOverlapEcal1min=0.0;
  float totalUnrecoverableOverlapEcal2min=0.0;
  float totalUnrecoverableOverlapHcalmin=0.0;
  float noPointerEcal=0.0;
  float noPointerHcal=0.0;

 typedef std::map <MCParticle*,ClusterImpl*> mapMCP2Clust;
    
   mapMCP2Clust p2clust;

  for (unsigned int i=0 ; i < _caloCollections.size(); ++i) 
    {    

      cout<<"name ="<<_caloCollections[i]<<endl;
       const LCCollection * relcol ;
       relcol = evt->getCollection("RelationCaloHit");
	 LCCollection * col = evt->getCollection(_caloCollections[i].c_str());
          unsigned int nhits = col->getNumberOfElements();
  
     
      double gearRMax = Global::GEAR->getEcalBarrelParameters().getExtent()[0];
    
      double zmax= Global::GEAR->getEcalEndcapParameters().getExtent()[2];
   
 LCRelationNavigator navigate(relcol); 

	    for (unsigned int j=0; j < nhits; ++j) 
             {
	   CalorimeterHit * calhit = dynamic_cast<CalorimeterHit*>(col->getElementAt(j));
	LCObjectVec objectVec = navigate.getRelatedToObjects(calhit);
		if (objectVec.size() > 0) {
          
           SimCalorimeterHit * simhit=dynamic_cast<SimCalorimeterHit*>( objectVec[0]);
	   unsigned int ncontrib=simhit->getNMCContributions();
	   int cellid=calhit->getCellID0();
           int layer=(cellid & MASK_K) >> SHIFT_K;  
	   int module=(cellid & MASK_M) >> SHIFT_M;  
	   if ( ncontrib>1)
	     {   // so there is an overlap  
	       float maxmix=0.0;
               unsigned int ttmax=0;
               unsigned int ttmin=0;
	       float minmix=simhit->getEnergyCont(0);
	       for ( unsigned int tt=0 ; tt<ncontrib ; ++tt)
		 {
		   if (simhit->getEnergyCont(tt)> maxmix)
		     {
		       maxmix=simhit->getEnergyCont(tt);
		       ttmax=tt;
		     }
                   if (simhit->getEnergyCont(tt)< minmix)
		     {
		       minmix=simhit->getEnergyCont(tt);
		       ttmin=tt;
		     }
		 }
                 MCParticle * par =simhit->getParticleCont(ttmax);
		 if ( par !=NULL)
		     {

                      first_jump:
		     
		       if (p2clust[par] == NULL) 
			 {
                           
			   float x=par->getVertex()[0];
			   float y=par->getVertex()[1];
			   float z=par->getVertex()[2];
			   float radius=sqrt(x*x+y*y);
        	
			 if( (radius<gearRMax && fabs(z)<zmax) ) 
			     {
                	       ClusterImpl * cluster = new ClusterImpl();
			       p2clust[par]=cluster;
			       cluster->addHit(calhit,(float)1.0);
			     }else{
			        par=par->getParents()[0];
				goto first_jump;
  			     }
			 }else{
			      ClusterImpl * cluster = p2clust[par];
			      cluster->addHit(calhit,(float)1.0);
		         }         
	            
                     }else{// par !=0
		        noPointerEcal+=simhit->getEnergyCont(ttmax);
 
		     }//par !=0    

                        totalUnrecoverableOverlapEcal1max+=maxmix;
		        totalUnrecoverableOverlapEcal1min+=minmix;
	      }else{ // contributions
                    MCParticle * par =simhit->getParticleCont(0);
		 if ( par !=NULL)
		     {

                       second_jump:
		    
		       if (p2clust[par] == NULL) 
			 {
                           
			   float x=par->getVertex()[0];
			   float y=par->getVertex()[1];
			   float z=par->getVertex()[2];
			   float radius=sqrt(x*x+y*y);
        
		        if( (radius<gearRMax && fabs(z)<zmax) ) 
			     {
                	       ClusterImpl * cluster = new ClusterImpl();
			       p2clust[par]=cluster;
			       cluster->addHit(calhit,(float)1.0);
			     }else{
			        par=par->getParents()[0];
				goto second_jump;
  			     }
			 }else{
			      ClusterImpl * cluster = p2clust[par];
			      cluster->addHit(calhit,(float)1.0);
		         }         
	            
                     }else{// par !=0
		       
			noPointerEcal+=simhit->getEnergyCont(0);
		     }//par !=0    
                    	           
	          }// ncontribution
		}
	     }// over nhits
       
    }   

    		
  LCCollectionVec * clscol = new LCCollectionVec(LCIO::CLUSTER);
  LCCollectionVec * relationcol = new LCCollectionVec(LCIO::LCRELATION);

    mapMCP2Clust::iterator pos;

    for (pos = p2clust.begin(); pos != p2clust.end(); ++pos) 
    {
      MCParticle * mcp = pos->first;
   
      ClusterImpl * cluster = pos->second;
      
      if(cluster != NULL){
      CalorimeterHitVec calohitvec = cluster->getCalorimeterHits();
      int nhcl = (int)calohitvec.size();
    
	float totene = 0.0;
	float totecal = 0.0;
	float tothcal = 0.0;
	for (int i=0; i< nhcl; ++i) 
	  {
	  CalorimeterHit * calhit = calohitvec[i];
	  totene +=  calhit->getEnergy();	
	  }	
	cluster->setEnergy(totene);
	clscol->addElement(cluster);
   
	LCRelationImpl * rel = new LCRelationImpl(cluster,mcp,(float)1.0);
	relationcol->addElement( rel ); 

     }   
    }
    evt->addCollection(clscol,_trueClustCollection.c_str());
    evt->addCollection(relationcol,_trueClustToMCP.c_str());


    _nEvt++;
    
}


void ClusterCheater::check( LCEvent * evt ) { }
  
void ClusterCheater::end(){ } 

