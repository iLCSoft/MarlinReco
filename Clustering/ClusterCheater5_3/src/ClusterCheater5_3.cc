#include "ClusterCheater5_3.h"
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
#include <stack>
#include <math.h>
// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/CalorimeterParameters.h>
//#include "MarlinCED.h"
//#include "DrowUtil.h"
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
bool to_be_saved(const MCParticle * parm , const float cut , const CalorimeterHit* calhit); 
ClusterCheater5_3 aClusterCheater5_3 ;


ClusterCheater5_3::ClusterCheater5_3() : Processor("ClusterCheater5_3") {

  _description = "Creates true clusters..." ;

  registerOutputCollection( LCIO::CLUSTER,
			    "TrueClusterCollection",
			    "Collection of True Clusters",
			    _trueClustCollection ,
			    std::string("TrueClusters"));

  registerOutputCollection( LCIO::LCRELATION,
			    "TrueClusterToMCPCollection",
			    "Relation Collection Cluster to MCP",
			    _trueClustToMCP,
			    std::string("TrueClusterToMCP"));

  std::vector<std::string> caloCollections;

  caloCollections.push_back(std::string("ECAL"));
  caloCollections.push_back(std::string("HCAL"));
  

  registerInputCollections( LCIO::CALORIMETERHIT,
			    "CaloCollections",
			    "Calorimeter Collection Names",
			    _caloCollections ,
			    caloCollections);

  registerInputCollection( LCIO::LCRELATION,
			   "RelCollection",
			   "SimCaloHit to CaloHit Relations Collection Name",
			   _relCollection ,
			   std::string("RelationCaloHit"));

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection",
			   "Calorimeter Collection Names",
			   _MCcollection ,
			   std::string("MCParticle"));

  registerProcessorParameter("CutBackscatter",
			     "Not to connect hist from  backscatter",
			     _backcut ,
			     (int)1 );

 registerProcessorParameter("MinHitsInCluster",
			     "Minimal number of hits in cluster",
			     _Nmin ,
			     (int)0 );

}

void ClusterCheater5_3::init() {
    _nRun = -1;
    _nEvt = 0;
    // MarlinCED::init(this) ;
    _nlost=0;
  // const gear::CalorimeterParameters& pHcalBarrel = Global::GEAR->getHcalBarrelParameters();
 const gear::TPCParameters& pTPC = Global::GEAR->getTPCParameters();
 gearRMax =pTPC.getDoubleVal("tpcOuterRadius");
 zmax= pTPC.getMaxDriftLength();
}


void ClusterCheater5_3::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 

void ClusterCheater5_3::processEvent( LCEvent * evt ) { 
 
  //float cut =50.0;
  //float totalUnrecoverableOverlapEcal1max=0.0;
  //float totalUnrecoverableOverlapEcal2max=0.0;
  //float totalUnrecoverableOverlapHcalmax=0.0;
  //float totalUnrecoverableOverlapEcal1min=0.0;
  //float totalUnrecoverableOverlapEcal2min=0.0;
  //float totalUnrecoverableOverlapHcalmin=0.0;
  float noPointerEcal=0.0;
  //float noPointerHcal=0.0;
 
 typedef std::map <MCParticle*,ClusterImpl*> mapMCP2Clust;
     
   mapMCP2Clust p2clust;

  for (unsigned int i=0 ; i < _caloCollections.size(); ++i) 
    {    

      //  cout<<"name ="<<_caloCollections[i]<<endl;
       const LCCollection * relcol ;
       relcol = evt->getCollection(_relCollection.c_str());
	 LCCollection * col = evt->getCollection(_caloCollections[i].c_str());
          unsigned int nhits = col->getNumberOfElements();
  
     
     
    
     
	  LCRelationNavigator navigate(relcol); 

	    for (unsigned int j=0; j < nhits; ++j) 
             {
	   CalorimeterHit * calhit = dynamic_cast<CalorimeterHit*>(col->getElementAt(j));
	   LCObjectVec objectVec = navigate.getRelatedToObjects(calhit);

	   if (objectVec.size() > 0) {
          
	     vector<MCParticle*> vpar;
             vector<float>contrib;
	     SimCalorimeterHit * simhit=0;
	     unsigned int ncontrib;
	     vector<MCParticle*>::iterator testp;
	     for( unsigned int k=0;k<objectVec.size();++k)
	       {
		simhit=dynamic_cast<SimCalorimeterHit*>( objectVec[k]);
		ncontrib=simhit->getNMCContributions();
                
                for ( unsigned int kk=0; kk<ncontrib;++kk)
		  {
		    MCParticle * par=simhit->getParticleCont(kk);
                 
                    int countp=0;
                    testp=find(vpar.begin(),vpar.end(),par);
		    countp=testp-vpar.begin();

		    if( testp == vpar.end()) 
			{
			  vpar.push_back(par);
			  contrib.push_back(simhit->getEnergyCont(kk));
			}else{
			 contrib[countp]+=simhit->getEnergyCont(kk);
		        }
		  }
	       }
	     vector<float>::iterator naj;
	     naj=max_element(contrib.begin(),contrib.end());
	    
              
             int max=naj-contrib.begin();
	    
	     MCParticle * par =vpar[max];
	   
		 if ( par !=NULL)
		     {

                       second_jump:
		    
		       if (p2clust[par] == NULL) 
			 {
                           
			   float x=par->getVertex()[0];
			   float y=par->getVertex()[1];
			   float z=par->getVertex()[2];
			   float radius=sqrt(x*x+y*y);
			         x=par->getEndpoint()[0];
			         y=par->getEndpoint()[1];
			   float z1=par->getEndpoint()[2];
			   float radius1=sqrt(x*x+y*y);


			   if( ((radius<gearRMax && fabs(z)<zmax) && 
			       (radius1>gearRMax || fabs(z1)>zmax )) || 
			       (par->isBackscatter()&&_backcut) ) 
 			     {  
		
                	       ClusterImpl * cluster = new ClusterImpl();
			       p2clust[par]=cluster;
			       //  if( to_be_saved(par , cut,calhit )) 
			       cluster->addHit(calhit,(float)1.0);
			      }else{ 
			             if(par->getParents().size()!=0)
				       {
					 par=par->getParents()[0];
					 goto second_jump;
				       }else{
				       ++_nlost;
				       }
  			     }
			 }else{
			      ClusterImpl * cluster = p2clust[par];
			      //  if( to_be_saved(par , cut,calhit )) 
			      cluster->addHit(calhit,(float)1.0);
		         }         
	            
                     }else{// par !=0
		       
			noPointerEcal+=simhit->getEnergyCont(0);
		     }//par !=0    
                    	           
	     
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
// 	float totecal = 0.0;
// 	float tothcal = 0.0;
	if( nhcl>_Nmin)
	  {
	    float center[3];
            center[0]=0.0;    center[1]=0.0;    center[2]=0.0;
	    for (int i=0; i< nhcl; ++i) 
	      {
		CalorimeterHit * calhit = calohitvec[i];
		float energy=calhit->getEnergy();
		center[0]+= (calhit->getPosition()[0])*energy;
		center[1]+= (calhit->getPosition()[1])*energy;
		center[2]+= (calhit->getPosition()[2])*energy;
		totene += energy;	
	      }	

	    center[0]=center[0]/totene;
	    center[1]=center[1]/totene;
	    center[2]=center[2]/totene;
	    cluster->setPosition(center);
	    cluster->setEnergy(totene);
	    clscol->addElement(cluster);

	    LCRelationImpl * rel = new LCRelationImpl(cluster,mcp,(float)1.0);
	    relationcol->addElement( rel ); 
	  }
     }   
    }
    evt->addCollection(clscol,_trueClustCollection.c_str());
    evt->addCollection(relationcol,_trueClustToMCP.c_str());

    // cout << " n clusters at the end " <<  clscol->getNumberOfElements() << endl;
    _nEvt++;
    
}


void ClusterCheater5_3::check( LCEvent * evt ) { }
  
void ClusterCheater5_3::end()
{  

  cout << "cluster cheater  n hits lost = " << _nlost << endl;
} 

bool to_be_saved(const MCParticle * parm , const float cut,const CalorimeterHit * calhit ) 
{
   

 vector<MCParticle*> particlesToRec ;
 
 for( unsigned int pp=0;pp<parm->getDaughters().size();++pp)
    {
        MCParticle * par =parm->getDaughters()[pp];
	
	 stack <MCParticle*> krk;
           
     if(1) 
 
	   {
            int testz=count(particlesToRec.begin(),particlesToRec.end(),par);
	     if (testz==0)
	       {
		
		 // particlesToRec.push_back(par);
	         krk.push(par);
	       }
	   }
	     while(!krk.empty())
	       {
                
		     par=krk.top();    
                      krk.pop();
		   
		     if(1)
		 
		      {
                           int testz=count(particlesToRec.begin(),particlesToRec.end(),par);

			   if (testz==0)
			     {		 
			       particlesToRec.push_back(par);			        
			     }
			   for( unsigned int i=0;i<par->getDaughters().size();++i)
			     {
			       krk.push(par->getDaughters()[i]);
			       particlesToRec.push_back(par->getDaughters()[i]);
			     }
		      }


	       }//while
	  
    }

       vector<MCParticle*> particlesToBack ;
     	 stack <MCParticle*> krk;
 for ( unsigned int i=0; i< particlesToRec.size();++i)
   {
     if( particlesToRec[i]->isBackscatter() )
       {
	 particlesToBack.push_back( particlesToRec[i]);
	 krk.push(particlesToRec[i]);
       }
   }

          while(!krk.empty())
	       {
                
		  MCParticle*  par=krk.top();    
		    krk.pop();
		   
		    
                           int testz=count(particlesToBack.begin(),particlesToBack.end(),par);

			   if (testz==0)
			     {		 
			       particlesToBack.push_back(par);			        
			     }
			   for( unsigned int i=0;i<par->getDaughters().size();++i)
			     {
			       krk.push(par->getDaughters()[i]);
			       particlesToBack.push_back(par->getDaughters()[i]);
			     }
		   


	       }//while

	  // cout << "N BCKSCATTER = "<<particlesToBack.size()<< endl;
 bool test=true;
 float x1=calhit->getPosition()[0];
 float y1=calhit->getPosition()[1];
 float z1=calhit->getPosition()[2];
	 float px=parm->getMomentum()[0];
         float py=parm->getMomentum()[1];
         float pz=parm->getMomentum()[2];
  for ( unsigned int i=0; i< particlesToBack.size();++i)
   {   
   	 float x=particlesToBack[i]->getEndpoint()[0];
         float y=particlesToBack[i]->getEndpoint()[1];
         float z=particlesToBack[i]->getEndpoint()[2];



	 float dd=px*x1+py*y1+pz*z1;
	 float d=sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1));

         if ( d< cut || dd<0.0) 
	   {
             test=false;
	  }                    
   }
  return test;

}
