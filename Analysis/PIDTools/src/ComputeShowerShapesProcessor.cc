#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <vector>

#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/CalorimeterParameters.h>

#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ClusterImpl.h"
#include <lcio.h>

#include "CalorimeterHitType.h"

//#include "Api/PandoraApi.h"

#include "ComputeShowerShapesProcessor.hh"

ComputeShowerShapesProcessor aComputeShowerShapesProcessor ;

ComputeShowerShapesProcessor::ComputeShowerShapesProcessor()
  : Processor("ComputeShowerShapesProcessor") {
  
  // Processor description
  _description = "Cluster Shower Profile extraction using Fitting" ;
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "PFOCollection",
                           "PFO collection name",
                           _PfoCollection,
                           std::string("PandoraPFOs"));
  
  registerProcessorParameter("ClusterCollectionName",
			     "Cluster collection name",
			     _ClusterCollection,
			     std::string("PandoraClusters"));
  
  registerProcessorParameter("RadiationLength_Ecal",
			     "RadiationLength of Absorbers",
			     _X01,
			     float(3.50));
  
  registerProcessorParameter("RadiationLength_Hcal",
			     "RadiationLength of Absorbers",
			     _X02,
			     float(17.57));
  
  registerProcessorParameter("MoliereRadius_Ecal",
			     "Moliere radius of Absorbers",
			     _Rm1,
			     float(9.00));
  
  registerProcessorParameter("MoliereRadius_Hcal",
			     "Moliere radius of Absorbers",
			     _Rm2,
			     float(17.19));
  
  
  
} 

void ComputeShowerShapesProcessor::init( LCEvent * evt ) { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  
  _PFOCol = evt->getCollection( _ClusterCollection ) ;
  //_myShowerShapes = new ComputeShowerShapes(); 
  //define variable names
  std::vector<std::string> ClusterShapeNames;
  ClusterShapeNames.push_back("chi2");
  ClusterShapeNames.push_back("max_Ed");
  ClusterShapeNames.push_back("showerMax");
  ClusterShapeNames.push_back("absorption_Length");
  ClusterShapeNames.push_back("showerMax_photon");
  ClusterShapeNames.push_back("showerMax_ratio");
  ClusterShapeNames.push_back("Rm");
  ClusterShapeNames.push_back("showerstart");
  ClusterShapeNames.push_back("functionstart");
  ClusterShapeNames.push_back("a");
  ClusterShapeNames.push_back("b");
  ClusterShapeNames.push_back("c");
  ClusterShapeNames.push_back("d");
  ClusterShapeNames.push_back("max_Ed_hit");
  ClusterShapeNames.push_back("showerMax_hit");
  ClusterShapeNames.push_back("xt90");
  ClusterShapeNames.push_back("xl20");
  
  _PFOCol->parameters().setValues("ClusterShapeParameters",ClusterShapeNames);
  
  printParameters();
  
}

void ComputeShowerShapesProcessor::processRunHeader( LCRunHeader* run) { 
} 

void ComputeShowerShapesProcessor::processEvent( LCEvent * evt ) { 
  _PFOCol = evt->getCollection( _PfoCollection ) ;
  int nClusters = _PFOCol->getNumberOfElements();
  
  const unsigned int ecal_Index(0) ;
  const unsigned int hcal_Index(1) ;
  const unsigned int yoke_Index(2) ;
  const unsigned int lcal_Index(3) ;
  const unsigned int lhcal_Index(4);
  const unsigned int bcal_Index(5) ;
  
  for (int iClu=0;iClu<nClusters;++iClu) {
    ReconstructedParticle* part=(ReconstructedParticle*) _PFOCol->getElementAt( iClu );

    for(unsigned int jClu=0;jClu<part->getClusters().size();jClu++){
      ClusterImpl* pCluster = (ClusterImpl*) part->getClusters()[jClu];  //only use first cluster shape

      const unsigned int nHitsInCluster(pCluster->getCalorimeterHits().size());

      float clusterEnergy(0.);
      float *pHitE = new float[nHitsInCluster];
      float *pHitX = new float[nHitsInCluster];
      float *pHitY = new float[nHitsInCluster];
      float *pHitZ = new float[nHitsInCluster];
      int *typ = new int[nHitsInCluster];

      for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit)
	{
	  EVENT::CalorimeterHit *pCalorimeterHit = (CalorimeterHit*)(pCluster->getCalorimeterHits()[iHit]);

	  const float caloHitEnergy(pCalorimeterHit->getEnergy());
	  
	  pHitE[iHit] = caloHitEnergy;
	  pHitX[iHit] = pCalorimeterHit->getPosition()[0];
	  pHitY[iHit] = pCalorimeterHit->getPosition()[1];
	  pHitZ[iHit] = pCalorimeterHit->getPosition()[2];
	  clusterEnergy += caloHitEnergy;

	  switch (CHT(pCalorimeterHit->getType()).caloID())
	    {
	    case CHT::ecal:  typ[iHit]=ecal_Index; break;
	    case CHT::hcal:  typ[iHit]=hcal_Index; break;
	    case CHT::yoke:  typ[iHit]=yoke_Index; break;
	    case CHT::lcal:  typ[iHit]=lcal_Index; break;
	    case CHT::lhcal: typ[iHit]=lhcal_Index; break;
	    case CHT::bcal:  typ[iHit]=bcal_Index; break;
	    default: streamlog_out(DEBUG) << " no subdetector found for hit with type: " << pCalorimeterHit->getType() << std::endl;
	    }
	}
    
      pClusterShapes = new ClusterShapes(nHitsInCluster, pHitE, pHitX, pHitY, pHitZ);
      pClusterShapes->setHitTypes(typ);   //set hit types

      //here is cluster shape study - cluster transverse & longitudinal information
      //define variables
      float chi2,a,b,c,d,xl0,CoG[3],xStart[3]; //for fitting parameters
      float X0[2]={0,0};  //in mm. //this is the exact value of tangsten and iron
      float Rm[2]={0,0};  //in mm. need to change to estimate correctly times 2
      X0[0]=_X01;
      X0[1]=_X02;
      Rm[0]=_Rm1;
      Rm[1]=_Rm2;

      //get barrel detector surfce
      const gear::CalorimeterParameters& eCalDet = marlin::Global::GEAR->getEcalBarrelParameters(); 
      const EVENT::DoubleVec& ecal_ext = eCalDet.getExtent();
      float ecalrad=(float)ecal_ext[0];    //1.847415655e+03;   //in mm 

      //get endcap detector surfce
      const gear::CalorimeterParameters& pCalDet = marlin::Global::GEAR->getEcalEndcapParameters(); 
      const EVENT::DoubleVec& pcal_ext = pCalDet.getExtent();
      float plugz=(float)pcal_ext[2];    //2.450000000e+03;   //in mm 
      
      //looking for the hit which corresponds to the nearest hit from IP in the direction of the center of gravity
      int index_xStart=0;
      float lCoG=0.0,tmpcos=0.0,tmpsin=0.0,detsurface=0.0;
      CoG[0]=pClusterShapes->getEigenVecInertia()[0];
      CoG[1]=pClusterShapes->getEigenVecInertia()[1];
      CoG[2]=pClusterShapes->getEigenVecInertia()[2];
      //CoG2[0]=pCluster->getPosition()[0];
      //CoG2[1]=pCluster->getPosition()[1];
      //CoG2[2]=pCluster->getPosition()[2];
      
      lCoG=sqrt(CoG[0]*CoG[0]+CoG[1]*CoG[1]+CoG[2]*CoG[2]);
      tmpcos=CoG[2]/lCoG;
      tmpsin=sqrt(CoG[0]*CoG[0]+CoG[1]*CoG[1])/lCoG;
      pClusterShapes->fit3DProfile(chi2,a,b,c,d,xl0,xStart,index_xStart,X0,Rm);  //is this good??
      float lxstart=sqrt(xStart[0]*xStart[0]+xStart[1]*xStart[1]);
      //calculate detector surface
      if(fabs(xStart[2])<plugz){   //if in the barrel
	detsurface=(lxstart-ecalrad)/tmpsin;
      }else{  //if in plug
	detsurface=(fabs(xStart[2])-plugz)/fabs(tmpcos);
      }
      
      //float maxed=a*pow(b/c,b)*exp(-b);   //for simple fitting
      float maxed = a*c*gsl_sf_gammainv(b)*pow(b-1,b-1)*exp(-b+1);  //for advanced multiply with fabs(d) to avoid NaN
      float maxlength_pho=(1.0*std::log(clusterEnergy/(X0[0] * 0.021/Rm[0]))-0.5);  //this definition, +0.5 if gamma
      EVENT::FloatVec shapes;
      //these variables are fit based variables
      shapes.push_back(chi2);
      shapes.push_back(maxed);
      shapes.push_back(((b-1.0)*X0[0]/c+xl0+detsurface)/(2.0*X0[0]));
      shapes.push_back(1/fabs(d));
      shapes.push_back(maxlength_pho);
      shapes.push_back(((b-1.0)/c)/maxlength_pho);
      shapes.push_back(Rm[0]*2.0);
      shapes.push_back(detsurface);
      shapes.push_back(xl0);
      shapes.push_back(a);
      shapes.push_back(b);
      shapes.push_back(c);
      shapes.push_back(d);
      //these variables are detector based variables
      shapes.push_back(pClusterShapes->getEmax(xStart,index_xStart,X0,Rm));
      shapes.push_back(pClusterShapes->getsmax(xStart,index_xStart,X0,Rm));
      shapes.push_back(pClusterShapes->getxt90(xStart,index_xStart,X0,Rm));
      shapes.push_back(pClusterShapes->getxl20(xStart,index_xStart,X0,Rm));
      
      //add shower shapes
      pCluster->setShape(shapes);

      delete pClusterShapes;
      delete[] pHitE; delete[] pHitX; delete[] pHitY; delete[] pHitZ;
    }
  }  
}

void ComputeShowerShapesProcessor::check( LCEvent * evt ) { 
}

void ComputeShowerShapesProcessor::end() { 
  //delete _myShowerShapes;
}
 
