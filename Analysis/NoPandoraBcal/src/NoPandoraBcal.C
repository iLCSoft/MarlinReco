#include <math.h>
#include "NoPandoraBcal.h"

#include <EVENT/LCCollection.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/ParticleIDImpl.h>
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/BField.h>

using namespace lcio ;
using namespace marlin ;


NoPandoraBcal aNoPandoraBcal ;


NoPandoraBcal::NoPandoraBcal() : Processor("NoPandoraBcal") {
  
  // modify processor description
  _description = "Throws out PFOs reconstructed in the BCAL" ;
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "PFOCollectionName" ,
                           "Name of the input PFO collection"  ,
                           _PFOName ,
                           std::string("PandoraPFOs") ) ;
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "NewRecoCollectionName" ,
                           "Name of the new RECO collection"  ,
                           _newRECOName ,
                           std::string("PandoraPFOsnoBcal") ) ;
  registerInputCollection( LCIO::CLUSTER,
                           "ClusterCollectionName" ,
                           "Name of the input PFO collection"  ,
                           _CluName ,
                           std::string("PandoraClusters") ) ;
  registerOutputCollection( LCIO::CLUSTER,
                           "NewClusterCollectionName" ,
                           "Name of the new Cluster collection"  ,
                           _newCluName ,
                           std::string("PandoraClustersnoBcal") ) ;
}

void NoPandoraBcal::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;
  //number_of_events = 0;
  _nRun = 0 ;
  _nEvt = 0 ;
  std::cout << "Init done" << std::endl;
}

void NoPandoraBcal::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

//--------------------------------------------------------------------------------------------------------------------------------------

void NoPandoraBcal::processEvent( LCEvent * evt ) { 

  double pi ;
  pi = atan(1.)*4. ;  
  LCCollection* col_PFOs = evt->getCollection(_PFOName); 
  LCCollectionVec* NEWRECO = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  ReconstructedParticle* p=0;
  NEWRECO->setSubset(true);

  LCCollection* col_Clus = evt->getCollection(_CluName); 
  LCCollectionVec* NEWCLU = new LCCollectionVec(LCIO::CLUSTER);
  Cluster* clu=0;
  NEWCLU->setSubset(true);

  bool keep = true;
  FloatVec pe;
  int nPFOs=col_PFOs ->getNumberOfElements(); 
  for (int i=0; i<nPFOs; ++i)
    {
      keep = true;
      p = dynamic_cast<ReconstructedParticle*> ( col_PFOs -> getElementAt( i ) );
      const ClusterVec & Clusters = p->getClusters();
      streamlog_out(DEBUG1) << " pfo nb " << i << " energy " << p->getEnergy() << 
        " no of clusters " <<  Clusters.size() << std::endl ;
      if ( Clusters.size() > 0 ) 
	{
	  for (unsigned short j=0; j< Clusters.size(); ++j)
	    {
	      FloatVec pe = Clusters[j]->getSubdetectorEnergies();
	      streamlog_out(DEBUG1) << "     cluster nb " << j << " energies " 
                << pe[0] << "" 
                << pe[1] << " "
                << pe[2] << " "
                << pe[3] << " "
                << pe[4] << " "
		<< pe[5] << std::endl ;
	      if ( (pe[3] != 0.0 || pe[5] !=0.0 ) && keep)
		{
	          CalorimeterHitVec hv = Clusters[j]->getCalorimeterHits();
                  if ( hv.size() > 0 ) {
                    for ( unsigned jjj=0 ; jjj< hv.size() ; jjj++ ) {
                      double hitz=fabs(hv[jjj]->getPosition()[2]) ;
                      double hitr=sqrt(hv[jjj]->getPosition()[0]*hv[jjj]->getPosition()[0]+hv[jjj]->getPosition()[1]*hv[jjj]->getPosition()[1]);
                      streamlog_out(DEBUG1) << " hitz = " << hitz <<  " hitr = " << hitr << std::endl;
                      if ( (hitz> 2505. && hitz< 2635) &&  (hitr> 65 && hitr<215) ) {
                        streamlog_out(DEBUG2) << " hit nb: " << jjj << " hitz = " << hitz <<  " hitr = " << hitr << std::endl;
                        streamlog_out(DEBUG2) << "     SKIP ! " << std::endl;
                        keep = false;
                        break;
                      }
                    }
                  } else {
  		  if ( fabs(Clusters[j]->getPosition()[2]) > 2300. || pe[5] !=0.0 )
		    {
		      keep = false;
         	      streamlog_out(DEBUG2) << "     cluster nb " << j << " energies " 
                        << pe[0] << "" 
                        << pe[1] << " "
                        << pe[2] << " "
                        << pe[3] << " "
                        << pe[4] << " "
		        << pe[5] << std::endl ;
                     streamlog_out(DEBUG2) << "     SKIP ! " << std::endl;
		    }
                  }
		}
	    }
	}
      
      if (keep) 
	{
	  NEWRECO->addElement(p);
          for ( unsigned jjj=0 ; jjj<Clusters.size() ; jjj++ ) {
	    clu=Clusters[jjj] ;
  	    NEWCLU->addElement(clu);
            
          }    
	}
    }
  if (col_PFOs ->getNumberOfElements() != NEWRECO ->getNumberOfElements()) {
    streamlog_out(DEBUG3) << "Removed Particles from the PFO Collection: Before " << col_PFOs ->getNumberOfElements() << " after: " <<
      NEWRECO ->getNumberOfElements() << std::endl;
  }
  if (col_Clus ->getNumberOfElements() != NEWCLU ->getNumberOfElements()) {
    streamlog_out(DEBUG3) << "Removed Clusters from the Cluster Collection: Before " << col_Clus ->getNumberOfElements() << " after: " <<
      NEWCLU ->getNumberOfElements() << std::endl;
  }
  evt->addCollection(NEWRECO,_newRECOName) ;
  evt->addCollection(NEWCLU,_newCluName) ;
  _nEvt ++ ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void NoPandoraBcal::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NoPandoraBcal::end()
{ 

}
