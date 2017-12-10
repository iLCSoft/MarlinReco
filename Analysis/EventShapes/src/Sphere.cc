#include "Sphere.h"
#include <iostream>
#include <vector>
#include "jama_eig.h"
#include "tnt_math_utils.h"
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include "marlin/Exceptions.h"
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif
#include "IMPL/LCEventImpl.h" 
#include "IMPL/LCCollectionVec.h"
#include "IMPL/SimCalorimeterHitImpl.h"
#include "IMPL/CalorimeterHitImpl.h"
#include "IMPL/MCParticleImpl.h" 
#include "IMPL/TrackerHitImpl.h" 
#include "IMPL/TrackImpl.h" 
#include "IMPL/ClusterImpl.h" 
#include "IMPL/ReconstructedParticleImpl.h" 
#include "IMPL/ParticleIDImpl.h" 
#include "IMPL/LCFlagImpl.h" 

#include "IMPL/LCRelationImpl.h"


#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>

using namespace lcio ;
using namespace marlin ;
using namespace std;
using namespace JAMMA;
using namespace TNT;
Sphere  aSphere ;


Sphere::Sphere() : Processor("Sphere") {
  
  // modify processor description
  _description = "Sphere calculates eigenvalues of sphericity tensor" ;
  
   
  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			   "CollectionName" , 
			   "Name of the ReconstructedParticle collection"  ,
			   _colName ,
			   std::string("RecontructedParticle") ) ;
  
  registerProcessorParameter("r_value",
                             " exponent in sphericity tensor use 2.0 for classical 1.0 for C,D",
			     _r,
                             float(2.0));
  registerProcessorParameter("eigenvalues_parameter_name",
                             "name of parameter to store the results ",
                             _dumpobjectname,
			     std::string("sphericity_tensor_eigenvalues"));

}


void Sphere::init() { 

  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void Sphere::processRunHeader( LCRunHeader*  /*run*/) { 

  _nRun++ ;
} 

void Sphere::processEvent( LCEvent * evt ) { 

   // this gets called for every event 
  // usually the working horse ...
 
//   float Pduzi; 
 double Pduzi; 
  float sp[3][3];
  float norm=0.0;
//   float *pp;
  const double *pp ;
//   float dvojka=2.0;
  double dvojka=2.0;

  ReconstructedParticle* p;
  LCCollection* col ; 
  
  try{ col = evt -> getCollection( _colName ) ; }
  catch(EVENT::DataNotAvailableException){
    streamlog_out(DEBUG)  << "Cannot find PFO Collection in event/run  " << evt->getEventNumber() <<" / "<< evt->getRunNumber() <<std::endl;
   streamlog_out(DEBUG) << "Skipping this event!" << std::endl;
   throw marlin::SkipEventException(this);
  }
  
  if (col->getNumberOfElements() == 0){
    streamlog_out(DEBUG) << "PFO Collection is empty in event/run  " << evt->getEventNumber() <<" / "<< evt->getRunNumber() <<std::endl;
    return;
  }



     for(int i=0;i<3;i++)
       {
    for (int j=0;j<3;j++)       
       {     
         sp[j][i]=0.0;
       }
       }

  if( col != 0 ){
    
    int nRecP = col->getNumberOfElements()  ;
    
    for(int kk=0; kk< nRecP ; kk++){
      p=dynamic_cast<ReconstructedParticle*>( col->getElementAt( kk ) );

//       pp=const_cast<float*>(p->getMomentum());
      pp = p->getMomentum() ;

     Pduzi=sqrt(pow(pp[0],dvojka)+pow(pp[1],dvojka)+pow(pp[2],dvojka));

     norm=norm+pow(Pduzi,(double) _r); 
        for(int j=0;j<3;j++)      
	   {
	for(int i=0;i<3;i++)     
	   {
	  sp[j][i]=sp[j][i]+pp[i]*pp[j]*pow(Pduzi,(_r-dvojka));
           }	
	   }
         
    } 
  }

      for(int j=0;j<3;j++)      
	   {
	for(int i=0;i<3;i++)     
	   {
	     sp[j][i]=sp[j][i]/norm;
	     //  cout << "sp tenzor "<< sp[j][i] << endl;
           }	
	   }
        float lre[3];
	Eigenvalue test(sp);
	test.getRealEigenvalues(lre);
	// cout << lre[0] << "   " << lre[1] << "  " << lre[2] << endl;  
        float sphericity;
        float aplanarity;
        string out_name;
         sphericity=1.5*(lre[0]+lre[1]);
         aplanarity=1.5*lre[0];
         FloatVec udri_me_do_zore;
         udri_me_do_zore.push_back(lre[0]);
	 udri_me_do_zore.push_back(lre[1]);
	 udri_me_do_zore.push_back(lre[2]);
	 if (_r==2.0){
         col->parameters().setValue("sphericity",sphericity);
	 col->parameters().setValue("aplanarity",aplanarity);
	 col->parameters().setValues(_dumpobjectname.c_str(),udri_me_do_zore);
         }
	 if (_r==1.0){
             float ccc;
             float ddd;
	     ccc=3.0*(lre[0]*lre[1]+lre[0]*lre[2]+lre[2]*lre[1]);
	      ddd=27.0*lre[0]*lre[1]*lre[2];
         col->parameters().setValue("C",ccc);
	 col->parameters().setValue("D",ddd);
	 col->parameters().setValues(_dumpobjectname.c_str(),udri_me_do_zore);
         }
         else 
         {
	     col->parameters().setValues(_dumpobjectname.c_str(),udri_me_do_zore);
         }
  _nEvt ++ ;
}



void Sphere::check( LCEvent *  /*evt*/ ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void Sphere::end(){ 
  
//   std::cout << "MyProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

