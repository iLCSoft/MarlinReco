#include "PrepareSmearRECParticles.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;


#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCRelation.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/TrackImpl.h>
#include "IMPL/LCFlagImpl.h" 
#include "UTIL/LCRelationNavigator.h"

#include <gear/GEAR.h>
#include <gear/BField.h>
#include <marlin/Global.h>
#include "HelixClass.h"

#include <gsl/gsl_randist.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#define coutEv -1
#define coutUpToEv 10

using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

PrepareSmearRECParticles aPrepareSmearRECParticles ;

template<typename T>
T sgn(T n)
{
if (n < 0) return -1;
if (n > 0) return 1;
return 0;
}

PrepareSmearRECParticles::PrepareSmearRECParticles() : Processor("PrepareSmearRECParticles") 
{
  // modify processor description
  _description = "PrepareSmearRECParticles converts MCParticles to ReconstructedParticles  and smeares impact parameter resolution and puts them into a new collection making sure all the information which is needed to run the TauFinder is provided. ";

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCCollectionName" , 
			   "Name of the MCParticle collection"  ,
			   _colNameMC ,
			   std::string("MCParticle") ) ;
                                        
  registerProcessorParameter( "outputColMC" ,
                              "Name of the output Collection of refilled information"  ,
                              _outcolMC ,
                              std::string("MCParticles_tau")) ;
   
  registerProcessorParameter( "D0res_a" ,
                              "Impact parameter resolution value for a in microns"  ,
                              _D0res_a ,
                              (double)2.) ;

  registerProcessorParameter( "D0res_b" ,
                              "Impact parameter resolution value for b in microns"  ,
                              _D0res_b ,
                              (double)10.) ;
  
  registerProcessorParameter( "momres" ,
                              "Momentum resolution delta(1/pt)"  ,
                              _momres ,
                              (double)5e-5) ;

  registerProcessorParameter( "Eres" ,
                              "Energy resolution in calorimeter"  ,
                              _Eres ,
                              (double)0.6) ;

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RecCollection",
                            "Collection of Rec Particles for TauFinder",
                            _outcolMC ,
                            std::string("MCParticles_tau"));
    
  registerOutputCollection( LCIO::LCRELATION,
			   "MCRECLinkCollectionName" , 
			   "Name of the MC Truth ReconstructedParticle collection"  ,
			   _colNameMCTruth ,
			   std::string("MCRecLink") ) ;

 registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "outputColTracks",
                            "Collection of Rec Particles for TauFinder",
                            _outcolTracks ,
                            std::string("Tracks_tau"));
}


void PrepareSmearRECParticles::init() 
{ 
  std::cout << "INIT CALLED" << std::endl;
  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  // usually a good idea to
  printParameters() ;
  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
  _nRun = 0 ;
  _nEvt = 0 ;
  //intialise random number generator 
  _random = gsl_rng_alloc(gsl_rng_ranlxs2);

  std::cout << "INIT IS DONE" << std::endl;
}

void PrepareSmearRECParticles::processRunHeader( LCRunHeader* run) 
{ 
  _nRun++ ;
} 

void PrepareSmearRECParticles::processEvent( LCEvent * evt ) 
{ 
  // this gets called for every event 
  // usually the working horse ...
  
  LCCollection *colMC;
  try {
    colMC = evt->getCollection( _colNameMC ) ;
  } catch (Exception e) {
    colMC = 0;
  }

  _nEvt = evt->getEventNumber();  
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"EVENT "<<_nEvt<<endl;
  
  LCCollectionVec *reccol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec *trackcol = new LCCollectionVec(LCIO::TRACK);
  //LCRelation stuff
  LCCollectionVec *mc_relationcol = new LCCollectionVec(LCIO::LCRELATION);
  mc_relationcol->parameters().setValue(std::string("FromType"),LCIO::RECONSTRUCTEDPARTICLE);
  mc_relationcol->parameters().setValue(std::string("ToType"),LCIO::MCPARTICLE);
  
  HelixClass *helix = new HelixClass();
  //convert MCPARTICLES 
  if( colMC != 0 ) 
    {
      int nMCP = colMC->getNumberOfElements();
      for(int k=0; k < nMCP; k++) 
	{
	  MCParticle *mc = dynamic_cast<MCParticle*> (colMC->getElementAt(k) );
	  if(mc->getGeneratorStatus()==1)//only stable ones
	    {
	       //filter out the neutrinos and other invisibles
	      if(  fabs(mc->getPDG())==12 || fabs(mc->getPDG())==14 || fabs(mc->getPDG())==16
		   || fabs(mc->getPDG())==1000022)
		continue;
	      //detector acceptance
	      double Cos_T  = fabs(mc->getMomentum()[2]) / sqrt(pow(mc->getMomentum()[0],2)+pow(mc->getMomentum()[1],2) + pow(mc->getMomentum()[2],2));
	      double p_t  = sqrt(pow(mc->getMomentum()[0],2) + pow(mc->getMomentum()[1],2));
	      
	      if (p_t   > 0.2 && Cos_T <  0.99)
		{
		  ReconstructedParticleImpl *rec = new ReconstructedParticleImpl();
		  //copy values from MC to REC
		  //smear momentum resolution, leave z untouched
		  float mommc[3];
		  for(int i=0;i<3;i++)
		    mommc[i]=mc->getMomentum()[i];
		  double pt=sqrt(mommc[0]*mommc[0]+mommc[1]*mommc[1]);
		  double phi = atan2( mommc[1],mommc[0] );
		  double sigma_pt=_momres*pt*pt;
		  double pt_smear = pt+gsl_ran_gaussian( _random,sigma_pt ) ; 
		  mommc[0] = pt_smear*cos(phi) ; 
		  mommc[1]=  pt_smear*sin(phi);
		  //cout<<pt<<" "<<sigma_pt<<" "<<sqrt(mommc[0]*mommc[0]+mommc[1]*mommc[1])<<endl;	      
		  rec->setMomentum(mommc);
		  rec->setType(mc->getPDG());
		  rec->setMass(mc->getMass());
		  rec->setCharge(mc->getCharge());
		  //add the track if charged, so that information for D0 is present
		  if(mc->getCharge())
		    {
		      
		      TrackImpl *track=new TrackImpl();
		      float ver[3];
		      float mom[3];
		      float rp[3];
		      for(int i=0;i<3;i++){
			ver[i]=mc->getVertex()[i];
			mom[i]=rec->getMomentum()[i];
		      }
		      helix->Initialize_VP(ver,mom,mc->getCharge(),_bField);
		      for(int i=0;i<3;i++)
			rp[i]=helix->getReferencePoint()[i];
		      track->setReferencePoint(rp);
		      track->setPhi (helix->getPhi0());
		      track->setOmega (helix->getOmega());
		      track->setTanLambda (helix->getTanLambda());
		 
		      double p=sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
		      //smear energy, in case of charged particle the same as momentum resolution
		      rec->setEnergy(sqrt(p*p+mc->getMass()*mc->getMass()));
		      //smearing of D0, try 2 < a < 6 micron and 10 < b < 20 micron
// 		      double pt=sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
// 		      double theta=atan(pt/fabs(mom[2]));
// 		      double sigma_D0 = sqrt(_D0res_a*_D0res_a + _D0res_b*_D0res_b/(p*p*pow(sin(theta),3))); 
// 		      double mean_D0=fabs(helix->getD0());
// 		      double D0=mean_D0+gsl_ran_gaussian(_random,sigma_D0)/1000.;
// 		      track->setD0(fabs(D0));
// 		      //cout<<mean_D0<<" "<<sigma_D0<<" "<<fabs(D0)<<" "<<pt<<"  "<<p<<endl;
// 		      double sigma_z = sqrt(_D0res_a*_D0res_a + _D0res_b*_D0res_b/(p*pow(sin(theta),5)));
// 		      double mean_z=helix->getZ0();
// 		      track->setZ0(mean_z+gsl_ran_gaussian(_random,sigma_z)/1000);
		      //no D0 smearing
		      track->setD0(fabs(helix->getD0()));
		      track->setZ0 (helix->getZ0());
		      rec->addTrack(track);
		      trackcol->addElement(track);
		    
		    }
		  else
		    {
		      //smear energy, in case of neutral particle given by calorimter
		      double sigma_E=_Eres*mc->getEnergy()/sqrt(mc->getEnergy());
		      double Esmeared=mc->getEnergy()+gsl_ran_gaussian(_random,sigma_E);
		      if(Esmeared<0)
			Esmeared=0;
		      rec->setEnergy(Esmeared);
		    }
		
		  //cout<<mc->getCharge()<<" "<<mc->getEnergy()<<" "<<rec->getEnergy()<<endl;
		  reccol->addElement( rec );
		  LCRelationImpl *rel = new LCRelationImpl(rec,mc);
		  mc_relationcol->addElement( rel );
		}
	    }
	}
    }
  delete helix;
  
  
  evt->addCollection(reccol,_outcolMC);
  evt->addCollection(trackcol,_outcolTracks);
  evt->addCollection(mc_relationcol,_colNameMCTruth);
  
   
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"--------------------------------------------------------------------------------------------"<<endl;
   
  _nEvt ++ ;
  
}



void PrepareSmearRECParticles::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void PrepareSmearRECParticles::end(){ 
  
   gsl_rng_free(_random);

  std::cout << "PrepareSmearRECParticles::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    << std::endl ;
  
  //Close File here
  

}


