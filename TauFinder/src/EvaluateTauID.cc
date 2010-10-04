#include "EvaluateTauID.h"
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
#include <EVENT/LCObject.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/VertexImpl.h>
#include "UTIL/LCRelationNavigator.h"

#include <gear/GEAR.h>
#include <gear/BField.h>
#include <marlin/Global.h>

#include "HelixClass.h"
#include "SimpleLine.h"
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#define coutEv -1
#define coutUpToEv 100

using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

EvaluateTauID aEvaluateTauID ;

EvaluateTauID::EvaluateTauID() : Processor("EvaluateTauID") 
{
  // modify processor description
  _description = "EvaluateTauID checks performance of TauID and writes output to root file." ;
  
  // register steering parameters: name, description, class-variable, default value
  
  registerProcessorParameter( "inputCol" ,
                              "Name of the input Collection"  ,
                              _incol ,
                              std::string("TauRec_PFA")) ;
  
  registerProcessorParameter( "relCol" ,
                              "Name of the LCRelation otput Collection"  ,
                              _colNameMCRecLink ,
                              std::string("MCRecLink")) ;

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCCollectionName           " , 
			   "Name of the MCParticle collection"  ,
			   _colNameMC ,
			   std::string("MCParticlesSkimmed") ) ;
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "RECOCollectionName" , 
			   "Name of the ReconstructedParticle collection"  ,
			   _colNameRECO ,
			   std::string("PandoraPFOs") ) ;
  
  registerInputCollection( LCIO::LCRELATION,
			   "RECOMCTRUTHCollectionName" , 
			   "Name of the MC Truth PFA collection"  ,
			   _colNameMCTruth ,
			   std::string("RecoMCTruthLink") ) ;
  
  registerInputCollection( LCIO::LCRELATION,
			   "RECOMCTRUTHCollectionName" , 
			   "Name of the MC Truth PFA collection"  ,
			   _colNameMCRecLink ,
			   std::string("MCRecLink") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "TauRecCollection",
                            "Collection of Tau Candidates",
                            _incol ,
                            std::string("TauRec_PFA"));
 
  
}


void EvaluateTauID::init() 
{ 
  std::cout << "INIT CALLED" << std::endl;
  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
 
  evtuple=new TNtuple("evtuple","evtuple","EvID:Ntaus_mc:Ntaus_rec");
  taumatchtuple=new TNtuple("taumatch","taumatch","EvID:mcEfull:mcE:mcpt:mcnQ:mcnN:recE:recpt:recD0:Eseed:recnQ:recnN");
  mctuple=new TNtuple("mcmiss","mcmiss","EvID:E:Evis:ptvis:nQ:nN:D1:D2");
  faketuple =new TNtuple("fake","fake","EvID:parentpdg:D1:D2:recE:recpt:recD0");

  
  std::cout << "INIT IS DONE" << std::endl;
}

void EvaluateTauID::processRunHeader( LCRunHeader* run) 
{ 
  _nRun++ ;
} 

void EvaluateTauID::processEvent( LCEvent * evt ) 
{ 
  // this gets called for every event 
  // usually the working horse ...
  
  LCCollection *colMC,*colMCTruth, *colTau,*colMCRecLink;

  try {
    colMC = evt->getCollection( _colNameMC ) ;
  } catch (Exception e) {
    colMC = 0;
  }
 
  try {
    colTau = evt->getCollection( _incol ) ;
  } catch (Exception e) {
    colTau = 0;
  }
  
  try {
    colMCTruth = evt->getCollection( _colNameMCTruth ) ;
  } catch (Exception e) {
    colMCTruth = 0;
  }
 
   try {
    colMCRecLink = evt->getCollection( _colNameMCRecLink ) ;
  } catch (Exception e) {
    colMCRecLink = 0;
  }
 
  _nEvt = evt->getEventNumber();  
  
  int ntau_mc=0;
  int ntau_rec=0;
 
  LCRelationNavigator* relationNavigatorPFOMC = 0;
  LCRelationNavigator* relationNavigatorMC = 0;
  if( colMCTruth != 0)
    relationNavigatorPFOMC = new LCRelationNavigator( colMCTruth );
  if( colMCRecLink != 0)
    relationNavigatorMC = new LCRelationNavigator( colMCRecLink );
 
  bool isfake=false;
  int nfakes=0;

  if( colTau != 0)
    {
      int nT = colTau->getNumberOfElements();
      ntau_rec=nT;
      
      if(_nEvt<coutUpToEv || _nEvt==coutEv)
	cout<<"EVENT "<<_nEvt<<" with "<<nT<<" taus"<<endl;
     
      for(int k=0; k < nT; k++) 
	{
	  ReconstructedParticle *tau = dynamic_cast <ReconstructedParticle*>( colTau->getElementAt( k ) );
	  const double *pvec=tau->getMomentum();
	  double pt=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
	  std::vector< ReconstructedParticle * > tauvec=tau->getParticles();
	  double Eseed=0, D0=0;
	  int nQ=0,nN=0;
	  bool istau=false;
	  bool contaminated=false;
	  MCParticle *mctau=NULL;
	  //loop over content of tau
	  for(unsigned int o=0;o<tauvec.size();o++)
	    {
	      //find seed track for D0
	      if(tauvec[o]->getCharge()!=0)
		{
		  nQ++;
		  if(tauvec[o]->getEnergy()>Eseed && tauvec[o]->getTracks().size()!=0 )
		    {
		      D0=(float)tauvec[o]->getTracks()[0]->getD0();
		      Eseed=tauvec[o]->getEnergy();
		    }
		}
	      else
		nN++;
	
	      //follow the chain back to mc truth
	      //for filled MC particles
	    
	      if(relationNavigatorMC)
		{
		  EVENT::LCObjectVec relobj = relationNavigatorMC->getRelatedToObjects(tauvec[o]);
		  for(unsigned int m=0;m<relobj.size();m++)
		    {
		      MCParticle *mc=dynamic_cast <MCParticle*>(relobj[m]);
		      //check whether particles parent is really a tau:
		      MCParticle *dummy=mc;
		      MCParticle *parent=mc;
		      int size=mc->getParents().size();
		      while(size!=0)
			{
			  dummy=parent->getParents()[0];
			  size=dummy->getParents().size();
			  parent=dummy;
			  if(fabs(parent->getPDG())==15)
			    size=0;
			}
		      if(fabs(parent->getPDG())==15)
			{
			  istau=true;
			  mctau=parent;
			}
		      else
			contaminated=true;
		    }
		  
		}
	      
	      //for PFA objects
	      if(relationNavigatorPFOMC)
		{
		  EVENT::LCObjectVec relobjMC = relationNavigatorPFOMC->getRelatedToObjects(tauvec[o]);
		  for(unsigned int m=0;m<relobjMC.size();m++)
		    {
		      MCParticle *mc=dynamic_cast <MCParticle*>(relobjMC[m]);
		      //check whether particles parent is really a tau:
		      MCParticle *dummy=mc;
		      MCParticle *parent=mc;
		      int size=mc->getParents().size();
		      while(size!=0)
			{
			  dummy=parent->getParents()[0];
			  size=dummy->getParents().size();
			  parent=dummy;
			  if(fabs(parent->getPDG())==15)
			    size=0;
			}
		      if(fabs(parent->getPDG())==15)
			{
			  istau=true;
			  mctau=parent;
			}
		      else
			contaminated=true;
		    }
		}
	      
	    }//tau content
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    cout<<"REC: "<<tau->getEnergy()<<" "<<nQ<<" "<<nN<<" "<<D0<<endl;
	  
	  //compare tau with mc truth
	  if(mctau)
	    {
	      double Evis=0,ptvis=0;
	      int mcnQ=0,mcnN=0;
	      
	      LoopDaughters(mctau,Evis,ptvis,mcnQ,mcnN);

	      taumatchtuple->Fill(_nEvt,mctau->getEnergy(),Evis,ptvis,nQ,nN,tau->getEnergy(),pt,D0,Eseed,nQ,nN);
	    }
	  else //fake
	    { 
	      int d1=0,d2=0,pdg=0;
	      for(unsigned int o=0;o<tauvec.size();o++)
		{
		  if(relationNavigatorMC)
		    {
		      EVENT::LCObjectVec relobj = relationNavigatorMC->getRelatedToObjects(tauvec[o]);
		      for(unsigned int m=0;m<relobj.size();m++)
			{
			  MCParticle *mc=dynamic_cast <MCParticle*>(relobj[m]);
			  if(mc->getPDG()==22)
			    continue;
			  MCParticle *dummy=mc;
			  MCParticle *parent=mc;
			  int size=mc->getParents().size();
			  while(size!=0)
			    {
			      dummy=parent->getParents()[0];
			      size=dummy->getParents().size();
			      parent=dummy;
			    }
			  pdg=parent->getPDG();
			  if(parent->getDaughters().size())
			    d1=parent->getDaughters()[0]->getPDG();
			  if(parent->getDaughters().size()>1)
			    d2=parent->getDaughters()[1]->getPDG();
			}
		    }
		}
	      faketuple->Fill(_nEvt,pdg,d1,d2,tau->getEnergy(),pt,D0);
	      isfake=true;
	      nfakes++;
	    }
	}
    }//colTau
  
  if( colMC != 0 ) 
    {
      int nMCP = colMC->getNumberOfElements();
      if(_nEvt<coutUpToEv || _nEvt==coutEv)
	cout<<"MCTRUTH: "<<endl;
      for(int k=0; k < nMCP; k++) 
	{
	  MCParticle *particle = dynamic_cast<MCParticle*> (colMC->getElementAt(k) );
	  if(particle->getGeneratorStatus()!=3 && fabs(particle->getPDG())==15)
	    {
	      ntau_mc++;
	      double Evis=0,ptvis=0;
	      int nQ=0,nN=0;
	      LoopDaughters(particle,Evis,ptvis,nQ,nN);
	    
	      if(_nEvt<coutUpToEv || _nEvt==coutEv)
		cout<<particle->getPDG()<<" "<<Evis<<" "<<ptvis<<" "<<nQ<<" "<<nN<<endl;
	      
	      mctuple->Fill(_nEvt,particle->getEnergy(),Evis,ptvis,nQ,nN,particle->getDaughters()[0]->getPDG(),particle->getDaughters()[1]->getPDG());
	    	      
	    }//tau
	}
    }
  
  //filling the tuple
  evtuple->Fill(_nEvt,ntau_mc,ntau_rec);
   
  //cleanup
  delete relationNavigatorMC;
  delete relationNavigatorPFOMC;
  
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"--------------------------------------------------------------------------------------------"<<endl;
   
  _nEvt ++ ;
   
}

void EvaluateTauID::LoopDaughters(MCParticle *particle,double &Evis,double &ptvis,int &nQ, int &nN)
{
  for(unsigned int d=0;d<particle->getDaughters().size();d++)
    { 
      MCParticle *daughter=particle->getDaughters()[d];
      //only particles visible in the detector
      if(daughter->hasLeftDetector()==false || fabs(daughter->getPDG())==13)
	{
	  if (daughter->getGeneratorStatus()==1)
	    {
	      Evis+=daughter->getEnergy();
	      const double *mc_pvec=daughter->getMomentum();
	      ptvis+=sqrt(mc_pvec[0]*mc_pvec[0]+mc_pvec[1]*mc_pvec[1]);
	    }
	  if(daughter->getCharge()==1 || daughter->getCharge()==-1)
	    nQ++;
	  else
	    nN++;
	}
      if(daughter->getDaughters().size())
	LoopDaughters(daughter,Evis,ptvis,nQ,nN);
    }
}
 

void EvaluateTauID::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EvaluateTauID::end(){ 
  

  std::string end(".root");
  std::string fname = _incol+end;
  const char *filename=fname.c_str();
  rootfile = new TFile(filename,"RECREATE");
  
  std::cout << "EvaluateTauID::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "<<std::endl;
  

  evtuple->Write();
  mctuple->Write();
  taumatchtuple->Write();
  faketuple->Write();
    
  //Close File here
  rootfile->Write();
  rootfile->Close();
}


