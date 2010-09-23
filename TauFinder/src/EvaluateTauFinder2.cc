#include "EvaluateTauFinder2.h"
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

#define coutEv 0
#define coutUpToEv 10

using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

EvaluateTauFinder2 aEvaluateTauFinder2 ;



EvaluateTauFinder2::EvaluateTauFinder2() : Processor("EvaluateTauFinder2") 
{
  // modify processor description
  _description = "EvaluateTauFinder2 checks performance of TauFinder and writes output to root file." ;
  
  // register steering parameters: name, description, class-variable, default value
  
  registerProcessorParameter( "inputCol" ,
                              "Name of the input Collection"  ,
                              _incol ,
                              std::string("TauRec_MC")) ;
  
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCCollectionName           " , 
			   "Name of the MCParticle collection"  ,
			   _colNameMC ,
			   std::string("MCParticlesSkimmed") ) ;
  
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "TauRecCollection",
			   "Collection of Tau Candidates",
			   _incol ,
			   std::string("TauRec_PFA"));
  
  
  registerInputCollection( LCIO::LCRELATION,
			   "MCTauLinkCollectionName" , 
			   "Name of the MC Truth Reconstructed Tau collection"  ,
			   _colNameMCRecLink ,
			   std::string("MCRecLink") ) ;
   
}


void EvaluateTauFinder2::init() 
{ 
  std::cout << "INIT CALLED" << std::endl;
  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  _ntau_correct= 0;
  _dEsum= 0;
  _dEsumsq= 0;
  _ndE= 0;
  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
 
  evtuple=new TNtuple("evtuple","evtuple","EvID:Ntaus_mc:Ntaus_rec:missed");
  tautuple=new TNtuple("tautuple","tautuple","EvID:recE:recp:recpt:recD0:recQ:recN:Q");
  mctautuple=new TNtuple("mctautuple","mctautuple","EvID:Evis:ptvis:pvis:mcQ:mcN");
  mcmisstuple=new TNtuple("mcmiss","mcmiss","EvID:E:D0:D1:D2");
  taumatchtuple=new TNtuple("taumatch","taumatch","EvID:recE:recp:recpt:recD0:recQ:recN:mcE:mcp:mcpt:mcQ:mcN:mcEfull");
  faketuple =new TNtuple("fake","fake","EvID:parentpdg:D1:D2:recE:recp:recD0");
  
  
  std::cout << "INIT IS DONE" << std::endl;
}

void EvaluateTauFinder2::processRunHeader( LCRunHeader* run) 
{ 
  _nRun++ ;
} 

void EvaluateTauFinder2::processEvent( LCEvent * evt ) 
{ 
  // this gets called for every event 
  // usually the working horse ...
  
  LCCollection *colMC, *colTau;
  LCCollection *colMCRecLink;
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
    colMCRecLink = evt->getCollection( _colNameMCRecLink ) ;
  } catch (Exception e) {
    colMCRecLink = 0;
  }
  
  LCRelationNavigator* relationNavigatorMC = 0;
  if( colMCRecLink != 0)
    relationNavigatorMC = new LCRelationNavigator( colMCRecLink );
  
  _nEvt = evt->getEventNumber();  
  
  int ntau_mc=0;
  int ntau_rec=0;
  int missed=0;

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
	  double p=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]+pvec[2]*pvec[2]);
	  //double phi=180./TMath::Pi()*atan(pvec[1]/pvec[0]);
	  //double theta=180./TMath::Pi()*atan(pt/fabs(pvec[2])); 
	  //const ReconstructedParticleVec &tauvec=tau->getParticles();
	  std::vector< ReconstructedParticle * > tauvec=tau->getParticles();
	  double Eseed=0,D0=0;
	  int Q=0,N=0;
	  bool istau=false;
	  MCParticle *mctau=NULL;
	  for(unsigned int o=0;o<tauvec.size();o++)
	    {
	      //find seed track for D0
	      if(tauvec[o]->getCharge()!=0)
		{
		  Q++;
		  if(tauvec[o]->getEnergy()>Eseed && tauvec[o]->getTracks().size()!=0 )
		    {
		      D0=(float)tauvec[o]->getTracks()[0]->getD0();
		      Eseed=tauvec[o]->getEnergy();
		    }
		}
	      else
		N++;
	      
	      if(relationNavigatorMC)
		{		
		  EVENT::LCObjectVec relobj = relationNavigatorMC->getRelatedToObjects(tauvec[o]);
		  for(unsigned int m=0;m<relobj.size();m++)
		    {
		      MCParticle *mc=dynamic_cast <MCParticle*>(relobj[m]);
		      MCParticle *dummy=mc;
		      MCParticle *parent=mc;
		      int size=mc->getParents().size();
		      while(size!=0)
			{
			  dummy=parent->getParents()[0];
			  size=dummy->getParents().size();
			  parent=dummy;
			  if(fabs(parent->getPDG())==15)
			    {
			      size=0;
			      istau=true;
			      mctau=parent;
			    }
			}
		    }
		}
	    }
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    cout<<"Tau: "<<k<<" "<<tau->getEnergy()<<" "<<p<<" "<<D0
		<<" "<<Q<<" "<<N<<" "<<tauvec.size()<<endl;
	  if(Q>4)
	    {
	      cout<<_nEvt<<" "<<tauvec.size()<<" "<<Q<<" "<<N<<endl;
	      for(unsigned int o=0;o<tauvec.size();o++)
		cout<<k<<" "<<tauvec[o]->getCharge()<<" "<<tauvec[o]->getType()<<" "<<tauvec[o]->getEnergy()<<endl;
	    }

	  tautuple->Fill(_nEvt,tau->getEnergy(),p,pt,D0,Q,N,tau->getCharge());
	  if(mctau)
	    {
	      double Evis=0,ptvis=0,pvis=0;
	      int ntr=0,qtr=0;
	      MCParticle *highEtrack=NULL;
	      LoopDaughters(highEtrack,mctau,Evis,ptvis,pvis,ntr,qtr);
	      taumatchtuple->Fill(_nEvt,tau->getEnergy(),p,pt,D0,Q,N,Evis,pvis,ptvis,qtr,ntr,mctau->getEnergy());
	    }
	}
    }
  
  if( colMC != 0 ) 
    {
      int nMCP = colMC->getNumberOfElements();
      if(_nEvt<coutUpToEv || _nEvt==coutEv)
	cout<<"MCTRUTH: "<<endl;
      HelixClass * helix = new HelixClass();
      for(int k=0; k < nMCP; k++) 
	{
	  MCParticle *particle = dynamic_cast<MCParticle*> (colMC->getElementAt(k) );
	  if(particle->getGeneratorStatus()!=3 && fabs(particle->getPDG())==15)
	    {
	      ntau_mc++;
	      double Evis=0,ptvis=0,pvis=0,D0=0;
	      int ntr=0,qtr=0;
	      MCParticle *highEtrack=NULL;
	      LoopDaughters(highEtrack,particle,Evis,ptvis,pvis,ntr,qtr);
	      mctautuple->Fill(_nEvt,Evis,ptvis,pvis,qtr,ntr);
	      //find out which mc taus do not have a link to the rec
	      if(relationNavigatorMC)
		{
		  bool hasRel=false;
		  LoopDaughtersRelation(particle,relationNavigatorMC ,hasRel);
		  if(!hasRel)
		    {		    
		      missed++;
		      mcmisstuple->Fill(_nEvt,Evis,D0,particle->getDaughters()[0]->getPDG(),particle->getDaughters()[1]->getPDG());
		      if(_nEvt<coutUpToEv || _nEvt==coutEv)
			cout<<"Missed: "<<Evis<<" "<<D0<<" "<<particle->getDaughters()[0]->getPDG()<<" "<<particle->getDaughters()[1]->getPDG()<<endl;
		    }
		}
	      if(_nEvt<coutUpToEv || _nEvt==coutEv)
		cout<<"MC: "<<Evis<<" "<<qtr<<" "<<ntr<<endl;
	      
	    }//tau
	 	   
	}
      delete helix;
    }
  
  //filling the tuple
  evtuple->Fill(_nEvt,ntau_mc,ntau_rec,missed);
 
   if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"--------------------------------------------------------------------------------------------"<<endl;
  
  
  _nEvt ++ ;
  //cleanup
  
  delete relationNavigatorMC;
    
}

void EvaluateTauFinder2::LoopDaughters(MCParticle *hEt,MCParticle *particle,double &Evis,
				       double &ptvis,double &pvis,int &ntr,int &qtr)
{
  double trackE=0;
  for(unsigned int d=0;d<particle->getDaughters().size();d++)
    { 
      MCParticle *daughter=particle->getDaughters()[d];
      //only particles visible in the detector
      if(daughter->hasLeftDetector()==false || fabs(daughter->getPDG())==13
	 && fabs(daughter->getPDG())!=1000022)
	{
	  if (daughter->getGeneratorStatus()==1)
	    {
	      if(daughter->getCharge() && daughter->getEnergy()>trackE)
		{
		  trackE=daughter->getEnergy();
		  hEt=daughter;
		}
	      Evis+=daughter->getEnergy();
	      const double *mc_pvec=daughter->getMomentum();
	      ptvis+=sqrt(mc_pvec[0]*mc_pvec[0]+mc_pvec[1]*mc_pvec[1]);
	      pvis+=sqrt(mc_pvec[0]*mc_pvec[0]+mc_pvec[1]*mc_pvec[1]+mc_pvec[2]*mc_pvec[2]);
	      if(daughter->getCharge()!=0)
		qtr++;
	      else
		ntr++;
	    }
	}
      if(daughter->getDaughters().size())
	LoopDaughters(hEt,daughter,Evis,ptvis,pvis,ntr,qtr);
    }
}
 
void EvaluateTauFinder2::LoopDaughtersRelation(MCParticle *particle, LCRelationNavigator* relationNavigatorMC ,bool &relToTau)
{
  for(unsigned int d=0;d<particle->getDaughters().size();d++)
    { 
      MCParticle *daughter=particle->getDaughters()[d];
      //only particles visible in the detector
      if(daughter->hasLeftDetector()==false || fabs(daughter->getPDG())==13)
	{
	  if (daughter->getGeneratorStatus()==1)
	    {
	      //relation to the filled reconstructed particle
	      EVENT::LCObjectVec relobjTO = relationNavigatorMC->getRelatedFromObjects(daughter);
	      if(relobjTO.size())
		relToTau=true;
	      
	    }
	}
      if(relToTau)
	break;
      if(daughter->getDaughters().size())
	LoopDaughtersRelation(daughter,relationNavigatorMC,relToTau);
    }
}


void EvaluateTauFinder2::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EvaluateTauFinder2::end(){ 
  

  std::string end(".root");
  std::string fname = _incol+end;
  const char *filename=fname.c_str();
  rootfile = new TFile(filename,"RECREATE");
  
  std::cout << "EvaluateTauFinder2::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "<<std::endl;
  

  evtuple->Write();
  tautuple->Write();
  mctautuple->Write();
  mcmisstuple->Write();
  taumatchtuple->Write();
  faketuple->Write();
  
  //Close File here
  rootfile->Write();
  rootfile->Close();
}


