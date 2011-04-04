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

   registerProcessorParameter( "kappa" ,
                              "kappa for charge weighting"  ,
                              _kappa ,
			       (double) 0.4) ;

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
  taumatchtuple=new TNtuple("taumatch","taumatch","EvID:mcEfull:mcE:mcpt:mcnQ:mcnN:recE:recpt:Eseed:recnQ:recnN:D1:D2:D3:D4");
  mctuple=new TNtuple("mctuple","mctuple","EvID:E:Evis:ptvis:nQ:nN:D1:D2:D3:D4:Qw:rec");
  faketuple =new TNtuple("fake","fake","EvID:parentpdg:D1:D2:recE:recpt:recD0");
  Qtuple =new TNtuple("Q","Q","EvID:Q:nTr:Qw:trueTau");
  
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
  bool print=false;

  if( colMC != 0 ) 
    {
      int nMCP = colMC->getNumberOfElements();
      for(int k=0; k < nMCP; k++) 
	{
	  MCParticle *particle = dynamic_cast<MCParticle*> (colMC->getElementAt(k) );
	  if(particle->getGeneratorStatus()!=3 && fabs(particle->getPDG())==15)
	    {
	      
	      ntau_mc++;
	      double Evis=0,ptvis=0,Qw=0,pQ=0;
	      int nQ=0,nN=0;
	      LoopDaughters(particle,Evis,ptvis,nQ,nN,Qw,pQ);
	      int D[4]={0,0,0,0};
	      int nvisD=0;
	      for(unsigned int d=0;d< particle->getDaughters().size();d++)
		{
		  if(particle->getDaughters()[d]->getGeneratorStatus()==1 || particle->getDaughters()[d]->getPDG()==111)
		    {
		      int pdg=particle->getDaughters()[d]->getPDG();
		      if(abs(pdg)!=16 && abs(pdg)!=12 && abs(pdg)!=14)
			{
			  D[nvisD]=pdg;
			  nvisD++;
			}
		    }
		  else
		    {
		      MCParticle *daughter=particle->getDaughters()[d];
		      for(unsigned int d1=0;d1< daughter->getDaughters().size();d1++)
			{
			  if(daughter->getDaughters()[d1]->getGeneratorStatus()==1)
			    {
			      int pdg=daughter->getDaughters()[d1]->getPDG();
			      if(abs(pdg)!=16 && abs(pdg)!=12 && abs(pdg)!=14)
				{
				  D[nvisD]=pdg;
				  nvisD++;
				}
			    }
			}
		    }
		}
	      //Try to find a corresponding reconstructed tau
	      bool found= false;
	      if( colTau != 0)
		{
		  int nT = colTau->getNumberOfElements();
		  for(int k=0; k < nT; k++) 
		    {
		      ReconstructedParticle *tau = dynamic_cast <ReconstructedParticle*>( colTau->getElementAt( k ) );
		      std::vector< ReconstructedParticle * > tauvec=tau->getParticles();
		      for(unsigned int o=0;o<tauvec.size();o++)
			{
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
				      // cout<<parent<<" "<<particle<<endl;
				      if(parent==particle)
					found=true;
				    }
				}
			    }
			}
		    }
		}
	      
	      if(_nEvt<coutUpToEv || _nEvt==coutEv || print)
		{
		  if(found)
		    cout<<"MCTRUTH: "<<found<<" "<<particle->getPDG()<<" "<<Evis<<" "<<ptvis<<" "<<nQ<<" "<<nN<<
		      " "<<atan(sqrt(particle->getMomentum()[1]*particle->getMomentum()[1]+particle->getMomentum()[0]*particle->getMomentum()[0])/particle->getMomentum()[2])<<" "<<D[0]<<" "<<D[1]<<" "<<D[2]<<" "<<D[3]<<endl;
		  else
		    cout<<"MCMISSED: "<<found<<" "<<particle->getPDG()<<" "<<Evis<<" "<<ptvis<<" "<<nQ<<" "<<nN<<
		      " "<<atan(sqrt(particle->getMomentum()[1]*particle->getMomentum()[1]+particle->getMomentum()[0]*particle->getMomentum()[0])/particle->getMomentum()[2])<<" "<<D[0]<<" "<<D[1]<<" "<<D[2]<<" "<<D[3]<<endl;
		}
	      mctuple->Fill(_nEvt,particle->getEnergy(),Evis,ptvis,nQ,nN,D[0],D[1],D[2],D[3],Qw/pQ,found);
	      if(found==false)
		print=true;
	    }//tau
	}
    }

  

  if( colTau != 0)
    {
      int nT = colTau->getNumberOfElements();
      ntau_rec=nT;
      
      if(_nEvt<coutUpToEv || _nEvt==coutEv || print)
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
	  double Qweighted=0;
	  double ptsum=0;
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
		  double pt=sqrt(tauvec[o]->getMomentum()[0]*tauvec[o]->getMomentum()[0]
				 +tauvec[o]->getMomentum()[1]*tauvec[o]->getMomentum()[1]);
		  Qweighted+=tauvec[o]->getCharge()*pow(tauvec[o]->getEnergy(),_kappa);
		  ptsum+=pow(tauvec[o]->getEnergy(),_kappa);
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
	  float trueTau=0;
	  if(mctau)
	    trueTau=1;
	  Qtuple->Fill(_nEvt,tau->getCharge(),nQ,Qweighted/ptsum,trueTau);
	   
	  //compare tau with mc truth
	  if(mctau)
	    {
	      double Evis=0,ptvis=0,Qw=0,pQ=0;
	      int mcnQ=0,mcnN=0;
	      
	      LoopDaughters(mctau,Evis,ptvis,mcnQ,mcnN,Qw,pQ);
	      int D[4]={0,0,0,0};
	      int nvisD=0;
	      for(unsigned int d=0;d< mctau->getDaughters().size();d++)
		{
		  if(mctau->getDaughters()[d]->getGeneratorStatus()==1 || mctau->getDaughters()[d]->getPDG()==111)
		    {
		      int pdg=mctau->getDaughters()[d]->getPDG();
		      if(abs(pdg)!=16 && abs(pdg)!=12 && abs(pdg)!=14)
			{
			  D[nvisD]=pdg;
			  nvisD++;
			}
		    }
		  else
		    {
		      MCParticle *daughter=mctau->getDaughters()[d];
		      for(unsigned int d1=0;d1< daughter->getDaughters().size();d1++)
			{
			  if(daughter->getDaughters()[d1]->getGeneratorStatus()==1)
			    {
			      int pdg=daughter->getDaughters()[d1]->getPDG();
			      if(abs(pdg)!=16 && abs(pdg)!=12 && abs(pdg)!=14)
				{
				  D[nvisD]=pdg;
				  nvisD++;
				}
			    }
			}
		    }
		}
	      taumatchtuple->Fill(_nEvt,mctau->getEnergy(),Evis,ptvis,nQ,nN,tau->getEnergy(),pt,Eseed,nQ,nN,D[0],D[1],D[2],D[3]);
	      if(_nEvt<coutUpToEv || _nEvt==coutEv || print)
		cout<<"REC MATCHED: "<<tau->getEnergy()<<" "<<nQ<<" "<<nN<<" "<<D0<<endl;
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
	      if(_nEvt<coutUpToEv || _nEvt==coutEv || print)
		cout<<"REC FAKE: "<<tau->getEnergy()<<" "<<nQ<<" "<<nN<<" "<<D0<<endl;
	      isfake=true;
	      nfakes++;
	    }
	}
    }//colTau
  
 
  
  //filling the tuple
  evtuple->Fill(_nEvt,ntau_mc,ntau_rec);
   
  //cleanup
  delete relationNavigatorMC;
  delete relationNavigatorPFOMC;
  
  if(_nEvt<coutUpToEv || _nEvt==coutEv || print)
    cout<<"--------------------------------------------------------------------------------------------"<<endl;
   
  _nEvt ++ ;
   
}

void EvaluateTauID::LoopDaughters(MCParticle *particle,double &Evis,double &ptvis,int &nQ, int &nN, double &Qw, double &pQ)
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
	      if(daughter->getCharge()==1 || daughter->getCharge()==-1)
		{ 
		  //Qw+=daughter->getCharge()*pow(sqrt(mc_pvec[0]*mc_pvec[0]+mc_pvec[1]*mc_pvec[1]),_kappa);
		  Qw+=daughter->getCharge()*pow(daughter->getEnergy(),_kappa);
		  //pQ+=pow(sqrt(mc_pvec[0]*mc_pvec[0]+mc_pvec[1]*mc_pvec[1]),_kappa);
		  pQ+=pow(daughter->getEnergy(),_kappa);
		  nQ++;
		}
	      else
		nN++;
	    }
	}
      if(daughter->getDaughters().size())
	LoopDaughters(daughter,Evis,ptvis,nQ,nN,Qw,pQ);
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
  Qtuple->Write(); 
  //Close File here
  rootfile->Write();
  rootfile->Close();
}


