#include "EvaluateTauFinder.h"
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

#define coutUpToEv 10


using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

EvaluateTauFinder aEvaluateTauFinder ;

struct TAU {   // Declare  struct type
  double E,phi,theta,D0,NQ,N;   // Declare member types
};   

bool MyAngleSort( TAU p1, TAU p2)
{
  return fabs(p1.phi) > fabs(p2.phi);
}

EvaluateTauFinder::EvaluateTauFinder() : Processor("EvaluateTauFinder") 
{
  // modify processor description
  _description = "EvaluateTauFinder checks performance of TauFinder and writes output to root file." ;
  
  // register steering parameters: name, description, class-variable, default value
  
 registerProcessorParameter( "FileName_Signal",
			      "Name of the Signal output file "  ,
			      _OutputFile_Signal ,
			      std::string("Signal.root") ) ;

  registerProcessorParameter( "inputCol" ,
                              "Name of the input Collection"  ,
                              _incol ,
                              std::string("TauRec_PFA")) ;
  registerProcessorParameter( "relCol" ,
                              "Name of the LCRelation otput Collection"  ,
                              _colNameTauRecLink ,
                              std::string("TauRecLink_PFO")) ;

   registerInputCollection( LCIO::MCPARTICLE,
			   "MCCollectionName" , 
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

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "TauRecCollection",
                            "Collection of Tau Candidates",
                            _incol ,
                            std::string("TauRec_PFA"));
 
  registerInputCollection( LCIO::LCRELATION,
			   "MCTauLinkCollectionName" , 
			   "Name of the MC Truth Reconstructed Tau collection"  ,
			   _colNameTauRecLink ,
			   std::string("TauRecLink_PFO") ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "MCTauLinkCollectionName" , 
			   "Name of the MC Truth Reconstructed Tau collection"  ,
			   _colNameMCRecLink ,
			   std::string("MCRecLink") ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "TracksTauLinkCollectionName" , 
			   "Name of the Track Truth Reconstructed Tau collection"  ,
			   _colNameTracksRecLink ,
			   std::string("TracksRecLink") ) ;

  
}


void EvaluateTauFinder::init() 
{ 
  std::cout << "INIT CALLED" << std::endl;
  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  _ntot_rec= 0;
  _ntot_mc= 0;
  _ntau_correct= 0;
  _dEsum= 0;
  _dEsumsq= 0;
  _ndE= 0;
  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
 

  rootfile = new TFile((_OutputFile_Signal).c_str (),"RECREATE");

  evtuple=new TNtuple("evtuple","evtuple","EvID:Ntaus_mc:Ntaus_rec:missed:WpD1:WpD2:WmD1:WmD2:LD");
  tautuple=new TNtuple("tautuple","tautuple","EvID:mcE:mcPhi:mcTheta:mcD0:recE:recPhi:recTheta:recD0:recNQ:recN:LD");
  mcmisstuple=new TNtuple("mcmiss","mcmiss","EvID:E:pt:theta:D0:D1:D2");
  taumatchtuple=new TNtuple("taumatch","taumatch","EvID:E:mcE:mcp:mcpt:mcPhi:mcTheta:mcD0:recE:recp:recpt:recPhi:recTheta:recD0:LD");
  tauexacttuple=new TNtuple("tauexact","tauexact","EvID:E:mcE:mcp:mcpt:mcD0:recE:recp:recpt:recD0:ED0seed:LD");
  faketuple =new TNtuple("fake","fake","EvID:parentpdg:D1:D2:recE:recp:recD0:pt:theta:N");
  topofaketuple =new TNtuple("topofake","topofake","EvID:nfake:WpD1:WpD2:WmD1:WmD2");
  leptons=new TNtuple("Leptons","Leptons","EvID:PDG");

  std::cout << "INIT IS DONE" << std::endl;
}

void EvaluateTauFinder::processRunHeader( LCRunHeader* run) 
{ 
  _nRun++ ;
} 

void EvaluateTauFinder::processEvent( LCEvent * evt ) 
{ 
  // this gets called for every event 
  // usually the working horse ...
  
  LCCollection *colMC, *colRECO, *colMCTruth, *colTau;
  LCCollection *colTauRecLink, *colMCRecLink, *colTracksRecLink;
  try {
    colMC = evt->getCollection( _colNameMC ) ;
  } catch (Exception e) {
    colMC = 0;
  }
  
  try {
    colRECO = evt->getCollection( _colNameRECO ) ;
  } catch (Exception e) {
    colRECO = 0;
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
    colTauRecLink = evt->getCollection( _colNameTauRecLink ) ;
  } catch (Exception e) {
    colTauRecLink = 0;
  }
  try {
    colMCRecLink = evt->getCollection( _colNameMCRecLink ) ;
  } catch (Exception e) {
    colMCRecLink = 0;
  }
  
  try {
    colTracksRecLink = evt->getCollection( _colNameTracksRecLink ) ;
  } catch (Exception e) {
    colTracksRecLink = 0;
  }
  
  
  _nEvt = evt->getEventNumber();  
  
  int ntau_mc=0;
  int ntau_rec=0;
  int missed=0;

  LCRelationNavigator* relationNavigatorTau = 0;
  LCRelationNavigator* relationNavigatorMC = 0;
  LCRelationNavigator* relationNavigatorTracks = 0;
  LCRelationNavigator* relationNavigatorPFOMC = 0;
  if( colTauRecLink != 0)
    relationNavigatorTau = new LCRelationNavigator( colTauRecLink );
  if( colMCRecLink != 0)
    relationNavigatorMC = new LCRelationNavigator( colMCRecLink );
  if( colTracksRecLink != 0)
    relationNavigatorTracks = new LCRelationNavigator( colTracksRecLink );
  if( colMCTruth != 0)
    relationNavigatorPFOMC = new LCRelationNavigator( colMCTruth );
    
  bool isfake=false;
  int nfakes=0;
  int LDev=0;
  std::vector<TAU> rectauvec;
  if( colTau != 0)
    {
      int nT = colTau->getNumberOfElements();
      ntau_rec=nT;
      int LD=0;
      if(_nEvt<coutUpToEv || _nEvt==coutEv)
	cout<<"EVENT "<<_nEvt<<" with "<<nT<<" taus"<<endl;
      HelixClass *helix = new HelixClass();
      HelixClass *mc_helix = new HelixClass();
      for(int k=0; k < nT; k++) 
	{
	  ReconstructedParticle *tau = dynamic_cast <ReconstructedParticle*>( colTau->getElementAt( k ) );
	  const double *pvec=tau->getMomentum();
	  double pt=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
	  double p=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]+pvec[2]*pvec[2]);
	  double phi=180./TMath::Pi()*atan(pvec[1]/pvec[0]);
	  double theta=180./TMath::Pi()*atan(pt/fabs(pvec[2])); 
	  std::vector< ReconstructedParticle * > tauvec=tau->getParticles();
	  double Eseed=0,D0=0;
	  
	  if(tauvec.size()==1)
	    leptons->Fill(_nEvt,tauvec[0]->getType());
	  int NQ=0;
	  for(unsigned int o=0;o<tauvec.size();o++)
	    {
	      //find seed track for D0
	      if(tauvec[o]->getCharge()!=0)
		{
		  NQ++;
		  if(tauvec[o]->getEnergy()>Eseed && tauvec[o]->getTracks().size()!=0 )
		    {
		      D0=(float)tauvec[o]->getTracks()[0]->getD0();
		      Eseed=tauvec[o]->getEnergy();
		    }
		}
	      //check for leptonic decay
	      if(tauvec[o]->getType()==11 ||tauvec[o]->getType()==13)
		LD=1;
	    }
	  // float mom[3];
// 	  float ver[3];
	  
// 	  for (int icomp=0; icomp<3; ++icomp) {
// 	    mom[icomp]=(float)tau->getMomentum()[icomp];
// 	    VertexImpl *vtx=dynamic_cast<VertexImpl*>(tau->getStartVertex());
// 	    if(vtx)
// 	      {
// 		const float *vpos=vtx->getPosition();
// 		ver[icomp]=vpos[icomp];
// 	      }
// 	    else
// 	      ver[icomp]=0;
// 	  }
	  
// 	  float charge = tau->getCharge(); 
// 	  helix->Initialize_VP(ver,mom,charge,_bField);
	  //double D0=fabs(helix->getD0());
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    cout<<tau->getEnergy()<<" "<<phi<<" "<<theta<<" "<<D0<<" "<<tauvec.size()<<endl;
	  
	  TAU rtau;
	  rtau.E=tau->getEnergy();
	  rtau.phi=phi;
	  rtau.theta=theta;
	  rtau.D0=D0;
	  rtau.N=tauvec.size();	 
	  rtau.NQ=NQ;	 
	  rectauvec.push_back(rtau);
	  _ntot_rec++;
	  
	  //follow the chain back to mc truth
	  if(relationNavigatorTau)
	    {
	      bool istau=false;
	      bool contaminated=false;
	      MCParticle *mctau=NULL;
	      EVENT::LCObjectVec relobjFROM = relationNavigatorTau->getRelatedToObjects(tau);
	      for(unsigned int o=0;o<relobjFROM.size();o++)
		{
		  ReconstructedParticle *rec=dynamic_cast <ReconstructedParticle*>(relobjFROM[o]);
		  if(relationNavigatorMC)
		    {
		      EVENT::LCObjectVec relobj = relationNavigatorMC->getRelatedToObjects(rec);
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
		  if(relationNavigatorPFOMC)
		    {
		      EVENT::LCObjectVec relobjMC = relationNavigatorPFOMC->getRelatedToObjects(rec);
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
		}
			
	      //compare tau with mc truth
	      if(mctau)
		{
		  float mc_mom[3];
		  float mc_ver[3];
		  const double *mc_pvec=mctau->getMomentum();
		  double mc_pt=sqrt(mc_pvec[0]*mc_pvec[0]+mc_pvec[1]*mc_pvec[1]);
		  double mc_phi=180./TMath::Pi()*atan(mc_pvec[1]/mc_pvec[0]);
		  double  mc_theta=180./TMath::Pi()*atan(mc_pt/fabs(mc_pvec[2]));
		  
		  for (int icomp=0; icomp<3; ++icomp) {
		    mc_mom[icomp]=(float)mctau->getMomentum()[icomp];
		    mc_ver[icomp]=(float)mctau->getDaughters()[0]->getVertex()[icomp];
		  }
		  float mc_charge = mctau->getCharge(); 
		  mc_helix->Initialize_VP(mc_ver,mc_mom,mc_charge,_bField);
		  double mc_D0=fabs(mc_helix->getD0());
		  double Evis=0,ptvis=0,pvis=0;
		  
		  LoopDaughters(mctau,Evis,ptvis,pvis);
		  
		  taumatchtuple->Fill(_nEvt,mctau->getEnergy(),Evis,pvis,ptvis,mc_phi,mc_theta,mc_D0,tau->getEnergy(),p,pt,phi,theta,D0,LD);
		  if(!contaminated)
		    tauexacttuple->Fill(_nEvt,mctau->getEnergy(),Evis,pvis,ptvis,mc_D0,tau->getEnergy(),p,pt,D0,Eseed,LD);
		  _dEsum+=Evis-tau->getEnergy();
		  _dEsumsq+=(Evis-tau->getEnergy())*(Evis-tau->getEnergy());
		  _ndE++;
		}
	      if(istau)
		_ntau_correct++;
	      else
		{ 
		  int d1=0,d2=0,pdg=0;
		  for(unsigned int o=0;o<relobjFROM.size();o++)
		    {
		      ReconstructedParticle *rec=dynamic_cast <ReconstructedParticle*>(relobjFROM[o]);
		      if(relationNavigatorMC)
			{
			  EVENT::LCObjectVec relobj = relationNavigatorMC->getRelatedToObjects(rec);
			  for(unsigned int m=0;m<relobj.size();m++)
			    {
			      MCParticle *mc=dynamic_cast <MCParticle*>(relobj[m]);
			      if(mc->getCharge()==0)
				continue;
			      MCParticle *dummy=mc;
			      MCParticle *parent=mc;
			      int size=mc->getParents().size();
			      while(size!=0)
				{
				  dummy=parent->getParents()[0];
				  size=dummy->getParents().size();
				  parent=dummy;
				  if(parent->getGeneratorStatus()==2)
				    break;
				}
			      pdg=parent->getPDG();
			      if(parent->getDaughters().size())
				d1=parent->getDaughters()[0]->getPDG();
			      if(parent->getDaughters().size()>1)
				d2=parent->getDaughters()[1]->getPDG();
			    }
			}
		    }
		  const double *tpvec=tau->getMomentum();
		  double tpt=sqrt(tpvec[0]*tpvec[0]+tpvec[1]*tpvec[1]);
		  double ttheta=180./TMath::Pi()*atan(tpt/fabs(tpvec[2]));
		  faketuple->Fill(_nEvt,pdg,d1,d2,tau->getEnergy(),p,D0,tpt,ttheta,tau->getParticles().size());
		  isfake=true;
		  nfakes++;
		}
	    }//relNavTau
	}

      delete helix;
      delete mc_helix;
    }
  
  int D1=0,D2=0,D3=0,D4=0;
  std::vector<TAU> mctauvec;
  if( colMC != 0 ) 
    {
      int nMCP = colMC->getNumberOfElements();
      if(_nEvt<coutUpToEv || _nEvt==coutEv)
	cout<<"MCTRUTH: "<<endl;
      HelixClass * helix = new HelixClass();
      for(int k=0; k < nMCP; k++) 
	{
	  MCParticle *particle = dynamic_cast<MCParticle*> (colMC->getElementAt(k) );

	  if(particle->getGeneratorStatus()==2 && fabs(particle->getPDG())==15 && particle->getDaughters().size()>1 )
	    {
	      if(fabs(particle->getDaughters()[0]->getPDG())==15)
		continue;
	      ntau_mc++;
	      _ntot_mc++;
	      const double *pvec=particle->getMomentum();
	      double pt=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
	      double phi=180./TMath::Pi()*atan(pvec[1]/pvec[0]);
	      double theta=180./TMath::Pi()*atan(pt/fabs(pvec[2])); 
	      float mom[3];
	      float ver[3];
	      for (int icomp=0; icomp<3; ++icomp) {
		mom[icomp]=(float)particle->getMomentum()[icomp];
		ver[icomp]=(float)particle->getDaughters()[0]->getVertex()[icomp];
	      }
	      float charge = particle->getCharge(); 
	      helix->Initialize_VP(ver,mom,charge,_bField);
	      double D0=fabs(helix->getD0());
	      double Evis=0,ptvis=0,pvis=0;
	      LoopDaughters(particle,Evis,ptvis,pvis);
	      
	      TAU mctau;
	      mctau.E=Evis;
	      mctau.phi=phi;
	      mctau.theta=theta;
	      mctau.D0=D0;
	      mctauvec.push_back(mctau);
	      if(_nEvt<coutUpToEv || _nEvt==coutEv)
		cout<<Evis<<" "<<phi<<" "<<theta<<" "<<D0<<endl;
	      
	      //find out which mc taus do not have a link to the rec
	      if(relationNavigatorMC && relationNavigatorTau )
		{
		  bool hasRel=false;
		  LoopDaughtersRelation(particle,relationNavigatorTau ,relationNavigatorMC ,hasRel);
		  if(!hasRel)
		    {		    
		      missed++;
		      int d1=0,d2=0;
		      if(particle->getDaughters().size()==2)
			{
			  d1=particle->getDaughters()[0]->getPDG();
			  d2=particle->getDaughters()[1]->getPDG();
			}
		      mcmisstuple->Fill(_nEvt,Evis,pt,theta,D0,d1,d2);
		      if(_nEvt<coutUpToEv || _nEvt==coutEv)
			cout<<"Missed: "<<Evis<<" "<<D0<<" "<<d1<<" "<<d2<<endl;
		    }
		}
	    }//tau
	  if(particle->getGeneratorStatus()<3 && fabs(particle->getPDG())==24)
	    {
	      if(particle->getPDG()==24 && particle->getDaughters().size()==2)
		{
		  D1=particle->getDaughters()[0]->getPDG();
		  D2=particle->getDaughters()[1]->getPDG();
		}
	      if(particle->getPDG()==-24 && particle->getDaughters().size()==2)
		{
		  D3=particle->getDaughters()[0]->getPDG();
		  D4=particle->getDaughters()[1]->getPDG();
		}
	    }
	}
      delete helix;
    }
  if(isfake)
    topofaketuple->Fill(_nEvt,nfakes,D1,D2,D3,D4);
  //filling the tuple
  evtuple->Fill(_nEvt,ntau_mc,ntau_rec,missed,D1,D2,D3,D4,LDev);
  //sort the mc t and rec taus for comparison
  std::sort(mctauvec.begin(), mctauvec.end(), MyAngleSort);
  std::sort(rectauvec.begin(), rectauvec.end(), MyAngleSort);
  
  unsigned int common=mctauvec.size();
  if(mctauvec.size()>rectauvec.size())
    common=rectauvec.size();
  for(unsigned int p=0;p<common;p++)
    tautuple->Fill(_nEvt,mctauvec[p].E,mctauvec[p].phi,mctauvec[p].theta,mctauvec[p].D0,rectauvec[p].E,rectauvec[p].phi,rectauvec[p].theta,rectauvec[p].D0,rectauvec[p].NQ,rectauvec[p].N);
  if(mctauvec.size()>rectauvec.size())
    {
      for(unsigned int p=common;p<mctauvec.size();p++)
	tautuple->Fill(_nEvt,mctauvec[p].E,mctauvec[p].phi,mctauvec[p].theta,mctauvec[p].D0,0,0,0,0,0,0);
    }
  if(mctauvec.size()<rectauvec.size())
    {
      for(unsigned int p=common;p<rectauvec.size();p++)
	tautuple->Fill(_nEvt,0,0,0,0,rectauvec[p].E,rectauvec[p].phi,rectauvec[p].theta,rectauvec[p].D0,rectauvec[p].NQ,rectauvec[p].N);
    }
  
  
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"--------------------------------------------------------------------------------------------"<<endl;
  
  
  _nEvt ++ ;
  //cleanup
  delete relationNavigatorTau;
  delete relationNavigatorMC;
  delete relationNavigatorTracks;
  delete relationNavigatorPFOMC;
  
}

void EvaluateTauFinder::LoopDaughters(MCParticle *particle,double &Evis,double &ptvis,double &pvis)
{
  for(unsigned int d=0;d<particle->getDaughters().size();d++)
    { 
      MCParticle *daughter=particle->getDaughters()[d];
      //only particles visible in the detector
      if(daughter->hasLeftDetector()==false || fabs(daughter->getPDG())==13)
	{
	  if (daughter->getGeneratorStatus()==1)
	    {
	       //filter out the neutrinos and other invisibles
	      if(!(  fabs(daughter->getPDG())==12 || fabs(daughter->getPDG())==14 || fabs(daughter->getPDG())==16
		     || fabs(daughter->getPDG())==1000022))
		{
		  //		  if(_nEvt<coutUpToEv || _nEvt==coutEv)
		  //cout<<"D vis "<<d<<" "<<daughter->getPDG()<<" "<<daughter->getEnergy()<<endl;
		  Evis+=daughter->getEnergy();
		  const double *mc_pvec=daughter->getMomentum();
		  ptvis+=sqrt(mc_pvec[0]*mc_pvec[0]+mc_pvec[1]*mc_pvec[1]);
		  pvis+=sqrt(mc_pvec[0]*mc_pvec[0]+mc_pvec[1]*mc_pvec[1]+mc_pvec[2]*mc_pvec[2]);
		}
	    }
	  
	}
      if(daughter->getDaughters().size())
	LoopDaughters(daughter,Evis,ptvis,pvis);
    }
}
 
void EvaluateTauFinder::LoopDaughtersRelation(MCParticle *particle,LCRelationNavigator* relationNavigatorTau ,
					      LCRelationNavigator* relationNavigatorMC ,bool &relToTau)
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
	      for(unsigned int o=0;o<relobjTO.size();o++)
		{
		  //relation to the reconstructed tau
		  ReconstructedParticle *rec=dynamic_cast <ReconstructedParticle*>(relobjTO[o]);
		  EVENT::LCObjectVec relobj = relationNavigatorTau->getRelatedFromObjects(rec);
		  if(relobj.size())
		    relToTau=true;
		}
	    }
	}
      if(relToTau)
	break;
      if(daughter->getDaughters().size())
	LoopDaughtersRelation(daughter,relationNavigatorTau,relationNavigatorMC,relToTau);
    }
}


void EvaluateTauFinder::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EvaluateTauFinder::end(){ 
  

  
  std::cout << "EvaluateTauFinder::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "<<std::endl;
  

//   evtuple->Write();
//   tautuple->Write();
//   mcmisstuple->Write();
//   taumatchtuple->Write();
//   tauexacttuple->Write();
//   faketuple->Write();
//   topofaketuple->Write();
  
  //Close File here
  rootfile->Write();
  delete rootfile;
  // rootfile->Close();
}


