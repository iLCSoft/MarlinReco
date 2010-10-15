#include "TauFinder.h"
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
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include "UTIL/LCRelationNavigator.h"

#include <gear/GEAR.h>
#include <gear/BField.h>
#include <marlin/Global.h>

#include "HelixClass.h"
#include "SimpleLine.h"
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#define coutEv -1
#define coutUpToEv 0

using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

TauFinder aTauFinder ;


bool MyEnergySort( ReconstructedParticle *p1, ReconstructedParticle *p2)
{
  return fabs(p1->getEnergy()) > fabs(p2->getEnergy());
}

TauFinder::TauFinder() : Processor("TauFinder") 
{
  // modify processor description
  _description = "TauFinder writes tau candidates as ReconstructedParticles into collection. It runs on a collection of ReconstructedParticels, if you want  to run on MCParticles you have to convert them before hand (use e.g. PrepareRECParticles processor)" ;
  
  // register steering parameters: name, description, class-variable, default value
  
  registerProcessorParameter( "inputCol" ,
                              "Name of the input Collection"  ,
                              _incol ,
                              std::string("PandoraPFOs")) ;

  registerProcessorParameter( "outputCol" ,
                              "Name of the output Collection"  ,
                              _outcol ,
                              std::string("TauRec_PFA")) ;
 
  //the link between the reconstructed tau and the input particles used for it
  registerProcessorParameter( "relCol" ,
                               "Name of the LCRelation output Collection"  ,
                               _colNameTauRecLink ,
			      std::string("TauRecLink_PFA")) ;

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "TauRecCollection",
                            "Collection of Tau Candidates",
                            _outcol ,
                            std::string("TauRec_PFA"));
  
  registerOutputCollection( LCIO::LCRELATION,
			    "TauRecLinkCollectionName" , 
			    "Name of the Tau link to ReconstructedParticle collection"  ,
			    _colNameTauRecLink ,
			    std::string("TauRecLink_PFA") ) ;
  
  registerProcessorParameter( "pt_cut" ,
                              "Cut on pt to suppress background"  ,
                              _ptcut ,
                              (float)0.) ;
  
  registerProcessorParameter( "searchConeAngle" ,
                              "Opening angle of the search cone for tau jet in rad"  ,
                              _coneAngle ,
                              (float)0.05) ;

  registerProcessorParameter( "isolationConeAngle" ,
                              "Outer isolation cone around search cone of tau jet in rad (relativ to cone angle)"  ,
                              _isoAngle ,
                              (float)0.02) ;
  
  registerProcessorParameter( "isolationEnergy" ,
                              "Energy allowed within isolation cone region"  ,
                              _isoE ,
                              (float)5.0) ;
  
  registerProcessorParameter( "D0seedmax" ,
                              "Limit on D0 for the track seeding the tau jet"  ,
                              _D0seedmax ,
                              (float)0.5) ;

  registerProcessorParameter( "D0seedmin" ,
                              "Limit on D0 for the track seeding the tau jet"  ,
                              _D0seedmin ,
                              (float)1e-5) ;

  registerProcessorParameter( "ptseed" ,
                              "Minimum tranverse momentum of tau seed"  ,
                              _ptseed ,
                              (float)5.0) ;
  
  registerProcessorParameter( "invariant_mass" ,
                              "Upper limit on invariant mass of tau candidate"  ,
                              _minv ,
                              (float)2.0) ;
}


void TauFinder::init() 
{ 
  std::cout << "INIT CALLED" << std::endl;
  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
  _fail_minv=0;
  _fail_minv_neg=0;
  _fail_Qtr=0;
  _fail_isoE=0;
  std::cout << "INIT IS DONE" << std::endl;
}

void TauFinder::processRunHeader( LCRunHeader* run) 
{ 
  _nRun++ ;
} 

void TauFinder::processEvent( LCEvent * evt ) 
{ 
  // this gets called for every event 
  // usually the working horse ...
  
  LCCollection *colRECO;
  
  try {
    colRECO = evt->getCollection( _incol ) ;
  } catch (Exception e) {
    colRECO = 0;
  }
 
  //LCRelation: to keep information from which particles the tau was made
  LCCollectionVec *relationcol = new LCCollectionVec(LCIO::LCRELATION);
  relationcol->parameters().setValue(std::string("FromType"),LCIO::RECONSTRUCTEDPARTICLE);
  relationcol->parameters().setValue(std::string("ToType"),LCIO::RECONSTRUCTEDPARTICLE);
 
  _nEvt = evt->getEventNumber();  
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"------EVENT "<<_nEvt<<"---"<<endl;
  
  LCCollectionVec * reccol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
 
  //sort all input particles into charged and neutral
  std::vector<ReconstructedParticle*> Avector;//all particles
  std::vector<ReconstructedParticle*> Qvector;//charged particles
  std::vector<ReconstructedParticle*> Nvector;//neutral particles


  if( colRECO != 0)
    {
      int nRCP = colRECO->getNumberOfElements();
      for (int i = 0; i < nRCP; i++) 
	{
	  ReconstructedParticle *particle = dynamic_cast <ReconstructedParticle*>( colRECO->getElementAt( i ) );
	  double pt=sqrt(particle->getMomentum()[0]*particle->getMomentum()[0]
			 +particle->getMomentum()[1]*particle->getMomentum()[1]);
	
	  if(pt<_ptcut)   
	    continue;
	  Avector.push_back(particle);
	  if(particle->getCharge()==1 || particle->getCharge()==-1)
	    Qvector.push_back(particle);
	  else
	    Nvector.push_back(particle);
	}
    }//colRECO


  //sort mc vec according to energy
  std::sort(Qvector.begin(), Qvector.end(), MyEnergySort);
  std::sort(Nvector.begin(), Nvector.end(), MyEnergySort);
  
  //vector to hold tau candidates
  std::vector<std::vector<ReconstructedParticle*> > tauvec;
  bool finding_done=false;
  while(Qvector.size() && !finding_done)
    finding_done= FindTau(Qvector,Nvector,tauvec);
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"------------------------"<<endl;

 //combine associated particles to tau
  std::vector<std::vector<ReconstructedParticle*> >::iterator iterT=tauvec.begin();
  std::vector<ReconstructedParticleImpl* > tauRecvec;
  //remember number of charged tracks in each tau for possible merging later
  std::vector<int> QTvec;
  std::vector<int> NTvec;
  for(unsigned int p=0;p<tauvec.size();p++)
    {
      ReconstructedParticleImpl *taurec=new ReconstructedParticleImpl();
      double E=0;
      double charge=0;
      double mom[3]={0,0,0};
      int chargedtracks=0;
      int neutraltracks=0;
      std::vector<ReconstructedParticle*> tau=tauvec[p];
     
      for(unsigned int tp=0;tp<tau.size();tp++)
	{
	  //add up energy and momentum
	  E+=tau[tp]->getEnergy();
	  charge+=tau[tp]->getCharge();
	  mom[0]+=tau[tp]->getMomentum()[0];
	  mom[1]+=tau[tp]->getMomentum()[1];
	  mom[2]+=tau[tp]->getMomentum()[2];
	  if(tau[tp]->getCharge())
	    chargedtracks++;
	  else
	    neutraltracks++;
	  LCRelationImpl *rel = new LCRelationImpl(taurec,tau[tp]);
	  relationcol->addElement( rel );
	  taurec->addParticle(tau[tp]);
	}
      
      double pt_tau=sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
      double psquare=pt_tau*pt_tau+mom[2]*mom[2];
      double mass_inv=0;
      
      if(E*E<psquare)
	mass_inv=E-sqrt(psquare);
      else
	mass_inv=sqrt(E*E-psquare);
     
      //check for invariant mass
      if(mass_inv>_minv || mass_inv<-0.001 || chargedtracks>4 || chargedtracks==0)
	{
	  if(mass_inv>_minv)
	    _fail_minv++;
	  if(mass_inv<-0.001)
	    _fail_minv_neg++;
	  if(chargedtracks>4 || chargedtracks==0)
	    _fail_Qtr++;

	  iterT=tauvec.erase(iterT);
	  p--;
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    {
	      double phi=180./TMath::Pi()*atan(mom[1]/mom[0]);
	      double theta=180./TMath::Pi()*atan(pt_tau/fabs(mom[2])); 
	      cout<<"Tau candidate failed: minv="<<mass_inv<<"   pt="<<pt_tau<<" "<<E<<" Q trks:"<<chargedtracks<<" N trks:"<<neutraltracks<<" "<<phi<<" "<<theta<<endl;
	    }
	  continue;
	}
      else
	++iterT;
      
      int pdg=15;
      if(charge<0)
	pdg=-15;
	
      taurec->setEnergy(E);
      taurec->setCharge(charge);
      taurec->setMomentum(mom);
      taurec->setType(pdg);
      if(_nEvt<coutUpToEv || _nEvt==coutEv)
	{
	  double phi=180./TMath::Pi()*atan(mom[1]/mom[0]);
	  double theta=180./TMath::Pi()*atan(pt_tau/fabs(mom[2])); 
	  cout<<"Tau candidate "<<p<<": "<<E<<" Q trks:"<<chargedtracks<<" N trks:"<<neutraltracks<<" "<<phi<<" "<<theta<<endl;
	}
      QTvec.push_back(chargedtracks);
      NTvec.push_back(neutraltracks);
      tauRecvec.push_back(taurec);
    }
  //merge taus that are very close together, because they are likely to be from 1 tau that got split in algorithm
  LCRelationNavigator *relationNavigator = new LCRelationNavigator( relationcol );
  if(tauRecvec.size()>1)
    {
      std::vector<ReconstructedParticleImpl*>::iterator iterC=tauRecvec.begin();
      std::vector<ReconstructedParticleImpl*>::iterator iterF=tauRecvec.begin();
      int erasecount=0;
      for ( unsigned int t=0; t<tauRecvec.size() ; t++ )
	{
	  ReconstructedParticleImpl *tau=dynamic_cast<ReconstructedParticleImpl*>(tauRecvec[t]);
	  const double *mom=tau->getMomentum();
	  double pt_tau=sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
	  double phi=180./TMath::Pi()*atan(mom[1]/mom[0]);
	  double theta=180./TMath::Pi()*atan(pt_tau/fabs(mom[2])); 
	  double angle=10000000;
	 
	  for ( unsigned int t2=t+1; t2<tauRecvec.size() ; t2++ )
	    {
	      iterC=tauRecvec.begin()+t2;
	      ReconstructedParticleImpl *taun=dynamic_cast<ReconstructedParticleImpl*>(tauRecvec[t2]);
	      
	      const double *momn=taun->getMomentum();
	      angle=acos((mom[0]*momn[0]+mom[1]*momn[1]+mom[2]*momn[2])/
			 (sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])*
			  sqrt(momn[0]*momn[0]+momn[1]*momn[1]+momn[2]*momn[2])));
	      
	      if(angle<_coneAngle)
		{
		  double E=tau->getEnergy();
		  double En=E+taun->getEnergy();
		  tau->setEnergy(En);
		  double newp[3]={mom[0]+momn[0],mom[1]+momn[1],mom[2]+momn[2]};
		  tau->setMomentum(newp);		  
		  tau->setCharge(tau->getCharge()+taun->getCharge());

		  if(_nEvt<coutUpToEv || _nEvt==coutEv)
		    {
		      cout<<" Tau Merging: "<<endl;
		      cout<<t<<" "<<E<<" "<<phi<<" "<<theta<<endl;
		      cout<<t2<<" "<<taun->getEnergy()<<" "<<angle<<" -> "<<En<<" "<<QTvec[t]+QTvec[t2]<<endl;
		    }
		  
		  //check for invariant mass and number of tracks
		  double ptn=sqrt(newp[0]*newp[0]+newp[1]*newp[1]);
		  double psquaren=ptn*ptn+newp[2]*newp[2];
		  double mass_inv=0;
		  if(En*En<psquaren)
		    mass_inv=En*En-psquaren;
		  else
		    mass_inv= sqrt(En*En-psquaren);
		  //failed to merge
		  if(mass_inv>_minv || mass_inv<-0.001 ||  QTvec[t+erasecount]+QTvec[t2+erasecount]>4)
		    {
		      if(mass_inv>_minv)
			_fail_minv++;
		      if(mass_inv<-0.001)
			_fail_minv_neg++;
		      if(QTvec[t+erasecount]+QTvec[t2+erasecount]>4)
			_fail_Qtr++;
		      
		      delete *iterC;
		      tauRecvec.erase(iterC);
		      erasecount++;
		      t2--;
		      delete *iterF;
		      tauRecvec.erase(iterF);
		      erasecount++;
		      if(tauRecvec.size()>t)
			{
			  tau=*iterF;
			  mom=tau->getMomentum();
			  pt_tau=sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
			  phi=180./TMath::Pi()*atan(mom[1]/mom[0]);
			  theta=180./TMath::Pi()*atan(pt_tau/fabs(mom[2])); 
			  angle=10000000;
			}
		    }
		  else //merge
		    {
		      //set the relations and add particles from one tau to the other
		      std::vector< ReconstructedParticle * > mergetaus=taun->getParticles();
		      for(unsigned int p=0;p<mergetaus.size();p++)
			tau->addParticle(mergetaus[p]);
		      EVENT::LCObjectVec relobjFROM = relationNavigator->getRelatedToObjects(taun);
		      for(unsigned int o=0;o<relobjFROM.size();o++)
			{
			  ReconstructedParticle *rec=dynamic_cast <ReconstructedParticle*>(relobjFROM[o]);
			  LCRelationImpl *rel = new LCRelationImpl(tau,rec);
			  relationcol->addElement( rel );
			}
		      delete *iterC;
		      tauRecvec.erase(iterC);
		      erasecount++;
		      t2--;
		    }
		}
	    }
	  iterF++;
	}
    }
  delete relationNavigator;
  //test for isolation and too many tracks
  std::vector<ReconstructedParticleImpl*>::iterator iter=tauRecvec.begin();
  int erasecount=0;
  for ( unsigned int t=0; t<tauRecvec.size() ; t++ )
    {
      ReconstructedParticleImpl *tau=dynamic_cast<ReconstructedParticleImpl*>(tauRecvec[t]);
      double E_iso=0;
      int nparticles=0;
      const double *pvec_tau=tau->getMomentum();
      //too many particles in tau 
      if(QTvec[t+erasecount]+NTvec[t+erasecount]>10 || QTvec[t+erasecount]>4)
	{
	  _fail_Qtr++;
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    cout<<"Tau "<<tau->getEnergy()<<": too many particles: "<<QTvec[t+erasecount]<<" "<<NTvec[t+erasecount]<<endl;
	  delete *iter;
	  tauRecvec.erase(iter);
	  erasecount++;
	  t--;
	  continue;
	}
      //isolation
      for ( unsigned int s=0; s<Avector.size() ; s++ )
	{
	  ReconstructedParticle *track=dynamic_cast<ReconstructedParticle*>(Avector[s]);
	  const double *pvec=track->getMomentum();
	  double angle=acos((pvec[0]*pvec_tau[0]+pvec[1]*pvec_tau[1]+pvec[2]*pvec_tau[2])/
			    (sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]+pvec[2]*pvec[2])*
			     sqrt(pvec_tau[0]*pvec_tau[0]+pvec_tau[1]*pvec_tau[1]+pvec_tau[2]*pvec_tau[2])));
	  if(angle>_coneAngle && angle<_isoAngle+_coneAngle)
	    {
	      nparticles++;
	      E_iso+=track->getEnergy();
	    }
	}
      
      if(E_iso>_isoE)
	{
	  _fail_isoE++;
	   if(_nEvt<coutUpToEv || _nEvt==coutEv)
	     cout<<"Tau "<<tau->getEnergy()<<": Isolation Energy: "<<E_iso<<" in "<<nparticles<<" particles"<<endl;
	  delete *iter;
	  tauRecvec.erase(iter);
	  erasecount++;
	  t--;
	}
      else
	{
	  reccol->addElement(tau);  
	  if(QTvec[t+erasecount]>4)
	    cout<<"Tau "<<tau->getEnergy()<<" "<<QTvec[t+erasecount]<<" "<<NTvec[t+erasecount]<<endl;
	  iter++;
	}
    }
  
 
  evt->addCollection(reccol,_outcol);
  evt->addCollection(relationcol,_colNameTauRecLink);
  
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"--------------------------------------------------------------------------------------------"<<endl;
  
  _nEvt ++ ;
  
}

bool TauFinder::FindTau(std::vector<ReconstructedParticle*> &Qvec,std::vector<ReconstructedParticle*> &Nvec,
			std::vector<std::vector<ReconstructedParticle*> > &tauvec)
{
  std::vector<ReconstructedParticle*>  tau;
  if(Qvec.size()==0)
    {
      if(_nEvt<coutUpToEv || _nEvt==coutEv)
	cout<<"No charged particle in event!"<<endl;
      return true;
    }
  double OpAngleMax=0,D0seed=0;
  //find a good tauseed, check impact parameter 
  ReconstructedParticle *tauseed;
  std::vector<ReconstructedParticle*>::iterator iterS=Qvec.begin();
  HelixClass *helix = new HelixClass();
  for ( unsigned int s=0; s<Qvec.size() ; s++ )
    {
      tauseed=dynamic_cast<ReconstructedParticle*>(Qvec[s]);
      const EVENT::TrackVec &tv=dynamic_cast<const EVENT::TrackVec &>(tauseed->getTracks());
      float momtr[3];
      double pt=0;
      if(tv.size())
	{
	  double p=0;
	  //take the track with the highest momentum to get impact parameter
	  //(Note: definition of impact parameter depends on the track model 
	  for(unsigned int s=0;s<tv.size();s++)
	    {
	      //momentum of track assuming B along z
	      double pt=fabs(_bField/tv[s]->getOmega())*3e-4;
	      double mom=fabs(pt/cos(atan(tv[s]->getTanLambda())));	  
	      if(mom>p)
		D0seed=tv[s]->getD0();
	      p=mom;
	    }
	  for (int icomp=0; icomp<3; ++icomp) 
	    momtr[icomp]=(float)tauseed->getMomentum()[icomp];

	  pt=sqrt(momtr[0]*momtr[0]+momtr[1]*momtr[1]);
	}
      else
	D0seed=0;
      if(D0seed>_D0seedmin && D0seed<_D0seedmax && pt>_ptseed)
     	break;
      else
	{
	  iterS++;
	  tauseed=NULL;
	}
    }
  delete helix;
  if(!tauseed)
    {
      if(_nEvt<coutUpToEv || _nEvt==coutEv)
	cout<<"no further tau seed! D0="<<D0seed<<endl;
      return true;
    }
  
  double  Etau=tauseed->getEnergy();
 
  tau.push_back(tauseed);
  
  const double *pvec=tauseed->getMomentum();
  double pt=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
  double p=sqrt(pt*pt+pvec[2]*pvec[2]);
  double phi=180./TMath::Pi()*atan(pvec[1]/pvec[0]);
  double theta=180./TMath::Pi()*atan(pt/fabs(pvec[2])); 
 
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"seeding: "<<tauseed->getType()<<"\t"<<tauseed->getEnergy()<<"\t"<<p<<"\t"<<theta<<"\t"<<phi<<"\t"<<D0seed<<endl;
 
  Qvec.erase(iterS);
  double pvec_tau[3]={0,0,0};
  pvec_tau[0]=tauseed->getMomentum()[0];
  pvec_tau[1]=tauseed->getMomentum()[1];
  pvec_tau[2]=tauseed->getMomentum()[2];

  //assign charged particles to tau candidate
  std::vector<ReconstructedParticle*>::iterator iterQ=Qvec.begin();
  for (unsigned int s=0; s<Qvec.size() ; s++ )
    {
      ReconstructedParticle *track=dynamic_cast<ReconstructedParticle*>(Qvec[s]);

      const double *pvec=track->getMomentum();
      double angle=acos((pvec[0]*pvec_tau[0]+pvec[1]*pvec_tau[1]+pvec[2]*pvec_tau[2])/
			(sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]+pvec[2]*pvec[2])*
			 sqrt(pvec_tau[0]*pvec_tau[0]+pvec_tau[1]*pvec_tau[1]+pvec_tau[2]*pvec_tau[2])));
      double pt=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
      double p=sqrt(pt*pt+pvec[2]*pvec[2]);
      double phi=180./TMath::Pi()*atan(pvec[1]/pvec[0]);
      double theta=180./TMath::Pi()*atan(pt/fabs(pvec[2]));
      
      if(angle<_coneAngle)
	{
	  if(angle>OpAngleMax)
	    OpAngleMax=angle;
	  tau.push_back(Qvec[s]);
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    std::cout<<"Adding Q: "<<track->getType()<<"\t"<<track->getEnergy()<<"\t"<<p<<"\t"<<theta<<"\t"<<phi<<std::endl;
	  Etau+=Qvec[s]->getEnergy();
	  //combine to new momentum
	  for(int i=0;i<3;i++){
	    pvec_tau[i]=pvec_tau[i]+track->getMomentum()[i];
	  }
	  Qvec.erase(iterQ);
	  s--;
	}
      else
	iterQ++;
    }
  //assign neutral particles to tau candidate
  std::vector<ReconstructedParticle*>::iterator iterN=Nvec.begin();
  for (unsigned int s=0; s<Nvec.size() ; s++ )
    {
      ReconstructedParticle *track=Nvec[s];
 
      double *pvec=(double*)track->getMomentum();
      double angle=acos((pvec[0]*pvec_tau[0]+pvec[1]*pvec_tau[1]+pvec[2]*pvec_tau[2])/
			(sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]+pvec[2]*pvec[2])*
			 sqrt(pvec_tau[0]*pvec_tau[0]+pvec_tau[1]*pvec_tau[1]+pvec_tau[2]*pvec_tau[2])));
      double pt=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
      double p=sqrt(pt*pt+pvec[2]*pvec[2]);
      double phi=180./TMath::Pi()*atan(pvec[1]/pvec[0]);
      double theta=180./TMath::Pi()*atan(pt/fabs(pvec[2]));
      
      if(angle<_coneAngle)
	{
	  if(angle>OpAngleMax)
	    OpAngleMax=angle;
	  tau.push_back(Nvec[s]);
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    std::cout<<"Adding N: "<<track->getType()<<"\t"<<track->getEnergy()<<"\t"<<p<<"\t"<<theta<<"\t"<<phi<<std::endl;
	  Etau+=Nvec[s]->getEnergy();
	  //combine to new momentum
	  for(int i=0;i<3;i++){
	    pvec_tau[i]=pvec_tau[i]+track->getMomentum()[i];
	  }
	  s--;
	  Nvec.erase(iterN);
	}
      else 
	iterN++;
    }
  tauvec.push_back(tau);
  return false;
}

void TauFinder::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TauFinder::end(){ 
  
  
  std::cout << "TauFinder::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    << std::endl ;
  std::cout << "Reasons for Failure:   " <<std::endl;
  std::cout << "High inverse mass:     " << _fail_minv<< std::endl ;
  std::cout << "Negative inverse mass: " << _fail_minv_neg<< std::endl ;
  std::cout << "No or to many tracks:  " << _fail_Qtr<< std::endl ;
  std::cout << "No isolation        :  " << _fail_isoE<< std::endl ;
 

}


