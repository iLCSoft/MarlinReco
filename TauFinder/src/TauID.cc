#include "TauID.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCRelation.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <gear/GEAR.h>
#include <gear/BField.h>
#include <marlin/Global.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#define coutEv -1
#define coutUpToEv 0

using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

TauID aTauID ;


bool MyEnergySort_TauID( ReconstructedParticle *p1, ReconstructedParticle *p2)
{
  return fabs(p1->getEnergy()) > fabs(p2->getEnergy());
}

TauID::TauID() : Processor("TauID") 
{
  // modify processor description
  _description = "TauID writes tau candidates as ReconstructedParticles into collection. It runs on a collection of ReconstructedParticels, if you want  to run on MCParticles you have to convert them before hand (use e.g. PrepareRECParticles processor)" ;
  
  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter( "inputCol" ,
                              "Name of the input Collection"  ,
                              _incol ,
                              std::string("PandoraPFOs")) ;

  registerProcessorParameter( "outputCol" ,
                              "Name of the output Collection"  ,
                              _outcol ,
                              std::string("TauID_PFA")) ;
 
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "RECOCollectionName" , 
			   "Name of the ReconstructedParticle collection"  ,
			   _incol ,
			   std::string("PandoraPFOs") ) ;
  
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "TauRecCollection",
                            "Collection of Tau Candidates",
                            _outcol ,
                            std::string("TauRec_PFA"));
    
  registerProcessorParameter( "pt_cut" ,
                              "Cut on pt to suppress background"  ,
                              _ptcut ,
                              (float)0.) ;
  
  registerProcessorParameter( "searchConeAngle" ,
                              "Opening angle of the search cone for tau jet in theta-phi"  ,
                              _coneAngle ,
                              (float)0.7) ;

  registerProcessorParameter( "isolationConeAngle" ,
                              "Outer isolation cone around search cone of tau jet in theta-phi(relativ to cone angle)"  ,
                              _isoAngle ,
                              (float)0.2) ;
  
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
  

}


void TauID::init() 
{ 
  std::cout << "INIT CALLED" << std::endl;
  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  // usually a good idea to
  printParameters() ;

  _fail_Qtr=0;
  _fail_Q=0;
  _fail_isoE=0;
  _nRun = 0 ;
  _nEvt = 0 ;
  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
 
  std::cout << "INIT IS DONE" << std::endl;
}

void TauID::processRunHeader( LCRunHeader* run) 
{ 
  _nRun++ ;
} 

void TauID::processEvent( LCEvent * evt ) 
{ 
  //input collection
  LCCollection *colRECO;
  
  try {
    colRECO = evt->getCollection( _incol ) ;
  } catch (Exception e) {
    colRECO = 0;
  }
 
  _nEvt = evt->getEventNumber();  
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"------EVENT "<<_nEvt<<"---"<<endl;
  
  //output collection
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
	  //user pt cut, to supress background or low energy particles
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
  std::sort(Qvector.begin(), Qvector.end(), MyEnergySort_TauID);
  std::sort(Nvector.begin(), Nvector.end(), MyEnergySort_TauID);
  
  //vector to hold tau candidates
  std::vector<ReconstructedParticleImpl*>  tauvec;
  bool finding_done=false;
  while(Qvector.size() && !finding_done)
    finding_done= FindTau(Qvector,Nvector,tauvec);
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"------------------------"<<endl;
  
  //combine associated particles to tau
  for(unsigned int p=0;p<tauvec.size();p++)
    {
      ReconstructedParticleImpl* tau=tauvec[p];
      double E=0;
      double charge=0;
      double mom[3]={0,0,0};
      
      for(unsigned int tp=0;tp<tau->getParticles().size();tp++)
	{
	  //add up energy and momentum
	  E+=tau->getParticles()[tp]->getEnergy();
	  charge+=tau->getParticles()[tp]->getCharge();
	  mom[0]+=tau->getParticles()[tp]->getMomentum()[0];
	  mom[1]+=tau->getParticles()[tp]->getMomentum()[1];
	  mom[2]+=tau->getParticles()[tp]->getMomentum()[2];
	}
      
      int pdg=15;
      if(charge<0)
	pdg=-15;
	
      tau->setEnergy(E);
      tau->setCharge(charge);
      tau->setMomentum(mom);
      tau->setType(pdg);
    }
  
  //merge taus that are very close together, because they are likely to be from 1 tau that got split in algorithm
  if(tauvec.size()>1)
    {      
      std::vector<ReconstructedParticleImpl*>::iterator iterC=tauvec.begin();
     
      for ( unsigned int t=0; t<tauvec.size() ; t++ )
	{
	  ReconstructedParticleImpl *tau=dynamic_cast<ReconstructedParticleImpl*>(tauvec[t]);
	  const double *mom=tau->getMomentum();
	  double p=sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
	  double angle=100000;
	  
	  for ( unsigned int t2=t+1; t2<tauvec.size() ; t2++ )
	    {
	      iterC=tauvec.begin()+t2;
	      ReconstructedParticleImpl *taun=dynamic_cast<ReconstructedParticleImpl*>(tauvec[t2]);
	      
	      const double *momn=taun->getMomentum();
	      double pn=sqrt(momn[0]*momn[0]+momn[1]*momn[1]+momn[2]*momn[2]);

	      double phin=atan(momn[1]/momn[0]);
	      double etan=0.5*std::log((pn+momn[2])/(pn-momn[2]));
	      //angle=sqrt((phi-phin)*(phi-phin)+(eta-etan)*(eta-etan));
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
		  for(unsigned int i=0;i<taun->getParticles().size();i++)
		    tau->addParticle(taun->getParticles()[i]);
		  
		  // if(_nEvt<coutUpToEv || _nEvt==coutEv)
// 		    {
// 		      cout<<" Tau Merging: "<<endl;
// 		      cout<<t<<" "<<E<<" "<<phi<<" "<<eta<<endl;
// 		      cout<<t2<<" "<<taun->getEnergy()<<" "<<phin<<" "<<etan<<endl;
// 		      cout<<angle<<" -> "<<En<<" "<<endl;
		  //  }
		  delete *iterC;
		  tauvec.erase(iterC);
		  t2--;
		}
	    }
	}
    }
  //test tau candidates for isolation, charge and too many tracks
  std::vector<ReconstructedParticleImpl*>::iterator iter=tauvec.begin();
  for ( unsigned int t=0; t<tauvec.size() ; t++ )
    {
      ReconstructedParticleImpl *tau=dynamic_cast<ReconstructedParticleImpl*>(tauvec[t]);
      //charge consistency
      double kappa=0.1;
      double Qweighted=0,ptsum=0;
      double E=tau->getEnergy();
      double mom[3]={0,0,0};
      mom[0]=tau->getMomentum()[0];
      mom[1]=tau->getMomentum()[1];
      mom[2]=tau->getMomentum()[2];
      double p=sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);

      double phi=atan(mom[1]/mom[0]);
      double eta=0.5*std::log((p+mom[2])/(p-mom[2]));
      
      double theta=atan(sqrt(mom[0]*mom[0]+mom[1]*mom[1])/mom[2]);

      for(unsigned int i=0;i<tau->getParticles().size();i++)
	{
	  //  double pt=sqrt(tau->getParticles()[i]->getMomentum()[0]*tau->getParticles()[i]->getMomentum()[0]
	  //	 +tau->getParticles()[i]->getMomentum()[1]*tau->getParticles()[i]->getMomentum()[1]);
	  Qweighted+=tau->getParticles()[i]->getCharge()*pow(tau->getParticles()[i]->getEnergy(),kappa);
	  ptsum+=pow(tau->getParticles()[i]->getEnergy(),kappa);
	}
      if(abs(Qweighted/ptsum)<0.12 )
	{
	  _fail_Q++;
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    cout<<"Tau "<<"("<<tau->getEnergy()<<","<<theta<<") failed charge: "<<tau->getCharge()<<" "<<abs(Qweighted/ptsum)<<endl;
	  delete *iter;
	  tauvec.erase(iter);
	  t--;
	  continue;
	}


      //too many particles in tau 
      int nQparticles=0, nNparticles=0;
      for(unsigned int i=0;i<tau->getParticles().size();i++)
	{
	  if(tau->getParticles()[i]->getCharge()==1 || tau->getParticles()[i]->getCharge()==-1)
	    nQparticles++;
	  else
	    nNparticles++; 
	}
      if(tau->getParticles().size()>10 || nQparticles>4 || nQparticles==0 )
	{
	  _fail_Qtr++;
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    cout<<"Tau "<<"("<<tau->getEnergy()<<","<<theta<<"): has too many particles: "<<tau->getParticles().size()<<" "<<nQparticles<<endl;
	  delete *iter;
	  tauvec.erase(iter);
	  t--;
	  continue;
	}
      
      //isolation
      double E_iso=0;
      for ( unsigned int s=0; s<Avector.size() ; s++ )
	{
	  ReconstructedParticle *track=dynamic_cast<ReconstructedParticle*>(Avector[s]);
	  const double *momn=track->getMomentum();
	  double pn=sqrt(momn[0]*momn[0]+momn[1]*momn[1]+momn[2]*momn[2]);
	  double angle=acos((mom[0]*momn[0]+mom[1]*momn[1]+mom[2]*momn[2])/
			 (sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])*
			  sqrt(momn[0]*momn[0]+momn[1]*momn[1]+momn[2]*momn[2])));
	  if(angle>_coneAngle && angle<_isoAngle+_coneAngle && E_iso<track->getEnergy())
	    E_iso=track->getEnergy();
	}
      if(E_iso>_isoE)
	{
	  _fail_isoE++;
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    cout<<"Tau "<<"("<<tau->getEnergy()<<","<<theta<<") failed Isolation Energy: "<<E_iso<<endl;
	  delete *iter;
	  tauvec.erase(iter);
	  t--;
	}
      else
	{
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    cout<<"Tau found:"<<tau->getEnergy()<<" "<<nQparticles<<" "<<nNparticles<<" "<<theta<<endl;
	  reccol->addElement(tau);  
	  iter++;
	}
    }
  
 
  evt->addCollection(reccol,_outcol);
    
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"--------------------------------------------------------------------------------------------"<<endl;
  
  _nEvt ++ ;
  
}

bool TauID::FindTau(std::vector<ReconstructedParticle*> &Qvec,std::vector<ReconstructedParticle*> &Nvec,
			std::vector<ReconstructedParticleImpl* > &tauvec)
{
  ReconstructedParticleImpl *tau=new ReconstructedParticleImpl();
  if(Qvec.size()==0)
    {
      // if(_nEvt<coutUpToEv || _nEvt==coutEv)
      //	cout<<"No charged particle in event!"<<endl;
      return true;
    }
 
  //find a good tauseed, check impact parameter 
  ReconstructedParticle *tauseed;
  double D0seed=0;
  std::vector<ReconstructedParticle*>::iterator iterS=Qvec.begin();
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
  
  if(!tauseed)
    {
     //  if(_nEvt<coutUpToEv || _nEvt==coutEv)
// 	cout<<"no further tau seed! D0="<<D0seed<<endl;
      return true;
    }
  
  tau->addParticle(tauseed);

  double  Etau=tauseed->getEnergy();  
  const double *pvec=tauseed->getMomentum();
  double pt=sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
  double p=sqrt(pt*pt+pvec[2]*pvec[2]);
  double theta=atan(pt/pvec[2]);

 
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"seeding: "<<tauseed->getType()<<"\t"<<tauseed->getEnergy()<<"\t"<<theta<<"\t"<<D0seed<<endl;
 
  Qvec.erase(iterS);
  double pvec_tau[3]={0,0,0};
  pvec_tau[0]=tauseed->getMomentum()[0];
  pvec_tau[1]=tauseed->getMomentum()[1];
  pvec_tau[2]=tauseed->getMomentum()[2];

  //assign charged particles to tau candidate
  std::vector<ReconstructedParticle*>::iterator iterQ=Qvec.begin();
  for (unsigned int s=0; s<Qvec.size() ; s++ )
    {
      ReconstructedParticle *track=Qvec[s];

      const double *pvec2=track->getMomentum();
      double pt2=sqrt(pvec2[0]*pvec2[0]+pvec2[1]*pvec2[1]);
      double p2=sqrt(pt2*pt2+pvec2[2]*pvec2[2]);

      double phi2=atan(pvec2[1]/pvec2[0]);
      double eta2=0.5*std::log((p2+pvec2[2])/(p2-pvec2[2]));
      //double angle=sqrt((phi-phi2)*(phi-phi2)+(eta-eta2)*(eta-eta2));
      double angle=acos((pvec2[0]*pvec_tau[0]+pvec2[1]*pvec_tau[1]+pvec2[2]*pvec_tau[2])/
			(sqrt(pvec2[0]*pvec2[0]+pvec2[1]*pvec2[1]+pvec2[2]*pvec2[2])*
			 sqrt(pvec_tau[0]*pvec_tau[0]+pvec_tau[1]*pvec_tau[1]+pvec_tau[2]*pvec_tau[2])));
      if(angle<_coneAngle)
	{
	  tau->addParticle(Qvec[s]);
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    std::cout<<"Adding Q: "<<track->getType()<<"\t"<<track->getEnergy()<<"\t"<<angle<<std::endl;
	  Etau+=Qvec[s]->getEnergy();
	  //combine to new momentum
	  for(int i=0;i<3;i++){
	    pvec_tau[i]=pvec_tau[i]+track->getMomentum()[i];
	  }
	  //adjust cone axis
	  pt=sqrt(pvec_tau[0]*pvec_tau[0]+pvec_tau[1]*pvec_tau[1]);
	  p=sqrt(pt*pt+pvec_tau[2]*pvec_tau[2]);
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
 
      const double *pvec2=track->getMomentum();
      double pt2=sqrt(pvec2[0]*pvec2[0]+pvec2[1]*pvec2[1]);
      double p2=sqrt(pt2*pt2+pvec2[2]*pvec2[2]);

      double phi2=atan(pvec2[1]/pvec2[0]);
      double eta2=0.5*std::log((p2+pvec2[2])/(p2-pvec2[2]));
      //double angle=sqrt((phi-phi2)*(phi-phi2)+(eta-eta2)*(eta-eta2));
      double angle=acos((pvec2[0]*pvec_tau[0]+pvec2[1]*pvec_tau[1]+pvec2[2]*pvec_tau[2])/
			(sqrt(pvec2[0]*pvec2[0]+pvec2[1]*pvec2[1]+pvec2[2]*pvec2[2])*
			 sqrt(pvec_tau[0]*pvec_tau[0]+pvec_tau[1]*pvec_tau[1]+pvec_tau[2]*pvec_tau[2])));

      double theta2=atan(pt2/pvec2[2]);
      if(angle<_coneAngle)
	{
	  tau->addParticle(Nvec[s]);
	  if(_nEvt<coutUpToEv || _nEvt==coutEv)
	    std::cout<<"Adding N: "<<track->getType()<<"\t"<<track->getEnergy()<<"\t"<<theta2<<"\t"<<angle<<std::endl;
	  Etau+=Nvec[s]->getEnergy();
	  //combine to new momentum
	  for(int i=0;i<3;i++){
	    pvec_tau[i]=pvec_tau[i]+track->getMomentum()[i];
	  }
	  //adjust cone axis
	  pt=sqrt(pvec_tau[0]*pvec_tau[0]+pvec_tau[1]*pvec_tau[1]);
	  p=sqrt(pt*pt+pvec_tau[2]*pvec_tau[2]);
	  s--;
	  Nvec.erase(iterN);
	}
      else 
	iterN++;
    }
  tauvec.push_back(tau);
  return false;
}

void TauID::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TauID::end(){ 
  
  
  std::cout << "TauID::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    << std::endl ;
  std::cout << "Reasons for Failure:   " <<std::endl;
  std::cout << "No or to many tracks:  " << _fail_Qtr<< std::endl ;
  std::cout << "No isolation:          " << _fail_isoE<< std::endl ;
  std::cout << "Charge inconsistency:  " << _fail_Q<< std::endl ;
 

}


