#include "PrepareRECParticles.h"
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
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/VertexImpl.h>
#include "IMPL/LCFlagImpl.h" 
#include "UTIL/LCRelationNavigator.h"

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

PrepareRECParticles aPrepareRECParticles ;

template<typename T>
T sgn(T n)
{
if (n < 0) return -1;
if (n > 0) return 1;
return 0;
}


PrepareRECParticles::PrepareRECParticles() : Processor("PrepareRECParticles") 
{
  // modify processor description
  _description = "PrepareRECParticles converts input to ReconstructedParticles and puts them into a new collection making sure all the information which is needed to run the TauFinder is provided. ";

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCCollectionName           " , 
			   "Name of the MCParticle collection"  ,
			   _colNameMC ,
			   std::string("MCParticlesSkimmed") ) ;
                                        
  registerInputCollection( LCIO::TRACK,
			   "TrackCollectionName" , 
			   "Name of the Track collection"  ,
			   _colNameTrack ,
			   std::string("LDCTracks") ) ;

  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "PandoraPFOCollectionName" , 
			   "Name of the PandoraPFO collection"  ,
			   _colNamePFO ,
			   std::string("PandoraPFOs") ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "RECOMCTRUTHCollectionName" , 
			   "Name of the MC Truth PFA collection"  ,
			   _colNamePFOMCTruth ,
			   std::string("RecoMCTruthLink") ) ;

  registerProcessorParameter( "outputColMC" ,
                              "Name of the output Collection of refilled information"  ,
                              _outcolMC ,
                              std::string("MCParticles_tau")) ;
  
  registerProcessorParameter( "outputColTracks" ,
                              "Name of the output Collection of refilled information"  ,
                              _outcolTracks ,
                              std::string("Tracks_tau")) ;
 
  registerProcessorParameter( "outputColPFO" ,
                              "Name of the output Collection of refilled information"  ,
                              _outcolPFO ,
                              std::string("PFO_tau")) ;

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RecCollection",
                            "Collection of Rec Particles for TauFinder",
                            _outcolMC ,
                            std::string("MCParticles_tau"));
  
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RecCollection",
                            "Collection of Rec Particles for TauFinder",
                            _outcolTracks ,
                            std::string("Tracks_tau"));
  
   registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RecCollection",
                            "Collection of Rec Particles for TauFinder",
                            _outcolPFO ,
                            std::string("PFO_tau"));
  
  registerOutputCollection( LCIO::LCRELATION,
			   "MCRECLinkCollectionName" , 
			   "Name of the MC Truth ReconstructedParticle collection"  ,
			   _colNameMCTruth ,
			   std::string("MCRecLink") ) ;

  registerOutputCollection( LCIO::LCRELATION,
			    "TrackRECLinkCollectionName" , 
			    "Name of the Track Truth ReconstructedParticle collection"  ,
			    _colNameTrackTruth ,
			    std::string("TracksRecLink") ) ;

  registerOutputCollection( LCIO::LCRELATION,
			    "PFORECLinkCollectionName" , 
			    "Name of the Track Truth ReconstructedParticle collection"  ,
			    _colNamePFOTruth ,
			    std::string("PFORecLink") ) ;
}


void PrepareRECParticles::init() 
{ 
  std::cout << "INIT CALLED" << std::endl;
  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  // usually a good idea to
  printParameters() ;
  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
  _nRun = 0 ;
  _nEvt = 0 ;
   
  std::cout << "INIT IS DONE" << std::endl;
}

void PrepareRECParticles::processRunHeader( LCRunHeader* run) 
{ 
  _nRun++ ;
} 

void PrepareRECParticles::processEvent( LCEvent * evt ) 
{ 
  // this gets called for every event 
  // usually the working horse ...
  
  LCCollection *colMC, *colTrack, *colPFO,*colMCTruth;
  try {
    colMC = evt->getCollection( _colNameMC ) ;
  } catch (Exception e) {
    colMC = 0;
  }
  
  try {
    colTrack = evt->getCollection( _colNameTrack ) ;
  } catch (Exception e) {
    colTrack = 0;
  }
  
  try {
    colPFO = evt->getCollection( _colNamePFO ) ;
  } catch (Exception e) {
    colPFO = 0;
  }
  
  try {
    colMCTruth = evt->getCollection( _colNamePFOMCTruth ) ;
  } catch (Exception e) {
    colMCTruth = 0;
  }

  _nEvt = evt->getEventNumber();  
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"EVENT "<<_nEvt<<endl;
  
  LCCollectionVec *reccol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec *trackcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec *pfocol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  //need to store the vertices because the RecontructeParticle only holds the pointer
  LCCollectionVec *vtxcol = new LCCollectionVec(LCIO::VERTEX);
  //LCRelation stuff
  LCCollectionVec *mc_relationcol = new LCCollectionVec(LCIO::LCRELATION);
  mc_relationcol->parameters().setValue(std::string("FromType"),LCIO::RECONSTRUCTEDPARTICLE);
  mc_relationcol->parameters().setValue(std::string("ToType"),LCIO::MCPARTICLE);
  LCCollectionVec *track_relationcol = new LCCollectionVec(LCIO::LCRELATION);
  track_relationcol->parameters().setValue(std::string("FromType"),LCIO::RECONSTRUCTEDPARTICLE);
  track_relationcol->parameters().setValue(std::string("ToType"),LCIO::TRACK);
  LCCollectionVec *pfo_relationcol = new LCCollectionVec(LCIO::LCRELATION);
  pfo_relationcol->parameters().setValue(std::string("FromType"),LCIO::RECONSTRUCTEDPARTICLE);
  pfo_relationcol->parameters().setValue(std::string("ToType"),LCIO::RECONSTRUCTEDPARTICLE);
   
  LCRelationNavigator* relationNavigatorPFOMC = 0;
  if( colMCTruth != 0)
    relationNavigatorPFOMC = new LCRelationNavigator( colMCTruth );
  
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
	      if( (mc->getMass()==0 && fabs(mc->getPDG())!=22)
		  || mc->getPDG()==1000022)
		continue;
	      ReconstructedParticleImpl *rec = new ReconstructedParticleImpl();
	      //copy values from MC to REC
	      rec->setMomentum(mc->getMomentum());
	      rec->setType(mc->getPDG());
	      rec->setEnergy(mc->getEnergy());
	      rec->setMass(mc->getMass());
	      rec->setCharge(mc->getCharge());
	      VertexImpl *vtx=new VertexImpl();
	      vtx->setPosition(mc->getVertex()[0],mc->getVertex()[1],mc->getVertex()[2]);
	      vtxcol->addElement(vtx);
	      rec->setStartVertex(vtx);
	      reccol->addElement( rec );
	      LCRelationImpl *rel = new LCRelationImpl(rec,mc);
	      mc_relationcol->addElement( rel );
	    }
	}
    }
   
  //convert TRACKS
  if( colTrack != 0 ) 
    {
      int nt=colTrack->getNumberOfElements();
      for(int n=0;n<nt;n++)
	{
	  Track *tr=dynamic_cast <Track*>(colTrack->getElementAt(n));
	  ReconstructedParticleImpl *trec = new ReconstructedParticleImpl();
	  //momentum of track assuming B along z
	  double pt=fabs(_bField/tr->getOmega())*3e-4;
	  double p=fabs(pt/cos(atan(tr->getTanLambda())));	  
	  double mom[3];
	  mom[0]=pt*cos(tr->getPhi());
	  mom[1]=pt*sin(tr->getPhi());
	  mom[2]=pt*tr->getTanLambda();
	  double charge = sgn(_bField/tr->getOmega());
	  trec->setMomentum(mom);
	  trec->setType(-1);
	  trec->setEnergy(p);
	  trec->setCharge(charge);
	  VertexImpl *vtx=new VertexImpl();
	  vtx->setPosition(tr->getReferencePoint()[0],tr->getReferencePoint()[1],tr->getReferencePoint()[2]);
	  vtxcol->addElement(vtx);
	  trec->setStartVertex(vtx);
	  trackcol->addElement( trec );
	  LCRelationImpl *rel = new LCRelationImpl(trec,tr);
	  track_relationcol->addElement( rel );
	}
    }
  
  //convert PANDORAPFOs
  //here we also have to make sure a StartVertex is set correctly
  if( colPFO != 0 ) 
    {
      int nt=colPFO->getNumberOfElements();
      for(int n=0;n<nt;n++)
	{
	  ReconstructedParticle *pfo=dynamic_cast < ReconstructedParticle*>(colPFO->getElementAt(n));
	  ReconstructedParticleImpl *prec = new ReconstructedParticleImpl(); 
	  //copy the pfo object
	  prec->setType (pfo->getType());
	  prec->setMomentum (pfo->getMomentum());
	  prec->setEnergy(pfo->getEnergy());
	  prec->setCharge (pfo->getCharge());
	  prec->setCovMatrix(pfo->getCovMatrix());
	  prec->setMass(pfo->getMass());
	  prec->setReferencePoint (pfo->getReferencePoint());
	  prec->setGoodnessOfPID(pfo->getGoodnessOfPID());
	  //copying of ParticleId vector causes problems because it is also deleted in ~ReconstructedParticle
	  //left it out because it is not needed and if so can be retrieved from the original PFO via Relations
	  EVENT::ParticleID *pid=dynamic_cast<EVENT::ParticleID *>(pfo->getParticleIDUsed());
	  prec->setParticleIDUsed(pid);
	  const EVENT::ReconstructedParticleVec &recpv=dynamic_cast<const EVENT::ReconstructedParticleVec &>(pfo->getParticles());
	  for(unsigned int s=0;s<recpv.size();s++)
	    prec->addParticle(recpv[s]);
	  const EVENT::ClusterVec &cv=pfo->getClusters();
	  for(unsigned int s=0;s<cv.size();s++)
	    prec->addCluster(cv[s]);
	  //check if a StartVertex is set. If not compute one from tracks
	  const EVENT::TrackVec &tv=dynamic_cast<const EVENT::TrackVec &>(pfo->getTracks());
	  VertexImpl *vtx=new VertexImpl();
	  bool vtxset=false;
	  Vertex *vtx2=dynamic_cast<Vertex*>(pfo->getStartVertex());
	  const float *pos=0;
	  if(vtx2)
	    {
	      pos=vtx2->getPosition();
	      vtx->setPosition(pos[0],pos[1],pos[2]);
	      vtxset=true;
	    }
	  else
	    vtx->setPosition(0.,0.,0.);
	  
	  if(tv.size() && !vtxset)
	    {
	      double p=0;
	      //take the track with the highest momentum to set Vertex
	      for(unsigned int s=0;s<tv.size();s++)
		{
		  prec->addTrack(tv[s]);				
		  //momentum of track assuming B along z
		  double pt=fabs(_bField/tv[s]->getOmega())*3e-4;
		  double mom=fabs(pt/cos(atan(tv[s]->getTanLambda())));	  
		  if(mom>p)
		    vtx->setPosition(tv[s]->getReferencePoint()[0],tv[s]->getReferencePoint()[1],tv[s]->getReferencePoint()[2]);
		  p=mom;
		}
	    }
	  
	  vtxcol->addElement(vtx);
	  prec->setStartVertex(vtx);
	  pfocol->addElement( prec );
	  LCRelationImpl *rel = new LCRelationImpl(prec,pfo);
	  pfo_relationcol->addElement( rel );
	}

    }

  
  evt->addCollection(reccol,_outcolMC);
  evt->addCollection(mc_relationcol,_colNameMCTruth);
  evt->addCollection(trackcol,_outcolTracks);
  evt->addCollection(track_relationcol,_colNameTrackTruth);
  evt->addCollection(pfocol,_outcolPFO);
  evt->addCollection(pfo_relationcol,_colNamePFOTruth);
  evt->addCollection(vtxcol,"VTXCollection");
  
  
  if(_nEvt<coutUpToEv || _nEvt==coutEv)
    cout<<"--------------------------------------------------------------------------------------------"<<endl;
   
  _nEvt ++ ;
  
}



void PrepareRECParticles::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void PrepareRECParticles::end(){ 
  
  
  std::cout << "PrepareRECParticles::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    << std::endl ;
  
  //Close File here
  

}


