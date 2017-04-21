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

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

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

  registerProcessorParameter( "outputColMC" ,
                              "Name of the output Collection of refilled information"  ,
                              _outcolMC ,
                              std::string("MCParticles_tau")) ;
  
  registerProcessorParameter( "outputColTracks" ,
                              "Name of the output Collection of refilled information"  ,
                              _outcolTracks ,
                              std::string("Tracks_tau")) ;
 
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

  
}


void PrepareRECParticles::init() 
{ 
  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  // usually a good idea to
  printParameters() ;
  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
  _nRun = 0 ;
  _nEvt = 0 ;
   
}

void PrepareRECParticles::processRunHeader( LCRunHeader* )
{ 
  _nRun++ ;
} 

void PrepareRECParticles::processEvent( LCEvent * evt ) 
{ 
  // this gets called for every event 
  // usually the working horse ...
  
  LCCollection *colMC, *colTrack;
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
  
  LCCollectionVec *reccol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec *trackcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  //LCRelation stuff
  LCCollectionVec *mc_relationcol = new LCCollectionVec(LCIO::LCRELATION);
  mc_relationcol->parameters().setValue(std::string("FromType"),LCIO::RECONSTRUCTEDPARTICLE);
  mc_relationcol->parameters().setValue(std::string("ToType"),LCIO::MCPARTICLE);
  LCCollectionVec *track_relationcol = new LCCollectionVec(LCIO::LCRELATION);
  track_relationcol->parameters().setValue(std::string("FromType"),LCIO::RECONSTRUCTEDPARTICLE);
  track_relationcol->parameters().setValue(std::string("ToType"),LCIO::TRACK);
  
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
	      //add the track if charged, so that information for D0 is present
	      if(mc->getCharge())
		{
		  TrackImpl *track=new TrackImpl();
		  float ver[3];
		  float mom[3];
		  float rp[3];
		  for(int i=0;i<3;i++){
		    ver[i]=mc->getVertex()[i];
		    mom[i]=mc->getMomentum()[i];
		  }
		  helix->Initialize_VP(ver,mom,mc->getCharge(),_bField);
		  for(int i=0;i<3;i++)
		    rp[i]=helix->getReferencePoint()[i];
		  track->setReferencePoint(rp);
		  track->setD0(fabs(helix->getD0()));
		  track->setPhi (helix->getPhi0());
		  track->setOmega (helix->getOmega());
		  track->setZ0 (helix->getZ0());
		  track->setTanLambda (helix->getTanLambda());
		  rec->addTrack(track);
		}
	      reccol->addElement( rec );
	      LCRelationImpl *rel = new LCRelationImpl(rec,mc);
	      mc_relationcol->addElement( rel );
	    }
	}
    }
  delete helix;
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
	  trec->addTrack(tr);
	  trackcol->addElement( trec );
	  LCRelationImpl *rel = new LCRelationImpl(trec,tr);
	  track_relationcol->addElement( rel );
	}
    }
  
  
  evt->addCollection(reccol,_outcolMC);
  evt->addCollection(mc_relationcol,_colNameMCTruth);
  evt->addCollection(trackcol,_outcolTracks);
  evt->addCollection(track_relationcol,_colNameTrackTruth);

  _nEvt++;

}



void PrepareRECParticles::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void PrepareRECParticles::end(){ 
  
  
  streamlog_out( DEBUG ) << "PrepareRECParticles::end()  " << name()
                         << " processed " << _nEvt << " events in " << _nRun << " runs "
                         << std::endl ;
  
  //Close File here
  

}


