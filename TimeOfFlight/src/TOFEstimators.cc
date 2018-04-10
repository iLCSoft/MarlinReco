#include "TOFEstimators.h"
#include "TOFUtils.h"

#include <iostream>
#include <cmath>
#include <functional>


#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/Operators.h>

#include "DDRec/Vector3D.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "marlin/ProcessorEventSeeder.h"
#include <marlin/AIDAProcessor.h>
#include <marlin/Global.h>


#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/ICloud2D.h>
#include <AIDA/IHistogram2D.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//---- ROOT -----
#include "TH1F.h" 

using namespace lcio ;
using namespace marlin ;
using namespace TOFUtils ;

TOFEstimators aTOFEstimators ;



TOFEstimators::TOFEstimators() : Processor("TOFEstimators") {

  // modify processor description
  _description = "TOFEstimators compute some estimators for the time of flight from calorimeter hits" ;

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "ReconstructedParticleCollection" , 
			   "Name of the ReconstructedParticle collection"  ,
			   _colNamePFO ,
			   std::string("PandoraPFOs")
    );


    registerProcessorParameter("MaxLayerNumber",
			     "Use only calorimeter hits up to MaxLayerNumber in TOF estimators",
			     _maxLayerNum,
			     int(100) );


  registerProcessorParameter("TimeResolution",
			     "Assumed time resolution per hit in ps",
			     _resolution,
			     float(0.) );

}



void TOFEstimators::init() { 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);

  Global::EVENTSEEDER->registerProcessor(this);
}


void TOFEstimators::processRunHeader( LCRunHeader* ) {

  _nRun++ ;
} 



void TOFEstimators::processEvent( LCEvent * evt ) { 


  // use the global Marlin random seed for this processor
  gsl_rng_set( _rng, Global::EVENTSEEDER->getSeed(this) ) ;   

  streamlog_out( DEBUG ) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;
  
  

  // get the PFO collection from the event if it exists
  LCCollection* colPFO = nullptr ;
    
  try{
    colPFO = evt->getCollection( _colNamePFO ) ;
  }
  catch(lcio::Exception){
    streamlog_out( DEBUG6 ) << " collection " << _colNamePFO
			    << " not found in event - nothing to do ... " << std::endl ;
  }

  if( colPFO->getTypeName() != LCIO::RECONSTRUCTEDPARTICLE ) {

    streamlog_out( ERROR ) << " collection " << _colNamePFO
			   << " not of type LCIO::RECONSTRUCTEDPARTICLE " << std::endl ;

    colPFO = nullptr ;
  }
    
  if( colPFO != nullptr ){
    
    int nPFO = colPFO->getNumberOfElements()  ;

    // split into charged and neutral PFOs

    std::vector<ReconstructedParticle*> chargedPFOs ;
    std::vector<ReconstructedParticle*> neutralPFOs ;

    chargedPFOs.reserve( nPFO ) ;
    neutralPFOs.reserve( nPFO ) ;

    for(int i=0; i< nPFO ; ++i){ 

      ReconstructedParticle* p = static_cast<ReconstructedParticle*>(  colPFO->getElementAt( i ) ) ;

      if(  p->getClusters().size() != 1 ){
	
	streamlog_out( DEBUG1 ) << " ignore particle w/ cluster number other than one:  " <<  *p << std::endl ; 
	continue ;
      }

      if( std::fabs( p->getCharge() ) < 0.1  && p->getTracks().size() == 0  ) {

	neutralPFOs.push_back( p ) ;
      }
      else if ( std::fabs( p->getCharge() ) > 0.1  && p->getTracks().size() == 1  ) {

	chargedPFOs.push_back( p ) ;
      }
      else {
	streamlog_out( DEBUG1 ) << " ignore particle w/ track number other than zero or one:  " <<  *p << std::endl ; 
      }
    }

    streamlog_out( DEBUG2 ) << "  --- will compute TOF estimators for " << chargedPFOs.size()
			    << " charged and " << neutralPFOs.size()
			    << " neutral particles " << std::endl ; 


    for( auto* p : chargedPFOs ) {
    

      streamlog_out( DEBUG ) << " -----  compute TOF estimators for charged particle : " << *p << std::endl ;


      const double* mom = p->getMomentum() ;
      double momentum = sqrt( mom[0] * mom[0] +  mom[1] * mom[1] +  mom[2] * mom[2] ) ;
      double energy =  p->getEnergy() ;
      
      Cluster* clu = p->getClusters()[0] ;

      Track* trk =  p->getTracks()[0] ;


      const TrackState* tscalo = trk->getTrackState( TrackState::AtCalorimeter ) ; 	

      float x_ref  = tscalo->getReferencePoint()[0] ;
      float y_ref  = tscalo->getReferencePoint()[1] ;
      float z_ref  = tscalo->getReferencePoint()[2] ;



      // =======  use only Ecal hits  (requires the CalorimeterHitType to be set in the digitizer )
      //          with time information ( > 1 ps) and layer <= max layer

      unsigned maxLayerNum = unsigned(  _maxLayerNum ) ;
      std::function<bool(CalorimeterHit*)> selectHits =  [maxLayerNum](CalorimeterHit* h){

	return ( isEcal( h )            &&
		 h->getTime() > 1.e-3   &&
		 layer( h ) <=  maxLayerNum  ) ;
      } ;

      // -------------------------------------------------------------------------------------------

      const CalorimeterHitVec& cluhv = clu->getCalorimeterHits() ;
      
      // create a vector of extended handle objects for relevant calorimeter hits
      
      CaloHitUPtrVec uniqueVec ;
      uniqueVec.reserve( cluhv.size()  ) ;
      
      //-------------------------------------------------------------------------------
      
      CaloHitLayerMap layerMap ;

      for( auto* clh : cluhv ){
	
	if( selectHits( clh ) ){

	  uniqueVec.push_back(  std::unique_ptr<CaloHitData>( new CaloHitData( clh) )  ) ; 

	  CaloHitData* ch = uniqueVec.back().get() ;

	  ch->layer = layer( ch->lcioHit ) ;
	  ch->timeResolution = _resolution ; 

	  ch->smearedTime  = ( _resolution > 0. ?
			       gsl_ran_gaussian( _rng, _resolution / 1000. ) : // convert ps to ns 
			       ch->lcioHit->getTime() ) ;

	  ch->distanceFromIP = dd4hep::rec::Vector3D( clh->getPosition()[0],
						      clh->getPosition()[1],
						      clh->getPosition()[2] ).r() ; 




	  layerMap[ ch->layer ].push_back( ch ) ;
	}
      }
      
      streamlog_out( DEBUG ) << " --- map with hits per layer : " << std::endl ;
      for( auto m : layerMap ){
	streamlog_out( DEBUG ) << "  ----- layer " << m.first << " : " << std::endl ;
	for( auto ch : m.second )
	  streamlog_out( DEBUG ) << "            " << caloTypeStr( ch->lcioHit ) << std::endl ; 		
      }

      
      //=========================  loop over all relevant hits in this particle's cluster ================

      streamlog_out( DEBUG ) << " -----  calorimeter hits used for estimators : " << std::endl ;

      for( auto it = uniqueVec.begin() ; it != uniqueVec.end() ; ++it ){

	CaloHitData* ch = it->get() ;
	CalorimeterHit* calohit = ch->lcioHit ; 

	streamlog_out( DEBUG ) <<  "     ----- " << caloTypeStr( calohit )
			       <<  " -- time : " <<  calohit->getTime()  <<   std::endl ;



	// compute distance of hit from the IP
	float x  = calohit->getPosition()[0] ;
	float y  = calohit->getPosition()[1] ;
	float z  = calohit->getPosition()[2] ;
	      
	float distFromIP = std::sqrt( x*x + y*y + z*z ) ;

	float distFromRef = std::sqrt( (x-x_ref)*(x-x_ref) 
				       + (y-y_ref)*(y-y_ref) 
				       + (z-z_ref)*(z-z_ref)  ) ;
	      
      }
   //=========================================================================================
    
    }
  }


  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() << std::endl ;



  _nEvt ++ ;
}



void TOFEstimators::check( LCEvent *) {}



void TOFEstimators::end(){ 

  gsl_rng_free( _rng );
  

  streamlog_out(MESSAGE) << "TOFEstimators::end()  " << name() 
			 << " processed " << _nEvt << " events in " << _nRun << " runs "
			 << std::endl ;

}


//*******************************************************************************************************

