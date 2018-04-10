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
#include <UTIL/PIDHandler.h>

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
#include "TH2F.h" 

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


  _TOFNames = { "TOFFirstHit",
		"TOFClosestHits", "TOFClosestHitsError",
		"TOFCluster",     "TOFClusterError" } ;

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
    


    PIDHandler pidh( colPFO );
    int algoID = pidh.addAlgorithm( name()  , _TOFNames);


    int nPFO = colPFO->getNumberOfElements()  ;

    // split into charged and neutral PFOs

    std::vector<ReconstructedParticle*> chargedPFOs ;
    std::vector<ReconstructedParticle*> neutralPFOs ;

    chargedPFOs.reserve( nPFO ) ;
    neutralPFOs.reserve( nPFO ) ;


    for(int i=0; i< nPFO ; ++i){ 

      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(  colPFO->getElementAt( i ) ) ;

      bool isCharged = false ;

      if(  pfo->getClusters().size() != 1 ){
	
	streamlog_out( DEBUG1 ) << " ignore particle w/ cluster number other than one:  " <<  *pfo << std::endl ; 
	continue ;
      }

      if( std::fabs( pfo->getCharge() ) < 0.1  && pfo->getTracks().size() == 0  ) {

	neutralPFOs.push_back( pfo ) ;
      }

      else if ( std::fabs( pfo->getCharge() ) > 0.1  && pfo->getTracks().size() == 1  ) {

	chargedPFOs.push_back( pfo ) ;

	isCharged = true ;

      } else {

	streamlog_out( DEBUG1 ) << " ignore particle w/ track number other than zero or one:  " <<  *pfo << std::endl ; 
	continue ;
      }

      streamlog_out( DEBUG1 ) << " --- compute TOF estimators for particle : " << *pfo << std::endl ;


      
      // =======  use only Ecal hits  (requires the CalorimeterHitType to be set in the digitizer )
      //          with time information ( > 1 ps) and layer <= max layer

      unsigned maxLayerNum = unsigned(  _maxLayerNum ) ;
      std::function<bool(CalorimeterHit*)> selectHits =  [maxLayerNum](CalorimeterHit* h){

	return ( isEcal( h )            &&
		 h->getTime() > 1.e-3   &&
		 layer( h ) <=  maxLayerNum  ) ;
      } ;

      // -------------------------------------------------------------------------------------------

      Cluster* clu = pfo->getClusters()[0] ;
      const CalorimeterHitVec& cluhv = clu->getCalorimeterHits() ;
      
      // create vectors of extended handle objects for relevant calorimeter hits
      // one w/ unique_ptr for memory handling
      
      CaloHitUPtrVec uniqueVec ;
      uniqueVec.reserve( cluhv.size()  ) ;
      
      CaloHitDataVec caloHitVec ;
      caloHitVec.reserve( cluhv.size()  ) ;
      //-------------------------------------------------------------------------------
      
      CaloHitLayerMap layerMap ;

      for( auto* clh : cluhv ){
	
	if( selectHits( clh ) ){

	  uniqueVec.push_back(  std::unique_ptr<CaloHitData>( new CaloHitData( clh) )  ) ; 

	  CaloHitData* ch = uniqueVec.back().get() ;

	  caloHitVec.push_back( ch ) ;

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
      
      
      if( layerMap.empty() ) {

	streamlog_out( DEBUG1 ) << " --- not suitable Ecal hits found for particle " << std::endl ;
	continue ;
      }

      // streamlog_out( DEBUG ) << " --- map with hits per layer : " << std::endl ;
      // for( auto m : layerMap ){
      // 	streamlog_out( DEBUG ) << "  ----- layer " << m.first << " : " << std::endl ;
      // 	for( auto ch : m.second )
      // 	  streamlog_out( DEBUG ) << "            " << caloTypeStr( ch->lcioHit ) << std::endl ; 		
      // }
      
      
      // --- define reference point: track state at calo for charged - hit closest to IP for neutral
      //     and direction of straight line - either from IP or from track state at calo

      dd4hep::rec::Vector3D refPoint ;
      dd4hep::rec::Vector3D unitDir ;
      
      if( isCharged ){

	Track* trk =  pfo->getTracks()[0] ;
	const TrackState* tscalo = trk->getTrackState( TrackState::AtCalorimeter ) ; 	
	
	refPoint = { tscalo->getReferencePoint()[0],
		     tscalo->getReferencePoint()[1],
		     tscalo->getReferencePoint()[2] } ;
		
	float tanL = tscalo->getTanLambda() ;
	float theta = std::atan( 1. / tanL ) ;

	unitDir = dd4hep::rec::Vector3D( 1. ,  tscalo->getPhi() , theta , dd4hep::rec::Vector3D::spherical ) ;

	
      } else {  // neutral particle

	CaloHitDataVec chv =  layerMap.begin()->second ; // only look in first layer w/ hits   

	CaloHitData* closestHit =
	  *std::min_element( chv.begin() , chv.end () ,
			    [](CaloHitData* c0, CaloHitData* c1 ){ return c0->distanceFromIP < c1->distanceFromIP  ; }
	    ) ; 
	  
	refPoint = { closestHit->lcioHit->getPosition()[0],
		     closestHit->lcioHit->getPosition()[1],
		     closestHit->lcioHit->getPosition()[2]  } ; 

	dd4hep::rec::Vector3D cluPos( clu->getPosition()[0], 
				      clu->getPosition()[1], 
				      clu->getPosition()[2] ) ;

	unitDir = cluPos.unit() ;
	
      } 
      
      streamlog_out( DEBUG2 ) << " ----- use reference point for TOF : " << refPoint << std::endl ;

      streamlog_out( DEBUG ) << " -----  calorimeter hits considered for estimators : " << std::endl ;


      // ------  loop again over hits and fill missing data

      for( auto ch : caloHitVec ){

	CalorimeterHit* calohit = ch->lcioHit ; 

	dd4hep::rec::Vector3D pos = { ch->lcioHit->getPosition()[0], 
				      ch->lcioHit->getPosition()[1], 
				      ch->lcioHit->getPosition()[2] } ;

	ch->distanceFromReferencePoint = ( pos - refPoint ).r()   ; 


	ch->distancefromStraightline = computeDistanceFromLine( calohit, refPoint, unitDir ) ;

	  
	streamlog_out( DEBUG ) <<  "     ----- " << caloTypeStr( calohit )
			       <<  " --  " <<  ch->toString()   <<   std::endl ;

      }

      // ---- now get hits that are closest to the extrapolated line
      CaloHitDataVec tofHits = findHitsClosestToLine( layerMap ) ;


      streamlog_out( DEBUG )   <<  " ***** hits used for the TOF estimator : " << std::endl ;
      for( auto ch : tofHits ){
	streamlog_out( DEBUG ) <<  "     ----- " << ch->toString() << std::endl ;
      }

      auto t_dt     = computeTOFEstimator( tofHits ) ;
      auto t_dt_clu = computeTOFEstimator( caloHitVec ) ; 
      
      streamlog_out( DEBUG2 ) << "  #### tof ( first ) : " <<  tofHits[0]->smearedTime << " +/- " << 0 << std::endl ; 
      streamlog_out( DEBUG2 ) << "  #### tof ( straight line ) : " << t_dt.first << " +/- " << t_dt.second << std::endl ; 
      streamlog_out( DEBUG2 ) << "  #### tof ( cluster ) : " << t_dt_clu.first << " +/- " << t_dt_clu.second << std::endl ; 

      FloatVec TOF_params = { tofHits[0]->smearedTime,
			      t_dt.first, t_dt.second,
			      t_dt_clu.first, t_dt_clu.second } ; 

      pidh.setParticleID(pfo , 0, 0 , 0.0 , algoID, TOF_params );

													      
   //=========================================================================================
    
    }

    streamlog_out( DEBUG2 ) << "  --- will compute TOF estimators for " << chargedPFOs.size()
			    << " charged and " << neutralPFOs.size()
			    << " neutral particles " << std::endl ; 
    


   for( auto* pfo : chargedPFOs ) {
      
      
      streamlog_out( DEBUG ) << " -----  compute TOF estimators for charged particle : " << *pfo << std::endl ;
      
      
      const double* mom = pfo->getMomentum() ;
      double momentum = sqrt( mom[0] * mom[0] +  mom[1] * mom[1] +  mom[2] * mom[2] ) ;
      double energy =  pfo->getEnergy() ;
      
      Cluster* clu = pfo->getClusters()[0] ;
      
      Track* trk =  pfo->getTracks()[0] ;
      
      
      const TrackState* tscalo = trk->getTrackState( TrackState::AtCalorimeter ) ; 	
      
      float x_ref  = tscalo->getReferencePoint()[0] ;
      float y_ref  = tscalo->getReferencePoint()[1] ;
      float z_ref  = tscalo->getReferencePoint()[2] ;
      
      
    }
    
  }

  streamlog_out(DEBUG2) << "   processed event: " << evt->getEventNumber() 
			<< "   in run:  " << evt->getRunNumber() << std::endl ;



  _nEvt ++ ;
}



void TOFEstimators::check( LCEvent *evt) {

  streamlog_out( DEBUG ) << " ****  check called !!! " << std::endl ; 


  if( isFirstEvent() ){
    _h.resize(5) ;
    int nBins = 100 ;
    _h[0] = new TH2F( "hbetaFirstHit", "beta vs momentum - first hit ", nBins, 1. , 10., nBins, 0.93 , 1.03 ) ; 
    _h[1] = new TH2F( "hbetaCloseHit", "beta vs momentum - closest hits ", nBins, 1. , 10., nBins, 0.93 , 1.03 ) ; 
    _h[2] = new TH2F( "hbetaCluster",  "beta vs momentum - cluster hits ", nBins, 1. , 10., nBins, 0.93 , 1.03 ) ; 
  }

  // get the PFO collection from the event if it exists
  LCCollection* colPFO = nullptr ;
  try{ colPFO = evt->getCollection( _colNamePFO ) ; } catch(lcio::Exception){}

  if( colPFO != nullptr ){
  
    PIDHandler pidh( colPFO );
    int algoID = pidh.getAlgorithmID( name() );
    int fh_idx = pidh.getParameterIndex(algoID,"TOFFirstHit") ;
    int ch_idx = pidh.getParameterIndex(algoID,"TOFClosestHits") ;
    int cl_idx = pidh.getParameterIndex(algoID,"TOFCluster") ;
    
    
    int nPFO = colPFO->getNumberOfElements()  ;

    for( int i=0 ; i< nPFO ; ++i){

      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(  colPFO->getElementAt( i ) ) ;

      const ParticleID& tofPID = pidh.getParticleID( pfo , algoID ) ;
      
      const FloatVec& tofParams = tofPID.getParameters() ;

      streamlog_out( DEBUG ) << " ****  found TOF parameters for pfo w/ size " << tofParams.size() << std::endl ;

      if( !tofParams.empty() && pfo->getTracks().size() == 1 ){

	const double* mom = pfo->getMomentum() ;
	double momentum = sqrt( mom[0] * mom[0] +  mom[1] * mom[1] +  mom[2] * mom[2] ) ;

	Track* trk =  pfo->getTracks()[0] ;
	float length = computeFlightLength( trk ) ;

	double beta_fh = ( length / tofParams[ fh_idx  ] ) / 299.8 ;
	double beta_ch = ( length / tofParams[ ch_idx  ] ) / 299.8 ;
	double beta_cl = ( length / tofParams[ cl_idx  ] ) / 299.8 ;

	_h[0]->Fill( momentum , beta_fh );
	_h[1]->Fill( momentum , beta_ch );
	_h[2]->Fill( momentum , beta_cl );


      } 

    }
  }
    
}




void TOFEstimators::end(){ 

  gsl_rng_free( _rng );
  

  streamlog_out(MESSAGE) << "TOFEstimators::end()  " << name() 
			 << " processed " << _nEvt << " events in " << _nRun << " runs "
			 << std::endl ;

}


//*******************************************************************************************************

