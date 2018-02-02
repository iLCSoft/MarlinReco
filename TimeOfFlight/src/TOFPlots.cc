#include "TOFPlots.h"
#include <iostream>
#include <math.h>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

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


TOFPlots aTOFPlots ;


// struct with cluster timing parameters
struct CluTime{

  float clutime ;
  float time05perc ;
  float time10perc ;
  float time20perc ;
  float timeoffastesthit ;
  float time05hits ;
  float time10hits ;
  float time20hits ;
  float cor_ref_time05perc ;
  float cor_ref_time10perc ;
  float cor_ref_time20perc ;
  float cor_ref_time05hits ;
  float cor_ref_time10hits ;
  float cor_ref_time20hits ;

} ;


struct CorRefTime : lcrtrel::LCFloatExtension<CorRefTime> {} ; 
struct SmearedTime : lcrtrel::LCFloatExtension<SmearedTime> {} ; 


CluTime computeClusterTimes(EVENT::Cluster* clu, const  float* refPoint, float timeResolution,  gsl_rng* rng ) ;


// compute the path length along the track from the IP to the entry point into the calorimeter
float computeFlightLength( EVENT::Track* trk ) ;


// helper enum defining histogram index in vector 
namespace TOFHistos {
  enum index{
    hMCPEnergy, hPFOEnergy, hCluEnergy, hHitEnergy, hCorrectedTime, hCorrectedTimeFromRefPoint, 
    htimerefpoint100ps, htimerefpoint50ps, htimerefpoint10ps, haveragetime, haveragecorrectedtime, 
    haveragecorrectedreftime, 
    clutime ,
    time05perc ,
    time10perc ,
    time20perc ,timeoffastesthit,
    time05hits ,
    time10hits ,
    time20hits ,
    cor_ref_time05perc ,
    cor_ref_time10perc ,
    cor_ref_time20perc ,
    cor_ref_time05hits ,
    cor_ref_time10hits ,
    cor_ref_time20hits ,
    flightlength ,
    h_beta_05hits ,
    h_beta_05perc ,
    h_beta_10hits ,
    h_beta_10perc ,
    h_beta_20hits ,
    h_beta_20perc ,
    hBeta_5hitsvsMomentum ,
    hBeta_5percvsMomentum ,
    hBeta_10hitsvsMomentum ,
    hBeta_10percvsMomentum ,
    hBeta_20hitsvsMomentum ,
    hBeta_20percvsMomentum ,
 //-----  keep Size as last :
    Size   
  };

  
  class Histograms{
  public:
    Histograms(std::vector<TH1*>& v) : _h(&v) {}
  
    void create(int idx, const char* n, int nBin=100, double min=0., double max=0. ){
      create( idx , n , n , nBin, min , max ) ; 
    }

    void create(int idx, const char* n, const char* t,  int nBin=100, double min=0., double max=0. ){

      _h->at( idx ) = new TH1D( n, t , nBin , min, max ) ;

      streamlog_out( DEBUG ) << " create histo " <<  n << " at index " << idx << std::endl ;
    }

    void create(int idx, const char* n, const char* t,  int nBin , double* bins ){

      _h->at( idx ) = new TH1D( n, t , nBin , bins ) ;

      streamlog_out( DEBUG ) << " create histo " <<  n << " at index " << idx << std::endl ;
    }

    void fill( int idx , double val, double weight=1.0 ){  _h->at( idx )->Fill( val , weight ) ; }

  protected:

    std::vector<TH1*>* _h;

  };
}

using namespace TOFHistos ;

TOFPlots::TOFPlots() : Processor("TOFPlots") {

  // modify processor description
  _description = "TOFPlots creates plots related to potential time of flight measurements with the calorimter" ;


  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCPCollectionName" , 
			   "Name of the MCParticle collection"  ,
			   _colNameMCP ,
			   std::string("MCParticle")
			   );

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "PFOCollectionName" , 
			   "Name of the ReconstructedParticle collection"  ,
			   _colNamePFO ,
			   std::string("PandoraPFOs")
			   );
}



void TOFPlots::init() { 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);

  Global::EVENTSEEDER->registerProcessor(this);
}


void TOFPlots::processRunHeader( LCRunHeader* ) {

  _nRun++ ;
} 



void TOFPlots::processEvent( LCEvent * evt ) { 


  gsl_rng_set( _rng, Global::EVENTSEEDER->getSeed(this) ) ;   
  streamlog_out( DEBUG4 ) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;
  
  // this gets called for every event 
  // usually the working horse ...
  
  
  // define a histogram pointer

  static AIDA::ICloud2D* hHitDistvsTime ;   
  static AIDA::ICloud2D* hHitDistRefvsTime ;   
  static AIDA::IHistogram2D* hBeta_5hitsvsMomentum ;   
  static AIDA::IHistogram2D* hBeta_5percvsMomentum ;
  static AIDA::IHistogram2D* hBeta_10hitsvsMomentum ;
  static AIDA::IHistogram2D* hBeta_10percvsMomentum ;
  static AIDA::IHistogram2D* hBeta_20hitsvsMomentum ;
  static AIDA::IHistogram2D* hBeta_20percvsMomentum ;
 
  Histograms h(_h) ;
    
  if( isFirstEvent() ){
      
    // this creates a directory for this processor ....
    AIDAProcessor::histogramFactory( this ) ;
      
    _h.resize( TOFHistos::Size ) ;
      
    const int nBins = 100 ;

    h.create(hMCPEnergy  , "hMCPEnergy ", "energy of the MCParticles", nBins ) ; 
    h.create(hPFOEnergy  , "hPFOEnergy ", "energy of the ReconstructedParticles", nBins ) ; 
    h.create(hCluEnergy  , "hCluEnergy ", "energy of the cluster", nBins ) ; 
    h.create(hHitEnergy  , "hHitEnergy ", "energy of the Hit", nBins ) ; 
    h.create(hCorrectedTime  , "hCorrectedTime", "hittime compared to c", nBins ) ; 
    h.create(hCorrectedTimeFromRefPoint  , "hCorrectedTimeFromRefPoint", "hittime from Refpoint compared to c", nBins ) ; 
    h.create(htimerefpoint100ps  , "htimerefpoint100ps", "hittime from Refpoint compared to c with gaus100ps", nBins ) ; 
    h.create(htimerefpoint50ps  , "htimerefpoint50ps", "hittime from Refpoint compared to c with gaus50ps", nBins ) ; 
    h.create(htimerefpoint10ps  , "htimerefpoint10ps", "hittime from Refpoint compared to c with gaus10ps", nBins ) ; 
    h.create(haveragetime  , "haveragetime", "average hittime per cluster", nBins , 7.2 , 8.5 ) ; 

    h.create(clutime ,   "clutime_0ps" , "average cluster time (all hits) smeared with 0ps",    nBins , 7.3 , 8. ) ; 
    h.create(time05perc ,"time05perc_0ps" , "average cluster time ( 5% of hits) smeared with 0ps", nBins , 7.0 , 7.6 ) ; 
    h.create(time10perc ,"time10perc_0ps" , "average cluster time (10% of hits) smeared with 0ps", nBins , 7.1 , 7.7 ) ; 
    h.create(time20perc ,"time20perc_0ps" , "average cluster time (20% of hits) smeared with 0ps", nBins , 7.1 , 7.7 ) ; 
    h.create(timeoffastesthit ,"timeoffastesthit_0ps" , "time of fastest hit smeared with 0ps", nBins , 7. ,7.5 ) ; 
    
    h.create(time05hits ,"time05hits_0ps" , "average cluster time ( 5 fastest hits) smeared with 0ps", nBins , 7.1 , 7.6 ) ; 
    h.create(time10hits ,"time10hits_0ps" , "average cluster time (10 fastest hits) smeared with 0ps", nBins , 7.1 , 7.6 ) ; 
    h.create(time20hits ,"time20hits_0ps" , "average cluster time (20 fastest hits) smeared with 0ps", nBins , 7.2 , 7.7 ) ; 
    
    h.create(cor_ref_time05perc ,"cor_ref_time05perc_0ps" , "average cluster cor_ref_time ( 5% of hits) smeared with 0ps", nBins , 7. ,7.5 ) ; 
    h.create(cor_ref_time10perc ,"cor_ref_time10perc_0ps" , "average cluster cor_ref_time (10% of hits) smeared with 0ps", nBins , 7. ,7.5 ) ; 
    h.create(cor_ref_time20perc ,"cor_ref_time20perc_0ps" , "average cluster cor_ref_time (20% of hits) smeared with 0ps", nBins , 7. ,7.5 ) ; 
    
    h.create(cor_ref_time05hits ,"cor_ref_time05hits_0ps" , "average cluster cor_ref_time ( 5 fastest hits) smeared with 0ps", nBins , 7. ,7.5 ) ; 
    h.create(cor_ref_time10hits ,"cor_ref_time10hits_0ps" , "average cluster cor_ref_time (10 fastest hits) smeared with 0ps", nBins , 7. ,7.5 ) ; 
    h.create(cor_ref_time20hits ,"cor_ref_time20hits_0ps" , "average cluster cor_ref_time (20 fastest hits) smeared with 0ps", nBins , 7. ,7.5 ) ; 

    h.create(haveragecorrectedtime  , "haveragecorrectedtime", "average corrected hittime per cluster", nBins , 0. , 1. ) ; 
    h.create(haveragecorrectedreftime  , "haveragecorrectedreftime", "average corrected hittime from refpoint per cluster", 
	     nBins , 7. , 8. ) ; 

    h.create(flightlength , "flightlength" , "total lenght of flight from the IP to the RefPoint", nBins , 1800 , 2200 ) ; 

    h.create(h_beta_05hits , "beta_05hits_10ps" , "beta factor with cor_ref_time05hits smeared with 10ps", nBins , 0.98 , 1.05 ) ; 
    h.create(h_beta_05perc , "beta_05perc_10ps" , "beta factor with cor_ref_time05perc smeared with 10ps", nBins , 0.98 , 1.05 ) ; 
    h.create(h_beta_10hits , "beta_10hits_10ps" , "beta factor with cor_ref_time10hits smeared with 10ps", nBins , 0.98 , 1.05 ) ; 
    h.create(h_beta_10perc , "beta_10perc_10ps" , "beta factor with cor_ref_time10perc smeared with 10ps", nBins , 0.98 , 1.05 ) ; 
    h.create(h_beta_20hits , "beta_20hits_10ps" , "beta factor with cor_ref_time20hits smeared with 10ps", nBins , 0.98 , 1.05 ) ; 
    h.create(h_beta_20perc , "beta_20perc_10ps" , "beta factor with cor_ref_time20perc smeared with 10ps", nBins , 0.98 , 1.05 ) ; 

    hHitDistvsTime  = AIDAProcessor::histogramFactory(this)->
      createCloud2D( "hHitDistvsTime", "global r of Hit vs time", nBins ) ; 
    
    hHitDistRefvsTime  = AIDAProcessor::histogramFactory(this)->
      createCloud2D( "hHitDistRefvsTime", "dist of Hit from calo entry point vs time", nBins ) ; 

    hBeta_5hitsvsMomentum  = AIDAProcessor::histogramFactory(this)->
      createHistogram2D( "hBeta_5hitsvsMomentum_0ps", "Particles speed in c (beta_05hits) vs momentum smeared with 0ps", nBins, 1. , 10., nBins, 0.93 , 1.03 ) ; 

    hBeta_5percvsMomentum  = AIDAProcessor::histogramFactory(this)->
      createHistogram2D( "hBeta_5percvsMomentum_0ps", "Particles speed in c (beta_05perc) vs momentum smeared with 0ps", nBins, 1. , 10., nBins, 0.93 , 1.03 ) ; 

    hBeta_10hitsvsMomentum  = AIDAProcessor::histogramFactory(this)->
      createHistogram2D( "hBeta_10hitsvsMomentum_0ps", "Particles speed in c (beta_10hits) vs momentum smeared with 0ps", nBins, 1. , 10., nBins, 0.93 , 1.03 ) ; 

    hBeta_10percvsMomentum  = AIDAProcessor::histogramFactory(this)->
      createHistogram2D( "hBeta_10percvsMomentum_0ps", "Particles speed in c (beta_10perc) vs momentum smeared with 0ps", nBins, 1. , 10., nBins, 0.93 , 1.03 ) ; 
    
    hBeta_20hitsvsMomentum  = AIDAProcessor::histogramFactory(this)->
      createHistogram2D( "hBeta_20hitsvsMomentum_0ps", "Particles speed in c (beta_20hits) vs momentum smeared with 0ps", nBins, 1. , 10., nBins, 0.93 , 1.03 ) ; 
    
    hBeta_20percvsMomentum  = AIDAProcessor::histogramFactory(this)->
      createHistogram2D( "hBeta_20percvsMomentum_0ps", "Particles speed in c (beta_20perc) vs momentum smeared with 0ps", nBins, 1. , 10., nBins, 0.93 , 1.03 ) ; 
    
  }

  LCCollection* colMCP=0 ;
  LCCollection* colPFO=0 ;
    
  try{ colMCP = evt->getCollection( _colNameMCP ) ; } catch(lcio::Exception){}
  try{ colPFO = evt->getCollection( _colNamePFO ) ; } catch(lcio::Exception){}


  if( colMCP != NULL ){
    int nMCP = colMCP->getNumberOfElements()  ;

    for(int i=0; i< nMCP ; i++){
      MCParticle* p = dynamic_cast<MCParticle*>( colMCP->getElementAt( i ) ) ;

      // only look at true generator particles
      if( p->getGeneratorStatus() != 1 )
	continue ;
	    
      // fill histogram from LCIO data :
      h.fill( hMCPEnergy,  p->getEnergy() ) ;
    } 
  }



  if( colPFO != NULL ){

    int nPFO = colPFO->getNumberOfElements()  ;

    
    // loop over all PFOs ----------------------------------------
    for(int i=0; i< nPFO ; i++){ 
      ReconstructedParticle* p = dynamic_cast<ReconstructedParticle*>( colPFO->getElementAt( i ) ) ;



      const double* mom = p->getMomentum() ;

      double momentum = sqrt( mom[0] * mom[0] +  mom[1] * mom[1] +  mom[2] * mom[2] ) ;

      // only interested in charged particles
      if( std::fabs( p->getCharge() ) < 0.5 )
	continue ;

      // fill histogram from LCIO data :
      h.fill( hPFOEnergy,  p->getEnergy() ) ;

      // get clusters
      const ClusterVec& cv = p->getClusters() ; 
	    
      // take only the frst 
      Cluster* clu = ( !cv.empty() ? cv[0] : 0 ) ;

      if( clu == 0 ){
	streamlog_out(DEBUG7) << " no cluster in PFO ...." << std::endl ; 
	continue ;
      }


      // get tracks
      const TrackVec& tv = p->getTracks() ; 
	    
      // take only the frst 
      Track* trk = ( !tv.empty() ? tv[0] : 0 ) ;

      if( trk == 0 ){
	streamlog_out(DEBUG7) << " no track in PFO ...." << std::endl ; 
	continue ;
      }

      const TrackState* tscalo = trk->getTrackState( TrackState::AtCalorimeter ) ; 	
      float x_ref  = tscalo->getReferencePoint()[0] ;
      float y_ref  = tscalo->getReferencePoint()[1] ;
      float z_ref  = tscalo->getReferencePoint()[2] ;




      h.fill( hCluEnergy, clu->getEnergy() ) ;
	   
      // getCalorimeterHits
      const CalorimeterHitVec& chv = clu->getCalorimeterHits() ;

      float total_time = 0 ;
      float total_correctedtime = 0 ;
      float total_correctedreftime = 0 ;
      int number_of_hits = 0 ;

      //=========================  loop over all hits in this particle's cluster ================
      for( CalorimeterHit* calohit : chv){ 

	h.fill( hHitEnergy, calohit->getEnergy() ) ;

	float time = calohit->getTime() ;

	if( time < 1.e-3 ){ // ignore hits if time 0 or less than 1 ps
	  continue ;
	}


	long cellID = calohit->getCellID0() ;
	UTIL::BitField64 bf("system:5") ;
	bf.setValue(cellID) ;
	int systemID = bf["system"] ;

	if( systemID != UTIL::ILDDetID::ECAL ){
	  continue ;
	}

	// compute distance of hit from the IP
	float x  = calohit->getPosition()[0] ;
	float y  = calohit->getPosition()[1] ;
	float z  = calohit->getPosition()[2] ;
	      
	float distFromIP = std::sqrt( x*x + y*y + z*z ) ;

	float distFromRef = std::sqrt( (x-x_ref)*(x-x_ref) 
				       + (y-y_ref)*(y-y_ref) 
				       + (z-z_ref)*(z-z_ref)  ) ;
	      
     
	if( time < 10. ){
	  hHitDistvsTime->fill( time, distFromIP ) ;
	  hHitDistRefvsTime->fill( time, distFromRef ) ;
	}

	// compute calohittime from Refpoint compared to a particle travelling with speed of light

	float hittimewithcfromrefpoint = distFromRef / 299.8 ;

	float correctedtimefromrefpoint = time - hittimewithcfromrefpoint ;
	      

	// get hit times assuming different resolutions
	// time is given in ns - assume 10, 50 and 100 ps time resolution

	float smear_100ps = gsl_ran_gaussian( _rng, .1 ) ;
	float timerefpoint100ps = correctedtimefromrefpoint + smear_100ps ;
	float smear_50ps = gsl_ran_gaussian( _rng, .05 ) ;
	float timerefpoint50ps = correctedtimefromrefpoint + smear_50ps ;
	float smear_10ps = gsl_ran_gaussian( _rng, .01 ) ;
	float timerefpoint10ps = correctedtimefromrefpoint + smear_10ps ;

	// compute calohittime compared to a particle travelling with speed of light

	float hittimewithc = std::sqrt( x*x + y*y + z*z )/299.8 ;

	float correctedtime = time - hittimewithc ;

	total_time += time ;
	total_correctedtime += correctedtime ;
	total_correctedreftime += correctedtimefromrefpoint ;
	++number_of_hits ;

	if( correctedtime < 0.3 ){
	  h.fill( hCorrectedTime, correctedtime ) ;
	  h.fill( hCorrectedTimeFromRefPoint, correctedtimefromrefpoint ) ;
	}
	      
	if( timerefpoint100ps < 8.5 ){
	  h.fill( htimerefpoint100ps, timerefpoint100ps ) ;
	}

	if( timerefpoint50ps < 8.5 ){
	  h.fill( htimerefpoint50ps, timerefpoint50ps ) ;
	}

	if( timerefpoint10ps < 8.5 ){
	  h.fill( htimerefpoint10ps, timerefpoint10ps ) ;
	}
      }

      float averagetime = total_time / (float) number_of_hits ;
      float averagecorrectedtime = total_correctedtime / (float) number_of_hits ;
      float averagecorrectedreftime = total_correctedreftime / (float) number_of_hits ;


      double timeRes = 0.0 ; // time resolution in ns
      CluTime ct = computeClusterTimes( clu,  tscalo->getReferencePoint(), timeRes , _rng ) ;

      h.fill( clutime,  ct.clutime ) ;
      h.fill( time05perc,  ct.time05perc ) ;
      h.fill( time10perc,  ct.time10perc ) ;
      h.fill( time20perc,  ct.time20perc ) ;
      h.fill( time05hits,  ct.time05hits ) ;
      h.fill( time10hits,  ct.time10hits ) ;
      h.fill( time20hits,  ct.time20hits ) ;

      h.fill( cor_ref_time05perc,  ct.cor_ref_time05perc ) ;
      h.fill( cor_ref_time10perc,  ct.cor_ref_time10perc ) ;
      h.fill( cor_ref_time20perc,  ct.cor_ref_time20perc ) ;
      h.fill( cor_ref_time05hits,  ct.cor_ref_time05hits ) ;
      h.fill( cor_ref_time10hits,  ct.cor_ref_time10hits ) ;
      h.fill( cor_ref_time20hits,  ct.cor_ref_time20hits ) ;
      h.fill( timeoffastesthit,    ct.timeoffastesthit ) ;


      float length = computeFlightLength( trk) ;

      h.fill( flightlength , length ) ;

      float beta_05hits = ( length / ct.cor_ref_time05hits ) / 299.8 ;
      float beta_05perc = ( length / ct.cor_ref_time05perc ) / 299.8 ;
      float beta_10hits = ( length / ct.cor_ref_time10hits ) / 299.8 ;
      float beta_10perc = ( length / ct.cor_ref_time10perc ) / 299.8 ;
      float beta_20hits = ( length / ct.cor_ref_time20hits ) / 299.8 ;
      float beta_20perc = ( length / ct.cor_ref_time20perc ) / 299.8 ;

      h.fill( h_beta_05hits , beta_05hits ) ;
      h.fill( h_beta_05perc , beta_05perc ) ;
      h.fill( h_beta_10hits , beta_10hits ) ;
      h.fill( h_beta_10perc , beta_10perc ) ;
      h.fill( h_beta_20hits , beta_20hits ) ;
      h.fill( h_beta_20perc , beta_20perc ) ;

      hBeta_5hitsvsMomentum->fill( momentum, beta_05hits ) ;
      hBeta_5percvsMomentum->fill( momentum, beta_05perc ) ;
      hBeta_10hitsvsMomentum->fill( momentum, beta_10hits ) ;
      hBeta_10percvsMomentum->fill( momentum, beta_10perc ) ;
      hBeta_20hitsvsMomentum->fill( momentum, beta_20hits ) ;
      hBeta_20percvsMomentum->fill( momentum, beta_20perc ) ;
      

      streamlog_out(DEBUG) << "   averagetime : " << averagetime 
			   << " and from  computeClusterTimes: " <<  ct.cor_ref_time05perc << std::endl ;


      h.fill( haveragetime, averagetime ) ;
      h.fill( haveragecorrectedtime, averagecorrectedtime ) ;
      h.fill( haveragecorrectedreftime, averagecorrectedreftime ) ;

      //=========================================================================================
    } 
  }



  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() << std::endl ;



  _nEvt ++ ;
}



void TOFPlots::check( LCEvent *) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void TOFPlots::end(){ 

  gsl_rng_free( _rng );
  

  streamlog_out(MESSAGE) << "TOFPlots::end()  " << name() 
			 << " processed " << _nEvt << " events in " << _nRun << " runs "
			 << std::endl ;

}


//*******************************************************************************************************


CluTime computeClusterTimes(EVENT::Cluster* clu, const float* refPoint, float timeResolution,  gsl_rng* rng){

  struct CluTime ct ;

  // getCalorimeterHits
  const CalorimeterHitVec& cluHv = clu->getCalorimeterHits() ;
  std::vector< EVENT::CalorimeterHit* > chv ;
  chv.reserve(  cluHv.size() ) ;

  float total_time = 0 ;
  int number_of_hits = 0 ;



  //=========================  loop over all hits in this particle's cluster ================
  for( CalorimeterHit* calohit : cluHv){ 

    float time = calohit->getTime() ;

    // here we apply a random time shift to simulate the resolution
    float rndTimeShift = gsl_ran_gaussian( rng, timeResolution ) ;
    time += rndTimeShift ;

    calohit->ext<SmearedTime>() = time ;


    if( time < 1.e-3 ){ // ignore hits if time 0 or less than 1 ps
      continue ;
    }
    long cellID = calohit->getCellID0() ;
    UTIL::BitField64 bf("system:5") ;
    bf.setValue(cellID) ;
    int systemID = bf["system"] ;
    
    if( systemID != UTIL::ILDDetID::ECAL ){ // only look at Ecal (barrel hits)
      continue ;
    }

    ++number_of_hits ;
    
    // copy the hit into the hit vector
    chv.push_back( calohit ) ;

    total_time += time ;
  }
  ct.clutime = total_time / number_of_hits ;
  
  // sort the hit vector wrt time using a lambda expression 
  std::sort( chv.begin(), chv.end(), []( EVENT::CalorimeterHit* h1, EVENT::CalorimeterHit* h2) {
      return h1->ext<SmearedTime>() < h2->ext<SmearedTime>() ;   
    });
  

  float total_time05perc = 0 ;
  int number_of_hits05perc = 0 ;

  float total_time10perc = 0 ;
  int number_of_hits10perc = 0 ;

  float total_time20perc = 0 ;
  int number_of_hits20perc = 0 ;

  float total_time05hits = 0 ;
  float total_time10hits = 0 ;
  float total_time20hits = 0 ;

  //  for( CalorimeterHit* calohit : chv){ 

  for(size_t i=0 ; i < chv.size() ; ++i){

    EVENT::CalorimeterHit* calohit = chv.at(i);

    float time = calohit->ext<SmearedTime>() ;

    // compute time corrected for time of flight from reference point
    float x  = calohit->getPosition()[0] ;
    float y  = calohit->getPosition()[1] ;
    float z  = calohit->getPosition()[2] ;

    float x_ref  = refPoint[0] ;
    float y_ref  = refPoint[1] ;
    float z_ref  = refPoint[2] ;

    float distFromRef = std::sqrt( (x-x_ref)*(x-x_ref) 
				   + (y-y_ref)*(y-y_ref) 
				   + (z-z_ref)*(z-z_ref)  ) ;
    
    calohit->ext<CorRefTime>() = time - distFromRef / 299.8  ;

    if( i < (0.05*chv.size () ) ){
      ++number_of_hits05perc ;
      total_time05perc += time ;
    }

     if( i < (0.10*chv.size () ) ){
      ++number_of_hits10perc ;
      total_time10perc += time ;
    }

    if( i < (0.20*chv.size () ) ){
      ++number_of_hits20perc ;
      total_time20perc += time ;
    }

    if( i < 5 ){
      total_time05hits += time ;
    }

    if( i < 10 ){
      total_time10hits += time ;
    }
 
    if( i < 20 ){
      total_time20hits += time ;
    }
   
  }

  ct.time05perc = total_time05perc / number_of_hits05perc ;
  ct.time10perc = total_time10perc / number_of_hits10perc ;
  ct.time20perc = total_time20perc / number_of_hits20perc ;
  ct.time05hits = total_time05hits / 5 ;
  ct.time10hits = total_time10hits / 10 ;
  ct.time20hits = total_time20hits / 20 ;

  ct.timeoffastesthit = ( chv.empty() ? 0. : chv[0]->getTime() ) ;


  // sort the hit vector again - now wrt to corrected time using a lambda expression 
  std::sort( chv.begin(), chv.end(), []( EVENT::CalorimeterHit* h1, EVENT::CalorimeterHit* h2) {
      return h1->ext<CorRefTime>() < h2->ext<CorRefTime>() ;   
    });

  float total_cor_ref_time05perc = 0 ;
  int number_of_cor_ref_hits05perc = 0 ;

  float total_cor_ref_time10perc = 0 ;
  int number_of_cor_ref_hits10perc = 0 ;

  float total_cor_ref_time20perc = 0 ;
  int number_of_cor_ref_hits20perc = 0 ;

  float total_cor_ref_time05hits = 0 ;
  float total_cor_ref_time10hits = 0 ;
  float total_cor_ref_time20hits = 0 ;

  for(size_t i=0 ; i < chv.size() ; ++i){

    EVENT::CalorimeterHit* calohit = chv.at(i);

    double time = calohit->ext<CorRefTime>() ;


    if( i < (0.05*chv.size () ) ){
      ++number_of_cor_ref_hits05perc ;
      total_cor_ref_time05perc += time ;
    }

     if( i < (0.10*chv.size () ) ){
      ++number_of_cor_ref_hits10perc ;
      total_cor_ref_time10perc += time ;
    }

    if( i < (0.20*chv.size () ) ){
      ++number_of_cor_ref_hits20perc ;
      total_cor_ref_time20perc += time ;
    }

    if( i < 5 ){
      total_cor_ref_time05hits += time ;
    }

    if( i < 10 ){
      total_cor_ref_time10hits += time ;
    }
 
    if( i < 20 ){
      total_cor_ref_time20hits += time ;
    }
   
  }

  ct.cor_ref_time05perc = total_cor_ref_time05perc / number_of_cor_ref_hits05perc ;
  ct.cor_ref_time10perc = total_cor_ref_time10perc / number_of_cor_ref_hits10perc ;
  ct.cor_ref_time20perc = total_cor_ref_time20perc / number_of_cor_ref_hits20perc ;
  ct.cor_ref_time05hits = total_cor_ref_time05hits / 5 ;
  ct.cor_ref_time10hits = total_cor_ref_time10hits / 10 ;
  ct.cor_ref_time20hits = total_cor_ref_time20hits / 20 ;




  return std::move( ct ) ;
}



float computeFlightLength( EVENT::Track* trk){

  const TrackState* tsIP = trk->getTrackState( TrackState::AtIP ) ; 	
  const TrackState* tscalo = trk->getTrackState( TrackState::AtCalorimeter ) ; 	

  float Phicalo = tscalo->getPhi() ;
  float PhiIP = tsIP->getPhi() ;

  float Omega = tsIP->getOmega()  ;
  float tanL = tsIP->getTanLambda() ;
  

  float length = (PhiIP-Phicalo)*(1/Omega) * sqrt( 1 + tanL * tanL ) ;


  return length ;
}
