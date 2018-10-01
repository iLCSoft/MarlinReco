#include "TrackZVertexGrouping.h"

// #include <iostream>
// #include <cmath>
// #include <functional>
// #include <algorithm>


#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/VertexImpl.h>

// #include <UTIL/BitField64.h>
// #include <UTIL/ILDConf.h>
#include <UTIL/Operators.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include <marlin/AIDAProcessor.h>
#include <marlin/Global.h>

#include <AIDA/IHistogramFactory.h>

//---- ROOT -----
#include "TH1F.h"

using namespace lcio ;
using namespace marlin ;

TrackZVertexGrouping aTrackZVertexGrouping ;



TrackZVertexGrouping::TrackZVertexGrouping() : Processor("TrackZVertexGrouping") {

  // modify processor description
  _description = "TrackZVertexGrouping: group tracks based on their Z0 significance";

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::TRACK,
			   "TrackCollection" , 
			   "Name of the Track input collection"  ,
			   _colNameTracks,
			   std::string("MarlinTrkTracks")
    );


  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "TrackGroupPFOs" , 
			    "Name of the ReconstructedParticle output collection"  ,
			    _colNameTrkGroupPFOs ,
			    std::string("TrackGroupPFOs")
    );
  
  registerOutputCollection( LCIO::VERTEX,
			    "TrackGroupVertices",
			    "Name of the Vertex output collection"  ,
			    _colNameTrkGroupVertices ,
			    std::string("TrackGroupPFOs")
    );
  
  registerProcessorParameter("Z0SignificanceCut",
			     "Cut for merging tracks groups in Z0 significance",
			     _z0SignificanceCut,
			     float(1.7) );
  
}



void TrackZVertexGrouping::init() { 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

}


void TrackZVertexGrouping::processRunHeader( LCRunHeader* ) {

  _nRun++ ;
} 



void TrackZVertexGrouping::processEvent( LCEvent * evt ) { 

  streamlog_out(DEBUG2) << "   process event: " << evt->getEventNumber() 
			<< "   in run:  " << evt->getRunNumber() << std::endl ;

  // get the Track collection from the event if it exists
  LCCollection* colTrk = nullptr ;
   
  try{
    colTrk = evt->getCollection( _colNameTracks ) ;
  }
  catch(lcio::Exception){
    streamlog_out( DEBUG6 ) << " collection " << _colNameTracks
			    << " not found in event - nothing to do ... " << std::endl ;
  }

  if( colTrk->getTypeName() != LCIO::TRACK ) {

    streamlog_out( ERROR ) << " collection " << _colNameTracks
			   << " not of type LCIO::TRACK " << std::endl ;

    colTrk = nullptr ;
  }
    
  if( colTrk == nullptr )
    return ;
    


  int nTrk = colTrk->getNumberOfElements()  ;
    
    
  for(int i=0; i< nTrk ; ++i){ 
    
    Track* trk = dynamic_cast<Track*>(  colTrk->getElementAt( i ) ) ;
    
  }
      
													      
  //=========================================================================================
    
  
  

  _nEvt ++ ;

}



void TrackZVertexGrouping::check( LCEvent *evt) {

  streamlog_out( DEBUG ) << " --- check called ! " << std::endl ; 


  // create some histograms with beta vs momentum for charged particles
  
  if( isFirstEvent() ){

    // this creates a directory for this processor ....
    AIDAProcessor::histogramFactory( this ) ;

    // _h.resize(5) ;
    // int nBins = 100 ;
    // _h[0] = new TH2F( "hbetaFirstHitsChrg", "beta vs momentum - first hit - charged",    nBins, .1 , 10., nBins, 0.93 , 1.03 ) ; 
    // _h[1] = new TH2F( "hbetaCloseHitsChrg", "beta vs momentum - closest hits - charged", nBins, .1 , 10., nBins, 0.93 , 1.03 ) ; 

    // _h[2] = new TH2F( "hbetaFirstHitsNeut", "beta vs momentum - first hit - neutral",    nBins, .1 , 10., nBins, 0.93 , 1.03 ) ; 
    // _h[3] = new TH2F( "hbetaCloseHitsNeut", "beta vs momentum - closest hits - neutral", nBins, .1 , 10., nBins, 0.93 , 1.03 ) ; 

    // _h[4] = new TH2F( "hbetaLastTrkHit",    "beta vs momentum - last tracker hit", nBins, .1 , 10., nBins, 0.93 , 1.03 ) ; 

  }
  
}

void TrackZVertexGrouping::end(){ 


  streamlog_out(MESSAGE) << "TrackZVertexGrouping::end()  " << name() 
			 << " processed " << _nEvt << " events in " << _nRun << " runs "
			 << std::endl ;

}


//*******************************************************************************************************

