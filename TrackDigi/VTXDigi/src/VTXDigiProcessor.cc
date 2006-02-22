/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "VTXDigiProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>
#include <EVENT/MCParticle.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

VTXDigiProcessor aVTXDigiProcessor ;


VTXDigiProcessor::VTXDigiProcessor() : Processor("VTXDigiProcessor") {
  
  // modify processor description
  _description = "VTXDigiProcessor should create VTX TrackerHits from SimTrackerHits" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "CollectionName" , 
			      "Name of the SimTrackerHit collection"  ,
			      _colName ,
			      std::string("vxd00_VXD") ) ;
  
  registerProcessorParameter( "SITCollectionName" , 
			      "Name of the SimTrackerHit collection"  ,
			      _colNameSIT ,
			      std::string("sit00_SIT") ) ;

}


void VTXDigiProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void VTXDigiProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void VTXDigiProcessor::processEvent( LCEvent * evt ) { 

  LCCollection* STHcol = 0 ;
  try{
    STHcol = evt->getCollection( _colName ) ;
  }
  catch(DataNotAvailableException &e){
  }

  
    if( STHcol != 0 ){    
    
      LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT )  ;

      int nSimHits = STHcol->getNumberOfElements()  ;

      for(int i=0; i< nSimHits; i++){
      
        SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;

        const int celId = SimTHit->getCellID() ;

         const double *pos ;
        pos =  SimTHit->getPosition() ;  
        
        double xSmear = 0.0 ;
        double ySmear = 0.0 ;
        double zSmear = 0.0 ;

        double newPos[3] ;
        newPos[0] = pos[0] + xSmear;
        newPos[1] = pos[1] + ySmear;
        newPos[2] = pos[2] + zSmear;        

        float de_dx ;
        float dedxSmear = 0.0 ;
        de_dx = SimTHit->getdEdx() ;

        de_dx = de_dx + dedxSmear ; 

        MCParticle *mcp ;
        mcp = SimTHit->getMCParticle() ;

        //store hit variables
        TrackerHitImpl* trkHit = new TrackerHitImpl ;

        //FIXME: SJA: this is a temporary work around the set'er should take a const double * 
        trkHit->setPosition( newPos ) ;

        trkHit->setdEdx( de_dx ) ;
        trkHit->setType( 100+celId ) ;
      
        // 	  push back the SimTHit for this TrackerHit
        trkHit->rawHits().push_back( SimTHit ) ;
        trkhitVec->addElement( trkHit ) ; 
      

      }
      evt->addCollection( trkhitVec , "VTXTrackerHits") ;
    }
    
  try{
    STHcol = evt->getCollection( _colNameSIT ) ;
  }
  catch(DataNotAvailableException &e){
  }
  
  if( STHcol != 0 ){    
    
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT )  ;
    
    int nSimHits = STHcol->getNumberOfElements()  ;
    
    for(int i=0; i< nSimHits; i++){
      
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      
      const int celId = SimTHit->getCellID() ;

      const double *pos ;
      pos =  SimTHit->getPosition() ;  
      
      double xSmear = 0.0 ;
      double ySmear = 0.0 ;
      double zSmear = 0.0 ;
      
      double newPos[3] ;
      newPos[0] = pos[0] + xSmear;
      newPos[1] = pos[1] + ySmear;
      newPos[2] = pos[2] + zSmear;        
      
      float de_dx ;
      float dedxSmear = 0.0 ;
      de_dx = SimTHit->getdEdx() ;
      
      de_dx = de_dx + dedxSmear ; 
      
      MCParticle *mcp ;
      mcp = SimTHit->getMCParticle() ;
      
      //store hit variables
      TrackerHitImpl* trkHit = new TrackerHitImpl ;
      
      //FIXME: SJA: this is a temporary work around the set'er should take a const double * 
      trkHit->setPosition( newPos ) ;
      
      trkHit->setdEdx( de_dx ) ;
      trkHit->setType( 400+celId ) ;
      
      // 	  push back the SimTHit for this TrackerHit
      trkHit->rawHits().push_back( SimTHit ) ;
      trkhitVec->addElement( trkHit ) ; 
      
      
    }
    evt->addCollection( trkhitVec , "SITTrackerHits") ;
  }
  

    
  _nEvt ++ ;

}



  void VTXDigiProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void VTXDigiProcessor::end(){ 
  
//   std::cout << "VTXDigiProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

