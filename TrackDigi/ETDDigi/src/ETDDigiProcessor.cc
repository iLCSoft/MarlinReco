/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "ETDDigiProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>
#include <EVENT/MCParticle.h>
//#include <math.h>
#include <cmath>
#include <algorithm>

#include <gsl/gsl_randist.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

ETDDigiProcessor aETDDigiProcessor ;


ETDDigiProcessor::ETDDigiProcessor() : Processor("ETDDigiProcessor") {
  
  // processor description
  _description = "ETDDigiProcessor creates ETD TrackerHits from SimTrackerHits" ;
  
  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "CollectionName" , 
                           "Name of the SimTrackerHit collection"  ,
                           _colName ,
                           std::string("ETDTrackerHits") ) ;
  
  registerProcessorParameter( "PointResolution" ,
                              "Point Resolution in ETD"  ,
                              _pointReso ,
                              (float)0.010) ;
  
  registerProcessorParameter( "MomentumCutForDRays",
                              "Momentum Cut For D Rays (GeV)",
                              _momCut ,
                              float(10.0));
  
  registerProcessorParameter( "RemoveDrays",
                              "Remove D rays?",
                              _removeDrays ,
                              int(0));

  std::vector<int> activeLayers ;
  activeLayers.push_back( -2 ) ;
  activeLayers.push_back(  2 ) ;
  
  registerProcessorParameter( "ActiveLayers",
                              "only hits from active layers are digitized (mimic stereo layers)",
                              _activeLayers,
                              activeLayers );

  registerOutputCollection( LCIO::TRACKERHIT,
                            "OutputCollectionName" , 
                            "Name of the TrackerHit output collection"  ,
                            _outColName ,
                            std::string("ETDTrackerHits") ) ;

}


void ETDDigiProcessor::init() { 

  // usually a good idea to
  printParameters() ;
  _nRun = 0 ;
  _nEvt = 0 ;

  //intialise random number generator 
  r = gsl_rng_alloc(gsl_rng_ranlxs2);
  
}

void ETDDigiProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void ETDDigiProcessor::processEvent( LCEvent * evt ) { 

  streamlog_out( DEBUG ) << " processing evt : " << evt->getEventNumber() << std::endl ;


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
      
      
      // fg: the layer number is the cellID (for Mokka at least) 
      int cellID = SimTHit->getCellID() ;
      
      if( find( _activeLayers.begin() ,_activeLayers.end() , cellID ) == _activeLayers.end() ) {
        
        continue ;   // ----------------- ignore hit 
      } 
      
      
      //         if (_momCut > 0.) {
      //           float mom = 0;
      //           for (int i=0;i<3;++i) 
      //             mom += SimTHit->getMomentum()[i]*SimTHit->getMomentum()[i];
      //           mom = sqrt(mom);
      //           bool accept = true; 
      //           if (_removeDrays==1 && mom < _momCut) 
      //             accept = false;
      //           if (accept) {
      
      // fg: logic should be - if (_removeDrays) check for momentum cut ...
      if (_removeDrays==1 && _momCut > 0. ){
        
        float mom = SimTHit->getMomentum()[0]*SimTHit->getMomentum()[0]
          +         SimTHit->getMomentum()[1]*SimTHit->getMomentum()[1]
          +         SimTHit->getMomentum()[2]*SimTHit->getMomentum()[2] ;
        
        
        mom = sqrt(mom);
        
        if (mom < _momCut) 
          
          continue ;    // ----------------- ignore hit 
      }
      
      const double *pos ;
      pos =  SimTHit->getPosition() ;  
      
      double xSmear = gsl_ran_gaussian(r,_pointReso);
      double ySmear = gsl_ran_gaussian(r,_pointReso);
      
      double newPos[3] ;
      newPos[0] = pos[0] + xSmear;
      newPos[1] = pos[1] + ySmear;
      // No semaring of Z coordinate
      // position of ETD layer is fixed along Z axis
      newPos[2] = pos[2] ;
      
      float de_dx ;
      float dedxSmear = 0.0 ;
      de_dx = SimTHit->getdEdx() ;
      de_dx = de_dx + dedxSmear ; 
      MCParticle *mcp ;
      mcp = SimTHit->getMCParticle() ;
      
      //store hit variables
      TrackerHitImpl* trkHit = new TrackerHitImpl ;
      
      //FIXME: SJA: this is a temporary work around the set'er should take a const double * 
      trkHit->setPosition(  newPos  ) ;
      
      trkHit->setdEdx( de_dx ) ;
      trkHit->setType( 200+abs(SimTHit->getCellID()));
      float covMat[TRKHITNCOVMATRIX]={_pointReso*_pointReso,0.,_pointReso*_pointReso,0.,0.,0.};
      trkHit->setCovMatrix(covMat);      
      // push back the SimTHit for this TrackerHit
      trkHit->rawHits().push_back( SimTHit ) ;
      trkhitVec->addElement( trkHit ) ; 
      
      
    }
    evt->addCollection( trkhitVec ,  _outColName ) ;
    
  }
  
  _nEvt ++ ;
}



  void ETDDigiProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ETDDigiProcessor::end(){ 

  gsl_rng_free(r);  
//   std::cout << "ETDDigiProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

