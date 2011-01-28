/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "FTDDigiProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>
#include <EVENT/MCParticle.h>
#include <math.h>

#include <marlin/Global.h>
#include "marlin/ProcessorEventSeeder.h"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/BField.h>
#include <cmath>

#include <gsl/gsl_randist.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

FTDDigiProcessor aFTDDigiProcessor ;



/** helper class to compute the layerID of FTD hits from the z-position */
class FTDLayer{
protected:
  
  //  std::map<int,unsigned> _layerMap  ;
  std::vector< std::pair< double, double > > _layerVec ;

public:
  FTDLayer( const gear::GearParameters& param ){
    
    const DoubleVec& zv = param.getDoubleVals( "FTDZCoordinate" ) ;
    
    
    _layerVec.resize(  zv.size() ) ;
    
    for(unsigned i=0 ; i < zv.size() ; ++i ) { 
      
      //       int zPos = int( std::floor( zv[i] + 0.5  ) ) ;
      //       // store layerId for +/- 1 mm
      //       _layerMap[ zPos - 1 ] = i ;
      //       _layerMap[ zPos     ] = i ;
      //       _layerMap[ zPos + 1 ] = i ;
      
      //       streamlog_out( DEBUG )  << " FTDLayer : " << i <<  "  at : " 
      //                               << zPos << " +/- 1 " << std::endl ;

      
      if( i != 0 ) {  // not first
        
        _layerVec[i].first =  ( zv[i-1] + zv[i] ) / 2. ;
        
      } else {
        
        _layerVec[i].first =  zv[i] - 100. ;
      }
      
      if( i+1 < zv.size() ) { // not last 
        
        _layerVec[i].second = ( zv[i] + zv[i+1] ) / 2. ;
        
      } else {
        
        _layerVec[i].second = zv[i] + 100. ;
        
      }

      streamlog_out( DEBUG )  << " FTDLayer : " << i <<  "  at : " 
                              << _layerVec[i].first << "  -  " <<  _layerVec[i].second 
                              << std::endl ;
    }

  }
  
  unsigned layerID( double zPos ) {
    
    //     std::map<int,unsigned>::iterator it = _layerMap.find( int( std::floor( zPos + 0.5  ) ) ) ;  
    
    //     if( it == _layerMap.end() ){
    
    //       stringstream em ;
    //       em << " FTDDigiProcessor::FTDLayer cannot find layer for z-position : " 
    //          << zPos << std::endl ;
    
    //       throw Exception( em.str() ) ;
    //     }
    
    //     return it->second ;
    //   }
    
    double zAbs = std::abs( zPos ) ;
    
    for(unsigned i=0 ; i < _layerVec.size() ; ++i ) { 
      
      if(   _layerVec[i].first  < zAbs  && zAbs <  _layerVec[i].second ) {
        
        return i ;
      }
    }


    stringstream em ;
    em << " FTDDigiProcessor::FTDLayer cannot find layer for z-position : " 
       << zPos << std::endl ;
    
    throw Exception( em.str() ) ;
  }

} ;


FTDDigiProcessor::FTDDigiProcessor() : Processor("FTDDigiProcessor") {
  
  // processor description
  _description = "FTDDigiProcessor creates FTD TrackerHits from SimTrackerHits" ;
  
  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "CollectionName" , 
                           "Name of the SimTrackerHit collection"  ,
                           _colName ,
                           std::string("ftd01_FTD") ) ;
  
  registerProcessorParameter( "PointResolution" ,
                              "Point Resolution in FTD"  ,
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

  registerOutputCollection( LCIO::TRACKERHIT,
                            "OutputCollectionName" , 
                            "Name of the TrackerHit output collection"  ,
                            _outColName ,
                            std::string("FTDTrackerHits") ) ;

}


void FTDDigiProcessor::init() { 

  // usually a good idea to
  printParameters() ;
  _nRun = 0 ;
  _nEvt = 0 ;

  //intialise random number generator 
  r = gsl_rng_alloc(gsl_rng_ranlxs2);
  Global::EVENTSEEDER->registerProcessor(this);

  const gear::GearParameters& pFTD = Global::GEAR->getGearParameters("FTD");

  _ftdLayers = new FTDLayer( pFTD )  ;

}

void FTDDigiProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void FTDDigiProcessor::processEvent( LCEvent * evt ) { 

  gsl_rng_set( r, Global::EVENTSEEDER->getSeed(this) ) ;   
  streamlog_out( DEBUG ) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;

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

        
        if (_momCut > 0.) {
          float mom = 0;
          for (int i=0;i<3;++i) 
            mom += SimTHit->getMomentum()[i]*SimTHit->getMomentum()[i];

          mom = sqrt(mom);

          bool accept = true; 
          if (_removeDrays==1 && mom < _momCut) 
            accept = false;

          if (accept) {
    
//             const int celId = SimTHit->getCellID() ;
          
            const double *pos ;
            pos =  SimTHit->getPosition() ;  
            
            double xSmear = gsl_ran_gaussian(r,_pointReso);
            double ySmear = gsl_ran_gaussian(r,_pointReso);
          
            double newPos[3] ;
            newPos[0] = pos[0] + xSmear;
            newPos[1] = pos[1] + ySmear;
            // No semaring of Z coordinate
            // position of FTD layer is fixed along Z axis
            newPos[2] = pos[2] ;
          
            float edep ;
            float dedxSmear = 0.0 ;
            edep = SimTHit->getEDep() ;
            edep = edep + dedxSmear ; 
            MCParticle *mcp ;
            mcp = SimTHit->getMCParticle() ;
            
            //store hit variables
            TrackerHitImpl* trkHit = new TrackerHitImpl ;

            //FIXME: SJA: this is a temporary work around the set'er should take a const double * 
            trkHit->setPosition(  newPos  ) ;

            trkHit->setEDep( edep ) ;

            //            trkHit->setType( 200+abs(SimTHit->getCellID()));
            //fg: get layerID from z-position
            trkHit->setType( 201 + _ftdLayers->layerID(  pos[2] ) ) ;
            // FIXME:  need proper way of decoding the layerID in the TrackerHit
            //         for now we use this 'workaround' 201 is expected as FTD ID in 
            //         SiliconTracking


            float covMat[TRKHITNCOVMATRIX]={_pointReso*_pointReso,0.,_pointReso*_pointReso,0.,0.,0.};
            trkHit->setCovMatrix(covMat);      
            // push back the SimTHit for this TrackerHit
            trkHit->rawHits().push_back( SimTHit ) ;
            trkhitVec->addElement( trkHit ) ; 
          }
        }
      }
      evt->addCollection( trkhitVec ,  _outColName ) ;
    }
    
  _nEvt ++ ;
}



  void FTDDigiProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void FTDDigiProcessor::end(){ 

  gsl_rng_free(r);  
//   std::cout << "FTDDigiProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

  delete _ftdLayers ;
}

