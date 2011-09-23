/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "NewFTDDigiProcessor.h"
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

NewFTDDigiProcessor aNewFTDDigiProcessor ;



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


NewFTDDigiProcessor::NewFTDDigiProcessor() : Processor("NewFTDDigiProcessor") {
  
  // processor description
  _description = "NewFTDDigiProcessor creates FTD TrackerHits from SimTrackerHits" ;
  
  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "CollectionName" , 
                           "Name of the SimTrackerHit collection"  ,
                           _colName ,
                           std::string("FTDCollection") ) ;
       	
  registerProcessorParameter( "PointResolutionPixel" ,
                              "Point Resolution in FTD Pixel Disks (mm)"  ,
                              _pointResoPixel ,
                              (float)0.005) ;
      	
  registerProcessorParameter( "PointResolutionStripR" ,
                              "Point Resolution in FTD Strip Disks (mm)"  ,
                              _pointResoStrip ,
                              (float)0.007) ;

  registerProcessorParameter( "PointResolutionStripPhi" ,
                              "Point Resolution in FTD Strip Disks (mm)"  ,
                              _pointResoStripPhi ,
                              (float)0.007) ;
	
  registerProcessorParameter( "MomentumCutForDRays",
                              "Momentum Cut For D Rays (GeV)",
                              _momCut ,
                              float(10.0));
      
  registerProcessorParameter( "minIdStrip",
                              "Layer Id of the first strip disk",
                              _minIdStrip ,
                              int(0));


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


void NewFTDDigiProcessor::init() { 

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

void NewFTDDigiProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void NewFTDDigiProcessor::processEvent( LCEvent * evt ) { 

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
    
//             const int celId = SimTHit->getCellID0() ;
          
            const double *pos ;
            pos =  SimTHit->getPosition() ;  
            
            //-----> simTHit->getCellID0() <= 200+3 --> Pixel
            //  Strip or pixel?
            const int LayerId = _ftdLayers->layerID(  pos[2] );//SimTHit->getCellID0() ;
            double newPos[3];
            float pointReso;
            // _minIdStrip CLIC-ILD = 3
          
            //Smearing in RPhi for Strips
            if ( abs(LayerId) >= _minIdStrip )
              {
                
                pointReso = _pointResoStrip;
                double R = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
                double Phi  = atan2(pos[1], pos[0]);
                
                streamlog_out(DEBUG) << "Strip "    << LayerId 
                                     << " z "       << pos[2]
                                     << " R res "   << pointReso 
                                     << " Phi res " << _pointResoStripPhi 
                                     << std::endl;						
                
                double _pointResoPhi = asin(_pointResoStripPhi/sqrt( _pointResoStripPhi*_pointResoStripPhi + R*R));
                double rSmear = gsl_ran_gaussian(r, _pointResoStrip);
                double phiSmear = gsl_ran_gaussian(r, _pointResoPhi);
                
                double newR = R + rSmear;
                double newPhi = Phi + phiSmear;
                streamlog_out(DEBUG) << "R " << R << " -> " << newR 
                                     << " Phi " << Phi << " --> " << newPhi << endl; 
                
                newPos[0] = newR*cos( newPhi );
                newPos[1] = newR*sin( newPhi );
                
              }//Smearing for Pixels
            else if (  abs(LayerId) < _minIdStrip )
              {
                pointReso = _pointResoPixel;
                
                streamlog_out(DEBUG) << "Pixel " << LayerId 
                                     << " z "    << pos[2]
                                     << " Res "  << _pointResoPixel 
                                     << std::endl;						
                
                double xSmear = gsl_ran_gaussian(r,_pointResoPixel);
                double ySmear = gsl_ran_gaussian(r,_pointResoPixel);
                newPos[0] = pos[0] + xSmear;
                newPos[1] = pos[1] + ySmear;
              }
            else {
              streamlog_out(ERROR) << "Unexpected Error! The cell Identification is not inside [1,7]: " << abs(LayerId) << std::endl;
              exit(1);
            }
            
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

            //            trkHit->setType( 200+abs(SimTHit->getCellID0()));
            //fg: get layerID from z-position
            trkHit->setType( 201 + _ftdLayers->layerID(  pos[2] ) ) ;
            // FIXME:  need proper way of decoding the layerID in the TrackerHit
            //         for now we use this 'workaround' 201 is expected as FTD ID in 
            //         SiliconTracking

            //pointReso depends on either Strip or Pixel!
            float covMat[TRKHITNCOVMATRIX]={pointReso*pointReso,0.,pointReso*pointReso,0.,0.,0.};
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



void NewFTDDigiProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void NewFTDDigiProcessor::end(){ 

  gsl_rng_free(r);  
//   std::cout << "NewFTDDigiProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

  delete _ftdLayers ;
}

