/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "VTXNoiseClusters.h"

#include "VXDGeometry.h"

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>

#include <iostream>


#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/SimTrackerHitImpl.h>

// #include <CLHEP/Random/RandGauss.h>
#include <gsl/gsl_randist.h>

#include <cmath>
#include <math.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

VTXNoiseClusters aVTXNoiseClusters ;


VTXNoiseClusters::VTXNoiseClusters() : Processor("VTXNoiseClusters") {
  
  // modify processor description
  _description = "VTXNoiseClusters adds SimTrackerHits with salt'n pepper noise clusters" ;
  
  
  // register steering parameters: name, description, class-variable, default value
  
  FloatVec densityDefault(6) ;
  for(int i=0 ;i<6; i++ ) 
    densityDefault[i] = 0.0 ;
  
  
  registerProcessorParameter( "HitDensityPerLayer_VTX" ,
                              "hit densities (hits/cm^2) per VXD layer"  ,
                              _densities ,
                              densityDefault ) ;
	
  
  // Input collection
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "VTXCollectionName" , 
                           "Name of the VTX SimTrackerHit collection"  ,
                           _colNameVTX ,
                           std::string("VXDCollection") ) ;
  
  
//   registerProcessorParameter( "PointResolutionRPhi_VTX" ,
//                               "R-Phi Resolution in VTX"  ,
//                               _pointResoRPhiVTX ,
//                               float(0.0027)) ;
//   registerProcessorParameter( "PointResolutionZ_VTX" , 
//                               "Z Resolution in VTX" ,
//                               _pointResoZVTX ,
//                               float(0.0027));
  
  registerProcessorParameter( "RandomSeed" , 
                              "random seed - default 42" ,
                              _ranSeed ,
                              int(42) ) ;

}


void VTXNoiseClusters::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  _vxdGeo = new VXDGeometry( Global::GEAR ) ;

  // initialize gsl random generator
  // ranlux algorithm of Lüscher, which produces 'luxury random numbers'
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2) ;

  gsl_rng_default_seed = _ranSeed ;

}


void VTXNoiseClusters::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void VTXNoiseClusters::processEvent( LCEvent * evt ) { 
}

void VTXNoiseClusters::modifyEvent( LCEvent * evt ) { 

  LCCollection* col = 0 ;
  
  try{
    
    col = evt->getCollection( _colNameVTX ) ;
    
  }    
  catch(DataNotAvailableException &e){
    
    
    streamlog_out( WARNING ) << " VTX collection " << _colNameVTX  
                             << " not found - do nothing ! " << std::endl ;
    
    
    return ;
    
  }
  
//   //fg ++++++++++++++ test and debug code ++++++++++++++++++++++
//   _vxdGeo->test() ;
  
//   if( col != 0 ){
    
//     for( int i=0 ; i < col->getNumberOfElements() ; ++i ){

//       SimTrackerHit* sth = dynamic_cast<SimTrackerHit*>( col->getElementAt(i) ) ; 
      
//       gear::Vector3D pos( sth->getPosition()[0], sth->getPosition()[1], sth->getPosition()[2]  ) ;

//       std::pair<int,int> id = _vxdGeo->getLadderID( pos ) ;

//       if( id.first < 0 ) {
        
//         streamlog_out( WARNING ) <<  " VTX hit outside sensitive : " 
//                                  <<  pos 
//                                  << " in gear: " 
//                                  << Global::GEAR->getVXDParameters().isPointInSensitive( pos ) 
//                                  << std::endl ;

        

//       } else {

//          gear::Vector3D lad = _vxdGeo->lab2Ladder( pos , id.first, id.second )  ;
//          streamlog_out( DEBUG4 ) << " %%% " << lad ;
//       }
//     } 
//   }


//   return ;
//   //fg ++++++++++++++ test and debug code ++++++++++++++++++++++
  

  //get VXD geometry info
  const gear::VXDParameters& gearVXD = Global::GEAR->getVXDParameters() ;
  const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout(); 
  

  if( (unsigned) layerVXD.getNLayers() !=   _densities.size()  ){
    
    
    streamlog_out( ERROR  ) << " *************************************************** " << std::endl 
                            << " wrong number of hit densities: " <<  _densities.size() 
                            << " for " << layerVXD.getNLayers() << " VXD layers "
                            << "  - do nothing ! " << std::endl 
                            << " *************************************************** " << std::endl  ;
    
    return ;
    
  }
              
  double hgap = 	gearVXD.getShellGap() / 2. ;

  for( int i=0 ; i <  layerVXD.getNLayers() ; i++ ) {
    
    double 	width = layerVXD.getSensitiveWidth (i) ;
    double 	len =   layerVXD.getSensitiveLength (i) ; 

    // area of one double ladder (+z and -z)
    double area = 2 * len * width / 100. ;  // mm^2 to cm^2
    
    
    int nLad  = layerVXD.getNLadders(i) ;
    
    for( int j=0 ; j < nLad ; j++ ) {

      int  nHit = gsl_ran_poisson( _rng ,  _densities[ i ] * area  ) ;

      streamlog_out( DEBUG ) << " layer: " << i << " - ladder: " << j 
                             << " density: " <<  _densities[ i ] 
                             << " area: " <<   area
                             << " #hits: " << nHit 
                             << std::endl ;

      for( int k=0 ; k< nHit ; k++){

        double  l1 = gsl_rng_uniform( _rng ) ;
        double  l2 = gsl_rng_uniform( _rng ) ;
       
        // ------ compute a point on the ladder 
        // - (x,y) origin is in the center of the ladder:
        double k1 = 0.5 - l1 ; 
        // - z origin is at z==0 in the lab

        double k2 =  ( (k % 2) ?  -1. : 1. ) * l2 ; 

        gear::Vector3D lad( 0 , k1 * width, k2 * len + hgap ) ;  
        
        gear::Vector3D lab = _vxdGeo->ladder2Lab( lad, i , j )  ;

        // double check if point is in sensitive
        if( ! gearVXD.isPointInSensitive( lab ) ){
          streamlog_out( WARNING ) << " created hit outside sensitve volume " << lab << std::endl ;
        }        
        
        SimTrackerHitImpl *hit = new SimTrackerHitImpl() ;

        double pos[3] = { lab[0] , lab[1] , lab[2]  } ;

        hit->setPosition( pos ) ;

        hit->setdEdx( 0. ) ; // FIXME: which dedx should be used for noise hits 

        //FIXME: encode a proper cellID
        hit->setCellID(  i + 1  ) ; // fg: here we'd like to have ladder id as well ....




        // TODO : now we need to add some cluster paramaters  to the hit :

        // ...

        // ...


        col->addElement( hit ) ; 

      }
    }
  }

  _nEvt ++ ;
}



  void VTXNoiseClusters::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void VTXNoiseClusters::end(){ 
  
  std::cout << "VTXNoiseClusters::end()  " << name() 
            << " processed " << _nEvt << " events in " << _nRun << " runs "
            << std::endl ;
}


