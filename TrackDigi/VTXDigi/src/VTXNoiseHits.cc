/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "VTXNoiseHits.h"

#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>

#include <IMPL/TrackerHitImpl.h>

// #include <CLHEP/Random/RandGauss.h>
#include <gsl/gsl_randist.h>

#include <cmath>
#include <math.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

VTXNoiseHits aVTXNoiseHits ;


VTXNoiseHits::VTXNoiseHits() : Processor("VTXNoiseHits") {
  
  // modify processor description
  _description = "VTXNoiseHits should create VTX TrackerHits from SimTrackerHits" ;
  
  
  // register steering parameters: name, description, class-variable, default value
  
  FloatVec densityDefault(5) ;
  for(int i=0 ;i<5; i++ ) 
    densityDefault[i] = 0.0 ;
  
  
  registerProcessorParameter( "HitDensityPerLayer_VTX" ,
                              "hit densities (hits/cm^2) per VXD layer"  ,
                              _densities ,
                              densityDefault ) ;
	
  
  // Input collection
  registerInputCollection( LCIO::TRACKERHIT,
                           "VTXCollectionName" , 
                           "Name of the VTX TrackerHit collection"  ,
                           _colNameVTX ,
                           std::string("VTXTrackerHits") ) ;
  
  
  registerProcessorParameter( "PointResolutionRPhi_VTX" ,
                              "R-Phi Resolution in VTX"  ,
                              _pointResoRPhiVTX ,
                              float(0.0027)) ;
	
  
  registerProcessorParameter( "PointResolutionZ_VTX" , 
                              "Z Resolution in VTX" ,
                              _pointResoZVTX ,
                              float(0.0027));
  
}


void VTXNoiseHits::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  // initialize gsl random generator
  // ranlux algorithm of Lüscher, which produces 'luxury random numbers'
  r = gsl_rng_alloc(gsl_rng_ranlxs2) ;

  //FIXME: do we want a seed from a processor parameter ?
}


void VTXNoiseHits::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void VTXNoiseHits::processEvent( LCEvent * evt ) { 

  LCCollection* col = 0 ;
  
  try{
    
    col = evt->getCollection( _colNameVTX ) ;
    
  }    
  catch(DataNotAvailableException &e){
    
    
    streamlog_out( WARNING ) << " VTX collection " << _colNameVTX  
                             << " not found - do nothing ! " << std::endl ;
    
    
    return ;
    
  }
  
  //get VXD geometry info
  const gear::VXDParameters& gearVXD = Global::GEAR->getVXDParameters() ;
  const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout(); 
  
  //gear::Vector3D hitvec(pos[0],pos[1],pos[2]);

  if( (unsigned) layerVXD.getNLayers() !=   _densities.size()  ){
    
    
    streamlog_out( ERROR  ) << " *************************************************** " << std::endl 
                            << " wrong number of hit densities: " <<  _densities.size() 
                            << " for " << layerVXD.getNLayers() << " VXD layers "
                            << "  - do nothing ! " << std::endl 
                            << " *************************************************** " << std::endl  ;
    
    return ;
    
  }
              
//  double gap = 	gearVXD.getShellGap () ;

  for( int i=0 ; i <  layerVXD.getNLayers() ; i++ ) {

    double 	phi0 = layerVXD.getPhi0 (i) ;    
    double 	dist = layerVXD.getSensitiveDistance (i) ;

    double 	thick = layerVXD.getSensitiveThickness (i) ;
    double 	offs =  layerVXD.getSensitiveOffset (i) ;
    double 	width = layerVXD.getSensitiveWidth (i) ;

    // -----fg:  gear length is half length really !!!!!
    double 	len =   2 * layerVXD.getSensitiveLength (i) ; 
      
 
    int nLad  = layerVXD.getNLadders(i) ;

    for( int j=0 ; j < nLad ; j++ ) {

      double phi = phi0 + j *  ( 2 * M_PI ) /  nLad  ; 

      // point in middle of sensitive ladder at z=-len/2
      gear::Vector3D p0( dist + thick / 2. ,  phi ,  - len / 2. ,   gear::Vector3D::cylindrical ) ;     

      // direction vector along ladder in rphi (negative) 
      gear::Vector3D v0( - (width / 2.  - offs ) ,  phi +  M_PI / 2. , 0. ,   gear::Vector3D::cylindrical ) ;     
      
      // point p1:   'lower left corner of sensitive surface' (seen from outside)
      gear::Vector3D p1 = p0 + v0 ;    
      
      // v1: direction along  rphi
      gear::Vector3D v1( width ,  phi +  M_PI / 2. , 0. ,   gear::Vector3D::cylindrical ) ;  

      // v2: direction along  ladder (z-axis)
      gear::Vector3D v2(  0. ,  0. ,  len   ) ;  
      

      double area = len * width / 100. ;  // mm^2 to cm^2


      int  nHit = gsl_ran_poisson( r,  _densities[ i ] * area  ) ;

      streamlog_out( DEBUG ) << " layer: " << i << " - ladder: " << j 
                             << " density: " <<  _densities[ i ] 
                             << " area: " <<   area
                             << " #hits: " << nHit 
                             << std::endl ;


      

      for( int k=0 ; k< nHit ; k++){

        double  l1 = gsl_rng_uniform( r ) ;
        double  l2 = gsl_rng_uniform( r ) ;
       
        gear::Vector3D  ph = p1 + l1 * v1 + l2 * v2 ;
        
        if( ! gearVXD.isPointInSensitive( ph ) ){
          
          streamlog_out( WARNING ) << " created hit outside sensitve volume " << ph << std::endl ;
        }        

        TrackerHitImpl *hit = new TrackerHitImpl() ;

        double pos[3] = { ph[0] , ph[1] , ph[2]  } ;

        hit->setPosition( pos ) ;

        hit->setEDep( 0. ) ; // FIXME: which dedx should be used for noise hits 

        hit->setType( 101 + i  ) ; // encoding used in VTXDigi: 100 + layernum + 1

        //fg: FIXME: this is not quite the proper error matrix for a ladder
        // and it is also not in cartesian coordiantes
        // however this is what is done in VTXDigiProcessor as well. ..
        float covMat[TRKHITNCOVMATRIX]={ 0., 0.,  _pointResoRPhiVTX * _pointResoRPhiVTX, 
                                         0., 0.,  _pointResoZVTX    * _pointResoZVTX } ;
        hit->setCovMatrix(covMat);      
        
//         const LCObjectVec& simHits  = hit->getRawHits() ;
//         streamlog_out( DEBUG ) << "   hit->getRawHits() " << simHits.size() << std::endl ;
//         for( LCObjectVec::const_iterator objIt = simHits.begin() ;
//              objIt != simHits.end() ; ++objIt ){
//         }
        
        col->addElement( hit ) ; 

      }
    }
  }

  _nEvt ++ ;
}



  void VTXNoiseHits::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void VTXNoiseHits::end(){ 
  
  std::cout << "VTXNoiseHits::end()  " << name() 
            << " processed " << _nEvt << " events in " << _nRun << " runs "
            << std::endl ;
}


