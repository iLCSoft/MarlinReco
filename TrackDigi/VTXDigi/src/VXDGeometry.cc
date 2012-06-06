/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "VXDGeometry.h"
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>

#include "streamlog/streamlog.h"
// #include <iostream>

#include <gsl/gsl_randist.h>

#include <cmath>
#include <math.h>

  
bool isEqual(double d0 , double d1 ) {
  
  static double EPSILON = 1e-12 ; 
  
  if( d1 == d0 ) 
    return true ;
  
  if( std::abs( d1 )  < EPSILON && std::abs( d0 ) < EPSILON  ) 
    return true ;
      
  if( std::abs( d0 - d1 )  < EPSILON ) // absolute difference small 
    return true ;
        
  // return relative difference smaller than epsilon 
  return ( 2. * std::abs ( d0 - d1 )  / (std::abs( d1 )  + std::abs( d0 )) ) < EPSILON ;


}
bool isEqual(gear::Vector3D v0 , gear::Vector3D v1 ) {
  
  return 
    isEqual(  v0[0] , v1[0] ) && 
    isEqual(  v0[1] , v1[1] ) && 
    isEqual(  v0[2] , v1[2] ) ;
}



VXDGeometry::VXDGeometry(gear::GearMgr* gearMgr) : 
  _gearMgr( gearMgr ) {
  
  init() ;

}

void VXDGeometry::init() { 

  // init VXD parameters ....
  const gear::VXDParameters& gearVXD = _gearMgr->getVXDParameters() ;
  const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout(); 
  
  
  streamlog_out( DEBUG0 ) << " initializing VXD ladder geometry from gear ... " 
                          << std::endl  ;
  
  _vxdLadders.resize( layerVXD.getNLayers() ) ; 
  _vxdLayers.resize(  layerVXD.getNLayers() ) ; 


  // shell gap  - same for all layers
  double  gap = gearVXD.getShellGap () ; 

  for( int i=0 ; i <  layerVXD.getNLayers() ; i++ ) {
    
    double 	phi0 = layerVXD.getPhi0 (i) ;    
    double 	dist = layerVXD.getSensitiveDistance (i) ;
    
    double 	thick = layerVXD.getSensitiveThickness (i) ;
    double 	offs  = layerVXD.getSensitiveOffset (i) ;
    double 	width = layerVXD.getSensitiveWidth (i) ;
    
    // -----fg:  gear length is half length really !!!!!
    double 	len =  layerVXD.getSensitiveLength (i) ; 

    
    int nLad  = layerVXD.getNLadders(i) ;
    
    //   double rMin ;
    //   double rMax ;
    //   double length ;
    //   double width ;
    //   double ladderArea ;
    //   int nLadders ;
    
    _vxdLayers[i].width = width ;
    _vxdLayers[i].length = layerVXD.getSensitiveLength (i) ;

    _vxdLayers[i].gap  = 	gap ;
    _vxdLayers[i].thickness  = 	thick ;

    _vxdLayers[i].ladderArea = 2 * width * len ;
    _vxdLayers[i].nLadders = nLad ;

    double locA = dist+thick ;
    double locB = width/2. + std::abs(offs) ;

    _vxdLayers[i].rMax = std::sqrt( locA*locA + locB*locB )  ;
    _vxdLayers[i].rMin = dist ;
    
    streamlog_out( DEBUG4 ) << " layer: " << i 
                            << " phi0 : " << phi0
                            << " offs : " << offs
                            << " width : " << width
                            << " length : " << _vxdLayers[i].length
                            << " ladderArea : " << _vxdLayers[i].ladderArea
                            << " rMin : " << _vxdLayers[i].rMin 
                            << " rMax : " << _vxdLayers[i].rMax 
                            << std::endl ;
    
    _vxdLadders[i].resize( nLad ) ;
    
    for( int j=0 ; j < nLad ; j++ ) {

      // rotation around z-axis
      double phi = phi0 + j *  ( 2 * M_PI ) /  nLad  ; 
      

      // translation vector
      gear::Vector3D t0( dist + thick / 2. , phi, 0., gear::Vector3D::cylindrical ) ;
      
      double sign = offs / std::abs( offs ) ;

      gear::Vector3D t1( std::abs( offs ),  
                         phi + sign * M_PI / 2.,
                         0,  // keep z coordinate
                         gear::Vector3D::cylindrical) ;
      
      VXDLadder& l = _vxdLadders[i][j] ;
      l.phi = phi ;
      l.trans = t0 + t1  ;

    }
  }
}


gear::Vector3D VXDGeometry::lab2LadderPos( gear::Vector3D labPos, int layerID, int ladderID) {
  
  gear::Vector3D u = labPos -  _vxdLadders[layerID][ladderID].trans ;
  
  return gear::Vector3D( u.rho(), 
                         u.phi() - _vxdLadders[layerID][ladderID].phi, 
                         u.z(),  
                         gear::Vector3D::cylindrical ) ;
}

gear::Vector3D VXDGeometry::ladder2LabPos( gear::Vector3D ladderPos, int layerID, int ladderID){

  gear::Vector3D r( ladderPos.rho(), 
                    ladderPos.phi() + _vxdLadders[layerID][ladderID].phi, 
                    ladderPos.z(),  
                    gear::Vector3D::cylindrical ) ;

  return r + _vxdLadders[layerID][ladderID].trans ;

}

gear::Vector3D VXDGeometry::lab2LadderDir( gear::Vector3D labDir, int layerID, int ladderID) {
  
  return gear::Vector3D( labDir.rho(), 
                         labDir.phi() - _vxdLadders[layerID][ladderID].phi, 
                         labDir.z(),  
                         gear::Vector3D::cylindrical ) ;
}

gear::Vector3D VXDGeometry::ladder2LabDir( gear::Vector3D ladderDir, int layerID, int ladderID){

  return gear::Vector3D( ladderDir.rho(), 
                         ladderDir.phi() + _vxdLadders[layerID][ladderID].phi, 
                         ladderDir.z(),  
                         gear::Vector3D::cylindrical ) ;
  
}
  

 
std::pair<int,int> VXDGeometry::getLadderID( gear::Vector3D labPos, int layerID ) {
  
  if( layerID < 0 ) { // need to search the layer first.... 

    double r = labPos.rho() ;
    unsigned n = _vxdLayers.size() ;
    
    for(unsigned i = 0 ; i < n ;  ++i ) { // loop over all layers
      
      if( _vxdLayers[i].rMin <= r && r <= _vxdLayers[i].rMax ) { // candidate layer
        
        //        bool isSensitive = false ;
        unsigned m = _vxdLadders[i].size() ;

        for(unsigned j = 0 ; j < m ;  ++j ) { // check 
          
          gear::Vector3D l = lab2LadderPos( labPos , i , j ) ;
          
          double z_abs = std::abs( l.z() ) - _vxdLayers[i].gap / 2. ;
          
          if( ( 0 < z_abs && z_abs < _vxdLayers[i].length         ) && 
              ( std::abs( l.y() )  < _vxdLayers[i].width / 2.     ) && 
              ( std::abs( l.x() )  < _vxdLayers[i].thickness / 2. )    ) {
            
            return std::make_pair(  i , j )  ;
          }
        }
      }
    }
    return std::make_pair(  -1, -1  )  ;
  }

  unsigned m = _vxdLadders[layerID].size() ;
  
  for(unsigned j = 0 ; j < m ;  ++j ) { // check 
    
    gear::Vector3D l = lab2LadderPos( labPos , layerID , j ) ;
    
    double z_abs = std::abs( l.z() ) - _vxdLayers[layerID].gap / 2. ;
    
    if( ( 0 < z_abs && z_abs < _vxdLayers[layerID].length         ) && 
        ( std::abs( l.y() )  < _vxdLayers[layerID].width / 2.     ) && 
        ( std::abs( l.x() )  < _vxdLayers[layerID].thickness / 2. )    ) {
      
      return std::make_pair( layerID , j )  ;
    }
  }
  return std::make_pair(  layerID , -1  )  ;
}
  



void VXDGeometry::test() {

  int nHit = 1000 ;

  // initialize gsl random generator
  // ranlux algorithm of Lüscher, which produces 'luxury random numbers'
  gsl_rng* r = gsl_rng_alloc(gsl_rng_ranlxs2) ;
  
  
  const gear::VXDParameters& gearVXD = _gearMgr->getVXDParameters() ;
  double hgap = gearVXD.getShellGap() / 2.  ;

  for(unsigned i=0 ; i< _vxdLayers.size() ; ++i) {
    
    for(unsigned j=0 ; j< _vxdLadders[i].size() ; ++j) {


      double width = _vxdLayers[i].width ;
      double len   = _vxdLayers[i].length ;
      
      // FIXME: what about the sensitive gap ?
      
      streamlog_out( MESSAGE ) << " testing  " << nHit << " random points for ladder " 
                               << j  << " of layer " << i << std::endl ;
      
      bool isOK = true  ;
      
      for( int k=0 ; k< nHit ; k++){
        
        double  l1 = gsl_rng_uniform( r ) ;
        double  l2 = gsl_rng_uniform( r ) ;
        
        // ------ compute a coordinate in the labframe 
        // - (x,y) origin is in the center of the ladder:
        double k1 = 0.5 - l1 ; 
        // - z origin is at z==0 in the lab
        double k2 =  ( -1. + (k % 2) * 2. ) * l2 ; // -1. k even , +1. k odd 

        
        gear::Vector3D ladRnd( 0 , k1 * width, k2 * len + hgap ) ;  
        
        gear::Vector3D lab = ladder2LabPos( ladRnd, i , j )  ;
        
        gear::Vector3D ladder = lab2LadderPos( lab, i , j )  ;
        


        if( ! isEqual( ladRnd , ladder )  ) {
          

          isOK = false ;

          streamlog_out( WARNING )   << " l1,l2: " << l1 << ", " << l2 << std::endl 
                                     << " ladRnd: " << ladRnd 
                                     << " lab   : " <<  lab
                                     << " ladder: " <<  ladder
                                     << std::endl ;
          

          streamlog_out( DEBUG4 ) <<  "isEqual(ladRnd[0],ladder[0]): "<< isEqual(ladRnd[0],ladder[0] ) 
                                  << std::endl  
                                  <<  "isEqual(ladRnd[1],ladder[1]): "<< isEqual(ladRnd[1],ladder[1] )
                                  << std::endl
                                  <<  "isEqual(ladRnd[2],ladder[2]): "<< isEqual(ladRnd[2],ladder[2] ) 
                                  << std::endl ;

        }
        
        
        if( ! gearVXD.isPointInSensitive( lab ) ){
          
          isOK = false ;
          
          streamlog_out( WARNING ) << " created hit outside sensitve volume " << lab << std::endl ;
        }        


        std::pair<int,int> id  =  getLadderID( lab ) ;
        if( (unsigned)id.first != i ){

          isOK = false ;
          
          streamlog_out( WARNING ) << "  getLadderID layer wrong : " << id.first
                                   << " instead of " << i << std::endl ;
        }        
        if( (unsigned)id.second != j ){

          isOK = false ;
          
          streamlog_out( WARNING ) << "  getLadderID ladder wrong : " << id.second
                                   << " instead of " << j << " layer " << i << std::endl ;
        }        
        


        if(  id != getLadderID( lab , i ) ){

          isOK = false ;
          
          streamlog_out( WARNING ) << "  getLadderID with ID=" <<i << " - ladder wrong : " 
                                   << id.second<< std::endl ;
        } 

        

        streamlog_out( DEBUG )  << " ladRnd: " << ladRnd 
                                << " lab   : " <<  lab
                                << " ladder: " <<  ladder
                                << std::endl ;
        
      }
      if (isOK ) 
        streamlog_out( MESSAGE ) << " tested  " << nHit << " random points for ladder " 
                                 << j  << " of layer " << i << " OK ! " << std::endl ;
      
    }
  } // ------ loop over layers

}
