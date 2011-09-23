/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "VTXBgClusters.h"
#include "VXDGeometry.h"
#include "VXDClusterParameters.h"

//#include <iostream>

#include <EVENT/LCCollection.h>
//#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
//#include <IMPL/TrackerHitImpl.h>
//#include <EVENT/MCParticle.h>
// #include "random.h"
// #include <CLHEP/Random/RandGauss.h>
//#include <gsl/gsl_randist.h>

#include <gear/GEAR.h>
#include <marlin/Global.h>
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#endif 

#include <cmath>
#include <math.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

VTXBgClusters aVTXBgClusters ;


VTXBgClusters::VTXBgClusters() : Processor("VTXBgClusters") {
  
  // modify processor description
  _description = "VTXBgClusters should create VTX TrackerHits from SimTrackerHits" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "RemoveDrays" ,
                              "Remove D-rays ?",
                              _removeDRays,
                              int(0));

  registerProcessorParameter( "MomentumCutForDRays" , 
                              "Momentum Cut For D Rays (MeV)",
                              _momCut ,
                              float(10.0));

  registerProcessorParameter( "Debug",
                              "Debugging option",
                              _debug,
                              int(0)); 

  registerProcessorParameter( "Modif",
                              "sensor modification option",
                              _mod,
                              int(0)); 
 
  // Input collections
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "VTXCollectionName" , 
                           "Name of the VTX SimTrackerHit collection"  ,
                           _colNameVTX ,
                           std::string("VXDCollection") ) ;

  //fg: not needed...  
//   // Output collections
//   registerOutputCollection( LCIO::TRACKERHIT,
//                            "VTXHitCollection" , 
//                            "Name of the vxd TrackerHit output collection"  ,
//                            _outColNameVTX ,
//                            std::string("VTXTrackerHits") ) ;

}


void VTXBgClusters::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  _vxdGeo = new VXDGeometry( Global::GEAR ) ;

  //fg: not needed
  //r = gsl_rng_alloc(gsl_rng_ranlxs2);

  //fg:  these should become processor parameters at some point:
  if (_mod==0){
    _epi=1;
    for(int i=0; i<6; i++) {
      _pitch[i]=25;
      _it[i]=200;
      _it[0]=50;
      _it[1]=50;
    }
    std::cout << "------------STANDARD mode for sensor--------------- " << std::endl;

  } else{
    _epi=0.015/0.05;
    for(int i=0; i<6; i++) {
      _pitch[i]=33;
      _pitch[0]=20;
      _pitch[1]=20;
      _it[i]=100;
      _it[0]=25;
      _it[1]=25;
    }
    std::cout << "------------MODIFIED mode for sensor--------------- " << std::endl;
  }
}

void VTXBgClusters::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void VTXBgClusters::processEvent( LCEvent * evt ) { 

  LCCollection* STHcol = 0 ;
  try{
    STHcol = evt->getCollection( _colNameVTX ) ;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG ) << "Collection " << _colNameVTX.c_str() << " is unavailable in event " 
                           << _nEvt << std::endl;

  }
  
  
  if( STHcol != 0 ){    
    
    //fg: not needed    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT )  ;
    
    int nSimHits = STHcol->getNumberOfElements()  ;
    
    for(int i=0; i< nSimHits; i++){
      
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      
      bool accept = 1;

      float totMomentum = 0;
      if (_removeDRays) { // check if hit originates from delta-electron 
        for (int i=0;i<3;++i) 
          totMomentum+=SimTHit->getMomentum()[i]*SimTHit->getMomentum()[i];
        totMomentum = sqrt(totMomentum);
        if (totMomentum < _momCut)
          accept = 0;
      }
      
      if (accept == 1) {

        /* the equation of the ellipse is
           y'**2/sma**2 + z'**2/SMA**2
           where y'-z' is the reference system which coincide with the axes of the ellipse 
           and with origin in the center of the ellipse.
           z' coincide with is the direction of the particle.
        */

        const int cellId = SimTHit->getCellID0() ;
        //fg: cellId is 'layer + 1' - we need layerId 
        const int layerId = ( cellId & 0xff ) - 1 ;

        //calculate the two axis of the ellipse that approximate the cluster
        double effpl = SimTHit->getPathLength()*_epi; //effective path length on the ladder surface

        double sma; //ellipse minor axis
        if (_mod==0)
          sma=3*_pitch[layerId];
        else sma=4*_pitch[layerId];

        double SMA = ((int)(effpl*1000/_pitch[layerId] -1/4) +4)*_pitch[layerId]; //ellipse major axis - effpl in mm
        sma*= (0.001 / 2.) ; //in mm and half axis
        SMA*= (0.001 / 2.) ; //in mm and half axis

 
        //calculate the ladder where the hit has taken place
        gear::Vector3D pos(SimTHit->getPosition()[0],SimTHit->getPosition()[1],SimTHit->getPosition()[2]) ;
        int ladderId = (_vxdGeo->getLadderID(pos,layerId)).second;

       streamlog_out( DEBUG ) << " eff. path length : " << effpl 
                               << " SMA " << SMA 
                               << " sma " << sma 
                               << " layerId :  " <<   layerId 
                               << " ladderId :  " <<   ladderId 
                               << std::endl ;


        //calculate the momentum of the hit in the ladder ref.syst.
        gear::Vector3D mom(SimTHit->getMomentum()[0],SimTHit->getMomentum()[1],SimTHit->getMomentum()[2]); 
        gear::Vector3D laddermom = _vxdGeo->lab2LadderDir(mom,layerId,ladderId) ;
        double modladdermom = sqrt(laddermom[1]*laddermom[1]+laddermom[2]*laddermom[2]);

        //calculate the direction of the SMA (which lays on the y-z plane)
        gear::Vector3D dSMAladder(0,SMA*laddermom[1]/modladdermom,SMA*laddermom[2]/modladdermom);

        //calculate the direction of the sma (orthogonal to SMA on the y-z plane)
        gear::Vector3D dsmaladder(0,sma*laddermom[2]/modladdermom,sma*laddermom[1]/modladdermom);

//         //calculate the direction of the SMA in the lab frame
//         gear::Vector3D dSMA = _vxdGeo->ladder2LabDir(dSMAladder,layerId,ladderId);
//         //calculate the direction of the sma in the lab frame
//         gear::Vector3D dsma = _vxdGeo->ladder2LabDir(dsmaladder,layerId,ladderId);
//        SimTHit->ext< ClusterParams >() = new VXDClusterParameters( dSMA , dsma , layerId, ladderId)  ;

        // compute the hit position in ladder coordinates:
        gear::Vector3D locPos = _vxdGeo->lab2LadderPos(pos,layerId,ladderId);

        SimTHit->ext< ClusterParams >() = new VXDClusterParameters(locPos, 
                                                                   dSMAladder, dsmaladder, 
                                                                   layerId,ladderId )  ;
        
      }      
      
    }
  }

  _nEvt ++ ;
}



void VTXBgClusters::check( LCEvent * evt ) { 

  // create some histograms with VXD hit densities and cluster sizes...

#ifdef MARLIN_USE_AIDA
  struct H1D{
    enum { 
      hitsLayer1,
      hitsLayer2,
      hitsLayer3,
      hitsLayer4,
      hitsLayer5,
      hitsLayer6,
      size 
    }  ;
  };

  if( isFirstEvent() ) {
    
    _hist1DVec.resize( H1D::size )   ;
    
    float  hitMax =  1000. ;
    _hist1DVec[ H1D::hitsLayer1 ] = AIDAProcessor::histogramFactory(this)->createHistogram1D( "hitsLayer1", 
											      "hits Layer 1", 
											      100, 0. ,hitMax ) ; 
    _hist1DVec[ H1D::hitsLayer2 ] = AIDAProcessor::histogramFactory(this)->createHistogram1D( "hitsLayer2", 
											      "hits Layer 2", 
											      100, 0. , hitMax ) ; 
    _hist1DVec[ H1D::hitsLayer3 ] = AIDAProcessor::histogramFactory(this)->createHistogram1D( "hitsLayer3", 
											      "hits Layer 3", 
											      100, 0. , hitMax ) ; 
    _hist1DVec[ H1D::hitsLayer4 ] = AIDAProcessor::histogramFactory(this)->createHistogram1D( "hitsLayer4", 
											      "hits Layer 4", 
											      100, 0. ,hitMax ) ; 
    _hist1DVec[ H1D::hitsLayer5 ] = AIDAProcessor::histogramFactory(this)->createHistogram1D( "hitsLayer5", 
											      "hits Layer 5", 
											      100, 0. , hitMax ) ; 
    _hist1DVec[ H1D::hitsLayer6 ] = AIDAProcessor::histogramFactory(this)->createHistogram1D( "hitsLayer6", 
											      "hits Layer 6", 
											      100, 0. , hitMax ) ; 


    _hist2DVec.resize( H1D::size )   ;
    
    for(int i=0 ; i < H1D::size ; ++i ){

      std::stringstream name ; 
      name << "clusSizLayer"  << i ;
      
       std::stringstream comment ;
       comment << "cluster sizes in layer " << i ;
      

      int nx = 200 ;
      float xMin = 0. ;
      float xMax = 4. ; 

      int ny =  200 ;
      float yMin = 0. ;
      float yMax = 4. ;
      
      
      _hist2DVec[i] = AIDAProcessor::histogramFactory(this)
        ->createHistogram2D( name.str() , 
                             comment.str() , 
                             nx, xMin, xMax ,
                             ny, yMin, yMax ) ; 
      
    }
    
    

  }
#endif



  LCCollection* vxdCol = 0 ; 

  int nhit = 0 ;
  int nHitL1 = 0 ;
  int nHitL2 = 0 ;
  int nHitL3 = 0 ;
  int nHitL4 = 0 ;
  int nHitL5 = 0 ;
  int nHitL6 = 0 ;

  
  try { 
    vxdCol = evt->getCollection( _colNameVTX ) ;

    int nH = vxdCol->getNumberOfElements() ;
    

    for(int i=0; i<nH ; ++i){
      
      SimTrackerHit* sth = dynamic_cast<SimTrackerHit*>(  vxdCol->getElementAt(i) ) ;
      
      int layer = ( sth->getCellID0() & 0xff ) - 1 ;
     
      // --- cluster axes in lab frame
      VXDClusterParameters* cluP = sth->ext< ClusterParams >() ;
      


      if( cluP != 0 ) { // physics hits have no cluster parameters

        // now get cluster axis vector:
        
        gear::Vector3D axisAlad = cluP->getClusterAxisA() ; 
        gear::Vector3D axisBlad = cluP->getClusterAxisB() ; 

        // fg: fixme this is not the exact projection of the rectangle....
        double za = std::abs( axisAlad.z() ) ;
        double zb = std::abs( axisBlad.z() ) ;
        double ya = std::abs( axisAlad.y() ) ;
        double yb = std::abs( axisBlad.y() ) ;

        double zPro = ( za > zb ?  za : zb ) ;
        double yPro = ( ya > yb ?  ya : yb ) ;

#ifdef MARLIN_USE_AIDA
        _hist2DVec[ layer ]->fill( zPro ,  yPro  ) ; 
#endif         
        


        //====== uncomment for  testing  the isPointInClusterEllipse method ============================
        /*
          gear::Vector3D pos(  sth->getPosition()[0] , sth->getPosition()[1] , sth->getPosition()[2] ) ; 
          gear::Vector3D ladPos  =  _vxdGeo->lab2LadderPos(pos , layer , cluP->getLadderId() ) ;
          
          if( ! cluP->isPointInClusterEllipse( ladPos ) )
          
          streamlog_out( ERROR ) << " point " << ladPos << "   is not in ellipse " 
          << " a: " << axisAlad << " b: " << axisBlad 
          << " at " << cluP->getClusterPosition()
          << std::endl ;
          
          if( cluP->isPointInClusterEllipse( ladPos + 1.01 * axisAlad ) )
          
          streamlog_out( ERROR ) << " point " << ladPos + 1.01 * axisAlad << "   is in ellipse " 
          << " a: " << axisAlad << " b: " << axisBlad 
          << " at " << cluP->getClusterPosition()
          << std::endl ;
        */ 
        //====== uncomment for  testing  the isPointInClusterEllipse method ============================
        

      }
      
      switch (layer) {
        
      case 0 :
        nHitL1++ ; 
        break ;      
      case 1 :
        nHitL2++ ; 
        break ;      
      case 2 :
        nHitL3++ ; 
        break ;      
      case 3 :
        nHitL4++ ; 
        break ;      
      case 4 :
        nHitL5++ ; 
        break ;      
      case 5 :
        nHitL6++ ; 
        break ;      
      }
      
    }
  } catch( DataNotAvailableException& e) {}
  
#ifdef MARLIN_USE_AIDA
  _hist1DVec[ H1D::hitsLayer1 ]->fill( nHitL1 ) ;
  _hist1DVec[ H1D::hitsLayer2 ]->fill( nHitL2 ) ;
  _hist1DVec[ H1D::hitsLayer3 ]->fill( nHitL3 ) ;
  _hist1DVec[ H1D::hitsLayer4 ]->fill( nHitL4 ) ;
  _hist1DVec[ H1D::hitsLayer5 ]->fill( nHitL5 ) ;
  _hist1DVec[ H1D::hitsLayer6 ]->fill( nHitL6 ) ;
#endif
}


void VTXBgClusters::end(){ 
  
  streamlog_out( DEBUG4 )  << "VTXBgClusters::end()  " << name() 
                           << " processed " << _nEvt << " events in " << _nRun << " runs "
                           << std::endl ;
}

