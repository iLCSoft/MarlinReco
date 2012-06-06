/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
//#ifdef USE_ROOT
#include "VTXNoiseClusters.h"
#include "VXDClusterParameters.h"

#include "VXDGeometry.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#endif 

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>

#include <iostream>
#include <sstream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/SimTrackerHitImpl.h>

// #include <CLHEP/Random/RandGauss.h>
#include <gsl/gsl_randist.h>

#include <cmath>
#include <math.h>

#include "TFile.h" 

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
	
  

  StringVec rootDefault(7) ;
  rootDefault[0] = "sap2VXD05_1.root" ;
  rootDefault[1] = "hst1" ;
  rootDefault[2] = "hst2" ;
  rootDefault[3] = "hst3" ;
  rootDefault[4] = "hst4" ;
  rootDefault[5] = "hst5" ;
  rootDefault[6] = "hst6" ;

  registerProcessorParameter( "RootHistograms" ,
                              "root file name and histogram names (one per layer)"  ,
                              _rootNames ,
                              rootDefault ) ;
	

  // Input collection
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "VTXCollectionName" , 
                           "Name of the VTX SimTrackerHit collection"  ,
                           _colNameVTX ,
                           std::string("VXDCollection") ) ;
  
  
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




  const gear::VXDParameters& gearVXD = Global::GEAR->getVXDParameters() ;
  const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout(); 
  
  if( (unsigned) layerVXD.getNLayers() != _rootNames.size() - 1 ){
    
    streamlog_out( ERROR  ) << " *************************************************** " << std::endl 
                            << " wrong number of histograms: " <<  _rootNames.size() 
                            << " for " << layerVXD.getNLayers() << " VXD layers "
                            << "  - do nothing ! " << std::endl 
                            << " *************************************************** " << std::endl  ;

  }


  const std::string& filename = _rootNames[0] ;
  
  _hfile = new TFile( filename.c_str() );

  if( ! _hfile->IsOpen() ) {
    
    std::string mess(" can't open root file :") ;
    mess +=  _rootNames[0] ;
    
    throw Exception( mess ) ;
  }
  
  streamlog_out(DEBUG4) << "Read histograms from file : " << filename << endl;
  
  
  //  _hfile->ls()  ;
  
  _hist.resize( layerVXD.getNLayers() ) ;
  
  for(int i=0;i<layerVXD.getNLayers(); ++i ){
    
    
    streamlog_out(DEBUG4) <<  "    " << _rootNames[i+1].c_str() << ":  .... "  ;
    
    _hist[i] = dynamic_cast<TH2F*> (  _hfile->Get( _rootNames[i+1].c_str()  ) ) ;
    
    if( !_hist[i] ) {
      std::string mess(" can't read  histogram :") ;
      mess +=  _rootNames[i+1] ;
      throw Exception( mess ) ;
    }
    streamlog_out(DEBUG4) <<  " OK  at :"  <<  _hist[i]  << std::endl ;

    streamlog_out(DEBUG) <<  " 2D histo - bins (x,y): "  
                         << _hist[i]->GetNbinsX() <<   ", "  
                         << _hist[i]->GetXaxis()->GetBinLowEdge(1) << ", "
                         << _hist[i]->GetXaxis()->GetBinLowEdge( _hist[i]->GetNbinsX() ) 
      +  _hist[i]->GetXaxis()->GetBinWidth(  _hist[i]->GetNbinsX() )  << " , " 
                         << _hist[i]->GetNbinsY() <<   ", "  
                         << _hist[i]->GetYaxis()->GetBinLowEdge(1) << ", "
                         << _hist[i]->GetYaxis()->GetBinLowEdge( _hist[i]->GetNbinsY() ) 
      +  _hist[i]->GetYaxis()->GetBinWidth(  _hist[i]->GetNbinsY() ) 
                               << std::endl ;
    
  }

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
        
        gear::Vector3D lab = _vxdGeo->ladder2LabPos( lad, i , j )  ;

        // double check if point is in sensitive
        if( ! gearVXD.isPointInSensitive( lab ) ){
          streamlog_out( WARNING ) << " created hit outside sensitve volume " << lab << std::endl ;
        }        
        
        SimTrackerHitImpl *hit = new SimTrackerHitImpl() ;

        double pos[3] = { lab[0] , lab[1] , lab[2]  } ;

        hit->setPosition( pos ) ;

        hit->setEDep( 0. ) ; // FIXME: which dedx should be used for noise hits 

        //FIXME: encode a proper cellID
        hit->setCellID(  i + 1  ) ; // fg: here we'd like to have ladder id as well ....


        // now we need to add some cluster parameters  to the hit :
        double cluZ, cluRPhi ;


       _hist[i]->GetRandom2( cluZ, cluRPhi ) ;
        

       // --- cluster axes in ladder frame
       gear::Vector3D axisAlad( 0, cluRPhi/2. , 0    ) ;
       gear::Vector3D axisBlad( 0,   0     , cluZ/2. ) ;
       
       // --- cluster axes in lab frame
       //        gear::Vector3D axisAlab = _vxdGeo->ladder2LabDir( axisAlad , i , j )  ;
       //        gear::Vector3D axisBlab = _vxdGeo->ladder2LabDir( axisBlad, i , j )  ;
       
       // hit position in ladder frame:

       
       hit->ext< ClusterParams >() = new VXDClusterParameters( lad, axisAlad , axisBlad , i , j )  ;


       col->addElement( hit ) ; 
       
      }
    }
  }


  _nEvt ++ ;
}



  void VTXNoiseClusters::check( LCEvent * evt ) { 


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
    
    float  hitMax =  100000. ;
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

      std::stringstream name("clusSizLayer") ;
      name << i ;
      
       std::stringstream comment("cluster sizes in layer ") ;
      comment << i ;
      

      int nx =  _hist[i]->GetNbinsX() ;
      float xMin = _hist[i]->GetXaxis()->GetBinLowEdge(1)  ;
      float xMax = _hist[i]->GetXaxis()->GetBinLowEdge( _hist[i]->GetNbinsX() ) 
        +  _hist[i]->GetXaxis()->GetBinWidth(  _hist[i]->GetNbinsX() ) ;

      int ny =  _hist[i]->GetNbinsY() ;
      float yMin = _hist[i]->GetYaxis()->GetBinLowEdge(1)  ;
      float yMax = _hist[i]->GetYaxis()->GetBinLowEdge( _hist[i]->GetNbinsY() ) 
        +  _hist[i]->GetYaxis()->GetBinWidth(  _hist[i]->GetNbinsY() ) ;
      
      
      _hist2DVec[i] = AIDAProcessor::histogramFactory(this)
        ->createHistogram2D( name.str() , 
                             comment.str() , 
                             nx, xMin, xMax ,
                             ny, yMin, yMax ) ; 
      
    }
    
    

  }
#endif



  LCCollection* vxdCol = 0 ; 

//  int nhit = 0 ;
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
      
      int layer = ( sth->getCellID0() & 0xff )  ;
     
      // --- cluster axes in lab frame
      VXDClusterParameters* cluP = sth->ext< ClusterParams >() ;
      

      //      gear::Vector3D pos(  sth->getPosition()[0] , sth->getPosition()[1] , sth->getPosition()[2] ) ; 

      if( cluP != 0 ) { // physics hits have no cluster parameters

        // now get cluster axis vector:

        gear::Vector3D axisAlab = cluP->getClusterAxisA() ; //- pos ;
        gear::Vector3D axisBlab = cluP->getClusterAxisB() ; //- pos ;
        
        double cluA = axisAlab.r() ;
        double cluB = axisBlab.r() ;
        
        //         streamlog_out( DEBUG ) << " axisAlab " << axisAlab 
        //                                << " axisBlab " << axisBlab  << std::endl ;
        
#ifdef MARLIN_USE_AIDA
        _hist2DVec[ layer-1 ]->fill( cluB *2. , cluA * 2. ) ;
#endif
      }
      
      switch (layer) {
        
      case 1 :
        nHitL1++ ; 
        break ;      
      case 2 :
        nHitL2++ ; 
        break ;      
      case 3 :
        nHitL3++ ; 
        break ;      
      case 4 :
        nHitL4++ ; 
        break ;      
      case 5 :
        nHitL5++ ; 
        break ;      
      case 6 :
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


void VTXNoiseClusters::end(){ 
  
  std::cout << "VTXNoiseClusters::end()  " << name() 
            << " processed " << _nEvt << " events in " << _nRun << " runs "
            << std::endl ;
}


//#endif // USE_ROOT
