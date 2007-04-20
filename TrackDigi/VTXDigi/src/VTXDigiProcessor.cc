/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "VTXDigiProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>
#include <EVENT/MCParticle.h>
// #include "random.h"
// #include <CLHEP/Random/RandGauss.h>
#include <gsl/gsl_randist.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

VTXDigiProcessor aVTXDigiProcessor ;


VTXDigiProcessor::VTXDigiProcessor() : Processor("VTXDigiProcessor") {
  
  // modify processor description
  _description = "VTXDigiProcessor should create VTX TrackerHits from SimTrackerHits" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "PointResolutionRPhi_VTX" ,
                              "R-Phi Resolution in VTX"  ,
                              _pointResoRPhi_VTX ,
                               float(0.0040)) ;
	
  registerProcessorParameter( "PointResolutionZ_VTX" , 
                              "Z Resolution in VTX" ,
                              _pointResoZ_VTX ,
                              float(0.0040));

  registerProcessorParameter( "PointResolutionRPhi_SIT" ,
                              "R-Phi Resolution in SIT"  ,
                              _pointResoRPhi_SIT ,
                               float(0.010)) ;
	
  registerProcessorParameter( "PointResolutionZ_SIT" , 
                              "Z Resolution in SIT" ,
                              _pointResoZ_SIT ,
                              float(0.010));

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
 
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "VTXCollectionName" , 
                           "Name of the SimTrackerHit collection"  ,
                           _colNameVTX ,
                           std::string("vxd00_VXD") ) ;
  
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "SITCollectionName" , 
                           "Name of the SimTrackerHit collection"  ,
                           _colNameSIT ,
                           std::string("sit00_SIT") ) ;
  
  registerOutputCollection( LCIO::TRACKERHIT,
                           "VTXHitCollection" , 
                           "Name of the vxd TrackerHit output collection"  ,
                           _outColNameVTX ,
                           std::string("VTXTrackerHits") ) ;


  registerOutputCollection( LCIO::TRACKERHIT,
                            "SITHitCollection" , 
                            "Name of the sit TrackerHit output collection"  ,
                            _outColNameSIT ,
                            std::string("SITTrackerHits") ) ;
  
}


void VTXDigiProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  r = gsl_rng_alloc(gsl_rng_ranlxs2);
  
}

void VTXDigiProcessor::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void VTXDigiProcessor::processEvent( LCEvent * evt ) { 

  for (int iColl=0;iColl<2;++iColl) {

    LCCollection* STHcol = 0 ;
    try{
      if (iColl==0)
        STHcol = evt->getCollection( _colNameVTX ) ;
      else 
        STHcol = evt->getCollection( _colNameSIT ) ;
    }
    catch(DataNotAvailableException &e){
      if (_debug == 1) {
        if (iColl==0)
          std::cout << "Collection " << _colNameVTX.c_str() << " is unavailable in event " << _nEvt << std::endl;
        else 
          std::cout << "Collection " << _colNameSIT.c_str() << " is unavailable in event " << _nEvt << std::endl;
      }
    }

  
    if( STHcol != 0 ){    
    
      LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT )  ;
    
      int nSimHits = STHcol->getNumberOfElements()  ;
    
      for(int i=0; i< nSimHits; i++){
      
        SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;

        bool accept = 1;
        if (_removeDRays) { // check if hit originates from delta-electron 
          float totMomentum = 0;
          for (int i=0;i<3;++i) 
            totMomentum+=SimTHit->getMomentum()[i]*SimTHit->getMomentum()[i];
          totMomentum = sqrt(totMomentum);
          if (totMomentum < _momCut)
            accept = 0;
        }
        
        if (accept == 1) {
          
          const int celId = SimTHit->getCellID() ;
          
          const double *pos ;
          pos =  SimTHit->getPosition() ;  

          double xSmear;
          double zSmear;
          
          if (iColl==0) {        
            xSmear = gsl_ran_gaussian(r,_pointResoRPhi_VTX);
            zSmear = gsl_ran_gaussian(r,_pointResoZ_VTX);
            _pointResoRPhi = _pointResoRPhi_VTX;
            _pointResoZ = _pointResoZ_VTX;
          }
          else {
            xSmear = gsl_ran_gaussian(r,_pointResoRPhi_SIT);
            zSmear = gsl_ran_gaussian(r,_pointResoZ_SIT);
            _pointResoRPhi = _pointResoRPhi_SIT;
            _pointResoZ = _pointResoZ_SIT; 
          }


          double newPos[3] ;
          double phi = atan2(pos[1],pos[0]);
          double rad = sqrt(pos[1]*pos[1]+pos[0]*pos[0]);
          double phi_new = phi + xSmear/rad;
          newPos[0] = rad*cos(phi_new);
          newPos[1] = rad*sin(phi_new);
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
          if (iColl==0) 
            trkHit->setType(100+celId ); 
          else
            trkHit->setType(400+celId);
          float covMat[TRKHITNCOVMATRIX]={0.,0.,_pointResoRPhi*_pointResoRPhi,0.,0.,_pointResoZ*_pointResoZ};
          trkHit->setCovMatrix(covMat);      
          // 	  push back the SimTHit for this TrackerHit
          trkHit->rawHits().push_back( SimTHit ) ;
          trkhitVec->addElement( trkHit ) ; 
        }      

      }
      if (iColl==0) 
         evt->addCollection( trkhitVec , _outColNameVTX ) ;
      else 
        evt->addCollection( trkhitVec , _outColNameSIT) ;
    }
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

