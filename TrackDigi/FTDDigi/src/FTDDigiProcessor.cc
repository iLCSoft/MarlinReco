/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "FTDDigiProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>
#include <EVENT/MCParticle.h>
#include <math.h>

#include <gsl/gsl_randist.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

FTDDigiProcessor aFTDDigiProcessor ;


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
  
}

void FTDDigiProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void FTDDigiProcessor::processEvent( LCEvent * evt ) { 

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

          if (mom > _momCut) {
    
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
            trkHit->setType( SimTHit->getCellID());
            // FIXME: SJA this needs to set to FTD code
            // trkHit->setType( 100+celId ) ;
            float covMat[TRKHITNCOVMATRIX]={_pointReso,0.,0.,_pointReso,0.,0.};
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

}

