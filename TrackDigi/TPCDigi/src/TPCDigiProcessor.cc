/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "TPCDigiProcessor.h"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include <gsl/gsl_randist.h>



#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>

//stl exception handler
#include <stdexcept>
#include "constants.h"
#include "voxel.h"

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
//

using namespace lcio ;
using namespace marlin ;
using namespace constants ;
using namespace std ;

TPCDigiProcessor aTPCDigiProcessor ;

// put this here as I dont't know how to get it into the tpc class or if that is even the write place for it
bool compare_phi( Voxel_tpc * a, Voxel_tpc * b){return ( a->getPhiIndex() < b->getPhiIndex() );}

TPCDigiProcessor::TPCDigiProcessor() : Processor("TPCDigiProcessor") 
{
  
  // modify processor description
  _description = "Produces TPC TrackerHit collection from SimTrackerHit collection, smeared in RPhi and Z" ;
  
  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "CollectionName" , 
                           "Name of the SimTrackerHit collection"  ,
                           _colName ,
                           std::string("STpc01_TPC") ) ;

  registerOutputCollection( LCIO::TRACKERHIT,
                            "TPCTrackerHitsCol" , 
                            "Name of the digitized TrackerHit collection"  ,
                           _TPCTrackerHitsCol ,
                            std::string("TPCTrackerHits") ) ;


  registerProcessorParameter( "PointResolutionRPhi" ,
                              "R-Phi Resolution constant in TPC"  ,
                              _pointResoRPhi ,
                               (float)0.160) ;

 registerProcessorParameter( "DiffusionCoeffRPhi" ,
                              "R-Phi Diffusion Co-efficent in TPC"  ,
                              _diffRPhi ,
                               (float)0.0) ;

  registerProcessorParameter( "PointResolutionZ" ,
                              "Z Resolution constant in TPC"  ,
                              _pointResoZ ,
                               (float)0.5) ;

 registerProcessorParameter( "PixZ" ,
                              "Defines spatial slice in Z"  ,
                              _pixZ ,
                               (float)1.4) ;

 registerProcessorParameter( "PixRP" ,
                              "Defines spatial slice in RP"  ,
                              _pixRP ,
                               (float)1.0) ;




}


void TPCDigiProcessor::init() 
{ 
  
  // From GNU documentation:
  // A replacement for the standard terminate_handler which prints 
  // more information about the terminating exception (if any) on stderr. Call ...
  std::set_terminate (__gnu_cxx::__verbose_terminate_handler);

  // usually a good idea to
  printParameters() ;

  //intialise random number generator 
  _random = gsl_rng_alloc(gsl_rng_ranlxs2);
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void TPCDigiProcessor::processRunHeader( LCRunHeader* run) 
{ 

  _nRun++ ;
} 

void TPCDigiProcessor::processEvent( LCEvent * evt ) 
{ 

  static bool firstEvent = true ;
  
  // this gets called for every event 
  // usually the working horse ...

//     int skipToEvent = 0;
//     if(_nEvt<skipToEvent) {
//       cout << "skipping event " << _nEvt << endl;
//      ++_nEvt;
//       return;
//     }

  
  if(firstEvent==true) cout << "TPCDigiProcessor called for first event" << endl;

  firstEvent = false ;

  LCCollection* STHcol = 0 ;
  try{
    STHcol = evt->getCollection( _colName ) ;
  }
  catch(DataNotAvailableException &e){
  }
  
  if( STHcol != 0 ){
  
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT )  ;
    
    int n_sim_hits = STHcol->getNumberOfElements()  ;

    //    cout << "number of SimHits = " << n_sim_hits << endl;
    // Assume initialy that there is no merging 
    
    //
    const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
    const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

    vector< vector <Voxel_tpc *> > tpcRowHits;

    // set size of row_hits to hold (n_rows) vectors
    tpcRowHits.resize(padLayout.getNRows());

    map< Voxel_tpc *,SimTrackerHit *> tpcHitMap;
  
    for(int i=0; i< n_sim_hits; i++){
      
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      
      double *pos;
      float de_dx;
      EVENT::MCParticle *mcp;
      
      pos = (double*) SimTHit->getPosition(); 
      de_dx = SimTHit->getdEdx();

      mcp = SimTHit->getMCParticle() ;

      //      cout << "x position for this hit is " << pos[0] << endl; 
      //      cout << "y position for this hit is " << pos[1] << endl; 
      //      cout << "z position for this hit is " << pos[2] << endl; 
      //      cout << "de/dx for this hit is " << de_dx << endl; 
      //      cout << "MCParticle PID for this hit is " << mcp->getPDG() << endl; 
      //      cout << "x =  " << x << endl; 

      //  SMEARING

       
      double aReso =_pointResoRPhi*_pointResoRPhi;
      double driftLenght = gearTPC.getMaxDriftLength() - fabs(pos[2]);
      if (driftLenght <0) { 
        std::cout << " TPCDigiProcessor : Warning! driftLenght < 0 " << driftLenght << " --> Check out your GEAR file!!!!" << std::endl; 
        std::cout << "Setting driftLenght to 0.1" << std::endl;
        driftLenght = 0.10;
      }
      double bReso = _diffRPhi*_diffRPhi;
      double tpcRPhiRes = sqrt(aReso + bReso*driftLenght);
 

//       std::cout << "_pointResoRPhi = " <<_pointResoRPhi << std::endl;
//       std::cout << "_diffRPhi = " << _diffRPhi << std::endl;
//       std::cout << "tpcRPhiRes = " << tpcRPhiRes << std::endl;
//       std::cout << "_pointResoZ = " << _pointResoZ << std::endl;

      double randrp = gsl_ran_gaussian(_random,tpcRPhiRes);
      double randz =  gsl_ran_gaussian(_random,_pointResoZ);

      // Make sure that the radius is equal to a pad radius
      // Get current hit radius

      double rad = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);

      //
      //      cout << "x position before smearing " << pos[0] << endl;     
      //      cout << "y position before smearing " << pos[1] << endl;     
      //      cout << "z position before smearing " << pos[2] << endl;     
      //

      pos[0] = pos[0] - randrp * pos[1]/rad;
      pos[1] = pos[1] + randrp * pos[0]/rad;
      pos[2] = pos[2] + randz;
      
      //
      //      cout << "x position after smearing " << pos[0] << endl;     
      //      cout << "y position after smearing " << pos[1] << endl;     
      //      cout << "z position after smearing " << pos[2] << endl;     
      //
      
      // At this point hits are mearly smeared now they must be digitised to trackerhits

      rad = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);

      //get phi of current hit
      float phi = atan2(pos[1],pos[0]);
      if (phi<0.) phi=phi+twopi;


      //get row index of current hit
      
      // modified to ues GEAR



      int padIndex = padLayout.getNearestPad(rad,phi);

      //      cout << "padIndex = " << padIndex << endl;

      const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;

      double gearRMin = planeExt[0] ;
      double gearRMax = planeExt[1] ;

      if(fabs(pos[2])>gearTPC.getMaxDriftLength()) pos[2] = (fabs(pos[2])/pos[2])*gearTPC.getMaxDriftLength();

      if(rad>gearRMax){
        pos[0] = gearRMax*cos(phi);
        pos[1] = gearRMax*sin(phi);
        cout << "Radius is greater than TPC RMax " << endl;
        cout << "rad = " << rad << endl; 
        cout << "the tpc OuterRadius = " << gearRMax << endl; 
      }

      if(rad<gearRMin){
        pos[0] = gearRMin*cos(phi);
        pos[1] = gearRMin*sin(phi);
        cout << "Radius is less than TPC RMin" << endl;
        cout << "rad = " << rad << endl; 
        cout << "the tpc InnerRadius = " << gearRMin << endl; 
      }
      
      int iRowHit = padLayout.getRowNumber(padIndex);


//       cout << "padLayout.getPadWidth(0) = " <<padLayout.getPadWidth(0) << endl;  
//       cout << "padLayout.getPadWidth(padIndex) = " <<padLayout.getPadWidth(padIndex) << endl;  
//       cout << "padLayout.getPadHeight(0) = " <<padLayout.getPadHeight(0) << endl;  
//       cout << "padLayout.getPadHeight(padIndex) = " <<padLayout.getPadHeight(padIndex) << endl;  

      //je: commented out next line as proposed by Kristian Harder
      //gear::Point2D padCoord = padLayout.getPadCenter(padIndex);

      //get phi index of current hit

      int iPhiHit = padLayout.getPadNumber(padIndex);

      int NumberOfTimeSlices =  (int) ((2.0 * gearTPC.getMaxDriftLength()) / _pixZ);

      //get z index of current hit

      int iZHit = (int) ( (float) NumberOfTimeSlices * 
                         ( gearTPC.getMaxDriftLength() + pos[2] ) / ( 2.0 * gearTPC.getMaxDriftLength() ) ) ;


      if(iZHit<0) iZHit=0;
      if(iZHit>NumberOfTimeSlices) iZHit=NumberOfTimeSlices;

      
      Voxel_tpc * atpcVoxel = new Voxel_tpc(iRowHit,iPhiHit,iZHit, pos, de_dx);

      tpcRowHits.at(iRowHit).push_back(atpcVoxel);
      
      tpcHitMap[atpcVoxel] = SimTHit; 


      //      cout << "a voxel hit "<< i << " has been added to row " << iRowHit << endl;  

    }    

    cout << "finished looping over simhits" << endl;
    // Add background hits here
    
    vector <Voxel_tpc *> row_hits;

    
    //      cout << "get the row hits" << endl;
    
    for (int i = 0; i<padLayout.getNRows(); ++i){
      
      row_hits = tpcRowHits[i];

      //      sort(row_hits.begin(), row_hits.end(), compare_phi );

//       cout << "row = " << i << "  row_hits.size() = " << row_hits.size() << endl;

      for (unsigned int j = 0; j<row_hits.size(); ++j){
        
        //        cout << "got the row hit j = " << j << "  row_hits.size() = " << row_hits.size() << endl;
        
        for (unsigned int k = j+1; k<row_hits.size(); k++){

          //          cout << "got the row hits k = " << k << endl;
          
          if(row_hits[j]->getPhiIndex() > row_hits[k]->getPhiIndex()){
            //            cout << "phi index of hit ["<<j<<"] = " << row_hits[j]->getPhiIndex() << endl; 
            //            cout << "phi index of hit ["<<k<<"] = " << row_hits[k]->getPhiIndex() << endl; 
          }
          

          
          if(abs(row_hits[j]->getZIndex()-row_hits[k]->getZIndex())<=1){

              if(abs(row_hits[j]->getPhiIndex()-row_hits[k]->getPhiIndex())<=1||abs(row_hits[j]->getPhiIndex()-row_hits[k]->getPhiIndex())==(int)padLayout.getPadsInRow(i).size()){

                /*              
                 cout << "*&^*&*&*&*&*&*&*&**&*&*&*&*&&*&*&**&*&*&" << endl;
                 cout << "double hit candidate found in row: " << i <<  endl;

                 cout << "row_hits[j]->getZ() " << row_hits[j]->getZ() << endl;
                 cout << "row_hits[k]->getZ() " << row_hits[k]->getZ() << endl;
                

                 cout << "row_hits dX^2 + dY^2 " << 
                   ( ( ( row_hits[j]->getX() - row_hits[k]->getX() ) * 
                       ( row_hits[j]->getX() - row_hits[k]->getX() ) ) +
                     ( ( row_hits[j]->getY() - row_hits[k]->getY() ) * 
                       ( row_hits[j]->getY() - row_hits[k]->getY() ) ) )
                      << endl;
                 cout << "padLayout.getPadWidth(0)^2 " << padLayout.getPadWidth(0) * padLayout.getPadWidth(0)
                      << endl;
                */

                if(fabs( row_hits[j]->getZ() - row_hits[k]->getZ() ) < _pixZ ) {

                  if((((row_hits[j]->getX()-row_hits[k]->getX())*(row_hits[j]->getX()-row_hits[k]->getX()))
                      +((row_hits[j]->getY()-row_hits[k]->getY())*(row_hits[j]->getY()-row_hits[k]->getY())))
//                      <  padLayout.getPadWidth(0) *  padLayout.getPadWidth(0) ){
// FIXME: SJA: the function getPadWidth does not return 2.2mm as set in gear_ldc.xml so use hard coded number for now
                     <  2.2 *  2.2 ){


                    //SimTrackerHit* Hit1 = tpcHitMap[row_hits[j]];
                    //SimTrackerHit* Hit2 = tpcHitMap[row_hits[k]];


                    /*
                    
                    cout << "double hit found in row: " << i 
                         << "    MCP1 : " << Hit1->getMCParticle()->id() 
                         << " (" <<  Hit1->getMCParticle()->getPDG() << ")   " 
                         << "Energy : " << Hit1->getMCParticle()->getEnergy() << "     " 
                         << "MCP2 : " << Hit2->getMCParticle()->id() 
                         << " (" <<  Hit1->getMCParticle()->getPDG() << ")   " 
                         << "Energy : " << Hit2->getMCParticle()->getEnergy()
                         << endl;
                    */


                    /*
                    float xd = (float)row_hits[j]->getX();
                    float yd = (float)row_hits[j]->getY();
                    float zd = (float)row_hits[j]->getZ();
	    
                    float xd2 = (float)row_hits[k]->getX();
                    float yd2 = (float)row_hits[k]->getY();
                    float zd2 = (float)row_hits[k]->getZ();
                    */

                    row_hits[j]->setAdjacent(row_hits[k]);
                    row_hits[k]->setAdjacent(row_hits[j]);
                  }		  
                }
              }
          }
        }
        
        // FIXME:SJA: 
        // 	At this point the double hits have been identified and at present they are not added to the
        // 	tracker hit collection. What should be done with them will be decided later.
        
        if(row_hits[j]->getNumberOfAdjacent()==0){
          //store hit variables
          TrackerHitImpl* trkHit = new TrackerHitImpl ;
          double pos[3] = {row_hits[j]->getX(),row_hits[j]->getY(),row_hits[j]->getZ()}; 
          trkHit->setPosition(pos);
          trkHit->setdEdx(row_hits[j]->getdEdx());
          trkHit->setType( 500 );
         
          
          double aReso =_pointResoRPhi*_pointResoRPhi;
          double driftLenght = gearTPC.getMaxDriftLength() - fabs(pos[2]);
          if (driftLenght <0) { 
            std::cout << " TPCDigiProcessor : Warning! driftLenght < 0 "  
                      << driftLenght << " --> Check out your GEAR file!!!!" << std::endl; 
            std::cout << "Setting driftLenght to 0.1" << std::endl;
            driftLenght = 0.10;
          }
          double bReso = _diffRPhi*_diffRPhi;
          double tpcRPhiRes = sqrt(aReso + bReso*driftLenght);
          float covMat[TRKHITNCOVMATRIX]={0.,0.,float(tpcRPhiRes*tpcRPhiRes),0.,0.,float(_pointResoZ*_pointResoZ)};
          trkHit->setCovMatrix(covMat);      
          
          if(pos[0]*pos[0]+pos[1]*pos[1]>0.0 * 0.0){
            // 	  push back the SimTHit for this TrackerHit
            trkHit->rawHits().push_back( tpcHitMap[row_hits[j]] );                        
            trkhitVec->addElement( trkHit ); 
          }

        }
      }
    }
    
    cout << "finished row hits" << endl;
    
    // set the parameters to decode the type information in the collection
    // for the time being this has to be done manually
    // in the future we should provide a more convenient mechanism to 
    // decode this sort of meta information
    StringVec typeNames ;
    IntVec typeValues ;
    typeNames.push_back( LCIO::TPCHIT ) ;
    typeValues.push_back( 1 ) ;
    trkhitVec->parameters().setValues("TrackerHitTypeNames" , typeNames ) ;
    trkhitVec->parameters().setValues("TrackerHitTypeValues" , typeValues ) ;
    
    evt->addCollection( trkhitVec , _TPCTrackerHitsCol ) ;
    
    
    for (int i = 0; i<padLayout.getNRows(); ++i){
      
      vector <Voxel_tpc *> current_row = tpcRowHits[i];
      
      for (unsigned int j = 0; j<current_row.size(); ++j){
        
        delete current_row[j];
      }
    }    
  }
  _nEvt++;
  
}



void TPCDigiProcessor::check( LCEvent * evt ) 
{ 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TPCDigiProcessor::end()
{ 

  gsl_rng_free(_random);
  //   cout << "TPCDigiProcessor::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << endl ;
  
}
