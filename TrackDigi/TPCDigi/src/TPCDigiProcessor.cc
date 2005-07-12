/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "TPCDigiProcessor.h"
#include <iostream>

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
//romans random class for normal distribution
//fixme the random number should have the default value set as 1
#include "random.h"
#include "constants.h"
#include "voxel.h"
#include "tpc.h"

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
//

using namespace lcio ;
using namespace marlin ;
using namespace constants ;

TPCDigiProcessor aTPCDigiProcessor ;

// put this here as I dont't know how to get it into the tpc class or if that is even the write place for it
bool compare_phi( Voxel_tpc * a, Voxel_tpc * b){return ( a->getPhiIndex() < b->getPhiIndex() );};

TPCDigiProcessor::TPCDigiProcessor() : Processor("TPCDigiProcessor") 
{
  
  // modify processor description
  _description = "Produces TPC TrackerHit collection from SimTrackerHit collection, smeared in RPhi and Z" ;
  
  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "CollectionName" , 
			      "Name of the SimTrackerHit collection"  ,
			      _colName ,
			      std::string("tpc04_TPC") ) ;
}


void TPCDigiProcessor::init() 
{ 
  
  // From GNU documentation:
  // A replacement for the standard terminate_handler which prints 
  // more information about the terminating exception (if any) on stderr. Call ...
  std::set_terminate (__gnu_cxx::__verbose_terminate_handler);

  // usually a good idea to
  printParameters() ;
  
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
  
  if(firstEvent==true) std::cout << "TPCDigiProcessor called for first event" << std::endl;

  firstEvent = false ;
  
  LCCollection* STHcol = evt->getCollection( _colName ) ;
  
  if( STHcol != 0 ){
  
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT )  ;
    
    int n_sim_hits = STHcol->getNumberOfElements()  ;

    //    cout << "number of SimHits = " << n_sim_hits << endl;
    // Assume initialy that there is no merging 
    
    //
    const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;


  
    for(int i=0; i< n_sim_hits; i++){
      
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      
      double *pos;
      float de_dx;
      EVENT::MCParticle *mcp;
      
      pos = (double*) SimTHit->getPosition(); 
      de_dx = SimTHit->getdEdx();

      mcp = SimTHit->getMCParticle() ;

      //      std::cout << "x position for this hit is " << pos[0] << std::endl; 
      //      std::cout << "y position for this hit is " << pos[1] << std::endl; 
      //      std::cout << "z position for this hit is " << pos[2] << std::endl; 
      //      std::cout << "de/dx for this hit is " << de_dx << std::endl; 
      //      std::cout << "MCParticle PID for this hit is " << mcp->getPDG() << std::endl; 
      //      std::cout << "x =  " << x << std::endl; 

      //  SMEARING


      double tpcRPhiResMax = gearTPC.getDoubleVal("tpcRPhiResMax");
      double tpcRPhiRes = tpcRPhiResMax-fabs(pos[2])/gearTPC.getMaxDriftLength()*0.10;
      double tpcZRes = gearTPC.getDoubleVal("tpcZRes");

      RandomNumberGenerator RandomNumber;

      //      float randrp = the_tpc->getTpcRphiResMax() * (*(RandomNumber.Gauss(1.0)));
      //      float randz = the_tpc->getTpcZRes() * (*(RandomNumber.Gauss(1.0)));
      double randrp = tpcRPhiRes * (*(RandomNumber.Gauss(1.0)));
      double randz = tpcZRes * (*(RandomNumber.Gauss(1.0)));

      // Make sure that the radius is equal to a pad radius
      // Get current hit radius

      double rad = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);

      //
      //      std::cout << "x position before smearing " << pos[0] << std::endl;     
      //      std::cout << "y position before smearing " << pos[1] << std::endl;     
      //      std::cout << "z position before smearing " << pos[2] << std::endl;     
      //

      pos[0] = pos[0] - randrp * pos[1]/rad;
      pos[1] = pos[1] + randrp * pos[0]/rad;
      pos[2] = pos[2] + randz;
      
      //
      //      std::cout << "x position after smearing " << pos[0] << std::endl;     
      //      std::cout << "y position after smearing " << pos[1] << std::endl;     
      //      std::cout << "z position after smearing " << pos[2] << std::endl;     
      //
      
      // At this point hits are mearly smeared now they must be digitised to trackerhits

      rad = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);

      //get phi of current hit
      float phi = atan2(pos[1],pos[0]);
      if (phi<0.) phi=phi+twopi;


      //get row index of current hit
      
      //      int irow_hit = (int)((rad-the_tpc->getInnerRadius())/the_tpc->getPadRPitch());

      // modified to ues GEAR

      const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

      const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;

      double gearRMin = planeExt[0] ;
      double gearRMax = planeExt[1] ;

      // this needs to be looked at within the context of GEAR i.e. getPadNumber()
      int irow_hit = (int)((rad- ( gearRMin ))/the_tpc->getPadRPitch());


      //

      if(irow_hit<0.) {
	cout << "row index of hit less than zero : irow = " << irow_hit << endl;
	cout << "rad = " << rad << endl; 
	cout << "the_tpc->getInnerRadius() = " << the_tpc->getInnerRadius() << endl; 
	cout << "the_tpc->getPadRPitch() = " << the_tpc->getPadRPitch() << endl; 

	irow_hit = 0;

      }

      if(irow_hit>the_tpc->getNumberOfRows()-1){
	cout << "row index of hit greater than number of pad rows: irow = " << irow_hit << endl;
	cout << "rad = " << rad << endl; 
	cout << "the_tpc->getInnerRadius() = " << the_tpc->getInnerRadius() << endl; 
	cout << "the_tpc->getPadRPitch() = " << the_tpc->getPadRPitch() << endl; 

	irow_hit = the_tpc->getNumberOfRows()-1;
      }

      //get row radius needed to calculate the number of pixels in row
      //      float row_radius = twopi*((the_tpc->getPadRPitch()/2.)+(float)irow_hit*the_tpc->getPadRPitch());
      //      int number_of_pixels_in_row = int(row_radius / the_tpc->getPixelPhiWidth());

      //      cout << "about to get number of pixels for row = " << irow_hit << endl;
      int number_of_pixels_in_row = the_tpc->getNumberOfPhiSegments(irow_hit);
      //      cout << "got number of pixels rows" << endl;

      //get phi index of current hit
      int iphi_hit = (int)((number_of_pixels_in_row * phi) / twopi);

      //get z index of current hit
      int iz_hit = (int)((float)(the_tpc->getNumberOfTimeSlices())*
			 ((the_tpc->getHalfLength()+pos[2])/(2.*the_tpc->getHalfLength())));
      
      if(iz_hit<0) iz_hit=0;
      if(iz_hit>the_tpc->getNumberOfTimeSlices()) iz_hit=the_tpc->getNumberOfTimeSlices();

      
      Voxel_tpc * a_tpc_voxel = new Voxel_tpc(irow_hit,iphi_hit,iz_hit, pos, de_dx);
      
      the_tpc->putRowHit(irow_hit, a_tpc_voxel);
      
      the_tpc->putSimTrackerHit(a_tpc_voxel , SimTHit);

      //      cout << "a voxel hit "<< i << " has been added to row " << irow_hit << endl;  

    }    

    cout << "finished looping over simhits" << endl;

    // Add background hits here
    
    vector <Voxel_tpc *> row_hits;
    
    //      cout << "get the row hits" << endl;
    
    for (int j = 0; j<the_tpc->getNumberOfRows(); j++){
      
      row_hits = the_tpc->getRowHits(j);

      //      sort(row_hits.begin(), row_hits.end(), compare_phi );

      //      cout << "row = " << j << "  row_hits.size() = " << row_hits.size() << endl;

      //       if(row_hits.size()==1){
      //       //store hit variables
      //       TrackerHitImpl* trkHit = new TrackerHitImpl ;
      //       double pos[3] = {row_hits[0]->getX(),row_hits[0]->getY(),row_hits[0]->getZ()}; 
      //       trkHit->setPosition(pos);
      //       trkHit->setdEdx(row_hits[0]->getdEdx());
      //       trkHit->setType( 500 );
      
      //       // 	  push back the SimTHit for this TrackerHit
      //       trkHit->rawHits().push_back( the_tpc->getSimTrackerHit(row_hits[0]) );

      //       trkhitVec->addElement( trkHit ); 
      //       }

      for (int i = 0; i<row_hits.size(); i++){
	
	//	cout << "got the row hits i = " << i << "  row_hits.size() = " << row_hits.size() << endl;
	
	for (int k = i+1; k<row_hits.size(); k++){

	  //	  if(row_hits[i]->getPhiIndex() > row_hits[k]->getPhiIndex()){
	  //	    cout << "phi index of hit ["<<i<<"] = " << row_hits[i]->getPhiIndex() << endl; 
	  //	    cout << "phi index of hit ["<<k<<"] = " << row_hits[k]->getPhiIndex() << endl; 
	  //	  }

	  //	  cout << "got the row hits k = " << k << endl;
	  
	  if(abs(row_hits[i]->getZIndex()-row_hits[k]->getZIndex())<=1){
	    if(abs(row_hits[i]->getPhiIndex()-row_hits[k]->getPhiIndex())<=1||abs(row_hits[i]->getPhiIndex()-row_hits[k]->getPhiIndex())==the_tpc->getNumberOfPhiSegments(j)){

	      
	      //	      cout << "*&^*&*&*&*&*&*&*&**&*&*&*&*&&*&*&**&*&*&" << endl;
	      //	      cout << "double hit candidate found in row: " << j <<  endl;
	      
	      if(fabs(row_hits[i]->getZ()-row_hits[k]->getZ())<the_tpc->getTimeWidth()){
		if((((row_hits[i]->getX()-row_hits[k]->getX())*(row_hits[i]->getX()-row_hits[k]->getX()))
		    +((row_hits[i]->getY()-row_hits[k]->getY())*(row_hits[i]->getY()-row_hits[k]->getY())))
		   <the_tpc->getPixelPhiWidth()*the_tpc->getPixelPhiWidth()){
		  
		  cout << "double hit found in row: " << j <<  endl;
		  row_hits[i]->setAdjacent(row_hits[k]);
		  row_hits[k]->setAdjacent(row_hits[i]);
		}		  
	      }
	    }
	  }
	}

	// FIXME:SJA: 
	// 	At this point the double hits have been identified and at present they are not added to the
	// 	tracker hit collection. What should be done with them will be decided later.

	if(row_hits[i]->getNumberOfAdjacent()==0){
	  //store hit variables
	  TrackerHitImpl* trkHit = new TrackerHitImpl ;
	  double pos[3] = {row_hits[i]->getX(),row_hits[i]->getY(),row_hits[i]->getZ()}; 
	  trkHit->setPosition(pos);
	  trkHit->setdEdx(row_hits[i]->getdEdx());
	  trkHit->setType( 500 );
	  
	  // 	  push back the SimTHit for this TrackerHit
	  trkHit->rawHits().push_back( the_tpc->getSimTrackerHit(row_hits[i]) );

	  trkhitVec->addElement( trkHit ); 

	}
      }
    }
    
    cout << "finished row hits" << endl;
    
    //    for (int j = 0; j<the_tpc->getNumberOfRows(); j++){
      
    //      cout << "row number = " << j << endl;

    //      row_hits = the_tpc->getRowHits(j);
      
    //      for (int i = 0; i<row_hits.size(); i++){
	
    //	cout << "number of adjacents = " << row_hits[i]->getNumberOfAdjacent() << endl;	  

    //      }
    //    }

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
    
    evt->addCollection( trkhitVec , "TPCTrackerHits") ;
    
    
    for (int j = 0; j<the_tpc->getNumberOfRows(); j++){
      
      vector <Voxel_tpc *> current_row = the_tpc->getRowHits(j);
      
      for (int i = 0; i<current_row.size(); i++){
	
	delete current_row[i];
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
  
  //   std::cout << "TPCDigiProcessor::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
  
}

