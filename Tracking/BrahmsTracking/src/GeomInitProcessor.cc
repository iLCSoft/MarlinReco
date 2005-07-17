#include "GeomInitProcessor.h"
#include <iostream>

#include "tpc.h"
#include "marlin_tpcgeom.h"

//FIXME: SJA: get rid of using namespace and use e.g. tpcgeom::pixzres
using namespace tpcgeom;

GeomInitProcessor aGeomInitProcessor ;

GeomInitProcessor::GeomInitProcessor() : Processor("GeomInitProcessor") {
  
  // modify processor description
  _description = "Sets up geometry" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "" , 
			      ""  ,
			      _colName ,
			      std::string("") ) ;
}


void GeomInitProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void GeomInitProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void GeomInitProcessor::processEvent( LCEvent * evt ) { 

  static bool firstEvent = true ;
    
  // this gets called for every event 
  // usually the working horse ...

  if(firstEvent==true) std::cout << "GeomInitProcessor called for first event" << std::endl;

  firstEvent = false ;
  
  // set up the tpc
    
  // carefull here the last two values are the pixel sizes
  //  the_tpc = new Tpc(tpc_outerAR,tpc_innerAR,tpc_length,tpc_padrows,tpc_pitch,tpc_PixRPhi,tpc_PixZ);
  the_tpc = new Tpc(tpcacro,tpcacri,zdrift,nrtpc,tpcpadr,pix_rp,pix_z,tpc_rphi_res_max,tpc_z_res);
  
  
  _nEvt ++ ;
  
}



void GeomInitProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void GeomInitProcessor::end(){ 
  
//   std::cout << "GeomInitProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

