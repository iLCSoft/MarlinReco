#include "GeomEndProcessor.h"
#include <iostream>

#include "tpc.h"


GeomEndProcessor aGeomEndProcessor ;

GeomEndProcessor::GeomEndProcessor() : Processor("GeomEndProcessor") {
  
  // modify processor description
  _description = "Deletes geometry objects" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "" , 
			      ""  ,
			      _colName ,
			      std::string("") ) ;
}

void GeomEndProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void GeomEndProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void GeomEndProcessor::processEvent( LCEvent * evt ) { 

  static bool firstEvent = true ;
    
  // this gets called for every event 
  // usually the working horse ...

  if(firstEvent==true) std::cout << "GeomEndProcessor called for first event" << std::endl;

  firstEvent = false ;

  delete the_tpc;

  _nEvt ++ ;
  
}



void GeomEndProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void GeomEndProcessor::end(){ 
  
//   std::cout << "GeomEndProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

