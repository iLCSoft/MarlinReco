#include "BrahmsEndProcessor.h"
#include <iostream>

#include"tkmcbank.h"
#include"tpchitbank.h"
#include"tkhitbank.h"
#include "tktebank.h"



BrahmsEndProcessor aBrahmsEndProcessor ;

BrahmsEndProcessor::BrahmsEndProcessor() : Processor("BrahmsEndProcessor") {
  
  // modify processor description
  _description = "Deletes tkbanks " ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "" , 
			      ""  ,
			      _colName ,
			      std::string("") ) ;
}

void BrahmsEndProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void BrahmsEndProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void BrahmsEndProcessor::processEvent( LCEvent * evt ) { 

  static bool firstEvent = true ;
    
  // this gets called for every event 
  // usually the working horse ...

  if(firstEvent==true) std::cout << "BrahmsEndProcessor called for first event" << std::endl;

  firstEvent = false ;

  delete TkMCBank;
  delete TPCHitBank;
  delete TkHitBank;
  delete TkTeBank;

  _nEvt ++ ;
  
}



void BrahmsEndProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void BrahmsEndProcessor::end(){ 
  
//   std::cout << "BrahmsEndProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

