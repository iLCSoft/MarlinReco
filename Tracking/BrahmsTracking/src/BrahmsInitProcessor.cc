#include "BrahmsInitProcessor.h"
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


#include "tkmcbank.h"
#include "tpchitbank.h"
#include "tkhitbank.h"
#include "tktebank.h"

BrahmsInitProcessor aBrahmsInitProcessor ;

BrahmsInitProcessor::BrahmsInitProcessor() : Processor("BrahmsInitProcessor") {
  
  // modify processor description
  _description = "Sets up the tkbanks for the fortran tracking" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "" , 
			      ""  ,
			      _colName ,
			      std::string("") ) ;
}


void BrahmsInitProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void BrahmsInitProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void BrahmsInitProcessor::processEvent( LCEvent * evt ) { 

  static bool firstEvent = true ;
    
  // this gets called for every event 
  // usually the working horse ...

  if(firstEvent==true) std::cout << "BrahmsInitProcessor called for first event" << std::endl;

  firstEvent = false ;
  
  TkMCBank = new Tk_MC_Bank;
  TPCHitBank = new TPC_Hit_Bank;  
  TkHitBank = new Tk_Hit_Bank;  
  TkTeBank = new Tk_Te_Bank;  
  
  _nEvt ++ ;
  
}



void BrahmsInitProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void BrahmsInitProcessor::end(){ 
  
//   std::cout << "BrahmsInitProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

