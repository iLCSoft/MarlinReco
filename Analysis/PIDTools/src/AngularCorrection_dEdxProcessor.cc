#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/AIDAProcessor.h>

#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/TrackImpl.h>

#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>
#include <DD4hep/DD4hepUnits.h>

#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"

// MarlinUtil                                                                                                                                                                                             
#include "LCGeometryTypes.h"
#include "SimpleHelix.h"
#include "HelixClass.h"

#include "signal.h"

#include "AngularCorrection_dEdxProcessor.hh"

AngularCorrection_dEdxProcessor aAngularCorrection_dEdxProcessor ;

AngularCorrection_dEdxProcessor::AngularCorrection_dEdxProcessor()
  : Processor("AngularCorrection_dEdxProcessor") {
  
  // Processor description
  _description = "Correct_AngularCorrection_dEdxProcessor: Makes a hard angular-based correction of dEdx for all the Tracks in the event. ATTENTION: this processor rewrites the MarlinTrk Collection and it is to be used only for simulations produced with ILCSoft from v02 to v02-02-01" ;
  
  registerInputCollection(LCIO::TRACK,
			  "LDCTrackCollection",
			  "LDC track collection name",
			  _LDCTrackCollection,
			  std::string("MarlinTrkTracks"));

  registerProcessorParameter( "Write_dEdx",
			      "If set, the calculated dEdx value will be rewritten in the its corresponding track (default: false).",
			      _writedEdx,
			      bool(false));

  std::vector<float> _newpar = {  0.970205,
				  0.0007506,  
				  4.41781e-8,
				  5.8222e-8};

  // ***************************************
  // Fit of NormLamdaFullAll_1 (dEdxAnalyser) using single particle samples recosntructed with v02-02-01
  // 2021/04
  // (including the default previous angular correction)
  //Minimizer is Linear
  //Chi2                      =      807.398
  //NDf                       =           27
  //p0                        =     0.970205   +/-   0.000127468 
  //p1                        =   0.00075065   +/-   1.55853e-05 
  //p2                        =  4.41781e-08   +/-   5.09662e-07 
  //p3                        =   5.8222e-08   +/-   4.71913e-09 
  //
  registerProcessorParameter( "AngularCorrectionParameters",
			      "parameter for new angular correction dedx= uncorrected_dedx  / f, with f= pol3(lambda)",
			      _par,
			      _newpar);



} 

void AngularCorrection_dEdxProcessor::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  // it's usually a good idea to
  printParameters();
  
}

void AngularCorrection_dEdxProcessor::processRunHeader( LCRunHeader* ) { 
} 

void AngularCorrection_dEdxProcessor::processEvent( LCEvent * evt ) { 

  //fill values
  if (_writedEdx)
    {
      _LDCCol = evt->getCollection( _LDCTrackCollection ) ;
      int nTrkCand = _LDCCol->getNumberOfElements();
      
      for (int iTRK=0;iTRK<nTrkCand;++iTRK) {
	
	TrackImpl * trkCand = (TrackImpl*) _LDCCol->getElementAt( iTRK );
	
	float dedx=trkCand->getdEdx();
	float dedx_error=trkCand->getdEdxError();
	float trklambda = trkCand->getTanLambda();

	float lambda = fabs(atan(trklambda)*180./M_PI);
	double  f3 = 1 / (_par[0] + _par[1] * lambda  + _par[2] * pow(lambda,2) + _par[3] * pow(lambda,3) );
	
	double new_dedx = dedx*f3;
	
	streamlog_out(DEBUG) << "Original dEdx: " <<dedx <<" Error: "<<dedx_error <<std::endl;
	streamlog_out(DEBUG) << "NeW dEdx: " <<new_dedx <<" Error: "<<dedx_error <<std::endl;
	
	//fill values
	trkCand->setdEdx(new_dedx);
	trkCand->setdEdxError(dedx_error);
      }
    }   else 
    {
      streamlog_out(ERROR) << " Why do you use this processor and not re-write dEdx ?? Check your steering file." <<std::endl;
      raise(SIGSEGV);
    }
}

void AngularCorrection_dEdxProcessor::check( LCEvent * ) { 
}

void AngularCorrection_dEdxProcessor::end() { 
}


