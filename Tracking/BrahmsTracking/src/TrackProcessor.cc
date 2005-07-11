#include "TrackProcessor.h"
#include <iostream>
#include <vector>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogram1D.h>
#endif

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCRelation.h>
#include <UTIL/LCRelationNavigator.h>


#include"constants.h"

using namespace std;
using namespace constants;

TrackProcessor aTrackProcessor ;

TrackProcessor::TrackProcessor() : Processor("TrackProcessor") {
  
  // modify processor description
  _description = "Plots the delta inverse momentum of tracks reconstructed in the TPC" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "CollectionName" , 
			      "Name of the Track collection"  ,
			      _colName ,
			      std::string("TPC_Tracks") ) ;

#ifdef MARLIN_USE_ROOT
   registerProcessorParameter( "ROOTFile" ,
 			      " name of output file"  ,
 			      _rootfilename,
 			      std::string("masshist.root") ) ;
  registerProcessorParameter( "Histoname" ,
			      " name of the histogram" ,
			      _roothistname,
			      std::string("InvMass"));
  registerProcessorParameter( "nBins" ,
                                  " number of bins" ,
			      _rootnbins, 50 );
  registerProcessorParameter( "HistLow" ,
			      " low edge of first bin" ,
			      _roothistlow, (float) 0. );
  registerProcessorParameter( "HistHigh" ,
			      " upper edge of last bin" ,
			      _roothisthigh, (float) 200. );
  
#endif

}

void TrackProcessor::init() { 

  // usually a good idea to
  printParameters() ;



  
#ifdef MARLIN_USE_ROOT
  // open ROOT-File
  //  std::cout << "ROOT Filename" << _rootfilename.c_str() << std::endl;
   _histoFile=new TFile(_rootfilename.c_str(),"RECREATE");

    // create histogram
  _massHisto=new TH1F(_roothistname.c_str(),_roothistname.c_str()
  		      ,_rootnbins,_roothistlow,_roothisthigh);
  
#endif
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void TrackProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void TrackProcessor::processEvent( LCEvent * evt ) { 

//   static bool firstEvent = true ;

//   // this gets called for every event 
//   // usually the working horse ...

//   if(firstEvent==true) std::cout << "TrackProcessor called for first event" << std::endl;

  if( isFirstEvent() ) std::cout << "TrackProcessor called for first event" << std::endl;

#ifdef MARLIN_USE_AIDA
      // define a histogram pointer
      static AIDA::ICloud1D* invdpcloud ;

      if(  isFirstEvent() ) { 
    
	invdpcloud = 
	  AIDAProcessor::histogramFactory(this)->
	  createCloud1D( "invdp", "1/dp of tracks", 40 ) ; 
      }

#endif

  LCCollection* Tcol = evt->getCollection( _colName ) ;
  LCCollection* LCRcol = evt->getCollection(  "MC_Track_Relations" ) ;
  
  if( Tcol != 0 ){

    int n_tracks = Tcol->getNumberOfElements()  ;   
    
    std::cout << "the number of tracks in this event = " << n_tracks << std::endl; 

    for(int i=0; i< n_tracks; i++){
      
      Track* track = dynamic_cast<Track*>( Tcol->getElementAt( i ) ) ;

      //      std::cout << "track curvature = " << track->getOmega() << std::endl;

// 	    // transformation from 1/p to 1/R = CONSB * (1/p) / sin(theta)
// 	    // CONSB is given by 1/R = (c*B)/(pt*10^9) where B is in T and pt in GeV  
	    
      float CONSB = (2.99792458*4.)/(10*1000.);    // divide by 1000 m->mm

      float momentum1 = fabs(1/( (sin(twopi/4.-track->getTanLambda()) * track->getOmega())/CONSB));


      if(LCRcol !=0 )
	{
	  LCRelationNavigator* nav = new LCRelationNavigator(LCRcol);
	  
	  const LCObjectVec& LCRelVec = nav->getRelatedToObjects(track); 
	  const FloatVec& LCRelWeightsVec = nav->getRelatedToWeights(track);
	  
	  
	  
	  int greatestweight=0;
	  float oldweight=0;
	  
	  for(int j =0; j<LCRelWeightsVec.size(); j++){
	    
	    if( LCRelWeightsVec[j] > oldweight) greatestweight=j;
	    
	  }
	  
	  MCParticle* mcp = dynamic_cast<MCParticle*>(LCRelVec.at(greatestweight));	  

// 	  const float * mom;
// LCIO v01-05 needs double !
	  const double * mom;

	  if(mcp != NULL) {
	    
	    mom = mcp->getMomentum();
	    
	    float p = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+ mom[2]*mom[2]);
	    
	    //	  std::cout << "mc particle momemtum = " << p << std::endl;
	    
#ifdef MARLIN_USE_ROOT
	    _massHisto->Fill( (1./momentum1) - (1./p) ) ;
#endif
	    
#ifdef MARLIN_USE_AIDA
	    invdpcloud->fill( (1./momentum1) - (1./p) ) ;
#endif
	  }
	  
	}
      
      



#ifdef MARLIN_USE_AIDA

//     // define a histogram pointer
//     static AIDA::ICloud1D* hZmass ;    
    
//     if( firstEvent ) { 
      
//       hZmass = 
// 	AIDAProcessor::histogramFactory(this)->
// 	createCloud1D( "hZmass", "invariant mass of tracks", 40 ) ; 
//     }
    
//     // fill histogram from LCIO data :
    
//     for(int i=0; i<(n_tracks-1); i++)
//       {
// 	for(int j=i+1;j<n_tracks; j++)
// 	  {

// 	    Track* track1 = dynamic_cast<Track*>( Tcol->getElementAt( i ) ) ;
// 	    Track* track2 = dynamic_cast<Track*>( Tcol->getElementAt( j ) ) ;	    

	    
// 	    // transformation from 1/p to 1/R = CONSB * (1/p) / sin(theta)
// 	    // CONSB is given by 1/R = (c*B)/(pt*10^9) where B is in T and pt in GeV  
	    
// 	    float CONSB = (2.99792458*4.)/(10*1000.);    // divide by 1000 m->mm

// 	    float momentum1 = fabs(1/( (sin(twopi/4.-track1->getTanLambda()) * track1->getOmega())/CONSB));
	    
// 	    float phi1 = track1->getPhi();
// 	    float theta1 = twopi/4.-track1->getTanLambda(); 
	    
// 	    float pz1 = cos(theta1)*momentum1;
//             float py1 = sin(phi1)*sin(theta1)*momentum1;
//             float px1 = cos(phi1)*sin(theta1)*momentum1;

// 	    float momentum2 = fabs(1/( (sin(twopi/4.-track2->getTanLambda()) * track2->getOmega())/CONSB));
	    
// 	    float phi2 = track2->getPhi();
// 	    float theta2 = twopi/4.-track2->getTanLambda(); 
	    
// 	    float pz2 = cos(theta2)*momentum2;
//             float py2 = sin(phi2)*sin(theta2)*momentum2;
//             float px2 = cos(phi2)*sin(theta2)*momentum2;

// 	    float Mass_mu = 0.1057;
// 	    float Mass_Z;
// 	    //    Double_t P_mu = (p[0] + p[1])/2.;
	
// 	    float E_1 = sqrt(momentum1*momentum1 + Mass_mu*Mass_mu);
// 	    float E_2 = sqrt(momentum2*momentum2 + Mass_mu*Mass_mu);
	
// 	    Mass_Z = sqrt((Mass_mu*Mass_mu+Mass_mu*Mass_mu)+(2*E_1*E_2)-2.*(px1*px2+py1*py2+pz1*pz2));

// 	    if(Mass_Z>80. && Mass_Z<100.) hZmass->fill( Mass_Z ) ;

// 	  }
	

	
#endif
    }
  }
  
  //   firstEvent = false 
  _nEvt ++ ;

}



void TrackProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TrackProcessor::end(){ 

#ifdef MARLIN_USE_ROOT
  // close root file
    _histoFile->Write();
    _histoFile->Close();
#endif

//   std::cout << "TrackProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

