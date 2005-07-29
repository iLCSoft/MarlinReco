#include "TrackingPerformanceProcessor.h"
#include <iostream>
#include <cfortran.h>
#include <map>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include <UTIL/LCRelationNavigator.h>

#include"tkmcbank.h"
#include"marlin_tpcgeom.h"
#include"constants.h"

#include "HelixClass.h"

using namespace constants;
using namespace std;

TrackingPerformanceProcessor aTrackingPerformanceProcessor ;

TrackingPerformanceProcessor::TrackingPerformanceProcessor() : Processor("TrackingPerformanceProcessor")
{
  // modify processor description
  _description = "Plots tracking performace histogrammes" ;

  registerProcessorParameter( "CollectionName" ,
			      "Name of the MC collection" ,
			      _colName ,
			      std::string("MCParticle") ) ;

}

void TrackingPerformanceProcessor::init()
{

  // usually a good idea to 
  printParameters();

#ifdef MARLIN_USE_ROOT
  // open ROOT-File
   _histoFile=new TFile("trackperf.root","RECREATE");
   // create histogram
  _efficiencyHisto=new TH1I("efficiency","efficiency all"
  		      ,3,0,3);
    // create histogram
  _efficiencyHisto1=new TH1I("efficiency1","efficiency 1-5GeV"
  		      ,3,0,3);
     // create histogram
  _efficiencyHisto2=new TH1I("efficiency2","efficiency 5-15GeV"
  		      ,3,0,3);
     // create histogram
  _efficiencyHisto3=new TH1I("efficiency3","efficiency 15-70GeV"
  		      ,3,0,3);
     // create histogram
  _efficiencyHisto4=new TH1I("efficiency4","efficiency > 70GeV"
  		      ,3,0,3);
   
      // create histogram
  _invpHisto1=new TH1F("invp1","invp 1-5GeV"
  		      ,100,-0.006,0.006);
     // create histogram
  _invpHisto2=new TH1F("invp2","invp 5-15GeV"
  		      ,100,-0.006,0.006);
     // create histogram
  _invpHisto3=new TH1F("invp3","invp 15-70GeV"
  		      ,100,-0.006,0.006);
     // create histogram
  _invpHisto4=new TH1F("invp4","invp > 70GeV"
  		      ,100,-0.006,0.006);
    // create histogram
  _ghostHisto=new TH1I("ghostHisto","number of ghosts in the TPC"
  		      ,20,0,21);
#endif
 
  _nRun = 0 ;
  _nEvt = 0 ;

}

void TrackingPerformanceProcessor::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
}

void TrackingPerformanceProcessor::processEvent( LCEvent * evt )
{
  static bool firstEvent = true ;

  // this gets called for every event
  // usually the working horse ...

  if(firstEvent==true) std::cout << "TrackingPerformanceProcessor called for first event" << std::endl;
 
  firstEvent = false ;
  
  LCCollection* MCcol = evt->getCollection( _colName ) ;
  LCCollection* STHcol = evt->getCollection( "tpc04_TPC" ) ;
  LCCollection* LCRcol = evt->getCollection(  "MCTrackRelations" ) ;
  LCRelationNavigator* nav = new LCRelationNavigator(LCRcol);
  LCCollection* Tcol = evt->getCollection( "TPCTracks" ) ;

  std::map <const MCParticle* , int > mc_hit_map;

  float number_of_ghosts = 0;

  std::cout << "Looking for TPC Tracks in Event" << std::endl;
  if( Tcol != 0) {

    int n_tracks = Tcol->getNumberOfElements();
    if( n_tracks==0 ) {
      std::cout << "No TPC Tracks in Event" << std::endl;
      //      getchar();
    }
    else {

      if(LCRcol !=0 )
	{
	  for(int i=0; i< n_tracks; i++){
	    Track* track = dynamic_cast<Track*>( Tcol->getElementAt( i ) ) ;

	    const LCObjectVec& LCRelVec = nav->getRelatedToObjects(track); 
	    const FloatVec& LCRelWeightsVec = nav->getRelatedToWeights(track);

	    float greatestweight=0;
	  
	    for(int j =0; j<LCRelWeightsVec.size(); j++){
	      
	      if( LCRelWeightsVec[j] > greatestweight) greatestweight=LCRelWeightsVec[j];
	    
	    }
	    
	    //	    std::cout << "greatest weight = " << greatestweight  << std::endl; 
	    if(greatestweight<0.75) number_of_ghosts++;

	  }
	  //	  number_of_ghosts = number_of_ghosts/n_tracks;
	  std::cout << "number of ghosts = " << number_of_ghosts  << std::endl; 
#ifdef MARLIN_USE_ROOT
	  _ghostHisto->Fill(number_of_ghosts);
#endif
	}
	  

      
      if( MCcol != 0 ){
	
	MCParticle* mcparticle;
	
	if( STHcol != 0 ){
	  
	  int n_sim_hits = STHcol->getNumberOfElements()  ;
	  
	  // create MAP which stores the number of Sim Trackerhits produced by MCParticle
	  
	  for(int i=0; i< n_sim_hits; i++){
	    
	    SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
	    mcparticle = SimTHit->getMCParticle();
	    // FIXME: SJA: check whether values for non existing keyword are really intialised to zero 
	    if (mcparticle!=NULL) mc_hit_map[mcparticle] = mc_hit_map[mcparticle] + 1;
	  }
	}
	
	int n_particles = MCcol->getNumberOfElements()  ;
	
	for(int i=0; i<n_particles; i++){
	  
	  mcparticle = dynamic_cast<MCParticle*>( MCcol->getElementAt( i ) ) ;
	  
	  
	//	std::cout << "MCParticle PID for this particle " << mcparticle->getPDG() << std::endl; 
	  
	  	  
	  //--------------------------------------------------
	  // tracking analysis
	  //	  std::cout << "the number of hits for this MC particle " << i << " is " << mc_hit_map[mcparticle] << std::endl;

	  //	  if(!mcparticle->isCreatedInSimulation()){	  
	  if (1) {		      

	    const double * mom;
	      mom = mcparticle->getMomentum();
	      
	      double mc_p = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+ mom[2]*mom[2]);	      
	      double mc_pt = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
	      double mc_theta = asin(mc_pt/mc_p);

	      const MCParticleVec& mcdaughtersvec = mcparticle->getDaughters();
	      
	      bool daughter_in_tpc = false;

	      for(int mcd=0; mcd < mcdaughtersvec.size(); mcd++){

		const double * mcdvetex = mcdaughtersvec.at(mcd)->getVertex();

		if(mcdvetex[2]*mcdvetex[2]<tpcgeom::zdrift*tpcgeom::zdrift){

		  double mcdrad2 = mcdvetex[0]*mcdvetex[0] + mcdvetex[1]*mcdvetex[1];

		  if(mcdrad2<(tpcgeom::tpcacro*tpcgeom::tpcacro)){

		    //		    if(mc_hit_map[mcdaughtersvec.at(mcd)] > 0) daughter_in_tpc = true;
		    if(sqrt(mcdrad2)>tpcgeom::tpcacri){
		      if(mc_hit_map[mcdaughtersvec.at(mcd)] > 0) daughter_in_tpc = true;
		      if(daughter_in_tpc){
			std::cout << "the PID of this daughter is " << mcdaughtersvec.at(mcd)->getPDG()<< std::endl;
			std::cout << "the radius of daugher = " << sqrt(mcdrad2) << std::endl;
			std::cout << "the z of daughter = " << sqrt(mcdvetex[2]*mcdvetex[2]) << std::endl;
			std::cout << "number of hits = " << mc_hit_map[mcdaughtersvec.at(mcd)] << std::endl;
			std::cout << "==============================" << std::endl;
		      }
		    }
		  }
		}
		

	      }
	      
	      //	      if (daughter_with_tracker_hits) std::cout << "this MCParticle has daughters with hits in the TPC" << std::endl; 

	      if(mc_hit_map[mcparticle]>2){
		if(mc_p>1.){

#ifdef MARLIN_USE_ROOT
 		  _efficiencyHisto->Fill(0) ;
		  if(mc_p>1.&&mc_p<5.) _efficiencyHisto1->Fill(0) ;
		  if(mc_p>5.&&mc_p<15.) _efficiencyHisto2->Fill(0) ;
		  if(mc_p>15.&&mc_p<70.) _efficiencyHisto3->Fill(0) ;
		  if(mc_p>70.) _efficiencyHisto4->Fill(0) ;
#endif		  
	    
		  if(LCRcol !=0 ) {

		    const LCObjectVec& LCRelVec = nav->getRelatedFromObjects(mcparticle);
		    const FloatVec& LCRelWeightsVec = nav->getRelatedFromWeights(mcparticle);

		    //		  std::cout << "the number of hits for this MCParticle is " << mc_hit_map[mcparticle] << std::endl;
		    //		  std::cout << "the numebr of tracks for this MCParticle is " << LCRelVec.size() << std::endl;

		    //		    if (!daughter_in_tpc) {
		    if (1) {		      

		      if(LCRelVec.size() == 1){ // LCRelVec.size() == 1
		      //		      if (1) {		      
		
			// test for d0

			float vertex[3];
			float momentum[3];
			
			for (int i=0; i<3; ++i) {
			  vertex[i]   = (float)mcparticle->getVertex()[i];
			  momentum[i] = (float)mcparticle->getMomentum()[i];
			}
			
			float charge = mcparticle->getCharge();
			// Magnetic field in Tesla
			float B = 4.0 ;
			
			HelixClass * helix = new HelixClass();
			helix->Initialize_VP(vertex, momentum, charge, B);
			
			float d0 = helix->getD0();
			float z0 = helix->getZ0();
			float phi = helix->getPhi0();
			if (phi<0.) phi = phi + 2*3.1415926535898 ;
			float omega = helix->getOmega();
			float tanLambda = helix->getTanLambda();
			
			cout << "d0 from MC = " << d0 ;
			cout << "  Z0 from MC = " << z0 ;
			cout << "  omega from MC = " << omega ;
			cout << "  tanLamda from MC = " << tanLambda ;
			cout << "  phi0 from MC = " << phi << endl;   
			//end of test for d0
			
			
			// 	    // transformation from 1/p to 1/R = CONSB * (1/p) / sin(theta)
			// 	    // CONSB is given by 1/R = (c*B)/(pt*10^9) where B is in T and pt in GeV  
			
			double CONSB = (2.99792458*4.)/(10*1000.);    // divide by 1000 m->mm
			Track* track = dynamic_cast<Track*>( LCRelVec.at(0) );

			cout << "d0 from LEP = " << track->getD0() ;
			cout << "  Z0 from LEP = " << track->getZ0() ;
			cout << "  omega from LEP = " << track->getOmega() ;
			cout << "  tanLamda from LEP = " << track->getTanLambda() ;
			cout << "  phi0 from LEP = " << track->getPhi() << endl;   
			cout << endl;
			
			double pt = fabs(1/(track->getOmega()/CONSB));
			
			double lambda = atan(track->getTanLambda());
			double invp = 1/(pt/sin(twopi/4.-lambda));
			
			//			double invp = 1/fabs(1/( (sin(twopi/4.-lambda) * track->getOmega())/CONSB));
			// 			std::cout << "#######################" << std::endl;
			// 			std::cout << "theta = "<< sin(twopi/4.-lambda) << std::endl;
			// 			std::cout << "mctheta = "<< mc_theta << std::endl;
			// 			std::cout << "LCRelWeightsVec[j] = " << LCRelWeightsVec.at(0) << std::endl; 
			// 			std::cout << "pt = " << pt << std::endl; 
			// 			std::cout << "1/invp = " <<1/invp << std::endl; 
			// 			std::cout << "mc_p = " << mc_p << std::endl; 
			// 			std::cout << "(1/mc_p) - invp = " <<(1/mc_p)-invp << std::endl; 
			// 			std::cout << "#######################" << std::endl;			


#ifdef MARLIN_USE_ROOT
 			_efficiencyHisto->Fill(1) ;
			if(mc_p>1.&&mc_p<5.) _efficiencyHisto1->Fill(1) ;
			if(mc_p>5.&&mc_p<15.) _efficiencyHisto2->Fill(1) ;
			if(mc_p>15.&&mc_p<70.) _efficiencyHisto3->Fill(1) ;
			if(mc_p>70.) _efficiencyHisto4->Fill(1) ;
			
			if(mc_p>1.&&mc_p<5.) _invpHisto1->Fill(invp-(1/mc_p));
			if(mc_p>5.&&mc_p<15.) _invpHisto2->Fill(invp-(1/mc_p));
			if(mc_p>15.&&mc_p<70.) _invpHisto3->Fill(invp-(1/mc_p));
			if(mc_p>70.) _invpHisto4->Fill((1/mc_p)-invp);			
#endif			
			
 		      }
		      
		      if(LCRelVec.size() > 1){
			if((mom[0]*mom[0]+mom[1]*mom[1])>4.){//pt > 4GeV
			  std::cout << "mom[0]*mom[0]+mom[1]*mom[1] = "  << mom[0]*mom[0]+mom[1]*mom[1]<< std::endl;
#ifdef MARLIN_USE_ROOT
			  _efficiencyHisto->Fill(2) ;
			  if(mc_p>1.&&mc_p<5.) _efficiencyHisto1->Fill(2) ;
			  if(mc_p>5.&&mc_p<15.) _efficiencyHisto2->Fill(2) ;
			  if(mc_p>15.&&mc_p<70.) _efficiencyHisto3->Fill(2) ;
			  if(mc_p>70.) _efficiencyHisto4->Fill(2) ;
#endif
			  std::cout << "this MCParticle has no daughters with hits in the TPC but has " << 
			  LCRelVec.size() << " tracks"<< std::endl; 
			  std::cout << "px = " << mom[0] << " py = " << mom[1]<< " pz = " << mom[2] << std::endl;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  
	  
	  //--------------------------------------------------
	  
	}
      }
    }
  }
  std::cout << "TrackingPerformanceProcessor: number of events processed = " << _nEvt << std::endl;
  delete nav;  
  _nEvt++;
}

void TrackingPerformanceProcessor::check( LCEvent * evt )
{
  // nothing to check here
}

void TrackingPerformanceProcessor::end()

{
#ifdef MARLIN_USE_ROOT
  // close root file
    _histoFile->Write();
    _histoFile->Close();
#endif
  std::cout << "TrackingPerformanceProcessor::end() " << std::endl;
}
