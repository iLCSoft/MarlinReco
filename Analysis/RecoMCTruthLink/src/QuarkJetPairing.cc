#include "QuarkJetPairing.h"

#include <vector>
#include <iostream>
#include <string>

#include <EVENT/LCCollection.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/Operators.h>

#define _USE_MATH_DEFINES

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include "TNtupleD.h"
#include <TProfile.h>
#include "TMVA/Reader.h"
#include <TMath.h>
#include <TLorentzVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <TVector3.h>
#include <cmath>
#include <math.h>

#include "UTIL/LCRelationNavigator.h"
#include "EVENT/LCIO.h"
#include "marlin/Exceptions.h"

#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/ReconstructedParticle.h>
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/LCCollectionVec.h" 
//#include <root/TLorentzVector.h>
#include "IMPL/LCCollectionVec.h"


using namespace lcio ;
using namespace marlin ;
using namespace CLHEP ;
using namespace std;
using namespace TMVA;

QuarkJetPairing aQuarkJetPairing ;

QuarkJetPairing::QuarkJetPairing() : Processor("QuarkJetPairing") {
  // modify processor description
  _description = "" ;
  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "InputJetCollection" ,
                           "Input collection of jets",
                           _inputJetCollection,
                           std::string("NewJets"));


  registerInputCollection( LCIO::MCPARTICLE,
                           "InputMCPCollection" ,
                           "Input collection of MCParticles",
                           _inputMCPCollection,
                           std::string("MCParticlesSkimmed")); 
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void QuarkJetPairing::init() {
  streamlog_out(DEBUG) << " " << std::endl ;
  printParameters() ;

  _nRun = 0;
  _nEvt = 0;
  
  /////////////////////////creates outputfile//////////////////////////
  //  TFile *outRootFile = new  TFile ("output.root", "RECREATE");
  /////////////////////////////////////////////////////////////////////



}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void QuarkJetPairing::processRunHeader( LCRunHeader*  /*run*/) {
  
  _nRun++;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void QuarkJetPairing::processEvent( LCEvent * evt ) { 

 

  _nEvt++;
  ///////////////////////////////////frequency distribution////////////////////////////////////
  static int counter_events=0;
  static int counter_alpha = 0;
  static int counter_energy = 0;
  static int counter_theta = 0;
  
  counter_events++;
  //////////////////////////////////////////////////set variables to zero/////////////////////////////////////////////////

  float energy = 0;
 
  float jet_ene[4];
  float quark_ene[4];
     
  float jet_theta[4];
  float quark_theta[4]; 
     
  float  jet_px[4];
  float quark_px[4];
     
  float  jet_py[4];
  float quark_py[4]; 
     
  float jet_pz[4];
  float quark_pz[4];
     
  float jet_ptot[4];
  float quark_ptot[4];
     
  int quark[4][24]; 
  int jet[4][24]; 
     
  for(int i=0; i<4; i++){
       
    jet_ene[i]=0;
    quark_ene[i]=0;
       
    jet_theta[i]=0;
    quark_theta[i]=0;
       
    jet_px[i]=0;
    quark_px[i]=0;
    jet_py[i]=0;
    quark_py[i]=0;
    jet_pz[i]=0;
    quark_pz[i]=0;
    jet_ptot[i]=0;
    quark_ptot[i]=0;
  } 

     
  ///////////////////////////////JET COLLECTION//////////////////////////////////
     
  LCCollection* jetcol = evt->getCollection( _inputJetCollection ) ; 

   
  if (jetcol != 0) {
       
    int nJETS = jetcol->getNumberOfElements()  ;
       
    if (nJETS != 4) return; 

    for(int i=0; i< nJETS ; i++){
         
 
      ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>( jetcol->getElementAt( i ) ) ;
               
      if (j) {
	    
	jet_px[i]=j->getMomentum()[0]; 
	jet_py[i]=j->getMomentum()[1];
	jet_pz[i]=j->getMomentum()[2];
	jet_ene[i]=j->getEnergy();
 	    
      }  
      jet_ptot[i]=sqrt(pow(jet_px[i], 2)+pow(jet_py[i], 2)+pow(jet_pz[i], 2));
    }
  }  

  ////////////////////////////MCParticle Collection/////////////////////////////
 
   LCCollection* mcpcol = evt->getCollection( _inputMCPCollection ) ;
  if (mcpcol != 0) {
         
    for(int e=6; e<10 ; e++){
         
      MCParticle* q = dynamic_cast<MCParticle*>( mcpcol->getElementAt( e ) ) ;
               
      if (q) {
	    
  	quark_px[e-6]=q->getMomentum()[0];
  	quark_py[e-6]=q->getMomentum()[1];
  	quark_pz[e-6]=q->getMomentum()[2];
  	quark_ene[e-6]=q->getEnergy();
      }
      quark_ptot[e-6]=sqrt(pow(quark_px[e-6], 2)+pow(quark_py[e-6], 2)+pow(quark_pz[e-6], 2));
    } 
  }
      
  /////////////////////////////////////////////set variables////////////////////////////////////////////////
 
  ipair = 0; 
  iperm = 0; 
 
  int counter_jets = 0; 
   
  Float_t theta[4];
  
  Float_t alpha[24];
  Float_t alpha_min = 9999; 
  int iperm_min = 0; 
   
  //////////////////////////////////////////permutation algorithm///////////////////////////////////////////
 
  for(int i=0; i<4; i++){
    energy+=quark_ene[i];
  }
 
    if(energy >495){
      counter_energy++;
  
      for (int index1=0; index1<4; index1++) {
	for (int index2=0; index2<4; index2++) {
	  if (index2 != index1) {
	    for (int index3=0; index3<4; index3++) {
	      if ((index3 != index2) && (index3 != index1)) {
		for (int index4=0; index4<4; index4++) {
		  if ((index4 != index3) && (index4 != index2) && (index4 != index1)) {
		    quark[0][iperm]=index1;
		    quark[1][iperm]=index2;
		    quark[2][iperm]=index3;
		    quark[3][iperm]=index4;
		    jet[0][iperm]=0;
		    jet[1][iperm]=1;
		    jet[2][iperm]=2;
		    jet[3][iperm]=3;
		    iperm++;
		  }
		}  
	      }
	    }
	  }	     
	}
      }
    
      for(int n = 0; n < iperm; n++){
	alpha[n] = 0;
	for (int ijet=0; ijet<4; ijet++) { 
	  alpha[n]+= acos(((jet_px[jet[ijet][n]]*quark_px[quark[ijet][n]])+(jet_py[jet[ijet][n]]*quark_py[quark[ijet][n]])+(jet_pz[jet[ijet][n]]*quark_pz[quark[ijet][n]]))/(jet_ptot[jet[ijet][n]]*quark_ptot[quark[ijet][n]]));
	}
	if(alpha[n] < alpha_min){
	  alpha_min = alpha[n];
	  iperm_min=n;
	}
	alpha[n]=0;
      }
    
      // message<DEBUG>( log() << "alpha_min: " << alpha_min ) ;
       
      if(alpha_min < 1){
	counter_alpha++;
      
	for(int ijet = 0; ijet<4; ijet++){
	
	  //////////////////////////////////////////// Determination of theta //////////////////////////
	
	  theta[ijet]=0;
	  quark_theta[ijet] = atan((sqrt((pow(quark_px[quark[ijet][iperm_min]],2))+(pow(quark_py[quark[ijet][iperm_min]],2))))/quark_pz[quark[ijet][iperm_min]]);
	  jet_theta[ijet] = atan((sqrt((pow(jet_px[jet[ijet][iperm_min]],2))+(pow(jet_py[jet[ijet][iperm_min]],2))))/jet_pz[jet[ijet][iperm_min]]);
	  if(quark_theta[ijet] < 0){
	    quark_theta[ijet] = M_PI+quark_theta[ijet];	
	  }
	  if(jet_theta[ijet] < 0){
	    jet_theta[ijet] = M_PI+jet_theta[ijet];	
	  }
	
	  if(((quark_theta[ijet] > 0.14) && (quark_theta[ijet]<3)) &&((jet_theta[ijet] > 0.14) && (jet_theta[ijet]<3)) ){
	    counter_jets++;
	 
	  } 
	  else{  
	    throw marlin::SkipEventException(this); 
	    counter_jets = 0;
	  }
	}
      }
      else{  
	throw marlin::SkipEventException(this);
      }
    }
    else{
      throw marlin::SkipEventException(this);
    }

    if(counter_jets == 4){

    counter_theta++;
  } 

  message<DEBUG>( log() << "Counter_events: " << counter_events ) ; 
  message<DEBUG>( log() << "Counter_energy: " << counter_energy ) ; 
  message<DEBUG>( log() << "Counter_alpha: " << counter_alpha ) ; 
  message<DEBUG>( log() << "Counter_theta: " << counter_theta ) ; 
      
  counter_jets = 0;
  message<DEBUG>( log() << "Event: " << evt->getEventNumber() ) ;
   
 
} 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void QuarkJetPairing::check( LCEvent *  /*evt*/ ) {


}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void QuarkJetPairing::end(){

  std::cout << "Finish run" <<std::endl;
}
