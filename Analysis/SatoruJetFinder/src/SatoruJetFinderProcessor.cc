#include "SatoruJetFinderProcessor.h"
#include <iostream>
#include "EVENT/LCIO.h"
#include "EVENT/LCCollection.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/LCCollectionVec.h" 

using namespace std;

SatoruJetFinderProcessor aSatoruJetFinderProcessor ;


SatoruJetFinderProcessor::SatoruJetFinderProcessor() 
    : Processor("SatoruJetFinderProcessor") {
}
 
void SatoruJetFinderProcessor::init(){

  std::cout << "SatoruJetFinderProcessor::init()  " << name() 
	    << std::endl 
	    << "  parameters: " << std::endl 
	    << *parameters()  ;
  InputCollection=
      parameters()->getStringVal("InputCollection");
  OutputCollection=
      parameters()->getStringVal("OutputCollection");
  JetFindingMode= parameters()->getStringVal("Mode");
  Debug = parameters()->getIntVal("Debug");  
 
  if (JetFindingMode != "Manual"){
      if (JetFindingMode=="DurhamNjet"){
	  GlobalMode="0B";
	  PrimaryJetFindingMode=5;
	  NJetRequested= parameters()->getIntVal("NJetRequested");
	  cout << "GlobalMode "<<GlobalMode << " NJetRequested " << NJetRequested << endl;
      }else if(JetFindingMode=="DurhamYCut"){
	  GlobalMode="0A";
	  PrimaryJetFindingMode=5;
	  YCut[0]= parameters()->getFloatVal("YCut");
      }else if(JetFindingMode=="ConeBlanka"){
	  NJetRequested = 0;
	  Threshold = 0.7; 
	  PrimaryJetFindingMode = 12;
	  YCut[0]=  0.2;
	  YCut[1]= 0.7;
	  GlobalMode="0A";
      }else if(JetFindingMode=="Satoru"){
	  Threshold = 0.01;
	  NJetRequested= parameters()->getIntVal("NJetRequested");
	  PrimaryJetFindingMode = 5;
	  YCut[0] = 0.0;
	  MergingMode  = 2;
	  MergingThreshold= 0.;
	  SecondJetFindingMode =4;
	  GlobalMode = "2C";
      }
  }else{
      GlobalMode=
	  parameters()->getStringVal("GlobalMode");
      NJetRequested=
	  parameters()->getIntVal("NJetRequested");
      Threshold=
	  parameters()->getFloatVal("Threshold");
      PrimaryJetFindingMode=
	  parameters()->getIntVal("PrimaryJetFindingMode");
      YCut[0]=
	  parameters()->getFloatVal("YCut");
      MergingMode=
	  parameters()->getIntVal("MerginMode");
      MergingThreshold=
	  parameters()->getFloatVal("MergingThreshold");
      SecondJetFindingMode=
	  parameters()->getIntVal("SecondJetFindingMode");
  }
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void SatoruJetFinderProcessor::processRunHeader( LCRunHeader* run) { 
  std::cout << "SatoruJetFinderProcessor::processRun()  " << name() 
	    << " in run " << run->getRunNumber() 
	    << std::endl ;  
  _nRun++ ;
} 

void SatoruJetFinderProcessor::processEvent( LCEvent * evt ) { 
     LCCollectionVec* JetsCol= new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    // just a simple example    
    //first write the name of all collection included ...
    cout << " start" << endl;
    GoSatoru(evt,JetsCol);
    cout << " after" << endl;
  _nEvt ++ ;
}

void SatoruJetFinderProcessor::end(){
  std::cout << "SatoruJetFinderProcessor::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    << std::endl ;
};


/* *********************************************************************** */




void SatoruJetFinderProcessor::GoSatoru(LCEvent * evt,LCCollection* JetsCol){
    PutPartons(evt);
    //    WritePartons();
    CallSatoru(evt);
    GetJets(evt,JetsCol);
    //GetPointerFromPartonToJet();
};


void SatoruJetFinderProcessor::PutPartons(LCEvent * evt){
    LCCollection* enflowcol=evt->getCollection(InputCollection);
    int nenflow =  enflowcol->getNumberOfElements(); 
    PartonsWorkArray.NumberOfPartons=nenflow; 
    for ( int ienflow=0; ienflow<nenflow ; ienflow++){
	  ReconstructedParticle* enflow = 
	      dynamic_cast<ReconstructedParticle*>
	      (enflowcol->getElementAt( ienflow ));
	  for(int i=0;i<3;i++){
	      PartonsWorkArray.Momentum[ienflow*4+i]=(enflow->getMomentum())[i];	  }
	  PartonsWorkArray.Momentum[ienflow*4+3]=enflow->getEnergy();
    }
};



void SatoruJetFinderProcessor::WritePartons(){
    for (int iparton=0;iparton<PartonsWorkArray.NumberOfPartons;iparton++){
	cout << "Px,Py,Pz,E: "<< 
	    PartonsWorkArray.Momentum[iparton*4+0] << ", " <<
	    PartonsWorkArray.Momentum[iparton*4+1] << ", " <<
	    PartonsWorkArray.Momentum[iparton*4+2] << ", " <<
	    PartonsWorkArray.Momentum[iparton*4+3] << endl;    
    }
};


void SatoruJetFinderProcessor::CallSatoru(LCEvent * evt){
     int DimensionOfInputArray=4;
     int DimensionOfOutputArray=4;
     int MaximalNumberOfJets=20;
     float YMinus,YPlus;
     int IError;
     int GlobalModusLenght=GlobalMode.length()-1;

     syjkrn_(GlobalMode.c_str(),
	     NJetRequested,Threshold,
	     PrimaryJetFindingMode,YCut,
	     MergingMode,MergingThreshold,
	     SecondJetFindingMode,
	     PartonsWorkArray.NumberOfPartons,
	     DimensionOfInputArray,PartonsWorkArray.Momentum,
	     MaximalNumberOfJets,
	     JetsWorkArray.NumberOfJets,PartonsWorkArray.PointerParticleToJets,
	     DimensionOfOutputArray,JetsWorkArray.Momentum,
	     YMinus,YPlus,IError,GlobalModusLenght);
};

void SatoruJetFinderProcessor::GetJets(LCEvent * evt,LCCollection* JetsCol){
    LCCollection* enflowcol=evt->getCollection(InputCollection);
    cout << " " << endl;
    cout << " number of jets found: " <<JetsWorkArray.NumberOfJets << endl;
    for( int ijets=0;ijets<JetsWorkArray.NumberOfJets; ijets++){
	ReconstructedParticleImpl* Jets = new ReconstructedParticleImpl;
	float momentum[3],energy;
	momentum[0]=JetsWorkArray.Momentum[4*ijets+0];
	momentum[1]=JetsWorkArray.Momentum[4*ijets+1];
	momentum[2]=JetsWorkArray.Momentum[4*ijets+2];
	energy=JetsWorkArray.Momentum[4*ijets+3];
	Jets->setMomentum(momentum);
	Jets->setEnergy(energy);
	for( int iobj=0;iobj<PartonsWorkArray.NumberOfPartons;iobj++){
	    if(PartonsWorkArray.PointerParticleToJets[iobj]==(ijets+1)){
		ReconstructedParticle* enflow = 
		    dynamic_cast<ReconstructedParticle*> 
		    (enflowcol->getElementAt(iobj));
		Jets->addParticle(enflow);
	    }
	}
	JetsCol->addElement(Jets);
	
    }  
    evt->addCollection(JetsCol ,OutputCollection) ; 
    cout << "bla" << endl;
    cout << "Next Event" << endl;
};


void SatoruJetFinderProcessor::WriteJets(){};
void SatoruJetFinderProcessor::WriteParameters(){};
