//-------------------------------------------------------------------------------------------------------------
//   The YThresh module calculates the yThresh variable, which is the crossover value of the yCut jet finding  
//   variable from NMin to NMin+1 jets found using durhamycut.  For example, if NMin=2 the yThresh variable is 
//   the value of yCut above which durhamycut returns 2 jets, below which it returns 3 jets. The value of      
//   yThresh will be stored as a parameter in the ReconstructedParticle collection with name y[NMin][NMin+1],  
//   ie y23 for NMin=2. SatoruJetFinder must be installed for this package to run.
//   For questions/comments please email Ben Hooberman at benhooberman@berkeley.edu
//-------------------------------------------------------------------------------------------------------------

#include "YThresh.h"
#include <iostream>
#include <math.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Cluster.h>
#include <EVENT/Track.h>
#include <UTIL/LCTOOLS.h>
#include <EVENT/ReconstructedParticle.h>

using namespace lcio ;
using namespace marlin ;
using namespace std;

YThresh aYThresh ;

YThresh::YThresh() : Processor("YThresh") {
  
  // modify processor description
  _description = "YThresh finds the crossover value of the yCut variable from NMin to NMin+1 jets found using durhamycut";
  
  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                 "RecoParticleCollection",
			     "Name of the input ReconstructedParticle collection",
			     _inputCollection,
			     std::string("RecoParticles"));

  registerProcessorParameter("NMin",
			     "min number of jets, ie. NMin=2 will give y23 variable",
			     _nMin,
			     int(2));

  registerProcessorParameter("PrintOutput",
			     "toggle print text output",
			     _output,
			     int(0));

  registerProcessorParameter("NIterations",
			     "number of iterations",
			     _nIter,
			     int(20));

  registerProcessorParameter("NMinParticles",
			     "min num particles for ythresh calculation",
			     _nMinParticles,
			     int(3));

  registerProcessorParameter("YStart",
			     "starting value for yCut",
			     _yStart,
			     float(0.1));

}

void YThresh::init() { 

  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  _globalMode="0A";
  _primaryJetFindingMode=5;

}

void YThresh::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void YThresh::processEvent( LCEvent * evt ) { 

  try{
    
    LCCollection * recoColl=evt->getCollection(_inputCollection.c_str());
    _nRecoParticles=recoColl->getNumberOfElements();

    float ymax=0;

    if(_nRecoParticles>=_nMinParticles){
      
      putPartons(evt);

      //Phase 1: find yCut value for which NJets<=NMin
      //start with yCut=YStart, increase yCut until this condition is satisfied
      _yCut[0]=_yStart;
      
      do{	
	callSatoru(evt);
      	if(_jetsWorkArray.NumberOfJets<=_nMin) break;
	else _yCut[0]=_yCut[0]*2;
      }while(true);

      //Phase 2: Start with yCut value determined in Phase 1 (ymax)
      //Find yThresh value between 0,ymax satisfying NJets=NMin for yCut>yThresh, NJets=NMin+1 for yCut<yThresh
      ymax=(float)_yCut[0];
      _yCut[0]=_yCut[0]/2;
      
      for(int i=1;i<=_nIter;i++){
	callSatoru(evt);	
	if(_jetsWorkArray.NumberOfJets>_nMin)  _yCut[0]+=ymax/(pow(2.,i+1));
	else                                   _yCut[0]-=ymax/(pow(2.,i+1));
      }
    }
    else{
      cout<<"YThresh Error "<<_nRecoParticles<<" particles<"<<_nMinParticles<<endl;
      _yCut[0]=-1.;
    }

    //Store yThresh value with name y[NMin][NMin+1] (ie y23 for NMin=2) as parameter in RecoParticleCollection
    stringstream tempStream;      
    tempStream<<'y'<<_nMin<<_nMin+1;      
    string yname=tempStream.str();      
    recoColl->parameters().setValue(yname.c_str(),_yCut[0]);
    
    if(_output>0)cout<<"YThresh Event "<<_nEvt<<" "<<yname<<"="<<_yCut[0]<<endl;

    if(_yCut[0]>=0){

      //Check to make sure that NJets=NMin for yCut=yThresh+epsilon, NJets=NMin+1 for yCut=yThresh-epsilon
      float epsilon=ymax/(pow(2.,_nIter));
      
      _yCut[0]+=epsilon;
      callSatoru(evt);
      int n1=_jetsWorkArray.NumberOfJets;
      
      _yCut[0]-=2*epsilon;
      callSatoru(evt);
      int n2=_jetsWorkArray.NumberOfJets;
      
      if(!(n1==_nMin && n2==_nMin+1)) cout<<"ERROR INCORRECT YTHRESH VALUE OBTAINED"<<endl;

    }

  }catch(DataNotAvailableException &e){
    cout<<"YThresh: cannot find collection "<<_inputCollection.c_str()<<endl;
  }
    
    _nEvt ++ ;
}

void YThresh::putPartons(LCEvent * evt){

  LCCollection* enflowcol=evt->getCollection(_inputCollection);
  int nenflow =  enflowcol->getNumberOfElements();
  _partonsWorkArray.NumberOfPartons=nenflow;
  for ( int ienflow=0; ienflow<nenflow ; ienflow++){
      ReconstructedParticle* enflow =
        dynamic_cast<ReconstructedParticle*>
        (enflowcol->getElementAt( ienflow ));
      for(int i=0;i<3;i++){
        _partonsWorkArray.Momentum[ienflow*4+i]=(enflow->getMomentum())[i]; }
      _partonsWorkArray.Momentum[ienflow*4+3]=enflow->getEnergy();
  }
};

void YThresh::callSatoru(LCEvent * evt){

  int DimensionOfInputArray=4;
  int DimensionOfOutputArray=4;
  int MaximalNumberOfJets=20;
  float YMinus,YPlus;
  int IError;
  int GlobalModusLength=_globalMode.length();

  syjkrn_(_globalMode.c_str(),
	  _nJetRequested,
	  _threshold,
	  _primaryJetFindingMode,
	  _yCut,
	  _mergingMode,
	  _mergingThreshold,
	  _secondJetFindingMode,
	  _partonsWorkArray.NumberOfPartons,
	  DimensionOfInputArray,
	  _partonsWorkArray.Momentum,
	  MaximalNumberOfJets,
	  _jetsWorkArray.NumberOfJets,
	  _partonsWorkArray.PointerParticleToJets,
	  DimensionOfOutputArray,
	  _jetsWorkArray.Momentum,
	  YMinus,
	  YPlus,
	  IError,
	  GlobalModusLength);


};

void YThresh::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void YThresh::end(){   
  
  std::cout << "YThresh::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    << std::endl ;
  
}
