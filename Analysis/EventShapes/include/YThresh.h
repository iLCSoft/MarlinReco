//-------------------------------------------------------------------------------------------------------------
//   The YThresh module calculates the yThresh variable, which is the crossover value of the yCut jet finding  
//   variable from NMin to NMin+1 jets found using durhamycut.  For example, if NMin=2 the yThresh variable is 
//   the value of yCut above which durhamycut returns 2 jets, below which it returns 3 jets. The value of      
//   yThresh will be stored as a parameter in the ReconstructedParticle collection with name y[NMin][NMin+1],  
//   ie y23 for NMin=2. SatoruJetFinder must be installed for this package to run.
//   For questions/comments please email Ben Hooberman at benhooberman@berkeley.edu
//-------------------------------------------------------------------------------------------------------------


#ifndef YThresh_h
#define YThresh_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
//#include "TH1.h"
#include "IMPL/LCCollectionVec.h"

using namespace lcio ;

struct SatoruPartonArray{
  int NumberOfPartons;
  float Momentum[1200];
  int PointerParticleToJets[450];
};
struct SatoruJetsArray{
  int NumberOfJets;
  float Momentum[80];
};

extern "C" {
  void syjkrn_( const char* GlobalModus_,
		// Primary jet finding                                                                       
		int &NJetRequested_,float &Threshold_,
		int &PrimaryJetFindingMode_,float* YCut_,
		//  Reassosiation after first jet finding                                                    
		int &MergingMode_,float &MergingThreshold_,
		// Second Jet finding mode                                                                   
		int &SecondJetFindingMode_,
		// input array--> array of                                                                   
		// PPartons(DimensionOfInputArray(4-6),NumberOfPartons)                                      
		int &NumberOfPartons_,
		int &DimensionOfInputArray_,float *PPartons_,
		// Output                                                                                    
		int &MaximalNumberOfJets_,
		int &NJetsFound_,int *PointerParicleToJet_,
		int &DimensionOfOutputArray_,float *PJets_,
		float &YMinus_,float &YPlus_,
		// error flag                                                                                
		int &IError_,int GlobalModusLenght_)
       ;}

using namespace marlin ;
using namespace std;

class YThresh : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new YThresh ; }
   
  YThresh() ;

  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ;   
  virtual void check( LCEvent * evt ) ;   
  void putPartons(LCEvent * evt);
  void callSatoru(LCEvent * evt);
  virtual void end() ;
  
 protected:

  int _nRun,_nEvt;
  int _nRecoParticles,_printOutput;
  int _nMin,_nIter;
  std::string _inputCollection;
  std::string _jetFindingMode;
  SatoruPartonArray _partonsWorkArray;
  SatoruJetsArray _jetsWorkArray;                                                                        
  std::string _globalMode;
  int _nMinParticles;
  int _output;
  int _nJetRequested;
  float _threshold;
  int _primaryJetFindingMode;
  float _yCutParam;
  float _rConeParam;
  float _epsConeParam;
  float _yCut[2];
  int _mergingMode;
  float _mergingThreshold;
  int _secondJetFindingMode;
  float _yStart;
} ;

#endif



