#ifndef LikelihoodPIDProcessor_hh
#define LikelihoodPIDProcessor_hh 1

#include <string>
#include <vector>
#include <marlin/Processor.h>
#include <iostream>
#include <cstring>
#include <fstream>
#include <signal.h>

#include <EVENT/LCCollection.h>

using namespace lcio ;
using namespace marlin ;

class LikelihoodPID;
class LowMomentumMuPiSeparationPID_BDTG;

class LikelihoodPIDProcessor : public Processor{
public:
  virtual Processor*  newProcessor() { return new LikelihoodPIDProcessor ; }
  LikelihoodPIDProcessor();
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run);
  virtual void processEvent( LCEvent * evt );
  virtual void check( LCEvent * evt );
  virtual void end();
 
private:
  void createParticleIDClass(int parttype, ReconstructedParticle *part, PIDHandler &pidh, int algoID, float MVAoutput);
  
  LikelihoodPID *_myPID{};
  std::string _description{};
  std::string _inputPFOsCollection{};
  std::string _PDFName{};
  std::vector<std::string> _weightFileName{};
 
  std::vector<std::string> _methodstorun{};
  std::string _methodstorun_version{};

  EVENT::FloatVec _energyBoundary{};
  LCCollection* _pfoCol{};
  std::vector<int> _pdgTable{};
  std::vector<std::string> _particleNames{};
  std::vector<std::string> _dEdxNames{};

  std::vector<float> _dEdxParamsElectron{};
  std::vector<float> _dEdxParamsMuon{};
  std::vector<float> _dEdxParamsPion{};
  std::vector<float> _dEdxParamsKaon{};
  std::vector<float> _dEdxParamsProton{};
  std::vector<float> _cost{};

  LowMomentumMuPiSeparationPID_BDTG *_mupiPID{};

  bool _basicFlg{}, _dEdxFlg{}, _showerShapesFlg{};
  int _UseBayes{};
  bool _UseMVA{};
  float _dEdxNormalization{}, _dEdxErrorFactor{};
};

#endif 
