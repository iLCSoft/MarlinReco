#ifndef CreatePDFProcessor_hh
#define CreatePDFProcessor_hh 1

#include <string>
#include <vector>
#include <marlin/Processor.h>

#include <EVENT/LCCollection.h>

#include "LikelihoodPID.hh"

#include "TFile.h"
#include "TH1.h"

using namespace lcio ;
using namespace marlin ;

//class LikelihoodPID;

class CreatePDFProcessor : public Processor{
public:
  virtual Processor*  newProcessor() { return new CreatePDFProcessor ; }
  CreatePDFProcessor();
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run);
  virtual void processEvent( LCEvent * evt );
  virtual void check( LCEvent * evt );
  virtual void end();
 
private:
  void CalculateDeltaPosition(float charge, TVector3 &p, const float* caylpos, float* delpos);

  std::string _description = {};
  std::string _PfoCollection = {};
  std::string _LinkCollection = {};

  LCCollection* _PFOCol = NULL;
  LCCollection* _LinkCol = NULL;

  std::vector<float> _dEdxParamsElectron = {};
  std::vector<float> _dEdxParamsMuon = {};
  std::vector<float> _dEdxParamsPion = {};
  std::vector<float> _dEdxParamsKaon = {};
  std::vector<float> _dEdxParamsProton = {};
  float _dEdxNormalization{}, _dEdxErrorFactor{}, _bfield{};

  LikelihoodPID *_myPID{};
  TFile* _fpdf{};
  TH1F* pidvariable[6][21]{};
  std::string _filename{};

  std::string itos(int i)
  {
    std::stringstream s;
    s << i;
    return s.str();
  }
};

#endif 
