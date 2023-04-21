#ifndef ComprehensivePIDProcessor_hh
#define ComprehensivePIDProcessor_hh 1

#include <marlin/MarlinConfig.h>
#include <marlin/Processor.h>
//#include <LCIOSTLTypes.h>
#include <InputAlgorithm.h>
#include <string>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMVA/Reader.h>

#ifdef MARLIN_AIDA //AIDA
#include <marlin/AIDAProcessor.h>
#endif

using namespace lcio;
using namespace marlin;


class ComprehensivePIDProcessor : public Processor{

public:

  virtual Processor*  newProcessor() { return new ComprehensivePIDProcessor; }

  ComprehensivePIDProcessor(const ComprehensivePIDProcessor&) = delete;
  ComprehensivePIDProcessor& operator=(const ComprehensivePIDProcessor&) = delete;
  ComprehensivePIDProcessor();
  virtual ~ComprehensivePIDProcessor() = default;

  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run);
  virtual void processEvent( LCEvent* evt );
  virtual void check( LCEvent* evt );
  virtual void end();

  //virtual void RegisterBranch();
  void RegisterBranchD(std::vector<TTree*>& treeVec, double* var, const char* name);
  void RegisterBranchF(std::vector<TTree*>& treeVec, float* var, const char* name);
  void RegisterBranchI(std::vector<TTree*>& treeVec, int* var, const char* name);
 
private:

  int _mode=0;
  int _mode2nBins=0;
  std::string _weightsfolder{};
  std::string _reffile{};

  std::string _PFOColName{};
  std::string _RecoMCTruthLinkName{};

  int _nEvt=0, _nRun=0;
  long long int _nPFO=0;
  std::vector<int> _signalPDGs{};
  std::vector<int> _backgroundPDGs{};
  std::vector<std::string> _inputAlgoSpecs{};
  std::vector<std::string> _inputAlgoTypes{};
  std::vector<std::string> _inputAlgoNames{};

  std::vector<float> _inputAlgo1ParamsF{};
  std::vector<std::string> _inputAlgo1ParamsS{};

  std::vector<cpid::InputAlgorithm*> _inputAlgorithms{};

  std::vector<std::string> _observablesNames{};
  std::vector<float> _observablesValues{};
  TTree* _observablesTree{};

  std::vector<std::pair<float,float> > _weightFiles{};

  int _nAlgos=0, _nObs=0, _nObs_base=0;
  std::vector<int> _nObsKum{};

  //TTree* _observablesTree_pions{};
  //TTree* _observablesTree_kaons{};

  std::vector<TTree*> _observablesTrees{};
  TMVA::Reader* _TReader{};

  float _momabs=0, _lambda=0;
  float _d0=0, _z0=0;
  float _PDG=0, _nTracks=0;

  const static int _nReason = 6;
  double _nRejectedPFOs[_nReason];
  TH1D* _rejectedPFOs{};

  const static int _nPart = 5;
  double _nWrongMCPDG[_nPart];
  TH1D* _wrongMCPDG{};


  //std::string _outputRootFile{};

};

#endif 
