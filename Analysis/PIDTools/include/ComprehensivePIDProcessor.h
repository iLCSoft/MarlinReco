#ifndef ComprehensivePIDProcessor_hh
#define ComprehensivePIDProcessor_hh 1

#include <marlin/MarlinConfig.h>
#include <marlin/Processor.h>
#include <InputAlgorithm.h>
#include <TrainingModel.h>
#include <string>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2I.h>
#include <TCanvas.h>
//#include <TImage.h>
#include <TStyle.h>
#include <TMVA/Reader.h>

//#ifdef MARLIN_AIDA //AIDA
//#include <marlin/AIDAProcessor.h>
//#endif

using namespace lcio;
using namespace marlin;

/* Comprehensive Particle Identification (CPID) Processor
 *
 * The CPIDProcessor gathers PID-related observables, trains a model to optimise PID determination and infers from the model to data.
 *
 */

class ComprehensivePIDProcessor : public Processor{

public:

  virtual Processor*  newProcessor() { return new ComprehensivePIDProcessor; }

  ComprehensivePIDProcessor(const ComprehensivePIDProcessor&) = delete;
  ComprehensivePIDProcessor& operator=(const ComprehensivePIDProcessor&) = delete;
  ComprehensivePIDProcessor();
  virtual ~ComprehensivePIDProcessor() = default;

  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run);
  virtual void processEvent( LCEvent* );
  virtual void check( LCEvent* evt );
  virtual void end();

  // cross check if training observables are among the available extracted ones
  bool CheckTrainingObservables(const std::vector<std::string>& trainObs, const std::vector<std::string>& compObs);

  // standard creator and reader of training reference file; created after TrainingModel::initTrain, read before TrainingModel::initInfer
  void CreateReferenceFile(int n);
  void ReadReferenceFile(int n);
  std::string ReferenceFile(int n);

  void PlotTH2(TCanvas* can, TH2* hist);
 
private:

  // processor parameters

  std::string _PFOColName{};
  std::string _RecoMCTruthLinkName{};

  bool _modeExtract=false;
  bool _modeTrain=false;
  bool _modeInfer=false;

  std::string _TTreeFileName{};
  std::vector<std::string> _inputAlgoSpecs{};
  std::vector<std::string> _trainModelSpecs{};
  std::vector<std::string> _reffile{};

  std::vector<int> _signalPDGs{};
  std::vector<int> _backgroundPDGs{};

  float _momMin=1, _momMax=100;
  bool _momLog = true;
  int _momNBins = 12;

  float _cutD0=0, _cutZ0=0;
  float _cutLamMin=0, _cutLamMax=0;
  int _cutNTracksMin=0, _cutNTracksMax=0;

  // processor-wide variables

  int _nEvt=0, _nRun=0;
  long long int _nPFO=0;

  TFile* _TTreeFile{};

  std::vector<std::string> _inputAlgoTypes{};
  std::vector<std::string> _inputAlgoNames{};
  std::vector<std::string> _trainModelTypes{};
  std::vector<std::string> _trainModelNames{};

  std::vector<cpid::InputAlgorithm*> _inputAlgorithms{};
  std::vector<cpid::TrainingModel*>  _trainingModels{};

  std::vector<std::string> _observablesNames{};
  std::vector<float> _observablesValues{};
  TTree* _observablesTree{};

  std::vector<std::string> _trainingObservables{};
  std::vector<std::vector<std::string> >_weightFiles{};
  std::vector<std::vector<std::pair<float,float> > >_weightFileBrackets{};

  int _nAlgos=0, _nModels=0, _nObs=0, _nObs_base=0;
  std::vector<int> _nObsKum{};

  std::vector<float> _momBins{};

  std::vector<int> _nMomBins{};
  std::vector<float*> _inferObsValues{};
  std::vector<int> _SnB{};

  float _momabs=0, _lambda=0;
  float _d0=0, _z0=0;
  float _PDG=0, _nTracks=0;

  const static int _nReason = 6;
  unsigned int _nRejectedPFOs[_nReason];
  TH1D* _rejectedPFOs{};

  std::string _plotFolder=".", _fileFormat=".png";
  TH2I* _PDGCheck{};

};

#endif 
