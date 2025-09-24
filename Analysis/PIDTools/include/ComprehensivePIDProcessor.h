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
 * The processors uses modules to achieve this: InputAlgorithms to collect PID information, TrainingModels to train and infer.
 * The modules are made available via dynamic loading, using a base class and a Mgr.
 * They are then steered in the processor steering file.
 * For this, the modules to be used need to be specified in the _inputAlgoSpecs and _trainModelSpecs, respectively.
 * Specify [type]:[name] of a module, or only [type] in which case name=type.
 * The processor then checks if float and/or string vectors are specified as processor parameters called [name].F and [name].S, respectively.
 * The function of these parameters is specified in each individual module description.
 * The modules are created and initialised in init.
 *
 * The modes in which the processor can run are steered via the _mode flags.
 * It extracts the information from the PFOs via the InputAlgorithms in processEvent.
 * These are put in a TTree and can then be stored in a root file.
 * If specified, training is run in end. Training can be run with the previously generated root file only, i.e. without extraction.
 * If specified, inference is run in processEvent. It needs extraction.
 * Training and inference cannot be done simultaneously.
 *
 * The training creates a reference file in which the used observables and the used signal and background PDGs are stored, as well as the location of the model's weight files together with the corresponding momentum bin.
 * The inference reads the reference file and loads the TrainingModel based on the information therein.
 *
 * For each observable in each InputAlgorithm one branch is created in the TTree called [algorithm_name]_[observable_name].
 * In addition, the processor also collects for each PFO: momentum, lambda angle (polar angle relative to z=0), the PDG of the likeliest connected MCParticle, the number of associated tracks as well as d0 and z0 of the first track (if any).
 * The training observables can be specified. If not done so, the defaults are all observables associated with the specified InputAlgorithms + momentum and lambda.
 *
 * A number of PFO acceptance cuts can be applied, based on its properties or that of its first track.
 * Signal and background PDGs can be defined. Only PFOs which originate from either of these are used.
 * The consequence of belonging to either of these groups lies in the details of the individual TrainingModels.
 * The PFOs are divided into momentum bins, for each of which a separate TrainingModel is created and trained/inferred from.
 * If inference is specified, confusion matrix plots of all signal PDGs can be created in the specified folder during end.
 *
 * @param _PFOColName - Name of the PFO input collection (ReconstructedParticleImpl).
 *    string, default: "PandoraPFOs".
 * @param _RecoMCTruthLinkName - Name of the link from PFOs to MCParticles input collection (LCRelation).
 *    string, default: "RecoMCTruthLink".
 *
 * @param _modeExtract - Set true to extract PID observables via the specified InputAlgorithms.
 *    bool, default: false.
 * @param _modeTrain - Set true to train via the specified TrainingModels. Cannot be true at the same time as _modeInfer.
 *    bool, default: false.
 * @param _modeInfer - Set true to infer PID from the specified TrainingModels. Cannot be true at the same time as _modeTrain.
 *    bool, default: false.
 *
 * @param _TTreeFileName - Name of the root file in which the TTree with all observables is stored; in case of extraction it is an optional output with no output if left empty, otherwise it is a necessary input.
 *    string, default: "TTreeFile.root".
 * @param _inputAlgoSpecs - List of input algorithms; for each specify type:name or only type (then name=type).
 *    string vector, default: {}.
 * @param _trainModelSpecs - List of training models; for each specify type:name or only type (then name=type).
 *    string vector, default: {}.
 * @param _reffile - Reference file(s). If only one file but several training models are specified the reference files are auto-numbered.
 *    string vector, default: {"Ref.txt"}.
 * @param _trainingObservables - List of observables that should be used for traning. If empty, all observables from the specified algorithms + momabs + lambda are used.
 *    string vector, default: {}.
 *
 * @param _signalPDGs - List of PDG numbers that are considered signal.
 *    int vector, default: {11,13,211,321,2212}.
 * @param _backgroundPDGs - List of PDG numbers that are considered background.
 *    int vector, default: {}.
 *
 * @param _plotFolder - Folder in which the automatic confusion matrix plots of inference will be put, is created if not already existing; if empty, no plots are created.
 *    string, default: "CPID_Plots".
 * @param _fileFormat - File format of the automatic confusion matrix plots of inference.
 *    string, default: ".png".
 *
 * @param _momMin - For momentum bins: minimum momentum / GeV.
 *    float, default: 1.
 * @param _momMax - For momentum bins: maximum momentum / GeV.
 *    float, default: 100
 * @param _momLog - For momentum bins: should the momentum bins be logarithmic.
 *    bool, default: true.
 * @param _momNBins - For momentum bins: number of momentum bins.
 *    int, default: 12.
 *
 * @param _cutD0 - PFOs whose first track have a d0 larger than the given value will be ignored; set to 0 to accept all particles.
 *    float, default: 0.
 * @param _cutZ0 - PFOs whose first track have a z0 larger than the given value will be ignored; set to 0 to accept all particles.
 *    float, default: 0.
 * @param _cutLamMin - PFOs whose first track have an angle lambda (relative to the cathode) smaller than the given value will be ignored; set to 0 to accept all particles.
 *    float, default: 0.
 * @param _cutLamMax - PFOs whose first track have an angle lambda (relative to the cathode) larger than the given value will be ignored; set to 0 to accept all particles.
 *    float, default: 0.
 * @param _cutNTracksMin - PFOs with fewer (<) tracks than the given value are ignored; set to -1 to accept all PFOs.
 *    int, default: -1.
 * @param _cutNTracksMax - PFOs with more (>) tracks than the given value are ignored; set to -1 to accept all PFOs.
 *    int, default: -1.
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

  void PlotTH2(TCanvas* can, TH2* hist, int effpur=0);
 
private:

  // processor parameters

  std::string _PFOColName{};
  std::string _RecoMCTruthLinkName{};

  bool _modeExtract=false;
  bool _modeTrain=false;
  bool _modeInfer=false;
  bool _addMCPID=false;

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
  bool _writeTTreeFile = false;

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
  std::vector<TH2I*> _PDGCheck{};
  std::vector<std::vector<TH2I*> > _PDGChecks{};

};

#endif 
