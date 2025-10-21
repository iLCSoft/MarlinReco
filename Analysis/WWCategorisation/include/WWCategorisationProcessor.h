#ifndef WWCategorisationProcessor_h
#define WWCategorisationProcessor_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"

using namespace lcio;
using namespace marlin;

class WWCategorisationProcessor : public Processor {

  /*
   * The WWCategorisationProcessor categorises an event into a true category as well as a basic and an advanced
   * reconstructed category. These classify W(W) (i.e. including 'single-W') events into 6 categories: 0 - hadronic
   * (both W decay hadronically) 1 - invisible semileptonic (one W decays hadronically, but the lepton is not within
   * detector acceptance/detected) 2 - semileptonic electron 3 - semileptonic muon 4 - semileptonic tauon 5 - leptonic
   * (both W decay leptonically)
   *
   * The categories as well as the 4 observables used to derive the reconstructed categories
   * (invariant mass, number of non-isolated PFOs, missing pT, missing energy) are stored as event parameters:
   * WWCategorisation.TrueCat, .RecoCatBasic, .RecoCatAdvanced, .nPFO (all int), .mInv, .misspT, .missE (all float).
   * The angular cut, above which a semileptonic electron is considered outside of the detector acceptance, can be
   * adjusted via the _MCElCat_CosThCut parameter.
   *
   * The true category does not give relevant results for event that are not W(W).
   * The basic reconstruction categorises via number and flavour of isolated leptons + one invariant mass cut for the
   * hadronic channel. The advanced reconstruction uses a combination of cuts on the 4 observables. For details see
   * Master thesis by A. Silva (work in progress).
   *
   * Input collections are MCParticles as well as the PFOs split into isolated particles (e, mu, tau, gamma) and
   * non-isolated PFOs. In addition, if corresponding file names are provided, a TTree of the event observables is
   * written into a TFile and the cofusion matrix (advanced reco vs. true) is printed.
   *
   * @author A. Silva and U. Einhaus, DESY
   * @date 05/2024
   */

public:
  virtual Processor* newProcessor() { return new WWCategorisationProcessor; }

  WWCategorisationProcessor();

  WWCategorisationProcessor(const WWCategorisationProcessor&) = delete;
  WWCategorisationProcessor& operator=(const WWCategorisationProcessor&) = delete;

  virtual void init();
  virtual void processRunHeader(LCRunHeader*);
  virtual void processEvent(LCEvent* evt);
  virtual void check(LCEvent*);
  virtual void end();

  // counting leptons in an array
  int CountingLeptons(std::vector<int>& daughters_pdg_vector);

  // create a plot of the confusion matrix in the pwd
  virtual void PlotConfusionMatrix(TCanvas* canvas, TH2* histogram);

  static bool IsLepton(int pdg);

protected:
  // collections
  std::string _MCParColName{};
  std::string _ElectronColName{};
  std::string _MuonColName{};
  std::string _TauColName{};
  std::string _PhotonColName{};
  std::string _PFOsMinusIsolatedObjetcs{};

  // processor parameters
  float _MCElCat_CosThCut{};
  std::string _TTreeFileName{};
  std::string _ConfusionMatrixFileName{};

  // flags for creation of files for TTree and ConfusionMatrix
  bool _doTT = false, _doCM = false;

  // root objects
  TFile* _TTreeFile{};
  TTree* _observablesTree{};
  TH2F* _confusion_matrix{};
  TH2F* _confusion_matrix_ef{};

  // other parameters

  static const int _n_cat = 7;
  std::map<int, int> _sl_subcat = {{11, 2}, {13, 3}, {15, 4}};

  long long int _nTrue[_n_cat];
  long long int _nReco[_n_cat];

  float _mInv = 0, _misspT = 0, _missE = 0;
  int _nPFO = 0;
  int _true_cat = 0, _reco_cat_basic = 0, _reco_cat_advanced = 0;

  int _nRun{};
  int _nEvt{};
};

#endif
