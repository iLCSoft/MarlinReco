#ifndef EvaluateTauFinder_h
#define EvaluateTauFinder_h 1

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TTree.h"
#include "UTIL/LCRelationNavigator.h"
#include "lcio.h"
#include "marlin/Processor.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <string>

using namespace lcio;
using namespace marlin;

struct MyParticle;

/**  Evaluation processor for TauFinder
 *
 * @author A. Muennich, CERN
 */

class EvaluateTauFinder : public Processor {

public:
  virtual Processor* newProcessor() { return new EvaluateTauFinder; }

  EvaluateTauFinder();
  EvaluateTauFinder(const EvaluateTauFinder&) = delete;
  EvaluateTauFinder& operator=(const EvaluateTauFinder&) = delete;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();
  virtual void LoopDaughters(MCParticle* particle, double& Evis, double& ptvis, double& pvis);
  virtual void LoopDaughtersRelation(MCParticle* particle, LCRelationNavigator* relationNavigatorTau,
                                     LCRelationNavigator* relationNavigatorMC, bool& ralToTau);

protected:
  /** Input collection name.
   */
  std::string _colNameMC{}, _colNameTrack{};
  std::string _colNameMCTruth{}, _incol{}, _colNamePFORecLink{};
  std::string _colNameMCRecLink{}, _colNameTracksRecLink{}, _colNameTauRecLink{};
  std::string _OutputFile_Signal{};
  int _nRun = -1;
  int _nEvt = -1;

  double _ntot_rec = 0.0;
  double _ntot_mc = 0.0;
  double _ntau_correct = 0.0;
  double _dEsum = 0.0;
  double _dEsumsq = 0.0;
  int _ndE = 0.0;

  float _bField = 0.0;
  TFile* rootfile = NULL;
  TNtuple* leptons = NULL;
  TNtuple* tautuple = NULL;
  TNtuple* mcmisstuple = NULL;
  TNtuple* taumatchtuple = NULL;
  TNtuple* tauexacttuple = NULL;
  TNtuple* evtuple = NULL;
  TNtuple* faketuple = NULL;
  TNtuple* topofaketuple = NULL;
};

#endif
