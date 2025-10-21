#ifndef LeptonIDProcessor_h
#define LeptonIDProcessor_h 1

#include <lcio.h>
#include <map>
#include <marlin/EventModifier.h>
#include <marlin/Processor.h>
#include <string>

#include <EVENT/Cluster.h>
#include <EVENT/MCParticle.h>

#include <marlinutil/WeightedPoints3D.h>

#include <TMVA/Reader.h>
#include <TTree.h>

using namespace lcio;
using namespace marlin;

/**  Example processor for marlin.
 *
 *  If compiled with MARLIN_USE_AIDA
 *  it creates a histogram (cloud) of the MCParticle energies.
 *
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4>
 *  A histogram.
 *
 * @param CollectionName Name of the MCParticle collection
 *
 * @author F. Gaede, DESY
 * @version $Id: MyProcessor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $
 */

class LeptonIDProcessor : public Processor //, public EventModifier
{

public:
  virtual Processor* newProcessor() { return new LeptonIDProcessor; }

  LeptonIDProcessor();

  LeptonIDProcessor(const LeptonIDProcessor&) = delete;
  LeptonIDProcessor& operator=(const LeptonIDProcessor&) = delete;

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

  virtual const std::string& name() const { return Processor::name(); }

  virtual void check(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

protected:
  void clear();
  void getShowerShapes(Cluster* clu);
  void getClusterShapes(WeightedPoints3D& wgpt);
  void setupTree();
  void setupMVAReader();

  /** Input collection name.
   */
  std::string _PandoraPFOsColName{};
  std::string _InputPFOsColName{};
  std::string _RecoMCTruthLinkName{};
  std::string _MCTruthRecoLinkName{};
  std::string _PandoraClustersColName{};

  std::string _weightfile{};
  std::string _mvaname{};

  std::string _pidname{};
  std::string _dEdxname{};

  std::vector<std::string> _usedVariables{};

  bool _buildTree = true;
  bool _evalMVA = false;

  float _MCTruthRecoCweightCut{};
  float _MCTruthRecoTweightCut{};
  float _RecoMCTruthCweightCut{};
  float _RecoMCTruthTweightCut{};

  TTree* _tree{};
  TMVA::Reader _mvareader;

  int _truePDG{};
  int _parentPDG{};
  int _pandoraPDG{};
  int _likelihoodPDG{};

  std::map<std::string, float_t> _vars = {
      {"seenP", 0.0},
      {"seenLogP", 0.0},
      {"seenCosTheta", 0.0},
      {"seenE", 0.0},
      {"seenEcalDep", 0.0},
      {"seenHcalDep", 0.0},
      {"seenYokeDep", 0.0},
      {"e_over_p", 0.0},
      {"ecal_share", 0.0},
      {"dEdxDist_e", 0.0},
      {"shape0", 0.0},
      {"shape1", 0.0},
      {"shape2", 0.0},
      {"shape3", 0.0},
      {"shape4", 0.0},
      {"shape5", 0.0},
      {"shape6", 0.0},
      {"shape7", 0.0},
      {"shape8", 0.0},
      {"shape9", 0.0},
      {"shape10", 0.0},
      {"shape11", 0.0},
      {"shape12", 0.0},
      {"shape13", 0.0},
      {"shape14", 0.0},
      {"shape15", 0.0},
      {"shape16", 0.0},
      {"shape17", 0.0},
      {"shape18", 0.0},
      {"shape19", 0.0},
      {"shape20", 0.0},
      {"shape21", 0.0},
      {"shape22", 0.0},
      {"shape23", 0.0},
      {"cluEllipsoid_r1", 0.0},
      {"cluEllipsoid_r2", 0.0},
      {"cluEllipsoid_r3", 0.0},
      {"cluEllipsoid_vol", 0.0},
      {"cluEllipsoid_r_ave", 0.0},
      {"cluEllipsoid_density", 0.0},
      {"cluEllipsoid_eccentricity_T", 0.0},
      {"cluEllipsoid_eccentricity_L", 0.0},
  };
};

#endif
