#ifndef CLICPfoSelectorAnalysis_h
#define CLICPfoSelectorAnalysis_h 1

#include "PfoUtilities.h"
#include "marlin/Processor.h"
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCRelation.h>
#include <UTIL/LCRelationNavigator.h>
#include "lcio.h"
#include <string>

#include "TTree.h"
#include "TGraph.h"
#include "TH1F.h"

using namespace lcio ;
using namespace marlin ;

using namespace std;

/**  CLICPfoSelectorAnalysis processor
 * 
 *  Run on the PFO input collection and create:
 * - a TTree with the PFO variables used in the CLICPfoSelector
 * - cluster time vs pT graphs for each particle category and region
 * - PFO energy sum histos for each particle category and region
 * Possibility to detect if the PFO belongs to signal/overlay
 * Possibility to check if the track and the cluster belonging to
 * the same PFO were produced by at least one common MCParticle
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of ReconstructedParticles.
 *  Needs the collection of MCParticles.
 *  Needs the collection of LCRelations - to do the match track/cluster.
 *
 *  <h4>Output</h4> 
 *  A TTree.
 *  Time vs pT graphs.
 *  Energy histos.
 * 
 */

class CLICPfoSelectorAnalysis : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CLICPfoSelectorAnalysis ; }
  
  CLICPfoSelectorAnalysis() ;
  
  //initializing the variables in the TTree  
  virtual void init() ;

  virtual void processRunHeader( LCRunHeader* run ) ;

  virtual void processEvent( LCEvent * evt ) ; 

  virtual void end() ;

  //filling the TTree  
  void fillTree(LCEvent * evt, string collName); 

  //filling the graphs and histos  
  void fillPlots();

 protected:

  // Input collection name
  string colNamePFOs{};
  string treeName{};
  float cutCosTheta = 0.975;
  int minECalHits = 5;
  int minHcalEndcapHits = 5;
  float forwardCosThetaForHighEnergyNeutralHadrons = 0.95, forwardHighEnergyNeutralHadronsEnergy = 10.0;
  bool analyzePhotons = true, analyzeChargedPfos = true, analyzeNeutralHadrons = true;
  bool analyzeAll = true, analyzeSignal = true, analyzeOverlay = true;

  int _nRun{};
  int _nEvt{};

  //Variables in the TTree
  TTree *pfo_tree = NULL;
  int eventNumber = 0, runNumber = 0;
  int type = 0;
  double p = 0.0, px = 0.0, py = 0.0, pz = 0.0, pT = 0.0;
  double costheta = 0.0, energy = 0.0, mass = 0.0, charge = 0.0;
  int nTracks = 0, nClusters = 0;
  double clusterTime = 0.0, clusterTimeEcal = 0.0, clusterTimeHcalEndcap = 0.0;
  int nCaloHits = 0, nEcalHits = 0, nHcalEndCapHits = 0;
  int trk_clu_sameMCPart = 0, atLeastOneSignal = 0;

  //List of plots
  vector<string> particleCategories{};
  vector<string> generationCategories{};
  map<string,TGraph*> g_timeVsPt{};
  map<string,TGraph*> g_timeVsPt_central{};
  map<string,TGraph*> g_timeVsPt_forward{};
  TH1F* h_energy_tot{};
  TH1F* h_energy_tot_signal{};
  TH1F* h_energy_tot_background{};
  map<string,TH1F*> h_energy{};
  map<string,TH1F*> h_energy_central{};
  map<string,TH1F*> h_energy_forward{};
  map<string,double> energy_tot{};
  map<string,double> energy_tot_central{};
  map<string,double> energy_tot_forward{};
  float en_min = 0.0, en_max = 500;

  //MC particles collections
  string m_inputPhysicsParticleCollection{};
  string m_recoMCTruthLink{};
  string m_SiTracksMCTruthLink{};
  string m_ClusterMCTruthLink{};
  vector<MCParticle*> physicsParticles{};

} ;

#endif



