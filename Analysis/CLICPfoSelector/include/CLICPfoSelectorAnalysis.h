#ifndef CLICPfoSelectorAnalysis_h
#define CLICPfoSelectorAnalysis_h 1

#include "marlin/Processor.h"
#include <EVENT/ReconstructedParticle.h>
#include "lcio.h"
#include <string>

#include "TTree.h"
#include "TGraph.h"
#include "TH1F.h"

using namespace lcio ;
using namespace marlin ;

using namespace std;

#define FORMATTED_OUTPUT_TRACK_CLUSTER(out, N1, E1,E2,E3,N2,E4,N3,E5,E6,E7) \
    out <<                                                                                      \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E3        <<                                   \
    std::right << std::setw(widthInt)      <<    N2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E4        <<                                   \
    std::right << std::setw(widthInt  )    <<    N3        <<                                   \
    std::right << std::setw(widthFloat)    <<    E5        <<                                   \
    std::right << std::setw(widthFloat)    <<    E6        <<                                   \
    std::right << std::setw(widthFloat)    <<    E7  << std::endl

/**  CLICPfoSelectorAnalysis processor
 * 
 *  It creates a TTree with the PFO variables used in the CLIC selection.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of ReconstructedParticles.
 *
 *  <h4>Output</h4> 
 *  A TTree.
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

  void fillScatterPlots();

 protected:

  // Input collection name
  string colNamePFOs{};
  string treeName{};
  float cutCosTheta;
  bool analyzePhotons, analyzeChargedPfos, analyzeNeutralHadrons;

  int _nRun{};
  int _nEvt{};

  //Variables in the TTree
  TTree *pfo_tree = NULL;
  int type = 0;
  double p = 0.0, px = 0.0, py = 0.0, pz = 0.0, pT = 0.0;
  double costheta = 0.0, energy = 0.0, mass = 0.0, charge = 0.0;
  int nTracks = 0, nClusters = 0;
  double clusterTime = 0.0, clusterTimeEcal = 0.0, clusterTimeHcalEndcap = 0.0;
  int nCaloHits = 0, nEcalHits = 0, nHcalEndCapHits = 0;

  int eventNumber = 0, nPartMC = 0, nPartPFO = 0.0;

  //List of scatter plots
  vector<string> particleCategories;
  map<string,TGraph*> timeVsPt;
  map<string,TGraph*> timeVsPt_barrel;
  map<string,TGraph*> timeVsPt_endcap;

} ;

#endif



