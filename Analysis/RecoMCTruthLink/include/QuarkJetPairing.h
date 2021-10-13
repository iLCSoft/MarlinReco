#ifndef QuarkJetPairing_h
#define QuarkJetPairing_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>

#include "IMPL/LCCollectionVec.h"

using namespace lcio ;
using namespace marlin ;

/** Processor for marlin to find the smallest deviation between reconstructed jets and jet initialized quarks. Here on basis of WW->4f hadronic decays at 500 GeV from the DBD database
 * The processor looks for the smallest sum of angle (scalar product) between quarks and jets and declare it as the best permutation. Additional a theta cut is implemented to fullfill detector constraints
 *
 * @author B.Hermberg, DESY
 */

class QuarkJetPairing : public Processor {
 public:
  virtual Processor*  newProcessor() { return (new QuarkJetPairing) ; }
  QuarkJetPairing() ;
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ;
  virtual void check( LCEvent * evt ) ;
  virtual void end() ;

 protected:
  std::string _inputJetCollection{};
  std::string _inputMCPCollection{};

  int ipair{}; 
  int iperm{};
  int index[4]{};

 
  int jet_array[4]{};
  int quark_array[4]{};

  enum{nJETS=4};
  enum{NPERM=6};

  int permutation[nJETS]{};


  // float jet_ptot[4]{};
  // float jet_px[4]{};
  // float jet_py[4]{};
  // float jet_pz[4]{};
  // float jet_ene[4]{};

  // float quark_ptot[4]{};
  // float quark_px[4]{};
  // float quark_py[4]{};
  // float quark_pz[4]{};
  // float quark_ene[4]{};

  float phi_jet_quark[4][4]{};
  float momentum[3]{}, jetenergy{};
  int _nRun{};
  int _nEvt{};

} ;

#endif
