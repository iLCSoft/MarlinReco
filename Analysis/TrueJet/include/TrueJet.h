#ifndef TrueJet_h
#define TrueJet_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/LCRelation.h>
#include <UTIL/LCRelationNavigator.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/ParticleIDImpl.h>
#include <string>


using namespace lcio ;
using namespace marlin ;


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
 * @version $Id: TrueJet.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class TrueJet : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TrueJet ; }
  
  
  TrueJet() ;

  TrueJet& operator=(const TrueJet&) = delete ;
  TrueJet(const TrueJet&) = delete ;  
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt) ; 
  

  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  void readback( LCEvent * evt);


 protected:

  /** Input collection name.
   */
  std::string  _MCParticleColllectionName{};
  std::string _recoParticleCollectionName{};
  std::string _recoMCTruthLink{};
 /**  ouput collection names */
  std::string _trueJetCollectionName{};
  std::string _finalColourNeutralCollectionName{};
  std::string _initialColourNeutralCollectionName{};
  std::string _trueJetPFOLink{};
  std::string _trueJetMCParticleLink{};
  std::string _finalElementonLink{};
  std::string _initialElementonLink{};
  std::string _finalColourNeutralLink{};
  std::string _initialColourNeutralLink{};
  int _nRun{};
  int _nEvt{};

 
private:

  void getPyjets(LCCollection* mcpcol);

  void true_lepton();
  void cluster();
  void string();
  void assign_jet(int jet1,int jet2,int this_fafp);
  void first_parton(int this_partic,int this_jet,int& first_partic,int& last_94_parent,int& nfsr,int& info,int& info2);
  int flavour(int k2) ;
  void fix94() ;
  void isr() ;
  void grouping() ;
  void fix_top();
  LCEvent * evt{};
  MCParticleVec  mcp_pyjets{};

    int jet[4001]{};
    int companion[4001]{};
    double p[4001][6]{};
    int k[4001][7]{};
    int nlund{};
    int i_jet{};
    int first_line{};

    // TYPE jets_summary_t
    int njet{};
    int n_djb{};
    int n_dje{};
    int nstr{};  // *2
    int n_hard_lepton{};
    int nboson{};
    int nclu{};  // *2
    int nisr{};
    int n_mixed{};
    int n_jetless{};
    int n_2jet_clu{};
    int n_0_E_jets{};
    int n_beam_jet{};

    // TYPE jet_t
    int fafp_last[26]{};  // fafp = fermion-antifermion pair. _last is the last created = the first if there was no gluon splitting
                        // or from a W from eg a top decay
    int elementon[26]{};  // I just invented that word (collective for quarks, leptons and elementary bosons) 
    int boson[26]{};
    int fafp_boson[26]{};  // fafp of the boson creating this jet. 
    int fafp[26]{};       // first fafp of the jet = either fafp_last (jet not from boson) or fafp_boson (jet from boson)
    int nfsr[26]{};
    int type[26]{};
    int group[26][26]{};  // group this jet belongs to. 
    int dijet_begining[26]{}; // the initial di-jet of the jet. the fafp:s of all jets in this group is the initial fafp ie. ie the boson
                             // kids (colour singlet!)
    int dijet_end[26]{};       // idem for the final di-jet (could be the same as the initial, or could be gluon/W induced)
    // pjet(5) =
    double tmom[26][3]{};
    double tE[26]{};

    //TYPE dijet_t (begin)
    int jets_begin[26][26]{};
    //singlet_four_p(5)
    //psum_four_p(5)

    //TYPE dijet_t (end)
    int jets_end[26][26]{};
    //singlet_four_p(5)
    //psum_four_p(5)

} ;

#endif



