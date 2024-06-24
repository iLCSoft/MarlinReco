#ifndef TrueJet_h
#define TrueJet_h 1

#include "marlin/Processor.h"
#include "lcio.h"

#include <EVENT/MCParticle.h>
#include <EVENT/LCCollection.h>

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
 /**  ouput collection names. The corresponding collections are filled from the internal data structues at the end of processEvent */
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
  void stdhep_reader_bug_workaround( int line94 );
  void true_lepton();
  void cluster();
  void string();
  void assign_jet(int jet1,int jet2,int this_fafp);
  void first_parton(int this_partic,int this_jet,int& first_partic,int& last_94_parent,int& nfsr,int& info,int& info2);
  int flavour(int k2) ;
  void isr() ;
  void photon();
  void grouping() ;
  LCEvent * evt{};

                             // The information in the MCParticleCollection, recast into "pyjets"
                             // format. This is the structure the processor works with internally.
                             // Note that line-numbering starts at 1, not 0

                             // getting from one to the other:
                             //  mcp_pyjets[k_py]=mcp  gets the MCParticle corresponding to a pythia-line
                             //  i_py= mcp->ext<MCPpyjet>()  gets the pythia-line of an MCParticle.

    double p[4001][6]{};     // p,E and M - direct copy from the MCParticle  
    int k[4001][9]{};        // status, pid, and history information:
                             //  1:status , 2:pdg, 3:(first) mother, 4&5 first and last daugther,
                             //  6: last mother, 7&8; colour&anti-colour code

                             // history codes (3,4,5,6) refers to the line (first index) of p and k.
                             // The  first and last daugther codes always indicate a consecutive
                             // set of lines with the same mother, while between the first and last mother
                             // lines, there can be lines *not* having the daughter in question, so
                             // a loop-n-check between thes two is always needed!

                             // status codes: 0:created in simulation, 1: stable, 11:decaying, 
                             //   21:documentation. Any of these +30: particle is overlay

                             // pdg: either the particle pdg code, or one of Pythia's internal
                             // codes: 91; cluster, i.e. a hadron created *during* the parton-shower
                             //  92: string, i.e. the object that is the parent of all hadrons at the
                             //  end of the parton-shower. It's direct parent is a group of lines,
                             //  starting and ending with quarks, with a number of colour-conected
                             //  gluons inbetween. 94: a many-to-many object used to allow to redistribute
                             //  four-momentum between it's components (total 4-mom conserved)

                             // Note than k is *not* a direct copy of the information in the MCParticle.
                             // In particular, many idiosyncarcies are corrected for - mostly concerning
                             // 94 objects. Also the (inaptly named) colour-code is sometimes also used
                             // for leptons, to indicate how these should be grouped, in case the generator
                             // information contains hints of this.

    int nlund{};             //  number of entries in this event
    MCParticleVec  mcp_pyjets{};  // Cross-reference

                             //  final results of the processor at particle level:

    int jet[4001]{};         //  Jet number of the particle 
    int companion[4001]{};   //  Associated jet (i.e. the jet that creates a colour-neutral
                             //  object together with the jet this particle belongs to)

    int current_jet{};
    int first_line{};

                             //  final results of the processor at jet level:
    // TYPE jets_summary_t
    int njet{};         // number of jets found
    int n_djb{};        // number of initial di-jets (=colour neutrals)
    int n_dje{};        // number of final di-jets (=colour neutrals), If there are no gluon-splitting in the P.S.
                        // n_djb = n_dje.
    int nstr{};         // number of strings   
    int n_hard_lepton{};// number of hard (ME) leptons
    int nboson{};       // number of bosons (normally gluons) yielding ffbar pairs *during* the P.S.
    int nclu{};         // number of "clusters", i.e. hadrons created *during* the P.S. They have two quark parents,
                        // so even if there is only one hadron created (usually, but not always the case), there will be
                        // an empty jet.
    int nisr{};         // number of ISRs
    int nphot{};        // number of M.E. photons (not counting ISR)
    int n_mixed{};      // number of "mixed" jets - jets with final particles from different true jets. should be 0
    int n_jetless{};    // number of jets with no assigned identy. should be 0
    int n_2jet_clu{};   // number of clusters actually giving two non empty jets (see above)
    int n_0_E_jets{};   // number of empty (zero energy) jets - typically from clusters (see above)
    int n_beam_jet{};   // number of beam (overlay) jets. Either 1 or 0 (is there was no overlay at all in the event)


    // TYPE jet_t         //    Index in the following is the jet-number, and when the variable refers to a particle,
                          //    the line in pyjets is implied

    int fafp_last[26]{};  // fafp = fermion-antifermion pair. _last is the last created = the first if there was no gluon splitting
                          // or from a W from eg a top decay
    int elementon[26]{};  // I just invented that word (collective for quarks, leptons and elementary bosons) 
                          // The elemeton is the particle directly before the final hadrons or leptons
                          // for hadrons, it would be the two last-generation quarks, for leptons usually simply the
                          // initial leptons.
    int boson[26]{};      // boson (gluons) inside the P.S. giving qqbars
    int fafp_boson[26]{}; // fafp of the boson creating this jet. 
    int fafp[26]{};       // first fafp of the jet = either fafp_last (jet not from boson) or fafp_boson (jet from boson)
    int nfsr[26]{};       // number of FSR photons
    int type[26]{};       // type of the jet: 1:string, 2:lepton, 3:cluster, 4:ISR, 5:overlay, 6:ME photon
    int group[26][26]{};  // group this jet belongs to. 
    int dijet_begining[26]{}; // the initial di-jet of the jet. the fafp:s of all jets in this group is the initial fafp ie. 
                              // ie the fermion anti-fermion pair the boson directly decayed to. A colour singlet.
    int dijet_end[26]{};      // idem for the final di-jet (could be the same as the initial, or could be gluon/W induced)
    double tmom[26][3]{}; // true momentum
    double tE[26]{};      // true energy.

    //TYPE dijet_t (begin)
    int jets_begin[26][26]{};    // list of jets the initial colour singlet gives rise to

    //TYPE dijet_t (end)
    int jets_end[26][26]{};     // list of jets the final colour singlet gives rise to

    bool _whiz1 = false ;
    bool _top_event = false ;
    bool _higgs_to_glue_glue = false ;
} ;

#endif



