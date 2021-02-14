#ifndef TrueJet_Parser_h
#define TrueJet_Parser_h 1

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

struct MCPpyjet : LCIntExtension<MCPpyjet> {};
struct JetIndex : LCIntExtension<JetIndex> {};
struct IcnIndex : LCIntExtension<IcnIndex> {};
struct FcnIndex : LCIntExtension<FcnIndex> {};
// LCRelationNavigator* reltrue_tj =0;

class TrueJet_Parser {
  
 public:
  
  
  TrueJet_Parser()   ;
  virtual ~TrueJet_Parser() ;

  // These two lines avoid frequent compiler warnings when using -Weffc++
  TrueJet_Parser( const TrueJet_Parser& ) = delete;
  TrueJet_Parser& operator=( const TrueJet_Parser& ) = delete;

  virtual    std::string get_recoMCTruthLink(){ return _recoMCTruthLink  ;};
 
  ReconstructedParticleVec* getJets();       // get all true jets in event
  ReconstructedParticleVec* getFinalcn();    // get all final colour neutrals in event
  ReconstructedParticleVec* getInitialcn();  // get all initial colour neutrals in event

  void getall(LCEvent* event)  ;  // Convenient method to get all collections needed: no need to call
                                  // the three above individually. Also all navigators are set up
                                  // with this call - See below.

  void delall() ;                 // Tidy up, to be called at end of each event.

  int njets() { return tjcol->getNumberOfElements(); };  
                                  // Get the total number of true jets in the event

  int jetindex(ReconstructedParticle* jet ) { return jet->ext<JetIndex>() ; };
                                  // Get the index of jet in various arrays. That is, ijet in the
                                  // methods below.

  const ReconstructedParticle* jet(int ijet) {return jets->at(ijet); }
                                  // Get the jet object that is jet-number ijet.

  double Eseen(int ijet) { return jets->at(ijet)->getEnergy(); };
                                  // Seen energy of jet ijet

  double Mseen(int ijet) { return jets->at(ijet)->getMass(); };
                                  // Seen mass of jet ijet

  const double* pseen(int ijet) { return jets->at(ijet)->getMomentum(); };
                                  // Seen 3-momentum of jet ijet

  const double* p4seen(int ijet) ;
                                  // Seen 4-momentum of jet ijet (component 0 is E)

  const ReconstructedParticleVec& seen_partics(int ijet) { return jets->at(ijet)->getParticles(); };
                                  // Get all PFOs in jet ijet. 

  const MCParticleVec& true_partics(int ijet) ;
                                  // Get all MCPs in jet ijet.

  int type_jet(int ijet) { return jets->at(ijet)->getParticleIDs()[0]->getType() ; };
                                  // Get the type jet ijet: 1 = hadronic, from string
                                  //                        2 = leptonic
                                  //                        3 = hadronic, from cluster
                                  //                        4 = ISR
                                  //                        5 = overlay
                                  //                        6 = M.E. photon
                                  // If the jet is completely unseen, the type is negated.
                                  // If the jet originates from a gluon-splitting, the type is
                                  // modified to type + (jet# radiating the gluon)*100

  double Etrue(int ijet) ; 
                                  // True energy of jet ijet

  double Mtrue(int ijet) ;
                                  // True mass of jet ijet

  const double* ptrue(int ijet) ;
                                  // True 3-momentum of jet ijet

  const double* p4true(int ijet) ;
                                  // True 4-momentum of jet ijet (component 0 is E)

  double Equark(int ijet) {return ( final_elementon(ijet) != NULL ? final_elementon(ijet)->getEnergy() : 0.);}; 
                                  // Energy of the last quark/lepton before "hadronisation" jet ijet

  double Mquark(int ijet) {return ( final_elementon(ijet) != NULL ? final_elementon(ijet)->getMass() : 0 );};
                                 // Mass of the last quark/lepton before "hadronisation" jet ijet

  const double* pquark(int ijet);
                                 // 3-momentum of the last quark/lepton before "hadronisation" jet ijet

  const double* p4quark(int ijet) ;
                                 // 4-momentum of the last quark/lepton before "hadronisation" jet ijet 

  double Etrueseen(int ijet) ; 
                                  // Sum of true energy of all the particles of jet ijet, that were detected

  double Mtrueseen(int ijet) ;
                                  // Sum of true mass of all the particles of jet ijet, that were detected

  const double* ptrueseen(int ijet) ;
                                  // Sum of true 3-momentum of all the particles of jet ijet, that were detected

  const double* p4trueseen(int ijet) ;
                                  // Sum of true 4-momentum of all the particles of jet ijet, that were detected 
                                  // (component 0 is E)

  const IntVec& final_siblings( int ijet ); 
                                  // list of jets grouped with jet ijet at the end of the parton-shower

  const IntVec& initial_siblings( int ijet ); 
                                  // list of jets grouped with jet ijet seen from the  beginning of the parton-shower
                                  // (i.e. those that, together with jet ijet, belong to the same initial colour neutral)

  int mcpjet( MCParticle* mcp);
  int mcpicn( MCParticle* mcp);
  int mcpfcn( MCParticle* mcp);
                                  // Get jetnumber, initial or final colour neutral of true particle
                                  // mcp

  int recojet( ReconstructedParticle* reco);
  int recoicn( ReconstructedParticle* reco);
  int recofcn( ReconstructedParticle* reco);
                                  // Get jetnumber, initial or final colour neutral of reconstructed 
                                  // particle reco.


  int final_cn( int ijet ); 
                                  // the id of the final colour-neutral that jet ijet comes form

  int initial_cn( int ijet ); 
                                  // the id of the initial colour-neutral that jet ijet comes form

  const IntVec& jets_of_final_cn( int ifcn ); 
                                  // the list of the jets final colour-neutral ifcn gives rise to

  const IntVec& jets_of_initial_cn( int iicn ); 
                                  // the list of the jets initial colour-neutral iicn gives rise to

  int nicn() { return icncol->getNumberOfElements(); };
                                  // Number of initial colour neutrals

  int type_icn_parent(int iicn) { return initialcns->at(iicn)->getParticleIDs()[0]->getType() ; };
                                 // type of initial colour neutral iicn:  1 or 3 : c.n. is a quark pair
                                 //                                            2 : c.n. is a lepton pair
                                 //                                            4 : c.n. is an ISR
                                 //                                            6 : c.n. is a ME photon
                                 

  const IntVec& type_icn_comps(int iicn) ;
                                 // types of each of the quarks/leptons constituting this
                                 // initial colour neutral = type of the corresponding jet

  int pdg_icn_parent(int iicn) { return initialcns->at(iicn)->getParticleIDs()[0]->getPDG() ; };
                                 // PDG of parent of initial colour neutral iicn, i.e. the boson (23=Z, 24=W, 25=H ...)
  
  const IntVec& pdg_icn_comps(int iicn) ;
                                 // PDGs of each of the quarks/leptons constituting this
                                 // initial colour neutral

  double E_icn(int iicn) { return initialcns->at(iicn)->getEnergy();};
                                 // Energy of initial colour neutral iicn.

  double M_icn(int iicn) {; return initialcns->at(iicn)->getMass();};
                                 // Mass of initial colour neutral iicn.

  const double* p_icn(int iicn){ return initialcns->at(iicn)->getMomentum();} ;
                                 // 3-momentum of initial colour neutral iicn.

  const double* p4_icn(int iicn) ;
                                 // 4-momentum  of initial colour neutral iicn (component 0 is E).

  int nfcn() { return fcncol->getNumberOfElements(); };
                                  // Number of final colour neutrals

  int type_fcn_parent(int ifcn) { return finalcns->at(ifcn)->getParticleIDs()[0]->getType() ; };
                                 // type of parent of final colour neutral ifcn. Same as that of the
                                 // two jets constituting the c.n. (bare type, i.e. always positive,
                                 // no gluon-splitting indication)

  const IntVec& type_fcn_comps(int ifcn) ;
                                 // types of the quarks/leptons constituting this final colour neutral.
                                 // Equal to the type of the corresponding jet.

  int pdg_fcn_parent(int ifcn) { return finalcns->at(ifcn)->getParticleIDs()[0]->getPDG() ; };
                                 // PDG of parent of final colour neutral ifcn:  (92=string, 91=cluster,
                                 // 22=photon (ISR or ME), any lepton PDG=lepton pair)

  const IntVec& pdg_fcn_comps(int ifcn) ;
                                 // PDGs of the quarks/leptons constituting this final colour neutral

  double E_fcn(int ifcn){return finalcns->at(ifcn)->getEnergy(); } ; 
                                 // Energy of final colour neutral ifcn.

  double M_fcn(int ifcn){return finalcns->at(ifcn)->getMass();} ;
                                 // Mass of final colour neutral ifcn.

  const double* p_fcn(int ifcn){return finalcns->at(ifcn)->getMomentum(); } ;
                                 // 3-momentum of final colour neutral ifcn.

  const double* p4_fcn(int ifcn) ;
                                 // 4-momentum  of final colour neutral ifcn (component 0 is E). 

  MCParticle* initial_elementon(int ijet) ;
                                 // The MCParticle that is at the beginning of the parton shower
                                 // leading to jet ijet

  MCParticle* final_elementon(int ijet) ;
                                 // The MCParticle that is at the end of the parton shower
                                 // leading to jet ijet

  const MCParticleVec& elementons_final_cn(int ifcn) ;
                                 // list of MCParticles constituting the final colour neutral

  const MCParticleVec& elementons_initial_cn(int ifcn) ;
                                 // list of MCParticles constituting the initial colour neutral

                                 // All the following navigators are set up by a call to getall ...
                                 // Note that most, if not all, information one can get using
                                 // these navigators can already be obtained in a simpler way using
                                 // the method above!

                                 // We list what they connect like so:
                                 // "a xxx -> yyys" means that navigator relzzz gives
                                 //
                                 // yyyvec = relzzz->getRelated T o Objects( a xxx ) , 
                                 //
                                 // and, consequently, that the relation in the other direction is gotten by
                                 //
                                 // xxxvec = relzzz->getRelated F r o m Objects( a yyy ) 
                                 //
                                 // In [] the LCIO type of the vector-elements are given. a truejet is a
                                 // ReconstructedParticle, and  jets->at(ijet) gives the truejet object
                                 // at index ijet (i.e. the identifier used in all the other methods),
                                 // while in the other way jetindex( truejet-object ) gives the index.


    LCRelationNavigator* relfcn{};      // a truejet to final colour neutral(s)    [ ReconstructedParticle:s ]
    LCRelationNavigator* relicn{};      // a truejet to initial colour neutral(s)   [ ReconstructedParticle:s ]

    LCRelationNavigator* relfp{};       // a truejet to its final quarks/leptons   [ MCParticle:s ]
    LCRelationNavigator* relip{};       // a truejet to its initial quarks/leptons [ MCParticle:s ]


    LCRelationNavigator* reltjreco{};   // a truejet to all seen particles in it    [ ReconstructedParticle:s ]
    LCRelationNavigator* reltjmcp{};    // a truejet to all true particles in it   [ MCParticle:s ]


                                 // Example:

//       ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( pfocol->getElementAt( j ) ) ;
//       LCObjectVec jetvec = reltjreco->getRelatedFromObjects( pfo );
//       int ijet = jetindex( dynamic_cast<ReconstructedParticle*>(jetvec[0]));

                                // will give you the jetnumber (ijet) that PFO j belongs to
                                // (Note that there is a method above for the opposite operation:
                                //  seen_partics(ijet) returns a ReconstructedParticlVec of all
                                //  PFOs belonging to jet ijet)

    LCRelationNavigator* reltrue_tj{};  // = RecoMCTruthLink

    LCCollection* tjcol{};
    LCCollection* fcncol{};
    LCCollection* icncol{};
    ReconstructedParticleVec* jets{};
    ReconstructedParticleVec* finalcns{};
    ReconstructedParticleVec* initialcns{};

 protected:
 /**  input collection names */
  std::string _trueJetCollectionName{};
  std::string _finalColourNeutralCollectionName{};
  std::string _initialColourNeutralCollectionName{};
  std::string _trueJetPFOLink{};
  std::string _trueJetMCParticleLink{};
  std::string _finalElementonLink{};
  std::string _initialElementonLink{};
  std::string _finalColourNeutralLink{};
  std::string _initialColourNeutralLink{};
  std::string _recoMCTruthLink{};
  int _COUNT_FSR{};
private:
  LCEvent* evt{};
  double p4[4]{};
  double p3[3]{}; 
  IntVec* intvec{};
  MCParticleVec* mcpartvec{};
} ;
#endif
