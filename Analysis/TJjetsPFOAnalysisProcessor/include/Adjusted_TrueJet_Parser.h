#ifndef Adjusted_TrueJet_Parser_h
#define Adjusted_TrueJet_Parser_h 1

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

struct MCPpyjet : LCIntExtension<MCPpyjet> {} ;
struct JetIndex : LCIntExtension<JetIndex> {} ;
struct IcnIndex : LCIntExtension<IcnIndex> {} ;
struct FcnIndex : LCIntExtension<FcnIndex> {} ;
//LCRelationNavigator* reltrue_tj =0;

class Adjusted_TrueJet_Parser {
  
 public:
  
  
  Adjusted_TrueJet_Parser()   ;
  ~Adjusted_TrueJet_Parser() ;
  
  
  virtual    std::string get_recoMCTruthLink(){ return _recoMCTruthLink  ;};
 
  ReconstructedParticleVec* getJets();
  ReconstructedParticleVec* getFinalcn();
  ReconstructedParticleVec* getInitialcn();

  void getall(LCEvent* event)  ;

  void delall() ;
  int njets() { return tjcol->getNumberOfElements(); };
  int jetindex(ReconstructedParticle* jet ) { return jet->ext<JetIndex>() ; };

  double Eseen(int ijet) { return jets->at(ijet)->getEnergy(); };
  double Mseen(int ijet) { return jets->at(ijet)->getMass(); };
  const double* pseen(int ijet) { return jets->at(ijet)->getMomentum(); };
  const double* p4seen(int ijet) ;
  const ReconstructedParticleVec& seen_partics(int ijet) { return jets->at(ijet)->getParticles(); };

  int type_jet(int ijet) { return jets->at(ijet)->getParticleIDs()[0]->getType() ; };

  double Etrue(int ijet) ; 
  double Mtrue(int ijet) ;
  const double* ptrue(int ijet) ;
  const double* p4true(int ijet) ;

  double Equark(int ijet) {return ( final_elementon(ijet) != NULL ? final_elementon(ijet)->getEnergy() : 0.);}; 
  double Mquark(int ijet) {return ( final_elementon(ijet) != NULL ? final_elementon(ijet)->getMass() : 0 );};
  const double* pquark(int ijet);
  const double* p4quark(int ijet) ;

  double Etrueseen(int ijet) ; 
  double Mtrueseen(int ijet) ;
  const double* ptrueseen(int ijet) ;
  const double* p4trueseen(int ijet) ;

  const IntVec& final_siblings( int ijet ); 
  const IntVec& initial_siblings( int ijet ); 

  int final_cn( int ijet ); 
  int initial_cn( int ijet ); 

  const IntVec& jets_of_final_cn( int ifcn ); 
  const IntVec& jets_of_initial_cn( int iicn ); 

  int nicn() { return icncol->getNumberOfElements(); };
  int type_icn_parent(int iicn) { return initialcns->at(iicn)->getParticleIDs()[0]->getType() ; };
  const IntVec& type_icn_comps(int iicn) ;
  int pdg_icn_parent(int iicn) { return initialcns->at(iicn)->getParticleIDs()[0]->getPDG() ; };
  const IntVec& pdg_icn_comps(int iicn) ;
  double E_icn(int iicn) { return initialcns->at(iicn)->getEnergy();}; 
  double M_icn(int iicn) {; return initialcns->at(iicn)->getMass();};
  const double* p_icn(int iicn){ return initialcns->at(iicn)->getMomentum();} ;
  const double* p4_icn(int iicn) ;

  int nfcn() { return fcncol->getNumberOfElements(); };
  int type_fcn_parent(int ifcn) { return finalcns->at(ifcn)->getParticleIDs()[0]->getType() ; };
  const IntVec& type_fcn_comps(int ifcn) ;
  int pdg_fcn_parent(int ifcn) { return finalcns->at(ifcn)->getParticleIDs()[0]->getPDG() ; };
  const IntVec& pdg_fcn_comps(int ifcn) ;
  double E_fcn(int ifcn){return finalcns->at(ifcn)->getEnergy(); } ; 
  double M_fcn(int ifcn){return finalcns->at(ifcn)->getMass();} ;
  const double* p_fcn(int ifcn){return finalcns->at(ifcn)->getMomentum(); } ;
  const double* p4_fcn(int ifcn) ;

  MCParticle* initial_elementon(int ijet) ;
  MCParticle* final_elementon(int ijet) ;

  const MCParticleVec& elementons_final_cn(int ifcn) ;
  const MCParticleVec& elementons_initial_cn(int ifcn) ;


    LCRelationNavigator* relfcn ;
    LCRelationNavigator* relicn ;
    LCRelationNavigator* relfp ;
    LCRelationNavigator* relip ;
    LCRelationNavigator* reltjreco ;
    LCRelationNavigator* reltjmcp ;
    LCRelationNavigator* reltrue_tj ;
    LCCollection* tjcol ;
    LCCollection* fcncol ;
    LCCollection* icncol ;
    ReconstructedParticleVec* jets;
    ReconstructedParticleVec* finalcns;
    ReconstructedParticleVec* initialcns;

 protected:
 /**  input collection names */
  std::string _trueJetCollectionName ;
  std::string _finalColourNeutralCollectionName;
  std::string _initialColourNeutralCollectionName;
  std::string _trueJetPFOLink ;
  std::string _trueJetMCParticleLink ;
  std::string _finalElementonLink ;
  std::string _initialElementonLink ;
  std::string _finalColourNeutralLink ;
  std::string _initialColourNeutralLink;
  std::string _recoMCTruthLink;
  int _COUNT_FSR;
private:
  LCEvent* evt ;
  double p4[4];
  double p3[3]; 
  IntVec* intvec;
  MCParticleVec* mcpartvec;
} ;
#endif


