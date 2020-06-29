#ifndef UTILITIES_H
#define UTILITIES_H

// *******************************************************
// some useful functions
// *******************************************************

#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include "TROOT.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace lcio ;
using namespace std;

namespace isolep{

  Int_t getMCSerial(MCParticle *mcPart, LCCollection *colMCP);
  //MCParticle* getLinkedMCParticle(ReconstructedParticle *recPart, LCCollection *colMCTL, Double_t &weight, Int_t &nMCTL);
  MCParticle *getMCParticle(ReconstructedParticle *recPart, LCCollection *colMCTL);
  MCParticle *getMCParticle(ReconstructedParticle *recPart, LCCollection *colMCTL, Double_t &weight);
  MCParticle *getMCParticle(ReconstructedParticle *recPart, LCCollection *colMCTL, Double_t &weight, Int_t &nMCTL);
  Int_t getLinkedMCParticle(ReconstructedParticle *recPart, LCCollection *colMCTL, Double_t &weight, Int_t &nMCTL);
  Int_t getOriginalPDG(MCParticle *mcPart, Bool_t iHiggs = 0);
  Int_t getOriginalPDGForIsoLep(MCParticle *mcPart, LCCollection *colMC);
  Int_t getOriginalPDGForIsoLep(MCParticle *mcPart);
  Int_t getOriginalSerial(MCParticle *mcPart, LCCollection *colMCP, Bool_t iHiggs = 0);
  Int_t getOriginalSerial(ReconstructedParticle *recPart, LCCollection *colMCTL, LCCollection *colMCP, Bool_t iHiggs = 0);
  Int_t getOriginalSerialForZHH(MCParticle *mcPart, LCCollection *colMCP);
  Double_t getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone);
  Double_t getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, Int_t mode);
  Double_t getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, 
			 std::vector<lcio::ReconstructedParticle*> &conePFOs);
  Double_t getInvariantMass(ReconstructedParticle *recPart1, ReconstructedParticle *recPart2);
  Double_t getInvariantMass(ReconstructedParticle *recPart1, ReconstructedParticle *recPart2, ReconstructedParticle *recPart3);
  
  Int_t getLeptonID(ReconstructedParticle *recPart);
  Bool_t getFSRTag(ReconstructedParticle *motherPart, ReconstructedParticle *recPart, Double_t fCosFSRCut = 0.999);
  Bool_t getSplitTag(ReconstructedParticle *motherPart, ReconstructedParticle *recPart); 
//TVector3 getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, Bool_t woFSR);
//void getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, Bool_t woFSR, TVector3 coneEnergy0);
  void getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, Bool_t woFSR, Double_t coneEnergy[3], Double_t pFSR[4]);
  void getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, Bool_t woFSR, Double_t coneEnergy[3], Double_t pFSR[4], 
		     Double_t conCone2, Double_t pCone2[4], Int_t &nConePhoton);
  TLorentzVector getFSRMomentum(ReconstructedParticle *recPart, LCCollection *colPFO);
  Int_t isSelectedByFastJet( ReconstructedParticle *pfo, LCCollection *colFastJet, Double_t &ratioEPartEJet, Double_t &ratioPTMJet);
  void doPhotonRecovery(ReconstructedParticle *electron, LCCollection *colPFO, ReconstructedParticleImpl *recoElectron, Double_t fCosFSRCut = 0.999);   
  void doPhotonRecovery(ReconstructedParticle *electron, LCCollection *colPFO, ReconstructedParticleImpl *recoElectron, Double_t fCosFSRCut, Int_t lepType, std::vector<lcio::ReconstructedParticle*> &photons);
  Bool_t isFoundInVector(ReconstructedParticle *pfo, std::vector<lcio::ReconstructedParticle*> &pfos);
  
  Double_t jetFunction( TLorentzVector lortz, float beta, float power);
  Double_t jetFunction( TLorentzVector lortz, float virtual_scale=20); // GeV
  
  Int_t calculateEnergyComponents(LCCollection *colMC, LCCollection *colMCTL, ReconstructedParticle *jet, Double_t energyComponents[4]);
  Int_t calculateEnergyComponents(LCCollection *colMC, LCCollection *colMCTL, ReconstructedParticle *jet, Double_t energyComponents[4],
				  Int_t iColor1, Int_t iColor2, Int_t iColor3);

  Double_t getLikelihood(TString fname, TString hist, Double_t mass);
  
  void listMCParticles(LCCollection *colMC);
  
  Bool_t isOverlay(ReconstructedParticle *pfo, LCCollection *colMCTL);
  
  void dumpJetParticles(ReconstructedParticle *jet, LCCollection *colMC, LCCollection *colMCTL);
  
  Int_t getVertexComponents(LCCollection *colMC, LCCollection *colMCTL, ReconstructedParticle *vetex, Double_t energyComponents[2], Int_t nparticles[2]);

  Int_t isVertexFromOverlay( ReconstructedParticle *vertex, LCCollection *colMC, LCCollection *colMCTL);
  
  Double_t getJetDistance( ReconstructedParticle *i, ReconstructedParticle *j, TString algorithm, Double_t R);
  
  Int_t isFromVertex( ReconstructedParticle *recPart, LCCollection *colVertex);

  std::vector<Int_t> getHiggsDecayModes(LCCollection *colMC);
  std::vector<Int_t> getNumberOfOverlayEvents(Double_t fEcm, LCCollection *colMC);

  void mcDebug(LCCollection *colMC);
  TLorentzVector getLorentzEcm(Double_t fEcm);
  TLorentzVector getLorentzEcm(Double_t fEcm,Bool_t isCrossingAngle);
  Double_t getRecoilMass(TLorentzVector lortzEcm, TLorentzVector lortzZ);
  Double_t getAcoPlanarity(TLorentzVector lortz1, TLorentzVector lortz2);
  TLorentzVector getLorentzVector(ReconstructedParticle *pfo);
  Double_t getCosTheta(ReconstructedParticle *part1, ReconstructedParticle *part2);
  Double_t getCosTheta(TLorentzVector part1, TLorentzVector part2);

}

#endif
