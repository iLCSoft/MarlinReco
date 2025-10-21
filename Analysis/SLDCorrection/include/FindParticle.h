#ifndef FindParticle_h_1
#define FindParticle_h_1

#include "DDMarlinCED.h"
#include "SLDCorrectionTypes.h"
#include "TVector3.h"
#include "UTIL/LCRelationNavigator.h"
#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>

//	getTrueDecayProducts finds true charged (trueChargedDecayProducts) and neutral (trueNeutralDecayProducts) decay
// products of a hadron (parentMCP) that decays semi-leptonically to a charged lepton (SLDLepton)
void getTrueDecayProducts(const MCP& parentMCP, const MCP& SLDLepton, const MCP& trueNeutrino,
                          MCPVector& trueNeutralDecayProducts, MCPVector& trueChargedDecayProducts);

//	getTruePVADecayProducts follows the decay chain of a semi-leptonic decay of a hadron (parentMCP) from
// MCParticles that decays to a charged lepton (true:SLDLepton, linkes reconstructed particle: linkedRecoLepton with
// weights MCParticle <--> ReconstructedParticle: weightMCLeptoRecoLep & weightRecoLeptoMCLep) and corresponding
// neutrino (trueNeutrino) 	getTruePVADecayProducts finds the reconstructed charged (truePVAChargedDecayProducts)
// and neutral (truePVANeutralDecayProducts) decay products of the semi-leptonic decay using the LCRelationNavigators
// from MCParticle to ReconstrucedParticles and vice-versa
void getTruePVADecayProducts(const MCP& parentMCP, const MCP& SLDLepton, const MCP& trueNeutrino,
                             RecoParticle& linkedRecoLepton, float& weightRecoLeptoMCLep, float& weightMCLeptoRecoLep,
                             PFOVector& truePVANeutralDecayProducts, PFOVector& truePVAChargedDecayProducts,
                             const LCRelationNavigator& RecoMCParticleNav,
                             const LCRelationNavigator& MCParticleRecoNav);

//	getLinkedRecoLepton finds the reconstructed particle (linkedRecoLepton) that is linked to the charged lepton
//(SLDLepton) from the semi-leptonic decay using LCRelationNavigators (MCParticleRecoNav & RecoMCParticleNav) and finds
// the weight of links in both direction SLDLepton <--> linkedRecoLepton: weightMCLeptoRecoLep & weightRecoLeptoMCLep
void getLinkedRecoLepton(const MCP& SLDLepton, RecoParticle& linkedRecoLepton,
                         const LCRelationNavigator& RecoMCParticleNav, const LCRelationNavigator& MCParticleRecoNav,
                         float& weightRecoLeptoMCLep, float& weightMCLeptoRecoLep);

//	getLinkedPFO returns a reconstructed particle that is linked to mcParticle with highest weight using
// LCRelationNavigators in both directions (MCParticleRecoNav & RecoMCParticleNav) 	if getChargedPFO = true: looks
// for the reconstructed particle with highest weight of its track to mcParticle 	if getNeutralPFO = true: looks
// for the reconstructed particle with highest weight of its cluster to mcParticle 	weightMCPtoPFO & weightPFOtoMCP
// are the highest weight of returned Reconstructed Particle to mcParticle 	if a Reconstructed Particle is found,
// foundPFO is set to true, else is set to false
RecoParticle getLinkedPFO(const MCP& mcParticle, const LCRelationNavigator& RecoMCParticleNav,
                          const LCRelationNavigator& MCParticleRecoNav, const bool& getChargedPFO,
                          const bool& getNeutralPFO, float& weightPFOtoMCP, float& weightMCPtoPFO, bool& foundPFO);

#endif
