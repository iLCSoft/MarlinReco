#ifndef FindParticle_h_1
#define FindParticle_h_1


#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include "DDMarlinCED.h"
#include "TVector3.h"
#include "SLDCorrectionTypes.h"

void getTrueDecayProducts( const MCP &parentMCP , const MCP &SLDLepton , const MCP &trueNeutrino , MCPVector &trueNeutralDecayProducts , MCPVector &trueChargedDecayProducts );

void getTruePVADecayProducts( const MCP &parentMCP , const MCP &SLDLepton , const MCP &trueNeutrino , RecoParticle &linkedRecoLepton , float &weightRecoLeptoMCLep , float &weightMCLeptoRecoLep , PFOVector &truePVANeutralDecayProducts , PFOVector &truePVAChargedDecayProducts , const LCRelationNavigator &RecoMCParticleNav , const LCRelationNavigator &MCParticleRecoNav );

void getLinkedRecoLepton( const MCP &mcParticle , RecoParticle &linkedRecoLepton , const LCRelationNavigator &RecoMCParticleNav , const LCRelationNavigator &MCParticleRecoNav , float &weightRecoLeptoMCLep , float &weightMCLeptoRecoLep );

RecoParticle getLinkedPFO( const MCP &mcParticle , const LCRelationNavigator &RecoMCParticleNav , const LCRelationNavigator &MCParticleRecoNav , const bool &getChargedPFO , const bool &getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO , bool &foundPFO );

#endif
