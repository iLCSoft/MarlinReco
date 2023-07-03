#ifndef FindParticle_h_1
#define FindParticle_h_1


#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include "DDMarlinCED.h"
#include "TVector3.h"

typedef EVENT::ReconstructedParticle* pfo;

typedef std::vector<EVENT::ReconstructedParticle*> pfoVector;

typedef std::vector<EVENT::MCParticle*> mcpVector;

void getTrueDecayProducts( EVENT::MCParticle *parentMCP , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *trueNeutrino , mcpVector &trueNeutralDecayProducts , mcpVector &trueChargedDecayProducts );

void getTruePVADecayProducts( EVENT::MCParticle *parentMCP , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *trueNeutrino , pfo &linkedRecoLepton , float &weightRecoLeptoMCLep , float &weightMCLeptoRecoLep , pfoVector &truePVANeutralDecayProducts , pfoVector &truePVAChargedDecayProducts , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav );

void getLinkedRecoLepton( EVENT::MCParticle *mcParticle , pfo &linkedRecoLepton , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , float &weightRecoLeptoMCLep , float &weightMCLeptoRecoLep );

EVENT::ReconstructedParticle* getLinkedPFO( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedPFO , bool getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO , bool &foundPFO );

#endif
