#include "FindParticle.h"

void getTrueDecayProducts(const MCP& parentMCP, const MCP& SLDLepton, const MCP& trueNeutrino,
                          MCPVector& trueNeutralDecayProducts, MCPVector& trueChargedDecayProducts) {
  for (long unsigned int i_daughter = 0; i_daughter < (parentMCP->getDaughters()).size(); ++i_daughter) {
    MCP daughter = parentMCP->getDaughters()[i_daughter];
    if (daughter == SLDLepton || daughter == trueNeutrino)
      continue;
    bool trueParticleExist = false;
    if (daughter->getGeneratorStatus() == 1) {
      if (fabs(daughter->getCharge()) < 0.1) {
        for (unsigned int i_par = 0; i_par < trueNeutralDecayProducts.size(); ++i_par) {
          if (daughter == trueNeutralDecayProducts[i_par])
            trueParticleExist = true;
        }
        if (!trueParticleExist)
          trueNeutralDecayProducts.push_back(daughter);
      } else {
        for (unsigned int i_par = 0; i_par < trueChargedDecayProducts.size(); ++i_par) {
          if (daughter == trueChargedDecayProducts[i_par])
            trueParticleExist = true;
        }
        if (!trueParticleExist)
          trueChargedDecayProducts.push_back(daughter);
      }
    } else {
      getTrueDecayProducts(daughter, SLDLepton, trueNeutrino, trueNeutralDecayProducts, trueChargedDecayProducts);
    }
  }
}

void getTruePVADecayProducts(const MCP& parentMCP, const MCP& SLDLepton, const MCP& trueNeutrino,
                             RecoParticle& linkedRecoLepton, float& weightRecoLeptoMCLep, float& weightMCLeptoRecoLep,
                             PFOVector& truePVANeutralDecayProducts, PFOVector& truePVAChargedDecayProducts,
                             const LCRelationNavigator& RecoMCParticleNav,
                             const LCRelationNavigator& MCParticleRecoNav) {
  for (long unsigned int i_daughter = 0; i_daughter < (parentMCP->getDaughters()).size(); ++i_daughter) {
    MCP daughter = parentMCP->getDaughters()[i_daughter];
    ReconstructedParticle* linkedChargedPFO{};
    ReconstructedParticle* linkedNeutralPFO{};
    bool foundLinkedChargedPFO = false;
    bool foundLinkedNeutralPFO = false;
    float weightChargedPFOtoMCP = 0.0;
    float weightChargedMCPtoPFO = 0.0;
    float weightNeutralPFOtoMCP = 0.0;
    float weightNeutralMCPtoPFO = 0.0;
    bool recoParticleExist = false;
    if (daughter == trueNeutrino) {
      continue;
    } else if (daughter == SLDLepton) {
      getLinkedRecoLepton(SLDLepton, linkedRecoLepton, RecoMCParticleNav, MCParticleRecoNav, weightRecoLeptoMCLep,
                          weightMCLeptoRecoLep);
      getTruePVADecayProducts(SLDLepton, SLDLepton, trueNeutrino, linkedRecoLepton, weightNeutralPFOtoMCP,
                              weightNeutralMCPtoPFO, truePVANeutralDecayProducts, truePVAChargedDecayProducts,
                              RecoMCParticleNav, MCParticleRecoNav);
    } else {
      linkedChargedPFO = getLinkedPFO(daughter, RecoMCParticleNav, MCParticleRecoNav, true, false,
                                      weightChargedPFOtoMCP, weightChargedMCPtoPFO, foundLinkedChargedPFO);
      linkedNeutralPFO = getLinkedPFO(daughter, RecoMCParticleNav, MCParticleRecoNav, false, true,
                                      weightNeutralPFOtoMCP, weightNeutralMCPtoPFO, foundLinkedNeutralPFO);
    }
    if (foundLinkedChargedPFO) {
      for (unsigned int i_par = 0; i_par < truePVAChargedDecayProducts.size(); ++i_par) {
        if (linkedChargedPFO == truePVAChargedDecayProducts[i_par])
          recoParticleExist = true;
      }
      for (unsigned int i_par = 0; i_par < truePVANeutralDecayProducts.size(); ++i_par) {
        if (linkedChargedPFO == truePVANeutralDecayProducts[i_par])
          recoParticleExist = true;
      }
      if (!recoParticleExist && linkedChargedPFO != linkedRecoLepton)
        truePVAChargedDecayProducts.push_back(linkedChargedPFO);
    } else if (foundLinkedNeutralPFO) {
      for (unsigned int i_par = 0; i_par < truePVAChargedDecayProducts.size(); ++i_par) {
        if (linkedNeutralPFO == truePVAChargedDecayProducts[i_par])
          recoParticleExist = true;
      }
      for (unsigned int i_par = 0; i_par < truePVANeutralDecayProducts.size(); ++i_par) {
        if (linkedNeutralPFO == truePVANeutralDecayProducts[i_par])
          recoParticleExist = true;
      }
      if (!recoParticleExist && linkedNeutralPFO != linkedRecoLepton) {
        if ((linkedNeutralPFO->getTracks()).size() == 0) {
          truePVANeutralDecayProducts.push_back(linkedNeutralPFO);
        } else {
          truePVAChargedDecayProducts.push_back(linkedNeutralPFO);
        }
      }
    } else {
      getTruePVADecayProducts(daughter, SLDLepton, trueNeutrino, linkedRecoLepton, weightChargedPFOtoMCP,
                              weightChargedMCPtoPFO, truePVANeutralDecayProducts, truePVAChargedDecayProducts,
                              RecoMCParticleNav, MCParticleRecoNav);
    }
  }
}

void getLinkedRecoLepton(const MCP& SLDLepton, RecoParticle& linkedRecoLepton,
                         const LCRelationNavigator& RecoMCParticleNav, const LCRelationNavigator& MCParticleRecoNav,
                         float& weightRecoLeptoMCLep, float& weightMCLeptoRecoLep) {
  ReconstructedParticle* linkedRecoLep{};
  bool foundLinkedRecoLepton = false;
  linkedRecoLep = getLinkedPFO(SLDLepton, RecoMCParticleNav, MCParticleRecoNav, true, false, weightRecoLeptoMCLep,
                               weightMCLeptoRecoLep, foundLinkedRecoLepton);
  if (!foundLinkedRecoLepton) {
    for (long unsigned int i_daughter = 0; i_daughter < (SLDLepton->getDaughters()).size(); ++i_daughter) {
      MCP daughter = SLDLepton->getDaughters()[i_daughter];
      getLinkedRecoLepton(daughter, linkedRecoLep, RecoMCParticleNav, MCParticleRecoNav, weightRecoLeptoMCLep,
                          weightMCLeptoRecoLep);
    }
  }
  linkedRecoLepton = linkedRecoLep;
}

RecoParticle getLinkedPFO(const MCP& mcParticle, const LCRelationNavigator& RecoMCParticleNav,
                          const LCRelationNavigator& MCParticleRecoNav, const bool& getChargedPFO,
                          const bool& getNeutralPFO, float& weightPFOtoMCP, float& weightMCPtoPFO,
                          bool& foundlinkedPFO) {
  streamlog_out(DEBUG1) << "" << std::endl;
  streamlog_out(DEBUG1) << "	Look for PFO linked to visible MCParticle:" << std::endl;
  streamlog_out(DEBUG1) << *mcParticle << std::endl;
  ReconstructedParticle* linkedPFO{};
  foundlinkedPFO = false;
  const EVENT::LCObjectVec& PFOvec = MCParticleRecoNav.getRelatedToObjects(mcParticle);
  const EVENT::FloatVec& PFOweightvec = MCParticleRecoNav.getRelatedToWeights(mcParticle);
  streamlog_out(DEBUG0) << "	Visible MCParticle is linked to " << PFOvec.size() << " PFO(s)" << std::endl;
  weightPFOtoMCP = 0.0;
  weightMCPtoPFO = 0.0;
  double maxweightPFOtoMCP = 0.;
  double maxweightMCPtoPFO = 0.;
  int iPFOtoMCPmax = -1;
  int iMCPtoPFOmax = -1;
  for (unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++) {
    double pfo_weight = 0.0;
    double trackWeight = (int(PFOweightvec.at(i_pfo)) % 10000) / 1000.0;
    double clusterWeight = (int(PFOweightvec.at(i_pfo)) / 10000) / 1000.0;
    if (getChargedPFO && !getNeutralPFO) {
      pfo_weight = trackWeight;
    } else if (getNeutralPFO && !getChargedPFO) {
      pfo_weight = clusterWeight;
    } else {
      pfo_weight = (trackWeight > clusterWeight ? trackWeight : clusterWeight);
    }
    streamlog_out(DEBUG0) << "	Visible MCParticle linkWeight to PFO: " << PFOweightvec.at(i_pfo)
                          << " (Track: " << trackWeight << " , Cluster: " << clusterWeight << ")" << std::endl;
    ReconstructedParticle* testPFO = (ReconstructedParticle*)PFOvec.at(i_pfo);
    if (pfo_weight > maxweightMCPtoPFO) //&& track_weight >= m_MinWeightTrackMCTruthLink )
    {
      maxweightMCPtoPFO = pfo_weight;
      iMCPtoPFOmax = i_pfo;
      streamlog_out(DEBUG0) << "	PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType()
                            << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
    }
  }
  if (getChargedPFO && maxweightMCPtoPFO < 0.5) {
    streamlog_out(DEBUG1) << "	MCParticle has link weight lower than 0.8 ( " << maxweightMCPtoPFO
                          << " ), looking for linked PFO in clusters" << std::endl;
    maxweightMCPtoPFO = 0.0;
    for (unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++) {
      double pfo_weight = (int(PFOweightvec.at(i_pfo)) / 10000) / 1000.0;
      streamlog_out(DEBUG0) << "	Visible MCParticle linkWeight to PFO: " << PFOweightvec.at(i_pfo)
                            << " (Track: " << (int(PFOweightvec.at(i_pfo)) % 10000) / 1000.0
                            << " , Cluster: " << (int(PFOweightvec.at(i_pfo)) / 10000) / 1000.0 << ")" << std::endl;
      ReconstructedParticle* testPFO = (ReconstructedParticle*)PFOvec.at(i_pfo);
      if (pfo_weight > maxweightMCPtoPFO) //&& track_weight >= m_MinWeightTrackMCTruthLink )
      {
        maxweightMCPtoPFO = pfo_weight;
        iMCPtoPFOmax = i_pfo;
        streamlog_out(DEBUG0) << "	PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType()
                              << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
      }
    }
  }
  if (iMCPtoPFOmax != -1) {
    ReconstructedParticle* testPFO = (ReconstructedParticle*)PFOvec.at(iMCPtoPFOmax);
    const EVENT::LCObjectVec& MCPvec = RecoMCParticleNav.getRelatedToObjects(testPFO);
    const EVENT::FloatVec& MCPweightvec = RecoMCParticleNav.getRelatedToWeights(testPFO);
    for (unsigned int i_mcp = 0; i_mcp < MCPvec.size(); i_mcp++) {
      double mcp_weight = 0.0;
      double trackWeight = (int(MCPweightvec.at(i_mcp)) % 10000) / 1000.0;
      double clusterWeight = (int(MCPweightvec.at(i_mcp)) / 10000) / 1000.0;
      if (getChargedPFO && !getNeutralPFO) {
        mcp_weight = trackWeight;
      } else if (getNeutralPFO && !getChargedPFO) {
        mcp_weight = clusterWeight;
      } else {
        mcp_weight = (trackWeight > clusterWeight ? trackWeight : clusterWeight);
      }
      MCParticle* testMCP = (MCParticle*)MCPvec.at(i_mcp);
      if (mcp_weight > maxweightPFOtoMCP) //&& mcp_weight >= m_MinWeightTrackMCTruthLink )
      {
        maxweightPFOtoMCP = mcp_weight;
        iPFOtoMCPmax = i_mcp;
        streamlog_out(DEBUG0) << "	MCParticle at index: " << testMCP->id() << " has PDG: " << testMCP->getPDG()
                              << " and PFO to MCParticle link weight is " << mcp_weight << std::endl;
      }
    }
    if (iPFOtoMCPmax != -1) {
      if (MCPvec.at(iPFOtoMCPmax) == mcParticle) {
        linkedPFO = testPFO;
        foundlinkedPFO = true;
      }
    }
  }

  if (foundlinkedPFO) {
    streamlog_out(DEBUG1) << "	Linked PFO to MCParticle found successfully " << std::endl;
    streamlog_out(DEBUG1) << *linkedPFO << std::endl;
    weightPFOtoMCP = maxweightPFOtoMCP;
    weightMCPtoPFO = maxweightMCPtoPFO;
    return linkedPFO;
  } else {
    streamlog_out(DEBUG1) << "	Couldn't Find a PFO linked to MCParticle" << std::endl;
    return NULL;
  }
}
