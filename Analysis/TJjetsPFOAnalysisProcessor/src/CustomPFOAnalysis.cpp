#include "TJjetsPFOAnalysisProcessor.h"

void TJjetsPFOAnalysisProcessor::Clear()
{
    m_pfoVector.clear();
    m_pfoTargetVector.clear();

    m_nPfosTotal = 0;
    m_nPfosNeutralHadrons = 0;
    m_nPfosPhotons = 0;
    m_nPfosTracks = 0;
    m_pfoEnergyTotal = 0.f;
    m_pfoEnergyNeutralHadrons = 0.f;
    m_pfoEnergyPhotons = 0.f;

    m_pfoEnergyTracks = 0.f;
    m_pfoECalToEmEnergy = 0.f;
    m_pfoECalToHadEnergy = 0.f;
    m_pfoHCalToEmEnergy = 0.f;
    m_pfoHCalToHadEnergy = 0.f;
    m_pfoMuonToEnergy = 0.f;
    m_pfoOtherEnergy = 0.f;
    m_pfoMassTotal = 0.f;

    m_pfoEnergies.clear();
    m_pfoPx.clear();
    m_pfoPy.clear();
    m_pfoPz.clear();
    m_pfoCosTheta.clear();

    m_pfoTargetEnergies.clear();
    m_pfoTargetPx.clear();
    m_pfoTargetPy.clear();
    m_pfoTargetPz.clear();
    m_pfoTargetCosTheta.clear();

    m_pfoPdgCodes.clear();
    m_pfoTargetPdgCodes.clear();

    m_nPfoTargetsTotal = 0;
    m_nPfoTargetsNeutralHadrons = 0;
    m_nPfoTargetsPhotons = 0;
    m_nPfoTargetsTracks = 0;

    m_pfoTargetsEnergyTotal = 0.f;
    m_pfoTargetsEnergyNeutralHadrons = 0.f;
    m_pfoTargetsEnergyPhotons = 0.f;
    m_pfoTargetsEnergyTracks = 0.f;

    m_mcEnergyENu = 0.f;
    m_mcEnergyFwd = 0.f;
    m_eQQ = -99.f;
    m_eQ1 = -99.f;
    m_eQ2 = -99.f;
    m_costQQ = -99.f;
    m_costQ1 = -99.f;
    m_costQ2 = -99.f;
    m_mQQ = -99.f;
    m_thrust = -99.f;
    m_qPdg = -99;

}

//------------------------------------------------------------------------------------------------------------------------------------------

void TJjetsPFOAnalysisProcessor::ExtractCollections(JetContentPair* jet_content)
{
    // Extract reconstructed pfo collection
    try
    {

        ReconstructedParticleVec recos = jet_content->second;

        for (unsigned int i = 0, nElements = recos.size(); i < nElements; ++i)
        {
            const EVENT::ReconstructedParticle *pReconstructedParticle = recos[i];

            if (NULL == pReconstructedParticle)
                throw EVENT::Exception("Collection type mismatch");

            m_pfoVector.push_back(pReconstructedParticle);
        }
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Could not extract input particle collection: " << m_inputPfoCollection << std::endl;
    }

    // Extract mc pfo collection
    MCParticleList preliminary_mcPfoCandidates, mcPfoCandidates;

    try
    {

        MCParticleVec mcs = jet_content->first;

        for (unsigned int i = 0, nElements = mcs.size(); i < nElements; ++i)
        {
            MCParticle* pMCParticle = mcs[i];

            if (pMCParticle == NULL) {
              throw EVENT::Exception("Collection type mismatch");
            }

            // if (!pMCParticle->getParents().empty()) {
            //   continue;
            // }

            this->ApplyPfoSelectionRules(pMCParticle, preliminary_mcPfoCandidates);
        }
        streamlog_out(DEBUG) << "N leftover mcs: " << preliminary_mcPfoCandidates.size() << " N input: " << mcs.size() <<  std::endl;


        std::set<EVENT::MCParticle*>::iterator it;
        for (it = preliminary_mcPfoCandidates.begin(); it != preliminary_mcPfoCandidates.end(); it++) {
            if ( hasSomeParentsInMCList( *it, preliminary_mcPfoCandidates ) ) {
              continue;
            }
            mcPfoCandidates.insert(*it);
        }
    }
    catch (...)
    {
        streamlog_out(WARNING) << "Could not extract mc particle collection " << m_mcParticleCollection << std::endl;
    }

    m_pfoTargetVector.insert(m_pfoTargetVector.begin(), mcPfoCandidates.begin(), mcPfoCandidates.end());
    std::sort(m_pfoTargetVector.begin(), m_pfoTargetVector.end(), TJjetsPFOAnalysisProcessor::SortPfoTargetsByEnergy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TJjetsPFOAnalysisProcessor::hasSomeParentsInMCList(EVENT::MCParticle *pMCParticle, MCParticleList &mcs) const {
  // Check if the MCparticle has some parents/grandparents/grandgrandparents/... in the given MCParticleList
  bool found_parent = false;
  MCParticleVec pMCParticleParents (pMCParticle->getParents());
  MCParticleVec mcs_vector (mcs.begin(), mcs.end());

  if ( ! areDisjointVectors( pMCParticleParents, mcs_vector ) ) {
    found_parent = true;
  } else { // Check next higher generation
    for (unsigned int i=0; i<pMCParticleParents.size(); i++) {
      if ( hasSomeParentsInMCList(pMCParticleParents[i], mcs) ) {
        found_parent = true;
      }
    }
  }

  return found_parent;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void TJjetsPFOAnalysisProcessor::ApplyPfoSelectionRules(EVENT::MCParticle *pMCParticle, MCParticleList &mcPfoCandidates) const
{
    const float innerRadius(std::sqrt(pMCParticle->getVertex()[0] * pMCParticle->getVertex()[0] + pMCParticle->getVertex()[1] * pMCParticle->getVertex()[1] + pMCParticle->getVertex()[2] * pMCParticle->getVertex()[2]));
    const float outerRadius(std::sqrt(pMCParticle->getEndpoint()[0] * pMCParticle->getEndpoint()[0] + pMCParticle->getEndpoint()[1] * pMCParticle->getEndpoint()[1] + pMCParticle->getEndpoint()[2] * pMCParticle->getEndpoint()[2]));
    const float momentum(std::sqrt(pMCParticle->getMomentum()[0] * pMCParticle->getMomentum()[0] + pMCParticle->getMomentum()[1] * pMCParticle->getMomentum()[1] + pMCParticle->getMomentum()[2] * pMCParticle->getMomentum()[2]));

    if ((mcPfoCandidates.find(pMCParticle) == mcPfoCandidates.end()) &&
        (outerRadius > m_mcPfoSelectionRadius) && (innerRadius <= m_mcPfoSelectionRadius) && (momentum > m_mcPfoSelectionMomentum) &&
        !((pMCParticle->getPDG() == 2212 || pMCParticle->getPDG() == 2112) && (pMCParticle->getEnergy() < m_mcPfoSelectionLowEnergyNPCutOff)))
    {
        mcPfoCandidates.insert(pMCParticle);
    }
    // else
    // {
    //     const EVENT::MCParticleVec &daughterVector(pMCParticle->getDaughters());
    //
    //     for (EVENT::MCParticleVec::const_iterator iter = daughterVector.begin(), iterEnd = daughterVector.end(); iter != iterEnd; ++iter)
    //     {
    //         this->ApplyPfoSelectionRules(*iter, mcPfoCandidates);
    //     }
    // }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TJjetsPFOAnalysisProcessor::MakeQuarkVariables(JetContentPair* jet_content)
{
    MCParticleVector mcQuarkVector;

///////////////////////////////////////////////////////////////////
    try
    {

        MCParticleVec mcs = jet_content->first;
        MCParticleList mcs_list =  MCParticleList(mcs.begin(),mcs.end());

        for (unsigned int i = 0, nElements = mcs.size(); i < nElements; ++i)
        {
            EVENT::MCParticle *pMCParticle = mcs[i];

            if (NULL == pMCParticle)
                throw EVENT::Exception("Collection type mismatch");

            const int absPdgCode(std::abs(pMCParticle->getPDG()));

            // By default, the primary quarks are the ones without parents included in the jet-mc
            if (!m_lookForQuarksWithMotherZ)
            {
                if ((absPdgCode >= 1) && (absPdgCode <= 6) && !hasSomeParentsInMCList(pMCParticle, mcs_list) )// pMCParticle->getParents().empty()) //!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    mcQuarkVector.push_back(pMCParticle);
            }
            else
            {
                // For MC files generated in the SLIC environment, the primary quarks have parents; the mother should be the Z-boson
                if ((absPdgCode >= 1) && (absPdgCode <= 6))
                {
                    if ((pMCParticle->getParents().size() == 1) && ((pMCParticle->getParents())[0]->getPDG() == 23))
                        mcQuarkVector.push_back(pMCParticle);
                }
            }
        }
    }
    catch (...)
    {
        streamlog_out(WARNING) << "Could not extract mc quark information" << std::endl;
    }
///////////////////////////////////////////////////////////////////

    if (!mcQuarkVector.empty())
    {
        m_qPdg = std::abs(mcQuarkVector[0]->getPDG());
        float energyTot(0.f);
        float costTot(0.f);

        for (unsigned int i = 0; i < mcQuarkVector.size(); ++i)
        {
            const float px(mcQuarkVector[i]->getMomentum()[0]);
            const float py(mcQuarkVector[i]->getMomentum()[1]);
            const float pz(mcQuarkVector[i]->getMomentum()[2]);
            const float energy(mcQuarkVector[i]->getEnergy());
            const float p(std::sqrt(px * px + py * py + pz * pz));
            const float cost(std::fabs(pz) / p);
            energyTot += energy;
            costTot += cost * energy;
        }

        m_thrust = costTot / energyTot;
    }

    if (mcQuarkVector.size() == 2)
    {
        const float pQ1x = mcQuarkVector[0]->getMomentum()[0];
        const float pQ1y = mcQuarkVector[0]->getMomentum()[1];
        const float pQ1z = mcQuarkVector[0]->getMomentum()[2];

        const float pQ2x = mcQuarkVector[1]->getMomentum()[0];
        const float pQ2y = mcQuarkVector[1]->getMomentum()[1];
        const float pQ2z = mcQuarkVector[1]->getMomentum()[2];

        const float pQ1[3] = {pQ1x, pQ1y, pQ1z};
        const float pQ2[3] = {pQ2x, pQ2y, pQ2z};
        const float pQQ[3] = {pQ1[0] + pQ2[0], pQ1[1] + pQ2[1], pQ1[2] + pQ2[2]};

        const TLorentzVector q1(pQ1[0], pQ1[1], pQ1[2], mcQuarkVector[0]->getEnergy());
        const TLorentzVector q2(pQ2[0], pQ2[1], pQ2[2], mcQuarkVector[1]->getEnergy());
        const TLorentzVector qq = q1 + q2;

        m_mQQ = qq.M();
        m_eQQ = mcQuarkVector[0]->getEnergy() + mcQuarkVector[1]->getEnergy();

        const float pQ1Tot(std::sqrt(pQ1[0] * pQ1[0] + pQ1[1] * pQ1[1] + pQ1[2] * pQ1[2]));
        const float pQ2Tot(std::sqrt(pQ2[0] * pQ2[0] + pQ2[1] * pQ2[1] + pQ2[2] * pQ2[2]));
        const float pQQTot(std::sqrt(pQQ[0] * pQQ[0] + pQQ[1] * pQQ[1] + pQQ[2] * pQQ[2]));

        if (std::fabs(pQQTot) > std::numeric_limits<float>::epsilon())
            m_costQQ = pQQ[2] / pQQTot;

        if (std::fabs(pQ1Tot) > std::numeric_limits<float>::epsilon())
            m_costQ1 = pQ1[2] / pQ1Tot;

        if (std::fabs(pQ2Tot) > std::numeric_limits<float>::epsilon())
            m_costQ2 = pQ2[2] / pQ2Tot;

        m_eQ1 = pQ1Tot;
        m_eQ2 = pQ2Tot;

        if (m_printing > 0)
        {
            std::cout << " eQQ    = " << m_eQQ << std::endl
                      << " eQ1    = " << m_eQ1 << std::endl
                      << " eQ2    = " << m_eQ2 << std::endl
                      << " costQQ = " << m_costQQ << std::endl
                      << " costQ1 = " << m_costQ1 << std::endl
                      << " costQ2 = " << m_costQ2 << std::endl
                      << " mQQ    = " << m_mQQ << std::endl
                      << " Thrust = " << m_thrust << std::endl
                      << " QPDG   = " << m_qPdg << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TJjetsPFOAnalysisProcessor::PerformPfoAnalysis()
{
    float momTot[3] = {0.f, 0.f, 0.f};

    // Extract quantities relating to reconstructed pfos
    for (ParticleVector::const_iterator iter = m_pfoVector.begin(), iterEnd = m_pfoVector.end(); iter != iterEnd; ++iter)
    {
        const EVENT::ReconstructedParticle *pPfo = *iter;

        ++m_nPfosTotal;
        m_pfoEnergyTotal += pPfo->getEnergy();
        m_pfoPdgCodes.push_back(pPfo->getType());
        m_pfoEnergies.push_back(pPfo->getEnergy());

        m_pfoPx.push_back(pPfo->getMomentum()[0]);
        m_pfoPy.push_back(pPfo->getMomentum()[1]);
        m_pfoPz.push_back(pPfo->getMomentum()[2]);

        const float momentum(std::sqrt(pPfo->getMomentum()[0] * pPfo->getMomentum()[0] + pPfo->getMomentum()[1] * pPfo->getMomentum()[1] + pPfo->getMomentum()[2] * pPfo->getMomentum()[2]));
        const float cosTheta((momentum > std::numeric_limits<float>::epsilon()) ? pPfo->getMomentum()[2] / momentum : -999.f);
        m_pfoCosTheta.push_back(cosTheta);

        if (!pPfo->getTracks().empty())
        {
            // Charged pfos
            ++m_nPfosTracks;
            m_pfoEnergyTracks += pPfo->getEnergy();
        }
        else
        {
            // Neutral pfos
            float cellEnergySum(0.f);
            const EVENT::ClusterVec &clusterVec(pPfo->getClusters());

            for (EVENT::ClusterVec::const_iterator clustIter = clusterVec.begin(), clustIterEnd = clusterVec.end(); clustIter != clustIterEnd; ++clustIter)
            {
                const EVENT::CalorimeterHitVec &calorimeterHitVec((*clustIter)->getCalorimeterHits());

                for (EVENT::CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter)
                {
                    const EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;
                    cellEnergySum += pCalorimeterHit->getEnergy();
                }
            }

            const float correctionFactor((cellEnergySum < std::numeric_limits<float>::epsilon()) ? 0.f : pPfo->getEnergy() / cellEnergySum);

            if (22 == pPfo->getType())
            {
                // Photons
                ++m_nPfosPhotons;
                m_pfoEnergyPhotons += pPfo->getEnergy();

                for (EVENT::ClusterVec::const_iterator clustIter = clusterVec.begin(), clustIterEnd = clusterVec.end(); clustIter != clustIterEnd; ++clustIter)
                {
                    const EVENT::CalorimeterHitVec &calorimeterHitVec((*clustIter)->getCalorimeterHits());

                    for (EVENT::CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter)
                    {
                        const EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;

                        const float hitEnergy(correctionFactor * pCalorimeterHit->getEnergy());
                        const CHT cht(pCalorimeterHit->getType());

                        if (cht.is(CHT::ecal))
                        {
                            m_pfoECalToEmEnergy += hitEnergy;
                        }
                        else if (cht.is(CHT::hcal))
                        {
                            m_pfoHCalToEmEnergy += hitEnergy;
                        }
                        else if (cht.is(CHT::muon))
                        {
                            m_pfoMuonToEnergy += hitEnergy;
                        }
                        else
                        {
                            m_pfoOtherEnergy += hitEnergy;
                        }
                    }
                }
            }
            else
            {
                // Neutral hadrons
                ++m_nPfosNeutralHadrons;
                m_pfoEnergyNeutralHadrons += pPfo->getEnergy();

                for (EVENT::ClusterVec::const_iterator clustIter = clusterVec.begin(), clustIterEnd = clusterVec.end(); clustIter != clustIterEnd; ++clustIter)
                {
                    const EVENT::CalorimeterHitVec &calorimeterHitVec((*clustIter)->getCalorimeterHits());

                    for (EVENT::CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter)
                    {
                        const EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;

                        const float hitEnergy(correctionFactor * pCalorimeterHit->getEnergy());
                        const CHT cht(pCalorimeterHit->getType());

                        if (cht.is(CHT::ecal))
                        {
                            m_pfoECalToHadEnergy += hitEnergy;
                        }
                        else if (cht.is(CHT::hcal))
                        {
                            m_pfoHCalToHadEnergy += hitEnergy;
                        }
                        else if (cht.is(CHT::muon))
                        {
                            m_pfoMuonToEnergy += hitEnergy;
                        }
                        else
                        {
                            m_pfoOtherEnergy += hitEnergy;
                        }
                    }
                }
            }
        }

        momTot[0] += pPfo->getMomentum()[0];
        momTot[1] += pPfo->getMomentum()[1];
        momTot[2] += pPfo->getMomentum()[2];

        if (m_printing > 0)
        {
            std::cout << " RECO PFO, pdg: " << pPfo->getType() << ", E: " << pPfo->getEnergy() << ", nTracks: " << pPfo->getTracks().size()
                      << ", nClusters: " << pPfo->getClusters().size() << ", charge: " << pPfo->getCharge() << std::endl;
        }
    }

    m_pfoMassTotal = std::sqrt(m_pfoEnergyTotal * m_pfoEnergyTotal - momTot[0] * momTot[0] - momTot[1] * momTot[1] - momTot[2] * momTot[2]);

    // Extract quantities relating to pfo targets
    for (MCParticleVector::const_iterator iter = m_pfoTargetVector.begin(), iterEnd = m_pfoTargetVector.end(); iter != iterEnd; ++iter)
    {
        const EVENT::MCParticle *pMCParticle = *iter;

        ++m_nPfoTargetsTotal;
        m_pfoTargetsEnergyTotal += pMCParticle->getEnergy();
        m_pfoTargetPdgCodes.push_back(pMCParticle->getPDG());
        m_pfoTargetEnergies.push_back(pMCParticle->getEnergy());

        m_pfoTargetPx.push_back(pMCParticle->getMomentum()[0]);
        m_pfoTargetPy.push_back(pMCParticle->getMomentum()[1]);
        m_pfoTargetPz.push_back(pMCParticle->getMomentum()[2]);

        const float momentum(std::sqrt(pMCParticle->getMomentum()[0] * pMCParticle->getMomentum()[0] + pMCParticle->getMomentum()[1] * pMCParticle->getMomentum()[1] + pMCParticle->getMomentum()[2] * pMCParticle->getMomentum()[2]));
        const float cosTheta((momentum > std::numeric_limits<float>::epsilon()) ? pMCParticle->getMomentum()[2] / momentum : -999.f);
        m_pfoTargetCosTheta.push_back(cosTheta);

        if (std::fabs(cosTheta) > 0.98f)
            m_mcEnergyFwd += pMCParticle->getEnergy();

        if ((std::abs(pMCParticle->getPDG()) == 12) || (std::abs(pMCParticle->getPDG()) == 14) || (std::abs(pMCParticle->getPDG()) == 16))
            m_mcEnergyENu += pMCParticle->getEnergy();

        if (22 == pMCParticle->getPDG())
        {
            ++m_nPfoTargetsPhotons;
            m_pfoTargetsEnergyPhotons += pMCParticle->getEnergy();
        }
        else if ((11 == std::abs(pMCParticle->getPDG())) || (13 == std::abs(pMCParticle->getPDG())) || (211 == std::abs(pMCParticle->getPDG()))) // TODO, more options here?
        {
            ++m_nPfoTargetsTracks;
            m_pfoTargetsEnergyTracks += pMCParticle->getEnergy();
        }
        else
        {
            ++m_nPfoTargetsNeutralHadrons;
            m_pfoTargetsEnergyNeutralHadrons += pMCParticle->getEnergy();
        }
    }

    if (m_printing > 0)
    {
        std::cout << " EVENT                : " << m_nEvt << std::endl
                  << " NPFOs                : " << m_nPfosTotal << " (" << m_nPfosTracks << " + " << m_nPfosPhotons + m_nPfosNeutralHadrons << ")" << std::endl
                  << " RECONSTRUCTED ENERGY : " << m_pfoEnergyTotal << std::endl
                  << " RECO ENERGY + eNu    : " << m_pfoEnergyTotal + m_mcEnergyENu << std::endl;
    }

    // Fill basic histograms
    m_hPfoEnergySum->Fill(m_pfoEnergyTotal, 1.);

    if ((m_qPdg >= 1) && (m_qPdg <= 3) && (m_thrust <= 0.7f))
        m_hPfoEnergySumL7A->Fill(m_pfoEnergyTotal + m_mcEnergyENu, 1.);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TJjetsPFOAnalysisProcessor::SortPfoTargetsByEnergy(const EVENT::MCParticle *const pLhs, const EVENT::MCParticle *const pRhs)
{
    return (pLhs->getEnergy() > pRhs->getEnergy());
}
