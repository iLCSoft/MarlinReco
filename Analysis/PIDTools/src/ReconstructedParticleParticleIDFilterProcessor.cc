#include "ReconstructedParticleParticleIDFilterProcessor.h"

#include <EVENT/LCEvent.h>
#include <EVENT/ReconstructedParticle.h>
#include <Exceptions.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCParametersImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/PIDHandler.h>

#include <algorithm>

ReconstructedParticleParticleIDFilterProcessor aReconstructedParticleParticleIDFilterProcessor;

struct LCParametersFilterAccessor : public IMPL::LCParametersImpl {
  void purgeAlgoMetaData(const std::string& algoName) {
    streamlog_out(DEBUG) << "Purging ParticleID meta information for algorithm " << algoName << '\n';

    auto paramNamesHandle = _stringMap.extract("ParameterNames_" + algoName);
    if (paramNamesHandle.empty()) {
      streamlog_out(DEBUG) << "Could not find ParameterNames_" << algoName << " to purge\n";
    } else {
      streamlog_out(DEBUG) << "Purged ParameterNames_" << algoName << '\n';
    }
  }
};

struct PIDHandlerFilterAccessor : public UTIL::PIDHandler {
  void purgeAlgoMetaData(const std::string& algoName) {
    auto& algoIDMap    = _cpm.map();
    auto  nameIDHandle = algoIDMap.extract(algoName);
    if (nameIDHandle.empty()) {
      streamlog_out(DEBUG) << "Could not find the meta information mapping the ID to the name for " << algoName
                           << " cannot purge it\n";
      return;
    } else {
      streamlog_out(DEBUG) << "Purged meta information for algorithm " << algoName << " from metadata\n";
    }

    // Take out the parameter names meta information as well, otherwise the
    // PIDHandler will re-insert them into the collection parameters
    auto& paramNameMap    = _pNames;
    auto  paramNameHandle = paramNameMap.extract(nameIDHandle.mapped());
    if (paramNameHandle.empty()) {
      streamlog_out(DEBUG) << "Could not find parameter names for " << algoName << " cannot purge them\n";
    } else {
      streamlog_out(DEBUG) << "Purged parameter names for " << algoName << " from metadata\n";
    }
  }
};

struct RecoParticleIDFilterAccessor : public IMPL::ReconstructedParticleImpl {
  void purgeParticleIDAlgorithm(const int algoId) {
    auto& pids   = _pid;  // Better debugging experience with a local variable
    auto  lastIt = std::remove_if(pids.begin(), pids.end(), [algoId](auto pid) {
      if (pid->getAlgorithmType() == algoId) {
        // We need to cleanup the ParticleID object here, because this is the
        // last time that we really have access to it, because remove_if does
        // not guarantee anything for the objects it "removes". Hence, a
        // followup loop to delete the pointers before we erase them does not
        // work, because we will be cleaning up the wrong things.
        delete pid;
        return true;
      }
      return false;
    });

    pids.erase(lastIt, pids.end());
  }
};

ReconstructedParticleParticleIDFilterProcessor::ReconstructedParticleParticleIDFilterProcessor()
    : marlin::Processor("ReconstructedParticleParticleIDFilterProcessor") {
  _description = "ReconstructedParticleParticleIDFilterProcessor: Filters ParticleID objects from ReconstructedParticles";

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "RecoParticleCollection", "collection to filter", m_inputCollName,
                          std::string("PandoraPFOs"));

  registerProcessorParameter("FilterPIDAlgos", "PID algorithm names to filter", m_filterPidAlgos,
                             std::vector<std::string>{});
}

void ReconstructedParticleParticleIDFilterProcessor::processEvent(LCEvent* event) {
  try {
    auto  recoColl      = dynamic_cast<IMPL::LCCollectionVec*>(event->getCollection(m_inputCollName));
    auto  pidHandler    = UTIL::PIDHandler(recoColl);
    auto& collParams    = static_cast<LCParametersFilterAccessor&>(recoColl->parameters());
    auto& filterHandler = static_cast<PIDHandlerFilterAccessor&>(pidHandler);

    for (const auto& algoName : m_filterPidAlgos) {
      streamlog_out(DEBUG) << "Now filtering algorithm " << algoName << std::endl;
      try {
        const auto pidAlgoId = pidHandler.getAlgorithmID(algoName);
        collParams.purgeAlgoMetaData(algoName);
        filterHandler.purgeAlgoMetaData(algoName);

        for (int i = 0; i < recoColl->getNumberOfElements(); ++i) {
          auto reco = static_cast<RecoParticleIDFilterAccessor*>(recoColl->getElementAt(i));
          reco->purgeParticleIDAlgorithm(pidAlgoId);
        }
      } catch (UnknownAlgorithm&) {
        streamlog_out(WARNING) << "Cannot purge algorithm " << algoName << " from " << m_inputCollName
                               << " because it is not known to be set on this collection" << std::endl;
      }
    }

  } catch (DataNotAvailableException& e) {
    streamlog_out(WARNING) << "Input collection " << m_inputCollName << " is not available" << std::endl;
  }
}
