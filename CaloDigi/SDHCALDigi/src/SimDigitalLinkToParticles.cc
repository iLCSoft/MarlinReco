#include "SimDigitalLinkToParticles.h"

#include <EVENT/CalorimeterHit.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>

#include <cassert>

SimDigitalLinkToParticles aSimDigitalLinkToParticles;

SimDigitalLinkToParticles::SimDigitalLinkToParticles() : Processor("SimDigitalLinkToParticles") {
  _description = "This processor links calorimeter hits to MC particles";

  std::vector<std::string> inputCollections = {};
  registerInputCollections(LCIO::CALORIMETERHIT, "inputHitCollections", "Sim Calorimeter Hit Collections",
                           _inputCollections, inputCollections);

  std::vector<std::string> inputRelCollections = {};
  registerInputCollections(LCIO::LCCOLLECTION, "inputRelationCollections",
                           "Collections of generic objects containing additional step informations",
                           _inputRelCollections, inputRelCollections);

  std::vector<std::string> outputRelCollections = {};
  registerProcessorParameter("outputRelationCollections", "output hit relation Collection Names", _outputRelCollections,
                             outputRelCollections);
}

void SimDigitalLinkToParticles::init() {
  assert(_inputCollections.size() == _inputRelCollections.size() &&
         _inputCollections.size() == _outputRelCollections.size());
}

void SimDigitalLinkToParticles::processEvent(LCEvent* evt) {
  for (unsigned int i(0); i < _inputCollections.size(); ++i) {
    try {
      std::string inputColName = _inputCollections.at(i);
      std::string inputRelColName = _inputRelCollections.at(i);
      std::string outputRelColName = _outputRelCollections.at(i);

      LCCollection* inputCol = evt->getCollection(inputColName.c_str());
      LCCollection* inputRelCol = evt->getCollection(inputRelColName.c_str());

      LCCollectionVec* outputRelCol = processCollection(inputCol, inputRelCol);

      evt->addCollection(outputRelCol, outputRelColName.c_str());
    } catch (DataNotAvailableException&) {
    }
  }
}

LCCollectionVec* SimDigitalLinkToParticles::processCollection(LCCollection* inputCol, LCCollection* inputRelCol) {
  LCFlagImpl flag;
  flag.setBit(LCIO::CHBIT_LONG);

  LCCollectionVec* outputRelCol = new LCCollectionVec(LCIO::LCRELATION);
  outputRelCol->setFlag(flag.getFlag());

  LCRelationNavigator navi(inputRelCol);

  int nHitsInCol = inputCol->getNumberOfElements();

  for (int i = 0; i < nHitsInCol; ++i) {
    CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>(inputCol->getElementAt(i));

    if (navi.getRelatedToObjects(hit).size() > 0) {
      SimCalorimeterHit* simHit = dynamic_cast<SimCalorimeterHit*>(navi.getRelatedToObjects(hit)[0]);

      std::map<MCParticle*, unsigned int> linkMap;
      unsigned int total = 0U;
      for (int j = 0; j < simHit->getNMCContributions(); ++j) {
        MCParticle* particle = simHit->getParticleCont(j);
        linkMap[particle]++;
        total++;
      }

      for (const auto& it : linkMap)
        outputRelCol->addElement(new LCRelationImpl(it.first, hit, 1.f * it.second / total));

    } else {
      streamlog_out(WARNING) << "could not find relation to sim calo hit !" << std::endl;
    }
  }

  return outputRelCol;
}
