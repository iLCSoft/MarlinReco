#ifndef DIGITIZER_REALISTICCALORECO_H
#define DIGITIZER_REALISTICCALORECO_H 1

#include "lcio.h"
#include "marlin/Processor.h"

#include <EVENT/CalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>

#include <UTIL/CellIDDecoder.h>

#include <IMPL/LCFlagImpl.h>

#include "CalorimeterHitType.h"

#include <string>
#include <vector>

using namespace lcio;
using namespace marlin;

/** === RealisticCaloReco Processor === <br>
    realistic reconstruction of calorimeter hits
    e.g. apply sampling fraction correction
    virtual class, technology indenpendent
    D.Jeans 02/2016.

    24 March 2016: removed gap corrections - to be put into separate processor
    changed relations: now keep relation between reconstructed and simulated hits.
*/

class RealisticCaloReco : virtual public Processor {

public:
  RealisticCaloReco();
  RealisticCaloReco(const RealisticCaloReco&) = delete;
  RealisticCaloReco& operator=(const RealisticCaloReco&) = delete;

  virtual void init();
  virtual void processRunHeader(LCRunHeader* run);
  virtual void processEvent(LCEvent* evt);
  virtual void check(LCEvent* evt);
  virtual void end();

protected:
  float getLayerCalib(int ilayer);
  virtual float reconstructEnergy(const CalorimeterHit* hit) = 0; // to be overloaded, technology-specific

  // parameters
  std::vector<std::string> _inputHitCollections{};
  std::vector<std::string> _inputRelCollections{};
  std::vector<std::string> _outputHitCollections{};
  std::vector<std::string> _outputRelCollections{};
  std::vector<float> _calibrCoeff{};
  std::vector<int> _calLayers{};

  std::string _cellIDLayerString{};

  // internal variables
  LCFlagImpl _flag{};
  LCFlagImpl _flag_rel{};

  CellIDDecoder<CalorimeterHit>* _idDecoder{};
};

#endif
