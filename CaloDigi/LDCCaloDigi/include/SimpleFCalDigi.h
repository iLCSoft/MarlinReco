#ifndef SimpleFCalDigi_H
#define SimpleFCalDigi_H 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <string>
#include <vector>

using namespace lcio;
using namespace marlin;

/** === SimpleFCalDigi Processor === <br>
 *  Simple calorimeter digitizer for the LCal Processor. <br>
 *  Converts SimCalorimeterHit collections to one
 *  CalorimeterHit collection applying a threshold and an calibration constant...
 */

class SimpleFCalDigi : public Processor {

public:
  virtual Processor* newProcessor() { return new SimpleFCalDigi; }

  SimpleFCalDigi();

  virtual void init();

  virtual void processRunHeader(LCRunHeader* run);

  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  virtual void end();

protected:
  int _nRun{};
  int _nEvt{};

  std::vector<std::string> _fcalCollections{};

  std::string _outputFcalCollection{};
  std::string _outputRelCollection{};

  std::string _cellIDLayerString{};

  std::string _caloLayout{};
  std::string _caloID{};
  std::string _caloType{};

  float _thresholdFcal{};
  float _calibrCoeffFcal{};
  bool _fixLCalHits{};

  float xing_angle{};
  float zMin{};
  float dZ{};
  float rMin{};
  float cellDimR{};
  float cellDimPhi{};
  float WThickness{};
};

#endif
