#ifndef SimDigital_HHH
#define SimDigital_HHH

#include "lcio.h"
#include "marlin/Processor.h"
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>

#include <marlin/Global.h>

#include <TF2.h>
#include <TFile.h>
#include <TH1.h>

#include "CalorimeterHitType.h" //in MarlinUtil
#include "marlinutil/LCGeometryTypes.h"

#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogramFactory.h>

#include "ChargeInducer.h"
#include "ChargeSpreader.h"
#include "EfficiencyManager.h"
#include "SimDigitalGeom.h"

class TH1F;
class TF1;
class TTree;

namespace AIDA {
class ITuple;
}

using namespace lcio;
using namespace marlin;

/** Digitization for the SDHcal - based on NewLDCCaloDigi.
 *
 *  @author  G.Grenier, IPNL
 *  @author  R.Han, IPNL
 *	@author  A.Steen, IPNL
 *	@author  G.Garillot, IPNL
 *	@author  B.Li, IPNL
 *  @version $Id$
 */

struct StepAndCharge {
  StepAndCharge() : step{} {}
  StepAndCharge(LCVector3D vec, float _length, float _time) : step{vec}, stepLength{_length}, time{_time} {}
  LCVector3D step;
  float charge = 0;
  float stepLength = 0;
  float time = 0;
};

struct AsicKey {
  AsicKey(int l, int aI = -1, int aJ = -1) : layerID(l), asicI(aI), asicJ(aJ) {}
  int layerID;
  int asicI;
  int asicJ;

  bool operator<(const AsicKey& b) const {
    return std::tie(layerID, asicI, asicJ) < std::tie(b.layerID, b.asicI, b.asicJ);
  }
  bool operator==(const AsicKey& b) const {
    return std::tie(layerID, asicI, asicJ) == std::tie(b.layerID, b.asicI, b.asicJ);
  }
};

class SimDigital : public Processor {
public:
  virtual Processor* newProcessor() { return new SimDigital; }
  SimDigital();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent* evt);

  SimDigital(const SimDigital& toCopy) = delete;
  void operator=(const SimDigital& toCopy) = delete;

private:
  // intermediate storage class
  struct hitMemory {
    hitMemory() : ahit(nullptr), relatedHits(), maxEnergydueToHit(-1), rawHit(-1) {}

    std::unique_ptr<CalorimeterHitImpl> ahit = nullptr;

    std::set<int> relatedHits{};
    float maxEnergydueToHit = -1;
    int rawHit = -1;

    hitMemory(const hitMemory& other) = delete;
    hitMemory& operator=(const hitMemory& other) = delete;
  };

  typedef std::map<dd4hep::CellID, hitMemory> cellIDHitMap;

  void processCollection(LCCollection* inputCol, LCCollectionVec*& outputCol, LCCollectionVec*& outputRelCol,
                         CHT::Layout layout);
  cellIDHitMap createPotentialOutputHits(LCCollection* col, SimDigitalGeomCellId* aGeomCellId);

  void removeAdjacentStep(std::vector<StepAndCharge>& vec);
  void fillTupleStep(const std::vector<StepAndCharge>& vec, int level);
  void removeHitsBelowThreshold(cellIDHitMap& myHitMap, float threshold);
  void applyThresholds(cellIDHitMap& myHitMap);

  std::vector<std::string> _inputCollections{};

  std::vector<std::string> _outputCollections{};
  std::vector<std::string> _outputRelCollections{};

  LCFlagImpl flag{};
  LCFlagImpl flagRel{};

  std::map<std::string, int> _counters{};
  std::vector<float> _thresholdHcal{};

  std::vector<double> _hitCharge = {};

  std::map<dd4hep::CellID, std::vector<LCGenericObject*>> geneMap = {};

  float _cellSize = 0;
  float _gasGapWidth = 1.2f;

  // charge spreader
  std::string chargeSpreaderOption = "Uniform";
  std::string spreaderMapFile = "";
  ChargeSpreaderParameters chargeSpreaderParameters;
  ChargeSpreader* chargeSpreader = nullptr;

  std::string polyaOption = "Uniform";
  std::string polyaMapFile = "";
  float polyaQbar = 0.0f;
  float polyaTheta = 0.0f;
  ChargeInducer* chargeInducer = nullptr;
  int _polyaRandomSeed = 1;

  float _angleCorrPow = 0.4f;

  double timeCut = std::numeric_limits<double>::max();
  double stepLengthCut = -1.0;

  bool _linkSteps = false;
  bool _doThresholds = true;

  std::string efficiencyOption = "Uniform";
  std::string effMapFile = "";
  EfficiencyManager* efficiency = nullptr;
  float _constEffMapValue = 0.97f;

  float _absZstepFilter = 0.0005f;
  float _minXYdistanceBetweenStep = 0.5f;
  bool _keepAtLeastOneStep = true;

  AIDA::ITuple* _debugTupleStepFilter = nullptr;
  AIDA::ITuple* _tupleStepFilter = nullptr;
  AIDA::ITuple* _tupleCollection = nullptr;

  AIDA::IHistogram1D* _histoCellCharge = nullptr;

  std::string _encodingType = "LCGEO";
};

#endif
