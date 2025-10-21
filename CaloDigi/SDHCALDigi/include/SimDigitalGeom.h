#ifndef SimDigitalGeom_h
#define SimDigitalGeom_h

#include <utility>

#include <IMPL/LCCollectionVec.h>
#include <marlin/Processor.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

#include <EVENT/LCGenericObject.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>

#include "CalorimeterHitType.h" //in MarlinUtil
#include "marlinutil/LCGeometryTypes.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Factories.h"
#include "DDRec/DDGear.h"
#include "DDRec/DetectorData.h"
#include "DDRec/DetectorSurfaces.h"
#include "DDRec/MaterialManager.h"

namespace AIDA {
class ITuple;
}

struct StepAndCharge;

struct PotentialSameTrackID {
  PotentialSameTrackID(int _pdgStep, int _pdgParent) : PDGStep{_pdgStep}, PDGParent{_pdgParent} {}

  int PDGStep;
  int PDGParent;

  bool operator<(const PotentialSameTrackID& b) const {
    return std::tie(PDGStep, PDGParent) < std::tie(b.PDGStep, b.PDGParent);
  }
  bool operator==(const PotentialSameTrackID& b) const {
    return std::tie(PDGStep, PDGParent) == std::tie(b.PDGStep, b.PDGParent);
  }
};

class SimDigitalGeomCellId {
public:
  SimDigitalGeomCellId(LCCollection* inputCol, LCCollectionVec* outputCol);
  virtual ~SimDigitalGeomCellId();

  void setCellSize(float size) { _cellSize = size; }

  virtual float getCellSize() = 0;
  virtual void setLayerLayout(CHT::Layout layout) = 0;

  std::vector<StepAndCharge> decode(SimCalorimeterHit* hit, bool link);

protected:
  virtual void processGeometry(SimCalorimeterHit* hit) = 0;
  void createStepAndChargeVec(SimCalorimeterHit* hit, std::vector<StepAndCharge>& vec, bool link);

  void linkSteps(std::vector<StepAndCharge>& vec);

public:
  virtual std::unique_ptr<CalorimeterHitImpl> encode(int delta_I, int delta_J) = 0;

  int I() const { return _Iy; }
  int J() const { return _Jz; }
  int K() const { return _trueLayer; }
  int stave() const { return _stave; }
  int module() const { return _module; }
  int tower() const { return _tower; }

  const LCVector3D& normalToRPCPlane() const { return _normal; }
  const LCVector3D& Iaxis() const { return _Iaxis; }
  const LCVector3D& Jaxis() const { return _Jaxis; }

  SimDigitalGeomCellId(const SimDigitalGeomCellId& toCopy) = delete;
  void operator=(const SimDigitalGeomCellId& toCopy) = delete;

protected:
  CHT::Layout _currentHCALCollectionCaloLayout = CHT::any;

  dd4hep::CellID _cellIDvalue = 0;
  CellIDDecoder<SimCalorimeterHit> _decoder;
  CellIDEncoder<CalorimeterHitImpl> _encoder;

  float _cellSize = 0.0f;

  int _trueLayer = -999;
  int _stave = -999;
  int _module = -999;
  int _tower = -999;
  int _Iy = -999;
  int _Jz = -999;

  LCVector3D _normal;
  LCVector3D _Iaxis;
  LCVector3D _Jaxis;

  const float* _hitPosition = nullptr;

  std::string _cellIDEncodingString = "";

  // geometry debug tuples
public:
  static void bookTuples(const marlin::Processor* proc);

protected:
  void fillDebugTupleGeometryHit();
  void fillDebugTupleGeometryStep(SimCalorimeterHit* hit, const std::vector<StepAndCharge>& stepsInIJZcoord);

  static AIDA::ITuple* _tupleHit;
  enum {
    TH_CHTLAYOUT,
    TH_MODULE,
    TH_TOWER,
    TH_STAVE,
    TH_LAYER,
    TH_I,
    TH_J,
    TH_X,
    TH_Y,
    TH_Z,
    TH_NORMALX,
    TH_NORMALY,
    TH_NORMALZ,
    TH_IX,
    TH_IY,
    TH_IZ,
    TH_JX,
    TH_JY,
    TH_JZ
  };
  static AIDA::ITuple* _tupleStep;
  enum {
    TS_CHTLAYOUT,
    TS_HITCELLID,
    TS_NSTEP,
    TS_HITX,
    TS_HITY,
    TS_HITZ,
    TS_STEPX,
    TS_STEPY,
    TS_STEPZ,
    TS_DELTAI,
    TS_DELTAJ,
    TS_DELTALAYER,
    TS_TIME
  };
};

class SimDigitalGeomCellIdLCGEO : public SimDigitalGeomCellId {
public:
  SimDigitalGeomCellIdLCGEO(LCCollection* inputCol, LCCollectionVec* outputCol);
  virtual ~SimDigitalGeomCellIdLCGEO();

  virtual float getCellSize();
  virtual void setLayerLayout(CHT::Layout layout);

  virtual std::unique_ptr<CalorimeterHitImpl> encode(int delta_I, int delta_J);

  SimDigitalGeomCellIdLCGEO(const SimDigitalGeomCellIdLCGEO& toCopy) = delete;
  void operator=(const SimDigitalGeomCellIdLCGEO& toCopy) = delete;

protected:
  virtual void processGeometry(SimCalorimeterHit* hit);

  std::vector<std::string> _encodingString = {"layer", "stave", "module", "tower", "x", "y"};

  dd4hep::rec::LayeredCalorimeterData* _caloData = nullptr;
};

class SimDigitalGeomCellIdPROTO : public SimDigitalGeomCellId {
public:
  SimDigitalGeomCellIdPROTO(LCCollection* inputCol, LCCollectionVec* outputCol);
  virtual ~SimDigitalGeomCellIdPROTO();

  void setCellSize(float size) { _cellSize = size; }
  virtual float getCellSize() { return _cellSize; }
  virtual void setLayerLayout(CHT::Layout layout);

  virtual std::unique_ptr<CalorimeterHitImpl> encode(int delta_I, int delta_J);

  SimDigitalGeomCellIdPROTO(const SimDigitalGeomCellIdPROTO& toCopy) = delete;
  void operator=(const SimDigitalGeomCellIdPROTO& toCopy) = delete;

protected:
  virtual void processGeometry(SimCalorimeterHit* hit);

  std::vector<std::string> _encodingString = {"K-1", "", "", "", "I", "J"};
};

#endif // SimDigitalGeom_h
