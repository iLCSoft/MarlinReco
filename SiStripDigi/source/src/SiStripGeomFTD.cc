
#include "SiStripGeomFTD.h"

// Include standard header files
#include "Colours.h"
#include "PhysicalConstants.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>

// Include CLHEP header files
#include <CLHEP/Vector/EulerAngles.h>
#include <CLHEP/Vector/Rotation.h>

// Include Gear header files
#include <gear/BField.h>
// #include <gear/GearParameters.h>
#include <gearimpl/FTDLayerLayoutImpl.h>
#include <gearimpl/FTDParametersImpl.h>
#include <gearimpl/Vector3D.h>

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Include Encoder stuff
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>

namespace sistrip {

SiStripGeomFTD::SiStripGeomFTD(const std::string& detector, const double& pitchFront, const double& pitchRear)
    : SiStripGeom(detector), _ftdParams(0), _ftdLayer(0), _pitchFront(pitchFront), _pitchRear(pitchRear) {
  if (_gearType != "FTD") {
    streamlog_out(ERROR) << "Wrong instantiation of the class. "
                         << "See SiStripBuilder class, should have an incoherence!" << std::endl;
    exit(0);
  }
}

SiStripGeomFTD::~SiStripGeomFTD() {}

// Method initializing this class - reads Gear parameters from XML file
void SiStripGeomFTD::initGearParams() {

  // BField: needed for every tracker detector
  try {
    const gear::BField& bField = marlin::Global::GEAR->getBField();

    gear::Vector3D auxvec(0.0, 0.0, 0.0);

    _magField.setX((bField.at(auxvec)).x() * T);
    _magField.setY((bField.at(auxvec)).y() * T);
    _magField.setZ((bField.at(auxvec)).z() * T);
  } catch (gear::UnknownParameterException& e) {
    std::cout << "No magnetic field found in gear file!" << std::endl;
  }

  //------Get the geometry from the gear file-----//
  try {
    _ftdParams = const_cast<gear::FTDParameters*>(&(marlin::Global::GEAR->getFTDParameters()));
  } catch (gear::UnknownParameterException& e) {
    std::cout << "No FTD found in gear file" << std::endl;
    exit(-1);
  }

  try {
    _ftdLayer = const_cast<gear::FTDLayerLayout*>(&_ftdParams->getFTDLayerLayout());
  } catch (gear::UnknownParameterException& e) {
    std::cout << "Unexpected Error! No FTDLayerLayout found in"
              << " FTDParameters class" << std::endl;
    exit(-1);
  }

  // Initializing the data members, extracting from gear

  _numberOfLayers = _ftdLayer->getNLayers();
  // Layers goes from 0,..., 2N-1
  _layerType.reserve(2 * _numberOfLayers);
  _layerZ.reserve(2 * _numberOfLayers);
  _layerRadius.reserve(2 * _numberOfLayers);
  _layerPhi0.reserve(2 * _numberOfLayers);
  _layerHalfPhi.reserve(2 * _numberOfLayers);
  _layerPetalOpAngle.reserve(2 * _numberOfLayers);
  _layerOuterRadius.reserve(2 * _numberOfLayers);
  _sensorThick.reserve(2 * _numberOfLayers);
  _sensorWidth.reserve(2 * _numberOfLayers);
  _sensorWidth2.reserve(2 * _numberOfLayers);
  _sensorLength.reserve(2 * _numberOfLayers);
  _layerRealID.reserve(2 * _numberOfLayers);
  _ladderThick.reserve(2 * _numberOfLayers);
  _ladderOffsetZ.reserve(2 * _numberOfLayers);
  _ladderZOffsetSign0.reserve(2 * _numberOfLayers);
  _numberOfLadders.reserve(2 * _numberOfLayers);
  _ladderLength.reserve(2 * _numberOfLayers);
  for (int i = 0; i < _numberOfLayers; ++i) {
    _layerType.push_back(_ftdLayer->getSensorType(i));
    _layerZ.push_back(_ftdLayer->getZposition(i) * mm);
    _layerRadius.push_back(_ftdLayer->getSensitiveRinner(i) * mm);
    _layerPhi0.push_back(_ftdLayer->getPhi0(i));
    _layerHalfPhi.push_back(_ftdLayer->getPhiHalfDistance(i));
    _layerPetalOpAngle.push_back(_ftdLayer->getPhiHalfDistance(i));
    _layerOuterRadius.push_back(_ftdLayer->getSensitiveRinner(i) + _sensorLength[i] * mm);
    _layerRealID.push_back(i + 1);

    _sensorThick.push_back(_ftdLayer->getSensitiveThickness(i) * mm);
    _sensorWidth.push_back(_ftdLayer->getSensitiveLengthMax(i) * mm);  // x-direction
    _sensorWidth2.push_back(_ftdLayer->getSensitiveLengthMin(i) * mm); // x-direction
    _sensorLength.push_back(_ftdLayer->getSensitiveWidth(i) * mm);     // y-direction
    _ladderThick.push_back(_ftdLayer->getSupportThickness(i) * mm);
    _ladderZOffsetSign0.push_back(_ftdLayer->getZoffsetSign0(i));
    _ladderOffsetZ.push_back(_ftdLayer->getZoffset(i) * mm);
    _numberOfLadders.push_back(_ftdLayer->getNPetals(i));
    _ladderLength.push_back(_ftdLayer->getSupportWidth(i) * mm);
  }
  // Negative layers
  const unsigned int ridsize = _layerRealID.size();
  for (unsigned int i = 0; i < ridsize; ++i) {
    _layerRealID.push_back(-1 * _layerRealID.at(i));
  }
  _layerType.insert(_layerType.end(), _layerType.begin(), _layerType.end());
  _layerZ.insert(_layerZ.end(), _layerZ.begin(), _layerZ.end());
  _layerRadius.insert(_layerRadius.end(), _layerRadius.begin(), _layerRadius.end());
  _layerPhi0.insert(_layerPhi0.end(), _layerPhi0.begin(), _layerPhi0.end());
  _layerHalfPhi.insert(_layerHalfPhi.end(), _layerHalfPhi.begin(), _layerHalfPhi.end());
  _layerPetalOpAngle.insert(_layerPetalOpAngle.end(), _layerPetalOpAngle.begin(), _layerPetalOpAngle.end());
  _layerOuterRadius.insert(_layerOuterRadius.end(), _layerOuterRadius.begin(), _layerOuterRadius.end());
  _sensorThick.insert(_sensorThick.end(), _sensorThick.begin(), _sensorThick.end());
  _sensorWidth.insert(_sensorWidth.end(), _sensorWidth.begin(), _sensorWidth.end());
  _sensorWidth2.insert(_sensorWidth2.end(), _sensorWidth2.begin(), _sensorWidth2.end());
  _sensorLength.insert(_sensorLength.end(), _sensorLength.begin(), _sensorLength.end());
  _ladderThick.insert(_ladderThick.end(), _ladderThick.begin(), _ladderThick.end());
  _ladderOffsetZ.insert(_ladderOffsetZ.end(), _ladderOffsetZ.begin(), _ladderOffsetZ.end());
  _ladderZOffsetSign0.insert(_ladderZOffsetSign0.end(), _ladderZOffsetSign0.begin(), _ladderZOffsetSign0.end());
  _numberOfLadders.insert(_numberOfLadders.end(), _numberOfLadders.begin(), _numberOfLadders.end());
  _ladderLength.insert(_ladderLength.end(), _ladderLength.begin(), _ladderLength.end());

  _sensorPitchInFront = std::vector<double>(2 * _numberOfLayers, _pitchFront); //*um);
  _sensorNStripsInFront.reserve(2 * _numberOfLayers);
  _sensorPitchInRear = std::vector<double>(2 * _numberOfLayers, _pitchRear); //*um);
  _sensorNStripsInRear.reserve(2 * _numberOfLayers);
  for (int i = 0; i < _numberOfLayers; ++i) {
    // Pitch defined in the middle
    const double xmaxsensorup = _ftdLayer->getSensitiveLengthMax(i) * mm;
    const double ymiddlesensorup = _ftdLayer->getSensitiveWidth(i) * mm / 2.0;
    const double widthmiddle = xmaxsensorup - (2.0 * tan(_ftdLayer->getPhiHalfDistance(i)) * ymiddlesensorup);

    _sensorNStripsInFront.push_back((int)(widthmiddle / _sensorPitchInFront[i]) + 1);
    _sensorNStripsInRear.push_back((int)(widthmiddle / _sensorPitchInRear[i]) + 1);
  }
  _sensorNStripsInFront.insert(_sensorNStripsInFront.end(), _sensorNStripsInFront.begin(), _sensorNStripsInFront.end());
  _sensorNStripsInRear.insert(_sensorNStripsInRear.end(), _sensorNStripsInRear.begin(), _sensorNStripsInRear.end());
}

void SiStripGeomFTD::updateCanonicalCellID(const int& cellID, const int& stripType, const int& stripID,
                                           UTIL::BitField64* cellEnc) {
  std::map<std::string, int> bfmap = decodeCellID(cellID);
  (*cellEnc)["subdet"] = ILDDetID::FTD;
  short int layer = bfmap["layer"];
  short int module = bfmap["module"];
  short int sensor = bfmap["sensor"];

  // decodeCellID(layer,module,sensor,cellID);

  int realLayer = getLayerRealID(layer);

  (*cellEnc)["side"] = abs(realLayer) / realLayer;
  (*cellEnc)["layer"] = abs(realLayer);
  (*cellEnc)["module"] = module + 1; // Note that we want to keep the geant4 style
  (*cellEnc)["sensor"] = sensor;
  (*cellEnc)["stripType"] = stripType;
  (*cellEnc)["stripID"] = stripID;
}

//
//! Returns the input cellID0 where the field sensor is put to 0
int SiStripGeomFTD::cellID0withSensor0(const int& cellID0) const {
  UTIL::BitField64 cellDec(LCTrackerCellID::encoding_string());
  cellDec.setValue((lcio::long64)cellID0);

  cellDec["sensor"] = 0;

  return cellDec.lowWord();
}

//
// Decode Strip type and Strip ID (CellID1) (int version)
std::pair<StripType, int> SiStripGeomFTD::decodeStripID(const int& cellID1) const {
  UTIL::BitField64 cellDec("stripType:2,stripID:11");
  cellDec.setValue((lcio::long64)cellID1);

  StripType stripType = static_cast<StripType>(cellDec["stripType"].value());
  int stripID = cellDec["stripID"];

  return std::pair<StripType, int>(stripType, stripID);
}

//
// Decode Strip type and Strip ID (CellID1) (BitField version)
std::pair<StripType, int> SiStripGeomFTD::decodeStripID(const UTIL::BitField64& cellDec) const {
  if (cellDec["subdet"] != ILDDetID::FTD) {
    streamlog_out(ERROR) << "SiStripGeomFTD::decodeStripID "
                         << " - subdetector is not FTD!!" << std::endl;
    exit(-1);
  }

  StripType stripType = static_cast<StripType>(cellDec["stripType"].value());
  int stripID = cellDec["stripID"];

  return std::pair<StripType, int>(stripType, stripID);
}

// FIXME: TO BE PUT IN GEAR ??
// Returns the layer, module and sensor id. It returns in C-type indexs:
//---- The notation used in this code:
//     Layer:  0,...., 6 (Positive Z)
//             7,....,13 (Negative Z)
//     Petal:  0,....,15
//     Sensor: 1 2 3 4
// Argument codified as string
std::map<std::string, int> SiStripGeomFTD::decodeCellID(const int& cellID) const {
  UTIL::BitField64 cellDec(LCTrackerCellID::encoding_string());
  cellDec.setValue((lcio::long64)cellID);

  return decodeCellID(cellDec);
}

// FIXME: TO BE PUT IN GEAR ??
// Returns the layer, module and sensor id. It returns in C-type indexs:
//---- The notation used in this code:
//     Layer:  0,...., 6 (Positive Z)
//             7,....,13 (Negative Z)
//     Petal:  0,....,15
//     Sensor: 1 2 3 4
// Argument codified in the BitField64 class
std::map<std::string, int> SiStripGeomFTD::decodeCellID(const UTIL::BitField64& cellDec) const {
  // FIXME: Quizas incluir un decodeCellID0 en gear??

  // Checking we are using the canonical codification
  // FIXME:: los strings son diferentes!! pensar como solucionar
  /*if(cellDec.fielDescription() != LCTrackerCellID::encoding_string())
  {
          streamlog_out(ERROR) << "SiStripGeomFTD::decodeCellID0 "
                  << " - the codification used is not the canonical.\n "
                  << " Canonical codification string: " << cellDec.fieldDescription()
                  << " " << LCTrackerCellID::encoding_string()
                  << std::endl;
          exit(-1);
  }*/

  if (cellDec["subdet"] != ILDDetID::FTD) {
    streamlog_out(ERROR) << "SiStripGeomFTD::decodeCellID0 "
                         << " - subdetector is not FTD!!" << std::endl;
    exit(-1);
  }

  // Recall: codification in Mokka/G4:
  //     Layer: 1,...7; Side:  -1,+1; Petal: 1...16; Sensor: 1,2,3,4
  std::map<std::string, int> codmap;
  codmap["layer"] = getLayerIDCTypeNo(cellDec["side"] * cellDec["layer"]);
  codmap["module"] = cellDec["module"] - 1;
  codmap["sensor"] = cellDec["sensor"];

  return codmap;
}

// Added
double SiStripGeomFTD::getLadderOffsetX(const short int& layerID) const {

  if (_sensorWidth.size() < (unsigned int)layerID) {
    streamlog_out(ERROR) << "SiStripGeomFTD::getLadderOffsetX - layerID: " << layerID << " out of range!!!"
                         << std::endl;
  }

  return _sensorWidth[layerID] / 2.0;
}

// TRANSFORMATION METHODS - GLOBAL REF. SYSTEM

// The Local Reference Frame system (LRF) is defined in the way that:
//  - is defined positive
//  - the electric drift field is oriented in the negative direction
//    of the X-axis
//  - all the hits inside a sensors have positive local coordinates
//  - the strips collecting holes (Single Side as Double side sensors)
//    are faced to the IP
//  - the strips collecting holes are in the sensor plane X=0 (the holes
//    are collected by the strips oriented in Z-direction).
// These requerimemts make placing the LRF:
//
//                             --------
//                             |      |              Face near IP
//    1 and 2 sensors           |    |
//                               |__| X <---
//
//
//                             --------
//                             |      |
//    3 and 4 sensors           |    |               Face far IP
//                        ---> X |__|
//

//
// Method transforming given point from global ref. system to local ref. system
// (parameters: diskID, petalID, sensorID and space point in global ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformPointToLocal(short int diskID, short int petalID, short int sensorID,
                                                        const CLHEP::Hep3Vector& globalPoint) {
  // Initialize local point
  CLHEP::Hep3Vector localPoint(globalPoint);

  // Phi Angle of the Petal (on the centroid of the petal with respect the x-axis)
  const double phi0 = _ftdLayer->getPhiPetalCd(diskID, petalID);

  //
  // RECALL: ftd gear have the convention mm (distance), here we have cm
  //

  // Extract the Z of the sensor: diferent position depending the sensor
  const double zsensorCd = _ftdLayer->getSensitiveZposition(diskID, petalID, sensorID) * mm;
  const int zsign = (int)(zsensorCd / fabs(zsensorCd));
  const double sensorthickness = _ftdLayer->getSensitiveThickness(diskID) * mm;
  // Sensor 3 and 4: Displacing to the trapezoid farest the IP
  double zsensor = zsensorCd + zsign * sensorthickness / 2.0;
  // Sensor 1 and 2: Displacing to the trapezoid facing the IP
  if (sensorID < 3) {
    zsensor = zsensorCd - zsign * sensorthickness / 2.0;
  }

  // And the X and Y CentroiD position: over the smallest side of the trapezoid
  const int xsign = (int)(fabs(globalPoint.getX()) / globalPoint.getX());
  const double xsensorCd = xsign * (_ftdLayer->getSensitiveRinner(diskID) * mm) * fabs(cos(phi0));

  const int ysign = (int)(fabs(globalPoint.getY()) / globalPoint.getY());
  const double ysensorCd = ysign * (_ftdLayer->getSensitiveRinner(diskID) * mm) * fabs(sin(phi0));

  // The (0,0,0) position of LFR
  CLHEP::Hep3Vector localOrigin(xsensorCd, ysensorCd, zsensor);

  // Translating the globalPoint to the local Origin
  localPoint -= localOrigin;

  // Perform rotation to get the system local:
  //  Z-positives sensors 1,2 |
  //  Z-negatives sensors 3,4 |-- Same transformation
  //
  //  Z-positives sensors 3,4 | ---> La antigua
  //  Z-negatives sensors 1,2 |-- Same transformation
  double rotZangle = -phi0;
  double rotYangle = -M_PI / 2.0;
  if ((sensorID < 3 && zsign > 0) || (sensorID > 2 && zsign < 0)) {
    rotZangle = M_PI - phi0;
    rotYangle = M_PI / 2.0;
  }
  localPoint.rotateZ(rotZangle);
  localPoint.rotateY(rotYangle);

  // Avoiding X,Y and Z negative values--> displacing from the
  // centroid to the edge
  const double longtrapezoidedge = _ftdLayer->getSensitiveLengthMax(diskID) * mm;
  localPoint += CLHEP::Hep3Vector(0.0, longtrapezoidedge / 2., 0.0);

  streamlog_out(DEBUG4) << "============================================\n"
                        << "SiStripGeomFTD::transformPointToLocal \n"
                        << " Sensor ID (C-vector style): \n "
                        << "   DISK: " << diskID << " PETAL: " << petalID << " SENSOR: " << sensorID << "\n"
                        << "   Petal Phi: " << phi0 * 180.0 / M_PI << "\n"
                        << std::setprecision(3) << " Origen of Local frame [mm]: (" << xsensorCd / mm << ","
                        << ysensorCd / mm << "," << zsensor / mm << ")\n"
                        << " Hit (Global ref. frame) [mm]:" << globalPoint / mm << "\n"
                        << " Hit (Local ref. frame)  [mm]:" << localPoint / mm << "\n"
                        << " Maximum dimensions: x < " << sensorthickness / mm << ", y < " << longtrapezoidedge / mm
                        << ", z < " << _ftdLayer->getSensitiveWidth(diskID) * mm / mm << " [mm]\n"
                        << "============================================" << std::endl;

  // Return space point in local ref. system
  return localPoint;
}

//
// Method transforming given vector from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and vector in global ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformVecToLocal(short int diskID, short int petalID, short int sensorID,
                                                      const CLHEP::Hep3Vector& globalVec) {
  // Initialize local vector
  CLHEP::Hep3Vector localVec(globalVec);

  const double phi0 = _ftdLayer->getPhiPetalCd(diskID, petalID);
  // Extract z-sign and
  const int realLayerID = getLayerRealID(diskID);
  const int zsign = abs(realLayerID) / realLayerID;
  // Sensors 3,4 z-positive and 1,2 z-negative
  double rotZangle = -phi0;
  double rotYangle = -M_PI / 2.0;
  // Sensors 1,2 z-positive and 3,4 z-negative
  if ((sensorID < 3 && zsign > 0) || (sensorID > 2 && zsign < 0)) {
    rotZangle = M_PI - phi0;
    rotYangle = M_PI / 2.0;
  }
  localVec.rotateZ(rotZangle);
  localVec.rotateY(rotYangle);

  // Return vector in local ref. system
  return localVec;
}

//
// Method transforming given matrix 3x3 from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and matrix in global ref. system)
//
CLHEP::HepMatrix SiStripGeomFTD::transformMatxToLocal(short int layerID, short int ladderID,
                                                      const CLHEP::HepMatrix& globalMatrix) {
  // FIXME: NOT DONE YET
  std::cout << "SiStripGeomFTD::transformMatxToLocal -- METHOD NOT IMPLEMENTED" << std::endl;
  exit(0);
  // Initialize local matrix 3x3 to zero values
  CLHEP::HepMatrix localMatrix(3, 3, 0);

  // Initialize rotation matrices: R, R^T (transposition)
  CLHEP::HepMatrix rotMatrix(3, 3, 0);
  CLHEP::HepMatrix rotMatrixT(3, 3, 0);

  // Gear type: VXD
  if (_gearType == "VXD") {

    // Calculate rotation angles
    double theta = getLadderTheta(layerID); // ALWAYS M_PI/2
    double phi = getLadderPhi(layerID, ladderID);

    // Calculate rotation matrices - help
    CLHEP::HepRotation rotMatrixZ(CLHEP::Hep3Vector(0, 0, 1), -phi);
    CLHEP::HepRotation rotMatrixY(CLHEP::Hep3Vector(0, 1, 0), +theta);
    CLHEP::HepRotation rotMatrixHelp(rotMatrixY * rotMatrixZ);

    CLHEP::HepRotation rotMatrixHelpT(rotMatrixHelp.inverse());

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        rotMatrix[i][j] = rotMatrixHelp[i][j];
        rotMatrixT[i][j] = rotMatrixHelpT[i][j];
      }
    }
  }

  // Gear type: unknown - error
  else {
    streamlog_out(ERROR) << "Unknown gear type!" << std::endl;

    exit(0);
  }

  // Transform given matrix - rotation wrt global ref. system
  localMatrix = rotMatrix * globalMatrix * rotMatrixT;

  // Return matrix in local ref. system
  return localMatrix;
}

// TRANSFORMATION METHODS - LOCAL REF. SYSTEM

//
// Method transforming given point from local ref. system to global ref. system
// (parameters: layerID, ladderID, sensorID and space point in local ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformPointToGlobal(short int diskID, short int petalID, short int sensorID,
                                                         const CLHEP::Hep3Vector& localPoint) {
  // Initialize global point
  CLHEP::Hep3Vector globalPoint(localPoint);

  // Calculate rotation angles
  double theta = getLadderTheta(diskID); // ALWAYS M_PI/2.0
  double phi = getLadderPhi(diskID, petalID);

  // Taking account of the z-side
  const int realLayerID = getLayerRealID(diskID);
  const int zsign = abs(realLayerID) / realLayerID;
  double rotZangle = -phi;
  double rotYangle = -theta;
  if ((sensorID < 3 && zsign > 0) || (sensorID > 2 && zsign < 0)) {
    rotZangle = M_PI - phi;
    rotYangle = theta;
  }

  // Perform translation - to the center of a petal (in local frame)
  // Note that the center of the petal we decide defining in the back face
  // of the sensor
  const double xlocalCd = 0.0; // Already in the back face_ftdLayer->getSensitiveThickness(diskID)*mm;
  const double ylocalCd = _ftdLayer->getSensitiveLengthMax(diskID) / 2.0 * mm;
  const double zlocalCd = 0.0; // already in the low edge of the petal (see ysensorCd)

  globalPoint -= CLHEP::Hep3Vector(xlocalCd, ylocalCd, zlocalCd);

  // Perform rotation - with respect to the local centre of a ladder
  globalPoint.rotateY(-rotYangle);
  globalPoint.rotateZ(-rotZangle);
  // Now we have the vector in global coordinates from the center of the petal
  // to the hit

  // Find (0,0,0) position of local coordinate system (in global coord.)
  // Extract the Z of the sensor
  const double zsensorCd = _ftdLayer->getSensitiveZposition(diskID, petalID, sensorID) * mm;
  const double sensorthickness = _ftdLayer->getSensitiveThickness(diskID) * mm;
  // Sensor 3 and 4: Displacing to the trapezoid farest the IP
  double zsensor = zsensorCd + zsign * sensorthickness / 2.0;
  // Sensor 1 and 2: Displacing to the trapezoid facing the IP
  if (sensorID < 3) {
    zsensor = zsensorCd - zsign * sensorthickness / 2.0;
  }
  /* OLD
  // Sensors 3 and 4: Displacing to the backed the IP face
  double zsensorback = zsensor+zsign*sensorthickness/2.0;
  if( sensorID < 3 )
  {
          zsensorback = zsensor-zsign*sensorthickness/2.0;
  } END OLD*/
  // And the X and Y CentroiD position: over the smallest side of the trapezoid
  const double xsensorCd = _ftdLayer->getSensitiveRinner(diskID) * mm * cos(phi);
  const double ysensorCd = _ftdLayer->getSensitiveRinner(diskID) * mm * sin(phi);
  CLHEP::Hep3Vector localOrigin(xsensorCd, ysensorCd, zsensor);

  // Perform translation - to the global system
  globalPoint += localOrigin;

  streamlog_out(DEBUG4) << "============================================\n"
                        << "SiStripGeomFTD::transformPointToGlobal \n"
                        << " Sensor ID (C-vector style): \n "
                        << "   DISK: " << diskID << " PETAL: " << petalID << " SENSOR: " << sensorID << "\n"
                        << "   Petal Phi: " << phi * 180.0 / M_PI << "\n"
                        << std::setprecision(3) << " Origen of Local frame [mm]: (" << xsensorCd / mm << ","
                        << ysensorCd / mm << "," << zsensor / mm << ")\n"
                        << " Hit (Global ref. frame) [mm]:" << globalPoint / mm << "\n"
                        << " Hit (Local ref. frame)  [mm]:" << localPoint / mm << "\n"
                        << "============================================" << std::endl;

  // Return space point in global ref. system
  return globalPoint;
}

//
// Method transforming given vector from local ref. system to global ref. system
// (parameters: layerID, ladderID, sensorID and vector in local ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformVecToGlobal(short int diskID, short int petalID, short int sensorID,
                                                       const CLHEP::Hep3Vector& localVec) {
  // Initialize global vector
  CLHEP::Hep3Vector globalVec(localVec);

  // Calculate rotation angles
  double theta = getLadderTheta(diskID); // ALWAYS PI/2
  double phi = getLadderPhi(diskID, petalID);

  const int realLayerID = getLayerRealID(diskID);
  const int zsign = abs(realLayerID) / realLayerID;
  double rotZangle = -phi;
  double rotYangle = -theta;
  if ((sensorID < 3 && zsign > 0) || (sensorID > 2 && zsign < 0)) {
    rotZangle = M_PI - phi;
    rotYangle = theta;
  }
  globalVec.rotateY(-rotYangle);
  globalVec.rotateZ(-rotZangle);

  // Return vector in global ref. system
  return globalVec;
}

//
// Method transforming given matrix 3x3 from local ref. system to global ref. system
// (parameters: layerID, ladderID, sensorID and matrix in local ref. system)
//
CLHEP::HepMatrix SiStripGeomFTD::transformMatxToGlobal(short int layerID, short int ladderID,
                                                       const CLHEP::HepMatrix& localMatrix) {
  // FIXME: NOT DONE YET
  std::cout << "SiStripGeomFTD::transformMatxToGlobal -- METHOD NOT IMPLEMENTED" << std::endl;
  exit(0);
  // Initialize local matrix 3x3 to zero values
  CLHEP::HepMatrix globalMatrix(3, 3, 0);

  // Initialize rotation matrices: R, R^T (transposition)
  CLHEP::HepMatrix rotMatrix(3, 3, 0);
  CLHEP::HepMatrix rotMatrixT(3, 3, 0);

  // Gear type: VXD
  if (_gearType == "VXD") {

    // Calculate rotation angles
    double theta = getLadderTheta(layerID);
    double phi = getLadderPhi(layerID, ladderID);

    // Calculate rotation matrices - help
    CLHEP::HepRotation rotMatrixY(CLHEP::Hep3Vector(0, 1, 0), -theta);
    CLHEP::HepRotation rotMatrixZ(CLHEP::Hep3Vector(0, 0, 1), +phi);
    CLHEP::HepRotation rotMatrixHelp(rotMatrixZ * rotMatrixY);

    CLHEP::HepRotation rotMatrixHelpT(rotMatrixHelp.inverse());

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        rotMatrix[i][j] = rotMatrixHelp[i][j];
        rotMatrixT[i][j] = rotMatrixHelpT[i][j];
      }
    }
  }

  // Gear type: unknown - error
  else {
    streamlog_out(ERROR) << "Unknown gear type!" << std::endl;

    exit(0);
  }

  // Transform given matrix - rotation wrt local ref. system
  globalMatrix = rotMatrix * localMatrix * rotMatrixT;

  // Return matrix in global ref. system
  return globalMatrix;
}

// OTHER METHODS - GLOBAL REF. SYSTEM

//
// Get info whether the given point is inside of Si sensor (parameters: layerID,
// space point in local ref. system)
//
bool SiStripGeomFTD::isPointInsideSensor(short int diskID, short int petalID, short int sensorID,
                                         const CLHEP::Hep3Vector& point) const {
  // Angle defining the petals
  const double phi0 = _ftdLayer->getPhiHalfDistance(diskID);

  // Calculate the maximum of each axis
  const double xmax = _ftdLayer->getSensitiveThickness(diskID) * mm;
  const double ygap = (_ftdLayer->getSensitiveWidth(diskID) * mm - point.getZ()) * tan(phi0);
  const double ymax = _ftdLayer->getSensitiveLengthMax(diskID) * mm - ygap;
  const double zmax = _ftdLayer->getSensitiveWidth(diskID) * mm;

  bool isIn = true;
  // Boundary set +- epsilon
  if ((point.getX() > xmax + EPS * um) || (point.getX() < -EPS * um)) {
    streamlog_out(ERROR) << std::setprecision(3) << " position out of sensor! X [um] = " << point.getX() / um
                         << " (Maximum x: " << xmax / um << ")" << std::endl;
    isIn = false;
  }

  if ((point.getY() > ymax + EPS * um) || (point.getY() < ygap - EPS * um)) {
    // FIXME: HAY UN PUNTO EN simHitFTD.slcio que sale del detector: investigar
    streamlog_out(ERROR) << std::setprecision(3) << " position out of sensor! Y [mm] = " << point.getY() / mm
                         << " (Maximum y: " << ymax / mm << ", Minimum y:" << ygap / mm << ")" << std::endl;
    isIn = false;
  }

  if ((point.getZ() > zmax + EPS * um) || (point.getZ() < -EPS * um)) {
    streamlog_out(ERROR) << std::setprecision(3) << " position out of sensor! Z [mm] = " << point.getZ() / mm
                         << " (Maximum z: " << zmax / mm << ")" << std::endl;
    isIn = false;
  }

  // Return if out or not
  return isIn;
}

// OTHER METHODS - LOCAL REF. SYSTEM

//
// Get info whether the given point is out of Si sensor (parameters: layerID,
// space point in local ref. system)
//
bool SiStripGeomFTD::isPointOutOfSensor(short int layerID, const CLHEP::Hep3Vector& point) const {

  bool isOut = false;

  // Boundary set +- epsilon
  double tanAlpha = (getSensorWidthMax(layerID) - getSensorWidthMin(layerID)) / 2. / getSensorLength(layerID);
  double actualSensorWidth = (getSensorWidthMax(layerID) - 2 * tanAlpha * point.getZ());

  // Recalculate point into "local" local system (width depends on posZ)
  CLHEP::Hep3Vector recalcPoint(point.getX(), point.getY() - getSensorWidthMax(layerID) / 2. + actualSensorWidth / 2.,
                                point.getZ());

  if ((recalcPoint.getX() > (getSensorThick(layerID) + EPS * um)) || (recalcPoint.getX() < (-EPS * um)) ||
      (recalcPoint.getY() > (actualSensorWidth + EPS * um)) || (recalcPoint.getY() < (-EPS * um)) ||
      (recalcPoint.getZ() > (getSensorLength(layerID) + EPS * um)) || (recalcPoint.getZ() < (-EPS * um))) {
    isOut = true;
  }

  // Return if out or not
  return isOut;
}

//
// Get number of strips in sensors
//
int SiStripGeomFTD::getSensorNStrips(const int& diskID, const int& sensorID) const {
  std::vector<int> const* sensorNStrips = 0;

  if (sensorID == 1 || sensorID == 2) {
    sensorNStrips = &_sensorNStripsInFront;
  } else if (sensorID == 3 || sensorID == 4) {
    sensorNStrips = &_sensorNStripsInRear;
  } else {
    streamlog_out(ERROR) << "SiStripGeomFTD::getSensorNStrips "
                         << " - Inconsistent ID of sensor: '" << sensorID << "'" << std::endl;
    exit(-1);
  }

  if (sensorNStrips->size() > (unsigned short int)diskID) {
    return (*sensorNStrips)[diskID];
  } else {
    return 0;
  }
}

//
// Get strip ID for RPhi  (perpendicular to the beam axis),
// points are given in local ref. system.
//
// NOTES: Necesitem la posRPhi per designar quin  strip estic utilitzant,
//       la pos Z es necesaria per definir a quina zona del trapezoide
//       estem, recorda que la Pitch= Pitch(Z)
/*int SiStripGeomFTD::getStripIDInRPhi(short int diskID, double posRPhi, double posZ ) const
{
        // Get pitch
        double sensPitch = getSensorPitchInRPhi(diskID,posZ);
        if(sensPitch == 0)
        {
                streamlog_out(ERROR) << "SiStripGeomFTD::getStripIDInRPhi "
                        << "- division by zero (sensPitch is zero)!!!"
                        << std::endl;
                exit(-1);
        }
        // Get number of strips
        int sensNStrips = getSensorNStripsInRPhi(diskID);

        int stripID;
        // Calculate stripID
        if(posRPhi <= 0.)
        {
                stripID = 0;
        }
        else
        {
                stripID = floor(posRPhi/sensPitch);
                if(stripID >= sensNStrips)
                {
                        stripID = sensNStrips - 1;
                }
        }
        // Error
        if(stripID >= sensNStrips)
        {
                streamlog_out(ERROR) << "SiStripGeom::getStripIDInRPhi "
                        << "- stripID in RPhi greater than number of strips!!!"
                        << std::endl;
                exit(-1);
        }
std::cout << " SiStripGeomFTD::getStripIDInRPhi" << std::endl;
std::cout << "----> " << posRPhi << std::endl;
std::cout << "----> " << sensPitch << std::endl;
std::cout << "----> " << sensNStrips << std::endl;
std::cout << "----> " << stripID << " (en zlocal=" << posZ << ")" << std::endl;
std::cout << "----> " << sensPitch*stripID << " (Y-position)" << std::endl;
std::cout << "END SiStripGeomFTD::getStripIDInRPhi" << std::endl;
        return stripID;
}*/
//! Get director vector of a strip (inside the local Ref. system)
//! The vector is defined to describe the strip from z_local=0
CLHEP::Hep3Vector SiStripGeomFTD::getStripUnitVector(const int& diskID, const int& sensorID, const int& stripID) const {
  // Extract y in z=0 and in z=L
  const double yatz0 = getStripPosY(diskID, sensorID, stripID, 0.0);
  const double yatzL = getStripPosY(diskID, sensorID, stripID, getSensorLength(diskID));
  // Defining positive u-component if yatz0 > Width/2.0
  const double Dy = yatzL - yatz0;
  const int u_sign = fabs(Dy) / Dy;

  const double alpha = fabs(atan(Dy / getSensorLength(diskID)));
  // FIXME: Control that alpha<pi/2

  return CLHEP::Hep3Vector(0.0, ((double)u_sign) * sin(alpha), cos(alpha));
}

//
// Get the point crossing two strips. The returned point is placed in the middle
// of the Front sensor (x=getSensorThick/2.0)
//
CLHEP::Hep3Vector SiStripGeomFTD::getCrossLinePoint(const int& diskID, const int& petalID, const int& stripIDFront,
                                                    const int& stripIDRear) const {
  // Extract y in z=0 in front and rear
  const double yatz0Front = getStripPosY(diskID, 1, stripIDFront, 0.0);
  const double yatz0Rear = getStripPosY(diskID, 3, stripIDRear, 0.0);
  // And putting them in the same reference system (the front ref. local system)
  const double yatz0RearPrime = getSensorWidthMax(diskID) - yatz0Rear;

  // Director vectors of the strips
  CLHEP::Hep3Vector uFront = getStripUnitVector(diskID, 1, stripIDFront);
  CLHEP::Hep3Vector uRear = getStripUnitVector(diskID, 3, stripIDRear);
  // Reference system front
  uRear.rotateZ(-M_PI);

  // Strip line parametrized with t => r = (yatz0,0)+t(uz,uy)
  const double extProd = uFront.getY() * uRear.getZ() - uFront.getZ() * uRear.getY();
  const double tR = (yatz0RearPrime - yatz0Front) * uFront.getZ() / extProd;

  const double y = yatz0RearPrime + tR * uRear.getY();
  const double z = tR * uRear.getZ();

  return CLHEP::Hep3Vector(getSensorThick(diskID) / 2.0, y, z);
}

//
// Get the stereo angle given a sensor
//
double SiStripGeomFTD::getStereoAngle(const int& diskID, const int& sensorID) const {
  // Stereo angle (if proceed)
  double stAngle = 0.0;

  if (sensorID == 3 || sensorID == 4) {
    stAngle = 50e-3;
  }

  return stAngle;
}

// double SiStripGeomFTD::get

//
// Get strip ID, point is given in local ref. system
//
//
int SiStripGeomFTD::getStripID(const int& diskID, const int& sensorID, const double& posRPhi,
                               const double& posZ) const {
  // FIXME: BORRA
  double zRot = 0.0;
  double posRPhiRot = 0.;
  double yorigen = 0.0;
  // Get pitch
  double sensPitch = getSensorPitch(diskID, sensorID, posZ);
  if (sensPitch < 1e-40) {
    streamlog_out(ERROR) << "SiStripGeomFTD::getStripID "
                         << "- division by zero (sensPitch is zero)!!!" << std::endl;
    exit(-1);
  }
  // Get number of strips
  int sensNStrips = getSensorNStrips(diskID, sensorID);

  int stripID;
  // Calculate stripID
  if (posRPhi <= 0.) {
    stripID = 0;
  } else {
    CLHEP::Hep3Vector pointRot = transformPointToRotatedLocal(diskID, sensorID, CLHEP::Hep3Vector(0.0, posRPhi, posZ));
    posRPhiRot = pointRot.getY();
    zRot = pointRot.getZ();
    // Put the strip 0 in the edge of the petal (recall
    // the y=0 of the local system begins in the long edge)
    const double tanPhi = tan(_ftdLayer->getPhiHalfDistance(diskID));
    yorigen = (getSensorLength(diskID) - pointRot.getZ()) * tanPhi;

    stripID = floor((posRPhiRot - yorigen) / sensPitch);
    if (stripID >= sensNStrips) {
      stripID = sensNStrips - 1;
    }
  }
  // Error
  if (stripID >= sensNStrips) {
    streamlog_out(ERROR) << "SiStripGeom::getStripID "
                         << "- stripID in RPhi greater than number of strips!!!" << std::endl;
    exit(-1);
  }

  /*std::cout << "-----+ SiStripGeomFTD::getStripID" << std::endl;
  std::cout << "----> Y   =" << posRPhi << " (y/Pitch):" << (posRPhi-yorigen)/sensPitch << std::endl;
  std::cout << "----> YRot =" << posRPhiRot << " (y/Pitch):" << (posRPhiRot-yorigen)/sensPitch << std::endl;
  std::cout << "----> Pitch=" << sensPitch << std::endl;
  std::cout << "----> Number Strips:" << sensNStrips << std::endl;
  std::cout << "----> StripID:" << stripID << std::endl;
  std::cout <<  "----> yOrigenRot:" << yorigen << std::endl;
  std::cout <<  "----> zRot:" << zRot << std::endl;
  std::cout << "----> NS*sID=" << (yorigen+sensPitch*stripID) << " (Y-position)" << std::endl;
  std::cout << "-----+END SiStripGeomFTD::getStripIDInRPhi +------" << std::endl;*/
  return stripID;
}

// METHODS TO IDENTIFY

//
// Get sensor pitch for RPhi strips (perpendicular to the beam axis).
// Currently it means the front face to the IP sensors of the petals.
// The posZ is needed because the trapezoid shape of the sensors therefore pitch=pitch(Z)
//
/*double SiStripGeomFTD::getSensorPitchInRPhi(short int diskID, double posZ) const
{
        if( (getLayerType(diskID) == strip) &&
                        (_sensorPitchInFront.size() > (unsigned int)diskID) )
        {
                const double tanAlpha = tan(_ftdLayer->getPhiHalfDistance(diskID));
                if( (posZ*um>-EPS*um) && (posZ<(getSensorLength(diskID) + EPS*um)) ) //FIXME::::
                {
                        return ((getSensorWidth(diskID) - 2.0*tanAlpha*posZ)
                                        /getSensorNStripsInRPhi(diskID));
                }
                else {
        std::cout << getSensorLength(diskID) << std::endl;
        std::cout << posZ << "   -EPS:" << -EPS*um << " Cumple?: " << (posZ*um>-EPS*um) << std::endl;
                        streamlog_out(ERROR)
                                << std::setiosflags(std::ios::fixed | std::ios::internal )
                                << std::setprecision(2)
                                << "SiStripGeomFTD::getSensorPitchInRPhi - posZ: "
                                << std::setw(4) << posZ << " out of range!!!"
                                << std::resetiosflags(std::ios::showpos)
                                << std::setprecision(0)
                                << std::endl;
                        exit(0);;
                }
        }
        else
        {
                return 0.;
        }
}*/

//
// Get sensor pitch
// The posZ is needed because the trapezoid shape of the sensors therefore pitch=pitch(Z)
// The posZ is in the local ref. system
//
double SiStripGeomFTD::getSensorPitch(const int& diskID, const int& sensorID, const double& posZ) const {
  if ((posZ < -EPS * um) && (posZ > (getSensorLength(diskID) + EPS * um))) {
    streamlog_out(ERROR) << std::setiosflags(std::ios::fixed | std::ios::internal) << std::setprecision(2)
                         << "SiStripGeomFTD::getSensorPitch - posZ: " << std::setw(4) << posZ << " out of range!!!"
                         << std::resetiosflags(std::ios::showpos) << std::setprecision(0) << std::endl;
    exit(0);
  }

  // Changing reference system: rotate the petal in its centroid point
  // to obtain the strips with the stereo angle
  const double tanPhi = tan(_ftdLayer->getPhiHalfDistance(diskID));
  CLHEP::Hep3Vector point = transformPointToRotatedLocal(diskID, sensorID, CLHEP::Hep3Vector(0., 0.0, posZ));

  // Now we can deal the petal in the usual way.
  const double width_z = 2.0 * (getSensorLength(diskID) - point.getZ()) * tanPhi;

  double totalwidth = getSensorWidthMax(diskID);

  return ((totalwidth - width_z) / getSensorNStrips(diskID, sensorID));
}

//
// Get the Y position (in the local ref. frame) for the given strip, posZ
// is the point in the local reference system of the petal
//
double SiStripGeomFTD::getStripPosY(const int& diskID, const int& sensorID, const int& stripID,
                                    const double& posZ) const {
  // Can't be stripID = 0
  if (stripID <= 0) {
    streamlog_out(ERROR) << "SoStripGeomFTD::getStripPosY "
                         << "- incoherent stripID=" << stripID << "!!" << std::endl;
    exit(-1);
  }

  // Extract zPos position in the local rotated frame (if proceed)
  const double zPosPrime = transformPointToRotatedLocal(diskID, sensorID, CLHEP::Hep3Vector(0.0, 0.0, posZ)).getZ();
  // Extract origen of the sensitive part
  const double tanPhi = tan(_ftdLayer->getPhiHalfDistance(diskID));
  const double yorigen = (getSensorLength(diskID) - zPosPrime) * tanPhi;

  // Get pitch (da el pitch teniendo en cuenta lo que tiene que tener)
  double sensPitch = getSensorPitch(diskID, sensorID, posZ);
  if (sensPitch < 1e-40) {
    streamlog_out(ERROR) << "SiStripGeomFTD::getStripPosY "
                         << "- division by zero (sensPitch is zero)!!!" << std::endl;
    exit(-1);
  }

  // Calculate position (Center of the strip)
  double posRPhi = yorigen + sensPitch * (stripID + 0.5);
  // Error
  if ((posRPhi < 0.) || (posRPhi > getSensorWidthMax(diskID))) {
    streamlog_out(ERROR) << "SiStripGeomFTD::getStripPosY - position out of sensor!!!" << std::endl;
    exit(-1);
  }

  // Return R-Phi position of given strip in local ref. system
  CLHEP::Hep3Vector pointInSRL =
      transformPointFromRotatedLocal(diskID, sensorID, CLHEP::Hep3Vector(0.0, posRPhi, zPosPrime));

  return pointInSRL.getY();
}

//
// Transforming a given point to the local ref. frame of a petal which is rotated
// around its center an angle stAngle
//

CLHEP::Hep3Vector SiStripGeomFTD::transformPointToRotatedLocal(const int& diskID, const int& sensorID,
                                                               const CLHEP::Hep3Vector& point) const {
  CLHEP::Hep3Vector pointPrime(point);

  const double stAngle = getStereoAngle(diskID, sensorID);

  // Changing reference system: rotate the petal in its centroid point
  // to obtain the strips with the stereo angle
  double ycentre = getSensorWidthMax(diskID) / 2.0;
  const double zcentre = getSensorLength(diskID) / 2.0;

  // The centre of the rotation
  CLHEP::Hep3Vector rotcentre(0.0, ycentre, zcentre);
  // Translating Zpos to the new ref. system
  pointPrime -= rotcentre;

  // Now we have the hit point from the system placed in the centre
  // of the sensor and rotated the stereo angle

  // We need to return to the local ref. frame (defined along this code)
  // placed in a petal 'rotated'
  rotcentre.rotateX(-stAngle);
  // And extract the point in that system
  pointPrime += rotcentre;

  // Rotating the vector to the new coordinate system--->NOOO
  // pointPrime.rotateX(-stAngle);
  // Consistency Check
  /*CLHEP::Hep3Vector Orig(0.,0.,0.);
  Orig = CLHEP::Hep3Vector(0.0,ycentre,zcentre)-rotcentre;
  std::cout << "To-- pointPrime+O =? P :" << pointPrime + Orig << " =? " << point << std::endl;*/

  return pointPrime;
}

//
// Inverse transformation of transformPointToRotatedLocal (see this function)
//

CLHEP::Hep3Vector SiStripGeomFTD::transformPointFromRotatedLocal(const int& diskID, const int& sensorID,
                                                                 const CLHEP::Hep3Vector& pointPrime) const {
  CLHEP::Hep3Vector point(pointPrime);

  const double stAngle = getStereoAngle(diskID, sensorID);

  // center of the petal
  double ycentre = getSensorWidthMax(diskID) / 2.0;
  const double zcentre = getSensorLength(diskID) / 2.0;

  CLHEP::Hep3Vector rotcentre(0.0, ycentre, zcentre);
  rotcentre.rotateX(-stAngle);

  // Positioning in the center of the rotating frame
  point -= rotcentre;

  // And positioning in the origin
  point += CLHEP::Hep3Vector(0.0, ycentre, zcentre);

  /*	// Consistency Check
          CLHEP::Hep3Vector Orig(0.,0.,0.);
          Orig = CLHEP::Hep3Vector(0.0,ycentre,zcentre)-rotcentre;
          std::cout << " pointPrime+O =? P :" << pointPrime + Orig << " =? " << point << std::endl;*/

  return point;
}

// PRINT METHODS

//
// Method printing general Gear parameters
//

void SiStripGeomFTD::printGearParams() const {
  streamlog_out(MESSAGE3) << std::endl
                          << " " << DUNDERL << DBLUE << "Gear parameters:" << ENDCOLOR << " " << std::endl
                          << std::endl;

  // Gear type: BField
  streamlog_out(MESSAGE3) << "  B field [T]:       " << _magField / T << std::endl << std::endl;

  // Gear type: VXD
  if (_gearType == "VXD") {

    // Print general info
    for (int i = 0; i < _numberOfLayers; ++i) {

      streamlog_out(MESSAGE2) << std::endl;
      streamlog_out(MESSAGE2) << "  Layer: " << _layerRealID[i] << std::endl;
      if (_layerType[i] == pixel) {
        streamlog_out(MESSAGE2) << "  LayerType:       " << "pixel" << std::endl;
      } else {
        streamlog_out(MESSAGE2) << "  LayerType:       " << "strip" << std::endl;
      }
      streamlog_out(MESSAGE2) << "  NumberOfLadders: " << _numberOfLadders[i] << std::endl;
      streamlog_out(MESSAGE2) << "  Radius[mm]:      " << _layerRadius[i] / mm << std::endl;
      streamlog_out(MESSAGE2) << "  Width[mm]:       " << _ladderWidth[i] / mm << std::endl;
      streamlog_out(MESSAGE2) << "  Length[mm]:      " << _ladderLength[i] / mm << std::endl;
      streamlog_out(MESSAGE2) << "  Petal SemiAngle: " << _layerHalfPhi[i] << std::endl;
      streamlog_out(MESSAGE2) << "  Phi0:            " << _layerPhi0[i] << std::endl;
      streamlog_out(MESSAGE2) << "  Theta:           " << _layerTheta[i] << std::endl;
      streamlog_out(MESSAGE2) << "  OffsetY[mm]:     " << _ladderOffsetY[i] / mm << std::endl;
      //        streamlog_out(MESSAGE2) << "  OffsetZ[mm]:     " << _ladderOffsetZ[i]/mm << std::endl; OUT//
    }

  }
  // Gear type: unknown - error
  else {
    streamlog_out(ERROR) << "Unknown gear type!" << std::endl;

    exit(0);
  }
}
//
// Method printing sensor Gear parameters (parameters: sensorID)
//
void SiStripGeomFTD::printSensorParams(short int layerID) const {
  // Print sensor parameters
  streamlog_out(MESSAGE2) << "  * Parameters: " << _sensorThick[layerID] / um << "um thick, "
                          << "with " << _sensorPitchInRear[layerID] / um << "um pitch "
                          << "and " << _sensorNStripsInRear[layerID] << " strips in Z"
                          << ", resp. " << _sensorPitchInFront[layerID] / um << "um pitch "
                          << "and " << _sensorNStripsInFront[layerID] << " strips in R-Phi." << std::endl;
}

} // namespace sistrip
