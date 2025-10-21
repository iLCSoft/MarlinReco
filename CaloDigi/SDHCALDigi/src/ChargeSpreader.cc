#include "ChargeSpreader.h"
#include "SimDigital.h"

#include <cmath>

#include <TFile.h>
#include <TTree.h>

#include <marlin/Exceptions.h>
#include <marlin/VerbosityLevels.h>

using namespace marlin;

ChargeSpreader::ChargeSpreader() : chargeMap(), parameters() {}

ChargeSpreader::~ChargeSpreader() {}

void ChargeSpreader::addCharge(float charge, float posI, float posJ, SimDigitalGeomCellId*) {
  if (parameters.padSeparation > parameters.cellSize)
    return;

  float chargeTotCheck = 0;

  int minCellI = static_cast<int>(std::round(posI - parameters.range) / parameters.cellSize);
  int maxCellI = static_cast<int>(std::round(posI + parameters.range) / parameters.cellSize);

  int minCellJ = static_cast<int>(std::round(posJ - parameters.range) / parameters.cellSize);
  int maxCellJ = static_cast<int>(std::round(posJ + parameters.range) / parameters.cellSize);

  for (int I = minCellI; I <= maxCellI; ++I) {
    float minI = (I - 0.5f) * parameters.cellSize - posI + 0.5f * parameters.padSeparation;
    float maxI = (I + 0.5f) * parameters.cellSize - posI - 0.5f * parameters.padSeparation;

    for (int J = minCellJ; J <= maxCellJ; ++J) {
      float minJ = (J - 0.5f) * parameters.cellSize - posJ + 0.5f * parameters.padSeparation;
      float maxJ = (J + 0.5f) * parameters.cellSize - posJ - 0.5f * parameters.padSeparation;

      float integralResult = computeIntegral(minI, maxI, minJ, maxJ);

      chargeMap[I_J_Coordinates(I, J)] += charge * integralResult * normalisation;

      if (chargeMap[I_J_Coordinates(I, J)] < 0)
        streamlog_out(MESSAGE) << "!!!!!!!!!!Negative Charge!!!!!!!!!!" << std::endl
                               << " X " << posJ << " " << minJ << " " << maxJ << std::endl
                               << " Y " << posI << " " << minI << " " << maxI << std::endl;

      chargeTotCheck += chargeMap[I_J_Coordinates(I, J)];
    }
  }
  streamlog_out(DEBUG) << " Charge = " << charge << " ; total splitted charge = " << chargeTotCheck << std::endl;
}

GaussianSpreader::GaussianSpreader() : ChargeSpreader() { normalisation = 1.0f; }

GaussianSpreader::~GaussianSpreader() {}

void GaussianSpreader::init() {
  assert(parameters.erfWidth.size() == parameters.erfWeigth.size());

  float _normalisation = 0;
  for (unsigned int i = 0; i < parameters.erfWidth.size(); i++) {
    streamlog_out(DEBUG) << "Erf function parameters " << i + 1 << " : " << parameters.erfWidth[i] << ", "
                         << parameters.erfWeigth[i] << std::endl;
    _normalisation += parameters.erfWeigth.at(i) * parameters.erfWidth.at(i) * parameters.erfWidth.at(i) * M_PI;
  }

  normalisation = 1.f / _normalisation;

  streamlog_out(DEBUG) << "Charge splitter normalisation factor: " << normalisation << std::endl;
  streamlog_out(DEBUG) << "range : " << parameters.range << " ; padseparation : " << parameters.padSeparation
                       << std::endl;
}

float GaussianSpreader::computeIntegral(float x1, float x2, float y1, float y2) const {
  float integralResult = 0.0f;

  for (unsigned int n = 0; n < parameters.erfWidth.size(); n++) {
    integralResult += fabs(std::erf(x2 / parameters.erfWidth.at(n)) - std::erf(x1 / parameters.erfWidth.at(n))) *
                      fabs(std::erf(y2 / parameters.erfWidth.at(n)) - std::erf(y1 / parameters.erfWidth.at(n))) *
                      parameters.erfWeigth.at(n) * M_PI * parameters.erfWidth.at(n) * parameters.erfWidth.at(n) / 4;
  }
  return integralResult;
}

ExactSpreader::ExactSpreader() : ChargeSpreader() {}

ExactSpreader::~ExactSpreader() {}

void ExactSpreader::init() {
  float _normalisation = static_cast<float>(2 * M_PI);
  normalisation = 1.0f / _normalisation;

  streamlog_out(DEBUG) << "Charge splitter normalisation factor: " << normalisation << std::endl;
  streamlog_out(DEBUG) << "range : " << parameters.range << " ; padseparation : " << parameters.padSeparation
                       << std::endl;
}

float ExactSpreader::computeIntegral(float x1, float x2, float y1, float y2) const {
  float term1 = std::atan(y2 * x2 / (parameters.d * std::sqrt(parameters.d * parameters.d + y2 * y2 + x2 * x2)));
  float term2 = std::atan(y1 * x2 / (parameters.d * std::sqrt(parameters.d * parameters.d + y1 * y1 + x2 * x2)));
  float term3 = std::atan(y2 * x1 / (parameters.d * std::sqrt(parameters.d * parameters.d + y2 * y2 + x1 * x1)));
  float term4 = std::atan(y1 * x1 / (parameters.d * std::sqrt(parameters.d * parameters.d + y1 * y1 + x1 * x1)));

  return (term1 - term2 - term3 + term4);
}

ExactSpreaderPerAsic::ExactSpreaderPerAsic(std::string fileName) : ExactSpreader(), dMap() { readFile(fileName); }

ExactSpreaderPerAsic::~ExactSpreaderPerAsic() {}

void ExactSpreaderPerAsic::readFile(std::string fileName) {
  TFile* file = TFile::Open(fileName.c_str(), "READ");
  if (!file) {
    std::cerr << "ERROR : file " << fileName << " not found for MultiplicityChargeSplitterExactPerAsic::readFile"
              << std::endl;
    throw;
  }

  TTree* tree = dynamic_cast<TTree*>(file->Get("tree"));
  if (!tree) {
    std::cerr << "ERROR : tree not present in file " << fileName << std::endl;
    throw;
  }

  int asicID;
  int layerID;
  float dAsic;
  std::vector<double>* position = NULL;

  tree->SetBranchAddress("LayerID", &layerID);
  tree->SetBranchAddress("AsicID", &asicID);
  tree->SetBranchAddress("d", &dAsic);
  tree->SetBranchAddress("Position", &position);

  int iEntry = 0;
  while (tree->GetEntry(iEntry++)) {
    int iAsic = static_cast<int>((position->at(0) - 10.408) / (8 * 10.408));
    int jAsic = static_cast<int>((position->at(1) - 10.408) / (8 * 10.408));
    int K = static_cast<int>((position->at(2) - 26.131) / 26.131 + 0.5);

    if (asicID == -1 && layerID != -1) // global value for layer
      dMap.insert(std::make_pair(AsicKey(K), dAsic));
    else
      dMap.insert(std::make_pair(AsicKey(K, iAsic, jAsic), dAsic));
  }
  file->Close();
}

void ExactSpreaderPerAsic::addCharge(float charge, float posI, float posJ, SimDigitalGeomCellId* cellID) {
  AsicKey asicKey(cellID->K(), (cellID->I() - 1) / 8, (cellID->J() - 1) / 8);

  std::map<AsicKey, float>::iterator it = dMap.find(asicKey);

  if (it == dMap.end()) {
    it = dMap.find(AsicKey(cellID->K())); // else search for layer mul
    if (it == dMap.end())
      parameters.d = dGlobal;
    else
      parameters.d = it->second;
  } else
    parameters.d = it->second;

  ExactSpreader::addCharge(charge, posI, posJ, cellID);
}
