#include "EfficiencyManager.h"
#include "SimDigital.h"

#include <marlin/VerbosityLevels.h>

#include <TFile.h>
#include <TTree.h>

EfficiencyManager::EfficiencyManager() {}
EfficiencyManager::~EfficiencyManager() {}

UniformEfficiency::UniformEfficiency(float val) : EfficiencyManager(), value(val) {}
UniformEfficiency::~UniformEfficiency() {}

float UniformEfficiency::getEfficiency(SimDigitalGeomCellId*) { return value; }

AsicEfficiency::AsicEfficiency(std::string fileName, float globalVal) : UniformEfficiency(globalVal), effMap() {
  readFile(fileName);
}

AsicEfficiency::~AsicEfficiency() {}

void AsicEfficiency::readFile(std::string fileName) {
  TFile* file = TFile::Open(fileName.c_str(), "READ");
  if (!file) {
    std::cerr << "ERROR : file " << fileName << " not found for AsicEfficiency::readFile" << std::endl;
    throw;
  }

  TTree* tree = dynamic_cast<TTree*>(file->Get("tree"));
  if (!tree) {
    std::cerr << "ERROR : tree not present in file " << fileName << std::endl;
    throw;
  }

  int padID;
  int asicID;
  int layerID;
  std::vector<double>* efficiencies = NULL;
  std::vector<double>* position = NULL;

  tree->SetBranchAddress("LayerID", &layerID);
  tree->SetBranchAddress("AsicID", &asicID);
  tree->SetBranchAddress("PadID", &padID);
  tree->SetBranchAddress("Efficiencies", &efficiencies);
  tree->SetBranchAddress("Position", &position);

  int iEntry = 0;
  while (tree->GetEntry(iEntry++)) {
    if (padID != -1)
      continue;

    int iAsic = static_cast<int>((position->at(0) - 10.408) / (8 * 10.408));
    int jAsic = static_cast<int>((position->at(1) - 10.408) / (8 * 10.408));
    int K = static_cast<int>((position->at(2) - 26.131) / 26.131 + 0.5);

    if (asicID == -1 && layerID != -1) // global value for layer
      effMap.insert(std::make_pair(AsicKey(K), efficiencies->at(0)));
    else
      effMap.insert(std::make_pair(AsicKey(K, iAsic, jAsic), efficiencies->at(0)));
  }

  file->Close();
}

float AsicEfficiency::getEfficiency(SimDigitalGeomCellId* cellID) {
  AsicKey asicKey(cellID->K(), (cellID->I() - 1) / 8, (cellID->J() - 1) / 8);

  std::map<AsicKey, float>::const_iterator it = effMap.find(asicKey);

  if (it == effMap.end()) {
    it = effMap.find(AsicKey(cellID->K())); // else search for layer mul
    if (it == effMap.end())
      return value;
    else
      return it->second;
  } else
    return it->second;
}
