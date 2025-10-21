#include "ChargeInducer.h"
#include "SimDigital.h"

#include <sstream>

#include <marlin/VerbosityLevels.h>

#include <TFile.h>
#include <TTree.h>

ChargeInducer::ChargeInducer() : generator() {}

ChargeInducer::~ChargeInducer() {}
void ChargeInducer::setSeed(unsigned int value) { generator.seed(value); }

UniformPolya::UniformPolya(float _qbar, float _theta)
    : ChargeInducer(), gammadist(std::gamma_distribution<float>(_qbar / _theta, _theta)) {}

UniformPolya::~UniformPolya() {}

float UniformPolya::getCharge(SimDigitalGeomCellId*) { return gammadist(generator); }

AsicPolya::AsicPolya(float _qbar, float _theta, std::string fileName) : UniformPolya(_qbar, _theta), polyaMap() {
  readFile(fileName);
}

AsicPolya::~AsicPolya() {}

void AsicPolya::readFile(std::string fileName) {
  TFile* file = TFile::Open(fileName.c_str(), "READ");
  if (!file) {
    std::cerr << "ERROR : file " << fileName << " not found for AsicPolya::readFile" << std::endl;
    throw;
  }

  TTree* tree = dynamic_cast<TTree*>(file->Get("tree"));
  if (!tree) {
    std::cerr << "ERROR : tree not present in file " << fileName << std::endl;
    throw;
  }

  int asicID;
  int layerID;
  float qbarAsic;
  float deltaAsic;
  std::vector<double>* position = NULL;

  tree->SetBranchAddress("LayerID", &layerID);
  tree->SetBranchAddress("AsicID", &asicID);
  tree->SetBranchAddress("qbar", &qbarAsic);
  tree->SetBranchAddress("delta", &deltaAsic);
  tree->SetBranchAddress("Position", &position);

  int iEntry = 0;
  while (tree->GetEntry(iEntry++)) {
    if (qbarAsic <= 0 || deltaAsic <= 0)
      continue;

    float alpha = qbarAsic / deltaAsic;
    float delta = deltaAsic;

    std::gamma_distribution<float> localDist(alpha, delta);

    int iAsic = static_cast<int>((position->at(0) - 10.408) / (8 * 10.408));
    int jAsic = static_cast<int>((position->at(1) - 10.408) / (8 * 10.408));
    int K = static_cast<int>((position->at(2) - 26.131) / 26.131 + 0.5);

    if (asicID == -1 && layerID != -1) // global value for layer
      polyaMap.insert(std::make_pair(AsicKey(K), localDist));
    else
      polyaMap.insert(std::make_pair(AsicKey(K, iAsic, jAsic), localDist));
  }
  file->Close();
}

float AsicPolya::getCharge(SimDigitalGeomCellId* cellID) {
  //	int asicKey = (cellID.I()-1)/8 + ((cellID.J()-1)/8)*12 + cellID.K()*1000 ;
  AsicKey asicKey(cellID->K(), (cellID->I() - 1) / 8, (cellID->J() - 1) / 8);

  std::map<AsicKey, std::gamma_distribution<float>>::iterator it = polyaMap.find(asicKey);

  if (it == polyaMap.end()) {
    it = polyaMap.find(AsicKey(cellID->K())); // else search for layer polya
    if (it == polyaMap.end())
      return gammadist(generator);
    else
      return it->second(generator);
  } else
    return it->second(generator);
}
