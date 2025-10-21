#include <AlgorithmMgr.h>
#include <InputAlgorithm.h>
// #include <InputAlgorithm_dEdx.hh>

#include <lcio.h>

namespace cpid {

InputAlgorithm::InputAlgorithm(const std::string& typeName) {
  _description = "description not set by author";
  _typeName = typeName;

  // register processor in map
  AlgorithmMgr::instance()->registerAlgorithm(this);
}

std::vector<std::string> InputAlgorithm::init(const std::vector<float>&, const std::vector<std::string>&) {
  std::vector<std::string> obsNames{};
  return obsNames;
}

std::vector<std::pair<float, float>> InputAlgorithm::extractObservables(ReconstructedParticleImpl*, LCCollection*,
                                                                        int) {
  std::vector<std::pair<float, float>> obsValues{};
  std::cout << " IA::exObs called " << std::endl;

  // if (_algoName == "dEdx") obsValues = cpid::extract_dEdx(pfo);

  std::cout << " IA::exObs returning... " << std::endl;
  return obsValues;
}

} // namespace cpid
