#include <ModelMgr.h>
#include <TrainingModel.h>
// #include <InputAlgorithm_dEdx.hh>

#include <lcio.h>

namespace cpid {

TrainingModel::TrainingModel(const std::string& typeName) {
  _description = "description not set by author";
  _typeName = typeName;

  ModelMgr::instance()->registerModel(this);
}

std::vector<std::string> TrainingModel::initTraining(TrainingModelInterface&) {
  std::vector<std::string> weightfiles{};
  return weightfiles;
}

void TrainingModel::initInference(TrainingModelInterface&) {}

void TrainingModel::runTraining(TTree*) {}

const std::vector<float> TrainingModel::runInference(int) {
  const std::vector<float> result{};
  return result;
}

} // namespace cpid
