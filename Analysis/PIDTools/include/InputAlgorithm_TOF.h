#ifndef InputAlgorithm_TOF_h
#define InputAlgorithm_TOF_h 1

#include <InputAlgorithm.h>

using namespace lcio;
namespace cpid {

class InputAlgorithm_TOF : public InputAlgorithm {

public:
  InputAlgorithm_TOF(const InputAlgorithm_TOF&) = delete;

  InputAlgorithm_TOF& operator=(const InputAlgorithm_TOF&) = delete;

  virtual ~InputAlgorithm_TOF() = default;

  InputAlgorithm_TOF();

  virtual InputAlgorithm* newAlgorithm() { return new InputAlgorithm_TOF; }

  virtual std::vector<std::string> init(const std::vector<float>& inparF, const std::vector<std::string>& inparS);

  // The work horse
  virtual std::vector<std::pair<float, float>> extractObservables(ReconstructedParticleImpl* pfo, LCCollection* pfo_col,
                                                                  int PDG);

private:
  std::string _TOFAlgoName{};
};

} // end namespace cpid

#endif
