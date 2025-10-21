#ifndef InputAlgorithm_TOF223_h
#define InputAlgorithm_TOF223_h 1

#include <InputAlgorithm.h>

using namespace lcio;
namespace cpid {

class InputAlgorithm_TOF223 : public InputAlgorithm {

public:
  InputAlgorithm_TOF223(const InputAlgorithm_TOF223&) = delete;

  InputAlgorithm_TOF223& operator=(const InputAlgorithm_TOF223&) = delete;

  virtual ~InputAlgorithm_TOF223() = default;

  InputAlgorithm_TOF223();

  virtual InputAlgorithm_TOF223* newAlgorithm() { return new InputAlgorithm_TOF223; }

  virtual std::vector<std::string> init(const std::vector<float>& inparF, const std::vector<std::string>& inparS);

  // The work horse
  virtual std::vector<std::pair<float, float>> extractObservables(ReconstructedParticleImpl* pfo, LCCollection* pfo_col,
                                                                  int PDG);

private:
  std::string _TOFAlgoName{};
};

} // end namespace cpid

#endif
