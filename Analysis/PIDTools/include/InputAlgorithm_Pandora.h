#ifndef InputAlgorithm_Pandora_h
#define InputAlgorithm_Pandora_h 1

#include <InputAlgorithm.h>

using namespace lcio;
namespace cpid {

class InputAlgorithm_Pandora : public InputAlgorithm {

public:
  InputAlgorithm_Pandora(const InputAlgorithm_Pandora&) = delete;

  InputAlgorithm_Pandora& operator=(const InputAlgorithm_Pandora&) = delete;

  virtual ~InputAlgorithm_Pandora() = default;

  InputAlgorithm_Pandora();

  virtual InputAlgorithm* newAlgorithm() { return new InputAlgorithm_Pandora; }

  virtual std::vector<std::string> init(const std::vector<float>& inparF, const std::vector<std::string>& inparS);

  // The work horse
  virtual std::vector<std::pair<float, float>> extractObservables(ReconstructedParticleImpl* pfo, LCCollection* pfo_col,
                                                                  int PDG);

private:
};

} // end namespace cpid

#endif
