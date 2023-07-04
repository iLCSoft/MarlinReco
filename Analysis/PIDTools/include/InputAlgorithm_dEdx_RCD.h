#ifndef InputAlgorithm_dEdx_RCD_h
#define InputAlgorithm_dEdx_RCD_h 1

#include <InputAlgorithm.h>

using namespace lcio;
namespace cpid {

  class InputAlgorithm_dEdx_RCD : public InputAlgorithm {

  public:

    InputAlgorithm_dEdx_RCD(const InputAlgorithm_dEdx_RCD&) = delete;

    InputAlgorithm_dEdx_RCD& operator=(const InputAlgorithm_dEdx_RCD&) = delete;

    virtual ~InputAlgorithm_dEdx_RCD() = default;

    InputAlgorithm_dEdx_RCD();

    virtual InputAlgorithm* newAlgorithm() {return new InputAlgorithm_dEdx_RCD;}

    virtual std::vector<std::string> init(const std::vector<float>& inparF, const std::vector<std::string>& inparS);

    // The work horse
    virtual std::vector<std::pair<float,float> > extractObservables(ReconstructedParticleImpl* pfo, LCCollection* pfo_col, int PDG);

    double BB_curve(double mass, double mom, std::vector<double>& pars);

  private:

    std::vector<std::vector<double> > _RC_parameters{};
    float _scaling=1;
  };

} // end namespace cpid

#endif
