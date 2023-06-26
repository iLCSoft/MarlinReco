#ifndef InputAlgorithm_dEdx_h
#define InputAlgorithm_dEdx_h 1

#include <InputAlgorithm.h>

using namespace lcio;
namespace cpid {

  class InputAlgorithm_dEdx : public InputAlgorithm {

  public:

    InputAlgorithm_dEdx(const InputAlgorithm_dEdx&) = delete;

    InputAlgorithm_dEdx& operator=(const InputAlgorithm_dEdx&) = delete;

    virtual ~InputAlgorithm_dEdx() = default;

    InputAlgorithm_dEdx();

    virtual InputAlgorithm* newAlgorithm() {return new InputAlgorithm_dEdx;}

    virtual std::vector<std::string> init( std::vector<float> inparF, std::vector<std::string> inparS);
  
    // The work horse
    virtual std::vector<std::pair<float,float> > extractObservables(ReconstructedParticleImpl* pfo, LCCollection* pfo_col, int PDG);

  private:

  };

} // end namespace cpid

#endif 
