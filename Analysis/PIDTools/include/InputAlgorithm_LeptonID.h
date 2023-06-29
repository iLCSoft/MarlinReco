#ifndef InputAlgorithm_LeptonID_h
#define InputAlgorithm_LeptonID_h 1

#include <InputAlgorithm.h>

using namespace lcio;
namespace cpid {

  class InputAlgorithm_LeptonID : public InputAlgorithm {

  public:

    InputAlgorithm_LeptonID(const InputAlgorithm_LeptonID&) = delete;

    InputAlgorithm_LeptonID& operator=(const InputAlgorithm_LeptonID&) = delete;

    virtual ~InputAlgorithm_LeptonID() = default;

    InputAlgorithm_LeptonID();

    virtual InputAlgorithm* newAlgorithm() {return new InputAlgorithm_LeptonID;}

    virtual std::vector<std::string> init(const std::vector<float>& inparF, const std::vector<std::string>& inparS);
  
    // The work horse
    virtual std::vector<std::pair<float,float> > extractObservables(ReconstructedParticleImpl* pfo, LCCollection* pfo_col, int PDG);

  private:

    std::string _LeptonIDAlgoName{};

  };

} // end namespace cpid

#endif 
