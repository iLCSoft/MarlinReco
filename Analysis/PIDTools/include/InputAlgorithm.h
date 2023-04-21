#ifndef InputAlgorithm_h
#define InputAlgorithm_h 1

#include <string>
#include <vector>
#include <lcio.h>
#include <streamlog/streamlog.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/LCCollection.h>

#define sloM streamlog_out(MESSAGE)
#define sloE streamlog_out(ERROR)

using namespace lcio;
namespace cpid {

//class LikelihoodPID;
//class LowMomentumMuPiSeparationPID_BDTG;

  class InputAlgorithm {

  public:

    InputAlgorithm(const InputAlgorithm&) = delete;

    InputAlgorithm& operator=(const InputAlgorithm&) = delete;

    virtual ~InputAlgorithm() = default;

    InputAlgorithm(const std::string& typeName);

    virtual InputAlgorithm* newAlgorithm() = 0;

    virtual const std::string & type() const {return _typeName;}

    virtual const std::string & name() const {return _algoName;}

    virtual void setName(std::string name) {_algoName = name;}

    virtual std::vector<std::string> init(std::vector<float>, std::vector<std::string> );

    void print_init(std::vector<float> inparF, std::vector<std::string> inparS)
    {
      sloM << _algoName << " init float parameters: ";
      for (unsigned int i=0; i<inparF.size(); ++i) {sloM << inparF[i] << " ";}
      sloM << std::endl;
      sloM << _algoName << " init string parameters: ";
      for (unsigned int i=0; i<inparS.size(); ++i) {sloM << inparS[i] << " ";}
      sloM << std::endl;
    }

    // The work horse
    virtual std::vector<std::pair<float,float> > extractObservables(ReconstructedParticleImpl* pfo, LCCollection* pfo_col);

  protected:

    std::string _typeName{};
    std::string _algoName{};
    std::string _description{};

  };

} // end namespace cpid

#endif 
