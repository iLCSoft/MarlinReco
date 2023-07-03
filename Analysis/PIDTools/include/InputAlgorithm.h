#ifndef InputAlgorithm_h
#define InputAlgorithm_h 1

#include <string>
#include <vector>
#include <lcio.h>
#include <streamlog/streamlog.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/LCCollection.h>

#define sloD streamlog_out(DEBUG)
#define sloM streamlog_out(MESSAGE)
#define sloE streamlog_out(ERROR)

using namespace lcio;
namespace cpid {


  class InputAlgorithm {

  public:

    InputAlgorithm(const InputAlgorithm&) = delete;

    InputAlgorithm& operator=(const InputAlgorithm&) = delete;

    virtual ~InputAlgorithm() = default;

    InputAlgorithm(const std::string& typeName);

    virtual InputAlgorithm* newAlgorithm() = 0;

    const std::string & type() const {return _typeName;}

    const std::string & name() const {return _algoName;}

    void setName(std::string name) {_algoName = name;}

    virtual std::vector<std::string> init(const std::vector<float>&, const std::vector<std::string>&);

    void print_init(const std::vector<float>& inparF, const std::vector<std::string>& inparS)
    {
      sloM << _algoName << " init float parameters: ";
      for (unsigned int i=0; i<inparF.size(); ++i) {sloM << inparF[i] << " ";}
      sloM << std::endl;
      sloM << _algoName << " init string parameters: ";
      for (unsigned int i=0; i<inparS.size(); ++i) {sloM << inparS[i] << " ";}
      sloM << std::endl;
    }

    // The work horse
    virtual std::vector<std::pair<float,float> > extractObservables(ReconstructedParticleImpl* pfo, LCCollection* pfo_col, int PDG);

  protected:

    std::string _typeName{};
    std::string _algoName{};
    std::string _description{};

  };

  using InputAlgorithmPtr = std::unique_ptr<InputAlgorithm>;

} // end namespace cpid

#endif 
