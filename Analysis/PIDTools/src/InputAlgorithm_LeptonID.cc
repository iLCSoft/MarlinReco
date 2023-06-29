#include <InputAlgorithm_LeptonID.h>
#include <UTIL/PIDHandler.h>

using namespace lcio;

namespace cpid {

  InputAlgorithm_LeptonID InputAlgorithm_LeptonID;

  InputAlgorithm_LeptonID::InputAlgorithm_LeptonID(): InputAlgorithm("LeptonID")
  {
    _description = "Returns LeptonID values of PFO";
  }

  std::vector<std::string> InputAlgorithm_LeptonID::init(const std::vector<float>& inparF, const std::vector<std::string>& inparS)
  {
    print_init(inparF, inparS);

    if (inparS.size()>0) _LeptonIDAlgoName = inparS[0];
    else _LeptonIDAlgoName = "LeptonID";

    std::vector<std::string> obsNames{"elness","muness","piness"};
    return obsNames;
  }

  std::vector<std::pair<float,float> > InputAlgorithm_LeptonID::extractObservables(ReconstructedParticleImpl* pfo, LCCollection* col_pfo, int)
  {
    std::vector<std::pair<float,float> > obsValues;
    for (int i=0; i<3; ++i) obsValues.push_back(std::pair<float,float>{-1,-1});

    PIDHandler PIDHan(col_pfo);

    int LeptonIDAlgoID  = PIDHan.getAlgorithmID(_LeptonIDAlgoName);
    const ParticleID& LeptonID = PIDHan.getParticleID(pfo,LeptonIDAlgoID);
    const FloatVec& LeptonIDParams = LeptonID.getParameters();

    for (unsigned int i=0; i<LeptonIDParams.size()&&i<3; ++i) obsValues[i]=std::pair<float,float>{LeptonIDParams[i], 0};

    return obsValues;

  }

}
