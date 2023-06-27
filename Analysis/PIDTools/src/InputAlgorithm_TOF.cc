#include <InputAlgorithm_TOF.h>
#include <UTIL/PIDHandler.h>
#include <math.h>

using namespace lcio;

namespace cpid {

  InputAlgorithm_TOF anInputAlgorithm_TOF;

  InputAlgorithm_TOF::InputAlgorithm_TOF() : InputAlgorithm("TOF")
  {
    _description = "Returns TOF beta of PFO";
  }

  std::vector<std::string> InputAlgorithm_TOF::init(const std::vector<float>& inparF, const std::vector<std::string>& inparS)
  {
    print_init(inparF, inparS);

    if (inparS.size()>0)
      _TOFAlgoName = inparS[0];
    else
    {
      sloE << _algoName << ": too few (<1) string parameters! TOFAlgoName not specified in TOF::init." << std::endl;
      throw std::runtime_error("parameters error");
    }

    std::vector<std::string> obsNames{"beta"};
    return obsNames;
  }

  std::vector<std::pair<float,float> > InputAlgorithm_TOF::extractObservables( ReconstructedParticleImpl* pfo, LCCollection* col_pfo, int)
  {
    std::vector<std::pair<float,float> > obsValues;

    PIDHandler PIDHan(col_pfo);

    int TOFAlgoID  = PIDHan.getAlgorithmID(_TOFAlgoName);
    int TOFParaID_closest  = PIDHan.getParameterIndex(TOFAlgoID,"TOFClosestHits");
    int TOFParaID_length   = PIDHan.getParameterIndex(TOFAlgoID,"TOFFlightLength");

    const ParticleID& TOFPID = PIDHan.getParticleID(pfo,TOFAlgoID);
    const FloatVec& TOFParams = TOFPID.getParameters();

    if (int(TOFParams.size())>TOFParaID_length)
    {
      float TOFlength = TOFParams[TOFParaID_length];
      float TOFbeta_ch = (TOFlength/TOFParams[TOFParaID_closest]) / 299.8;
      if (isnan(TOFbeta_ch) || isinf(TOFbeta_ch)) TOFbeta_ch = -1;
      obsValues.push_back(std::pair<float,float>{TOFbeta_ch, 0});
    }
    else
      obsValues.push_back(std::pair<float,float>{-1,-1} );

    return obsValues;
  }

}
