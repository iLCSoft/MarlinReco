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

  std::vector<std::string> InputAlgorithm_TOF::init(std::vector<float> inparF, std::vector<std::string> inparS){

//    sloM << "InputAlgorithm_TOF init float parameters: ";
//    for (unsigned int i=0; i<inparF.size(); ++i) {sloM << inparF[i] << " ";}
//    sloM << std::endl;
//
//    sloM << "InputAlgorithm_TOF init string parameters: ";
//    for (unsigned int i=0; i<inparS.size(); ++i) {sloM << inparS[i] << " ";}
//    sloM << std::endl;

    print_init(inparF, inparS);

    if (inparS.size()>0)
      _TOFAlgoName = inparS[0];
    else
    {
      sloE << "TOFAlgoName not specified in TOF::init!" << std::endl;
      //_TOFAlgoName = inparS[0];
    }

    std::vector<std::string> obsNames{"beta"};
    return obsNames;
  }

  std::vector<std::pair<float,float> > InputAlgorithm_TOF::extractObservables( ReconstructedParticleImpl* pfo, LCCollection* col_pfo, int){

    //sloM << " IA_TOF::exObs called " << std::endl;

    std::vector<std::pair<float,float> > obsValues;

    PIDHandler PIDHan(col_pfo);
    //sloM << " IA_TOF::exObs - PID handled " << std::endl;

    int TOFAlgoID  = PIDHan.getAlgorithmID(_TOFAlgoName);
    //sloM << " IA_TOF::exObs - TOF Algo ID read" << std::endl;
    int TOFParaID_closest  = PIDHan.getParameterIndex(TOFAlgoID,"TOFClosestHits");
    int TOFParaID_length   = PIDHan.getParameterIndex(TOFAlgoID,"TOFFlightLength");
    //sloM << " IA_TOF::exObs - TOF IDs set " << std::endl;

    //ReconstructedParticleImpl* mypfo = dynamic_cast<ReconstructedParticleImpl*>(pfo);

    const ParticleID& TOFPID = PIDHan.getParticleID(pfo,TOFAlgoID);
    const FloatVec& TOFParams = TOFPID.getParameters();
    //sloM << " IA_TOF::exObs - PID read " << std::endl;

    if (int(TOFParams.size())>TOFParaID_length)
    {
      float TOFlength = TOFParams[TOFParaID_length];
      float TOFbeta_ch = (TOFlength/TOFParams[TOFParaID_closest]) / 299.8;
      if (isnan(TOFbeta_ch) || isinf(TOFbeta_ch)) TOFbeta_ch = -1;
      obsValues.push_back(std::pair<float,float>{TOFbeta_ch, 0});
    }
    else
      obsValues.push_back(std::pair<float,float>{-1,-1} );

    //sloM << " IA_TOF::exObs obsValues: " << obsValues[0].first << " " << obsValues[0].second << std::endl;

    return obsValues;

  }

}
