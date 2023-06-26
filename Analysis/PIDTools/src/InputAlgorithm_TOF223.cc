#include <InputAlgorithm_TOF223.h>
#include <UTIL/PIDHandler.h>
#include <math.h>

using namespace lcio;

namespace cpid {

  InputAlgorithm_TOF223 anInputAlgorithm_TOF223;

  InputAlgorithm_TOF223::InputAlgorithm_TOF223() : InputAlgorithm("TOF223")
  {
    _description = "Returns TOF mass of PFO";
  }

  std::vector<std::string> InputAlgorithm_TOF223::init(std::vector<float> inparF, std::vector<std::string> inparS){

//    sloM << "InputAlgorithm_TOF223 init float parameters: ";
//    for (unsigned int i=0; i<inparF.size(); ++i) {sloM << inparF[i] << " ";}
//    sloM << std::endl;
//
//    sloM << "InputAlgorithm_TOF223 init string parameters: ";
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

    std::vector<std::string> obsNames{"mass"};
    return obsNames;
  }

  std::vector<std::pair<float,float> > InputAlgorithm_TOF223::extractObservables( ReconstructedParticleImpl* pfo, LCCollection* col_pfo, int){

    //sloM << " IA_TOF::exObs called " << std::endl;

    std::vector<std::pair<float,float> > obsValues;

    PIDHandler PIDHan(col_pfo);
    //sloM << " IA_TOF::exObs - PID handled " << std::endl;

    int TOFAlgoID  = PIDHan.getAlgorithmID(_TOFAlgoName);
    //sloM << " IA_TOF::exObs - TOF Algo ID read" << std::endl;
    int TOFParaID_time   = PIDHan.getParameterIndex(TOFAlgoID,"timeOfFlight");
    int TOFParaID_length = PIDHan.getParameterIndex(TOFAlgoID,"trackLength");
    int TOFParaID_mom    = PIDHan.getParameterIndex(TOFAlgoID,"momentumHM");
    //sloM << " IA_TOF::exObs - TOF IDs set " << std::endl;

    //ReconstructedParticleImpl* mypfo = dynamic_cast<ReconstructedParticleImpl*>(pfo);

    const ParticleID& TOFPID = PIDHan.getParticleID(pfo,TOFAlgoID);
    const FloatVec& TOFParams = TOFPID.getParameters();
    //sloM << " IA_TOF::exObs - PID read " << std::endl;

    if (int(TOFParams.size())>TOFParaID_length)
    {
      float TOFlength = TOFParams[TOFParaID_length];
      float TOFbeta = (TOFlength/TOFParams[TOFParaID_time]) / 299.8;
      float TOFmom = TOFParams[TOFParaID_mom];
      float TOFmass = (1/TOFbeta/TOFbeta>1) ? sqrt(TOFmom*TOFmom*(1/TOFbeta/TOFbeta-1)) : -1;
      if (isnan(TOFmass) || isinf(TOFmass)) TOFmass = -1;
      obsValues.push_back(std::pair<float,float>{TOFmass, 0});
    }
    else
      obsValues.push_back(std::pair<float,float>{-1,-1} );

    //sloM << " IA_TOF::exObs obsValues: " << obsValues[0].first << " " << obsValues[0].second << std::endl;

    return obsValues;

  }

}
