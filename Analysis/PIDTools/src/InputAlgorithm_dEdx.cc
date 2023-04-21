#include <InputAlgorithm_dEdx.h>
#include <IMPL/TrackImpl.h>

using namespace lcio;

namespace cpid {

  InputAlgorithm_dEdx anInputAlgorithm_dEdx;

  InputAlgorithm_dEdx::InputAlgorithm_dEdx(): InputAlgorithm("dEdx")
  {
    _description = "Returns dE/dx value of PFO track";
  }

  std::vector<std::string> InputAlgorithm_dEdx::init(std::vector<float> inparF, std::vector<std::string> inparS)
  {
//    sloM << "InputAlgorithm_dEdx init float parameters: ";
//    for (unsigned int i=0; i<inparF.size(); ++i) {sloM << inparF[i] << " ";}
//    sloM << std::endl;
//
//    sloM << "InputAlgorithm_dEdx init string parameters: ";
//    for (unsigned int i=0; i<inparS.size(); ++i) {sloM << inparS[i] << " ";}
//    sloM << std::endl;

    print_init(inparF, inparS);

    std::vector<std::string> obsNames{"dEdx"};
    return obsNames;
  }

  std::vector<std::pair<float,float> > InputAlgorithm_dEdx::extractObservables(ReconstructedParticleImpl* pfo, LCCollection* )
  {
    std::vector<std::pair<float,float> > obsValues;

    //ReconstructedParticleImpl* mypfo = dynamic_cast<ReconstructedParticleImpl*>(pfo);
    TrackVec tracks = pfo->getTracks();
    if (tracks.size() > 0)
    {
      Track* track = tracks[0];
      std::pair<float,float> dEdxValue {track->getdEdx(),track->getdEdxError()};
      obsValues.push_back(dEdxValue);
    }
    else
      obsValues.push_back(std::pair<float,float>{-1,-1});

    //sloM << " IA_dEdx::exObs obsValues: " << obsValues[0].first << " " << obsValues[0].second << std::endl;

    return obsValues;

  }

}
