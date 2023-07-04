#include <InputAlgorithm_dEdx.h>
#include <IMPL/TrackImpl.h>

using namespace lcio;

namespace cpid {

  InputAlgorithm_dEdx anInputAlgorithm_dEdx;

  InputAlgorithm_dEdx::InputAlgorithm_dEdx(): InputAlgorithm("dEdx")
  {
    _description = "Returns dE/dx value of PFO track";
  }

  std::vector<std::string> InputAlgorithm_dEdx::init(const std::vector<float>& inparF, const std::vector<std::string>& inparS)
  {
    print_init(inparF, inparS);

    std::vector<std::string> obsNames{"dEdx"};
    return obsNames;
  }

  std::vector<std::pair<float,float> > InputAlgorithm_dEdx::extractObservables(ReconstructedParticleImpl* pfo, LCCollection* , int)
  {
    std::vector<std::pair<float,float> > obsValues;

    TrackVec tracks = pfo->getTracks();
    if (tracks.size() > 0)
    {
      Track* track = tracks[0];
      std::pair<float,float> dEdxValue {track->getdEdx(),track->getdEdxError()};
      obsValues.push_back(dEdxValue);
    }
    else
      obsValues.push_back(std::pair<float,float>{-1,-1});

    return obsValues;

  }

}
