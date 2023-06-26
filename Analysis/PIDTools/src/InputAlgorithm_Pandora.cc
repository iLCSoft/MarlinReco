#include <InputAlgorithm_Pandora.h>
#include <IMPL/TrackImpl.h>

using namespace lcio;

namespace cpid {

  InputAlgorithm_Pandora anInputAlgorithm_Pandora;

  InputAlgorithm_Pandora::InputAlgorithm_Pandora(): InputAlgorithm("Pandora")
  {
    _description = "Returns dE/dx value of PFO track";
  }

  std::vector<std::string> InputAlgorithm_Pandora::init(std::vector<float> inparF, std::vector<std::string> inparS)
  {
//    sloM << "InputAlgorithm_dEdx init float parameters: ";
//    for (unsigned int i=0; i<inparF.size(); ++i) {sloM << inparF[i] << " ";}
//    sloM << std::endl;
//
//    sloM << "InputAlgorithm_dEdx init string parameters: ";
//    for (unsigned int i=0; i<inparS.size(); ++i) {sloM << inparS[i] << " ";}
//    sloM << std::endl;

    print_init(inparF, inparS);

    std::vector<std::string> obsNames{"PandoraPDG"};
    return obsNames;
  }

  std::vector<std::pair<float,float> > InputAlgorithm_Pandora::extractObservables(ReconstructedParticleImpl* pfo, LCCollection* , int)
  {
    std::vector<std::pair<float,float> > obsValues;

    ReconstructedParticleImpl* mypfo = dynamic_cast<ReconstructedParticleImpl*>(pfo);
    obsValues.push_back(std::pair<float,float>{abs(mypfo->getType()),0});

    return obsValues;

  }

}
