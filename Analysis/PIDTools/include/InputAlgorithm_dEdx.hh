//#include <IMPL/ReconstructedParticleImpl.h>
//#include <IMPL/TrackImpl.h>
#include <InputAlgorithm.h>

namespace cpid {

  class InputAlgorithm_dEdx : InputAlgorithm {

//  public:
//
//    InputAlgorithm(const InputAlgorithm&) = delete;
//
//    InputAlgorithm& operator=(const InputAlgorithm&) = delete;
//
//    virtual ~InputAlgorithm() = default;
//
//    InputAlgorithm();
//
//    virtual std::vector<std::string> init( std::vector<float> inparF, std::vector<std::string> inparS, std::string obs);
//
//    // The work horse
//    virtual std::vector<std::pair<double,double> > extractObservables( ReconstructedParticleImpl* pfo );
//
//  private:
//
//    std::string _algoName{};
//
  };


//  std::vector<std::string> init_dEdx(std::vector<float> inparF, std::vector<std::string> inparS){
//
//    std::vector<std::string> obsNames{"dEdx"};
//    std::cout << "init: " << obsNames[0] << std::endl;
//    return obsNames;
//
//  }
//
//  std::vector<std::pair<double,double> > extract_dEdx( ReconstructedParticleImpl* pfo ){
//
//    std::vector<std::pair<double,double> > obsValues;
//    obsValues.push_back(std::pair<double,double>(0,0));
//    std::cout << " IA:ex_dEdx called " << std::endl;
//    //ReconstructedParticleImpl* mypfo = dynamic_cast<ReconstructedParticleImpl*>(pfo);
//    //std::cout << " dc mypfo done " << std::endl;
//    if (pfo->getTracks().size()==0) return obsValues;
//    Track* track = pfo->getTracks()[0];
//    //std::cout << " track received " << std::endl;
//
//    //std::vector<std::pair<double,double> > obsValues;
//    std::pair<double,double> dEdxValue = {track->getdEdx(),track->getdEdxError()};
//    //dEdxValue.first = track->getdEdx();
//    //dEdxValue.second = track->getdEdxError();
//    //std::cout << " dEdx values found " << std::endl;
//    obsValues[0] = dEdxValue;
//    //std::cout << " dEdx values stored " << std::endl;
//    return obsValues;
//
//  }

}
