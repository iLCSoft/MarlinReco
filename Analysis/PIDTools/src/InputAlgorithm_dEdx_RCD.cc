#include <InputAlgorithm_dEdx_RCD.h>
#include <IMPL/TrackImpl.h>
#include <math.h>

using namespace lcio;

namespace cpid {

  InputAlgorithm_dEdx_RCD anInputAlgorithm_dEdx_RCD;

  InputAlgorithm_dEdx_RCD::InputAlgorithm_dEdx_RCD(): InputAlgorithm("dEdx_RCD")
  {
    _description = "Returns dE/dx value of PFO track";
  }

  std::vector<std::string> InputAlgorithm_dEdx_RCD::init(std::vector<float> inparF, std::vector<std::string> inparS)
  {
//    sloM << "InputAlgorithm_dEdx_RCD init float parameters: ";
//    for (unsigned int i=0; i<inparF.size(); ++i) {sloM << inparF[i] << " ";}
//    sloM << std::endl;
//
//    sloM << "InputAlgorithm_dEdx_RCD init string parameters: ";
//    for (unsigned int i=0; i<inparS.size(); ++i) {sloM << inparS[i] << " ";}
//    sloM << std::endl;

    print_init(inparF, inparS);

    if (inparF.size()<25) {sloE << "InputAlgorithm dEdx_RCD received too few (<25) parameters for reference curves!" << std::endl;}

    auto it = inparF.begin();
    for (int i=0; i<5; ++i) _RC_parameters.push_back(std::vector<double>(it+i*5,it+i*5+5));

    std::vector<std::string> obsNames{"elDis","muDis","piDis","kaDis","prDis"};
    return obsNames;
  }

  std::vector<std::pair<float,float> > InputAlgorithm_dEdx_RCD::extractObservables(ReconstructedParticleImpl* pfo, LCCollection* )
  {
    std::vector<std::pair<float,float> > obsValues;
    for (int i=0; i<5; ++i) obsValues.push_back(std::pair<float,float>{-1000,-1});

    const double* mom = pfo->getMomentum();
    double momabs = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);

    double mass[] = {0.000510998, 0.105658, 0.139570, 0.493677, 0.938272};

    TrackVec tracks = pfo->getTracks();
    if (tracks.size() > 0)
    {
      Track* track = tracks[0];
      std::pair<float,float> dEdx {track->getdEdx(),track->getdEdxError()};
      //sloM << dEdx.first << " " << dEdx.second << std::endl;
      if (dEdx.first>0 && dEdx.second>0)
        for (int i=0; i<5; ++i)
        {
          double BB = BB_curve(mass[i], momabs, _RC_parameters[i]);
          std::pair<float,float> ret{-1000,-1};
          if (fabs((dEdx.first-BB)/dEdx.second)<1000)
          {
            ret.first = (dEdx.first-BB)/dEdx.second;
            ret.second = dEdx.second;
          }
          else
          {sloM << "  " <<  mass[i] << " " << momabs << " " << BB << " " << dEdx.first << " " << dEdx.second << std::endl;}
          obsValues[i]=ret;
        }
    }

    return obsValues;

  }

  double InputAlgorithm_dEdx_RCD::BB_curve(double mass, double mom, std::vector<double> pars)
  {
    if (pars.size()<5) return 0;

    double bg = mom/mass;
    double bb  = sqrt(bg*bg/(1.+bg*bg));
    double tmax = pars[2]*pow(bg,2.);
    return 1.35e-7*(.5*pars[0]*log(pars[1]*pow(bg,2.)*tmax) - pars[3]*bb*bb - pars[4]*bg/2.)/(bb*bb);
  }

}
