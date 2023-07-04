#ifndef InputAlgorithm_dNdx_h
#define InputAlgorithm_dNdx_h 1

#include <InputAlgorithm.h>
#include <HelixClass.h>

using namespace lcio;
namespace cpid {

  class InputAlgorithm_dNdx : public InputAlgorithm {

  public:

    InputAlgorithm_dNdx(const InputAlgorithm_dNdx&) = delete;

    InputAlgorithm_dNdx& operator=(const InputAlgorithm_dNdx&) = delete;

    virtual ~InputAlgorithm_dNdx() = default;

    InputAlgorithm_dNdx();

    virtual InputAlgorithm* newAlgorithm() {return new InputAlgorithm_dNdx;}

    virtual std::vector<std::string> init(const std::vector<float>& inparF, const std::vector<std::string>& inparS);
  
    // The work horse
    virtual std::vector<std::pair<float,float> > extractObservables(ReconstructedParticleImpl* pfo, LCCollection* pfo_col, int PDG);

    float nClusters(float bg, int gas);

    float computePointInTracker(HelixClass& helix, float* ref, float radius, float zmax);

  private:

    int _gas=0;
    float _scaling=1, _bField=-3.5;
    float _innRad=0.5, _outRad=2, _zLen=2;

  };

} // end namespace cpid

#endif 
