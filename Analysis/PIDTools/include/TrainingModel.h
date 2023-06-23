#ifndef TrainingModel_h
#define TrainingModel_h 1

#include <string>
#include <vector>
#include <lcio.h>
#include <InputAlgorithm.h>
#include <streamlog/streamlog.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/LCCollection.h>
#include <TTree.h>

using namespace lcio;
namespace cpid {

  struct TrainingModelInterface
    {
      std::vector<float> inparF{};
      std::vector<std::string> inparS{};

      std::vector<int> signalPDGs{};
      std::vector<int> backgroundPDGs{};
      std::vector<float> momBins{};
      std::vector<float*> obsValues{};
      std::vector<std::string> obsNames{};
      std::vector<std::string> weightFiles{};
    };


  class TrainingModel {

  public:

    TrainingModel(const TrainingModel&) = delete;

    TrainingModel& operator=(const TrainingModel&) = delete;

    virtual ~TrainingModel() = default;

    TrainingModel(const std::string& typeName);

    virtual TrainingModel* newModel() = 0;

    virtual const std::string & type() const {return _typeName;}

    virtual const std::string & name() const {return _modelName;}

    virtual void setName(std::string name) {_modelName = name;}


    virtual std::vector<std::string> initTraining(TrainingModelInterface&);

    virtual void initInference(TrainingModelInterface&);

    void print_init(std::vector<float> inparF, std::vector<std::string> inparS)
    {
      sloM << _modelName << " init float parameters: ";
      for (unsigned int i=0; i<inparF.size(); ++i) {sloM << inparF[i] << " ";}
      sloM << std::endl;
      sloM << _modelName << " init string parameters: ";
      for (unsigned int i=0; i<inparS.size(); ++i) {sloM << inparS[i] << " ";}
      sloM << std::endl;
    }

    virtual void runTraining(TTree*);

    virtual const std::vector<float> runInference(int);


  protected:

    std::string _typeName{};
    std::string _modelName{};
    std::string _description{};


  };

} // end namespace cpid

#endif 
