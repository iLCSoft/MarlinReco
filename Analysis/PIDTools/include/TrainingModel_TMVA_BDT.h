#ifndef TrainingModel_TMVA_BDT_h
#define TrainingModel_TMVA_BDT_h 1

#include <TrainingModel.h>
#include <TMVA/Tools.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Reader.h>

using namespace lcio;
namespace cpid {

  class TrainingModel_TMVA_BDT : public TrainingModel {

  public:

    TrainingModel_TMVA_BDT(const TrainingModel_TMVA_BDT&) = delete;

    TrainingModel_TMVA_BDT& operator=(const TrainingModel_TMVA_BDT&) = delete;

    virtual ~TrainingModel_TMVA_BDT() = default;

    TrainingModel_TMVA_BDT();

    virtual TrainingModel* newModel() {return new TrainingModel_TMVA_BDT;}


    virtual std::vector<std::string> initTraining(TrainingModelInterface& tmi);

    virtual void initInference(TrainingModelInterface& tmi);

    virtual void runTraining(TTree* inTree);

    virtual const std::vector<float> runInference(int momBracket);



  private:

    TrainingModelInterface _tmi{};
    std::string _weightsfolder{};
    std::string _facOpt{};  // options when creating the factory
    std::string _facLod{};  // options when loader does PrepareAndTestTree
    std::string _facMet{};  // options when factory does BookMethod
    std::string _facCut{};  // general cuts; in addition to momentum bins and sig/bkg PDG cuts
    TMVA::Reader* _TReader{};

    int _nEvtMin4Train = 10;
  };

} // end namespace cpid

#endif 
