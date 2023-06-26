#include <TrainingModel_TMVA_BDT_MC.h>
#include <ModelMgr.h>
//#include <InputAlgorithm_dEdx.hh>

//#include <filesystem>
#include <sys/stat.h>
#include <lcio.h>
#include <TFile.h>
#include <TCut.h>
#include <TMVA/Factory.h>

namespace cpid {

  TrainingModel_TMVA_BDT_MC aTrainingModel_TMVA_BDT_MC;

  TrainingModel_TMVA_BDT_MC::TrainingModel_TMVA_BDT_MC(): TrainingModel("TMVA_BDT_MC")
  {
    _description = "Runs Root TMVA Multiclass BDT training";
  }

  std::vector<std::string> TrainingModel_TMVA_BDT_MC::initTraining(TrainingModelInterface& tmi)
  {
    _tmi = tmi;
    print_init(tmi.inparF, tmi.inparS);
    _weightsfolder = _modelName + std::string("/weights");

    if (tmi.inparS.size()<4)
    {
      sloE << "TrainingModel TMVA_BDT_MC received too few (<4) string parameters for factory options!" << std::endl;
      throw std::runtime_error("parameters error");
    }
    _facOpt = tmi.inparS[0];
    _facLod = tmi.inparS[1];
    _facMet = tmi.inparS[2];
    _facCut = tmi.inparS[3];

    struct stat meta;
    if (stat(_modelName.c_str(), &meta)) mkdir(_modelName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    for (unsigned int i=0; i<_tmi.momBins.size()-1; ++i)
    {
      std::stringstream s; s << "/TMVA_BDT_MC_" << i << ".weights.xml";
      _tmi.weightFiles.push_back(_weightsfolder + s.str());
    }
    return _tmi.weightFiles;
  }

    void TrainingModel_TMVA_BDT_MC::initInference(TrainingModelInterface& tmi)
  {
    _tmi = tmi;
    if (_tmi.obsNames.size()!=_tmi.obsValues.size())
    {
      sloE << "Observables names (" << _tmi.obsNames.size() << ") and values (" << _tmi.obsValues.size() << ") have different sizes!" << std::endl;
      throw std::runtime_error("observables error");
    }
    _TReader = new TMVA::Reader("!Color:!Silent");
    for (unsigned int i=0; i<_tmi.obsNames.size(); ++i)
    {
      _TReader->AddVariable(_tmi.obsNames[i], _tmi.obsValues[i]);
      sloM << "Variable " << _tmi.obsNames[i] << " added to TReader" << std::endl;
    }

    for (unsigned int i=0; i<_tmi.weightFiles.size(); ++i)
    {
      std::stringstream nm; nm << "MomBin_" << i;
      _TReader->BookMVA(nm.str(), _tmi.weightFiles[i]);
    }
  }

  void TrainingModel_TMVA_BDT_MC::runTraining(TTree* inTree)
  {
    for (unsigned int i=0; i<_tmi.momBins.size()-1; ++i)
    {
      std::stringstream m; m << "momabs > " << _tmi.momBins[i] << " && momabs < " << _tmi.momBins[i+1];
      TCut cutMom = m.str().c_str();
      TCut cutExt = _facCut.c_str();
      TCut cutSig = cutMom && cutExt;

      TMVA::DataLoader loader(_modelName.c_str());
      for (unsigned int v=0; v<_tmi.obsNames.size(); ++v) loader.AddVariable(_tmi.obsNames[v]);
      std::stringstream o; o << _modelName << "/TMVA_BDT_MC_Out_" << i << ".root";
      TFile* outfile = TFile::Open(o.str().c_str(), "RECREATE");
      TMVA::Factory factory("TMVA", outfile, _facOpt);

      for (int pdg : _tmi.signalPDGs)
      {
        std::stringstream sb; sb << "PDG == " << pdg << " || PDG == " << -pdg;
        TCut cutPDG = cutSig && sb.str().c_str();
        std::stringstream nm; nm << "PDG=" << pdg;
        loader.AddTree(inTree, nm.str(), 1, cutPDG);
      }

      loader.PrepareTrainingAndTestTree("", _facLod);
      std::stringstream s; s << "BDT_MC_" << i;
      factory.BookMethod(&loader, TMVA::Types::kBDT, s.str().c_str(), _facMet);

      factory.TrainAllMethods();
      factory.TestAllMethods();
      factory.EvaluateAllMethods();
    }
  }

  const std::vector<float> TrainingModel_TMVA_BDT_MC::runInference(int momBracket)
  {
    std::stringstream nm; nm << "MomBin_" << momBracket;
    const std::vector<float> r = _TReader->EvaluateMulticlass(nm.str());

    std::vector<float> eval_out = r;
    for (unsigned i=0; i<_tmi.backgroundPDGs.size(); ++i) eval_out.push_back(0);

    return eval_out;
  }

}
