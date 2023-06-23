#include <TrainingModel_TMVA_BDT.h>
#include <ModelMgr.h>
//#include <InputAlgorithm_dEdx.hh>

//#include <filesystem>
#include <sys/stat.h>
#include <lcio.h>
#include <TFile.h>
#include <TCut.h>
#include <TMVA/Factory.h>

namespace cpid {

  TrainingModel_TMVA_BDT aTrainingModel_TMVA_BDT;

  TrainingModel_TMVA_BDT::TrainingModel_TMVA_BDT(): TrainingModel("TMVA_BDT")
  {
    _description = "Runs Root TMVA BDT training";
  }

  std::vector<std::string> TrainingModel_TMVA_BDT::initTraining(TrainingModelInterface& tmi)
  {
    _tmi = tmi;
    print_init(tmi.inparF, tmi.inparS);
    _weightsfolder = _modelName + std::string("/weights");

    if (tmi.inparS.size()<4) {sloE << "TrainingModel TMVA_BDT received too few (<4) string parameters for factory options!" << std::endl;}
    _facOpt = tmi.inparS[0];
    _facLod = tmi.inparS[1];
    _facMet = tmi.inparS[2];
    _facCut = tmi.inparS[3];

    struct stat meta;
    if (stat(_modelName.c_str(), &meta)) mkdir(_modelName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    for (unsigned int i=0; i<_tmi.momBins.size()-1; ++i)
    {
      std::stringstream s; s << "/TMVA_BDT_" << i << ".weights.xml";
      _tmi.weightFiles.push_back(_weightsfolder + s.str());
    }
    return _tmi.weightFiles;
  }

    void TrainingModel_TMVA_BDT::initInference(TrainingModelInterface& tmi)
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

  void TrainingModel_TMVA_BDT::runTraining(TTree* inTree)
  {
    for (unsigned int i=0; i<_tmi.momBins.size()-1; ++i)
    {
      std::stringstream m; m << "momabs > " << _tmi.momBins[i] << " && momabs < " << _tmi.momBins[i+1];
      TCut cutMom = m.str().c_str();
      TCut cutExt = _facCut.c_str();
      TCut cutSig = cutMom && cutExt;
      TCut cutBkg = cutMom && cutExt;

      TCut cutPDG = "";
      for (int pdg : _tmi.signalPDGs)
      {
        std::stringstream sb; sb << "PDG == " << pdg << " || PDG == " << -pdg;
        cutPDG = cutPDG || sb.str().c_str();
      }
      cutSig = cutSig && cutPDG;

      cutPDG = "";
      for (int pdg : _tmi.backgroundPDGs)
      {
        std::stringstream sb; sb << "PDG == " << pdg << " || PDG == " << -pdg;
        cutPDG = cutPDG || sb.str().c_str();
      }
      cutBkg = cutBkg && cutPDG;

      TMVA::DataLoader loader(_modelName.c_str());
      for (unsigned int v=0; v<_tmi.obsNames.size(); ++v) loader.AddVariable(_tmi.obsNames[v]);
      std::stringstream o; o << _modelName << "/TMVA_BDT_Out_" << i << ".root";
      TFile* outfile = TFile::Open(o.str().c_str(), "RECREATE");
      TMVA::Factory factory("TMVA", outfile, _facOpt);
      loader.AddSignalTree(inTree);
      loader.AddBackgroundTree(inTree);
      loader.PrepareTrainingAndTestTree(cutSig,cutBkg, _facLod);
      std::stringstream s; s << "BDT_" << i;
      factory.BookMethod(&loader, TMVA::Types::kBDT, s.str().c_str(), _facMet);

      factory.TrainAllMethods();
      factory.TestAllMethods();
      factory.EvaluateAllMethods();
    }
  }

  const std::vector<float> TrainingModel_TMVA_BDT::runInference(int momBracket)
  {
    std::stringstream nm; nm << "MomBin_" << momBracket;
    float r = _TReader->EvaluateMVA(nm.str());

    std::vector<float> eval_out{};
    for (unsigned i=0; i<_tmi.signalPDGs.size(); ++i) eval_out.push_back(r);
    for (unsigned i=0; i<_tmi.backgroundPDGs.size(); ++i) eval_out.push_back(-r);

    return eval_out;
  }

}
