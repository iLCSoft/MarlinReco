#include <ComprehensivePIDProcessor.h>

#include <AlgorithmMgr.h>
#include <EVENT/LCCollection.h>
#include <IMPL/MCParticleImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <InputAlgorithm.h>
#include <ModelMgr.h>
#include <TrainingModel.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>

#include <TClass.h>
#include <TCut.h>
#include <TText.h>
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

using namespace lcio;
using namespace marlin;
using namespace cpid;

ComprehensivePIDProcessor aComprehensivePIDProcessor;

ComprehensivePIDProcessor::ComprehensivePIDProcessor() : Processor("ComprehensivePIDProcessor") {

  _description = "ComprehensivePIDProcessor: Determines and evaluates PID";

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "PFOCollection", "PFO input collection", _PFOColName,
                          std::string("PandoraPFOs"));

  registerInputCollection(LCIO::LCRELATION, "RecoMCTruthLink",
                          "LCRelation from the reconstructed PFOs to the MC particles", _RecoMCTruthLinkName,
                          std::string("RecoMCTruthLink"));

  registerProcessorParameter("modeExtract",
                             "Set true to extract the specified observables from the PFOs; if false you need to set "
                             "the TTreeFileName instead; default: true.",
                             _modeExtract, true);

  registerProcessorParameter(
      "modeTrain", "Set true to internally train an MVA with the specified observables from the PFOs; default: false.",
      _modeTrain, false);

  registerProcessorParameter("modeInfer",
                             "Set true to infer the trained MVA with the specified observables to the PFOs; if true "
                             "you need to provide the reference and weight files; default: false.",
                             _modeInfer, false);

  registerProcessorParameter(
      "TTreeFileName",
      "Name of the root file in which the TTree with all observables is stored; in case of extraction it is an "
      "optional output with no output if left empty, otherwise it is a necessary input; default: TTreeFile.root",
      _TTreeFileName, std::string("TTreeFile.root"));

  registerProcessorParameter(
      "inputAlgoSpecs",
      "List of input algorithms; for each specify type:name or only type (then name=type); default: {}",
      _inputAlgoSpecs, std::vector<std::string>());

  registerProcessorParameter(
      "trainModelSpecs",
      "List of training models; for each specify type:name or only type (then name=type); default: {}",
      _trainModelSpecs, std::vector<std::string>());

  std::vector<std::string> ref{std::string("Ref.txt")};
  registerProcessorParameter("reffile",
                             "Reference file(s), if only one file but several training models are specified the "
                             "reference files are auto-numbered; default: {Ref.txt}.",
                             _reffile, ref);

  registerProcessorParameter("trainingObservables",
                             "List of observables that should be used for traning; if empty, all observables from the "
                             "specified algorithms + momabs + lambda are used; default: {}",
                             _trainingObservables, std::vector<std::string>());

  std::vector<int> signalPDGs = {11, 13, 211, 321, 2212};
  registerProcessorParameter("signalPDGs",
                             "List of PDG numbers that are considered signal; default: {11,13,211,321,2212}.",
                             _signalPDGs, signalPDGs);

  registerProcessorParameter("backgroundPDGs", "List of PDG numbers that are considered background; default: {}.",
                             _backgroundPDGs, std::vector<int>());

  registerProcessorParameter(
      "plotFolder",
      "Folder in which the automatic confusion matrix plots of inference will be put, is created if not already "
      "existing; if empty, no plots are created; default: CPID_Plots  [current working directory]",
      _plotFolder, std::string("CPID_Plots"));

  registerProcessorParameter("fileFormat",
                             "File format of the automatic confusion matrix plots of inference; default: .png",
                             _fileFormat, std::string(".png"));

  registerProcessorParameter("momMin", "For momentum bins: minimum momentum cut / GeV; default: 1.", _momMin, float(1));

  registerProcessorParameter("momMax", "For momentum bins: maximum momentum cut / GeV; default: 100.", _momMax,
                             float(100));

  registerProcessorParameter("momLog", "For momentum bins: should the momentum bins be logarihtmic; default: true.",
                             _momLog, true);

  registerProcessorParameter("momNBins", "For momentum bins: number of momentum bins; default: 12.", _momNBins,
                             int(12));

  registerProcessorParameter("cutD0",
                             "PFOs whose first track have a d0 larger than the given value will be ignored; set to 0 "
                             "to accept all particles; default: 0.",
                             _cutD0, float(0));

  registerProcessorParameter("cutZ0",
                             "PFOs whose first track have a z0 larger than the given value will be ignored; set to 0 "
                             "to accept all particles; default: 0.",
                             _cutZ0, float(0));

  registerProcessorParameter("cutLamMin",
                             "PFOs whose first track have an angle lambda (relative to the cathode) smaller than the "
                             "given value will be ignored; set to 0 to accept all particles; default: 0.",
                             _cutLamMin, float(0));

  registerProcessorParameter("cutLamMax",
                             "PFOs whose first track have an angle lambda (relative to the cathode) larger than the "
                             "given value will be ignored; set to 0 to accept all particles; default: 0.",
                             _cutLamMax, float(0));

  registerProcessorParameter(
      "cutNTracksMin",
      "PFOs with fewer (<) tracks than the given value are ignored; set to -1 to accept all PFOs; default: -1.",
      _cutNTracksMin, int(-1));

  registerProcessorParameter(
      "cutNTracksMax",
      "PFOs with more (>) tracks the given value are ignored; set to -1 to accept all PFOs; default: -1.",
      _cutNTracksMax, int(-1));
}

void ComprehensivePIDProcessor::init() {
  // usually a good idea to
  printParameters();

  if (_modeTrain && _modeInfer) {
    sloE << "I cannot train and infer in the same process!" << std::endl;
    throw std::runtime_error("mode error");
  }

  if (_TTreeFileName == "") {
    sloM << "TTreeFileName is empty, no root file is created" << std::endl;
  }
  if (_plotFolder == "") {
    sloM << "plotFolder is empty, no plots are created" << std::endl;
  }

  _nEvt = 0;
  _nPFO = 0;

  _nAlgos = _inputAlgoSpecs.size();
  _nModels = _trainModelSpecs.size();

  if (_modeExtract) {
    if (_TTreeFileName != "") {
      _writeTTreeFile = true;
      _TTreeFile = new TFile(_TTreeFileName.c_str(), "RECREATE");
      _TTreeFile->cd();
      _observablesTree = new TTree("observablesTree", "Tree of all observable values");
    }
  } else {
    _TTreeFile = new TFile(_TTreeFileName.c_str());
    _observablesTree = (TTree*)_TTreeFile->Get("observablesTree");
  }

  _observablesNames.push_back("momabs");
  _observablesNames.push_back("lambda");
  _observablesNames.push_back("d0");
  _observablesNames.push_back("z0");
  _observablesNames.push_back("PDG");
  _observablesNames.push_back("nTracks");
  _nObs_base = _observablesNames.size();

  // set up momentum bins
  _momBins.push_back(_momMin);
  sloM << "momBins: " << _momBins[0];
  for (int i = 0; i < _momNBins; ++i) {
    if (_momLog)
      _momBins.push_back(pow(10, log10(_momMin) + (log10(_momMax) - log10(_momMin)) * (i + 1) / (float)(_momNBins)));
    else
      _momBins.push_back(_momMin + (_momMax - _momMin) * (i + 1) / (float)(_momNBins));
    sloM << " " << _momBins[i + 1];
  }
  sloM << std::endl << "-------------------------------------------------" << std::endl;

  // prepare optional parameters of modules
  std::shared_ptr<StringParameters> pars = parameters();
  std::vector<float> inparF{};
  std::vector<std::string> inparS{};

  // set up and initialise input algorithms
  AlgorithmMgr::instance()->printAvailableAlgorithmTypes();

  for (int i = 0; i < _nAlgos; ++i) {
    // algorithms can be specified as type:name or just type (then name = type)
    std::size_t p = _inputAlgoSpecs[i].find(":");
    if (p == std::string::npos) {
      _inputAlgoTypes.push_back(_inputAlgoSpecs[i]);
      _inputAlgoNames.push_back(_inputAlgoSpecs[i]);
    } else {
      _inputAlgoTypes.push_back(_inputAlgoSpecs[i].substr(0, p));
      _inputAlgoNames.push_back(_inputAlgoSpecs[i].substr(p + 1));
    }

    sloD << "creating new algorithm: " << _inputAlgoTypes[i] << ":" << _inputAlgoNames[i] << std::endl;

    _inputAlgorithms.push_back(AlgorithmMgr::instance()->createAlgorithm(_inputAlgoTypes[i]));
    _inputAlgorithms[i]->setName(_inputAlgoNames[i]);

    inparF.clear();
    inparS.clear();
    std::string sF = _inputAlgoNames[i] + ".F";
    std::string sS = _inputAlgoNames[i] + ".S";
    pars->getFloatVals(sF, inparF);
    pars->getStringVals(sS, inparS);

    std::vector<std::string> obsNames = _inputAlgorithms[i]->init(inparF, inparS);
    for (unsigned int j = 0; j < obsNames.size(); ++j) {
      std::string s = _inputAlgoNames[i];
      s.append("_");
      _observablesNames.push_back(s.append(obsNames[j]));
    }

    sloM << "New algorithm created: " << _inputAlgorithms[i]->type() << ":" << _inputAlgorithms[i]->name() << "\n"
         << std::endl;
  }

  sloM << "Input algorithm inits done" << std::endl;
  sloM << "-------------------------------------------------" << std::endl;

  _nObs = _observablesNames.size();
  _observablesValues.resize(_nObs, 0);

  // set up and initialise training/inference models
  ModelMgr::instance()->printAvailableModelTypes();

  if (_modeTrain || _modeInfer) {
    for (int i = 0; i < _nModels; ++i) {
      // models can be specified as type:name or just type (then name = type)
      std::size_t p = _trainModelSpecs[i].find(":");
      if (p == std::string::npos) {
        _trainModelTypes.push_back(_trainModelSpecs[i]);
        _trainModelNames.push_back(_trainModelSpecs[i]);
      } else {
        _trainModelTypes.push_back(_trainModelSpecs[i].substr(0, p));
        _trainModelNames.push_back(_trainModelSpecs[i].substr(p + 1));
      }

      _trainingModels.push_back(ModelMgr::instance()->createModel(_trainModelTypes[i]));
      _trainingModels[i]->setName(_trainModelNames[i]);

      sloM << "New model created: " << _trainingModels[i]->type() << ":" << _trainingModels[i]->name() << std::endl;

      if (_modeInfer) {
        _nMomBins.push_back(0);
        std::vector<std::pair<float, float>> none{};
        _weightFileBrackets.push_back(none);
        std::vector<std::string> empty{};
        _weightFiles.push_back(empty);
        ReadReferenceFile(i);
      }
      sloM << std::endl;
    }

    // if training observables are not specified in the steering file use the extracted observables + momabs and lambda
    // (without d0, z0, MCPDG and nTrack)
    if (_trainingObservables.size() == 0) {
      _trainingObservables = _observablesNames;
      _trainingObservables.erase(_trainingObservables.begin() + 2, _trainingObservables.begin() + 6);
    }

    sloM << "Training/inference observables: " << _trainingObservables.size() << " -";
    for (unsigned int o = 0; o < _trainingObservables.size(); ++o) {
      sloM << " " << _trainingObservables[o];
    }
    sloM << std::endl << "Extracted observables: " << _observablesNames.size() << " -";
    for (unsigned int o = 0; o < _observablesNames.size(); ++o) {
      sloM << " " << _observablesNames[o];
    }
    sloM << std::endl;

    bool obs_ok = true;
    if (_modeExtract) {
      obs_ok = CheckTrainingObservables(_trainingObservables, _observablesNames);
      if (!obs_ok) {
        sloE << "Target training/inference observables not (completely) among extracted observables!" << std::endl;
        throw std::runtime_error("oberservables error");
      }
    }
    // todo: else: check if target obs are in root file

    if (_modeInfer) {
      for (std::string tObs : _trainingObservables)
        for (int p = 0; p < _nObs; ++p)
          if (tObs == _observablesNames[p])
            _inferObsValues.push_back(&(_observablesValues[p]));
    }

    TrainingModelInterface tmi;
    tmi.signalPDGs = _signalPDGs;
    tmi.backgroundPDGs = _backgroundPDGs;
    tmi.obsNames = _trainingObservables;
    tmi.momBins = _momBins;

    for (int i = 0; i < _nModels; ++i) {
      tmi.inparF.clear();
      tmi.inparS.clear();
      std::string sF = _trainModelNames[i] + ".F";
      std::string sS = _trainModelNames[i] + ".S";
      pars->getFloatVals(sF, tmi.inparF);
      pars->getStringVals(sS, tmi.inparS);

      if (_modeTrain) {
        _weightFiles.push_back(_trainingModels[i]->initTraining(tmi));
        CreateReferenceFile(i);
      }

      if (_modeInfer) {
        tmi.obsValues = _inferObsValues;
        tmi.weightFiles = _weightFiles[i];
        _trainingModels[i]->initInference(tmi);

        int nSig = _signalPDGs.size();
        std::stringstream h1;
        h1 << "PDGCheck_" << i;
        std::stringstream h2;
        h2 << "PDG Confusion Matrix, " << _trainingModels[i]->name();
        _PDGCheck.push_back(new TH2I(h1.str().c_str(), h2.str().c_str(), nSig, -.5, nSig - .5, nSig, -.5, nSig - .5));
        _PDGCheck[i]->SetXTitle("MCTruth PDG");
        _PDGCheck[i]->SetYTitle("PFO PDG");
        for (int p = 0; p < nSig; p++) {
          std::stringstream s;
          s << _signalPDGs[p];
          _PDGCheck[i]->GetXaxis()->SetBinLabel(p + 1, s.str().c_str());
          _PDGCheck[i]->GetYaxis()->SetBinLabel(p + 1, s.str().c_str());
        }
        //_PDGCheck->GetXaxis()->SetBinLabel(nSig+1,"other");
        //_PDGCheck->GetYaxis()->SetBinLabel(nSig+1,"other");
        sloD << _PDGCheck[i]->GetName() << std::endl;

        std::vector<TH2I*> tempV{};
        for (int m = 0; m < _nMomBins[i]; ++m) {
          std::stringstream t1;
          t1 << "PDGCheck_" << i << "_" << m;
          std::stringstream t2;
          t2 << "PDG Confusion Matrix, " << _trainingModels[i]->name() << ", " << _weightFileBrackets[i][m].first
             << "<p/GeV<" << _weightFileBrackets[i][m].second;
          tempV.push_back(new TH2I(t1.str().c_str(), t2.str().c_str(), nSig, -.5, nSig - .5, nSig, -.5, nSig - .5));
          tempV[m]->SetXTitle("MCTruth PDG");
          tempV[m]->SetYTitle("PFO PDG");
          for (int p = 0; p < nSig; p++) {
            std::stringstream s;
            s << _signalPDGs[p];
            tempV[m]->GetXaxis()->SetBinLabel(p + 1, s.str().c_str());
            tempV[m]->GetYaxis()->SetBinLabel(p + 1, s.str().c_str());
          }
        }
        _PDGChecks.push_back(tempV);
        for (int m = 0; m < _nMomBins[i]; ++m) {
          sloD << _PDGChecks[i][m]->GetName() << std::endl;
        }
      }
    }

    sloM << "Model inits done" << std::endl;
    sloM << "-------------------------------------------------" << std::endl;
  }

  if (_modeExtract && _writeTTreeFile) {
    for (int iObs = 0; iObs < _nObs; ++iObs) {
      std::stringstream s;
      s << _observablesNames[iObs].c_str() << "/F";
      _observablesTree->Branch(_observablesNames[iObs].c_str(), &(_observablesValues[iObs]), s.str().c_str());
    }
  }

  sloM << "Branches done" << std::endl;
  sloM << "-------------------------------------------------" << std::endl;

  for (int i = 0; i < _nReason; ++i)
    _nRejectedPFOs[i] = 0;
  _rejectedPFOs = new TH1D("rejectedPFOs", "Number of rejected PFOs by reason", _nReason, .5, _nReason + .5);
  _rejectedPFOs->GetXaxis()->SetBinLabel(1, "no MCPDG found");
  _rejectedPFOs->GetXaxis()->SetBinLabel(2, "MCPDG not among signal PDGs");
  _rejectedPFOs->GetXaxis()->SetBinLabel(3, "too few tracks");
  _rejectedPFOs->GetXaxis()->SetBinLabel(4, "too many tracks");
  _rejectedPFOs->GetXaxis()->SetBinLabel(5, "too large d0");
  _rejectedPFOs->GetXaxis()->SetBinLabel(6, "too large z0");
  _rejectedPFOs->SetYTitle("abundance");
}

void ComprehensivePIDProcessor::processRunHeader(LCRunHeader*) {
  _nRun++;
  //_nEvt = 0;
}

void ComprehensivePIDProcessor::processEvent(LCEvent* evt) {

  // sloM << _nEvt << std::endl;

  if (_nAlgos == 0)
    return;

  if (_modeExtract || _modeInfer) {
    LCCollection *col_pfo{}, *col_pfo2mc{};

    try {
      col_pfo = evt->getCollection(_PFOColName);
      col_pfo2mc = evt->getCollection(_RecoMCTruthLinkName);
    } catch (DataNotAvailableException& e) {
      sloM << "Input PFO or RecoMCTruthLink collection not found - skipping event " << _nEvt << std::endl;
      return;
    }

    int n_pfo = col_pfo->getNumberOfElements();
    LCRelationNavigator rel_pfo2mc(col_pfo2mc);

    PIDHandler pidh(col_pfo);
    std::vector<int> algoID{};
    std::vector<int> allPDGs{};
    if (_modeInfer) {
      std::vector<std::string> PDGness{};
      for (int pdg : _signalPDGs) {
        std::stringstream pn;
        pn << pdg << "-ness";
        PDGness.push_back(pn.str());
        allPDGs.push_back(pdg);
      }
      for (int pdg : _backgroundPDGs) {
        std::stringstream pn;
        pn << pdg << "-ness";
        PDGness.push_back(pn.str());
        allPDGs.push_back(pdg);
      }
      for (int m = 0; m < _nModels; ++m)
        algoID.push_back(pidh.addAlgorithm(_trainModelNames[m], PDGness));
    }

    for (int i = 0; i < n_pfo; ++i) {
      ReconstructedParticleImpl* pfo = dynamic_cast<ReconstructedParticleImpl*>(col_pfo->getElementAt(i));

      // get the true PDG  (same as in dEdxAnalyser)
      const LCObjectVec& mcparVec = rel_pfo2mc.getRelatedToObjects(pfo);
      const FloatVec& mcparWei = rel_pfo2mc.getRelatedToWeights(pfo);
      int MCPDG = 0;
      float bestwei = 0;
      int bestk = -1;
      for (unsigned k = 0; k < mcparVec.size(); ++k) {
        MCParticleImpl* mcpar = dynamic_cast<MCParticleImpl*>(mcparVec[k]);
        float trkwei = (int(mcparWei[k]) % 10000) / 1000.;
        float cluwei = (int(mcparWei[k]) / 10000) / 1000.;
        float mcparwei = trkwei > cluwei ? trkwei : cluwei;
        if (mcparwei > bestwei) {
          MCPDG = mcpar->getPDG();
          bestwei = mcparwei;
          bestk = k;
        }
      }

      // apply and register cuts
      if (_modeTrain) {
        if (bestk == -1) {
          ++_nRejectedPFOs[0];
          for (unsigned k = 0; k < mcparVec.size(); ++k) {
            sloD << mcparWei[k] << "|" << (int(mcparWei[k]) % 10000) / 1000. << "  ";
          }
          sloD << std::endl;
          continue;
        }
        if (std::find(_signalPDGs.begin(), _signalPDGs.end(), abs(MCPDG)) == _signalPDGs.end() &&
            std::find(_backgroundPDGs.begin(), _backgroundPDGs.end(), abs(MCPDG)) == _backgroundPDGs.end()) {
          ++_nRejectedPFOs[1];
          sloD << MCPDG << std::endl;
          continue;
        }
      }

      int nTracks = pfo->getTracks().size();
      if (_cutNTracksMin > -1 && nTracks < _cutNTracksMin) {
        ++_nRejectedPFOs[2];
        continue;
      }
      if (_cutNTracksMax > -1 && nTracks > _cutNTracksMax) {
        ++_nRejectedPFOs[3];
        continue;
      }
      if (nTracks > 0) {
        if (_cutD0 && pfo->getTracks()[0]->getD0() > _cutD0) {
          ++_nRejectedPFOs[4];
          continue;
        }
        if (_cutZ0 && pfo->getTracks()[0]->getZ0() > _cutZ0) {
          ++_nRejectedPFOs[5];
          continue;
        }
      }

      ++_nPFO;

      if (_modeExtract) {
        sloD << " extracting... " << std::endl;

        const double* mom = pfo->getMomentum();
        _momabs = sqrt(mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2]);
        _observablesValues[0] = _momabs;

        double tanl = mom[2] / sqrt(mom[0] * mom[0] + mom[1] * mom[1]);
        _observablesValues[1] = fabs(atan(tanl)) * 180 / M_PI;

        if (nTracks > 0) {
          _observablesValues[2] = pfo->getTracks()[0]->getD0();
          _observablesValues[3] = pfo->getTracks()[0]->getZ0();
        }

        _observablesValues[4] = MCPDG;
        _observablesValues[5] = nTracks;

        int iObs = _nObs_base;

        for (int a = 0; a < _nAlgos; ++a) {
          // sloM << "using algorithm " << _inputAlgorithms[a]->type() << std::endl;
          std::vector<std::pair<float, float>> obs = _inputAlgorithms[a]->extractObservables(pfo, col_pfo, MCPDG);
          for (unsigned int j = 0; j < obs.size(); ++j) {
            _observablesValues[iObs] = obs[j].first;
            ++iObs;
          }
        }
        sloD << "observables extracted " << std::endl;

        if (_writeTTreeFile)
          _observablesTree->Fill();
        sloD << "trees filled" << std::endl;
      }

      if (_modeInfer) {
        sloD << "start inference" << std::endl;
        for (int m = 0; m < _nModels; ++m)
          for (int j = 0; j < _nMomBins[m]; ++j) {
            if (_momabs < _weightFileBrackets[m][j].first || _momabs > _weightFileBrackets[m][j].second)
              continue;
            const std::vector<float> eval_out = _trainingModels[m]->runInference(j);

            auto eval_max = max_element(eval_out.begin(), eval_out.end());
            int index = std::distance(eval_out.begin(), eval_max);

            for (unsigned int k = 0; k < eval_out.size(); ++k) {
              sloD << " " << eval_out[k];
            }
            sloD << " | " << index << " " << allPDGs[index] << " " << abs(MCPDG) << std::endl;

            pidh.setParticleID(pfo, 0, allPDGs[index], 0, algoID[m], eval_out);
            for (int ipdg = 0; ipdg < (int)_signalPDGs.size(); ++ipdg)
              if (_signalPDGs[ipdg] == abs(MCPDG)) {
                _PDGCheck[m]->Fill(ipdg, index);
                _PDGChecks[m][j]->Fill(ipdg, index);
              }
            // sloM << MCPDG << " " << index << " " << allPDGs[index] << " " << _PDGCheck->GetEntries() << std::endl;
            break;
          }
      }
    } // for pfos
  }

  _nEvt++;
}

void ComprehensivePIDProcessor::check(LCEvent*) {}
void ComprehensivePIDProcessor::end() {
  sloM << "-------------------------------------------------" << std::endl << "end()" << std::endl;
  for (int i = 0; i < _nReason; ++i) {
    _rejectedPFOs->SetBinContent(i + 1, _nRejectedPFOs[i]);
    sloM << _rejectedPFOs->GetXaxis()->GetBinLabel(i + 1) << " " << _nRejectedPFOs[i] << "\n";
  }
  sloM << "PFOs passed: " << _nPFO << std::endl;

  if (_writeTTreeFile) {
    _TTreeFile->cd();
    _observablesTree->Write();
    _rejectedPFOs->Write();
  }

  if (_modeTrain) {
    sloM << "Begin training..." << std::endl << std::endl;
    _observablesTree->Print();

    for (int i = 0; i < _nModels; ++i)
      _trainingModels[i]->runTraining(_observablesTree);
    sloM << "Training done" << std::endl << std::endl;
  }

  if (_modeInfer && _plotFolder != "") {
    bool cont = true;
    if (!std::filesystem::exists(_plotFolder)) {
      try {
        std::filesystem::create_directory(std::filesystem::path(_plotFolder));
      } catch (...) {
        sloM << "plotFolder does not exist and could not be created, plotting is skipped.";
        cont = false;
      }
    }
    if (cont) {
      TCanvas* can = new TCanvas("Canvas2D", "Canvas2D", 500, 500);
      gStyle->SetPalette(kBird);
      can->SetLogz();
      can->SetGrid(0, 0);
      gStyle->SetOptStat(0);

      for (int m = 0; m < _nModels; ++m) {
        PlotTH2(can, _PDGCheck[m]);
        PlotTH2(can, _PDGCheck[m], 1);
        if (_writeTTreeFile)
          _PDGCheck[m]->Write();
        for (int j = 0; j < _nMomBins[m]; ++j) {
          PlotTH2(can, _PDGChecks[m][j]);
          if (_writeTTreeFile)
            _PDGChecks[m][j]->Write();
        }
      }
    }
  }

  if (_writeTTreeFile)
    _TTreeFile->Close();

  for (int i = 0; i < _nAlgos; ++i)
    delete _inputAlgorithms[i];
  for (int i = 0; i < _nModels; ++i)
    delete _trainingModels[i];

  sloM << "Numer of Events: " << _nEvt << "    Numer of PFOs: " << _nPFO << std::endl;
  sloM << "-------------------------------------------------" << std::endl;
}

bool ComprehensivePIDProcessor::CheckTrainingObservables(const std::vector<std::string>& trainObs,
                                                         const std::vector<std::string>& compObs) {
  for (const std::string& tObs : trainObs) {
    if (std::find(compObs.begin(), compObs.end(), tObs) == compObs.end())
      return false;
  }
  return true;
}

void ComprehensivePIDProcessor::CreateReferenceFile(int n) {
  sloM << "Creating reference file..." << std::endl;
  std::string refname = ReferenceFile(n);
  std::ofstream reffile;
  reffile.open(refname);

  for (std::string obs : _trainingObservables) {
    reffile << obs << " ";
  }
  reffile << "\n";

  for (int pdg : _signalPDGs) {
    reffile << pdg << " ";
  }
  reffile << "| ";
  for (int pdg : _backgroundPDGs) {
    reffile << pdg << " ";
  }
  reffile << "\n\n";

  for (int i = 0; i < _momNBins; ++i) {
    reffile << _momBins[i] << " " << _momBins[i + 1] << " " << _weightFiles[n][i] << "\n";
  }
  reffile.close();
  sloM << "Reference file " << refname << " created." << std::endl;
}

void ComprehensivePIDProcessor::ReadReferenceFile(int n) {
  std::string refname = ReferenceFile(n);
  std::vector<std::string> lineV;
  std::string line;
  std::ifstream reffile(refname);
  if (reffile.is_open()) {
    getline(reffile, line);
    boost::split(lineV, line, [](char c) { return c == ' '; });
    if (lineV.back() == "")
      lineV.pop_back();
    if (n == 0)
      _trainingObservables = lineV;
    if (n > 0 and lineV.size() != _trainingObservables.size()) {
      sloE << "Number of target observables not consistent among reference files!" << std::endl;
      throw std::runtime_error("oberservables error");
    }
    if (n > 0 and not CheckTrainingObservables(_trainingObservables, lineV)) {
      sloE << "Specified target observables not consistent among reference files!" << std::endl;
      throw std::runtime_error("oberservables error");
    }

    getline(reffile, line);
    if (n == 0) {
      boost::split(lineV, line, [](char c) { return c == ' '; });
      for (std::string l : lineV)
        if (l == "|")
          _SnB.push_back(0);
        else if (l != "")
          _SnB.push_back(stoi(l));
    }
    // todo: make also SnB safe against incoherent reference files

    while (getline(reffile, line)) {
      boost::split(lineV, line, [](char c) { return c == ' '; });
      if (lineV.size() <= 1)
        continue;
      std::pair<float, float> mombracket(std::stof(lineV[0]), std::stof(lineV[1]));
      _weightFileBrackets[n].push_back(mombracket);
      _weightFiles[n].push_back(lineV[2]);
    }
    reffile.close();
    _nMomBins[n] = _weightFileBrackets[n].size();
  } else {
    sloE << "Reference file " << refname << " could not be opened!" << std::endl;
    throw std::runtime_error("reffile error");
  }
  sloM << "Reference file successfully read." << std::endl;
}

std::string ComprehensivePIDProcessor::ReferenceFile(int n) {
  std::string refname{};
  if (_nModels > int(_reffile.size())) {
    refname = _reffile[0];
    std::size_t p = refname.rfind(".");
    refname.insert(p, std::string("_") + std::to_string(n));
  } else
    refname = _reffile[n];
  return refname;
}

void ComprehensivePIDProcessor::PlotTH2(TCanvas* can, TH2* hist, int effpur) {
  hist->Draw("colz");
  std::stringstream s;
  s << _plotFolder << "/" << hist->GetName();

  std::vector<double> sumX{}, sumY{};
  if (effpur) // calculate row and column sums
  {
    int nX = hist->GetNbinsX();
    int nY = hist->GetNbinsY();
    sumX.resize(nY);
    sumY.resize(nX);
    std::fill(sumX.begin(), sumX.end(), 0);
    std::fill(sumY.begin(), sumY.end(), 0);
    for (int i = 0; i < nX; ++i)
      for (int j = 0; j < nY; ++j) {
        sumX[j] = sumX[j] + hist->GetBinContent(i + 1, j + 1);
        sumY[i] = sumY[i] + hist->GetBinContent(i + 1, j + 1);
      }
  }
  if (effpur == 1 || effpur == 3) // print efficiency
  {
    for (int i = 0; i < hist->GetNbinsX(); ++i) {
      std::stringstream e;
      e << hist->GetBinContent(i + 1, i + 1) / sumY[i];
      std::string eff = e.str();
      eff.resize(4);
      TText* t = new TText(i - .45, i + .05, eff.c_str());
      t->Draw();
    }
    s << "_eff";
  }
  if (effpur == 2 || effpur == 3) // print purity
  {
    for (int j = 0; j < hist->GetNbinsY(); ++j) {
      std::stringstream p;
      p << hist->GetBinContent(j + 1, j + 1) / sumX[j];
      std::string pur = p.str();
      pur.resize(4);
      TText* t = new TText(j - .2, j - .45, pur.c_str());
      t->Draw();
    }
    s << "_pur";
  }

  can->Update();
  s << _fileFormat;
  can->Print(s.str().c_str());
}
