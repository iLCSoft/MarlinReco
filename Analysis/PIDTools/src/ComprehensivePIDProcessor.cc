#include <ComprehensivePIDProcessor.h>

#include <EVENT/LCCollection.h>
#include <UTIL/LCRelationNavigator.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/MCParticleImpl.h>
#include <InputAlgorithm.h>
#include <AlgorithmMgr.h>
#include <TrainingModel.h>
#include <ModelMgr.h>
#include <UTIL/PIDHandler.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <TClass.h>

using namespace lcio;
using namespace marlin;
using namespace cpid;

ComprehensivePIDProcessor aComprehensivePIDProcessor;

ComprehensivePIDProcessor::ComprehensivePIDProcessor() : Processor("ComprehensivePIDProcessor") {

  _description = "ComprehensivePIDProcessor: Determines and evaluates PID" ;

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
               "PFOCollection",
               "Pandora PFOs",
               _PFOColName,
               std::string("PandoraPFOs"));

  registerInputCollection(LCIO::LCRELATION,
           "RecoMCTruthLink",
           "LCRelation from the reconstructed PFOs to the MC particles",
           _RecoMCTruthLinkName,
           std::string("RecoMCTruthLink"));

  registerProcessorParameter("modeExtract",
                 "Set true to extract the specified observables from the PFOs; if false you need to set the TTreeFileName instead; default: true.",
                 _modeExtract,
                 true);

  registerProcessorParameter("modeTrain",
                 "Set true to internally train an MVA with the specified observables from the PFOs; default: false.",
                 _modeTrain,
                 false);

  registerProcessorParameter("modeInfer",
                 "Set true to infer the trained MVA with the specified observables to the PFOs; if true you need to provide the weight files; default: false.",
                 _modeInfer,
                 false);

  registerProcessorParameter("TTreeFileName",
                 "Name of the file in which the TTree with all observables is stored; optional output in case of extraction, otherwise necessary input; default: TTreeFile.root",
                 _TTreeFileName,
                 std::string("TTreeFile.root"));

  registerProcessorParameter("inputAlgoSpecs",
                 "List of input algorithms; for each specify only type (name is then type) or type:name; default: {}",
                 _inputAlgoSpecs,
                 std::vector<std::string>());

  registerProcessorParameter("trainModelSpecs",
                 "List of training models; for each specify only type (name is then type) or type:name; default: {}",
                 _trainModelSpecs,
                 std::vector<std::string>());

  registerProcessorParameter("trainingObservables",
                 "List of observables that should be used for traning; if empty, all observables from the specified algorithms + momabs + lambda are used; default: {}",
                 _trainingObservables,
                 std::vector<std::string>());

  std::vector<std::string> ref{std::string("Ref.txt")};
  registerProcessorParameter("reffile",
                 "Reference file(s), if only one file but several training models are specified the reference files are auto-numbered; default: {Ref.txt}.",
                 _reffile,
                 ref);

  std::vector<int> signalPDGs = {11,13,211,321,2122};
  registerProcessorParameter("signalPDGs",
                 "Vector of PDG numbers that are considered signal; default: {11,13,211,321,2212}.",
                 _signalPDGs,
                 signalPDGs);

  registerProcessorParameter("backgroundPDGs",
                 "Vector of PDG numbers that are considered signal; default: {}.",
                 _backgroundPDGs,
                 std::vector<int>());


  registerProcessorParameter("momMin",
                 "For training: minimum momentum cut / GeV; default: 1.",
                 _momMin,
                 float(1));

  registerProcessorParameter("momMax",
                 "For training: maximum momentum cut / GeV; default: 100.",
                 _momMax,
                 float(100));

  registerProcessorParameter("momLog",
                 "For training: should the momentum bins be logarihtmic; default: true.",
                 _momLog,
                 true);

  registerProcessorParameter("momNBins",
                 "For training: number of momentum bins; default: 12.",
                 _momNBins,
                 int(12));

  registerProcessorParameter("cutD0",
                 "Tracks with a d0 larger than the given value will be ignored. Set to 0 to accept all particles.",
                 _cutD0,
                 float(0));

  registerProcessorParameter("cutZ0",
                 "Tracks with a z0 larger than the given value will be ignored. Set to 0 to accept all particles.",
                 _cutZ0,
                 float(0));

  registerProcessorParameter("cutLamMin",
                 "Tracks with an angle lambda (relative to the cathode) smaller than the given value will be ignored. Set to 0 to accept all particles.",
                 _cutLamMin,
                 float(0));

  registerProcessorParameter("cutLamMax",
                 "Tracks with an angle lambda (relative to the cathode) larger than the given value will be ignored. Set to 0 to accept all particles.",
                 _cutLamMax,
                 float(0));

}

void ComprehensivePIDProcessor::init() {
  // usually a good idea to
  printParameters();

  if (_modeTrain && _modeInfer)
  {
    sloE << "I cannot train and infer in the same process!" << std::endl;
    throw std::runtime_error("mode error");
  }

  _nEvt = 0;
  _nPFO = 0;

  _nAlgos  = _inputAlgoSpecs.size();
  _nModels = _trainModelSpecs.size();


  //AIDAProcessor::tree(this);
  if (_modeExtract)
  {
    if (_TTreeFileName!="")
    {
      _TTreeFile = new TFile(_TTreeFileName.c_str(), "RECREATE");
      _TTreeFile->cd();
    }
    _observablesTree = new TTree("observablesTree","Tree of all observable values");
  }
  else
  {
    _TTreeFile = new TFile(_TTreeFileName.c_str());
    _observablesTree = (TTree*)_TTreeFile->Get("observablesTree");
    //_observablesTree->Print();
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
  for (int i=0; i<_momNBins; ++i)
  {
    if (_momLog) _momBins.push_back(pow( 10, log10(_momMin) + (log10(_momMax)-log10(_momMin))*(i+1)/(float)(_momNBins) ));
    else _momBins.push_back(_momMin + (_momMax-_momMin)*(i+1)/(float)(_momNBins));
    sloM << " " << _momBins[i+1];
  }
  sloM << std::endl << "-------------------------------------------------" << std::endl;

  // prepare optional parameters of modules
  std::shared_ptr<StringParameters> pars = parameters();
  std::vector<float> inparF{};
  std::vector<std::string> inparS{};

  // set up and initialise input algorithms
  AlgorithmMgr::instance()->printAvailableAlgorithmTypes();

  for (int i=0; i<_nAlgos; ++i)
  {
    std::size_t p = _inputAlgoSpecs[i].find(":");
    if (p==std::string::npos)
    {
      _inputAlgoTypes.push_back(_inputAlgoSpecs[i]);
      _inputAlgoNames.push_back(_inputAlgoSpecs[i]);
    }
    else
    {
      _inputAlgoTypes.push_back(_inputAlgoSpecs[i].substr(0,p));
      _inputAlgoNames.push_back(_inputAlgoSpecs[i].substr(p+1));
    }

    //sloM << "creating new algorithm: " << _inputAlgoTypes[i] << ":" << _inputAlgoNames[i] << std::endl;

    _inputAlgorithms.push_back(AlgorithmMgr::instance()->createAlgorithm(_inputAlgoTypes[i]));
    _inputAlgorithms[i]->setName(_inputAlgoNames[i]);

    inparF.clear();
    inparS.clear();
    std::string sF = _inputAlgoNames[i]; sF.append(".F");
    std::string sS = _inputAlgoNames[i]; sS.append(".S");
    pars->getFloatVals(sF,inparF);
    pars->getStringVals(sS,inparS);

    std::vector<std::string> obsNames = _inputAlgorithms[i]->init(inparF, inparS);
    for (unsigned int j=0; j<obsNames.size(); ++j)
    {
      std::string s = _inputAlgoNames[i];
      s.append("_");
      _observablesNames.push_back(s.append(obsNames[j]));
    }

    sloM << "New algorithm created: " << _inputAlgorithms[i]->type() << ":" << _inputAlgorithms[i]->name() << "\n" <<  std::endl;
  }

  sloM << "Input algorithm inits done" << std::endl;
  sloM << "-------------------------------------------------" << std::endl;

  _nObs = _observablesNames.size();
  _observablesValues.resize(_nObs,0);

  // set up and initialise training/inference models
  ModelMgr::instance()->printAvailableModelTypes();

  if (_modeTrain || _modeInfer)
  {
    for (int i=0; i<_nModels; ++i)
    {
      std::size_t p = _trainModelSpecs[i].find(":");
      if (p==std::string::npos)
      {
        _trainModelTypes.push_back(_trainModelSpecs[i]);
        _trainModelNames.push_back(_trainModelSpecs[i]);
      }
      else
      {
        _trainModelTypes.push_back(_trainModelSpecs[i].substr(0,p));
        _trainModelNames.push_back(_trainModelSpecs[i].substr(p+1));
      }

      _trainingModels.push_back(ModelMgr::instance()->createModel(_trainModelTypes[i]));
      _trainingModels[i]->setName(_trainModelNames[i]);

      sloM << "New model created: " << _trainingModels[i]->type() << ":" << _trainingModels[i]->name() << "\n" <<  std::endl;

      if (_modeInfer)
      {
        _nMomBins.push_back(0);
        std::vector<std::pair<float,float> > none{};
        _weightFileBrackets.push_back(none);
        std::vector<std::string> empty{};
        _weightFiles.push_back(empty);
        ReadReferenceFile(i);
      }
    }

    if (_trainingObservables.size()==0)
    {
      _trainingObservables = _observablesNames;
      _trainingObservables.erase(_trainingObservables.begin()+2,_trainingObservables.begin()+6);
    }

    sloM << "Training/inference observables: " << _trainingObservables.size() << " -";
    for (unsigned int o=0; o<_trainingObservables.size(); ++o) {sloM << " " << _trainingObservables[o];}
    sloM << std::endl << "Extracted observables: " << _observablesNames.size() << " -";
    for (unsigned int o=0; o<_observablesNames.size(); ++o) {sloM << " " << _observablesNames[o];}
    sloM << std::endl;

    bool obs_ok = true;
    if (_modeExtract)
    {
      obs_ok = CheckTrainingObservables(_trainingObservables, _observablesNames);
      if (!obs_ok)
      {
        sloE << "Target training/inference observables not (completely) among extracted observables!" << std::endl;
        throw std::runtime_error("oberservables error");
      }
    }
    // todo: else: check if target obs are in root file

    if (_modeInfer)
    {
       for (std::string tObs : _trainingObservables)
        for (int p=0; p<_nObs; ++p)
          if (tObs == _observablesNames[p]) _inferObsValues.push_back(&(_observablesValues[p]));
    }

    TrainingModelInterface tmi;
    tmi.signalPDGs = _signalPDGs;
    tmi.backgroundPDGs = _backgroundPDGs;
    tmi.obsNames = _trainingObservables;
    tmi.momBins = _momBins;

    for (int i=0; i<_nModels; ++i)
    {
      inparF.clear();
      inparS.clear();
      std::string sF = _trainModelNames[i]; sF.append(".F");
      std::string sS = _trainModelNames[i]; sS.append(".S");
      pars->getFloatVals(sF,inparF);
      pars->getStringVals(sS,inparS);

      tmi.inparF = inparF;
      tmi.inparS = inparS;

      if (_modeTrain)
      {
        _weightFiles.push_back(_trainingModels[i]->initTraining(tmi));
        CreateReferenceFile(i);
      }

      if (_modeInfer)
      {
        tmi.obsValues = _inferObsValues;
        tmi.weightFiles = _weightFiles[i];
        _trainingModels[i]->initInference(tmi);
      }
    }

    sloM << "Model inits done" << std::endl;
    sloM << "-------------------------------------------------" << std::endl;
  }

  if (_modeExtract)
  {
    for (int iObs=0; iObs<_nObs; ++iObs)
    {
      std::stringstream s; s << _observablesNames[iObs].c_str() << "/F";
      _observablesTree-> Branch(_observablesNames[iObs].c_str(), &(_observablesValues[iObs]), s.str().c_str());
    }
    _observablesTree->Write();
  }

  sloM << "Branches done" << std::endl;
  sloM << "-------------------------------------------------" << std::endl;

  for (int i=0; i<_nReason; ++i) _nRejectedPFOs[i] = 0;
  _rejectedPFOs = new TH1D("rejectedPFOs","Number of rejected PFOs by reason",_nReason,.5,_nReason+.5);
  _rejectedPFOs->GetXaxis()->SetBinLabel(1,"no MCPDG found");
  _rejectedPFOs->GetXaxis()->SetBinLabel(2,"MCPDG not among signal PDGs");
  _rejectedPFOs->GetXaxis()->SetBinLabel(3,"no track found");
  _rejectedPFOs->GetXaxis()->SetBinLabel(4,"more than one track found");
  _rejectedPFOs->GetXaxis()->SetBinLabel(5,"too large d0");
  _rejectedPFOs->GetXaxis()->SetBinLabel(6,"too large z0");
  _rejectedPFOs->SetYTitle("abundance");

  for (int i=0; i<_nPart; ++i) _nWrongMCPDG[i] = 0;
  _wrongMCPDG = new TH1D("wrongMCPDG","Number of MCPDG by origin file",_nPart,.5,_nPart+.5);
  _wrongMCPDG->GetXaxis()->SetBinLabel(1,"electron file");
  _wrongMCPDG->GetXaxis()->SetBinLabel(2,"muon file");
  _wrongMCPDG->GetXaxis()->SetBinLabel(3,"pion file");
  _wrongMCPDG->GetXaxis()->SetBinLabel(4,"kaon file");
  _wrongMCPDG->GetXaxis()->SetBinLabel(5,"proton file");
  _wrongMCPDG->SetYTitle("abundance");

  int nSig = _signalPDGs.size();
  _PDGCheck = new TH2I("PDGCheck", "PDG cross check PFO vs. MCTruth",nSig+1,-.5,nSig+.5,nSig+1,-.5,nSig+.5);
  _PDGCheck->SetXTitle("MCTruth PDG");
  _PDGCheck->SetYTitle("PFO PDG");
  for (int i=0; i<nSig; i++)
  {
    std::stringstream s; s << _signalPDGs[i];
    _PDGCheck->GetXaxis()->SetBinLabel(i+1,s.str().c_str());
    _PDGCheck->GetYaxis()->SetBinLabel(i+1,s.str().c_str());
  }
  _PDGCheck->GetXaxis()->SetBinLabel(nSig+1,"other");
  _PDGCheck->GetYaxis()->SetBinLabel(nSig+1,"other");

  sloM << _PDGCheck->GetName() << std::endl;
}


void ComprehensivePIDProcessor::processRunHeader( LCRunHeader* ) {
  _nRun++ ;
  //_nEvt = 0;
}

void ComprehensivePIDProcessor::processEvent(LCEvent* evt) {

  //sloM << _nEvt << std::endl;

  if (_nAlgos==0) return;


  if (_modeExtract || _modeInfer)
  {
    LCCollection *col_pfo{}, *col_pfo2mc{};

    try
      {
        col_pfo = evt->getCollection( _PFOColName );
        col_pfo2mc = evt->getCollection( _RecoMCTruthLinkName );
      }
      catch(DataNotAvailableException &e)
      {
        sloM << "Input PFO or RecoMCTruthLink collection not found - skipping event " << _nEvt << std::endl;
        return;
      }

    int n_pfo = col_pfo->getNumberOfElements();
    LCRelationNavigator rel_pfo2mc(col_pfo2mc);

    PIDHandler pidh(col_pfo);
    std::vector<int> algoID{};
    std::vector<int> allPDGs{};
    if (_modeInfer)
    {
      std::vector<std::string> PDGness{};
      for (int pdg : _signalPDGs)
      {
        std::stringstream pn; pn << pdg << "-ness";
        PDGness.push_back(pn.str());
        allPDGs.push_back(pdg);
      }
      for (int pdg : _backgroundPDGs)
      {
        std::stringstream pn; pn << pdg << "-ness";
        PDGness.push_back(pn.str());
        allPDGs.push_back(pdg);
      }
      for (int m=0; m<_nModels; ++m) algoID.push_back(pidh.addAlgorithm(_trainModelNames[m], PDGness));
    }

    for (int i=0; i<n_pfo; ++i)
    {
      ReconstructedParticleImpl* pfo = dynamic_cast<ReconstructedParticleImpl*>(col_pfo->getElementAt(i));

      // get the true PDG  (same as in dEdxAnalyser)
      const LCObjectVec& mcparVec = rel_pfo2mc.getRelatedToObjects(pfo);
      const FloatVec& mcparWei = rel_pfo2mc.getRelatedToWeights(pfo);
      int MCPDG = 0;
      float bestwei = 0;
      int bestk = -1;
      for (unsigned k=0; k<mcparVec.size(); ++k)
      {
        MCParticleImpl* mcpar = dynamic_cast<MCParticleImpl*>(mcparVec[k]);
        float mcparwei = (int(mcparWei[k])%10000)/1000.;
        if (mcparwei > bestwei)
        {
          MCPDG = mcpar->getPDG();
          bestwei = mcparwei;
          bestk = k;
        }
      }

      if (bestk==-1) {++_nRejectedPFOs[0]; continue;}
      if (std::find(_signalPDGs.begin(),_signalPDGs.end(),abs(MCPDG)) == _signalPDGs.end() && std::find(_backgroundPDGs.begin(),_backgroundPDGs.end(),abs(MCPDG)) == _backgroundPDGs.end())
      {
        ++_nRejectedPFOs[1];
        ++_nWrongMCPDG[_nEvt/100000];
        continue;
      }
      if (pfo->getTracks().size()<1) {++_nRejectedPFOs[2]; continue;}
      //if (pfo->getTracks().size()>1) {++_nRejectedPFOs[3]; continue;}
      if (pfo->getTracks()[0]->getD0()>100) {++_nRejectedPFOs[4]; continue;}
      if (pfo->getTracks()[0]->getZ0()>100) {++_nRejectedPFOs[5]; continue;}

      ++_nPFO;

      if (_modeExtract)
      {

        const double* mom = pfo->getMomentum();
        _momabs = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);

        double tanl = mom[2]/sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
        _lambda = fabs(atan(tanl)) *180/M_PI;

        _d0 = pfo->getTracks()[0]->getD0();
        _z0 = pfo->getTracks()[0]->getZ0();

        _PDG = MCPDG;
        _nTracks = pfo->getTracks().size();

        _observablesValues[0]=_momabs;
        _observablesValues[1]=_lambda;
        _observablesValues[2]=_d0;
        _observablesValues[3]=_z0;
        _observablesValues[4]=_PDG;
        _observablesValues[5]=_nTracks;

        int iObs=_nObs_base;

        for (int a=0; a<_nAlgos; ++a)
        {
          //sloM << "using algorithm " << _inputAlgorithms[a]->type() << std::endl;
          std::vector<std::pair<float,float> > obs = _inputAlgorithms[a]->extractObservables(pfo, col_pfo);
          for (unsigned int j=0; j<obs.size(); ++j)
          {
            _observablesValues[iObs] = obs[j].first;
            ++iObs;
          }
        }
        sloD << "observables extracted " << std::endl;

        _observablesTree->Fill();
        sloD << "trees filled" << std::endl;
      }

      if (_modeInfer)
      {
        for (int m=0; m<_nModels; ++m) for (int j=0; j<_nMomBins[m]; ++j)
        {
          if (_momabs<_weightFileBrackets[m][j].first || _momabs>_weightFileBrackets[m][j].second) continue;
          const std::vector<float> eval_out = _trainingModels[m]->runInference(j);
          sloD << " > inference was run" << std::endl;

          auto eval_max = max_element(eval_out.begin(),eval_out.end());
          int index = std::distance(eval_out.begin(), eval_max);

          pidh.setParticleID(pfo, 0, allPDGs[index], 0, algoID[m], eval_out);
          for (int ipdg=0; ipdg<(int)_signalPDGs.size(); ++ipdg) if (_signalPDGs[ipdg]==abs(MCPDG)) _PDGCheck->Fill(ipdg, index);
          //sloM << MCPDG << " " << index << " " << allPDGs[index] << " " << _PDGCheck->GetEntries() << std::endl;
          break;
        }
      }
    } // for pfos

  }

  _nEvt++;

}

void ComprehensivePIDProcessor::check( LCEvent* ) {}
void ComprehensivePIDProcessor::end()
{
  sloM << "-------------------------------------------------" << std::endl << "end()" << std::endl;
  for (int i=0; i<_nReason; ++i) _rejectedPFOs->SetBinContent(i+1,_nRejectedPFOs[i]);
  for (int i=0; i<_nPart; ++i) _wrongMCPDG->SetBinContent(i+1,_nWrongMCPDG[i]);

  if (_modeExtract)
  {
    _TTreeFile->cd();
    //_observablesTree->Print();
    _observablesTree->Write();
  }

  if (_modeTrain)
  {
    sloM << "Begin training..." << std::endl << std::endl;
    _observablesTree->Print();
    for (int i=0; i<_nModels; ++i) _trainingModels[i]->runTraining(_observablesTree);
    sloM << "Training done" << std::endl << std::endl;
  }

  TCanvas* can = new TCanvas;
  //TImage* img = TImage::Create();
  gStyle->SetPalette(kBird);
  //can->SetGrid();
  //can->SetLogx();
  can->SetLogz();
  can->SetGrid(0,0);
  gStyle->SetOptStat(0);

  sloM << _PDGCheck->GetName() << std::endl;
  PlotTH2(can, _PDGCheck);

  _TTreeFile->Close();
  sloM << "Numer of Events: " << _nEvt << "    Numer of PFOs: " << _nPFO << std::endl;
}


bool ComprehensivePIDProcessor::CheckTrainingObservables(std::vector<std::string>& trainObs, std::vector<std::string>& compObs)
{
  bool check_ok = true;
  for (std::string tObs : trainObs)
  {
    bool found = false;
    for (std::string cObs : compObs)
      if (tObs == cObs) {found = true; break;}

    if (!found) {check_ok = false; break;}
  }
  return check_ok;
}

void ComprehensivePIDProcessor::CreateReferenceFile(int n)
{
  sloM << "Creating reference file..." << std::endl;
  std::string refname = ReferenceFile(n);
  std::ofstream reffile;
  reffile.open(refname);

  for (std::string obs : _trainingObservables) {reffile << obs << " ";}
  reffile << "\n";

  for (int pdg : _signalPDGs) {reffile << pdg << " ";}
  reffile << "| ";
  for (int pdg : _backgroundPDGs) {reffile << pdg << " ";}
  reffile << "\n\n";

  for (int i=0; i<_momNBins; ++i) {reffile << _momBins[i] << " " << _momBins[i+1] << " " << _weightFiles[n][i] << "\n";}
  reffile.close();
  sloM << "Reference file " << refname << " created." << std::endl;
}

void ComprehensivePIDProcessor::ReadReferenceFile(int n)
{
  std::string refname = ReferenceFile(n);
  std::vector<std::string> lineV;
  std::string line;
  std::ifstream reffile (refname);
  if (reffile.is_open())
  {
    getline(reffile,line);
    boost::split(lineV, line, [](char c){return c == ' ';});
    if (lineV.back()=="") lineV.pop_back();
    if (n==0) _trainingObservables = lineV;
    if (n>0 and lineV.size()!=_trainingObservables.size())
    {
      sloE << "Number of target observables not consistent among reference files!" << std::endl;
      throw std::runtime_error("oberservables error");
    }
    if (n>0 and not CheckTrainingObservables(_trainingObservables, lineV))
    {
      sloE << "Specified target observables not consistent among reference files!" << std::endl;
      throw std::runtime_error("oberservables error");
    }

    getline(reffile,line);
    if (n==0)
    {
      boost::split(lineV, line, [](char c){return c == ' ';});
      for (std::string l : lineV)
        if (l=="|") _SnB.push_back(0);
        else if (l!="") _SnB.push_back(stoi(l));

    }
    // todo: make also SnB safe against incoherent reference files

    while (getline(reffile,line))
    {
      boost::split(lineV, line, [](char c){return c==' ';});
      if (lineV.size()<=1) continue;
      std::pair<float,float> mombracket(std::stof(lineV[0]),std::stof(lineV[1]));
      _weightFileBrackets[n].push_back(mombracket);
      _weightFiles[n].push_back(lineV[2]);
    }
    reffile.close();
    _nMomBins[n] = _weightFileBrackets[n].size();
  }
  else
  {
    sloE << "Reference file " << refname << " could not be opened!" << std::endl;
    throw std::runtime_error("reffile error");
  }
  sloM << "Reference file successfully read." << std::endl;
}

std::string ComprehensivePIDProcessor::ReferenceFile(int n)
{
  std::string refname{};
  if (_nModels > int(_reffile.size()))
  {
    refname = _reffile[0];
    std::size_t p = refname.rfind(".");
    refname.insert(p, std::string("_")+std::to_string(n));
  }
  else refname = _reffile[n];
  return refname;
}

void ComprehensivePIDProcessor::PlotTH2(TCanvas* can, TH2* hist)
{
  hist->Draw("colz");
  can->Update();
  std::stringstream s; s << _plotFolder << "/" << hist->GetName() << _fileFormat;
  can->Print(s.str().c_str());
}
