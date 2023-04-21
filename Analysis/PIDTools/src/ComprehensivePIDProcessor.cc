#include <ComprehensivePIDProcessor.h>

#include <EVENT/LCCollection.h>
#include <UTIL/LCRelationNavigator.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/MCParticleImpl.h>
#include <InputAlgorithm.h>
//#include <InputAlgorithm_TOF.h>
#include <AlgorithmMgr.h>
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

  registerProcessorParameter("operationMode",
                 "Operation mode of processor; 1 = extract observables, 2 = apply PID to PFOs; default: 1.",
                 _mode,
                 int(1));

  registerProcessorParameter("mode2nBins",
                 "Number of momentum bins in operation mode 2; default: 12.",
                 _mode2nBins,
                 int(12));

  registerProcessorParameter("weightsfolder",
                 "Folder were weights are or should be stored; default: ./dataset/weights.",
                 _weightsfolder,
                 std::string("./dataset/weights"));

  registerProcessorParameter("reffile",
                 "Reference file; default: Ref.txt.",
                 _reffile,
                 std::string("Ref.txt"));

  std::vector<int> signalPDGs = {11,13,211,321,2122};
  registerProcessorParameter("signalPDGs",
                 "Vector of PDG numbers that are considered signal; default: {11,13,211,321,2212}.",
                 _signalPDGs,
                 signalPDGs);

  registerProcessorParameter("backgroundPDGs",
                 "Vector of PDG numbers that are considered signal; default: {}.",
                 _backgroundPDGs,
                 std::vector<int>());

  //std::vector<std::string > inputalgos{};
  registerProcessorParameter("inputAlgoSpecs",
                 "List of input algorithms; for each specify only type (name is then type) or type:name; default: {}",
                 _inputAlgoSpecs,
                 std::vector<std::string>());

  registerOptionalParameter("inputAlgo1ParamsD",
                 "Vector of doubles that are input parameters for InputAlgorithm 1; default: {}.",
                 _inputAlgo1ParamsF,
                 std::vector<float>());

  //std::vector<std::string> inputAlgo1ParamsS = {};
  registerOptionalParameter("inputAlgo1ParamsS",
                 "Vector of strings that are input parameters for InputAlgorithm 1; default: {}.",
                 _inputAlgo1ParamsS,
                 std::vector<std::string>());

}

void ComprehensivePIDProcessor::init() {
  // usually a good idea to
  printParameters();

  _nEvt = 0;
  _nPFO = 0;

  _nAlgos = _inputAlgoSpecs.size();
  if (_mode==1 && _nAlgos==0) {
    streamlog_out(WARNING) << "No input algorithms given, will not do anything!" << std::endl;
    return;
  }

  AlgorithmMgr::instance()->printAvailableAlgorithmTypes();
  sloM << "===================" << std::endl;

  AIDAProcessor::tree(this);

  _observablesTree = new TTree("observablesTree","Tree of all observable values");
  _observablesTrees.push_back(_observablesTree);
  _nObs = 0;

  _observablesNames.push_back("momabs");
  _observablesNames.push_back("lambda");
  _observablesNames.push_back("d0");
  _observablesNames.push_back("z0");
  _observablesNames.push_back("PDG");
  _observablesNames.push_back("nTracks");
  _nObs_base = 6;

  std::shared_ptr<StringParameters> pars = parameters();
  std::vector<float> inparF{};
  std::vector<std::string> inparS{};

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

    std::vector<std::string> obsNames = _inputAlgorithms[i]->init(inparF,inparS);
    for (unsigned int j=0; j<obsNames.size(); ++j)
    {
      std::string s = _inputAlgoNames[i];
      s.append("_");
      _observablesNames.push_back(s.append(obsNames[j]));
    }


    sloM << "new algorithm created: " << _inputAlgorithms[i]->type() << ":" << _inputAlgorithms[i]->name() << "\n" <<  std::endl;
  }
  sloM << "inits done" << std::endl;

  _nObs = _observablesNames.size();
  _observablesValues.resize(_nObs,0);

  for (int iObs=0; iObs<_nObs; ++iObs) RegisterBranchF(_observablesTrees, &(_observablesValues[iObs]), _observablesNames[iObs].c_str());

  sloM << "branches done" << std::endl;
  sloM << "===================" << std::endl;

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

  if (_mode==2)
  {
    _TReader = new TMVA::Reader( "!Color:!Silent" );

    std::vector<std::string> obslist;
    std::vector<std::string> lineV;
    std::string line;
    std::ifstream reffile (_reffile.c_str());
    if (reffile.is_open())
    {
      getline(reffile,line);
      boost::split(obslist, line, [](char c){return c == ' ';});

      for (unsigned int o=0; o<obslist.size(); ++o)
        for (int p=0; p<_nObs; ++p)
          if (obslist[o] == _observablesNames[p]) _TReader->AddVariable(_observablesNames[p], &(_observablesValues[p]));

      int i=0;
      while (getline(reffile,line))
      {
        boost::split(lineV, line, [](char c){return c==' ';});
        if (lineV.size()<=1) continue;
        std::pair<float,float> mombracket(std::stof(lineV[0]),std::stof(lineV[1]));
        _weightFiles.push_back(mombracket);
        std::stringstream nm; nm << "MomBin_" << i;
        _TReader->BookMVA(nm.str(), lineV[2]);
        ++i;
      }
      reffile.close();
      _mode2nBins = i;
    }
  }

}

void ComprehensivePIDProcessor::processRunHeader( LCRunHeader* ) {
  _nRun++ ;
  //_nEvt = 0;
}

void ComprehensivePIDProcessor::processEvent(LCEvent* evt) {

  //sloM << _nEvt << std::endl;

  if (_nAlgos==0) return;

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
  int algoID=0;
  if (_mode==2)
  {
    std::vector<std::string> PDGness{};
    for (unsigned int s=0; s<_signalPDGs.size(); ++s)
    {
      std::stringstream ss; ss << _signalPDGs[s] << "-ness";
      PDGness.push_back(ss.str());
    }
    algoID = pidh.addAlgorithm("dEdx_RCD__Pandora", PDGness);
  }

  for (int i=0; i<n_pfo; ++i)
  {
    //sloM << " i " << i << std::endl;
    ReconstructedParticleImpl* pfo = dynamic_cast<ReconstructedParticleImpl*>(col_pfo->getElementAt(i));
    //sloM << " dc 1 done " << i << std::endl;

    // get the true PDG
    const LCObjectVec& mcparVec = rel_pfo2mc.getRelatedToObjects(pfo);
    const FloatVec& mcparWei = rel_pfo2mc.getRelatedToWeights(pfo);
    int MCPDG = 0;
    float bestwei = 0;
    int bestk = -1;
    for (unsigned k=0; k<mcparVec.size(); ++k)
    {
      MCParticleImpl* mcpar = dynamic_cast<MCParticleImpl*>(mcparVec[k]);
      float mcparwei = mcparWei[k];
      if (mcparwei > bestwei)
      {
        MCPDG = mcpar->getPDG();
        bestwei = mcparwei;
        bestk = k;
      }
    }

    if (bestk==-1) {++_nRejectedPFOs[0]; continue;}
    if (std::find(_signalPDGs.begin(),_signalPDGs.end(),abs(MCPDG)) == _signalPDGs.end())
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

    //sloM << "observables extracted " << std::endl;

    if (_mode==1)
    {
      _observablesTree->Fill();
      //sloM << "trees filled" << std::endl;
    }

    if (_mode==2)
    {
      for (int j=0; j<_mode2nBins; ++j)
      {

        if (_momabs<_weightFiles[j].first || _momabs>_weightFiles[j].second) continue;
        std::stringstream nm; nm << "MomBin_" << j;
        const std::vector<float> eval_out = _TReader->EvaluateMulticlass(nm.str());

        auto eval_max = max_element(eval_out.begin(),eval_out.end());
        int index = std::distance(eval_out.begin(), eval_max);

//        int x = 0;
//        if (index != 0)
//        {
//          sloM << x << " "; ++x;
//          sloM << index << " ";
//          sloM << _signalPDGs[index] << " == ";
//            for (unsigned int e=0; e<eval_out.size(); ++e)
//          {
//            sloM << eval_out[e] << " ";
//          }
//          //sloM << eval_max << " ";
//
//          sloM << std::endl;
//        }


        pidh.setParticleID(pfo, 0, _signalPDGs[index], 0, algoID, eval_out);
        break;
      }
      //else if (_momabs>1) {sloM << bin << " " << _momabs << std::endl;}
    }
    //sloM << " ex done " << i << std::endl;
    //sloM << _nEvt << " " << i << " " << obs1[0].first << " " << obs1[0].second << std::endl;
    //sloM << _nEvt << " " << i << " " << obs2[0].first << " " << obs2[0].second << std::endl;
  }

  _nEvt++;

}

void ComprehensivePIDProcessor::check( LCEvent* evt ) {if (evt){}}
void ComprehensivePIDProcessor::end()
{
  for (int i=0; i<_nReason; ++i) _rejectedPFOs->SetBinContent(i+1,_nRejectedPFOs[i]);
  for (int i=0; i<_nPart; ++i) _wrongMCPDG->SetBinContent(i+1,_nWrongMCPDG[i]);

  sloM << "Numer of Events: " << _nEvt << "    Numer of PFOs: " << _nPFO << std::endl;
}

void ComprehensivePIDProcessor::RegisterBranchD(std::vector<TTree*>& treeVec, double* var, const char* name)
{
  std::stringstream s; s << name << "/D";
  for (unsigned int i=0; i<treeVec.size(); ++i) treeVec[i]->Branch(name, var, s.str().c_str());
}

void ComprehensivePIDProcessor::RegisterBranchF(std::vector<TTree*>& treeVec, float* var, const char* name)
{
  std::stringstream s; s << name << "/F";
  for (unsigned int i=0; i<treeVec.size(); ++i) treeVec[i]->Branch(name, var, s.str().c_str());
}

void ComprehensivePIDProcessor::RegisterBranchI(std::vector<TTree*>& treeVec, int* var, const char* name)
{
  std::stringstream s; s << name << "/I";
  for (unsigned int i=0; i<treeVec.size(); ++i) treeVec[i]->Branch(name, var, s.str().c_str());
}




