#include <string>
#include <vector>

#include <EVENT/LCCollection.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/PIDHandler.h>

#include "TLorentzVector.h"

#include "LikelihoodPID.hh"
#include "LikelihoodPIDProcessor.hh"
#include "LowMomentumMuPiSeparationPID_BDTG.hh"

LikelihoodPIDProcessor aLikelihoodPIDProcessor;

LikelihoodPIDProcessor::LikelihoodPIDProcessor() : Processor("LikelihoodPIDProcessor") {

  // Processor description
  _description = "Particle ID code using Bayesian Classifier";

  std::vector<std::string> xmlfiles;
  xmlfiles.push_back("TMVAClassification_BDTG_02GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_03GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_04GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_05GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_06GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_07GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_08GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_09GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_10GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_11GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_12GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_13GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_14GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_15GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_16GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_17GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_18GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_19GeVP_clusterinfo.weights.xml");
  xmlfiles.push_back("TMVAClassification_BDTG_20GeVP_clusterinfo.weights.xml");

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "RecoParticleCollection",
                          "Input collection of Reconstrcuted Particle", _inputPFOsCollection,
                          std::string("PandoraPFOs"));

  registerProcessorParameter("EnergyBoundaries", "Boundaries for energy", _energyBoundary, EVENT::FloatVec(0, 1.0e+07));

  registerProcessorParameter("FilePDFName", "rootfile of PDF", _PDFName, std::string("pdf_ParticleID_ok.root"));

  registerProcessorParameter("FileWeightFormupiSeparationName", "weight file for low momentum mu pi separation",
                             _weightFileName, xmlfiles);

  std::vector<float> parselectron, parsmuon, parspion, parskaon, parsproton;
  parselectron.push_back(-2.40638e-03);
  parselectron.push_back(7.10337e-01);
  parselectron.push_back(2.87718e-01);
  parselectron.push_back(-1.71591e+00);
  parselectron.push_back(0.0);
  registerProcessorParameter("dEdxParameter_electron", "dEdx Parameters for electron", _dEdxParamsElectron,
                             parselectron);

  parsmuon.push_back(8.11408e-02);
  parsmuon.push_back(9.92207e-01);
  parsmuon.push_back(7.58509e+05);
  parsmuon.push_back(-1.70167e-01);
  parsmuon.push_back(4.63670e-04);
  registerProcessorParameter("dEdxParameter_muon", "dEdx Parameters for muon", _dEdxParamsMuon, parsmuon);

  parspion.push_back(8.10756e-02);
  parspion.push_back(-1.45051e+06);
  parspion.push_back(-3.09843e+04);
  parspion.push_back(2.84056e-01);
  parspion.push_back(3.38131e-04);
  registerProcessorParameter("dEdxParameter_pion", "dEdx Parameters for pion", _dEdxParamsPion, parspion);

  parskaon.push_back(7.96117e-02);
  parskaon.push_back(4.13335e+03);
  parskaon.push_back(1.13577e+06);
  parskaon.push_back(1.80555e-01);
  parskaon.push_back(-3.15083e-04);
  registerProcessorParameter("dEdxParameter_kaon", "dEdx Parameters for Kaon", _dEdxParamsKaon, parskaon);

  parsproton.push_back(7.78772e-02);
  parsproton.push_back(4.49300e+04);
  parsproton.push_back(9.13778e+04);
  parsproton.push_back(1.50088e-01);
  parsproton.push_back(-6.64184e-04);
  registerProcessorParameter("dEdxParameter_proton", "dEdx Parameters for proton", _dEdxParamsProton, parsproton);

  registerProcessorParameter("dEdxNormalization", "dEdx Normalization", _dEdxNormalization, float(1.350e-7));

  registerProcessorParameter("dEdxErrorFactor", "dEdx Error factor", _dEdxErrorFactor, float(7.55));

  registerProcessorParameter("UseBayesian", "PID is based on Bayesian or Likelihood", _UseBayes, int(2));

  std::vector<float> costmat;
  costmat.push_back(0.0);
  costmat.push_back(1.0);
  costmat.push_back(1.0);
  costmat.push_back(1.0);
  costmat.push_back(1.0);

  costmat.push_back(1.0);
  costmat.push_back(0.0);
  costmat.push_back(9.0);
  costmat.push_back(1.0);
  costmat.push_back(1.0);

  costmat.push_back(1.0);
  costmat.push_back(1.0);
  costmat.push_back(0.0);
  costmat.push_back(1.0);
  costmat.push_back(1.0);

  costmat.push_back(1.0);
  costmat.push_back(1.0);
  costmat.push_back(30.0);
  costmat.push_back(0.0);
  costmat.push_back(1.0);

  costmat.push_back(1.0);
  costmat.push_back(1.0);
  costmat.push_back(1.0);
  costmat.push_back(10.0);
  costmat.push_back(0.0);

  registerProcessorParameter("CostMatrix", "cost matrix for risk based bayesian classifier", _cost, costmat);

  registerProcessorParameter("UseLowMomentumMuPiSeparation", "MVA mu/pi separation should be used or not", _UseMVA,
                             bool(true));

  std::vector<std::string> methods_to_run;
  methods_to_run.push_back("BasicVariablePID");
  methods_to_run.push_back("LowMomMuID");
  methods_to_run.push_back("ShowerShapesPID");
  methods_to_run.push_back("dEdxPID");
  methods_to_run.push_back("LikelihoodPID");

  registerProcessorParameter(
      "PIDMethodsToRun",
      "methods to be run, default: BasicVariablePID LowMomMuID ShowerShapesPID dEdxPID LikelihoodPID", _methodstorun,
      methods_to_run);

  std::string methods_to_run_version = "";
  registerProcessorParameter("PIDMethodsToRun_version", "version of the methods (for rerunning purposes)",
                             _methodstorun_version, methods_to_run_version);
}

void LikelihoodPIDProcessor::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;

  // dEdx parameters
  double pars[28];
  for (int i = 0; i < 5; i++)
    pars[i] = _dEdxParamsElectron[i];
  for (int i = 0; i < 5; i++)
    pars[i + 5] = _dEdxParamsMuon[i];
  for (int i = 0; i < 5; i++)
    pars[i + 10] = _dEdxParamsPion[i];
  for (int i = 0; i < 5; i++)
    pars[i + 15] = _dEdxParamsKaon[i];
  for (int i = 0; i < 5; i++)
    pars[i + 20] = _dEdxParamsProton[i];
  pars[25] = (double)_UseBayes;
  pars[26] = (double)_dEdxNormalization;
  pars[27] = (double)_dEdxErrorFactor;

  _pdgTable.push_back(11);
  _pdgTable.push_back(13);
  _pdgTable.push_back(211);
  _pdgTable.push_back(321);
  _pdgTable.push_back(2212);

  // for parameters except for dEdxPID
  _particleNames.push_back("electronLikelihood");
  _particleNames.push_back("muonLikelihood");
  _particleNames.push_back("pionLikelihood");
  _particleNames.push_back("kaonLikelihood");
  _particleNames.push_back("protonLikelihood");
  _particleNames.push_back("hadronLikelihood");
  _particleNames.push_back("MVAOutput_mupiSeparation");
  _particleNames.push_back("electronProbability");
  _particleNames.push_back("muonProbability");
  _particleNames.push_back("pionProbability");
  _particleNames.push_back("kaonProbability");
  _particleNames.push_back("protonProbability");
  _particleNames.push_back("hadronProbability");
  _particleNames.push_back("electron_dEdxdistance");
  _particleNames.push_back("muon_dEdxdistance");
  _particleNames.push_back("pion_dEdxdistance");
  _particleNames.push_back("kaon_dEdxdistance");
  _particleNames.push_back("proton_dEdxdistance");

  // for parameters of dEdxPID
  _dEdxNames.push_back("electronLikelihood");
  _dEdxNames.push_back("muonLikelihood");
  _dEdxNames.push_back("pionLikelihood");
  _dEdxNames.push_back("kaonLikelihood");
  _dEdxNames.push_back("protonLikelihood");
  _dEdxNames.push_back("hadronLikelihood");
  _dEdxNames.push_back("MVAOutput_mupiSeparation");
  _dEdxNames.push_back("electronProbability");
  _dEdxNames.push_back("muonProbability");
  _dEdxNames.push_back("pionProbability");
  _dEdxNames.push_back("kaonProbability");
  _dEdxNames.push_back("protonProbability");
  _dEdxNames.push_back("hadronProbability");
  _dEdxNames.push_back("electron_dEdxdistance");
  _dEdxNames.push_back("muon_dEdxdistance");
  _dEdxNames.push_back("pion_dEdxdistance");
  _dEdxNames.push_back("kaon_dEdxdistance");
  _dEdxNames.push_back("proton_dEdxdistance");

  _myPID = new LikelihoodPID(_PDFName, pars, _cost);

  // mupi separation class
  _mupiPID = new LowMomentumMuPiSeparationPID_BDTG(_weightFileName);

  printParameters();

  // CHECK that init settings are correct
  bool allnamecorrect = false;
  for (size_t i = 0; i < _methodstorun.size(); i++) {
    {
      if (_methodstorun.at(i).compare("dEdxPID") == 0 || _methodstorun.at(i).compare("LowMomMuID") == 0 ||
          _methodstorun.at(i).compare("BasicVariablePID") == 0 || _methodstorun.at(i).compare("ShowerShapesPID") == 0 ||
          _methodstorun.at(i).compare("LikelihoodPID") == 0)
        allnamecorrect = true;

      if (allnamecorrect == false) {
        throw EVENT::Exception(_methodstorun.at(i) +
                               std::string(" is not in the list of valid methods: BasicVariablePID LowMomMuID "
                                           "ShowerShapesPID dEdxPID LikelihoodPID"));
      }
      allnamecorrect = false;
    }
  }
}

void LikelihoodPIDProcessor::processRunHeader(LCRunHeader* /*run*/) {}

void LikelihoodPIDProcessor::processEvent(LCEvent* evt) {
  _pfoCol = evt->getCollection(_inputPFOsCollection);

  int npfo = _pfoCol->getNumberOfElements();
  PIDHandler pidh(_pfoCol); // BasicPID

  int algoID[10];
  for (size_t i = 0; i < _methodstorun.size(); i++) {
    if (_methodstorun.at(i).compare("dEdxPID") == 0)
      algoID[i] = pidh.addAlgorithm(_methodstorun.at(i) + _methodstorun_version, _dEdxNames);
    else
      algoID[i] = pidh.addAlgorithm(_methodstorun.at(i) + _methodstorun_version, _particleNames);
  }

  for (int i = 0; i < npfo; i++) {
    ReconstructedParticleImpl* part = dynamic_cast<ReconstructedParticleImpl*>(_pfoCol->getElementAt(i));

    if (part->getCharge() == 0)
      continue; // avoid neutral particle

    EVENT::ClusterVec clu = part->getClusters();
    lcio::Track* trk = part->getTracks()[0];
    TLorentzVector pp(part->getMomentum()[0], part->getMomentum()[1], part->getMomentum()[2], part->getEnergy());

    Int_t parttype = -1;
    //////////////////////////////////////////////////////////////////////////////////
    for (size_t ialgo = 0; ialgo < _methodstorun.size(); ialgo++) {
      if (_methodstorun.at(ialgo) == "LowMomMuID") {
        // Low momentum Muon identification
        //  (from 0.2 GeV until 2 GeV)
        if (_UseMVA) {
          Float_t MVAoutput = -100.0;
          parttype = -1;
          if (pp.P() < 2.5) {
            parttype = _mupiPID->MuPiSeparation(pp, trk, clu);
            MVAoutput = _mupiPID->getMVAOutput();
          }
          // create PIDHandler
          createParticleIDClass(parttype, part, pidh, algoID[ialgo], MVAoutput);
        }
      } else {
        // several partivle IDs performed
        // use just basic variables
        if (_methodstorun.at(ialgo) == "BasicVariablePID") {
          _myPID->setBasicFlg(true);
          _myPID->setdEdxFlg(false);
          _myPID->setShowerShapesFlg(false);
        }

        if (_methodstorun.at(ialgo).compare("dEdxPID") == 0) {
          _myPID->setBasicFlg(false);
          _myPID->setdEdxFlg(true);
          _myPID->setShowerShapesFlg(false);
        }

        if (_methodstorun.at(ialgo).compare("ShowerShapesPID") == 0) {
          _myPID->setBasicFlg(false);
          _myPID->setdEdxFlg(false);
          _myPID->setShowerShapesFlg(true);
        }

        if (_methodstorun.at(ialgo).compare("LikelihoodPID") == 0) {
          _myPID->setBasicFlg(true);
          _myPID->setdEdxFlg(true);
          _myPID->setShowerShapesFlg(true);
        }

        parttype = _myPID->Classification(pp, trk, clu);
        if (parttype < 0)
          parttype = 2;

        // mu-pi Separation for very low momentum tracks (from 0.2 GeV until 2 GeV)
        Float_t MVAoutput = -100.0;
        if ((parttype == 1 || parttype == 2) && (_UseMVA && pp.P() < 2.0 && clu.size() != 0)) {
          int tmpparttype = _mupiPID->MuPiSeparation(pp, trk, clu);
          if (_mupiPID->isValid())
            parttype = tmpparttype;
          MVAoutput = _mupiPID->getMVAOutput();
        }

        // create PIDHandler
        createParticleIDClass(parttype, part, pidh, algoID[ialgo], MVAoutput);
      }
    }

    /*ParticleIDImpl *PIDImpl=new ParticleIDImpl();
    if(parttype==0){
      PIDImpl->setType(11);
      PIDImpl->setLikelihood((float) posterior[0]);
    }
    if(parttype==1){
      PIDImpl->setType(13);
      PIDImpl->setLikelihood((float) posterior[1]);
    }
    if(parttype==2){
      PIDImpl->setType(211);
      PIDImpl->setLikelihood((float) posterior[2]);
    }
    if(parttype==3){
      PIDImpl->setType(321);
      PIDImpl->setLikelihood((float) posterior[3]);
    }
    if(parttype==4){
      PIDImpl->setType(2212);
      PIDImpl->setLikelihood((float) posterior[4]);
    }

    //posterior probability
    PIDImpl->addParameter((float) posterior[0]);   //electron
    PIDImpl->addParameter((float) posterior[1]);   //muon
    PIDImpl->addParameter((float) posterior[2]);   //pion
    PIDImpl->addParameter((float) posterior[3]);   //kaon
    PIDImpl->addParameter((float) posterior[4]);   //proton

    //add particle ID
    part->addParticleID(PIDImpl);
    //FloatVec chk=part->getParticleIDs()[0]->getParameters();
    //std::cout << chk[0] << std::endl;
    //add to PFO collection(so far, this is only PID)
    part->setParticleIDUsed(PIDImpl);
    part->setGoodnessOfPID(PIDImpl->getLikelihood());*/
  }
}

void LikelihoodPIDProcessor::check(LCEvent* /*evt*/) {}

void LikelihoodPIDProcessor::end() {
  delete _myPID;
  delete _mupiPID;
}

void LikelihoodPIDProcessor::createParticleIDClass(int parttype, ReconstructedParticle* part, PIDHandler& pidh,
                                                   int algoID, float MVAoutput) {

  Double_t* posterior = _myPID->GetPosterior();
  Double_t* likelihood = _myPID->GetLikelihood();
  std::vector<float> likelihoodProb;

  if (pidh.getAlgorithmName(algoID).find("LowMomMuID") == std::string::npos) {
    for (int j = 0; j < 6; j++)
      likelihoodProb.push_back(likelihood[j]);
    likelihoodProb.push_back(MVAoutput);
    for (int j = 0; j < 6; j++)
      likelihoodProb.push_back(posterior[j]);
  } else {
    for (int j = 0; j < 6; j++)
      likelihoodProb.push_back(999.0);
    likelihoodProb.push_back(MVAoutput);
    for (int j = 0; j < 6; j++)
      likelihoodProb.push_back(0.0);
    for (int j = 0; j < 6; j++) {
      likelihood[j] = 999.0;
      posterior[j] = 0.0;
    }
  }

  // for dEdx PID
  if (pidh.getAlgorithmName(algoID).find("dEdxPID") != std::string::npos ||
      pidh.getAlgorithmName(algoID).find("LikelihoodPID") != std::string::npos) {
    likelihoodProb.push_back((float)_myPID->get_dEdxDist(0)); // electron hypothesis
    likelihoodProb.push_back((float)_myPID->get_dEdxDist(1)); // muon hypothesis
    likelihoodProb.push_back((float)_myPID->get_dEdxDist(2)); // pion hypothesis
    likelihoodProb.push_back((float)_myPID->get_dEdxDist(3)); // kaon hypothesis
    likelihoodProb.push_back((float)_myPID->get_dEdxDist(4)); // proton hypothesis

    // std::cout << "check dedx: " << parttype << " " << algoID << " "
    //	      << likelihoodProb[11] << " " << likelihoodProb[12] << " " << likelihoodProb[13] << " " <<
    // likelihoodProb[14] << " " << likelihoodProb[15] << std::endl;
  } else {
    likelihoodProb.push_back((float)0.0); // electron hypothesis
    likelihoodProb.push_back((float)0.0); // muon hypothesis
    likelihoodProb.push_back((float)0.0); // pion hypothesis
    likelihoodProb.push_back((float)0.0); // kaon hypothesis
    likelihoodProb.push_back((float)0.0); // proton hypothesis
  }

  // std::cout << "check posterior: " << posterior[0] << " "
  //	    << posterior[1] << " "
  //	    << posterior[2] << " "
  //	    << posterior[3] << " "
  //	    << posterior[4] << " " << std::endl;

  // set pid results
  // set each hadron type
  if (pidh.getAlgorithmName(algoID).find("dEdxPID") != std::string::npos ||
      pidh.getAlgorithmName(algoID).find("LikelihoodPID") != std::string::npos)
    pidh.setParticleID(part, 0, _pdgTable[parttype], (float)likelihood[parttype], algoID, likelihoodProb);
  else { // ignore each hadron type
    int tmppart = parttype;
    parttype = std::min(2, parttype);
    if (parttype == 2)
      tmppart = 5;
    pidh.setParticleID(part, 0, _pdgTable[parttype], (float)likelihood[tmppart], algoID, likelihoodProb);
  }

  return;
}
