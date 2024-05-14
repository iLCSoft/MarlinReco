#include "WWCategorisationProcessor.h"

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <array> 

#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/LCEventImpl.h" 
#include "IMPL/LCCollectionVec.h"
#include "IMPL/MCParticleImpl.h" 
#include "IMPL/ReconstructedParticleImpl.h"

#include "TText.h"
#include "TStyle.h"
#include "TLorentzVector.h"

using namespace lcio;
using namespace marlin;


WWCategorisationProcessor  aWWCategorisationProcessor ;

WWCategorisationProcessor::WWCategorisationProcessor() : Processor("WWCategorisationProcessor") {
  
  _description = "Categorisation of all WW decay channels" ;
  
  registerInputCollection(LCIO::MCPARTICLE,
			   "MCParticles",
			   "Name of the MCParticle collection",
			   _MCParColName,
			   std::string("MCParticlesSkimmed"));
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
         "IsolatedElectrons",
         "Name of the ReconstructedParticle collection",
         _ElectronColName,
         std::string("IsolatedElectrons"));
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			   "IsolatedMuons",
			   "Name of the ReconstructedParticle collection",
			   _MuonColName,
			   std::string("IsolatedMuons"));
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			   "IsolatedTaus",
			   "Name of the ReconstructedParticle collection",
			   _TauColName,
			   std::string("IsolatedTaus"));
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "IsolatedPhotons",
			   "Name of the ReconstructedParticle collection",
			   _PhotonColName,
			   std::string("IsolatedPhotons"));
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			   "PFOsminusphoton",
			   "Name of the ReconstructedParticle collection",
			   _PFOsMinusIsolatedObjetcs,
			   std::string("PFOsminusphoton"));

  registerProcessorParameter("TTreeFileName",
          "Name of the root file in which the TTree with the 4 event observables is stored; if left empty no root file is created; default: TTreeFile.root",
          _TTreeFileName,
          std::string("TTreeFile.root"));
  registerProcessorParameter("ConfusionMatrixFileName",
          "Name of the png file in which the confusion matrix is stored; if left empty no file is created; default: ConfusionMatrix.png",
          _ConfusionMatrixFileName,
          std::string("ConfusionMatrix.png"));
}

void WWCategorisationProcessor::init(){
  // usually a good idea to
  printParameters() ;

  if (_TTreeFileName != ""){
    _doTT = true;
    _TTreeFile = new TFile(_TTreeFileName.c_str(), "RECREATE");
    _TTreeFile->cd();

    _observablesTree = new TTree("observablesTree","Tree of all observable values");

    _observablesTree->Branch("mInv",&_mInv,"mInv/F");
    _observablesTree->Branch("nPFO",&_nPFO,"nPFO/I");
    _observablesTree->Branch("misspT",&_misspT,"misspT/F");
    _observablesTree->Branch("missE",&_missE,"missE/F");

    _observablesTree->Branch("true_cat",&_true_cat,"true_cat_basic/I");
    _observablesTree->Branch("reco_cat_basic",&_reco_cat_basic,"reco_cat_basic/I");
    _observablesTree->Branch("reco_cat_advanced",&_reco_cat_advanced,"reco_cat_advancedc/I");
  }

  if (_ConfusionMatrixFileName != "")
  {
    _doCM = true;
    _confusion_matrix = new TH2F("Attempt_ConfusionMatrix", " ", _n_cat, -0.5, _n_cat-.5, _n_cat, -0.5, _n_cat-.5);
    _confusion_matrix_ef = new TH2F("ConfusionMatrix_efficiency", " ", _n_cat, -0.5, _n_cat-.5, _n_cat, -0.5, _n_cat-.5);
    _confusion_matrix_ef->SetXTitle("True Event Category");
    _confusion_matrix_ef->SetYTitle("Reco Event Category");
    _confusion_matrix_ef->GetXaxis()->SetBinLabel(1,"hadronic");
    _confusion_matrix_ef->GetXaxis()->SetBinLabel(2,"#splitline{true e#nu qq}{#theta < 6.27 deg} "); // #splitline{first line}{second line} # use latex with #
    _confusion_matrix_ef->GetXaxis()->SetBinLabel(3,"#splitline{true e#nu qq}{#theta > 6.27 deg} ");
    _confusion_matrix_ef->GetXaxis()->SetBinLabel(4,"true #mu#nu qq");
    _confusion_matrix_ef->GetXaxis()->SetBinLabel(5,"true #tau#nu qq");
    _confusion_matrix_ef->GetXaxis()->SetBinLabel(6,"leptonic");
    _confusion_matrix_ef->GetXaxis()->SetBinLabel(7,"other");

    _confusion_matrix_ef->GetYaxis()->SetBinLabel(1,"hadronic");
    _confusion_matrix_ef->GetYaxis()->SetBinLabel(2,"#splitline{semileptonic}{invisible}");
    _confusion_matrix_ef->GetYaxis()->SetBinLabel(3,"semileptonic e");
    _confusion_matrix_ef->GetYaxis()->SetBinLabel(4,"semileptonic #mu");
    _confusion_matrix_ef->GetYaxis()->SetBinLabel(5,"semileptonic #tau");
    _confusion_matrix_ef->GetYaxis()->SetBinLabel(6,"leptonic");
    _confusion_matrix_ef->GetYaxis()->SetBinLabel(7,"other");
  }

 for (int i=0; i < _n_cat; ++i){  // number of events in each category
    _nTrue.push_back(0);
    _nReco.push_back(0);
  }
  
  _nRun = 0 ;
  _nEvt = 0 ;
}

void WWCategorisationProcessor::processRunHeader(LCRunHeader* ){
  _nRun++ ;
}

void WWCategorisationProcessor::processEvent( LCEvent*  evt) {

  streamlog_out(DEBUG) <<" processing event " << evt->getEventNumber() << "  in run " << evt->getRunNumber() << std::endl ;
  
  LCCollection *col_muon{}, *col_electron{}, *col_tau{}, *col_mcparticles{}, *col_pfominusisolatedthings{}, *col_photons{};

  try{
    col_mcparticles            = evt->getCollection(_MCParColName);
    col_electron               = evt->getCollection(_ElectronColName);
    col_muon                   = evt->getCollection(_MuonColName);
    col_tau                    = evt->getCollection(_TauColName);
    col_photons                = evt->getCollection(_PhotonColName);
    col_pfominusisolatedthings = evt->getCollection(_PFOsMinusIsolatedObjetcs);
  }
  catch(DataNotAvailableException &e){
    streamlog_out(MESSAGE) << "Input collections not found - skipping event " << _nEvt << std::endl;
    return;
  }

  // default category value: "other", works as cross check that nothing is forgotten
  int true_cat = _n_cat-1;
  int reco_cat = _n_cat-1;
  
  // get number of isolated leptons and photons -> used for basic categorisation
  int n_electrons = col_electron->getNumberOfElements();
  int n_muons = col_muon->getNumberOfElements();
  int n_taus = col_tau->getNumberOfElements();
  int n_photons = col_photons->getNumberOfElements();
  int n_leptons = n_electrons + n_muons + n_taus;


  // MCTruth categorisation

  // particle 5 is - supposed to be! - the incoming electron/positron after ISR
  MCParticle* mc = dynamic_cast <MCParticle*> (col_mcparticles->getElementAt(5));
  int incoming_pdg = abs(mc->getPDG()); // get PDG of incoming e- or e+ after ISR

  MCParticle* myelectron{}; // space for electron in sl channel for later check

  // If the incoming particle is an electron/positron, check size of the daughters array; 4,3 or 1. If anything else: print
  if(incoming_pdg == 11){
    std::vector<MCParticle*> daughters = mc->getDaughters();
    std::vector<MCParticle*> general_mcdaughters{}; // space for MC daughters
    std::vector<int> general_daughters_pdgs = {}; // space for their PDGs

    if (daughters.size() <= 3){ //if daughters are of length 3 or 2 
      for (unsigned int i = 0; i < daughters.size(); ++i){ //go through all daughters and get their PDG
        int daughters_pdg = abs(daughters[i]->getPDG());
        
        if (daughters_pdg > 20 && daughters_pdg < 30){ // make sure the incoming particle daughter is a boson
          std::vector<MCParticle*> boson_daughters = daughters[i]->getDaughters();

          for (unsigned int j = 0; j < boson_daughters.size(); ++j){ // go through the boson daughters and add them to daughters and pdg lists
            general_mcdaughters.push_back(boson_daughters[j]);
            general_daughters_pdgs.push_back(abs(boson_daughters.at(j)->getPDG()));
          }
        }
        else{ //including remaining daughters that are not from W or Z
          general_mcdaughters.push_back(daughters[i]); //including mc information of remaining daughters that are not from W or Z
          general_daughters_pdgs.push_back(daughters_pdg);  //including PDGs of remaining daughters that are not from W or Z 
        } //by the end we get all 4f daughters PDGs and all 4f daughters mc info
      }
    }
    else if(daughters.size() == 4){ //go through all daughters of e-/e+ and save their PDG, then use the function again to count and categorise the events
      for(unsigned int i = 0; i < daughters.size(); ++i){
        general_mcdaughters.push_back(daughters[i]);
        general_daughters_pdgs.push_back(abs(daughters[i]->getPDG()));
      }
    }else{ //PRINT OUT ACTUAL SIZE IF NOT ANY OF THE ONES MENTIONED
      streamlog_out(DEBUG) << "atual size of daughters: " << daughters.size() << std::endl;
    }

    int channel = CountingLeptons(general_daughters_pdgs); //which channel -> based of how many final state leptons we found
    if (channel == 0) true_cat = 0; // hadronic
    else if (channel == 2) true_cat = 5; // leptonic
    else if (channel == 1){ // semileptonic, do sub-categorisation
      streamlog_out(DEBUG) << "MC daughters: [ ";
      for (unsigned int k = 0; k < general_daughters_pdgs.size(); ++k){
        streamlog_out(DEBUG) << general_daughters_pdgs[k] << ", ";

        if (IsLepton(general_daughters_pdgs[k])){
        true_cat = _sl_subcat[general_daughters_pdgs[k]];
        if (general_daughters_pdgs[k] == 11) myelectron = general_mcdaughters[k];
        break;
        }
      }
    }
    streamlog_out(DEBUG) << "]" << std::endl;
    streamlog_out(DEBUG) << "How many leptons: " << channel << std::endl;
  }
  else{ // particle 5 is not e+/-, print for debugging
      streamlog_out(MESSAGE) << "MCParticle no. 5 PDG: " << incoming_pdg << std::endl;
  }

  // very forward (i.e. invisible) electrons in semileptonic channel from 'single-W' contribution go to separate true category
  if (true_cat == 2){
    std::vector<MCParticle*> deepest_electron{}; // sl electron after FSR -> getting daughters of the electron until we have a FSR electron (radiated a photon)
    bool last = false;
    while(!last){
      deepest_electron = myelectron->getDaughters();
      last = true;
      for(unsigned int w = 0; w < deepest_electron.size(); w++){
        int deepest_electron_pdg = abs(deepest_electron[w]->getPDG());
        if(((deepest_electron_pdg > 90) && (deepest_electron_pdg < 99)) || (deepest_electron_pdg == 11)){
          myelectron = deepest_electron[w];
          last = false;
          break;
        }
      }
    }
    TLorentzVector lv = TLorentzVector(myelectron->getMomentum(), myelectron->getEnergy());
    double mctheta = lv.Theta();
    double cosine_mctheta = fabs(cos(mctheta));
    streamlog_out(DEBUG) << "cosine(#theta): " << cosine_mctheta << std::endl;
    if(cosine_mctheta > 0.994) true_cat = 1;
  }

  // true category established
  _true_cat = true_cat;

  // collect properties for reconstructed event category

  LCCollection* cols[] = {col_electron, col_muon, col_tau, col_photons};
  int n_parts[] = {n_electrons, n_muons, n_taus, n_photons};
  double momtot_max[] = {0, 0, 0, 0};  // highest total momentum for each isolated species
  TLorentzVector lv_IsoCumul{};  // cumulative 4-momentum for all isolated particles

  for (int i = 0; i < 4; ++i){  // loop through isolated species
    for (int j = 0; j < n_parts[i]; ++j){ // loop through all isolated particles
    EVENT::ReconstructedParticle* par = dynamic_cast <EVENT::ReconstructedParticle*>(cols[i]->getElementAt(j));
    lv_IsoCumul += TLorentzVector(par->getMomentum(), par->getEnergy()); // sum 4-mom for all particles
    auto mom = par->getMomentum();
    double momtot = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
    if (momtot > momtot_max[i]) momtot_max[i] = momtot; // assign highest momentum for that species
    }
  }

  // number of non-isolated PFOs before and after overlay cut
  int n_nonIsoPFOs = col_pfominusisolatedthings->getNumberOfElements();
  _nPFO = 0;

  // cumulative 4-vectors with and without overlay
  TLorentzVector lv_nonIsoPFOs_overlay{};
  TLorentzVector lv_nonIsoPFOs_nooverlay{};

  for(int i = 0; i < n_nonIsoPFOs; i++){ // loop over non-isolated PFOs
    EVENT::ReconstructedParticle* remainingpfo = dynamic_cast <EVENT::ReconstructedParticle*>(col_pfominusisolatedthings->getElementAt(i));
    TLorentzVector lv_remainingpfo = TLorentzVector(remainingpfo->getMomentum(), remainingpfo->getEnergy());
    double logpT = log10(lv_remainingpfo.Perp());
    double sinTheta = sin(fabs(lv_remainingpfo.Theta()));

    // summing all 4-vectors
    lv_nonIsoPFOs_overlay += TLorentzVector(remainingpfo->getMomentum(), remainingpfo->getEnergy());

    if(logpT > -2*sinTheta + 0.2){ // get all PFOs above this line -> triangle cut
      lv_nonIsoPFOs_nooverlay += TLorentzVector(remainingpfo->getMomentum(), remainingpfo->getEnergy());
      _nPFO ++;
    }
  }

  // invariant mass after triangle cut
  _mInv = lv_nonIsoPFOs_nooverlay.M();

  // cross check: 4Vec of remaining PFOs
  streamlog_out(DEBUG) << "Remaining PFOs 4v: [ "; 
  for (double k = 0; k < 4; k++){
    streamlog_out(DEBUG) << lv_nonIsoPFOs_nooverlay[k] << ", ";
  }
  streamlog_out(DEBUG) << "]" << std::endl;

  // add isolated PFOs to calcuated missing pT and energy
  TLorentzVector full_event4v_general = lv_nonIsoPFOs_overlay + lv_IsoCumul;
  TLorentzVector full_event4v_trianglecut = lv_nonIsoPFOs_nooverlay + lv_IsoCumul;

  // take into account beam crossing angle, get missing momentum and energy
  TLorentzVector cm4v_initialconditions = TLorentzVector(1.75,0,0,250);  
  TLorentzVector missing_neutrino4v_general = full_event4v_general - cm4v_initialconditions; // without overlay cut
  TLorentzVector missing_neutrino4v_trianglecut = full_event4v_trianglecut - cm4v_initialconditions; // with triangle overlay cut

  // missing pT and energy after triangle overlay cut
  _misspT = missing_neutrino4v_trianglecut.Perp();
  _missE = -missing_neutrino4v_trianglecut[3];

  // cross check: 4Vec of missing energy/momentum (neutrinos)
  streamlog_out(DEBUG) << "neutrino4v_general: [ "; 
  for (double k = 0; k < 4; k++){
    streamlog_out(DEBUG) << missing_neutrino4v_general[k] << ", ";
  }
  streamlog_out(DEBUG) << "]" << std::endl;
  

  // basic reco categorisation
  if(n_leptons >= 2) reco_cat = 5;  // leptonic
  else if(n_leptons == 1){ // semileptonic
    if (n_electrons == 1) reco_cat = 2;
    if (n_muons == 1) reco_cat = 3;
    if (n_taus == 1) reco_cat = 4;
  }
  else if(n_leptons == 0){
    reco_cat = 1;  // invisible
    if (_mInv >= 140) reco_cat = 0;  // hadronic
  }

  _reco_cat_basic = reco_cat;

  // advanced reco categorisation
  if(reco_cat > 1 && (_nPFO <= 20 || _mInv <= 40)) reco_cat = 5; // any non-hadronic to leptonic
  if(_nPFO >= 100 || _mInv >= 140) reco_cat = 0; // any to hadronic

  if(reco_cat == 5 && (40 <= _mInv && _mInv <= 140)){  // leptonic to semileptonic
    if (momtot_max[0] > momtot_max[1] && momtot_max[0] > momtot_max[2]) reco_cat = 2; // sl electron
    if (momtot_max[1] > momtot_max[0] && momtot_max[1] > momtot_max[2]) reco_cat = 3; // sl muon
    if (momtot_max[2] > momtot_max[0] && momtot_max[2] > momtot_max[1]) reco_cat = 4; // sl tau
  }

  if(reco_cat == 0){ // hadronic to ...
    if(_nPFO <= 20) reco_cat = 5; // leptonic
    if((20 < _nPFO && _nPFO <= 40) || _misspT >= 40) reco_cat = 1; // semileptonic invisible
  }

  _reco_cat_advanced = reco_cat;


  // store in event parameters
  LCParameters& pars = evt->parameters();
  pars.setValue("WWCategorisation.TrueCat", _true_cat);
  pars.setValue("WWCategorisation.RecoCatBasic", _reco_cat_basic);
  pars.setValue("WWCategorisation.RecoCatAdvanced", _reco_cat_advanced);
  pars.setValue("WWCategorisation.mInv", _mInv);
  pars.setValue("WWCategorisation.nPFO", _nPFO);
  pars.setValue("WWCategorisation.misspT", _misspT);
  pars.setValue("WWCategorisation.missE", _missE);

  // store in TTree and/or confusion matrix
  if (_doTT) _observablesTree->Fill();
  if (_doCM) _confusion_matrix->Fill(true_cat, reco_cat);

  _nTrue[true_cat]++;
  _nReco[reco_cat]++;

  _nEvt ++;
}

void WWCategorisationProcessor::check(LCEvent *) {}

void WWCategorisationProcessor::end(){
  // efficiency confusion matrix
  for(unsigned int reco_idx = 1; reco_idx < _nReco.size()+1; reco_idx++){
    for(unsigned int true_idx = 1; true_idx < _nTrue.size()+1; true_idx++){
      double val = _nTrue[true_idx-1] != 0 ? _confusion_matrix->GetBinContent(true_idx,reco_idx)/_nTrue[true_idx-1] : 0;
      _confusion_matrix_ef->SetBinContent(true_idx, reco_idx, val);
    }
  }

  if (_doTT){
    _TTreeFile->cd();
    _observablesTree->Write();
  }

  if (_doCM) {
    TCanvas* can = new TCanvas;
    PlotConfusionMatrix(can, _confusion_matrix_ef);
    if (_doTT){
      _confusion_matrix->Write();
      _confusion_matrix_ef->Write();
      can->Write();
    }
  }

  if (_doTT) _TTreeFile->Close();
}

int WWCategorisationProcessor::CountingLeptons(std::vector<int> daughters_pdg_vector){
  int numberof_leptons = 0;
  for(long unsigned int i = 0; i < daughters_pdg_vector.size(); i++){
    if(( daughters_pdg_vector[i] == 11) || (daughters_pdg_vector[i] == 13) || (daughters_pdg_vector[i] == 15)){
      numberof_leptons++;
    }
  }
  return numberof_leptons;
}

bool WWCategorisationProcessor::IsLepton(int pdg){
  return (pdg==11 || pdg==13 || pdg==15);
}

void WWCategorisationProcessor::PlotConfusionMatrix(TCanvas* can, TH2* histogram){
  histogram->Draw("text colz");

  for  (long unsigned int i=1; i<_nTrue.size()+1; ++i){
    std::stringstream t; t << _nTrue.at(i-1);
    std::string t1 = t.str(); t1.resize(6);
    TText* text2 = new TText(i -1.15 , _n_cat-.45 ,t1.c_str());
    text2->SetTextSize(0.035);
    text2->Draw();
  }

  gStyle->SetOptStat(0);
  histogram->SetMinimum(0);
  histogram->SetMaximum(1);
  can->Update();
  can->Print(_ConfusionMatrixFileName.c_str());
}
