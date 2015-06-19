#ifndef LIKELIHOODPID_hh
#define LIKELIHOODPID_hh 1

#include <string>

#include "TLorentzVector.h"

#include "EVENT/Cluster.h"
#include "EVENT/Track.h"

class TFile;
class TH1F;

class LikelihoodPID{
public:
  LikelihoodPID(std::string fname);
  ~LikelihoodPID();
   
  //Bool_t Class_electron(Int_t trkid, jetdata data);
  //Bool_t Class_muon(Int_t trkid, jetdata data);
  Int_t Classification(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  Double_t *GetPosterior();
  Double_t getCorrEnergy(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  Double_t getCorrEnergy(TLorentzVector pp, Int_t parttype);

  //for leptonIDlikelihood
  Double_t get_dEdxChi2(Int_t parttype, TVector3 p, Float_t hit, Double_t dEdx);

private:
  Double_t get_Norm(Double_t dedx, Float_t hit, Double_t trkcos);
  Double_t BetheBloch(Double_t x, Double_t mass, Double_t *pars);

  Int_t Class_electron(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  Int_t Class_muon(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  Int_t Class_hadron(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  Double_t getValue(Int_t type, Int_t valtype, Double_t value);
  Double_t getPenalty(Int_t ptype, Int_t hypothesis, Double_t p);

  Double_t par[5][5];
  TFile* fpdf;
  TH1F* pdf[6][20];

  //define masses
  Double_t emass;
  Double_t mmass;
  Double_t pimass;
  Double_t kmass;
  Double_t pmass;

  //threshold
  Double_t threshold[5];
  //penality
  Double_t penalty[5][5];
  Double_t prior[5];

  Double_t fact[5][5];

  Double_t _posterior[5];

  EVENT::FloatVec shapes;
};

#endif 
