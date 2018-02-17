#ifndef LIKELIHOODPID_hh
#define LIKELIHOODPID_hh 1

#include <string>
#include <sstream>
#include "TLorentzVector.h"
#include "EVENT/Cluster.h"
#include "EVENT/Track.h"
#include "TFile.h"
#include "TH1.h"

class LikelihoodPID{
public:
  LikelihoodPID(){};
  LikelihoodPID(std::string fname, double *pars, std::vector<float> cost);
  LikelihoodPID(double *pars);
  ~LikelihoodPID();
   
  //bool Class_electron(int trkid, jetdata data);
  //bool Class_muon(int trkid, jetdata data);
  int Classification(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  double *GetPosterior();
  double *GetLikelihood();
  double getCorrEnergy(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  double getCorrEnergy(TLorentzVector pp, int parttype);

  //for leptonIDlikelihood
  double get_dEdxDist(int parttype);

  //for each kind of PID flg
  void setBasicFlg(bool flg) {_basicFlg = flg;}
  void setdEdxFlg(bool flg) {_dEdxFlg = flg;}
  void setShowerShapesFlg(bool flg) {_showerShapesFlg = flg;}

  double get_dEdxChi2(int parttype, TVector3 p, float hit,  double dEdx);
  double get_dEdxFactor(int parttype, TVector3 p, float hit,  double dEdx);

  void CalculateDeltaPosition(float charge, TVector3 p, const float* calpos);

private:
  double get_Norm( double dedx, float hit,  double trkcos);
  double BetheBloch( double x,  double mass,  double *pars);
  
  int Class_electron(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  int Class_muon(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  int Class_hadron(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  const  double getValue(int type, int valtype,  double value);
  double getPenalty(int ptype, int hypothesis,  double p);
  
  double par[5][5] = {};
  TFile* fpdf = NULL;
  TH1F* pdf[6][20] ={};
  
  //define masses
  double emass = 0.0;
  double mmass = 0.0;
  double pimass = 0.0;
  double kmass = 0.0;
  double pmass = 0.0;
  
  //threshold
  double threshold[5] = {};
  //penality
  double penalty[5][5] = {};
  double prior[5] = {};
  
  double fact[5][5] = {};
  
  //weights for hadron likelihood calculation
  double _weights[6][20] = {};

  //posterior
  double _posterior[6] = {};   //add ahdron type
  double _likelihood[6] = {};  //add hadron type
  
  //distance from bethe bloch line with each particle hypothesis
  double _dEdxDist[5] = {};
  
  //for shower profile
  EVENT::FloatVec shapes;

  bool _basicFlg{}, _dEdxFlg{}, _showerShapesFlg{};
  int _usebayes{}, _usecorr{};
  float _dEdxnorm{}, _dEdxerrfact{}, _bfield{};
  double _delpos[3]{};
};

#endif 
