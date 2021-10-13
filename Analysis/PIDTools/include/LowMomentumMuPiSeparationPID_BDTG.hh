#ifndef LowMomentumMuPiSeparationPID_BDTG_hh
#define LowMomentumMuPiSeparationPID_BDTG_hh 1

#include <string>

#include "TLorentzVector.h"
#include <EVENT/LCCollection.h>

#include "EVENT/Cluster.h"
#include "EVENT/Track.h"
#include "TMVA/Reader.h"

class LowMomentumMuPiSeparationPID_BDTG{
public:

  LowMomentumMuPiSeparationPID_BDTG(const LowMomentumMuPiSeparationPID_BDTG&) = delete;
  LowMomentumMuPiSeparationPID_BDTG& operator=(const LowMomentumMuPiSeparationPID_BDTG&) = delete;

  LowMomentumMuPiSeparationPID_BDTG(std::vector< std::string > fname);
  
  ~LowMomentumMuPiSeparationPID_BDTG();
  
  TMVA::Reader *reader{};
  
  Int_t MuPiSeparation(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  Float_t Dclus{}, EclOvPtr{}, Rmean{}, Rrms{};
  TString weightfile{};
  
  Float_t getMVAOutput();
  bool isValid();
  
  Float_t mvaout{};
  bool _isValid{};
  EVENT::FloatVec shapes{};
    
};

#endif 
