/*
--------------- implement Particle Identification for low momentum muon pion separation  -------------

This code includes class for Particle ID.
 
Boosted Decision Tree with Gradient Boosting is used for separation.
In the training, pure muon and pion samples generated using Particle Gun are used.
In total 19 different low  momenta are investigated starting from 0.2 GeV to 2 GeV.

The samples are generated as following:
- shoot the particles directly to ECAL (keep Z position fixed) and 
- change the position of the gun with a defined step size (changes depending on the number of particles generated) and 
- smear them over the polar and azimuthal angles uniformly

============== G4macro for particle gun:=====================
/generator/generator particleGun
/gun/position 400 400 2450 mm            
/gun/positionStep 0.075 0.075 0.0 mm
/gun/direction 0 0 1
/gun/directionSmearingMode uniform
/gun/phiSmearing 360 deg
/gun/thetaSmearing 90 deg
/gun/particle mu-
/gun/momentum 1.0 GeV
/run/beamOn 20000
/gun/info
=============================================================

Four discriminative varibales are used in the separation;
- depth of the cluster
- mean value of the radius of hits wrt cluster cog 
- RMS value of the radius of hits wrt cluster cog
- cluster energy over track momentum

For any comments or questions: hale.sert@desy.de 
*/

#include <vector>
#include <string>

#include <EVENT/LCCollection.h>

#include <sstream>
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"

#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ClusterImpl.h"
#include <lcio.h>

#include "LowMomentumMuPiSeparationPID_BDTG.hh"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


using namespace std;

LowMomentumMuPiSeparationPID_BDTG::LowMomentumMuPiSeparationPID_BDTG(std::vector<std::string> fname){

    reader = new TMVA::Reader("Silent");

   //Add variables to Reader
   reader->AddVariable( "Dclus", &Dclus);
   reader->AddVariable( "EclOvPtr", &EclOvPtr);
   reader->AddVariable( "Rmean",    &Rmean);
   reader->AddVariable( "Rrms",     &Rrms);

   // Book method(s)  with weight file                                                                           
   std::vector<std::string> Pvalue;
    Pvalue.push_back("02GeVP");
    Pvalue.push_back("03GeVP");
    Pvalue.push_back("04GeVP");
    Pvalue.push_back("05GeVP");
    Pvalue.push_back("06GeVP");
    Pvalue.push_back("07GeVP");
    Pvalue.push_back("08GeVP");
    Pvalue.push_back("09GeVP");
    Pvalue.push_back("10GeVP");
    Pvalue.push_back("11GeVP");
    Pvalue.push_back("12GeVP");
    Pvalue.push_back("13GeVP");
    Pvalue.push_back("14GeVP");
    Pvalue.push_back("15GeVP");
    Pvalue.push_back("16GeVP");
    Pvalue.push_back("17GeVP");
    Pvalue.push_back("18GeVP");
    Pvalue.push_back("19GeVP");
    Pvalue.push_back("20GeVP");


    TString myMethod;
    
    for(int i=0; i < 19; i++){
        myMethod = "BDTG_"+ Pvalue[i]+ "_clusterinfo";   
        weightfile = fname[i];
        reader->BookMVA( myMethod, weightfile );     
      }
}

Int_t LowMomentumMuPiSeparationPID_BDTG::MuPiSeparation(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec){
   
    double tmpid=-1;
    //  std::cout << "==> start ClassApplication" << std::endl;
    
    Dclus=0;
    EclOvPtr=0;
    Rmean=0;
    Rrms=0;
    
    // BDTG cut values obtained after training
   Float_t cut02 = -0.0093;  // S/sqrt(S+B) = 87.8137
   Float_t cut03 = -0.1344;  // S/sqrt(S+B) = 86.6146
   Float_t cut04 = -0.0118;  // S/sqrt(S+B) = 87.5220
   Float_t cut05 = -0.0965 ; // S/sqrt(S+B) = 87.4062
   Float_t cut06 = -0.2278;  // S/sqrt(S+B) = 89.0582
   Float_t cut07 =  0.0030;  // S/sqrt(S+B) = 90.5962
   Float_t cut08 = -0.1708;  // S/sqrt(S+B) = 91.3660
   Float_t cut09 =  0.0060;  // S/sqrt(S+B) = 92.1172
   Float_t cut10 = -0.1957;  // S/sqrt(S+B) = 93.3845
   Float_t cut11 = -0.3255;  // S/sqrt(S+B) = 93.5227
   Float_t cut12 = -0.1523;  // S/sqrt(S+B) = 94.2114
   Float_t cut13 = -0.1290;  // S/sqrt(S+B) = 94.0152
   Float_t cut14 = -0.1301;  // S/sqrt(S+B) = 94.5915
   Float_t cut15 = -0.1729;  // S/sqrt(S+B) = 95.0980
   Float_t cut16 = -0.1704;  // S/sqrt(S+B) = 95.6708
   Float_t cut17 = -0.1978;  // S/sqrt(S+B) = 95.7158
   Float_t cut18 = -0.0110;  // S/sqrt(S+B) = 95.7112
   Float_t cut19 = -0.1871;  // S/sqrt(S+B) = 95.8019
   Float_t cut20 = -0.2328;  // S/sqrt(S+B) = 95.7793    
   
   
   // get variables used low momentum muon and pion separation --
   // cluster shape parameters
    if(cluvec.size()!=0){      
        shapes=cluvec[0]->getShape();  
    }

   float clene=0;
   float clpox=0;
   float clpoy=0;
   float clpoz=0;
   // float chene[10000]={0};
   float chpox[10000]={0};
   float chpoy[10000]={0};
   // float chpoz[10000]={0};

   // int nhit=0;

  float Rhits_clu = 0;
  float Rhits2_clu = 0;
  float rsum_2cluster =0;
  float r2sum_2cluster =0;
  float isum_2cluster =0;
  int m_clu=0;
  float cosalphacluster=0;
  unsigned nch=0;
 
  // const EVENT::TrackStateVec & tss = trk->getTrackStates() ;
  
  const lcio::TrackState* ts ;
  
  ts = trk->getTrackState( lcio::TrackState::AtCalorimeter )  ;
  
  float tscpx = ts->getReferencePoint()[0] ;
  float tscpy = ts->getReferencePoint()[1] ;
  float tscpz = ts->getReferencePoint()[2] ;
  
  for(size_t iclu=0; iclu < cluvec.size(); iclu++){
      
      float Eclu=0;
      if(cluvec[iclu]->getEnergy() > Eclu){
          Eclu=cluvec[iclu]->getEnergy();
          clene =  cluvec[iclu]->getEnergy();
          clpox =  cluvec[iclu]->getPosition()[0];
          clpoy =  cluvec[iclu]->getPosition()[1];
          clpoz =  cluvec[iclu]->getPosition()[2];
          
          EVENT::CalorimeterHitVec tempvec= cluvec[iclu]->getCalorimeterHits();
          
          nch =tempvec.size();
          for( unsigned jhit=0; jhit < nch ; ++jhit ) {
              // nhit = tempvec.size();
              // chene[jhit] =  tempvec[jhit]->getEnergy();
              chpox[jhit] =  tempvec[jhit]->getPosition()[0];
              chpoy[jhit] =  tempvec[jhit]->getPosition()[1];
              // chpoz[jhit] =  tempvec[jhit]->getPosition()[2];
          }
      }
  }
  for(size_t iclu=0; iclu < cluvec.size(); iclu++){
      for( unsigned jhit=0; jhit < nch ; ++jhit ) {
          
          Rhits_clu=sqrt(pow(chpox[jhit]-clpox,2)+pow(chpoy[jhit]-clpoy,2));
          Rhits2_clu=pow(chpox[jhit]-clpox,2)+pow(chpoy[jhit]-clpoy,2);
          
          rsum_2cluster +=Rhits_clu;
          r2sum_2cluster += Rhits2_clu;
          isum_2cluster = m_clu++;
      }
      
      // if cluster cog is at endcap
      _isValid=false;
      if(abs(clpoz) > 2450){
	_isValid=true;
          cosalphacluster = (clpoz - tscpz)/(sqrt(pow(clpox -tscpx,2)+pow(clpoy - tscpy,2)+pow(clpoz- tscpz,2)));
          Dclus = (clpoz - tscpz)/cosalphacluster;
          EclOvPtr = clene/pp.P(); //truemom; 
          Rmean = rsum_2cluster/isum_2cluster;
          Rrms = sqrt(r2sum_2cluster/isum_2cluster);
      }
      else {
          Dclus=0;
          EclOvPtr=0;
          Rmean=0;
          Rrms=0;
      } 
  }
  
  if(shapes.size()!=0){
      if(Dclus!=0 && EclOvPtr!=0 && Rmean !=0 &&  Rrms!=0){ 
          if(0.15< pp.P() && pp.P()<= 0.25){
              mvaout = reader->EvaluateMVA("BDTG_02GeVP_clusterinfo");
              if(mvaout > cut02) tmpid=1; 
              else tmpid=2;
          }
          else if(0.25< pp.P() && pp.P()<=0.35){
              mvaout = reader->EvaluateMVA("BDTG_03GeVP_clusterinfo");
              if(mvaout > cut03) tmpid=1;
              else tmpid=2;
          }
          else if(0.35< pp.P() && pp.P()<=0.45){
              mvaout = reader->EvaluateMVA("BDTG_04GeVP_clusterinfo");
              if(mvaout > cut04) tmpid=1;
              else tmpid=2;
          }
          else if(0.45< pp.P() && pp.P()<=0.55){
              mvaout = reader->EvaluateMVA("BDTG_05GeVP_clusterinfo");
              if(mvaout > cut05) tmpid=1;
              else tmpid=2;
          }
          else if(0.55< pp.P() && pp.P()<=0.65){
              mvaout = reader->EvaluateMVA("BDTG_06GeVP_clusterinfo");
              if(mvaout > cut06) tmpid=1;
              else tmpid=2;
          }
          else if(0.65< pp.P() && pp.P()<=0.75){
              mvaout = reader->EvaluateMVA("BDTG_07GeVP_clusterinfo");
              if(mvaout > cut07) tmpid=1;
              else tmpid=2;
          }
          else if(0.75< pp.P() && pp.P()<=0.85){
              mvaout = reader->EvaluateMVA("BDTG_08GeVP_clusterinfo");
              if(mvaout > cut08) tmpid=1;
              else tmpid=2;
          }
          else if(0.85< pp.P() && pp.P()<=0.95){
              mvaout = reader->EvaluateMVA("BDTG_09GeVP_clusterinfo");
              if(mvaout > cut09) tmpid=1;
              else tmpid=2;
          }
          else if(0.95< pp.P() && pp.P()<=1.05){
              mvaout = reader->EvaluateMVA("BDTG_10GeVP_clusterinfo");
              if(mvaout > cut10) tmpid=1;
              else tmpid=2;
          }
          else if(1.05< pp.P() && pp.P()<=1.15){
              mvaout = reader->EvaluateMVA("BDTG_11GeVP_clusterinfo");
              if(mvaout > cut11) tmpid=1;
              else tmpid=2;
          }
          else if(1.15< pp.P() && pp.P()<=1.25){
              mvaout = reader->EvaluateMVA("BDTG_12GeVP_clusterinfo");
              if(mvaout > cut12) tmpid=1;
              else tmpid=2;
          }
          else if(1.25< pp.P() && pp.P()<=1.35){
              mvaout = reader->EvaluateMVA("BDTG_13GeVP_clusterinfo");
              if(mvaout > cut13) tmpid=1;
              else tmpid=2;
          }
          else if(1.35< pp.P() && pp.P()<=1.45){
              mvaout = reader->EvaluateMVA("BDTG_14GeVP_clusterinfo");
              if(mvaout > cut14) tmpid=1;
              else tmpid=2;
          }
          else if(1.45< pp.P() && pp.P()<=1.55){
              mvaout = reader->EvaluateMVA("BDTG_15GeVP_clusterinfo");
              if(mvaout > cut15) tmpid=1;
              else tmpid=2;
          }
          else if(1.55< pp.P() && pp.P()<=1.65){
              mvaout = reader->EvaluateMVA("BDTG_16GeVP_clusterinfo");
              if(mvaout > cut16) tmpid=1;
              else tmpid=2;
          }
          else if(1.65< pp.P() && pp.P()<=1.75){
              mvaout = reader->EvaluateMVA("BDTG_17GeVP_clusterinfo");
              if(mvaout > cut17) tmpid=1;
              else tmpid=2;
          }
          else if(1.75< pp.P() && pp.P()<=1.85){
              mvaout = reader->EvaluateMVA("BDTG_18GeVP_clusterinfo");
              if(mvaout > cut18) tmpid=1;
              else tmpid=2;
          }
          else if(1.85< pp.P() && pp.P()<=1.95){
              mvaout = reader->EvaluateMVA("BDTG_19GeVP_clusterinfo");
              if(mvaout > cut19) tmpid=1;
              else tmpid=2;
          }
          else if(1.95< pp.P() && pp.P()<=2.05){
              mvaout = reader->EvaluateMVA("BDTG_20GeVP_clusterinfo");
              if(mvaout > cut20) tmpid=1;
              else tmpid=2;
          }
          else mvaout = -100.0;
      }
  } // shape size
  
return tmpid;
   
}// end of loop


Float_t LowMomentumMuPiSeparationPID_BDTG::getMVAOutput(){
    return mvaout;
}

bool LowMomentumMuPiSeparationPID_BDTG::isValid(){
    return _isValid;
}

LowMomentumMuPiSeparationPID_BDTG::~LowMomentumMuPiSeparationPID_BDTG(){
    delete reader;
}
