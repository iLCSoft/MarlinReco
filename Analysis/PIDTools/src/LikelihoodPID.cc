/*
--------------- inplement Particle Identofication using Global likelihood -------------

This code includes class for Particle ID. The method is Naive Bayes Classfier
1. First, check whether a particle is electron or not. 
   If so, the particle is identified as electron and return electron status code.
2. And then, check whether a particle a particle is muonn or not. 
   If so, the particle is identified as muon and return muon status code.
3. If a particle is neither electron nor muon, 
   try to classify a particle to each Hadron type(Pion, Kaon, or Proton)

TODO:
Bethe-Bloch parameters should be moved to steer file
Threshold parameters should be moved to steer file

make flags to choose Particle ID method 

risk-minimization and MAP
 */

#include <string>
#include <sstream>
#include <TFile.h>
#include <TH1F.h>
#include "LikelihoodPID.hh"

using namespace std;

std::string itos(int i)  
{
  std::stringstream s;
  s << i;
  return s.str();
}

LikelihoodPID::LikelihoodPID(string fname){
  //set parameters for bethebloch
  //electron                                                                                                                   
  par[0][0]=-2.40638e-03;   //3.36491e-04;                                            
  par[0][1]=7.10337e-01;    //1.20930e-11;                                            
  par[0][2]=2.87718e-01;    //1.28546e-01;                                                    
  par[0][3]=-1.71591e+00;    //-1.65870e+00;                                               
  par[0][4]=0.0;  //do not fit this parameter to avoid abnormal termination                          

  //muon                                                                                                            
  par[1][0]=8.11408e-02;
  par[1][1]=9.92207e-01;
  par[1][2]=7.58509e+05;
  par[1][3]=-1.70167e-01;
  par[1][4]=4.63670e-04;

  //pion                                                                                            
  par[2][0]=8.10756e-02;
  par[2][1]=-1.45051e+06;
  par[2][2]=-3.09843e+04;
  par[2][3]=2.84056e-01;
  par[2][4]=3.38131e-04;

  //kaon                                                                             
  par[3][0]=7.96117e-02;
  par[3][1]=4.13335e+03;
  par[3][2]=1.13577e+06;
  par[3][3]=1.80555e-01;
  par[3][4]=-3.15083e-04;

  //proton                                                                            
  par[4][0]=7.78772e-02;
  par[4][1]=4.49300e+04;
  par[4][2]=9.13778e+04;
  par[4][3]=1.50088e-01;
  par[4][4]=-6.64184e-04;

  //set mass
  emass=0.000510998;
  mmass=0.105658;
  pimass=0.139570;
  kmass=0.493677;
  pmass=0.938272;

  //get pdf plots
  string ffstr1 = fname;        //"pdf/pdf_ParticleID_ok.root";
  fpdf=new TFile(ffstr1.c_str());

  string hname,hname2;
  for(Int_t i=0;i<6;i++){
    //ep
    hname="hep" + itos(i+1);
    hname2="hep" + itos(i+1) + "_2";
    pdf[i][0]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //ehad
    hname="hehad" + itos(i+1);
    hname2="hehad" + itos(i+1) + "_2";
    pdf[i][1]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //mucal
    hname="hmucal" + itos(i+1);
    hname2="hmucal" + itos(i+1) + "_2";
    pdf[i][2]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //chi2
    hname="hchi2" + itos(i+1);
    hname2="hchi2" + itos(i+1) + "_2";
    pdf[i][3]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //showermax/exp.showermax
    hname="hldiscrepancy" + itos(i+1);
    hname2="hldiscrepancy" + itos(i+1) + "_2";
    pdf[i][4]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //absorption length
    hname="htdiscrepancy" + itos(i+1);
    hname2="htdiscrepancy" + itos(i+1) + "_2";
    pdf[i][5]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //xl20
    hname="hxl20" + itos(i+1);
    hname2="hxl20" + itos(i+1) + "_2";
    pdf[i][6]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //likeliele
    hname="hlikeliele" + itos(i+1);
    hname2="hlikeliele" + itos(i+1) + "_2";
    pdf[i][7]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //likelimuo
    hname="hlikelimuo" + itos(i+1);
    hname2="hlikelimuo" + itos(i+1) + "_2";
    pdf[i][8]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //likelipi
    hname="hlikelipi" + itos(i+1);
    hname2="hlikelipi" + itos(i+1) + "_2";
    pdf[i][9]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //likelik
    hname="hlikelik" + itos(i+1);
    hname2="hlikelik" + itos(i+1) + "_2";
    pdf[i][10]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //likelip
    hname="hlikelip" + itos(i+1);
    hname2="hlikelip" + itos(i+1) + "_2";
    pdf[i][11]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //deltax
    hname="hdeltax" + itos(i+1);
    hname2="hdeltax" + itos(i+1) + "_2";
    pdf[i][12]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    //deltaz
    hname="hdeltaz" + itos(i+1);
    hname2="hdeltaz" + itos(i+1) + "_2";
    pdf[i][13]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
  }
  
  //normalize histograms
  Double_t weight=1.0;
  for(Int_t i=0;i<6;i++){
    for(Int_t j=0;j<14;j++){
      //normalize histograms
      weight=pdf[i][j]->Integral(0,pdf[i][j]->GetNbinsX()+1,"");
      pdf[i][j]->Scale(1.0/weight);
    }
  }

  //set threshold
  //this is original
  threshold[0]=TMath::Exp(-0.55);
  threshold[1]=TMath::Exp(-0.7);
  threshold[2]=0.0; 
  threshold[3]=TMath::Exp(-0.6); 
  threshold[4]=TMath::Exp(-0.6);
  
  //set factor
  fact[0][0]=1.0;
  fact[0][1]=1.0;
  fact[0][2]=1.0;
  fact[0][3]=1.0;
  fact[0][4]=1.0;
  fact[1][0]=1.0;
  fact[1][1]=1.0;
  fact[1][2]=1.0;
  fact[1][3]=1.0;
  fact[1][4]=1.0;
  fact[2][0]=1.0;
  fact[2][1]=1.0;
  fact[2][2]=1.0/1.02239;
  fact[2][3]=1.0/1.00616;
  fact[2][4]=1.0/0.965142;
  fact[3][0]=1.0;
  fact[3][1]=1.0;
  fact[3][2]=1.0/1.01857;
  fact[3][3]=1.0/1.00038;
  fact[3][4]=1.0/0.986345;
  fact[4][0]=1.0;
  fact[4][1]=1.0;
  fact[4][2]=1.0/1.09023;
  fact[4][3]=1.0/1.00211;
  fact[4][4]=1.0/1.00129;
  return;
}

LikelihoodPID::~LikelihoodPID(){

  fpdf->Close();

  return;
}

//public
Int_t LikelihoodPID::Classification(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec){
  Int_t tmpid=-1;

  //set prior first
  for(Int_t i=0;i<5;i++){
    prior[i]=0.2;
    _posterior[i]=prior[i];
  }

  //get dedx distance
  for(Int_t i=0;i<5;i++) _dEdxDist[i]=get_dEdxChi2(i, pp.Vect(), trk->getdEdxError(),trk->getdEdx());

  //first check whether electron or not
  for(Int_t i=0;i<5;i++){
    _likelihood[i]=0.0;
    if(_basicFlg==false && _dEdxFlg==true && _showerShapesFlg==false) prior[i]=0.2;
  }
  tmpid=Class_electron(pp, trk, cluvec);
  if(tmpid==0) return tmpid;

  //second check whether muon or not
  for(Int_t i=0;i<5;i++){
    _likelihood[i]=0.0;
    if(_basicFlg==false && _dEdxFlg==true && _showerShapesFlg==false) prior[i]=0.2;
  }
  tmpid=Class_muon(pp, trk, cluvec);
  if(tmpid==1) return tmpid;

  //third classify hadrons
  for(Int_t i=0;i<5;i++){
    _likelihood[i]=0.0;
    if(_basicFlg==false && _dEdxFlg==true && _showerShapesFlg==false) prior[i]=0.2;
  }
  tmpid=Class_hadron(pp, trk, cluvec);
  return tmpid;
}

Double_t *LikelihoodPID::GetPosterior(){

  return _posterior;
}

Double_t *LikelihoodPID::GetLikelihood(){

  return _likelihood;
}

//track energy correction
Double_t LikelihoodPID::getCorrEnergy(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec){
  Int_t parttype=Classification(pp, trk, cluvec);
  if(parttype==-1) parttype=2;  //move to pion class

  return getCorrEnergy(pp, parttype);
}

Double_t LikelihoodPID::getCorrEnergy(TLorentzVector pp, Int_t parttype){
  if(parttype==-1) parttype=2;  //move to pion class

  Double_t tmpmass=0.0;
  if(parttype==0) tmpmass=emass;
  if(parttype==1) tmpmass=mmass;
  if(parttype==2) tmpmass=pimass;
  if(parttype==3) tmpmass=kmass;
  if(parttype==4) tmpmass=pmass;

  return sqrt(pp.P()*pp.P()+tmpmass*tmpmass);
}

//private
Double_t LikelihoodPID::get_dEdxChi2(Int_t parttype, TVector3 p, Float_t hit, Double_t dEdx){
  //get parameters for chi2
  Double_t tmppar[5],tmpmass=0.0;
  for(Int_t i=0;i<5;i++) tmppar[i]=par[parttype][i];
  //getmass
  switch(parttype){
  case 0:
    tmpmass=emass;
    break;
  case 1:
    tmpmass=mmass;
    break;
  case 2:
    tmpmass=pimass;
    break;
  case 3:
    tmpmass=kmass;
    break;
  case 4:
    tmpmass=pmass;
    break;
  }

  //cal. polar angle
  Double_t trkcos=p.CosTheta();

  //get nomalized dEdx
  Double_t dEdx_Norm=get_Norm(dEdx, hit, trkcos);

  //get expected dE/dx
  Double_t ExpdEdx=BetheBloch(p.Mag(),tmpmass,tmppar);

  //cout << "check " << emass << " " << tmpmass << " " << trkcos << " " << dEdx_Norm << " " << ExpdEdx << endl;

  //get chi2!!(so far 5% error imposed. conservative)
  Double_t chi2=TMath::Power((dEdx_Norm-ExpdEdx)/(0.05*dEdx_Norm),2.0);
  if(dEdx_Norm-ExpdEdx<0.0) chi2=-chi2;    //get signed chi2
  return chi2;
}

//public
Double_t LikelihoodPID::get_dEdxDist(Int_t parttype){
  //get parameters for dedxdist
  Double_t dedxdist=0.0;

  if(_dEdxDist[parttype]!=0.0) dedxdist = sqrt(fabs(_dEdxDist[parttype]))*_dEdxDist[parttype]/fabs(_dEdxDist[parttype]);
  if(dedxdist!=dedxdist) dedxdist=0.0;
  return dedxdist;
}

Double_t LikelihoodPID::get_Norm(Double_t dedx, Float_t hit, Double_t trkcos){
  //cal. hit dep.
  Double_t f1=1.0;   //1.0+TMath::Exp(-hit/1.468);   //already corrected
  //cal. polar angle dep.
  //Double_t c=1.0/sqrt(1.0-trkcos*trkcos);
  Double_t f2=1.0;   //1.0/(1.0-0.08887*TMath::Log(c));    //already corrected
  
  return dedx*f1*f2/(1.350e-7);
}

Double_t LikelihoodPID::BetheBloch(Double_t x, Double_t mass, Double_t *pars){
  Double_t bg=x/mass;
  Double_t b=sqrt(bg*bg/(1.0+bg*bg));
  //Double_t g=bg/b;
  Double_t tmax=pars[2]*TMath::Power(bg,2.0);   ///(1.0+pars[3]*g+pars[4]);

  return (0.5*pars[0]*TMath::Log(pars[1]*TMath::Power(bg,2.0)*tmax)-pars[3]*b*b-pars[4]*bg/2.0)/(b*b);
}
 
Int_t LikelihoodPID::Class_electron(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec){
  //Bayesian approach to classify electron!
  
  //get deposit energy
  Double_t ecal=0.0,hcal=0.0,mucal=0.0;
  if(cluvec.size()!=0){
    for(unsigned int i=0;i<cluvec.size();i++){
      ecal+=cluvec[i]->getSubdetectorEnergies()[0];
      hcal+=cluvec[i]->getSubdetectorEnergies()[1];
      mucal+=cluvec[i]->getSubdetectorEnergies()[2];
    }
  }

  //get track shape parameters
  if(cluvec.size()!=0) shapes=cluvec[0]->getShape();

  //get variables
  Double_t var[11]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //11 variales used(is it OK?)
  //for track variables
  var[0]=(ecal+hcal)/pp.P();
  if(ecal+hcal!=0.0) var[1]=ecal/(ecal+hcal);
  if(shapes.size()!=0){
    var[2]=shapes[0];
    var[3]=shapes[5];
    var[4]=fabs(shapes[3])/(shapes[6]);
    var[5]=shapes[15]/(2.0*3.50);

  }else var[2]=-1.0;

  var[6]=get_dEdxChi2(0,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[7]=get_dEdxChi2(1,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[8]=get_dEdxChi2(2,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[9]=get_dEdxChi2(3,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[10]=get_dEdxChi2(4,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[6]=-0.5*fabs(var[6]);
  var[7]=-0.5*fabs(var[7]);
  var[8]=-0.5*fabs(var[8]);
  var[9]=-0.5*fabs(var[9]);
  var[10]=-0.5*fabs(var[10]);

  //get likelihood for each class
  Double_t posterior[5]={0.0,0.0,0.0,0.0,0.0};
  //Double_t risk[5]={0.0,0.0,0.0,0.0,0.0};
  Double_t okval[5]={0.0,0.0,0.0,0.0,0.0};
  Double_t total=0.0;
  Double_t priorprob[5]={0.0,0.0,0.0,0.0,0.0};
  Int_t valtype=0;
  for(Int_t j=0;j<11;j++){   //variables
    //avoid very low momentum tracks (no energy deposit in the cal.)
    if(var[0]==0.0 && j<=5) continue;

    //cout << "value: " << j << " " << var[j] << endl;
    //first, get likelihood
    for(Int_t i=0;i<5;i++){   //particle type
      if(j==0) valtype=0;
      if(j==1) valtype=1;
      if(j==2) valtype=3;
      if(j==3) valtype=4;
      if(j==4) valtype=5;
      if(j==5) valtype=6;
      if(j==6) valtype=7;
      if(j==7) valtype=8;
      if(j==8) valtype=9;
      if(j==9) valtype=10;
      if(j==10) valtype=11;
      
      okval[i]=getValue(i,valtype,var[j]);   //likelihood
      if(mucal==0.0 && i==1) okval[i]=getValue(5,valtype,var[j]);   //likelihood
      //cout << "var: " << i << " " << j << " " << var[j] << " " << okval[i] << endl;
    }
    
    //for basic variables flg
    if(!_basicFlg && j<2) continue;
    //for cluster shape flg
    if(!_showerShapesFlg && j>=2 && j<=5) continue;
    //for dEdx flg
    if(!_dEdxFlg && j>=6) continue; 
   
    //to avoid strange value for dE/dx
    if(j>=6 && var[j]<=-50.0) continue;
    
    //cal. prior probability
    for(Int_t i=0;i<5;i++)   //particle type
      priorprob[i]=prior[i]/
	(prior[0]+prior[1]+prior[2]+prior[3]+prior[4]);
        
    //cal total probability for each particle type
    total=0.0;
    for(Int_t i=0;i<5;i++) total+=priorprob[i]*okval[i];
    
    //cal. log posterior probability!
    for(Int_t i=0;i<5;i++){   //particle type
      ///20150708 change it from posterior probability to likelihood
      posterior[i]=TMath::Log(okval[i])+TMath::Log(priorprob[i])-TMath::Log(total);
      _likelihood[i]+=TMath::Log(okval[i]);

    }

    //if(_basicFlg==false && _dEdxFlg==true && _showerShapesFlg==false)
    // cout << "check posterior: " << TMath::Exp(posterior[0]) << " " << TMath::Exp(posterior[1]) << " " <<  TMath::Exp(posterior[2]) << endl;
    
    //bayesian updating
    for(Int_t i=0;i<5;i++) prior[i]=TMath::Exp(posterior[i]);
  }
  
  //cal. risk
  //first. get penalty
  /*for(Int_t i=0;i<5;i++){   //particle type
    for(Int_t j=0;j<5;j++){
      penalty[i][j]=1.0;  
      if(i==j) penalty[i][j]=1.0e-50;
    }
  }

  for(Int_t i=0;i<5;i++){   //particle type
    risk[i]=penalty[i][0]*TMath::Exp(posterior[0])+penalty[i][1]*TMath::Exp(posterior[1])
      +penalty[i][2]*TMath::Exp(posterior[2])+penalty[i][3]*TMath::Exp(posterior[3])+penalty[i][4]*TMath::Exp(posterior[4]);
      }*/

  //check the rule
  Int_t okflg=-1;
  Double_t tmppp=-1.0e+100;
  //---------- here is likelihood based determination rule -----------
  for(Int_t i=0;i<5;i++){
    tmppp=TMath::Max(posterior[i], tmppp);
    //if(posterior[0]<posterior[i]){
    //  okflg=false;
    //  break;
    //}
  }
  
  for(Int_t i=0;i<5;i++){
    if(fabs(tmppp-posterior[i])<1.0e-6) okflg=i;
  }
  //------------------------------------------------------------------

  //---------- here is risk based determination rule ----------- 
  /*tmppp=1.0e+100;
  for(Int_t i=0;i<5;i++){
    tmppp=TMath::Min(risk[i], tmppp);
  }
  
  for(Int_t i=0;i<5;i++){
    if(fabs(tmppp-risk[i])<1.0e-8) okflg=i;
    }*/
  //------------------------------------------------------------------
  
  //threshold is set.
  //if(okflg==0 && posterior[okflg]<TMath::Log(threshold[okflg])) okflg=-1;

  //save posterior
  _posterior[0]=TMath::Exp(posterior[0]);
  _posterior[1]=TMath::Exp(posterior[1]);
  _posterior[2]=TMath::Exp(posterior[2]);
  _posterior[3]=TMath::Exp(posterior[3]);
  _posterior[4]=TMath::Exp(posterior[4]);

  return okflg;
}

Int_t LikelihoodPID::Class_muon(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec){
  //Bayesian approach to classify muon!
  
  //get deposit energy
  Double_t ecal=0.0,hcal=0.0,mucal=0.0;
  if(cluvec.size()!=0){
    for(unsigned int i=0;i<cluvec.size();i++){
      ecal+=cluvec[i]->getSubdetectorEnergies()[0];
      hcal+=cluvec[i]->getSubdetectorEnergies()[1];
      mucal+=cluvec[i]->getSubdetectorEnergies()[2];
    }
  }

  //get track shape parameters
  if(cluvec.size()!=0) shapes=cluvec[0]->getShape();

  //get variables
  Double_t var[12]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //12 variales used(is it OK?)
  var[0]=(ecal+hcal)/pp.P();
  if(ecal+hcal!=0.0) var[1]=ecal/(ecal+hcal);
  var[2]=mucal;
  if(shapes.size()!=0){
    var[3]=shapes[0];
    var[4]=shapes[5];
    var[5]=fabs(shapes[3])/(shapes[6]);
    var[6]=shapes[15]/(2.0*3.50);
  }else var[3]=-1.0;
  
  var[7]=get_dEdxChi2(0,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[8]=get_dEdxChi2(1,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[9]=get_dEdxChi2(2,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[10]=get_dEdxChi2(3,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[11]=get_dEdxChi2(4,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[7]=-0.5*fabs(var[7]);
  var[8]=-0.5*fabs(var[8]);
  var[9]=-0.5*fabs(var[9]);
  var[10]=-0.5*fabs(var[10]);
  var[11]=-0.5*fabs(var[11]);

  //get likelihood for each class
  Double_t posterior[5]={0.0,0.0,0.0,0.0,0.0};
  //Double_t risk[5]={0.0,0.0,0.0,0.0,0.0};
  Int_t valtype=0;
  Double_t okval[5]={0.0,0.0,0.0,0.0,0.0};
  Double_t total=0.0;     //[5]={0.0,0.0,0.0,0.0,0.0};
  Double_t priorprob[5]={0.0,0.0,0.0,0.0,0.0};
  for(Int_t j=0;j<12;j++){   //variables

    if(var[0]>0.0){
      //avoid pion misID when energy deposit to mucal is zero
      if(var[2]==0.0 && (j==2 || j==8)) continue;
      //don't use when ep>0.0
      else if(var[2]>0.0 && j>=8) continue;
    }else{
      //avoid very low momentum tracks (no energy deposit in the cal.)
      if(j<=6 || j>=8) continue;
    }

    //first, get likelihood
    for(Int_t i=0;i<5;i++){   //particle type(electron isn't be checked)
      if(j==0) valtype=0;
      if(j==1) valtype=1;
      if(j==2) valtype=2;
      if(j==3) valtype=3;
      if(j==4) valtype=4;
      if(j==5) valtype=5;
      if(j==6) valtype=6;
      if(j==7) valtype=7;
      if(j==8) valtype=8;
      if(j==9) valtype=9;
      if(j==10) valtype=10;
      if(j==11) valtype=11;
      
      okval[i]=getValue(i,valtype,var[j]);   //likelihood
      if(var[2]==0.0 && i==1) okval[i]=getValue(5,valtype,var[j]);   //likelihood 
    }

    //for basic variables flg
    if(!_basicFlg && j<=2) continue;
    //for cluster shape flg
    if(!_showerShapesFlg && j>2 && j<=6) continue;
    //for dEdx flg
    if(!_dEdxFlg && j>6) continue; 
    
    //to avoid strange value for dE/dx
    if(j>=7 && j<=11 && var[j]<=-50.0) continue;
    
    //cal. prior probability
    for(Int_t i=0;i<5;i++)   //particle type
      priorprob[i]=prior[i]/
	(prior[0]+prior[1]+prior[2]+prior[3]+prior[4]);

    //cal total probability for each particle type
    total=0.0;
    for(Int_t i=0;i<5;i++) total+=priorprob[i]*okval[i];
    
    //cal. log likelihood first
    for(Int_t i=0;i<5;i++){   //particle type
      //20150708 change it from posterior probability to simple likelihood
      posterior[i]=TMath::Log(okval[i])+TMath::Log(priorprob[i])-TMath::Log(total); 
      _likelihood[i]+=TMath::Log(okval[i]);
    }
    
    //bayesian updating
    for(Int_t i=0;i<5;i++) prior[i]=TMath::Exp(posterior[i]);
  }

  //cal. risk
  //first. get penalty
  /*for(Int_t i=0;i<5;i++){   //particle type
    for(Int_t j=0;j<5;j++){
      penalty[i][j]=1.0; 
      if(i==j) penalty[i][j]=1.0e-50; 
    }
  }

  for(Int_t i=0;i<5;i++){   //particle type
    risk[i]=penalty[i][0]*TMath::Exp(posterior[0])+penalty[i][1]*TMath::Exp(posterior[1])
      +penalty[i][2]*TMath::Exp(posterior[2])+penalty[i][3]*TMath::Exp(posterior[3])+penalty[i][4]*TMath::Exp(posterior[4]);
      }*/

  //check the rule
  Int_t okflg=-1;
  Double_t tmppp=-1.0e+100;
  //---------- here is likelihood based determination rule -----------
  for(Int_t i=0;i<5;i++){
    tmppp=TMath::Max(posterior[i], tmppp);
    //if(posterior[0]<posterior[i]){
    //  okflg=false;
    //  break;
    //}
  }
  
  for(Int_t i=0;i<5;i++){
    if(fabs(tmppp-posterior[i])<1.0e-6) okflg=i;
  }
  //------------------------------------------------------------------

  //---------- here is risk based determination rule ----------- 
  /*tmppp=1.0e+100;
  for(Int_t i=0;i<5;i++){
    tmppp=TMath::Min(risk[i], tmppp);
  }
  
  for(Int_t i=0;i<5;i++){
    if(fabs(tmppp-risk[i])<1.0e-8) okflg=i;
    }*/
  //------------------------------------------------------------------
  
  //if(okflg==1 && posterior[okflg]<TMath::Log(threshold[okflg])) okflg=-1;   //threshold is set.

  //save posterior
  _posterior[0]=TMath::Exp(posterior[0]);
  _posterior[1]=TMath::Exp(posterior[1]);
  _posterior[2]=TMath::Exp(posterior[2]);
  _posterior[3]=TMath::Exp(posterior[3]);
  _posterior[4]=TMath::Exp(posterior[4]);

  return okflg;  
}

Int_t LikelihoodPID::Class_hadron(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec){
  //Bayesian approach to classify hadrons!

  //get deposit energy
  Double_t ecal=0.0,hcal=0.0,mucal=0.0;
  if(cluvec.size()!=0){
    for(unsigned int i=0;i<cluvec.size();i++){
      ecal+=cluvec[i]->getSubdetectorEnergies()[0];
      hcal+=cluvec[i]->getSubdetectorEnergies()[1];
      mucal+=cluvec[i]->getSubdetectorEnergies()[2];
    }
  }

  //get track shape parameters
  if(cluvec.size()!=0) shapes=cluvec[0]->getShape();

  //get variables
  Double_t var[10]; //10 variales used(is it OK?)
  var[0]=(ecal+hcal)/pp.P();
  if(ecal+hcal!=0.0)var[1]=ecal/(ecal+hcal);
  if(shapes.size()!=0){
    var[2]=shapes[0];
    var[3]=fabs(shapes[3])/(shapes[6]);
    var[4]=shapes[15]/(2.0*3.50);
  }else var[2]=-1.0;
 
  var[5]=get_dEdxChi2(0,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[6]=get_dEdxChi2(1,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[7]=get_dEdxChi2(2,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[8]=get_dEdxChi2(3,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[9]=get_dEdxChi2(4,pp.Vect(),trk->getdEdxError(),trk->getdEdx());
  var[5]=-0.5*fabs(var[5]);
  var[6]=-0.5*fabs(var[6]);
  var[7]=-0.5*fabs(var[7]);
  var[8]=-0.5*fabs(var[8]);
  var[9]=-0.5*fabs(var[9]);

  //get likelihood for each class
  Double_t posterior[5]={TMath::Log(0.20),TMath::Log(0.20),TMath::Log(0.20),TMath::Log(0.20),TMath::Log(0.20)};
  //Double_t risk[5]={0.0,0.0,0.0,0.0,0.0};
  Int_t valtype=0;
  Double_t okval[5]={0.0,0.0,0.0,0.0,0.0};
  Double_t total[5]={0.0,0.0,0.0,0.0,0.0};
  Double_t priorprob[5]={0.0,0.0,0.0,0.0,0.0};
  for(Int_t j=0;j<10;j++){   //variables
    //avoid very low momentum tracks (no energy deposit in the cal.)
    if(var[0]==0.0 && j<=4) continue;

    //first, get likelihood
    for(Int_t i=0;i<5;i++){   //particle type(electron&muon aren't be checked)
      if(j==0) valtype=0;
      if(j==1) valtype=1;
      if(j==2) valtype=3;
      if(j==3) valtype=5;
      if(j==4) valtype=6;
      if(j==5) valtype=7;
      if(j==6) valtype=8;
      if(j==7) valtype=9;
      if(j==8) valtype=10;
      if(j==9) valtype=11;
      
      okval[i]=getValue(i,valtype,var[j]);   //likelihood
    }

    //for basic variables flg
    if(!_basicFlg && j<2) continue;
    //for cluster shape flg
    if(!_showerShapesFlg && j>=2 && j<=4) continue;
    //for dEdx flg
    if(!_dEdxFlg && j>=5) continue; 
    
    //to avoid strange value for dE/dx
    if(j>=5 && var[j]<=-50.0) continue;
    
    //cal. prior probability
    for(Int_t i=0;i<5;i++)   //particle type
      priorprob[i]=prior[i]/
	(prior[0]+prior[1]+prior[2]+prior[3]+prior[4]);
        
    //cal total probability
    for(Int_t i=0;i<5;i++){   //particle type
      total[i]=priorprob[0]*okval[0]+priorprob[1]*okval[1]+priorprob[2]*okval[2]+priorprob[3]*okval[3]+priorprob[4]*okval[4];
      if(total[i]==0.0) total[i]=1.0e-100;
    }

    //cal. log posterior probability!
    for(Int_t i=0;i<5;i++){   //particle type
      //20150708 change it from posterior probabity to simple likelihood
      posterior[i]=TMath::Log(okval[i])+TMath::Log(priorprob[i])-TMath::Log(total[i]);
      _likelihood[i]+=TMath::Log(okval[i]);
    }
    
    //bayesian updating
    for(Int_t i=0;i<5;i++) prior[i]=TMath::Exp(posterior[i]);
  }
    
  //cal. risk
  //first. get penalty
  /*for(Int_t i=2;i<5;i++){   //particle type
    for(Int_t j=2;j<5;j++){
      penalty[i][j]=getPenalty(i,j,pp.P());
    }
    Double_t tmppenalty=penalty[i][i];
    for(Int_t j=2;j<5;j++){
      penalty[i][j]=(fabs(penalty[i][j]-tmppenalty)+tmppenalty)/tmppenalty;
      penalty[i][j]=penalty[i][j]-1.0;
      
    }
    }

  for(Int_t i=0;i<5;i++){   //particle type
    risk[i]=penalty[i][0]*TMath::Exp(posterior[0])+penalty[i][1]*TMath::Exp(posterior[1])
      +penalty[i][2]*TMath::Exp(posterior[2])+penalty[i][3]*TMath::Exp(posterior[3])+penalty[i][4]*TMath::Exp(posterior[4]);
      }*/

  //check the rule
  Int_t okflg=-1;
  Double_t tmppp=-1.0e+100;
  //---------- here is likelihood based determination rule -----------
  for(Int_t i=0;i<5;i++){
    tmppp=TMath::Max(posterior[i], tmppp);
    //if(posterior[0]<posterior[i]){
    //  okflg=false;
    //  break;
    //}
  }
  
  for(Int_t i=0;i<5;i++){
    if(fabs(tmppp-posterior[i])<1.0e-6) okflg=i;
  }

  //move to pion when very bad
  if(TMath::Exp(tmppp)<0.21) okflg=2;
  //------------------------------------------------------------------

  //---------- here is risk based determination rule ----------- 
  /*tmppp=1.0e+100;
  for(Int_t i=0;i<5;i++){
    tmppp=TMath::Min(risk[i], tmppp);
  }
  
  for(Int_t i=0;i<5;i++){
    if(fabs(tmppp-risk[i])<1.0e-6) okflg=i;
  }
  if(posterior[okflg]<=TMath::Log(0.2)) okflg=TMath::Max(2,okflg);  //do not go to leptons*/
  //------------------------------------------------------------------
  
  //check penalty -1 means undefined
  //if(okflg>=2 && posterior[okflg]<TMath::Log(threshold[okflg])) okflg=-1;   //threshold is set.

  //save posterior
  _posterior[0]=TMath::Exp(posterior[0]);
  _posterior[1]=TMath::Exp(posterior[1]);
  _posterior[2]=TMath::Exp(posterior[2]);
  _posterior[3]=TMath::Exp(posterior[3]);
  _posterior[4]=TMath::Exp(posterior[4]);

  return okflg;  
}

const Double_t LikelihoodPID::getValue(Int_t type, Int_t valtype, Double_t value){

  /*Int_t nbins=pdf[type][valtype]->GetNbinsX();
  Double_t interval=pdf[type][valtype]->GetBinWidth(1);
  Double_t minval=pdf[type][valtype]->GetBinLowEdge(1);

  Double_t val=1.0e-30;
  Int_t bin=0;
  for(Int_t i=0;i<nbins;i++){
    if(valtype!=1 && valtype<7){
      if(value<minval){
	bin=0;
	break;
      }
      if(value>=minval+i*interval && value<minval+(i+1)*interval){
	bin=i+1;
	break;
      }
      if(value>=minval+nbins*interval){
	bin=nbins+1;
	break;
      }
    }else if(valtype>=7){   //for hadem distribution!(because 1.0 is the max value!)
      if(value<=minval){
	bin=0;
	break;
      }
      if(value>minval+i*interval && value<=minval+(i+1)*interval){
	bin=i+1;
	break;
      }
      if(value>minval+nbins*interval){
	bin=nbins+1;
	break;
      }
    }else if(valtype==1){    //ehad very strange graph!
      if(value<=minval){
	bin=1;
	break;
      }
      if(value>minval+i*interval && value<=minval+(i+1)*interval){
	bin=i+1;
	break;
      }
      if(value>minval+nbins*interval){
	bin=nbins+1;
	break;
      }
    }
    }*/

  
  Double_t val=1.0e-30;
  Int_t bin=0;

  bin=pdf[type][valtype]->GetXaxis()->FindBin(value);

  //get probability
  val=pdf[type][valtype]->GetBinContent(bin);

  return val;
}

Double_t LikelihoodPID::getPenalty(Int_t ptype, Int_t hypothesis, Double_t p){
  Double_t par[3]={0.0,0.0,0.0};
  //set parameters
  switch(ptype){
  case 0:  //electron
    if(hypothesis==0){
      par[0]=0.0;
      par[1]=1.0;
      par[2]=1.00178;

    }else if(hypothesis==1){
      par[0]=0.101945;
      par[1]=2.18719;
      par[2]=0.990000;

    }else if(hypothesis==2){
      par[0]=0.104521;
      par[1]=1.98490;
      par[2]=0.990000;

    }else if(hypothesis==3){
      par[0]=0.234299;
      par[1]=0.276835;
      par[2]=0.965973;

    }else{
      par[0]=0.414613;
      par[1]=-0.132530;
      par[2]=0.949679;

    }

    break;
  case 1:   //muon
    if(hypothesis==0){
      par[0]=-0.00512190;
      par[1]=-0.211835;
      par[2]=1.00024;

    }else if(hypothesis==1){
      par[0]=0.0;
      par[1]=1.00;
      par[2]=0.999684;

    }else if(hypothesis==2){
      par[0]=0.00833283;
      par[1]=9.99995;
      par[2]=0.999027;

    }else if(hypothesis==3){
      par[0]=0.0964021;
      par[1]=-0.214469;
      par[2]=0.989688;

    }else{
      par[0]=0.318674;
      par[1]=-0.197755;
      par[2]=0.968436;
    }

    break;
  case 2:   //pion
    if(hypothesis==0){
      par[0]=-0.0123577;
      par[1]=-0.141521;
      par[2]=1.00273;

    }else if(hypothesis==1){
      par[0]=-0.00558462;
      par[1]=-0.136941;
      par[2]=1.00135;

    }else if(hypothesis==2){
      par[0]=0.0;
      par[1]=1.0;
      par[2]=1.00001;

    }else if(hypothesis==3){
      par[0]=0.122083;
      par[1]=-0.1333923;
      par[2]=0.976863;

    }else{
      par[0]=0.401111;
      par[1]=-0.116807;
      par[2]=0.930906;
    }

    break;
  case 3:    //kaon
    if(hypothesis==0){
      par[0]=-0.102300;
      par[1]=-0.139570;
      par[2]=1.01362;

    }else if(hypothesis==1){
      par[0]=-0.0973257;
      par[1]=-0.138932;
      par[2]=1.01293;

    }else if(hypothesis==2){
      par[0]=-0.0936450;
      par[1]=-0.138469;
      par[2]=1.01242;

    }else if(hypothesis==3){
      par[0]=0.0;
      par[1]=1.0;
      par[2]=0.999865;

    }else{
      par[0]=0.223317;
      par[1]=-0.101273;
      par[2]=0.973484;
    }

    break;
  case 4:   //proton
    if(hypothesis==0){
      par[0]=-0.260150;
      par[1]=-0.0990612;
      par[2]=1.02854;

    }else if(hypothesis==1){
      par[0]=-0.256503;
      par[1]=-0.0976728;
      par[2]=1.02811;

    }else if(hypothesis==2){
      par[0]=-0.253788;
      par[1]=-0.0966732;
      par[2]=1.02779;

    }else if(hypothesis==3){
      par[0]=-0.183031;
      par[1]=-0.0742123;
      par[2]=1.01965;

    }else{
      par[0]=0.0;
      par[1]=1.0;
      par[2]=0.999791;
    }

    break;
  }

  return par[0]/sqrt(p*p+par[1])+par[2]; 
}

