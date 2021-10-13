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

#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"
#include "DD4hep/DD4hepUnits.h"

#include "LikelihoodPID.hh"

using namespace std;

std::string itos(int i)  
{
  std::stringstream s;
  s << i;
  return s.str();
}

LikelihoodPID::LikelihoodPID(double *pars){
  //set parameters for bethebloch
  //electron                                                                                                                   
  for(int i=0;i<5;i++) par[0][i]=pars[i];
  //muon                                                                                                                   
  for(int i=0;i<5;i++) par[1][i]=pars[i+5];
  //pion                                                                                                                   
  for(int i=0;i<5;i++) par[2][i]=pars[i+10];
  //kaon                                                                                                                   
  for(int i=0;i<5;i++) par[3][i]=pars[i+15];
  //proton                                                                                                                   
  for(int i=0;i<5;i++) par[4][i]=pars[i+20];
                                                
  //choose method
  _usebayes=(int)pars[25];
  _dEdxnorm=(float)pars[26];
  _dEdxerrfact=pars[27];

  //set mass
  emass=0.000510998;
  mmass=0.105658;
  pimass=0.139570;
  kmass=0.493677;
  pmass=0.938272;

  //get B field
  double bfield[3]={};
  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
  lcdd.field().magneticField({0.0,0.0,0.0}, bfield);
  _bfield = (float)bfield[2]/dd4hep::tesla;
  
  return;
}

LikelihoodPID::LikelihoodPID(string fname, double *pars, std::vector<float> cost){
  //set parameters for bethebloch
  //electron                                                                                                                   
  for(int i=0;i<5;i++) par[0][i]=pars[i];
  //muon                                                                                                                   
  for(int i=0;i<5;i++) par[1][i]=pars[i+5];
  //pion                                                                                                                   
  for(int i=0;i<5;i++) par[2][i]=pars[i+10];
  //kaon                                                                                                                   
  for(int i=0;i<5;i++) par[3][i]=pars[i+15];
  //proton                                                                                                                   
  for(int i=0;i<5;i++) par[4][i]=pars[i+20];
                                                
  //choose method
  _usebayes=(int)pars[25];
  _dEdxnorm=(float)pars[26];
  _dEdxerrfact=pars[27];

  //set mass
  emass=0.000510998;
  mmass=0.105658;
  pimass=0.139570;
  kmass=0.493677;
  pmass=0.938272;

  //get pdf plots
  string ffstr1 = fname;   // "./pdf_standard_s_05_v01.root";
  fpdf=new TFile(ffstr1.c_str());

  string hname,hname2;
  //for(int i=0;i<6;i++) pdf[i] =new TH1F()[14];
  for(int i=0;i<6;i++){
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

    TObject *obj;
    //cluster depth
    hname="hcludepth" + itos(i+1);
    fpdf->GetObject(hname.c_str(), obj);
    if(obj){
      hname2="hcludepth" + itos(i+1) + "_2";
      pdf[i][14]=(TH1F*)obj->Clone(hname2.c_str());

      //deltax
      hname="hdeltax" + itos(i+1);
      fpdf->GetObject(hname.c_str(), obj);
      hname2="hdeltax" + itos(i+1) + "_2";
      pdf[i][12]=(TH1F*)obj->Clone(hname2.c_str());

      //deltaz
      hname="hdeltaz" + itos(i+1);
      fpdf->GetObject(hname.c_str(), obj);
      hname2="hdeltaz" + itos(i+1) + "_2";
      pdf[i][13]=(TH1F*)obj->Clone(hname2.c_str());

      //hitmean
      hname="hhitmean" + itos(i+1);
      fpdf->GetObject(hname.c_str(), obj);
      hname2="hhitmean" + itos(i+1) + "_2";
      pdf[i][15]=(TH1F*)obj->Clone(hname2.c_str());

      //hitrms
      hname="hhitrms" + itos(i+1);
      fpdf->GetObject(hname.c_str(), obj);
      hname2="hhitrms" + itos(i+1) + "_2";
      pdf[i][16]=(TH1F*)obj->Clone(hname2.c_str());
    }else{
      pdf[i][12]=NULL;  //backward compatibility
      pdf[i][13]=NULL;  //backward compatibility
      pdf[i][14]=NULL;  //backward compatibility
      pdf[i][15]=NULL;  //backward compatibility
      pdf[i][16]=NULL;  //backward compatibility
    }
    
    //deltaz
    /*hname="hdeltaz" + itos(i+1);
      hname2="hdeltaz" + itos(i+1) + "_2";
      pdf[i][13]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
      //cluster depth
      hname="hcludepth" + itos(i+1);
      hname2="hcludepth" + itos(i+1) + "_2";
      pdf[i][14]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
      //hitmean
      hname="hhitmean" + itos(i+1);
      hname2="hhitmean" + itos(i+1) + "_2";
      pdf[i][15]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
      //hitrms
      hname="hhitrms" + itos(i+1);
      hname2="hhitrms" + itos(i+1) + "_2";
      pdf[i][16]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());*/
  }

  //normalize histograms
  double weight=1.0;
  for(int i=0;i<6;i++){
    for(int j=0;j<17;j++){
      //normalize histograms
      if(pdf[i][j]!=NULL){
	weight=pdf[i][j]->Integral(0,pdf[i][j]->GetNbinsX()+1,"");
	pdf[i][j]->Scale(1.0/weight);
	_weights[i][j]=weight;  //for hadron likelihood calculation
      }
    }
  }
    
  //normalize weights
  for(int j=0;j<17;j++){
    if(pdf[2][j]!=NULL && pdf[3][j]!=NULL && pdf[4][j]!=NULL){
      double denom = _weights[2][j]+_weights[3][j]+_weights[4][j];
      for(int i=2;i<5;i++) _weights[i][j] = _weights[i][j]/denom;
    }else
      for(int i=2;i<5;i++) _weights[i][j] = 1.0 / 3.0;  //same weights
  }

  //set threshold
  //this is original
  threshold[0]=0.0;  //TMath::Exp(-0.55);
  threshold[1]=0.0;  //TMath::Exp(-0.7);
  threshold[2]=0.0; 
  threshold[3]=0.0;  //TMath::Exp(-0.6); 
  threshold[4]=0.0;  //TMath::Exp(-0.6);
  
  //set cost
  if(_usebayes==1){
    for(int i=0;i<5;i++){
      for(int j=0;j<5;j++){
	if(i==j) penalty[i][j]=1.0e-50;
	else penalty[i][j]=1.0;
      }
    }
  }else if(_usebayes==2){
    for(int i=0;i<5;i++){ 
      for(int j=0;j<5;j++){
	penalty[i][j]=cost[i*5+j];
      }
    }
  }

  //get z component of B field
  double bfield[3]={};
  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
  lcdd.field().magneticField({0.0,0.0,0.0}, bfield);
  _bfield = (float)bfield[2]/dd4hep::tesla;
  //cout << "magnetic field: " << _bfield << endl;

  return;
}

LikelihoodPID::~LikelihoodPID(){

  fpdf->Close();

  return;
}

//public
int LikelihoodPID::Classification(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec){
  int tmpid=-1;

  //set prior first
  for(int i=0;i<5;i++){
    prior[i]=0.2;
    _posterior[i]=prior[i];
    _likelihood[i]=0.0;
  }
  _posterior[5]=prior[2]+prior[3]+prior[4];
  _likelihood[5]=0.0;

  //get dedx distance
  for(int i=0;i<5;i++) _dEdxDist[i]=get_dEdxChi2(i, pp.Vect(), trk->getdEdxError(),trk->getdEdx());

  //first check whether electron or not
  for(int i=0;i<5;i++){
    _likelihood[i]=0.0;
    if(_basicFlg==false && _dEdxFlg==true && _showerShapesFlg==false) prior[i]=0.2;
  }
  _likelihood[5]=0.0;  //hadron type
  
  tmpid=Class_electron(pp, trk, cluvec);
  if(tmpid==0) return tmpid;

  //second check whether muon or not
  for(int i=0;i<5;i++){
    _likelihood[i]=0.0;
    if(_basicFlg==false && _dEdxFlg==true && _showerShapesFlg==false) prior[i]=0.2;
  }
  _likelihood[5]=0.0;  //hadron type

  tmpid=Class_muon(pp, trk, cluvec);
  if(tmpid==1) return tmpid;

  //third classify hadrons
  for(int i=0;i<5;i++){
    _likelihood[i]=0.0;
    if(_basicFlg==false && _dEdxFlg==true && _showerShapesFlg==false) prior[i]=0.2;
  }
  _likelihood[5]=0.0;  //hadron type

  tmpid=Class_hadron(pp, trk, cluvec);

  //avoid strange value
  if(_posterior[0]!=_posterior[0]){  //cannot estimate likelihood and probability
    for(int i=0;i<6;i++){
      _likelihood[i] = 999.0;
      _posterior[i] = 0.0;
    }
  }else if(_likelihood[0]==0.0 &&_likelihood[1]==0.0 &&_likelihood[2]==0.0 && _likelihood[3]==0.0 && _likelihood[4]==0.0){
    for(int i=0;i<6;i++){
      _likelihood[i] = 999.0;
      _posterior[i] = 0.0;
    }
  }

  return tmpid;
}

double *LikelihoodPID::GetPosterior(){

  return _posterior;
}

double *LikelihoodPID::GetLikelihood(){

  return _likelihood;
}

//track energy correction
double LikelihoodPID::getCorrEnergy(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec){
  int parttype=Classification(pp, trk, cluvec);
  if(parttype==-1) parttype=2;  //move to pion class

  return getCorrEnergy(pp, parttype);
}

double LikelihoodPID::getCorrEnergy(TLorentzVector pp, int parttype){
  if(parttype==-1) parttype=2;  //move to pion class

  double tmpmass=0.0;
  if(parttype==0) tmpmass=emass;
  if(parttype==1) tmpmass=mmass;
  if(parttype==2) tmpmass=pimass;
  if(parttype==3) tmpmass=kmass;
  if(parttype==4) tmpmass=pmass;

  return sqrt(pp.P()*pp.P()+tmpmass*tmpmass);
}

//private
double LikelihoodPID::get_dEdxChi2(int parttype, TVector3 p, float hit, double dEdx){
  //get parameters for chi2
  double tmppar[5],tmpmass=0.0;
  for(int i=0;i<5;i++) tmppar[i]=par[parttype][i];
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
  // double trkcos=p.CosTheta();

  //get nomalized dEdx
  double dEdx_Norm=get_Norm(dEdx);

  //get expected dE/dx
  double ExpdEdx=BetheBloch(p.Mag(),tmpmass,tmppar);

  //cout << "check " << emass << " " << tmpmass << " " << trkcos << " " << dEdx_Norm << " " << ExpdEdx << endl;
  double dEdx_Error = _dEdxerrfact * dEdx_Norm * hit/dEdx;    //change 20151218

  //get chi2!!
  double chi2=TMath::Power((dEdx_Norm-ExpdEdx)/dEdx_Error,2.0);
  if(dEdx_Norm-ExpdEdx<0.0) chi2=-chi2;    //get signed chi2

  return chi2;
}

double LikelihoodPID::get_dEdxFactor(int /*parttype*/, TVector3 /*p*/, float hit, double dEdx){
  //get parameters for chi2
  // double tmppar[5],tmpmass=0.0;
  // for(int i=0;i<5;i++) tmppar[i]=par[parttype][i];
  // //getmass
  // switch(parttype){
  // case 0:
  //   tmpmass=emass;
  //   break;
  // case 1:
  //   tmpmass=mmass;
  //   break;
  // case 2:
  //   tmpmass=pimass;
  //   break;
  // case 3:
  //   tmpmass=kmass;
  //   break;
  // case 4:
  //   tmpmass=pmass;
  //   break;
  // }

  //cal. polar angle
  // double trkcos=p.CosTheta();

  //get nomalized dEdx
  double dEdx_Norm=get_Norm(dEdx);

  //cout << "check " << emass << " " << tmpmass << " " << trkcos << " " << dEdx_Norm << " " << ExpdEdx << endl;
  double dEdx_Error = _dEdxerrfact * dEdx_Norm * hit/dEdx;    //change 20151218

  //get likelihood factor
  double factor=-0.5*TMath::Log(2.0*TMath::Pi()*dEdx_Error*dEdx_Error);

  return factor;
}

//public
double LikelihoodPID::get_dEdxDist(int parttype){
  //get parameters for dedxdist
  double dedxdist=0.0;

  if(_dEdxDist[parttype]!=0.0) dedxdist = sqrt(fabs(_dEdxDist[parttype]))*_dEdxDist[parttype]/fabs(_dEdxDist[parttype]);
  if(dedxdist!=dedxdist) dedxdist=0.0;
  return dedxdist;
}

double LikelihoodPID::get_Norm(double dedx){
  return dedx/_dEdxnorm;
}

double LikelihoodPID::BetheBloch(double x, double mass, double *pars){
  double bg=x/mass;
  double b=sqrt(bg*bg/(1.0+bg*bg));
  //double g=bg/b;
  double tmax=pars[2]*TMath::Power(bg,2.0);   ///(1.0+pars[3]*g+pars[4]);

  return (0.5*pars[0]*TMath::Log(pars[1]*TMath::Power(bg,2.0)*tmax)-pars[3]*b*b-pars[4]*bg/2.0)/(b*b);
}
 
int LikelihoodPID::Class_electron(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec){
  //Bayesian approach to classify electron!
  
  //get deposit energy
  double ecal=0.0,hcal=0.0,mucal=0.0;
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
  double var[11]={}; //11
  //for track variables
  var[0]=(ecal+hcal)/pp.P();
  if(ecal+hcal!=0.0) var[1]=ecal/(ecal+hcal);
  if(shapes.size()!=0){
    var[2]=shapes[0+4];
    var[3]=shapes[5+4];
    var[4]=fabs(shapes[3+4])/(shapes[6+4]);
    var[5]=shapes[15+4]/(2.0*3.50);
  }else var[2]=-1.0;

  var[6]=get_dEdxChi2(0,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[7]=get_dEdxChi2(1,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[8]=get_dEdxChi2(2,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[9]=get_dEdxChi2(3,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[10]=get_dEdxChi2(4,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2

  var[6]=-0.5*fabs(var[6])+get_dEdxFactor(0,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[7]=-0.5*fabs(var[7])+get_dEdxFactor(1,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[8]=-0.5*fabs(var[8])+get_dEdxFactor(2,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[9]=-0.5*fabs(var[9])+get_dEdxFactor(3,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[10]=-0.5*fabs(var[10])+get_dEdxFactor(4,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood

  //get likelihood for each class
  double posterior[5]={0.0,0.0,0.0,0.0,0.0};
  double risk[5]={0.0,0.0,0.0,0.0,0.0};
  double okval[5]={0.0,0.0,0.0,0.0,0.0};
  double total=0.0;
  double priorprob[5]={0.0,0.0,0.0,0.0,0.0};
  int valtype=0;
  for(int j=0;j<7;j++){   //variables
    //avoid very low momentum tracks (no energy deposit in the cal.)
    if(var[0]==0.0 && j<=5) continue;

    //first, get likelihood
    for(int i=0;i<5;i++){   //particle type
      if(j==0) valtype=0;
      if(j==1) valtype=1;
      if(j==2) valtype=3;
      if(j==3) valtype=4;
      if(j==4) valtype=5;
      if(j==5) valtype=6;
      if(j==6) valtype=7;
      
      if(j<6){
	okval[i]=getValue(i,valtype,var[j]);   //likelihood
	if(mucal==0.0 && i==1) okval[i]=getValue(5,valtype,var[j]);   //likelihood
      }else if(j==6)
	okval[i]=std::max(std::exp(var[6+i]), 1.0e-300);   //just use dE/dx likelihood 
    }
    
    //for basic variables flg
    if(!_basicFlg && j<2) continue;
    //for cluster shape flg
    if(!_showerShapesFlg && (j>=2 && j<=5)) continue;
    //for dEdx flg
    if(!_dEdxFlg && j==6) continue; 
    //don't use some dEdxin the case of LikelihoodPID
    if(_basicFlg && _showerShapesFlg && _dEdxFlg && !(j<=6)) continue;

    //to avoid strange value for dE/dx
    //if(j>=6 && var[j]<=-50.0) continue;
    if(j>=6 && var[j]!=var[j]) continue;
    
    //cal. prior probability
    for(int i=0;i<5;i++)   //particle type
      priorprob[i]=prior[i]/
	(prior[0]+prior[1]+prior[2]+prior[3]+prior[4]);
        
    //cal total probability for each particle type
    total=0.0;
    for(int i=0;i<5;i++) total+=priorprob[i]*okval[i];
    
    //cal. log posterior probability!
    for(int i=0;i<5;i++){   //particle type
      ///20150708 change it from posterior probability to likelihood
      posterior[i]=TMath::Log(okval[i])+TMath::Log(priorprob[i])-TMath::Log(total);
      _likelihood[i]+=TMath::Log(okval[i]);
    }
    //hadron type likelihood
    _likelihood[5] += TMath::Log(_weights[2][valtype]*okval[2] + _weights[3][valtype]*okval[3] + _weights[4][valtype]*okval[4]);
    
    //bayesian updating
    for(int i=0;i<5;i++) prior[i]=TMath::Exp(posterior[i]);

  }

  //cal. risk
  //first. get penalty
  // for(int i=0;i<5;i++){   //particle type
  //   for(int j=0;j<5;j++){
  //     penalty[i][j]=1.0;  
  //     if(i==j) penalty[i][j]=1.0e-50;
  //   }
  // }

  for(int i=0;i<5;i++){   //particle type
    risk[i]=penalty[i][0]*TMath::Exp(posterior[0])+penalty[i][1]*TMath::Exp(posterior[1])
      +penalty[i][2]*TMath::Exp(posterior[2])+penalty[i][3]*TMath::Exp(posterior[3])+penalty[i][4]*TMath::Exp(posterior[4]);
  }

  //check the rule
  int okflg=-1;
  double tmppp=-1.0e+100;
  if(_usebayes==0){
    //---------- here is likelihood based determination rule -----------
    for(int i=0;i<5;i++){
      tmppp=TMath::Max(_likelihood[i], tmppp);
      //if(posterior[0]<posterior[i]){
      //  okflg=false;
      //  break;
      //}
    }
    
    for(int i=0;i<5;i++){
      //if(fabs(tmppp-posterior[i])<1.0e-6) okflg=i;
      if(fabs(tmppp-_likelihood[i])<1.0e-6) okflg=i;
    }
    //------------------------------------------------------------------
  }else{
    //---------- here is risk based determination rule ----------- 
    tmppp=1.0e+100;
    for(int i=0;i<5;i++){
      tmppp=TMath::Min(risk[i], tmppp);
    }
    
    for(int i=0;i<5;i++){
      if(fabs(tmppp-risk[i])<1.0e-8) okflg=i;
    }

    //threshold is set.
    if(okflg==0 && posterior[okflg]<TMath::Log(threshold[okflg])) okflg=-1;
    //------------------------------------------------------------------
  }

  //save posterior
  _posterior[0]=TMath::Exp(posterior[0]);
  _posterior[1]=TMath::Exp(posterior[1]);
  _posterior[2]=TMath::Exp(posterior[2]);
  _posterior[3]=TMath::Exp(posterior[3]);
  _posterior[4]=TMath::Exp(posterior[4]);
  //hadron type
  _posterior[5]=_posterior[2]+_posterior[3]+_posterior[4];
  
  return okflg;
}

int LikelihoodPID::Class_muon(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec){
  //Bayesian approach to classify muon!
  
  //get deposit energy
  double ecal=0.0,hcal=0.0,mucal=0.0;
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
  double var[17]={}; //17 variales used(is it OK?)
  var[0]=(ecal+hcal)/pp.P();
  if(ecal+hcal!=0.0) var[1]=ecal/(ecal+hcal);
  var[2]=mucal;

  if(shapes.size()!=0){
    var[3]=shapes[0+4];
    var[4]=shapes[5+4];
    var[5]=fabs(shapes[3+4])/(shapes[6+4]);
    var[6]=shapes[15+4]/(2.0*3.50);

    //for low momentum mu/pi separation
    var[14]= shapes[17+4];
    var[15]= shapes[18+4];
    var[16]= shapes[19+4];
  }else var[3]=-1.0;
  
  var[7]=get_dEdxChi2(0,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[8]=get_dEdxChi2(1,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[9]=get_dEdxChi2(2,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[10]=get_dEdxChi2(3,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[11]=get_dEdxChi2(4,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2

  var[7]=-0.5*fabs(var[7])+get_dEdxFactor(0,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[8]=-0.5*fabs(var[8])+get_dEdxFactor(1,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[9]=-0.5*fabs(var[9])+get_dEdxFactor(2,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[10]=-0.5*fabs(var[10])+get_dEdxFactor(3,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[11]=-0.5*fabs(var[11])+get_dEdxFactor(4,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood

  var[12] = _delpos[0];
  var[13] = _delpos[2];

  //get likelihood for each class
  double posterior[5]={0.0,0.0,0.0,0.0,0.0};
  double risk[5]={0.0,0.0,0.0,0.0,0.0};
  int valtype=0;
  double okval[5]={0.0,0.0,0.0,0.0,0.0};
  double total=0.0;     //[5]={0.0,0.0,0.0,0.0,0.0};
  double priorprob[5]={0.0,0.0,0.0,0.0,0.0};
  for(int j=0;j<13;j++){   //variables

    if(var[0]>0.0){
      //avoid pion misID when energy deposit to mucal is zero
      if(var[2]==0.0 && j==2) continue;   //(j==2 || j==7 || j>=9)) continue;
      //don't use when ep>0.0
      //else if(var[2]>0.0 && j>=9) continue;
    }else{
      //avoid very low momentum tracks (no energy deposit in the cal.)
      if(j<=6) continue;
      if(j>7) continue;
    }

    //first, get likelihood
    for(int i=0;i<5;i++){   //particle type(electron isn't be checked)
      if(j==0) valtype=0;
      if(j==1) valtype=1;
      if(j==2) valtype=2;
      if(j==3) valtype=3;
      if(j==4) valtype=4;
      if(j==5) valtype=5;
      if(j==6) valtype=6;
      if(j==7) valtype=7;
      if(j==8) valtype=11;
      if(j==9) valtype=12;
      if(j==10) valtype=13;
      if(j==11) valtype=14;
      if(j==12) valtype=15;
      
      if(j<7){
	okval[i]=getValue(i,valtype,var[j]);   //likelihood
	if(var[2]==0.0 && i==1) okval[i]=getValue(5,valtype,var[j]);   //likelihood 
      }else if(j==7){
	okval[i]=std::max(std::exp(var[7+i]), 1.0e-300);   //just use dE/dx likelihood 
      }else if(j==8 || j==9){
	okval[i]=getValue(i,valtype,var[j+4]);   //likelihood
	if(var[2]==0.0 && i==1) okval[i]=getValue(5,valtype,var[j+4]);   //likelihood 
      }else if(j>=10 && pp.P()>0.3 && pp.P()<6.0){  //for low momentum mu/pi
	okval[i]=getValue(i,valtype,var[j+4]);   //likelihood
	if(var[2]==0.0 && i==1) okval[i]=getValue(5,valtype,var[j+4]);   //likelihood 
      }
    }

    //for basic variables flg
    if(!_basicFlg && j<=2) continue;
    //for cluster shape flg
    if(!_showerShapesFlg && ((j>2 && j<=6) || j>=8)) continue;
    //for dEdx flg
    if(!_dEdxFlg && j==7) continue; 
    //don't use some dEdxin the case of LikelihoodPID
    if(_basicFlg && _showerShapesFlg && _dEdxFlg && !(j<=9 || (j>=10 && pp.P()>0.3 && pp.P()<6.0))) continue;
    if(mucal!=0.0 && j==10) continue;

    //to avoid strange value for dE/dx
    //if(j>=7 && j<=11 && var[j]<=-50.0) continue;
    if(j>=7 && j<=11 && var[j]!=var[j]) continue;
    
    //cal. prior probability
    for(int i=0;i<5;i++)   //particle type
      priorprob[i]=prior[i]/
	(prior[0]+prior[1]+prior[2]+prior[3]+prior[4]);

    //cal total probability for each particle type
    total=0.0;
    for(int i=0;i<5;i++) total+=priorprob[i]*okval[i];
    
    //cal. log likelihood first
    for(int i=0;i<5;i++){   //particle type
      //20150708 change it from posterior probability to simple likelihood
      posterior[i]=TMath::Log(okval[i])+TMath::Log(priorprob[i])-TMath::Log(total); 
      _likelihood[i]+=TMath::Log(okval[i]);
    }
    //hadron type likelihood
    _likelihood[5] += TMath::Log(_weights[2][valtype]*okval[2] + _weights[3][valtype]*okval[3] + _weights[4][valtype]*okval[4]);
    
    //bayesian updating
    for(int i=0;i<5;i++) prior[i]=TMath::Exp(posterior[i]);
  }

  //cal. risk
  //first. get penalty
  // for(int i=0;i<5;i++){   //particle type
  //   for(int j=0;j<5;j++){
  //     penalty[i][j]=1.0; 
  //     if(i==j) penalty[i][j]=1.0e-50; 
  //   }
  // }
  // penalty[1][2] = 9.0;   //TMath::Min(1.0, 9.0*pp.P(); Todo: momentum depepndence  

  for(int i=0;i<5;i++){   //particle type
    risk[i]=penalty[i][0]*TMath::Exp(posterior[0])+penalty[i][1]*TMath::Exp(posterior[1])
      +penalty[i][2]*TMath::Exp(posterior[2])+penalty[i][3]*TMath::Exp(posterior[3])+penalty[i][4]*TMath::Exp(posterior[4]);
  }

  //check the rule
  int okflg=-1;
  double tmppp=-1.0e+100;
  if(_usebayes==0){
    //---------- here is likelihood based determination rule -----------
    for(int i=0;i<5;i++){
      tmppp=TMath::Max(_likelihood[i], tmppp);
      //if(posterior[0]<posterior[i]){
      //  okflg=false;
      //  break;
      //}
    }
    
    for(int i=0;i<5;i++){
      //if(fabs(tmppp-posterior[i])<1.0e-6) okflg=i;
      if(fabs(tmppp-_likelihood[i])<1.0e-6) okflg=i;
    }
    //------------------------------------------------------------------
  }else{
    //---------- here is risk based determination rule ----------- 
    tmppp=1.0e+100;
    for(int i=0;i<5;i++){
      tmppp=TMath::Min(risk[i], tmppp);
    }
    
    for(int i=0;i<5;i++){
      if(fabs(tmppp-risk[i])<1.0e-8) okflg=i;
    }

    if(okflg==1 && posterior[okflg]<TMath::Log(threshold[okflg])) okflg=-1;   //threshold is set.
    //------------------------------------------------------------------
  }

  //save posterior
  _posterior[0]=TMath::Exp(posterior[0]);
  _posterior[1]=TMath::Exp(posterior[1]);
  _posterior[2]=TMath::Exp(posterior[2]);
  _posterior[3]=TMath::Exp(posterior[3]);
  _posterior[4]=TMath::Exp(posterior[4]);
  //hadron type
  _posterior[5]=_posterior[2]+_posterior[3]+_posterior[4];

  return okflg;  
}

int LikelihoodPID::Class_hadron(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec){
  //Bayesian approach to classify hadrons!

  //get deposit energy
  double ecal=0.0,hcal=0.0,mucal=0.0;
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
  double var[10]; //10 variales used(is it OK?)
  var[0]=(ecal+hcal)/pp.P();
  if(ecal+hcal!=0.0)var[1]=ecal/(ecal+hcal);
  if(shapes.size()!=0){
    var[2]=shapes[0+4];
    var[3]=fabs(shapes[3+4])/(shapes[6+4]);
    var[4]=shapes[15+4]/(2.0*3.50);
  }else var[2]=-1.0;
 
  var[5]=get_dEdxChi2(0,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[6]=get_dEdxChi2(1,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[7]=get_dEdxChi2(2,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[8]=get_dEdxChi2(3,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2
  var[9]=get_dEdxChi2(4,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //chi2

  var[5]=-0.5*fabs(var[5])+get_dEdxFactor(0,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[6]=-0.5*fabs(var[6])+get_dEdxFactor(1,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[7]=-0.5*fabs(var[7])+get_dEdxFactor(2,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[8]=-0.5*fabs(var[8])+get_dEdxFactor(3,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood
  var[9]=-0.5*fabs(var[9])+get_dEdxFactor(4,pp.Vect(),trk->getdEdxError(),trk->getdEdx());   //log likelihood

  //get likelihood for each class
  double posterior[5]={TMath::Log(0.20),TMath::Log(0.20),TMath::Log(0.20),TMath::Log(0.20),TMath::Log(0.20)};
  double risk[5]={0.0,0.0,0.0,0.0,0.0};
  int valtype=0;
  double okval[5]={0.0,0.0,0.0,0.0,0.0};
  double total[5]={0.0,0.0,0.0,0.0,0.0};
  double priorprob[5]={0.0,0.0,0.0,0.0,0.0};
  for(int j=0;j<6;j++){   //variables
    //avoid very low momentum tracks (no energy deposit in the cal.)
    if(var[0]==0.0 && j<=4) continue;

    //first, get likelihood
    for(int i=0;i<5;i++){   //particle type(electron&muon aren't be checked)
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
      
      if(j<5){
	okval[i]=getValue(i,valtype,var[j]);   //likelihood
      }else if(j>=5)
	okval[i]=std::max(std::exp(var[5+i]), 1.0e-300);   //just use dE/dx likelihood 
    }

    //for basic variables flg
    if(!_basicFlg && j<2) continue;
    //for cluster shape flg
    if(!_showerShapesFlg && j>=2 && j<=4) continue;
    //for dEdx flg
    if(!_dEdxFlg && j>=5) continue; 
    //don't use some dEdxin the case of LikelihoodPID
    if(_basicFlg && _showerShapesFlg && _dEdxFlg && !(j<=5)) continue;

    //to avoid strange value for dE/dx
    //if(j>=5 && var[j]<=-50.0) continue;
    if(j>=5 && var[5]!=var[5]) continue;
    
    //cal. prior probability
    for(int i=0;i<5;i++)   //particle type
      priorprob[i]=prior[i]/
	(prior[0]+prior[1]+prior[2]+prior[3]+prior[4]);
        
    //cal total probability
    for(int i=0;i<5;i++){   //particle type
      total[i]=priorprob[0]*okval[0]+priorprob[1]*okval[1]+priorprob[2]*okval[2]+priorprob[3]*okval[3]+priorprob[4]*okval[4];
      if(total[i]==0.0) total[i]=1.0e-100;
    }

    //cal. log posterior probability!
    for(int i=0;i<5;i++){   //particle type
      //20150708 change it from posterior probabity to simple likelihood
      posterior[i]=TMath::Log(okval[i])+TMath::Log(priorprob[i])-TMath::Log(total[i]);
      _likelihood[i]+=TMath::Log(okval[i]);
    }
    //hadron type likelihood
    _likelihood[5] += TMath::Log(_weights[2][valtype]*okval[2] + _weights[3][valtype]*okval[3] + _weights[4][valtype]*okval[4]);
    
    //bayesian updating
    for(int i=0;i<5;i++) prior[i]=TMath::Exp(posterior[i]);
  }
    
  //cal. risk
  //normalize
  // for(int i=0;i<5;i++){ 
  //   for(int j=0;j<5;j++){
  //     penalty[i][j]=1.0; 
  //     if(i==j) penalty[i][j]=1.0e-50;
  //   }
  // }
  //first. get penalty
  // for(int i=2;i<5;i++){   //particle type
  //   for(int j=2;j<5;j++){
  //     penalty[i][j]=getPenalty(i,j,pp.P());
  //   }
  //   double tmppenalty=penalty[i][i];
  //   for(int j=2;j<5;j++){
  //     penalty[i][j]=(fabs(penalty[i][j]-tmppenalty)+tmppenalty)/tmppenalty;
  //     penalty[i][j]=penalty[i][j]-1.0;
      
  //   }
  // }

  //weight of the cost
  //reserve default value
  double res43=penalty[4][3];
  double res32=penalty[3][2];
  if(_usebayes==2){
    if(pp.P()>1.0){
      penalty[4][3] = max(1.0, penalty[4][3]*(1.0-TMath::Exp(-0.772*pp.P())));  //10.0
      //penalty[3][2] = max(1.0, penalty[3][2]*TMath::Exp(0.0268*pp.P()));  //30.0
    }else{
      penalty[4][3] = 1.0;
      penalty[3][2] = 1.0;
    }

    if(fabs(pp.CosTheta())<1.0){
      penalty[4][3] = max(1.0, penalty[4][3]*(0.439*pp.CosTheta()*pp.CosTheta()-0.940*fabs(pp.CosTheta())+0.495));
      penalty[3][2] = max(1.0, penalty[3][2]*TMath::Exp(-1.46*fabs(pp.CosTheta())));  //30.0
    }

  }
  
  //penalty[4][2] = 5.0;

  for(int i=0;i<5;i++){   //particle type
    risk[i]=penalty[i][0]*TMath::Exp(posterior[0])+penalty[i][1]*TMath::Exp(posterior[1])
      +penalty[i][2]*TMath::Exp(posterior[2])+penalty[i][3]*TMath::Exp(posterior[3])+penalty[i][4]*TMath::Exp(posterior[4]);
  }

  penalty[4][3]=res43;
  penalty[3][2]=res32;

  //check the rule
  int okflg=-1;
  double tmppp=-1.0e+100;
  if(_usebayes==0){
    //---------- here is likelihood based determination rule -----------
    for(int i=0;i<5;i++){
      tmppp=TMath::Max(_likelihood[i], tmppp);
      //if(posterior[0]<posterior[i]){
      //  okflg=false;
      //  break;
      //}
    }
    
    for(int i=0;i<5;i++){
      //if(fabs(tmppp-posterior[i])<1.0e-6) okflg=i;
      if(fabs(tmppp-_likelihood[i])<1.0e-6) okflg=i;
    }

    //move to pion when very bad
    //if(TMath::Exp(tmppp)<0.21) okflg=2;
    if(tmppp==0.0) okflg=2;
    //------------------------------------------------------------------
  }else{
    //---------- here is risk based determination rule ----------- 
    tmppp=1.0e+100;
    for(int i=0;i<5;i++){
      tmppp=TMath::Min(risk[i], tmppp);
    }
    
    for(int i=0;i<5;i++){
      if(fabs(tmppp-risk[i])<1.0e-6) okflg=i;
    }
    if(okflg<2) okflg=2;
    if(posterior[okflg]<=TMath::Log(0.21)) okflg=TMath::Max(2,okflg);  //do not go to leptons

    //check penalty -1 means undefined
    if(okflg>=2 && posterior[okflg]<TMath::Log(threshold[okflg])) okflg=-1;   //threshold is set.
    //------------------------------------------------------------------
  }

  //save posterior
  _posterior[0]=TMath::Exp(posterior[0]);
  _posterior[1]=TMath::Exp(posterior[1]);
  _posterior[2]=TMath::Exp(posterior[2]);
  _posterior[3]=TMath::Exp(posterior[3]);
  _posterior[4]=TMath::Exp(posterior[4]);
  //hadron type
  _posterior[5]=_posterior[2]+_posterior[3]+_posterior[4];

  return okflg;  
}

double LikelihoodPID::getValue(int type, int valtype, double value){

  /*int nbins=pdf[type][valtype]->GetNbinsX();
  double interval=pdf[type][valtype]->GetBinWidth(1);
  double minval=pdf[type][valtype]->GetBinLowEdge(1);

  double val=1.0e-30;
  int bin=0;
  for(int i=0;i<nbins;i++){
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

  
  double val=1.0e-100;
  int bin=0;

  if(pdf[type][valtype]!=NULL){
    bin=pdf[type][valtype]->GetXaxis()->FindBin(value);
    
    //get probability
    val=pdf[type][valtype]->GetBinContent(bin);
    if(val==0.0) val= 1.0e-100;
  }else
    val = 1.0;

  return val;
}

double LikelihoodPID::getPenalty(int ptype, int hypothesis, double p){
  double param[3]={0.0,0.0,0.0};
  //set parameters
  switch(ptype){
  case 0:  //electron
    if(hypothesis==0){
      param[0]=0.0;
      param[1]=1.0;
      param[2]=1.00178;

    }else if(hypothesis==1){
      param[0]=0.101945;
      param[1]=2.18719;
      param[2]=0.990000;

    }else if(hypothesis==2){
      param[0]=0.104521;
      param[1]=1.98490;
      param[2]=0.990000;

    }else if(hypothesis==3){
      param[0]=0.234299;
      param[1]=0.276835;
      param[2]=0.965973;

    }else{
      param[0]=0.414613;
      param[1]=-0.132530;
      param[2]=0.949679;

    }

    break;
  case 1:   //muon
    if(hypothesis==0){
      param[0]=-0.00512190;
      param[1]=-0.211835;
      param[2]=1.00024;

    }else if(hypothesis==1){
      param[0]=0.0;
      param[1]=1.00;
      param[2]=0.999684;

    }else if(hypothesis==2){
      param[0]=0.00833283;
      param[1]=9.99995;
      param[2]=0.999027;

    }else if(hypothesis==3){
      param[0]=0.0964021;
      param[1]=-0.214469;
      param[2]=0.989688;

    }else{
      param[0]=0.318674;
      param[1]=-0.197755;
      param[2]=0.968436;
    }

    break;
  case 2:   //pion
    if(hypothesis==0){
      param[0]=-0.0123577;
      param[1]=-0.141521;
      param[2]=1.00273;

    }else if(hypothesis==1){
      param[0]=-0.00558462;
      param[1]=-0.136941;
      param[2]=1.00135;

    }else if(hypothesis==2){
      param[0]=0.0;
      param[1]=1.0;
      param[2]=1.00001;

    }else if(hypothesis==3){
      param[0]=0.122083;
      param[1]=-0.1333923;
      param[2]=0.976863;

    }else{
      param[0]=0.401111;
      param[1]=-0.116807;
      param[2]=0.930906;
    }

    break;
  case 3:    //kaon
    if(hypothesis==0){
      param[0]=-0.102300;
      param[1]=-0.139570;
      param[2]=1.01362;

    }else if(hypothesis==1){
      param[0]=-0.0973257;
      param[1]=-0.138932;
      param[2]=1.01293;

    }else if(hypothesis==2){
      param[0]=-0.0936450;
      param[1]=-0.138469;
      param[2]=1.01242;

    }else if(hypothesis==3){
      param[0]=0.0;
      param[1]=1.0;
      param[2]=0.999865;

    }else{
      param[0]=0.223317;
      param[1]=-0.101273;
      param[2]=0.973484;
    }

    break;
  case 4:   //proton
    if(hypothesis==0){
      param[0]=-0.260150;
      param[1]=-0.0990612;
      param[2]=1.02854;

    }else if(hypothesis==1){
      param[0]=-0.256503;
      param[1]=-0.0976728;
      param[2]=1.02811;

    }else if(hypothesis==2){
      param[0]=-0.253788;
      param[1]=-0.0966732;
      param[2]=1.02779;

    }else if(hypothesis==3){
      param[0]=-0.183031;
      param[1]=-0.0742123;
      param[2]=1.01965;

    }else{
      param[0]=0.0;
      param[1]=1.0;
      param[2]=0.999791;
    }

    break;
  }

  return param[0]/sqrt(p*p+param[1])+param[2]; 
}

void LikelihoodPID::CalculateDeltaPosition(float charge, TVector3 p, const float* calpos){
  //calculate some variables for preparation
  //transverse momentum
  double prphi=sqrt(pow(p[0],2)+pow(p[1],2));
  //momentum
  double pp=p.Mag();
  //radius in r-phi plane
  double radius=prphi/(0.3*_bfield);
  //tangent lambda
  double tanlam=p[2]/prphi;
  //change to unit vector
  p[0]=p[0]/prphi;
  p[1]=p[1]/prphi;
  p[2]=p[2]/pp;

  //chenge position vector to meter
  float calcalpos[3];
  calcalpos[0]=calpos[0]/1000.0;
  calcalpos[1]=calpos[1]/1000.0;
  calcalpos[2]=calpos[2]/1000.0;


  //radius at the tracker end
  double rradius=sqrt(pow(calcalpos[0],2)+pow(calcalpos[1],2));

  //cout << "check val. " << prphi << " " << radius << " " << tanlam
  //     << " " << rradius << endl;

  //cal. the position of the center of track circle
  TVector3 cc;
  cc[0]=-charge*p[1]*prphi/(0.3*_bfield);
  cc[1]=charge*p[0]*prphi/(0.3*_bfield);
  cc[2]=radius*tanlam;

  //cal. sign and cosine 2theta
  double sintheta=charge*rradius/(2*radius);
  double costheta=sqrt(1-sintheta*sintheta);
  double sin2theta=2*sintheta*costheta;
  double cos2theta=costheta*costheta-sintheta*sintheta;

  TVector3 calpos2;
  calpos2[2]=rradius*tanlam;
  calpos2[0]=cc[0]*(1-cos2theta)+cc[1]*sin2theta;
  calpos2[1]=-cc[0]*sin2theta+cc[1]*(1-cos2theta);

  //cal. difference of the position
  //float delpos[3];
  _delpos[0]=charge*fabs(calcalpos[0]-calpos2[0]);
  _delpos[1]=charge*fabs(calcalpos[1]-calpos2[1]);
  _delpos[2]=fabs(calcalpos[2]-calpos2[2]);


  //cout << "calpos: " << calcalpos[0] << " " << calcalpos[1] << " " << calcalpos[2] << endl; 
  //cout << "calpos2: " << calpos2[0] << " " << calpos2[1] << " " << calpos2[2] << endl; 
  //cout << "result: " << delpos[0] << " " << delpos[1] << " " << delpos[2] << endl; 
  return;
}
