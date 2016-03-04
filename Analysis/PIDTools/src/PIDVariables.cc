#include <PIDVariables.hh>

#include "TRandom3.h"


/*******************************************************
 *
 *   Implementation of PIDVariable_base and its derived classes
 *   Only the Update(...) methods must be redefined
 *
 ******************************************************/

TRandom3* PIDVariable_base::varRand = NULL;

int PIDVariable_base::Update(EVENT::ReconstructedParticle* particle)
{
  EVENT::ClusterVec cluvec=particle->getClusters();
  EVENT::TrackVec trk = particle->getTracks();
  TVector3 p3(particle->getMomentum());

  return Update(cluvec, trk, p3);
}

double PIDVariable_base::BetheBloch(const PIDParticles::PIDParticle_base* hypothesis, const float p) {

  float bg=p/hypothesis->mass;
  float b=sqrt(bg*bg/(1.0+bg*bg));
  //Double_t g=bg/b;
  float tmax=hypothesis->GetBBpars()[2]*TMath::Power(bg,2.0);   ///(1.0+pars[3]*g+pars[4]);

  return (0.5*hypothesis->GetBBpars()[0]*TMath::Log(hypothesis->GetBBpars()[1]*TMath::Power(bg,2.0)*tmax)
          - hypothesis->GetBBpars()[3]*b*b - hypothesis->GetBBpars()[4]*bg/2.0)/(b*b);
}



/***   (ECAL+HCAL)/p   ***/

// 1/10 of the minimum reconstructible pT at ILD
const float PID_CaloTotal::pCut = .01;

PID_CaloTotal::PID_CaloTotal() :
    PIDVariable_base("CaloTotal", "(ECAL+HCAL)/p", "")
{}

int PID_CaloTotal::Update(const EVENT::ClusterVec cluvec, const EVENT::TrackVec trax, const TVector3 p3)
{
  float p = p3.Mag();
  if(p < pCut) {
    SetOutOfRange();
    return MASK_InvalidMomentum;
  }

  float ecal=0., hcal=0.;
  if(cluvec.size()>0){
    for(unsigned int i=0; i<cluvec.size(); i++){
      FloatVec sde = cluvec[i]->getSubdetectorEnergies();
      ecal += sde[0];
      hcal += sde[1];
    }
    _value = (ecal +hcal) / p;
    return 0;
  }
  else {
    _value = 0.;
    return MASK_EmptyClusters;
  }
}


/***   ECAL/(ECAL+HCAL)   ***/

// If total calorimetric deposit smaller than 1 keV consider it as empty clusters
const float PID_CaloEFrac::caloCut = 1.e-6;

PID_CaloEFrac::PID_CaloEFrac() :
    PIDVariable_base("CaloEFrac", "ECAL/(ECAL+HCAL)", "")
{}

int PID_CaloEFrac::Update(const EVENT::ClusterVec cluvec, const EVENT::TrackVec trax, const TVector3 p3)
{
  float ecal=0., hcal=0.;
  if(cluvec.size()>0){
    for(unsigned int i=0; i<cluvec.size(); i++){
      FloatVec sde = cluvec[i]->getSubdetectorEnergies();
      ecal += sde[0];
      hcal += sde[1];
    }

    if ( ecal+hcal < caloCut ) {
      SetOutOfRange();
      return MASK_EmptyClusters;
    }
    else {
      _value = ecal/(ecal+hcal) ;
      return 0;
    }
  }
  else {
    SetOutOfRange();
    return MASK_EmptyClusters;
  }
}


/***   Muon System deposit   ***/

// If muon-system calorimetric deposit smaller than 1 keV consider it as empty clusters
const float PID_CaloMuSys::muSysCut = 1.e-6;

PID_CaloMuSys::PID_CaloMuSys() :
    PIDVariable_base("CaloMuSys", "E_{#mu system}", "GeV")
{}

int PID_CaloMuSys::Update(const EVENT::ClusterVec cluvec, const EVENT::TrackVec trax, const TVector3 p3)
{
  float mucal=0.;
  if(cluvec.size()>0){
    for(unsigned int i=0; i<cluvec.size(); i++){
      FloatVec sde = cluvec[i]->getSubdetectorEnergies();
      mucal += sde[2];
    }
    _value = mucal;
    if(_value < muSysCut) {
      if(varRand) {
        _value += varRand->Gaus(0., 1.e-6);
      }
    }
    return 0;
  }
  else {
    SetOutOfRange();
    return MASK_EmptyClusters;
  }
}


/***   Cluster shapes Masakazu   ***/

PID_CluShapeChi2::PID_CluShapeChi2() :
    PIDVariable_base("CluShapeChi2", "Cluster shape #chi^{2}", "")
{}

int PID_CluShapeChi2::Update(const EVENT::ClusterVec cluvec, const EVENT::TrackVec trax, const TVector3 p3)
{
  if (cluvec.size() < 1) {
    SetOutOfRange();
    return MASK_EmptyClusters;
  }

  FloatVec shapes=cluvec[0]->getShape();
  if(shapes.size()!=0){
    _value = shapes[0];
    if (_value < -.1 ) {
      if(varRand) {
        _value += varRand->Gaus(0., 1.e-6);
      }
    }
    return 0;
  }
  else {
    SetOutOfRange();
    return MASK_EmptyShapes;
  }
}

PID_CluShapeLDiscr::PID_CluShapeLDiscr() :
    PIDVariable_base("DiscrepancyL", "d_{Shower max} / d_{EM shower max}", "")
{}

int PID_CluShapeLDiscr::Update(const EVENT::ClusterVec cluvec, const EVENT::TrackVec trax, const TVector3 p3)
{
  if (cluvec.size() < 1) {
    SetOutOfRange();
    return MASK_EmptyClusters;
  }

  FloatVec shapes=cluvec[0]->getShape();
  if(shapes.size()!=0){
    _value = shapes[5];
    return 0;
  }
  else {
    SetOutOfRange();
    return MASK_EmptyShapes;
  }
}

PID_CluShapeTDiscr::PID_CluShapeTDiscr() :
    PIDVariable_base("DiscrepancyT", "Absorption length", "R_{m}")
{}

int PID_CluShapeTDiscr::Update(const EVENT::ClusterVec cluvec, const EVENT::TrackVec trax, const TVector3 p3)
{
  if (cluvec.size() < 1) {
    SetOutOfRange();
    return MASK_EmptyClusters;
  }

  FloatVec shapes=cluvec[0]->getShape();
  if(shapes.size()!=0){
    _value = shapes[3]/shapes[6];
    return 0;
  }
  else {
    SetOutOfRange();
    return MASK_EmptyShapes;
  }
}

PID_CluShapeXl20::PID_CluShapeXl20() :
    PIDVariable_base("Xl20", "xl20", "?")
{}

int PID_CluShapeXl20::Update(const EVENT::ClusterVec cluvec, const EVENT::TrackVec trax, const TVector3 p3)
{
  if (cluvec.size() < 1) {
    SetOutOfRange();
    return MASK_EmptyClusters;
  }

  FloatVec shapes=cluvec[0]->getShape();
  if(shapes.size()!=0){
    _value = shapes[15]/(2.0*3.50);
    return 0;
  }
  else {
    SetOutOfRange();
    return MASK_EmptyShapes;
  }
}


/***   dE/dx - Chi2 and Log(Chi2)   ***/

PID_dEdxChi2::PID_dEdxChi2(const PID_dEdxChi2& ref) :
    PIDVariable_base(ref.Name(), ref.Description(), ref.Unit()),
    _hypothesis(ref._hypothesis), _dEdx_MIP(ref._dEdx_MIP)
    {}

PID_dEdxChi2::PID_dEdxChi2(const PIDParticles::PIDParticle_base* hypothesis, const float dEdx_MIP) :
    PIDVariable_base(TString::Format("dEdx_chi2_%s", hypothesis->Name()).Data(),
        TString::Format("#chi2_{dE/dx %s}", hypothesis->Name()).Data(), ""),
        _hypothesis(hypothesis), _dEdx_MIP(dEdx_MIP)
{}

PID_dEdxChi2::~PID_dEdxChi2()
{
  delete _hypothesis;
}

int PID_dEdxChi2::Update(const EVENT::ClusterVec cluvec,
    const EVENT::TrackVec trax, const TVector3 p3)
{
  int result = 0;
  float ExpdEdx = BetheBloch(_hypothesis, p3.Mag());
  float dEdx;
  if(trax.size() > 0) { dEdx = trax.at(0)->getdEdx(); }
  else { dEdx = -_dEdx_MIP; result |= MASK_EmptyTracks; } // Keep an eye on get_dEdxChi2()...

  // Normalise dEdx to MIP
  dEdx /= _dEdx_MIP;

  if(dEdx>1.e-15) {
    //get chi2!!(so far 5% error assumed. conservative)
    double normdev = (dEdx-ExpdEdx)/(0.05*dEdx);
    _value = TMath::Sign(TMath::Power(normdev, 2), normdev); // Signed chi2
  }
  else {
    _value = -FLT_MAX;
    result |= MASK_ZerodEdx;
  }

  return result;

}

PID_dEdxLogChi2::PID_dEdxLogChi2(const PID_dEdxLogChi2 &ref) :
    PIDVariable_base(ref.Name(), ref.Description(), ref.Unit()),
    _hypothesis(ref._hypothesis), _dEdx_MIP(ref._dEdx_MIP)
    {}

PID_dEdxLogChi2::PID_dEdxLogChi2(const PIDParticles::PIDParticle_base* hypothesis, const float dEdx_MIP) :
    PIDVariable_base(Form("dEdx_LogChi2_%s", hypothesis->Name()),
        Form("Log(#chi2_{dE/dx %s})", hypothesis->Name()), ""),
        _hypothesis(hypothesis), _dEdx_MIP(dEdx_MIP)
{/*std::cout << _name << "; " << _description << "; " << _unit << std::endl;*/}

PID_dEdxLogChi2::~PID_dEdxLogChi2()
{
  delete _hypothesis;
}

int PID_dEdxLogChi2::Update(const EVENT::ClusterVec cluvec,
    const EVENT::TrackVec trax, const TVector3 p3)
{
  int result = 0;
  float ExpdEdx = BetheBloch(_hypothesis, p3.Mag());
  float dEdx;
  if(trax.size() > 0) { dEdx = trax.at(0)->getdEdx(); }
  else { dEdx = -_dEdx_MIP; result |= MASK_EmptyTracks; } // Keep an eye on get_dEdxChi2()...

  // Normalise dEdx to MIP
  dEdx /= _dEdx_MIP;

  if(dEdx>1.e-15) {
    //get chi2!!(so far 5% error assumed. conservative)
    float normdev = (dEdx-ExpdEdx)/(0.05*dEdx);
    _value = TMath::Sign(float( 2.*TMath::Log(fabs(normdev)+FLT_MIN) ), normdev);
  }
  else {
    _value = -FLT_MAX;
    result |= MASK_ZerodEdx;
  }

  return result;

}


/*******************************************************
 *
 *   Implementation of the PIDVariables base and derived
 *   classes
 *
 ******************************************************/



PIDVariables_base::PIDVariables_base() :
  _p(0.)
{}

PIDVariables_base::PIDVariables_base(EVENT::ReconstructedParticle* particle) :
  _p(0.)
{}

PIDVariables_base::~PIDVariables_base() {
  _varVec.clear();
}


int PIDVariables_base::Update(EVENT::ReconstructedParticle* _particle) {
  EVENT::ClusterVec cluvec=_particle->getClusters();
  EVENT::TrackVec trk = _particle->getTracks();
  TVector3 p3(_particle->getMomentum());

  return Update(cluvec, trk, p3);
}


int PIDVariables_base::Update(const EVENT::ClusterVec cluvec,
    const EVENT::TrackVec trax, const TVector3 p3){

  int result = 0;
  for (VarVec::iterator vit=_varVec.begin(); vit!=_varVec.end(); vit++) {
    result |= (*vit)->Update(cluvec, trax, p3);
  }
  _p = p3.Mag();

  return result;
}


void PIDVariables_base::SetOutOfRange() {
  for (VarVec::iterator vit=_varVec.begin(); vit!=_varVec.end(); vit++)
    { (*vit)->SetOutOfRange(); }
}

void PIDVariables_base::ClearVars() {
/*  for(unsigned int i=0; i<_varVec.size(); i++) {
    delete _varVec.at(i);
  }*/
  _varVec.clear();
}



/***  PIDVariables for the Likelihood PID processor ***/

PIDVariables_LLPID::PIDVariables_LLPID()
{
  Populate();
}


PIDVariables_LLPID::PIDVariables_LLPID(EVENT::ReconstructedParticle* particle)
{
  Populate();
  Update(particle);
}


PIDVariables_LLPID::~PIDVariables_LLPID()
{}

void PIDVariables_LLPID::Populate() {

  _varVec.push_back(new PID_CaloTotal);
  _varVec.push_back(new PID_CaloEFrac);
  _varVec.push_back(new PID_CaloMuSys);

  _varVec.push_back(new PID_CluShapeChi2);
  _varVec.push_back(new PID_CluShapeLDiscr);
  _varVec.push_back(new PID_CluShapeTDiscr);
  _varVec.push_back(new PID_CluShapeXl20);

  _varVec.push_back(new PID_dEdxChi2(&PIDParticles::electronProperties));
  _varVec.push_back(new PID_dEdxChi2(&PIDParticles::muonProperties));
  _varVec.push_back(new PID_dEdxChi2(&PIDParticles::pionProperties));
  _varVec.push_back(new PID_dEdxChi2(&PIDParticles::kaonProperties));
  _varVec.push_back(new PID_dEdxChi2(&PIDParticles::protonProperties));
}


PIDVariables_MvaPid::PIDVariables_MvaPid()
{
  Populate();
}

PIDVariables_MvaPid::PIDVariables_MvaPid(EVENT::ReconstructedParticle* particle)
{
  Populate();
  Update(particle);
}


PIDVariables_MvaPid::~PIDVariables_MvaPid()
{}

int PIDVariables_MvaPid::Update(EVENT::ReconstructedParticle* particle)
{
  int result = PIDVariables_base::Update(particle);
  RefreshMvaVars();
  return result;
}

void PIDVariables_MvaPid::SetOutOfRange()
{
  PIDVariables_base::SetOutOfRange();
  RefreshMvaVars();
}

void PIDVariables_MvaPid::RefreshMvaVars()
{
  for(unsigned int i=0; i<_varVec.size(); i++) {
    _mvaVars.at(i) = _varVec.at(i)->Value();
  }
}

void PIDVariables_MvaPid::Populate() {

  _varVec.push_back(new PID_CaloTotal);
  _varVec.push_back(new PID_CaloEFrac);
  _varVec.push_back(new PID_CaloMuSys);

  _varVec.push_back(new PID_CluShapeChi2);
  _varVec.push_back(new PID_CluShapeLDiscr);
  _varVec.push_back(new PID_CluShapeTDiscr);
  _varVec.push_back(new PID_CluShapeXl20);

  _varVec.push_back(new PID_dEdxLogChi2(&PIDParticles::electronProperties));
  _varVec.push_back(new PID_dEdxLogChi2(&PIDParticles::muonProperties));
  _varVec.push_back(new PID_dEdxLogChi2(&PIDParticles::pionProperties));
  _varVec.push_back(new PID_dEdxLogChi2(&PIDParticles::kaonProperties));
  _varVec.push_back(new PID_dEdxLogChi2(&PIDParticles::protonProperties));

  for(unsigned int i=0; i<_varVec.size(); i++) {
    _mvaVars.push_back(0.);
  }
}
