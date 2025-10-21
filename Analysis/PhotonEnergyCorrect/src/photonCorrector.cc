#include "photonCorrector.h"

#include <cassert>
#include <iostream>
#include <math.h>
#include <streamlog/streamlog.h>

using std::endl;

void photonCorrector::set_energyCorr_linearise(std::vector<float> pars) {
  if (pars.size() != 2) {
    streamlog_out(ERROR) << "Wrong number of parameters for energyCorr_linearise ! expecting 2 got " << pars.size()
                         << std::endl;
    assert(0);
  } else {
    _energyLin_const = pars[0];
    _energyLin_logen = pars[1];
  }
  return;
}
void photonCorrector::set_energyCorr_barrelPhi(std::vector<float> pars) {
  if (pars.size() != 5) {
    streamlog_out(ERROR) << "Wrong number of parameters for energyCorr_barrelPhi ! expecting 5 got " << pars.size()
                         << std::endl;
    assert(0);
  } else {
    _phiBarrelCorr_pos_const = pars[0];
    _phiBarrelCorr_pos_logen = pars[1];
    _phiBarrelCorr_depth = pars[2];
    _phiBarrelCorr_width1 = pars[3];
    _phiBarrelCorr_width2 = pars[4];
  }
}
void photonCorrector::set_energyCorr_costheta(std::vector<float> pars) {
  if (pars.size() != 12) {
    streamlog_out(ERROR) << "Wrong number of parameters for energyCorr_costheta ! expecting 12 got " << pars.size()
                         << std::endl;
    assert(0);
  } else {
    _costhCorr_gaus1_norm_const = pars[0];
    _costhCorr_gaus1_norm_logen = pars[1];
    _costhCorr_gaus1_mean = pars[2];
    _costhCorr_gaus1_sigm = pars[3];
    _costhCorr_gaus2_norm_const = pars[4];
    _costhCorr_gaus2_norm_logen = pars[5];
    _costhCorr_gaus2_mean = pars[6];
    _costhCorr_gaus2_sigm = pars[7];
    _costhCorr_gaus3_norm = pars[8];
    _costhCorr_gaus3_mean = pars[9];
    _costhCorr_gaus3_sigm = pars[10];
    _costhCorr_endcap_scale = pars[11];
  }
}
void photonCorrector::set_energyCorr_endcap(std::vector<float> pars) {
  if (pars.size() != 6) {
    streamlog_out(ERROR) << "Wrong number of parameters for energyCorr_endcap ! expecting 6, got " << pars.size()
                         << std::endl;
    assert(0);
  } else {
    _endcap_gaus1_norm = pars[0];
    _endcap_gaus1_mean = pars[1];
    _endcap_gaus1_sigm = pars[2];
    _endcap_gaus2_norm = pars[3];
    _endcap_gaus2_mean = pars[4];
    _endcap_gaus2_sigm = pars[5];
  }
}
void photonCorrector::set_phiCorr_barrel(std::vector<float> pars) {
  if (pars.size() != 18) {
    streamlog_out(ERROR) << "Wrong number of parameters for phiCorr_barrel ! expecting 18 got " << pars.size()
                         << std::endl;
    assert(0);
  } else {
    _phiBias_barrel_p0_1 = pars[0];
    _phiBias_barrel_p0_2 = pars[1];
    _phiBias_barrel_p0_3 = pars[2];
    _phiBias_barrel_p0_4 = pars[3];
    _phiBias_barrel_p1_1 = pars[4];
    _phiBias_barrel_p1_2 = pars[5];
    _phiBias_barrel_p2_1 = pars[6];
    _phiBias_barrel_p2_2 = pars[7];
    _phiBias_barrel_p3_1 = pars[8];
    _phiBias_barrel_p4_1 = pars[9];
    _phiBias_barrel_p5_1 = pars[10];
    _phiBias_barrel_p5_2 = pars[11];
    _phiBias_barrel_p5_3 = pars[12];
    _phiBias_barrel_p6_1 = pars[13];
    _phiBias_barrel_p6_2 = pars[14];
    _phiBias_barrel_p7_1 = pars[15];
    _phiBias_barrel_p7_2 = pars[16];
    _phiBias_barrel_p7_3 = pars[17];
  }
}

void photonCorrector::set_thetaCorr_barrel(std::vector<float> pars) {
  if (pars.size() != 4) {
    streamlog_out(ERROR) << "Wrong number of parameters for thetaCorr_barrel ! expecting 4 got " << pars.size()
                         << std::endl;
    assert(0);
  } else {
    _thetaBias_barrel_p0_1 = pars[0];
    _thetaBias_barrel_p0_2 = pars[1];
    _thetaBias_barrel_p1_1 = pars[2];
    _thetaBias_barrel_p1_2 = pars[3];
  }
}

void photonCorrector::set_thetaCorr_endcap(std::vector<float> pars) {
  if (pars.size() != 6) {
    streamlog_out(ERROR) << "Wrong number of parameters for thetaCorr_endcap ! expecting 6 got " << pars.size()
                         << std::endl;
    assert(0);
  } else {
    _thetaBias_endcap_p0_1 = pars[0];
    _thetaBias_endcap_p0_2 = pars[1];
    _thetaBias_endcap_p1_1 = pars[2];
    _thetaBias_endcap_p1_2 = pars[3];
    _thetaBias_endcap_p2_1 = pars[4];
    _thetaBias_endcap_p2_2 = pars[5];
  }
}

float photonCorrector::photonEnergyCorrection(EVENT::ReconstructedParticle* rp) {
  // first correct the energy for gaps
  float gapCorrEn = gapCompensatedEnergy(rp);
  // then correct for overall non-linearity
  return energyLinearise(gapCorrEn);
}

float photonCorrector::energyLinearise(float en) {
  float logen = log10(en);
  float energyCorrectionFactor = _energyLin_const + _energyLin_logen * logen;
  return en / energyCorrectionFactor;
}

float photonCorrector::gapCompensatedEnergy(EVENT::ReconstructedParticle* rp) {
  // returns corrected energy for the input PFO
  float gapcorrectionFactor(1.0);
  float en = rp->getEnergy();

  if (rp->getType() != 22) { // check that it's a photon-like PFO
    streamlog_out(WARNING) << "gapCompensate designed only for photons! not applying correction" << endl;
  } else {
    float cosTheta = rp->getMomentum()[2] /
                     sqrt(pow(rp->getMomentum()[0], 2) + pow(rp->getMomentum()[1], 2) + pow(rp->getMomentum()[2], 2));
    float phi = atan2(rp->getMomentum()[1], rp->getMomentum()[0]);
    gapcorrectionFactor *= gapCompensate_theta(en, cosTheta);
    if (fabs(cosTheta) < _barrelendcap_costhlimit) { // barrel
      gapcorrectionFactor *= gapCompensate_barrelPhi(en, phi);
    } else { // endcap
      float xAcross = getDistanceAcrossEndcapQuadrant(cosTheta, phi);
      gapcorrectionFactor *= gapCompensate_endcap(xAcross);
    }
  }
  return en * gapcorrectionFactor;
}

float photonCorrector::gapCompensate_barrelPhi(float en, float phi) {
  // correction for phi gaps in barrel
  float logen = log10(en);
  float foldphi = getBarrelFoldedPhi(phi);
  float par_gap_pos = _phiBarrelCorr_pos_const + _phiBarrelCorr_pos_logen * logen;
  float width = foldphi < par_gap_pos ? _phiBarrelCorr_width1 : _phiBarrelCorr_width2;
  float fitval = 1.;
  fitval += _phiBarrelCorr_depth * exp(-pow(foldphi - par_gap_pos, 2) / (2 * pow(width, 2)));
  return 1. / fitval;
}

float photonCorrector::gapCompensate_theta(float en, float costh) {
  // correct for gaps in theta
  float logen = log10(en);
  float absCosTh = fabs(costh);
  float par_gaus1_norm = _costhCorr_gaus1_norm_const + _costhCorr_gaus1_norm_logen * logen;
  float par_gaus2_norm = _costhCorr_gaus2_norm_const + _costhCorr_gaus2_norm_logen * logen;

  float fitval = 1.;
  fitval += par_gaus1_norm * exp(-0.5 * pow((absCosTh - _costhCorr_gaus1_mean) / _costhCorr_gaus1_sigm, 2));
  fitval += par_gaus2_norm * exp(-0.5 * pow((absCosTh - _costhCorr_gaus2_mean) / _costhCorr_gaus2_sigm, 2));
  fitval += _costhCorr_gaus3_norm * exp(-0.5 * pow((absCosTh - _costhCorr_gaus3_mean) / _costhCorr_gaus3_sigm, 2));
  if (absCosTh > _barrelendcap_costhlimit)
    fitval *= _costhCorr_endcap_scale;

  return 1. / fitval;
}

float photonCorrector::gapCompensate_endcap(float xAcross) {
  // compensate energy for cracks in endcap. xAcross is the local coord across the endcap quadrant
  float fitval = 1.;
  fitval += _endcap_gaus1_norm * exp(-0.5 * pow((xAcross - _endcap_gaus1_mean) / _endcap_gaus1_sigm, 2));
  fitval += _endcap_gaus2_norm * exp(-0.5 * pow((xAcross - _endcap_gaus2_mean) / _endcap_gaus2_sigm, 2));
  return 1. / fitval;
}

float photonCorrector::getDistanceAcrossEndcapQuadrant(float costh, float phi) {
  // this calculates the distance across an endcap quadrant (ie perpendicular to slab direction),
  // from the inner edge of quadrant
  assert(fabs(costh) > _barrelendcap_costhlimit);
  if (costh < 0)
    _assumed_endZ *= -1;
  float endX = _assumed_endZ * sin(acos(costh)) * cos(phi);
  float endY = _assumed_endZ * sin(acos(costh)) * sin(phi);
  // which module [quadrant] is the photon pointing at?
  int quad = -1;
  if (costh > 0) {
    if (endX > -_assumed_boxsize && endY > _assumed_boxsize)
      quad = 0;
    else if (endX > _assumed_boxsize && endY < _assumed_boxsize)
      quad = 1;
    else if (endX < _assumed_boxsize && endY < -_assumed_boxsize)
      quad = 2;
    else if (endX < -_assumed_boxsize && endY > -_assumed_boxsize)
      quad = 3;
  } else {
    if (endX < _assumed_boxsize && endY > _assumed_boxsize)
      quad = 0;
    else if (endX > _assumed_boxsize && endY > -_assumed_boxsize)
      quad = 1;
    else if (endX > -_assumed_boxsize && endY < -_assumed_boxsize)
      quad = 2;
    else if (endX < -_assumed_boxsize && endY < _assumed_boxsize)
      quad = 3;
  }
  float foldY(0);
  if (quad >= 0) { // not in the center box
    float foldPhi = phi + quad * acos(-1) / 2.;
    foldY = _assumed_endZ * sin(acos(costh)) * sin(foldPhi);
  }
  return foldY;
}

float photonCorrector::getBarrelFoldedPhi(float phi) {
  // fold the barrel phi into a single octant
  if (phi < 0)
    phi += 2 * acos(-1);
  int sector = phi / (acos(-1) / 4);
  float foldedPhi = phi - sector * acos(-1) / 4.;
  return foldedPhi;
}

void photonCorrector::photonDirectionCorrection(EVENT::ReconstructedParticle* rp, float& cor_theta, float& cor_phi) {
  // returns corrected direction (theta, phi) of photon PFOs
  float origEn = rp->getEnergy();
  float p(0);
  for (int i = 0; i < 3; i++) {
    p += pow(rp->getMomentum()[i], 2);
  }
  p = sqrt(p);
  float origCosth = rp->getMomentum()[2] / p;
  float origTheta = acos(origCosth);
  float origPhi = atan2(rp->getMomentum()[1], rp->getMomentum()[0]);
  cor_theta = getCorrectedTheta(origEn, origTheta);
  cor_phi = getCorrectedPhi(origEn, origCosth, origPhi);
  return;
}

float photonCorrector::getBarrelCorrectedPhi(float en, float phi) {
  // phi dorection corrention in the barrel
  float logen = log10(en);
  float phifold = getBarrelFoldedPhi(phi);
  float par[8] = {0};

  //    # energy dependence of function parameters
  par[0] =
      _phiBias_barrel_p0_1 + _phiBias_barrel_p0_2 / (1. + exp(_phiBias_barrel_p0_3 * (logen + _phiBias_barrel_p0_4)));
  par[1] = _phiBias_barrel_p1_1 + logen * _phiBias_barrel_p1_2;
  if (par[1] < 0.)
    par[1] = 0.;
  par[2] = _phiBias_barrel_p2_1 + _phiBias_barrel_p2_2 * logen;
  par[3] = _phiBias_barrel_p3_1;
  par[4] = _phiBias_barrel_p4_1;
  par[5] = _phiBias_barrel_p5_1 + logen * _phiBias_barrel_p5_2 + pow(logen, 2) * _phiBias_barrel_p5_3;
  par[6] = _phiBias_barrel_p6_1 + logen * _phiBias_barrel_p6_2;
  par[7] = _phiBias_barrel_p7_1 + logen * _phiBias_barrel_p7_2 + pow(logen, 2) * _phiBias_barrel_p7_3;

  // # and the function
  float y = par[0];                                                 // # constant
  y = y + par[1] * exp(-0.5 * pow((phifold - par[2]) / par[3], 2)); // # gaussian
  y = y + par[4] * sin(4 * phifold);                                // # sinusoidal
  y = y + par[5] * sin(8 * phifold);
  y = y + par[6] * sin(12 * phifold);
  y = y + par[7] * sin(16 * phifold);
  return phi - y;
}

float photonCorrector::getCorrectedPhi(float en, float costh, float phi) {
  float corPhi(phi);
  if (fabs(costh) < _barrelendcap_costhlimit) {
    corPhi = getBarrelCorrectedPhi(en, phi);
  }
  return corPhi;
}

float photonCorrector::getBarrelCorrectedTheta(float en, float theta) {
  // correct theta in the barrel
  float costh = cos(theta);
  float newtheta(theta);
  if (abs(costh) < _barrelendcap_costhlimit) {
    float logen = log10(en);
    float par[2] = {0};
    // energy dependence of function parameters
    par[0] = _thetaBias_barrel_p0_1 + logen * _thetaBias_barrel_p0_2;
    par[1] = _thetaBias_barrel_p1_1 + logen * _thetaBias_barrel_p1_2;
    // and the function
    float y = par[0] * costh + par[1] * (4 * pow(costh, 3) - 3 * costh);
    newtheta = theta - y;
  }
  return newtheta;
}

float photonCorrector::getEndcapCorrectedTheta(float en, float theta) {
  // correct theta in the endcap
  float logen = log10(en);
  float foldTheta = theta;
  if (theta > acos(-1) / 2.) {
    foldTheta = acos(-1.) - theta;
  }
  float par[3] = {0};
  par[0] = _thetaBias_endcap_p0_1 + _thetaBias_endcap_p0_2 * logen;
  par[1] = _thetaBias_endcap_p1_1 + _thetaBias_endcap_p1_2 * logen;
  par[2] = _thetaBias_endcap_p2_1 + _thetaBias_endcap_p2_2 * logen;
  float cor = par[0] + par[1] * foldTheta + par[2] * pow(foldTheta, 2.);
  if (theta > acos(-1) / 2.) {
    cor *= -1;
  }
  return theta - cor;
}

float photonCorrector::getCorrectedTheta(float en, float theta) {
  // returns corrected theta
  float corTheta(theta);
  float abscth = fabs(cos(theta));
  if (abscth < _barrelendcap_costhlimit - 0.05) { // don't correct overlap region for now
    corTheta = getBarrelCorrectedTheta(en, theta);
  } else if (abscth > _barrelendcap_costhlimit && abscth < cos(atan(_assumed_boxsize / _assumed_endZ))) {
    corTheta = getEndcapCorrectedTheta(en, theta);
  }
  return corTheta;
}

void photonCorrector::printParams() {

  streamlog_out(MESSAGE) << "photonCorrector::printParams" << endl;
  streamlog_out(MESSAGE) << "barrelendcap_limit          " << get_barrelendcap_limit() << endl;
  streamlog_out(MESSAGE) << "energyLin_const             " << get_energyLin_const() << endl;
  streamlog_out(MESSAGE) << "energyLin_logen             " << get_energyLin_logen() << endl;
  streamlog_out(MESSAGE) << "phiBarrelCorr_pos_const     " << get_phiBarrelCorr_pos_const() << endl;
  streamlog_out(MESSAGE) << "phiBarrelCorr_pos_logen     " << get_phiBarrelCorr_pos_logen() << endl;
  streamlog_out(MESSAGE) << "phiBarrelCorr_depth         " << get_phiBarrelCorr_depth() << endl;
  streamlog_out(MESSAGE) << "phiBarrelCorr_width1        " << get_phiBarrelCorr_width1() << endl;
  streamlog_out(MESSAGE) << "phiBarrelCorr_width2        " << get_phiBarrelCorr_width2() << endl;
  streamlog_out(MESSAGE) << "costhCorr_gaus1_norm_const  " << get_costhCorr_gaus1_norm_const() << endl;
  streamlog_out(MESSAGE) << "costhCorr_gaus1_norm_logen  " << get_costhCorr_gaus1_norm_logen() << endl;
  streamlog_out(MESSAGE) << "costhCorr_gaus1_mean        " << get_costhCorr_gaus1_mean() << endl;
  streamlog_out(MESSAGE) << "costhCorr_gaus1_sigm        " << get_costhCorr_gaus1_sigm() << endl;
  streamlog_out(MESSAGE) << "costhCorr_gaus2_norm_const  " << get_costhCorr_gaus2_norm_const() << endl;
  streamlog_out(MESSAGE) << "costhCorr_gaus2_norm_logen  " << get_costhCorr_gaus2_norm_logen() << endl;
  streamlog_out(MESSAGE) << "costhCorr_gaus2_mean        " << get_costhCorr_gaus2_mean() << endl;
  streamlog_out(MESSAGE) << "costhCorr_gaus2_sigm        " << get_costhCorr_gaus2_sigm() << endl;
  streamlog_out(MESSAGE) << "costhCorr_gaus3_norm        " << get_costhCorr_gaus3_norm() << endl;
  streamlog_out(MESSAGE) << "costhCorr_gaus3_mean        " << get_costhCorr_gaus3_mean() << endl;
  streamlog_out(MESSAGE) << "costhCorr_gaus3_sigm        " << get_costhCorr_gaus3_sigm() << endl;
  streamlog_out(MESSAGE) << "costhCorr_endcap_scale      " << get_costhCorr_endcap_scale() << endl;
  streamlog_out(MESSAGE) << "endcap_gaus1_norm           " << get_endcap_gaus1_norm() << endl;
  streamlog_out(MESSAGE) << "endcap_gaus1_mean           " << get_endcap_gaus1_mean() << endl;
  streamlog_out(MESSAGE) << "endcap_gaus1_sigm           " << get_endcap_gaus1_sigm() << endl;
  streamlog_out(MESSAGE) << "endcap_gaus2_norm           " << get_endcap_gaus2_norm() << endl;
  streamlog_out(MESSAGE) << "endcap_gaus2_mean           " << get_endcap_gaus2_mean() << endl;
  streamlog_out(MESSAGE) << "endcap_gaus2_sigm           " << get_endcap_gaus2_sigm() << endl;
  streamlog_out(MESSAGE) << "assumed_boxsize             " << get_assumed_boxsize() << endl;
  streamlog_out(MESSAGE) << "assumed_endZ                " << get_assumed_endZ() << endl;

  streamlog_out(MESSAGE) << "phiBias_barrel_p0_1         " << get_phiBias_barrel_p0_1() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p0_2         " << get_phiBias_barrel_p0_2() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p0_3         " << get_phiBias_barrel_p0_3() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p0_4         " << get_phiBias_barrel_p0_4() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p1_1         " << get_phiBias_barrel_p1_1() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p1_2         " << get_phiBias_barrel_p1_2() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p2_1         " << get_phiBias_barrel_p2_1() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p2_2         " << get_phiBias_barrel_p2_2() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p3_1         " << get_phiBias_barrel_p3_1() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p4_1         " << get_phiBias_barrel_p4_1() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p5_1         " << get_phiBias_barrel_p5_1() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p5_2         " << get_phiBias_barrel_p5_2() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p5_3         " << get_phiBias_barrel_p5_3() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p6_1         " << get_phiBias_barrel_p6_1() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p6_2         " << get_phiBias_barrel_p6_2() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p7_1         " << get_phiBias_barrel_p7_1() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p7_2         " << get_phiBias_barrel_p7_2() << endl;
  streamlog_out(MESSAGE) << "phiBias_barrel_p7_3         " << get_phiBias_barrel_p7_3() << endl;

  streamlog_out(MESSAGE) << "thetaBias_barrel_p0_1       " << get_thetaBias_barrel_p0_1() << endl;
  streamlog_out(MESSAGE) << "thetaBias_barrel_p0_2       " << get_thetaBias_barrel_p0_2() << endl;
  streamlog_out(MESSAGE) << "thetaBias_barrel_p1_1       " << get_thetaBias_barrel_p1_1() << endl;
  streamlog_out(MESSAGE) << "thetaBias_barrel_p1_2       " << get_thetaBias_barrel_p1_2() << endl;

  streamlog_out(MESSAGE) << "thetaBias_endcap_p0_1       " << get_thetaBias_endcap_p0_1() << endl;
  streamlog_out(MESSAGE) << "thetaBias_endcap_p0_2       " << get_thetaBias_endcap_p0_2() << endl;
  streamlog_out(MESSAGE) << "thetaBias_endcap_p1_1       " << get_thetaBias_endcap_p1_1() << endl;
  streamlog_out(MESSAGE) << "thetaBias_endcap_p1_2       " << get_thetaBias_endcap_p1_2() << endl;
  streamlog_out(MESSAGE) << "thetaBias_endcap_p2_1       " << get_thetaBias_endcap_p2_1() << endl;
  streamlog_out(MESSAGE) << "thetaBias_endcap_p2_2       " << get_thetaBias_endcap_p2_2() << endl;

  return;
}
