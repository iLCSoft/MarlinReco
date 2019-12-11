#include "photonCorrector.h"

#include <iostream>
#include <math.h>
#include <cassert>
#include <streamlog/streamlog.h>

using std::endl;

float photonCorrector::photonEnergyCorrection( EVENT::ReconstructedParticle* rp ) {
  // first correct for gaps
  float gapCorrEn = gapCompensatedEnergy( rp );
  // then correct for non-linearity
  return energyLinearise( gapCorrEn );
}


float photonCorrector::energyLinearise( float en ) {
  float logen = log10(en);
  float energyCorrectionFactor = _energyLin_const + _energyLin_logen*logen;
  return en/energyCorrectionFactor;
}

float photonCorrector::gapCompensatedEnergy( EVENT::ReconstructedParticle* rp ) {
  float gapcorrectionFactor(1.0);
  float en = rp->getEnergy();

  if ( rp->getType() != 22 ) {  // check that it's a photon-like PFO
    streamlog_out (WARNING) << "gapCompensate designed only for photons! not applying correction" << endl;
  } else {
    float cosTheta = rp->getMomentum()[2]/sqrt( pow(rp->getMomentum()[0],2)+pow(rp->getMomentum()[1],2)+pow(rp->getMomentum()[2],2) );
    float phi = atan2( rp->getMomentum()[1], rp->getMomentum()[0] );
    gapcorrectionFactor*=gapCompensate_theta( en, cosTheta );
    if ( fabs(cosTheta) < _barrelendcap_costhlimit ) { // barrel
      gapcorrectionFactor*=gapCompensate_barrelPhi( en, phi );
    } else { // endcap
      float xAcross = getDistanceAcrossEndcapQuadrant( cosTheta, phi );
      gapcorrectionFactor*=gapCompensate_endcap( xAcross );
    }
  }
  return en*gapcorrectionFactor;
}

float photonCorrector::gapCompensate_barrelPhi( float en, float phi ) {
  // correcion for phi gaps in barrel [new, from 2019/11/20]
  float logen = log10(en);
  float foldphi = getBarrelFoldedPhi(phi);
  float par_gap_pos = _phiBarrelCorr_pos_const + _phiBarrelCorr_pos_logen*logen;
  float width = foldphi < par_gap_pos ? _phiBarrelCorr_width1 : _phiBarrelCorr_width2 ;
  float  fitval = 1.;
  fitval += _phiBarrelCorr_depth * exp ( - pow(foldphi - par_gap_pos, 2)/(2*pow(width,2)) );
  return 1./fitval;
}


float photonCorrector::gapCompensate_theta( float en, float costh ) {
  float logen = log10(en);
  float absCosTh = fabs(costh);
  float par_gaus1_norm = _costhCorr_gaus1_norm_const + _costhCorr_gaus1_norm_logen*logen;
  float par_gaus2_norm = _costhCorr_gaus2_norm_const + _costhCorr_gaus2_norm_logen*logen;

  float  fitval = 1.;
  fitval += par_gaus1_norm*exp( -0.5*pow( (absCosTh-_costhCorr_gaus1_mean)/_costhCorr_gaus1_sigm , 2 ) );
  fitval += par_gaus2_norm*exp( -0.5*pow( (absCosTh-_costhCorr_gaus2_mean)/_costhCorr_gaus2_sigm , 2 ) );
  fitval += _costhCorr_gaus3_norm*exp( -0.5*pow( (absCosTh-_costhCorr_gaus3_mean)/_costhCorr_gaus3_sigm , 2 ) );
  if ( absCosTh>_barrelendcap_costhlimit ) fitval*=_costhCorr_endcap_scale;

  return 1./fitval;
}

float photonCorrector::gapCompensate_endcap( float xAcross ) {

  float fitval = 1.;
  fitval += _endcap_gaus1_norm*exp( -0.5*pow( (xAcross-_endcap_gaus1_mean)/_endcap_gaus1_sigm , 2 ) );
  fitval += _endcap_gaus2_norm*exp( -0.5*pow( (xAcross-_endcap_gaus2_mean)/_endcap_gaus2_sigm , 2 ) );

  return 1./fitval;
}


float photonCorrector::getDistanceAcrossEndcapQuadrant( float costh, float phi ) {
  // this calculates the distance across an endcap quadrant (ie perpendicular to slab direction),
  // from the inner edge of quadrant
  assert( fabs(costh) > _barrelendcap_costhlimit );
  if ( costh<0 ) _assumed_endZ*=-1;
  float endX = _assumed_endZ*sin(acos(costh))* cos(phi);
  float endY = _assumed_endZ*sin(acos(costh))* sin(phi);
  // which module [quadrant] is the photon pointing at?
  int quad = -1;
  if ( costh>0 ) {
    if      ( endX > -_assumed_boxsize && endY >  _assumed_boxsize ) quad=0;
    else if ( endX >  _assumed_boxsize && endY <  _assumed_boxsize ) quad=1;
    else if ( endX <  _assumed_boxsize && endY < -_assumed_boxsize ) quad=2;
    else if ( endX < -_assumed_boxsize && endY > -_assumed_boxsize ) quad=3;
  } else {
    if      ( endX <  _assumed_boxsize && endY >  _assumed_boxsize ) quad=0;
    else if ( endX >  _assumed_boxsize && endY > -_assumed_boxsize ) quad=1;
    else if ( endX > -_assumed_boxsize && endY < -_assumed_boxsize ) quad=2;
    else if ( endX < -_assumed_boxsize && endY <  _assumed_boxsize ) quad=3;
  }
  //  float foldX(0);
  float foldY(0);
  if ( quad>=0 ) { // not center box
    float foldPhi = phi + quad*acos(-1)/2.;
    // foldX = fabs(_assumed_endZ)*sin(acos(costh))* cos(foldPhi); // this should flip the -ve side endcap...
    foldY = _assumed_endZ*sin(acos(costh))* sin(foldPhi);
  }
  return foldY;
}

float photonCorrector::getBarrelFoldedPhi( float phi) {
  if ( phi<0 ) phi+=2*acos(-1);
  int sector = phi/(acos(-1)/4);
  float foldedPhi = phi - sector*acos(-1)/4.;
  return foldedPhi;
}

void photonCorrector::printParams() {

  streamlog_out (MESSAGE) << "photonCorrector::printParams" << endl;
  streamlog_out (MESSAGE) << "barrelendcap_limit          " << get_barrelendcap_limit         () << endl;
  streamlog_out (MESSAGE) << "energyLin_const             " << get_energyLin_const            () << endl;
  streamlog_out (MESSAGE) << "energyLin_logen             " << get_energyLin_logen            () << endl;
  streamlog_out (MESSAGE) << "phiBarrelCorr_pos_const     " << get_phiBarrelCorr_pos_const    () << endl;
  streamlog_out (MESSAGE) << "phiBarrelCorr_pos_logen     " << get_phiBarrelCorr_pos_logen    () << endl;
  streamlog_out (MESSAGE) << "phiBarrelCorr_depth         " << get_phiBarrelCorr_depth        () << endl;
  streamlog_out (MESSAGE) << "phiBarrelCorr_width1        " << get_phiBarrelCorr_width1       () << endl;
  streamlog_out (MESSAGE) << "phiBarrelCorr_width2        " << get_phiBarrelCorr_width2       () << endl;
  streamlog_out (MESSAGE) << "costhCorr_gaus1_norm_const  " << get_costhCorr_gaus1_norm_const () << endl;
  streamlog_out (MESSAGE) << "costhCorr_gaus1_norm_logen  " << get_costhCorr_gaus1_norm_logen () << endl;
  streamlog_out (MESSAGE) << "costhCorr_gaus1_mean        " << get_costhCorr_gaus1_mean       () << endl;
  streamlog_out (MESSAGE) << "costhCorr_gaus1_sigm        " << get_costhCorr_gaus1_sigm       () << endl;
  streamlog_out (MESSAGE) << "costhCorr_gaus2_norm_const  " << get_costhCorr_gaus2_norm_const () << endl;
  streamlog_out (MESSAGE) << "costhCorr_gaus2_norm_logen  " << get_costhCorr_gaus2_norm_logen () << endl;
  streamlog_out (MESSAGE) << "costhCorr_gaus2_mean        " << get_costhCorr_gaus2_mean       () << endl;
  streamlog_out (MESSAGE) << "costhCorr_gaus2_sigm        " << get_costhCorr_gaus2_sigm       () << endl;
  streamlog_out (MESSAGE) << "costhCorr_gaus3_norm        " << get_costhCorr_gaus3_norm       () << endl;
  streamlog_out (MESSAGE) << "costhCorr_gaus3_mean        " << get_costhCorr_gaus3_mean       () << endl;
  streamlog_out (MESSAGE) << "costhCorr_gaus3_sigm        " << get_costhCorr_gaus3_sigm       () << endl;
  streamlog_out (MESSAGE) << "costhCorr_endcap_scale      " << get_costhCorr_endcap_scale     () << endl;
  streamlog_out (MESSAGE) << "endcap_gaus1_norm           " << get_endcap_gaus1_norm          () << endl;
  streamlog_out (MESSAGE) << "endcap_gaus1_mean           " << get_endcap_gaus1_mean          () << endl;
  streamlog_out (MESSAGE) << "endcap_gaus1_sigm           " << get_endcap_gaus1_sigm          () << endl;
  streamlog_out (MESSAGE) << "endcap_gaus2_norm           " << get_endcap_gaus2_norm          () << endl;
  streamlog_out (MESSAGE) << "endcap_gaus2_mean           " << get_endcap_gaus2_mean          () << endl;
  streamlog_out (MESSAGE) << "endcap_gaus2_sigm           " << get_endcap_gaus2_sigm          () << endl;
  streamlog_out (MESSAGE) << "assumed_boxsize             " << get_assumed_boxsize            () << endl;
  streamlog_out (MESSAGE) << "assumed_endZ                " << get_assumed_endZ               () << endl;

  return;

}
