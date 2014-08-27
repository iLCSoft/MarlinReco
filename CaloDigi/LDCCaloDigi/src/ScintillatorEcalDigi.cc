#include "ScintillatorEcalDigi.h"
#include <assert.h>
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGauss.h"
#include <iostream>
using std::cout;
using std::endl;

// this applies some extra digitisation to scintillator+PPD hits (PPD=SiPM, MPPC)
// - poisson fluctuates the number of photo-electrons according to #PEs/MIP
// - applies PPD saturation according to #pixels
// Daniel Jeans, Jan/Feb 2014. 
// (split off from ILDCaloDigi Aug'14.)

ScintillatorEcalDigi::ScintillatorEcalDigi() {
  _pe_per_mip=-99;
  _calib_mip=-99;
  _npix=-99;
  _misCalibNpix=0;
  _pixSpread=0;
  _elecNoise=0;
}

float ScintillatorEcalDigi::getDigitisedEnergy(float energy) {

  float correctedEnergy(energy);

  if (_pe_per_mip<=0 || _calib_mip<=0 || _npix<=0) {
    cout << "ERROR, crazy parameters for ScintillatorEcalDigi: PE/MIP=" << _pe_per_mip << ", MIP calib=" << _calib_mip << ", #pixels=" << _npix << endl;
    cout << "you must specify at least the #PE/MIP, MIP calibration, and #pixels for realistic ECAL scintillator digitisation!!" << endl;
    cout << "refusing to proceed!" << endl;
    assert(0);
  }

  // 1. convert energy to expected # photoelectrons (via MIPs)
  float npe = _pe_per_mip*energy/_calib_mip;

  // 2. smear according to poisson (PE statistics)
  npe = CLHEP::RandPoisson::shoot( npe );

  if (_npix>0) {
    // 3. apply average sipm saturation behaviour
    npe = _npix*(1.0 - exp( -npe/_npix ) );

    // 4. apply intrinsic sipm fluctuations (see e.g arXiv:0706.0746 [physics.ins-det])
    float alpha = npe/_npix; // frac of hit pixels
    float width = sqrt( _npix*exp(-alpha)*( 1. - (1.+alpha)*exp(-alpha) ) );

    // make an integer correction
    int dpix = int( floor( CLHEP::RandGauss::shoot(0, width) + 0.5 ) );

    npe += dpix;
  }

  if (_pixSpread>0) {
    // variations in pixel capacitance
    npe *= CLHEP::RandGauss::shoot(1, _pixSpread/sqrt(npe) );
  }

  if (_elecNoise>0) {
    // add electronics noise
    npe += CLHEP::RandGauss::shoot(0, _elecNoise*_pe_per_mip);
  }

  if (_npix>0) {
    // 4. unfold the saturation
    // - miscalibration of npix
    float smearedNpix = _misCalibNpix>0 ? _npix*CLHEP::RandGauss::shoot( 1.0, _misCalibNpix ) : _npix;
    // - this is to deal with case when #pe is larger than #pixels (would mean taking log of negative number)
    float epsilon=1; // any small number...
    if ( npe>smearedNpix-epsilon ) npe=smearedNpix-epsilon;
    // - unfold saturation
    npe = -smearedNpix * std::log ( 1. - ( npe / smearedNpix ) );
  }

  // convert back to energy
  correctedEnergy = _calib_mip*npe/_pe_per_mip;

  return correctedEnergy;
}

