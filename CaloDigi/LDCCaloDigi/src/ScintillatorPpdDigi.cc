#include "ScintillatorPpdDigi.h"
#include <assert.h>
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandBinomial.h"
#include <iostream>
using std::cout;
using std::endl;

// this applies some extra digitisation to scintillator+PPD hits (PPD=SiPM, MPPC)
// - poisson fluctuates the number of photo-electrons according to #PEs/MIP
// - applies PPD saturation according to #pixels
// Daniel Jeans, Jan/Feb 2014. 
// (split off from ILDCaloDigi Aug'14.)

ScintillatorPpdDigi::ScintillatorPpdDigi() {
  _pe_per_mip=-99;
  _calib_mip=-99;
  _npix=-99;
  _misCalibNpix=0;
  _pixSpread=0;
  _elecNoise=0;
}

void ScintillatorPpdDigi::printParameters() {
  cout << "--------------------------------" << endl;
  cout << "ScintillatorPpdDigi parameters" << endl;
  cout << " pe_per_mip   = " <<  _pe_per_mip << endl;
  cout << " calib_mip    = " <<  _calib_mip << endl;
  cout << " npix         = " <<  _npix << endl;
  cout << " misCalibNpix = " <<  _misCalibNpix << endl;
  cout << " pixSpread    = " <<  _pixSpread << endl;
  cout << " elecNoise    = " <<  _elecNoise << endl;    
  cout << "--------------------------------" << endl;
  return;
}

float ScintillatorPpdDigi::getDigitisedEnergy(float energy) {

  float correctedEnergy(energy);

  if (_pe_per_mip<=0 || _calib_mip<=0 || _npix<=0) {
    cout << "ERROR, crazy parameters for ScintillatorPpdDigi: PE/MIP=" << _pe_per_mip << ", MIP calib=" << _calib_mip << ", #pixels=" << _npix << endl;
    cout << "you must specify at least the #PE/MIP, MIP calibration, and #pixels for realistic scintillator digitisation!!" << endl;
    cout << "refusing to proceed!" << endl;
    assert(0);
  }

  // 1. convert energy to expected # photoelectrons (via MIPs)
  float npe = _pe_per_mip*energy/_calib_mip;

  //oh: commented out Daniel's digitisation model. (npe -> poisson -> saturation -> stoykov smearing).
  // Stoykov smearing used with Gaussian shape for lack of better model.
  /*
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
*/
  
  //AHCAL TB style digitisation: npe -> saturation -> binomial smearing
  //shown to be mathematically equivalent to Daniel's model above, but slightly faster and automatically generates correct shape instead of Gaussian approximation
  
  if (_npix>0){
    // apply average sipm saturation behaviour
    npe = _npix*(1.0 - exp( -npe/_npix ) );
    
    //apply binomial smearing
    float p = npe/_npix; // fraction of hit pixels on SiPM
    npe = CLHEP::RandBinomial::shoot(_npix, p); //npe now quantised to integer pixels
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
    
    //oh: commented out daniel's implmentation of dealing with hits>smearedNpix. using linearisation of saturation-reconstruction for high amplitude hits instead.
    /*
    // - this is to deal with case when #pe is larger than #pixels (would mean taking log of negative number)
    float epsilon=1; // any small number...
    if ( npe>smearedNpix-epsilon ) npe=smearedNpix-epsilon;
    // - unfold saturation
    npe = -smearedNpix * std::log ( 1. - ( npe / smearedNpix ) );
    */
    
    const float r = 0.95; //this is the fraction of SiPM pixels fired above which a linear continuation of the saturation-reconstruction function is used. 0.95 of nPixel corresponds to a energy correction of factor ~3.

    if (npe < r*smearedNpix){ //current hit below linearisation threshold, reconstruct energy normally:
      npe = -smearedNpix * std::log ( 1. - ( npe / smearedNpix ) );
    } else { //current hit is aove linearisation threshold, reconstruct using linear continuation function:
      npe = 1/(1-r)*(npe-r*smearedNpix)-smearedNpix*std::log(1-r);
    }
  }

  // convert back to energy
  correctedEnergy = _calib_mip*npe/_pe_per_mip;

  return correctedEnergy;
}

