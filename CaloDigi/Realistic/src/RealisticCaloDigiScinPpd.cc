// Calorimeter digitiser for the IDC ECAL and HCAL
// For other detectors/models SimpleCaloDigi should be used

#include "RealisticCaloDigiScinPpd.h"

#include <marlin/Global.h>

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string>

#include "CLHEP/Random/RandBinomial.h"
#include "CLHEP/Random/RandGauss.h"

using namespace std;
// using namespace lcio ;
using namespace marlin;

RealisticCaloDigiScinPpd aRealisticCaloDigiScinPpd;

RealisticCaloDigiScinPpd::RealisticCaloDigiScinPpd() : RealisticCaloDigi::Processor("RealisticCaloDigiScinPpd") {

  _description = "Performs digitization of sim calo hits...";

  registerProcessorParameter("ppd_mipPe", "# Photo-electrons per MIP (scintillator): used to poisson smear #PEs if >0",
                             _PPD_pe_per_mip, (float)10.);

  registerProcessorParameter("ppd_npix", "total number of MPPC/SiPM pixels for implementation of saturation effect",
                             _PPD_n_pixels, (int)10000);

  registerProcessorParameter("ppd_npix_uncert", "fractional uncertainty of effective total number of MPPC/SiPM pixels",
                             _misCalibNpix, float(0.05));

  registerProcessorParameter("ppd_pix_spread", "variation of PPD pixel signal (as a fraction: 0.01=1%)", _pixSpread,
                             float(0.05));
}

float RealisticCaloDigiScinPpd::convertEnergy(float energy,
                                              int inUnit) const { // convert energy from input to output scale (NPE)
  if (inUnit == NPE)
    return energy;
  else if (inUnit == MIP)
    return _PPD_pe_per_mip * energy;
  else if (inUnit == GEVDEP)
    return _PPD_pe_per_mip * energy / _calib_mip;

  throw std::runtime_error("RealisticCaloDigiScinPpd::convertEnergy - unknown unit " + std::to_string(inUnit));
}

float RealisticCaloDigiScinPpd::digitiseDetectorEnergy(float energy) const {
  // input energy in deposited GeV
  // output in npe
  float npe = energy * _PPD_pe_per_mip / _calib_mip; // convert to pe scale

  if (_PPD_n_pixels > 0) {
    // apply average sipm saturation behaviour
    npe = _PPD_n_pixels * (1.0 - exp(-npe / _PPD_n_pixels));
    // apply binomial smearing
    float p = npe / _PPD_n_pixels;                      // fraction of hit pixels on SiPM
    npe = CLHEP::RandBinomial::shoot(_PPD_n_pixels, p); // npe now quantised to integer pixels

    if (_pixSpread > 0) {
      // variations in pixel capacitance
      npe *= CLHEP::RandGauss::shoot(1, _pixSpread / sqrt(npe));
    }
  }

  return npe;
}
