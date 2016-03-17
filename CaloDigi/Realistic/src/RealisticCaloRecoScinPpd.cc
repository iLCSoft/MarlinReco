#include "RealisticCaloRecoScinPpd.h"
#include <cmath>
#include <iostream>

RealisticCaloRecoScinPpd aRealisticCaloRecoScinPpd;

RealisticCaloRecoScinPpd::RealisticCaloRecoScinPpd() : RealisticCaloReco::Processor("RealisticCaloRecoScinPpd") {

  _description = "Performs fist reconstruction of scintillator calo hits";

  registerProcessorParameter("ppd_mipPe" ,
                             "# Photo-electrons per MIP (scintillator): used to poisson smear #PEs if >0" ,
                             _PPD_pe_per_mip,
                             (float)10.);

  registerProcessorParameter("ppd_npix" ,
                             "total number of MPPC/SiPM pixels for implementation of saturation effect" ,
                             _PPD_n_pixels,
                             (int)10000);
}

float RealisticCaloRecoScinPpd::reconstructEnergy(const CalorimeterHit* hit) {
  // here the input energy should be in NPE
  float energy = hit->getEnergy();

  // first de-saturate PPD response
  // this is the fraction of SiPM pixels fired above which a linear continuation of the saturation-reconstruction function is used. 
  // 0.95 of nPixel corresponds to a energy correction of factor ~3.
  const float r = 0.95;
  if (energy < r*_PPD_n_pixels){ //current hit below linearisation threshold, reconstruct energy normally:
    energy = -_PPD_n_pixels * std::log ( 1. - ( energy / _PPD_n_pixels ) );
  } else { //current hit is aove linearisation threshold, reconstruct using linear continuation function:
    energy = 1/(1-r)*(energy-r*_PPD_n_pixels)-_PPD_n_pixels*std::log(1-r);
  }
  // then go back to MIP scale
  energy/=_PPD_pe_per_mip;

  // what layer is this hit in?
  int layer   = (*_idDecoder) (hit)[_cellIDLayerString];
  // now correct for sampling fraction (calibration from MIP -> shower GeV)
  energy *= getLayerCalib( layer );

  return energy;
}


