// Calorimeter digitiser for the IDC ECAL and HCAL
// For other detectors/models SimpleCaloDigi should be used
#include "RealisticCaloDigiSilicon.h"

#include <marlin/Global.h>

#include <iostream>
#include <string>
#include <algorithm>
#include <assert.h>

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace lcio ;
using namespace marlin ;

RealisticCaloDigiSilicon aRealisticCaloDigiSilicon ;

RealisticCaloDigiSilicon::RealisticCaloDigiSilicon() : RealisticCaloDigi::Processor("RealisticCaloDigiSilicon") {
  _description = "Digitization of silicon simcalo hits" ;
  registerProcessorParameter("silicon_pairEnergy" ,
                             "energy required to create e-h pair in silicon (in eV)" ,
                             _ehEnergy,
                             (float)3.6);
}

float RealisticCaloDigiSilicon::convertEnergy( float energy, int inUnit ) const { // convert energy from input to output scale (MIP)
  // converts input energy to MIP scale
  if      ( inUnit==MIP ) return energy;
  else if ( inUnit==GEVDEP ) return energy/_calib_mip;

  throw std::runtime_error("RealisticCaloDigiSilicon::convertEnergy - unknown unit " + std::to_string(inUnit));
}


float RealisticCaloDigiSilicon::digitiseDetectorEnergy(float energy) const {
  // applies extra digitisation to silicon hits
  //  input energy in deposited GeV
  //  output is MIP scale
  float smeared_energy(energy);
  if ( _ehEnergy>0 ) {
    // calculate #e-h pairs
    float nehpairs = 1e9*energy/_ehEnergy; // check units of energy! _ehEnergy is in eV, energy in GeV
    // fluctuate it by Poisson (actually an overestimate: Fano factor actually makes it smaller, however even this overstimated effect is tiny for our purposes)
    smeared_energy *= CLHEP::RandPoisson::shoot( nehpairs )/nehpairs;
  }

  return smeared_energy/_calib_mip; // convert to MIP units
}


