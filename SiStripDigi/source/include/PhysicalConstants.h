#ifndef PHYSICALCONSTANTS_H
#define PHYSICALCONSTANTS_H 1

namespace sistrip {

//
// System of units
//
//  Basic units:
//
//    charge        in electrons         [e]
//    distance      in centimeters       [cm]
//    energy        in kiloelectronvolts [keV]
//    mag. filed    in Tesla             [T]
//    temperature   in Kelvin            [K]
//    time          in seconds           [s]
//    voltage       in volts             [V]

// Elementary charge
//   static const float       e = 1.602176462E-19;
   static const float       e = 1.;
   static const float   ePlus = 1.;
   static const float  eMinus =-1.;

// Charge
   static const float       C = 1/1.602176462E-19*e;
   static const float      fC = C / 1.E15;

// Distance
   static const float      cm = 1.;
   static const float       m = cm * 100;
   static const float      mm = cm / 10;
   static const float      um = cm / 1.E4;

// Energy
   static const float      eV = 1.;
   static const float     keV = eV * 1.E3;
   static const float     MeV = eV * 1.E6;
   static const float     GeV = eV * 1.E9;

// Temperature
   static const float       K = 1.;

// Time
   static const float       s = 1.;
   static const float      ms = s / 1.E3;
   static const float      us = s / 1.E6;
   static const float      ns = s / 1.E9;

// Voltage
   static const float       V = 1.;

// Magnetic field
   static const float       T = 1.*V*s/m/m;

//
//  Basic constants:
//
//    Boltzmann constant [eV/K]
//    Energy needed for creation of 1 e-h pair [eV]
//    Pi                 [1]

// Pi
   static const double     pi = 3.14159265358979323846;
   static const double piHalf = pi/2;

// Boltzmann constant in eV/K
   static const float       k = 8.617343 * 1.E-5 * eV/K;

// Energy needed for e-h pair creation
   static const float     Eeh = 3.65  * eV;

// Particle physics
   static const double e_mass = 0.510999 * MeV;
   static const double fine_str_const = 1./137.036;

} // Namespace

#endif // PHYSICALCONSTANTS_H
