//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// -------------------------------------------------------------------
//
// This product includes software developed by Members of the Geant4
// Collaboration ( http://cern.ch/geant4 ).
//
// GEANT4 Classes utilized:
//
// File names:    G4UniversalFluctuation (by Vladimir Ivanchenko)
//                G4MollerBhabhaModel (by Vladimir Ivanchenko)
//                G4MuBetheBlochModel (by Vladimir Ivanchenko)
//                G4BetheBlochModel (by Vladimir Ivanchenko)
//
// -------------------------------------------------------------------

#ifndef SIENERGYFLUCT_H
#define SIENERGYFLUCT_H

// Include LCIO classes
#include <lcio.h>
#include <EVENT/MCParticle.h>

using namespace lcio ;

namespace sistrip {

//!
//! Special class providing particle energy loss fluctuations in Si material (Landau fluctuations). The main method
//! simply follows the strategy taken in Geant4 class G4UniversalFluctuation by V. Ivanchenko. As the fluctuation is
//! strongly dependent on particle type and it's energy, detailed calculations of mean ionisation losses have been
//! implemented as well. The differ for hadrons (follows Geant4 class G4BetheBlochModel), muons (follows Geant4
//! class G4MuBetheBlochModel) and electrons & positrons (follows Geant4 class G4MollerBhabhaModel) ... All the details
//! about physics used can be found in http://cern.ch/geant4/UserDocumentation/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
//!
//! @author Z. Drasal, Charles University, Prague
//!
class SiEnergyFluct {

 public:

//!Constructor
   SiEnergyFluct(double cutOnDeltaRays);

//!Destructor
  ~SiEnergyFluct();

//!Method providing energy loss fluctuations, the model used to get the
//!fluctuations is essentially the same as in Glandz in Geant3 (Cern program
//!library W5013, phys332). L. Urban et al. NIM A362, p.416 (1995) and Geant4
//!Physics Reference Manual
   double SampleFluctuations(const MCParticle * part, const double length);


protected:

private:

//!Method calculating actual dEdx for hadrons - based on ComputeDEDXPerVolume method from G4BetheBlochModel Geant4 class
   double getHadronDEDX(const MCParticle * part);

//!Method calculating actual dEdx for muons - based on ComputeDEDXPerVolume method from G4MuBetheBlochModel Geant4 class
   double getMuonDEDX(const MCParticle * part);

//!Method calculating actual dEdx for electrons & positrons - based on ComputeDEDXPerVolume method G4MollerBhabhaModel from Geant4 class
   double getElectronDEDX(const MCParticle * part);

// Pointer to MCParticle given as a parameter during last call of SampleFluctuations method
   const MCParticle * _prevMCPart;

// Mean dE/dx calculated during last call of SampleFluctuations method
   double _prevMeanLoss;

// Cut on secondary electrons
   double _cutOnDeltaRays; //!< Cut on secondary electrons - must be the same as in Geant4

// Constants related to dE/dx
   double _twoln10;

// Constants related to Si material
   double _eexc;
   double _eexc2;
   double _KBethe;
   double _Zeff;
   double _th;

// Constants related to Si material & dEdx -> density effect
   double _aden;
   double _cden;
   double _mden;
   double _x0den;
   double _x1den;

   double _xgi[8];
   double _wgi[8];

   double _limitKinEnergy;
   double _logLimitKinEnergy;
   double _alphaPrime;

// Constants related to Si material & universal fluctuations
   double _minLoss;
   double _minNumberInteractionsBohr;
   double _nmaxCont1;
   double _nmaxCont2;

   double _facwidth;
   double _f1Fluct;
   double _f2Fluct;
   double _e1Fluct;
   double _e2Fluct;
   double _e1LogFluct;
   double _e2LogFluct;
   double _ipotFluct;
   double _ipotLogFluct;
   double _e0;

}; // Class

} // Namespace

#endif // SIENERGYFLUCT_H
