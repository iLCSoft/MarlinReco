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

#include "SiEnergyFluct.h"
#include "PhysicalConstants.h"
#include <algorithm>
#include <cmath>

// Include CLHEP classes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoisson.h"

// Include LCIO classes
#include <EVENT/MCParticle.h>

// Namespaces
using namespace CLHEP;

namespace sistrip {

//
// Constructor setups various constants for dE/dx calculations & universal fluctuations (for Si only)
//
SiEnergyFluct::SiEnergyFluct(double cutOnDeltaRays) : _cutOnDeltaRays(cutOnDeltaRays) {
  // Variables used to reduce recalculation of mean dE/dx
  _prevMCPart = 0;
  _prevMeanLoss = 0.;

  // Cut on secondary electrons
  _cutOnDeltaRays = cutOnDeltaRays;

  // Constants related to dE/dx
  _twoln10 = 2.0 * log(10.0);

  // Constants related to Si material & dE/dx
  _eexc = 173 * eV; // material->GetIonisation()->GetMeanExcitationEnergy();
  _eexc2 = _eexc * _eexc;
  _Zeff = 14;           // material->GetElectronDensity()/material->GetTotNbOfAtomsPerVolume();
  _th = 0.935414 * keV; // 0.25*sqrt(Zeff)

  // Constants related to Si material & dEdx -> density effect
  _aden = 0.160077;              // material->GetIonisation()->GetAdensity();
  _cden = 4.43505;               // material->GetIonisation()->GetCdensity();
  _mden = 3;                     // material->GetIonisation()->GetMdensity();
  _x0den = 0.2;                  // material->GetIonisation()->GetX0density();
  _x1den = 3;                    // material->GetIonisation()->GetX1density();
  _KBethe = 0.178328 * MeV / cm; // material->GetElectronDensity()*twopi_mc2_rcl2;

  // Constants related to dEdx in Si material & muons
  _xgi[0] = 0.0199;
  _xgi[1] = 0.1017;
  _xgi[2] = 0.2372;
  _xgi[3] = 0.4083;
  _xgi[4] = 0.5917;
  _xgi[5] = 0.7628;
  _xgi[6] = 0.8983;
  _xgi[7] = 0.9801;
  _wgi[0] = 0.0506;
  _wgi[1] = 0.1112;
  _wgi[2] = 0.1569;
  _wgi[3] = 0.1813;
  _wgi[4] = 0.1813;
  _wgi[5] = 0.1569;
  _wgi[6] = 0.1112;
  _wgi[7] = 0.0506;

  _limitKinEnergy = 100. * keV;
  _logLimitKinEnergy = log(_limitKinEnergy);
  _alphaPrime = fine_str_const / 2. / pi;

  // Constants related to universal fluctuations
  _minLoss = 10 * eV;
  _minNumberInteractionsBohr = 10.0;
  _nmaxCont1 = 4.;
  _nmaxCont2 = 16.;

  _facwidth = 1. / keV;
  _f1Fluct = 0.857143;       // material->GetIonisation()->GetF1fluct();
  _f2Fluct = 0.142857;       // material->GetIonisation()->GetF2fluct();
  _e1Fluct = 0.115437 * keV; // material->GetIonisation()->GetEnergy1fluct();
  _e2Fluct = 1.96 * keV;     // material->GetIonisation()->GetEnergy2fluct();
  _e1LogFluct = log(_e1Fluct);
  _e2LogFluct = log(_e2Fluct);
  _ipotFluct = 173. * eV; // material->GetIonisation()->GetLogMeanExcEnergy();
  _ipotLogFluct = log(_ipotFluct);
  _e0 = 10 * eV; // material->GetIonisation()->GetEnergy0fluct();
}

//
// Destructor
//
SiEnergyFluct::~SiEnergyFluct() {}

//
// Main method providing energy loss fluctuations, the model used to get the
// fluctuations is essentially the same as in Glandz in Geant3 (Cern program
// library W5013, phys332). L. Urban et al. NIM A362, p.416 (1995) and Geant4
// Physics Reference Manual
//
double SiEnergyFluct::SampleFluctuations(const MCParticle* part, const double length) {
  // Calculate particle related quantities
  double mass = part->getMass() * GeV;
  double momentum2 =
      (part->getMomentum()[0] * part->getMomentum()[0] + part->getMomentum()[1] * part->getMomentum()[1] +
       part->getMomentum()[2] * part->getMomentum()[2]) *
      GeV * GeV;
  double kineticEnergy = sqrt(momentum2 + mass * mass) - mass;
  double ratio = e_mass / mass;
  double tau = kineticEnergy / mass;
  double gam = tau + 1.0;
  double gam2 = gam * gam;
  double bg2 = tau * (tau + 2.0);
  double beta2 = bg2 / (gam2);
  int pdg = part->getPDG();
  double chargeSquare = part->getCharge() * part->getCharge();

  double maxT = 2.0 * e_mass * tau * (tau + 2.) / (1. + 2.0 * (tau + 1.) * ratio + ratio * ratio);
  maxT = std::min(_cutOnDeltaRays, maxT);

  // Recalculate mean loss - mean dE/dx, if necessary
  double meanLoss = 0.;

  if (part != _prevMCPart) {

    // Electron or positron
    if (pdg == 11)
      meanLoss = getElectronDEDX(part);

    // Muon
    else if (pdg == 13)
      meanLoss = getMuonDEDX(part);

    // Hadron
    else
      meanLoss = getHadronDEDX(part);
  }
  // not needed - same particle
  else {

    meanLoss = _prevMeanLoss;
  }

  // Save values for next call
  _prevMCPart = part;
  _prevMeanLoss = meanLoss;

  // Get mean loss in MeV
  meanLoss *= length;

  //
  // Start calculation
  double loss = 0.;
  double sigma = 0.;

  // Trick for very very small loss ( out of model validity)
  if (meanLoss < _minLoss) {
    return meanLoss;
  }

  //
  // Gaussian regime for heavy particles only
  if ((mass > e_mass) && (meanLoss >= _minNumberInteractionsBohr * maxT)) {

    double massrate = e_mass / mass;
    double tmaxkine = 2. * e_mass * beta2 * gam2 / (1. + massrate * (2. * gam + massrate));

    if (tmaxkine <= 2. * maxT) {

      // Sigma
      sigma = (1.0 / beta2 - 0.5) * _KBethe * maxT * length * chargeSquare;
      sigma = sqrt(sigma);

      double twomeanLoss = meanLoss + meanLoss;

      if (twomeanLoss < sigma) {

        double x;
        do {

          loss = twomeanLoss * RandFlat::shoot();
          x = (loss - meanLoss) / sigma;

        } while (1.0 - 0.5 * x * x < RandFlat::shoot());

      } else {

        do {

          loss = RandGaussQ::shoot(meanLoss, sigma);

        } while (loss < 0. || loss > twomeanLoss);
      }
      return loss;
    }
  }

  //
  // Glandz regime : initialisation

  // Trick for very small step or low-density material
  if (maxT <= _e0)
    return meanLoss;

  double a1 = 0.;
  double a2 = 0.;
  double a3 = 0.;

  // Correction to get better width even using stepmax
  //  if(abs(meanLoss- oldloss) < 1.*eV)
  //    samestep += 1;
  //  else
  //    samestep = 1.;
  //  oldloss = meanLoss;

  double width = 1. + _facwidth * meanLoss;
  if (width > 4.50)
    width = 4.50;
  double e1 = width * _e1Fluct;
  double e2 = width * _e2Fluct;

  // Cut and material dependent rate
  double rate = 1.0;
  if (maxT > _ipotFluct) {

    double w2 = log(2. * e_mass * beta2 * gam2) - beta2;

    if (w2 > _ipotLogFluct && w2 > _e2LogFluct) {

      rate = 0.03 + 0.23 * log(log(maxT / _ipotFluct));
      double C = meanLoss * (1. - rate) / (w2 - _ipotLogFluct);

      a1 = C * _f1Fluct * (w2 - _e1LogFluct) / e1;
      a2 = C * _f2Fluct * (w2 - _e2LogFluct) / e2;
    }
  }

  double w1 = maxT / _e0;
  if (maxT > _e0)
    a3 = rate * meanLoss * (maxT - _e0) / (_e0 * maxT * log(w1));

  // 'Nearly' Gaussian fluctuation if a1>nmaxCont2&&a2>nmaxCont2&&a3>nmaxCont2
  double emean = 0.;
  double sig2e = 0.;
  double sige = 0.;
  double p1 = 0.;
  double p2 = 0.;
  double p3 = 0.;

  // Excitation of type 1
  if (a1 > _nmaxCont2) {

    emean += a1 * e1;
    sig2e += a1 * e1 * e1;
  } else if (a1 > 0.) {

    p1 = double(RandPoisson::shoot(a1));
    loss += p1 * e1;

    if (p1 > 0.)
      loss += (1. - 2. * RandFlat::shoot()) * e1;
  }

  // Excitation of type 2
  if (a2 > _nmaxCont2) {

    emean += a2 * e2;
    sig2e += a2 * e2 * e2;
  } else if (a2 > 0.) {

    p2 = double(RandPoisson::shoot(a2));
    loss += p2 * e2;

    if (p2 > 0.)
      loss += (1. - 2. * RandFlat::shoot()) * e2;
  }

  // Ionisation
  double lossc = 0.;

  if (a3 > 0.) {

    p3 = a3;
    double alfa = 1.;

    if (a3 > _nmaxCont2) {

      alfa = w1 * (_nmaxCont2 + a3) / (w1 * _nmaxCont2 + a3);
      double alfa1 = alfa * log(alfa) / (alfa - 1.);
      double namean = a3 * w1 * (alfa - 1.) / ((w1 - 1.) * alfa);
      emean += namean * _e0 * alfa1;
      sig2e += _e0 * _e0 * namean * (alfa - alfa1 * alfa1);
      p3 = a3 - namean;
    }

    double w2 = alfa * _e0;
    double w = (maxT - w2) / maxT;

    int nb = RandPoisson::shoot(p3);

    if (nb > 0)
      for (int k = 0; k < nb; k++)
        lossc += w2 / (1. - w * RandFlat::shoot());
  }

  if (emean > 0.) {

    sige = sqrt(sig2e);
    loss += std::max(0., RandGaussQ::shoot(emean, sige));
  }

  loss += lossc;

  return loss;
}

//
// Method calculating actual dE/dx for hadrons
//
double SiEnergyFluct::getHadronDEDX(const MCParticle* part) {
  // Calculate particle related quantities
  double mass = part->getMass() * GeV;
  double momentum2 =
      (part->getMomentum()[0] * part->getMomentum()[0] + part->getMomentum()[1] * part->getMomentum()[1] +
       part->getMomentum()[2] * part->getMomentum()[2]) *
      GeV * GeV;
  double kineticEnergy = sqrt(momentum2 + mass * mass) - mass;
  double spin = 0.5;
  double charge2 = part->getCharge() * part->getCharge();
  double ratio = e_mass / mass;
  double tau = kineticEnergy / mass;
  double gam = tau + 1.0;
  double bg2 = tau * (tau + 2.0);
  double beta2 = bg2 / (gam * gam);

  double maxT = 2.0 * e_mass * tau * (tau + 2.) / (1. + 2.0 * (tau + 1.) * ratio + ratio * ratio);
  double cutEnergy = std::min(_cutOnDeltaRays, maxT);

  // Start with dE/dx
  double dedx = log(2.0 * e_mass * bg2 * cutEnergy / _eexc2) - (1.0 + cutEnergy / maxT) * beta2;

  // Spin 0.5 particles
  if (0.5 == spin) {

    double del = 0.5 * cutEnergy / (kineticEnergy + mass);
    dedx += del * del;
  }

  // Density correction
  double x = log(bg2) / _twoln10;

  if (x >= _x0den) {

    dedx -= _twoln10 * x - _cden;
    if (x < _x1den)
      dedx -= _aden * pow((_x1den - x), _mden);
  }

  // Shell correction --> not used
  // dedx -= 2.0*corr->ShellCorrection(p,material,kineticEnergy);

  // Now compute the total ionization loss
  if (dedx < 0.0)
    dedx = 0.0;

  dedx *= _KBethe * charge2 / beta2;

  // High order correction only for hadrons --> not used
  // if(!isIon) dedx += corr->HighOrderCorrections(p,material,kineticEnergy);

  return dedx;
}

//
// Method calculating actual dE/dx for muons
//
double SiEnergyFluct::getMuonDEDX(const MCParticle* part) {
  // Calculate particle related quantities
  double mass = part->getMass() * GeV;
  double momentum2 =
      (part->getMomentum()[0] * part->getMomentum()[0] + part->getMomentum()[1] * part->getMomentum()[1] +
       part->getMomentum()[2] * part->getMomentum()[2]) *
      GeV * GeV;
  double kineticEnergy = sqrt(momentum2 + mass * mass) - mass;
  double totEnergy = kineticEnergy + mass;
  double ratio = e_mass / mass;
  double tau = kineticEnergy / mass;
  double gam = tau + 1.0;
  double bg2 = tau * (tau + 2.0);
  double beta2 = bg2 / (gam * gam);

  double maxT = 2.0 * e_mass * tau * (tau + 2.) / (1. + 2.0 * (tau + 1.) * ratio + ratio * ratio);
  double cutEnergy = std::min(_cutOnDeltaRays, maxT);

  // Start with dE/dx
  double dedx = log(2.0 * e_mass * bg2 * cutEnergy / _eexc2) - (1.0 + cutEnergy / maxT) * beta2;

  double del = 0.5 * cutEnergy / totEnergy;
  dedx += del * del;

  // Density correction
  double x = log(bg2) / _twoln10;

  if (x >= _x0den) {

    dedx -= _twoln10 * x - _cden;
    if (x < _x1den)
      dedx -= _aden * pow((_x1den - x), _mden);
  }

  // Shell correction --> not used
  // dedx -= 2.0*corr->ShellCorrection(p,material,kineticEnergy);

  // Now compute the total ionization loss
  if (dedx < 0.0)
    dedx = 0.0;

  // Use radiative corrections of R. Kokoulin
  if (cutEnergy > _limitKinEnergy) {

    double logtmax = log(cutEnergy);
    double logstep = logtmax - _logLimitKinEnergy;
    double dloss = 0.0;
    double ftot2 = 0.5 / (totEnergy * totEnergy);

    for (int ll = 0; ll < 8; ll++) {
      double ep = exp(_logLimitKinEnergy + _xgi[ll] * logstep);
      double a1 = log(1.0 + 2.0 * ep / e_mass);
      double a3 = log(4.0 * totEnergy * (totEnergy - ep) / mass / mass);
      dloss += _wgi[ll] * (1.0 - beta2 * ep / maxT + ep * ep * ftot2) * a1 * (a3 - a1);
    }
    dedx += dloss * logstep * _alphaPrime;
  }

  dedx *= _KBethe / beta2;

  // High order corrections --> not used
  // dedx += corr->HighOrderCorrections(p,material,kineticEnergy);

  return dedx;
}

//
// Method calculating actual dE/dx for electrons & positrons
//
double SiEnergyFluct::getElectronDEDX(const MCParticle* part) {
  // Calculate particle related quantities
  double mass = part->getMass() * GeV;
  double momentum = (part->getMomentum()[0] * part->getMomentum()[0] + part->getMomentum()[1] * part->getMomentum()[1] +
                     part->getMomentum()[2] * part->getMomentum()[2]) *
                    GeV * GeV;
  double kineticEnergy = sqrt(momentum * momentum + mass * mass) - mass;
  double charge = part->getCharge();
  double tau = kineticEnergy / mass;
  double gam = tau + 1.0;
  double gamma2 = gam * gam;
  double bg2 = tau * (tau + 2.0);
  double beta2 = bg2 / (gamma2);

  // Calculate the dE/dx due to the ionization by Seltzer-Berger formula
  bool lowEnergy = false;
  double tkin = kineticEnergy;

  if (kineticEnergy < _th) {
    tkin = _th;
    lowEnergy = true;
  }
  double lowLimit = 0.2 * keV;

  double maxT = kineticEnergy;

  if (charge < 0.)
    maxT *= 0.5;

  double eexc = _eexc / e_mass;
  double eexc2 = eexc * eexc;

  double d = std::min(_cutOnDeltaRays, maxT) / e_mass;
  double dedx;

  // Electron
  if (charge < 0.) {

    dedx = log(2.0 * (tau + 2.0) / eexc2) - 1.0 - beta2 + log((tau - d) * d) + tau / (tau - d) +
           (0.5 * d * d + (2.0 * tau + 1.) * log(1. - d / tau)) / gamma2;
  }
  // Positron
  else {

    double d2 = d * d * 0.5;
    double d3 = d2 * d / 1.5;
    double d4 = d3 * d * 3.75;
    double y = 1.0 / (1.0 + gam);
    dedx = log(2.0 * (tau + 2.0) / eexc2) + log(tau * d) -
           beta2 * (tau + 2.0 * d - y * (3.0 * d2 + y * (d - d3 + y * (d2 - tau * d3 + d4)))) / tau;
  }

  // Density correction
  double x = log(bg2) / _twoln10;

  if (x >= _x0den) {

    dedx -= _twoln10 * x - _cden;
    if (x < _x1den)
      dedx -= _aden * pow(_x1den - x, _mden);
  }

  // Now compute the total ionization loss
  dedx *= _KBethe / beta2;
  if (dedx < 0.0)
    dedx = 0.0;

  // Lowenergy extrapolation

  if (lowEnergy) {

    if (kineticEnergy >= lowLimit)
      dedx *= sqrt(tkin / kineticEnergy);
    else
      dedx *= sqrt(tkin * kineticEnergy) / lowLimit;
  }

  return dedx;
}

} // namespace sistrip
