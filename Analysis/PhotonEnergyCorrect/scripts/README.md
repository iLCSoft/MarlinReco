# code to derive photon PFO energy and angular corrections as a fn of energy, position

Daniel Jeans, 2022-2025

Calculate corrections to energy and direction of photon-like PFOs.
For use by the MarlinReco/Analysis/PhotonEnergyCorrect/photonCorrectionProcessor

## input samples

- reconstructed single photon calibration samples of various fixed energies (lcio format, REC or DST)
- suggested range of photon energies ~1 GeV -> ~300 GeV, with at least 10 distinct energy points, approx equidistant in log(E)
- at least 100k photons per energy point, flat in cos(theta), with no overlay or crossing angle applied

## step 1: prepare root trees from lcio files

prepareInput.py :
- Adjust first part of prepareInput.py which defines the "filenames" dictionary:
  this should have as value an array of filenames (with path) of the input lcio files for a particular photon energy,
  and as key an energy-dependent label.
  The current implementation finds the input samples by matchng filenames in some directory.
- Adjust geometrical parameters (inner radius of ecal etc) if using different design than ILD_l5_v02
- Adjust PFO and MCParticle collection names, if necessary
- run it: python3 prepareInput.py
- from the input lcio files, this code prepares a root tree containing photon PFO/Cluster/MC properties

## step 2: calculate optimal parameters for predetermined correction functions

extractCorrection.py :
- set "inname" to the root file produced in the previous step
- set useClus : to correct PFO energies or Cluster energies
  probably usually should use PFO, unless they have already been corrected and
  you want to use the original uncorrected energies (which can be found in the cluster energies)
- set useMC (derive corrections as fn of truth or reco angles)
- run it: python3 extractCorrection.py
- input is the root file produced by prepareInput.py
- output: parameters to pass to photonCorrectionProcessor,
          and a pdf file showing the various dstributions and fits used
- check the output pdf looks "reasonable"
- paste the xml snippet printed at the end of the job into the steering for photonCorrectionProcessor
