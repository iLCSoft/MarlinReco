# code to derive photon PFO energy and angular corrections as a fn of energy, position

Daniel Jeans, 2022-2025

calculate corrections to energy and direction of photon-like PFOs.

## input samples

- single photon calibration samples of various fixed energies (lcio format)
- suggested range of photon energies ~1 GeV -> ~300 GeV, with of order 10 distinct energy points
- at least 100k photons per energy point, flat in cos(theta), with no overlay or crossing angle applied

## step 1

- python3 prepareInput.py
- from lcio files, makes root tree with photon PFO/Cluster/MC properties

## step 2

- python3 extractCorrection.py
- derive the corrections
- input is the root file produced by prepareInput.py
- output is the parameters to pass to photonCorrectionProcessor