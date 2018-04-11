
# TimeOfFlight

This package deals with computing and displaying time of flight parameters
based on the timing information stored in CalorimeterHits.


## Processors

### TOFEstimators

author: F.Gaede, DESY, 2018

Compute various estimators for time of flight and add these to the 
ReconstrucedParticles as PIDParameters.


### TOFPlots

author: N.Weinhold (DESY, internship 2017)

Computes various estimators for the time of flight from CalorimeterHits.
Creates ROOT histograms for these parameters.


## Examples

See [./scripts/tofestimators.xml](./scripts/tofestimators.xml) for and example 
steering file that computes TOF estimators for 0,10 and 50 ps single hit 
time resolution.
Use [./scripts/draw_tofestimators.C](./scripts/draw_tofestimators.C) to create 
plots of beta vs. momentum for the diffefrent estimators.

