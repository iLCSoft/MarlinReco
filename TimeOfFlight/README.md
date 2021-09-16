
# TimeOfFlight

This package deals with computing and displaying time of flight parameters
based on the timing information stored in CalorimeterHits.


## Processors

### TOFEstimators

author: N.Weinhold (DESY, internship 2017)
author: F.Gaede, DESY, 2018
author: B.Dudar, DESY, 2021
Compute various estimators for time of flight and add these to the
ReconstrucedParticles as PIDParameters.


## Example steering file

[./xml/steer.xml](./xml/steer.xml) is an example of a steering file that runs TOFEstimators processor. It computes momentum harmonic mean, time-of-flight and track length of the PandoraPFOs and stores this information inside PIDHandler. In the example time of flight calculated using closest Ecal hit assuming perfect time resolution. Then new slcio file "test.slcio" produced with LCIOOutputProcessor processor.
