## Description

TrackLength is the Marlin processor that calculates track length and square root of the harmonic mean of the squared momentum of charged particles.
Previously has been a part of TOFEstimators processor.

The detailed description of the source code is available in the [doxygen documentation](https://www.desy.de/~dudarboh/track_length_doc/html/index.html). <br>

If formulas are barely visible please switch to the light theme of the github or check the main page in the [doxygen documentation](https://www.desy.de/~dudarboh/track_length_doc/html/index.html).

## Requirements and limitations


TrackLength processor works only for PFOs with **exactly one Track and exactly one Cluster**. In other cases we write all output as 0.<br>

TrackLength processor **works with REC files but does not work with DST** or mini-DST files as it requires individual hit information from the tracker.<br>

You need to **run InitDD4hep processor before** running TrackLength processor. It depends on detector geometry to get magnetic field information.<br>

## Steering parameters

This processor has one steering parameters:

+ ReconstructedParticleCollection | default: "PandoraPFOs"<br>
  Collection of `ReconstructedParticle` objects to analyze.<br>
  The final results will be written in the PIDHandler inside this collection.

## Output parameters

This processor has four output parameters:

+ "momentumHMToSET" / "momentumHMToEcal"  (GeV/c)<br>
  A squared root of the harmonic mean of the squared momentum: <img src="https://render.githubusercontent.com/render/math?math=\sqrt{\langle p^{2} \rangle_{HM}}= \sqrt{ \sum_{i=0}^{n} \ell_{i} \bigg/ \sum_{i=0}^{n} \frac{\ell_{i}}{p_{i}^{2}} }">

  If momentum is constant, then it is just a momentum:  <img src="https://render.githubusercontent.com/render/math?math=p=\sqrt{\langle p^{2} \rangle_{HM}}=const"><br>
  If momentum changes, then using <img src="https://render.githubusercontent.com/render/math?math=\sqrt{\langle p^{2} \rangle_{HM}}"> instead of <img src="https://render.githubusercontent.com/render/math?math=p_{\mathrm{ECal}}"> or <img src="https://render.githubusercontent.com/render/math?math=p_{\mathrm{IP}}"> is mathematically rigorous choice for the relativistic particle assumption <img src="https://render.githubusercontent.com/render/math?math=p \gg m"> <br>
  See [Winfried A. Mitaroff paper](https://arxiv.org/abs/2107.02031) for the details.

  Sum over i up to n includes track length segments between all tracker hits for "momentumHMToSET" and includes one additional segment to the extrapolated track position at the Ecal surface for "momentumHMToEcal".

+ "trackLengthToSET" / "trackLengthToEcal" (mm)<br>
  A track length calculated between the point of the closest approach (PCA) to the IP `(0,0,0)`<br>
  and the last tracker hit for "trackLengthToSET" and ECal surface for "trackLengthToEcal".

  A total track length is calculated as a sum of track segments <img src="https://render.githubusercontent.com/render/math?math=\ell = \sum \ell_{i}"> obtained from iterating over track states defined at the IP, every tracker hit and extrapolated track position at the ECal for "trackLengthToEcal". <br>

  For every pair of track states we check an approximate number of revolutions of the helix with  <img src="https://render.githubusercontent.com/render/math?math=N_{turns} = \frac{\left |z_{i %2B 1} - z_{i}\right |}{|\tan{\lambda}|} \bigg / (\frac{2 \pi}{|\Omega|})  ">.

  In case <img src="https://render.githubusercontent.com/render/math?math=N_{\mathrm{turns}} < 0.5"> the track segment length is calculated between neighboring track states as:
  <img src="https://render.githubusercontent.com/render/math?math=\ell_{i} = \sqrt{\left( \frac{\varphi_{i %2B 1} - \varphi_{i}}{\Omega}\right)^2 %2B \left( z_{i %2B 1} - z_{i} \right)^2 }">.<br>
  In case <img src="https://render.githubusercontent.com/render/math?math=N_{\mathrm{turns}} > 0.5"> the track segment length is calculated between neighboring track states as:
   <img src="https://render.githubusercontent.com/render/math?math=\ell_{\mathrm{last}} = \frac{\left |z_{i %2B 1} - z_{i}\right |}{| \tan{\lambda} |} \sqrt{1 %2B \tan^2{\lambda} }">.<br>

  If helix has more than half revolution between two tracker states then we are unable to use the formula with phi angle, thus we use different formula that does not rely on the azimuthal angle but only on z coordinate and the dip of the helix. Later formula is less precise, so we use it only in these exception cases. More than half turn situation generally should not happen between TPC hits as neighboring TPC hits usually close to each other. But more than half turn often can happen between last TPC hit and the extrapolated track state at the ECal endcap.


## Steering file example

[./xml/steer.xml](./xml/steer.xml) is a steering file example that runs three TrackLength processor. It writes output parameters in the PIDHandlers of the **PandoraPFOs** collection which is then saved in the new **output.slcio** file using LCIOOutputProcessor.

To run this example one needs to setup iLCSoft environment.
If the reader has a NAF account he/she can setup iLCSoft environment with:<br>

    source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh

and run the example steering file with:<br>

    Marlin ./xml/steer.xml

Then one can look at the output.sclio file with, e.g.:<br>

    dumpevent output.slcio 1 | less

One can find output for PandoraPFOs which has new TrackLength algorithm attached...


    collection name : PandoraPFOs
    parameters:
    --------------- print out of ReconstructedParticle collection ---------------
                                        . . .
    parameter ParameterNames_MyTrackLengthProcessor [string]: trackLengthToSET, trackLengthToEcal, momentumHMToSET, momentumHMToEcal,
                                        . . .


... and see final results for each individual PFO


    ------------ detailed PID info: ---
      algorithms :
      [id: 9]   MyTrackLengthProcessor - params:  trackLengthToSET trackLengthToEcal momentumHMToSET momentumHMToEcal

      [particle] |  PDG   | likelihood |  type  |  algoId  | parameters :
                 |        |            |        |          |
      [00000071]                        . . .
                 |      0 | 0.0000e+00 | 000000 |        9 | [ trackLengthToSET : 2.62e+03, trackLengthToEcal : 2.85e+03, momentumHMToSET : 3.84e+00, momentumHMToEcal : 3.84e+00,]

## Analysis

After you run **TOFEstimators** and **TrackLength** processors you might want to run you analysis processor to e.g. calculate the mass of particles using time-of-flight information.

Here is the code example how to do that:
```cpp
    float YourAmazingAnalysisProcessor::getParameterFromPID(ReconstructedParticle* pfo, PIDHandler& pidHandler, std::string algorithmName, std::string parameterName){
        int algorithmID = pidHandler.getAlgorithmID(algorithmName);
        const ParticleID& pfoPID = pidHandler.getParticleID(pfo, algorithmID);
        const std::vector<float>& parameters = pfoPID.getParameters();
        int parIdx = pidHandler.getParameterIndex(algorithmID, parameterName);
        return parameters[parIdx];
    }

    void YourAmazingAnalysisProcessor::processEvent(LCEvent* event){
        LCCollection* pfos = event->getCollection("PandoraPFOs");
        PIDHandler pidHandler(pfos);

        for(int i=0; i < pfos->getNumberOfElements(); ++i){
            ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );

            float momentum = getParameterFromPID(pfo, pidHandler, "MyTrackLengthProcessor", "momentumHMToEcal"); // in GeV
            float trackLength = getParameterFromPID(pfo, pidHandler, "MyTrackLengthProcessor", "trackLengthToEcal"); // in mm
            float tof = getParameterFromPID(pfo, pidHandler, "MyTofClosest0ps", "timeOfFlight"); // in ns

            //calculate mass in GeV using relativistic momentum formula
            double mass = momentum * std::sqrt( std::pow(tof*CLHEP::c_light/trackLength, 2) - 1 );
        }
    }
```
## Authors
- B.Dudar, DESY, 2022<br>
