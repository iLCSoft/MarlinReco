## Description

TOFEstimators is the Marlin processor that calculates time-of-flight for charged particles that can be used for the mass reconstruction and particle identification.<br>

The detailed description of the source code is available in the [doxygen documentation](https://www.desy.de/~dudarboh/tof_doc/html/index.html). <br>

If formulas are barely visible please switch to the light theme of the GitHub or check the main page in the [doxygen documentation](https://www.desy.de/~dudarboh/tof_doc/html/index.html).

## Requirements and limitations

Time-of-flight is calculated only for PFOs with **exactly one Track and exactly one Cluster**. In other cases we write time-of-flight as 0.<br>

TOFEstimators **works with REC files but does not work with DST** or mini-DST files as it requires individual hit information from the calorimeter.<br>

You need to **run InitDD4hep processor before** running TOFEstimator. TOFEstimators is dependent on detector geometry to identify SET hits.<br>

## Steering parameters

TOFEstimators has five steering parameters:

+ ReconstructedParticleCollection | default: "PandoraPFOs"<br>
  Collection of `ReconstructedParticle` objects to analyze.<br>
  The final results will be written in the PIDHandler inside this collection.

+ ExtrapolateToEcal | default: true<br>
  If `true` time-of-flight is calculated using ECal hits via one of the methods chosen in **TofMethod** steering parameter.

  If `false` time-of-flight is calculated using the last tracker hit. Which is most likely SET hit for the barrel region.<br>
  If the last tracker hit is the SET hit time-of-flight is set to the average time of two SET strips.<br>
  If the last tracker hit is not the SET hit time-of-flight is set to `0` (endcap case).

+ TofMethod | options: "closest", "frankAvg", "frankFit" | default: "closest"<br>
  If **ExtrapolateToEcal** is set to `false` then this steering parameter is ignored.

  If **ExtrapolateToEcal** is set to `true` then it defines how to use ECal hits to calculate the time-of-flight to the ECal surface.

  - **closest** uses the closest ECal hit to the extrapolated track position at the ECal surface.<br>
  Time-of-flight is defined as: <img src="https://render.githubusercontent.com/render/math?math=\mathrm{TOF} = t_{\mathrm{closest}} - \frac{\left| \vec{r}_{\mathrm{track}} - \vec{r}_{\mathrm{closest}} \right|}{c}"> <br>
  If no ECal hits are found returns `0.0`.

  - **frankAvg** Select the closest ECal hit to the linearly extrapolated track line inside the ECal in each of the first *MaxEcalLayer* ECal layers. Define time-of-flight as an average of their time corrected for the distance to the track position at the ECal surface assuming speed of flight is the speed of light.
  <img src="https://render.githubusercontent.com/render/math?math=\mathrm{TOF} = \frac{1}{\mathrm{MaxEcalLayer}}\sum_{i}^{\mathrm{MaxEcalLayer}} \left( t_{i} - \frac{\left|\vec{r}_{\mathrm{track}} - \vec{r}_{i} \right|}{c} \right)"> <br>
  If no ECal hits are found returns `0.0`.

  - **frankFit** Select the closest ECal hit to the linearly extrapolated track line inside the ECal in each of the first *MaxEcalLayer* ECal layers. Use linear fit of the hits time as a function of the distance to the track position at the ECal surface <img src="https://render.githubusercontent.com/render/math?math=t=f(|\vec{r}_{\mathrm{track}} - \vec{r}_{\mathrm{hit}} |)"> to define time-of-flight as an extrapolated time at the extrapolated track position at the ECal surface <img src="https://render.githubusercontent.com/render/math?math=\mathrm{TOF}=f(|\vec{r}_{\mathrm{track}} - \vec{r}_{\mathrm{hit}} |=0)"> <br>
  If no ECal hits are found returns `0.0`.<br>
  If only one ECal hit is found, which is not enough to perform the linear fit, then returns the same as *closest* and *frankAvg*.

  for more illustrative explanations of the methods above see slides 10, 12, 13 from the following [LCWS2021 talk]((https://indico.cern.ch/event/995633/contributions/4259659/attachments/2209010/3738157/Bohdan_TOF_LCWS2021.pdf)).


+ TimeResolution | default: 0.0 (ps)<br>
  If **ExtrapolateToEcal** is set to `true` it defines times resolution of individual ECal hits in ps.

  If **ExtrapolateToEcal** is set to `false` it defines times resolution of individual SET strips in ps.

+ MaxEcalLayer | default: 10<br>
  Defines number of inner ECal layers to use for the **frankAvg** and **frankFit** algorithms.

  If **frankAvg** and **frankFit** are unused then this steering parameter is ignored.


## Output parameters

This processor has single output parameter - time-of-flight.

+ timeOfFlight (ns)<br>
  Time-of-flight of the particle assuming <img src="https://render.githubusercontent.com/render/math?math=t_{\mathrm{IP}}=0"> . For the detailed description of different time-of-flight estimators see the **TofMethod** bullet point in the Steering parameters section.

## Steering file example

[./xml/steer.xml](./xml/steer.xml) is a steering file example that runs three TOFEstimators processors: *MyTofClosest0ps*, *MyTofSET10ps* and *MyTofFrankAvg50ps*. They write output parameters in the PIDHandlers of the *PandoraPFOs* collection. Results are saved in the new *output.slcio* file using LCIOOutputProcessor.

To run this example one needs to setup iLCSoft environment.
If the reader has a NAF account he/she can setup iLCSoft environment with:<br>
`source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh`

and run the example steering file with:<br>
`Marlin ./xml/steer.xml`

Then one can look at the output.sclio file with, e.g.:<br>
`dumpevent output.slcio 1 | less`

One can find output for PandoraPFOs which has new TOF algorithms attached...


    collection name : PandoraPFOs
    parameters:
    --------------- print out of ReconstructedParticle collection ---------------
    parameter ParameterNames_MyTofClosest0ps [string]: timeOfFlight,
    parameter ParameterNames_MyTofFrankAvg50ps [string]: timeOfFlight,
    parameter ParameterNames_MyTofSET10ps [string]: timeOfFlight,


... and see final results for each individual PFO


    ------------ detailed PID info: ---
      algorithms :                                        
      [id: 9]   MyTofClosest0ps - params:  timeOfFlight
      [id: 11]   MyTofFrankAvg50ps - params:  timeOfFlight
      [id: 10]   MyTofSET10ps - params:  timeOfFlight

      [particle] |  PDG   | likelihood |  type  |  algoId  | parameters :
                 |        |            |        |          |              
      [00000073]                        . . .
                 |      0 | 0.0000e+00 | 000000 |        9 | [ timeOfFlight : 9.58e+00,]
                 |      0 | 0.0000e+00 | 000000 |       10 | [ timeOfFlight : 0.00e+00,]
                 |      0 | 0.0000e+00 | 000000 |       11 | [ timeOfFlight : 9.57e+00,]

## Analysis

After you run **TOFEstimators** and **TrackLength** processors you might want to run you analysis processor to e.g. calculate the mass of particles using time-of-flight information.

Here is the code example how to do that:

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
            float tof = getParameterFromPID(pfo, pidHandler, "MyTofClosest0ps", "TimeOfFlight"); // in ns

            //calculate mass in GeV using relativistic momentum formula
            double mass = momentum * std::sqrt( std::pow(tof*CLHEP::c_light/trackLength, 2) - 1 );
        }
    }


## Authors
- N.Weinhold, DESY, internship 2017<br>
- F.Gaede, DESY, 2018<br>
- B.Dudar, DESY, 2022<br>
