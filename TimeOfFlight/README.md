# Description

TOFEstimators is the Marlin processor that calculates track length and time-of-flight for charged particles that can be used for the mass reconstruction and particle identification.<br>
In addition it calculates square root of the harmonic mean of the momentum squared. One can use it as a better substitute for the particle momentum in the mass reconstruction.

[Doxygen documentation](https://www.desy.de/~dudarboh/tof_doc/html/index.html) contains a detailed description of the source code. <br>
In case reader uses dark theme on the github and finds hard to read the formulas he can find this README as the main page there.

## Requirements

TOFEstimators is dependent on InitDD4hep processor to extract detector geometry details. <br>
Make sure to run InitDD4hep before this processor is executed.

In addition, input sclio file must contain a lot of various tracker and calorimeter hit collections.<br>
Thus, input slcio file can be e.g. *REC* format, but **cannot** be *DST* or *mini-DST* format.

## Steering parameters

1. **ReconstructedParticleCollection**  / default: `"PandoraPFOs"`<br>
  Collection of `ReconstructedParticle` objects to analyze.<br>
  The final results will be written in the PIDHandler inside this collection.

2. **ExtrapolateToEcal** / default: `true`<br>
  If `true` the track length is calculated to the extrapolated track position at the ECal surface. Time-of-flight then calculated using ECal hits via one of the methods chosen in *TofMethod* steering parameter (see below).

  If `false` the track length is calculated to the last tracker hit. Time-of-flight then calculated using the last tracker hit. Which is SET hit for the barrel region and *undefined* for the endcap region.<br>
  If the last tracker hit is the SET hit time-of-flight is set to the average time of two SET strips.<br>
  If the last tracker hit is the TPC hit time-of-flight is set to `0`.

3. **TofMethod** / options: `"closest"`, `"frankAvg"`, `"frankFit"` / default: `"closest"`<br>
  If *ExtrapolateToEcal* is set to `true` then it defines how to calculate the time-of-flight at the ECal surface.

  - `"closest"` uses the closest ECal hit to the extrapolated track position at the ECal surface.<br>
  Time-of-flight is defined as: <img src="https://render.githubusercontent.com/render/math?math=\mathrm{TOF} = t_{\mathrm{closest}} - \frac{\left| \vec{r}_{\mathrm{track}} - \vec{r}_{\mathrm{closest}} \right|}{c}"> .

  - `"frankAvg"` Select the closest ECal hit to the linearly extrapolated track line inside the ECal in each of the first *MaxEcalLayer* ECal layers. Define time-of-flight as an average of their time corrected for the distance to the track position at the ECal surface assuming speed of flight is the speed of light.
  <img src="https://render.githubusercontent.com/render/math?math=\mathrm{TOF} = \frac{1}{\mathrm{MaxEcalLayer}}\sum_{i}^{\mathrm{MaxEcalLayer}} \left( t_{i} - \frac{\left|\vec{r}_{\mathrm{track}} - \vec{r}_{i} \right|}{c} \right)"> .

  - `"frankFit"` Select the closest ECal hit to the linearly extrapolated track line inside the ECal in each of the first *MaxEcalLayer* ECal layers. Use linear fit of the hits time as a function of the distance to the track position at the ECal surface <img src="https://render.githubusercontent.com/render/math?math=t=f(|\vec{r}_{\mathrm{track}} - \vec{r}_{\mathrm{hit}} |)"> to define time-of-flight as an extrapolated time at the ECal surface <img src="https://render.githubusercontent.com/render/math?math=\mathrm{TOF}=f(|\vec{r}_{\mathrm{track}} - \vec{r}_{\mathrm{hit}} |=0)"> .

  for more illustrative explanations of the methods above see slides 10, 12, 13 from the following [LCWS2021 talk]((https://indico.cern.ch/event/995633/contributions/4259659/attachments/2209010/3738157/Bohdan_TOF_LCWS2021.pdf)).

  If *ExtrapolateToEcal* is set to `false` then this steering parameter is ignored.


4. **TimeResolution** / default: `0.0`<br>
    If *ExtrapolateToEcal* is set to `true` it defines times resolution of individual ECal hits in ps.

    If *ExtrapolateToEcal* is set to `false` it defines times resolution of individual SET strips in ps.

5. **MaxEcalLayer** / default: `10`<br>
    Defines number of inner ECal layers to use for the *frankAvg* and *frankFit* algorithms.

    If *frankAvg* and *frankFit* are unused then this steering parameter is ignored.


## Output parameters

1. `momentumHM`

  A squared root of the harmonic mean of the squared momentum: <img src="https://render.githubusercontent.com/render/math?math=\sqrt{\langle p^{2} \rangle_{HM}}= \sqrt{ \sum_{i=0}^{n} \ell_{i} \bigg/ \sum_{i=0}^{n} \frac{\ell_{i}}{p_{i}^{2}} }">.

  If momentum is constant, then it is just a momentum:  <img src="https://render.githubusercontent.com/render/math?math=p=\sqrt{\langle p^{2} \rangle_{HM}}=const"> .<br>
  If momentum changes, then using <img src="https://render.githubusercontent.com/render/math?math=\sqrt{\langle p^{2} \rangle_{HM}}"> instead of <img src="https://render.githubusercontent.com/render/math?math=p_{\mathrm{ECal}}"> or <img src="https://render.githubusercontent.com/render/math?math=p_{\mathrm{IP}}"> is mathematically rigorous choice for the relativistic particle assumption <img src="https://render.githubusercontent.com/render/math?math=p \gg m"> . <br>
  See [Winfried A. Mitaroff paper](https://arxiv.org/abs/2107.02031) for the details.

2. `trackLength`

  A track length calculated between the point of the closest approach (PCA) to the IP `(0,0,0)`<br>
  and the ECal surface if *ExtrapolateToEcal* is `true`<br>
  or the last tracker hit if *ExtrapolateToEcal* is `false`.

  A total track length is calculated as a sum of track segments <img src="https://render.githubusercontent.com/render/math?math=\ell = \sum \ell_{i}"> obtained from iterating over track states defined at the IP, every tracker hit and the ECal surface if *ExtrapolateToEcal* is `true`. <br>
  A track segment length is calculated between neighboring track states as:
  <img src="https://render.githubusercontent.com/render/math?math=\ell_{i} = \sqrt{\left( \frac{\varphi_{i %2B 1} - \varphi_{i}}{\Omega}\right)^2 %2B \left( z_{i %2B 1} - z_{i} \right)^2 }"> .

  If *ExtrapolateToEcal* is `true` then we additionally check for the number of helix turns between  the last tracker hit and the extrapolated track position at the ECal surface. If <img src="https://render.githubusercontent.com/render/math?math=N_{\mathrm{turns}} > 0.5"> we are unable to use the formula above to calculate the length of the last segment. Thus we use different formula that does not rely on the azimuthal angle:

  <img src="https://render.githubusercontent.com/render/math?math=\ell_{\mathrm{last}} = \frac{\left |z_{i %2B 1} - z_{i}\right |}{\tan{\lambda}} \sqrt{1 %2B \tan^2{\lambda} }"> .


3. `timeOfFlight`

  Time-of-flight of the particle assuming <img src="https://render.githubusercontent.com/render/math?math=t_{\mathrm{IP}}=0"> . For the detailed description of different time-of-flight estimators see the **TofMethod** bullet point in the Steering parameters section.




## Steering file example

[./xml/steer.xml](./xml/steer.xml) is a steering file example that runs three TOFEstimators processors: *MyTofClosest0ps*, *MyTofSET10ps* and *MyTofFrankAvg50ps* and then writes their output parameters in the PIDHandlers of the *PandoraPFOs* collection in the new *output.slcio* file using LCIOOutputProcessor.

To run this example one needs to setup iLCSoft environment.
If the reader has a NAF account he/she can setup iLCSoft environment with:<br>
`source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-02/init_ilcsoft.sh`

and run the example steering file with:<br>
`Marlin ./xml/steer.xml`

Then one can look at the output.sclio file with, e.g.:<br>
`dumpevent output.slcio 1 | less`

Scrolling carefully one can find these lines which indicate that the new file contains our parameters for all three TOFEstimators:

```
collection name : PandoraPFOs
parameters:

--------------- print out of ReconstructedParticle collection ---------------

 flag:  0x0
. . .
parameter ParameterNames_MyTofClosest0ps [string]: momentumHM, trackLength, timeOfFlight,
parameter ParameterNames_MyTofFrankAvg50ps [string]: momentumHM, trackLength, timeOfFlight,
parameter ParameterNames_MyTofSET10ps [string]: momentumHM, trackLength, timeOfFlight,
. . .
```

Scrolling more down, one can see calculated parameters for each individual particle:<br>
```
------------ detailed PID info: ---

  algorithms :
  [id: 0]   BasicVariablePID - params:  electronLikelihood muonLikelihood pionLikelihood kaonLikelihood protonLikelihood hadronLikelihood M
  [id: 3]   LikelihoodPID - params:  electronLikelihood muonLikelihood pionLikelihood kaonLikelihood protonLikelihood hadronLikelihood MVAO
  [id: 4]   LowMomMuID - params:  electronLikelihood muonLikelihood pionLikelihood kaonLikelihood protonLikelihood hadronLikelihood MVAOutp
  [id: 9]   MyTofClosest0ps - params:  momentumHM trackLength timeOfFlight
  [id: 11]   MyTofFrankAvg50ps - params:  momentumHM trackLength timeOfFlight
  [id: 10]   MyTofSET10ps - params:  momentumHM trackLength timeOfFlight
  [id: 2]   ShowerShapesPID - params:  electronLikelihood muonLikelihood pionLikelihood kaonLikelihood protonLikelihood hadronLikelihood MV
  [id: 5]   TOFEstimators0ps - params:  TOFFirstHit TOFClosestHits TOFClosestHitsError TOFFlightLength TOFLastTrkHit TOFLastTrkHitFlightLen
  [id: 8]   TOFEstimators100ps - params:  TOFFirstHit TOFClosestHits TOFClosestHitsError TOFFlightLength TOFLastTrkHit TOFLastTrkHitFlightL
  [id: 6]   TOFEstimators10ps - params:  TOFFirstHit TOFClosestHits TOFClosestHitsError TOFFlightLength TOFLastTrkHit TOFLastTrkHitFlightLe
  [id: 7]   TOFEstimators50ps - params:  TOFFirstHit TOFClosestHits TOFClosestHitsError TOFFlightLength TOFLastTrkHit TOFLastTrkHitFlightLe
  [id: 1]   dEdxPID - params:  electronLikelihood muonLikelihood pionLikelihood kaonLikelihood protonLikelihood hadronLikelihood MVAOutput_


  [particle] |  PDG   | likelihood |  type  |  algoId  | parameters :
             |        |            |        |          |              
. . .
  [00000073] |   2212 | 2.0802e+00 | 000000 |        1 | [ electronLikelihood : -1.64e+02, muonLikelihood : -3.47e+01, pionLikelihood : -2.
             |      0 | 0.0000e+00 | 000000 |        4 | [ electronLikelihood : 9.99e+02, muonLikelihood : 9.99e+02, pionLikelihood : 9.99e
             |      0 | 0.0000e+00 | 000000 |        5 | [ TOFFirstHit : 9.58e+00, TOFClosestHits : 9.58e+00, TOFClosestHitsError : 2.71e-0
             |      0 | 0.0000e+00 | 000000 |        6 | [ TOFFirstHit : 9.57e+00, TOFClosestHits : 9.58e+00, TOFClosestHitsError : 7.16e-0
             |      0 | 0.0000e+00 | 000000 |        7 | [ TOFFirstHit : 9.54e+00, TOFClosestHits : 9.58e+00, TOFClosestHitsError : 6.32e-0
             |      0 | 0.0000e+00 | 000000 |        8 | [ TOFFirstHit : 9.66e+00, TOFClosestHits : 9.64e+00, TOFClosestHitsError : 8.34e-0
             |      0 | 0.0000e+00 | 000000 |        9 | [ momentumHM : 3.84e+00, trackLength : 2.85e+03, timeOfFlight : 9.58e+00,]
             |      0 | 0.0000e+00 | 000000 |       10 | [ momentumHM : 3.84e+00, trackLength : 2.62e+03, timeOfFlight : 0.00e+00,]
             |      0 | 0.0000e+00 | 000000 |       11 | [ momentumHM : 3.84e+00, trackLength : 2.85e+03, timeOfFlight : 9.57e+00,]
. . .
```
## Extracting output parameters

Here is the pseudocode example on how one could extract PIDHandler parameters of *MyTofClosest0ps* results from the slcio file in your own Marlin processor and then use them to e.g. reconstruct the mass:<br>
```
void YourAmazingProcessor::processEvent(LCEvent* event){

    LCCollection* pfos = event->getCollection("PandoraPFOs");
    PIDHandler pidHandler(pfos);

    int algorithmID = pidHandler.getAlgorithmID("MyTofClosest0ps");
    int momParIdx = pidHandler.getParameterIndex(algorithmID, "momentumHM");
    int lenParIdx = pidHandler.getParameterIndex(algorithmID, "trackLength");
    int tofParIdx = pidHandler.getParameterIndex(algorithmID, "timeOfFlight");


    for(int i=0; i < pfos.getNumberOfElements(); ++i){
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );

        const ParticleID& pfoID = pidHandler.getParticleID(pfo, algorithmID);
        const std::vector<float>& pars = pfoID.getParameters();

        double momentum = pars[momParIdx]; // in GeV
        double trackLength = pars[lenParIdx]; // in mm
        double tof = pars[tofParIdx]; // in ns

        //in GeV
        double mass = momentum * std::sqrt( std::pow(tof*CLHEP::c_light/trackLength, 2) - 1 );
    }
}
```


## Authors
- N.Weinhold, DESY, internship 2017<br>
- F.Gaede, DESY, 2018<br>
- B.Dudar, DESY, 2021<br>
