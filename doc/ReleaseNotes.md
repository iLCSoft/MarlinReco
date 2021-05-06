# v01-30

* 2021-03-03 Remi Ete ([PR#83](https://github.com/iLCSoft/MarlinReco/pull/83))
  - RealisticCaloDigi: Added new option for energy/charge and time integration
    - deal with slow/fast shaper of ROC chips and time estimate
    - new processor parameters:
       - **integrationMethod**: "Standard" for old implementation (default) or "ROC" emulating the behavior of ROC chip.
       - **fastShaper**: fast shaper time, unit in ns. Only for ROC method
       - **slowShaper**: slow shaper time, unit in ns. Only for ROC method
       - **timingResolution**: optional time resolution gaussian smearing to apply, unit . Only apply is > 0. Default is 0 (no smearing)

# v01-29

* 2021-02-26 Bohdan Dudar ([PR#88](https://github.com/iLCSoft/MarlinReco/pull/88))
  - TOFEstimators processor: adding new algorithms to calculate the ToF, track length and momentum.
    - New output parameters:
      - TOFClosest --- based on the closest ECAL hit
      - TOFFastest --- based on the fastest ECAL hit
      - TOFCylFit --- extrapolated time from the ECAL hits within a cylinder of a shower core
      - TOFClosestFit --- extrapolated time from the ECAL hits closest to the linear continuation of the track inside the ECAL
      - FlightLength --- the helix track length based on the track parameters from the TrackState at the ECAL
      - MomAtCal --- the momentum based on the track parameters from the TrackState at the ECAL
    - New steering parameters:
      - ProcessorVersion --- changes output between the idr version (before this patch) and the dev (after this patch). Default: idr
      - CylRadius --- the radius within which to select hits for the TOFCylFit method. Default: 5 mm

* 2021-02-22 Carl Mikael Berggren ([PR#89](https://github.com/iLCSoft/MarlinReco/pull/89))
  Heavily reworked version of TrueJet
  ===========================
  
  From the outside, little has changed:
  New jet-type (6) for non-isr photon from the hard interaction.
  Does not require recoparticles (nor recomctruthlink), so also works
  on generator-output lcio-files.
  Boolean steering flag to indicate Whizard1 (GDE) or Whizard2 (LCC) input
  (false by default, i.e. input is Whizard2)
  
  Internally, heavily re-worked.
  
  - removed methods fix94 and fix_top
  - New method (stdhep_reader_bug_workaround) called if _whiz1 is true  
  - New method - photon - added to treat the new jet type 6 (M.E. photon).
  - New boolean data-members _whiz1 (steering flag),  _higgs_to_glue_glue
  and _top_event (the latter two set at run-time event by event)
  - Handle the case if no reconstructed information (generator input). 
  - In top events, define the initial colour-neutral as  the top 
  - Default collectionnames changed to the mc2020 ones.
  
  Changes and enhancemts to TrueJet_Parser and Use_TrueJet to go with new TrueJet
  =================================================================
  
  TrueJet_Parser.
  
  - added  MCPseen : LCIntExtension<MCPseen>, used to avoid
  double counting of *true* particles (decay-in-flight ...) for true-of-seen.
  
  - New methods: mcpjet, mcpicn,mcpfcn, recojet, recoicn and recofcn to 
  returns the  corresponding jet, icn or fcn number of MCPs or PFOs.
  true_partics, to return the list of all mcps in a jet (reco_particle
  already existed).
  
  Use_TrueJet:
  
  - Gracefully handle case with no recoMCTruthLink. 
  
  - Exercise the new methods in TrueJet_Parser.

* 2021-02-10 Yasser Radkhorrami ([PR#85](https://github.com/iLCSoft/MarlinReco/pull/85))
  Improved the `ErrorFlow` processor:
  1.  An option added to enable/disable confusion term in the jet energy error
  2. For scaling angular uncertainties of jets, CovMat of the jet is scaled by scaling factors (as processor parameters)
  3. So far, charged PFOs used to be identified by charge, so the energy of neutral PFOs with tracks was added to photons/neutral hadrons energy. Now, neutral PFOs are classified as charged PFO, since the energy and momentum are obtained from tracks.
  4. An option added to propagate confusion term to **all** covariance matrix elements.
  
  - New processor parameters:
        -   "EnableConfusionTerm", "Enable/disable confusion term to be added to covariance matrix"
        -   "CovMatFactorPhotons", "A correction factor to be multiplied to angular uncertainties of photons"
        -   "CovMatFactorNeutralHadrons", "A correction factor to be multiplied to angular uncertainties of Neutral Hadrons"
        -   "PropagateConfusion2Mom", "Enable/disable Propagating uncertainty due to confusion to the Momentum components/All CovMat elements"

* 2021-02-10 Yasser Radkhorrami ([PR#85](https://github.com/iLCSoft/MarlinReco/pull/85))
  ...

* 2021-02-01 tmadlener ([PR#86](https://github.com/iLCSoft/MarlinReco/pull/86))
  - Migrate CI from travis to github actions

# v01-28

* 2020-09-02 Carl Mikael Berggren ([PR#81](https://github.com/iLCSoft/MarlinReco/pull/81))
  - Fix wrong signs in Jacobian in the transformation of the covariance matrix of cluster CoG and Energy to neutral PFO (E,px,py,pz), and double declaration of local variable Eerror.  Update example  CMakeLists.txt and AddClusterProperties.xml to work in the present world.

* 2020-09-02 Junping Tian ([PR#79](https://github.com/iLCSoft/MarlinReco/pull/79))
  - added a new processor for finding isolated photon
  - the default option for isolated muon finder is changed back to use Yoke energy

* 2020-08-31 Yasser Radkhorrami ([PR#80](https://github.com/iLCSoft/MarlinReco/pull/80))
  - An option is added to include full CovMat of neutral PFOs in jet error

# v01-27

* 2020-07-01 Junping Tian ([PR#78](https://github.com/iLCSoft/MarlinReco/pull/78))
  - add a new processor which can be used to obtain the input variables for isolated lepton training
  - add a MVA Classification macro for training

* 2020-07-01 Daniel Jeans ([PR#74](https://github.com/iLCSoft/MarlinReco/pull/74))
  - DDStripSplitter:
    - now run separately for barrel and endcap
    - create hit relations for split hits
    - some cleaning up of code (remove some histograms; fix compiler warnings, etc)
  - RecoMCTruthLinker
    - adapt to work with split hits
    - fix some compiler warnings
  - PhotonCorrectionProcessor
    - check that corrections requested before calculating them
    - add some debug printouts

* 2020-06-29 Junping Tian ([PR#76](https://github.com/iLCSoft/MarlinReco/pull/76))
  - added a processor which can be used to obtain the needed input variables for MVA training of IsolatedLeptonTagging
  - added a root macro of MVA Classification for the training.

* 2020-06-15 JennyListDESY ([PR#75](https://github.com/iLCSoft/MarlinReco/pull/75))
  - adding four-momentum covariance matrix calculation for V0s 
    based on covariance matrices of the two tracks

* 2020-05-13 Daniel Jeans ([PR#73](https://github.com/iLCSoft/MarlinReco/pull/73))
  - add photonPFO direction correction to "PhotonEnergyCorrect" package

* 2020-04-17 Daniel Jeans ([PR#72](https://github.com/iLCSoft/MarlinReco/pull/72))
  - bug fix in photonCorrectionProcessor:
         - adjust the magnitude of the *momentum* of photon PFOs, not just their energies.
            (this bug resulted in photon PFOs with inconsistent energy and momentum.)

# v01-26

* 2019-12-11 Daniel Jeans ([PR#71](https://github.com/iLCSoft/MarlinReco/pull/71))
  - remove hard-coded default correction; control everything from processor parameters.
  - add getter and print functionality to PhotonCorrector class

* 2019-12-03 Daniel Jeans ([PR#70](https://github.com/iLCSoft/MarlinReco/pull/70))
  BruteForceEcalGapFiller makes new hits both within and between ECAL modules to estimate energy lost in cracks. In this update:
  - allow switching off of correction in gaps between modules by processor parameter (in fact, set it off by default)
  - Fix a couple of compiler warnings
  - add photonCorrectionProcessor to correct photon PFO energies

* 2019-10-06 Remi Ete ([PR#69](https://github.com/iLCSoft/MarlinReco/pull/69))
  - Replaced raw usage of streamlog by pre-processor macros from streamlog package

* 2019-09-03 beyerja ([PR#68](https://github.com/iLCSoft/MarlinReco/pull/68))
  - Remove compiler warnings for ErrorFlow processor (default initialization in header, using shared_ptr instead of raw pointer).
  - Add ErrorFlow processor to list of processors to be compiled by default in MarlinReco.

* 2019-08-12 Ete Remi ([PR#67](https://github.com/iLCSoft/MarlinReco/pull/67))
  - SimpleFcalDigi processor:
      - Missing cellid 1 flag in fcal output collection, causing issues in the BeamCal reconstruction

* 2019-04-08 Ete Remi ([PR#66](https://github.com/iLCSoft/MarlinReco/pull/66))
  - ReconstructedParticleImpl_CopyProcessor:
      - Added relation collection to link old PFOs to new PFOs

* 2019-03-29 Matthias Artur Weber ([PR#65](https://github.com/iLCSoft/MarlinReco/pull/65))
  - IsolatedLeptonFinder: save track and cluster information for isolated lepton candidates

* 2019-01-30 Carl Mikael Berggren ([PR#64](https://github.com/iLCSoft/MarlinReco/pull/64))
  - TrueJet: 
    - Much improved assignment of PFOs to jets, and hence of seen jet energies
    - Get rid of all -Wall warnings

* 2019-01-16 Carl Mikael Berggren ([PR#63](https://github.com/iLCSoft/MarlinReco/pull/63))
  - bug fix in TrueJet_Parser
        - fix in TrueJet_Parser::final_cn, replace math.h by cmath in Use_TrueJet

* 2019-01-15 Junping Tian ([PR#62](https://github.com/iLCSoft/MarlinReco/pull/62))
  - added fix to IsolatedLeptonfinder
       - added a new feature for supporting applications not using impact parameters
       - added corresponding weight files trained for the new feature

* 2018-12-19 Ulrich Einhaus ([PR#60](https://github.com/iLCSoft/MarlinReco/pull/60))
  - Add ReconstructedParticleImpl_CopyProcessor

* 2018-11-09 Daniel Jeans ([PR#59](https://github.com/iLCSoft/MarlinReco/pull/59))
  - add new Marlin processor DDStripSplitter, which implements the Strip Splitting Algorithm for eg a scintillator strip ECAL.

* 2018-10-29 Ete Remi ([PR#58](https://github.com/iLCSoft/MarlinReco/pull/58))
  - IsolatedLeptonTaggingProcessor processor:
     - Fixed processor logic when PFO or Vertex input collections don't exist or are empty

* 2018-10-19 Erica Brondolin ([PR#57](https://github.com/iLCSoft/MarlinReco/pull/57))
  - CLICPfoSelectorAnalysis: Fix entries of time vs pt

* 2018-10-02 Frank Gaede ([PR#55](https://github.com/iLCSoft/MarlinReco/pull/55))
  - add sub-package Analysis/GammaGammaHadronRemoval
    - first processor *TrackZVertexGrouping* implementation of track grouping based on Z0 significance
    - algorithm developed by S.Sasikumar, DESY

# v01-25

* 2018-07-05 Jakob Beyer ([PR#52](https://github.com/ilcsoft/MarlinReco/pull/52))
  - Adding new analysis toolTJjetsPFOAnalysisProcessor: 
        - Combined the PFOAnalysis processor with the jet analysis power of the TrueJet/TrueJet_Parser tools to gain insight into individual jet behaviour and reconstruction.

* 2018-08-06 Erica Brondolin ([PR#54](https://github.com/ilcsoft/MarlinReco/pull/54))
  - Introduce CLICPfoSelectorAnalysis which runs on the PFO input collection and creates: a TTree with the PFO variables used in the CLICPfoSelector, cluster time vs pT graphs for each particle category and region, and PFO energy sum histos for each particle category and region
  - CLICPfoSelectorAnalysis has the possibility to detect if the PFO belongs to signal/overlay
  - CLICPfoSelectorAnalysis has the possibility to check if the track and the cluster belonging to the same PFO were produced by at least one common MCParticle

* 2018-07-18 Junping Tian ([PR#53](https://github.com/ilcsoft/MarlinReco/pull/53))
  - fixed IsolatedLeptonTagger
        - fixed the problem about track impact parameters in the new samples where interaction point is     smeared
        - some minor updates about pre-cut values and symmetric treatment for d0/z0 significance
        - new weights trained for new samples are provided

# v01-24-01

* 2018-04-18 Ete Remi ([PR#51](https://github.com/ilcsoft/MarlinReco/pull/51))
  - RecoMCThruthLinker processor
      - Turned WARNING message to DEBUG9 to avoid log pollution
      - this warning occurs only for the SDHcal case where one SimCalorimeterHit can create more than CalorimeterHit

# v01-24

* 2018-04-10 Guillaume ([PR#48](https://github.com/ilcsoft/MarlinReco/pull/48))
  - SDHCAL digitizer : 
     - Switch order of the LCRelation collection between SDHCAL SimCalorimeterHits and Digitized CalorimeterHits (now from CalorimeterHit to SimCalorimeterHit)

* 2018-04-18 Carl Mikael Berggren ([PR#50](https://github.com/ilcsoft/MarlinReco/pull/50))
  - improved TrueJet processor:
       - Fixed crash due to index out-of-range
       - Remove all compiler warnings except local shadow (cheked to be OK)
       -  id of initial ColourNeutrals fixed (should be a boson (W,Z,H))
       - MCParticle collection does not need to start with the beam-particles, back-tracking now also
  gracefully stops if the first entry is reached. This should allow for usage also for the DBD-250 samples, in which the initial beam-particles are missing in the MCParticle collections.
       - nitty-gritty special cases in history fixed.
       - Now also works for higgs-samples *except for h->gluon gluon*, which will need a completely different treatment fro back-track through the parton shower, as there is no quark-line to follow...

* 2018-04-17 Frank Gaede ([PR#49](https://github.com/ilcsoft/MarlinReco/pull/49))
  - add new package TimeOfFlight
       - use TOFEstimators processor to compute TOF parameters
       - will be added as PID object to the ReconstructedParticles

# v01-23

* 2018-01-31 Strahinja Lukic ([PR#37](https://github.com/iLCSoft/MarlinReco/pull/37))
  Updates of SiTracker_dEdxProcessor:
  
  - Cleaned up unnecessary code.
  - Added runtime protection against failure of `TrackImpl::getTrackState()`.
  - Replaced `MarlinTrk::IMarlinTrack::propagate()` which was making the processor ~1000X slower than necessary with `MarlinTrk::IMarlinTrack::extrapolate()`.
  - Removed unnecessary parameters.
  
  Updates of AnalyseSidEdxProcessor:
  
  -   Added minor runtime protections.

* 2018-01-09 Strahinja Lukic ([PR#35](https://github.com/iLCSoft/MarlinReco/pull/35))
  - `SiTracker_dEdxProcessor` was adapted to determine the barrel/endcap type of tracker detector by checking the layering extension, rather than the type flag as before. 
  - A bug was corrected in `SiTracker_dEdxProcessor` that caused miscalculation of total sensor thickness for some of the available dEdx estimators.

* 2018-02-25 KURATA Masakazu ([PR#42](https://github.com/iLCSoft/MarlinReco/pull/42))
  - Focus on low momentum mu/pi separation
  - Correct the corresponding change for createPDF processor
   
  - update: correct bad lines:
    1. use DD4hepAPI to get b field
    2. correct array initialization
    3. prevent memory leak when using ROOT

* 2018-02-25 KURATA Masakazu ([PR#42](https://github.com/iLCSoft/MarlinReco/pull/42))
  - Thank you for writing the text to appear in the release notes. It will show up
    exactly as it appears between the two bold lines
  - ...

* 2018-02-08 Daniel Jeans ([PR#38](https://github.com/iLCSoft/MarlinReco/pull/38))
  - Change order of hit relations: now from reco/digi -> sim [to be consistent with past practice and all other processors]
  - Add to/from type info to the relation collections
  - Fix warnings from the compiler

* 2018-02-08 Andreas Alexander Maier ([PR#33](https://github.com/iLCSoft/MarlinReco/pull/33))
  - This package is an extension to the IsolatedLeptonFinderProcessor. The default functionality is untouched.
    - Optionally, it dresses leptons with close-by particles. By default it dresses electrons and muons with photons.
    - The algorithm starts with the highest energy lepton and adds all photons (and, optionally, electrons) in a cone of a given size around to it. As the original, it creates a collection with the dressed leptons and another collections with all remaining particles, except the ones that were dressed into the leptons.
    - All compiler warnings are fixed

* 2018-03-13 Frank Gaede ([PR#43](https://github.com/iLCSoft/MarlinReco/pull/43))
  -  Fix for iLCSoft/LCIO#35

* 2018-03-28 Frank Gaede ([PR#46](https://github.com/iLCSoft/MarlinReco/pull/46))
  - Fix for the removal of DDSurfaces which have been merged into DDRec 
    -  includes from `DDSurfaces` -> `DDRec`
    - namespace `DDSurfaces` -> `dd4hep::rec`

* 2017-12-12 Frank Gaede ([PR#34](https://github.com/iLCSoft/MarlinReco/pull/34))
  - Remove all warnings of type "should be initialized in the member initialization list [-Weffc++]"
  - Remove all warnings of type "unused parameter 'run'" for processRunHeader( LCRunHeader*  /*run*/)
  - Remove all warnings of type "unused parameter 'evt'" for "check( LCEvent *  /*evt*/ )"

* 2018-03-23 Ulrich Einhaus ([PR#44](https://github.com/iLCSoft/MarlinReco/pull/44))
  - Compute_dEdXProcessor:
    - geometry issue: gear to DD4hep unit adaption, fixed low momentum problem
    - added dx calculation strategies
    - added various processor options, default are all old version
    - added documentation

# v01-22

* 2017-11-15 Ete Remi ([PR#30](https://github.com/ilcsoft/MarlinReco/pull/30))
  - New SDHCAL digitizer version from Lyon group (ggarillot)
    - step linking and 'Angular Correction' 
    - two processors, one for applying the threshold, one for the threshold energy factors

* 2017-11-15 Shaojun Lu ([PR#29](https://github.com/ilcsoft/MarlinReco/pull/29))
  - replace gear with DD4hep in FourMomentumCovMat and PIDTools
        - use DD4hep for accessing BField and LayeredCalorimeterData extension to replace Gear.
         - with this the ILD standard reconstruction does no longer need a  gear file

# v01-21-01

* 2017-11-10 Ete Remi ([PR#28](https://github.com/iLCSoft/MarlinReco/pull/28))
  - Missing memory allocation and delete for random engine
  - Added delete specification for copy constructor and assignment operator to avoid warning on compilation
  - Missing delete for two arrays causing memory leaks

# v01-21

* 2017-10-26 KURATA Masakazu ([PR#26](https://github.com/ilcsoft/MarlinReco/pull/26))
  - improved PIDTools
        - added additional smearing functon for dE/dx resolution.
        - and corrected the strange behavior for PID. That is most of the muons are identified as pions.

# v01-20

* 2017-07-20 Andre Sailer ([PR#17](https://github.com/iLCSoft/MarlinReco/pull/17))
  - TauFinder: fix memory leak (few kB/event)

* 2017-09-13 Frank Gaede ([PR#22](https://github.com/iLCSoft/MarlinReco/pull/22))
  - add a default value for parameter inputHitCollections in RealisticCaloDigi.cc
  - create package EventShapes_Fortran, Fixes #20 
         - moved YThres from EventShapes to EventShapes_Fortran
  - build EventShapes for C++ now

* 2017-09-27 libo929 ([PR#24](https://github.com/iLCSoft/MarlinReco/pull/24))
  - improve SDHCALDigi/src/SimDigital.cc
       - Fixed the exception of unknown name y
        - add optional parameter HCALCellSize to overwrite the one from dd4hep

* 2017-07-29 libo929 ([PR#18](https://github.com/iLCSoft/MarlinReco/pull/18))
  - Updated the SimDigital processor to make it compatible with lcgeo simulation
  - Removed the digitization for ECAL 
  - Fixed the warnings

* 2017-08-14 Ete Remi ([PR#19](https://github.com/iLCSoft/MarlinReco/pull/19))
  - BruteForceEcalGapFiller : Add missing information on calo hit type

* 2017-09-15 Daniel Jeans ([PR#23](https://github.com/iLCSoft/MarlinReco/pull/23))
  - modify algorithm to fill ECAL gaps in BruteForceEcalGapFiller:
        - overall linear suppression gap hit energies [as was done in DBD]
        - additional logarithmic suppression to reduce large energy gap hits
        - separate parameters for (a) gaps between modules, and (b) gaps within modules
        - default values of new parameters seem reasonable for ILD_l4_v02 model

* 2017-10-06 Andre Sailer ([PR#25](https://github.com/iLCSoft/MarlinReco/pull/25))
  - Drop unused and no longer existing header includes AidaSoft/DD4hep#241

# v01-19-01

* 2017-07-04 Daniel Jeans ([PR#16](https://github.com/iLCSoft/MarlinReco/pull/16))
  - update resolution formula in TPC digitisation (from Dimitra Tsionou)

# v01-19

* 2017-05-09 Andre Sailer ([PR#8](https://github.com/iLCSoft/MarlinReco/pull/8))
  - CLICPfoSelector: remove dependency on GEAR File; simplify calculation of TimeAtECal by using trackState at Calorimeter
  - CLICPfoSelector: fix coverity defect
  - TauFinder: fix coverity defects; use streamlog instead of cout

* 2017-06-23 Andre Sailer ([PR#13](https://github.com/iLCSoft/MarlinReco/pull/13))
  - TauFinder:: PrepareRecParticles: change second parameter with same name RecCollection --> RecCollection_Tracks
  - V0Finder: Change second parameter to its correct Name RxyCutGamma --> RxyCutLambda

* 2017-06-20 Andre Sailer ([PR#11](https://github.com/iLCSoft/MarlinReco/pull/11))
  - Adapt to changes in namespaces and LCDD -->  Detector

* 2017-06-20 Shaojun Lu ([PR#10](https://github.com/iLCSoft/MarlinReco/pull/10))
  - Set default from 'ON' to 'OFF' to build MarlinReco fortran sources.
  - It could be switch on in the release-ilcsoft.cfg if you want.
  - https://github.com/iLCSoft/iLCInstall/blob/master/releases/HEAD/release-ilcsoft.cfg#L116

* 2017-06-26 Andre Sailer ([PR#15](https://github.com/iLCSoft/MarlinReco/pull/15))
  - RecoMCTruthLinker: clean collections that are only created for linking. Move cleanup to end of process event instead of trackLinker so things are also cleaned up if there are no tracks. Fixes memory leaks

* 2017-06-26 Frank Gaede ([PR#14](https://github.com/iLCSoft/MarlinReco/pull/14))
  -  replace Gear with DDRec and fix warnings in
          - Compute_dEdxProcessor
          - V0Finder
          - KinkFinder
  - remove obsolete sub-packages:
         - BrahmsTracking
         - FullLDCTracking
         - SiliconTracking
         - TrackCheater
         - TrackbasedPflow
         - Wolf
         - BCalReco
         - ETDDigi
         - FTDDigi

* 2017-06-01 Matthias ([PR#9](https://github.com/iLCSoft/MarlinReco/pull/9))
  - read magnetic field from dd4hep in taufinder particle preparator
  - remove leftover condition for former printout of the taufinder, which prevents the desired increment of the iterator, and thus deletes wrong taucandidates when filling the reconstructed tau's

* 2017-04-22 Andre Sailer ([PR#7](https://github.com/iLCSoft/MarlinReco/pull/7))
  - Re-enable compilation of CLICPfoSelector, CLICDstChecker
  - Fix Warnings in  CLICPfoSelector, CLICDstChecker

* 2017-04-22 Andre Sailer ([PR#6](https://github.com/iLCSoft/MarlinReco/pull/6))
  - Added additional TauFinder processor from Astrid Munnich, see doc/TauFinder/TauFinderLCDNote.pdf for description and performance
  
  - Ignore Warnings from external headers in MarlinReco

# v01-18

* 2017-03-30 Andre Sailer ([PR#4](https://github.com/iLCSoft/MarlinReco/pull/4))
  - RecoMCTruthLinker: move encoder initialisation from constructor to function body to avoid accessing encoding_string too early

* 2017-03-24 Emilia Leogrande ([PR#2](https://github.com/iLCSoft/MarlinReco/pull/2))
  - Replace ILDCellID0 to LCTrackerCellID

* 2017-04-12 Shaojun Lu ([PR#5](https://github.com/iLCSoft/MarlinReco/pull/5))
  - update to latest changes in LCIO
         - removed deprecated SimTrackerHit::setCellID
         - don't use assignment for CellIDDecoder

# v01-17

S. Lukic
  - Adding SiPID package to MarlinReco. 
    For the moment it contains the processor which calculates dE/dx for 
    a silicon tracker and writes it to the lcio file.
    Primarily intended for CLICdet but generally applicable. 
    DD4hep is used to access geometry.
    CMakeLists.txt in MarlinReco was adjusted. The package is added only
    if DD4hep and MarlinTrk are found.

M.Petric
  - fix in TruthVertexFinder
 
D. Jeans
  - copy calodigi documentation into Realistic package (+ a few changes to help latex compile)

S. Lu
  - comment out 'ADD_MARLINRECO_PKG( ./Clustering/BCalReco)' for solving
    '*** Break *** segmentation violation' at the end of Marlin run.
    ToDo: To ckeck which ilcsoft package use the same resource (maybe root resource). 



# v01-16

A.Ebrahimi
  - add package ErrorFlow
    - the ErrorFlow processor computes jet-specific energy uncertainty
      by summing up the uncertainty of individual particles clustered
      in a jet object.


# v01-15

M. Kurata
   - Include dEdx angular correction, which is a feedback from Sviatoslav
   - add PDF creation processor and steering file
   - add hadron type likelihood and probability
   - PIDTools: Adding PIDVariables and PIDParticles to trunk
   - assign 999 to likelihood when cannot perform PID

F. Gaede
   - introduced anonymous namespace to compile w/ -pedantic
   -  made compatible with c++11
   - removed -ansi -pedantic -Wno-long-long
   - build sub-package Relaistic only if DD4hep is found
   - added implementation of getExtension() - copied from DDMarlinPandora
  
D. Jeans
   - introduction of RealisticCaloDigi/Reco processors: refactorisation and cleaning up of ILDCaloDigi
   - gap corrections split off into BruteForceEcalGapFiller processor; change to setup of reco hit relations
   - fixed logic in determination of calibration by layer

M. Chera
   - fixed bug such that the processor now skips events with empty PFO collections without crashing
   

# v01-14

M. Berggren
   - Much ameliorated linking for clusters, taking back-scattering and
     decays in flight into account in a more correct way when finding
     what MCP to connect to a sim-hit.
     Work-around for the inverted meaning of the vertexIsNotEndpointOfParent
     method in Mokka. Set processor parameter  InvertedNonDestructiveInteractionLogic
     to activate the work-around (default is false)
     More consistent debugging. All levels from 0 to 9 used.

   - improved RecoMCTruthLinker
     Changed weight in hit<->mcparticle relation to the CaloHit energy (instead of the
     SimHit one). Also means the sim-hits not digitised (ie. without related CaloHit)
     are ignored. This also feeds through to the weights in the Cluster<->MCparticle
     relation, which should now be more correct: correct calibration factors applied
     for each hit.
     When doing the final relinking in problematic cases, explicitly check *only* back-scatters.
     Added MCP->cluster link for the case with the fixup of the missing truth-link for
     the LCAL in the DBD samples. Before only cluster->MCP was filled in this case.
   - Correted the true-particle to calorimeter-hit link for the case that
     the original true particle has been squeezed out by the cluster-linker,
     typically because it was created before the calorimeter, but after the
     last sensitive layer in the tracker.


  - AddClusterProperties:
     Use the PFO energy (not the cluster) when assign the PFO 4-momentum.
     Also add the sum of weights (=hit energy) to the cluster shape parameters.
     (In the DBD production, this is the same as the cluster energy, but it
     is not guaranteed to be the case).
   - Better/safer use of collection-parameters for indexing in sub-detector energy and shape-parameter arrays.
   - Different/simpler(?) calculation of Cov(4mom) in PFO.

F.Gaede 
 - changed log level for "processing event" from MESSAGE to DEBUG8
   in TPCDigiProcessor.cc
   -> use Statusmonitor instead



G.Wilson
   - added GammaGammaCandidateTruthFilter
    New processor to select those GammaGammaCandidates (fitted pairs of Pandora PFO photons) based 
    on MC information indicating actual parentage from a true pi0, eta or eta' meson decaying to gamma gamma.
    Output collections are called TrueGammaGammaPi0s etc.
    Cuts are made on the parent meson being prompt (production vertex distance < 10 cm), and the amount of the PFO 
    energy coming from the MC photon following similar criteria developed in ILDPerformance/pi0.
    This aspect will need some further investigation/tuning.


J.List
 - improve TrueJet:
 - First version that runs also over Higgs samples. No crash/obvious problem in several 1000s of events, 
   but more sophisticated checks pending. When running with full RecoMCTruthLink, estimate of of seen of true 
   energy needs fix, both on standard Whizard and Higgs samples


 M. Kurata
    - change overall structure for particle ID in LikelihoodPIDProcessor
    - change dEdx error estimation

 
M.Berggren:

    Added AddClusterProperties processor, a processor that adds cluster
    position and direction (w/uncertainties) to clusters, energy error 
    (if not already set) to clusters and 4-mom (w/ covariance)
    to neutral PFOs. Depends on a new MarlinUtil class (WeightedPoints3D),
    so for now, the processor is out-commented in CMakeLists.txt, to
    avoid problems if MarlinUtil isnot up-to-date !


G.Wilson 
  - update GammaGammaCandidateFinder
    - compute and fill ReconstructedParticle covariance matrix
  - Add GammaGammaSolutionFinder and DistilledPFOCreator processors
  - Add some steering options - and fit quality checks
  - Add 1st implementation of graph based solution finding in GammaGammaSolutionFinder - depends on BOOST

T.Calancha
  - The correlations of the energy (e) with x,y,z,e were stored in the wrong position of the matrix
    in FourMomentumCovMat. Thx to Ali Ebrahimi.

D.Jeans 
  -  electronics dynamic range applied before unfolding of PPD saturation
    in ILDCaloDigi.cc /  ScintillatorPpdDigi.cc

H.Sert
  - a separate algorithm was created for low P mu-pi ID in LikelihoodPIDProcessor.cc

 J.Tian
 - leave skipping or not an event without isolep to user in IsolatedLeptonTaggingProcessor.cc

T.Suehara
  -  bugfix on track charge in TrackToRecoParticleConverter

 J.List
  - added dummy error functions to ClusterShapesMR for MarlinPandora development, to be filled with life once proper calculations 
    are available
  - made SimpleFCalDigi use MarlinUtil:caloTypeFromString(_caloType)


# v01-13

  F.Gaede
   - renamed SimpleLHCalDigi to SimpleFCalDigi
   - deleted SimpleLCalDigi
   => old steering files need to be updated !!

  M. Berggren 
   - added  parameter FixLCalHits SimpleLHCalDigi in order to:
    - fix of the wrong xyz for hits in LCal. LCal CaloHits (collection LCAL) 
      are now within a few micron of their correct position. LCal SimCaloHits 
      (collection LumiCalCollection) still have the flawed xyz.

  J.List
  - updated SimpleLHCalDigi to include functionality fromm SimpleLCalDigi
   -> can be used for all fwd calos

  J.Tian
   - added new package:
     Analysis/OverlayRemoval

  H.Sert
  - new processor low momentum pi/mu separation:
    Analysis/PIDTools/LowMomentumMuPiSeparationPID_BDTG

  G.Wilson
  - new photon finder package: 
    Analysis/GammaGammaCandidateFinder
    - provides a 4-vector w/  mass-constrained fit (using MarlinKinfit  - new dependency )

  M Berggren
  - fixed links in RecoMCTruthLink/src/RecoMCTruthLinker.cc
  - full links in both directions (true->seen, seen->true) for PFOs, Tracks and Clusters.



# v01-12

  J.Tian
   - added new Analysis/IsolatedLeptonTagging
      MVA based isolated lepton finder

  T.Calancha
    - fixes for Analysis/FourMomentumCovMat
       - printout and read-only exception

  M.Kurata
  - update PIDTools for several kinds of likelihood calculation

  S.Bilokin
   - added Analysis/VertexChargeRecovery

  T.Suehara:
   - Analysis/TauFinder added

  F.Gaede

  - added copy of ClusterShapes from MarlinUtil 
    -> use "include "ClusterShapesMR.hh"
       and marlinreco::ClusterShapes
  - applied to ComputeShowerShapesProcessor.cc
  - eventuall changes should be merged back to MarlinUtil
  
  - fixed warning  -Wc++11-narrowing in TPCDigitizer




# v01-11

     LDCCaloDigi/ILDCaloDigi: D.Jeans & O. Hartbrich
     - improved realism of silicon and sintillator/PPD calorimeter hit digitisation
     - ScintillatorPpdDigi class defines PPD response model
     - default behaviour is no realism: must be switched on in steering file
     - added documentation (also described in LC-TOOL-2014-011)
     - can regroup the virtual cells of a scintillator strip
     
     Clustering/hybridEcalSplitter: D. Jeans
     - can deal with scintillator strips with >1 virtual cells

    - S.Bilokin: new package TruthVertexFinder added 
      ( see ./Analysis/TruthVertexFinder )

    - M.Berggren: added new packages TrueJet and TrueJet_Parser 
      (./Analysis/TrueJet, ./Analysis/TrueJet_Parser)

    - F.Gaede: add parameter DontEncodeSide to TPCDigiProcessor
        - allows to process old Mokka and new DD4hep based simulation 
       
    - F.Gaede: modified calo digitizers in CaloDigi/LDCCaloDigi:
      ILDCaloDigi, NewLDCCaloDigi, SimpleLCalDigi, SimpleLHCalDigi, SimpleMuonDigi
      added parameters
        CellIDLayerString, CellIDModuleString, CellIDStaveString
      to allow to switch between old cellId encoding string
      using K/K-1,M,S and new one using layer,module,sensor
      ( additionally added CellIDIndexIString, CellIDIndexJString to ILDCaloDigi)
  
    - C.Calancha: new package Analysis/FourMomentumCovMat
 
    - M.Kurata: new package Analysis/PIDTools: TPC track dEdx calculation, 
      Cluster shower profile extraction, and Likelihood Particle Identification

 
# v01-10
    - no release notes, use 
       svn -v log 
      to see changes ... 


# v01-09

     - improved SDHCAL digitizer (G.Grenier):
        charge simulation : simplifies Polya function. 
	Charge dispatching on pads : fix integration if using functions (slows the process)
	and add new dispatching using sum of Erf functions (recommended option)


# v01-08

 - John Marshall:
			Update fillECALGaps function, so that it remains identical to that implemented in NewLDCCaloDigi processor (i.e. copy the NewLDCCaloDigi updates from 9th Jan 2013).
			Apply HCAL "other" calibration constant to both ring and plug calorimeter hit types.


# v01-07
   
  - Manqi Ruan: added new clas G2D:
     Include G2CD, a general digitizer for gaseous hadron calorimeter. It takes 1mm simulated Hcal hits as input,
     and output digitized hits with tunable sell size. Efficiency and Multiplicity effects are also taken into 
     account, which can be adjusted in the steering file
  
  - J.List: fixed BCalTagEfficiency:
      adapted efficiency parametrisation for ECM=1TeV

  - Katsushige Coterra: updated NewLDCCaloDigi.cc
    for the Scintillator Strip Ecal reconstruction


# v01-06
   - updated IsolatedLeptonfinder (T.Tanabe)
     - LAL Lepton Finder included in IsolatedLeptonFinderProcessor

   - updated BCalTagEfficiency (J.List)
     - new option in BCalTagEfficiency: allows topick correct MCParticles from BCALMCTruthLink collection written by SGV

  - updated FPCCDClustering (D.Kamai)
     -  t.mori enable to difine pixel size for each layer.


  - fixed memory leak in BCalReco - fix provided by D.Jeans 
      - (loop variable in BCalReconstruction::Free3DArray )

# v01-05

   - BCalReco (A.Rosca)
     - in BCalReco.h reduce memory use: #define MAXNENTR 80000	      
     - introduced new variable EdepErr.
     - small changes to read new variable, EdepErr, from the background map.
     - fixed bug in SearchTowers().
     - write also empty BCal collections
     - fixed several programming flaws/reimplemented parts of code

   - BCalTagEfficiency (J.List)
     - removed local copy of TDR background map, updated README and bcal_ild_05_v05.xml
     - some minor clean-up of print statements, comments etc

# v01-04-02

    - BCalTagEfficiency:
        - threw out interpolation between BeamCal cells, introduced optional write-out of background map
        - uses LCRelationNavigator now, cleaned up example steering
        - write out also empty collections
        - added 3.5T background map and up to date example steering for ILD00

    - BCalReco:
        - improved calculation of Cartesian coordinates from cell position given
            in terms of ring/pad number.

    - RecoMCTruthLink:
        - exchanged ttrlcol and ctrlcol in makeSkim()
            ( will have no effect except code sanity )

    - LDCCaloDigi:
        - Bug fix to get output hits stored in output collection+increase speed of integration

    - FPCCDDigi:
        - t.mori enable to difine pixel size for each layer.


# v01-04-01
     - fixed TPCDigiProcessor :
        - added bfield correction factor to point resolution
          ( 4.0 / bField )

     - BCalRec (A.Rosca)
       - cleaned up code (removed unnecessary code)
       - introduced new default bg map file
       - enclosed code with  
             if( nHits > 0 ) {...}
         in order to prevent warning about missing CellIDEncoder string

       - fixed indentation of code (fg, using emacs defaults) (r4088)
         No code change, i.e. same code as r4087 (use this for diffs)

     - BCalTagEfficiency (J.List) 
       - veryfied unchanged performance on LoI samples (modulo a bug fix); steering flag for
         map format added


# v01-04

   - SiliconTracker, SiliconTrack_CLIC, TruthTracker: Corrected use of getPointOnCircle from HelixClass (MarlinUtil) 
     documentaion wrong said that a float[3] is required, when in fact a float[6] is needed. 

   - RecoMCTruthLinker: Improved Debug output for track mcparticle relations.
     			Fix for the case were no tracker collections are present.  

   - FPCCDDigitizer:    Bug fixed - process for the point on ladder edge.

   - BCalReco: 		Removed more uses of cout and cerr. Replaced with streamlog.
     			Corrected float instead of integer value for cluster position


# v01-03

   - added option to enable/disable output of root file containing histograms to BCalReco.cc
     (T.Tanabe)

   - fixed SIT parameters (new ZPlanarParameters )  in KinkFinder.cc (F.Gaede)

   - use streamlog for messages in EventShapes and SatoruJetFinder packages (T.Tanabe)

    - Analysis/RecoMCTruthLink: Corrected track <-> mcparticle weight regarding space points, nhits+=2. Reduce verbosity levels on MESSA
    GE.
    - Analysis/CLICPfoSelector: corrected logic to stop ihit90 index becoming negative.
    - removed/fixed several compiler/linker warnings


# v01-02

    - new package Analysis/IsolatedLeptonFinderProcessor
      (R. Yonamine, KEK, T. Tanabe, U.Tokiyo)
    - TrackDigi/FPCCDDigi/src/FPCCDDigitizer: bug fixed - unified the layer number in each function.
    - TrackDigi/TPCDigi/src/TPCDigiProcessor: use lcio::ILDDetID::barrel instead of 0
    - made compatible w/ clang++ compiler


# v01-01-01

    - bug fixes in Analysis/RecoMCTruthLink/src/RecoMCTruthLinker.cc:

        - The LCRelation Collections for TrackerHit to SimTrackerHits have been merged into a list of collections. Internally a single collection of LCRelations is created.
        - create track and cluster relations in any case as they are needed for the final recomctruth-link
	- The LCRelation Collections for TrackerHit to SimTrackerHits have been merged into a list of collections. 
	  Internally a single temporary collection of LCRelations is created.


# v01-01
      
      - RecoMCTruthLinker
        -  added mcTruthTrackLink: inverse relation from MCParticles to Tracks
           the weight will be the fraction of all sim-hits from the MCParticle
           that contributed to this track
        - changed weight definition for trackMCTruthLink:
          - weight is the fraction of all hits on the the track that have 
            contributions from this MCParticle
          - computed as sumSimHits_from_this_MCParticle/total_number_of_hits_on_track
             -> weight can be large than 1.0 (if more than one hit from an MCParticle
                would be used in a given layer (delta ray))
             -> this definition of weights is more accurate for defining the 'fake hit rate'
                as (1.-weight) in the case of many merged sim hits for a given track
         - changed steering logic, such that it can be more easily called for
           certain sub-tasks, eg. Track-MCTruth-Link only:
         - made the following output collections optional:
           trackMCTruthLink, clusterMCTruthLink, RecoMCTruthLink, calo-hit MCTruthLink, skimmed MCParticle
           -> collections won't be created if empty name specified (default)
         - removed options:
           OutputClusterRelation, OutputCalohitRelation, OutputTrackRelation
           they are now also based on name empty/not empty 
         - changed default for UseTrackerHitRelations to 'true'
            - print WARNING of set to false
           -- F.Gaede

      - LDCCaloDigi
          Implement use of LCIO calorimeter Hit step position for SDHCAL digitizer 
          (multiplicity simulation) and add parameters to the processor to better 
          control the digitization -- Gerald Grenier IPNL

      - TPCDigi
          Bug fix: added smearing of merged TPC hits
          (using rather large errors - which may have to be iterated on.)

       - BCalReconstruction
       	  Improved calibration (A.Rosca)

# v01-00-01


	- FPCCDDigi 
	  - Fixed setting of the cov matrix values.
	  - Modified to make link to mcp only for single particle event. 
	  - Chaged the algorithm to find pixels.

	- BCalReco
	  - Added reconstructed particle collection and check histograms.
	  - New parameter for bg map file: BackgroundFilename

# v01-00

	- General
	  - New release for preparation for the DBD.

	- Added hybridEcalSplitter (K. Kotera)

	- SiStripDigi: Processor to digitize and clusterize hits for the FTD subdtector
                       First release of the code, still missing full validation; (J. Duarte)

	- VTXDigiProcessor: corrected smearing to be along ladder, previously it was incorrectly perpendicular to the ladder, 
                            due to the use of the ladders phi angle instead of it angle of inclination

	- FPCCDDigitizer: Fixed bug in helix approximation.- in previous version the intersection between particle and ladder 
                          could be rotated to the opposite side.

        - RecoMCTruthLinker: Bug fix: protect against undefined mother pointer (T.Tanabe )
	  		     Modified to be able to use TrackerHit relations, though not yet the default. 
			     Enable by setting UseTrackerHitRelations true in steering file

	- TPCDigiProcessor: Corrected problem where hits from alternative z halves could be considered to be adjacent
	  		    Removed SimHits from RawHits and Used Relations. The deprecated feature to store them in the 
			    rawhits can be enabled by setting UseRawHitsToStoreSimhitPointer true. 

# v00-30

	- General
	  - Enabled build without cernlib/fortran (option MARLINRECO_FORTRAN=off) - e.g. for macos-64bit
	    sub packages (BrahmsTracking, Satoru, MarlinKinfit,...) that need cernlib will not be built

	- Clustering   
          - removed obsolete sub package ClusterCheater (ClusterCheater5_3 is newer)

          - deprecated TrackwiseClustering
            ( code still in repository, however not built by default )

          - added new sub package BCalReco 
            beam cal reconstruction (A.Rosca, DESY)


	- Pflow  
	  - deprecated Pflow package Wolf and TrackBasedPflow
            ( code still in repository, however not built by default )

	- CaloDigi
	  - added SimDigital processor to LDCCaloDigi 
            SDHCAL digitization (G.Grenier/R.Han, INPL)    

	- Analysis
	  - Removed MarlinKinfit processor

	  - RecoMCTruthLinker 
	    - Modification to allow Bremsstrahlung photons to be written to Skimmed list (off by default)

	- Tracking
	  - Use ILDCellID0 from lcio Util

	  - Removed deprecated setIsReferencePointPCA

	  - New processor FPCCDClustering.cc was added ( Daisuke Kamai )
        (FPCCDData and FPCCDPixelHit classes are in MarlinUtil )

	  - TPCDigiProcessor 
	    - Changed exception to warning for the case of points with the same x-y coordinates being passed to getPadTheta and getPadPhi. 
	      Corrected the check function which was previously using float == float, to use fabs(float-float) < tolerance. 
	      In this case tolerance has been set to a tenth of a micron


# v00-20

        - Tracking
          - MarlinTrackFit
           - Use double precision for the conversion from Tanagra to LC Track Parameters. Use of the HelixClass has been removed from the conversion.
          - FullLDCTracking
           - Use faster HelixClass::getDistanceToPoint where possible.
           - Improved dynamic mememory allocation
           - General clean-up and removal of redundant code.
          - LEPTracking 
           - Improved dynamic mememory allocation
           - Upper limit on track momentum now 2 TeV.
          - KinkFinder
           - Improved dynamic mememory allocation
           - Make sure prongDaughters[i] has entries.           
          - SiliconTracking
           - Protect against missing connection from SimTrackerHit to MCParticle.
           - Reducing maximum number of allowed FTD combinations
           - In case of many potential FTD hit combinations, iteratively adjust phi sectors for triplet formation to avoid long processing times 

        - TrackDigi
           - Improved dynamic mememory allocation

        - Analysis
          - MarlinKinfit
           - Added new test classes IterationScanner and ParameterScanner
           - Added NewFitterGSL.
          - CLICPfoSelector
           - Put in possibility of applying time of flight corrections to hit times: CorrectHitTimesForTimeOfFlight
           - Add one final cut to CLICPFOSelector to remove high energy neutral hadrons in the far forward region which are formed from mulitple background particles.
           - Improved monitoring of background/physics pfos.
           - Added simple analysis processor to loop over pfos separated in to "physics event" and background.
           - Create Selected PFO Collections as a subset of All PFOs.
          - RecoMCTruthLink
           - Use push_back instead of directly assigning entries to the vector, in case we run out of reserved space.
           - Output cluster truth linked collections for Clusters and Hits can only be created if there are clusters/hits.

# v00-19

 	  - added new sub package FPCCDDigi (Daisuke Kamai)

	  - Analysis/MarlinKinfit  (M.Beckmann) : 
            several updates, new RootTracer, name changes (PhotonFitObject ->
	    SimplePhotonFitObject, PhotonFitObjectPxyg -> ISRPhotonFitObject)

       - TPCDigitizer (S.Aplin) 
         - protected against accessing NULL pointer from getMCParticle
   	    - clean up unused variables and unneeded ifdefs
	    - fixed unsigned int == int comparison warnings
         - hit smearing (M.Thomson) 

       - NewFTDDigiProcessor
         - Add NewFTDTrackDigi with smearing for strips or Pixels, depending LayerID

       - ILDCaloDigi.cc ( M.Thomson )
         - Added new digitiser for ILD calorimeters with timing cuts and ability to 
           treat barrel and endcap differently

       - SimpleMuonDigi( M.Thomson) 
         - improve calibration of muon hits and set a maximum energy for a single hit (factor
	      two improvement in energy resolution). The defaults are appropriat e for ILD00       

       - MaterialDB
         - new version for CLIC CDR studies CLICCDRMaterialDB
         - made gear parameters for SIT/SET (support) layer thicknesses DoubleVec
         - increased the max number of tracking surfaces in F77, see fkparm.inc

       - BrahmsTracking
         - removed tpc Ionisation Potential in F77
         - increased the max number of tracking surfaces in F77, see fkparm.inc

       - FullLDCTracking ( M.Thomson )
         - modifications to reduce the number of split or "ghost" tracks

       - V0Finder (M.Thomson) 
         - improved checking to avoid false positives

 	    - BCalTagEfficiency.cc: fixed memory leak bug (C. Bartels)

        - RecoMCTruthLinker (M.Berggren)
         - added additional links (LCRelations) between ReconstructedParticles and Tracks and Clusters
         
          - code fixes (J.Engels) 
            - made gcc 4.4 compliant
	    - fixed various issues with gfortran and g2c
            - renamed get/setdEdx() to get/setEDep() in TrackDigitizers 
            - fixed cmake issues          


# v00-18
      - added new package Tracking/KinkFinder  (M.Thomson, J.Marshall)

      - moved sub-package BCalTagEfficiency from MarlinAna to ./Analysis

      - Tracking/V0Finder : fix for Lambda0bar (M.Thomson)

      - improved compatibility for gcc 4.x
        - g2c/gfortran in sl5 with gcc 4.1.1
        - fixed CaloDigi/LDCCaloDigi/src/CHT_helper.cc for gcc 4.4.3

     bug fixes:

      - CaloDigi/SimpleLHCalDigi.cc, CaloDigi/SimpleMuonDigi.cc,
         added missing output relation collection to the event

      - bug fixes in SimpleMuonDigi, (New)LDCCaloDigi wrt. to
        ecoding of calorimeter layout in CalorimeterHit::type 
        (was allways CHT::any)

      - fixed bug in TPCDigitizer with hit-MCParticle association

      - PFOID/src/PFOID.cc: addedd missing information
	that order of pdfs must be consistent with parameter
	EnergyBoundaries in parameter description   

      - TrackDigi/TPCDigi/src/TPCDigiProcessor.cc: fixed a bug in the
	logic when checking for hits from the same MCParticle 
        (last hit compared with itself)


# v00-17-02
     - made compatible with MacOS
     bug fixes: 
     - incorrect library version numbers
      


# v00-17
        removed dangerous revision of impact parameters 
        d0 and z0 in MarlinTrackFit.cc

# v00-16

     *	MarlinKinfit (boehmej, beckmann & K. Fujii)

      -	Global covariance matrix now also filled by NewtonFitter. 
      -	Added Tracer classes and access to global covariance matrix in
	BaseFitter, plus some bug fix in PhotonFitObjectPxyg and
	PhotonFitObject 

      -	Various bug fixes.

      -	Patches applied for compiling on mac osx.

	
     *	TrackDigi (gaede, aplin, raspereza, deMasi & fujii)

      -	VTXDigi, cluster parameter class for user extension of
	SimTrackerHits. 
      -	Introduced coordinate transform for direction
	vectors in VXDGeometry. New class for coordinate transform between
	ladder and lab frame. 
      -	Salt'n pepper noise hits with cluster
	parameters. 
      -	Only keep link to SimTrackerHit if MCParticle pointer
	is not null. 
      -	Evolved version of TPCDigi that provides additional functionality
	to deal with background. Coupled to the Mokka Sensitive Detector
	Driver. TPCSD03.cc. Plus considerable clean up.
      -	VTXDigitizer.cc, added possibility to store for each digitized
	tracker hit information on the fired pixels.
      -	Initial version of VTXBgClusters.
      -	ETDDigiProcessor, VTXDigiProcessor, patches applied for compiling 
	on mac osx. 

	
     * 	CaloDigi (thomson)
      -	Updates for digital and semi-digital HCAL readout of DHCAL 

	
     *	PFOID (aplin)

      -	Update the hcal calo type determintation to use Helper class from
	MarlinUtil for decoding/encoding lcio::CalorimeterHit types for
	the ILD. 


     *	ZFinder (thomson)

      -	New package for finding Z->l+l-.


     *	BrahmsTracking (aplin)

      -	Removed hardcoded cut of chi^2<300.


     * 	RecoMCTruthLink (gaede)

      -	Protect against MCParticle null pointers.

      -> see Changelog for details 

# v00-15
     * PFOID (A.Raspereza) 
       - Updated version of PFOID processor. Energy dependent scheme of 
         particle identification is implemented.
       - bug fix in Histogram.cc (array out of bounds)
       -  New pdf files for energy dependent particle ID procedure

# v00-14
   * modified calorimeter digitizers to encode CalorimeterHit types
     with CalorimeterHitType class from MarlinUtil ( need v00-13 or higher )

# v00-13
  * TrackDigi (F.Gaede)
    - ETDDigiProcessor/VTXDigiProcessor: 
      added active(SET)Layers parameter to mimic stereo layer readout by ignoring some layers

  * LDCCaloDigi  (M.Thomson)
     - SimpleMuonDigi:  fix default cell id encoding
     - SimpleLHCalDigi: new digitizer
     - NewLDCCaloDigi: improved Calo
       digitiser for ECAL and HCAL which writes out a single collection
       for each input collection to allow for possible CellID encoding

  * MarlinKinfit (M.Beckmann)
     - Replaced OPALFitter and NewtonFitter by OPALFitterGSL and
       NewtonFitterGSL. New PhotonFitObject to describe ISR photons.
       Several Bug Fixes.

  * bug fixes (J.Engels,S.Aplin)
    - removed compiler warnings in BrahmsTracking (cfortran.h,fkexts.inc)
    - CMake: bug fix: libg2c in 64 bit, added -m32 to global link flags
    - BrahmsTracking: fixed of bounds checking, buffer sizes


# v00-12-01
     - bug fix release:

      BrahmsTracking (S.Aplin)
       tflnks.F:  fixed a bug where a
	  number of TPC padrows between N*32+1 and N*32+4 caused division
	  by 0 when looking for hits in the last four padrows

       FullLDCTracking corrected the GEAR definition for the disk z positions

# v00-12-00

   * improved PFOID  (S.Aplin, A.Rasperezza) 
    - fixed to compile w/ gcc4.x 
    - output particle collection dropped - ParticleID object is assigned
      directly to ReconstructedParticle
    - improved documentation 

   TrackDigi/VTXDigi:
   * fixed bug in VTXNoiseHits (F.Gaede)

   * new processor CCDDigitizer (K.Harder)
     development version of vertex digitizier specific for CCD
     technology. written by Stefan Uebelacker (RAL) based on
     VTXDigitizer and on CCD digitizer code by Nick Sinev. not
     recommended for general use yet.

   * FullLDCTracking: New option for track fit included

   * made cmake 2.6 compliant
   * added 32 bit compatibility build option
   * added LIBRARY_DIRS

                                    
# v00-11-00

     * FullLDCTracking: (A.Raspereza)
        added cut on #Silicon hits if no TPC hits present - default: >=4

     * added code for PFOID  (M.Ohlerich et. al)
       particle id package
       -> experimental - does not compile and is not built into library

# v00-10-04
      
      * RecoMCTruthLinker (F.Gaede)
	  added paramater daughtersECutMeV for supressing low energy deltas od decays in flight

      * bug fixes and improvements in BrahmsTracking FullLDCTracking (S.Aplin)
        - splitted non-looping tracks based on hit-helix proximity 
        - corrected conversion of covariance matrix
        - added steering parameter AlwaysRunCurlKiller in BrahmsTracking
        - bug fix in  clear function of tkhitbank

        -> see Changelog for details 

# v00-10-03

	* bug fix SiliconTracking/LEPTrackingProcessor (S. Aplin)
          subdetectorHitNumbers  for ETD and SET hits

        * improved V0Finder (A.Raspereza)
          added cut on the mass of the candidate vertices

# v00-10-01

    * Tracking/FullLDCTracking/src/FullLDCTracking.cc: Si and TPC
      tracks are allowed to be merged if fit error = 0 ( proposed by
      Mark ).

# v00-10

    * new package Tracking/V0Finder  (A.Raspereza)
      processor to perform identification of neutral vertices,
      originating from photon conversions and decays of K0S and
      Lambda0

# v00-09-01
      fixed bug in RecoMCTruthLinker
	  and recover missing MCParticle link for Lcal clusters

# v00-09

    * RecoMCTruthLink:   (M.Beckmann, F.Gaede)
      - new processor MCTruthJetEnergy: adds colleciton parameter 'MCTruthJetEnergies'

    * TrackDigi/VTXDigitizer  (A.Raspereza)
      - bug fix for phi0 of gear::VXDParameters

    * BrahmsTracking/LEPTrackingProcessor 
      - bug fix w/ wrong TPCHit collection
      - changed default binning for the emoval of hits ala CurlKiller
      - bug fix for iteration of CurlKiller

# v00-08
     
    * TrackDigi: 
       - new processor VTXNoiseHits.cc   (F.Gaede)
       -> creates random noise hits with layer dependent density
          
       - added  ETDDigiProcessor (copy of FTDDigiProcessor w/simple
	  gaussian smearing in x-y )

    * improved LEPTracking (S.Aplin)	
      - error flag in LCIO event if patrec fails
      - produce only TPCTracks collection
      - functionality of CurlKiller moved into LEPTracking
      - improved log messages (streamlog/debug flags f77)
      - increased f77 array sizes for hits

    * improved FullLDCTracking (S.Aplin, M.Thomson)
       - correction for non looping split tracks
       - bug fixes

    
    * additional bug fixes

       -> see Changelog for details
 
# v00-07

     * improved FullLDCTracking (A.Raspereza)
       include SET and ETD in pattern recognition and fit

     * improved SiliconTracking (A.Raspereza, S.Aplin))
        - fixed/improved FTD geometry code 
        - improved efficiency in patrec
 
     * added LDCCaloDigi package (M.Thomson)
       calo digi code with gap corrections


     * additional bug fixes and code cleaning

       -> see Changelog for details

# v00-06
	 
    * improved LEPTracking and TPCDigi(S.Aplin)
	- reduced printout and moved cout to debug streamlog_out
        - corrected covariance matrix definition of TrackerHits
        - removed VTX and SIT Hits from LEPTracking

    * removed obsolete package VertexTracking (use SiliconTracking)

    * new package Analysis/RecoMCTruthLink (F.Gaede)
       - create link between ReconstructedParticles and MCParticles 
       - create skimmed MCparticle list

    * improved SiliconTracking (A.Raspereza)
       - reduced printout	
       - improved Silicon track finding efficiency (hit triplet fitting)
       - bug fix for correct assignment of left-over Silicon hits


     * additional bug fixes

       -> see Changelog for details

# v00-05-03

       bug fix release 

	* TrackDigi/TPCDigi/TPCDigiProcessor:    (S.Aplin)
	  - fixed several memory leaks 

	* Tracking/SiliconTracking/include/MaterialDB (A.Raspereza)
	  - added support structures for SIT
	  - included treatment of the first layer of heavy-weight Si for FTD
     	    detector


# v00-05-02
  
       bug fix release 

          !!! requires MarlinUtil v00-06 !!!!


	* TrackDigi/VTXDigi/VTXDigiProcessor: (C. Lynch)

	  - vertex hits are now smeared on ladders 
 
	* TrackDigi/TPCDigi/TPCDigiProcessor:    (S.Aplin)
          - new parametrization of r-phi and z resolution smearing (from LC-TPC group) 
          - bug fixes 
          - added expert check histograms (enable w/ -D EXPERTCHECKPLOTS)
          - covariance matrix of TrackerHit are now x,y,z (used to be cylindrical)

	* Analysis/MarlinKinfit:  (J.List) 
	  - added soft constraints
          - added more example fitters
          - bug fixes

# v00-05-01
       bug fix release

        - MarlinKinfit improved and new example (J.List) 

        - fixed LEPTracking for large #pad rows in TPC (S.Aplin)
          ( currently 512 pad rows are allowed - increase by
            setting the value of N32BITREG in padrow.inc )

        - dvips for documentation does not send to printer (J.Engels)

        - added simple digitizers SimpleLCalDigi and SimpleMuonDigi (by M.Thomson)

        - ...

# v00-05
       - new processor YThresh, calculates the yThresh variable(B. Hooberman)
         * in Analysis/EventShapes
         * see README and doxygen documentation
 
       - modified TPC digitization for new TPC driver (S.Aplin)
       - removed silicon tracking from BrahmsTracking

       - new package MarlinKinFit (J.List)
	 * Analysis/MarlinKinfit
         * kinematic fits with constraints
         * see doxygen documentation       

       - improved TrackBasedPFlow (O.Wendt)
	 * example steering file for the
	   Track-Based Particle Flow processor
	 * bug fixes and improved performance
         * LEPTracking: use TrackerHit::getCovMatrix for resolution

       - make consistent use of GEAR (K.Harder, C.Lynch) 
         * TrackDigi/VTXDigi: introduce gear
         * use global GEAR B field setting instead of TPC parameters setting in:
           Analysis/CheckPlots, Pflow/TrackBasedPFlow, Pflow/Wolf, Tracking/BrahmsTracking
         * updated example gear files (Bfield)
	 * Tracking/SiliconTracking: adjustments for new
	   GEAR file format as of Mokka v06-04-p02

       - improved FullLDCTracking (A.Raspereza)
         * some bugs fixed, e.g. wrong reference point
         * detailed users manual: doc/FullLDCTracking_Manual.pdf
         * consistent use of gear::BField 
       	 * TrackDigi/TPCDigi: Modified class for digitization of the
	   TPC hits. Spatial resolutions as well as parameters tpcPixRP and
	   tpcPixZ, representing the physical pad width and time binning
	   respectively, are passed to the processor as external parameters
	   (implemented by Clare Lynch)
         * example steering file
         * introduced various fit options
         * improved performance
         * digitization of laddered VXD detectors
  
       - improved code and cmake build files (J.Engels)
         * removed unused variables (compiler warnings)
         * install target for header files
         * fortran compiler flags and libg2c
            
       - improved documentation (F.Gaede)
	 * automatic generation of API doc (doxygen)


    -> see ChangeLog for details 

# v00-04
   - cmake is now default build tool:  (J.Engels)
        # edit BuildSetup.cmake as needed
        mkdir build ; cd build
        cmake -C ../BuildSetup.cmake ..
        make install
     -> creates plugin library $MarlinReco/lib/libMarlinReco.so
       

   - improved Tracking  (A.Raspereza)
       Extended version of TrackCheater, FullLDCTracking and SiliconTracking.
       The numbers of hits from different subdetectors are stored using
       TRACK:subdetectorHitNumbers. trackfit.F has been modified to
       account for curling tracks.

   -  Clustering/PhotonFinderKit/EMShowerFinder  (O.Wendt)
	  Processor to find electro-magnetic showers. 
          It is based on the KIT package and takes
	  only ECAL hits into account.	The output is a collection of
	  clusters in which the electro-magnetic shower are stored

   - TrackBasedPFlow (O.Wendt)
       + debugging and improvement of performance
       + included EMShowerFinder for improved efficiency to find
	 electromagnetic showers

   - bug fixes
       + made compliant with SL4 (gcc3.4) 
       +...

    -> see ChangeLog for details 

# v00-03
  - improved Full Tracking  (A.Raspereza)

  - initial version of TrackBasedPFA (O.Wendt)

  - PhotonFinderKit (P. Krstonosic)

  - added cmake support (epxerimental)
 

 for details see the ChangeLog file


