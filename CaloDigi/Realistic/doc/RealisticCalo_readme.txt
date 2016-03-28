RealisticCalo* is a re-orgainsation of the calorimeter digitisation and reconstruction code, which was previously in ILDCaloDigi.
It's also designed to work only with DD4hep geometry information.

The task has been split into several tasks: 

1) hit digitisation

conversion of geant4 energy deposit to some calibrated electronics signal.
e.g. for silicon, converted to the MIP scale, for scintillator to the SiPM pixel scale.
effects such as miscalibration, dead cells, pe statistics are applied at this stage

2) hit reconstruction

the conversion of these hits to the EM energy scale.
This includes unfolding of SiPM response, and correction for the calorimeter sampling ratio.

** note that the calibration factors given here are for converting from the MIP scale to the total deposited energy scale.
** this is different to LDC/ILDCaloDigi, where the conversion was from the energy deposited in the active material to total deposited energy.

3) gap correction

If hits are situated directly across a dead area of the detector, their energies and the distance between them is used to guess how much energy was lost.
"Gap correction hits" are now written into a new calorimeterhit collection (previously, the energy of existing hits at wafer boundaries was increased).
For now, a BruteForceEcalGapFiller has been written, which uses the distance between hits to decide if they are situated across a dead region.
This is close to what was done previously in e.g. ILDCaloDigi.
In principle, it should be possible (and probably more efficient) to use cell indices to decide whether a hit is at the edge of a wafer.

----------------------

One instance of each processor should be used for each collection of hits (or set of collections with identical properties).

There are virtual base classes RealisticCaloDigi and RealisticCaloReco. 
Derived classes for a particular technology inherit from these: RealisticCaloDigiSilicon, RealisticCaloDigiScinPpd and RealisticCaloRecoSilicon, RealisticCaloRecoSilicon.

example of steering:

 <!-- digitisation for silicon ECAL barrel hits -->
 <processor name="ecalbarrelDigi" type="RealisticCaloDigiSilicon">
   <parameter name="inputHitCollections"> EcalBarrelCollection </parameter>
   <parameter name="outputHitCollections"> EcalBarrelCollectionDigi </parameter>
   <parameter name="outputRelationCollections"> EcalBarrelCollectionDigiRelation </parameter>
   <parameter name="threshold"> 0.5 </parameter>
   <parameter name="thresholdUnit"> MIP </parameter>
   <parameter name="timingCut"> 1 </parameter>
   <parameter name="timingCorrectForPropagation"> 1 </parameter>
   <parameter name="timingWindowMin"> -10 </parameter>
   <parameter name="timingWindowMax"> 100 </parameter>
   <parameter name="calibration_mip"> 1.7e-4 </parameter>
   <parameter name="elec_noise_mip"> 0.1 </parameter>
   <parameter name="elec_range_mip"> 2000 </parameter>
   <parameter name="CellIDLayerString"> layer </parameter>
   <parameter name="silicon_pairEnergy"> 3.6 </parameter>
 </processor>

 <!-- reconstruction for silicon ECAL barrel hits -->
 <processor name="ecalbarrelReco" type="RealisticCaloRecoSilicon">
   <parameter name="inputHitCollections"> EcalBarrelCollectionDigi </parameter>
   <parameter name="inputRelCollections"> EcalBarrelCollectionDigiRelation </parameter>
   <parameter name="outputHitCollections"> EcalBarrelCollectionReco </parameter>
   <parameter name="outputRelCollections"> EcalBarrelCollectionRecoRelation </parameter>
   <parameter name="calibration_layergroups"> 20 100 </parameter>
   <parameter name="calibration_factorsMipGev"> 6.58e-3 13.16e-3 </parameter>
   <parameter name="gap_correction"> 1 </parameter>
   <parameter name="CellIDLayerString"> layer </parameter>
 </processor>

 <!-- digitisation for scintillator HCAL barrel hits -->
 <processor name="hcalbarrelDigi" type="RealisticCaloDigiScinPpd">
   <parameter name="inputHitCollections"> HcalBarrelRegCollection </parameter>
   <parameter name="outputHitCollections"> HcalBarrelRegCollectionDigi </parameter>
   <parameter name="outputRelationCollections"> HcalBarrelRegCollectionDigiRelation </parameter>
   <parameter name="threshold"> 0.5 </parameter>
   <parameter name="thresholdUnit"> MIP </parameter>
   <parameter name="timingCut"> 1 </parameter>
   <parameter name="timingCorrectForPropagation"> 1 </parameter>
   <parameter name="timingWindowMin"> -10 </parameter>
   <parameter name="timingWindowMax"> 100 </parameter>
   <parameter name="calibration_mip"> 5.3e-4 </parameter>
   <parameter name="elec_noise_mip"> 0.1 </parameter>
   <parameter name="elec_range_mip"> 500 </parameter>
   <parameter name="CellIDLayerString"> layer </parameter>
   <parameter name="ppd_mipPe"> 15 </parameter>
   <parameter name="ppd_npix"> 1600 </parameter>
   <parameter name="ppd_npix_uncert"> 0 </parameter>
   <parameter name="ppd_pix_spread"> 0 </parameter>
 </processor>

 <!-- reconstruction for scintillator HCAL barrel hits -->
 <processor name="hcalbarrelReco" type="RealisticCaloRecoScinPpd">
   <parameter name="inputHitCollections"> HcalBarrelRegCollectionDigi </parameter>
   <parameter name="inputRelCollections"> HcalBarrelRegCollectionDigiRelation </parameter>
   <parameter name="outputHitCollections"> HcalBarrelRegCollectionReco </parameter>
   <parameter name="outputRelCollections"> HcalBarrelRegCollectionRecoRelation </parameter>
   <parameter name="calibration_layergroups"> 100 </parameter>
   <parameter name="calibration_factorsMipGev"> 3.0e-2 </parameter>
   <parameter name="ppd_mipPe"> 15 </parameter>
   <parameter name="ppd_npix"> 1600 </parameter>
   <parameter name="CellIDLayerString"> layer </parameter>
 </processor>

 <processor name="myBarrelGapFiller" type="BruteForceEcalGapFiller">
   <parameter name="inputHitCollection"> EcalBarrelCollectionReco </parameter>
   <parameter name="outputHitCollection"> EcalBarrelCollectionGapHits </parameter>
   <parameter name="CellIDLayerString"> layer </parameter>
   <parameter name="CellIDModuleString"> module </parameter>
   <parameter name="CellIDStaveString"> stave </parameter>
 </processor>
