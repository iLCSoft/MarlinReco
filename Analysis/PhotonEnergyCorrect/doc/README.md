## photonCorrectionProcessor

This is a processor to make energy corrections to photon PFOs.

Questions to daniel.jeans@kek.jp

Overview:
1. correct barrel overlap regions in phi
2. correct for barrel module boundaries in z
3. correct for module boundaries in endcaps
4. overall energy correction to restore linearity

corrections have been determined for the ILD_l5_o1_v02 model (for the cases in which the BruteForceEcalGapFiller both has and has not attempted to fill inter-module gaps).

DJ recommends that BruteForceEcalGapFiller is set to correct for gaps within modules, but not between them.

# extra explanation of processor parameters

```
 <processor name="MyphotonCorrectionProcessor" type="photonCorrectionProcessor">
 <!--photonCorrectionProcessor applies an energy correction to photon-like PFOs-->
  <!--verbosity level of this processor ("DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT")-->
  <!--parameter name="Verbosity" type="string">DEBUG </parameter-->
```
The name of the input PFO collection, and whether to adjust their energies
```
  <!--name of input PFO collection-->
  <parameter name="inputCollection" type="string">PandoraPFOs </parameter>
  <!--apply the corrected energies to the PFOs-->
  <parameter name="modifyPFOenergies" type="bool">true </parameter>
```
produce validation histograms (designed for single photon events)? and typical energy to set the histogram range
```
  <!--produce validation plots-->
  <parameter name="validationPlots" type="bool">false </parameter>
  <!--nominal photon energy (for validation plots)-->
  <parameter name="nominalEnergy" type="float">200 </parameter>
```
Which correction set to use? Some default sets are defined in photonCorrector (for ILD_l5_o1_v02 with and without gap hits between modules)...
```
  <!--use defaults correction parameters (<0: no ; >=0 : yes, as defined in photonCorrector)-->
  <parameter name="useCorrectorDefaultSet" type="int">1 </parameter>
```
...or you can set the various parameters by hand
In the barrel region the phi dependence is modelled as an asymmetric gaussian, whose mean position depends logarithmically on energy
```
  <!--barrel phi correction: gaussian depth-->
  <parameter name="phiBarrelCorr_depth" type="float">-999 </parameter>
  <!--barrel phi correction: central position (constant)-->
  <parameter name="phiBarrelCorr_pos_const" type="float">-999 </parameter>
  <!--barrel phi correction: central position (log(e) coeff)-->
  <parameter name="phiBarrelCorr_pos_logen" type="float">-999 </parameter>
  <!--barrel phi correction: gaussian width (left side)-->
  <parameter name="phiBarrelCorr_width1" type="float">-999 </parameter>
  <!--barrel phi correction: gaussian width (right side)-->
  <parameter name="phiBarrelCorr_width2" type="float">-999 </parameter>
```
the cos(theta) dependence is modelled as 3 Gaussians: the 2 intermodule gaps, and the barrel/endcap transition
```
  <!--barrel cos(theta) correction: gaus1: mean-->
  <parameter name="costhCorr_gaus1_mean" type="float">-999 </parameter>
  <!--barrel cos(theta) correction: gaus1: norm (constant)-->
  <parameter name="costhCorr_gaus1_norm_const" type="float">-999 </parameter>
  <!--barrel cos(theta) correction: gaus1: norm (log(e) coeff)-->
  <parameter name="costhCorr_gaus1_norm_logen" type="float">-999 </parameter>
  <!--barrel cos(theta) correction: gaus1: sigma-->
  <parameter name="costhCorr_gaus1_sigm" type="float">-999 </parameter>
  <!--barrel cos(theta) correction: gaus2: mean-->
  <parameter name="costhCorr_gaus2_mean" type="float">-999 </parameter>
  <!--barrel cos(theta) correction: gaus2: norm (constant)-->
  <parameter name="costhCorr_gaus2_norm_const" type="float">-999 </parameter>
  <!--barrel cos(theta) correction: gaus2: norm (log(e) coeff)-->
  <parameter name="costhCorr_gaus2_norm_logen" type="float">-999 </parameter>
  <!--barrel cos(theta) correction: gaus2: sigm-->
  <parameter name="costhCorr_gaus2_sigm" type="float">-999 </parameter>
  <!--barrel cos(theta) correction: gaus3: mean-->
  <parameter name="costhCorr_gaus3_mean" type="float">-999 </parameter>
  <!--barrel cos(theta) correction: gaus3: norm (constant)-->
  <parameter name="costhCorr_gaus3_norm" type="float">-999 </parameter>
  <!--barrel cos(theta) correction: gaus3: sigm-->
  <parameter name="costhCorr_gaus3_sigm" type="float">-999 </parameter>
```
In the endcap, we model the inter-module gaps within a quadrant as two gaussians:
```
  <!--across endcap module correction: gaus1 mean-->
  <parameter name="endcap_gaus1_mean" type="float">-999 </parameter>
  <!--across endcap module correction: gaus1 norm-->
  <parameter name="endcap_gaus1_norm" type="float">-999 </parameter>
  <!--across endcap module correction: gaus1 sigma-->
  <parameter name="endcap_gaus1_sigm" type="float">-999 </parameter>
  <!--across endcap module correction: gaus2 mean-->
  <parameter name="endcap_gaus2_mean" type="float">-999 </parameter>
  <!--across endcap module correction: gaus2 norm-->
  <parameter name="endcap_gaus2_norm" type="float">-999 </parameter>
  <!--across endcap module correction: gaus2 sigma-->
  <parameter name="endcap_gaus2_sigm" type="float">-999 </parameter>
```
extra energy correction for all PFOs in endcap region
```
  <!--extra correction factor for endcap-->
  <parameter name="costhCorr_endcap_scale" type="float">-999 </parameter>
```
Finally, an overall non-linear energy correction to give energy linearity
```
  <!--overall energy correction: constant term-->
  <parameter name="energyLin_const" type="float">-999 </parameter>
  <!--overall energy correction: log(e) coefficient-->
  <parameter name="energyLin_logen" type="float">-999 </parameter>
</processor>
```