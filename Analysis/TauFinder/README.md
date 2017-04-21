The Taufinder folder contains two processors for Tau finding, the
_TaJetClustering_ and the _TauFinder_ processors

See this README for TaJetClustering, see MarlinReco/doc/TauFinder for descriptions of the _TauFinder_ processor.


The processor _PrepareRECParticles_ processor can be used to create REC particles directly from MCParticles for testing tau reconstruction on MCTruth.
The processor _EvaluateTauFinder_ can be used to check the performance of the _TauFinder_ processor.


# TaJetClustering: processor for Tau finder

-- Simple explanation of TauFinder --
150708 Taikan Suehara / Kyushu University <suehara@phys.kyushu-u.ac.jp>

mainly targetted to obtain isolated tau from jet environment

Three steps:
1. Clustering
2. Primary cut
3. Cone cut

Parameters and default values:

0. Collections
  ```
    <!--Input PFO collection-->
    <parameter name="PFOCollection" type="string" lcioInType="ReconstructedParticle">PandoraPFOs </parameter>
    <!--Tau output collection-->
    <parameter name="OutputTauCollection" type="string" lcioOutType="ReconstructedParticle">TaJets </parameter>
    <!--Remained PFO collection not clustered-->
    <parameter name="RemainPFOCollection" type="string" lcioOutType="ReconstructedParticle">RemainPFOs </parameter>
  ```
1. Clustering
  ```
    <!--Tau mass for tau clustering [GeV]-->
    <parameter name="TauMass" type="double">2 </parameter>
    <!--Allowed cosine angle to be clustered-->
    <parameter name="TauCosAngle" type="double">0.98 </parameter>
  ```
2. Primary cut
  ```
    <!-- Skip ANY Primary and Cone cuts if true: should be only used in lepton-only final states! -->
    <parameter name="NoSelection" type="int">0 </parameter>

    <!--Primary cut include IMPLICIT selection of accepting only 1 or 3 tracks in jets:
              this loosen the counting of low energy tracks-->
    <parameter name="AcceptFlexibleLowEnergyTrack" type="int">1 </parameter>
    <!--Minimum jet energy to be accepted as taus-->
    <parameter name="MinimumJetEnergy" type="double">3 </parameter>
    <!--Minimum track energy to be accepted as taus-->
    <parameter name="MinimumTrackEnergy" type="double">2 </parameter>
    <!--Minimum track energy to be counted-->
    <parameter name="MinimumTrackEnergyAssoc" type="double">2 </parameter>
  ```
3. Cone cut: currently simple 1D cut
  ```
    <!-- No cone selection if true -->
    <parameter name="NoSelection" type="int">0 </parameter>

    <!--Minimum cosine angle for cone-->
    <parameter name="ConeMinCosAngle" type="double">0.9 </parameter>
    <!--Maximum cosine angle for cone-->
    <parameter name="ConeMaxCosAngle" type="double">1 </parameter>
    <!--Energy fraction of cone compared to central-->
    <parameter name="ConeMaxEnergyFrac" type="double">0.1 </parameter>
  ```
