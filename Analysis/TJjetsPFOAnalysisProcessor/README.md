TJjetsPFOAnalysisProcessor
====================
A processor that combines basic analysis of the PFOAnalysis processor with the jet analysis power of the TrueJet processor.

**Basic elements:**
  -  *PFOAnalysis* is part of LCPandoraAnalysis: https://github.com/PandoraPFA/LCPandoraAnalysis/
  - *TrueJet* and *TrueJet_Parser* are included in the Analysis part of MarlinReco: https://github.com/iLCSoft/MarlinReco

**Authors:**
  - The majority of this code is copy-pasted and slightly adjusted from the aforementioned basic elements, so first credit should be paid to their authors!
  - Jakob Beyer  <<jakob.beyer@desy.de>>

Usage:
---------------------
This processor can only be used in combination with the *TrueJet* processor. In the Marlin steering file first run *TrueJet*, then use its output to run this processor. The output will be a root-file with the event-information.


Description:
---------------------

The core idea of this code is to get a quick overview of how individual jets within the given sample behave and how they are reconstructed.

To do so one must first find the individual jets. For this the *TrueJet* processor must be applied to the events. Its results are interpreted with the use of the *TrueJet_Parser* tool which was integrated into this processor.

Then the basic analysis of the reconstructed and Monte-Carlo particles from *PFOAnalysis* is applied to the jets found by *TrueJet*.

  - **WARNING for working with the MCParticles:** The way the "reconstructable" MC particles are found had to be changed slightly from PFOAnalysis in order to avoid overthrowing the assignement to the TrueJet-jets.
