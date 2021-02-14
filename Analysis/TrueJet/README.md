#  TrueJet: A Marlin processor to group particles into jets, using true information.

Author:

  Mikael Berggren <mikael.berggren@desy.de>

## Brief description:

  Use this Marlin processor to create collections and navigators
to find jets actually coming from separate initial systems.
Both true, seen, and true-of-seen jet-energies are accessible.
The correct total energy can be calculated. Correct di-jet
pairings giving W or Z are provided. Jet-flavour is indicated.
Inner brems-strahlung and gluon induced jets are flagged.
Particles from overlaid gamma-gamma to hadrons are identified
The only exception is that currently _initial_ gluon-jets (from H->g g)
cannot be treated.


To use the the output from the processor in further analysis,
it is suggested to use the helper-class TrueJet_Parser.

The only processor-parameters used are collection names for output and
input (mcparticles(skimmed), pfos , and the mc-reco link.
The latter two are gracefully ignored if absent, as is the case
if the input LCIO file is the output from the event generator.
There is also one boolean flag to switch between Whizard1 and
and Whizard2 format of the MCParticle (i.e. between mc2020 and IDR DSTs).
The defaults should work when running on mc2020 DSTs,
so normally no parameters are needed to be specified. For IDR DSTs,
and for generator-output LCIO files, the name of the MCParticle
collection needs to be adjusted, and for IDR DSTs, the flag Whizard1
must be set to `true`.

See the examples folder for suggested steering-files. Play with the
Verbosity level. Put to MESSAGE or WARNING to make it almost shut up.

## Detailed Description:

  The processor uses the MCParticle collection, and the
RecoMCTruthLink navigator as input. It recreates the Pythia event
record from the MCParticles, corrects in-consistencies (both those 
already present in  the input file and those created by
simulation). It identifies the _initial_ and _final_ colour-neutral
systems, meaning those at the beginning and end of the parton-shower,
respectively. It identifies the decay-products of these systems and
the particles giving rise to them. In each such system, the final
products are assigned to one of two jets, by looking at the angle to
either of the two particles creating the colour-neutral system.
All decay-products - also those created in simulation - of the
primary hadrons are assigned to the same jet as their ancestor.
The jet-identity is also propagated backwards through the parton-
shower, until either the initial, hard particles from the matrix-
element are reached, or it is found that the particles came from
a radiated gluon or IVB. In the process, inner-bremsstrahlung photons
are also found, and added to the jet of the particle that radiated
them.

The  _final_ colour-neutral systems always give rise to two jets.
This is also true if the colour-neutral system is a "cluster"
i.e. a bound state of initial quarks (normally the  colour-neutral 
system is a "string" i.e.. a system of two quarks and a chain of gluons).
Usually, a "cluster" materialises to a single hadron, so one of the
two jets will be _empty_ . However, both jets will be present earlier
in the parton-shower.

The _initial_ colour-neutral system  might make more than two jets, 
due to gluon radiation. In the this case, TrueJet keeps track of 
which jets ultimately came from a given _initial_ colour-neutral system. 
The _initial_  colour-neutral system is the one that is expected to come 
from a single IVB-fermion-antifermion vertex in the initial hard process.

The same logic also works for leptons: They, too, are grouped together
two-by-two into groups associated with the same initial boson.

ISR photons are assigned to separate (un-grouped) jets, as are photons
explicitly present in the hard interaction, i.e. the final state in the
Whizard process definition. This includes gammas from H->gammagamma or
H->Zgamma.

All post-parton shower MCParticles that were not assigned to
any jet by this process are from overlaid gamma-gamma->hadrons
events, or beam-strahlung pairs, and are grouped together into one single jet.

For further explanations of the method , see the commented listing in 
the examples folder and  [this talk](http://agenda.linearcollider.org/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=6526).
(note that some collection-names have changed wrt. that talk).

**Process-parameters**:

<u>Input collection names</u>


|   parameter name         |          type         |            default |
|-----                     |-------                | ------- |
|  MCParticleCollection    | MCParticle            |   MCParticles |
|  RecoParticleCollection  | ReconstructedParticle |  PandoraPFOs |
|  RecoMCTruthLink         | LCRelation            |  RecoMCTruthLink |

<u>Output collection names</u>


|   parameter name         |          type         |            default |
|-----                     |-------                | ------- |
|  TrueJets                | ReconstructedParticle  | TrueJets  |
|  FinalColourNeutrals     | ReconstructedParticle  | FinalColourNeutrals |
|  InitialColourNeutrals   | ReconstructedParticle  | InitialColourNeutrals |
| | | |
|  TrueJetPFOLink          | LCRelation             | TrueJetPFOLink |
|  TrueJetMCParticleLink   | LCRelation             | TrueJetMCParticleLink |
|  FinalElementonLink      | LCRelation             | FinalElementonLink |
|  InitialElementonLink    | LCRelation             | InitialElementonLink |
|  FinalColourNeutralLink  | LCRelation             | FinalColourNeutralLink |
|  InitialColourNeutralLink| LCRelation             | InitialColourNeutralLink |


** Marlin steering sections **:

```

 <execute>
           .
           .
           .
  <processor name="mytruejet"/>
           .
           .
           .
 (optionally:)
  <processor name="DSTOutput"/>
 </execute>
           .
           .
           .
<processor name="mytruejet" type="TrueJet">
 ... (add process parameters if for some strange reason the default collection names are not OK)
</processor>
```

If DSTOutput called: need to add the TrueJet outputs to the parameters of 
processor name="DSTOutput", by adding


-       TrueJets
-       FinalColourNeutrals
-       InitialColourNeutrals
-       TrueJetPFOLink 
-       TrueJetMCParticleLink 
-       FinalElementonLink 
-       InitialElementonLink 
-       FinalColourNeutralLink 
-       InitialColourNeutralLink

to the   <parameter name="KeepCollectionNames" type="StringVec"> 
section of the DSTOutput processor parameters. If full functionality
of TrueJet_Parser is wanted, at least

-       MCParticles
-       RecoMCTruthLink
-       MCTruthRecoLink

should also be in the list, but they would be in any useful DST-output anyhow....

** Contents of the output collections **:

  <u>TrueJets (ReconstructedParticle):</u>

   `getEnergy, getMass, getMomentum, getCharge` returns those quantities, calculated from the values
    of all PFOs connected to the true jet.
    
   `getParticles` returns the list of all PFOs in the jet.
   `getParticleIDs()[0]->getType` returns the jet type as:
   
       1  : jet from string
       2  : jet is lepton
       3  : jet from cluster
       4  : jet is ISR
       5  : jet is overlay
       6  : jet is a photon from the matrix-element
       


if the jet came from a quark from gluon-splitting, [(jet radiating the gluon)+1]*100 is added to the type.
`getParticleIDs()[0]->getPDG()` returns the jet flavour as the PDG of the elementon (quark/lepton/photon) it
is associated to.

No other information is filled.

<u>  FinalColourNeutrals (ReconstructedParticle):</u>

`getEnergy`, `getMass`, `getMomentum`,  returns those quantities, as given by the values
of the string, or the particle produced by the cluster, or as the sum of the two
leptons from the initial boson, or as the sum of true energies of all stable particles
(for jet type 1,3,2,5), or simply as the true value of each photon (type 4 or 6)
`getParticles` returns the TrueJets connected to the system (always two, except for ISR)

   `getParticleIDs()[0]->getType` returns the  Colour Neutral type as:
   
       1  : c.n. is string
       2  : c.n. is lepton-pair
       3  : c.n. is cluster
       4  : c.n. is single ISR
       6  : c.n. is single photon from the M.E.

`getParticleIDs()[1:2]->getType` is the same as the type of jet 1:2 emerging from the
Colour Neutral system.

   `getParticleIDs()[0:2]->getPDG` returns the  Colour Neutral PDG as:
   
       0 : PDG of the system itself (92=string, 91=cluster, 22=photon (ISR or ME), any lepton PDG=lepton pair)
       1 : PDG of first elementon (quark, lepton or photon)
       2 : PDG of second elementon (quark or  lepton)

No other information is filled.

<u>  InitialColourNeutrals (ReconstructedParticle):</u>

`getEnergy`, `getMass`, `getMomentum`, returns those quantities, as given by the values
of the CMcluster, or the particle produced by the boson, in case there was no CMcluster
(sometimes happens for quarks, almost always for leptons)
`getParticles` returns the _TrueJets_ connected to the system (at least two, possibly more
if there was hard gluon radiation or (top) quark-decay during the parton shower)

   `getParticleIDs()[0]->getType` returns the Colour Neutral type as:
   
       1,3  : c.n. is hadronic
       2    : c.n. lepton-pair
       
`getParticleIDs()[1:n]->getType` is the same as the type of jet 1:n emerging from the
 Colour Neutral system.
   
   `getParticleIDs()[0:n]->getPDG` returns the  Colour Neutral PDG as:

       0   : PDG of the boson (23=Z, 24=W, 25=H ... )
       1:n : PDG of elementon (quark, lepton) 1:n

Note that the PDG of the boson is the _last_ boson in the case of a chain, specifically for
H->WW*, ZZ* or Zgamma, the boson is *not* the Higgs, but rather the W, Z or gamma. In H->ffbar, on
the other hand, the boson *is* the Higgs.

No other information is filled.
 
** Navigators: **

|   Name                      |  Object   |  Related objects |
| ----------                  | --------- | --------- |
|  TrueJetPFOLink             | TrueJet-> | all PFO of jet |
|  TrueJetMCParticleLink      | TrueJet-> | all MCParticles of jet |
|  FinalElementonLink         | TrueJet-> | the MCParticle that is the Final elementon (quark/lepton/photon) of the jet  |
|  InitialElementonLink       | TrueJet-> | the MCParticle that is the Initial elementon (quark/lepton/photon) of the jet  |
|  FinalColourNeutralLink     | TrueJet-> | FinalColourNeutrals of the jet |
|  InitialColourNeutralLink   | TrueJet-> | InitialColourNeutrals of the jet |


To make this clear: If `reltjmcp` is created by

```
    tjmcplcol  = evt->getCollection(  _trueJetMCParticleLink );
    reltjmcp = new LCRelationNavigator( tjmcplcol );
```

then


`reltjmcp->getRelated`*To*`Objects( aTrueJet )`
                     

returns the list of mcparticles in jet aTrueJet, while


`reltjmcp->getRelated`*From*`Objects( anMCparticle )`

returns the jet an MCParticle belongs to (It is in principle a list, but off course
it only contains one element!)

Differently put, the order of the elements in the table is the same as the order of
the arguments of the addRelation method of an `LCRelationNavigator`.

Except for _TrueJetMCParticle_, the weight of the relation has no particular meaning.
In _TrueJetMCParticle_ the meaning of the weight (as from `www =  reltjmcp->getRelated`*To*`Weights(  aTrueJet )` )
is:

      1   : stable particle from generator, or decayed particle where the sum of 4-mom of the daughters is 
            different from the 4-mom of the particle (which indicates that the detector simulation has modified
	    the particle before decay)  
     -1   : generator inner brems photon.

      11  : generator decayed particle.
      0   : created in simulation, or generator particle where the parent has been modified before decay,
            so that the 4-mom is inconsistent.

Hence, adding all MCP:s related with any jets of types 1,2,3 or 6 with `|weight| == 1` gives the total
4-mom of the hard reaction. Also adding jets of type 4 adds the ISR to the sum and finally
including type 5 jets also include the overlay events.
