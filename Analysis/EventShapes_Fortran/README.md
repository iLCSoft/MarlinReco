# YThres

The YThresh module calculates the yThresh variable, which is the crossover value of the yCut jet finding
variable from NMin to NMin+1 jets found using durhamycut.  For example, if NMin=2 the yThresh variable is
the value of yCut above which durhamycut returns 2 jets, below which it returns 3 jets. The value of
yThresh will be stored as a parameter in the ReconstructedParticle collection with name y[NMin][NMin+1],
ie y23 for NMin=2. SatoruJetFinder must be installed for this package to run.

For questions/comments please email Ben Hooberman at benhooberman@berkeley.edu

## Processor Parameters

```
string RecoParticleCollection: Name of the input ReconstructedParticle collection.
int NMin: Min number of jets. For example, NMin=2 will give y23, the crossover value from 2 to 3 jets.
int PrintOutput: toggle text output stating event number and value of yThresh.
int NIterations: the number of times the durhamycut algorithm is performed. The more iterations, the lower the error in yThresh.
int NMinParticles: minimum number of particles for yThresh calculation. This should usually remain at 3.
int YStart: the starting value of yCut for the yThresh calculation.  This value should be significantly
	larger than typical values of yThresh.
```

## XML

The following should be inserted into your xml file and called after the ReconstructedParticle collection is made.

```
<processor name="MyYThresh" type="YThresh">
  <!--Name of the ReconstructedParticle collection-->
  <parameter name="RecoParticleCollection" type="string">  RecoParticles  </parameter>
  <!--Min number of jets. ie. NMin=2 will return y23-->
  <parameter name="NMin" type="int">  2  </parameter>
  <!--Toggle text output of yThresh value found==>
  <parameter name="PrintOutput" type="int">  1  </parameter>
  <!--Number of times to perform durhamycut algorithm-->
  <parameter name="NIterations" type="int">  20  </parameter>
  <!--Starting value of yCut.  Should be larger than typical yThresh values-->
  <parameter name="YStart" type="float"> 0.1  </parameter>
</processor>
```
