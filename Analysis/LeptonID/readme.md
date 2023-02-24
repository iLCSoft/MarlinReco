# LeptonIDProcessor

This processor can be used to distinguish between electrons, muons and charged pions. This distinction also works relatively well at low momenta and in dense environments (jets). More information regarding the performance of the LeptonID and the motivation for building it can be found in [this presentation](doc/ma_reichenbach_2022-12-01.pdf).

## TL; DR
1. set which of the available variables to use in the steering file
2. set input file(s) to be used for training in the steering file
3. run the processor once in tree output mode
4. adapt and run the [training script](training/PID_tmva_multi_jet_dEdx_50.py) to create the weight file
5. set weight and input file(s) to be used for evaluation in the steering file
6. run the processor in eval mode

## Implementation and usage
The LeptonIDProcessor utilizes a multi-class BDT from [ROOT TMVA](https://root.cern/manual/tmva/) to classify given particles into either of the three categories ($e$, $\mu$, hadrons). The input variables available for the BDT consist of the currently available [cluster shapes](link). Some further variables are taken from the [`WeightedPoints3D` utility](link). Additionally, the dEdx distance to the electron curve can also be used as an input.

The BDT depends on a weight file that needs to be created first. For this, the processor offers the possibility to turn the BDT evaluation off and output a ROOT tree necessary for training. There is an additional [training script](training/PID_tmva_multi_jet_dEdx_50.py) to generate the weight file. Additionally, the training folder also contains [a Jupyter notebook](training/PID-var-comp.ipynb) that can be used to compare the distributions between the different classes for the variables in the training tree.

The weight file produced during the training then needs to be set in the steering file for the processor. Furthermore, the list of variables used in the training must be set accurately in the steering file.