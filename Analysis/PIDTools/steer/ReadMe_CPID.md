# Comprehensive Particle Identification (CPID) Processor

The CPIDProcessor gathers PID observables, trains a model to optimise PID determination and infers from the model to data.  
The processors uses modules to achieve this: InputAlgorithms to collect PID information, TrainingModels to train and infer.  
The modules are made available via dynamic loading, using a base class and a Mgr.  
They are then steered in the processor steering file.  
For this, the modules to be used need to be specified in the _inputAlgoSpecs and _trainModelSpecs, respectively.  
Specify [type]:[name] of a module, or only [type] in which case name=type.  
The processor then checks if float and/or string vectors are specified as processor parameters called [name].F and [name].S.  
The function of these parameters is specified in each individual module description.  
The modules are created and initialised in init.

The modes in which the processor can run are steered via the _mode flags.  
It extracts the information from the PFOs via the InputAlgorithms in processEvent.  
These are put in a TTree and can then be stored in a root file.  
If specified, training is run in end. Training can be run with the previously generated root file only, i.e. without extraction.  
If specified, inference is run in processEvent. It needs extraction.  
Training and inference cannot be done simultaneously.  

The training creates a reference file in which the used observables and the used signal and background PDGs are stored, as well as the location of the model's weight files together with the corresponding momentum bin.  
The inference reads the reference file and loads the TrainingModel based on the information therein.

For each observable in each InputAlgorithm one branch is created in the TTree called [algorithm_name]_[observable_name].  
In addition, the processor also collects for each PFO: momentum, lambda angle (polar angle relative to z=0), the PDG of the likeliest connected MCParticle, the number of associated tracks as well as d0 and z0 of the first track (if any).  
The training observables can be specified. If not done so, the defaults are all observables associated with the specified InputAlgorithms + momentum and lambda.  

A number of PFO acceptance cuts can be applied, based on its properties or that of its first track.  
Signal and background PDGs can be defined. Only PFOs which originate from either of these are used.  
The consequence of belonging to either of these groups lies in the details of the individual TrainingModels.  
The PFOs are divided into momentum bins, for each of which a separate TrainingModel is created and trained/inferred from.  
If inference is specified, a plot of the confusion matrix of all signal PDGs is created in the current working directory during end.  

Collection names:  

- **_PFOColName** - Name of the PFO input collection (ReconstructedParticleImpl).  
     string, default: PandoraPFOs.
- **_RecoMCTruthLinkName** - Name of the link from PFOs to MCParticles input collection (LCRelation).  
     string, default: RecoMCTruthLink.

Mode selection:

- **_modeExtract** - Set true to extract PID observables via the specified InputAlgorithms.  
     bool, default: false.
- **_modeTrain** - Set true to train via the specified TrainingModels. Cannot be true at the same time as _modeInfer.  
     bool, default: false.
- **_modeInfer** - Set true to infer PID from the specified TrainingModels. Cannot be true at the same time as _modeTrain.  
     bool, default: false.

Other parameters:

- **_TTreeFileName** - Name of the root file in which the TTree with all observables is stored; optional output in case of extraction, otherwise necessary input.  
     string, default: TTreeFile.root.
- **_inputAlgoSpecs** - List of input algorithms; for each specify type:name or only type (then name=type).  
     string vector, default: {}.
- **_trainModelSpecs** - List of training models; for each specify type:name or only type (then name=type).  
     string vector, default: {}.
- **_reffile** - Reference file(s). If only one file but several training models are specified the reference files are auto-numbered.  
     string vector, default: {Ref.txt}.
- **_trainingObservables** - List of observables that should be used for traning. If empty, all observables from the specified algorithms + momabs + lambda are used.  
     string vector, default: {}
