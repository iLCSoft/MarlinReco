 Comprehensive Particle Identification (CPID) Processor

 The CPIDProcessor gathers PID-related observables, trains a model to optimise PID determination and infers from the model to data.
 The processors uses modules to achieve this: InputAlgorithms to collect PID information, TrainingModels to train and infer.
 The modules are made available via dynamic loading, using a base class and a Mgr.
 They are then steered in the processor steering file.
 For this, the modules to be used need to be specified in the _inputAlgoSpecs and _trainModelSpecs, respectively.
 Specify [type]:[name] of a module, or only [type] in which case name=type.
 The processor then checks if float and/or string vectors are specified as processor parameters called [name].F and [name].S, respectively.
 The function of these parameters is specified in each individual module description.
 The modules are created and initialised in init.

 The modes in which the processor can run are steered via the _mode flags.
 It extracts the information from the PFOs via the InputAlgorithms in processEvent.
 These are put in a TTree and can then be stored in a root file.
 If specified, training is run in end. Training can be run with the previously generated root file only, i.e. without extraction.
 If specified, inference is run in processEvent. It needs extraction.
 Training and inference cannot be done simultaneously.

 The training creates a reference file in which the used observables and the used signal and background PDGs are stored, as well as the location of the model's weight files together with the corresponding momentum bracket.
 The inference reads the reference file and loads the TrainingModel based on the information therein.
