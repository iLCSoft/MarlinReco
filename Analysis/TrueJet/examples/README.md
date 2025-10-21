

## Examples of TrueJet use

`TrueJet` is part of the standard `libMarlinReco.so`, in your `MARLIN_DLL` . The steering file
_truejet_run_dstout.xml_ can be used to run a minimal job creating a DST with all the collections
`TrueJet` creates. Edit it to select a useful input DST (an example is in there, usable as-is
at DESY), modify the output name if wanted, etc. . Then

```

$MARLIN/bin/Marlin truejet_run_dstout.xml

```

will create _/tmp/newhej.slcio_ which can be further analysed.

To analyse this file, an example exercising all options that `TrueJet_Parser` gives is
given here, in _src/Use_TrueJet.cc_ . (`TrueJet_Parser` is in `MarlinUtil`, so `MarlinUtil` must be
required by `cmake`.)

The following puts the example in place:

```

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -C $ILCSOFT/ILCSoft.cmake ..
make install
cd ../lib
export MARLIN_DLL=${MARLIN_DLL}:`pwd`/libUse_TrueJet.so
cd ..

```

Then

```

$MARLIN/bin/Marlin truejet_run_ana_frdst.xml

```

runs the example, reading the first 10 events of _/tmp/newhej.slcio_, at the highest
debug-level.

Alternatively, an existing DST can be analysed on the fly using the steering-file _truejet_run_ana_onthefly.xml_ .
Also here, edit it to taste, notably selecting the input to use, and then

```

$MARLIN/bin/Marlin truejet_run_ana_onthefly.xml


```
will do the job, provided you did build `Use_TrueJet` as per above.

The files _event-865-run-82501-4f_WW_semileptonic-I250018.lis_ and _comments-on-event-865-run-82501.txt_
contains such a debug output, and comments on how to understand the information given.




## Further information on using TrueJet_Parser:

 `TrueJet_Parser` is class to help to interpret the collections created by `TrueJet`.



In the processor using it, three things are needed:

  - Make the class inherit `TrueJet_Parser`
  - Add the definitions of the collections names of TrueJet to the c'tor of your processor.
  - Implement the method `get_recoMCTruthLink`.
  - make sure _TrueJet_Parser.h_ is included (and is found in the path)


Specifically:

 In the header of your processor:

```

    .
    .
    .
#include "TrueJet_Parser.h"

    .
    .
    .

class yourprocessor : public Processor , public TrueJet_Parser {

 public

    .
    .
    .
  std::string get_recoMCTruthLink(){ return _recoMCTruthLink  ; } ;
    .
    .


```

(This assumes that your processor already has the name of the `RecoMCTruthLink`
collection as an input parameter and that it is stored in the
variable `_recoMCTruthLink` . I guess that that is the case for almost
any processor used for analysis !)

 Then in the source file:

Modify the c'tor this way:

```

yourprocessor::yourprocessor() : Processor("yourprocessor") {

    .
    .
    .
 registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "TrueJets" ,
                           "Name of the TrueJetCollection input collection"  ,
                           _trueJetCollectionName ,
                           std::string("TrueJets") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "FinalColourNeutrals" ,
                           "Name of the FinalColourNeutralCollection input collection"  ,
                           _finalColourNeutralCollectionName ,
                           std::string("FinalColourNeutrals") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "InitialColourNeutrals" ,
                           "Name of the InitialColourNeutralCollection input collection"  ,
                           _initialColourNeutralCollectionName ,
                           std::string("InitialColourNeutrals") ) ;


  registerInputCollection( LCIO::LCRELATION,
                            "TrueJetPFOLink" ,
                            "Name of the TrueJetPFOLink input collection"  ,
                            _trueJetPFOLink,
                            std::string("TrueJetPFOLink") ) ;

  registerInputCollection( LCIO::LCRELATION,
                            "TrueJetMCParticleLink" ,
                            "Name of the TrueJetMCParticleLink input collection"  ,
                            _trueJetMCParticleLink,
                            std::string("TrueJetMCParticleLink") ) ;

  registerInputCollection( LCIO::LCRELATION,
                            "FinalElementonLink rueJetMCParticleLink" ,
                            "Name of the  FinalElementonLink input collection"  ,
                            _finalElementonLink,
                            std::string("FinalElementonLink") ) ;

  registerInputCollection( LCIO::LCRELATION,
                            "InitialElementonLink" ,
                            "Name of the  InitialElementonLink input collection"  ,
                            _initialElementonLink,
                            std::string("InitialElementonLink") ) ;

  registerInputCollection( LCIO::LCRELATION,
                            "FinalColourNeutralLink" ,
                            "Name of the  FinalColourNeutralLink input collection"  ,
                            _finalColourNeutralLink,
                            std::string("FinalColourNeutralLink") ) ;

  registerInputCollection( LCIO::LCRELATION,
                            "InitialColourNeutralLink" ,
                            "Name of the  InitialColourNeutralLink input collection"  ,
                            _initialColourNeutralLink,
                            std::string("InitialColourNeutralLink") ) ;


    .
    .
    .
}

```

Then you are ready to use the methods of `TrueJet_Parser`:
somewhere at the beginning of `processEvent()`, add

```

    // tj is a pointer to a Truejet_Parser, with the data of this processor object:
    TrueJet_Parser* tj= this ;
    // this method gets all the collections needed + initialises a few convenient variables.
    tj->getall(evt);

```

That's it ! Build you processor as usual, but remember to require `MarlinUtil` in the _CMakeLists.txt_ .
For a working example, look in src/Use_TrueJet.cc !

## Developing TrueJet :

If you want to play with TrueJet, _truejet_CMakeLists.txt_ is a minimal cmake configuration that
can be used. Note, however, that since TrueJet is in libMarlinReco.so, you need to either
have your own installation of MarlinReco, and re-_make install_ it to get any modifications included
into the library, or give your version of the TrueJet class a different name.

