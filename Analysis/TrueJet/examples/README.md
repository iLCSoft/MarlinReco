##TrueJet_Parser:

class to help to interpret the collections created by `TrueJet`.
 
>  **Usage**:
>
>   Compile as any marlin-processor, add the created shared library to MARLIN_DLL.


In the processor using it, three things are needed:

  - Make the class inherit TrueJet_Parser
  - Add the definitions of the collections names of TrueJet to the c'tor of your processor.
  - Implement the method get_recoMCTruthLink.
  - make sure TrueJet_Parser.h is included (and is found in the path)


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

(This assumes that your processor already has the name of the _RecoMCTruthLink_
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

The you are ready to use the methods of TrueJet_Parser:
somewhere at the beginning of processEvent(), add

```

    // tj is a pointer to a Trujet_Parser, with the data of this processor object:
    TrueJet_Parser* tj= this ;
    // this method gets all the collections needed + initialises a few convenient variables.
    tj->getall(evt);

```

That's it !


For a working example, look in examples/Use_TrueJet !