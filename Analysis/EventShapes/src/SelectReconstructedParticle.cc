#include "SelectReconstructedParticle.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/LCIO.h>
#include <Exceptions.h>
// #include <EVENT/ReconstructedParticle.h>

// #include <EVENT/MCParticle.h>
// #include <IMPL/LCRunHeaderImpl.h>
// #include <IMPL/LCEventImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>


using namespace lcio ;
using namespace marlin ;
using namespace std ;


SelectReconstructedParticle aSelectReconstructedParticle ;


SelectReconstructedParticle::SelectReconstructedParticle() 
   : Processor("SelectReconstructedParticle") {
  
  // modify processor description
  _description = "SelectReconstructedParticle: Selects particles from all reconstructed particles to be used for the thrust finder" ;
  

  // register steering parameters: 
  // name, description, class-variable, default value

  registerProcessorParameter( "inputCollectionName" ,
      "Collection of reconstructed particles to chose from"  ,
      _inputCollectionName ,
      std::string(LCIO::RECONSTRUCTEDPARTICLE) ) ;
  registerProcessorParameter( "outputCollectionName" ,
      "Collection of selected reconstructed particles"  ,
      _outputCollectionName ,
      std::string("SelectedReconstructedParticle") ) ;
  registerProcessorParameter(
	"MinimumMomentum" ,
        "Minimum momentum a particle has to have to be used for the thrust calculation"  ,
	_minimumMomentum ,
	float(0) ) ;

}

void SelectReconstructedParticle::init() { 

  // usually a good idea to
  printParameters() ;

}

void SelectReconstructedParticle::processRunHeader( LCRunHeader* run) { 
  //   cout << "processing runheader" << endl;
} 

void SelectReconstructedParticle::processEvent( LCEvent * evt ) { 
  int good = 0;
  std::stringstream errorMsg;
  // cout << "processing event: " << endl;

  // get pointer to collection vec of input particles  
  LCCollection* inParVec = evt->getCollection(_inputCollectionName) ;
  if (inParVec->getTypeName()!=LCIO::RECONSTRUCTEDPARTICLE) {
    errorMsg << "\nProcessor: SelectReconstructedParticle \n" << 
      "Collection is of wrong type (" << inParVec->getTypeName() <<
      "). Processor requires collection tpye " << LCIO::RECONSTRUCTEDPARTICLE << "\n" ; 
    throw Exception(errorMsg.str()); 
  }

  // create collection of output particles 
  LCCollectionVec* outParVec = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE)  ;
  outParVec->setSubset(true);


  // loop over all Particles an do cuts 
  for (int i = 0; i < inParVec->getNumberOfElements() ; i++)
    {
      float totparmom;
      const double *pparmom;
      float x;
      pparmom = dynamic_cast<ReconstructedParticle*>(inParVec->getElementAt(i))->getMomentum();


      // (inParVec->getElementAt(i))->getMomentum();
      totparmom = sqrt( pparmom[0]*pparmom[0]
                      + pparmom[1]*pparmom[1]
                      + pparmom[2]*pparmom[2]); 
      if (totparmom >= _minimumMomentum &&
          totparmom <= 5000 &&
          !dynamic_cast<ReconstructedParticle*>(inParVec->getElementAt(i))->isCompound()) { 
	good++;
	if (good <3000) 
	outParVec->addElement(inParVec->getElementAt(i));
      }
    }

  // add collection to event
  evt->addCollection( outParVec, _outputCollectionName )  ;

  /* cout << "run: " << evt->getRunNumber() 
   << " Event no.: " << evt->getEventNumber() <<
    " Teilchen: " << inParVec->getNumberOfElements()
    << outParVec->getNumberOfElements() << endl; */ 

}

void SelectReconstructedParticle::end(){ 
  
//   std::cout << "MyProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}
