// *****************************************************
// Processor to save particles from jets to a PFO collection
//                        ----Junping
// *****************************************************
#include "JetPFOsProcessor.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/LCTypedVector.h>
#include <marlin/Exceptions.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;

JetPFOsProcessor aJetPFOsProcessor ;


JetPFOsProcessor::JetPFOsProcessor() : Processor("JetPFOsProcessor") {
  
  // modify processor description
  _description = "JetPFOsProcessor does whatever it does ..." ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputJetCollection" , 
			   "Name of the input jet collection"  ,
			   _colJet ,
			   std::string("JetsFromKt") ) ;

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			    "OutputPFOsCollection",
			    "Name of the output PFOs collection",
			    _colPFOsFromJet,
			    std::string("PFOsFromKtJet") );

}

void JetPFOsProcessor::init() { 

  streamlog_out(DEBUG) << "   JetPFOs Processor init called  " 
		       << std::endl ;
  
}

void JetPFOsProcessor::processRunHeader( LCRunHeader* run) { 
} 

void JetPFOsProcessor::processEvent( LCEvent * evt ) { 

  LCCollection *colJet = evt->getCollection(_colJet);
  if (!colJet) {
    std::cerr << "No Input Jet Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }
  LCCollectionVec *pPFOsFromJetCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  pPFOsFromJetCollection->setSubset(true);
  int nJets = colJet->getNumberOfElements();
  for (int i=0;i<nJets;i++) {
    ReconstructedParticle *jet = dynamic_cast<ReconstructedParticle*>(colJet->getElementAt(i));
    std::vector<lcio::ReconstructedParticle*> partVec = jet->getParticles();
    for (std::vector<lcio::ReconstructedParticle*>::const_iterator iPart=partVec.begin();iPart!=partVec.end();++iPart) {
      pPFOsFromJetCollection->addElement((*iPart));
    }
  }
  evt->addCollection(pPFOsFromJetCollection,_colPFOsFromJet.c_str());

}



void JetPFOsProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void JetPFOsProcessor::end(){ 
}
