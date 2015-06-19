#include <vector>
#include <string>

#include <EVENT/LCCollection.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/PIDHandler.h>

#include "TLorentzVector.h"

#include "LikelihoodPIDProcessor.hh"
#include "LikelihoodPID.hh"

LikelihoodPIDProcessor aLikelihoodPIDProcessor ;

LikelihoodPIDProcessor::LikelihoodPIDProcessor()
  : Processor("LikelihoodPIDProcessor") {
  
  // Processor description
  _description = "Particle ID code using Bayesian Classifier" ;
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "RecoParticleCollection" ,
			   "Input collection of Reconstrcuted Particle",
			   _inputPFOsCollection,
			   std::string("PandoraPFOs"));
  
  registerProcessorParameter( "EnergyBoundaries" ,
			      "Boundaries for energy",
			      _energyBoundary,
			      EVENT::FloatVec(0,1.0e+07));
  
  registerProcessorParameter( "FilePDFName" ,
			      "rootfile of PDF",
			      _PDFName,
			      std::string("pdf_ParticleID_ok.root") );
} 

void LikelihoodPIDProcessor::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  
  _pdgTable.push_back(11);
  _pdgTable.push_back(13);
  _pdgTable.push_back(211);
  _pdgTable.push_back(321);
  _pdgTable.push_back(2212);
  
  _particleNames.push_back("electron");
  _particleNames.push_back("muon");
  _particleNames.push_back("pion");
  _particleNames.push_back("kaon");
  _particleNames.push_back("proton");

  _myPID = new LikelihoodPID(_PDFName); 
  
  printParameters();
  
}

void LikelihoodPIDProcessor::processRunHeader( LCRunHeader* run) { 
} 

void LikelihoodPIDProcessor::processEvent( LCEvent * evt ) { 
  _pfoCol = evt->getCollection( _inputPFOsCollection ) ;

  
  int npfo = _pfoCol->getNumberOfElements();
  PIDHandler pidh(_pfoCol);
  int algoID = pidh.addAlgorithm("LikelihoodPID", _particleNames);
  for (int i = 0; i < npfo; i++ ) {
    ReconstructedParticleImpl* part = dynamic_cast<ReconstructedParticleImpl*>( _pfoCol->getElementAt(i) );
    
    if(part->getCharge()==0) continue;  //avoid neutral particle
    
    EVENT::ClusterVec clu=part->getClusters();
    lcio::Track* trk = part->getTracks()[0];
    TLorentzVector pp(part->getMomentum()[0],
		      part->getMomentum()[1],
		      part->getMomentum()[2],
		      part->getEnergy());
    
    Int_t parttype=_myPID->Classification(pp, trk, clu);
    if(parttype<0) parttype=2;
    Double_t *posterior=_myPID->GetPosterior();
    std::vector<float> likelihoodProb;
    for(int j=0;j<5;j++) likelihoodProb.push_back(posterior[j]);

    //std::cout << "check posterior: " << posterior[0] << " " 
    //	    << posterior[1] << " "
    //	    << posterior[2] << " "
    //	    << posterior[3] << " "
    //	    << posterior[4] << " " << std::endl;

    //set pid results
    pidh.setParticleID(part, 0, _pdgTable[parttype], (float)posterior[parttype], algoID, likelihoodProb);

    /*ParticleIDImpl *PIDImpl=new ParticleIDImpl();
    if(parttype==0){
      PIDImpl->setType(11);
      PIDImpl->setLikelihood((float) posterior[0]);
    }
    if(parttype==1){
      PIDImpl->setType(13);
      PIDImpl->setLikelihood((float) posterior[1]);
    }
    if(parttype==2){
      PIDImpl->setType(211); 
      PIDImpl->setLikelihood((float) posterior[2]);
    }
    if(parttype==3){
      PIDImpl->setType(321);
      PIDImpl->setLikelihood((float) posterior[3]);
    }
    if(parttype==4){
      PIDImpl->setType(2212);
      PIDImpl->setLikelihood((float) posterior[4]);
    }
    
    //posterior probability
    PIDImpl->addParameter((float) posterior[0]);   //electron
    PIDImpl->addParameter((float) posterior[1]);   //muon
    PIDImpl->addParameter((float) posterior[2]);   //pion
    PIDImpl->addParameter((float) posterior[3]);   //kaon
    PIDImpl->addParameter((float) posterior[4]);   //proton
    
    //add particle ID
    part->addParticleID(PIDImpl);
    //FloatVec chk=part->getParticleIDs()[0]->getParameters();
    //std::cout << chk[0] << std::endl;
    //add to PFO collection(so far, this is only PID)
    part->setParticleIDUsed(PIDImpl);
    part->setGoodnessOfPID(PIDImpl->getLikelihood());*/

  }
  
}

void LikelihoodPIDProcessor::check( LCEvent * evt ) { 
}

void LikelihoodPIDProcessor::end() { 
  delete _myPID;
}
