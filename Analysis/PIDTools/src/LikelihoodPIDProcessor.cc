#include <vector>
#include <string>

#include <EVENT/LCCollection.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/PIDHandler.h>

#include "TLorentzVector.h"

#include "LikelihoodPIDProcessor.hh"
#include "LikelihoodPID.hh"
#include "LowMomentumMuPiSeparationPID_BDTG.hh"

LikelihoodPIDProcessor aLikelihoodPIDProcessor ;

LikelihoodPIDProcessor::LikelihoodPIDProcessor()
  : Processor("LikelihoodPIDProcessor") {
  
  // Processor description
  _description = "Particle ID code using Bayesian Classifier" ;
 
  std::vector< std::string > xmlfiles;
  xmlfiles.push_back( " /weightFiles_forLowMomentumMuPiSeparation/TMVAClassification_BDTG_02GeVP_clusterinfo.weights.xml" );  
  xmlfiles.push_back( "TMVAClassification_BDTG_03GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_04GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_05GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_06GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_07GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_08GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_09GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_10GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_11GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_12GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_13GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_14GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_15GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_16GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_17GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_18GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_19GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_20GeVP_clusterinfo.weights.xml" );
  
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
  
  registerProcessorParameter( "FileWeightFormupiSeparationName" ,
  			      "weight file for low momentum mu pi separation",
  			      _weightFileName,
  			      xmlfiles );
}  
  


void LikelihoodPIDProcessor::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  
  _pdgTable.push_back(11);
  _pdgTable.push_back(13);
  _pdgTable.push_back(211);
  _pdgTable.push_back(321);
  _pdgTable.push_back(2212);
  
  //for parameters except for dEdxPID
  _particleNames.push_back("electronLikelihood");
  _particleNames.push_back("muonLikelihood");
  _particleNames.push_back("pionLikelihood");
  _particleNames.push_back("kaonLikelihood");
  _particleNames.push_back("protonLikelihood");
  _particleNames.push_back("MVAOutput_mupiSeparation");
  _particleNames.push_back("electronProbability");
  _particleNames.push_back("muonProbability");
  _particleNames.push_back("pionProbability");
  _particleNames.push_back("kaonProbability");
  _particleNames.push_back("protonProbability");
 
  //for parameters of dEdxPID
  _dEdxNames.push_back("electronLikelihood");
  _dEdxNames.push_back("muonLikelihood");
  _dEdxNames.push_back("pionLikelihood");
  _dEdxNames.push_back("kaonLikelihood");
  _dEdxNames.push_back("protonLikelihood");
  _dEdxNames.push_back("MVAOutput_mupiSeparation");
  _dEdxNames.push_back("electronProbability");
  _dEdxNames.push_back("muonProbability");
  _dEdxNames.push_back("pionProbability");
  _dEdxNames.push_back("kaonProbability");
  _dEdxNames.push_back("protonProbability");
  _dEdxNames.push_back("electron_dEdxdistance");
  _dEdxNames.push_back("muon_dEdxdistance");
  _dEdxNames.push_back("pion_dEdxdistance");
  _dEdxNames.push_back("kaon_dEdxdistance");
  _dEdxNames.push_back("proton_dEdxdistance");

  _myPID = new LikelihoodPID(_PDFName); 

  //mupi separation class
  //_mupi = new Class(_weightFileName);
    _mupiPID = new LowMomentumMuPiSeparationPID_BDTG(_weightFileName);

  printParameters();
  
}

void LikelihoodPIDProcessor::processRunHeader( LCRunHeader* run) { 
} 

void LikelihoodPIDProcessor::processEvent( LCEvent * evt ) { 
  _pfoCol = evt->getCollection( _inputPFOsCollection ) ;
  
  int npfo = _pfoCol->getNumberOfElements();
  PIDHandler pidh(_pfoCol);   //BasicPID
  int algoID1 = pidh.addAlgorithm("BasicVariablePID", _particleNames);
  int algoID2 = pidh.addAlgorithm("dEdxPID", _dEdxNames);
  int algoID3 = pidh.addAlgorithm("ShowerShapesPID", _particleNames);
  int algoID4 = pidh.addAlgorithm("LikelihoodPID", _particleNames);
  for (int i = 0; i < npfo; i++ ) {
    ReconstructedParticleImpl* part = dynamic_cast<ReconstructedParticleImpl*>( _pfoCol->getElementAt(i) );
    
    if(part->getCharge()==0) continue;  //avoid neutral particle
    
    EVENT::ClusterVec clu=part->getClusters();
    lcio::Track* trk = part->getTracks()[0];
    TLorentzVector pp(part->getMomentum()[0],
		      part->getMomentum()[1],
		      part->getMomentum()[2],
		      part->getEnergy());
    
    Int_t parttype=-1;

    //several partivle IDs performed
    //use just basic variables
    _myPID->setBasicFlg(true);
    _myPID->setdEdxFlg(false);
    _myPID->setShowerShapesFlg(false);
    parttype = _myPID->Classification(pp, trk, clu);
    if(parttype<0) parttype=2;

    //mu-pi Separation for very low momentum tracks (from 0.2 GeV until 2 GeV)
    Float_t MVAoutput = -1.0;
    if((parttype==1 || parttype==2) && pp.P()<2.0){
        parttype=_mupiPID->MuPiSeparation(pp, trk, clu);
        MVAoutput = _mupiPID->getMVAOutput();   
    }
    
    //create PIDHandler
    createParticleIDClass(parttype, part, pidh, algoID1, MVAoutput);

    //use just dEdx variables
    _myPID->setBasicFlg(false);
    _myPID->setdEdxFlg(true);
    _myPID->setShowerShapesFlg(false);
    parttype = _myPID->Classification(pp, trk, clu);
    if(parttype<0) parttype=2;

     //mu-pi Separation for very low momentum tracks (from 0.2 GeV until 2 GeV)
    if((parttype==1 || parttype==2) && pp.P()<2.0){
        parttype=_mupiPID->MuPiSeparation(pp, trk, clu);
        MVAoutput = _mupiPID->getMVAOutput();   
    }
    

    //create PIDHandler
    createParticleIDClass(parttype, part, pidh, algoID2, MVAoutput);

    //use just Shower Profile variables
    _myPID->setBasicFlg(false);
    _myPID->setdEdxFlg(false);
    _myPID->setShowerShapesFlg(true);
    parttype = _myPID->Classification(pp, trk, clu);
    if(parttype<0) parttype=2;

     //mu-pi Separation for very low momentum tracks (from 0.2 GeV until 2 GeV)
    if((parttype==1 || parttype==2) && pp.P()<2.0){
        parttype=_mupiPID->MuPiSeparation(pp, trk, clu);
        MVAoutput = _mupiPID->getMVAOutput();   
    }
    

    //create PIDHandler
    createParticleIDClass(parttype, part, pidh, algoID3, MVAoutput);

    //calculate global likelihood
    _myPID->setBasicFlg(true);
    _myPID->setdEdxFlg(true);
    _myPID->setShowerShapesFlg(true);
    parttype = _myPID->Classification(pp, trk, clu);
    if(parttype<0) parttype=2;

    //mu-pi Separation for very low momentum tracks (from 0.2 GeV until 2 GeV)
    if((parttype==1 || parttype==2) && pp.P()<2.0){
        parttype=_mupiPID->MuPiSeparation(pp, trk, clu);
        MVAoutput = _mupiPID->getMVAOutput();   
    }
    
    //create PIDHandler
    createParticleIDClass(parttype, part, pidh, algoID4, MVAoutput);

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

void LikelihoodPIDProcessor::createParticleIDClass(int parttype, ReconstructedParticle *part, PIDHandler &pidh, int algoID, float MVAoutput){

  Double_t *posterior = _myPID->GetPosterior();
  Double_t *likelihood = _myPID->GetLikelihood();
  std::vector<float> likelihoodProb;
  for(int j=0;j<5;j++) likelihoodProb.push_back(likelihood[j]);
  likelihoodProb.push_back(MVAoutput);
  for(int j=0;j<5;j++) likelihoodProb.push_back(posterior[j]);

  //std::cout << "check likelihood: " << parttype << " " << algoID << " "
  //	    << likelihoodProb[0] << " " << likelihoodProb[1] << " " << likelihoodProb[2] << " " << likelihoodProb[3] << " " << likelihoodProb[4] << std::endl;
  
  //for dEdx PID
  if(pidh.getAlgorithmName(algoID)=="dEdxPID"){
    likelihoodProb.push_back((float)_myPID->get_dEdxDist(0));  //electron hypothesis
    likelihoodProb.push_back((float)_myPID->get_dEdxDist(1));  //muon hypothesis
    likelihoodProb.push_back((float)_myPID->get_dEdxDist(2));  //pion hypothesis
    likelihoodProb.push_back((float)_myPID->get_dEdxDist(3));  //kaon hypothesis
    likelihoodProb.push_back((float)_myPID->get_dEdxDist(4));  //proton hypothesis

    //std::cout << "check dedx: " << parttype << " " << algoID << " "
    //	      << likelihoodProb[11] << " " << likelihoodProb[12] << " " << likelihoodProb[13] << " " << likelihoodProb[14] << " " << likelihoodProb[15] << std::endl;
  }

  //std::cout << "check posterior: " << posterior[0] << " " 
  //	    << posterior[1] << " "
  //	    << posterior[2] << " "
  //	    << posterior[3] << " "
  //	    << posterior[4] << " " << std::endl;
  
  //set pid results
  pidh.setParticleID(part, 0, _pdgTable[parttype], (float)likelihood[parttype], algoID, likelihoodProb);
  
  return;
}
