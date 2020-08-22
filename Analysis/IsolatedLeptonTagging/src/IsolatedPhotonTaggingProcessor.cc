// *****************************************************
// Processor for isolated photon tagging
//                        ----originally developped by S.Kawada and J.Tian
//                        ----for miniDST analysis
// *****************************************************
#include "IsolatedPhotonTaggingProcessor.h"

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Cluster.h>
#include <UTIL/LCTypedVector.h>
#include <EVENT/Track.h>
#include <marlin/Exceptions.h>
#include <EVENT/Vertex.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

//root
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMVA/Reader.h"
//local
#include "Utilities.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;

using namespace TMVA;
using namespace isolep;

IsolatedPhotonTaggingProcessor aIsolatedPhotonTaggingProcessor ;


IsolatedPhotonTaggingProcessor::IsolatedPhotonTaggingProcessor() : Processor("IsolatedPhotonTaggingProcessor") {
  
  // modify processor description
  _description = "IsolatedPhotonTaggingProcessor does whatever it does ..." ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputPandoraPFOsCollection" , 
			   "Name of the PandoraPFOs collection"  ,
			   _colPFOs ,
			   std::string("PandoraPFOs") ) ;

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			    "OutputPFOsWithoutIsoPhotonCollection",
			    "Name of the new PFOs collection without isolated photon",
			    _colNewPFOs,
			    std::string("PandoraPFOsWithoutIsoPhoton") );
  
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			    "OutputIsoPhotonsCollection",
			    "Name of collection with the selected isolated photon",
			    _colPhotons,
			    std::string("ISOPhotons") );
  
  registerProcessorParameter("DirOfISOPhotonWeights",
			     "Directory of Weights for the Isolated Photon MVA Classification"  ,
			     _isolated_photon_weights ,
			     std::string("isolated_photon_weights") ) ;
  
  registerProcessorParameter("CutOnTheISOPhotonMVA",
			     "Cut on the mva output of isolated photon selection"  ,
			     _mvaCut,
			     float(0.5) ) ;
  
  registerProcessorParameter("IsSelectingOneIsoPhoton",
			     "flag to select one most like isolated photon"  ,
			     _is_one_isophoton ,
			     bool(false) ) ;
  
  registerProcessorParameter("MinEForPhoton",
			     "Minimum energy for photon"  ,
			     _minE,
			     float(5.) ) ;
  
  registerProcessorParameter("CosConeSmall",
			     "cosine of the smaller cone"  ,
			     _cosConeSmall ,
			     float(0.98) ) ;
  
  registerProcessorParameter("CosConeLarge",
			     "cosine of the larger cone"  ,
			     _cosConeLarge ,
			     float(0.95) ) ;

  registerProcessorParameter("RatioNeutralConeEnergy",
			     "Cut on the ratio of Neutral Cone Energy over Photon Energy"  ,
			     _coneNeutralRatio,
			     float(0.2) ) ;

  registerProcessorParameter("RatioChargedConeEnergy",
			     "Cut on the ratio of Charged Cone Energy over Photon Energy"  ,
			     _coneChargedRatio,
			     float(0.1) ) ;
}

void IsolatedPhotonTaggingProcessor::init() { 

  streamlog_out(DEBUG) << "IsolatedPhotonTagging   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;


}

void IsolatedPhotonTaggingProcessor::processRunHeader( LCRunHeader*  /*run*/) { 
} 

void IsolatedPhotonTaggingProcessor::processEvent( LCEvent * evt ) { 

    
  streamlog_out(DEBUG) << "Hello, Isolated Photon Tagging!" << endl;

  // -- get PFO collection --
  LCCollection *colPFO = nullptr;
  try {
    colPFO = evt->getCollection(_colPFOs);
  }
  catch(DataNotAvailableException &e) {
    std::cerr << "No PFO collection found (" << _colPFOs << ") !" << std::endl;
    return;
  }

  Int_t nPFOs = colPFO->getNumberOfElements();
  std::vector<lcio::ReconstructedParticle*> newPFOs;
  std::vector<lcio::ReconstructedParticle*> isoPhotons;
  FloatVec isoPhotonTagging;

  // loop all the PFOs
  float energy_photon_max = -1.;
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    // here in principle can add any precut on each PFO
    newPFOs.push_back(recPart);
    Double_t energy = recPart->getEnergy();
    Double_t charge = recPart->getCharge();
    TVector3 momentum = TVector3(recPart->getMomentum());
    TLorentzVector lortz = TLorentzVector(momentum,energy);
    Double_t momentumMagnitude = momentum.Mag();
    // get cone information
    // double cone based isolation algorithm
    Double_t coneEnergy0[3] = {0.,0.,0.};  // {total, neutral, charged} cone energy
    Double_t pLargeCone[4]  = {0.,0.,0.,0.}; // 4-momentum of all particles inside the larger cone
    getConeEnergy(recPart,colPFO,_cosConeSmall,coneEnergy0,_cosConeLarge,pLargeCone);
    Double_t coneEN     = coneEnergy0[1];
    Double_t coneEC     = coneEnergy0[2];
    TLorentzVector lortzLargeCone = TLorentzVector(pLargeCone[0],pLargeCone[1],pLargeCone[2],pLargeCone[3]);
    TVector3 momentumLargeCone = lortzLargeCone.Vect();
    Double_t cosThetaWithLargeCone = 1.;
    if (momentumLargeCone.Mag() > 0.0000001) {
      cosThetaWithLargeCone = momentum.Dot(momentumLargeCone)/momentumMagnitude/momentumLargeCone.Mag();
    }
    Double_t energyRatioWithLargeCone = energy/(energy+lortzLargeCone.E());
    Double_t pandoraID = recPart->getType();
    if (TMath::Abs(charge) < 0.5 && pandoraID == 22) {
      if (energy < _minE) continue; // cut on the minimum energy
      // default option to select all isolated photons
      if (coneEC/energy < _coneChargedRatio && coneEN/energy < _coneNeutralRatio) {  // simple charged and neutral cone energy ratio cut
	if (!_is_one_isophoton) {
	  isoPhotons.push_back(recPart);
	}
	else {
	  if (energy > energy_photon_max) {
	    energy_photon_max = energy;
	    isoPhotons.clear();
	    isoPhotons.push_back(recPart);
	  }
	}
      }
    }
  }

  LCCollectionVec *pPFOsWithoutIsoPhotonCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec *pIsoPhotonCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  pPFOsWithoutIsoPhotonCollection->setSubset(true);
  pIsoPhotonCollection->setSubset(true);
  // save the selected photon to a new collection
  for (std::vector<lcio::ReconstructedParticle*>::const_iterator iObj=isoPhotons.begin();iObj<isoPhotons.end();++iObj) {
    pIsoPhotonCollection->addElement(*iObj);
  }
  // save other PFOs to a new collection
  for (std::vector<lcio::ReconstructedParticle*>::const_iterator iObj=newPFOs.begin();iObj<newPFOs.end();++iObj) {
    Bool_t isLep=kFALSE;
    for (std::vector<lcio::ReconstructedParticle*>::const_iterator iLep=isoPhotons.begin();iLep<isoPhotons.end();++iLep) {
      if ((*iObj) == (*iLep)) isLep = kTRUE;
    }
    if (!isLep) pPFOsWithoutIsoPhotonCollection->addElement(*iObj);
  }

  // add new collections
  evt->addCollection(pPFOsWithoutIsoPhotonCollection,_colNewPFOs.c_str());
  evt->addCollection(pIsoPhotonCollection,_colPhotons.c_str());


}



void IsolatedPhotonTaggingProcessor::check( LCEvent *  /*evt*/ ) { 
}


void IsolatedPhotonTaggingProcessor::end(){ 

  for (std::vector<TMVA::Reader*>::const_iterator ireader=_readers.begin();ireader!=_readers.end();++ireader) {
    delete *ireader;
  }

}

