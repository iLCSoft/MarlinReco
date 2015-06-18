// *****************************************************
// Add a new LCCollection to the event with with a copy
// of PandoraPFOs including the covariance matrix
// in the momenta space for those particles with
// charge non null.
//------------------------------
// TODO: 
// The method should be implemented on the
// ReconstructedParticle's from PandoraPFO collection,
// without adding a new collection to the event.
//
//
// *****************************************************
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <marlin/Exceptions.h>

#include "momentumCov.h"

MomentumCov aMomentumCov ;


MomentumCov::MomentumCov() : Processor("MomentumCov") {

  // modify processor description
  _description = "MomentumCov calculate the covariance";
  _description += " matrix in the momenta space." ;


  // register steering parameters: name, description,
  // class-variable, default value
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                          "InputPandoraPFOsCollection" ,
                          "Name of the PandoraPFOs collection"  ,
                          _colPFOs ,
                          std::string("PandoraPFOs")) ;

  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                           "OutputNewPFOsCollection",
                           std::string("Name of the new PFOs collection")+=
                           " with covariance matrix",
                           _colNewPFOs,
                           std::string("Newpfos"));

}

void MomentumCov::init() {

  streamlog_out(DEBUG) << "   init called  "
                       << std::endl ;

  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
}

void MomentumCov::processRunHeader(LCRunHeader* run) {

  _nRun++ ;
}

void MomentumCov::processEvent(LCEvent * evt) {


  // -- Read out PFO information --
  LCCollection *colPFO = evt->getCollection(_colPFOs);
  if (!colPFO) {
      streamlog_out(ERROR) << "No PFO Collection Found!" << std::endl;
      throw marlin::SkipEventException(this);
  }

  LCCollectionVec *pNewPFOsCollection =
    new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  // pNewPFOsCollection->setTransient();
  // pNewPFOsCollection->setDefault();
  // std::cerr << " FLAG " << pNewPFOsCollection->getFlag() << std::endl;

  for (int i=0;i<colPFO->getNumberOfElements();++i) {

      ReconstructedParticle  *recPar =
        dynamic_cast<ReconstructedParticle  *>(colPFO->getElementAt(i));

      ReconstructedParticleImpl *recParImpl = new ReconstructedParticleImpl();

      copy_reconstructedParticle(recPar, recParImpl);

      // lower triangule covariance matrix elements in momenta space
      // (px,py,pz,E): total 10 elements
      FloatVec covarianceP(10, 0.0);
      // check if charged particle. In that case,
      // calculate the covariance matrix and stored it.
      if (TMath::Abs(recPar->getCharge()) > 0.5) {

          getCovMatrixMomenta(recPar, covarianceP);
          recParImpl->setCovMatrix (covarianceP);
      }
      pNewPFOsCollection->addElement(recParImpl);
  }
  // add new collection to the event
  evt->addCollection(pNewPFOsCollection,_colNewPFOs.c_str());

  _nEvt++;
  
}


void MomentumCov::check(LCEvent * evt) {
  // nothing to check here - could be used to fill
  // checkplots in reconstruction processor
}


void MomentumCov::end(){


  streamlog_out(MESSAGE4) << " processed "               << _nEvt
                          << " events in "               << _nRun
                          << " runs "                    << std::endl ;
}



void MomentumCov::copy_reconstructedParticle(ReconstructedParticle const *rec_orig,
                                             ReconstructedParticleImpl *rec_new){
// Copy all members from pfo from PandoraPFOs
// collection to a new ReconstructedParticle.
// I am not able to modify directly the PandoraPFOs: elements are read-only.

  rec_new->setType     (rec_orig->getType());
  rec_new->setMomentum (rec_orig->getMomentum());
  rec_new->setEnergy   (rec_orig->getEnergy());
  rec_new->setCovMatrix(rec_orig->getCovMatrix());
  rec_new->setMass     (rec_orig->getMass());
  rec_new->setCharge   (rec_orig->getCharge());

  rec_new->setReferencePoint (rec_orig->getReferencePoint()) ;
  ParticleIDVec::const_iterator it1 = rec_orig->getParticleIDs().begin();
  for (;it1!=rec_orig->getParticleIDs().end();++it1) {
      rec_new->addParticleID(*it1) ;
  }

  rec_new->setParticleIDUsed (rec_orig->getParticleIDUsed()) ;
  rec_new->setGoodnessOfPID  (rec_orig->getGoodnessOfPID()) ;

  ReconstructedParticleVec::const_iterator it2 = rec_orig->getParticles().begin();
  for (;it2!=rec_orig->getParticles().end();++it2) {
      rec_new->addParticle(*it2) ;
  }

  ClusterVec::const_iterator it3 = rec_orig->getClusters().begin();
  for (;it3!=rec_orig->getClusters().end();++it3) {
      rec_new->addCluster(*it3) ;
  }

  TrackVec::const_iterator it4 = rec_orig->getTracks().begin();
  for (;it4!=rec_orig->getTracks().end();++it4) {
      rec_new->addTrack(*it4) ;
  }

  rec_new->setStartVertex    (rec_orig->getStartVertex()) ;
}

