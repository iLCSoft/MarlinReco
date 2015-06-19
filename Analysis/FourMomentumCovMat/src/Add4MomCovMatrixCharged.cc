// *****************************************************
//
// Fill covariance matrix on (P,E) for a 
// charged ReconstructedParticle
//
// *****************************************************
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <marlin/Exceptions.h>

#include "Add4MomCovMatrixCharged.h"
#include "algebraImplementation.h"

#include <iomanip>

Add4MomCovMatrixCharged aAdd4MomCovMatrixCharged ;


Add4MomCovMatrixCharged::Add4MomCovMatrixCharged() : Processor("Add4MomCovMatrixCharged") {

  // modify processor description
  _description = "Set the convariance matrix in (P,E) for all charged pfos";
  _description += " in PandoraPFOs Collection";
  // register steering parameters: name, description,
  // class-variable, default value
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                          "InputPandoraPFOsCollection" ,
                          "Name of the PandoraPFOs collection"  ,
                          _colPFOs ,
                          std::string("PandoraPFOs")) ;
}

void Add4MomCovMatrixCharged::init() {

  streamlog_out(DEBUG) << "   init called  "
                       << std::endl ;

  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
}

void Add4MomCovMatrixCharged::processRunHeader(LCRunHeader* run) {

  _nRun++ ;
}

void Add4MomCovMatrixCharged::processEvent(LCEvent * evt) {


  // -- Read out PFO information --
  LCCollection *colPFO = evt->getCollection(_colPFOs);
  if (!colPFO) {
      streamlog_out(ERROR) << "No PFO Collection Found!" << std::endl;
      throw marlin::SkipEventException(this);
  }

  for (int i=0;i<colPFO->getNumberOfElements();++i) {

      ReconstructedParticleImpl  *recPar =
        dynamic_cast<ReconstructedParticleImpl  *>(colPFO->getElementAt(i));

      FloatVec covarianceP(10, 0.0);
      // Only for charged particles
      if (TMath::Abs(recPar->getCharge()) > 0.5) {
          getCovMatrixMomenta(recPar, covarianceP);
          recPar->setCovMatrix (covarianceP);
      }
      
      streamlog_out(MESSAGE3) << "Set covariance matrix for charged pfo with charge "
                              << (int) recPar->getCharge() << " and (px,py,pz,e) = ("
                              << recPar->getMomentum()[0] << ","
                              << recPar->getMomentum()[1] << ","
                              << recPar->getMomentum()[2] << ","
                              << recPar->getEnergy() << "):" << std::endl
                              << "(null charged pfo covariance matrix not implemented yet)" << std::endl
                              << std::setprecision(6)
                              << " cov1[px px] " <<  recPar->getCovMatrix()[0] << std::endl
                              << " cov1[py px] " <<  recPar->getCovMatrix()[1] << std::endl
                              << " cov1[py py] " <<  recPar->getCovMatrix()[2] << std::endl
                              << " cov1[pz px] " <<  recPar->getCovMatrix()[3] << std::endl
                              << " cov1[pz py] " <<  recPar->getCovMatrix()[4] << std::endl
                              << " cov1[pz pz] " <<  recPar->getCovMatrix()[5] << std::endl
                              << " cov1[e px] "  <<  recPar->getCovMatrix()[6] << std::endl
                              << " cov1[e py] "  <<  recPar->getCovMatrix()[7] << std::endl
                              << " cov1[e pz] "  <<  recPar->getCovMatrix()[8] << std::endl
                              << " cov1[e e] "   <<  recPar->getCovMatrix()[9] << std::endl
                              << std::endl;
  }
  

  _nEvt++;
  
}


void Add4MomCovMatrixCharged::check(LCEvent * evt) {
  // nothing to check here - could be used to fill
  // checkplots in reconstruction processor
}


void Add4MomCovMatrixCharged::end(){


  streamlog_out(MESSAGE4) << " processed "               << _nEvt
                          << " events in "               << _nRun
                          << " runs "                    << std::endl ;
}

