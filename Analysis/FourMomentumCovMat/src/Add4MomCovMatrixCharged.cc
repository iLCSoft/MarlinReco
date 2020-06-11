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

void Add4MomCovMatrixCharged::processRunHeader(LCRunHeader*  /*run*/) {

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
      bool newcov = false;
      // Only for charged particles
      if (TMath::Abs(recPar->getCharge()) > 0.5) {
          getCovMatrixMomenta(recPar, covarianceP);
          newcov = true;
      }  
      else if (TMath::Abs(recPar->getCharge())<0.001) {
         EVENT::TrackVec tracks = recPar->getTracks();
         int itype = recPar->getType();
         // only treat 2-track PFOs which are photon, K0 or Lambda
         if (tracks.size() == 2 && (itype == 22 || itype == 310 || itype == 3122)) {
            for (int itrk=0; itrk < 2; itrk++) { 
               Track *trk = tracks[i]; 
               // FIXME: a version of getCovMatrixMomenta which directly takes a track would be great here,
               // to avoid creation of a dummy PFO
               ReconstructedParticleImpl *trackPFO = new ReconstructedParticleImpl;
               trackPFO->addTrack(trk);
               FloatVec covariancePtrk(10, 0.0);
               // this function is in MarlinReco/Analysis/FourMomentumCovMat/src/algebraImplementation.cc
               getCovMatrixMomenta(trackPFO, covariancePtrk);
               for ( size_t j = 0; j < 10; j++) covarianceP[j] += covariancePtrk[j];
               delete trackPFO;
            }
            // now covarianceP contains the combined covariance matrix (hopefully)
            newcov = true;
         }
      } 
      if (newcov) {
          try
              {
              
                  recPar->setCovMatrix (covarianceP);
     
                  streamlog_out(DEBUG5) << "Set covariance matrix for charged pfo with charge "
                                        << (int) recPar->getCharge() << " and (px,py,pz,e) = ("
                                        << recPar->getMomentum()[0] << ","
                                        << recPar->getMomentum()[1] << ","
                                        << recPar->getMomentum()[2] << ","
                                        << recPar->getEnergy() << "):" << std::endl
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
          
          catch (EVENT::ReadOnlyException &e)
              {
                  std::cout << "WARNING::we cannot set member of pfo because "
                            << "PandoraPFO is read only "
                            << "'" << e.what() << "'" << std::endl;
              }
      }
        
  }
  
  _nEvt++;
  
}


void Add4MomCovMatrixCharged::check(LCEvent *  /*evt*/) {
  // nothing to check here - could be used to fill
  // checkplots in reconstruction processor
}


void Add4MomCovMatrixCharged::end(){


  streamlog_out(MESSAGE4) << " processed "               << _nEvt
                          << " events in "               << _nRun
                          << " runs "                    << std::endl ;
}

