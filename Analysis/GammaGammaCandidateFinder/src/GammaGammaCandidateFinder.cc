#include "EVENT/LCIO.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "GammaGammaCandidateFinder.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
typedef CLHEP::HepLorentzVector LorentzVector ;
typedef CLHEP::Hep3Vector Vector3D ;

// MarlinKinfit stuff
#include "JetFitObject.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "NewtonFitterGSL.h"
#include "MassConstraint.h"

// Marlin stuff
#include <marlin/Global.h>
// the event display

#include <cstdlib>

using namespace lcio;

GammaGammaCandidateFinder aGammaGammaCandidateFinder;

GammaGammaCandidateFinder::GammaGammaCandidateFinder() : marlin::Processor("GammaGammaCandidateFinder") {

  registerProcessorParameter( "Printing" , 
			      "Print certain messages"  ,
			      _printing,
			       (int)1 ) ;

  std::string ggResonanceName = "Pi0";
  registerProcessorParameter( "GammaGammaResonanceName" , 
			      "Particle decaying to Gamma Gamma"  ,
			      _ggResonanceName,
			      ggResonanceName) ;

  std::string inputParticleCollectionName = "PandoraPFOs";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "InputParticleCollectionName" , 
			     "Input Particle Collection Name "  ,
			     _inputParticleCollectionName,
			      inputParticleCollectionName);

  std::string outputParticleCollectionName = "GammaGammaCandidates";
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "OutputParticleCollectionName" , 
			     "Output Particle Collection Name "  ,
			     _outputParticleCollectionName,
			     outputParticleCollectionName);

  registerProcessorParameter( "GammaMomentumCut" , 
			      "Minimum Momentum Cut for Photon (GeV)"  ,
			      _gammaMomentumCut,
			       (float)0.5) ;

  registerProcessorParameter( "GammaGammaResonanceMass" , 
			      "Nominal Mass of Particle Decaying to Gamma Gamma (GeV)"  ,
			      _ggResonanceMass,
			       (float)0.135) ;

  registerProcessorParameter( "MaxDeltaMgg" , 
			      "Maximum difference between candidate mass and GammaGama Resonance mass (GeV)"  ,
			      _dmggcut,
			      (float)0.040);

  return;

}

//===================================================================================

void GammaGammaCandidateFinder::init() { 
  return;
}

//===================================================================================

void GammaGammaCandidateFinder::processRunHeader( LCRunHeader* run) { 
  return;
}

//===================================================================================

void GammaGammaCandidateFinder::processEvent( LCEvent * evt ) { 

  // Make a new vector of particles
  LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  // Access PFO collection
  bool found = this->FindPFOs(evt);
  if(found){
    if(_printing>1)std::cout << "Analysis : " << _ggResonanceName << std::endl; 
    this->FindGammaGammaCandidates(recparcol);
  }

  // Add new collection to event
  evt->addCollection( recparcol , _outputParticleCollectionName.c_str() );
  
  return;
  
}

//===================================================================================

void GammaGammaCandidateFinder::end(){ 
  return;
}

//===================================================================================

bool GammaGammaCandidateFinder::FindPFOs( LCEvent* evt ) {

  bool tf = false;

  // clear old vector
  _pfovec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){    
    if(*name==_inputParticleCollectionName){
      LCCollection* col = evt->getCollection(*name);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
	_pfovec.push_back(recoPart);
      }
    }
  }

  if(_printing>1)std::cout << "Find PFOs : " << tf << std::endl; 

  return tf;
  
}

//===================================================================================

void GammaGammaCandidateFinder::FindGammaGammaCandidates(LCCollectionVec * recparcol) {

// Take initial BestPi0 implementation and extend to all combinations consistent with pi0 mass

// TODO : Check whether these things should be floats or doubles.

  if(_printing>1)std::cout << "FindGammaGammaCandidates : (nPFOs = " << _pfovec.size() << " )" << std::endl; 

  // Look for photon candidates
  std::vector<LorentzVector>pgamma;
  std::vector<ReconstructedParticle*>pGammaPfoVec;

  for(unsigned int i=0;i<_pfovec.size();i++){
// Photon-ID
// Require charge = 0
    if( fabs( _pfovec[i]->getCharge() - 0.0) < 0.5 ){
      bool photon = true;
// Many neutral PFOs are not photons - so also explicitly reject if this is not identified as a photon
      if(_pfovec[i]->getType()!=22)photon=false;
// Kinematic reduction of combinatorics with energy cut on accepted photons - but will also loose efficiency
      if(_pfovec[i]->getEnergy()<_gammaMomentumCut)photon=false;
      if(photon){
        if(_printing>1)std::cout << "FindGammaGammaCandidates : Photon " << _pfovec[i]->getType() << std::endl;        
        LorentzVector pgam(_pfovec[i]->getMomentum()[0] ,_pfovec[i]->getMomentum()[1] ,_pfovec[i]->getMomentum()[2], _pfovec[i]->getEnergy() );
        pgamma.push_back(pgam);
        pGammaPfoVec.push_back(_pfovec[i]);
      }
    }
  }
  
  if(_printing>1)std::cout << "FindGammaGammaCandidates : (nphotons = " << pGammaPfoVec.size() << " " << pgamma.size() << " )" << std::endl;  
 
  // loop over pairs of candidate photons and keep all combinations consistent with the given parent mass assumption

  ReconstructedParticle* pGammai = NULL;
  ReconstructedParticle* pGammaj = NULL;

  if(pgamma.size() >=2 ){
    for(unsigned int i=0;i<pgamma.size()-1;i++){  
      for(unsigned int j=i+1;j<pgamma.size();j++){  
        LorentzVector pgg = pgamma[i]+pgamma[j];
        float mgg = pgg.m();
        if(_printing>2)std::cout << "Testing for " << _ggResonanceName << " " << i << " " << j << "  M = " << mgg << std::endl; 
// Select pairs consistent with resonance mass for further processing (storage and/or fitting)
        if( fabs(mgg - _ggResonanceMass) < _dmggcut){

	  pGammai = pGammaPfoVec[i];
          pGammaj = pGammaPfoVec[j];

// Set up mass constrained fit using code similar to WW5CFIT in MarlinKinfit

          double mass_constraint = double(_ggResonanceMass);
          MassConstraint mc(mass_constraint);

          LorentzVector pgi = pgamma[i];
          LorentzVector pgj = pgamma[j];

// Note - we currently use hard-wired approximate error estimates partly because the covariance matrix of the 
// photon 4-vectors is not yet properly filled on input 
// TODO: Also need to account for the different/better errors implicit in photon conversions ..

          JetFitObject j1(pgi.e(), pgi.theta(), pgi.phi(), 0.16*std::sqrt(pgi.e()), 0.001, 0.001, 0.0);
//          j1->setName("Photon1");
          JetFitObject j2(pgj.e(), pgj.theta(), pgj.phi(), 0.16*std::sqrt(pgj.e()), 0.001, 0.001, 0.0);
//          j2->setName("Photon2");

          mc.addToFOList(j1);
          mc.addToFOList(j2);

// Choose fit engine method TODO Add these as options to the steering file.
          int _ifitter = 0;
          BaseFitter *pfitter;
          if (_ifitter == 1) {
            pfitter = new NewFitterGSL();
          }
          else if (_ifitter == 2) {
            pfitter = new NewtonFitterGSL();
          }
          else {
            pfitter = new OPALFitterGSL();
          }
          BaseFitter &fitter = *pfitter;

          fitter.addFitObject(j1);
          fitter.addFitObject(j2);
          fitter.addConstraint(mc);

          double fit_probability = 0.0;
          fit_probability = fitter.fit();

          int nIterations = fitter.getIterations();

          int errorCode = fitter.getError();

          int cov_dim;
          double * cov = fitter.getGlobalCovarianceMatrix(cov_dim);  // 6x6

          if(_printing>3){
             std::cout << "Constrained fit results  RC: " << errorCode << std::endl;
             std::cout << "Measured mass = " << mgg << std::endl;
             std::cout << "No. of iterations " << nIterations << std::endl;
             std::cout << "Fit probability = " << fit_probability << std::endl;
             std::cout << "Covariance matrix dimension " << cov_dim << std::endl;
             if(cov_dim==6){
                for (unsigned int i=0; i<6*6; i++){
                    std::cout << "Covariance matrix element " << i << " " << cov[i] << std::endl;                   
                }
             }
          }

      // Make the reconstructed particle using the result of the constrained fit
      // likely need some minimal fit probability cut - configurable for each gamma-gamma hypothesis
          ReconstructedParticleImpl * recoPart = new ReconstructedParticleImpl();
          recoPart->addParticle(pGammai);
          recoPart->addParticle(pGammaj);
          double Energy;
          double Mom[3];
      // The 4-vector of the fitted gamma-gamma is the sum of the two fitted 4-vectors.
          Energy = j1.getE()  + j2.getE();
          Mom[0] = j1.getPx() + j2.getPx();
          Mom[1] = j1.getPy() + j2.getPy();
          Mom[2] = j1.getPz() + j2.getPz();
          recoPart->setEnergy( Energy );
          recoPart->setMomentum( Mom );
          recoPart->setMass( mass_constraint );
          recoPart->setCharge( 0.0 );
      // Also need eventually to fill other data-members 
      // Most important to store the fit probability and the 4-vector covariance matrix. 
      // May need some math to convert the 6*6 covariance matrix into the 4*4 one.
      // For now - store the fit probability in the goodnessofPiD variable
          float goodnessOfPID = float(fit_probability);
//          if( mgg  < _ggResonanceMass)goodnessOfPID = - goodnessOfPID;
          recoPart->setGoodnessOfPID( goodnessOfPID );

/*    // Make the reconstructed particle 
      // Was initially set to the sum of the measured photon vectors.
          double Energy;
          double Mom[3];
          ReconstructedParticleImpl * recoPart = new ReconstructedParticleImpl();
          recoPart->addParticle(pGammai);
          recoPart->addParticle(pGammaj);
          double pxi = pGammai->getMomentum()[0];
          double pyi = pGammai->getMomentum()[1];
          double pzi = pGammai->getMomentum()[2];
          double ei  = pGammai->getEnergy();
          double pxj = pGammaj->getMomentum()[0];
          double pyj = pGammaj->getMomentum()[1];
          double pzj = pGammaj->getMomentum()[2];
          double ej  = pGammaj->getEnergy();
          Energy = ei+ej;
          Mom[0] = pxi + pxj;
          Mom[1] = pyi + pyj;
          Mom[2] = pzi + pzj;
      // set reconstructed particle parameters
          float MassSquared = (Energy*Energy-Mom[0]*Mom[0]-Mom[1]*Mom[1]-Mom[2]*Mom[2]);
          float Mass = 0.0;
          if(MassSquared>0)Mass=sqrt(MassSquared);
          recoPart->setMomentum( Mom );
          recoPart->setEnergy( Energy );
          recoPart->setMass( Mass );
          recoPart->setCharge( 0. );
*/
      // Default choice is Pi0
          recoPart->setType( 111 );
          if(_ggResonanceName=="Eta")recoPart->setType( 221 );
          if(_ggResonanceName=="EtaPrime")recoPart->setType( 331 );

          if(_printing>1)std::cout << "Fitted " << _ggResonanceName << " gg candidate " 
                                   << mgg << " " << i << " " << j << " " << pgi.e() << " " << pgj.e() << 
                                   "Fitted: " << j1.getE() << " " << j2.getE() 
                                   << " Fit probability = " << fit_probability << std::endl; 

      // add it to the collection
          recparcol->addElement( recoPart );
        }
      }
    }
  }
// unclear whether this is ever needed/useful - likely does no harm
  pgamma.clear();
  pGammaPfoVec.clear();

  return;
}
