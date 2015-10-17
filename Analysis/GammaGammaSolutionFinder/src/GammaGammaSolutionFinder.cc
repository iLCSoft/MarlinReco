#include "EVENT/LCIO.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "GammaGammaSolutionFinder.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
typedef CLHEP::HepLorentzVector LorentzVector ;
typedef CLHEP::Hep3Vector Vector3D ;

// Marlin stuff
#include <marlin/Global.h>
// the event display

// ROOT stuff
#include "TMath.h"
#include "TMatrixD.h"

#include <cstdlib>

using namespace lcio;

GammaGammaSolutionFinder aGammaGammaSolutionFinder;

GammaGammaSolutionFinder::GammaGammaSolutionFinder() : marlin::Processor("GammaGammaSolutionFinder") {

  registerProcessorParameter( "Printing" , 
			      "Print certain messages"  ,
			      _printing,
			       (int)1 ) ;

  std::string inputParticleCollectionName1 = "GammaGammaCandidatePi0s";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "InputParticleCollectionName1" , 
			     "Input Particle Collection Name1 "  ,
			     _inputParticleCollectionName1,
			      inputParticleCollectionName1);

  std::string inputParticleCollectionName2 = "GammaGammaCandidateEtas";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "InputParticleCollectionName2" , 
			     "Input Particle Collection Name2 "  ,
			     _inputParticleCollectionName2,
			      inputParticleCollectionName2);

  std::string inputParticleCollectionName3 = "GammaGammaCandidateEtaPrimes";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "InputParticleCollectionName3" , 
			     "Input Particle Collection Name3 "  ,
			     _inputParticleCollectionName3,
			      inputParticleCollectionName3);

  std::string outputParticleCollectionName = "GammaGammaSolutions";
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "OutputParticleCollectionName" , 
			     "Output Particle Collection Name "  ,
			     _outputParticleCollectionName,
			     outputParticleCollectionName);

  return;

}

//===================================================================================

void GammaGammaSolutionFinder::init() { 
  return;
}

//===================================================================================

void GammaGammaSolutionFinder::processRunHeader( LCRunHeader* run) { 
  return;
}

//===================================================================================

void GammaGammaSolutionFinder::processEvent( LCEvent * evt ) { 

  // Make a new vector of particles
  LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  recparcol->setSubset(true);                      

  // Access PFO collection
  bool found = this->FindPFOs(evt);
  if(found){
    if(_printing>1)std::cout << "Find GammaGamma Solution: " << std::endl; 
    this->FindGammaGammaSolutions(recparcol);
  }

  // Add new collection to event
  evt->addCollection( recparcol , _outputParticleCollectionName.c_str() );
  
  return;
  
}

//===================================================================================

void GammaGammaSolutionFinder::end(){ 
  return;
}

//===================================================================================

bool GammaGammaSolutionFinder::FindPFOs( LCEvent* evt ) {

// Add all 3 GammaGammaCandidate collections to the PFO vector (pi0s , etas, eta's)

  bool tf = false;

  // clear old vector
  _pfovec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;


// pi0s
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){    
    if(*name==_inputParticleCollectionName1){
      LCCollection* col = evt->getCollection(*name);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
	_pfovec.push_back(recoPart);
      }
    }
  }

  if(_printing>1)std::cout << "FindPFOs : (nPFOs = " << _pfovec.size() << " )" << std::endl; 

// etas
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){    
    if(*name==_inputParticleCollectionName2){
      LCCollection* col = evt->getCollection(*name);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
	_pfovec.push_back(recoPart);
      }
    }
  }
  if(_printing>1)std::cout << "FindPFOs : (nPFOs = " << _pfovec.size() << " )" << std::endl; 

// eta's
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){    
    if(*name==_inputParticleCollectionName3){
      LCCollection* col = evt->getCollection(*name);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
	_pfovec.push_back(recoPart);
      }
    }
  }

  if(_printing>1)std::cout << "FindPFOs : (nPFOs = " << _pfovec.size() << " )" << std::endl; 

  if(_printing>1)std::cout << "Find PFOs : " << tf << std::endl; 

  return tf;
  
}

//===================================================================================

void GammaGammaSolutionFinder::FindGammaGammaSolutions(LCCollectionVec * recparcol) {

  if(_printing>1)std::cout << "FindGammaGammaSolutions : (nPFOs = " << _pfovec.size() << " )" << std::endl;

  typedef std::set<const ReconstructedParticle*> PfoSet;
  PfoSet particles_assigned;    // Set with pointers to the daughter particles (photons for now) that are already used in the solution

  // For convenience sort the GammaGammaCandidates by fit probability
  std::sort(_pfovec.begin(),_pfovec.end(),GammaGammaSolutionFinder::PfoProbabilitySortFunction);

  // Algorithm 1. Based on fit probability ordering, add GammaGammaCandidates to the output GammaGammaParticles collection 
  // if they do not contain already assigned constituent particles. For now, the constituents are always two photons  
  // - but may not hurt to keep this generic for when other decay modes (eg. Dalitz decay) are included
  for(unsigned int i = 0; i < _pfovec.size(); i++){
      const ReconstructedParticleVec particles = _pfovec[i]->getParticles();
      if(_printing>3)std::cout << "FindGammaGammaSolutions: (nparticles = " << particles.size() << " )" << std::endl;
      bool assignable = true;   
      for(unsigned int j = 0; j < particles.size(); j++){           // TODO - maybe double check these are indeed photons ..
          const ReconstructedParticle *particle = particles[j];
          if( particles_assigned.find(particle) != particles_assigned.end() )assignable = false; // If particle is already assigned then it is not still assignable
      }
      if(assignable){
         recparcol->addElement(_pfovec[i]);                    // Add meson to the output collection if all constituent particles are still assignable
         for(unsigned int j = 0; j < particles.size(); j++){
             const ReconstructedParticle *particle = particles[j];
             particles_assigned.insert(particle);              // Update set of particles already assigned in this solution
         }
      }         
  }
  return;
}

bool GammaGammaSolutionFinder::PfoProbabilitySortFunction(ReconstructedParticle* lhs,ReconstructedParticle* rhs){

  // Sort gamma gamma candidates by type and by decreasing fit probability

  //  true if lhs goes before

/*
   const double lhs_energy  = lhs->getEnergy();
   const double rhs_energy  = rhs->getEnergy();
*/
   const int lhs_particleType  = lhs->getType();
   const int rhs_particleType  = rhs->getType();

  // GoodnessOfPID currently filled with fit probability - may at some point be more convenient to fill with the fit chi-squared ...
   const float lhs_GoodnessOfPID = lhs->getGoodnessOfPID();
   const float rhs_GoodnessOfPID = rhs->getGoodnessOfPID();

   if(lhs_particleType==rhs_particleType)return (lhs_GoodnessOfPID>rhs_GoodnessOfPID);  // This rank makes sense when GoodnessOfPID has the fit probability
   return (lhs_particleType<rhs_particleType);                                          // Favor lower valued types (pi0s cf etas cf eta's)

}

