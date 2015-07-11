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

/*
  for(unsigned int i=0;i<_pfovec.size();i++){
// Photon-ID
// Require type of photon
    if(_pfovec[i]->getType()==22){
      bool photon = true;
// Many neutral PFOs are not photons - so also explicitly reject if this is not identified as a photon
//      if(_pfovec[i]->getType()!=22)photon=false;
// Kinematic reduction of combinatorics with energy cut on accepted photons - but will also loose efficiency
      if(_pfovec[i]->getEnergy()<_gammaMomentumCut)photon=false;
      if(photon){
        LorentzVector pgam(_pfovec[i]->getMomentum()[0] ,_pfovec[i]->getMomentum()[1] ,_pfovec[i]->getMomentum()[2], _pfovec[i]->getEnergy() );
        pgamma.push_back(pgam);
        pGammaPfoVec.push_back(_pfovec[i]);
      }
    }
  }
*/
  
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
        if( fabs(mgg - _ggResonanceMass) < _dmggcut){
	  pGammai = pGammaPfoVec[i];
          pGammaj = pGammaPfoVec[j];
      // Make the reconstructed particle
          float Energy = 0;
          float Mom[3] = {0.,0.,0.};
          ReconstructedParticleImpl * recoPart = new ReconstructedParticleImpl();
          recoPart->addParticle(pGammai);
          recoPart->addParticle(pGammaj);
          float pxi = pGammai->getMomentum()[0];
          float pyi = pGammai->getMomentum()[1];
          float pzi = pGammai->getMomentum()[2];
          float ei  = pGammai->getEnergy();
          float pxj = pGammaj->getMomentum()[0];
          float pyj = pGammaj->getMomentum()[1];
          float pzj = pGammaj->getMomentum()[2];
          float ej  = pGammaj->getEnergy();
          Energy = ei+ej;
          Mom[0] = Mom[0] + pxi + pxj;
          Mom[1] = Mom[1] + pyi + pyj;
          Mom[2] = Mom[2] + pzi + pzj;

      // set reconstructed particle parameters
          float MassSquared = (Energy*Energy-Mom[0]*Mom[0]-Mom[1]*Mom[1]-Mom[2]*Mom[2]);
          float Mass = 0.0;
          if(MassSquared>0)Mass=sqrt(MassSquared);
          recoPart->setMomentum( Mom );
          recoPart->setEnergy( Energy );
          recoPart->setMass( Mass );
          recoPart->setCharge( 0. );
      // Default choice is Pi0
          recoPart->setType( 111 );
          if(_ggResonanceName=="Eta")recoPart->setType( 221 );
          if(_ggResonanceName=="EtaPrime")recoPart->setType( 331 );

          if(_printing>1)std::cout << "Found " << _ggResonanceName << " gg candidate " << Mass << " " << i << " " << j << " " << ei << " " << ej << std::endl; 

      // add it to the collection
          recparcol->addElement( recoPart );
        }
      }
    }
  }

// Not clear this is needed - but seem to be having issues with memory leaks.
  pgamma.clear();
  pGammaPfoVec.clear();

  return;
}
