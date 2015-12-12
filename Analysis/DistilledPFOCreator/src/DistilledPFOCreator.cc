#include "EVENT/LCIO.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "DistilledPFOCreator.h"
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

DistilledPFOCreator aDistilledPFOCreator;

DistilledPFOCreator::DistilledPFOCreator() : marlin::Processor("DistilledPFOCreator") {

  registerProcessorParameter( "Printing" , 
			      "Print certain messages"  ,
			      _printing,
			       (int)1 ) ;

  std::string inputParticleCollectionName1 = "PandoraPFOs";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "InputParticleCollectionName1" , 
			     "Input Particle Collection Name1 "  ,
			     _inputParticleCollectionName1,
			      inputParticleCollectionName1);

  std::string inputParticleCollectionName2 = "GammaGammaParticles";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "InputParticleCollectionName2" , 
			     "Input Particle Collection Name2 "  ,
			     _inputParticleCollectionName2,
			      inputParticleCollectionName2);

  std::string outputParticleCollectionName = "DistilledPFOs";
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "OutputParticleCollectionName" , 
			     "Output Particle Collection Name "  ,
			     _outputParticleCollectionName,
			     outputParticleCollectionName);

  return;

}

//===================================================================================

void DistilledPFOCreator::init() {
  if(_printing>1)printParameters(); 
  return;
}

//===================================================================================

void DistilledPFOCreator::processRunHeader( LCRunHeader* run) { 
  return;
}

//===================================================================================

void DistilledPFOCreator::processEvent( LCEvent * evt ) { 

  // Make a new vector of particles
  LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  recparcol->setSubset(true);                      

  // Access PFO collection
  bool found = this->FindPFOs(evt);
  if(found){
    if(_printing>1)std::cout << "Create DistilledPFOs : " << std::endl; 
    this->CreateDistilledPFOs(recparcol);
  }

  // Add new collection to event
  evt->addCollection( recparcol , _outputParticleCollectionName.c_str() );
  
  return;
  
}

//===================================================================================

void DistilledPFOCreator::end(){ 
  return;
}

//===================================================================================

bool DistilledPFOCreator::FindPFOs( LCEvent* evt ) {

// Find the PandoraPFOs and the GammaGammaParticles as input PFOs and store in appropriate vectors

  bool tf = false;

  // clear old vectors
  _pfovec.clear();
  _ggpfovec.clear();
  _mypfovec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;

// PandoraPFOs
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
  if(_printing>1)std::cout << "FindPFOs : (nPFOs (PandoraPFOs) = " << _pfovec.size() << " )" << std::endl; 

// GammaGammaParticles
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){    
    if(*name==_inputParticleCollectionName2){
      LCCollection* col = evt->getCollection(*name);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
	_ggpfovec.push_back(recoPart);
      }
    }
  }
  if(_printing>1)std::cout << "FindPFOs : (nPFOs (GammaGammaParticles) = " << _ggpfovec.size() << " )" << std::endl; 

  if(_printing>1)std::cout << "Find PFOs : " << tf << std::endl; 

  return tf;
  
}

//===================================================================================

void DistilledPFOCreator::CreateDistilledPFOs(LCCollectionVec * recparcol) {

  // Merge PandoraPFOs and GammaGammaParticles into DistilledPFOs taking care to avoid double counting
  //
  // Graham W. Wilson,  15th October 2015

  // Keep track of various energy sums
  float esumgg = 0.0;                          // GammaGammaParticles
  float esump = 0.0;                           // PandoraPFOs
  float esumd;                                 // DistilledPFOs
  float esumch = 0.0;
  float esumph = 0.0;
  int nphotons = 0;
  int ncharged = 0;
  int nothers = 0;
  float esumothers = 0.0;
  float chargesum = 0.0;
  float esumNH = 0.0;
  int nNH=0;
  int npi0=0; int neta=0; int netap=0;
  float epi0=0.0; float eeta=0.0; float eetap=0.0;
  float Vpi0=0.0; float Veta=0.0; float Vetap=0.0;
  int nppfo=0;

  if(_printing>1)std::cout << "CreateDistilledPFOs : (nPFOs initially = " << _pfovec.size() << " " << _ggpfovec.size() << " )" << std::endl;

  typedef std::set<const ReconstructedParticle*> PfoSet;
  PfoSet particles_used;    // Set with pointers to the daughter particles (photons for now) that are already used in the solution

  // For convenience sort the PandoraPFOs by energy
  std::sort(_pfovec.begin(),_pfovec.end(),DistilledPFOCreator::PfoEnergySortFunction);

  // And also sort the GammaGammaParticles by energy
  std::sort(_ggpfovec.begin(),_ggpfovec.end(),DistilledPFOCreator::PfoEnergySortFunction);

  for(unsigned int i = 0; i < _ggpfovec.size(); i++){
      if(_ggpfovec[i]->getType()==111){
         npi0++;
         epi0 += _ggpfovec[i]->getEnergy();
         Vpi0 += _ggpfovec[i]->getCovMatrix()[9];      // E-E term
      }
      if(_ggpfovec[i]->getType()==221){
         neta++;
         eeta += _ggpfovec[i]->getEnergy();
         Veta += _ggpfovec[i]->getCovMatrix()[9];      // E-E term
      }
      if(_ggpfovec[i]->getType()==331){
         netap++;
         eetap += _ggpfovec[i]->getEnergy();
         Vetap += _ggpfovec[i]->getCovMatrix()[9];      // E-E term
      }
      esumgg += _ggpfovec[i]->getEnergy();
      const ReconstructedParticleVec particles = _ggpfovec[i]->getParticles();
      if(_printing>3)std::cout << "CreateDistilledPFOs: (nparticles = " << particles.size() << " )" << std::endl;
      for(unsigned int j = 0; j < particles.size(); j++){           // TODO - maybe double check these are indeed photons ..
          const ReconstructedParticle *particle = particles[j];
          particles_used.insert(particle);                          // Keep note of particles that are constituents of a GammaGammaParticle
      }
  // Add this PFO to the output
      recparcol->addElement(_ggpfovec[i]);
  }
  if(_printing>3)std::cout << "CreateDistilledPFOs: (nparticles_used = " << particles_used.size() << " )" << std::endl;
          
  esumd = esumgg;
  for(unsigned int i = 0; i < _pfovec.size(); i++){
      esump += _pfovec[i]->getEnergy();
      if(abs(_pfovec[i]->getCharge())>0.5){
         ncharged++;
         chargesum+=_pfovec[i]->getCharge();
         esumch+= _pfovec[i]->getEnergy();
      }
      else{
         if(_pfovec[i]->getType()==2112){
            nNH++;
            esumNH += _pfovec[i]->getEnergy();
         }
         else if(_pfovec[i]->getType()==22){
            nphotons++;
            esumph += _pfovec[i]->getEnergy();
         }
         else{
            nothers++;
            esumothers += _pfovec[i]->getEnergy();
         }
      }
      if( particles_used.find(_pfovec[i]) == particles_used.end() ){   // Condition is true if the i'th PandoraPFO is not in the particles_used set 
  // So can add this PFO to the output
          recparcol->addElement(_pfovec[i]);
          esumd += _pfovec[i]->getEnergy();
          nppfo++;
      }
  }
  float V = Vpi0 + Veta + Vetap;
  int ntot = npi0 + neta + netap;

  if(_printing>3)std::cout << "CreateDistilledPFOs : (nPFOs finally = " << recparcol->getNumberOfElements() << " )" << std::endl;

  if(_printing>3)std::cout << "CreateDistilledPFOs SUMMARY " << ntot << " " << nppfo << " " 
                                                             << ncharged << " " << chargesum << " " << esumch << " "
                                                             << nphotons << " " << esumph << " "
                                                             << nNH << " " << esumNH << " "
                                                             << nothers << " " << esumothers << " "
                                                             << esump << " " 
                                                             << esumd << " "
                                                             << esumgg << " " 
                                                             << V << " " 
                                                             << npi0 << " " << neta << " " << netap << " " 
                                                             << epi0 << " " << eeta << " " << eetap << " " 
                                                             << Vpi0 << " " << Veta << " " << Vetap << " "  
                                                             << std::endl;
  return;
}

bool DistilledPFOCreator::PfoEnergySortFunction(ReconstructedParticle* lhs, ReconstructedParticle* rhs){

  //  Sort the ReconstructedParticles by energy

  //  true if lhs goes before

   const double lhs_energy  = lhs->getEnergy();
   const double rhs_energy  = rhs->getEnergy();

   return (lhs_energy>rhs_energy);

}

