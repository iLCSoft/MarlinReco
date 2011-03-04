#include "CLICDstChecker.h"
#include <iostream>
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <iostream>
#include <math.h>
#include <map>
#include <marlin/Global.h>
#include <CalorimeterHitType.h>
#include "ClusterShapes.h"
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/VXDParameters.h>
#include <gear/GearParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/BField.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>

#include <limits>

using namespace lcio ;
using namespace marlin ;

const int precision = 2;
const int widthFloat = 7;
const int widthInt = 5;

CLICDstChecker aCLICDstChecker ;

CLICDstChecker::CLICDstChecker() : Processor("CLICDstChecker") {  
  _description = "Prints out information in DST in a useful format" ;  

  registerProcessorParameter("Debug",
			     "Activate debugging?",
			     m_debug,
			     int(0));

  registerProcessorParameter("Monitoring",
			     "Monitoring",
			     m_monitoring,
			     int(1));


  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                         "InputPfoCollection",
                          "Input PFO Collection",
                          m_inputPfoCollection,
                          std::string("PandoraPFANewPFOs"));

  registerProcessorParameter(
                         "PfoToMcRelationCollection",
                          "Input PFO to MC particle relation Collection",
                          m_inputPfoToMcRelationCollection,
                          std::string("RecoMCTruthLink"));

  registerInputCollection(LCIO::MCPARTICLE,
                         "McParticleCollection",
                          "Input MC particle collection",
                          m_inputMcParticleCollection,
                          std::string("MCParticlesSkimmed"));
   
}



void CLICDstChecker::init() { 

  printParameters();  
  m_nRun = -1 ;
  m_nEvt = 0 ;
 
}

void CLICDstChecker::processRunHeader( LCRunHeader* run) { 

  m_nRun++ ;
  m_nEvt = 0;
  std::cout << std::endl;
  std::cout << "CLICDstChecker ---> new run : run number = " << m_nRun << std::endl;

} 

void CLICDstChecker::processEvent( LCEvent * evt ) { 

  std::cout << " New Event " << std::endl;

  m_evt = evt;

  if (m_debug>=1) {
      std::cout << std::endl;
      std::cout << "CLICDstChecker -> run = " << m_nRun 
		<< "  event = " << m_nEvt << std::endl;
      std::cout << std::endl;
  }

  std::vector<ReconstructedParticle*>physicsPfos;
  std::vector<MCParticle*>physicsMatchedMcParticle;
  std::vector<ReconstructedParticle*>backgroundPfos;

 
  try {
    std::cout << " Reading collection " << m_inputPfoCollection.c_str() << std::endl;
    LCCollection * col = evt->getCollection(m_inputPfoCollection.c_str());
    int nelem = col->getNumberOfElements(); 
    for (int iPfo=0; iPfo<nelem; ++iPfo)m_pfoVector.push_back(dynamic_cast<ReconstructedParticle*>(col->getElementAt(iPfo)));
    std::sort(m_pfoVector.begin(),m_pfoVector.end(),CLICDstChecker::PfoSortFunction);
  }
  catch( DataNotAvailableException &e ) {
    std::cout << m_inputPfoCollection.c_str() << " collection is unavailable" << std::endl;
  };



  try {
    std::cout << " Reading collection " << m_inputMcParticleCollection.c_str() << std::endl;
    LCCollection * col = evt->getCollection(m_inputMcParticleCollection.c_str());
    int nelem = col->getNumberOfElements(); 
    for (int iMc=0; iMc<nelem; ++iMc)m_mcSet.insert(dynamic_cast<MCParticle*>(col->getElementAt(iMc)));
  }
  catch( DataNotAvailableException &e ) {
    std::cout << m_inputMcParticleCollection.c_str() << " collection is unavailable" << std::endl;
  };
 


  m_pfoToMcNavigator = NULL;
  try {
    std::cout << " Reading collection " << m_inputPfoToMcRelationCollection.c_str() << std::endl;
    m_pfoToMcNavigator = new LCRelationNavigator( evt->getCollection(m_inputPfoToMcRelationCollection.c_str()) );
  }
  catch( DataNotAvailableException &e ) {
    std::cout << m_inputPfoToMcRelationCollection.c_str() << " collection is unavailable" << std::endl;
  }
  

  for (unsigned int iPfo=0; iPfo<m_pfoVector.size(); ++iPfo) {
    ReconstructedParticle * pPfo = m_pfoVector[iPfo];
    LCObjectVec objectVec = m_pfoToMcNavigator->getRelatedToObjects(pPfo);
    if (objectVec.size() > 0) {
      for(unsigned int imc = 0; imc < objectVec.size(); imc++){
	MCParticle * pMC = dynamic_cast<MCParticle*>(objectVec[imc]);
	// since only saving skimmed set of MC particles, check pMC points to an existing object 
	if(m_mcSet.count(pMC)!=0){
	  physicsPfos.push_back(pPfo);
	  physicsMatchedMcParticle.push_back(pMC);
	}else{
	  backgroundPfos.push_back(pPfo);
	}
      }
    }
  }


  for(unsigned int ipass = 0; ipass<2;ipass++){
    if (m_monitoring&&ipass==0)std::cout << " ************** Physics    Pfos ****************" << std::endl;
    if (m_monitoring&&ipass==1)std::cout << " ************** Background Pfos ****************" << std::endl;
    int nPfos;
    (ipass==0) ? nPfos = physicsPfos.size() : nPfos = backgroundPfos.size();
    for (unsigned int iPfo=0; iPfo<nPfos; ++iPfo) {
      ReconstructedParticle * pPfo;
      (ipass==0) ? pPfo = physicsPfos[iPfo] : pPfo = backgroundPfos[iPfo];
      MCParticle * pMc = NULL;
      if(ipass==0)physicsMatchedMcParticle[iPfo];
      const int id = pPfo->id();
      const int type = pPfo->getType();
      const bool isCompound = pPfo->isCompound(); 
      float momentum[3];
      for(unsigned int i=0;i<3;i++)momentum[i] = pPfo->getMomentum()[i];
      const float pT = sqrt(momentum[0]*momentum[0]+momentum[1]*momentum[1]);
      const float  p = sqrt(pT*pT+momentum[2]*momentum[2]);
      const float cosTheta = fabs(momentum[2])/p; 
      const float energy  = pPfo->getEnergy();
      //eTotalInput+=energy;
      const float mass = pPfo->getMass();
      const float charge = pPfo->getCharge();
      const ParticleIDVec particleIDs = pPfo->getParticleIDs();
      ParticleID *particleIDUsed = pPfo->getParticleIDUsed();
      const float goodnessOfPID = pPfo->getGoodnessOfPID();
      const ReconstructedParticleVec particles = pPfo->getParticles();
      const ClusterVec clusters = pPfo->getClusters();
      const TrackVec   tracks   = pPfo->getTracks();
      
      for(unsigned int i = 0; i< tracks.size(); i++){
	Track *track = tracks[i];
	const int nHitsVTX = track->getSubdetectorHitNumbers()[6];
	const int nHitsFTD = track->getSubdetectorHitNumbers()[7];
	const int nHitsSIT = track->getSubdetectorHitNumbers()[8];
	const int nHitsTPC = track->getSubdetectorHitNumbers()[9];
	const int nHitsSET = track->getSubdetectorHitNumbers()[10];
	const int nHitsETD = track->getSubdetectorHitNumbers()[11];
	const float d0    = track->getD0();
	const float z0    = track->getZ0();
	const float omega = track->getOmega();
	const float tanL  = track->getTanLambda();
	const float phi0  = track->getPhi();
      }
      
      for(unsigned int i = 0; i< clusters.size(); i++){
	const Cluster *cluster = clusters[i];
      }

      if (m_monitoring) {
	std::cout << std::fixed;
	std::cout << std::setprecision(precision);
	if(clusters.size()==0)FORMATTED_OUTPUT_TRACK(type,energy,pT,cosTheta,tracks.size(),0);
	if(tracks.size()==0)FORMATTED_OUTPUT_CLUSTER(type,energy,pT,cosTheta,0,clusters.size());
	if(tracks.size()>0&&clusters.size()>0)FORMATTED_OUTPUT_TRACK_CLUSTER(type,energy,pT,cosTheta,tracks.size(),clusters.size());
      }
    }
  }
  
  if (m_monitoring)std::cout << " ************** Physics Pfos ****************" << std::endl;

  
  CleanUp();
  
  m_nEvt++;

}


void CLICDstChecker::CleanUp(){
  m_pfoVector.clear();
  m_mcSet.clear();
  delete m_pfoToMcNavigator;
  m_pfoToMcNavigator = NULL;
 
}

void CLICDstChecker::check(LCEvent * evt) { }

void CLICDstChecker::end() {

}

bool CLICDstChecker::PfoSortFunction(ReconstructedParticle* lhs,ReconstructedParticle* rhs){

  //  true if lhs goes before

   const float lhs_energy  = lhs->getEnergy();
   const float rhs_energy  = rhs->getEnergy();
   const ClusterVec lhs_clusters = lhs->getClusters();
   const TrackVec   lhs_tracks   = lhs->getTracks();
   const TrackVec   rhs_tracks   = rhs->getTracks();
   const ClusterVec rhs_clusters = rhs->getClusters();
   const int lhs_particleType  = abs(lhs->getType());
   const int rhs_particleType  = abs(rhs->getType());

   int lhs_type(0);
   int rhs_type(0);
   
   if(lhs_clusters.size()==0)lhs_type=1;
   if(rhs_clusters.size()==0)rhs_type=1;
   if(lhs_particleType==22)lhs_type=10;
   if(rhs_particleType==22)rhs_type=10;
   if(lhs_particleType==2112)lhs_type=20;
   if(rhs_particleType==2112)rhs_type=20;

   if(lhs_type==rhs_type)return (lhs_energy>rhs_energy);
   return (lhs_type<rhs_type);
   
}

