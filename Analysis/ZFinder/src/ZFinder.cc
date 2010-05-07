#include "EVENT/LCIO.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "ZFinder.h"
//fg: for now MarlinReco should not depend on ROOT - use CLHEP instead
//#include "TLorentzVector.h"
//#include "TVector3.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
typedef CLHEP::HepLorentzVector LorentzVector ;
typedef CLHEP::Hep3Vector Vector3D ;


// Marlin stuff
#include <marlin/Global.h>
// the event display

#include <cstdlib>

using namespace lcio;

ZFinder aZFinder;

ZFinder::ZFinder() : marlin::Processor("ZFinder") {

  registerProcessorParameter( "Printing" , 
			      "Print certain messages"  ,
			      _printing,
			       (int)1 ) ;


  std::string zdecay = "ee";
  registerProcessorParameter( "ZDecay" , 
			      "Z decay mode"  ,
			      _zdecay,
			      zdecay) ;


  std::string inputParticleCollectionName = "PandoraPFOs";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "InputParticleCollectionName" , 
			     "Input Particle Collection Name "  ,
			     _inputParticleCollectionName,
			      inputParticleCollectionName);


  std::string outputParticleCollectionName = "eeX";
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "OutputParticleCollectionName" , 
			     "Output Particle Collection Name "  ,
			     _outputParticleCollectionName,
			     outputParticleCollectionName);


  registerProcessorParameter( "AddPhotons" , 
			      "Include photons with reconstructed Z"  ,
			      _addPhotons,
			       (int)1 ) ;

  registerProcessorParameter( "CanUseClusterEnergyForElectrons" , 
			      "Allow possibilit of cluster not track used for electrons",
			      _canUseClusterEnergyForElectrons,
			      (int)1 ); 


  registerProcessorParameter( "FermionMomentumCut" , 
			      "Momentum Cut for fermion from Z decay"  ,
			      _momentumCut,
			       (float)10.) ;

  registerProcessorParameter( "CosTrackGammaCut" , 
			      "Minimum cosine of track-photon angle"  ,
			      _cosTrackGammaCut,
			       (float)0.999);

  registerProcessorParameter( "MaxDeltaMz" , 
			      "Maximum difference between candidate and Z mass"  ,
			      _dmzcut,
			      (float)50.);

  registerProcessorParameter( "MuonEcalEnergyCut" , 
			      "Cut on muon ECAL energy"  ,
			      _muonEcalEnergyCut,
			      (float)2.5);

  registerProcessorParameter( "MuonHcalEnergyCut" , 
			      "Cut on Muon HCAL energy"  ,
			      _muonHcalEnergyCut,
			      (float)10.);


  registerProcessorParameter( "MuonHcalEnergyCut1" , 
			      "Cut on Muon HCAL energy (alt)"  ,
			      _muonHcalEnergyCut1,
			      (float)5.);

  registerProcessorParameter( "ElectronEcalEnergyCut" , 
			      "Cut on electron ECAL energy"  ,
			      _electronEcalEnergyCut,
			      (float)10.);

  registerProcessorParameter( "ElectronHcalEnergyCut" , 
			      "Cut on Electron HCAL energy"  ,
			      _electronHcalEnergyCut,
			      (float)10.);


  registerProcessorParameter( "ElectronEoverPCutLow" , 
			      "Cut on Electron E/p cut"  ,
			      _electronEoPCutLow,
			      (float)0.5);

  registerProcessorParameter( "ElectronEoverPCutHigh" , 
			      "Cut on Electron E/p cut"  ,
			      _electronEoPCutHigh,
			      (float)1.3);


  return;

}

//===================================================================================

void ZFinder::init() { 
  return;
}

//===================================================================================

void ZFinder::processRunHeader( LCRunHeader* run) { 
  return;
}

//===================================================================================

void ZFinder::processEvent( LCEvent * evt ) { 

  // Make a new vector of particles
  LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  // Access PFO collection
  bool found = this->FindPFOs(evt);
  if(found){
    if(_printing>1)std::cout << "Analysis : " << _zdecay << std::endl; 
    if(_zdecay.find("ee")!=std::string::npos)this->FindZee(recparcol);
    if(_zdecay.find("mm")!=std::string::npos)this->FindZmumu(recparcol);
    if(_zdecay.find("mumu")!=std::string::npos)this->FindZmumu(recparcol);
    if(_zdecay.find("tautau")!=std::string::npos)streamlog_out(ERROR) << "ZFinder tautau decay mode not implemented" << std::endl; 
    if(_zdecay.find("qq")!=std::string::npos)streamlog_out(ERROR) << "ZFinder qq decay mode not implemented" << std::endl;
  }

  // Add new collection to event
  evt->addCollection( recparcol , _outputParticleCollectionName.c_str() );
  
  return;
  
}

//===================================================================================

void ZFinder::end(){ 
  return;
}

//===================================================================================

bool ZFinder::FindPFOs( LCEvent* evt ) {

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

void ZFinder::FindZmumu(LCCollectionVec * recparcol) {

  if(_printing>1)std::cout << "FindZmumu : " << _pfovec.size() << std::endl; 

  // Look for muon candidates (loose cuts)
  std::vector<LorentzVector>pmuplus;
  std::vector<LorentzVector>pmuminus;
  std::vector<ReconstructedParticle*>pMplusPfoVec;
  std::vector<ReconstructedParticle*>pMminusPfoVec;
  for(unsigned int i=0;i<_pfovec.size();i++){  
    if(_pfovec[i]->getCharge()!=0){
      ClusterVec c1 = _pfovec[i]->getClusters();
      float ecal1 = 0;
      float hcal1 = 0;
      float muon1 = 0;
      if(c1.size()==1){
	ecal1 = c1[0]->getSubdetectorEnergies()[0] + c1[0]->getSubdetectorEnergies()[3];
	hcal1 = c1[0]->getSubdetectorEnergies()[1] + c1[0]->getSubdetectorEnergies()[4];
	muon1 = c1[0]->getSubdetectorEnergies()[2];    
      }
      bool muon = true;
      if(_pfovec[i]->getEnergy()<_momentumCut)muon=false;
      if(ecal1>_muonEcalEnergyCut)muon=false;
      if(hcal1>_muonHcalEnergyCut)muon=false;
      if(hcal1>_muonHcalEnergyCut1&& muon1<0.01)muon=false;
      if(muon){
	LorentzVector pmu(_pfovec[i]->getMomentum()[0] ,_pfovec[i]->getMomentum()[1] ,_pfovec[i]->getMomentum()[2], _pfovec[i]->getEnergy()  );
	if(_pfovec[i]->getCharge()>0.5){
	  pmuplus.push_back(pmu);
	  pMplusPfoVec.push_back(_pfovec[i]);
	}
	if(_pfovec[i]->getCharge()<-0.5){
	  pmuminus.push_back(pmu);
	  pMminusPfoVec.push_back(_pfovec[i]);
	}
      }
    }
  }  
 
  // loop over pairs of candidate muons and look for best Z candidates

  ReconstructedParticle* pMplus  = NULL;
  ReconstructedParticle* pMminus = NULL;
  float bestdmz = 999.;
  unsigned int   besti=0;
  unsigned int   bestj=0;
  for(unsigned int i=0;i<pmuplus.size();i++){  
    for(unsigned int j=0;j<pmuminus.size();j++){  
      LorentzVector pmumu = pmuplus[i]+pmuminus[j];
      float massmumu = pmumu.m();
      if( fabs(massmumu-91.2) < bestdmz){
	bestdmz = fabs(massmumu-91.2);
	besti=i; bestj=j;
	pMplus  = pMplusPfoVec[i];
	pMminus = pMminusPfoVec[j];
      }
    }
  }

  // Make the reconstructed particle and add any photons
  float Energy = 0;
  float Mom[3] = {0.,0.,0.};
  if(bestdmz<_dmzcut && pMplus!=NULL && pMminus!=NULL){
    ReconstructedParticleImpl * recoPart = new ReconstructedParticleImpl();
    recoPart->addParticle(pMplus);
    recoPart->addParticle(pMminus);
    float pxem = pMminus->getMomentum()[0];
    float pyem = pMminus->getMomentum()[1];
    float pzem = pMminus->getMomentum()[2];
    float eem  = pMminus->getEnergy();
    float pxep = pMplus->getMomentum()[0];
    float pyep = pMplus->getMomentum()[1];
    float pzep = pMplus->getMomentum()[2];
    float eep  = pMplus->getEnergy();
    Energy = eep+eem;
    Mom[0] = Mom[0] + pxem+pxep;
    Mom[1] = Mom[1] + pyem+pyep;
    Mom[2] = Mom[2] + pzem+pzep;

    // look for FSR photon candidates 
    if(_addPhotons>0){
      for(unsigned int i=0;i<_pfovec.size();i++){  
	if(_pfovec[i]->getType()==22){
	  float pxg  = _pfovec[i]->getMomentum()[0];
	  float pyg  = _pfovec[i]->getMomentum()[1];
	  float pzg  = _pfovec[i]->getMomentum()[2];
	  float eg   = _pfovec[i]->getEnergy();
	  float cosm = 0.;
	  if(eg>0&&eem>0)cosm=(pxg*pxem+pyg*pyem+pzg*pzem)/eg/eem;
	  float cosp = 0.;
	  if(eg>0&&eep>0)cosp=(pxg*pxep+pyg*pyep+pzg*pzep)/eg/eep;
	  if(cosp>_cosTrackGammaCut||cosm>_cosTrackGammaCut){
	    recoPart->addParticle(_pfovec[i]);
	    Energy+=eg;
	    Mom[0]+=pxg;
	    Mom[1]+=pyg;
	    Mom[2]+=pzg;
	  }
	}
      }
    }

    // set reconstructed particle parameters
    float Mass = (Energy*Energy-Mom[0]*Mom[0]-Mom[1]*Mom[1]-Mom[2]*Mom[2]);
    if(Mass>0)Mass=sqrt(Mass);
    recoPart->setMomentum( Mom );
    recoPart->setEnergy( Energy );
    recoPart->setMass( Mass );
    recoPart->setCharge( 0. );
    recoPart->setType( 94 );
    if(_printing>1)std::cout << "Found mmX : " << Mass << std::endl; 

    // add it to the collection
    recparcol->addElement( recoPart );
  }

  return;
}

//===================================================================================

void ZFinder::FindZee(LCCollectionVec * recparcol) {

  if(_printing>1)std::cout << "FindZee : " << _pfovec.size() << std::endl; 

  // Look for electron candidates (loose cuts)
  std::vector<LorentzVector>peplus;
  std::vector<LorentzVector>peminus;
  std::vector<ReconstructedParticle*>pEplusPfoVec;
  std::vector<ReconstructedParticle*>pEminusPfoVec;
  for(unsigned int i=0;i<_pfovec.size();i++){  
    if(_pfovec[i]->getCharge()!=0 &&_pfovec[i]->getEnergy()>0){
      ClusterVec c1 = _pfovec[i]->getClusters();
      float ecal1 = 0;
      float hcal1 = 0;
      float muon1 = 0;
      if(c1.size()==1){
	ecal1 = c1[0]->getSubdetectorEnergies()[0] + c1[0]->getSubdetectorEnergies()[3];
	hcal1 = c1[0]->getSubdetectorEnergies()[1] + c1[0]->getSubdetectorEnergies()[4];
	muon1 = c1[0]->getSubdetectorEnergies()[2];    
      }
      bool electron = true;
      if(_pfovec[i]->getEnergy()<_momentumCut)electron=false;
      if(ecal1<_electronEcalEnergyCut)electron=false;
      if(hcal1>_electronHcalEnergyCut)electron=false;
      float eop = ecal1/_pfovec[i]->getEnergy();
      if(eop < _electronEoPCutLow)electron=false;
      if(eop > _electronEoPCutHigh)electron=false;
      if(muon1>0.01)electron=false;
      if(abs(_pfovec[i]->getType())==11 && _pfovec[i]->getEnergy()>_momentumCut/2.0)electron=true;
      if(electron){
	LorentzVector pe(_pfovec[i]->getMomentum()[0] ,_pfovec[i]->getMomentum()[1] ,_pfovec[i]->getMomentum()[2], _pfovec[i]->getEnergy()  );
	if(_pfovec[i]->getCharge()>0.5){
	  peplus.push_back(pe);
	  pEplusPfoVec.push_back(_pfovec[i]);
	}
	if(_pfovec[i]->getCharge()<-0.5){
	  peminus.push_back(pe);
	  pEminusPfoVec.push_back(_pfovec[i]);
	}
      }
    }
  }  
  
  // loop over pairs of candidate electrons and look for best Z candidates

  ReconstructedParticle* pEplus  = NULL;
  ReconstructedParticle* pEminus = NULL;
  float bestdmz = 999.;
  unsigned int   besti = 999;
  unsigned int   bestj = 999;
  for(unsigned int i=0;i<peplus.size();i++){  
    for(unsigned int j=0;j<peminus.size();j++){  
      LorentzVector pee = peplus[i]+peminus[j];
      float massee = pee.m();
      if( fabs(massee-91.2) < bestdmz){
	bestdmz = fabs(massee-91.2);
	besti=i; bestj=j;
	pEplus = pEplusPfoVec[i];
	pEminus = pEminusPfoVec[j];
      }
    }
  }
  
  // Make the reconstructed particle and add any photons
  float Energy = 0;
  float Mom[3] = {0.,0.,0.};
  if(bestdmz<_dmzcut && pEplus!=NULL && pEminus!=NULL){
    ReconstructedParticleImpl * recoPart = new ReconstructedParticleImpl();
    recoPart->addParticle(pEplus);
    recoPart->addParticle(pEminus);
    float pxem = pEminus->getMomentum()[0];
    float pyem = pEminus->getMomentum()[1];
    float pzem = pEminus->getMomentum()[2];
    float eem  = pEminus->getEnergy();
    float pxep = pEplus->getMomentum()[0];
    float pyep = pEplus->getMomentum()[1];
    float pzep = pEplus->getMomentum()[2];
    float eep  = pEplus->getEnergy();

    // look for FSR and Bremstahlung photon candidates 
    std::vector<ReconstructedParticle*>photons;
    float trackEPlus  = 0;
    float trackEMinus  = 0;
    float clusterEPlus  = 0;
    float clusterEMinus = 0;
    float sigmaEPlus  = 999.;
    float sigmaEMinus = 999.;
    float chiEPlus  = 999.;
    float chiEMinus = 999.;
    float sigpPlus  = 0.;
    float sigpMinus = 0.;


    Vector3D vecEPlus;
    Vector3D vecEMinus;
    const EVENT::ClusterVec ci   = pEplus->getClusters();
    const EVENT::TrackVec   ti   = pEplus->getTracks();
    trackEPlus   = pEplus->getEnergy();
    if(ti.size()==1)sigpPlus = trackEPlus*sqrt(ti[0]->getCovMatrix()[5])/fabs(ti[0]->getOmega());
    if(ci.size()>0){
      clusterEPlus = ci[0]->getEnergy();
      sigmaEPlus = 0.18*sqrt(clusterEPlus);
      chiEPlus = (trackEPlus-clusterEPlus)/sqrt(sigmaEPlus*sigmaEPlus+sigpPlus*sigpPlus);
      vecEPlus.set(ci[0]->getPosition()[0],ci[0]->getPosition()[1],ci[0]->getPosition()[2]);
    }
    const EVENT::ClusterVec cj   = pEminus->getClusters();
    const EVENT::TrackVec   tj   = pEminus->getTracks();
    trackEMinus   = pEminusPfoVec[bestj]->getEnergy();
    if(tj.size()==1)sigpMinus = trackEMinus*sqrt(tj[0]->getCovMatrix()[5])/fabs(tj[0]->getOmega());
    if(cj.size()>0){
      clusterEMinus = cj[0]->getEnergy();
      sigmaEMinus = 0.18*sqrt(clusterEMinus);
      chiEMinus = (trackEMinus-clusterEMinus)/sqrt(sigmaEMinus*sigmaEMinus+sigpMinus*sigpMinus);
      vecEMinus.set(cj[0]->getPosition()[0],cj[0]->getPosition()[1],cj[0]->getPosition()[2]);
    }
    
    if(_addPhotons>0){
      // loop over particles
      for(unsigned int i=0;i<_pfovec.size();i++){  
	if(_pfovec[i]->getType()==22){
	  float pxg  = _pfovec[i]->getMomentum()[0];
	  float pyg  = _pfovec[i]->getMomentum()[1];
	  float pzg  = _pfovec[i]->getMomentum()[2];
	  float eg   = _pfovec[i]->getEnergy();
	  float cosm = 0.;
	  if(eg>0&&eem>0)cosm=(pxg*pxem+pyg*pyem+pzg*pzem)/eg/eem;
	  float cosp = 0.;
	  if(eg>0&&eep>0)cosp=(pxg*pxep+pyg*pyep+pzg*pzep)/eg/eep;
	  // found a photon candidate near to electron/positron	  
	  if(cosp>_cosTrackGammaCut||cosm>_cosTrackGammaCut){
	    
	    // check that this isn't a split-off cluster
	    float drp = 999.;
	    float drm = 999.;
	    const EVENT::ClusterVec c   = _pfovec[i]->getClusters();
	    if(c.size()==1){
	      // distance between cluster and the electron/positron cluster
	      Vector3D vecg(c[0]->getPosition()[0],c[0]->getPosition()[1],c[0]->getPosition()[2]);
	      Vector3D v =  vecg.cross(vecEPlus);
	      float magg = vecg.mag();
	      if(magg>0)drp = v.mag()/magg;
	      v =  vecg.cross(vecEMinus);
	      if(magg>0)drm = v.mag()/vecg.mag();
	    }
	    if(cosp>cosm){
	      bool merge = false;
	      float chiNew = (trackEPlus-clusterEPlus-eg)/sigmaEPlus;
	      if(drp<20.0)merge = true;
	      // if fairly close merge if chi2 for matching improves greatly
	      if(drp<30.0 && chiEPlus>4.0 && fabs(chiNew)<chiEPlus)merge = true;
	      if(drp<40.0 && chiEPlus>5.0 && fabs(chiNew)<chiEPlus)merge = true;
	      if(drp<50.0 && chiEPlus>7.0 && fabs(chiNew)<chiEPlus)merge = true;
	      // sanity check
	      if(  fabs(chiEPlus)<2.0 && chiNew*chiNew>chiEPlus*chiEPlus+5.0)merge = false;
	      // always merge if very close - can't expect reconstruction to work
	      if(drp<10.0)merge = true;
	      // if we merge cluster do nothing except add photon to Z
	      if(merge){
		recoPart->addParticle(_pfovec[i]);
		clusterEPlus += eg;
		sigmaEPlus = 0.18*sqrt(clusterEPlus);
		chiEPlus = (trackEPlus-clusterEPlus)/sqrt(sigmaEPlus*sigmaEPlus+sigpPlus*sigpPlus);
	      }
	      // if not merged add photon and energy to Z
	      if(!merge){
		recoPart->addParticle(_pfovec[i]);
		photons.push_back(_pfovec[i]);
	      }
	    }else{
	      float chiNew = (trackEMinus-clusterEMinus-eg)/sigmaEMinus;
	      bool merge = false;
	      // always merge if very close 
	      if(drm<20.0)merge = true;
	      // if fairly close merge if chi2 for matching improves greatly
	      if(drm<30.0 && chiEMinus>4.0 && fabs(chiNew)<chiEMinus)merge = true;
	      if(drm<40.0 && chiEMinus>5.0 && fabs(chiNew)<chiEMinus)merge = true;
	      if(drm<50.0 && chiEMinus>7.0 && fabs(chiNew)<chiEMinus)merge = true;
	      // sanity check
	      if(  fabs(chiEMinus)<2.0 && chiNew*chiNew>chiEMinus*chiEMinus+5.0)merge = false;
	      // always merge if very close - can't expect reconstruction to work
	      if(drm<10.0)merge = true;
	      // if we merge cluster do nothing
	      if(merge){
		recoPart->addParticle(_pfovec[i]);
		clusterEMinus += eg;
		sigmaEMinus = 0.18*sqrt(clusterEMinus);
		chiEMinus = (trackEMinus-clusterEMinus)/sqrt(sigmaEMinus*sigmaEMinus+sigpMinus*sigpMinus);
	      }
	      // if not merged then add photon to Z 
	      if(!merge){
		recoPart->addParticle(_pfovec[i]);
		photons.push_back(_pfovec[i]);
	      }
	    }
	  }
	}
      } // end of loop over particles
    } // end of adding photons


    // first the electrons and positrons
    Energy = eep+eem;
    Mom[0] = Mom[0] + pxem+pxep;
    Mom[1] = Mom[1] + pyem+pyep;
    Mom[2] = Mom[2] + pzem+pzep;
    // any photons
    for(unsigned int i=0;i<photons.size();i++){ 
      Energy += photons[i]->getEnergy();
      Mom[0] += photons[i]->getMomentum()[0];
      Mom[1] += photons[i]->getMomentum()[1];
      Mom[2] += photons[i]->getMomentum()[2];
    }
    float Mass = sqrt(Energy*Energy-Mom[0]*Mom[0]-Mom[1]*Mom[1]-Mom[2]*Mom[2]);


    // Use cluster for electron/positron energy ? 
    if(_canUseClusterEnergyForElectrons && (chiEPlus<-3.5 || chiEMinus < -3.5)){
      // look for evidence of merged electrons/photons
      float MomX[3]={0.,0.,0.};
      float EnergyX=0;
      float scaleP = 1.0;
      if(chiEPlus<-3.5)scaleP  = clusterEPlus/trackEPlus;
      float scaleM = 1.0;
      if(chiEMinus<-3.5)scaleM = clusterEMinus/trackEMinus;
      EnergyX = Energy  + eep*(scaleP-1.) +eem*(scaleM-1.);
      MomX[0] = Mom[0] + pxem*(scaleM-1.)+pxep*(scaleP-1.);
      MomX[1] = Mom[1] + pyem*(scaleM-1.)+pyep*(scaleP-1.);
      MomX[2] = Mom[2] + pzem*(scaleM-1.)+pzep*(scaleP-1.);
      float MassX = sqrt(EnergyX*EnergyX-MomX[0]*MomX[0]-MomX[1]*MomX[1]-MomX[2]*MomX[2]);
      if(fabs(MassX-91.2)<fabs(Mass-91.2)){
	Energy = EnergyX;
	Mom[0] = MomX[0];
	Mom[1] = MomX[1];
	Mom[2] = MomX[2];
	Mass = sqrt(Energy*Energy-Mom[0]*Mom[0]-Mom[1]*Mom[1]-Mom[2]*Mom[2]);
      }
    }

    // Set parameters of reconstructed Z
    if(Mass>0)Mass=sqrt(Mass);
    recoPart->setMomentum( Mom );
    recoPart->setEnergy( Energy );
    recoPart->setMass( Mass );
    recoPart->setCharge( 0. );
    recoPart->setType( 94 );
    if(_printing>1)std::cout << "Found eeX : " << Mass << std::endl; 

    // add it to the collection
    recparcol->addElement( recoPart );
  }

  return;
}

