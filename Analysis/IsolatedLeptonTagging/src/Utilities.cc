// *******************************************************
// some useful functions
// *******************************************************
#include "Utilities.h"
#include "TVector3.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include <UTIL/LCRelationNavigator.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <iomanip>

namespace isolep{

Int_t getMCSerial(MCParticle *mcPart, LCCollection *colMCP) {
  // get the serial code of the particle in the MC particle list
  Int_t nMCP = colMCP->getNumberOfElements();
  Int_t serialID;
  for (Int_t i=0;i<nMCP;i++) {
    MCParticle *mcp = dynamic_cast<MCParticle*>(colMCP->getElementAt(i));
    if (mcp == mcPart) {
      serialID = i;
      return serialID;
    }
  }
  return -1;
}

  MCParticle *getMCParticle(ReconstructedParticle *recPart, LCCollection *colMCTL) {
    Double_t weight;
    Int_t nMCTL;
    return getMCParticle(recPart,colMCTL,weight,nMCTL);
  }

  MCParticle *getMCParticle(ReconstructedParticle *recPart, LCCollection *colMCTL, Double_t &weight) {
    Int_t nMCTL;
    return getMCParticle(recPart,colMCTL,weight,nMCTL);
  }
  
  MCParticle *getMCParticle(ReconstructedParticle *recPart, LCCollection *colMCTL, Double_t &weight, Int_t &nMCTL) {
    // get the corresponding MC particle of one reconstructed particle using the MCTruthLinker information
    MCParticle *mcLinkedParticle = NULL;
    //    Int_t iLink = -1;
    LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);
    LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(recPart);
    FloatVec vecWgtMCTL = navMCTL->getRelatedToWeights(recPart);
    weight = 0.;
    nMCTL = vecMCTL.size();
    for (Int_t i=0;i<nMCTL;i++) {
      MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[i]);
      if (vecWgtMCTL[i] > weight) {   // find the linked particle with largest weight as the mc truth particle
	weight = vecWgtMCTL[i];
	mcLinkedParticle = mcPart;
      }
    }
    delete navMCTL;
    return mcLinkedParticle;
  }

Int_t getLinkedMCParticle(ReconstructedParticle *recPart, LCCollection *colMCTL, Double_t &weight, Int_t &nMCTL) {
  // get the corresponding MC particle of one reconstructed particle using the MCTruthLinker information
  //  MCParticle *mcLinkedParticle = NULL;
  Int_t iLink = -1;
  LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);
  LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(recPart);
  FloatVec vecWgtMCTL = navMCTL->getRelatedToWeights(recPart);
  Double_t mcEnergyMax = -1.0;
  weight = 0.;
  nMCTL = vecMCTL.size();
  for (Int_t i=0;i<nMCTL;i++) {
    MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[i]);
    Double_t mcEnergy = mcPart->getEnergy();
    if (mcEnergy > mcEnergyMax) {   // find the linked particle with largest energy as the mc truth particle
      mcEnergyMax = mcEnergy;
      //      mcLinkedParticle = mcPart;
      weight = vecWgtMCTL[i];
      iLink = i;
    }
  }
  //  return mcLinkedParticle;
  return iLink;
}

Int_t getOriginalPDG(MCParticle *mcPart, Bool_t iHiggs) {
  // get the PDG of the original particle where the MCParticle comes from
  //  Int_t nParents = mcPart->getNumberOfParents();
  Int_t nParents = mcPart->getParents().size();
  while (nParents > 0) {
    //    MCParticle *mother = mcPart->getParent(0);
    MCParticle *mother = mcPart->getParents()[0];
    //    nParents = mother->getNumberOfParents();
    nParents = mother->getParents().size();
    mcPart = mother;
    Int_t pdg = mcPart->getPDG();
    if (!iHiggs) {
      if (abs(pdg) == 24 || abs(pdg) == 23 || abs(pdg) == 25) break;
    }
    else {
      if (abs(pdg) == 25) break;
    }
  }
  Int_t originalPDG = mcPart->getPDG();
  return originalPDG;
}

  Int_t getOriginalPDGForIsoLep(MCParticle *mcPart) {
    // get the PDG of the original particle where the MCParticle comes from
    //  Int_t nParents = mcPart->getNumberOfParents();
    Int_t nParents = mcPart->getParents().size();
    //    if (mcPart->isOverlay()) return -22;
    Double_t charge = mcPart->getCharge();
    while (nParents > 0) {
      MCParticle *mother = mcPart->getParents()[0];
      if (nParents > 1) {
	for (Int_t i=0;i<nParents;i++) {
	  MCParticle *temp = mcPart->getParents()[i];
	  if (TMath::Abs(temp->getCharge() - charge) < 0.1) {
	    mother = temp;
	    break;
	  }
	}
      }
      mcPart = mother;
      nParents = mother->getParents().size();
      Int_t pdg = mcPart->getPDG();
      if (abs(pdg) == 13 || abs(pdg) == 11 || abs(pdg) == 15) {
	if (mcPart->getParents().size() > 0) {
	  MCParticle *mmother = mcPart->getParents()[0];	
	  Int_t mpdg = mmother->getPDG();
	  if (abs(mpdg) == 11) {
	    //	  std::cerr << "Mother is an electron!" << endl;
	    break;
	  }
	}
      }
      else if (abs(pdg) == 24 || abs(pdg) == 23 || pdg == 25 || pdg == 21 || abs(pdg) == 15) {
	break;
      }
      else if (abs(pdg) <= 6 && abs(pdg) >= 1) {
	break;
      }

    }
    //    if (mcPart->isOverlay()) return -22;
    Int_t originalPDG = mcPart->getPDG();
    return originalPDG;
  }

  Int_t getOriginalPDGForIsoLep(MCParticle *mcPart, LCCollection *colMC) {
    // get the PDG of the original particle where the MCParticle comes from
    Int_t nParents = mcPart->getParents().size();
    Int_t serial = getMCSerial(mcPart,colMC);
    if (mcPart->isOverlay()) return -22;
    MCParticle *previous = mcPart;
    while (nParents > 0 && serial > 11) {
      Double_t charge = previous->getCharge();
      MCParticle *mother = mcPart->getParents()[0];
      if (nParents > 1) {
	for (Int_t i=0;i<nParents;i++) {
	  MCParticle *temp = mcPart->getParents()[i];
	  if (TMath::Abs(temp->getCharge() - charge) < 0.1) {
	    mother = temp;
	    break;
	  }
	}
      }
      previous = mcPart;
      mcPart = mother;
      if (mcPart->isOverlay()) break;
      nParents = mcPart->getParents().size();
      serial = getMCSerial(mcPart,colMC);
    }
    if (mcPart->isOverlay()) return -22;
    Int_t originalPDG = mcPart->getPDG();
    return originalPDG;
  }

Int_t getOriginalSerial(MCParticle *mcPart, LCCollection *colMCP, Bool_t iHiggs) {
  // get the serial number of the original particle where the MCParticle comes from
  //  Int_t nParents = mcPart->getNumberOfParents();
  Int_t nParents = mcPart->getParents().size();
  while (nParents > 0) {
    //    MCParticle *mother = mcPart->getParent(0);
    MCParticle *mother = mcPart->getParents()[0];
    //    nParents = mother->getNumberOfParents();
    nParents = mother->getParents().size();
    mcPart = mother;
    Int_t pdg = mcPart->getPDG();
    if (!iHiggs) {
      if (abs(pdg) == 24 || abs(pdg) == 23 || abs(pdg) == 25) break;
    }
    else {
      if (abs(pdg) == 25) break;
    }
  }
  Int_t originalSerial = getMCSerial(mcPart,colMCP);
  return originalSerial;
}

Int_t getOriginalSerial(ReconstructedParticle *recPart, LCCollection *colMCTL, LCCollection *colMCP, Bool_t iHiggs) {
  // get the serial number of the original particle where the PFO comes from

  Int_t originalSerial = -1;
  LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);
  LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(recPart);
  if (vecMCTL.size() > 0) {
    MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[0]);
    originalSerial = getOriginalSerialForZHH(mcPart,colMCP);
  }
  return originalSerial;
}

Int_t getOriginalSerialForZHH(MCParticle *mcPart, LCCollection *colMCP) {
  // get the serial number of the original particle where the MCParticle comes from
  //  Int_t nParents = mcPart->getNumberOfParents();
  Int_t nParents = mcPart->getParents().size();
  while (nParents > 0) {
    //    MCParticle *mother = mcPart->getParent(0);
    MCParticle *mother = mcPart->getParents()[0];
    //    nParents = mother->getNumberOfParents();
    nParents = mother->getParents().size();
    mcPart = mother;
    Int_t pdg = mcPart->getPDG();
    if (abs(pdg) == 25) {
      if (mcPart->getParents().size() > 0) {
	MCParticle *mmother = mcPart->getParents()[0];	
	Int_t mpdg = mmother->getPDG();
	if (mpdg == 11) break;
      }
    }
  }
  Int_t originalSerial = getMCSerial(mcPart,colMCP);
  return originalSerial;
}

Int_t getLeptonID(ReconstructedParticle *recPart) {
  // electron identification using ratios of energies deposited in ECal, HCal and Momentum
  Int_t iLeptonType = 0;
  Double_t fElectronCut1 = 0.5;      // lower edge of totalCalEnergy/momentum
  Double_t fElectronCut2 = 1.3;      // upper edge of totalCalEnergy/momentum
  Double_t fElectronCut3 = 0.9;     // lower edge of ecalEnergy/totalCalEnergy
  Double_t fMuonCut1 = 0.3;      // upper edge of totalCalEnergy/momentum
  //  Double_t fMuonCut2 = 0.5;      // upper edge of ecalEnergy/totalCalEnergy
  Double_t fMuonCut3 = 1.2;      // lower edge of yoke energy
  Double_t fPhotonCut1 = 0.7;      // lower edge of totalCalEnergy/momentum
  Double_t fPhotonCut2 = 1.3;      // upper edge of totalCalEnergy/momentum
  Double_t fPhotonCut3 = 0.9;     // lower edge of ecalEnergy/totalCalEnergy

  //  Double_t energy = recPart->getEnergy();
  Double_t charge = recPart->getCharge();
  Double_t ecalEnergy = 0;
  Double_t hcalEnergy = 0;
  Double_t yokeEnergy = 0;
  Double_t totalCalEnergy = 0;
  //  Int_t nHits = 0;
  std::vector<lcio::Cluster*> clusters = recPart->getClusters();
  for (std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();iCluster!=clusters.end();++iCluster) {
    ecalEnergy += (*iCluster)->getSubdetectorEnergies()[0];
    hcalEnergy += (*iCluster)->getSubdetectorEnergies()[1];
    yokeEnergy += (*iCluster)->getSubdetectorEnergies()[2];
    ecalEnergy += (*iCluster)->getSubdetectorEnergies()[3];
    hcalEnergy += (*iCluster)->getSubdetectorEnergies()[4];
    //    CalorimeterHitVec calHits = (*iCluster)->getCalorimeterHits();
    //    nHits = calHits.size();
  }
  totalCalEnergy = ecalEnergy + hcalEnergy;
  TVector3 momentum = TVector3(recPart->getMomentum());
  Double_t momentumMagnitude = momentum.Mag();
  Double_t fEpsilon = 1.E-10;
  if (totalCalEnergy/momentumMagnitude > fElectronCut1 && totalCalEnergy/momentumMagnitude < fElectronCut2 &&
      ecalEnergy/(totalCalEnergy + fEpsilon) > fElectronCut3 && TMath::Abs(charge) > 0.5) {
    iLeptonType = 11;
  }
  else if (totalCalEnergy/momentumMagnitude > fPhotonCut1 && totalCalEnergy/momentumMagnitude < fPhotonCut2 &&
      ecalEnergy/(totalCalEnergy + fEpsilon) > fPhotonCut3 && TMath::Abs(charge) < 0.5) {
    iLeptonType = 22;
  }
  //  else if (TMath::Abs(charge) > 0.5 && totalCalEnergy/momentumMagnitude < fMuonCut1 && 
  //	   ecalEnergy/(totalCalEnergy + fEpsilon) < fMuonCut2) {
  else if (TMath::Abs(charge) > 0.5 && totalCalEnergy/momentumMagnitude < fMuonCut1 && 
  	   yokeEnergy > fMuonCut3) {
      iLeptonType = 13;
  }

  return iLeptonType;
}

Bool_t getFSRTag(ReconstructedParticle *motherPart, ReconstructedParticle *recPart, Double_t fCosFSRCut) {
  // tag the particle recPart if it is from the Bremmstrahlung or Final-State-Radiation of another particle motherPart
  Bool_t isFSR = kFALSE;
  //  Double_t charge = motherPart->getCharge(); // mother particle should be charged
  TVector3 momentumMother = TVector3(motherPart->getMomentum());
  TVector3 momentum = TVector3(recPart->getMomentum());
  Double_t cosFSR = momentum.Dot(momentumMother)/momentum.Mag()/momentumMother.Mag();
  Int_t iType = getLeptonID(recPart);
  Int_t iTypeMother = getLeptonID(motherPart);
  //  if (TMath::Abs(charge) > 0.5 && cosFSR > fCosFSRCut) {
  //    if (iType == 11 || iType == 22) {
  if (cosFSR > fCosFSRCut && (iTypeMother == 11 || iTypeMother == 13)) {
    if (iType == 22) {
      isFSR = kTRUE;
    }
  }
  return isFSR;  
}

Bool_t getSplitTag(ReconstructedParticle *motherPart, ReconstructedParticle *recPart) {
  // developed from Mark Tomthon's ZFinder Processor
  // tag the particle recPart if it is the cluster splited from motherPart
  Bool_t isSplit = kFALSE;

  // information of mother particle
  Double_t trackMom=0.;
  Double_t clusterMom=0.;
  Double_t sigmaMom=999.;
  Double_t chiMom=999.;
  Double_t sigpMom=0.;
  TVector3 vecMom;
  const EVENT::ClusterVec ci   = motherPart->getClusters();
  const EVENT::TrackVec   ti   = motherPart->getTracks();
  trackMom = motherPart->getEnergy();
  if(ti.size()==1)sigpMom = trackMom*sqrt(ti[0]->getCovMatrix()[5])/TMath::Abs(ti[0]->getOmega());
  if(ci.size()>0){
    clusterMom = ci[0]->getEnergy();
    sigmaMom   = 0.18*sqrt(clusterMom);
    chiMom = (trackMom-clusterMom)/TMath::Sqrt(sigmaMom*sigmaMom+sigpMom*sigpMom);
    vecMom = TVector3(ci[0]->getPosition()[0],ci[0]->getPosition()[1],ci[0]->getPosition()[2]);
  }

  // calculate the distance between cluster and the mother particle cluster
  Double_t dr = 999.;
  const EVENT::ClusterVec c   = recPart->getClusters();
  if(c.size()==1){
    TVector3 vecg(c[0]->getPosition()[0],c[0]->getPosition()[1],c[0]->getPosition()[2]);
    TVector3 v =  vecg.Cross(vecMom);
    float magg = vecg.Mag();
    if(magg>0) dr = v.Mag()/magg;
  }
  
  // criteria for tag
  Double_t eg = recPart->getEnergy();
  float chiNew = (trackMom-clusterMom-eg)/sigmaMom;
  if(dr<20.0) isSplit = true;
  // if fairly close merge if chi2 for matching improves greatly
  if(dr<30.0 && chiMom>4.0 && TMath::Abs(chiNew)<chiMom) isSplit = true;
  if(dr<40.0 && chiMom>5.0 && TMath::Abs(chiNew)<chiMom) isSplit = true;
  if(dr<50.0 && chiMom>7.0 && TMath::Abs(chiNew)<chiMom) isSplit = true;
  // sanity check
  if(TMath::Abs(chiMom)<2.0 && chiNew*chiNew>chiMom*chiMom+5.0) isSplit = false;
  // always merge if very close - can't expect reconstruction to work
  if(dr<10.0) isSplit = true;

  // empirical correction
  Double_t costheta = getCosTheta(motherPart,recPart);
  if (costheta<0.99995 && eg/trackMom < 0.03) isSplit = true;

  return isSplit;  
}

Double_t getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone) { 
  // get the cone energy of the particle
  Int_t nPFOs = colPFO->getNumberOfElements();
  Double_t coneEnergy = 0.;
  TVector3 momentum0 = TVector3(recPart->getMomentum());
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    if (pfo != recPart) {
      TVector3 momentum = TVector3(pfo->getMomentum());
      Double_t cosTheta = momentum0.Dot(momentum)/momentum0.Mag()/momentum.Mag();
      if (cosTheta > cosCone) {
	coneEnergy += pfo->getEnergy();
      }
    }
  }
  return coneEnergy;
}

void getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, Bool_t woFSR, Double_t coneEnergy[3], Double_t pFSR[4]) { 
  // get the cone energy of the particle
  //  woFSR = kTRUE;
  TLorentzVector lortzFSR = TLorentzVector(0,  0,  0,  0);
  Int_t nPFOs = colPFO->getNumberOfElements();
  TVector3 momentum0 = TVector3(recPart->getMomentum());
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    if (pfo == recPart) continue;
    Double_t energy = pfo->getEnergy();
    TVector3 momentum = TVector3(pfo->getMomentum());
    if (woFSR) {
      Bool_t isFSR = getFSRTag(recPart,pfo);
      if (isFSR) {
	lortzFSR += TLorentzVector(momentum,energy);
	continue;
      }
    }
    Double_t charge = pfo->getCharge();
    Double_t cosTheta = momentum0.Dot(momentum)/momentum0.Mag()/momentum.Mag();
    if (cosTheta > cosCone) {
      coneEnergy[0] += energy;
      if (TMath::Abs(charge) < 0.5) {
	coneEnergy[1] += energy;
      }
      else {
	coneEnergy[2] += energy;
      }
    }
  }
  pFSR[0] = lortzFSR.Px();
  pFSR[1] = lortzFSR.Py();
  pFSR[2] = lortzFSR.Pz();
  pFSR[3] = lortzFSR.E();

  return;
}

void getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, Bool_t woFSR, Double_t coneEnergy[3], Double_t pFSR[4], 
		   Double_t cosCone2, Double_t pCone2[4], Int_t &nConePhoton ) { 
  // get the cone energy of the particle
  // add another larger cone
  //  woFSR = kTRUE;
  TLorentzVector lortzFSR = TLorentzVector(0,  0,  0,  0);
  TLorentzVector lortzCon = TLorentzVector(0,  0,  0,  0);
  Int_t nPFOs = colPFO->getNumberOfElements();
  TVector3 momentum0 = TVector3(recPart->getMomentum());
  nConePhoton = 0;
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    if (pfo == recPart) continue;
    Double_t energy = pfo->getEnergy();
    TVector3 momentum = TVector3(pfo->getMomentum());
    if (woFSR) {
      Bool_t isFSR = getFSRTag(recPart,pfo);
      if (isFSR) {
	lortzFSR += TLorentzVector(momentum,energy);
	Int_t iType = getLeptonID(pfo);
	if (iType == 22) nConePhoton++;
	continue;
      }
    }
    Double_t charge = pfo->getCharge();
    Double_t cosTheta = momentum0.Dot(momentum)/momentum0.Mag()/momentum.Mag();
    if (cosTheta > cosCone) {
      coneEnergy[0] += energy;
      if (TMath::Abs(charge) < 0.5) {
	coneEnergy[1] += energy;
      }
      else {
	coneEnergy[2] += energy;
      }
      Int_t iType = getLeptonID(pfo);
      if (iType == 22) nConePhoton++;
    }
    if (cosTheta > cosCone2) {
      lortzCon += TLorentzVector(momentum,energy);
    }
  }
  pFSR[0] = lortzFSR.Px();
  pFSR[1] = lortzFSR.Py();
  pFSR[2] = lortzFSR.Pz();
  pFSR[3] = lortzFSR.E();
  pCone2[0] = lortzCon.Px();
  pCone2[1] = lortzCon.Py();
  pCone2[2] = lortzCon.Pz();
  pCone2[3] = lortzCon.E();
  return;
}

TLorentzVector getFSRMomentum(ReconstructedParticle *recPart, LCCollection *colPFO) { 
  // get the cone energy of the particle
  TLorentzVector lortzFSR = TLorentzVector(0.,0.,0.,0.);
  Int_t nPFOs = colPFO->getNumberOfElements();
  TVector3 momentum0 = TVector3(recPart->getMomentum());
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    if (pfo == recPart) continue;
    Bool_t isFSR = getFSRTag(recPart,pfo);
    if (! isFSR) continue;
    Double_t energy = pfo->getEnergy();
    TVector3 momentum = TVector3(pfo->getMomentum());
    lortzFSR += TLorentzVector(momentum,energy);
  }

  return lortzFSR;
}

Double_t getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, Int_t mode) { 
  // get the cone energy of the particle
  Int_t nPFOs = colPFO->getNumberOfElements();
  Double_t coneEnergy = 0.,coneEnergyC = 0.,coneEnergyN = 0.;
  TVector3 momentum0 = TVector3(recPart->getMomentum());
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    if (pfo != recPart) {
      TVector3 momentum = TVector3(pfo->getMomentum());
      Int_t iCharge = pfo->getCharge();
      Double_t cosTheta = momentum0.Dot(momentum)/momentum0.Mag()/momentum.Mag();
      if (cosTheta > cosCone) {
	coneEnergy += pfo->getEnergy();
	if (iCharge == 0) {
	  coneEnergyN += pfo->getEnergy();
	}
	else {
	  coneEnergyC += pfo->getEnergy();
	}
      }
    }
  }
  if (mode == 0) {
    return coneEnergy;
  }
  else if (mode == 1) {
    return coneEnergyC;
  }
  else if (mode == 2) {
    return coneEnergyN;
  }
  else {
    return 99999.;
  }
}

Double_t getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, 
		       std::vector<lcio::ReconstructedParticle*> &conePFOs) {
  // get the cone energy of the particle
  Int_t nPFOs = colPFO->getNumberOfElements();
  Double_t coneEnergy = 0.;
  TVector3 momentum0 = TVector3(recPart->getMomentum());
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    if (pfo != recPart) {
      TVector3 momentum = TVector3(pfo->getMomentum());
      Double_t cosTheta = momentum0.Dot(momentum)/momentum0.Mag()/momentum.Mag();
      if (cosTheta > cosCone) {
	coneEnergy += pfo->getEnergy();
	conePFOs.push_back(pfo);
      }
    }
  }
  return coneEnergy;
}

Double_t getInvariantMass(ReconstructedParticle *recPart1, ReconstructedParticle *recPart2) {
  // get the invariant mass of two particles
  TVector3 momentum1 = TVector3(recPart1->getMomentum());
  TVector3 momentum2 = TVector3(recPart2->getMomentum());
  Double_t energy1 = recPart1->getEnergy();
  Double_t energy2 = recPart2->getEnergy();
  Double_t invariantMass = 0.;
  Double_t invariantMass2 = (energy1+energy2)*(energy1+energy2)-(momentum1+momentum2).Mag2();
  if (invariantMass2 > 0.) {
    invariantMass = sqrt(invariantMass2);
  }
  else {
    invariantMass = -sqrt(TMath::Abs(invariantMass2));
  }
  return invariantMass;
}

Double_t getInvariantMass(ReconstructedParticle *recPart1, ReconstructedParticle *recPart2, ReconstructedParticle *recPart3) {
  // get the invariant mass of three particles
  TVector3 momentum1 = TVector3(recPart1->getMomentum());
  TVector3 momentum2 = TVector3(recPart2->getMomentum());
  TVector3 momentum3 = TVector3(recPart3->getMomentum());
  Double_t energy1 = recPart1->getEnergy();
  Double_t energy2 = recPart2->getEnergy();
  Double_t energy3 = recPart3->getEnergy();
  Double_t invariantMass = 0.;
  Double_t invariantMass2 = (energy1+energy2+energy3)*(energy1+energy2+energy3)-(momentum1+momentum2+momentum3).Mag2();
  if (invariantMass2 > 0.) {
    invariantMass = sqrt(invariantMass2);
  }
  else {
    invariantMass = -sqrt(TMath::Abs(invariantMass2));
  }
  return invariantMass;
}

Int_t isSelectedByFastJet( ReconstructedParticle *pfo, LCCollection *colFastJet, Double_t &ratioEPartEJet, Double_t &ratioPTMJet ) { 
  // check the PFO if it is selectred by FastJet clustering
  Int_t iFastJet = 0;
  Int_t nJets = colFastJet->getNumberOfElements();
  Int_t iJet = -1;
  for (Int_t i=0;i<nJets;i++) {
    ReconstructedParticle *jet = dynamic_cast<ReconstructedParticle*>(colFastJet->getElementAt(i));
    std::vector<lcio::ReconstructedParticle*> partVec = jet->getParticles();
    for (std::vector<lcio::ReconstructedParticle*>::const_iterator iPart=partVec.begin();iPart!=partVec.end();++iPart) {
      if ((*iPart) == pfo) {
	iFastJet = 1;
	iJet = i;
	break;
      }
    }
  }
  if (iJet >= 0) {
    ReconstructedParticle *theJet = dynamic_cast<ReconstructedParticle*>(colFastJet->getElementAt(iJet)); // the jet where the pfo belongs
    // get the variables used by LAL Lepton Finder
    Double_t ePart = pfo->getEnergy();
    TVector3 pPart = pfo->getMomentum();
    Double_t eJet  = theJet->getEnergy();
    TVector3 pJet  = theJet->getMomentum();
    //    Double_t mJet  = theJet->getMass();
    Double_t mJet  = eJet*eJet > pJet.Mag2() ? TMath::Sqrt(eJet*eJet-pJet.Mag2()) : TMath::Sqrt(-eJet*eJet+pJet.Mag2());
    ratioEPartEJet = ePart/eJet;
    ratioPTMJet = pPart.Pt(pJet)/mJet;
  }

  return iFastJet;
}

void doPhotonRecovery(ReconstructedParticle *electron, LCCollection *colPFO, ReconstructedParticleImpl *recoElectron, Double_t fCosFSRCut) {
  // recover the BS and FSR photons
  TLorentzVector lortzElectron = TLorentzVector(electron->getMomentum(),electron->getEnergy());
  Int_t nPFOs = colPFO->getNumberOfElements();
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    if (recPart == electron) continue;
    Bool_t isFSR = getFSRTag(electron,recPart,fCosFSRCut);
    if (! isFSR) continue;
    recoElectron->addParticle(recPart);
    Bool_t isSplit = getSplitTag(electron,recPart);
    if (isSplit) continue;
    lortzElectron += TLorentzVector(recPart->getMomentum(),recPart->getEnergy());
  }
  Double_t energy = lortzElectron.E();
  Double_t mass   = lortzElectron.M();
  Double_t momentum[3] = {lortzElectron.Px(),lortzElectron.Py(),lortzElectron.Pz()};
  Double_t charge = electron->getCharge();
  recoElectron->setMomentum(momentum);
  recoElectron->setEnergy(energy);
  recoElectron->setMass(mass);
  recoElectron->setCharge(charge);
  recoElectron->setType(94);

}

void doPhotonRecovery(ReconstructedParticle *electron, LCCollection *colPFO, ReconstructedParticleImpl *recoElectron, Double_t fCosFSRCut, 
		      Int_t lepType, std::vector<lcio::ReconstructedParticle*> &photons) {
  // recover the BS and FSR photons
  TLorentzVector lortzElectron = TLorentzVector(electron->getMomentum(),electron->getEnergy());
  recoElectron->addParticle(electron);
  Int_t nPFOs = colPFO->getNumberOfElements();
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    if (recPart == electron) continue;
    if (isFoundInVector(recPart,photons)) continue;
    Bool_t isFSR = getFSRTag(electron,recPart,fCosFSRCut);
    if (! isFSR) continue;
    photons.push_back(recPart);
    recoElectron->addParticle(recPart);
    if (lepType == 11) {
      // do split algorithm only for electron
      Bool_t isSplit = getSplitTag(electron,recPart);
      if (isSplit) continue;
    }
    else if (lepType == 13) {
    }
    lortzElectron += TLorentzVector(recPart->getMomentum(),recPart->getEnergy());
  }
  Double_t energy = lortzElectron.E();
  Double_t mass   = lortzElectron.M();
  Double_t momentum[3] = {lortzElectron.Px(),lortzElectron.Py(),lortzElectron.Pz()};
  Double_t charge = electron->getCharge();
  recoElectron->setMomentum(momentum);
  recoElectron->setEnergy(energy);
  recoElectron->setMass(mass);
  recoElectron->setCharge(charge);
  recoElectron->setType(94);

}

Bool_t isFoundInVector(ReconstructedParticle *pfo, std::vector<lcio::ReconstructedParticle*> &pfos) {
  for (std::vector<lcio::ReconstructedParticle*>::const_iterator ipfo=pfos.begin();ipfo!=pfos.end();++ipfo) {
    if ((*ipfo) == pfo) {
      return kTRUE;
    }
  }
  return kFALSE;
}

Double_t jetFunction( TLorentzVector lortz, float beta, float power) {
  // extended Georgi Jet Function
  
  Double_t energy = lortz.E();
  Double_t momentum = lortz.P();

  if (!energy) {
    return -9999999.;
  }
  else {
    return TMath::Power(energy,power)*(1.-beta+beta*momentum*momentum/energy/energy);
  }

}

Double_t jetFunction( TLorentzVector lortz, float virtual_scale) {
  // Junping Jet Function
  
  Double_t energy = lortz.E();
  Double_t mass   = lortz.M();

  if (!energy) {
    return -9999999.;
  }
  else {
    return energy*(1. - mass/virtual_scale);
  }

}

Int_t calculateEnergyComponents(LCCollection *colMC, LCCollection *colMCTL, ReconstructedParticle *jet, Double_t energyComponents[4]) {
  // calculate energies from different color singlets in a jet
  // [0]: from first Higgs; [1]: from second Higgs; [2]: from other color-singlet (mostly Z); [3]: from overlay

  for (Int_t i=0;i<4;i++) {
    energyComponents[i] = 0.;
  }

  LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);

  Int_t np = jet->getParticles().size();
  //  cerr << "Debug: Npar in jet = " << np << endl;
  for (Int_t i=0;i<np;i++) {
    ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle *>(jet->getParticles()[i]);
    Double_t energy = pfo->getEnergy();
    LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(pfo);
    //    cerr << "Debug: energy = " << energy << " ; mclink = " << vecMCTL.size() << endl;
    if (vecMCTL.size() > 0) {
      MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[0]);
      Int_t originalSerial = getOriginalSerialForZHH(mcPart,colMC);
      //      cerr << "Debug: Original Serial = " << originalSerial << endl;
      if (originalSerial == 8) { // 1st H
	energyComponents[0] += energy;
      }
      else if (originalSerial == 9) { // 2nd H
	energyComponents[1] += energy;
      }
      else if (mcPart->isOverlay()) { // overlay
	energyComponents[3] += energy;
      }
      else {
	energyComponents[2] += energy;
      }
    }
  }

  // find the largest energy component
  Int_t iMax = -1;
  Double_t energyMax = -1.;
  for (Int_t i=0;i<4;i++) {
    if (energyComponents[i] > energyMax) {
      energyMax = energyComponents[i];
      iMax = i;
    }
    // fraction
    //    energyComponents[i] /= jet->getEnergy();
  }

  return iMax+1;
}

Int_t calculateEnergyComponents(LCCollection *colMC, LCCollection *colMCTL, ReconstructedParticle *jet, Double_t energyComponents[4],
				Int_t iColor1, Int_t iColor2, Int_t iColor3) {
  // calculate energies from different color singlets in a jet
  // [0]: from first Higgs; [1]: from second Higgs; [2]: from other color-singlet (mostly Z); [3]: from overlay
  // for qqH

  for (Int_t i=0;i<4;i++) {
    energyComponents[i] = 0.;
  }

  LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);

  Int_t np = jet->getParticles().size();
  //  cerr << "Debug: Npar in jet = " << np << endl;
  for (Int_t i=0;i<np;i++) {
    ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle *>(jet->getParticles()[i]);
    Double_t energy = pfo->getEnergy();
    LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(pfo);
    //    cerr << "Debug: energy = " << energy << " ; mclink = " << vecMCTL.size() << endl;
    if (vecMCTL.size() > 0) {
      MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[0]);
      Int_t originalSerial = getOriginalSerialForZHH(mcPart,colMC);
      //      cerr << "Debug: Original Serial = " << originalSerial << endl;
      if (originalSerial == iColor1) { // H
	energyComponents[0] += energy;
      }
      else if (originalSerial == iColor2 || originalSerial == iColor3) { // q or q-bar from Z
	energyComponents[1] += energy;
      }
      else if (mcPart->isOverlay()) { // overlay
	energyComponents[3] += energy;
      }
      else {
	energyComponents[2] += energy;
      }
    }
  }

  // find the largest energy component
  Int_t iMax = -1;
  Double_t energyMax = -1.;
  for (Int_t i=0;i<4;i++) {
    if (energyComponents[i] > energyMax) {
      energyMax = energyComponents[i];
      iMax = i;
    }
    // fraction
    //    energyComponents[i] /= jet->getEnergy();
  }

  return iMax+1;
}

Double_t getLikelihood(TString fname, TString hist, Double_t mass)
{

  TFile file(fname);
  TH1D *hmass = dynamic_cast<TH1D *>(file.Get(hist));
  Int_t nbin = hmass->GetNbinsX();
  Double_t ntot = hmass->GetEntries();
  Double_t epsilon = 0.001;
  Double_t prob = 0.;
  for (Int_t j=0;j<nbin;j++) {
    Double_t nj = hmass->GetBinContent(j+1);
    Double_t tj = hmass->GetBinCenter(j+1);
    Double_t delta = hmass->GetBinWidth(j+1);
    Double_t hj = TMath::Power(4./3/ntot,1./5)*delta*TMath::Sqrt(ntot/(nj+epsilon));
    Double_t density = 1.0*nj/ntot/TMath::Sqrt(2*3.1416)/hj*TMath::Exp(-(mass-tj)*(mass-tj)/2/hj/hj);
    //    cerr << "hj: " << hj << " nj: " << nj << " tj: " << tj 
    //	 << " delta: " << delta << " mass: " << mass << endl;
    //    cerr << "Prob Density: " << density << endl;
    prob += density;
  }

  file.Close();
  //  cerr << "Prob: " << prob <<endl; 
  return prob;
}

void listMCParticles(LCCollection *colMC) {
  // a detail look into MC particles and their decay chain
  if (!colMC) {
    std::cerr << "No MC Collection Found!" << std::endl;
    return ;
  }
  Int_t nMCP = colMC->getNumberOfElements();

  std::cout << setw(6) << "Index" << setw(7) << "PDG" << setw(7) << "Mother" << setw(10) << "Charge" << setw(10) << "Mass"	 
	    << setw(11) << "Energy" << setw(18) << "NumberOfDaughters" << setw(16) << "NumberOfParents" 
	    << setw(10) << "Original" << setw(8) << "Overlay" << setw(8) << "FromSim" << std::endl;

  for (Int_t i=0;i<nMCP;i++) {
    MCParticle *mcPart = dynamic_cast<MCParticle*>(colMC->getElementAt(i));
    Int_t pdg = mcPart->getPDG();
    Int_t nparents = mcPart->getParents().size();
    Int_t motherpdg = 0;
    if (nparents > 0) {
      MCParticle *mother = mcPart->getParents()[0];
      motherpdg = mother->getPDG();
    }
    Double_t charge = mcPart->getCharge();
    Double_t mass = mcPart->getMass();
    Double_t energy = mcPart->getEnergy();
    //    TVector3 pv = TVector3(mcPart->getMomentum());
    Int_t ndaughters = mcPart->getDaughters().size();
#if 0
    Int_t daughterpdg = 0;
    if (ndaughters > 0) {
      MCParticle *daughter = mcPart->getDaughters()[0];
      daughterpdg = daughter->getPDG();
    }
#endif    
    //    TLorentzVector lortz = TLorentzVector(pv,energy);
    Int_t originalPDG = getOriginalPDG(mcPart);
    Int_t ioverlay = mcPart->isOverlay()? 1 : 0;
    Int_t icreatedinsim = mcPart->isCreatedInSimulation()? 1 : 0;
    std::cout << setw(6) << i << setw(7) << pdg << setw(7) << motherpdg << setw(10) << setprecision(3) << charge << setw(10) << mass 
	      << setw(11) << energy << setw(18) << ndaughters << setw(16) << nparents 
	      << setw(10) << originalPDG << setw(8) << ioverlay << setw(8) << icreatedinsim << std::endl;
  }
}

  Bool_t isOverlay(ReconstructedParticle *pfo, LCCollection *colMCTL) {

  Bool_t iOverlay = false;

  LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);
  LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(pfo);
  if (vecMCTL.size() > 0) {
    MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[0]);
    if (mcPart->isOverlay()) iOverlay = true;
  }
  
  return iOverlay;
  
}

  void dumpJetParticles(ReconstructedParticle *jet, LCCollection *colMC, LCCollection *colMCTL) {
    // list the particles of jets

    LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);

    TLorentzVector lortz = TLorentzVector(jet->getMomentum(),jet->getEnergy());
    cerr << "----Particles of the Jet: (E,Px,Py,Pz) = (" << lortz.E() << "," << lortz.Px() << ","
	 << lortz.Py() << "," << lortz.Pz() << ")------" << endl;

    Int_t np = jet->getParticles().size();
    cerr << setw(6) << "Index" << setw(7) << "PDG" << setw(7) << "Mother" << setw(10) << "Charge" << setw(10) << "Mass" 
	 << setw(11) << "Energy" << setw(18) << "NumberOfDaughters" << setw(16) << "NumberOfParents" 
	 << setw(10) << "Original" << setw(8) << "Overlay" << setw(8) << "FromSim" << endl;
    for (Int_t i=0;i<np;i++) {
      ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle *>(jet->getParticles()[i]);
      LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(pfo);
      if (vecMCTL.size() > 0) {
	MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[0]);
	Int_t pdg = mcPart->getPDG();
	Int_t nparents = mcPart->getParents().size();
	Int_t motherpdg = 0;
	if (nparents > 0) {
	  MCParticle *mother = mcPart->getParents()[0];
	  motherpdg = mother->getPDG();
	}
	Double_t charge = mcPart->getCharge();
	Double_t mass = mcPart->getMass();
	Double_t energy = mcPart->getEnergy();
	TVector3 pv = TVector3(mcPart->getMomentum());
	Int_t ndaughters = mcPart->getDaughters().size();
#if 0
	Int_t daughterpdg = 0;
	if (ndaughters > 0) {
	  MCParticle *daughter = mcPart->getDaughters()[0];
	  daughterpdg = daughter->getPDG();
	}
#endif	
	Int_t originalPDG = getOriginalPDG(mcPart);
	Int_t ioverlay = mcPart->isOverlay()? 1 : 0;
	Int_t icreatedinsim = mcPart->isCreatedInSimulation()? 1 : 0;
	cerr << setw(6) << i << setw(7) << pdg << setw(7) << motherpdg << setw(10) << setprecision(3) << charge 
	     << setw(10) << mass << setw(11) << energy << setw(18) << ndaughters << setw(16) << nparents 
	     << setw(10) << originalPDG << setw(8) << ioverlay << setw(8) << icreatedinsim << endl;
      }
    }
  }

  // *******************************************************
  // *******************************************************
  Int_t getVertexComponents(LCCollection *colMC, LCCollection *colMCTL, ReconstructedParticle *vertex, Double_t energyComponents[2], Int_t nparticles[2]) {
  // calculate energies from signal particles or overlay in a vertex
  // [0]: from signal; [1]: from overlay

    for (Int_t i=0;i<2;i++) {
      energyComponents[i] = 0.;
      nparticles[i] = 0;
    }
    
    LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);
    
    Int_t np = vertex->getParticles().size();
    for (Int_t i=0;i<np;i++) {
      ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle *>(vertex->getParticles()[i]);
      Double_t energy = pfo->getEnergy();
      LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(pfo);
      if (vecMCTL.size() > 0) {
	MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[0]);
	if (mcPart->isOverlay()) { // overlay
	  energyComponents[1] += energy;
	  nparticles[1] += 1;
	}
	else {
	  energyComponents[0] += energy;
	  nparticles[0] += 1;
	}
      }
    }
    
    // find the largest energy component
    if (energyComponents[0] > energyComponents[1]) {
      return 0;
    }
    else {
      return 1;
    }
  }
  
  // *******************************************************
  // *******************************************************
  Int_t isVertexFromOverlay( ReconstructedParticle *vertex, LCCollection *colMC, LCCollection *colMCTL) {
    // 1: from overlay; 
    // 0: from signal;
    Double_t energy[2];
    Int_t np[2];
    return getVertexComponents(colMC,colMCTL,vertex,energy,np);

  }

  // *******************************************************
  // *******************************************************
  Double_t getJetDistance( ReconstructedParticle *i, ReconstructedParticle *j, TString algorithm, Double_t R) {
    // distance betweet jet i and j

    TLorentzVector lortz_i = TLorentzVector(i->getMomentum(),i->getEnergy());
    TLorentzVector lortz_j = TLorentzVector(j->getMomentum(),j->getEnergy());
    Double_t pt_i = lortz_i.Pt();
    Double_t pt_j = lortz_j.Pt();
    Double_t ptmin = pt_i < pt_j ? pt_i : pt_j;
    Double_t y = 999.;
    if (algorithm == "kt") {
      Double_t phi_i = lortz_i.Phi();
      Double_t phi_j = lortz_j.Phi();
      Double_t rapidity_i = lortz_i.Rapidity();
      Double_t rapidity_j = lortz_j.Rapidity();
      Double_t deltaR2 = (phi_i-phi_j)*(phi_i-phi_j)+(rapidity_i-rapidity_j)*(rapidity_i-rapidity_j);
      y = ptmin*ptmin*deltaR2/R/R;
    }
    return y;
  }

  // *******************************************************
  // *******************************************************
  Int_t isFromVertex( ReconstructedParticle *recPart, LCCollection *colVertex) {
    // check if the particle is from any reconstructed vertex
    
    Int_t nVtx = colVertex->getNumberOfElements();
    for (Int_t i=0;i<nVtx;i++) {
      ReconstructedParticle *vertex = dynamic_cast<ReconstructedParticle*>(colVertex->getElementAt(i));
      Int_t np = vertex->getParticles().size();
      for (Int_t j=0;j<np;j++) {
	ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle *>(vertex->getParticles()[j]);
	if (pfo == recPart) return 1;
      }
    }
    
    return 0;
  }

  // *******************************************************
  // *******************************************************
  std::vector<Int_t> getHiggsDecayModes(LCCollection *colMC) {
    // get the Higgs decay mode by MC information

    std::vector<Int_t> nHDecay;
    Int_t nHbb=0,nHWW=0,nWqq=0,nWlv=0,nHgg=0,nHZZ=0,nHcc=0;
    Int_t nHmumu=0,nHgamgam=0,nZvv=0,nZqq=0,nZee=0,nZmm=0,nHtautau=0;
    Int_t nMCP = colMC->getNumberOfElements();
    for (Int_t i=0;i<nMCP;i++) {
      MCParticle *mcPart = dynamic_cast<MCParticle*>(colMC->getElementAt(i));
      Int_t pdg = mcPart->getPDG();
      Int_t ndaughters = mcPart->getDaughters().size();
      Int_t daughterpdg = 0;
      if (ndaughters > 0) {
	MCParticle *daughter = mcPart->getDaughters()[0];
	daughterpdg = daughter->getPDG();
      }
      Int_t nparents = mcPart->getParents().size();
      Int_t motherpdg = 0;
      if (nparents > 0) {
	MCParticle *mother = mcPart->getParents()[0];
	motherpdg = mother->getPDG();
      }
      if (pdg == 25 && abs(daughterpdg) == 5) nHbb++;
      if (pdg == 25 && abs(daughterpdg) == 24)nHWW++;
      if (pdg == 25 && abs(daughterpdg) == 21) nHgg++;
      if (pdg == 25 && abs(daughterpdg) == 23) nHZZ++;
      if (pdg == 25 && abs(daughterpdg) == 4) nHcc++;
      if (pdg == 25 && abs(daughterpdg) == 13) nHmumu++;
      if (pdg == 25 && abs(daughterpdg) == 15) nHtautau++;
      if (pdg == 25 && abs(daughterpdg) == 22) nHgamgam++;
      if (abs(pdg) == 24 && abs(motherpdg) == 25 && abs(daughterpdg) > 10 && abs(daughterpdg) < 20) nWlv++;
      if (abs(pdg) == 24 && abs(motherpdg) == 25 && abs(daughterpdg) < 10 && abs(daughterpdg) > 0) nWqq++;
      if (abs(pdg) == 23 && abs(motherpdg) == 25 && (abs(daughterpdg) == 12 || abs(daughterpdg) == 14 || abs(daughterpdg) == 16)) nZvv++;
      if (abs(pdg) == 23 && abs(motherpdg) == 25 && abs(daughterpdg) < 10 && abs(daughterpdg) > 0) nZqq++;
      if (abs(pdg) == 23 && abs(motherpdg) == 25 && abs(daughterpdg) == 11) nZee++;
      if (abs(pdg) == 23 && abs(motherpdg) == 25 && abs(daughterpdg) == 13) nZmm++;
    }

    nHDecay.push_back(nHbb);
    nHDecay.push_back(nHWW);
    nHDecay.push_back(nHgg);
    nHDecay.push_back(nHtautau);
    nHDecay.push_back(nHcc);
    nHDecay.push_back(nHZZ);
    nHDecay.push_back(nHgamgam);
    nHDecay.push_back(nHmumu);
    nHDecay.push_back(nWlv);
    nHDecay.push_back(nWqq);
    nHDecay.push_back(nZvv);
    nHDecay.push_back(nZqq);
    nHDecay.push_back(nZee);
    nHDecay.push_back(nZmm);
    return nHDecay;
  }

  // *******************************************************
  // *******************************************************
  std::vector<Int_t> getNumberOfOverlayEvents(Double_t fEcm, LCCollection *colMC) {
    // return number of overlay event per signal event by MC

    std::vector<Int_t> nOvlEvents;
    Int_t nOvl=0,nBgam=0,nBelectron=0,nBpositron=0;
    Int_t nMCP = colMC->getNumberOfElements();
    for (Int_t i=0;i<nMCP;i++) {
      MCParticle *mcPart = dynamic_cast<MCParticle*>(colMC->getElementAt(i));
      Int_t pdg = mcPart->getPDG();
#if 0
      Int_t ndaughters = mcPart->getDaughters().size();
      Int_t daughterpdg = 0;
      if (ndaughters > 0) {
	MCParticle *daughter = mcPart->getDaughters()[0];
	daughterpdg = daughter->getPDG();
      }
      Int_t nparents = mcPart->getParents().size();
      Int_t motherpdg = 0;
      if (nparents > 0) {
	MCParticle *mother = mcPart->getParents()[0];
	motherpdg = mother->getPDG();
      }
#endif      
      Int_t ibeam = TMath::Abs(mcPart->getEnergy()-fEcm/2) < 0.1 ? 1 : 0;
      Int_t iovl = mcPart->isOverlay() ? 1 : 0;
      if (iovl == 1 && ibeam == 1 && pdg == 22)  nBgam +=1;
      if (iovl == 1 && ibeam == 1 && pdg == 11)  nBelectron +=1;
      if (iovl == 1 && ibeam == 1 && pdg == -11) nBpositron +=1;
    }

    nOvl = (nBgam + nBelectron + nBpositron)/2;
    nOvlEvents.push_back(nOvl);
    nOvlEvents.push_back(nBgam);
    nOvlEvents.push_back(nBelectron);
    nOvlEvents.push_back(nBpositron);

    return nOvlEvents;
  }
  
  // *******************************************************
  // *******************************************************
  void mcDebug(LCCollection *colMC) {
    // list the generator particles

    Int_t nMCP = colMC->getNumberOfElements();
    cerr << setw(6) << "Index" << setw(7) << "PDG" << setw(7) << "Mother" << setw(10) << "Charge" << setw(10) << "Mass"	 << setw(11) << "Energy" << setw(18) << "NumberOfDaughters" << setw(16) << "NumberOfParents" << setw(10) << "Original" << setw(8) << "Overlay" << setw(8) << "FromSim" << endl;
    for (Int_t i=0;i<nMCP;i++) {
      MCParticle *mcPart = dynamic_cast<MCParticle*>(colMC->getElementAt(i));
      Int_t pdg = mcPart->getPDG();
      Int_t nparents = mcPart->getParents().size();
      Int_t motherpdg = 0;
      if (nparents > 0) {
	MCParticle *mother = mcPart->getParents()[0];
	motherpdg = mother->getPDG();
      }
      Double_t charge = mcPart->getCharge();
      Double_t mass = mcPart->getMass();
      Double_t energy = mcPart->getEnergy();
      TVector3 pv = TVector3(mcPart->getMomentum());
      Int_t ndaughters = mcPart->getDaughters().size();
#if 0
      Int_t daughterpdg = 0;
      if (ndaughters > 0) {
	MCParticle *daughter = mcPart->getDaughters()[0];
	daughterpdg = daughter->getPDG();
      }
#endif      
      TLorentzVector lortz = TLorentzVector(pv,energy);
      Int_t originalPDG = getOriginalPDG(mcPart);
      Int_t ioverlay = mcPart->isOverlay()? 1 : 0;
      Int_t icreatedinsim = mcPart->isCreatedInSimulation()? 1 : 0;
      cerr << setw(6) << i << setw(7) << pdg << setw(7) << motherpdg << setw(10) << setprecision(3) << charge << setw(10) << mass << setw(11) << energy << setw(18) << ndaughters << setw(16) << nparents << setw(10) << originalPDG << setw(8) << ioverlay << setw(8) << icreatedinsim << endl;
    }

  }

  // *******************************************************
  // *******************************************************
  TLorentzVector getLorentzEcm(Double_t fEcm) {
    const Double_t fCrossAngle = 0.014;
    return TLorentzVector(fEcm*TMath::Sin(fCrossAngle/2),0.,0.,fEcm);
  }
  TLorentzVector getLorentzEcm(Double_t fEcm,Bool_t isCrossingAngle) {
    const Double_t fCrossAngle = 0.014;
    if (isCrossingAngle) {
      return TLorentzVector(fEcm*TMath::Sin(fCrossAngle/2),0.,0.,fEcm);
    }
    else {
      return TLorentzVector(0.,0.,0.,fEcm);
    }
  }

  // *******************************************************
  // *******************************************************
  Double_t getRecoilMass(TLorentzVector lortzEcm, TLorentzVector lortzZ) {
    TLorentzVector lortzRecoil = lortzEcm - lortzZ;
    return lortzRecoil.M();
  }

  // *******************************************************
  // *******************************************************
  Double_t getAcoPlanarity(TLorentzVector lortz1, TLorentzVector lortz2) {
    Double_t phi1 = lortz1.Phi();
    Double_t phi2 = lortz2.Phi();
    Double_t acop12 = TMath::Abs(phi1-phi2);
    return acop12 > TMath::Pi() ? TMath::Pi()*2. - acop12 : acop12;
  }

  // *******************************************************
  // *******************************************************
  TLorentzVector getLorentzVector(ReconstructedParticle *pfo) {
     return TLorentzVector(pfo->getMomentum(),pfo->getEnergy());
  }

  // *******************************************************
  // *******************************************************
  Double_t getCosTheta(ReconstructedParticle *part1, ReconstructedParticle *part2) {
    TVector3 momentum1 = TVector3(part1->getMomentum());
    TVector3 momentum2 = TVector3(part2->getMomentum());
    Double_t cosTheta = momentum1.Dot(momentum2)/momentum1.Mag()/momentum2.Mag();
    return cosTheta;
  }

  // *******************************************************
  // *******************************************************
  Double_t getCosTheta(TLorentzVector part1, TLorentzVector part2) {
    TVector3 momentum1 = TVector3(part1.Vect());
    TVector3 momentum2 = TVector3(part2.Vect());
    Double_t cosTheta = momentum1.Dot(momentum2)/momentum1.Mag()/momentum2.Mag();
    return cosTheta;
  }

} // namespace mylib
