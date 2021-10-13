// *****************************************************
// Processor for training isolated leptons selection
//                        ----Junping
// *****************************************************
#include "IsoLepTrainingProcessor.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/Cluster.h>
#include <UTIL/LCTypedVector.h>
#include <EVENT/Track.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/ParticleID.h>
#include <marlin/Exceptions.h>
#include <EVENT/Vertex.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNtupleD.h"
#include "TVector3.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "Utilities.h"

//#define __SINGLEP__

using namespace lcio ;
using namespace marlin ;
using namespace std;

using namespace isolep;

IsoLepTrainingProcessor aIsoLepTrainingProcessor ;


IsoLepTrainingProcessor::IsoLepTrainingProcessor() : Processor("IsoLepTrainingProcessor") {
  
  // modify processor description
  _description = "IsoLepTrainingProcessor does whatever it does ..." ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::MCPARTICLE,
			   "InputMCParticlesCollection" , 
			   "Name of the MCParticle collection"  ,
			   _colMCP ,
			   std::string("MCParticlesSkimmed") ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "InputMCTruthLinkCollection" , 
			   "Name of the MCTruthLink collection"  ,
			   _colMCTL ,
			   std::string("RecoMCTruthLink") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputPandoraPFOsCollection" , 
			   "Name of the PandoraPFOs collection"  ,
			   _colPFOs ,
			   std::string("PandoraPFOs") ) ;

  registerInputCollection( LCIO::VERTEX,
			   "InputPrimaryVertexCollection" , 
			   "Name of the Primary Vertex collection"  ,
			   _colPVtx ,
			   std::string("PrimaryVertex") ) ;
  
  registerProcessorParameter("IsLepTune",
			   "Is lepton tune?"  ,
			   _is_lep_tune ,
			   bool(true) ) ;

  registerProcessorParameter("MCDebugging" ,
			     "set true if you want to check generator information",
			     _mcdebug,
			     bool(false)
			     );

  registerProcessorParameter("IsForSignal",
			   "Is for signal?"  ,
			   _is_for_sig ,
			   bool(true) ) ;

  registerProcessorParameter("IsoLepType",
			   "Isolated Lepton Type"  ,
			   _iso_lep_type ,
			   int(13) ) ;
}

void IsoLepTrainingProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void IsoLepTrainingProcessor::processRunHeader( LCRunHeader* ) { 

  _nRun++ ;
} 

void IsoLepTrainingProcessor::processEvent( LCEvent * evt ) { 

    
  // this gets called for every event 
  // usually the working horse ...
  _nEvt++;

#if 1
  Double_t fEtrackCut = -1.;        // lower edge of each PFO energy 
#else
  Double_t fEtrackCut = 0.05;        // lower edge of each PFO energy 
#endif
  Double_t fCosConeCut = 0.98;   // the angle of cone around the direction of pfo
  Double_t fCosLargeConeCut = 0.95; // angel of large cone around the pfo

  TDirectory *last = gDirectory;
  gFile->cd("/");

  cerr << endl << "Hello, MVA Lepton Training! Event No. " << _nEvt << endl;

  static TNtupleD *hGen = 0;
  if (!hGen) {
    stringstream tupstr_gen;
    tupstr_gen << "nhbb:pdg:npvt:zipmc:xipmc:yipmc:zipvf:xipvf:yipvf:nlepmc:ievt:seriallep:iov:xerrip:yerrip:zerrip:ntrksip" << ":"
	       << "cosmc:phimc:pmc"
	       << ends;
    hGen = new TNtupleD("hGen","",tupstr_gen.str().data());
  }
  static TNtupleD *hLep = 0;
  if (!hLep) {
    stringstream tupstr_lep;
    tupstr_lep << "pdg:nlepmc:ievt:seriallep:iov:orig"
	       << ends;
    hLep = new TNtupleD("hLep","",tupstr_lep.str().data());
  }

  // -- Get the MCTruth Linker --
  LCCollection *colMCTL = evt->getCollection(_colMCTL);
  LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);

#ifndef __SINGLEP__
  // -- Read out MC information --  
  LCCollection *colMC = evt->getCollection(_colMCP);
  if (!colMC) {
    std::cerr << "No MC Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }
  // get the truth information
  Int_t nMCP = colMC->getNumberOfElements();
#endif
  Int_t pdgLepMC = -1,serialLepMC = -1,nMCLep=0,iOvlLepMC=-999;
  Double_t z_IPMC=999.,x_IPMC=999.,y_IPMC=999.;
  std::vector<Int_t> pdgLepMCv,serialLepMCv;
  TLorentzVector lortzLepMC;
  pdgLepMC = _iso_lep_type;
#ifndef __SINGLEP__
  for (Int_t i=0;i<nMCP;i++) {
    MCParticle *mcPart = dynamic_cast<MCParticle*>(colMC->getElementAt(i));
    Int_t pdg = mcPart->getPDG();
    Double_t energy = mcPart->getEnergy();
    TVector3 pv = TVector3(mcPart->getMomentum());
    TLorentzVector lortz = TLorentzVector(pv,energy);
    //    Int_t orig = getOriginalPDGForIsoLep(mcPart);
    Int_t orig = getOriginalPDGForIsoLep(mcPart,colMC);    
    // Int_t status = mcPart->getGeneratorStatus();
    if (i>=6 && i<=11) {
      if (TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13 || TMath::Abs(pdg) == 15) {
	pdgLepMC = pdg;
	pdgLepMCv.push_back(pdgLepMC);
      }
    }
    // Int_t origw = pdgLepMC > 0? -24 : 24;
    //    if ((TMath::Abs(pdg)==pdgLepMC) && (orig==24||orig==pdgLepMC) && status==1) {
    //    if (pdg==pdgLepMC && (orig==origw||orig==pdgLepMC) && status==1 && TMath::Abs(pdg)==_iso_lep_type && !(mcPart->isOverlay())) {
      //    if (status == 1 && TMath::Abs(pdg)==_iso_lep_type && TMath::Abs(orig) == 13 && !(mcPart->isOverlay())) {
    if (i==0) {
      serialLepMC=i;
      serialLepMCv.push_back(i);
      nMCLep++;
      iOvlLepMC = mcPart->isOverlay() ? 1 : 0;
      lortzLepMC = lortz;
      Double_t data_lep[20];
      data_lep[ 0] = pdg;
      data_lep[ 1] = nMCLep;
      data_lep[ 2] = _nEvt;
      data_lep[ 3] = serialLepMC;
      data_lep[ 4] = iOvlLepMC;
      data_lep[ 5] = orig;
      hLep->Fill(data_lep);
    }
    if (i == 0) {
      TVector3 vip = TVector3(mcPart->getVertex());
      z_IPMC = vip[2];
      x_IPMC = vip[0];
      y_IPMC = vip[1];
    }
  }

  if (_mcdebug) {
    mcDebug(colMC);
  }

  // primary vertex from MCTruth: z_IPMC

  // primary vertex from VertexFinder (LCFIPlus)
  LCCollection *colPVtx = evt->getCollection(_colPVtx);
  Int_t npvtx = colPVtx->getNumberOfElements();
  Vertex *pvtx = dynamic_cast<Vertex*>(colPVtx->getElementAt(0));
  TVector3 v_pvtx = TVector3(pvtx->getPosition());
  Double_t z_pvtx = v_pvtx[2];
  Double_t x_pvtx = v_pvtx[0];
  Double_t y_pvtx = v_pvtx[1];
  Double_t xerr_pvtx = TMath::Sqrt(pvtx->getCovMatrix()[0]);
  Double_t yerr_pvtx = TMath::Sqrt(pvtx->getCovMatrix()[2]);
  Double_t zerr_pvtx = TMath::Sqrt(pvtx->getCovMatrix()[5]);

  ReconstructedParticle *pvtx_rp = pvtx->getAssociatedParticle();
  Int_t ntrks_pvtx = pvtx_rp->getParticles().size();
  
  std::vector<Int_t> nHDecay;
  nHDecay = getHiggsDecayModes(colMC);
  Double_t nHbb = nHDecay[0]; // tag H---> b b
  Double_t data_gen[20];
  data_gen[ 0] = nHbb;
  data_gen[ 1] = pdgLepMC;
  data_gen[ 2] = npvtx;
  data_gen[ 3] = z_IPMC;
  data_gen[ 4] = x_IPMC;
  data_gen[ 5] = y_IPMC;
  data_gen[ 6] = z_pvtx;  
  data_gen[ 7] = x_pvtx;  
  data_gen[ 8] = y_pvtx;  
  data_gen[ 9] = nMCLep;  
  data_gen[10] = _nEvt;
  data_gen[11] = serialLepMC;
  data_gen[12] = iOvlLepMC;
  data_gen[13] = xerr_pvtx;
  data_gen[14] = yerr_pvtx;
  data_gen[15] = zerr_pvtx;
  data_gen[16] = ntrks_pvtx;
  data_gen[17] = lortzLepMC.CosTheta();
  data_gen[18] = lortzLepMC.Phi();
  data_gen[19] = lortzLepMC.P();
  hGen->Fill(data_gen);
#endif
  //  return;
  if (_is_for_sig && nMCLep == 0)  {
    std::cerr << "This is not an signal event with isolated muon or electron" << std::endl;
    throw marlin::SkipEventException(this);
  }

  
  // -- Read out PFO information --
  LCCollection *colPFO = evt->getCollection(_colPFOs);
  if (!colPFO) {
    std::cerr << "No PFO Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }
  Int_t nPFOs = colPFO->getNumberOfElements();

  static TNtupleD *hPfo = 0;
  if (!hPfo) {
    stringstream tupstr_pfo;
    tupstr_pfo << "ntracks:charge:mcpdg:motherpdg:deltae:mmotherpdg:ndaughters" << ":"
	       << "mcoriginal:energy:type:pid"                                  << ":"
	       << "totalcalenergy:momentum:ecalenergy:hcalenergy:coneenergy"    << ":"
	       << "nmctl:mcwgt:ievt:irun"                                       << ":"
	       << "nhits:ncones:nconechg:nconeneu:coneec:coneen:energylink"     << ":"
	       << "costheta:yokeenergy:energycor:momentumcor"                   << ":"
	       << "d0:z0:r0:deltad0:deltaz0:nsigd0:nsigz0:nsigr0:iov"           << ":"
	       << "coslarcon:energyratio:nphoton:ratioecal:ratiototcal"         << ":"
	       << "isim:pdgmc:orig:isorig:zipvf:zipmc"                          << ":"
	       << "x0:y0:xipmc:yipmc:xipvf:yipvf:mcserial:serialLepMC"          << ":"
	       << "cosmc:phimc:pmc"
	       << ends;
    hPfo = new TNtupleD("hPfo","",tupstr_pfo.str().data());
  }

  // loop all the PFOs
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(recPart);
    FloatVec vecWgtMCTL = navMCTL->getRelatedToWeights(recPart);
    Int_t mcpdg,motherpdg,mmotherpdg;
    Double_t mcwgt=0.;
    mcpdg = 0;
    motherpdg = -99999;
    mmotherpdg = -99999;
    Double_t deltaE = -99999.;
    Double_t energyLink = -99999.;
    Int_t mcoriginal = 0;
    Int_t mcoriginalIsoLep = 0;
    Int_t mcndaughters = 0;
    Int_t nMCTL = vecMCTL.size();
    Int_t iOverlay = 0;
    Int_t iCreatedInSim = 0;
    Int_t mcserial = -1;
#ifndef __SINGLEP__
    if (vecMCTL.size() > 0) {
      //      MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[0]);
      MCParticle *mcPart = getMCParticle(recPart,colMCTL,mcwgt);      
      if (!mcPart) continue;
      if (mcPart->isOverlay()) iOverlay = 1;
      if (mcPart->isCreatedInSimulation()) iCreatedInSim = 1;
      mcpdg = mcPart->getPDG();
      //      mcwgt = vecWgtMCTL[0];
      deltaE = mcPart->getEnergy()-recPart->getEnergy();
      energyLink = mcPart->getEnergy();
      mcoriginal = getOriginalPDG(mcPart);
      mcoriginalIsoLep = getOriginalPDGForIsoLep(mcPart);
      //      mcoriginalIsoLep = getOriginalPDGForIsoLep(mcPart,colMC);      
      motherpdg = 0;
      mcndaughters = mcPart->getDaughters().size();
      if (mcPart->getParents().size() != 0) {
	MCParticle *motherPart = mcPart->getParents()[0];
	motherpdg = motherPart->getPDG();
	mmotherpdg = 0;
	if (motherPart->getParents().size() != 0) {
	  MCParticle *mmotherPart = motherPart->getParents()[0];
	  mmotherpdg = mmotherPart->getPDG();
	}
      }
      mcserial = getMCSerial(mcPart,colMC);
    }
    //    Int_t mcserial = getMCSerial(recPart,colMCTL,colMC);
    //    Int_t isOrigLep = mcserial==serialLepMC? 1:0;
    Int_t isOrigLep = 0;
    for (size_t j=0;j<serialLepMCv.size();j++) {
      if (mcserial == serialLepMCv[j]) {
	isOrigLep = 1;
	break;
      }
    }
#else
    Int_t isOrigLep = 0;
#endif
    Double_t energy = recPart->getEnergy();
    Double_t charge = recPart->getCharge();
    Int_t itype = recPart->getType();
    Int_t pid = 0;
    TrackVec tckvec = recPart->getTracks();
    Int_t ntracks = tckvec.size();
    Double_t d0=0.,z0=0.,deltad0=0.,deltaz0=0.,nsigd0=0.,nsigz0=0.;
    Double_t phi0=0.,x0=0.,y0=0.;
    if (ntracks > 0) {
      d0 = tckvec[0]->getD0();
      z0 = tckvec[0]->getZ0();
      phi0 = tckvec[0]->getPhi();
      x0 = -d0*TMath::Sin(phi0);
      y0 = d0*TMath::Cos(phi0);
#ifndef __SINGLEP__
      z0 -= z_pvtx;
#endif      
      deltad0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[0]);
      deltaz0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[9]);
      nsigd0 = d0/deltad0;
      nsigz0 = z0/deltaz0;
    }
    Double_t r0 = TMath::Sqrt(d0*d0+z0*z0);
    Double_t nsigr0 = TMath::Sqrt(nsigd0*nsigd0+nsigz0*nsigz0);
    Double_t data[100];
    data[0] = ntracks;
    data[1] = charge;
    data[2] = mcpdg;
    data[3] = motherpdg;
    data[4] = deltaE;
    data[5] = mmotherpdg;
    data[6] = mcndaughters;
    data[7] = mcoriginal;
    data[8] = energy;
    data[9] = itype;
    data[10]= pid;
    if (energy > fEtrackCut) {
      Double_t ecalEnergy = 0;
      Double_t hcalEnergy = 0;
      Double_t yokeEnergy = 0;
      Double_t totalCalEnergy = 0;
      Int_t nHits = 0;
      std::vector<lcio::Cluster*> clusters = recPart->getClusters();
      for (std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();iCluster!=clusters.end();++iCluster) {
	ecalEnergy += (*iCluster)->getSubdetectorEnergies()[0];
	hcalEnergy += (*iCluster)->getSubdetectorEnergies()[1];
	yokeEnergy += (*iCluster)->getSubdetectorEnergies()[2];
	ecalEnergy += (*iCluster)->getSubdetectorEnergies()[3];
	hcalEnergy += (*iCluster)->getSubdetectorEnergies()[4];
	CalorimeterHitVec calHits = (*iCluster)->getCalorimeterHits();
	//	nHits += (*iCluster)->getCalorimeterHits().size();
	nHits = calHits.size();
      }
      totalCalEnergy = ecalEnergy + hcalEnergy;
      TVector3 momentum = TVector3(recPart->getMomentum());
      Double_t momentumMagnitude = momentum.Mag();
      Double_t cosTheta = momentum.CosTheta();
      //get cone information
      //      std::vector<lcio::ReconstructedParticle*> conePFOs;
      //      Double_t coneEnergy = getConeEnergy(recPart,colPFO,fCosConeCut,conePFOs);
      Bool_t woFSR = kTRUE;
      Double_t coneEnergy0[3] = {0.,0.,0.};
      Double_t pFSR[4] = {0.,0.,0.,0.};
      Double_t pLargeCone[4]  = {0.,0.,0.,0.};
      Int_t nConePhoton = 0;
      //      getConeEnergy(recPart,colPFO,fCosConeCut,woFSR,coneEnergy0,pFSR);
      getConeEnergy(recPart,colPFO,fCosConeCut,woFSR,coneEnergy0,pFSR,fCosLargeConeCut,pLargeCone,nConePhoton);
      Double_t coneEnergy = coneEnergy0[0];
      Double_t coneEN     = coneEnergy0[1];
      Double_t coneEC     = coneEnergy0[2];
      TLorentzVector lortzFSR = TLorentzVector(pFSR[0],pFSR[1],pFSR[2],pFSR[3]);
      TLorentzVector lortzLargeCone = TLorentzVector(pLargeCone[0],pLargeCone[1],pLargeCone[2],pLargeCone[3]);
      TVector3 momentumLargeCone = lortzLargeCone.Vect();
      Double_t cosThetaWithLargeCone = 1.;
      if (momentumLargeCone.Mag() > 0.0000001) {
	cosThetaWithLargeCone = momentum.Dot(momentumLargeCone)/momentumMagnitude/momentumLargeCone.Mag();
      }
      Double_t energyRatioWithLargeCone = energy/(energy+lortzLargeCone.E());
      Double_t energyCorr = energy + lortzFSR.E();
      TVector3 momentumCorr = momentum + TVector3(lortzFSR.Px(),lortzFSR.Py(),lortzFSR.Pz());
      Double_t momentumMagCorr = momentumCorr.Mag();
      Double_t ratioECal = 0., ratioTotalCal = 0.;
      if (ecalEnergy > 0.) ratioECal = ecalEnergy/totalCalEnergy;
      //      ratioTotalCal = totalCalEnergy/momentumMagnitude;
      ratioTotalCal = totalCalEnergy/(momentumMagnitude+0.00000001);
      //      Int_t nConePFOs    = conePFOs.size();
      Int_t nConePFOs    = 0;
      Int_t nConeCharged = 0;
      Int_t nConeNeutral = 0;
      // save the pfo information
      data[11] = totalCalEnergy;
      data[12] = momentumMagnitude;
      data[13] = ecalEnergy;
      data[14] = hcalEnergy;
      data[15] = coneEnergy;
      data[16] = nMCTL;
      data[17] = mcwgt;
      data[18] = _nEvt;
      data[19] = _nRun;
      data[20] = nHits;
      data[21] = nConePFOs;
      data[22] = nConeCharged;
      data[23] = nConeNeutral;
      data[24] = coneEC;
      data[25] = coneEN;
      data[26] = energyLink;
      data[27] = cosTheta;
      data[28] = yokeEnergy;
      data[29] = energyCorr;
      data[30] = momentumMagCorr;
      data[31] = d0;
      data[32] = z0;
      data[33] = r0;
      data[34] = deltad0;
      data[35] = deltaz0;
      data[36] = nsigd0;
      data[37] = nsigz0;
      data[38] = nsigr0;
      data[39] = iOverlay;
      data[40] = cosThetaWithLargeCone;
      data[41] = energyRatioWithLargeCone;
      data[42] = nConePhoton;
      data[43] = ratioECal;
      data[44] = ratioTotalCal;
      data[45] = iCreatedInSim;
      data[46] = pdgLepMC;
      data[47] = mcoriginalIsoLep;
      data[48] = isOrigLep;
      data[49] = z_pvtx;
      data[50] = z_IPMC;      
      data[51] = x0;      
      data[52] = y0;
      data[53] = x_IPMC;  
      data[54] = y_IPMC;  
      data[55] = x_pvtx;  
      data[56] = y_pvtx;  
      data[57] = mcserial;  
      data[58] = serialLepMC;  
      data[59] = lortzLepMC.CosTheta();
      data[60] = lortzLepMC.Phi();
      data[61] = lortzLepMC.P();
      if (_is_lep_tune) {
	if(_is_for_sig) {
	  if (isOrigLep==1) hPfo->Fill(data);
	}
	else {
	  hPfo->Fill(data);
	}
      }
    }
  }

  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
  		       << "   in run:  " << evt->getRunNumber() 
  		       << std::endl ;

  //  _nEvt ++ ;

  last->cd();
}



void IsoLepTrainingProcessor::check( LCEvent * ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void IsoLepTrainingProcessor::end(){ 

  cerr << "IsoLepTrainingProcessor::end()  " << name() 
       << " processed " << _nEvt << " events in " << _nRun << " runs "
       << endl ;
  
}
