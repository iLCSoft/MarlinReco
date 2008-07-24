#include "PFOID.h"
#include <iostream>
#include <fstream>
#include "math.h"

#include "HelixClass.h"
#include "ClusterShapes.h"
// GEAR include files
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/BField.h>

#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/Cluster.h>

#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/LCCollectionVec.h>
 
#include "PDF.h"

#ifdef root_out 
#include "TFile.h"
#include "TH1D.h"
#endif

using namespace lcio ;
using namespace marlin ;

PFOID aPFOID ;

PFOID::PFOID() : Processor("PFOID") {
  _description = "PFOID looks for reconstructed particles and checks wether there are muons" ;

  registerProcessorParameter("RecoParticleCollection","Name of the by Wolf reconstructed Particle collection",_recoCol,std::string("RecoParticles")) ;

  registerProcessorParameter("NewRecoParticleCollection","Name of the new reconstructed Particle collection",_newrecoCol,std::string("NewRecoParticles")) ;

  registerProcessorParameter("FilePDFName","Name of file containing pdfs",_filename,std::string("pdf.txt")) ;

  registerProcessorParameter("neutralFilePDFName","Name of file containing pdfs for neutral particle",_filename1,std::string("npdf.txt")) ;
};


void PFOID::init() {
  printParameters() ;
  _nRun = 0 ;
  _nEvt = 0 ;
  std::cout << " ************************** " << std::endl ;
  std::cout << "         PFOID            " << std::endl ;
  std::cout << " ************************** " << std::endl ;

  pdf = new PDF(_filename.c_str());
  npdf = new PDF(_filename1.c_str());

  noClusterParticle=0;

  // temporary
#ifdef root_out
  file = new TFile("likelihood.root","RECREATE");
  his_e = new TH1D("lh_e","Likelihood electron",50,0.,1.);
  his_mu = new TH1D("lh_mu","Likelihood muon",50,0.,1.);
  his_pi = new TH1D("lh_pi","Likelihood pion",50,0.,1.);
  his_g = new TH1D("lh_g","Likelihood gamma",50,0.,1.);
  his_nh = new TH1D("lh_nh","Likelihood neutral hadron",50,0.,1.);
#endif
};


void PFOID::processRunHeader( LCRunHeader * run ) {
  _nRun++ ;
};


void PFOID::processEvent( LCEvent * evt ) {
  std::cout << "%% PFOID: Event " << _nEvt << std::endl;
  const gear::BField& gearBField = Global::GEAR->getBField();
  _bField = float(gearBField.at(gear::Vector3D(0.0,0.0,0.0)).z()); //BField should be constant and point in z direction
  //const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  //_bField = float(gearTPC.getDoubleVal("BField"));

  LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  try{

    LCCollection *col = evt->getCollection(_recoCol.c_str());
    int nRecos = col->getNumberOfElements() ;
    for(int i=0; i<nRecos; i++){  // over all reco. particles
      ReconstructedParticle *recopart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));

      fill_info(i,recopart); // fills also pdf->VO or npdf->VO

      ReconstructedParticleImpl *rcpi = new ReconstructedParticleImpl();
      rcpi->setType(recopart->getType());
      rcpi->setMomentum(recopart->getMomentum());
      rcpi->setEnergy(recopart->getEnergy());
      rcpi->setCharge(recopart->getCharge());
      rcpi->setMass(recopart->getMass());
      rcpi->setCovMatrix(recopart->getCovMatrix());
      rcpi->setReferencePoint(recopart->getReferencePoint());

      ClusterVec clv = recopart->getClusters();
      for(unsigned int cls=0; cls<clv.size(); cls++)
	rcpi->addCluster(clv[cls]);
      TrackVec trv = recopart->getTracks();
      for(unsigned int trk=0; trk<trv.size(); trk++)
	rcpi->addTrack(trv[trk]);
      ReconstructedParticleVec rpv = recopart->getParticles();
      for(unsigned int rps=0; rps<rpv.size(); rps++)
	rcpi->addParticle(rpv[rps]);


      if(clv.size()>0){
	if(info.withTrack){

	  // potential for automizing
	  std::cout << " RecoParticle " << i << std::endl;
	  std::cout << "  likelihood(electron) : ";
	  std::cout << pdf->GetLikelihood("electron") << std::endl;
	  std::cout << "  likelihood(myon)     : ";
	  std::cout << pdf->GetLikelihood("muon") << std::endl;
	  std::cout << "  likelihood(pion)     : ";
	  std::cout << pdf->GetLikelihood("pion") << std::endl;
	  
	  
	  ParticleIDImpl *pIDe = new ParticleIDImpl();
	  
	  double maxv=0., lhd;
	  
	  // PID to recopart: electron
	  if(rcpi->getCharge()<0){
	    pIDe->setPDG(11);
	  }else{
	    pIDe->setPDG(-11);
	  }
	  pIDe->setType(0);
	  lhd=pdf->GetLikelihood("electron");
	  pIDe->setLikelihood((float)lhd);
	  pIDe->setAlgorithmType(0);
	  rcpi->addParticleID(pIDe);
	  maxv=lhd;
	  rcpi->setType(0);
	  rcpi->setMass(0.000511);
	  rcpi->setParticleIDUsed(pIDe);
	  rcpi->setGoodnessOfPID((float)lhd);
	  
	  ParticleIDImpl *pIDm = new ParticleIDImpl();
	  // PID to recopart: muon
	  if(rcpi->getCharge()<0){
	    pIDm->setPDG(13);
	  }else{
	    pIDm->setPDG(-13);
	  }
	  pIDm->setType(1);
	  lhd=pdf->GetLikelihood("muon");
	  pIDm->setLikelihood((float)lhd);
	  pIDm->setAlgorithmType(0);
	  rcpi->addParticleID(pIDm);
	  if(lhd>maxv){
	    maxv=lhd;
	    rcpi->setType(1);
	    rcpi->setMass(0.1056);
	    rcpi->setParticleIDUsed(pIDm);
	    rcpi->setGoodnessOfPID((float)lhd);
	  }
	  
	  ParticleIDImpl *pIDh = new ParticleIDImpl();
	  // PID to recopart: charged hadron
	  if(rcpi->getCharge()<0.){
	    pIDh->setPDG(-211);
	  }else{
	    pIDh->setPDG(211);
	  }
	  pIDh->setType(2);
	  lhd=pdf->GetLikelihood("pion");
	  pIDh->setLikelihood((float)lhd);
	  pIDh->setAlgorithmType(0);
	  rcpi->addParticleID(pIDh);
	  if(lhd>maxv){
	    rcpi->setMass(0.1396);
	    rcpi->setType(2);
	    maxv=lhd;
	    rcpi->setParticleIDUsed(pIDh);
	    rcpi->setGoodnessOfPID((float)lhd);	  
	  }
	  
#ifdef root_out
	  his_e->Fill(pdf->GetLikelihood("electron"));
	  his_mu->Fill(pdf->GetLikelihood("muon"));
	  his_pi->Fill(pdf->GetLikelihood("pion"));
#endif
	}else{

	  std::cout << " RecoParticle " << i << std::endl;
	  std::cout << "  likelihood(photon) : " ;
	  std::cout << npdf->GetLikelihood("photon") << std::endl;
	  std::cout << "  likelihood(kaon/neutron)     : " ;
	  std::cout << npdf->GetLikelihood("kaon/neutron") << std::endl;
	  
	  
	  ParticleIDImpl *pIDg = new ParticleIDImpl();
	  double maxv=0., lhd;
	  // PID to recopart: gamma
	  pIDg->setPDG(22);
	  pIDg->setType(3);
	  lhd=npdf->GetLikelihood("photon");
	  pIDg->setLikelihood((float)lhd);
	  pIDg->setAlgorithmType(0);
	  rcpi->addParticleID(pIDg);
	  maxv=lhd;
	  rcpi->setMass(0.0);
	  rcpi->setType(3);
	  rcpi->setParticleIDUsed(pIDg);
	  rcpi->setGoodnessOfPID((float)lhd);
	  
	  ParticleIDImpl *pIDnh = new ParticleIDImpl();
	  // PID to recopart: neutral hadron
	  pIDnh->setPDG(130);
	  pIDnh->setType(4);
	  lhd=npdf->GetLikelihood("kaon/neutron");
	  pIDnh->setLikelihood((float)lhd);
	  pIDnh->setAlgorithmType(0);
	  rcpi->addParticleID(pIDnh);
	  if(lhd>maxv){
	    maxv=lhd;
	    rcpi->setType(4);
	    rcpi->setMass(0.495);
	    rcpi->setParticleIDUsed(pIDnh);
	    rcpi->setGoodnessOfPID((float)lhd);
	  }

	  std::cout << "  EtoN_ecal : " << info.EtoN_ecal;
	  std::cout << "  EtoN_hcal : " << info.EtoN_hcal;
	  std::cout << "  dmean : " << info.dmean << std::endl;
	  std::cout << "  ex : " << info.ex;
	  std::cout << "  L1 : " << info.L1;
	  std::cout << "  L2 : " << info.L2;
	  std::cout << "  L3 : " << info.L3 << std::endl;
	  
#ifdef root_out
	  his_g->Fill(npdf->GetLikelihood("photon"));
	  his_nh->Fill(npdf->GetLikelihood("kaon/neutron"));
#endif
	}// if ( withTrack ...
      
      }else{

	ParticleIDImpl *pID = new ParticleIDImpl();
	pID->setPDG(666);
	pID->setType(5);
	double lhd=1;
	pID->setLikelihood((float)lhd);
	pID->setAlgorithmType(0);
	rcpi->addParticleID(pID);
	rcpi->setType(5);
	rcpi->setMass(0.0);
	rcpi->setParticleIDUsed(pID);
	rcpi->setGoodnessOfPID((float)lhd);
      }

      recparcol->addElement( rcpi );
    }// for i ...
  }catch(DataNotAvailableException &e){ }

  evt->addCollection( recparcol , _newrecoCol.c_str());

  std::cout << "%% PFOID: End of PFOID\n" << std::endl; 
  _nEvt++ ;

};




void PFOID::check( LCEvent * evt ) {};




void PFOID::end() {
#ifdef root_out
  file->Write();
  file->Close();
#endif
  delete pdf;
  delete npdf;
};




void PFOID::init_info(){
  info.px=0.0;
  info.py=0.0;
  info.pz=0.0;
  info.Eecal=0.0;
  info.Ehcal=0.0;
  info.ex=0.0;
  info.dmean=0.0;
  info.Edmean=0.0;
  info.EtoN_ecal = 0.0;
  info.EtoN_hcal = 0.0;
  info.EecalToEtot = 0.0;
  info.Necal=0;
  info.Nhcal=0;
  info.L1=-1;
  info.L2=-1;
  info.L3=-1;
  info.EL1=100;
  info.EL2=100;
  info.EL3=100;
  info.withTrack=false;
}



void PFOID::fill_info(int i, ReconstructedParticle *rp){
  bool noCluster=false;

  init_info();

  HelixClass * helix = new HelixClass();
  TrackVec tv = rp->getTracks();
  ClusterVec cv = rp->getClusters();
  std::cout << " RP: " << i << " #Tracks: " << tv.size() << " #Cluster: " << cv.size() << std::endl;

  // Cluster Layers and energies
  int nClusters = cv.size();

  float *ca, *cx ,*cy, *cz;
  int N, nHits;
  float Cpos[3] = {0., 0., 0.}, Esum=0;
  
  if(nClusters>0){
    for(int j=0; j<nClusters; j++){  // over each cluster of particle
      Cluster *cl = cv[j];
      info.Eecal += (double)cl->getSubdetectorEnergies()[0];
      info.Ehcal += (double)cl->getSubdetectorEnergies()[1];
      CalorimeterHitVec chv = cl->getCalorimeterHits();
      int nCalHits = chv.size();
      
      for(int k=0; k<nCalHits; k++){  // over each hit in cluster
	int cellid = chv[k]->getCellID0();
	int layer = cellid >> 24;
	if(chv[k]->getType()==1){
	  info.Nhcal +=1 ;
	  if(layer>=info.L1){
	    info.L3=info.L2 ;
	    info.L2=info.L1 ;
	    info.L1=(double)layer ;
	  }else{
	    if(layer>=info.L2){
	      info.L3=info.L2 ;
	      info.L2=(double)layer;
	    }else{
	      if(layer>info.L3)
		info.L3=(double)layer;
	    }
	  }
	}else{
	  info.Necal += 1;
	  if(layer<=info.EL1){
	    info.EL3=info.EL2 ;
	    info.EL2=info.EL1 ;
	    info.EL1=(double)layer ;
	  }else{
	    if(layer<=info.EL2){
	      info.EL3=info.EL2 ;
	      info.EL2=(double)layer;
	    }else{
	      if(layer<info.EL3)
		info.EL3=(double)layer;
	    }
	  }
	}
      }// for k ...
    }// for j ...
    
    // ClusterShape and postion
    N=(int)(info.Necal+info.Nhcal);
    ca = new float[N]; // for cluster shape
    cx = new float[N]; //      - " -
    cy = new float[N]; //      - " -
    cz = new float[N]; //      - " -
    nHits=0;  // counts all calorimeter hits (N)
        
    for(int j=0; j<nClusters; j++){  // over each cluster of particle
      Cluster *cl = cv[j];
      CalorimeterHitVec chv = cl->getCalorimeterHits();
      int nCalHits = chv.size();
      
      for(int kk=0; kk<3; kk++)
	Cpos[kk] += cl->getPosition()[kk]*cl->getEnergy() ;

      Esum += cl->getEnergy();
      
      for(int k=0; k<nCalHits; k++){  // over each hit in cluster
	if(nHits>=N){
	  std::cerr << "%% PFOID: Error in Calorimeter hits loop!" << std::endl;
	  break;
	}
	ca[nHits]=chv[k]->getEnergy();
	cx[nHits]=chv[k]->getPosition()[0];
	cy[nHits]=chv[k]->getPosition()[1];
	cz[nHits]=chv[k]->getPosition()[2];
	nHits++;
      }// for k ...
    }// for j ...
    for(int kk=0; kk<3; kk++) Cpos[kk] /= Esum;
    
    ClusterShapes *clsp = new ClusterShapes(N,ca,cx,cy,cz);
    
    float width = clsp->getWidth();
    float evi[3];
    for(int kk=0; kk<3; kk++)  evi[kk] = clsp->getEigenValInertia()[kk];
    
    float max = evi[0];
    if(evi[1]>max) max = evi[1];
    if(evi[2]>max) max = evi[2];
    info.ex = (double)width/(double)max; // excentricity
  }else{
    noCluster=true;
    std::cout << " WARNING!! This particle has no cluster" << std::endl;
    noClusterParticle++;
  }// if(nCluster>0)

  // Tracks
  float d0, fi, omega, z0, tanlambda;
  if(tv.size()>0){  // if tracks available
    info.withTrack=true;
    Track * track = NULL;
    unsigned int nOTH=0;
    for(unsigned int tri=0; tri<tv.size(); tri++)
      if(tv[tri]->getTrackerHits().size()>nOTH){
	nOTH=tv[tri]->getTrackerHits().size();
	track = tv[tri];   // take Track from TPC (if possible)
      }
    d0        = track->getD0() ;
    fi        = track->getPhi() ;
    omega     = track->getOmega() ;
    z0        = track->getZ0() ;
    tanlambda = track->getTanLambda() ;
    helix->Initialize_Canonical(fi, d0, z0, omega, tanlambda, _bField);
    info.px = (double)helix->getMomentum()[0];
    info.py = (double)helix->getMomentum()[1];
    info.pz = (double)helix->getMomentum()[2];
  }else{            // if only cluster available
    info.px = (double)rp->getMomentum()[0];
    info.py = (double)rp->getMomentum()[1];
    info.pz = (double)rp->getMomentum()[2];
  }// if (tv.size ...
  
  float Distance[3];
  nHits=0;
  
  for(int j=0; j<nClusters; j++){  // over each cluster of particle
    Cluster *cl = cv[j];
    CalorimeterHitVec chv = cl->getCalorimeterHits();
    int nCalHits = chv.size();
    
    for(int k=0; k<nCalHits; k++){  // over each hit in cluster
      if(nHits>=N){
	std::cerr << "%% Output2a: Error in Calorimeter hits loop!" << std::endl;
	break;
      }
      
      float xpoint[3];
      for (int icomp=0;icomp<3;++icomp)
	xpoint[icomp] = (float)chv[k]->getPosition()[icomp];
      double dist, Edist;
      if(tv.size()>0){
	helix->getDistanceToPoint(xpoint, Distance);
	dist=fabs(Distance[2]);
      }else{
	float tt=0, norm2=0;
	for(int ci=0; ci<3; ci++){
	  tt += Cpos[ci]*xpoint[ci];
	  norm2 += Cpos[ci]*Cpos[ci];
	}
	tt /= norm2;
	norm2=0;
	for(int ci=0; ci<3; ci++){
	  Distance[ci]=xpoint[ci]-tt*Cpos[ci];
	  norm2 += Distance[ci]*Distance[ci];
	}
	dist=sqrt(norm2);
      }
      Edist=dist*chv[k]->getEnergy();
      info.dmean += dist;   // 3D-distance
      info.Edmean += Edist; // 3D-distance energy weighted
      
      nHits++;
    }// for k ...
  }// for j ...
  info.dmean = info.dmean/(float)N ;
  info.Edmean = info.Edmean/Esum ;
  
  if(info.Necal>0){
    info.EtoN_ecal = info.Eecal/info.Necal;
  }else{
    info.EtoN_ecal = -1.;
  }
  if(info.Nhcal>0){
    info.EtoN_hcal = info.Ehcal/info.Nhcal;
  }else{
    info.EtoN_hcal = -1.;
  }
  info.EecalToEtot = info.Eecal/(info.Eecal+info.Ehcal) ;
  
  if(!noCluster){
    if(tv.size()>0){
      pdf->VO->SetValue("Eecal", info.Eecal) ;
      pdf->VO->SetValue("Ehcal", info.Ehcal) ;
      pdf->VO->SetValue("ex", info.ex) ;
      pdf->VO->SetValue("dmean", info.dmean) ;
      pdf->VO->SetValue("Edmean", info.Edmean) ;
      pdf->VO->SetValue("EtoN_ecal", info.EtoN_ecal) ;
      pdf->VO->SetValue("EtoN_hcal", info.EtoN_hcal) ;
      pdf->VO->SetValue("Necal", info.Necal) ;
      pdf->VO->SetValue("Nhcal", info.Nhcal) ;
      pdf->VO->SetValue("L1", info.L1) ;
      pdf->VO->SetValue("L2", info.L2) ;
      pdf->VO->SetValue("L3", info.L3) ;
      pdf->VO->SetValue("EL1", info.EL1) ;
      pdf->VO->SetValue("EL2", info.EL2) ;
      pdf->VO->SetValue("EL3", info.EL3) ;
      pdf->VO->SetValue("EecalToEtot", info.EecalToEtot) ;
    }else{
      npdf->VO->SetValue("Eecal", info.Eecal) ;
      npdf->VO->SetValue("Ehcal", info.Ehcal) ;
      npdf->VO->SetValue("ex", info.ex) ;
      npdf->VO->SetValue("dmean", info.dmean) ;
      npdf->VO->SetValue("Edmean", info.Edmean) ;
      npdf->VO->SetValue("EtoN_ecal", info.EtoN_ecal) ;
      npdf->VO->SetValue("EtoN_hcal", info.EtoN_hcal) ;
      npdf->VO->SetValue("Necal", info.Necal) ;
      npdf->VO->SetValue("Nhcal", info.Nhcal) ;
      npdf->VO->SetValue("L1", info.L1) ;
      npdf->VO->SetValue("L2", info.L2) ;
      npdf->VO->SetValue("L3", info.L3) ;
      npdf->VO->SetValue("EL1", info.EL1) ;
      npdf->VO->SetValue("EL2", info.EL2) ;
      npdf->VO->SetValue("EL3", info.EL3) ;
      npdf->VO->SetValue("EecalToEtot", info.EecalToEtot) ;    
    }
  }

};
