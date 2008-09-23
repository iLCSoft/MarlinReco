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
#include <EVENT/LCFloatVec.h>

#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/LCCollectionVec.h>

#include <UTIL/PIDHandler.h>
 
#include "PDF.h"

using namespace lcio ;
using namespace marlin ;

PFOID aPFOID ;

PFOID::PFOID() : Processor("PFOID") {
  _description = "Performs particle identification" ;

  registerProcessorParameter("RecoParticleCollection",
			     "Name of the PFO collection",
			     _recoCol,
			     std::string("RecoParticles")) ;

  registerProcessorParameter("FilePDFName",
			     "Name of file containing pdfs",
			     _filename,std::string("pdf_n.txt")) ;

  registerProcessorParameter("neutralFilePDFName",
			     "Name of file containing pdfs for neutral particle",
			     _filename1,
			     std::string("npdf_n.txt")) ;

  registerProcessorParameter("Debug",
			     "Debugging?",
			     _debug,
			     int(0));

}


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

}


void PFOID::processRunHeader( LCRunHeader * run ) {
  _nRun++ ;
}


void PFOID::processEvent( LCEvent * evt ) {

  if (_debug>0)
    std::cout << "%% PFOID: Event " << _nEvt << std::endl;

  const gear::BField& gearBField = Global::GEAR->getBField();
  _bField = float(gearBField.at(gear::Vector3D(0.0,0.0,0.0)).z()); //BField should be constant and point in z direction

  try{


    LCCollection *col = evt->getCollection(_recoCol.c_str());
    PIDHandler partID ( col );

    StringVec pNames;
    pNames.push_back("Eccentricity");
    pNames.push_back("DMean");
    pNames.push_back("EtoN_ecal");
    pNames.push_back("EtoN_hcal");
    pNames.push_back("EecalToEtotal");
    pNames.push_back("L1");
    pNames.push_back("L2");
    pNames.push_back("L3");
    pNames.push_back("EL1");
    pNames.push_back("EL2");
    pNames.push_back("EL3");

    int algo = partID.addAlgorithm("PFOID",pNames);
    int nRecos = col->getNumberOfElements() ;
    for(int i=0; i<nRecos; i++){  // over all reco. particles
      ReconstructedParticle *recopart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));


      fill_info(i,recopart); // fills also pdf->VO or npdf->VO
      ClusterVec clv = recopart->getClusters();

      LCFloatVec params;
      params.clear();

      if(clv.size()>0){
	if(info.withTrack){

	  if (_debug>0) {
	    std::cout << " RecoParticle " << i << std::endl;
	    std::cout << "  likelihood(electron) : ";
	    std::cout << pdf->GetLikelihood("electron") << std::endl;
	    std::cout << "  likelihood(muon)     : ";
	    std::cout << pdf->GetLikelihood("muon") << std::endl;
	    std::cout << "  likelihood(pion)     : ";
	    std::cout << pdf->GetLikelihood("pion") << std::endl;
	  }
	  
	  if (_debug>1) {
	    std::cout << "  Ecal/Etot : " << info.EecalToEtot;
	    std::cout << "  EtoN_ecal : " << info.EtoN_ecal;
	    std::cout << "  EtoN_hcal : " << info.EtoN_hcal;
	    std::cout << "  dmean : " << info.dmean << std::endl;
	    std::cout << "  ex : " << info.ex;
	    std::cout << "  L1 : " << info.L1;
	    std::cout << "  L2 : " << info.L2;
	    std::cout << "  L3 : " << info.L3 << std::endl;
	  }	  
	  
	  double maxv=0.;
	  double lhde, lhdm, lhdp;

	  int pdg = 0;
	  int type = 0;
	  const int pdge = -11;
	  const int pdgm = -13;
	  const int pdgp = 211;
	  

	  lhde = pdf->GetLikelihood("electron");
	  lhdm = pdf->GetLikelihood("muon");
	  lhdp = pdf->GetLikelihood("pion");

	  pdg  = pdge;
	  maxv = lhde;
	  type = 0;

	  if (lhdm>maxv) {
	    maxv = lhdm;
	    pdg = pdgm;
	    type = 1;
	  }

	  if (lhdp>maxv) {
	    maxv = lhdp;
	    pdg = pdgp;
	    type = 2;
	  }

	  if (recopart->getCharge()<0)
	    pdg = -pdg;

	  params.push_back( float(info.ex) ) ;
	  params.push_back( float(info.dmean) ) ;
	  params.push_back( float(info.EtoN_ecal) );
	  params.push_back( float(info.EtoN_hcal) );
	  params.push_back( float(info.EecalToEtot) );
	  params.push_back( float(info.L1) );
	  params.push_back( float(info.L2) );
	  params.push_back( float(info.L3) );
	  params.push_back( float(info.EL1) );
	  params.push_back( float(info.EL2) );
	  params.push_back( float(info.EL3) );	  
	  
	  partID.setParticleID(recopart, type, pdg, float(maxv), algo, params );


	}else{

	  if (_debug>0) {
	    std::cout << " RecoParticle " << i << std::endl;
	    std::cout << "  likelihood(photon) : " ;
	    std::cout << npdf->GetLikelihood("photon") << std::endl;
	    std::cout << "  likelihood(kaon/neutron)     : " ;
	    std::cout << npdf->GetLikelihood("kaon/neutron") << std::endl;
	  }
	  
	  if (_debug>1) {
	    std::cout << "  Ecal/Etot : " << info.EecalToEtot;
	    std::cout << "  EtoN_ecal : " << info.EtoN_ecal;
	    std::cout << "  EtoN_hcal : " << info.EtoN_hcal;
	    std::cout << "  dmean : " << info.dmean << std::endl;
	    std::cout << "  ex : " << info.ex;
	    std::cout << "  L1 : " << info.L1;
	    std::cout << "  L2 : " << info.L2;
	    std::cout << "  L3 : " << info.L3 << std::endl;
	  }	  

	  int pdg,type;

	  const int pdg_gamma = 22;
	  const int pdg_k0l   = 130;
	  
	  double maxv;
	  double lhd_gamma;
	  double lhd_k0l;

	  lhd_gamma = npdf->GetLikelihood("photon");
	  lhd_k0l   = npdf->GetLikelihood("kaon/neutron");

	  pdg = pdg_gamma;
	  type = 3;
	  maxv = lhd_gamma;

	  if (lhd_k0l>maxv) {
	    maxv = lhd_k0l;
	    pdg = pdg_k0l;
	    type = 4;
	  }

	  params.push_back( float(info.ex) ) ;
	  params.push_back( float(info.dmean) ) ;
	  params.push_back( float(info.EtoN_ecal) );
	  params.push_back( float(info.EtoN_hcal) );
	  params.push_back( float(info.EecalToEtot) );
	  params.push_back( float(info.L1) );
	  params.push_back( float(info.L2) );
	  params.push_back( float(info.L3) );
	  params.push_back( float(info.EL1) );
	  params.push_back( float(info.EL2) );
	  params.push_back( float(info.EL3) );	  

	  partID.setParticleID(recopart, type, pdg, float(maxv), algo, params );


	}// if ( withTrack ...
      
      }else{ // unknown particle (no cluster)

	params.push_back( 0. );
	params.push_back( 0. );
	params.push_back( 0. );
	params.push_back( 0. );
	params.push_back( 0. );
	params.push_back( 0. );
	params.push_back( 0. );
	params.push_back( 0. );
	params.push_back( 0. );
	params.push_back( 0. );
	params.push_back( 0. );

	int pdg  = 999999;
	float lhd = 0.0;
	int type = 5;

	partID.setParticleID(recopart, type, pdg, lhd, algo, params );

      }

    }// for i ...
  }catch(DataNotAvailableException &e){ }


  if (_debug>0) {
    std::cout << "%% PFOID: End of PFOID\n" << std::endl; 
    std::cout << std::endl;
  }
  _nEvt++ ;

}




void PFOID::check( LCEvent * evt ) {}




void PFOID::end() {

  delete pdf;
  delete npdf;
}




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
  if (_debug > 0)
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

}
