#include "CreatePDFs.h"
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

#include <EVENT/LCCollection.h>
#include <EVENT/LCObject.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/Cluster.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/MCParticle.h>

using namespace lcio ;
using namespace marlin ;

CreatePDFs aCreatePDFs ;


//****************************************************//
CreatePDFs::CreatePDFs() : Processor("CreatePDFs") {
  _description = "Creates from a gauge sample pdf files for a likelihood based particle identification";

  registerProcessorParameter("MCParticleCollection","Name of the MC Particle collection",_mcCol,std::string("MCParticle")) ;

  registerProcessorParameter("RecoParticleCollection","Name of a reconstructed Particle collection containing PFOs like Wolf, PandoraPFA, ...",_recoCol,std::string("RecoParticles")) ;

  registerProcessorParameter("ChargedPDFfilename","Name of output file for charged particles",_filename_c,std::string("pdf.txt")) ;

  registerProcessorParameter("NeutralPDFfilename","Name of output file for neutral particles",_filename_n,std::string("npdf.txt")) ;

  std::vector<int> chStart;
  registerProcessorParameter( "ChargedRunStart" , "Run numbers where a new charged category sample starts " , _chStart, chStart);

  std::vector<int> nStart;
  registerProcessorParameter( "NeutralRunStart" , "Run numbers where a new neutral category sample starts " , _nStart, nStart);

  std::vector<int> pidCol;
  registerProcessorParameter("PIDCollection","PIDs of the particles, that should be here", _pidCol, pidCol);
  
  std::vector<int> NoOfPID;
  registerProcessorParameter("PIDnumbers","number of pids per category",_NoOfPID,NoOfPID);

  registerProcessorParameter("RelCollection","Name of relation collection",_relCol,std::string("RelationCaloHit")) ;

  registerProcessorParameter("MCTracksRelCollection","Name of MC to Tracks Relation collection",_MCTracksRelCol,std::string("LDCTracksMCP"))
;}


//****************************************************//
void CreatePDFs::init() {
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  _nParticle = 0 ;

  std::cout << " ************************** " << std::endl ;
  std::cout << "        CreatePDFs          " << std::endl ;
  std::cout << " ************************** " << std::endl ;

  std::cout << "\n Please check that Configure.h is correctly set!" << std::endl ;
  
  if(_chStart.size()==0){
    std::cout << " [CreatePDFs : WARN] - No charged PDF file is created;";
    std::cout << " No category event start numbers are specified!\n";
  }

  if(_nStart.size()==0){
    std::cout << " [CreatePDFs : WARN] - No neutral PDF file is created;";
    std::cout << " No category event start numbers are specified!\n";
  }

  for(unsigned int i=1; i<_chStart.size(); i++){
    if(_chStart[i-1]>=_chStart[i]){
      std::cout << " [CreatePDFs : ERR] - Category event start numbers should be in increasing order! (charged)\n";
      exit(-2);
    }
  }
  
  for(unsigned int i=1; i<_nStart.size(); i++){
    if(_nStart[i-1]>=_nStart[i]){
      std::cout << " [CreatePDFs : ERR] - Category event start numbers should be in increasing order! (neutral)\n";
      exit(-2);
    }
  }
  
  if(_chStart.size()>0 && _nStart.size()>0){
    if(_chStart[_chStart.size()-1]>=_nStart[0]){
      std::cout << " [CreatePDFs : ERR] - Category of charged particle samples must start before Category of neutral particle samples!\n";
      exit(-3);
    }
    if(_NoOfPID.size()!=_chStart.size()+_nStart.size()){
      std::cout << " [CreatePDFs : ERR] - Total number of categories should councide with array size of PID Numbers!\n";
      exit(-4);
    }
    double sum=0;
    for(unsigned int i=0; i<_NoOfPID.size(); i++)
      sum+=_NoOfPID[i];
    if(sum!=_pidCol.size()){
      std::cout << " [CreatePDFs : ERR] - Size of PID array should coincide with sum of PID Numbers!\n";
      exit(-5);
    }
  }

  if(_chStart.size()>0){
    std::string catnames[CNOCATS]=CCATS;
    std::string varnames[NOVARS]=VARIABLES;
    std::string cvarnames[CNOVARS]=CVARIABLES;
    
    double ranges[3*CNOVARS]=CRANGES;

    curr_cat=catnames[0];

    pdf = new PDF(CPDFNAME,CNOCATS,catnames,CNOHISTS,NOVARS,varnames);
    
    for(int j=0; j<CNOHISTS; j++){
      int bins=(int)ranges[3*j];
      double start=ranges[3*j+2];
      double width=(ranges[3*j+1]-start)/ranges[3*j];
      pdf->InitHistogram(cvarnames[j],CHISTDIM,&cvarnames[j],&start,&width,&bins);
    }
  }else{
    pdf=NULL;
  }

  if(_nStart.size()>0){
    std::string catnames[NNOCATS]=NCATS;
    std::string varnames[NOVARS]=VARIABLES;
    std::string nvarnames[NNOVARS]=NVARIABLES;

    if(_chStart.size()==0)
      curr_cat=catnames[0];
    
    double ranges[3*NNOVARS]=NRANGES;

    npdf = new PDF(NPDFNAME,NNOCATS,catnames,NNOHISTS,NOVARS,varnames);
    
    for(int j=0; j<NNOHISTS; j++){
      int bins=(int)ranges[3*j];
      double start=ranges[3*j+2];
      double width=(ranges[3*j+1]-start)/ranges[3*j];
      npdf->InitHistogram(nvarnames[j],NHISTDIM,&nvarnames[j],&start,&width,&bins);
    }
  }else{
    npdf=NULL;
  }


  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  _bField = float(gearTPC.getDoubleVal("BField"));
  

  myindex=0;
  MyPidCol.clear();
  for(int j=0; j<_NoOfPID[0]; j++)
    MyPidCol.push_back(_pidCol[j]);

  noClusterParticle=0;  
}



//***********************************************************//
void CreatePDFs::processRunHeader( LCRunHeader * run ) {
  int oldIndex=myindex;

  if(_nStart.size()>0){
    if(_chStart.size()>0 && _nRun<_nStart[0])
      for(unsigned int i=0; i<_chStart.size(); i++)
	if(_nRun>=_chStart[i]) myindex=i;
    for(unsigned int i=0; i<_nStart.size(); i++)
      if(_nRun>=_nStart[i]) myindex=_chStart.size()+i;
  }else{
    for(unsigned int i=0; i<_chStart.size(); i++)
      if(_nRun>=_chStart[i]) myindex=i;
  }

  if(myindex!=oldIndex){
    std::cout << "%% CreatePDFs: " << _nParticle << " written to category " << curr_cat << std::endl;
    _nParticle=0;
    std::cerr << " Index changed " << myindex << std::endl;
    if(((unsigned int)myindex)<_chStart.size()){
      std::cerr << " charged " << std::endl;
      std::string catnames[CNOCATS]=CCATS;
      curr_cat=catnames[myindex];
    }else{
      std::cerr << " neutral " << std::endl;
      std::string catnames[NNOCATS]=NCATS;
      curr_cat=catnames[myindex-_chStart.size()];
    }
    std::cerr << " cat " << curr_cat << std::endl;
    MyPidCol.clear();
    int Istart=0;
    for(int i=0; i<myindex; i++)
      Istart+=_NoOfPID[i];
    std::cerr << " Istart " << Istart << std::endl;
    for(int j=0; j<_NoOfPID[myindex]; j++)
      MyPidCol.push_back(_pidCol[Istart+j]);
    std::cerr << "%% CreatePDFs: proceed with PDGs: ";
    for(unsigned int i=0; i<MyPidCol.size(); i++)
      std::cerr << MyPidCol[i] << "  ";
    std::cerr << std::endl;
  }

  std::cerr << "%% CreatePDFs: proceed with PDGs: ";
  for(unsigned int i=0; i<MyPidCol.size(); i++)
    std::cerr << MyPidCol[i] << "  ";
  std::cerr << std::endl;

  _nRun++ ;

}


//***********************************************************//
void CreatePDFs::processEvent( LCEvent * evt ) {
  
  std::cerr << "%% CreatePDFs: Event " << _nEvt << "  category : " << curr_cat << std::endl;

  try{
 
    // relation between hits and MC particles
    const LCCollection *relcol = evt->getCollection(_relCol.c_str());
    LCRelationNavigator navigate(relcol);

    const LCCollection *TrackToMCP = evt->getCollection(_MCTracksRelCol.c_str());
    LCRelationNavigator TrackToMCPNav(TrackToMCP);

    const LCCollection *mccol = evt->getCollection(_mcCol.c_str());
    int nMCs = mccol->getNumberOfElements() ;

    // collection of reconstructed particles
    LCCollection *col = evt->getCollection(_recoCol.c_str());
    int nRecos = col->getNumberOfElements() ;

    std::cout << " # reconstructed particles: " << nRecos << std::endl;

    for(int i=0; i<nRecos; i++){  // over all reco. particle
      bool noCluster=false;
      init_info();
      HelixClass * helix = new HelixClass();

      ReconstructedParticle *rp = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
      TrackVec tv = rp->getTracks();
      ClusterVec cv = rp->getClusters();
      std::cout << " RP: " << i << " #Tracks: " << tv.size() << " #Cluster: " << cv.size() << std::endl;

      // Cluster Layers and energies
      int nClusters = cv.size(), myPDG=0;
      SimCalorimeterHit *sd=NULL;
      bool EcalHit=false, charge =false;
      float *ca, *cx, *cy, *cz;
      ClusterShapes *clsp;

      if(nClusters>0){
	for(int j=0; j<nClusters; j++){   // over each cluster of particle
	  Cluster *cl = cv[j];
	  info.Eecal += (double)cl->getSubdetectorEnergies()[0];
	  info.Ehcal += (double)cl->getSubdetectorEnergies()[1];
	  CalorimeterHitVec chv = cl->getCalorimeterHits();
	  int nCalHits = chv.size();
	  
	  for(int k=0; k<nCalHits; k++){  // over each hit in cluster
	    LCObjectVec ov = navigate.getRelatedToObjects(chv[k]);
	    
	    int cellid = chv[k]->getCellID0();
	    int layer = cellid >> 24;     // layer number
	    if(chv[k]->getType()==1){     // HCAL layer
	      info.Nhcal +=1 ;
	      if(layer>=info.L1){
		info.L3=info.L2 ;
		info.L2=info.L1 ;
		info.L1=(double)layer ;
		if(!EcalHit) sd = dynamic_cast<SimCalorimeterHit*>(ov[0]);
	      }else{
		if(layer>=info.L2){
		  info.L3=info.L2 ;
		  info.L2=(double)layer;
		}else{
		  if(layer>info.L3)
		    info.L3=(double)layer;
		}
	      }
	    }else{                        // ECAL layer
	      info.Necal += 1;
	      EcalHit=true;
	      if(layer<=info.EL1){
		info.EL3=info.EL2 ;
		info.EL2=info.EL1 ;
		info.EL1=(double)layer ;
		sd = dynamic_cast<SimCalorimeterHit*>(ov[0]);
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
	
	
	// MC info via RelationCaloHit
	charge=false;
	myPDG=0;
	if(sd!=NULL){
	  unsigned int ncontrib = sd->getNMCContributions();
	  float maxmix=0.0;
	  unsigned int ttmax=0;
	  for(unsigned int tt=0; tt<ncontrib; tt++){
	    if(sd->getEnergyCont(tt)>maxmix){
	      maxmix=sd->getEnergyCont(tt);
	      ttmax=tt;
	    }
	  }
	  if(sd->getParticleCont(ttmax)!=0){
	    std::cout << " Cluster MC : could be " << sd->getParticleCont(ttmax)->getPDG() << " (# parents: " << sd->getParticleCont(ttmax)->getParents().size() << ")  E: " << sd->getParticleCont(ttmax)->getEnergy() << " GeV ";
	    std::cout << " ( Cluster energy: " << info.Eecal+info.Ehcal << " GeV )";
	    myPDG=sd->getParticleCont(ttmax)->getPDG();
	    charge=(fabs(sd->getParticleCont(ttmax)->getCharge())>=0.0001);
	  }else{
	    std::cout << " no cluster MC info";
	  }
	}
	std::cout << std::endl;
	
	
	// ClusterShape and postion
	int N=(int)(info.Necal+info.Nhcal);
	ca = new float[N]; // for cluster shape
	cx = new float[N]; //      - " -
	cy = new float[N]; //      - " -
	cz = new float[N]; //      - " -
	int nHits=0;  // counts all calorimeter hits (N)
	float Cpos[3] = {0., 0., 0.}, Esum=0; // total cluster position & energy
	
	for(int j=0; j<nClusters; j++){  // over each cluster of particle
	  Cluster *cl = cv[j];
	  CalorimeterHitVec chv = cl->getCalorimeterHits();
	  int nCalHits = chv.size();
	  
	  for(int kk=0; kk<3; kk++)
	    Cpos[kk] += cl->getPosition()[kk]*cl->getEnergy() ;

	  Esum += cl->getEnergy();
	  
	  for(int k=0; k<nCalHits; k++){  // over each hit in cluster
	    if(nHits>=N){
	      std::cerr << "%% CreatePDFs: Error in Calorimeter hits loop!" << std::endl;
	      break;
	    }
	    ca[nHits]=chv[k]->getEnergy();
	    cx[nHits]=chv[k]->getPosition()[0];
	    cy[nHits]=chv[k]->getPosition()[1];
	    cz[nHits]=chv[k]->getPosition()[2];
	    nHits++;
	  }// for k ...
	}// for j ...
	for(int kk=0; kk<3; kk++)	Cpos[kk] /= Esum;
	
	clsp = new ClusterShapes(N,ca,cx,cy,cz);
	
	float width = clsp->getWidth();
	float evi[3];
	for(int kk=0; kk<3; kk++)	evi[kk] = clsp->getEigenValInertia()[kk];
	
	float max = evi[0];
	if(evi[1]>max) max = evi[1];
	if(evi[2]>max) max = evi[2];
	info.ex = (double)width/(double)max; // excentricity
	
	// Tracks
	float d0, fi, omega, z0, tanlambda;
	if(tv.size()>0){  // if tracks available
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
	      std::cerr << "%% CreatePDFs: Error in Calorimeter hits loop!" << std::endl;
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
	      float tt=0, norm2=0.;
	      for(int ci=0; ci<3; ci++){
		tt += Cpos[ci]*xpoint[ci];   // scalar product
		norm2 += Cpos[ci]*Cpos[ci];  // norm
	      }
	      tt /= norm2;       // parameter on straight line with min-dist
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
      }else{
	noCluster=true;
	std::cout << " WARNING!! This particle has no cluster" << std::endl;
	noClusterParticle++;
      } // if(nClusters>0)

      // MC Track info
      if(tv.size()>0){
        for(int j=0; j<nMCs; j++){
          MCParticle *mcp=dynamic_cast<MCParticle*>(mccol->getElementAt(j));
          LCObjectVec ovt = TrackToMCPNav.getRelatedFromObjects(mcp);
          Track *track=NULL;
          FloatVec weight = TrackToMCPNav.getRelatedFromWeights(mcp);
          float weightmax=-1.0;
          for(unsigned int i=0; i<ovt.size(); i++){
            if(weight[i]>weightmax){
              weightmax=weight[i];
              track = dynamic_cast<Track*>(ovt[i]);
            }// if ...
          }// for( int i=0 ...
          if(track!=NULL){
            for(unsigned int k=0; k<tv.size(); k++){
              if(track==tv[k]){
		for(unsigned int kk=0; kk<MyPidCol.size(); kk++)
		  if(mcp->getPDG()==MyPidCol[kk]){
		    myPDG=mcp->getPDG();
		    charge=(fabs(mcp->getCharge())>=0.0001);
		    std::cout << " Track MC : could be " << mcp->getPDG() << " (# parents: " << mcp->getParents().size() << ")  E: " << mcp->getEnergy() << " GeV ";
		    std::cout << " ( Track energy: " << rp->getEnergy() << " GeV )" << std::endl;
		  }
	      }
            }
          }
        }
      }
      
      bool write=false;
      for(unsigned int kk=0; kk<MyPidCol.size(); kk++){
	std::cout << " PID: " << MyPidCol[kk] ;
	if(myPDG==MyPidCol[kk]){ 
	  std::cout << " x";
	  if(charge && tv.size()>0) write=true; 
	  if(!charge && tv.size()==0) write=true;
	}
	std::cout << std::endl;
      }

      if(write){
	if(charge){
	  if(!noCluster){
	    std::cout << " Start to fill histogram" << std::endl;
	  
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

	    pdf->FillHistograms(curr_cat);
	    std::cout << "%% CreatePDFs: Written to histogram (charged) ... " << std::endl;
	    _nParticle++;
	  }
	}else{ 
	  if(!noCluster){
	    std::cout << " Start to fill histogram" << std::endl;
	    
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
	  
	    npdf->FillHistograms(curr_cat);
	    std::cout << "%% CreatePDFs: Written to histogram (neutral) ... " << std::endl;
	    _nParticle++;
	  }
	}
      }
      
      if(!noCluster){
        delete helix;
        delete ca;
        delete cx;
        delete cy;
        delete cz;
        delete clsp;
      }
    }// for i ... Reco Particles
  }catch(DataNotAvailableException &e){ }

  std::cout << "%% CreatePDFs: End of CreatePDFs\n" << std::endl; 
  _nEvt++ ;

}

void CreatePDFs::check( LCEvent * evt ) {
}

void CreatePDFs::end() {
  std::cout << "%% CreatePDFs: " << _nParticle << " written to category " << curr_cat << std::endl;
  std::cout << "%% CreatePDFs: # particles without cluster: " << noClusterParticle << std::endl;
  if(pdf!=NULL){
    pdf->WritePDF(_filename_c);
    delete pdf;
  }
  if(npdf!=NULL){
    npdf->WritePDF(_filename_n);
    delete npdf;
  }
}


void CreatePDFs::init_info(){
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
  info.Necal=0;
  info.Nhcal=0;
  info.L1=-1;
  info.L2=-1;
  info.L3=-1;
  info.EL1=100;
  info.EL2=100;
  info.EL3=100;
  info.EecalToEtot=0.;
}
