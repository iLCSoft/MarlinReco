#include <cmath>
#include <vector>

#include <marlin/Global.h>
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"
#include "DD4hep/DD4hepUnits.h"

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCRelation.h"
#include "EVENT/MCParticle.h"
#include "IMPL/ClusterImpl.h"
#include <lcio.h>

#include "TVector3.h"

#include "CreatePDFProcessor.hh"

CreatePDFProcessor aCreatePDFProcessor ;

CreatePDFProcessor::CreatePDFProcessor()
  : Processor("CreatePDFProcessor") {
  
  // Processor description
  _description = "Processor for PDF templetes creation";

  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "OutputPDFFilename" ,
			   "Input collection of Reconstrcuted Particle",
			   _filename,
			   std::string("pdf_standard_v01.root"));

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "RecoParticleCollection" ,
			   "Input collection of Reconstrcuted Particle",
			   _PfoCollection,
			   std::string("PandoraPFOs"));
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "MCTruthLinkCollection" ,
			   "Input collection of Reconstrcuted Particle",
			   _LinkCollection,
			   std::string("RecoMCTruthLink"));

  std::vector< float > parselectron,parsmuon,parspion,parskaon,parsproton;
  parselectron.push_back(-2.40638e-03);
  parselectron.push_back(7.10337e-01);
  parselectron.push_back(2.87718e-01);
  parselectron.push_back(-1.71591e+00);
  parselectron.push_back(0.0);
  registerProcessorParameter( "dEdxParameter_electron" ,
  			      "dEdx Parameters for electron",
  			      _dEdxParamsElectron,
  			      parselectron);

  parsmuon.push_back(8.11408e-02);
  parsmuon.push_back(9.92207e-01);
  parsmuon.push_back(7.58509e+05);
  parsmuon.push_back(-1.70167e-01);
  parsmuon.push_back(4.63670e-04);
  registerProcessorParameter( "dEdxParameter_muon" ,
  			      "dEdx Parameters for muon",
  			      _dEdxParamsMuon,
  			      parsmuon);

  parspion.push_back(8.10756e-02);
  parspion.push_back(-1.45051e+06);
  parspion.push_back(-3.09843e+04);
  parspion.push_back(2.84056e-01);
  parspion.push_back(3.38131e-04);
  registerProcessorParameter( "dEdxParameter_pion" ,
  			      "dEdx Parameters for pion",
  			      _dEdxParamsPion,
  			      parspion);

  parskaon.push_back(7.96117e-02);
  parskaon.push_back(4.13335e+03);
  parskaon.push_back(1.13577e+06);
  parskaon.push_back(1.80555e-01);
  parskaon.push_back(-3.15083e-04);
  registerProcessorParameter( "dEdxParameter_kaon" ,
  			      "dEdx Parameters for Kaon",
  			      _dEdxParamsKaon,
  			      parskaon);

  parsproton.push_back(7.78772e-02);
  parsproton.push_back(4.49300e+04);
  parsproton.push_back(9.13778e+04);
  parsproton.push_back(1.50088e-01);
  parsproton.push_back(-6.64184e-04);
  registerProcessorParameter( "dEdxParameter_proton" ,
  			      "dEdx Parameters for proton",
  			      _dEdxParamsProton,
  			      parsproton);

  registerProcessorParameter( "dEdxNormalization" ,
  			      "dEdx Normalization",
  			      _dEdxNormalization,
  			      float(1.350e-7));

  registerProcessorParameter( "dEdxErrorFactor" ,
  			      "dEdx Error factor",
  			      _dEdxErrorFactor,
  			      float(7.55));
} 

void CreatePDFProcessor::init() { 
  //streamlog_out(DEBUG) << "   init called  " << std::endl ;
  //dEdx parameters                                     
  double pars[28];
  for(int i=0;i<5;i++) pars[i]=_dEdxParamsElectron[i];
  for(int i=0;i<5;i++) pars[i+5]=_dEdxParamsMuon[i];
  for(int i=0;i<5;i++) pars[i+10]=_dEdxParamsPion[i];
  for(int i=0;i<5;i++) pars[i+15]=_dEdxParamsKaon[i];
  for(int i=0;i<5;i++) pars[i+20]=_dEdxParamsProton[i];
  pars[25] =true;    //(double) _UseBayes;
  pars[26] =(double) _dEdxNormalization;
  pars[27] =(double) _dEdxErrorFactor;

  //create histograms
  std::string histtxt,labeltxt;
  for(int i=0;i<6;i++){
    for(int j=0;j<21;j++){
      switch(j){
      case 0:
        //momentum          
        histtxt="hmomentum" + itos(i+1);
        labeltxt="p(track)";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",50,0.0,250.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 1:
        //ep                                                           
        histtxt="hep" + itos(i+1);
        labeltxt="E/p";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",40,0.0,2.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 2:
        //ehad                    
        histtxt="hehad" + itos(i+1);
        labeltxt="ECAL/(ECAL+HCAL)";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",50,0.0,1.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 3:
        //photonEnergy(reserved)      
        histtxt="hphotonEnergy" + itos(i+1);
        labeltxt="E(photon) (GeV)";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",50,0.0,50.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 4:
        //mucal          
        histtxt="hmucal" + itos(i+1);
        labeltxt="MuCAL (GeV)";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",100,0.0,20.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 5:
        //chi2               
        histtxt="hchi2" + itos(i+1);
        labeltxt="#chi^{2}";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",104,-2.0,50.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 6:
        //shower max/expected shower max                    
        histtxt="hldiscrepancy" + itos(i+1);
        labeltxt="ShowerMax/EM ShowerMax";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",65,-30.0,100.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 7:
        //absorption length/expected shower max       
        histtxt="htdiscrepancy" + itos(i+1);
        labeltxt="absorption Length/Rm";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",100,0.0,2.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 8:
        //shower max using hits      
        histtxt="hxl20" + itos(i+1);
        labeltxt="xl20";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",100,0.0,100.0);
        //pidvariable[i][j]=new TH1F(histtxt.c_str(),"",80,0.0,40.0);    
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 9:
        //likelihood using electron hypothesis 
        histtxt="hlikeliele" + itos(i+1);
        labeltxt="Log(likelihood) (electron)";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",200,-100.0,0.0);
        //pidvariable[i][j]=new TH1F(histtxt.c_str(),"",80,0.0,40.0);  
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 10:
        //likelihood using muon hypothesis                       
        histtxt="hlikelimuo" + itos(i+1);
        labeltxt="Log(likelihood) (#mu)";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",200,-100.0,0.0);
        //pidvariable[i][j]=new TH1F(histtxt.c_str(),"",80,0.0,40.0);     
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 11:
        //likelihood using pion hypothesis  
        histtxt="hlikelipi" + itos(i+1);
        labeltxt="Log(likelihood) (#pi)";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",200,-100.0,0.0);
        //pidvariable[i][j]=new TH1F(histtxt.c_str(),"",80,0.0,40.0);  
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 12:
        //likelihood using kaon hypothesis 
        histtxt="hlikelik" + itos(i+1);
        labeltxt="Log(likelihood) (K)";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",200,-100.0,0.0);
        //pidvariable[i][j]=new TH1F(histtxt.c_str(),"",80,0.0,40.0);  
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 13:
        //likelihood using proton hypothesis    
        histtxt="hlikelip" + itos(i+1);
        labeltxt="Log(likelihood) (p)";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",200,-100.0,0.0);
        //pidvariable[i][j]=new TH1F(histtxt.c_str(),"",80,0.0,40.0);   
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 14:
        //likelihood ratio
        histtxt="hlikelikpi" + itos(i+1);
        labeltxt="Log(likelihood(K)/likelihood(#pi))";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",160,-50.0,30.0);
        //pidvariable[i][j]=new TH1F(histtxt.c_str(),"",80,0.0,40.0);   
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 15:
        //likelihood ratio                                                              
        histtxt="hlikelipk" + itos(i+1);
        labeltxt="Log(likelihood(p)/likelihood(K))";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",160,-50.0,30.0);
        //pidvariable[i][j]=new TH1F(histtxt.c_str(),"",80,0.0,40.0);    
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 16:
        //delta x                                                                
        histtxt="hdeltax" + itos(i+1);
        labeltxt="#Deltax#timesQ";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",200,-5.0,5.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 17:
        //delta z              
        histtxt="hdeltaz" + itos(i+1);
        labeltxt="#Deltaz";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",200,-5.0,5.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 18:
        //cluster depth            
        histtxt="hcludepth" + itos(i+1);
        labeltxt="cluster depth";
        pidvariable[i][j]=new TH1F(histtxt.c_str(),"",200,1500.0,3500.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 19:
        //hit mean             
        histtxt="hhitmean" + itos(i+1);
        labeltxt="hit mean";
        if(i==0) pidvariable[i][j]=new TH1F(histtxt.c_str(),"",100,0.0,2500.0);
        if(i!=0) pidvariable[i][j]=new TH1F(histtxt.c_str(),"",100,0.0,1000.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      case 20:
        //hitrms              
        histtxt="hhitrms" + itos(i+1);
        labeltxt="hitrms";
        if(i==0) pidvariable[i][j]=new TH1F(histtxt.c_str(),"",100,0.0,2500.0);
        if(i!=0) pidvariable[i][j]=new TH1F(histtxt.c_str(),"",100,0.0,1000.0);
        pidvariable[i][j]->SetMarkerStyle(20);
        pidvariable[i][j]->SetFillColor(0);
        pidvariable[i][j]->SetLineColor(1);
        pidvariable[i][j]->SetLineWidth(2.0);
        pidvariable[i][j]->GetYaxis()->SetTitle("Normalized Events");
        pidvariable[i][j]->GetXaxis()->SetTitle(labeltxt.c_str());
        pidvariable[i][j]->SetMinimum(0.0);
        break;
      }
      if(i==0) pidvariable[i][j]->SetLineColor(5);
      if(i==1) pidvariable[i][j]->SetLineColor(6);
      if(i==2) pidvariable[i][j]->SetLineColor(2);
      if(i==3) pidvariable[i][j]->SetLineColor(3);
      if(i==4) pidvariable[i][j]->SetLineColor(4);
      if(i==5) pidvariable[i][j]->SetLineColor(7);
    }
  }

  _myPID = new LikelihoodPID(pars);

  //get B field                                          
  double bfield[3]={};
  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
  lcdd.field().magneticField({0.0,0.0,0.0}, bfield);
  _bfield = (float)bfield[2]/dd4hep::tesla;

  printParameters();  
}

void CreatePDFProcessor::processRunHeader( LCRunHeader* /*run*/) { 
} 

void CreatePDFProcessor::processEvent( LCEvent * evt ) { 
  _PFOCol = evt->getCollection( _PfoCollection ) ;
  _LinkCol = evt->getCollection( _LinkCollection ) ;
  int npfo = _PFOCol->getNumberOfElements();
  int nln = _LinkCol->getNumberOfElements();

  std::map<int, int> pid;
  pid.insert(std::map<int, int>::value_type(11, 0));
  pid.insert(std::map<int, int>::value_type(13, 1));
  pid.insert(std::map<int, int>::value_type(211, 2));
  pid.insert(std::map<int, int>::value_type(321, 3));
  pid.insert(std::map<int, int>::value_type(2212, 4));

  for(int i=0;i<npfo;i++){
    ReconstructedParticle* part = (ReconstructedParticle*) _PFOCol->getElementAt(i);
    if(part->getCharge()==0) continue;  //neutral
    
    EVENT::Track* trk = part->getTracks()[0];
    EVENT::ClusterVec cluvec = part->getClusters();
    
    TVector3 pp = part->getMomentum();
    
    //get deposit energy
    double ecal=0.0,hcal=0.0,mucal=0.0;
    if(cluvec.size()!=0){
      for(unsigned int j=0;j<cluvec.size();j++){
	ecal+=cluvec[j]->getSubdetectorEnergies()[0];
	hcal+=cluvec[j]->getSubdetectorEnergies()[1];
	mucal+=cluvec[j]->getSubdetectorEnergies()[2];
      }
    }
    
    //check MC particle
    int pdfid=-1;
    for(int j=0;j< nln;j++){
      LCRelation* plc = (LCRelation*) _LinkCol->getElementAt( j );
      ReconstructedParticle* fp = (ReconstructedParticle*) plc->getFrom();
      
      if(part == fp){
	MCParticle* partmc = (MCParticle*) plc->getTo();
	int mcpdg = fabs(partmc->getPDG());
	if(mcpdg!=11 && mcpdg!=13 && mcpdg!=211 && mcpdg!=321 && mcpdg!=2212) continue;
	pdfid = pid[mcpdg];
	if(pdfid==1 && mucal==0.0) pdfid=5;
	break;
      }
    }  //done.
    
    if(pdfid<0) continue;  //do not match to MC particle
    
    //start to fill historams
    //momentum
    pidvariable[pdfid][0]->Fill(pp.Mag());
 
    if(cluvec.size()!=0){
      pidvariable[pdfid][1]->Fill((ecal+hcal)/pp.Mag());
      pidvariable[pdfid][2]->Fill(ecal/(ecal+hcal));
      pidvariable[pdfid][4]->Fill(mucal);
    
      //get shower profiles
      FloatVec shapes=cluvec[0]->getShape();      
      pidvariable[pdfid][5]->Fill(shapes[0+4]);
      pidvariable[pdfid][6]->Fill(shapes[5+4]);
      pidvariable[pdfid][7]->Fill(fabs(shapes[3+4])/shapes[6+4]);
      pidvariable[pdfid][8]->Fill(fabs(shapes[15+4])/7.00);

      //for mu/pi<5.0GeV?
      if(pp.Mag()>0.3 && pp.Mag()<6.0){
	pidvariable[pdfid][18]->Fill(shapes[17+4]);
	pidvariable[pdfid][19]->Fill(shapes[18+4]);
	pidvariable[pdfid][20]->Fill(shapes[19+4]);
      }
    }

    //get dedx info
    double dedx = trk->getdEdx();
    double dedxerr = trk->getdEdxError();
    double ehypo = -0.5*fabs(_myPID->get_dEdxChi2(0, pp, dedxerr, dedx))+_myPID->get_dEdxFactor(0,pp, dedxerr, dedx);
    double muhypo = -0.5*fabs(_myPID->get_dEdxChi2(1, pp, dedxerr, dedx))+_myPID->get_dEdxFactor(1, pp, dedxerr, dedx);
    double pihypo = -0.5*fabs(_myPID->get_dEdxChi2(2, pp, dedxerr, dedx))+_myPID->get_dEdxFactor(2, pp, dedxerr, dedx);
    double khypo = -0.5*fabs(_myPID->get_dEdxChi2(3, pp, dedxerr, dedx))+_myPID->get_dEdxFactor(3, pp, dedxerr, dedx);
    double phypo = -0.5*fabs(_myPID->get_dEdxChi2(4, pp, dedxerr, dedx))+_myPID->get_dEdxFactor(4, pp, dedxerr, dedx);
    
    pidvariable[pdfid][9]->Fill(ehypo);
    pidvariable[pdfid][10]->Fill(muhypo);
    pidvariable[pdfid][11]->Fill(pihypo);
    pidvariable[pdfid][12]->Fill(khypo);
    pidvariable[pdfid][13]->Fill(phypo);
    
    pidvariable[pdfid][14]->Fill(std::log(exp(khypo)/exp(pihypo)));
    pidvariable[pdfid][15]->Fill(std::log(exp(phypo)/exp(khypo)));
    
    //deltax and deltaz
    float delpos[3]={0.0,0.0,0.0};
    if(cluvec.size()!=0){   // && pp.Mag()<6.0){
      CalculateDeltaPosition(part->getCharge(), pp, cluvec[0]->getPosition(), delpos);
      pidvariable[pdfid][16]->Fill(delpos[0]);
      pidvariable[pdfid][17]->Fill(delpos[2]);
    }
  }  
}

void CreatePDFProcessor::check( LCEvent * /*evt*/ ) { 
}

void CreatePDFProcessor::end() { 
  //delete _myShowerShapes;
  //open tfile
  _fpdf = new TFile(_filename.c_str(), "RECREATE");

  for(int i=0;i<6;i++){
    for(int j=0;j<21;j++){
      pidvariable[i][j]->Write();
    }
  }

  _fpdf->Close();

  delete _fpdf;
  delete _myPID;
}

void CreatePDFProcessor::CalculateDeltaPosition(float charge, TVector3& p, const float* calpos, float* delpos){
  //calculate some variables for preparation
  //transverse momentum
  double prphi=sqrt(pow(p[0],2)+pow(p[1],2));
  //momentum
  double pp=p.Mag();
  //radius in r-phi plane
  double radius=prphi/(0.3*_bfield);
  //tangent lambda
  double tanlam=p[2]/prphi;
  //change to unit vector
  p[0]=p[0]/prphi;
  p[1]=p[1]/prphi;
  p[2]=p[2]/pp;

  //chenge position vector to meter
  float calcalpos[3];
  calcalpos[0]=calpos[0]/1000.0;
  calcalpos[1]=calpos[1]/1000.0;
  calcalpos[2]=calpos[2]/1000.0;


  //radius at the tracker end
  double rradius=sqrt(pow(calcalpos[0],2)+pow(calcalpos[1],2));

  //cout << "check val. " << prphi << " " << radius << " " << tanlam
  //     << " " << rradius << endl;

  //cal. the position of the center of track circle
  TVector3 cc;
  cc[0]=-charge*p[1]*prphi/(0.3*_bfield);
  cc[1]=charge*p[0]*prphi/(0.3*_bfield);
  cc[2]=radius*tanlam;

  //cal. sign and cosine 2theta
  double sintheta=charge*rradius/(2*radius);
  double costheta=sqrt(1-sintheta*sintheta);
  double sin2theta=2*sintheta*costheta;
  double cos2theta=costheta*costheta-sintheta*sintheta;

  TVector3 calpos2;
  calpos2[2]=rradius*tanlam;
  calpos2[0]=cc[0]*(1-cos2theta)+cc[1]*sin2theta;
  calpos2[1]=-cc[0]*sin2theta+cc[1]*(1-cos2theta);

  //cal. difference of the position
  //float delpos[3];
  delpos[0]=charge*fabs(calcalpos[0]-calpos2[0]);
  delpos[1]=charge*fabs(calcalpos[1]-calpos2[1]);
  delpos[2]=fabs(calcalpos[2]-calpos2[2]);


  //cout << "calpos: " << calcalpos[0] << " " << calcalpos[1] << " " << calcalpos[2] << endl; 
  //cout << "calpos2: " << calpos2[0] << " " << calpos2[1] << " " << calpos2[2] << endl; 
  //cout << "result: " << delpos[0] << " " << delpos[1] << " " << delpos[2] << endl; 
  return;
}
 
