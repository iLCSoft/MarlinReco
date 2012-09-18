
#include "BCalReco.h"
#include "BCalReconstruction.h"

#include <UTIL/CellIDDecoder.h>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogram3D.h>
#include <AIDA/IAxis.h>
#include <AIDA/ICloud2D.h>
#endif

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <TRandom3.h>

#define epsR 0.1
#define epsPhi 2e-05



using namespace lcio;
using namespace marlin;
using namespace std;


BCalReco aBCalReco;



BCalReco::BCalReco() : Processor("BCalReco"){
  _description = "Filling Histograms with deposited energy from beamstrahlung pairs in BeamCal";

  registerInputCollection( LCIO::MCPARTICLE,
			   "CollectionName" ,
			   "Name of the MCParticle collection" ,
			   _colName ,
			   std::string("MCParticle") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "RecoPartCollectionName" ,
			    "Name of the Beamcal reconstructed particles collection" ,
			    _BCALcolName ,
			    std::string("BCALParticles") );

  registerOutputCollection( LCIO::CLUSTER,
			    "BCALClusterName" ,
			    "Name of the Beamcal clusters collection" ,
			    _BCALClustersName ,
			    std::string("BCALClusters") );

  registerOutputCollection( LCIO::LCRELATION,
			    "BCALMCTruthLinkName" ,
			    "Name of the Beamcal cluster to mc-particle relation collection" ,
			    _BCALMCTruthLinkName ,
			    std::string("BCALMCTruthLink") );

  registerProcessorParameter( "BackgroundFilename"  ,
			      " Name of input file for background energy density" ,
			      backgroundfilename ,
			      std::string("bg_aver.sv01-14-p00.mILD_o1_v05.E1000-B1b_ws.PBeamstr-pairs.I210000.root") );

  registerProcessorParameter( "EventClustersHistoRange",
			      " Event Clusters Histgram X axis range [n] (10 by default)",
			      _eventClustersHistoRange,
			      10 );

  //registerProcessorParameter("BeamCalHitCol","Collection of SimCalorimeterHits in BeamCal",_HitCol,std::string("BeamCalCollection"));
  registerProcessorParameter("BeamCalHitCol","Collection of CalorimeterHits in BeamCal",_HitCol,std::string("BeamCalCollection"));

  registerProcessorParameter( "WriteHistograms",
			      " Enable output of histogram root file",
			      _writeHistograms,
			      false );
}




void BCalReco::init(){

  gStyle->SetOptStat(0001110);

  printParameters();

  _nRun = 0 ; //total number of runs per job
  totalEvt=0; //total number of events per job
  totalEnergy=0; //total deposited sensor energy per job
  nHits=0;


  kk = 0;

  //st=new TStyle("Default", "Default Style");


  //read parameters from gear file
  const gear::CalorimeterParameters& bcparam = Global::GEAR->getBeamCalParameters(); //in bcparam are read all parameters from the gearfile; the gear file is already opened by default by Marlin (the name is known from the Marlin steering file) 

  geometry.segmentation = bcparam.getDoubleVals("phi_segmentation");
  nRings=geometry.segmentation.size();
  const DoubleVec& ext = bcparam.getExtent() ;
  geometry.rMin = ext[0] ;
  geometry.rMax = ext[1] ;
  geometry.zMin = ext[2] ;
  geometry.zMax = ext[3] ;
  geometry.DArinner = bcparam.getDoubleVal("dead_area_outer_r");
  geometry.sPhi = bcparam.getDoubleVal("cylinder_starting_phi");
  geometry.dPhi = bcparam.getDoubleVal("cylinder_spanning_phi");
  const gear::LayerLayout & l = bcparam.getLayerLayout() ;
  nLayers = l.getNLayers() ;
  geometry.xAngle=bcparam.getDoubleVal("beam_crossing_angle");
  geometry.pairMon=bcparam.getDoubleVal("pairsMonitorZ");
  geometry.cellSize0=l.getCellSize0(0);
	
  //precompute limiting theta values for particle-calo hit checking
  _thetaMin = atan ( geometry.rMin / fabs(geometry.zMax));
  _thetaMax = atan ( geometry.rMax / fabs(geometry.zMin));
	

  std::cout.precision(15);

  streamlog_out(MESSAGE) << "PARAMETERS READ FROM GEAR FILE: " << std::endl<<" \ncrossing angle= "<<fixed<<geometry.xAngle<<" \nstarting Phi= "<<fixed<<geometry.sPhi<<" \ndelta(spanning) Phi= "<<fixed<<geometry.dPhi<<" \ndead area radius= "<<fixed<<geometry.DArinner<<" \nnumber of BeamCal Layers= "<<nLayers<<" \nnumber of Rings= "<<nRings<<" \ninner radius= "<<fixed<<geometry.rMin<<" \nouter radius= "<<fixed<<geometry.rMax<<" \nBeamCal starting Z= "<<fixed<<geometry.zMin<<" \nBeamCal ending Z= "<<fixed<<geometry.zMax<<" \npairMon Z = "<<fixed<<geometry.pairMon<<" \nfixed Pad Size (size along the radius) = "<<fixed<<geometry.cellSize0<<std::endl;

  streamlog_out(MESSAGE)<<"\n\nSegmentation:\n";

  for(int j=0;j<nRings;j++) streamlog_out(MESSAGE)<<"Ring "<<j<<" : "<<geometry.segmentation[j]<<"\n";

  streamlog_out(MESSAGE)<<"\n\nLayer Properties:\n";
  for(int i=1;i<=nLayers;i++){
    geometry.distance[i]=l.getDistance(i-1);
    streamlog_out(MESSAGE)<<"BeamCal LayerNumber="<<i<<" , z="<<geometry.distance[i]<<" , geometry.cellSize0="<<l.getCellSize0(i-1)<<std::endl;
  }

  streamlog_out(MESSAGE)<<std::endl<<"INITIALIZATION\n\nLAYER [1,"<<nLayers<<"] = BEAMCAL\tLAYER "<<nLayers+1<<" = PAIR MONITOR\n\n";

  TFile ftmp(backgroundfilename.c_str());

  if ( ftmp.IsZombie() ) {
    streamlog_out(ERROR) << "Could not read data file. Exit." << endl;
    exit(1);
  }

  //      cout << "TEST GetEnergyErr 1" << endl;

  //      ftmp.ls();

  TTree *ttmp = (TTree*) ftmp.Get("tBcDensAverage");

  //      cout << "TEST GetEnergyErr 2" << endl;


  if ( ttmp->IsZombie() ) {
    streamlog_out(ERROR) << "Could not find data tree \"tBcDensAverage\". Exit." << endl;
    exit(1);
  }



  //      cout << "TEST GetEnergyErr 3" << endl;

  Nentr = ttmp->GetEntries();
  ftmp.Close();


  //INITIALIZATION OF HISTOGRAMS

  coordhitsxyP = new TH2F("coordhitsxyP", "Hits XY plane", 150, -150., 150., 150, -150., 150.); //
  coordhitsxyN = new TH2F("coordhitsxyN", "Hits XY plane", 150, -150., 150., 150, -150., 150.); //

  coordhitszP = new TH1F("coordhitszP", "Hits YZ", 100, 3000., 4000.); //
  coordhitszN = new TH1F("coordhitszN", "Hits YZ", 100, -4000., -3000.); //
  coordclusterzP = new TH1F("coordclusterzP", "Clusters Coord Z", 100, 3000., 4000.); //
  coordclusterzN = new TH1F("coordclusterzN", "Clusters Coord Z", 100, -4000., -3000.); //


  coordclusterxyP = new TH2F("coordclusterxyP", "Clusters XY plane", 150, -150., 150., 150, -150., 150.); //
  coordclusterxyN = new TH2F("coordclusterxyN", "Clusters XY plane", 150, -150., 150., 150, -150., 150.); //
	
  clusterPosHistoN = new TH2F("clusterPosHistoN", "Clusters XY plane", 150, -150., 150., 150, -150., 150.); //
  clusterPosHisto = new TH2F("clusterPosHisto", "Clusters XY plane", 150, -150., 150., 150, -150., 150.); //
  clusterPos3DHistoN = new TH3F("clusterPos3DHistoN", "Clusters 3D", 150, -4000., 4000., 150, -150., 150., 150, -150., 150.); //
  clusterPos3DHisto = new TH3F("clusterPos3DHisto", "Clusters 3D", 150, -4000., 4000., 150, -150., 150., 150, -150., 150.); //
  clusterPos3DHistoAll = new TH3F("clusterPos3DHistoAll", "Clusters 3D", 150, -4000., 4000., 150, -150., 150., 150, -150., 150.); //

  particlePosHistoN = new TH2F("particlePosHistoN", "Particles XY plane", 150, -150., 150., 150, -150., 150.); //
  particlePosHisto = new TH2F("particlePosHisto", "Particles XY plane", 150, -150., 150., 150, -150., 150.); //
  particlePos3DHistoN = new TH3F("particlePos3DHistoN", "Particles 3D", 150, -4000., 4000., 150, -150., 150., 150, -150., 150.); //
  particlePos3DHisto = new TH3F("particlePos3DHisto", "Particles 3D", 150, -4000., 4000., 150, -150., 150., 150, -150., 150.); //
  particlePos3DHistoAll = new TH3F("particlePos3DHistoAll", "Particles 3D", 150, -4000., 4000., 150, -150., 150., 150, -150., 150.); //

  MCparticlePosHistoN = new TH2F("MCparticlePosHistoN", "Particles XY plane", 150, -150., 150., 150, -150., 150.); //
  MCparticlePosHisto = new TH2F("MCparticlePosHisto", "Particles XY plane", 150, -150., 150., 150, -150., 150.); //
  MCparticlePos3DHistoN = new TH3F("MCparticlePos3DHistoN", "Particles 3D", 150, -4000., 4000., 150, -150., 150., 150, -150., 150.); //
  MCparticlePos3DHisto = new TH3F("MCparticlePos3DHisto", "Particles 3D", 150, -4000., 4000., 150, -150., 150., 150, -150., 150.); //
  MCparticlePos3DHistoAll = new TH3F("MCparticlePos3DHistoAll", "Particles 3D", 150, -4000., 4000., 150, -150., 150., 150, -150., 150.); //


  numberclusters = new TH1F("numberclusters", "Clusters", 10, 0., 10.); // number of clusters in an event
  numberparticles = new TH1F("numberparticles", "Particles", 10, 0., 10.); // number of particles in an event



  energyFW = new TH1F("energyp", "Energy in FW", 50, 0., 300.); //
  energyBW = new TH1F("energyn", "Energy in BW", 50, 0., 300.); //
  clusterNegEnergyHisto = new TH1F("clusterNegEnergyHisto", "Cluster Energy", 50, 0., 300.); //
  clusterPosEnergyHisto = new TH1F("clusterPosEnergyHisto", "Cluster Energy", 50, 0., 300.); //
  particleNegEnergyHisto = new TH1F("particleNegEnergyHisto", "Particle Energy", 50, 0., 300.); //
  particlePosEnergyHisto = new TH1F("particlePosEnergyHisto", "Particle Energy", 50, 0., 300.); //
	
  particleNegMomeHisto = new TH1F("particleNegMomeHisto", "Particle Momentum", 50, 0., 300.); //
  particlePosMomeHisto = new TH1F("particlePosMomeHisto", "Particle Momentum", 50, 0., 300.); //
	

  SigPosMean  = new TH1F("SigPosMean","Energy signal hits FW",100,0.,5.);
  SigNegMean  = new TH1F("SigNegMean","Energy signal hits BW",100,0.,5.);
  SigGausPosMean  = new TH1F("SigGausPosMean","Energy signal hits + bkgd fluctuations FW",100,0.,5.);
  SigGausNegMean  = new TH1F("SigGausNegMean","Energy signal hits + bkgd fluctuations BW",100,0.,5.);

  coordhitsxyP->SetFillColor(33);
  coordhitsxyP->SetFillStyle(4003);
  coordhitsxyP->GetXaxis()->SetTitle("x, mm");
  coordhitsxyP->GetYaxis()->SetTitle("y, mm");

  coordhitsxyN->SetFillColor(5);  //eperlayerp->SetFillStyle(5);
  coordhitsxyN->GetXaxis()->SetTitle("x, mm");
  coordhitsxyN->GetYaxis()->SetTitle("y, mm");



  coordclusterxyP->SetFillColor(33);
  coordclusterxyP->SetFillStyle(4003);
  coordclusterxyP->GetXaxis()->SetTitle("x, mm");
  coordclusterxyP->GetYaxis()->SetTitle("y, mm");


  coordhitszP->SetFillColor(5);   //eperlayerp->SetFillStyle(5);
  coordhitszP->GetXaxis()->SetTitle("z, mm");
  coordhitszP->GetYaxis()->SetTitle("events");

  numberclusters->SetFillColor(5);        //eperlayerp->SetFillStyle(5);
  numberclusters->GetXaxis()->SetTitle("N");
  numberclusters->GetYaxis()->SetTitle("events");
	
  numberparticles->SetFillColor(5);        //eperlayerp->SetFillStyle(5);
  numberparticles->GetXaxis()->SetTitle("N");
  numberparticles->GetYaxis()->SetTitle("events");
	

  coordhitszN->SetFillColor(5);   //eperlayerp->SetFillStyle(5);
  coordhitszN->GetXaxis()->SetTitle("z, mm");
  coordhitszN->GetYaxis()->SetTitle("events");


  energyFW->SetFillColor(5);      //eperlayerp->SetFillStyle(5);
  energyFW->GetXaxis()->SetTitle("E, GeV");
  energyFW->GetYaxis()->SetTitle("events");

  clusterNegEnergyHisto->SetFillColor(5);    //eperlayerp->SetFillStyle(5);
  clusterNegEnergyHisto->GetXaxis()->SetTitle("E, GeV");
  clusterNegEnergyHisto->GetYaxis()->SetTitle("events");

  clusterPosEnergyHisto->SetFillColor(5);    //eperlayerp->SetFillStyle(5);
  clusterPosEnergyHisto->GetXaxis()->SetTitle("E, GeV");
  clusterPosEnergyHisto->GetYaxis()->SetTitle("events");
	
  particleNegEnergyHisto->SetFillColor(5);    //eperlayerp->SetFillStyle(5);
  particleNegEnergyHisto->GetXaxis()->SetTitle("E, GeV");
  particleNegEnergyHisto->GetYaxis()->SetTitle("events");

  particlePosEnergyHisto->SetFillColor(5);    //eperlayerp->SetFillStyle(5);
  particlePosEnergyHisto->GetXaxis()->SetTitle("E, GeV");
  particlePosEnergyHisto->GetYaxis()->SetTitle("events");	
	
  particleNegMomeHisto->SetFillColor(5);    //eperlayerp->SetFillStyle(5);
  particleNegMomeHisto->GetXaxis()->SetTitle("E, GeV");
  particleNegMomeHisto->GetYaxis()->SetTitle("events");

  particlePosMomeHisto->SetFillColor(5);    //eperlayerp->SetFillStyle(5);
  particlePosMomeHisto->GetXaxis()->SetTitle("E, GeV");
  particlePosMomeHisto->GetYaxis()->SetTitle("events");	

  energyBW->SetFillColor(5);      //eperlayerp->SetFillStyle(5);
  energyBW->GetXaxis()->SetTitle("E, GeV");
  energyBW->GetYaxis()->SetTitle("events");

  SigPosMean->SetFillColor(5);       //eperlayerp->SetFillStyle(5);
  SigPosMean->GetXaxis()->SetTitle("GeV");
  SigPosMean->GetYaxis()->SetTitle(" ");

  SigNegMean->SetFillColor(5);       //eperlayerp->SetFillStyle(5);
  SigNegMean->GetXaxis()->SetTitle("GeV");
  SigNegMean->GetYaxis()->SetTitle(" ");

  SigGausNegMean->SetFillColor(5);       //eperlayerp->SetFillStyle(5);
  SigGausNegMean->GetXaxis()->SetTitle("GeV");
  SigGausNegMean->GetYaxis()->SetTitle(" ");

  SigGausPosMean->SetFillColor(5);       //eperlayerp->SetFillStyle(5);
  SigGausPosMean->GetXaxis()->SetTitle("GeV");
  SigGausPosMean->GetYaxis()->SetTitle(" ");


  //number of cluster in an event
  _numClustersHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D("numClustersHisto",
									       "Number of clusters in event [n]", _eventClustersHistoRange, 1, _eventClustersHistoRange);
  assert (_numClustersHisto);


  //number of particles in an event
  _numParticHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D("numParticHisto",
									     "Number of particles in event [n]", _eventClustersHistoRange, 1, _eventClustersHistoRange);
  assert (_numParticHisto);

  streamlog_out(MESSAGE)<<"\n\n\nFINISHED INITIALIZATION\n\n\n"<<std::endl;

}



void BCalReco::processRunHeader( LCRunHeader * run ) {

  if(_nRun!=0)
    //                streamlog_out(MESSAGE)<<"RUN "<<_nRun<<" HAS "<<_nEvt<<" TOTAL EVENTS "<<runEnergy<< "GEV TOTAL DEPOSITED SENSOR ENERGY. "<<"ZNMIN="<<znmin<<"mm   ZNMAX="<<znmax<<"mm  ZPMIN="<<zpmin<<"mm  ZPMAX="<<zpmax<<"mm"<<std::endl;
    _nRun++;
  _nEvt = 0; //number of events for each run
  runEnergy=0; //total deposited sensor energy for each run
  znmin=0;
  znmax=-DBL_MAX;
  zpmin=DBL_MAX;
  zpmax=0;


}






void BCalReco::processEvent(LCEvent * evt){


  std::vector<ClusterImpl*> clustersNew;
  std::vector<ReconstructedParticleImpl*> particlesNew;


  //Create output collection
  LCCollectionVec* BCALClusters = new LCCollectionVec(LCIO::CLUSTER);
  LCCollectionVec* BCALCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  //aici

  int IP = 3595;
  float delta = 4.0; //layer thickness, mm
  double coordX = 0., coordY = 0.,coordZ = 0.;
  double radius = 0.;

  //Get elements for collection to study
  LCCollection* col = evt->getCollection( _colName ) ;

  if(col != 0) {

    int nMCP = col->getNumberOfElements();
	    
    //cout << "nMCP = " << nMCP << endl;

    for (int i = 0; i < nMCP; i++) {

      MCParticle* p = dynamic_cast<MCParticle*>( col->getElementAt(i) );

      if(!p->isCreatedInSimulation()) {
			   
	int pdg = p->getPDG();
				
	//cout << "partic ID = " << pdg << endl;
				   
	if(abs(pdg) == 11)  {

	  double* momentum = (double*) p->getMomentum();
	  double pTheta;

	  //calculate its angle from momentum
	  pTheta = atan ( sqrt(momentum[0]*momentum[0]+momentum[1]*momentum[1]) / fabs(momentum[2]));
	  //pTheta = atan2 ( sqrt(momentum[0]*momentum[0]+momentum[1]*momentum[1]) , momentum[2]);

	  //cout <<"Theta min = " << _thetaMin << endl;
	  //cout <<"Theta max = " << _thetaMax << endl;
	  //cout <<"Theta particle = " << pTheta << endl;
						


	  //if the particle hit the calorimeter
	  //    if ( pTheta < _thetaMax && pTheta > _thetaMin ) {

	  //check side
	  int MCpartIsNegative = p->getMomentum()[2] < 0;
					
									     				
	  //depending on side :
	  if (MCpartIsNegative)
	    {
					
					 			                
	      //get energy
	      double MCpartENeg = (double) p->getEnergy();
					
	      coordZ = - (double) IP - 5.*delta;
	      radius = (coordZ * MCpartENeg) / p->getMomentum()[2];
					 
	      coordX = radius * p->getMomentum()[0] / MCpartENeg ;
	      coordY = radius * p->getMomentum()[1] / MCpartENeg ;
					  
                                         
	      //fill position histograms
	      MCparticlePosHistoN->Fill(coordX, coordY,1);

	      MCparticlePos3DHistoN->Fill(coordZ, coordX, coordY,1);
	      MCparticlePos3DHistoAll->Fill(coordZ, coordX, coordY,1);



	    }
	  else //if (MCpartIsNegative)
	    {
	      //same as for negative side
						
	      //get energy
	      double MCpartEPos = (double) p->getEnergy();
					
	      coordZ = (double) IP + 5.*delta;
	      radius = (coordZ * MCpartEPos) / p->getMomentum()[2];
					 
	      coordX = radius * p->getMomentum()[0] / MCpartEPos ;
	      coordY = radius * p->getMomentum()[1] / MCpartEPos ;  
                                     	
						
	      //fill position histograms
	      MCparticlePosHisto->Fill(coordX, coordY,1);

	      MCparticlePos3DHisto->Fill(coordZ, coordX, coordY,1);
	      MCparticlePos3DHistoAll->Fill(coordZ, coordX, coordY,1);


	    } //else

	  //} // particle hit beamcal
	} // only electrons
      } //only source particles
    } // end of MCparticle loop

  } // end col




  //aici

  //CONSTRUCTION OF DYNAMIC ARRAY TO STORE THE INFORMATION FROM THE DETECTOR FOR THE RECONSTRUCTION CODE

  cells = new (nothrow) cell ** [nLayers+2];
  if(cells==0) {streamlog_out(ERROR)<<"No more space for Layer allocation"; exit(1);}

  for(int i=1;i<=nLayers+1;i++){
    cells[i]=new (nothrow) cell * [nRings];
    if(cells[i]==0) {streamlog_out(ERROR)<<"No more space for Ring allocation at layer "<<i; exit(1);}
  }

  for(int j=0;j<nRings;j++){

    int N=(int)(round(2*TMath::Pi()/geometry.segmentation[j]));

    //streamlog_out(DEBUG2)<<"\n\n\nRING "<<j<<"\nNumber of assigned cells: "<<N;

    for(int i=1;i<=nLayers+1;i++){

      cells[i][j]=new (nothrow) cell [N];
      if(cells[i][j]==0) {streamlog_out(ERROR)<<"\nNo more space for cell allocation at layer="<<i<<" ring="<<j; exit(1);}
    }

    int k=0;
    double phi=geometry.sPhi;

    //               streamlog_out(DEBUG2)<<"check sPhi "<< phi <<N;


    do{
      if(k>=N) {
	streamlog_out(DEBUG2)<<"\nNot enough cells allocated for ring="<<j<<" phi="<<phi<<" k="<<k<<"\n";
	fflush(stdout);
	exit(1);
      }

      if ((geometry.sPhi+geometry.dPhi-epsPhi<= phi) && (geometry.rMin+j*geometry.cellSize0 <= geometry.DArinner+epsR )) 
	{
	  //streamlog_out(DEBUG2)<< "\nArrived to Dead Area.\tRing "<<j<<" LastConsideredPhi = "<<phi-geometry.segmentation[j]<<" NextPhi =  "<<phi<< " geometry.sPhi+geometry.dPhi-NextPhi = "<<-phi+geometry.sPhi+geometry.dPhi<<" geometry.sPhi+geometry.dPhi-LastConsideredPhi = "<<-phi+geometry.sPhi+geometry.dPhi+geometry.segmentation[j];
	  break;
	}


      for(int i=1;i<=nLayers+1;i++){

	cells[i][j][k].sRin=geometry.rMin+j*geometry.cellSize0;
	cells[i][j][k].sRout=geometry.rMin+(j+1)*geometry.cellSize0;
	if(i==nLayers+1)
	  cells[i][j][k].sZstart=geometry.pairMon;
	else
	  cells[i][j][k].sZstart=geometry.distance[i];
	cells[i][j][k].sZend=0;
	cells[i][j][k].sEdepNeg=0.;
	cells[i][j][k].sEdepPos=0.;
	cells[i][j][k].sSphi=phi;
	if(phi>=2*TMath::Pi()) cells[i][j][k].sSphi-=2*TMath::Pi();
	cells[i][j][k].sDphi=geometry.segmentation[j];
	cells[i][j][k].sPos[0]=i;
	cells[i][j][k].sPos[1]=j;
	cells[i][j][k].sPos[2]=k;

      }

      phi=phi+geometry.segmentation[j];
      k++;

    }while(phi<geometry.sPhi+2*TMath::Pi()-epsPhi);

    //streamlog_out(DEBUG2)<<"\nOK\tRing "<<j<<" LastConsideredPhi = "<<phi-geometry.segmentation[j]<<" NextPhi = "<<phi<< " NextPhi-geometry.sPhi =  "<<phi-geometry.sPhi-2*TMath::Pi()<<" LastConsideredPhi-sPhi = "<<phi-geometry.sPhi-2*TMath::Pi()-geometry.segmentation[j];

    nbPhis[j]=k;
  }


  if(_nEvt == 1){
    streamlog_out(DEBUG2)<<"\n\n\nFINAL NUMBER OF PADS PER RING:\n";
    for(int j=0;j<nRings;j++)
      streamlog_out(DEBUG2)<<"ring "<<j<<" = "<<nbPhis[j]<<"\n";
  }

  try{//getCollection

    _nEvt ++;
    totalEvt++;

    LCCollection* inParVecSim = evt->getCollection(_HitCol); //save in inParVecSim all data related to the current event from the collection
    //CellIDDecoder<SimCalorimeterHit> mydecoder(inParVecSim);
    CellIDDecoder<CalorimeterHit> mydecoder(inParVecSim);

    //    streamlog_out(DEBUG2)<<"step 1 "<<"\n";

    try{//getElementAt
      evtEnergy=0; //total energy deposition for each event
      nHits=inParVecSim->getNumberOfElements(); //get total number of hits from current event
      for(int j=0; j<nHits; j++){
	//SimCalorimeterHit *hit = dynamic_cast<SimCalorimeterHit*>(inParVecSim->getElementAt(j)); // read each hit
				
	CalorimeterHit *hit = dynamic_cast<CalorimeterHit*>(inParVecSim->getElementAt(j)); // read each hit
	//CalorimeterHit *hit = dynamic_cast<CalorimeterHit*>(inParVecSim->getElementAt(j)); // read each hit
				
	if(hit){
                    
	  side    = mydecoder(hit)[ "S-1" ] ;
	  //if(side!=0) continue;
	  //side=0 -> pos, side=1 -> neg
	  cylind  = mydecoder(hit)[ "I" ] ; //r 0-16
	  sector  = mydecoder(hit)[ "J" ] ; //phi 0-max_nb_of_pads
	  layer   = mydecoder(hit)[ "K" ]; //z with the Mokka codification: nLayers=30; layer=[1,30] for BeamCal; layer=31 for PairMonitor
				
				

	  double pos[3] = {hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] } ; //x,y,z
							

	  E = hit->getEnergy();
				
	  //streamlog_out(DEBUG2) << "CHECK:::" << E << endl;
                                   

	  //IF ROTSATION IS NEEDED
	  /*
	    if( posTemp[2] < 0) geometry.xAngle=-TMath::Abs(geometry.xAngle);
	    else geometry.xAngle=TMath::Abs(geometry.xAngle);
	    pos[0] = TMath::Cos(geometry.xAngle/2)*posTemp[0]-TMath::Sin(geometry.xAngle/2)*posTemp[2];
	    pos[1] = posTemp[1];
	    pos[2] = TMath::Sin(geometry.xAngle/2)*posTemp[0]+TMath::Cos(geometry.xAngle/2)*posTemp[2];
	  */


	  if(cylind<0 || cylind>nRings-1 || layer<1 || layer>nLayers+1 || sector<0 || sector>nbPhis[cylind]-1) {
	    streamlog_out(DEBUG2)<<" WRONG  I or J or K\n";
	    continue;
	  }


	  if(side==1) cells[layer][cylind][sector].sEdepNeg+=E;
	  if(side==0) cells[layer][cylind][sector].sEdepPos+=E;



	  if( pos[2] < 0) geometry.xAngle=-TMath::Abs(geometry.xAngle);
	  else geometry.xAngle=TMath::Abs(geometry.xAngle);
	  pos[0] = TMath::Cos(geometry.xAngle/2000)*pos[0]-TMath::Sin(geometry.xAngle/2000)*pos[2];
	  pos[1] = pos[1];
	  pos[2] = TMath::Sin(geometry.xAngle/2000)*pos[0]+TMath::Cos(geometry.xAngle/2000)*pos[2];



	  if(side==0){
	    coordhitsxyP->Fill(pos[0],pos[1],1);
	    coordhitszP->Fill(pos[2],1);
	  }
	  else{
	    coordhitsxyN->Fill(pos[0],pos[1],1);
	    coordhitszN->Fill(pos[2],1);
	  }

	  fflush(stdout);


	  evtEnergy+=E;
	}//for
	// streamlog_out(DEBUG2) << "Run"<<_nRun<<" Event "<<_nEvt<<": Number of hits= "<< nHits<<" Energy deposition for event = " << evtEnergy << "  [GeV]"<<std::endl;

	runEnergy+=evtEnergy;
	totalEnergy+=evtEnergy;
      } // if (hit)	


    }//try getElementAt)
    catch(Exception& e){
      fflush(stdout);
      streamlog_out(DEBUG9) << "Exception in at least one getHit of Event "<<_nEvt <<": " << e.what() << std::endl;
    }



  }//try getCollection
  catch(Exception& e){
    fflush(stdout);
    streamlog_out(DEBUG9) <<"Exception in Event "<<_nEvt<<" with get BeamCalCollection: " << e.what() << endl;
  }


  TFile f(backgroundfilename.c_str());
  if ( f.IsZombie() ) {
    streamlog_out(ERROR) << "Could not read data file. Exit." << endl;
    exit(1);
  }


  t = (TTree*)f.Get("tBcDensAverage");


  if(t->IsZombie()) 
    streamlog_out(DEBUG9) <<"bgcoeff: Tree \"tBcDensAverage\" not found"<<endl;
    


  //      cout << "TEST GetEnergyErr 4" << endl;


  t->SetBranchAddress("sPos",    sPos);
  t->SetBranchAddress("sRin",   &sRin);  // cm
  t->SetBranchAddress("sRout",  &sRout); // cm
  t->SetBranchAddress("sZstart",&sZstart);
  t->SetBranchAddress("sZend",  &sZend);
  t->SetBranchAddress("sSphi",  &sSphi); // rad
  t->SetBranchAddress("sDphi",  &sDphi); // rad
  t->SetBranchAddress("sEdep",  &sEdep); // GeV
  t->SetBranchAddress("sCellArea", &sCellArea);
  t->SetBranchAddress("sEnDens", &sEnDens);
  t->SetBranchAddress("sEnDensErr", &sEnDensErr);


  for (int i=0; i<Nentr; i++){
    bgc[i].sPos[0] = 0;
    bgc[i].sPos[1] = 0;
    bgc[i].sPos[2] = 0;
    bgc[i].sRin    = 0.;
    bgc[i].sRout   = 0.;
    bgc[i].sZstart = 0.;
    bgc[i].sZend   = 0.;
    bgc[i].sSphi   = 0.;
    bgc[i].sDphi   = 0.;
    bgc[i].sEdep   = 0.;
    bgc[i].sCellArea = 0.;
    bgc[i].sEnDens = 0.;
    bgc[i].sEnDensErr = 0.;
  }

  sPos[0] = 0; sPos[1] = 0; sPos[2] = 0;
  sRin=0.;
  sRout=0.;
  sZstart=0.;
  sZend=0.;
  sSphi=0.; sDphi=0.;
  sEdep=0.;
  sCellArea=0.;sEnDens=0.;sEnDensErr=0.;

    
  for(int i=0; i<Nentr; i++){
    t->GetEntry(i);
    bgc[i].sPos[0] = sPos[0];
    bgc[i].sPos[1] = sPos[1];
    bgc[i].sPos[2] = sPos[2];
    bgc[i].sRin = sRin;
    bgc[i].sRout = sRout;
    bgc[i].sZstart = sZstart;
    bgc[i].sZend = sZend;
    bgc[i].sSphi = sSphi;
    bgc[i].sDphi = sDphi;
    bgc[i].sEdep = sEdep;
    bgc[i].cell_area = sCellArea;  // mm2
    bgc[i].sEnDens = sEnDens; // GeV/mm2
    bgc[i].sEnDensErr = sEnDensErr;
  }


  //      cout << "TEST 2 " << "Pad = " <<bgc[100].sPos[2] << endl;


  int index = 0;
  double EdepPos[36810],EdepNeg[36810];
  double EdepPosTmp[36810],EdepNegTmp[36810];
  double EdepPosNoiseTmp[36810],EdepNegNoiseTmp[36810];
  double EdepTotPos = 0.,EdepTotNeg = 0.;
  double EdepNoiseTotPos = 0.,EdepNoiseTotNeg = 0.;
  TRandom3 r;
  for(int i = 0; i < 36810; i++){                
    EdepPos[i] = 0.;
    EdepNeg[i] = 0.;
    EdepPosTmp[i] = 0.;
    EdepNegTmp[i] = 0.;
    EdepPosNoiseTmp[i] = 0.;
    EdepNegNoiseTmp[i] = 0.;
  }
			

  for (int i=1; i<=nLayers; ++i) {
    for (int j=0; j<nRings; ++j) {
      for (int k=0; k<nbPhis[j]; ++k){
	EdepPos[index] = cells[i][j][k].sEdepPos;
	EdepNeg[index] = cells[i][j][k].sEdepNeg;
	index++;
      }
    }
  }

  //		if(_nEvt == 2){
          
  //                cout << "Nentr = " <<Nentr << " " << "Index = " <<index << endl;
  //                cout << "nLayers = " <<nLayers << " " << "nRings = " <<nRings << endl;

  for(int i = 0; i < index - 1; i++)
    {
		    

      EdepPosNoiseTmp[i] = EdepPos[i] + r.Gaus(0,bgc[2*i].sEnDensErr);
      EdepNegNoiseTmp[i] = EdepNeg[i] + r.Gaus(0,bgc[2*i+1].sEnDensErr);
			


      if(EdepPosNoiseTmp[i] > bgc[2*i].sEnDensErr){
	EdepPosNoiseTmp[i] = EdepPosNoiseTmp[i];
	EdepPosTmp[i] = EdepPos[i];
      }
      else {
	EdepPosTmp[i] = 0.;
	EdepPosNoiseTmp[i] = 0.;
      }
      if(EdepNegNoiseTmp[i] > bgc[2*i+1].sEnDensErr){
	EdepNegNoiseTmp[i] = EdepNegNoiseTmp[i];
	EdepNegTmp[i] = EdepNeg[i];
      }
      else {
	EdepNegTmp[i] = 0.;
	EdepNegNoiseTmp[i] = 0.;
      }
      EdepTotPos += EdepPosTmp[i]; 
      EdepTotNeg += EdepNegTmp[i];
      EdepNoiseTotPos += EdepPosNoiseTmp[i]; 
      EdepNoiseTotNeg += EdepNegNoiseTmp[i];

    }
  //                }

  if(EdepTotPos > 0)
    SigPosMean->Fill(EdepTotPos,1);
  if(EdepTotNeg > 0)
    SigNegMean->Fill(EdepTotNeg,1);
  if(EdepNoiseTotPos > 0)
    SigGausPosMean->Fill(EdepNoiseTotPos,1);
  if(EdepNoiseTotNeg > 0)
    SigGausNegMean->Fill(EdepNoiseTotNeg,1);

  int index2 = 0;
  for (int i=1; i<=nLayers; ++i) {
    for (int j=0; j<nRings; ++j) {
      for (int k=0; k<nbPhis[j]; ++k){
	cells[i][j][k].sEdepPos = EdepPosNoiseTmp[index2];
	cells[i][j][k].sEdepNeg = EdepNegNoiseTmp[index2];
                                                   
	index2++;
      }
    }
  }



  //instantiate reconstruction class

  bc_en = new BCalReconstruction();


  //reconstruction

  celule = new BCalReconstruction::CellType** [nLayers];
  for (int i=1; i<nLayers; ++i) {
    celule[i]=new BCalReconstruction::CellType* [nRings];
    for (int j=0; j<nRings; ++j) {
      celule[i][j]=new BCalReconstruction::CellType [nbPhis[j]];
      for (int k=0; k<nbPhis[j]; ++k){
	celule[i][j][k].sRin=cells[i][j][k].sRin;
	celule[i][j][k].sRout=cells[i][j][k].sRout;
	celule[i][j][k].sZstart=cells[i][j][k].sZstart;
	celule[i][j][k].sZend=cells[i][j][k].sZend;
	celule[i][j][k].sEdepNeg=cells[i][j][k].sEdepNeg;
	celule[i][j][k].sEdepPos=cells[i][j][k].sEdepPos;
	celule[i][j][k].sSphi=cells[i][j][k].sSphi;
	celule[i][j][k].sDphi=cells[i][j][k].sDphi;
	celule[i][j][k].sPos[0]=cells[i][j][k].sPos[0];
	celule[i][j][k].sPos[1]=cells[i][j][k].sPos[1];
	celule[i][j][k].sPos[2]=cells[i][j][k].sPos[2];
      }
    }
  }


  obiect = bc_en->GetReconstrCoordinates(nLayers,nRings,nbPhis,celule);

  //cout << "Reconstructed side = " <<" " << obiect.side[0] << " " << " " << obiect.side[1] << endl;
		  
  float momPClu[3] = {0.,0.,0.}, energyClu = 0., thetaClu = 0., phiClu = 0., radClu = 0.;
  const float m = 0. ;
  const float q = 1e+19 ;
  float posit[3] = {0.,0.,0.};

  if(obiect.side[0] == 1){
    coordclusterxyP->Fill(obiect.CoordX[0],obiect.CoordY[0],1);
    coordclusterzP->Fill(obiect.CoordZ[0],1);
    energyFW->Fill(obiect.RecEne[0],1); 
		  
		  
    ReconstructedParticleImpl* particle = new ReconstructedParticleImpl;
    ClusterImpl* cluster = new ClusterImpl;

		
    //            cout << "Charge = " << q << endl;


    particle->setMass(m) ;
    particle->setCharge(q) ;

    cluster->setEnergy( obiect.RecEne[0] );
                
    posit[0] = obiect.CoordX[0];
    posit[1] = obiect.CoordY[0];
    posit[2] = obiect.CoordZ[0];

    cluster->setPosition( posit );

    radClu = sqrt(posit[0]*posit[0] + posit[1]*posit[1] + posit[2]*posit[2]) ;
		
    thetaClu = acos(posit[2]/radClu) ;
    phiClu = atan2 ( posit[1],posit[0] ) ;

    energyClu = obiect.RecEne[0] ;

    momPClu[0] = energyClu * sin ( thetaClu ) * cos ( phiClu ) ;
    momPClu[1] = energyClu * sin ( thetaClu ) * sin ( phiClu ) ;
    momPClu[2] = energyClu * cos ( thetaClu ) ;
		
    //cout << "FW::Check cosine phi =" << cos ( phiClu ) << endl;

    particle->setMomentum ( momPClu ) ;
    particle->setEnergy ( energyClu ) ;

    particle->addCluster( cluster ) ;


    BCALClusters->addElement(cluster);
    BCALCol->addElement(particle);
		  
		                   
  }	  
	
  if(obiect.side[1] == 1){
    coordclusterxyN->Fill(obiect.CoordX[1],obiect.CoordY[1],1);
    coordclusterzN->Fill(obiect.CoordZ[1],1);
    energyBW->Fill(obiect.RecEne[1],1);
		  
    ReconstructedParticleImpl* particle = new ReconstructedParticleImpl;
    ClusterImpl* cluster = new ClusterImpl;

		
    //            cout << "Charge = " << q << endl;



    particle->setMass(m) ;
    particle->setCharge(q) ;

    cluster->setEnergy( obiect.RecEne[1] );
                
    posit[0] = obiect.CoordX[1];
    posit[1] = obiect.CoordY[1];
    posit[2] = obiect.CoordZ[1];

    cluster->setPosition( posit );
		

    radClu = sqrt(posit[0]*posit[0] + posit[1]*posit[1] + posit[2]*posit[2]) ;
		
    thetaClu = acos(posit[2]/radClu) ;
    phiClu = atan2 ( posit[1],posit[0] ) ;

    energyClu = obiect.RecEne[1] ;

    momPClu[0] = energyClu * sin ( thetaClu ) * cos ( phiClu ) ;
		
		
    //cout << "BW::Check cosine phi =" <<  posit[1]/posit[0] << endl;
		
    momPClu[1] = energyClu * sin ( thetaClu ) * sin ( phiClu ) ;
    momPClu[2] = energyClu * cos ( thetaClu ) ;

    particle->setMomentum ( momPClu ) ;
    particle->setEnergy ( energyClu ) ;

    particle->addCluster( cluster ) ; 


    BCALClusters->addElement(cluster);
    BCALCol->addElement(particle);
		  
		                   		  
		  
		    
  }
	

  evt->addCollection(BCALClusters, _BCALClustersName);
  evt->addCollection(BCALCol, _BCALcolName);

  for (int i = 1; i < nLayers; ++i) {
    for (int j = 0; j < nRings; ++j)
      delete [] celule[i][j];
    delete [] celule[i];
  }
  delete [] celule;

  for (int i = 1; i < nLayers + 2; ++i) {
    for (int j = 0; j < nRings; ++j)
      delete [] cells[i][j];

    delete [] cells[i];
  }
  delete [] cells;


  delete bc_en;





  f.Close();


}

void BCalReco::check(LCEvent *evt){

  int j;
  int numElements;
	
	
  //cout << "BCalReco::check() : event# " << _nEvt << endl;


  //read clusters and fill basic histograms
  // split clusters into positive and negative sides
  //position of clusters
  //reconstructed energy
  try
    {
      LCCollection * cluster = evt->getCollection( _BCALClustersName );
      if (cluster)
	{
	  numElements = cluster->getNumberOfElements()  ;
			
	  //cout << "numElements = " << numElements << endl;

	  //first very basic histo
                 
	  numberclusters->Fill(numElements,1);


	  for (j = 0; j < numElements ; j++)
	    {
	      ClusterImpl * cl = dynamic_cast<ClusterImpl*>( cluster->getElementAt(j) );
	      if (cl)
		{
		  //check side
		  int clIsNegative = cl->getPosition()[2] < 0;

		  //depending on side :
		  if (clIsNegative)
		    {

		      //get energy and put it into a histo
		      double clENeg = (double) cl->getEnergy();

		      clusterNegEnergyHisto->Fill(clENeg,1);
                                         
		      //fill position histograms
		      clusterPosHistoN->Fill(cl->getPosition()[0], cl->getPosition()[1],1);

		      clusterPos3DHistoN->Fill(cl->getPosition()[2], cl->getPosition()[0], cl->getPosition()[1],1);
		      clusterPos3DHistoAll->Fill(cl->getPosition()[2], cl->getPosition()[0], cl->getPosition()[1],1);



		    }
		  else //if (clIsNegative)
		    {
		      //same as for negative side
						
		      //get energy and put it into a histo
		      double clEPos = (double) cl->getEnergy();

		      clusterPosEnergyHisto->Fill(clEPos,1);
                                         	
		      //fill position histograms
		      clusterPosHisto->Fill(cl->getPosition()[0], cl->getPosition()[1],1);

		      clusterPos3DHisto->Fill(cl->getPosition()[2], cl->getPosition()[0], cl->getPosition()[1],1);
		      clusterPos3DHistoAll->Fill(cl->getPosition()[2], cl->getPosition()[0], cl->getPosition()[1],1);

		    } //else
		} //if (cl)
	    } //for(j = 0; j < numElements ; j++)





	} //if (cluster)
    } //try

  catch(DataNotAvailableException &e)
    {
      //  _numClustersHisto->fill( 0 ); // no clusters found in event
    } //catch


  int numPartElements;

  //read particles and fill basic histograms
  // split particles into positive and negative sides
  //position of particles
  //reconstructed energy
  try
    {
      LCCollection * particle = evt->getCollection( _BCALcolName );
      if (particle)
	{
	  numPartElements = particle->getNumberOfElements()  ;

	  //first very basic histo
                       
	  numberparticles->Fill(numPartElements,1);
			 
			 
	  int IP = 3595;  
	  float delta = 4.0; //layer thickness, mm 
	  double coordX = 0., coordY = 0.,coordZ = 0.;
	  double radius = 0.;
			 
	  for (j = 0; j < numPartElements ; j++)
	    {
	      ReconstructedParticleImpl * part = dynamic_cast<ReconstructedParticleImpl*>( particle->getElementAt(j) );
	      if (part)
		{
		  //check side
		  int partIsNegative = part->getMomentum()[2] < 0;
					
									     				
		  //depending on side :
		  if (partIsNegative)
		    {
					
					 			                
		      //get energy and put it into a histo
		      double partENeg = (double) part->getEnergy();
					
		      //get momentum and put it into a histo
		      double partNegP = sqrt(part->getMomentum()[0] * part->getMomentum()[0] + part->getMomentum()[1] * part->getMomentum()[1] + part->getMomentum()[2] * part->getMomentum()[2]);

		      particleNegEnergyHisto->Fill(partENeg,1);
		      particleNegMomeHisto->Fill(partNegP,1);
					 
		      coordZ = - (double) IP - 5.*delta;
		      radius = (coordZ * partENeg) / part->getMomentum()[2];
					 
		      //cout << "BW::Px = " << part->getMomentum()[0] << endl;
		      //cout << "BW::Py = " << part->getMomentum()[1] << endl;
					 
		      coordX = radius * part->getMomentum()[0] / partENeg ;
		      coordY = radius * part->getMomentum()[1] / partENeg ;
					  
                                         
		      //fill position histograms
		      particlePosHistoN->Fill(coordX, coordY,1);

		      particlePos3DHistoN->Fill(coordZ, coordX, coordY,1);
		      particlePos3DHistoAll->Fill(coordZ, coordX, coordY,1);



		    }
		  else //if (partIsNegative)
		    {
		      //same as for negative side
						
		      //get energy and put it into a histo
		      double partEPos = (double) part->getEnergy();
					
		      //get momentum and put it into a histo
		      double partPosP = sqrt(part->getMomentum()[0] * part->getMomentum()[0] + part->getMomentum()[1] * part->getMomentum()[1] + part->getMomentum()[2] * part->getMomentum()[2]);

		      particlePosEnergyHisto->Fill(partEPos,1);
		      particlePosMomeHisto->Fill(partPosP,1);
					 
		      coordZ = (double) IP + 5.*delta;
		      radius = (coordZ * partEPos) / part->getMomentum()[2];
					 
		      //cout << "FW::Px = " << part->getMomentum()[0] << endl;
		      //cout << "FW::Py = " << part->getMomentum()[1] << endl;
					 
		      coordX = radius * part->getMomentum()[0] / partEPos ;
		      coordY = radius * part->getMomentum()[1] / partEPos ;  
                                     	
						
		      //fill position histograms
		      particlePosHisto->Fill(coordX, coordY,1);

		      particlePos3DHisto->Fill(coordZ, coordX, coordY,1);
		      particlePos3DHistoAll->Fill(coordZ, coordX, coordY,1);


		    } //else
		} //if (part)
	    } //for(j = 0; j < numPartElements ; j++)


	} //if (particle)
    } //try

  catch(DataNotAvailableException &e)
    {
      //   _numParticHisto->fill( 0 ); // no clusters found in event
    } //catch




}








void BCalReco::end(){

  //        streamlog_out(MESSAGE)<<"\n\n\nEND OF RUN "<<_nRun<<":\n"<<_nEvt<<" TOTAL EVENTS\nTOTAL DEPOSITED SENSOR ENERGY = "<<runEnergy<<" [GeV]\nZNMIN = "<<znmin<<" [mm]\nZNMAX = "<<znmax<<" [mm]\nZPMIN = "<<zpmin<<" [mm]\nZPMAX = "<<zpmax<<" [mm]"<<std::endl;


  streamlog_out(MESSAGE)<<"\n\n\nEND OF JOB:\n"<<totalEvt<<" EVENTS IN "<<_nRun << " RUNS\nTOTAL DEPOSITED SENSOR ENERGY = "<<totalEnergy<<" [GeV]"<<std::endl;

  // write the histos out

  if (_writeHistograms) {
    file  = TFile::Open("bcalrootfile.root","RECREATE");

    coordhitsxyP->Write();
    coordhitsxyN->Write();
    coordclusterxyP->Write();
    coordclusterzP->Write();
    coordclusterzN->Write();
    coordclusterxyN->Write();
    clusterPosHistoN->Write();
    clusterPosHisto->Write();
    clusterPos3DHistoN->Write();
    clusterPos3DHisto->Write();
    clusterPos3DHistoAll->Write();
    particlePosHistoN->Write();
    particlePosHisto->Write();
    particlePos3DHistoN->Write();
    particlePos3DHisto->Write();
    particlePos3DHistoAll->Write();

    MCparticlePosHistoN->Write();
    MCparticlePosHisto->Write();
    MCparticlePos3DHistoN->Write();
    MCparticlePos3DHisto->Write();
    MCparticlePos3DHistoAll->Write();

    coordhitszP->Write();
    coordhitszN->Write();
    energyFW->Write();
    clusterNegEnergyHisto->Write();
    clusterPosEnergyHisto->Write();
    particleNegEnergyHisto->Write();
    particlePosEnergyHisto->Write();
    particleNegMomeHisto->Write();
    particlePosMomeHisto->Write();
    energyBW->Write();
    numberclusters->Write();
    numberparticles->Write();
    SigPosMean->Write();
    SigNegMean->Write();
    SigGausPosMean->Write();
    SigGausNegMean->Write();

    file->Write();
    file->Close();
  }

}



