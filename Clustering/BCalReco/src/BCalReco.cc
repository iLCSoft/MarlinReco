
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
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimCalorimeterHit.h>

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

        registerProcessorParameter( "EventClustersHistoRange",
                                " Event Clusters Histgram X axis range [n] (10 by default)",
                                _eventClustersHistoRange,
                                10 );

        registerProcessorParameter("BeamCalHitCol","Collection of SimCalorimeterHits in BeamCal",_HitCol,std::string("BeamCalCollection"));
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
        geometry.sPhi = bcparam.getDoubleVal("cylinder_staring_phi");
        geometry.dPhi = bcparam.getDoubleVal("cylinder_spanning_phi");
        const gear::LayerLayout & l = bcparam.getLayerLayout() ;
        nLayers = l.getNLayers() ;
        geometry.xAngle=bcparam.getDoubleVal("beam_crossing_angle");
        geometry.pairMon=bcparam.getDoubleVal("pairsMonitorZ");
        geometry.cellSize0=l.getCellSize0(0);

        std::cout.precision(15);

        std::cout << "PARAMETERS READ FROM GEAR FILE: " << std::endl<<" \ncrossing angle= "<<fixed<<geometry.xAngle<<" \nstarting Phi= "<<fixed<<geometry.sPhi<<" \ndelta(spanning) Phi= "<<fixed<<geometry.dPhi<<" \ndead area radius= "<<fixed<<geometry.DArinner<<" \nnumber of BeamCal Layers= "<<nLayers<<" \nnumber of Rings= "<<nRings<<" \ninner radius= "<<fixed<<geometry.rMin<<" \nouter radius= "<<fixed<<geometry.rMax<<" \nBeamCal starting Z= "<<fixed<<geometry.zMin<<" \nBeamCal ending Z= "<<fixed<<geometry.zMax<<" \npairMon Z = "<<fixed<<geometry.pairMon<<" \nfixed Pad Size (size along the radius) = "<<fixed<<geometry.cellSize0<<std::endl;

        std::cout<<"\n\nSegmentation:\n";

        for(int j=0;j<nRings;j++) std::cout<<"Ring "<<j<<" : "<<geometry.segmentation[j]<<"\n";

        std::cout<<"\n\nLayer Properties:\n";
        for(int i=1;i<=nLayers;i++){
                geometry.distance[i]=l.getDistance(i-1);
                std::cout<<"BeamCal LayerNumber="<<i<<" , z="<<geometry.distance[i]<<" , geometry.cellSize0="<<l.getCellSize0(i-1)<<std::endl;
        }

        std::cout<<std::endl<<"INITIALIZATION\n\nLAYER [1,"<<nLayers<<"] = BEAMCAL\tLAYER "<<nLayers+1<<" = PAIR MONITOR\n\n";

   TFile ftmp("bg_aver_LDC_3.5T_14mrad_AntiDID_NominalBeamParam.root");

  if ( ftmp.IsZombie() ) {
      cerr << "Could not read data file. Exit." << endl;
      exit(1);
   }

//      cout << "TEST GetEnergyErr 1" << endl;

//      ftmp.ls();

   TTree *ttmp = (TTree*) ftmp.Get("tBcDensAverage");

//      cout << "TEST GetEnergyErr 2" << endl;


  if ( ttmp->IsZombie() ) {
      cerr << "Could not find data tree \"tBcDensAverage\". Exit." << endl;
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
        clusterPos3DHistoN = new TH3F("clusterPos3DHistoN", "Clusters 3D", 150, -150., 150., 150, -150., 150., 150, -150., 150.); //
        clusterPos3DHisto = new TH3F("clusterPos3DHisto", "Clusters 3D", 150, -150., 150., 150, -150., 150., 150, -150., 150.); //

        numberclusters = new TH1F("numberclusters", "Clusters", 10, 0., 10.); //


        energyFW = new TH1F("energyp", "Energy in FW", 50, 0., 300.); //
        energyBW = new TH1F("energyn", "Energy in BW", 50, 0., 300.); //
        clusterEnergyHisto = new TH1F("clusterEnergyHisto", "Cluster Energy", 50, 0., 300.); //
        momPxMC = new TH1F("momMC", "Momentum Px MC particles", 100, -50., 50.); //


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
        numberclusters->GetXaxis()->SetTitle("z, mm");
        numberclusters->GetYaxis()->SetTitle("events");

        coordhitszN->SetFillColor(5);   //eperlayerp->SetFillStyle(5);
        coordhitszN->GetXaxis()->SetTitle("z, mm");
        coordhitszN->GetYaxis()->SetTitle("events");


        energyFW->SetFillColor(5);      //eperlayerp->SetFillStyle(5);
        energyFW->GetXaxis()->SetTitle("E, GeV");
        energyFW->GetYaxis()->SetTitle("events");

        clusterEnergyHisto->SetFillColor(5);    //eperlayerp->SetFillStyle(5);
        clusterEnergyHisto->GetXaxis()->SetTitle("E, GeV");
        clusterEnergyHisto->GetYaxis()->SetTitle("events");


        energyBW->SetFillColor(5);      //eperlayerp->SetFillStyle(5);
        energyBW->GetXaxis()->SetTitle("E, GeV");
        energyBW->GetYaxis()->SetTitle("events");

        momPxMC->SetFillColor(5);       //eperlayerp->SetFillStyle(5);
        momPxMC->GetXaxis()->SetTitle("PxTotMC, GeV");
        momPxMC->GetYaxis()->SetTitle("events");

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


        std::cout<<"\n\n\nFINISHED INITIALIZATION\n\n\n"<<std::endl;

}



void BCalReco::processRunHeader( LCRunHeader * run ) {

        if(_nRun!=0)
//                std::cout<<"RUN "<<_nRun<<" HAS "<<_nEvt<<" TOTAL EVENTS "<<runEnergy<< "GEV TOTAL DEPOSITED SENSOR ENERGY. "<<"ZNMIN="<<znmin<<"mm   ZNMAX="<<znmax<<"mm  ZPMIN="<<zpmin<<"mm  ZPMAX="<<zpmax<<"mm"<<std::endl;
        _nRun++;
        _nEvt = 0; //number of events for each run
        runEnergy=0; //total deposited sensor energy for each run
        znmin=0;
        znmax=-DBL_MAX;
        zpmin=DBL_MAX;
        zpmax=0;


}






void BCalReco::processEvent(LCEvent * evt){

         kk += 1;

        if(kk < 4000) {
//        if(kk < 1517) {
//        if(kk < 100) {

        std::vector<ClusterImpl*> clustersNew;

        //Get elements for collection to study
        LCCollection* col = evt->getCollection( _colName ) ;


        //Create output collection
        LCCollectionVec* BCALClusters = new LCCollectionVec(LCIO::CLUSTER);


        if(col != 0) {

        int nMCP = col->getNumberOfElements();

        float pin[3]={0, 0, 0};

        mcp = 0;
        pxIPtot = 0.;

        for (int i = 0; i < nMCP; i++) {

          MCParticle* p = dynamic_cast<MCParticle*>( col->getElementAt(i) );

          pdg[mcp] = p->getPDG();

          if(abs(pdg[mcp]) == 11 || abs(pdg[mcp]) == 22) {


                //Get initial energy
                energy[mcp] = p->getEnergy();

                //Get momentum
                pin[0] = p->getMomentum() [0];
                pin[1] = p->getMomentum() [1];
                pin[2] = p->getMomentum() [2];

                pxIP[mcp] = pin[0];
                pyIP[mcp] = pin[1];
                pzIP[mcp] = pin[2];

                pxIPtot += pxIP[mcp];

          }
        } // end of MCparticle loop

   } // end col


                  momPxMC->Fill(pxIPtot,1);


        //CONSTRUCTION OF DYNAMIC ARRAY TO STORE THE INFORMATION FROM THE DETECTOR FOR THE RECONSTRUCTION CODE

        cells = new (nothrow) cell ** [nLayers+2];
        if(cells==0) {std::cerr<<"No more space for Layer allocation"; exit(1);}

        for(int i=1;i<=nLayers+1;i++){
                cells[i]=new (nothrow) cell * [nRings];
                if(cells[i]==0) {std::cerr<<"No more space for Ring allocation at layer "<<i; exit(1);}
        }

        for(int j=0;j<nRings;j++){

                int N=(int)(round(2*TMath::Pi()/geometry.segmentation[j]));

                //std::cout<<"\n\n\nRING "<<j<<"\nNumber of assigned cells: "<<N;

                for(int i=1;i<=nLayers+1;i++){

                        cells[i][j]=new (nothrow) cell [N];
                        if(cells[i][j]==0) {std::cerr<<"\nNo more space for cell allocation at layer="<<i<<" ring="<<j; exit(1);}
                }

                int k=0;
                double phi=geometry.sPhi;

//               std::cout<<"check sPhi "<< phi <<N;


                do{
                        if(k>=N) {
                                std::cout<<"\nNot enough cells allocated for ring="<<j<<" phi="<<phi<<" k="<<k<<"\n";
                                fflush(stdout);
                                exit(1);
                        }

                        if ((geometry.sPhi+geometry.dPhi-epsPhi<= phi) && (geometry.rMin+j*geometry.cellSize0 <= geometry.DArinner+epsR )) 
{
                                //std::cout<< "\nArrived to Dead Area.\tRing "<<j<<" LastConsideredPhi = "<<phi-geometry.segmentation[j]<<" NextPhi =  "<<phi<< " geometry.sPhi+geometry.dPhi-NextPhi = "<<-phi+geometry.sPhi+geometry.dPhi<<" geometry.sPhi+geometry.dPhi-LastConsideredPhi = "<<-phi+geometry.sPhi+geometry.dPhi+geometry.segmentation[j];
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

                //std::cout<<"\nOK\tRing "<<j<<" LastConsideredPhi = "<<phi-geometry.segmentation[j]<<" NextPhi = "<<phi<< " NextPhi-geometry.sPhi =  "<<phi-geometry.sPhi-2*TMath::Pi()<<" LastConsideredPhi-sPhi = "<<phi-geometry.sPhi-2*TMath::Pi()-geometry.segmentation[j];

                nbPhis[j]=k;
        }


             if(_nEvt == 1){
        std::cout<<"\n\n\nFINAL NUMBER OF PADS PER RING:\n";
        for(int j=0;j<nRings;j++)
                std::cout<<"ring "<<j<<" = "<<nbPhis[j]<<"\n";
              }

        try{//getCollection

                _nEvt ++;
                totalEvt++;

                LCCollection* inParVecSim = evt->getCollection(_HitCol); //save in inParVecSim all data related to the current event from the collection
                CellIDDecoder<SimCalorimeterHit> mydecoder(inParVecSim);

//    std::cout<<"step 1 "<<"\n";

                try{//getElementAt
                        evtEnergy=0; //total energy deposition for each event
                        nHits=inParVecSim->getNumberOfElements(); //get total number of hits from current event
                        for(int j=0; j<nHits; j++){
                                SimCalorimeterHit *hit = dynamic_cast<SimCalorimeterHit*>(inParVecSim->getElementAt(j)); // read each hit
                                side    = mydecoder(hit)[ "S-1" ] ;
                                //if(side!=0) continue;
                                //side=0 -> pos, side=1 -> neg
                                cylind  = mydecoder(hit)[ "I" ] ; //r 0-16
                                sector  = mydecoder(hit)[ "J" ] ; //phi 0-max_nb_of_pads
                                layer   = mydecoder(hit)[ "K" ]; //z with the Mokka codification: nLayers=30; layer=[1,30] for BeamCal; layer=31 for PairMonitor


                                double pos[3]={hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]}; //x,y,z

                                E=hit->getEnergy();


                                //IF ROTSATION IS NEEDED
                                /*
                                if( posTemp[2] < 0) geometry.xAngle=-TMath::Abs(geometry.xAngle);
                                else geometry.xAngle=TMath::Abs(geometry.xAngle);
                                pos[0] = TMath::Cos(geometry.xAngle/2)*posTemp[0]-TMath::Sin(geometry.xAngle/2)*posTemp[2];
                                pos[1] = posTemp[1];
                                pos[2] = TMath::Sin(geometry.xAngle/2)*posTemp[0]+TMath::Cos(geometry.xAngle/2)*posTemp[2];
                                */


                                if(cylind<0 || cylind>nRings-1 || layer<1 || layer>nLayers+1 || sector<0 || sector>nbPhis[cylind]-1) {
                                        std::cout<<" WRONG  I or J or K\n";
                                        continue;
                                }


                                if(side==1) cells[layer][cylind][sector].sEdepNeg+=E;
                                if(side==0) cells[layer][cylind][sector].sEdepPos+=E;


/*
                                if( pos[2] < 0) geometry.xAngle=-TMath::Abs(geometry.xAngle);
                                else geometry.xAngle=TMath::Abs(geometry.xAngle);
                                pos[0] = TMath::Cos(geometry.xAngle/2000)*pos[0]-TMath::Sin(geometry.xAngle/2000)*pos[2];
                                pos[1] = pos[1];
                                pos[2] = TMath::Sin(geometry.xAngle/2000)*pos[0]+TMath::Cos(geometry.xAngle/2000)*pos[2];
*/



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
                        std::cout << "Run"<<_nRun<<" Event "<<_nEvt<<": Number of hits= "<< nHits<<" Energy deposition for event = " << evtEnergy << "  [GeV]"<<std::endl;

                        runEnergy+=evtEnergy;
                        totalEnergy+=evtEnergy;



                }//try getElementAt)
                catch(Exception& e){
                        fflush(stdout);
                        std::cerr << "Exception in at least one getHit of Event "<<_nEvt <<": " << e.what() << std::endl;
                }



        }//try getCollection
        catch(Exception& e){
                fflush(stdout);
                std::cerr<<"Exception in Event "<<_nEvt<<" with get BeamCalCollection: " << e.what() << endl;
        }


//aici 1

//      cout << "TEST 1 " << "Nentr = " <<Nentr << endl;

//      cout << "TEST GetEnergyErr 5" << endl;

   TFile f("bg_aver_LDC_3.5T_14mrad_AntiDID_NominalBeamParam.root");
  if ( f.IsZombie() ) {
      cerr << "Could not read data file. Exit." << endl;
      exit(1);
   }


   t = (TTree*)f.Get("tBcDensAverage");


  if(t->IsZombie()) {
      cerr<<"bgcoeff: Tree \"tBcDensAverage\" not found"<<endl;
    }


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
//                      if(i == 0 || i == 1 || i == 2 || i == 3 || i == 4 || i == 5 || i == 6 || i == 7 || i == 8 || i == 9 || i == 10 || i == 11 )
//                            cout <<"Compare energy err " << " " << EdepPos[i] << " " << bgc[2*i].sEnDensErr*100 << endl;
                      if(EdepPos[i] > bgc[2*i].sEnDensErr*100){
                         EdepPosNoiseTmp[i] = EdepPos[i] + r.Gaus(0,bgc[2*i].sEnDensErr*100);
                         EdepPosTmp[i] = EdepPos[i];
                      }
                      else {
                         EdepPosTmp[i] = 0.;
			 EdepPosNoiseTmp[i] = 0.;
                      }
                      if(EdepNeg[i] > bgc[2*i+1].sEnDensErr*100){
                         EdepNegNoiseTmp[i] = EdepNeg[i] + r.Gaus(0,bgc[2*i+1].sEnDensErr*100);
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
////                       cout << " " << "Energy Pos = " << EdepPos[i] << " " << "RMS Map = " << bgc[2*i].sEnDensErr << endl;
////                       cout << " " << "Energy Neg = " << EdepNeg[i] << " " << "RMS Map = " << bgc[2*i+1].sEnDensErr << endl;
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



//aici 2


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

//                  cout << "Reconstructed side = " <<" " << obiect.side[0] << " " << " " << obiect.side[1] << endl;

        if(obiect.side[0] == 1 || obiect.side[1] == -1){
                  coordclusterxyP->Fill(obiect.CoordX[0],obiect.CoordY[0],1);
                  coordclusterzP->Fill(obiect.CoordZ[0],1);
                  energyFW->Fill(obiect.RecEne[0],1);
                  coordclusterxyN->Fill(obiect.CoordX[1],obiect.CoordY[1],1);
                  coordclusterzN->Fill(obiect.CoordZ[1],1);
                  energyBW->Fill(obiect.RecEne[1],1);
        }
        else if (obiect.side[0] == 0 && obiect.side[1] == 0)
                  cout << "Unreconstructed object in event = " <<" " <<_nEvt << endl;

         for (int i = 0; i < 2; ++i) {

            ClusterImpl* cluster = new ClusterImpl();

                cluster->setEnergy( obiect.RecEne[i] );

                float posit[3];
                posit[0] = obiect.CoordX[i];
                posit[1] = obiect.CoordY[i];
                posit[2] = obiect.CoordZ[i];

                cluster->setPosition( posit );


                clustersNew.push_back(cluster);

        BCALClusters->addElement(clustersNew[i]);

        }


        evt->addCollection(BCALClusters, _BCALClustersName);

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



//    cout << "End of processing event " << " " <<_nEvt << " " << "counter kk =" << " " << kk << endl;

  f.Close();

}


}





void BCalReco::check(LCEvent *evt){

        int j;
        int numElements;
        std::vector <ClusterImpl*> clustersPos, clustersNeg;
        double totalClENeg = 0, totalClEPos = 0;
        double _totalRecoEnergyNegative = 0;
        double _totalRecoEnergyPositive = 0;


//read clusters and fill basic histograms
        // position, hit position (energy weighted)
        // split clusters into positive and negative sides
        // sum up energy reconstructed in clusters
        try
        {
                LCCollection * clusters = evt->getCollection( _BCALClustersName );
                if (clusters)
                {
                        numElements = clusters->getNumberOfElements()  ;

                        //first very basic histo
                        _numClustersHisto->fill(numElements);
                         numberclusters->Fill(numElements,1);





                        for (j = 0; j < numElements ; j++)
                        {
                                ClusterImpl * cl = dynamic_cast<ClusterImpl*>( clusters->getElementAt(j) );
                                if (cl)
                                {
                                        //check side
                                        int clIsNegative = cl->getPosition()[2] < 0;

                                        //get energy and put it into a histo
                                        double clE = cl->getEnergy();

                                         clusterEnergyHisto->Fill(clE,1);

                                        //depending on side :
                                        if (clIsNegative)
                                        {
                                                //put into set for the sie
                                                clustersNeg.push_back( cl );

                                                //sum up resonstructed energy
                                                totalClENeg += clE;

                                                _totalRecoEnergyNegative += clE;

                                                //fill position histograms
                                                clusterPosHistoN->Fill(cl->getPosition()[0], cl->getPosition()[1],1);

                                                clusterPos3DHistoN->Fill(-cl->getPosition()[2], cl->getPosition()[0],
                                                        cl->getPosition()[1],1);



                                        }
                                        else //if (clIsNegative)
                                        {
                                                //same as for negative side
                                                clustersPos.push_back( cl );

                                                totalClEPos += clE;

                                                _totalRecoEnergyPositive += clE;

                                                clusterPosHisto->Fill(cl->getPosition()[0], cl->getPosition()[1],1);

                                                clusterPos3DHisto->Fill(cl->getPosition()[2], cl->getPosition()[0],
                                                        cl->getPosition()[1],1);

                                        } //else
                                } //if (cl)
                        } //for(j = 0; j < numElements ; j++)





                } //if (clusters)
        } //try

        catch(DataNotAvailableException &e)
       {
                _numClustersHisto->fill( 0 ); // no clusters found in event
        } //catch



        //celanup
        clustersPos.clear();
        clustersNeg.clear();

}








void BCalReco::end(){

//        std::cout<<"\n\n\nEND OF RUN "<<_nRun<<":\n"<<_nEvt<<" TOTAL EVENTS\nTOTAL DEPOSITED SENSOR ENERGY = "<<runEnergy<<" [GeV]\nZNMIN = "<<znmin<<" [mm]\nZNMAX = "<<znmax<<" [mm]\nZPMIN = "<<zpmin<<" [mm]\nZPMAX = "<<zpmax<<" [mm]"<<std::endl;


        std::cout<<"\n\n\nEND OF JOB:\n"<<totalEvt<<" EVENTS IN "<<_nRun << " RUNS\nTOTAL DEPOSITED SENSOR ENERGY = "<<totalEnergy<<" [GeV]"<<std::endl;

   // write the histos out

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
        coordhitszP->Write();
        coordhitszN->Write();
        energyFW->Write();
        clusterEnergyHisto->Write();
        energyBW->Write();
        momPxMC->Write();
        numberclusters->Write();
        SigPosMean->Write();
        SigNegMean->Write();
        SigGausPosMean->Write();
        SigGausNegMean->Write();

        file->Write();
        file->Close();

}



