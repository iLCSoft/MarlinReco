#ifndef BCalReco_h
#define BCalReco_h 

#define MAXLAYERS 35
#define MAXRINGS 20
#define MAXPHIS 200

#define MAXHITS 2000

// includes for the BCalReconstruction

#include "BCalReconstruction.h"


#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "lcio.h"
#include "TMath.h" //
#include "TH1D.h" //
#include "TH2F.h" //
#include "TH3F.h" //
#include "TFile.h" //
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGaxis.h"

#include "TROOT.h"

#include <string>

using namespace lcio ;
using namespace marlin ;

namespace AIDA
{
    class IHistogram1D;
    class IHistogram2D;
    class IHistogram3D;
}


#include "gear/GEAR.h"
#include "gear/GearParameters.h"
#include "gear/CalorimeterParameters.h"
#include "gear/LayerLayout.h"
#include "gearimpl/LayerLayoutImpl.h"


/** This processor reconstructs hits in the beamcal ....
  *  @author A.Rosca, DESY
  *  @version $Id$ 
  */


class BCalReco : public Processor {

	public:

		virtual Processor* newProcessor() { return new BCalReco ; }

		BCalReco() ;

		virtual void init() ;

		virtual void processRunHeader( LCRunHeader * run ) ;

		virtual void processEvent( LCEvent * evt ) ;

		virtual void check( LCEvent * evt ) ;

		virtual void end() ;

	protected:

		std::string _HitCol;
		std::string _MCCol;
		std::string _colName;
		std::string _BCALcolName;
		std::string _particle;
		std::string _BCALMCTruthLinkName;
		std::string _BCALClustersName;

	
		//cell tree
		typedef struct {
			double sRin,sRout,sZstart,sZend,sSphi,sDphi,sEdepNeg,sEdepPos;
			int sPos[3];
		} cell;

                //cluster_t tree
                typedef struct {
                        int nclu,in;
                        double sRin,sRout,sZstart,sZend,sSphi,sDphi,sEdepNeg,sEdepPos;
                        int sPos[3];
                } cluster_t;

                typedef struct {
                        int sPos[3];
                        double sRin;
                        double sRout;
                        double sZstart;
                        double sZend;
                        double sSphi;
                        double sDphi;
                        double sEdep;
                        double cell_area;
                        double sCellArea;
                        double sEnDens;
                        double sEnDensErr;
                } scoeff;

		cell myCell, *** cells;
                BCalReconstruction::CellType ***celule;
                BCalReconstruction::RecCorr obiect;             

                BCalReconstruction *bc_en;

		std::vector<cell> allInput;
                std::vector<cluster_t> clu;
                std::vector<cluster_t> emptyclu;
                std::vector<std::vector<cluster_t> > allClu; // vector of clusters
                
		
		//calorimeter parameters;
		int nLayers,nRings,nbPhis[MAXRINGS];
                typedef struct {
		      double xAngle,rMin,rMax,zMin,zMax,DArinner,sPhi,dPhi,pairMon,cellSize0,cellSize1,distance[MAXLAYERS+1];
		      std::vector<double> segmentation; 
                } beamcalgeom;

                 beamcalgeom geometry;

		//calculated cell parameters;
		double pos[3],E;
		int side,cylind,sector,layer;

		double coord[MAXHITS][5]; //for every hit: x,y,r,phi,E;
		int nHits; // number of hits for every event
		
		
		//calculated calorimeter parameters;
		int _nRun ; //total number of runs per job
		int totalEvt; //total number of events per job
		double totalEnergy; //total deposited sensor energy per job
		

		int _nEvt ; //number of events for each run
		double znmin,znmax,zpmin,zpmax; //z limits for each run
		double runEnergy; //total deposited sensor energy for each run			
		

		double evtEnergy; //total deposited sensor energy for each event			

                int kk ; //counter for the events

                int Nentr ; // number of entries in background map

                //per particle
		int mcp;
		enum {MCP = 100};
                int pdg[MCP];
                float energy[MCP], pxIP[MCP], pyIP[MCP], pzIP[MCP];
                float pxIPtot, pyIPtot, pzIPtot;

		TTree *t;
		TH1D *beamcalhitsn, *beamcalhitsp, *eperlayern, *eperlayerp; 
		TH2F *coordhitsxyP; 
		TH2F *coordhitsxyN; 
                TH1F *coordhitszP;
                TH1F *coordhitszN;
                TH1F *coordclusterzP;
                TH1F *coordclusterzN;
                TH1F *energyFW;
                TH1F *energyBW;
                TH1F *momPxMC;
                TH1F *numberclusters;
                TH1F *SigPosMean,*SigNegMean;
                TH1F *SigGausPosMean,*SigGausNegMean;
                TH1F *clusterEnergyHisto;
		TH2F *coordclusterxyP; 
		TH2F *coordclusterxyN; 
                TH2F *clusterPosHistoN,*clusterPosHisto;
                TH3F *clusterPos3DHistoN,*clusterPos3DHisto;
		TFile *file;
		TTree *tree;
		//TStyle *st;


  int sPos[3];
  double sRin, sRout;
  double sZstart, sZend;
  double sSphi, sDphi;
  double sEdep;
  double sCellArea,sEnDens,sEnDensErr;

  scoeff bgc[73620];



/**
* Number Of Clusters Histogram
*/
        AIDA::IHistogram1D* _numClustersHisto;

/**
* This stores the EventClusterHistoRange option
*/
        int _eventClustersHistoRange;




};


#endif
