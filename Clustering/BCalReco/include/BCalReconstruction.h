// BCalReconstruction class
// 15/03/2011

#ifndef BCalReconstruction_h 
#define BCalReconstruction_h 1

#include "TROOT.h"
#include<vector>

using namespace std;

class BCalReconstruction
{
public:
                typedef struct {
                        int side[2];
                        double RecEne[2];
                        double ErrEne[2];
                        double CoordX[2];
                        double CoordY[2];
                        double CoordZ[2];
                        double RecRad[2];
                        double RecPhi[2];
                } RecCorr;

                typedef struct {
                        double sRin,sRout,sZstart,sZend,sSphi,sDphi,sEdepNeg,sEdepPos;
                        int sPos[3];
                } CellType;

public:
   BCalReconstruction(){
      Init();
   }
   ~BCalReconstruction(){
      Destroy();
   };
   void Print(RecCorr obj);
   RecCorr GetReconstrCoordinates(int number_layers, int number_rings, int number_pads[], CellType ***info_detector);
   // returns reconstructed clusters in either the FW, or BW BeamCal
   // RecEne, reconstructed energy, provided with an old calibration curve
   // ErrEne, background fluctuation in the pad where the cluster was found, read from the backgound map, 
   //bg_aver_LDC_3.5T_14mrad_AntiDID_NominalBeamParam.root
   // RecRad, RecPhi, ring and pad where the cluster was found
   // CoordX, CoordY, CoordZ cartesian coordinates in the global ILD reference system, crossing angle 14 mrad
   // side, 0,1,-1 -> no, FW, BW reconstruction

private:
                static const int maxlayers = 35;
                static const int maxrings = 22;
                static const int maxphis = 200;


                //calorimeter parameters;
                int the_Layers,the_Rings;
                int the_Phis[maxrings];

                //information from detector
                CellType ***BcCells;

                typedef struct {
                        int nclu,in;
                        double sRin,sRout,sZstart,sZend,sSphi,sDphi,sEdepNeg,sEdepPos;
                        int sPos[3];
                } cluster;


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


private:
   void Init();
   void Destroy();
   vector<vector<int> > SearchTowers(int the_Chains[maxrings][maxphis][maxlayers]);
   RecCorr SearchClustersFW(CellType ***info_detector);
   RecCorr SearchClustersBW(CellType ***info_detector);
   double GetEnergyCalib(double energy);
   double GetCoordRotX(double ring, double pad, float IP, float angle);
   double GetCoordY(double ring, double pad);
   double GetCoordRotZ(double ring, double pad, float IP, float angle);
   double GetEnergyErr(double ring, double pad);
   vector<vector<int> > getVector(int rows, int cols);
   void Free2DArray(int **p2DArray);
   void Free3DArray(CellType ***p3DArray);
};


#endif
