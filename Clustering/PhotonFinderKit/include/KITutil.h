#ifndef KITUTIL_h
#define KITUTIL_h 1

#include "lcio.h"
#include "EVENT/LCIO.h"
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <vector>
#include <list>
#include <stack>
#include <string>
#include <UTIL/CellIDDecoder.h>
#include "Phys_Geom_Database.h" 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_pow_int.h>
#include <ANN/ANN.h>

using namespace lcio;
using namespace std;

#define MAXARRAYSIZE 100000

class Tmpcl2;
/**
 *    Basic hit class for reconstruction, contains the calorimeter hit plus     
 *    additional parameters
 *    @authors P.Krstonosic (DESY)
 */
class Superhit2
{
 public:

  /**
   *  Constructor , needs absolute MIP energy (calibrated ) , and pointer to 
   *  calorimeter hit
   */
  Superhit2(float E,CalorimeterHit* chitin);  
  /**
   *    Destructor
   */
  ~Superhit2();
  /**
   *   Returns position of the hit i=1 transformed position, i=2 true position
   */
  const float* getPosition(int i);
  /**
   *   Pointer to the LCIO calorimeter hit
   */
  CalorimeterHit* chit;
  
  bool connect;
  /**
   *   Transformed position of the hit
   */
  float point[3];  // transform position
  /**
   *   Real coordinate of the hit, copy of calorimeter hit data for direct access
   */
  float pointt[3]; // true position
  /**
   *   Energy of hit i terms of MIP
   */ 
  float mip;
  /**
   *   MIP value [GeV]
   */
  float mipE;
  /**
   *   Number of neighbours
   */
  int top;         // topological parameter
  /**
   *    Is hit assigned to a cluster or not 
   */
  bool is_assigned;  
  /**
   *    Vector of pointers to the neighbouring hits of Superhit2 type
   */
  vector <Superhit2*> neighbours;
  /**
   *   Pointer to the cluster to wich hit belong 
   */
  Tmpcl2 * cl;
  /**
   *  stove as coded by Mokka 
   */
  int S;
  /**
   *  layer as coded by Mokka 
   */
  int K;
  /**
   *  module as coded by Mokka
   */
  int M;
  /**
   *  0 by constructor  1 for ecal 2 for hcal 
   */ 
  int tip; 
 
};  
/**
 *    Basic cluster class for reconstruction
 *    @authors P.Krstonosic (DESY)
 */
class Tmpcl2{

 public : 
  /**
    * Constructor.
    */
   Tmpcl2();
  /**
    * Destructor. 
    */
  ~Tmpcl2();
  /**
   *  Calculates center of inerta for the objects in hits vector
   */
  void calcCenter();  
  /**
   *  returnes double[3] for the center as calculated with calcCenter method
   */ 
  double* getCenter();
  /**
   *  returnes energy of cluster in GeV.
   */
  double getEnergy();
  /**
   *  calculates eigenvalues and eigenvectors of inertia tensor
   */
  void findInertia();

  // vector<Superhit2*> konektori; 
  /**
   * hit in the cluster
   */
   vector<Superhit2*> hits;
   /**
    * pointers to clusters that are contained in this cluster
    */
  vector<Tmpcl2*> daughters;
  /**
   * pointers to clusters that contain this cluster
   */
  vector<Tmpcl2*> parents;
  /**
   * energy  in GeV (sum over hits in cluster) 
   */
  double energy;  
  /**
   *  position of cluster center (x,y,z)
   */
  double center[3];
  /**
   *  principal axes of inertia tensor , calculated in calcInertia
   */
  double direction[3];
  /**
   * normalized eigenvalues of inertia tensor
   */
  double inteigen[3]; 
  /**
   * eigenvectors of inertia tensor
   */       
  double inteigenvec[9];     

  /**
   *  Internal type of the cluster
   */
  int type ;      
};

class Photon2
{

 public:  
  /**
   * Constructor.
   */
  Photon2(double Ein,double* pravac, double* pocetak);
   /**
    * Destructor.
    */
  ~Photon2();
  
  void Prob(CalorimeterHit* ch,double cut,double* out);
  
  // data- stvari koje se racunaju jednom i gotovo 
  double z1;
  double z2;
  double k1;
  double k2;
  double k3;
  double k4;
  double p1;
  double p2;
  double p3;
  double y;
  double eprime;
  double Z;
  double x0eff;
  double sampling;
  double Eceff;
  double Rm;
  double Fs;
  double Thom;
  double Tsam;

  double alfahom;
  double alfasam; 
  double betasam;
  //  
  double Ee;
  double dir[3];
  double start[3];
};
typedef vector<Superhit2*>  Shitvec2;
typedef vector<Tmpcl2*>     Tmpclvec2;

/**
 * container for holding the numbers needed for energy estimation
 */
typedef struct{
  /**
   * level of cluster
   */
  int level;
  /**
   * nominal i.e. input energy in GeV
   */
  double Enom;
  /**
   *   Eestimate = a+b*Ecore ,  a parameter of this funcition
   */
  double a;
  /**
   *   Eestimate = a+b*Ecore ,  b parameter of this funcition
   */
  double b;
  /**
   *   lower validity range of energy estimation funciton in GeV
   */
  double Emin;
 /**
   *   upper validity range of energy estimation funciton in GeV
   */
  double Emax;
}  CoreCalib2;

/**
 * container for storing the EM shower core candidates
 */
typedef struct{
  /**
   * pointer to the cluster
   */
  Tmpcl2* cl;
  /**
   *  level of the cluster
   */
  int level;
  /**
   *  center of the candidate
   */
  double X[3];
  /**
   *  flag to activate deactivate 
   */
  bool active;
}PROTSEED2;

/**
 *  container for keeping the parameters of the core fineder together
 */
typedef struct{
  /**
   * fluctuation suppresion cut 
   */
  double Rcut;
  /**
   *  distance cut for core merging 
   */
  double Distcut;
  /**
   *  angular cut for core merging (value of the cosine is stored not the angle!)
   */
  double Coscut;
  /**
   * minimal number of hits needed for 0-th level cluster to be accepted as a core candidate
   */
  unsigned int MinHit0;
  /**
   * minimal number of hits in i-th level cluster to be accepted for splitting  of the core
   */
  unsigned int MinHitSplit;
}CoreCut2;


/**
 *  Creation of superhits, input is ECAL collection ,it's decoded and pointer to resulting container for
 *  superhits
 */
void CreateAllShits2(LCCollection* colt,CellIDDecoder<CalorimeterHit>& id,vector<Superhit2*>* calo);

/**
 *  Global precalculation function , iput is vector of superhits, ECAL decoder, number of hits, and 
 *  number of neighbors  cut for hit separation
 */
void TotalPrecalc2(vector<Superhit2*>* calo,CellIDDecoder<CalorimeterHit>& id,
		  unsigned int nelem, int Ncut);
/**
 *  Precalculation function used internaly by TotalPrecalc 
 */
void Precalc2(vector< Superhit2* >& shvec,double r, double z, double cell, double dist,bool tip,int stove,int module,CellIDDecoder<CalorimeterHit>& id);

/**
 *  Basic function for transformation of hit coordinates
 */
void GridTransform2( CalorimeterHit* clh,float& radius, float& halfz, float& cellsize,float*X,
		    bool tip,int stove,int module,CellIDDecoder<CalorimeterHit>& id);

/**
 *  Global EM core finding function , iput is vector of superhits, array of Tmpcl vectors bbb for internal
 *  computation, vector of EM core candidates prs2, and parameters of the algorithm
 */
void FindCores2(Shitvec2* secal1, Tmpclvec2* bbb , vector <PROTSEED2> * prs2,
		unsigned int N, vector<float> miipstep, CoreCut2 Ccut);
/**
 *  NN clustering 
 */
void cluster5( vector<Superhit2*>* shv, vector<Tmpcl2*>* clv);
/**
 * returns energy estimate for a given core , input  int level, double core energy in GeV, and 
 * calibration data table
 */
double  giveMeEEstimate2(int nivo,double Ecore, vector<CoreCalib2> cc);
/**
 * Example function for creation of the energy estimaiton table 
 */
void CreateCalibrationLDC00(vector<CoreCalib2>* cc);


void LineCaloIntersectD2( double* X1, double* dir,double&d,double&zmax, double*X);
void LineCaloIntersect2(double* X1, double* X2,double&d,double&zmax,  double*X);
double LinePointDistance2( double* X1, double* X2, double* X0);
void PointOnLine3(const double* X1,const double* X2,const float* X0,double* Xline);
void PointOnLine22(const double* Xstart,const double* dir,const float* X0,double* Xline);
void ModuleNormal2(double* X1,double& zmax, double* X0);
void ClusterInCluster2(Tmpcl2* cl, vector<Tmpcl2*>& clv);
double D_cl_cl2(Tmpcl2* cl1,Tmpcl2* cl2) ;
inline double Dot2(double* X1,double* X2);
void ClusterInCluster2(Tmpcl2* cl, vector<Tmpcl2*>& clv,vector<Tmpcl2*>& clout);

#endif



