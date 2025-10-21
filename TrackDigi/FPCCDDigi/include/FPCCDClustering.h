#ifndef FPCCDClustering_h
#define FPCCDClustering_h 1

#include "lcio.h"
#include "marlin/EventModifier.h"
#include "marlin/Processor.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <gsl/gsl_randist.h>

#include <gear/GEAR.h>
#include <marlin/Global.h>

/** ======= FPCCDClustering ========== <br>
 *
 * Some modifications are done by mori.
 * 1. The utility for making LCRelation between TrackerHit and SimTrackerHit is set.
 * 2. The utility for studying the position resolution of valiable clusters is set.
 * 3. The system of reading position resolution map for valiable clusters is set.
 * 4. Pair-BG cluster rejection algorithms are set.
 * @author Tatsuya Mori, Tohoku University: 2014-02-10
 *
 *
 *
 *
 *
 * Produces VXDTrackerHits collection from VTXPixelHits collections for FPCCD. <br>
 *
 * Parameters of this process
 *   Debug : default(0). : if 1, print debug information
 *   FPCCD_PixelSize : default(0.005) : FPCCD pixel size, which is used
 *                     digitization
 *   PointResolutionRPhi: default(0.001440) : resolution value assigned to
 *                     coveriance matric of TrackerHits.
 *   PointResolutionZ:    default(0.001440) : resolution value assigned to
 *                     coveriance matric of TrackerHits.
 *
 * <br>
 * @author Akiya Miyamoto, KEK: 2010-04-19
 *
 */

using namespace lcio;

class FPCCDData;
class FPCCDPixelHit;
class FPCCDDigitizer;
class TTree;
class TFile;

typedef std::pair<unsigned int, unsigned int> FPCCDHitLoc_t;
typedef std::map<FPCCDHitLoc_t, FPCCDPixelHit*> FPCCDLadderHit_t;
typedef std::vector<FPCCDPixelHit*> FPCCDCluster_t;
typedef std::vector<FPCCDCluster_t*> FPCCDClusterVec_t;

// =================================================================
class FPCCDClustering : public marlin::Processor, public marlin::EventModifier {

public:
  FPCCDClustering(const FPCCDClustering&) = delete;
  FPCCDClustering& operator=(const FPCCDClustering&) = delete;

  virtual Processor* newProcessor() { return new FPCCDClustering; }

  FPCCDClustering();

  virtual const std::string& name() const { return Processor::name(); }

  virtual void modifyEvent(LCEvent* evt);

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Called for every event - the working horse.
   */

  virtual void check(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

  // make TrackerHits from VTXPixelHits
  void makeTrackerHitVec(FPCCDData* pxHits, LCCollection* STHcol, LCCollectionVec* relCol, LCCollectionVec* trkHitvec);

  // Initialize Geometry data
  void InitGeometry();

protected:
  // Make clusters in a ladder
  void makeClustersInALadder(int layer, FPCCDLadderHit_t& ladderHit, FPCCDClusterVec_t& cvec);

  // Make TrackerHit from clusters
  void makeTrackerHit(LCCollection* STHcol, int layer, int ladder, FPCCDClusterVec_t& cvec,
                      std::multimap<std::pair<int, int>, SimTrackerHit*> relMap, LCCollectionVec* relCol,
                      LCCollectionVec* trkHitVec);

  void EnergyDigitizer(FPCCDPixelHit* aHit);

protected:
  std::string _colNameSTH{};
  std::string _colNameVTX{};
  std::string _outColNameVTX{};
  std::string _outRelColNameVTX{};

  int _nRun{};
  int _nEvt{};
  int _debug{};

  bool _new_tracking_system{};
  bool _remove_pixelhits_collection{};
  bool _makeRelation{}; // 2012_12_20 added to escape BG hit linking error.
  int _energyDigitization{};
  int _randomNoise{};
  float _pixelSize{};
  FloatVec _pixelSizeVec{};
  float _pixelheight{};
  float _pointResoRPhi{};
  float _pointResoZ{};

  double _electronsPerKeV{};
  double _threshold{};
  double _electronNoiseRate{};
  int _electronsPerStep{};
  int _nbitsForEdep{};

  int _ranSeed{};
  gsl_rng* _rng{};

  int _nLayer{};    // Number of layers
  int _maxLadder{}; // max number of ladders

  struct GeoData_t {
    int nladder{};
    double rmin{};                // distance of inner surface of sensitive region from IP
    double dphi{};                // azimuthal angle step of each ladder
    double phi0{};                // aximuthal angle offset
    std::vector<double> cosphi{}; // cos[phi_ladder], cos_phi of each ladder
    std::vector<double> sinphi{}; // sin[phi_ladder], sin_phi of each ladder
    std::vector<double> ladder_incline{};
    double sthick{};      // sensitive region thickness
    double sximin{};      // minimum xi of sensitive region.
    double sximax{};      // maximum xi of sensitive region
    double hlength{};     // ladder's half length in z
    int num_xi_pixel{};   // Number of xi pixel in this ladder
    int num_zeta_pixel{}; // Number of zeta pixel in this ladder
  };
  std::vector<GeoData_t> _geodata{};

  //  FPCCDCluster_t _cluster;  // Contains one ladder of clusters

  //*****************For Root Writing*******************//
  /*
  _tree->Branch("nlinks",&_link.nlink,"nlinks/I");
  325   _tree->Branch("weight",_link.weight,"weight[nlinks]/F");
  326   _tree->Branch("simthits_x",_simthits.x,"simthits_x[nlinks]/D");
  327   _tree->Branch("simthits_y",_simthits.y,"simthits_y[nlinks]/D");
  328   _tree->Branch("simthits_z",_simthits.z,"simthits_z[nlinks]/D");
  329   _tree->Branch("simthits_pdg",_simthits.pdg,"simthits_pdg[nlinks]/I");
  330
  331   _tree->Branch("trkhits_x",_trkhits.x,"trkhits_x[nlinks]/D");
  332   _tree->Branch("trkhits_y",_trkhits.y,"trkhits_y[nlinks]/D");
  333   _tree->Branch("trkhits_z",_trkhits.z,"trkhits_z[nlinks]/D");
  334   _tree->Branch("trkhits_CWidth_RPhi",_trkhits.CWidth_RPhi,"trkhits_CWidth_RPhi[nlinks]/D");
  335   _tree->Branch("trkhits_CWidth_Z",_trkhits.CWidth_Z,"trkhits_CWidth_Z[nlinks]/D");
  */

  TTree* _tree{};
  TFile* _rootf{};
#define MAX_LINK 5000
  struct {
    unsigned int nlink;
    float weight[MAX_LINK];
  } _link{};

  struct {
    double x[MAX_LINK];
    double y[MAX_LINK];
    double z[MAX_LINK];
    double R[MAX_LINK];
    double vx[MAX_LINK];
    double vy[MAX_LINK];
    double vz[MAX_LINK];
    double mcp_px[MAX_LINK];
    double mcp_py[MAX_LINK];
    double mcp_pz[MAX_LINK];
    double mcp_Pt[MAX_LINK];
    double mcp_energy[MAX_LINK];
    float mcp_time[MAX_LINK];
    int mcp_isCreatedInSimulation[MAX_LINK];
    int mcp_isBackscatter[MAX_LINK];
    int mcp_vertexIsNotEndpointOfParent[MAX_LINK];
    int mcp_isDecayedInTracker[MAX_LINK];
    int mcp_isDecayedInCalorimeter[MAX_LINK];
    int mcp_hasLeftDetector[MAX_LINK];
    int mcp_isStopped[MAX_LINK];
    float mcp_d0[MAX_LINK];
    float mcp_phi0[MAX_LINK];
    float mcp_omega[MAX_LINK];
    float mcp_z0[MAX_LINK];
    float mcp_tanL[MAX_LINK];
    double xi[MAX_LINK];
    double zeta[MAX_LINK];
    float edep[MAX_LINK];

    float px[MAX_LINK];
    float py[MAX_LINK];
    float pz[MAX_LINK];
    double pAbs[MAX_LINK];
    int pdg[MAX_LINK];
    double mass[MAX_LINK];
    double theta[MAX_LINK];
    double phi[MAX_LINK];
    double area_theta[MAX_LINK];
    double area_phi[MAX_LINK];
    int layer[MAX_LINK];
    int ladder[MAX_LINK];
  } _simthits{};

  struct {
    double x[MAX_LINK];
    double y[MAX_LINK];
    double z[MAX_LINK];
    double tposX[MAX_LINK];
    double tposY[MAX_LINK];
    double tposZ[MAX_LINK];
    double cx[MAX_LINK];
    double cy[MAX_LINK];
    double cz[MAX_LINK];
    double dot[MAX_LINK];
    double R[MAX_LINK];
    double area_theta[MAX_LINK];
    double area_phi[MAX_LINK];
    unsigned int CWidth_RPhi[MAX_LINK];
    unsigned int CWidth_Z[MAX_LINK];
    unsigned int tilt[MAX_LINK];
    unsigned int nPix[MAX_LINK];
    int layer[MAX_LINK];
    int ladder[MAX_LINK];
    double xi[MAX_LINK];
    double zeta[MAX_LINK];
    float edep[MAX_LINK];
  } _trkhits{};

  struct {
    double RPhi[MAX_LINK];
    double Z[MAX_LINK];
  } _diff{};

  std::string _rootFileName{};
  std::string _treeName{};

  bool _positionReso_ReadingFile_ON{};
  std::string _positionReso_ReadingFile{};
  typedef std::map<short int, float> ResoMap;
  ResoMap _resolutionMapRPhi{};
  ResoMap _resolutionMapZ{};
  std::ifstream _fin{};
  void calcTrackParameterOfMCP(MCParticle* pmcp, double* par);

  struct firstCut {
    bool isActive{};
    IntVec RPhiWidth{};
    IntVec ZWidth{};
    IntVec nPix{};
  } _firstCut{};

  struct Mori2ndCut {
    bool isActive{};
    FloatVec zpar{};
  } _m2Cut{};

  struct Kamai2ndCut {
    bool isActive{};
    FloatVec bpar{};
    IntVec minZWidth{};
  } _k2Cut{};
};

#endif
