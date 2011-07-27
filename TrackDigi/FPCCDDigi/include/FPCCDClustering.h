#ifndef FPCCDClustering_h
#define FPCCDClustering_h 1

#include "marlin/Processor.h"
#include "marlin/EventModifier.h"
#include "lcio.h"
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include <string>
#include <vector>

#include <gsl/gsl_randist.h>

#include <marlin/Global.h>
#include <gear/GEAR.h>

/** ======= FPCCDClustering ========== <br>
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
typedef std::pair<unsigned int, unsigned int> FPCCDHitLoc_t;
typedef std::map<FPCCDHitLoc_t, FPCCDPixelHit*>  FPCCDLadderHit_t;
typedef std::vector<FPCCDPixelHit*> FPCCDCluster_t;
typedef std::vector<FPCCDCluster_t*> FPCCDClusterVec_t;


// =================================================================
class FPCCDClustering : public marlin::Processor, public marlin::EventModifier {
  
 public:
  
  virtual Processor*  newProcessor() { return new FPCCDClustering ; }
  
  FPCCDClustering() ;

  virtual const std::string & name() const { return Processor::name() ; }

  virtual void modifyEvent( LCEvent * evt ) ; 
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  
  virtual void check( LCEvent * evt ) ;   
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

  // make TrackerHits from VTXPixelHits
  void makeTrackerHitVec(FPCCDData &pxHits, LCCollectionVec &colvec);

  // Initialize Geometry data
  void InitGeometry();

 protected:
  // Make clusters in a ladder
  void makeClustersInALadder( int layer, FPCCDLadderHit_t &ladderHit, FPCCDClusterVec_t &cvec);

  // Make TrackerHit from clusters
  void makeTrackerHit(int layer, int ladder, FPCCDClusterVec_t &cvec, LCCollectionVec &trkHitVec);

  void EnergyDigitizer(FPCCDPixelHit* aHit);
 protected:

  std::string _colNameVTX ;
  std::string _outColNameVTX ;

  int _nRun ;
  int _nEvt ;
  int _debug;

  bool _new_tracking_system ;

  int _energyDigitization;
  int _randomNoise;
  float _pixelSize;
  float _pixelheight;
  float _pointResoRPhi;
  float _pointResoZ;

  double _electronsPerKeV;
  double _threshold;
  double _electronNoiseRate;
  int _electronsPerStep;
  int _nbitsForEdep;

  int _ranSeed;
  gsl_rng* _rng;
  
  int _nLayer;  // Number of layers
  int _maxLadder; // max number of ladders
  struct GeoData_t {
    int nladder;
    double rmin;  // distance of inner surface of sensitive region from IP
    double dphi;  // azimuthal angle step of each ladder
    double phi0;  // aximuthal angle offset
    std::vector<double> cosphi;  // cos[phi_ladder], cos_phi of each ladder
    std::vector<double> sinphi;  // sin[phi_ladder], sin_phi of each ladder
    std::vector<double> ladder_incline;
    double sthick;  // sensitive region thickness
    double sximin;  // minimum xi of sensitive region.
    double sximax;  // maximum xi of sensitive region
    double hlength; // ladder's half length in z
    int num_xi_pixel;      // Number of xi pixel in this ladder
    int num_zeta_pixel;    // Number of zeta pixel in this ladder
  };
  std::vector<GeoData_t> _geodata;

//  FPCCDCluster_t _cluster;  // Contains one ladder of clusters

} ;

#endif
