#ifndef FPCCDDigitizer_h
#define FPCCDDigitizer_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gearimpl/Vector3D.h>
#include <gear/BField.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>

/** ======= FPCCDDigitizer ========== <br>
 * Produces VTXPixelHits collection from VXDCollection collections. <br> 
 * VXDCollection is a SimTrackerHit object and VTXPixelHits collection
 * is a LCGenericObject object, where FPCCD pixel hit information is stored.
 * Format of LCGenericObject is described in FPCCDData.cc
 * 
 * PixelHit ID is given by ( layer, ladder, xi, zeta ). xi and zeta are 
 * pixel address in a coordinate system local to a ladder. zeta is along 
 * Z axis, xi is in the ladder plain and eta is perpendicular to the plain 
 * of ladder. Further description will be found in the source code of 
 * FPCCDID::encodeFPCCDID(...)
 *
 * Conversion of SimTrackerHit to PixelHits is performed in 
 * FPCCDDigitizer::makePixelHits(...).  
 *
 * Parameters of this process
 *   Debug : default(0). : if 1, print debug information
 *   FPCCD_PixelSize : default(0.005) : FPCCD pixel size, which is used 
 *                     digitization
 *  
 * <br>
 * @author Akiya Miyamoto, KEK: 2010-04-19
 */

class FPCCDData;

namespace EVENT {
  class SimTrackerHit;
}

typedef struct {
  int layer;
  int ladder;
  int xi;
  int zeta;
}   FPCCDID_t;
  

// =================================================================
class FPCCDDigitizer : public marlin::Processor {

 public:
  
  virtual Processor*  newProcessor() { return new FPCCDDigitizer ; }
  
  
  FPCCDDigitizer() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ;   
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

  // make PixelHits from SimTrackerHits.
  void makePixelHits(const EVENT::SimTrackerHit *simHit, FPCCDData &pxHits);

  // Initialize Geometry data
  void InitGeometry();

  // Obtain FPCCDID data from layer number and hit position
  FPCCDID_t encodeFPCCDID(const int layer, const gear::Vector3D pos);

  void printBasicInfoToGlobal(const int nhits, const EVENT::SimTrackerHit *simthit);
  int getLadderID(const EVENT::SimTrackerHit* simthit);
  double getToLocalRad(const EVENT::SimTrackerHit* simthit,const int ladderID);
  gear::Vector3D* getLocalPos(const EVENT::SimTrackerHit* simthit,const double tolocalrad);
  gear::Vector3D* getMomAtLocalPos(const EVENT::SimTrackerHit* simthit,const double tolocalrad);
  gear::Vector3D* getTopBottomPosOnLadder(gear::Vector3D* pos,gear::Vector3D* mom,const char* topbottom);
  gear::Vector3D* getInOutPosOfHelixOnLadder(const EVENT::SimTrackerHit* simthit,gear::Vector3D* pos,gear::Vector3D* mom,gear::Vector3D* BField,float charge,const char* inout);
  void ModifyIntoLadder(gear::Vector3D* bemodifiedpos,const EVENT::SimTrackerHit* simthit,gear::Vector3D* pos,gear::Vector3D* mom);
  void setLoopRange(int* looprange);
  std::vector<gear::Vector3D*> getIntersectionOfTrkAndPix(const EVENT::SimTrackerHit* simthit,gear::Vector3D* pos,gear::Vector3D* mom,gear::Vector3D* top,gear::Vector3D* bottom);
  std::map<gear::Vector3D*, double> getLocalPixel(const EVENT::SimTrackerHit* simthit, std::vector<gear::Vector3D*> edgeofpixel);
  gear::Vector3D* FindPixel(gear::Vector3D* f_fst, gear::Vector3D* f_nxt, int f_layer, char* updown);


 protected:

  std::string _colNameVTX ;
  std::string _outColNameVTX ;

  int _nRun ;
  int _nEvt ;
  int _debug;
  
  float _momCut;
  float _pixelSize;

  float _pixelsizex;
  float _pixelsizey;
  float _pixelheight;

  double _electronsPerKeV;
  double _threshold;

  bool _isSignal;
  
  int _sloopx,_eloopx,_sloopy,_eloopy;
  
// Variables to store geometry information 
  int _nLayer;  // Number of layers
  int _maxLadder; // max no. of ladder in each layer

  struct GeoData_t {
    int nladder;
    double rmin;  // distance of inner surface of sensitive region from IP
    double dphi;  // azimuthal angle step of each ladder
    double phi0;  // aximuthal angle offset
    std::vector<double> cosphi;  // cos[phi_ladder], cos_phi of each ladder
    std::vector<double> sinphi;  // sin[phi_ladder], sin_phi of each ladder
    double sthick;  // sensitive region thickness
    double sximin;  // minimum xi of sensitive region.
    double sximax;  // maximum xi of sensitive region
    double hlength; // ladder's half length in z
  };
  std::vector<GeoData_t> _geodata;

} ;

#endif
