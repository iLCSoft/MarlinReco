#ifndef FPCCDDigitizer_h
#define FPCCDDigitizer_h 1

#include "FPCCDPixelHit.h"

#include "marlin/Processor.h"
#include "marlin/EventModifier.h"
#include "lcio.h"
#include <string>
#include <vector>

#include <marlin/Global.h>
#include <IMPL/LCEventImpl.h>
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
namespace IMPL {
  class SimTrackerHitImpl;
}

typedef struct {
  int layer;
  int ladder;
  int xi;
  int zeta;
}   FPCCDID_t;
  

// =================================================================
class FPCCDDigitizer : public marlin::Processor, public marlin::EventModifier {

 public:
  
  virtual Processor*  newProcessor() { return new FPCCDDigitizer ; }
  
  
  FPCCDDigitizer() ;

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
  //  virtual void processEvent( LCEvent * evt ) ;
  
  virtual void check( LCEvent * evt ) ;   
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

  // make PixelHits from SimTrackerHits.
  void makePixelHits(IMPL::SimTrackerHitImpl *simHit, FPCCDData &pxHits);

  // Initialize Geometry data
  void InitGeometry();

  void printBasicInfoToGlobal(const int nhits, IMPL::SimTrackerHitImpl *simthit);
  int getLadderID(const gear::Vector3D* pos, const int layer);
  gear::Vector3D* getLocalPos(const gear::Vector3D* pos, const int layer,const int ladder);
  void getInOutPosOnLadder(IMPL::SimTrackerHitImpl* simthit, gear::Vector3D* outpos, gear::Vector3D* inpos, gear::Vector3D* pos,gear::Vector3D* mom);
  void getInOutPosOfHelixOnLadder(IMPL::SimTrackerHitImpl* simthit,gear::Vector3D* outpos, gear::Vector3D* inpos, gear::Vector3D* pos,gear::Vector3D* mom,gear::Vector3D* BField,float charge);
  void ModifyIntoLadder(gear::Vector3D* bemodifiedpos,const int layer,gear::Vector3D* pos,gear::Vector3D* mom);
  void setLoopRange(int* looprange);
  std::vector<gear::Vector3D*> getIntersectionOfTrkAndPix(const int layer,gear::Vector3D* top,gear::Vector3D* bottom);
  std::map< std::pair< double, double>*, double> getLocalPixel(IMPL::SimTrackerHitImpl* simthit, std::vector<gear::Vector3D*> edgeofpixel);
  std::pair< double, double>* FindPixel(gear::Vector3D* f_fst, gear::Vector3D* f_nxt, int f_layer);
  void makeNewSimTHit(IMPL::SimTrackerHitImpl* simthit, gear::Vector3D* newpos, gear::Vector3D* newmom, int layer, int ladder, double newPathLength);
  int inSensitiveRegion( gear::Vector3D* pos, int layer);
  
 protected:

  std::string _colNameVTX ;
  std::string _colNameSIT ;
  std::string _colNameSET ;
  std::string _outColNameVTX ;
  std::string _outColNameSIT ;
  std::string _outColNameSET ;

  int _nRun ;
  int _nEvt ;
  int _debug;

  bool _modifySimTHit;
  bool _ladder_Number_encoded_in_cellID;
  bool _vtx_only;
  
  float _pixelSize;
  float _pixelheight;

  float _momCut;
  
  double _sigmaConst;

  float _pointResoRPhi, _pointResoZ;
  float _pointResoRPhi_SIT, _pointResoZ_SIT;
  float _pointResoRPhi_SET, _pointResoZ_SET;
  
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
