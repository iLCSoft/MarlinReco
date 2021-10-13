#ifndef anaPix_h
#define anaPix_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>

#include <marlin/Global.h>
#include <gear/GEAR.h>

/** ======= anaPix ========== <br>
 * anaPix processor is a utility for studying pixel occupancy.
 * Original code is written by Daisuke Kamai.
 * Some modifications are done by Mori.
 * @author Tatsuya Mori, Tohoku University: 2014-02-10
 */


using namespace lcio;

class FPCCDData;

// =================================================================
class anaPix : public marlin::Processor {
  
 public:
  anaPix(const anaPix&) = delete;
  anaPix& operator=(const anaPix&) = delete;
  
  virtual Processor*  newProcessor() { return new anaPix ; }
  
  
  anaPix() ;
  
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

  // make TrackerHits from VTXPixelHits
  void fillTTree(FPCCDData &pxHits);

  // Initialize Geometry data
  void InitGeometry();

 protected:

  std::string _colNameVTX{};
 
  int _nRun{};
  int _nEvt{};
  int _debug{};
  FloatVec _pixelSizeVec{};
  float _pixelSize{};
  float _pointResoRPhi{};
  float _pointResoZ{};
  std::string _rootFileName{};
  
  TFile* outroot{};
  TTree* hTreePix{};
  TTree* hTreeLocalPix{};
  
  int _nLayer{};  // Number of layers
  int _maxLadder{}; // max number of ladders
  struct GeoData_t {
    int nladder{};
    double rmin{};  // distance of inner surface of sensitive region from IP
    double dphi{};  // azimuthal angle step of each ladder
    double phi0{};  // aximuthal angle offset
    std::vector<double> cosphi{};  // cos[phi_ladder], cos_phi of each ladder
    std::vector<double> sinphi{};  // sin[phi_ladder], sin_phi of each ladder
    double sthick{};  // sensitive region thickness
    double sximin{};  // minimum xi of sensitive region.
    double sximax{};  // maximum xi of sensitive region
    double hlength{}; // ladder's half length in z
  };
  std::vector<GeoData_t> _geodata{};

} ;

#endif
