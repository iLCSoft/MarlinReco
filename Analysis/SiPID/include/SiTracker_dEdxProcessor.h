#ifndef SiTracker_dEdxProcessor_h
#define SiTracker_dEdxProcessor_h 1

#include <string>
#include <chrono>

#include "marlin/Processor.h"
#include "MarlinTrk/Factory.h"

#include "lcio.h"
#include <TTree.h>
#include <TFile.h>

#include <DDRec/SurfaceManager.h>
#include <DDRec/DetectorData.h>
#include <LayerFinder.h>

using namespace lcio ;
using namespace marlin ;



/**  SiTracker_dEdxProcessor for Marlin.
 *
 *  Calculates dEdx for planar silicon trackers and stores the information
 *  with the tracks in the lcio file.
 *
 * S. Lukic, Vinca, Belgrade
 * Dec 2016 - Jan 2018
 */


struct dEdxPoint{
public:
  dEdxPoint(const double _dE, const double _dx);
  dEdxPoint(const dEdxPoint&);

  double Get_dE() const { return dE; }
  double Get_dx() const { return dx; }
  double Get_dEdx() const { return dEdx; }

private:
  double dE;
  double dx;
  double dEdx;
} ;


class SiTracker_dEdxProcessor : public Processor {
  
 public:
  virtual Processor*  newProcessor() { return new SiTracker_dEdxProcessor ; }
  
  SiTracker_dEdxProcessor() ;
  SiTracker_dEdxProcessor(const SiTracker_dEdxProcessor &) = delete;
  
  virtual ~SiTracker_dEdxProcessor() ;

  SiTracker_dEdxProcessor & operator = (const SiTracker_dEdxProcessor &) = delete;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called at the end of every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void check( LCEvent * evt ) ; 
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  // Evaluation methods for dE/dx
  typedef std::vector<dEdxPoint> dEdxVec;
  static bool dEdxOrder(dEdxPoint p1, dEdxPoint p2) { return p1.Get_dEdx() < p2.Get_dEdx() ; }
  static double truncFractionUp;
  static double truncFractionLo;
  static double dEdxGeneralTruncMean(dEdxVec, double &dEdxError,
                              const double truncLo=0,
                              const double truncHi=0);
  static double dEdxMean(dEdxVec, double &dEdxError);
  static double dEdxMedian(dEdxVec, double &dEdxError);
  static double dEdxTruncMean(dEdxVec, double &dEdxError);
  static double dEdxHarmonic(dEdxVec, double &dEdxError);
  static double dEdxHarmonic2(dEdxVec, double &dEdxError);
  static double dEdxWgtHarmonic(dEdxVec, double &dEdxError);
  static double dEdxWgtHarmonic2(dEdxVec, double &dEdxError);


  typedef double (*evalChoice)(dEdxVec, double &dEdxError);

  protected:

  evalChoice dEdxEval{};

  /*** Steerable parameters ***/
  // Name of the track collection
  std::string m_trackCollName{};
  /* Tracker hit collection names.
   * Must be in the same order as tracker detector elements in LCDD.
   * (Check the order of tracker hit collections in the input LCIO file.)
   */
  StringVec m_trkHitCollNames{};

  // Shall we cheat the sensitive thicknesses?
  bool m_cheatSensorThicknesses{};
  // Cheat values for sensitive thicknesses. Default = -1. (no cheating)
  FloatVec m_sensThicknessCheatVals{};

  // Choice of estimator for dEdx
  std::string m_dEdxEstimator{};

  /*** Detector-related objects ***/
  const dd4hep::rec::SurfaceMap *surfMap;
  MarlinTrk::IMarlinTrkSystem *trkSystem;

  double _bField{};

  LayerFinder *layerFinder{};

  int lastRunHeaderProcessed{};

  unsigned nTimers = 8;
  std::vector<std::chrono::duration<double>> timers;
  std::chrono::high_resolution_clock::time_point lastTP;
  std::chrono::high_resolution_clock::time_point newTP;
  void addTime(int i) {
    newTP = std::chrono::high_resolution_clock::now();
    timers.at(i) += std::chrono::duration_cast<std::chrono::duration<double>>(newTP - lastTP);
    lastTP = newTP;
  }

} ;

#endif



