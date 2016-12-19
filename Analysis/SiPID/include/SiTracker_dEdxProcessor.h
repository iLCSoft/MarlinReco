#ifndef SiTracker_dEdxProcessor_h
#define SiTracker_dEdxProcessor_h 1

#include <string>

#include "marlin/Processor.h"
//#include <MarlinTrk/IMarlinTrkSystem.h>
#include "MarlinTrk/Factory.h"

#include "lcio.h"
#include <TTree.h>
#include <TFile.h>

#include <DDRec/SurfaceManager.h>
#include <DDRec/DetectorData.h>
#include <LayerFinder.h>

using namespace lcio ;
using namespace marlin ;
using namespace DD4hep;


/**  SiTracker_dEdxProcessor for marlin.
 *
 *  Calculates dEdx for planar silicon trackers and stores the information
 *  with the tracks in the lcio file.
 *
 * @author S. Lukic, Vinca, Belgrade
 * December 2016
 */

class SiTracker_dEdxProcessor : public Processor {
  
 public:
  virtual Processor*  newProcessor() { return new SiTracker_dEdxProcessor ; }
  
  SiTracker_dEdxProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   * Really?
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void check( LCEvent * evt ) ; 
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
 protected:

  /*** Steerable parameters ***/
  // Input collection names
  std::string m_trackColName;
  /* Tracker hit collection names.
   * Must be in the same order as tracker detector elements in LCDD.
   * (Check the order of tracker hit collections in the input LCIO file.)
   */
  StringVec m_trkHitCollNames;
  // Bit mask which tracker detector elements to use (respecting the order in LCDD)
  int m_elementMask;

  // Shall we cheat the sensitive thicknesses?
  bool m_cheatSensorThicknesses;
  // Cheat values for sensitive thicknesses. Default = -1. (no cheating)
  FloatVec m_sensThicknessCheatVals;

  const DDRec::SurfaceMap *surfMap;
  MarlinTrk::IMarlinTrkSystem *trkSystem;

  double _bField;

  LayerFinder *collFinder;

  int lastRunHeaderProcessed;
} ;

#endif



