#ifndef VTXBgClusters_h
#define VTXBgClusters_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <string>
#include <vector>

#ifdef MARLIN_USE_AIDA
#include <AIDA/AIDA.h>
typedef std::vector<AIDA::IHistogram1D*> Hist1DVec;
typedef std::vector<AIDA::IHistogram2D*> Hist2DVec;
#endif

// #include <gsl/gsl_rng.h>
//  // STUFF needed for GEAR
//  #include <marlin/Global.h>
//  #include <gear/GEAR.h>
//  #include <gear/VXDParameters.h>
//  #include <gear/VXDLayerLayout.h>

using namespace lcio;
using namespace marlin;

class VXDGeometry;

/** ======= VTXBgClusters ========== <br>
 * Add Cluster parameters to VXD hits, according to projection of path length into the ladder...
 * ...
 *
 * @param RemoveDrays When this flag is set to 1 hits produced by delta-electrons are removed
 * from output collections <br>
 * (default value 0) <br>
 * @param MomentumCutForDRays The upper cut on delta-electron momentum (in MeV) <br>
 * (default value 10) <br>
 * @param Debug When this flag is set to one, debugging regime is enabled with a lot of printouts <br>
 * (default value is 0) <br>
 * <br>
 *
 *
 * @author R. de Masi, IReS
 */
class VTXBgClusters : public Processor {

public:
  VTXBgClusters(const VTXBgClusters&) = delete;
  VTXBgClusters& operator=(const VTXBgClusters&) = delete;

  virtual Processor* newProcessor() { return new VTXBgClusters; }

  VTXBgClusters();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

  //   // find phi of point in projection 2D, taken from gear::VXDParameters
  //   double getPhiPoint( gear::Vector3D p ) const ;
  //   // find phi in correct range, taken from gear::VXDParameters
  //   double correctPhiRange( double Phi ) const ;

protected:
  std::string _colNameVTX{};
  //   std::string _outColNameVTX ;

  int _nRun{};
  int _nEvt{};
  int _debug{};
  int _mod{};
  int _removeDRays{};
  //   float _pointResoRPhi,_pointResoRPhi_VTX;
  //   float _pointResoZ,_pointResoZ_VTX;
  float _momCut{};
  float _epi{};      // epitaxial thickness in unit of 50 um (standard thickness)
  float _pitch[6]{}; // pitch in um
  float _it[6]{};    // integration time in ns

  //  gsl_rng * r ;
  VXDGeometry* _vxdGeo{};

#ifdef MARLIN_USE_AIDA
  Hist1DVec _hist1DVec{};
  Hist2DVec _hist2DVec{};
#endif
};

#endif
