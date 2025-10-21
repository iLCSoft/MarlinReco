// //////////////////////////////////////////////////////////////////////// //
//                                                                          //
// SiStripClus - Marlin Processor - provides clustering in Si strip sensors //
//                                                                          //
// //////////////////////////////////////////////////////////////////////// //

#ifndef SISTRIPCLUS_H
#define SISTRIPCLUS_H 1

// Define ROOT output if needed
// #define ROOT_OUTPUT

#include <queue>
#include <vector>

// Include CLHEP header files
#include <CLHEP/Vector/ThreeVector.h>

// Include Clus header files
#include "SiStripDigi.h"
#include "SiStripGeom.h"
#include "StripCluster.h"

// Include LCIO header files
#include <EVENT/LCCollection.h>
#include <IMPL/TrackerPulseImpl.h>
#include <UTIL/LCRelationNavigator.h>
#include <lcio.h>

// Include Marlin
#include <marlin/Global.h>
#include <marlin/Processor.h>

// ROOT classes
#ifdef ROOT_OUTPUT
#include <TFile.h>
#include <TTree.h>
#endif

// Namespaces
using namespace lcio;
using namespace marlin;

namespace sistrip {

// Typedefs
typedef const std::vector<std::string> ConstStringVec;
typedef std::vector<StripCluster*> ClsVec;
typedef std::vector<double*> DoubleVec;
typedef std::vector<LCCollection*> LCCollectionVec;
typedef std::vector<std::string> StringVec;
typedef std::queue<std::string> StringQue;

//! Marlin processor intended for cluster finding - uses digitized data worked out with SiStripDigi
//!
//! @author Z. Drasal, Charles University Prague
//!

class SiStripClus : public Processor {
public:
  //! Method that returns a new instance of this processor
  virtual Processor* newProcessor() { return new SiStripClus(); }

  //! Constructor - set processor description and register processor parameters
  SiStripClus();

  //! Method called at the beginning of data processing - used for initialization
  virtual void init();

  //! Method called for each run - used for run header processing
  virtual void processRunHeader(LCRunHeader* run);

  //! Method called for each event - used for event data processing
  virtual void processEvent(LCEvent* event);

  //! Method called after each event - used for data checking
  virtual void check(LCEvent* event);

  //! Method called after all data processing
  virtual void end();

protected:
  // MAIN CLUSTER METHOD

  //! Method searching for clusters - first, strips above _SNseed threshold, so-called seed
  //! strips, are defined. Then the strips adjacent to the seeds and above _SNadjacent
  //! threshold are extracted. Finally, clusters taken as Gaussian are calculated (if total
  //! charge is above _SNtotal threshold) and their mean positions and sigmas are saved in either
  //! R-Phi or Z. Finally, they are mixed into 3D cluster. (input parameter: sensor map
  //! of strips with total integrated charge, output parameter: sensor map of clusters
  //! found by this algorithm)
  ClsVec findClus(SensorStripMap& sensorMap);

  // OTHER METHODS
  //! Method calculating hits from given clusters
  void calcHits(ClsVec& clsVec, IMPL::LCCollectionVec* colOfTrkHits);

  //! Method calculating hit resolution, i.e. covariance matrix
  float* calcResolution(const int& layerID, const double& hitTheta, const double& posZ);

  //! Method to update and store the Sensor strip map
  void updateMap(TrackerPulseImpl* pulse, SensorStripMap& sensorMap);
  //! Method to release memory of the SensorStripMap
  void releaseMap(SensorStripMap& sensorMap);

  // PRINT METHODS
  //! Method printing processor parameters
  void printProcessorParams() const;

  //! Method printing hit info
  void printHitInfo(const StripCluster* pCluster) const;

  // VARIABLES

  // Collection names
  std::string _inColName;          //!< LCIO input collection name
  std::string _outColName;         //!< LCIO output collection name
  std::string _relColNamePlsToSim; //!< LCIO input relation collection name  - TrackerPulse <-> SimTrkHit

  // ClusterFinder parameters - set by users
  float _CMSnoise;   //!< Common mode subtracted noise, set in ENC
  float _SNseed;     //!< Signal to noise ratio cut for seed strips
  float _SNadjacent; //!< Signal to noise ratio cut for adjacent strips
  float _SNtotal;    //!< Signal to noise ratio cut for total cluster

  // Compensate Lorentz shift - set by users
  // float _TanOfAvgELorentzShift;          //!< Tangent of electrons' average Lorentz shift
  float _TanOfAvgHLorentzShift; //!< Tangent of holes' average Lorentz shift

  // Geometry parameters
  bool _floatStripsRPhi;  //!< Is every even strip floating in R-Phi?
  bool _floatStripsZ;     //!< Is every even strip floating in Z?
  SiStripGeom* _geometry; //!< All geometry information from Gear xml file

  // SVD resolution
  std::vector<float> _resSVDFirstInRPhi; //!< Mean strip resolution in RPhi - 1st layer; in [mm]
  std::vector<float> _resSVDOtherInRPhi; //!< Mean strip resolution in RPhi - other layers; in [mm]
  std::vector<float> _resSVDFirstInZ;    //!< Mean strip resolution in Z    - 1st layer; in [mm]
  std::vector<float> _resSVDOtherInZ;    //!< Mean strip resolution in Z    - other layers; in [mm]

  // Relation navigators
  LCRelationNavigator* _navigatorPls;

  // Pitch
  double _pitchFront; //!< Pitch in the middle of the front sensor
  double _pitchRear;  //!< Pitch in the middle of the rear sensor

  std::string _subdetector; //!< Name of the subdetector to be clusterize

  // Root output
#ifdef ROOT_OUTPUT

  TFile* _rootFile;
  TTree* _rootTree;

  int _rootEvtNum; //!< Event number to which the hit belongs

  int _rootLayerID;
  int _rootLadderID;
  int _rootSensorID;

  // In R-Phi
  double _rootSimRPhi;     //!< Simulated hit position in R-Phi (SimTrackerHit with heighest weight taken)
  double _rootRecRPhi;     //!< Reconstructed hit position in R-Phi
  double _rootClsSizeRPhi; //!< Cluster size in R-Phi
  int _rootMCPDGRPhi;      //!< PDG of particle which created current hit in R-Phi (MC Particle with highest weight)

  // Residuals
  double _rootResRPhi;   //!< Residual RPhi direction
  double _rootResR;      //!< Residual R direction
  double _rootResModule; //!< Residual Module (distance between
                         //!  simHit and recHit in the sensor plane)

  // In Z
  double _rootSimZ;     //!< Simulated hit position in Z (SimTrackerHit with heighest weight taken)
  double _rootRecZ;     //!< Reconstructed hit position in Z
  double _rootClsSizeZ; //!< Cluster size in Z
  int _rootMCPDGZ;      //!< PDG of particle which created current hit in Z (MC Particle with highest weight)
#endif

private:
  int _nRun;   //!< Run number
  int _nEvent; //!< Event number

  //
  //! Calculated and stored the hits adjacents
  template <class It>
  StripChargeMap& storeHitsAdjacents(StripChargeMap& clsStripsIn, It schmap, const It& endIt, StripChargeMap& cM);
};

} // namespace sistrip

#endif // SISTRIPCLUS_H
