#ifndef CLUSTERING_H
#define CLUSTERING_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include "CaloHitExtended.h"
#include "ClusterExtended.h"
#include "TrackExtended.h"
#include "EVENT/LCIO.h"

using namespace lcio ;
using namespace marlin ;


/** === TrackwiseClustering Processor === <br>
 *  Processor performs clustering of calorimeter hits. <br>
 *  No information on reconstructed tracks is used <br>
 *  in the clustering algorithm <br>
 *  implemented in the Processor. <br>
 *  Processor requires ECAL and <br>
 *  HCAL CalorimeterHit Collections. <br>
 *  These are specified with Processor Parameters <br>
 *  ECALCollections and HcalCollections. <br>
 *  TrackwiseClustering produces collection of Clusters <br>
 *  The name of the collection is specified with <br>
 *  Processor Parameter ClusterCollection. <br>
 *  The dependence of clustering procedure on the detector<br>
 *  geometry is minimized. 
 *  Processor needs the following calorimeter geometry <br>
 *  parameters : radius of ECAL barrel <br>
 *  n-fold symmetry of barrel <n = 8 for TESLA> <br>
 *  phi offset of the barrel stave w.r.t. <br>
 *  x-axis (0 for TESLA detector) <br>
 *  and +/- z coordinate of front face of ECAL endcaps <br> 
 *  All these parameters are passed to processor via <br>
 *  GEAR XML file <br>
 *  parameters : <br>
 *  DistanceForDirection, DistanceToTrackSeed, DistanceTrackBack, <br>
 *  StepTrackBack, ResolutionParameter, DistanceMergeForward, <br>
 *  NToDefineSP, NScanToMergeForward, TypeOfGenericDistance, <br>
 *  MinimalHitsInCluster, MaximalHitsToMerge, UseTracking, <br>
 *  DoMerging, ResolutionToMerge, WeightForResolution, WeightForDistance. <br>
 *  The meaning of these parameters as well as the overall <br>
 *  clustering procedure are described in <br> 
 *  <a href="http://www.desy.de/~rasp/Raspereza_pfa.pdf">
 *     http://www.desy.de/~rasp/Raspereza_pfa.pdf</a> <br>
 *  Please note, that the parameters steering clustering are <br>
 *  optimised for the TESLA detector on hadronic events <br>
 *  at Z pole (sqrts=91.2 GeV). <br>
 *  Printout of information on reconstructed clusters is <br>
 *  activated if processor parameter DisplayClusterInfo is set to 1. <br>
 *  @author A. Raspereza (DESY) <br>
 *  @version $Id$ <br>
 */


class TrackwiseClustering : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TrackwiseClustering ; }
  
  /** 
   * Constructor
   */  
  TrackwiseClustering() ;
  
  /** 
   * Method initialising Processor
   */
  virtual void init() ;
  
  /** Run header processor
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Event processor. Called for every event
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  /** Performs check at the end of every event   
   */
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Parameters steering clustering
   *  Adjusted for optimal clustering
   *  of Z0->hadrons events
   *  (Works both for digital RPC and Scintillator HCAL)
   */

  std::vector<float> _distanceTrackBack;
  std::vector<float> _stepTrackBack;
  std::vector<float> _resolutionParameter;
  std::vector<float> _distanceMergeForward;
  float _distanceToTrackSeed;
  float _distanceToDefineDirection;
  float _resolutionToMerge;
  int _nhit_merge_forward;
  int _use_tracks;
  int _nhit_minimal;
  int _typeOfGenericDistance;
  int _use_track;
  int _doMerging;
  int _doMergingForward;
  int _displayClusters;
  int _NDefineSP;
  int _nScanToMergeForward;

  /** Vector of strings specifying ECAL 
   *  collections of CalorimeterHits
   */
  std::vector<std::string> _ecalCollections;

  /** Vector of strings specifying HCAL 
   *  collections of CalorimeterHits
   */
  std::vector<std::string> _hcalCollections;

  /** Vector of strings specifying Reconstructed Track 
   *  collections 
   */
  std::vector<std::string> _trackCollections; 

  /** String specifying output Cluster 
   *  collection 
   */
  std::string _clusterCollection;

  ClusterExtendedVec _allSuperClusters;
  ClusterExtendedVec _allClusters;
  CaloHitExtendedVec _allHits;
  TrackExtendedVec _allTracks;

  /** Parameters specifying generic geometry of 
   *  calorimeter system
   */

  // z position of ECAL Endcap front face
  float _zofendcap;
  // radius of ECAL Barrel
  float _rofbarrel;
  // offset in Phi angle of Barrel 
  // (Phi = 0 for canonical LC detector)
  float _phiofbarrel;
  // Factor defining N_fold symmetry
  // N = 8 for canonical LC detector
  int _nsymmetry;
  // Theta of ENDCAP = atan(_rofbarrel/_zofendcap)
  float _thetaofendcap;

  float _weightForReso;
  float _weightForDist;

  int _nRun;
  int _nEvent;

  float _const_pi  ;
  float _const_2pi ;
  float _const_pi_n  ;
  float _const_2pi_n ;

  float _xmax_in_distance;
  float _xmin_in_distance;

  float _bField;


  void initialiseEvent( LCEvent * evt );
  float findResolutionParameter(CaloHitExtended * fromHit, CaloHitExtended * toHit);
  void CalculateGenericDistance(CaloHitExtended * calohit, float * dist); 
  void BubbleSort(CaloHitExtendedVec & input);  
  float DistanceBetweenPoints(float * x1, float * x2);
  void DisplayClusters(ClusterExtendedVec clusterVec);
  void GlobalSorting();
  void GlobalClustering();
  void CreateClusterCollection(LCEvent * evt, ClusterExtendedVec clusterVec);
  void mergeForward();
  void mergeLowMultiplicity();
  void calculateProperties(ClusterExtended * Cl);
  void propertiesForAll();
  void MergeTrackSegments(); 
  void CleanUp();


};


#endif
