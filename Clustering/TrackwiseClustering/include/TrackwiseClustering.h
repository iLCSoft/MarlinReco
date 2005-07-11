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


/** TrackwiseClustering Processor <br>
 *  Author : A.Raspereza   <br>
 *  VERY PRELIMINARY Version <br>
 *  Processor performs clustering of calorimeter hits. <br>
 *  No information on reconstructed tracks is used <br>
 *  in the clustering algorithm <br>
 *  implemented in the Processor. <br>
 *  Processor requires ECAL and <br>
 *  HCAL CalorimeterHit Collections. <br>
 *  These are specified with Processor Parameters <br>
 *  ECALCollections and HcalCollections. <br>
 *  TrackwiseClustering produces collection of Clusters <br>
 *  The name of collection is specified with <br>
 *  Processor Parameter ClusterCollection. <br>
 */
class TrackwiseClustering : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TrackwiseClustering ; }
  
  /** Constructor
   */  
  TrackwiseClustering() ;
  
  /** Method initialising Processor
   */
  virtual void init() ;
  
  /** Run header processor
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Event processor. Called for evergy event
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  /** Performs check at the end of evergy event   
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

  int _nRun;
  int _nEvent;

  float _const_pi  ;
  float _const_2pi ;
  float _const_pi_n  ;
  float _const_2pi_n ;

  float _xmax_in_distance;
  float _xmin_in_distance;


  void initialiseEvent( LCEvent * evt );
  float findResolutionParameter(CaloHitExtended * fromHit, CaloHitExtended * toHit);
  float CalculateGenericDistance(CaloHitExtended * calohit, int itype); 
  void BubbleSort(CaloHitExtendedVec & input);  
  float DistanceBetweenPoints(float * x1, float * x2);
  void DisplayClusters(ClusterExtendedVec clusterVec);
  void GlobalSorting();
  void GlobalClustering();
  void CreateClusterCollection(LCEvent * evt, ClusterExtendedVec clusterVec);
  void mergeForward();
  void mergeLowMultiplicity();
  void CleanUp();


};


#endif
