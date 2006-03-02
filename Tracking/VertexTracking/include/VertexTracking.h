#ifndef VERTEXTRACKING_H
#define VERTEXTRACKING_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <IMPL/TrackImpl.h>
#include "ClusterExtended.h"
#include "TrackExtended.h"
#include "TrackerHitExtended.h"
#include "HelixClass.h"

using namespace lcio ;
using namespace marlin ;


/** === Vertex Tracking Processor === <br>
 * Processor performing stand-alone pattern recognition
 * in the vertex detector. <br>
 * Algorithm starts with finding of hit triplets satisfying helix hypothesis 
 * in three different layers. To accelerate procedure, the 4-pi solid angle
 * is divided in NDivisionsInTheta and NDivisionsInPhi sectors in cosQ and Phi, respectively. 
 * Triplets are looked for in 2x2 window of adjacent sectors. Once triplet is found, attempt
 * is made to associate additional hits to track. Combinatin of hits is accepted for further analysis
 * if the Chi2 of the helix fit is less than certain predefined threshold. All accepted combinations
 * are sorted in ascending order of their Chi2. First track candidate in the sorted array is automatically accepted. 
 * The hits belonging to this track are marked as used, and track candidates sharing these hits are discarded.
 * The procedure proceeds with increasing index of track candidate in the sorted array until all track 
 * candidate have been output or discarded.
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collection of digitized vertex tracker hits. <br>
 * If such a collection with the user specified name does not exist <br>
 * processor takes no action.
 * <h4>Output</h4>
 * Processor produces an output collection of the Tracks. 
 * Collection has name "VTXTracks". <br>
 * @param HitCollectionName name of input TrackerHit collection <br>
 * (default parameter value : "VTXTrackerHits") <br>
 * @param NumberOfLayers number of layers in the vertex detector <br>
 * (default parameter value : 5 ) <br>
 * @param LayerCombinations combinations of layers used to search for hit triplets <br>
 * (default parameters : 4 3 2 4 3 1 4 2 1 3 2 1; this means that triplets are looked for 
 * first in layers 4 3 2, then in 4 3 1; 4 2 1 and finally 3 2 1. NOTE THAT LAYER INDEXING STARTS FROM 0.
 * LAYER 0 is the innermost layer ) <br>
 * @param NDivisionsInPhi Number of divisions in Phi <br>
 * (default value is 60) <br>
 * @param NDivisionsInTheta Number of divisions in cosQ <br>
 * (default value is 60) <br>
 * @param Chi2WRphiTriplet Cut on chi2 in R-Phi plane to accept track with 3 hits <br>
 * (default value is 10) <br>
 * @param Chi2WZTriplet Cut on chi2 in S-Z plane to accept track with 3 hits <br>
 * (default value is 10) <br>
 * @param Chi2WRphiQuartet Cut on chi2 in R-Phi plane to accept track with 4 hits <br>
 * (default value is 10) <br>
 * @param Chi2WZQuartet Cut on chi2 in S-Z plane to accept track with 4 hits <br>
 * (default value is 10) <br>
 * @param Chi2WRphiSeptet Cut on chi2 in R-Phi plane to accept track with 5 and more hits <br>
 * (default value is 10) <br>
 * @param Chi2WZSeptet Cut on chi2 in S-Z plane to accept track with 5 and more hits <br>
 * (default value is 10) <br>
 * <br> 
 * @author A. Raspereza (MPI Munich)<br>
 */
class VertexTracking : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new VertexTracking ; }
  
  
  VertexTracking() ;
  
  /**  
   * Initialization
   */
  virtual void init() ;
  
  /** Run header processor.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Event processor.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  int _nRun ;
  int _nEvt ;

  int _nDivisionsInPhi;
  int _nDivisionsInTheta;
  int _nLayers;

  std::string _VTXHitCollection;
  std::string _relCollection;
  
  std::vector<TrackerHitExtendedVec> _sectors;
  TrackExtendedVec _tracks5Hits;
  TrackExtendedVec _tracks4Hits;
  TrackExtendedVec _tracks3Hits;

  int Initialise(LCEvent * evt);
  void ProcessOneSector(int iSectorPhi, int iSectorTheta);
  void CleanUp();
  TrackExtended * TestTriplet(TrackerHitExtended * outerHit, 
			      TrackerHitExtended * middleHit,
			      TrackerHitExtended * innerHit,
			      HelixClass & helix);
  int BuildTrack(TrackerHitExtended * outerHit, 
		 TrackerHitExtended * middleHit,
		 TrackerHitExtended * innerHit,
		 HelixClass & helix, 
		 int innerlayer,
		 int iPhi,
		 int iTheta,
		 TrackExtended * trackAR);

  void Sorting( TrackExtendedVec & trackVec);
  void CreateVTXTrack(TrackExtended * trackAR );
  void AttachRemainingHits();
  int decode(int i, int j);

  float _bField;
  float _chi2WRPhiTriplet;
  float _chi2WRPhiQuartet;
  float _chi2WRPhiSeptet;
  float _chi2WZTriplet;
  float _chi2WZQuartet;
  float _chi2WZSeptet;
  float _minDistCutAttach;
  int _minimalLayerToAttach;

  double PI,TWOPI,PIOVER2;
  double _dPhi;
  double _dTheta;
  
  std::vector<int> _Combinations;

  float _resolutionRPhi;
  float _resolutionZ;
  float _distMinCut;
  
  float _phiCutForMerging;
  float _tanlambdaCutForMerging;

  int _print;

  float _distRPhi;
  float _distZ;

  TrackExtendedVec _trackImplVec;

} ;

#endif



