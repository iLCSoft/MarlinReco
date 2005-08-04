#ifndef RECOPARTICLESAR_H
#define RECOPARTICLESAR_H 1

#include "marlin/Processor.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCIO.h"
#include "lcio.h"
#include "ClusterExtended.h"
#include "TrackExtended.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;


/** === Wolf Processor === <br>
 *  Processor matches Clusters to Tracks,  <br>
 *  reconstructs particles and performs <br>
 *  particle identification on the basis of <br>
 *  the fraction of the energy deposited in ECAL. <br>
 *  If the fraction of the energy deposited in ECAL <br>
 *  greater than predefined threshold, particle  <br>
 *  is regarded to be EM one (photon or electron). <br>
 *  The corresponding threshold is specified with <br>
 *  processor parameter FractionEM. <br>  
 *  Particle kinematics is defined at point <br>
 *  of closest approach to IP. <br>
 *  The following convention is assumed <br>
 *  for particle types :               <br>
 *       1 = electron/positron         <br>
 *       2 = charged hadron or muon    <br>
 *       3 = photon                    <br>
 *       4 = neutral hadron            <br>
 *                                     <br>
 *  Processor requires Reconstructed Track Collection <br>
 *  and Reconstructed Cluster Collection. <br>
 *  The names of corresponding collections <br> 
 *  are specified through Processor Parameters <br>
 *  TrackCollection and ClusterCollection. <br>
 *  Processor produces an output in the form <br>
 *  of the ReconstructedParticle Collection. <br>
 *  Name of this collection is specified <br>
 *  with Processor Parameter ParticleCollection. <br> 
 *  Processor needs the following calorimeter geometry <br>
 *  parameters : radius of ECAL barrel <br>
 *  n-fold symmetry of barrel <n = 8 for TESLA> <br>
 *  phi offset of the barrel stave w.r.t. <br>
 *  x-axis (0 for TESLA detector) <br>
 *  and +/- z coordinate of front face of ECAL endcaps <br> 
 *  All these parameters are passed to processor via <br>
 *  GEAR XML file <br>
 *  Magnetic field (in units of Tesla) must be specified <br>
 *  with processor parameter BField <br>
 *  <br>
 *  <br>
 *  Processor has the following parameters : <br>
 *    - TrackCollection : name of the Track collection <br>
 *    - ClusterCollection : name of the Cluster collection <br>
 *    - ParticleCollection : name of ReconstructedParticle collection <br>
 *    - DistanceTrackToCluster : cut on distance between track <br>
 *      intersection point with inner boundary of calorimeter <br>
 *      and nearest calorimeter hit in a given cluster. <br>
 *      If distance between track intersection point and <br>
 *      nearest hit in cluster is less then this cut value, <br>
 *      then cluster is assigned to track <br>
 *    - FractionEM : If ratio of the ecal energy to the total energy <br>
 *      in cluster is greater that FractionEM, then particle <br>
 *      is considered to be electron/positron or photon. <br>
 *    - NativeTrackFitter : If 1, then internal track fitter is invoked <br>
 *      to define track parameters. If 0, then default track parameters <br>
 *      are used to calculate charge particle momentum and track intersecton <br>
 *      point with inner boundary of calorimeters. <br> 
 *    - MergeClusters : If 1, then attempt is made to merge additional <br>
 *      clusters to track if energy of already associated cluster is <br>
 *      still much less than track momentum. The condition to merge <br>
 *      additional clusters to track is quantified <br>
 *      with the following formula : <br>
 *      Ptrk - Ecluster > 3.0*HcalResolution*sqrt(Ptrk) <br>
 *      where Ptrk - track momentum, Ecluster - energy of associated cluster <br> 
 *      and HcalResolution is the energy resolution for hadrons <br>
 *      (this is assumed to be 0.5 for TESLA calorimeters) <br>
 *    - LowerMomentum : The above procedure of merging is applied if momentum <br>
 *      of track is greater than LowerMomentum. <br>
 *    - DistMergeCut : This parameter defines cut on distance of cluster <br>
 *      centre-of-gravity to helix for the cluster-to-track merging procedure <br>
 *      applied when Ptrk - Ecluster > 3.0*HcalResolution*sqrt(Ptrk) <br>
 *    A more detailed description of the algorithm can be found at <br>
 *    <a href="http://www.desy.de/~rasp/PFlowInMarlin.ps.gz">
 *     http://www.desy.de/~rasp/PFlowInMarlin.ps.gz</a> <br>
 *    @author A. Raspereza (DESY)<br>
 *    @version $ld: $<br>
 */
class Wolf : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new Wolf ; }
  
  
  Wolf() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
    
  virtual void check( LCEvent * evt ) ; 
    
  virtual void end() ;
  
  
 protected:

  int _nRun ;
  int _nEvt ;
  
  // Name of Track collection
  std::string _trackCollection;
  // Name of Cluster collection
  std::string _clusterCollection;
  // Name of Reconstructed Particle Collection
  std::string _particleCollection;

  // Parameters specifying detector geometry
  // (look at TrackwiseClustering Processor)
  float _zofendcap;
  float _rofbarrel;
  float _phiofbarrel;
  int _nSymmetry;
  int _trackFitter;

  float _distTrackToCluster;
  float _fractionEM;
  float _bField;

  float _rPhiCut;
  float _zCut;
  float _lowerMom;
  float _hcalReso;
  float _distMergeCut;
  int _mergeClusters;

  ClusterExtendedVec _clusterVec;
  TrackExtendedVec   _trackVec;

  void initialiseEvent( LCEvent * evt );
  void createPartCollection( LCEvent * evt);
  void ClusterTrackMatching();
  void defineIntersection( TrackExtended * track);    
  float DistanceBetweenPoints(float * x1, float * x2 );
  void MergeClustersToTracks();
  float angleVectors(float * vec1, float * vec2);
  void CleanUp();
} ;

#endif



