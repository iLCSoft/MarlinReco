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


/** Wolf Processor <br>
 *  Author : A.Raspereza <br>
 *  VERY PRELIMINARY Version <br>
 *  Processor matches Clusters to Tracks,  <br>
 *  reconstructs particles and performs very <br>
 *  simple particle ID on the basis of <br>
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
 *  are specified via Processor Parameters <br>
 *  TrackCollection and ClusterCollection. <br>
 *  Processor produces an output in the form <br>
 *  of the ReconstructedParticle Collection. <br>
 *  Name of this collection is specified <br>
 *  with Processor Parameter ParticleCollection. <br> 
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

  float _distTrackToCluster;
  float _fractionEM;
  float _bField;

  float _rPhiCut;
  float _zCut;

  ClusterExtendedVec _clusterVec;
  TrackExtendedVec   _trackVec;

  void initialiseEvent( LCEvent * evt );
  void createPartCollection( LCEvent * evt);
  void ClusterTrackMatching();
  void defineIntersection( TrackExtended * track);    
  float DistanceBetweenPoints(float * x1, float * x2 );
  void CleanUp();
} ;

#endif



