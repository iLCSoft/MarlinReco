#ifndef V0Finder_H
#define V0Finder_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include "TrackPair.h"

using namespace lcio ;
using namespace marlin ;


/** V0Finder Processor <br>
 *  V0Finder processor identify neutral vertices originating <br>
 *  from photon conversions and decays of K0S and Lamda0 <br>
 * <h4>Input collections and prerequisites</h4> 
 *  Processor requires collection of tracks. The name of the collection <br>
 *  is specified by the processor parameter "TrackCollection". <br> 
 *  If no collection with the specified name exist in event <br>
 *  processor takes no action <br>
 *  <h4>Output</h4>
 *  Processor produces LCIO collections of the reconstructed particles, <br> 
 *  and vertices, containing information on the reconstructed neutral vertices <br>
 *  Position of the vertex is accessed through the LCIO object VERTEX. <br>
 *  Four-vector of the vertex is stored in the object RECONSTRUCTEDPARTICLE <br>
 *  Type of the VERTEX is accessed through the method ReconstructedParticle::getType() <br>
 *  The convential codes are adopted to specify type of the neutral vertices <br>
 *  (22 - photon conversion, 310 - K0S, 3122 - Lamda0) <br>
 *  @param TrackCollection name of the input Track collection <br>
 *  (default value LDCTracks) <br>
 *  @param RecoParticleCollection name of the output collection of the ReconstructedParticles <br>
 *  (default value V0RecoParticles) <br>
 *  @param VertexCollection name of the output collection of Vertices <br>
 *  (default value V0Vertices) <br>
 *  @param CutOnRadius cut on the vertex radius vector. If vertex radius sqrt(x**2+y**2) <br>
 *  is less than cut value, the candidate vertex is discarded. In this case two tracks are assumed to <br>
 *  originate from the primary interaction vertex. <br>
 *  (default value 6.0 mm) <br>
 *  @param CutOnTrkDistance cut on the distance between tracks constituting neutral vertex <br>
 *  If the 3D distance between two tracks with opposite charge is less that cut value, <br>
 *  the two tracks are regarded as constituting neutral vertex <br>
 *  (default value 1.5 mm) <br>
 *  @param MassRangeGamma maximal allowed deviation in mass for gamma hypothesis <br>
 *  (default value 0.01 GeV) <br>
 *  @param MassRangeK0S maximal allowed deviation in mass for K0S hypothesis <br>
 *  (default value 0.01 GeV) <br>
 *  @param MassRangeL0 maximal allowed deviation in mass for L0 hypothesis <br>
 *  (default value 0.008 GeV) <br>
 *  @author A.Raspereza, DESY
 *  @version $Id: V0Finder.h,v 1.2 2008-07-10 11:36:24 rasp Exp $ 
 */
class V0Finder : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new V0Finder ; }
  
  
  V0Finder() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;

 protected:

  void Sorting( TrackPairVec & trkPairVec );
  float Rmin( Track* track );
  
  int _nRun ;
  int _nEvt ;
  
  std::string _trackColName;
  std::string _vertexColName;
  std::string _recoPartColName;
  
  float _rVertCut;
  float _dVertCut;

  float _bField;

  float MASSProton;
  float MASSPion;
  float MASSLambda0;
  float MASSK0S;
  float MASSGamma;

  float _deltaMassK0S;
  float _deltaMassL0;
  float _deltaMassGamma;

  float _rxyCutGamma;
  float _rxyCutK0S;
  float _rxyCutLambda;

  float _minTrackHitRatio;

} ;

#endif



