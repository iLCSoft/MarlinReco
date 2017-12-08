#ifndef KinkFinder_H
#define KinkFinder_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include "TrackPair.h"
#include "HelixClass.h"

using namespace lcio ;
using namespace marlin ;

/** KinkFinder Processor <br>
 *  KinkFinder processor identifes kinked tracks originating <br>
 *  from charged particle decays e.g. kaons, sigmas.  <br>
 *  KinkFinder also identifies prongs 1 -> many track matches and  <br>
 *  split tracks <br>
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
 *  The conventionial PDG codes are adopted to specify type of kink <br>
 *  Three collections (VERTEX and RECONSTRCUTED PARTICLE) can be written: Kinks, Prongs, Split Tracks <br>
 *  these are mutually exclusive.
 *  @param TrackCollection name of the input Track collection <br>
 *  (default value LDCTracks) <br>
 *  @param RecoParticleCollection name of the output collection of the Kink ReconstructedParticles <br>
 *  (default value KinkRecoParticles) <br>
 *  @param VertexCollection name of the output collection of Kink Vertices <br>
 *  (default value KinkVertices) <br>
 *  @param RecoParticleCollection name of the output collection of the Prong ReconstructedParticles <br>
 *  (default value ProngRecoParticles) <br>
 *  @param VertexCollection name of the output collection of Prong Vertices <br>
 *  (default value ProngVertices) <br>
 *  @param RecoParticleCollection name of the output collection of the Split Track ReconstructedParticles <br>
 *  (default value SplitRecoParticles) <br>
 *  @param VertexCollection name of the output collection of Split Track Vertices <br>
 *  (default value SplitVertices) <br>
 *  @param CutOnRadius minimum radius (mm) for identified kink <br>
 *  (default value 100) <br>
 *  @param DebugPrinting level of Debug printing <br>
 *  (default value 0) <br>
 *  @param HyperonDecayMassCut maximal allowed deviation in mass for Hyperon decay hypothesis <br>
 *  (default value 0.15 GeV) <br>
 *  @param HyperonTimeCut maximal cut on decay time/lifetime for Hyperon decay hypothesis <br>
 *  (default value 6 ) <br>
 *  @param SigmaDecayMassCut maximal allowed deviation in mass for Sigma decay hypothesis <br>
 *  (default value 0.15 GeV) <br>
 *  @param SigmaTimeCut maximal cut on decay time/lifetime for Sigma decay hypothesis <br>
 *  (default value 6 ) <br>
 *  @param KaonDecayMassCut maximal allowed deviation in mass for Kaon decay hypothesis <br>
 *  (default value 0.075 GeV) <br>
 *  @param PionDecayMassCut maximal allowed deviation in mass for Pion decay hypothesis <br>
 *  (default value 0.03 GeV) <br>
 *  @param KinkProjectionCutTPC maximum difference in projected track trajectories to be tagged a kink in the TPC  <br>
 *  (default value 20 mm) <br>
 *  @param TightKinkProjectionCutTPC tight cut projected track trajectories to be tagged a kink in the TPC  <br>
 *  (default value 5 mm) <br>
 *  @param VeryTightKinkProjectionCutTPC very tight cut projected track trajectories to be tagged a kink in the TPC  <br>
 *  (default value 1 mm) <br>
 *  @param KinkProjectionCutSIT maximum difference in projected track trajectories to be tagged a kink in the SIT  <br>
 *  (default value 10 mm) <br>
 *  @param LooseProjectionCutSIT loose cut in projected track trajectories to be tagged a kink in the SIT  <br>
 *  (default value 10 mm) <br>
 *  @param MaxDeltaTpcLayers maximum differenc in TPC pad rows between start-end of two associated tracks <br>
 *  (default value 10 ) <br>
 *  @param MinimumTrackHits minimum number of hits on tracks which are considered <br>
 *  (default value 5 ) <br>
 *  @param MinELambda minimum energy for Lammbda candidate <br>
 *  (default value 2 GeV) <br>
 *  @param SplitTrackMaxDeltaP maximum difference in momentum of split track segments <br>
 *  (default value 0.2 GeV) <br>
 *  @param SplitTrackMaxFracDeltaP maximum fractional difference in momentum of split track segments <br>
 *  (default value 0.02) <br>
 *  @author M. Thomson, DESY
 */

typedef struct {
  float x;
  float y;
  float z;
} vec3;

typedef struct {
  int   tracki;
  int   trackj;
  float vtx[3];
  float p[3];
  float mass;
  float distance;
  int   pdgCode;
} twoTrackIntersection_t;


/** KinkFinder Processor <br>
 *  KinkFinder processor identify kinked tracks <br>
 */
class KinkFinder : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new KinkFinder ; }
  
  
  KinkFinder() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;

 protected:

  void Sorting( TrackPairVec & trkPairVec );
  
  float kinkMass(HelixClass* parent, HelixClass* daughter, float daughterMass, float neutralMass);


  int _nRun{};
  int _nEvt{};
  
  std::string _trackColName{};

  std::string _kinkVertexColName{};
  std::string _prongVertexColName{};
  std::string _splitVertexColName{};

  std::string _kinkRecoPartColName{};
  std::string _prongRecoPartColName{};
  std::string _splitRecoPartColName{};
  
  float _rKinkCut{};
  int   _minTrackHits{};
  int   _maxDeltaTpcLayers{};
  int   _debugPrinting{};
  float _kaonDecayMassCut{};
  float _pionDecayMassCut{};
  float _sigmaDecayMassCut{};
  float _hyperonDecayMassCut{};
  float _sigmaTimeCut{};
  float _hyperonTimeCut{};
  float _tightDrCutTPC{};
  float _veryTightDrCutTPC{};
  float _drCutTPC{};
  float _drCutSIT{};
  float _looseDrCutSIT{};
  float _maxSplitTrackFracDeltaP{};
  float _maxSplitTrackDeltaP{};
  float _minELambda{};


  float _bField{};
  float _tpcInnerR{};
  float _tpcOuterR{};
  float _tpcZmax{};
  int   _tpcMaxRow{};
  int   _nLayersSIT{};
  int   _nLayersVTX{};
  std::vector<float> _rSIT{};
  std::vector<float> _rVTX{};


} ;

#endif



