#ifndef RecoMCTruthLinker_h
#define RecoMCTruthLinker_h 

#include "marlin/Processor.h"

#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include "lcio.h"

#include <set>


using namespace lcio ;


  
/** Creates four collections of LCRelations ("recoMCTruthLink", "trackMCTruthLink", "clusterMCTruthLink",
 *  ""clusterMCTruthLink" and "calohitMCTruthLink") with weighetd relations between true particles
 *  and reconstructed particles, tracks, clusters,  and calorimeter hits, respectively.
 *  The first is always created, the other three on request. By default, the first three are
 *  created.
 * 
 *  This relation is based on the number of hits for tracks, for hits weighted with the
 *  SimHit-energy for clusters and calorimeter hits. For tracks and clusters, the weight
 *  is the sum of hits from the considered true particle divided by the sum of all hits,
 *  for calorimter hits, it's simply the simHit energy of the hit.
 *  For the reconstructed particles, the relation could either be only to the true particle
 *  having the largest weight, or to all contributing true particles. In the former case,
 *  the pointer (and weight) will concern the true particle creating hits in the tracker,
 *  if there are any. Only for track-less seen particles would it point to the main
 *  contributor to the cluster. In the latter case, pointers are set up to all contributing
 *  true particles, and the weight is givean as (fractional contribution to track)+
 *  10000* (fractional contribution to cluster), with fractions given in permil (int).
 *  Hence: trackwgt = (int(wgt)%10000)/1000. and  clusterwgt = (int(wgt)/10000)/1000. 
 *  Which of the two is used is selected by the processor flag "FullRecoRelation" (default=true)
 *
 *  Setting the flag "FullRecoRelation" to false gives the same interpretation of "recoMCTruthLink" as 
 *  in the old version of the processor. If, in  addition "OutputTrackRelation" and  "OutputClusterRelation" 
 *  are also changed from the default true value to false, the created  output collections also agrees 
 *  with the old behaviour.
 *
 *  The calohitMCTruthLink LCRelation fixes errors in the "SimCalorimeterHitRelation", so that
 *  the originating true particle is always a particle entering the calorimeter : either it starts 
 *  outside a calorimeter, and ends inside, or back-scatters inside a calorimeter, then also ends 
 *  inside, but is not in the same cluster as it's pre-backscatter ancestors. Note that this is 
 *  a relation CalorimeterHit <-> MCParticle, not  SimCalorimeterHit <-> MCParticle as the 
 *  "SimCalorimeterHitRelation !
 *
 *  If a neutral particle with one cluster has no MC contribution assigned the MCParticle 
 *  pointing closest to the cluster is assigned and the weight is set to the negative scalar product
 *  of the MCParticle's momentum direction and the direction to the Cluster position. 
 *  (This fixes a bug in the Mokka LCal driver ov mokka-v06-06-p03).
 *
 *  <p>
 *  A skimmed MCParticle subset collection is created. It containes all particles created by the generator 
 *  program and all particles that have been reconstructed including all their parents.
 *  Additionally, the daughters of all decays in flight of particles specified in 'KeepDaughtersPDG' 
 *  (default: gamma, K0s and pi0) are kept in the skimmed list if the original particle 
 *  is in the skim ( either from the generator or from  reconstruction).
 *  
 * 
 *  <h4>Output</h4> 
 *  <ul>
 *  <li><b>trackMCTruthLink</b>:  holds LCRelations  that map the  tracks to the
 *                               corresponding MCParticle
 *  <li><b>clusterMCTruthLink</b>:  holds LCRelations  that map the  clusters to the
 *                               corresponding MCParticle
 *  <li><b>recoMCTruthLink</b>:  holds LCRelations  that map the reconstructed particles to the
 *                               corresponding MCParticle
 *  <li><b>calohitMCTruthLink</b>:  holds LCRelations  that map the calorimeter hits to the
 *                               corresponding MCParticle
 *  </li>
 *  <li><b>MCParticlesSkimmed</b>:  skimmed MCParticle collection - optional 
 *  </li>
 *  </ul>
 * 
 * @param MCParticleCollectionName      the MCParticle input collection
 * @param trackCollectionName           the ReconstructedParticles input collection
 * @param clusterCollectionName         the ReconstructedParticles input collection
 * @param SimTrackerHitRelation         relation betweeen simulated and digitized tracker hits
 * @param SimClusterHitRelation         relation betweeen simulated and digitized cluster hits
 * @param KeepDaughtersPDG              absolute PDG code of particles where daughter are to be kept (default: gamma,pi0,K0_S)
 * @param FullRecoRelation              Select which option to use for the reconstructed link ( default: full relation)
 * @param OutputTrackRelation           Output or not the track relation (default: output)
 * @param OutputClusterRelation         Output or not the cluster relation (default: output)
 * @param OutputCalohitRelation         Output or not the calohit relation (default: dont output)
 *
 * 
 * @param TrackMCTruthLinkName          name of output collection - default is "TrackMCTruthLink"
 * @param ClusterMCTruthLinkName        name of output collection - default is "ClusterMCTruthLink"
 * @param RecoMCTruthLinkName           name of output collection - default is "RecoMCTruthLink"
 * @param CalohitMCTruthLinkName        name of output collection - default is "CalohitMCTruthLink"
 * @param MCParticlesSkimmedName        skimmed MCParticle collection - default is "MCParticlesSkimmed"

 * 
 *  @author M. Berggren, DESY, based on RecoMCTruthLinker v 1.0 by F. Gaede, DESY. 
 *  @version $Id: RecoMCTruthLinker.h,v 2.0 2010/10/14 13:15:03 berggren Exp $ 
 */

class RecoMCTruthLinker : public marlin::Processor {
  

  typedef std::set< unsigned > PDGSet ;

public:
  
  virtual Processor*  newProcessor() { return new RecoMCTruthLinker ; }
  
  
  RecoMCTruthLinker() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void trackLinker(  LCCollection* mcpCol ,  LCCollection* trackCol,  LCCollection** ttrcol)  ;
  virtual void clusterLinker(  LCCollection* mcpCol ,  LCCollection* clusterCol, 
                               LCCollection* cHitRelCol , 
                               LCCollection** ctrcol, LCCollection** chittrlcol)  ;
  virtual void particleLinker(  LCCollection* particleCol ,  LCCollection* ttrcol, 
                               LCCollection* ctrlcol, LCCollection** ptrlcol )  ;
  virtual void check( LCEvent * evt ) ; 
  
  virtual void makeSkim(   LCCollection* mcpCol ,  LCCollection* ttrcol,  LCCollection* ctrcol ,LCCollectionVec**  skimVec) ; 
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
  
protected:
  
  void keepMCParticle( MCParticle* mcp ) ; 


  /**  input collection names */

  std::string _mcParticleCollectionName ;
  std::string _trackCollectionName ;
  std::string _clusterCollectionName ;
  std::string _recoParticleCollectionName ;
  std::string _trackHitRelationName;
  std::string _caloHitRelationName;

  /**  output collection names */
  std::string _trackMCTruthLinkName;
  std::string _clusterMCTruthLinkName;
  std::string _recoMCTruthLinkName;
  std::string _mcParticlesSkimmedName;
  std::string _calohitMCTruthLinkName;
  /**  output collection steering */
  bool  _FullRecoRelation;
  bool  _OutputTrackRelation;
  bool  _OutputClusterRelation;
  bool  _OutputCalohitRelation;
  float _eCutMeV ;
  
  IntVec _pdgVec ;

  int _nRun ;
  int _nEvt ;

  PDGSet _pdgSet ;

} ;



#endif
