#ifndef RecoMCTruthLinker_h
#define RecoMCTruthLinker_h 

#include "marlin/Processor.h"

#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include "lcio.h"
#include "EVENT/TrackerHit.h"
//#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>


#include <set>


namespace UTIL{
  class BitField64 ;
  class LCRelationNavigator ;
}

using namespace lcio ;


  
/** Optionally creates four collections of LCRelations ("recoMCTruthLink", "trackMCTruthLink", "mcTruthTrackLink", "clusterMCTruthLink",
 *  ""clusterMCTruthLink" and "calohitMCTruthLink") with weighetd relations between true particles
 *  and reconstructed particles, tracks, clusters,  and calorimeter hits, respectively.
 *  If the corresponding parameter for the output collection name is empty the collection is not created. 
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
 *  For example a weight of 0.95 for a charged particle link implies that 95 percent of the 
 *  SimTrackerHits used in the particles' track fit have been caused by the linked MCParticle.<br>
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
 *                                corresponding MCParticle - the weight is the fraction of all hits on the the track that have contributions from this MCParticle
 *  <li><b>mcTruthTrackLink</b>:  holds LCRelations  that map MCParticles to tracks - the weight is the fraction of all sim-hits from the MCParticle that contributed to this track
 *  <li><b>clusterMCTruthLink</b>: holds LCRelations  that map the  clusters to the
 *                                 corresponding MCParticle
 *  <li><b>recoMCTruthLink</b>:    holds LCRelations  that map the reconstructed particles to the
 *                                 corresponding MCParticle
 *  <li><b>calohitMCTruthLink</b>:  holds LCRelations  that map the calorimeter hits to the
 *                                  corresponding MCParticle
 *  </li>
 *  <li><b>MCParticlesSkimmed</b>:  skimmed MCParticle collection 
 *  </li>
 *  </ul>
 * 
 * @param MCParticleCollectionName      the MCParticle input collection
 * @param trackCollectionName           the ReconstructedParticles input collection
 * @param clusterCollectionName         the ReconstructedParticles input collection

 * @param UseTrackerHitRelations          use the rel collection for TrackerHits default true (false only for very old files)
 
 * @param TrackerHitsRelInputCollections the rel collection for TrackerHit collection 

 * @param SimClusterHitRelation         relation betweeen simulated and digitized cluster hits
 * @param KeepDaughtersPDG              absolute PDG code of particles where daughter are to be kept (default: gamma,pi0,K0_S)
 * @param FullRecoRelation              Select which option to use for the reconstructed link ( default: full relation)
 *
 * 
 * @param TrackMCTruthLinkName          name of output collection - default is ""
 * @param MCTruthTrackLinkName          name of output collection - default is ""
 * @param ClusterMCTruthLinkName        name of output collection - default is ""
 * @param RecoMCTruthLinkName           name of output collection - default is ""
 * @param CalohitMCTruthLinkName        name of output collection - default is ""
 * @param MCParticlesSkimmedName        skimmed MCParticle collection - default is ""
 *

 * @param UsingParticleGun              If Using Particle Gun Ignore Gen Stat - default is false
 
 * 
 *  @author M. Berggren, DESY, based on RecoMCTruthLinker v 1.0 by F. Gaede, DESY. 
 *          changelog F.Gaede, 02/2012  modified logic for optional output collections (collection name not empty) and added mcTruthTrackLink
 *  @version $Id$ 
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
  
  virtual void trackLinker(   LCEvent * evt, LCCollection* mcpCol ,  LCCollection* trackCol,  LCCollection** ttrcol ,  LCCollection** trtcol )  ;
  virtual void clusterLinker(  LCEvent * evt,   LCCollection* mcpCol ,  LCCollection* clusterCol, 
                                LCCollection** ctrcol, LCCollection** trccol,LCCollection** chittrlcol)  ;
  virtual void particleLinker(  LCCollection* mcpCol ,  LCCollection* particleCol ,  
                                LCCollection* ttrcol,  LCCollection* ctrlcol, 
                                LCCollection* trtlcol, LCCollection* trclcol,
                                LCCollection** ptrlcol , LCCollection** trplcol )  ;
  virtual void linkPrinter (  LCCollection* mcpCol, LCCollection* particleCol, LCCollection* ptrlcol, LCCollection* trplcol) ;
  virtual void check( LCEvent * evt ) ; 
  
  virtual void makeSkim(   LCCollection* mcpCol ,  LCCollection* ttrcol,  LCCollection* ctrcol ,LCCollectionVec**  skimVec) ; 
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
  
protected:
  
  virtual void mergeTrackerHitRelations(LCEvent * evt);
  virtual void mergeCaloHitRelations(LCEvent * evt);
  
  void keepMCParticle( MCParticle* mcp ) ; 

  const LCObjectVec* getSimHits( TrackerHit* trkhit, const FloatVec* weights = NULL);
  const LCObjectVec* getCaloHits( CalorimeterHit* calohit, const FloatVec* weights = NULL);
  
  UTIL::BitField64* _encoder;
  int getDetectorID(TrackerHit* hit) { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::ILDCellID0::subdet]; }


  /**  input collection names */

  std::string _mcParticleCollectionName ;
  std::string _trackCollectionName ;
  std::string _clusterCollectionName ;
  std::string _recoParticleCollectionName ;
 
  StringVec   _simTrkHitCollectionNames ;
  StringVec   _simCaloHitCollectionNames ;
 
  StringVec  _colNamesTrackerHitRelations ;
  StringVec   _caloHitRelationNames;
  
  LCCollection* _mergedTrackerHitRelCol ;
  LCRelationNavigator* _navMergedTrackerHitRel;
   LCCollection* _mergedCaloHitRelCol ;
  LCRelationNavigator* _navMergedCaloHitRel;
 
  bool _use_tracker_hit_relations;
  
  /**  output collection names */
  std::string _trackMCTruthLinkName;
  std::string _mCTruthTrackLinkName;
  std::string _clusterMCTruthLinkName;
  std::string _mCTruthClusterLinkName;
  std::string _recoMCTruthLinkName;
  std::string _mCTruthRecoLinkName;
  std::string _mcParticlesSkimmedName;
  std::string _calohitMCTruthLinkName;
  /**  output collection steering */
  bool  _FullRecoRelation;
  bool  _OutputTrackTruthRelation;
  bool  _OutputTruthTrackRelation;
  bool  _OutputClusterTruthRelation;
  bool  _OutputTruthClusterRelation;
  bool  _OutputCalohitRelation;
  bool  _OutputTruthRecoRelation;
  float _eCutMeV ;
  bool   _saveBremsstrahlungPhotons;
  float _bremsstrahlungEnergyCut;
  
  bool _using_particle_gun;
  
  IntVec _pdgVec ;

  int _nRun ;
  int _nEvt ;

  PDGSet _pdgSet ;

} ;



#endif
