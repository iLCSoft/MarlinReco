#ifndef RecoMCTruthLinker_h
#define RecoMCTruthLinker_h 

#include "marlin/Processor.h"

#include <EVENT/MCParticle.h>
#include "lcio.h"

#include <set>


using namespace lcio ;


  
/** Creates a collection of LCRelations ("RecoMCTruthLink")  with a weighted relation between the 
 *  ReconstructedParticles and their corresponding MCParticles. 
 *  This relation is based on the number of hits - energy weighted for neutral particles - 
 *  that have been used in creating the Track for charged and the Cluster for neutral 
 *  ReconstructedParticles.
 *  Every ReconstructedParticle is related to the MCParticle with the largest contribution
 *  stored in the weight of the Relation.
 *  For example a weight of 0.95 for a charged particle link implies that 95 percent of the 
 *  SimTrackerHits used in the particles' track fit have been caused by the linked MCParticle.<br>
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
 *  <li><b>RecoMCTruthLink</b>:  holds LCRelations  that map the  ReconstructedParticles to the
 *                               corresponding MCParticle
 *  </li>
 *  <li><b>MCParticlesSkimmed</b>:  skimmed MCParticle collection - optional 
 *  </li>
 *  </ul>
 * 
 * @param MCParticleCollectionName      the MCParticle input collection
 * @param RecoParticleCollectionName    the ReconstructedParticles input collection
 * @param SimTrackerHitRelation         relation betweeen simulated and digitized tracker hits
 * @param SimClusterHitRelation         relation betweeen simulated and digitized cluster hits
 * @param KeepDaughtersPDG              absolute PDG code of particles where daughter are to be kept
 *
 * 
 * @param RecoMCTruthLinkName        name of output collection - default is "RecoMCTruthLink"
 * @param MCParticlesSkimmedName     skimmed MCParticle collection - optional 
 * 
 *  @author F. Gaede, DESY
 *  @version $Id: RecoMCTruthLinker.h,v 1.2 2008-04-15 13:01:20 gaede Exp $ 
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
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
  
protected:
  
  void keepMCParticle( MCParticle* mcp ) ; 


  /**  input collection names */

  std::string _mcParticleCollectionName ;
  std::string _recoParticleCollectionName ;
  std::string _trackHitRelationName;
  std::string _caloHitRelationName;

  /**  ouput collection name */
  std::string _recoMCTruthLinkName;
  std::string _mcParticlesSkimmedName;
  
  IntVec _pdgVec ;

  int _nRun ;
  int _nEvt ;

  PDGSet _pdgSet ;

} ;



#endif
