#ifndef MCTruthJetEnergy_h
#define MCTruthJetEnergy_h 1

#include "marlin/Processor.h"
#include "lcio.h"

#ifdef MARLIN_USE_AIDA
 #include <vector>
 #include <AIDA/AIDA.h>
#endif

#include <string>

using namespace lcio ;
using namespace marlin ;


/** Calculates the true jet energy from MCParticles. Sums up the
 *  energy from all MCParticles that had the largest contribution
 *  to at least one of the ReconstructedParticles in the jet. 
 * 
 *  <h4>Input - Prerequisites</h4>
 *  jet collections
 *  MCParticle collection
 *  LCRelation(ReconstructedParticle,MCParticle)
 *
 *  <h4>Output</h4> 
 *  collection parameter 'MCTruthJetEnergies'
 * 
 * @param JetCollectionNames Name of the jet  collections
 * @param MCTruthRelationName of the mc truth relation
 * 
 * @author M.Beckmann, F.Gaede, DESY
 * @version $Id$
 */

class MCTruthJetEnergy : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new MCTruthJetEnergy ; }
  
  
  MCTruthJetEnergy() ;
  
  /** Called at the begin of the job before anything is read.
   *  Use to initialize the processor, e.g. book histograms.
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

  /** Input collection name.
   */
  StringVec _jetcolNames{}; 
  std::string  _relName{};

  int _nRun{};
  int _nEvt{};

//   float jetenergy, mcjetenergy;
//   int   njetparticles;
  
#ifdef MARLIN_USE_AIDA
  std::vector< AIDA::IHistogram1D* >  _jetEnergyHists{};    
  std::vector< AIDA::IHistogram1D* >  _jetEnergyResHists{};    
  AIDA::IHistogram2D* _jetEnergyTruthReco{};
#endif

} ;

#endif



