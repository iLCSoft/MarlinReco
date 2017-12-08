#include "marlin/Processor.h"
#include "EVENT/ReconstructedParticle.h"
#include "lcio.h"
#include <vector>
#include <set>
#include "IMPL/LCCollectionVec.h"

using namespace lcio ;

/** DistilledPFOCreator:<br>
 * 
 * @author Graham W. Wilson, University of Kansas
 */

class DistilledPFOCreator : public marlin::Processor {
  
 public:
  
 virtual marlin::Processor*  newProcessor() { return new DistilledPFOCreator ; }
  
  DistilledPFOCreator() ;

  /** Called at the beginning of the job before anything is read.
   *  Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 


  /** Called after data processing for clean up.
   */
  virtual void end() ;

  bool FindPFOs( LCEvent* evt );
  void CreateDistilledPFOs( LCCollectionVec *);

private:

  std::vector<ReconstructedParticle*>_pfovec{};
  std::vector<ReconstructedParticle*>_ggpfovec{};
  std::vector<ReconstructedParticle*>_mypfovec{};     // Not clear this is neeeded - but may be convenient to energy sort prior to finalization
  int   _printing{};
  std::string _inputParticleCollectionName1{};
  std::string _inputParticleCollectionName2{};
  std::string _outputParticleCollectionName{};
  static bool PfoEnergySortFunction(ReconstructedParticle* lhs,ReconstructedParticle* rhs);

protected:

} ;
