#include "marlin/Processor.h"
#include "EVENT/ReconstructedParticle.h"
#include "lcio.h"
#include <vector>
#include <set>
#include "IMPL/LCCollectionVec.h"

using namespace lcio ;

/** GammaGammaSolutionFinder:<br>
 * 
 * @author Graham W. Wilson, University of Kansas
 */

class GammaGammaSolutionFinder : public marlin::Processor {
  
 public:
  
 virtual marlin::Processor*  newProcessor() { return new GammaGammaSolutionFinder ; }
  
  GammaGammaSolutionFinder() ;

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
  void FindGammaGammaSolutions( LCCollectionVec *);

private:

  std::vector<ReconstructedParticle*>_pfovec;
  int   _printing;
  std::string _inputParticleCollectionName1;
  std::string _inputParticleCollectionName2;
  std::string _inputParticleCollectionName3;
  std::string _outputParticleCollectionName;
  static bool PfoProbabilitySortFunction(ReconstructedParticle* lhs,ReconstructedParticle* rhs);

protected:

} ;
