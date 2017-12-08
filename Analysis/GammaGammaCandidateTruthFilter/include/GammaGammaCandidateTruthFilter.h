#include <marlin/Processor.h>

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/LCIterator.h"
#include "UTIL/Operators.h"
#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitSet32.h>

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>
#include "gear/BField.h"

#include <iomanip>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

/** GammaGammaCandidateTruthFilter processor
 *  Checks which GammaGammaCandidates are correct
 * author:  Graham Wilson (based on ILDPerformance/pi0)
*/

class GammaGammaCandidateTruthFilter : public Processor {
  
 public:
 

  virtual Processor*  newProcessor() { return new GammaGammaCandidateTruthFilter ; }
  
  
  GammaGammaCandidateTruthFilter() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

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

  int nEvt{};

//  float _bField ;

 private:


  std::string _trueToReco{};
  std::string _recoToTrue{};
  std::string _mcParticleCollectionName{};
  std::string _gammaGammaParticleCollectionName{};
  std::string _ggResonanceName{};

  std::vector<ReconstructedParticle*>_pfovec{};
  int   _printing{};
  std::string _outputParticleCollectionName{};
  
} ;


