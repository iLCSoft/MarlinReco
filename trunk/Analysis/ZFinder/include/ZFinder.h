#include "marlin/Processor.h"
#include "EVENT/ReconstructedParticle.h"
#include "lcio.h"
#include <vector>
#include "IMPL/LCCollectionVec.h"

using namespace lcio ;


/** ZFinder:<br>
 * Returns the best Z->ee/Z->mm candidate in the event. Relatively
 * loose cuts are applied. For a physics analysis tighter cuts (applied
 * user analysis code) might be required. Bremstrahlung/FSR recovery is
 * on by default.
 * 
 * @author M.Thomson Cambridge
 */

class ZFinder : public marlin::Processor {
  
 public:
  
 virtual marlin::Processor*  newProcessor() { return new ZFinder ; }
  
  ZFinder() ;

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


  /** Called after data processing for clean up.
   */
  virtual void end() ;

  bool FindPFOs( LCEvent* evt );
  void FindZee( LCCollectionVec *);
  void FindZmumu( LCCollectionVec *);

private:

  std::vector<ReconstructedParticle*>_pfovec;
  int   _printing;
  std::string _zdecay;
  std::string _inputParticleCollectionName;
  std::string _outputParticleCollectionName;
  float _momentumCut;
  float _muonEcalEnergyCut;
  float _muonHcalEnergyCut;
  float _muonHcalEnergyCut1;
  float _electronEcalEnergyCut;
  float _electronHcalEnergyCut;
  float _electronEoPCutLow;
  float _electronEoPCutHigh;
  float _dmzcut;
  int   _addPhotons;
  int   _canUseClusterEnergyForElectrons;
  float _cosTrackGammaCut;

protected:

} ;

