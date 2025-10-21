#include "EVENT/ReconstructedParticle.h"
#include "IMPL/LCCollectionVec.h"
#include "lcio.h"
#include "marlin/Processor.h"
#include <vector>

using namespace lcio;

/** GammaGammaCandidateFinder:<br>
 *
 * (modelled after ZFinder processor)
 *
 * @author Graham W. Wilson, University of Kansas
 */

class GammaGammaCandidateFinder : public marlin::Processor {

public:
  virtual marlin::Processor* newProcessor() { return new GammaGammaCandidateFinder; }

  GammaGammaCandidateFinder();

  /** Called at the beginning of the job before anything is read.
   *  Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

  bool FindPFOs(LCEvent* evt);
  void FindGammaGammaCandidates(LCCollectionVec*);

private:
  std::vector<ReconstructedParticle*> _pfovec{};
  int _printing{};
  std::string _inputParticleCollectionName{};
  std::string _outputParticleCollectionName{};
  std::string _ggResonanceName{};
  float _gammaMomentumCut{};
  float _ggResonanceMass{};
  float _dmggcut{};
  double _fitProbabilityCut{};
  int _ifitter{};

protected:
};
