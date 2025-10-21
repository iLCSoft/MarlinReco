#ifndef IsolatedPhotonTaggingProcessor_h
#define IsolatedPhotonTaggingProcessor_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <string>
// #include <EVENT/LCCollection.h>
// #include <EVENT/ReconstructedParticle.h>

#include "TMVA/Reader.h"

/**  processor for isolated photon tagging.
 * @version $Id: IsolatedPhotonTaggingProcessor.h,v 1.0 2020-07-05 12:50:00 Junping & Shin-ichi Exp $
 */

class IsolatedPhotonTaggingProcessor : public marlin::Processor {

public:
  virtual marlin::Processor* newProcessor() { return new IsolatedPhotonTaggingProcessor; }

  IsolatedPhotonTaggingProcessor();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

protected:
  /** Input collection name.
   */
  std::string _colPFOs{};
  std::string _colNewPFOs{};
  std::string _colPhotons{};
  std::string _isolated_photon_weights{};

  bool _is_one_isophoton{};
  float _minE{};

  float _cosConeSmall{};
  float _cosConeLarge{};

  std::vector<TMVA::Reader*> _readers{};
  Float_t _coneec{}, _coneen{}, _momentum{}, _coslarcon{}, _energyratio{};
  Float_t _ratioecal{}, _ratiototcal{}, _nsigd0{}, _nsigz0{}, _yokeenergy{}, _totalcalenergy{};

  float _mvaCut{}, _coneNeutralRatio{}, _coneChargedRatio{};
};

#endif
