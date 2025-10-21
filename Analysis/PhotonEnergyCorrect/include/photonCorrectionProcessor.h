#ifndef photonCorrectionProcessor_h
#define photonCorrectionProcessor_h 1
#include "marlin/Processor.h"

#include "TFile.h"
#include "TH2F.h"
#include "photonCorrector.h"

// using namespace marlin ;

class photonCorrectionProcessor : public marlin::Processor {

public:
  virtual Processor* newProcessor() { return new photonCorrectionProcessor; }

  photonCorrectionProcessor();

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

  photonCorrectionProcessor(const photonCorrectionProcessor&) = delete;
  photonCorrectionProcessor& operator=(const photonCorrectionProcessor&) = delete;

protected:
  photonCorrector* _photonCorrector{};

  bool _modifyPFOenergies{};
  bool _modifyPFOdirections{};

  bool _validationPlots{};
  float _nominalEnergy{};

  std::string _inputCollection{};

  float _barrelendcap_limit_costh{};
  float _assumed_boxsize{};
  float _assumed_endZ{};

  std::vector<float> _energyCorr_linearise{};
  std::vector<float> _energyCorr_barrelPhi{};
  std::vector<float> _energyCorr_costheta{};
  std::vector<float> _energyCorr_endcap{};
  std::vector<float> _phiCorr_barrel{};
  std::vector<float> _thetaCorr_barrel{};
  std::vector<float> _thetaCorr_endcap{};
};

#endif
