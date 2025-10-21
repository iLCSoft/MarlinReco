#ifndef IsoLepTrainingProcessor_h
#define IsoLepTrainingProcessor_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <string>
// #include <iostream>
#include <sstream>

/**  Example processor for marlin.
 *
 *  If compiled with MARLIN_USE_AIDA
 *  it creates a histogram (cloud) of the MCParticle energies.
 *
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4>
 *  A histogram.
 *
 * @param CollectionName Name of the MCParticle collection
 *
 * @author F. Gaede, DESY
 * @version $Id: IsoLepTrainingProcessor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $
 */

class IsoLepTrainingProcessor : public marlin::Processor {

public:
  virtual Processor* newProcessor() { return new IsoLepTrainingProcessor; }

  IsoLepTrainingProcessor();

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
  std::string _colMCP{};
  std::string _colMCTL{};
  std::string _colPFOs{};
  std::string _colPVtx{};

  bool _mcdebug{};

  bool _is_lep_tune{};
  bool _is_for_sig{};

  int _nRun{};
  int _nEvt{};

  int _iso_lep_type{};
};

#endif
