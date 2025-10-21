#ifndef MyProcessor_h
#define MyProcessor_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <string>

using namespace lcio;
using namespace marlin;

/** Processor that calculates sphericity,aplanarity, C and D event parametres
 *   for detail explanation look
 *   <li> <a href="www.desy.de/~aplin/KP_spher.pdf">documentation</a></li>
 * @author P.K , DESY
 */
class Sphere : public Processor {

public:
  virtual Processor* newProcessor() { return new Sphere; }

  Sphere();

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
  std::string _colName{};
  /** Name of the parameter to store egenvalues of sphericity tensor
   */
  std::string _dumpobjectname{};
  float _r{};
  int _nRun{};
  int _nEvt{};
};

#endif
