#ifndef FORMTRUECLUSTERSKP_H
#define FORMTRUECLUSTERSKP_H 1

#include "HelixClass.h"
#include "lcio.h"
#include "marlin/Processor.h"
#include <EVENT/CalorimeterHit.h>
#include <EVENT/MCParticle.h>
#include <string>
#include <vector>

using namespace lcio;
using namespace marlin;

/** === Cluster Cheater 5_3 === <br>
 *  This processor constructs true clusters.<br>
 *  All the hits are collected. <br>
 *  Uses gear to get inner radius and z of ecal. <br>
 *    @author P. Krstonosic (DESY)<br>
 *    @version $ld: $<br>
 */
class ClusterCheater5_3 : public Processor {

public:
  virtual Processor* newProcessor() { return new ClusterCheater5_3; }

  ClusterCheater5_3();

  /** Initialization
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
  int _nRun{};
  int _nEvt{};

  double gearRMax{};
  double zmax{};
  int _nlost{};
  std::string _trueClustCollection{};
  std::vector<std::string> _caloCollections{};
  std::string _relCollection{};
  std::string _trueClustToMCP{};
  std::string _MCcollection{};

  int _backcut{};
  int _Nmin{};
};

#endif
