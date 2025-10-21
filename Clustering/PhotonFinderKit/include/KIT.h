#ifndef KIT_h
#define KIT_h 1
#include <iostream>
#ifdef MARLIN_USE_AIDA
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogramFactory.h>
#include <marlin/AIDAProcessor.h>
#endif
#include "KITutil.h"
#include "Phys_Geom_Database.h"
#include "lcio.h"
#include "marlin/Processor.h"
#include "math.h"
#include <EVENT/CalorimeterHit.h>
#include <EVENT/Cluster.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/MCParticleImpl.h>
#include <UTIL/LCTOOLS.h>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

using namespace lcio;
using namespace marlin;
using namespace std;
using namespace IMPL;
using namespace EVENT;

/** Example processor for marlin. If compiled with MARLIN_USE_AIDA
 *  it creates a histogram (cloud) of the MCParticle energies.
 */
class KIT : public Processor {

public:
  virtual Processor* newProcessor() { return new KIT; }

  KIT();

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
  int _options{};
  int _nRun{};
  int _nEvt{};

  std::string _Ecal_col{};
  std::string _CoreCollection{};

  string _ToClean{};
  int _CleanCut{};

  int _N{};
  vector<float> _miipstep{};

  int _MinHit0{};
  int _MinHitSplit{};
  double _Rcut{};
  double _Distcut{};
  double _Coscut{};
};

#endif
