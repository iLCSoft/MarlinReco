#ifndef CLICDSTCHECKER_H
#define CLICDSTCHECKER 1

#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include "PfoUtilities.h"
#include "lcio.h"
#include "marlin/Processor.h"

using namespace lcio;
using namespace marlin;

/** === CLICDstChecker Processor === <br>
 * Processor to check DST Pfos
 */

class CLICDstChecker : public Processor {
public:
  virtual Processor* newProcessor() { return new CLICDstChecker; }
  CLICDstChecker();
  CLICDstChecker(const CLICDstChecker&) = delete;
  CLICDstChecker& operator=(const CLICDstChecker&) = delete;
  virtual void    init();
  virtual void processRunHeader(LCRunHeader* run);
  virtual void processEvent(LCEvent* evt);
  virtual void check(LCEvent* evt);
  virtual void end();

protected:
  void CleanUp();

  int m_nRun = -1;
  int m_nEvt = -1;

  int         m_debug = 0;
  std::string m_inputPfoCollection{};  ///< Input PFO collection name
  std::string m_inputPfoToMcRelationCollection{};
  std::string m_inputMcParticleCollection{};
  int         m_monitoring     = 0;  ///< Whether to display monitoring information
  int         m_showBackground = 0;  ///< Whether to display background information

private:
  std::vector<ReconstructedParticle*> m_pfoVector{};
  std::set<MCParticle*>               m_mcSet{};
  LCRelationNavigator*                m_pfoToMcNavigator = NULL;
};

#endif
