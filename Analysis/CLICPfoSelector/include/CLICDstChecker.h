#ifndef CLICDSTCHECKER_H
#define CLICDSTCHECKER 1

#include "marlin/Processor.h"
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include "lcio.h"
#include <string>
#include <map>
#include <set>
#include <algorithm>

#define FORMATTED_OUTPUT_TRACK_CLUSTER(N1, E1,E2,E3,N2,N3)             		                \
    std::cout  <<                                                                               \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E3        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N2        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N3

#define FORMATTED_OUTPUT_TRACK(N1, E1,E2,E3,N2,N3)                                              \
    std::cout  <<                                                                               \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E3        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N2        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N3 

#define FORMATTED_OUTPUT_CLUSTER(N1, E1,E2,E3,N2,N3)                                            \
    std::cout  <<                                                                               \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E3        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N2        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N3 

#define FORMATTED_OUTPUT_MC(N1,E1)                                                              \
    std::cout  <<                                                                               \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1


using namespace lcio ;
using namespace marlin ;

/** === CLICDstChecker Processor === <br>
 * Processor to check DST Pfos
 */

class CLICDstChecker : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CLICDstChecker ; }  
  CLICDstChecker() ;
  CLICDstChecker(const CLICDstChecker&) = delete;
  CLICDstChecker& operator=(const CLICDstChecker&) = delete;
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;

 protected:

  void CleanUp(); 

  int m_nRun=-1;
  int m_nEvt=-1;

  int   m_debug=0;
  std::string m_inputPfoCollection{};                           ///< Input PFO collection name
  std::string m_inputPfoToMcRelationCollection{};
  std::string m_inputMcParticleCollection{};
  int         m_monitoring=0;                                   ///< Whether to display monitoring information
  int         m_showBackground=0;                               ///< Whether to display background information

 private:
  std::vector<ReconstructedParticle*> m_pfoVector{};
  std::set<MCParticle*> m_mcSet{};
  static bool PfoSortFunction(ReconstructedParticle* lhs,ReconstructedParticle* rhs); 
  LCRelationNavigator* m_pfoToMcNavigator=NULL;

} ;

#endif



