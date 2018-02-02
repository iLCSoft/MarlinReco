#ifndef MyProcessor_h
#define MyProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

using namespace lcio ;
using namespace marlin ;


class TH1 ;

/** Computes various estimators for the time of flight from CalorimeterHits.
 * Creates ROOT histograms for these parameters.
 * 
 * @author N. Weinhold, DESY internship 2017
 */

class TOFPlots : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TOFPlots ; }
  
  TOFPlots(const TOFPlots&) = delete;
  TOFPlots& operator=(const TOFPlots&) = delete;

  TOFPlots() ;
  
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
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection name.
   */
  std::string _colNameMCP{};
  std::string _colNamePFO{};

  int _nRun{};
  int _nEvt{};

  gsl_rng* _rng = nullptr ;

  std::vector<TH1*> _h{};

} ;

#endif



