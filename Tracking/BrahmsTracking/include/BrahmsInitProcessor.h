#ifndef BrahmsInitProcessor_h
#define BrahmsInitProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;


/** Dynamicly instaniates bank structures which replace the ZEBRA banks needed by the LEP 
* tracking routines used in FortranProcessor. The Banks are defined in tk*bank.h and accessed
* via a global pointer.
* @author S. Aplin, DESY
*/
class BrahmsInitProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new BrahmsInitProcessor ; }
  
  
  BrahmsInitProcessor() ;
  
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
  std::string _colName ;

  int _nRun ;
  int _nEvt ;
} ;

#endif



