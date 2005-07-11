#ifndef GeomInitProcessor_h
#define GeomInitProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;

/** Used in lieu of common Geometry API. Dynamicaly instansiates a TPC defined in marlin_tpcgeom.h
* which is accessed via a global pointer the_tpc.
* @author S. Aplin, DESY
*/

class GeomInitProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new GeomInitProcessor ; }
  
  
  GeomInitProcessor() ;
  
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



