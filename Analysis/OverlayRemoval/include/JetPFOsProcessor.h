#ifndef JetPFOsProcessor_h
#define JetPFOsProcessor_h 1

#include "marlin/Processor.h"
#include <string>

/**  A simple processor for to save particles from jets to a PFO collection
 * 
 * @author J. Tian, Junping
 */

class JetPFOsProcessor : public marlin::Processor {
  
 public:
  
  virtual marlin::Processor*  newProcessor() { return new JetPFOsProcessor ; }
  
  
  JetPFOsProcessor() ;
  
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
  std::string _colJet ;

  /** Output collection name.
   */
  std::string _colPFOsFromJet ;

} ;

#endif



