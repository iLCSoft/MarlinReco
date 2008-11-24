#ifndef TTBarExample_h
#define TTBarExample_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;


/**  An example processor for a kinematic fit
 *   
 *   ... testing a ttbar -> 6jets hypothesis
 *   with energy and momentum conservation
 *   and W mass constraints
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs two ReconstructedParticle collections, one with 4 light jets, one with 2 b-tagged jets
 *
 *  <h4>Output</h4> 
 *  Some histogram.
 * 
 * @param CollectionName Name of the ReconstructedParticle collection for light jets
 * @param CollectionName Name of the ReconstructedParticle collection for b-jets
 * 
 * @author J. List, DESY
 * @version $Id: TTBarExample.h,v 1.2 2008-11-24 11:01:01 beckmann Exp $ 
 */

class TTBarExample : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TTBarExample ; }
  
  
  TTBarExample() ;
  
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
  std::string _lightjetcolName ;
  std::string _bjetcolName ;
  /** Input parameter: center of mass energy.
   */
  float _ecm ;

  int _nRun ;
  int _nEvt ;
} ;

#endif



