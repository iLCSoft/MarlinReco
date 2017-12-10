#ifndef SelectReconstructedParticle_h
#define SelectReconstructedParticle_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;



/** Processor for marlin selecting reconstructed particle for further use. 
 */
class SelectReconstructedParticle : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new SelectReconstructedParticle;}
  
  
  SelectReconstructedParticle() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void modifyRunHeader( LCRunHeader* /*run*/ ) {}
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void modifyEvent( LCEvent * /*evt*/ ) {}
  
 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
protected:

  /** Input and output collection name.
   */

  std::string  _inputCollectionName{};
  std::string _outputCollectionName{};
  float _minimumMomentum{};

} ;

#endif
