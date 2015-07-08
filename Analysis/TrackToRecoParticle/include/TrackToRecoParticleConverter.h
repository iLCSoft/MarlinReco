#ifndef TrackToRecoParticleConverter_h
#define TrackToRecoParticleConverter_h 1

#include "marlin/Processor.h"
#include <string>

using namespace marlin ;

/**
 *  Simple creation of ReconstructedParticle collection encapsulating tracks
 * 
 * @author Tomohiko Tanabe, ICEPP, The University of Tokyo
 * @author Taikan Suehara, Dept. of Physics, Kyushu University
 * @version $Id$
 */

class TrackToRecoParticleConverter : public Processor
{
  
 public:
  
  virtual Processor*  newProcessor() { return new TrackToRecoParticleConverter ; }
  
  TrackToRecoParticleConverter() ;
  virtual ~TrackToRecoParticleConverter() ;
  
  /** Called at the begin of the job before anything is read.
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
  
 private:

  /** Input collection name.
   */
	std::string _inputTrackCollectionName;
	std::string _outputPFOCollectionName;
} ;

#endif



