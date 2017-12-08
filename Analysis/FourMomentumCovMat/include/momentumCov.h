#ifndef MOMENTUMCOV_H
#define MOMENTUMCOV_H 1

#include "marlin/Processor.h"

#include "algebraImplementation.h"


/*
   Process to obtain the covariance matrix in
   the momenta space.

   One new Collection is added to the event.
   The new collection is a copy of PandoraPFO
   with the inclusion of the covariance matrix
   in the momenta space for ReconstructedParticles
   with charge non null.

   Tue Apr  1 18:58:56 JST 2014 calancha

   ------------------------------
   Sat Apr  5 11:36:23 JST 2014
   Yesterday in the meeting i was pointed out that
   i can fill the covariance matrix directly in the
   Pandora PFOs. There is no need to create a new
   Collection. I am modifying the processor to
   follow that advice.

   ------------------------------
   Wed Jun 17 17:15:10 JST 2015
   Current version still create a new collection (hard copy
   of PandoraPFO) and implement the cov. matrix just on the
   particles belonging to this new collection.

 */

class MomentumCov : public marlin::Processor {

 public:

  virtual marlin::Processor*  newProcessor() { return new MomentumCov ; }


  MomentumCov() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( EVENT::LCRunHeader * ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( EVENT::LCEvent * ) ;


  virtual void check( EVENT::LCEvent * ) ;


  /** Called after data processing for clean up.
   */
  virtual void end() ;

  void copy_reconstructedParticle(ReconstructedParticle const *,
                                  ReconstructedParticleImpl *);


 protected:

  /** Input collection name.
   */
  std::string _colPFOs{};
  std::string _colNewPFOs{};
  int _nRun{};
  int _nEvt{};

} ;

#endif
