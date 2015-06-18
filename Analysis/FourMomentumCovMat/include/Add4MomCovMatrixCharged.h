#ifndef ADD4MOMCOVMATRIXCHARGED_H
#define ADD4MOMCOVMATRIXCHARGED_H 1

#include "marlin/Processor.h"


/*
   Process to obtain the covariance matrix in
   the momenta space.

   It read PandoraPFOs collection and for each
   charged pfo it fills the covariance matrix.

   Need to be run after Pandora to be able to
   modify the collection.

  C. Calancha <calancha@post.kek.jp>
  2015-06-18 22:02:19 
 */

class Add4MomCovMatrixCharged : public marlin::Processor {

 public:

  virtual marlin::Processor*  newProcessor() { return new Add4MomCovMatrixCharged ; }


  Add4MomCovMatrixCharged() ;

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


 protected:

  /** Input collection name.
   */
  std::string _colPFOs ;
  int _nRun  ;
  int _nEvt  ;

} ;

#endif
