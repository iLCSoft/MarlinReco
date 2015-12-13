#ifndef AddClusterProperties_h
#define AddClusterProperties_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <IMPL/ClusterImpl.h>


using namespace std ;
using namespace lcio ;
using namespace marlin ;


/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id: AddClusterProperties.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class AddClusterProperties : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new AddClusterProperties ; }
  
  
  AddClusterProperties() ;
  
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
  
  virtual void debuging(LCCollection* clucol ,ClusterImpl* clu,double* cog,double* cov,double* eval,double* eval_err,double* evp,double* evpe,double* evc,
			int np,double sum_wgt,double sum_wgtsqr,double sum_wgt4) ; 
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection name.
   */

  std::string _clusterCollectionName;
  std::string _PFOName ;

  int _nRun ;
  int _nEvt ;
} ;

#endif



