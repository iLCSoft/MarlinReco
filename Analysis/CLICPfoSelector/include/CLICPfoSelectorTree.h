#ifndef CLICPfoSelectorTree_h
#define CLICPfoSelectorTree_h 1

#include "marlin/Processor.h"
#include <EVENT/ReconstructedParticle.h>
#include "lcio.h"
#include <string>

#include "TH1F.h"
#include "TGraphErrors.h"
#include "TTree.h"

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
 * @version $Id: CLICPfoSelectorTree.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class CLICPfoSelectorTree : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CLICPfoSelectorTree ; }
  
  CLICPfoSelectorTree() ;
  
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

  void fillTree(LCEvent * evt, std::string collName); 

  //double profiling( std::string variable, std::string cut = "");

  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  // Input collection name
  std::string _colName_pfo{} ;

  int _nRun{} ;
  int _nEvt{} ;

  //Variables in the TTree
  TTree *m_tree = NULL;
  int type = 0;
  double p = 0.0, px = 0.0, py = 0.0, pz = 0.0, pT = 0.0;
  double costheta = 0.0;
  double energy = 0.0, mass = 0.0, charge = 0.0;
  int eventNumber = 0, nPartMC = 0, nPartPFO = 0.0;

} ;

#endif



