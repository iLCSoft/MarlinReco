#ifndef PrepareSmearRECParticles_h
#define PrepareSmearRECParticles_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include "TTree.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <gsl/gsl_rng.h>

using namespace lcio ;
using namespace marlin ;


/**  PreProcessor for TauFinder to provide necessary information
 *  and create universal input for TauFinder, so that the TauFinder
 *  can run on various input information
 *
 * @author A. Muennich, CERN
 *
 */

class PrepareSmearRECParticles : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new PrepareSmearRECParticles ; }
  
  
  PrepareSmearRECParticles() ;
  
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
  std::string _colNameMC;
  std::string _colNameMCTruth;
  std::string _outcolMC;
  double _bField;
  int _nRun ;
  int _nEvt ;
  double _D0res_a, _D0res_b,_momres,_Eres;
// gsl random number generator
  gsl_rng * _random ;

  
} ;

#endif



