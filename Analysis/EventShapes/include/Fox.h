#ifndef Fox_h
#define Fox_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;



/** Processor calculates Fox-Wolfram moments.
 *  input parameters are  name of the ReconstructedParticle collection
 *  and arbitrary number of moments entered just as integers 
 *  for example CalculateFoxWolframMoments 1 3 6 , moment of 0-th
 *  order is calculated by default (for later normalisation)
 *  so you don't need to enter this one.
 *  Output is stored as a parameter of te ReconstructedParticle collection
 *  with the name  FoxWolfram_moment(n) where n is actual number. 
 *  @author P.Krstonosic , DESY
 */
class Fox: public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new Fox ; }
  
  
   Fox() ;
  
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
  /** Vector of moments orders to be calculated
   */
  std::vector<int> _momentsToCalculate;
  int _nRun ;
  int _nEvt ;
} ;

#endif



