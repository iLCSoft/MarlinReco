#ifndef CurlKillerProcessor_h
#define CurlKillerProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;


/** === CurlKillerProcessor === <br>
 * Kills loopers in TPC to enable efficient pattern recognition
 * <h4>Input - Prerequisites</h4>
 * Input collection of all TrackerHits in TPC
 * <h4>Output</h4>
 * Output collection of TrackerHits in TPC where all hits
 * contributing to loopes are removed from
 * @param InputCollectionName Name of the input TrackerHit collection <br>
 * (default value is TPCTrackerHits) <br>
 * @param CutCollectionName Name of the output cut away TrackerHit collection <br>
 * (default value is cutTPCTrackeHits) <br>
 * @param RemainingCollectionName Name of the output collection of remaining TrackerHits <br>
 * (default value is remainingTPCTrackerHits) <br>
 * @param BinSize Bin size in square root of pad multiples <br>
 * (default value is 2) <br>
 * @param MultiplicityCut Cut for the number of hits allowed in one bin <br>
 * (default value is 4) <br>
 * @param PadHeight TPC PadHeight <br>
 * (default value is 6.2) <br>
 * @param PadWidth TPC PadWidth <br>
 * (default value is 2.2) <br>
 * @author S. Aplin, DESY
 * @version $Id: CurlKillerProcessor.h,v 1.4 2008-07-02 09:01:23 aplin Exp $ 
 * @deprecated As of ilcsoft v01-03-06 the functionality of CurlKiller has been moved into LEPTracking 
 */

class CurlKillerProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CurlKillerProcessor ; }
  
  
  CurlKillerProcessor() ;
  
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
  std::string _inputColName ;
  std::string _cutColName ;
  std::string _remainingColName ;
  int _binSize ;
  int _multiplicityCut ;

  int _nRun ;
  int _nEvt ;
} ;

#endif



