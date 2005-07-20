#ifndef LEPTrackingProcessor_h
#define LEPTrackingProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;



/** Produces Track collection from TPC TrackerHit collections using LEP tracking algorithms. 
* The Processors BrahmsInitProcessor and GeomInitProcessor must be placed ahead in the steering file
* in order that the TPC geomerty and TKBank structre are previously instantiated. The Geometric information
* for the tpc is defined both for the C++ and F77 implementations in marlin_tpcgeom.h . The LEP algorithms 
* are taken from Brahms but does not contain the DELPHI ambiguity resolver. The track finding is based on  
*  1) SUBROUTINE CIRCLE  (N.CHERNOV, G. OSOSKOV )  
*     REFERENCE:  COMPUTER PHYSICS COMMUNICATIONS VOL 33,P329   
*  2) 3-DIMENSIONAL ITERACTION  (MARTIN POPPE)     
*     REFERENCE:  ALEPH NOTE 87-102
* The final track fitting is done using a Kalman filter. At present the both the magnetic field 
* and the material desciption of the TPC is hard coded into the Fortran code, this will be addressd once 
* the common geometry API becomes available.
* @author S. Aplin, DESY    
*/
class LEPTrackingProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new LEPTrackingProcessor ; }
  
  
  LEPTrackingProcessor() ;
  
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

  /** Input collections name.
   */
  std::string _colNameTPC ;
  std::string _colNameVTX ;

  int _nRun ;
  int _nEvt ;
} ;

#endif



