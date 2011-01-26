#ifndef CLICTRACKSELECTOR_H
#define CLICTRACKSELECTOR 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include "TrackHitPair.h"
#include "HelixClass.h"
#include <map>
#include <set>

using namespace lcio ;
using namespace marlin ;

/** === CLICTrackSelector Processor === <br>
 * Processor to select good tracks based on timing
 */

class CLICTrackSelector : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CLICTrackSelector ; }  
  CLICTrackSelector() ;  
  float  TimeAtEcal(Track* pTrack);
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;

 protected:

  void CleanUp(); 

  int _nRun ;
  int _nEvt ;

  std::string _inputTrackCollection;
  std::string _inputTrackMCPCollName;
  std::string _selectedTrackCollection;
  std::string _selectedTrackMCPCollection;
  float PI, PIOVER2, TWOPI;
  float _bField;
  int _debug;
  int _createMap;
  int _cutOnTPCHits;
  int _cutOnSiHits;
  LCEvent * _evt;


} ;

#endif



