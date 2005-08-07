#ifndef FORMTRUETRACKSAR_H
#define FORMTRUETRACKSAR_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;


/** === Track Cheater Processor === <br>
 *  Constructs true tracks. True track is considered to comprise <br> 
 *  all the TrackerHits attributable to the same particle. <br>
 *  An user has to provide the names of TrackerHit collections <br>
 *  via processor parameter TrackerHitCollections. <br>
 *  An output collection of tracks is specified <br>
 *  with processor parameter TrueTrackCollection. <br>
 *  MCParticle - True Track relations are stored in the <br> 
 *  collection with name TrueTrackToMCP. True track is stored <br>
 *  in the corresponding collection if the number of <br>
 *  TrackerHits is greater than certain value. The latter is <br>
 *  specified with processor parameter MinimalHits. <br>
 *  Magnetic field is defined with processor parameter BField. <br>
 *  There is an option to perform a fit of TrackerHits assigned to <br>
 *  track. This is done if processor parameter FitTrueTrack is set <br>
 *  to 1. In this case parameters extracted from the fit are used <br>
 *  to define track parameters, namely D0, Z0, Phi, TanLambda and Omega. <br>
 *  If FitTrueTrack is set to 0 then true Monte Carlo information, namely <br>
 *  MCParticle momentum and vertex, is used to define <br>
 *  track parameters. Processor parameter HitToHelixInFit defines <br>
 *  the maximal allowed distance from TrackerHit to the helix derived <br>
 *  from true particle momentum and vertex. If distance between hit and <br>
 *  true helicoidal track is greater than HitToHelixInFit, hit is not assigned <br>
 *  to track. If the chi2 of the fit is greater that certain predefined <br>
 *  value (processor parameter Chi2Cut), then true MC information is used <br>
 *  to define track parameters. Only tracks with energy greater than <br>
 *  certain threshold, defined with processor parameter ECut, are retained. <br>
 *    @author A. Raspereza (DESY)<br>
 *    @version $Id: TrackCheater.h,v 1.4 2005-08-07 16:22:39 gaede Exp $<br>
 */
class TrackCheater : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TrackCheater ; }
  
  
  TrackCheater() ;
  
  /**  
   * Initialization
   */
  virtual void init() ;
  
  /** Run header processor.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Event processor.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  int _nRun ;
  int _nEvt ;

  int _fitTrueTrack;
  int _minimal_hits;
  float _chi2Cut;

  std::string _trueTracksCollection;
  std::vector<std::string> _trackerHitCollections;
  std::string _relCollection;
 
  float _bField;
  float _eCut;
  float _hitToHelixCut;
  float _hitToHelixInFit;
 
} ;

#endif



