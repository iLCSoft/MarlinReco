#ifndef FULLLDCTRACKING_H
#define FULLLDCTRACKING_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include "ClusterExtended.h"
#include "TrackExtended.h"
#include "TrackerHitExtended.h"
#include "HelixClass.h"
#include "GroupTracks.h"

using namespace lcio ;
using namespace marlin ;


class FullLDCTracking : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new FullLDCTracking ; }  
  FullLDCTracking() ;  
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  void prepareVectors( LCEvent * evt );
  void CleanUp();
  void MergeTPCandSiTracks();
  TrackExtended * CombineTracks(TrackExtended * tpcTrk, TrackExtended * siTrk);
  void Sorting(TrackExtendedVec & trackVec);
  void SelectCombinedTracks();
  void AddNotCombinedTracks(LCEvent * evt);

 protected:

  int _nRun ;
  int _nEvt ;

  int _simpleHelixFit;

  std::string _TPCTrackCollection;
  std::string _SiTrackCollection;
  std::string _VTXTrackerHitCollection;
  std::string _SITTrackerHitCollection;
  std::string _FTDTrackerHitCollection;
  std::string _TPCTrackerHitCollection;

  TrackExtendedVec _allSiTracks;
  TrackExtendedVec _allTPCTracks;
  TrackExtendedVec _allCombinedTracks;
  TrackExtendedVec _trkImplVec;
  TrackerHitExtendedVec _allTPCHits;
  TrackerHitExtendedVec _allVTXHits;
  TrackerHitExtendedVec _allFTDHits;
  TrackerHitExtendedVec _allSITHits;

  float _resolutionRPhi_TPC,_resolutionZ_TPC;
  float _resolutionRPhi_VTX,_resolutionZ_VTX;
  float _resolutionRPhi_FTD,_resolutionZ_FTD;
  float _resolutionRPhi_SIT,_resolutionZ_SIT;

  float PI, PIOVER2, TWOPI;

  float _deltaPhiForTracks,_deltaQForTracks;
  float _angleCutForMerging;

  float _bField;

  float _chi2PrefitCut;
  float _chi2FitCut;

  int _createMap;

} ;

#endif



