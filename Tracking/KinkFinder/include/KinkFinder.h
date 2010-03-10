#ifndef KinkFinder_H
#define KinkFinder_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include "TrackPair.h"
#include "HelixClass.h"

typedef struct {
  float x;
  float y;
  float z;
} vec3;

typedef struct {
  int   tracki;
  int   trackj;
  float vtx[3];
  float p[3];
  float mass;
  float distance;
  int   pdgCode;
} twoTrackIntersection_t;

using namespace lcio ;
using namespace marlin ;


/** KinkFinder Processor <br>
 *  KinkFinder processor identify kinked tracks <br>
 */
class KinkFinder : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new KinkFinder ; }
  
  
  KinkFinder() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;

 protected:

  void Sorting( TrackPairVec & trkPairVec );
  
  float kinkMass(HelixClass* parent, HelixClass* daughter, float daughterMass, float neutralMass);


  int _nRun ;
  int _nEvt ;
  
  std::string _trackColName;

  std::string _kinkVertexColName;
  std::string _prongVertexColName;
  std::string _splitVertexColName;

  std::string _kinkRecoPartColName;
  std::string _prongRecoPartColName;
  std::string _splitRecoPartColName;
  
  float _rKinkCut;
  int   _minTrackHits;
  int   _maxDeltaTpcLayers;
  int   _debugPrinting;
  float _kaonDecayMassCut;
  float _pionDecayMassCut;
  float _sigmaDecayMassCut;
  float _hyperonDecayMassCut;
  float _sigmaTimeCut;
  float _hyperonTimeCut;
  float _tightDrCutTPC;
  float _veryTightDrCutTPC;
  float _drCutTPC;
  float _drCutSIT;
  float _looseDrCutSIT;
  float _maxSplitTrackFracDeltaP;
  float _maxSplitTrackDeltaP;
  float _minELambda;


  float _bField;
  float _tpcInnerR;
  float _tpcOuterR;
  float _tpcZmax;
  int   _tpcMaxRow;
  int   _nLayersSIT;
  int   _nLayersVTX;
  std::vector<float> _rSIT;
  std::vector<float> _rVTX;


} ;

#endif



