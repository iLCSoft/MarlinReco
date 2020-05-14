#ifndef DDStripSplitter_h
#define DDStripSplitter_h 1
#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include "TVector3.h"
#include <vector>
#include <cmath>
#include "UTIL/CellIDDecoder.h"
#include "EVENT/CalorimeterHit.h"
#include "IMPL/LCCollectionVec.h"
#include <IMPL/LCFlagImpl.h>

#include "DDRec/DetectorData.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;

/*
implementation of Strip Splitting Algorithm, adapted for dd4hep based ILD models.
D. Jeans, Nov 2018.
 */

class DDStripSplitter : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new DDStripSplitter ; }
  
  DDStripSplitter() ;

  DDStripSplitter( const DDStripSplitter& ) = delete;
  DDStripSplitter& operator=(const DDStripSplitter&) = delete;
  
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
  
  virtual void setupGeometry();
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

 protected:

  /** Input collection name.
   */
  std::string _mcParticleCollectionName{};

  std::string _ecalCollectionEvenLayers{};
  std::string _ecalCollectionOddLayers{};

  std::string _inputRelationsColEven{};
  std::string _inputRelationsColOdd{};

  // output collection names
  std::string _stripIntersecCollName{};
  std::string _evenStripEndsCollName{};
  std::string _oddStripEndsCollName{};

  std::string _splitEcalCollection{};
  std::string _unsplitEcalCollection{};

  std::string _splitEcalRelCol{};

  std::pair < TVector3, TVector3 > getStripEnds(CalorimeterHit* hit, int orientation, bool barrel);
  TVector3 stripIntersect(CalorimeterHit* hit0, TVector3 axis0, CalorimeterHit* hit1, TVector3 axis1);
  std::vector <CalorimeterHit*> getVirtualHits(LCEvent* evt, CalorimeterHit* hit, int orientation, bool barrel );

  CellIDDecoder<CalorimeterHit>* _decoder{}; 
  CellIDDecoder<CalorimeterHit>* _decoder2{}; 

  float _stripLength{};
  float _stripWidth{};
  float _stripAspectRatio{};
  float _cellSize{};
  int   _symmetry{};
  int   _nVirtual{};
  int   _ecalStrip_default_nVirt{};
  bool _isBarrel{};

  bool _saveIntersections{};
  IMPL::LCCollectionVec* intersectionHits{};
  IMPL::LCCollectionVec* stripEndsEvenCol{};
  IMPL::LCCollectionVec* stripEndsOddCol{};

  int _evenIsTransverse{};
  enum {TRANSVERSE=0, LONGITUDINAL};

  dd4hep::rec::LayeredCalorimeterData* _caloGeomData{};

  LCFlagImpl _flag_rel{};
  std::string _cellIDLayerString{};
  std::string _cellIDModuleString{};
  std::string _cellIDStaveString{};

};


#endif



