#ifndef DDStripSplitter_h
#define DDStripSplitter_h 1
#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include <vector>
#include <cmath>
#include "UTIL/CellIDDecoder.h"
#include "EVENT/CalorimeterHit.h"
#include "IMPL/LCCollectionVec.h"

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
  std::vector <std::string> _ecalCollectionsEvenLayers{};
  std::vector <std::string> _ecalCollectionsOddLayers{};
  std::string _mcParticleCollectionName{};

  // output collection names
  std::string _stripIntersecCollName{};
  std::string _evenStripEndsCollName{};
  std::string _oddStripEndsCollName{};

  std::pair < TVector3, TVector3 > getStripEnds(CalorimeterHit* hit, int orientation, bool barrel);
  TVector3 stripIntersect(CalorimeterHit* hit0, TVector3 axis0, CalorimeterHit* hit1, TVector3 axis1);
  std::vector <CalorimeterHit*> getVirtualHits(LCEvent* evt, CalorimeterHit* hit, int orientation, bool barrel );

  CellIDDecoder<CalorimeterHit>* _decoder{}; 
  CellIDDecoder<CalorimeterHit>* _decoder2{}; 

  bool  _makePlots{};
  float _stripLength{};
  float _stripWidth{};
  float _stripAspectRatio{};
  float _cellSize{};
  int   _symmetry{};
  int   _nVirtual{};
  int   _ecalStrip_default_nVirt{};

  bool _saveIntersections{};
  IMPL::LCCollectionVec* intersectionHits{};
  IMPL::LCCollectionVec* stripEndsEvenCol{};
  IMPL::LCCollectionVec* stripEndsOddCol{};

  TFile* _fout{};
  TH2F* h_phiModuleCheck{};
  TH2F* h_phiThetaMC{};

  TH1F* h_stripDist_intercept{};
  TH1F* h_stripDist_nointercept{};

  TH2F* h_stavemodule[2]{};
  TH1F* h_layer[2]{};

  TH2F* h_staveX[2]{};
  TH2F* h_staveY[2]{};
  TH2F* h_staveZ[2]{};

  TH2F* h_moduleX[2]{};
  TH2F* h_moduleY[2]{};
  TH2F* h_moduleZ[2]{};

  TH2F* h_cth_phi[2][10][10]{};
  TH2F* h_XY[2][10][10]{};

  int _evenIsTransverse{};
  enum {TRANSVERSE=0, LONGITUDINAL};

  dd4hep::rec::LayeredCalorimeterData* _caloGeomData{};

  std::string _cellIDLayerString{};
  std::string _cellIDModuleString{};
  std::string _cellIDStaveString{};

};


#endif



