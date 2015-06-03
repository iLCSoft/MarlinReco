#ifndef hybridRecoProcessor_h
#define hybridRecoProcessor_h 1
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

using namespace lcio ;
using namespace marlin ;
using namespace std;


/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id$ 
 */

class hybridRecoProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new hybridRecoProcessor ; }
  
  hybridRecoProcessor() ;
  
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

  std::vector <std::string> _ecalCollectionsCells;
  std::vector <std::string> _ecalCollectionsTranStrips;
  std::vector <std::string> _ecalCollectionsLongStrips;

  std::pair < TVector3, TVector3 > getStripEnds(CalorimeterHit* hit, int orientation, bool barrel);
  TVector3 stripIntersect(CalorimeterHit* hit0, TVector3 axis0, CalorimeterHit* hit1, TVector3 axis1);
  std::vector <CalorimeterHit*> getVirtualHits(LCEvent* evt, CalorimeterHit* hit, int orientation, bool barrel );

  CellIDDecoder<CalorimeterHit>* _decoder; 
  CellIDDecoder<CalorimeterHit>* _decoder2; 

  bool  _makePlots;
  float _stripLength;
  float _stripWidth;
  float _stripAspectRatio;
  float _cellSize;
  int   _symmetry;
  int   _nVirtual;
  int _ecalStrip_default_nVirt;

  bool _saveIntersections;
  IMPL::LCCollectionVec* intersectionHits;
  IMPL::LCCollectionVec* stripEndsTransCol;
  IMPL::LCCollectionVec* stripEndsLongCol;

  TFile* _fout;
  TH2F* h_phiModuleCheck;
  TH2F* h_phiThetaMC;

  TH1F* h_stripDist_intercept;
  TH1F* h_stripDist_nointercept;

  TH2F* h_stavemodule[2];
  TH1F* h_layer[2];

  TH2F* h_staveX[2];
  TH2F* h_staveY[2];
  TH2F* h_staveZ[2];

  TH2F* h_moduleX[2];
  TH2F* h_moduleY[2];
  TH2F* h_moduleZ[2];


  TH2F* h_cth_phi[2][10][10];
  TH2F* h_XY[2][10][10];

  enum {TRANSVERSE, LONGITUDINAL};

};


#endif



