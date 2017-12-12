/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

/*
Evolved version of TPCDigi that provides additional functionality to deal with background. Couple to the Mokka Sensitive Detector Driver TPCSD03.cc

SJA:FIXME: Still needs to be tidied up for production release.

Three cases can be consider in the treatment of SimTrackerHits
i)   A clean isolated hit; this will be smeared according to the parametric point resolution and converted to a TrackerHit
ii)  Two or Three hits which are considered to be closer than the double hit resolution and which therefore cannot be viewed as seperable hits. These will be merged and be assigned a large associated measurement error.
iii) A continuous set of hits within one pad row which cannot be resolved as single hits, these are condidered to be charaterisable as background hits created by extremely low pt charged particles (pt < 10MeV) and therefore are removed from the hit collection.

The Driver has been modified to take an additional collection of SimTrackerHits which are produced by the Mokka TPC Sensitive Driver TPCSD03.cc. These hits are produced for particles which have very low pt and often do not move outside of the dimensions of a single pad row. These hits need to be treated differently as they do not cross any geometric boundaries in a Padrow based TPC Geometry. This negates the need to voxalise the TPC in Geant4 which has proved in the past to be prohibitive in terms of processing time due to the vastly increased number of geometric volumes. 

Steve Aplin 26 June 2009 (DESY)

*/

#ifndef TPCDigiProcessor_h
#define TPCDigiProcessor_h 1

#include <marlin/Processor.h>
#include <lcio.h>


#include <string>
#include <gsl/gsl_rng.h>

#ifdef MARLIN_USE_AIDA

#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>



//#define DIGIPLOTS



#ifdef DIGIPLOTS
// includes all AIDA header files
#include <AIDA/AIDA.h>
#endif

#endif

#include <vector>
#include <map>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDEncoder.h>

#include "CLHEP/Vector/TwoVector.h"
class Voxel_tpc;




using namespace lcio ;
using namespace marlin ;
#ifdef MARLIN_USE_AIDA
using namespace AIDA ;
#endif



/** ====== TPCDigiProcessor ====== <br>
 *
 * This Processor depends on Circle.h from MarlinUtil
 * 
 * Caution: This digitiser presently does not process space-point like SimTrackerHits which have been flagged with CellIDs set to the negetive row number. This must be implemented in future. 
 *Produces TPC TrackerHit collection from SimTrackerHit collection, smeared in r-phi and z. 
 * Double hits are identified but are currently not added to the collection. This may be change 
 * at a later date when criteria for their seperation is defined. The resolutions are defined in 
 * the GEAR stearing file. 
 *
 * Resolution in r-phi is calculated according to the formular <br>
 * sigma_{point}^2 = sigma_0^2 + Cd^2/N_{eff} * L_{drift}
 * Cd^2/N_{eff}} = 25^2/(22/sin(theta)*h/6mm)
 * Cd = 25 ( microns / cm^(1/2) )
 * (this is for B=4T, h is the pad height = pad-row pitch in mm,
 * theta is the polar angle)       
 *
 * At the moment resolution in z assumed to be independent of drift length. <br>
 *
 * The type of TPC TrackerHit is set to 500 via method TrackerHitImpl::setType(int type) <br>
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collections of SimTrackerHits in TPC <br>
 * <h4>Output</h4>
 * Processor produces collection of digitized TrackerHits in TPC <br>
 * @param CollectionName The name of input SimTrackerHit collection <br>
 * (default name STpc01_TPC)
 * @param RejectCellID0 Whether or not to reject SimTrackerHits with Cell ID 0. Mokka drivers
 * TPC00-TPC03 encoded the pad row number in the cell ID, which should always be non-zero anyway.
 * Drivers TPC04 and TPC05 do not simulate pad rows and thus have the cell ID set to zero for all hits.
 * You will need to set RejectCellID0 to 0 in order to use this processor with these drivers, but note
 * that the implications for track reconstruction are not strictly defined. Mokka driver TPC06 uses
 * a mixed approach with one hit per pad row having non-zero cell ID, extra hits having 0 cell ID.
 * Typically, unless you use TPC04 or TPC05, you should not touch this parameter. <br>
 * (default value 1)
 * @param TPCTrackerHitsCol The name of output collection of TrackerHits <br>
 * (default name TPCTrackerHits) <br>
 * <br>
 * @authors S. Aplin, DESY and A.Raspereza, MPI
 *
 * Changed 7/9/07 so that the const and diffusion resolution terms are taken as processor parameters rather than the gear file.
 * The parameters _pixZ and pixRP were also changed from gear parameters to processor parameters
 * clare.lynch@bristol.ac.uk
 *
 */
class TPCDigiProcessor : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new TPCDigiProcessor ; }
  
  
  TPCDigiProcessor() ;
  
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
  
  void writeVoxelToHit( Voxel_tpc* aVoxel) ;  
  void writeMergedVoxelsToHit( std::vector <Voxel_tpc*>* hitList ) ;  
  void plotHelixHitResidual(MCParticle *mcp, CLHEP::Hep3Vector *thisPointRPhi);
  double getPadPhi( CLHEP::Hep3Vector* thisPointRPhi, CLHEP::Hep3Vector* firstPointRPhi, CLHEP::Hep3Vector* middlePointRPhi, CLHEP::Hep3Vector* lastPointRPhi);
  double getPadTheta( CLHEP::Hep3Vector* firstPointRPhi, CLHEP::Hep3Vector* middlePointRPhi, CLHEP::Hep3Vector* lastPointRPhi );

protected:

  /** Input collection name.
   */
  std::string _padRowHitColName{};
  std::string _spacePointColName{};
  std::string _lowPtHitscolName{};

  
  /** Output collection name.
   */
  std::string _TPCTrackerHitsCol{};
  std::string _outRelColName{};

  bool _use_raw_hits_to_store_simhit_pointer{};
  
  int _rejectCellID0{};
  float _padWidth{};

  int _nRun{};
  int _nEvt{};

  EVENT::MCParticle* _mcp{};
  EVENT::MCParticle* _previousMCP{};
  EVENT::MCParticle* _nextMCP{};
  EVENT::MCParticle* _nMinus2MCP{};
  EVENT::MCParticle* _nPlus2MCP{};   

  SimTrackerHit* _SimTHit{};
  SimTrackerHit* _previousSimTHit{};
  SimTrackerHit* _nextSimTHit{};
  SimTrackerHit* _nPlus2SimHit{};
  SimTrackerHit* _nMinus2SimHit{};

  // gsl random number generator
  gsl_rng * _random{};

  bool _dontEncodeSide{};

  float _pointResoRPhi0{}; // Coefficient for RPhi point res independant of drift length 
  float _pointResoPadPhi{}; // Coefficient for the point res dependance on relative phi angle to the pad verticle 
  float _diffRPhi{}; // Coefficient for the rphi point res dependance on diffusion 
  int   _nEff{}; // number of effective electrons 


  float _pointResoZ0{}; // Coefficient Z point res independant of drift length 
  float _diffZ{}; // Coefficient for the Z point res dependance on diffusion 

  float _binningZ{};
  float _binningRPhi{};
  float _doubleHitResZ{};
  float _doubleHitResRPhi{};
  int _maxMerge{};

  int _nRechits{};

  std::vector< std::vector <Voxel_tpc *> > _tpcRowHits{};
  std::map< Voxel_tpc *,SimTrackerHit *> _tpcHitMap{};
  std::vector<float> _length{};
  int lenpos{};

  LCCollectionVec* _trkhitVec{};
  LCCollectionVec* _relCol{};  
  CellIDEncoder<TrackerHitImpl>* _cellid_encoder{};

  int  _NSimTPCHits{};
  int  _NBackgroundSimTPCHits{};
  int  _NPhysicsSimTPCHits{};
  int  _NPhysicsAbove02GeVSimTPCHits{};
  int  _NPhysicsAbove1GeVSimTPCHits{};
  int  _NRecTPCHits{};
  
  int  _NLostPhysicsTPCHits{};
  int  _NLostPhysicsAbove02GeVPtTPCHits{};
  int  _NLostPhysicsAbove1GeVPtTPCHits{};
  int  _NRevomedHits{};


#ifdef DIGIPLOTS
  IAnalysisFactory * _AF{};
  ITreeFactory * _TRF{};
  ITree * _TREE{};
  IHistogramFactory * _HF{};
  IHistogram1D * _phiDiffHisto{};
  IHistogram1D * _thetaDiffHisto{};
  IHistogram1D * _phiRelHisto{};
  IHistogram1D * _thetaRelHisto{};

  IHistogram1D * _phiDistHisto{};
  IHistogram1D * _rPhiPullHisto{};
  IHistogram1D * _rPhiDiffHisto{};
  IHistogram1D * _zDiffHisto{};
  IHistogram1D * _zPullHisto{};
  IHistogram2D * _zSigmaVsZHisto{};
  IHistogram1D * _zSigmaHisto{};
  IHistogram1D * _rPhiSigmaHisto{};
  IHistogram1D * _radiusCheckHisto{};
  IHistogram1D * _ResidualsRPhiHisto{};

  IHistogram1D * _NSimTPCHitsHisto{};
  IHistogram1D * _NBackgroundSimTPCHitsHisto{};
  IHistogram1D * _NPhysicsSimTPCHitsHisto{};
  IHistogram1D * _NPhysicsAbove02GeVSimTPCHitsHisto{};
  IHistogram1D * _NPhysicsAbove1GeVSimTPCHitsHisto{};
  IHistogram1D * _NRecTPCHitsHisto{};

  IHistogram1D * _NLostPhysicsTPCHitsHisto{};
  IHistogram1D * _NLostPhysicsAbove02GeVPtTPCHitsHisto{};
  IHistogram1D * _NLostPhysicsAbove1GeVPtTPCHitsHisto{};
  IHistogram1D * _NRevomedHitsHisto{};

  IHistogram1D * _NKeptPhysicsTPCHitsHistoPercent{};
  IHistogram1D * _NKeptPhysicsAbove02GeVPtTPCHitsHistoPercent{};
  IHistogram1D * _NKeptPhysicsAbove1GeVPtTPCHitsHistoPercent{};

#endif


} ;




#endif



