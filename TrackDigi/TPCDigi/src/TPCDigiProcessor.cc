/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

#include "TPCDigiProcessor.h"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#include <gsl/gsl_randist.h>
#include "marlin/VerbosityLevels.h"

#include "Circle.h"
#include "SimpleHelix.h"
#include"constants.h"
#include "LCCylinder.h"
#include <IMPL/LCFlagImpl.h>

//stl exception handler
#include <stdexcept>
#include "constants.h"
#include "voxel.h"

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/BField.h>
//

using namespace lcio ;
using namespace marlin ;
using namespace constants ;
using namespace std ;

#ifdef MARLIN_USE_AIDA
using namespace AIDA ;
#endif


TPCDigiProcessor aTPCDigiProcessor ;

bool compare_phi( Voxel_tpc* a, Voxel_tpc* b)  { 
  return ( a->getPhiIndex() < b->getPhiIndex() ) ; 
} 

bool compare_z( Voxel_tpc* a, Voxel_tpc* b) { 
  return ( a->getZIndex() < b->getZIndex() ) ; 
} 
  

TPCDigiProcessor::TPCDigiProcessor() : Processor("TPCDigiProcessor") 
{
  
  // modify processor description
  _description = "Produces TPC TrackerHit collection from SimTrackerHit collection, smeared in RPhi and Z. A search is made for adjacent hits on a pad row, if they are closer in z and r-phi than the steering parameters _doubleHitResRPhi (default value 2.0 mm) and _doubleHitResZ (default value 5.0 mm) they are considered to overlap. Clusters of hits smaller than _maxMerge (default value 3) are merged into a single tracker hit, with the position given as the average poision of the hits in phi and in z. Clusters which have _maxMerge hits or more are determined to be identifiable as multiple hits, and are not added to the tracker hit collection. This of course means that good hits caught up in a cluster of background hits will be lossed." ;
  
  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "TPCPadRowHitCollectionName" , 
                           "Name of the default pad-row based SimTrackerHit collection"  ,
                           _padRowHitColName ,
                           std::string("TPCCollection") ) ;

  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "TPCSpacePointCollectionName" , 
                           "Name of the additional space point collection which provides additional guide hits between pad row centers."  ,
                           _spacePointColName ,
                           std::string("TPCSpacePointCollection") ) ;

  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "TPCLowPtCollectionName" , 
                           "Name of the LowPt SimTrackerHit collection Produced by Mokka TPC Driver TPC0X"  ,
                           _lowPtHitscolName ,
                           std::string("TPCLowPtCollection") ) ;

  registerOutputCollection( LCIO::TRACKERHIT,
                            "TPCTrackerHitsCol" , 
                            "Name of the Output TrackerHit collection"  ,
                            _TPCTrackerHitsCol ,
                            std::string("TPCTrackerHits") ) ;

  registerProcessorParameter( "PointResolutionPadPhi" ,
                              "Pad Phi Resolution constant in TPC"  ,
                              _pointResoPadPhi ,
                              (float)0.900) ;

  registerProcessorParameter( "RejectCellID0" ,
                              "whether or not to use hits without proper cell ID (pad row)"  ,
                              _rejectCellID0 ,
                              (int)1) ;

  registerProcessorParameter( "PointResolutionRPhi" ,
                              "R-Phi Resolution constant in TPC"  ,
                              _pointResoRPhi0 ,
                              (float)0.050) ;

  registerProcessorParameter( "DiffusionCoeffRPhi" ,
                              "R-Phi Diffusion Coefficent in TPC"  ,
                              _diffRPhi ,
                              (float)0.025) ;

  registerProcessorParameter( "N_eff" ,
                              "Number of Effective electrons per pad in TPC"  ,
                              _nEff ,
                              (int)22) ;

  registerProcessorParameter( "PointResolutionZ" ,
                              "TPC Z Resolution Coefficent independent of diffusion"  ,
                              _pointResoZ0 ,
                              (float)0.4) ;

  registerProcessorParameter( "DiffusionCoeffZ" ,
                              "Z Diffusion Coefficent in TPC"  ,
                              _diffZ ,
                              (float)0.08) ;

  registerProcessorParameter( "HitSortingBinningZ" ,
                              "Defines spatial slice in Z"  ,
                              _binningZ ,
                              (float)5.0) ;

  registerProcessorParameter( "HitSortingBinningRPhi" ,
                              "Defines spatial slice in RP"  ,
                              _binningRPhi ,
                              (float)2.0) ;


  registerProcessorParameter( "DoubleHitResolutionZ" ,
                              "Defines the minimum distance for two seperable hits in Z"  ,
                              _doubleHitResZ ,
                              (float)5.0) ;

  registerProcessorParameter( "DoubleHitResolutionRPhi" ,
                              "Defines the minimum distance for two seperable hits in RPhi"  ,
                              _doubleHitResRPhi ,
                              (float)2.0) ;

  registerProcessorParameter( "MaxClusterSizeForMerge" ,
                              "Defines the maximum number of adjacent hits which can be merged"  ,
                              _maxMerge ,
                              (int)3) ;
}


void TPCDigiProcessor::init() 
{ 
  

  // From GNU documentation:
  // A replacement for the standard terminate_handler which prints 
  // more information about the terminating exception (if any) on stderr. Call ...
  std::set_terminate (__gnu_cxx::__verbose_terminate_handler);

#ifdef DIGIPLOTS
  //#ifdef EXPERTCHECKPLOTS
  /// Hook an AIDA implementation -----------------------------------------------

  // First create a pointer to the "IAnalysisFactory" of a specific AIDA
  // implementation. This factory can then be used to produce all other 
  // factories.
  _AF = AIDA_createAnalysisFactory();
  
  // Create a ITreeFactory. -----------------------------------------------------
  // A ITree can be used to store AIDA objects in memory or on disk.
  
  _TRF = _AF->createTreeFactory();
  
  /// Create a ITree object which is bound to a file. ---------------------------
  // You must always create a "ITree" object to create any other factory.
  /*
   * Creates a new tree and associates it with a store.
   * The store is assumed to be read/write.
   * The store will be created if it does not exist.
   * @param storeName The name of the store, if empty (""), the tree is 
   *                  created in memory and therefore will not be associated 
   *                  with a file.
   * @param storeType Implementation specific string, may control store type
   * @param readOnly If true the store is opened readonly, an exception if it 
   *                 does not exist
   * @param createNew If false the file must exist, if true the file will be 
   *                  created
   * @param options Other options, currently are not specified
   */
  // ITree * ITreeFactory::create(const std::string & storeName, 
  //                              const std::string & storeType = "", 
  //                              bool readOnly = false, 
  //                              bool createNew = false, 
  //                              const std::string & options = "") ;

  _TREE = _TRF->create("TPCDigi.root",
                     "root",
                     false,
                     true);

  /// Create an IHistogramFactory which is bound to the tree "*_TREE". -----------

  /*
   * Create an IHistogramFactory.
   * @param tree The ITree which created histograms will be associated to.
   * @return     The IHistogramFactory.
   */
  // IHistogramFactory * IAnalysisFactory::createHistogramFactory(ITree & tree);

  _HF = _AF->createHistogramFactory(*_TREE);
  
  _TREE->mkdir("Histograms");
 
  /*
   * Create a IHistogram1D.
   * @param path      The path of the created IHistogram. The path can either 
   *                  be a relative or full path.
   *                  ("/folder1/folder2/dataName" and 
   *                  "../folder/dataName" are valid paths).
   *                  All the directories in the path must exist. The 
   *                  characther `/` cannot be used in names; it is only 
   *                  used to delimit directories within paths.
   * @param title     The title of the IHistogram1D.
   * @param nBins     The number of bins of the x axis.
   * @param lowerEdge The lower edge of the x axis.
   * @param upperEdge The upper edge of the x axis.
   * @param options   The options for the IHistogram1D. The default is "".
   *                  "type=efficiency" for an efficiency IHistogram1D.
   * @return          The newly created IHistogram1D.
   */
  


  _phiDiffHisto = _HF->createHistogram1D("Histograms/phi_diff",
                                       "Calculated Phi - Track Phi",
                                       201, -0.05, 0.05);

  _thetaDiffHisto = _HF->createHistogram1D("Histograms/theta_diff",
                                         "Calculated Theta - Track Theta",
                                         201, -0.05, 0.05);
 
  _phiRelHisto = _HF->createHistogram1D("Histograms/padPhi",
                                      "Phi Relative to the Pad",
                                      201, 0.0, 6.3);

  _thetaRelHisto = _HF->createHistogram1D("Histograms/padtheta",
                                        "Theta Relative to the pad",
                                        201, 0.0, 6.3);

  _rPhiDiffHisto = _HF->createHistogram1D("Histograms/rPhiDiff",
                                     "rPhi_rec - rPhi_sim",
                                     201, -1.0, 1.0);

  _zDiffHisto = _HF->createHistogram1D("Histograms/zDiff",
                                     "Z_rec - Z_sim",
                                     201, -1.0, 1.0);

  _zPullHisto = _HF->createHistogram1D("Histograms/zPull",
                                       "(z_rec - z_sim) / Sigma_z",
                                       201, -10.0, 10.0);

  _phiDistHisto = _HF->createHistogram1D("Histograms/phiDist",
                                       "phi_rec - Phi_sim",
                                       201, -1.0, 1.0);

  _rPhiPullHisto = _HF->createHistogram1D("Histograms/rPhiPull",
                                       "(rPhi_rec - rPhi_sim) / Sigma_rPhi",
                                       201, -10.0, 10.0);

  _zSigmaVsZHisto = _HF->createHistogram2D("Histograms/zSigmaVsZ",
                                     "z Sigma vs Z ",
                                      3000, 0.0, 3000.0,
                                      201, -0.20, 5.20);

  _zSigmaHisto = _HF->createHistogram1D("Histograms/zSigma",
                                     "z Sigma ",
                                      201, -0.20, 5.20);

  _rPhiSigmaHisto = _HF->createHistogram1D("Histograms/rPhiSigma",
                                     "rPhi Sigma",
                                         201, -0.20, 0.20);

  _radiusCheckHisto = _HF->createHistogram1D("Histograms/radiusCheck",
                                     "R_hit - TPC Rmin - ((RowIndex + 0.5 )* padheight)",
                                         201, -0.20, 0.20);

  _ResidualsRPhiHisto = _HF->createHistogram1D("Histograms/ResidualsRPhi",
                                        "MC Track Phi - Hit Phi",
                                        50, -0.001, 0.001);

  _NSimTPCHitsHisto = _HF->createHistogram1D("Histograms/SimTPCHits",
                                        "Number of SimTPC Hits",
                                        100, 0.0, 1000000.0);

  _NBackgroundSimTPCHitsHisto = _HF->createHistogram1D("Histograms/NBackgroundSimTPCHits",
                                        "Number of Background SimTPC Hits",
                                        100, 0.0, 1000000.0);

  _NPhysicsSimTPCHitsHisto = _HF->createHistogram1D("Histograms/NPhysicsSimTPCHits",
                                        "Number of PhysicsSimTPC Hits",
                                        100, 0.0, 100000.0);

  _NPhysicsAbove02GeVSimTPCHitsHisto = _HF->createHistogram1D("Histograms/NPhysicsAbove02GeVTPCHits",
                                        "Number of PhysicsSimTPC Hits above 0.2GeV pt",
                                        100, 0.0, 100000.0);

  _NPhysicsAbove1GeVSimTPCHitsHisto = _HF->createHistogram1D("Histograms/NPhysicsAbove1GeVPtTPCHits",
                                        "Number of PhysicsSimTPC Hits above 1.0 GeV pt",
                                        100, 0.0, 100000.0);

  _NRecTPCHitsHisto = _HF->createHistogram1D("Histograms/NRecTPCHits",
                                        "Number of Rec TPC Hits",
                                        50, 0.0, 100000.0);

  _NLostPhysicsTPCHitsHisto = _HF->createHistogram1D("Histograms/NLostPhysicsTPCHits",
                                        "Number of PhysicsSimTPC Hits Lost",
                                        100, 0.0, 5000.0);

  _NLostPhysicsAbove02GeVPtTPCHitsHisto = _HF->createHistogram1D("Histograms/NLostPhysicsAbove02GeVPtTPCHits",
                                        "Number of PhysicsSimTPC Hits Lost above 0.2 GeV pt",
                                        100, 0.0, 5000.0);

  _NLostPhysicsAbove1GeVPtTPCHitsHisto = _HF->createHistogram1D("Histograms/NLostPhysicsAbove1GeVPtTPCHits",
                                        "Number of PhysicsSimTPC Hits Lost above 1.0 GeV pt",
                                        100, 0.0, 1000.0);

  _NRevomedHitsHisto = _HF->createHistogram1D("Histograms/NRevomedHits",
                                        "Number of Removed TPC hits",
                                        100, 0.0, 1000000.0);


  _NKeptPhysicsTPCHitsHistoPercent = _HF->createHistogram1D("Histograms/NKeptPhysicsTPCHitsPercent",
                                                     "Number of PhysicsSimTPC Hits Kept",
                                                     303, 0.0, 1.01);

  _NKeptPhysicsAbove02GeVPtTPCHitsHistoPercent = _HF->createHistogram1D("Histograms/NKeptPhysicsAbove02GeVPtTPCHitsPercent",
                                                                 "Number of PhysicsSimTPC Hits Kept above 0.2 GeV pt",
                                                                 303, 0.0, 1.01);

  _NKeptPhysicsAbove1GeVPtTPCHitsHistoPercent = _HF->createHistogram1D("Histograms/NKeptPhysicsAbove1GeVPtTPCHitsPercent",
                                                                "Number of PhysicsSimTPC Hits Kept above 1.0 GeV pt",
                                                                303, 0.0, 1.01);
  
#endif  


  printParameters() ;

  //intialise random number generator 
  _random = gsl_rng_alloc(gsl_rng_ranlxs2);

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void TPCDigiProcessor::processRunHeader( LCRunHeader* run) 
{ 
  _nRun++ ;
} 

void TPCDigiProcessor::processEvent( LCEvent * evt ) 
{ 

  int numberOfVoxelsCreated(0);

  _NSimTPCHits = 0;
  _NBackgroundSimTPCHits = 0;
  _NPhysicsSimTPCHits = 0;
  _NPhysicsAbove02GeVSimTPCHits = 0;
  _NPhysicsAbove1GeVSimTPCHits = 0;
  _NRecTPCHits = 0;

  _NLostPhysicsTPCHits = 0;
  _NLostPhysicsAbove02GeVPtTPCHits = 0;
  _NLostPhysicsAbove1GeVPtTPCHits = 0;
  _NRevomedHits = 0;

  static bool firstEvent = true;
  _tpcHitMap.clear();
  _tpcRowHits.clear();
  
  if(firstEvent==true) streamlog_out(MESSAGE) << "TPCDigiProcessor called for first event" << endl;

  firstEvent = false ;

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  // this gets the center of the first pad in the pad layout
  const gear::Vector2D padCoord = padLayout.getPadCenter(1) ;
  
  _trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT )  ;

  // first deal with the pad-row based hits from Mokka 
  LCCollection* STHcol = 0 ;
  try{
    STHcol = evt->getCollection( _padRowHitColName ) ;
  }
  catch(DataNotAvailableException &e){
  }
  
  if( STHcol != 0 ){
  
    int n_sim_hits = STHcol->getNumberOfElements()  ;

    LCFlagImpl colFlag( STHcol->getFlag() ) ;    

    _NSimTPCHits = n_sim_hits;
      
    streamlog_out(DEBUG4) << "number of Pad-Row based SimHits = " << n_sim_hits << std::endl;

    // this assumes that the pad width in r*phi is the same for all pads  
    _padWidth = padLayout.getPadWidth(0)*padCoord[0];
    // set size of row_hits to hold (n_rows) vectors
    _tpcRowHits.resize(padLayout.getNRows());

    _mcp=NULL;         
    _previousMCP=NULL; 
    _nextMCP=NULL;     
    _nMinus2MCP=NULL;  
    _nPlus2MCP=NULL;   
    
    _SimTHit=NULL;
    _previousSimTHit=NULL;
    _nextSimTHit=NULL;
    _nMinus2SimHit=NULL;
    _nPlus2SimHit=NULL;

    for(int i=0; i< n_sim_hits; i++){
      
      _SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;

      float de_dx;
      double padPhi(0.0);
      double padTheta (0.0);
      
      CLHEP::Hep3Vector *precedingPoint = (NULL);
      CLHEP::Hep3Vector *followingPoint = (NULL);
      CLHEP::Hep3Vector *nMinus2Point = (NULL);
      CLHEP::Hep3Vector *nPlus2Point = (NULL);
      
      CLHEP::Hep3Vector thisPoint(_SimTHit->getPosition()[0],_SimTHit->getPosition()[1],_SimTHit->getPosition()[2]);
      
      _mcp = _SimTHit->getMCParticle() ; 

      if(_mcp){ // SJA:FIX: the fact that it is a physics hit relies on the fact that for overlay the pointer to the mcp is set to NULL. This distinction may not always be true ...
        ++_NPhysicsSimTPCHits ;
        const float *mom= _SimTHit->getMomentum() ;
        double ptSQRD = mom[0]*mom[0]+mom[1]*mom[1] ; 
        if( ptSQRD > (0.2*0.2) ) ++_NPhysicsAbove02GeVSimTPCHits ;
        if( ptSQRD > 1.0 )  ++_NPhysicsAbove1GeVSimTPCHits ;
        
#ifdef EXPERTCHECKPLOTS
        if(_mcp) plotHelixHitResidual(_mcp, &thisPoint);
#endif  
        
      } else {
        ++_NBackgroundSimTPCHits;
      }
      
      if(colFlag.bitSet(LCIO::THBIT_MOMENTUM)) {
        
        const float * mcpMomentum = _SimTHit->getMomentum() ;
        
        CLHEP::Hep3Vector *mom = new CLHEP::Hep3Vector(mcpMomentum[0],mcpMomentum[1],mcpMomentum[2]);

        const double bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
        const double FCT = 2.99792458E-4;
        
        const double pt = mom->perp();
        const double radius = pt / (FCT*bField);
        
        const double tanLambda = mom->z()/pt;
        //        const double cosLambda = mom->z()/mom->r();
        //        const double sinLambda = mom->perp()/mom->r();
        
        padPhi = fabs(thisPoint.deltaPhi(*mom));
        padTheta = mom->theta();
        
      } 
      else {
        
        // as the momentum vector is not available at the hits perform the old gymnastics        
        if (!_mcp) { // then there is no way to know if this hit has consecutive hits from the same MCParticle
          // so just set nominal values theta=phi=90 
          padTheta = twopi/4.0 ;
          padPhi = twopi/4.0 ;    
        }
        else{
          // if there is at least one more hit after this one, set the pointer to the MCParticle for the next hit
          if (i < (n_sim_hits-1) ) {
            _nextSimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i+1 ) ) ;
            _nextMCP     = _nextSimTHit->getMCParticle() ;
          }
          // if there is at least two more hits after this one, set the pointer to the MCParticle for the next but one hit
          if (i < (n_sim_hits-2) ) {
            _nPlus2SimHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i+2 ));
            _nPlus2MCP   = _nPlus2SimHit->getMCParticle() ;            
          }
          
          if      ( _mcp==_previousMCP && _mcp==_nextMCP )    { // middle hit of 3 from the same MCParticle 
            // use a function which returns phi and takes arguments first, second, and third points to calculate the circle 
            precedingPoint = new CLHEP::Hep3Vector(_previousSimTHit->getPosition()[0],_previousSimTHit->getPosition()[1],_previousSimTHit->getPosition()[2]);
            followingPoint = new CLHEP::Hep3Vector(_nextSimTHit->getPosition()[0],_nextSimTHit->getPosition()[1],_nextSimTHit->getPosition()[2]);
            
            padPhi = getPadPhi( &thisPoint, precedingPoint, &thisPoint, followingPoint);
            padTheta = getPadTheta(precedingPoint, &thisPoint, followingPoint);

            delete precedingPoint;
            delete followingPoint;
          }
          else if ( _mcp==_nextMCP     && _mcp==_nPlus2MCP )  { // first  hit of 3 from the same MCParticle
            followingPoint = new CLHEP::Hep3Vector(_nextSimTHit->getPosition()[0],_nextSimTHit->getPosition()[1],_nextSimTHit->getPosition()[2]);
            nPlus2Point = new CLHEP::Hep3Vector(_nPlus2SimHit->getPosition()[0],_nPlus2SimHit->getPosition()[1],_nPlus2SimHit->getPosition()[2]);
            
            padPhi = getPadPhi( &thisPoint, &thisPoint, followingPoint, nPlus2Point);
            padTheta = getPadTheta(&thisPoint, followingPoint, nPlus2Point);

            delete followingPoint;
            delete nPlus2Point;
            
          }
          else if ( _mcp==_previousMCP && _mcp==_nMinus2MCP ) { // last   hit of 3 from the same MCParticle
            nMinus2Point = new CLHEP::Hep3Vector(_nMinus2SimHit->getPosition()[0],_nMinus2SimHit->getPosition()[1],_nMinus2SimHit->getPosition()[2]);
            precedingPoint = new CLHEP::Hep3Vector(_previousSimTHit->getPosition()[0],_previousSimTHit->getPosition()[1],_previousSimTHit->getPosition()[2]);
            
            padPhi = getPadPhi( &thisPoint, nMinus2Point, precedingPoint, &thisPoint);
            padTheta = getPadTheta(nMinus2Point, precedingPoint, &thisPoint);

            delete nMinus2Point;            
            delete precedingPoint;
            
          }
          else{ // the hit is isolated as either a single hit, or a pair of hits, from a single MCParticle  
            padTheta = twopi/4.0 ;
            padPhi = twopi/4.0 ;    
          }
        }
#ifdef EXPERTCHECKPLOTS

        if(colFlag.bitSet(LCIO::THBIT_MOMENTUM)) {
          
          const float * mcpMomentum = _SimTHit->getMomentum() ;
          
          CLHEP::Hep3Vector mom(mcpMomentum[0],mcpMomentum[1],mcpMomentum[2]);
          
          double trackPhi = mom.phi();
          
          if(trackPhi<0.0) trackPhi=trackPhi+twopi;
          if(trackPhi>twopi) trackPhi=trackPhi-twopi;
          if(trackPhi>twopi/2.0) trackPhi = trackPhi - twopi/2.0 ;
          
          double localPhi = thisPoint.phi() - padPhi;
          
          _phiRelHisto->fill(padPhi);
          _phiDiffHisto->fill((fabs(localPhi - trackPhi))/trackPhi);
          _thetaRelHisto->fill(padTheta);
          _thetaDiffHisto->fill( (sin(padTheta) - sin(mom.theta()))/sin(mom.theta()) );
          
          streamlog_out(DEBUG3) << "track Phi = " << trackPhi * (360.0 / twopi) << endl; 
          streamlog_out(DEBUG3) << "localPhi = " << localPhi * (360.0 / twopi) << endl; 
          streamlog_out(DEBUG3) << "pad Phi = " << padPhi * (360.0 / twopi) << endl; 
          streamlog_out(DEBUG3) << "pad Phi from track mom = " << ( thisPoint.phi() - trackPhi ) * (360.0 / twopi) << endl; 
          streamlog_out(DEBUG3) << "padTheta = " << padTheta * (360.0 / twopi) << endl; 
          streamlog_out(DEBUG3) << "padTheta from track mom = " << mom.theta() * (360.0 / twopi) << endl; 
          
        }
#endif        
        
      }
      
      //      int pad = padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi());
      int layerNumber = _SimTHit->getCellID();

      if(_rejectCellID0 && (layerNumber<1)) {
        continue;
      }

      de_dx = _SimTHit->getdEdx();

      //       cout << "x position for this hit is " << pos[0] << " - " << _SimTHit->getPosition()[0] << endl; 
      //       cout << "y position for this hit is " << pos[1] << " - " << _SimTHit->getPosition()[1] << endl; 
      //       cout << "z position for this hit is " << pos[2] << " - " << _SimTHit->getPosition()[2] << endl; 

      //       cout << "de/dx for this hit is " << de_dx << endl; 
      //       cout << "MCParticle PID for this hit is " << mcp->getPDG() << endl; 
      //       cout << "x =  " << x << endl; 
      
      // Calculate Point Resolutions according to Ron's Formula 

      // sigma_{RPhi}^2 = sigma_0^2 + Cd^2/N_{eff} * L_{drift}

      // sigma_0^2 = (50micron)^2 + (900micron*sin(phi))^2
      // Cd^2/N_{eff}} = 25^2/(22/sin(theta)*h/6mm)
      // Cd = 25 ( microns / cm^(1/2) )
      // (this is for B=4T, h is the pad height = pad-row pitch in mm,
      // theta is the polar angle)       

      // sigma_{z}^2 = (400microns)^2 + L_{drift}cm * (80micron/sqrt(cm))^2 
      
      double aReso =_pointResoRPhi0*_pointResoRPhi0 + (_pointResoPadPhi*_pointResoPadPhi * sin(padPhi)*sin(padPhi)) ;
      double driftLength = gearTPC.getMaxDriftLength() - (fabs(thisPoint.z()) - _cathode);

      if (driftLength <0) { 
        streamlog_out(DEBUG3) << " TPCDigiProcessor : Warning! driftLength < 0 " << driftLength << " --> Check out your GEAR file!!!!" << std::endl; 
        streamlog_out(DEBUG3) << "Setting driftLength to 0.1" << std::endl;
        streamlog_out(DEBUG3) << "gearTPC.getMaxDriftLength() = " << gearTPC.getMaxDriftLength() << std::endl; 
        driftLength = 0.10;
      }

      double padheight = padLayout.getPadHeight(padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi()));

      double bReso = ( (_diffRPhi * _diffRPhi) / _nEff ) * sin(padTheta) * ( 6.0 / (padheight) );

      double tpcRPhiRes = sqrt( aReso + bReso * (driftLength / 10.0) ); // driftLength in cm

      double tpcZRes  = sqrt(( _pointResoZ0 * _pointResoZ0 ) 
                             + 
                             ( _diffZ * _diffZ ) * (driftLength / 10.0) ); // driftLength in cm 

//             std::cout << "aReso = " << aReso << std::endl;
//             std::cout << "bReso = " << bReso << std::endl;
//             std::cout << "_pointResoPadPhi = " <<_pointResoPadPhi << std::endl;
//             std::cout << "_pointResoRPhi0 = " <<_pointResoRPhi0 << std::endl;
//             std::cout << "_pointResoZ0 = " << _pointResoZ0 << std::endl;
//             std::cout << "_diffRPhi = " << _diffRPhi << std::endl;
//             std::cout << "tpcRPhiRes = " << tpcRPhiRes << std::endl;
//             std::cout << "tpcZRes = " << tpcZRes << std::endl;


      int padIndex = padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi());

      const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;
      double TPCPadPlaneRMin = planeExt[0] ;
      double TPCPadPlaneRMax = planeExt[1] ;

      int iRowHit = padLayout.getRowNumber(padIndex);
      int iPhiHit = padLayout.getPadNumber(padIndex);
      int NBinsZ =  (int) ((2.0 * gearTPC.getMaxDriftLength()) / _binningZ);
      int iZHit = (int) ( (float) NBinsZ * ( gearTPC.getMaxDriftLength() + thisPoint.z() ) / ( 2.0 * gearTPC.getMaxDriftLength() ) ) ;

      if(iZHit<0) iZHit=0;
      if(iZHit>NBinsZ) iZHit=NBinsZ;

      // make sure that the hit lies at the middle of the pad ring
      thisPoint.setPerp(padLayout.getPadCenter(padIndex)[0]);

      if( (thisPoint.perp() < TPCPadPlaneRMin) || (thisPoint.perp() > TPCPadPlaneRMax) ) {
        streamlog_out(DEBUG3) << "Hit R not in TPC " << endl;
        streamlog_out(DEBUG3) << "R = " << thisPoint.perp() << endl; 
        streamlog_out(DEBUG3) << "the tpc InnerRadius = " << TPCPadPlaneRMin << endl;
        streamlog_out(DEBUG3) << "the tpc OuterRadius = " << TPCPadPlaneRMax << endl;
        streamlog_out(DEBUG3) << "Hit Dropped " << endl;
        continue;
      }

      if( (fabs(thisPoint.z()) > gearTPC.getMaxDriftLength()) ) {
        streamlog_out(DEBUG3) << "Hit Z not in TPC " << endl;
        streamlog_out(DEBUG3) << "Z = " << thisPoint.z() << endl; 
        streamlog_out(DEBUG3) << "the tpc Max Z = " << gearTPC.getMaxDriftLength() << endl;
        streamlog_out(DEBUG3) << "Hit Dropped " << endl;
        continue; 
      }
      
      //      Voxel_tpc * atpcVoxel = new Voxel_tpc(iRowHit,iPhiHit,iZHit, pos,  posRPhi, de_dx, tpcRPhiRes, tpcZRes);
      Voxel_tpc * atpcVoxel = new Voxel_tpc(iRowHit,iPhiHit,iZHit, thisPoint, de_dx, tpcRPhiRes, tpcZRes);
      
      _tpcRowHits.at(iRowHit).push_back(atpcVoxel);
      ++numberOfVoxelsCreated;
      
      _tpcHitMap[atpcVoxel] = _SimTHit; 

      //      cout << "a voxel hit "<< i << " has been added to row " << iRowHit << endl;  

      _nMinus2MCP = _previousMCP;
      _previousMCP = _mcp ;

      _nMinus2SimHit = _previousSimTHit;
      _previousSimTHit = _SimTHit;


    }
  }
  
  LCCollection* STHcolLowPt = 0 ;
  try{
    STHcolLowPt = evt->getCollection( _lowPtHitscolName ) ;
  }
  catch(DataNotAvailableException &e){
  }

  if(STHcolLowPt!=NULL){

    int n_sim_hitsLowPt = STHcolLowPt->getNumberOfElements()  ;

    _NBackgroundSimTPCHits += n_sim_hitsLowPt;
    _NSimTPCHits += n_sim_hitsLowPt;

    _padWidth = padLayout.getPadWidth(0)*padCoord[0];

    streamlog_out(DEBUG4) << "number of LowPt hits:" << n_sim_hitsLowPt << std::endl;

    for(int i=0; i< n_sim_hitsLowPt; i++){  

      _SimTHit = dynamic_cast<SimTrackerHit*>( STHcolLowPt->getElementAt( i ) ) ;
      
      CLHEP::Hep3Vector thisPoint(_SimTHit->getPosition()[0],_SimTHit->getPosition()[1],_SimTHit->getPosition()[2]);

      const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;
      double TPCPadPlaneRMin = planeExt[0] ;
      double TPCPadPlaneRMax = planeExt[1] ;

      int NBinsZ =  (int) ((2.0 * gearTPC.getMaxDriftLength()) / _binningZ);

      if( (thisPoint.perp() < TPCPadPlaneRMin) || (thisPoint.perp() > TPCPadPlaneRMax) ) {
        streamlog_out(DEBUG3) << "Hit R not in TPC " << endl;
        streamlog_out(DEBUG3) << "R = " << thisPoint.perp() << endl; 
        streamlog_out(DEBUG3) << "the tpc InnerRadius = " << TPCPadPlaneRMin << endl;
        streamlog_out(DEBUG3) << "the tpc OuterRadius = " << TPCPadPlaneRMax << endl;
        streamlog_out(DEBUG3) << "Hit Dropped " << endl;
        continue;
      }

      if( (fabs(thisPoint.z()) > gearTPC.getMaxDriftLength()) ) {
        streamlog_out(DEBUG3) << "Hit Z not in TPC " << endl;
        streamlog_out(DEBUG3) << "Z = " << thisPoint.z() << endl; 
        streamlog_out(DEBUG3) << "the tpc Max Z = " << gearTPC.getMaxDriftLength() << endl;
        streamlog_out(DEBUG3) << "Hit Dropped " << endl;
        continue; 
      }

      int padIndex = padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi());
      double padheight = padLayout.getPadHeight(padIndex);

      int iRowHit = padLayout.getRowNumber(padIndex);
      int iPhiHit = padLayout.getPadNumber(padIndex);
      int iZHit = (int) ( (float) NBinsZ * 
                          ( gearTPC.getMaxDriftLength() + thisPoint.z() ) / ( 2.0 * gearTPC.getMaxDriftLength() ) ) ;

      // shift the hit in r-phi to the nearest pad-row centre 
      thisPoint.setPerp(padLayout.getPadCenter(padIndex)[0]);

      double de_dx = _SimTHit->getdEdx();

      double tpcRPhiRes = _padWidth;
      double tpcZRes = _binningZ;

      Voxel_tpc * atpcVoxel = new Voxel_tpc(iRowHit,iPhiHit,iZHit, thisPoint, _SimTHit->getdEdx(), tpcRPhiRes, tpcZRes);
      
      _tpcRowHits.at(iRowHit).push_back(atpcVoxel);
      ++numberOfVoxelsCreated;      

      _tpcHitMap[atpcVoxel] = _SimTHit; 
     
    }
  }

  int number_of_adjacent_hits(0);

  streamlog_out(DEBUG4) << "finished looping over simhits, number of voxels = " << numberOfVoxelsCreated << endl;

  int numberOfhitsTreated(0);
    
  vector <Voxel_tpc *> row_hits;
    
  //      cout << "get the row hits" << endl;
    
  for (int i = 0; i<padLayout.getNRows(); ++i){
    
    row_hits = _tpcRowHits.at(i);
    
    //    cout << "got the row hits" << endl;

    //    std::sort(row_hits.begin(), row_hits.end(), compare_z );

    std::sort(row_hits.begin(), row_hits.end(), compare_phi );
    
    //    cout << "row = " << i << "  row_hits.size() = " << row_hits.size() << endl;
    
    for (unsigned int j = 0; j<row_hits.size(); ++j){


      ++numberOfhitsTreated;      
      //      cout << "got the row hit j = " << j << "  row_hits.size() = " << row_hits.size() << endl;
      
      for (unsigned int k = j+1; k<row_hits.size(); ++k){
        
        //        cout << "got the row hits k = " << k << endl;

//        // SJA:FIXME: the z separation here has to be dependant on theta, for highly dipped tracks the z-region affected will be considerably larger
//        if(row_hits[k]->getZIndex() > row_hits[k]->getZIndex()+1){ // this should be > row_hits[k]->getZIndex+1+(0.5*steplength*(costheta)/zbinsize)
//          break; // only compare hits in adjacent z bins
//        }

//        else if(row_hits[k]->getZIndex()==row_hits[j]->getZIndex() || fabs(row_hits[k]->getZ() - row_hits[j]->getZ()) < _doubleHitResZ ) {

        if(row_hits[k]->getPhiIndex() > (row_hits[j]->getPhiIndex())+2){ // here we need an AND to catch the wrap around 
          //cout <<  row_hits[j]->getPhiIndex()+2 << " " << row_hits[k]->getPhiIndex() << endl;
          break; // only compare hits in adjacent phi bins
        }
        
        else if( row_hits[k]->getPhiIndex()==row_hits[j]->getPhiIndex() 
                 || 
                 ( (fabs(row_hits[k]->getHep3Vector().deltaPhi(row_hits[j]->getHep3Vector()))) * row_hits[j]->getR()) < _doubleHitResRPhi ) {
          
          //          // compare phi
          //          double dPhi = fabs(row_hits[k]->getHep3Vector().deltaPhi(row_hits[j]->getHep3Vector()));
          //
          //          if( dPhi*row_hits[j]->getR() < _doubleHitResRPhi ){
          
          // compare z

          map <Voxel_tpc*,SimTrackerHit*> ::iterator it;
          
          SimTrackerHit* Hit1 = NULL;
          SimTrackerHit* Hit2 = NULL;
          
          it=_tpcHitMap.find(row_hits[j]);
          
          if(it!= _tpcHitMap.end()) {
            Hit1 = _tpcHitMap[row_hits[j]];
          }

          it=_tpcHitMap.find(row_hits[k]);

          if(it!= _tpcHitMap.end()) {
            Hit2 = _tpcHitMap[row_hits[k]];
          }
          
          double pathlengthZ1(0.0);
          double pathlengthZ2(0.0);
          
          if( Hit1 && Hit2 ){
            
            LCFlagImpl colFlag( STHcol->getFlag() );
            
            if( Hit1->getPathLength()!=0.0 && Hit2->getPathLength()!=0.0 && colFlag.bitSet(LCIO::THBIT_MOMENTUM)){
              
              const float * Momentum1 = Hit1->getMomentum() ;
              const float * Momentum2 = Hit2->getMomentum() ;
              
              CLHEP::Hep3Vector mom1(Momentum1[0],Momentum1[1],Momentum1[2]);
              CLHEP::Hep3Vector mom2(Momentum2[0],Momentum2[1],Momentum2[2]);
              
              pathlengthZ1 = fabs(Hit1->getPathLength() * (mom1.z() / mom1.perp()));
              pathlengthZ2 = fabs(Hit2->getPathLength() * (mom2.z() / mom2.perp()));
              
            } 
            else {
              pathlengthZ1 = (5.0) ; 
              pathlengthZ2 = (5.0) ; 
            }
            
            double dZ(0.0);
            double dZAlt(0.0);
            
            double absZ1 = fabs(row_hits[j]->getZ()); 
            double absZ2 = fabs(row_hits[k]->getZ());
            
            //          if( absZ1 < absZ2 ) {
            //            dZAlt = ( (absZ2 - 0.5*(5.0 + _binningZ)) - (absZ1 + 0.5*(5.0 + _binningZ)) );
            //          } else {
            //            dZAlt = ( (absZ1 - 0.5*(5.0 + _binningZ)) - (absZ2 + 0.5*(5.0 + _binningZ)) );
            //          }
            
            
            
            //          if(pathlengthZ1!=5.0 && pathlengthZ2 !=5.0 )
            //          cout << " pathlengthZ1 " << pathlengthZ1  <<  " pathlengthZ2 " << pathlengthZ2 << endl;
            //
            
            if ( absZ2 > absZ1 ) {
              dZAlt = absZ2 - absZ1;
            } else {
              dZAlt = absZ1 - absZ2;
            }
            
            //          
          
            //          if( dZ < _doubleHitResZ ){
           
            if( absZ2 > absZ1 ) {
              dZ = ( (absZ2 - 0.5*(pathlengthZ2 + _binningZ)) - (absZ1 + 0.5*(pathlengthZ1 + _binningZ)) );
            } else {
              dZ = ( (absZ1 - 0.5*(pathlengthZ1 + _binningZ)) - (absZ2 + 0.5*(pathlengthZ2 + _binningZ)) );
            }
            
        
            if( dZAlt < dZ ) {
              streamlog_out(DEBUG3) << "dZAlt " << dZAlt << endl;
              streamlog_out(DEBUG3) << "dZ " << dZ << endl;
              streamlog_out(DEBUG3) << "row_hits[j]->getZ() " << row_hits[j]->getZ() << endl;
              streamlog_out(DEBUG3) << "row_hits[k]->getZ() " << row_hits[k]->getZ() << endl;
              streamlog_out(DEBUG3) << "row_hits[j]->getPhi()*R " << row_hits[j]->getPhi()*row_hits[j]->getR() << endl;
              streamlog_out(DEBUG3) << "row_hits[k]->getPhi()*R " << row_hits[k]->getPhi()*row_hits[k]->getR() << endl;  
              streamlog_out(DEBUG3) << "pathlength 1 " << Hit1->getPathLength() << endl;  
              streamlog_out(DEBUG3) << "pathlength 2 " << Hit2->getPathLength() << endl;  
              streamlog_out(DEBUG3) << "pathlengthZ1 " << pathlengthZ1 << endl;  
              streamlog_out(DEBUG3) << "pathlengthZ2 " << pathlengthZ2 << endl;  
            }
            
            //if( dZAlt < _doubleHitResZ ){            
              if( dZ < _doubleHitResZ ){            
              
              
              //          if( true ){
              
              //              cout << "*&^*&*&*&*&*&*&*&**&*&*&*&*&&*&*&**&*&*&" << endl;
              //              cout << "double hit candidate found in row: " << i <<  endl;
              //              
              //
              //
              //              cout << "dZ " << dZ << endl;
              //              cout << "row_hits[j]->getZ() " << row_hits[j]->getZ() << endl;
              //              cout << "row_hits[k]->getZ() " << row_hits[k]->getZ() << endl;
              //              cout << "row_hits[j]->getPhi()*R " << row_hits[j]->getPhi()*row_hits[j]->getR() << endl;
              //              cout << "row_hits[k]->getPhi()*R " << row_hits[k]->getPhi()*row_hits[k]->getR() << endl;  
              //              cout << "pathlength 1 " << Hit1->getPathLength() << endl;  
              //              cout << "pathlength 2 " << Hit2->getPathLength() << endl;  
              //              cout << "pathlengthZ1 " << pathlengthZ1 << endl;  
              //              cout << "pathlengthZ2 " << pathlengthZ2 << endl;  
              //            cout << "costheta1 " << (mom1->z() / mom1->r()) << endl;  
              //            cout << "costheta2 " << (mom2->z() / mom2->r()) << endl;  
              //
              //            
              //            MCParticle* mcp1 = Hit1->getMCParticle();
              //            MCParticle* mcp2 = Hit2->getMCParticle();
              //
              //            cout << "double hit found in row: " << i  << "  "
              //                 << "P1 : " <<  mom1->r() << "     " 
              //                 << "P2 : " <<  mom2->r() << "     " ;
              //            
              //
              //            if( mcp1 && mcp2 ){
              //              cout << "    MCP1 : " << Hit1->getMCParticle()->id() 
              //                   << " (" <<  Hit1->getMCParticle()->getPDG() << ")   " 
              //                   << "MCP2 : " << Hit2->getMCParticle()->id() 
              //                   << " (" <<  Hit2->getMCParticle()->getPDG() << ")   " 
              //                   << endl;
              //            }
              //
              //           }
              row_hits[j]->setAdjacent(row_hits[k]);
              row_hits[k]->setAdjacent(row_hits[j]);
              ++number_of_adjacent_hits;
            }
          } else {
            streamlog_out(DEBUG3) << "Hit1=" << Hit1 << "Hit2=" << Hit2 << endl; 
          }
        }
      }
    }


    // now all hits have been checked for adjacent hits, go throught and write out the hits or merge

    for (unsigned int j = 0; j<row_hits.size(); ++j){

      Voxel_tpc* seed_hit = row_hits[j];

      if(seed_hit->IsMerged() || seed_hit->IsClusterHit()) { 
        continue;
      }

      if(seed_hit->getNumberOfAdjacent()==0){ // no adjacent hits so smear and write to hit collection
        writeVoxelToHit(seed_hit);        
      }

      else if(seed_hit->getNumberOfAdjacent() < (_maxMerge)){ // potential 3-hit cluster, can use simple average merge. 

        vector <Voxel_tpc*>* hitsToMerge = new vector <Voxel_tpc*>;

        int clusterSize = seed_hit->clusterFind(hitsToMerge);
        
        if( clusterSize <= _maxMerge ){ // merge cluster
          seed_hit->setIsMerged();
          //          std::cout << "merge cluster with " <<  clusterSize << " hits" << " hitsToMerge has size " << hitsToMerge->size()<< std::endl;
          writeMergedVoxelsToHit(hitsToMerge);  
        }
        delete hitsToMerge;
      } 
    } 
  }

  int numberOfHits(0);

  for (int i = 0; i<padLayout.getNRows(); ++i){
    row_hits = _tpcRowHits.at(i);
    for (unsigned int j = 0; j<row_hits.size(); ++j){
      numberOfHits++;
      Voxel_tpc* seed_hit = row_hits[j];
      if(seed_hit->IsMerged() || seed_hit->IsClusterHit() || seed_hit->getNumberOfAdjacent() > _maxMerge ) { 
        ++_NRevomedHits;
        if( (_tpcHitMap[ seed_hit ])->getMCParticle()!=NULL ) { 
          ++_NLostPhysicsTPCHits;
          const float *mom= (_tpcHitMap[ seed_hit ])->getMomentum() ;
          double ptSQRD = mom[0]*mom[0]+mom[1]*mom[1] ; 
          if( ptSQRD > (0.2*0.2) ) ++_NLostPhysicsAbove02GeVPtTPCHits ;
          if( ptSQRD > 1.0 )  ++_NLostPhysicsAbove1GeVPtTPCHits ;
        }
      }
    }
  }

  streamlog_out(DEBUG4) << "the number of adjacent hits is " <<  number_of_adjacent_hits << "  _doubleHitResZ " << _doubleHitResZ << endl;  

  streamlog_out(DEBUG4) << "number of rec_hits = "  << _NRecTPCHits << endl ;
  streamlog_out(DEBUG4) << "finished row hits " << numberOfHits << " " << numberOfhitsTreated << endl;    
  
  // set the parameters to decode the type information in the collection
  // for the time being this has to be done manually
  // in the future we should provide a more convenient mechanism to 
  // decode this sort of meta information

  StringVec typeNames ;
  IntVec typeValues ;
  typeNames.push_back( LCIO::TPCHIT ) ;
  typeValues.push_back( 1 ) ;
  _trkhitVec->parameters().setValues("TrackerHitTypeNames" , typeNames ) ;
  _trkhitVec->parameters().setValues("TrackerHitTypeValues" , typeValues ) ;
  
  evt->addCollection( _trkhitVec , _TPCTrackerHitsCol ) ;
  
  // delete voxels
  for (int i = 0; i<padLayout.getNRows(); ++i){
    vector <Voxel_tpc *>* current_row = &_tpcRowHits[i];  
    for (unsigned int j = 0; j<current_row->size(); ++j){
      delete current_row->at(j);
    }
  }

#ifdef DIGIPLOTS
  _NSimTPCHitsHisto->fill(_NSimTPCHits);
  _NBackgroundSimTPCHitsHisto->fill(_NBackgroundSimTPCHits);
  _NPhysicsSimTPCHitsHisto->fill(_NPhysicsSimTPCHits);
  _NPhysicsAbove02GeVSimTPCHitsHisto->fill(_NPhysicsAbove02GeVSimTPCHits);
  _NPhysicsAbove1GeVSimTPCHitsHisto->fill(_NPhysicsAbove1GeVSimTPCHits);
  _NRecTPCHitsHisto->fill(_NRecTPCHits);
  
  _NLostPhysicsTPCHitsHisto->fill(_NLostPhysicsTPCHits);
  _NLostPhysicsAbove02GeVPtTPCHitsHisto->fill(_NLostPhysicsAbove02GeVPtTPCHits);
  _NLostPhysicsAbove1GeVPtTPCHitsHisto->fill(_NLostPhysicsAbove1GeVPtTPCHits);
  _NRevomedHitsHisto->fill(_NRevomedHits);
  
  _NKeptPhysicsTPCHitsHistoPercent->fill( (float)(_NPhysicsSimTPCHits-_NLostPhysicsTPCHits) / (float)_NPhysicsSimTPCHits );
  _NKeptPhysicsAbove02GeVPtTPCHitsHistoPercent->fill( (float)(_NPhysicsAbove02GeVSimTPCHits-_NLostPhysicsAbove02GeVPtTPCHits) / (float)_NPhysicsAbove02GeVSimTPCHits );
  _NKeptPhysicsAbove1GeVPtTPCHitsHistoPercent->fill( (float)(_NPhysicsAbove1GeVSimTPCHits-_NLostPhysicsAbove1GeVPtTPCHits) / (float)_NPhysicsAbove1GeVSimTPCHits );
#endif

  streamlog_out(DEBUG4) << "_NSimTPCHits = " << _NSimTPCHits << endl;
  streamlog_out(DEBUG4) << "_NBackgroundSimTPCHits = " << _NBackgroundSimTPCHits << endl;
  streamlog_out(DEBUG4) << "_NPhysicsSimTPCHits = " << _NPhysicsSimTPCHits << endl;
  streamlog_out(DEBUG4) << "_NPhysicsAbove02GeVSimTPCHits = " << _NPhysicsAbove02GeVSimTPCHits << endl;
  streamlog_out(DEBUG4) << "_NPhysicsAbove1GeVSimTPCHits = " << _NPhysicsAbove1GeVSimTPCHits << endl;
  streamlog_out(DEBUG4) << "_NRecTPCHits = " << _NRecTPCHits<< endl;
  streamlog_out(DEBUG4) << "_NLostPhysicsTPCHits = " << _NLostPhysicsTPCHits << endl;
  streamlog_out(DEBUG4) << "_NLostPhysicsAbove02GeVPtTPCHits = " << _NLostPhysicsAbove02GeVPtTPCHits << endl;
  streamlog_out(DEBUG4) << "_NLostPhysicsAbove1GeVPtTPCHits = " << _NLostPhysicsAbove1GeVPtTPCHits << endl;
  streamlog_out(DEBUG4) << "_NRevomedHits = " << _NRevomedHits << endl;
  
 _nEvt++;  
}



void TPCDigiProcessor::check( LCEvent * evt ) 
{ 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TPCDigiProcessor::end()
{ 

  //#ifdef EXPERTCHECKPLOTS
#ifdef DIGIPLOTS
  _TREE->commit();
  _TREE->cd("/Histograms");
  _TREE->ls("..");

  _TREE->close();  
  cout << "EXPERTCHECKPLOTS Finished" << endl;
#endif

  gsl_rng_free(_random);
  cout << "TPCDigiProcessor::end()  " << name() 
       << " processed " << _nEvt << " events in " << _nRun << " runs "
       << endl ;
  //  
}

void TPCDigiProcessor::writeVoxelToHit( Voxel_tpc* aVoxel){

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  const gear::Vector2D padCoord = padLayout.getPadCenter(1) ;
  
  Voxel_tpc* seed_hit  = aVoxel;
    //store hit variables
  TrackerHitImpl* trkHit = new TrackerHitImpl ;
  //now the hit pos has to be smeared

  double tpcRPhiRes = seed_hit->getRPhiRes();
  double tpcZRes = seed_hit->getZRes();

  CLHEP::Hep3Vector* point = new CLHEP::Hep3Vector(seed_hit->getX(),seed_hit->getY(),seed_hit->getZ());

  double unsmearedPhi = point->phi();

  double randrp = gsl_ran_gaussian(_random,tpcRPhiRes);
  double randz =  gsl_ran_gaussian(_random,tpcZRes);

  point->setPhi( point->phi() + randrp/ point->perp() );
  point->setZ( point->z() + randz );

  // make sure the hit is not smeared beyond the TPC Max DriftLength
  if( fabs(point->z()) > gearTPC.getMaxDriftLength() ) point->setZ( (fabs(point->z()) / point->z() ) * gearTPC.getMaxDriftLength() );

  double pos[3] = {seed_hit->getX(),seed_hit->getY(),seed_hit->getZ()}; 
  trkHit->setPosition(pos);
  trkHit->setdEdx(seed_hit->getdEdx());
  trkHit->setType( 500 );
                
  //          cout << "row_hits->getY() = " << row_hits[j]->getY() << "  row_hits->getY() = " << row_hits[j]->getX() ;
  //          cout << "  phi = " <<  phi ;
  //          cout << "  tpcRPhiRes = " <<  tpcRPhiRes;
  //          cout << "  cos(phi)*cos(phi)*tpcRPhiRes*tpcRPhiRes = " << cos(phi)*cos(phi)*tpcRPhiRes*tpcRPhiRes ;
  //          cout << endl;
        
  // For no error in R
  float covMat[TRKHITNCOVMATRIX]={sin(unsmearedPhi)*sin(unsmearedPhi)*tpcRPhiRes*tpcRPhiRes,
                                  -cos(unsmearedPhi)*sin(unsmearedPhi)*tpcRPhiRes*tpcRPhiRes,
                                  cos(unsmearedPhi)*cos(unsmearedPhi)*tpcRPhiRes*tpcRPhiRes,
                                  0.,
                                  0.,
                                  float(tpcZRes*tpcZRes)};
        
  trkHit->setCovMatrix(covMat);      

  if( _tpcHitMap[seed_hit] == NULL ){
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: TPCDigiProcessor \n" 
             << "SimTracker Pointer is NULL throwing exception\n"
             << "\n" ;
    throw Exception(errorMsg.str());
  }

  if(pos[0]*pos[0]+pos[1]*pos[1]>0.0){ 
    // 	  push back the SimTHit for this TrackerHit
    trkHit->rawHits().push_back( _tpcHitMap[seed_hit] );                        
    _trkhitVec->addElement( trkHit ); 
    _NRecTPCHits++;
  }
        
#ifdef EXPERTCHECKPLOTS
  SimTrackerHit* theSimHit = _tpcHitMap[seed_hit];
  double rSimSqrd = theSimHit->getPosition()[0]*theSimHit->getPosition()[0] + theSimHit->getPosition()[1]*theSimHit->getPosition()[1];

  double phiSim = atan2(theSimHit->getPosition()[1],theSimHit->getPosition()[0]);

  double rPhiDiff = (point->phi() - phiSim)*sqrt(rSimSqrd);
  double rPhiPull = ((point->phi() - phiSim)*sqrt(rSimSqrd))/(sqrt((covMat[2])/(cos(point->phi())*cos(point->phi()))));
        
  double zDiff = point->getZ() - theSimHit->getPosition()[2];
  double zPull = zDiff/sqrt(covMat[5]);
        
        
  _rPhiDiffHisto->fill(rPhiDiff);
  _rPhiPullHisto->fill(rPhiPull);
  _phiDistHisto->fill(point->phi() - phiSim);
  _zDiffHisto->fill(zDiff);
  _zPullHisto->fill(zPull);
        
  _zSigmaVsZHisto->fill(seed_hit->getZ(),sqrt(covMat[5]));
  _rPhiSigmaHisto->fill(sqrt((covMat[2])/(cos(point->phi())*cos(point->phi()))));
  _zSigmaHisto->fill(sqrt(covMat[5]));
#endif
}

void TPCDigiProcessor::writeMergedVoxelsToHit( vector <Voxel_tpc*>* hitsToMerge){

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  const gear::Vector2D padCoord = padLayout.getPadCenter(1) ;

  TrackerHitImpl* trkHit = new TrackerHitImpl ;
  
  double sumZ = 0;
  double sumPhi = 0;
  double sumdEdx = 0;
  double R = 0;
  double lastR = 0;

  for(int ihitCluster = 0; ihitCluster < hitsToMerge->size(); ++ihitCluster){
    
    if( ihitCluster > 0 && (lastR - hitsToMerge->at(ihitCluster)->getR()) > 1e-10 ){
      //     cout << "Hits not at the same R" << endl;
      //     cout << "hit1 R = " << hitsToMerge->at(ihitCluster)->getR() << " hit2 R = " << lastR << "Rdiff = " << hitsToMerge->at(ihitCluster)->getR() - lastR << endl;
    }
    
    sumZ += hitsToMerge->at(ihitCluster)->getZ();
    sumPhi += hitsToMerge->at(ihitCluster)->getPhi();
    sumdEdx += hitsToMerge->at(ihitCluster)->getdEdx();
    hitsToMerge->at(ihitCluster)->setIsMerged();
//    cout << "hit"<< ihitCluster+2 << " Z = " << hitsToMerge->at(ihitCluster)->getZ() << endl;
//    cout << "hit"<< ihitCluster+2 << " RPhi = " << hitsToMerge->at(ihitCluster)->getR()*hitsToMerge->at(ihitCluster)->getPhi() << endl;
//    cout << "hit"<< ihitCluster+2 << " Phi = " << hitsToMerge->at(ihitCluster)->getPhi() << endl;
//    cout << "hit"<< ihitCluster+2 << " has " << hitsToMerge->at(ihitCluster)->getNumberOfAdjacent() << " adjacent hits" << endl;
    lastR = hitsToMerge->at(ihitCluster)->getR();
   
    trkHit->rawHits().push_back( _tpcHitMap[ hitsToMerge->at(ihitCluster) ] );                        

  }

  double avgZ = sumZ/(hitsToMerge->size());
  double avgPhi = sumPhi/(hitsToMerge->size());
 
  CLHEP::Hep3Vector* mergedPoint = new CLHEP::Hep3Vector(1.0,1.0,1.0);
  mergedPoint->setPerp(lastR);
  mergedPoint->setPhi(avgPhi);
  mergedPoint->setZ(avgZ);

  //store hit variables
  //now the hit pos has to be smeared
  double pos[3] = {mergedPoint->getX(),mergedPoint->getY(),mergedPoint->getZ()}; 
  trkHit->setPosition(pos);
  trkHit->setdEdx(sumdEdx);
  trkHit->setType( 500 );
        
  double phi = mergedPoint->getPhi();
  double tpcRPhiRes = _padWidth;
  double tpcZRes = _binningZ;
        
  //          cout << "row_hits->getY() = " << row_hits[j]->getY() << "  row_hits->getY() = " << row_hits[j]->getX() ;
  //          cout << "  phi = " <<  phi ;
  //          cout << "  tpcRPhiRes = " <<  tpcRPhiRes;
  //          cout << "  cos(phi)*cos(phi)*tpcRPhiRes*tpcRPhiRes = " << cos(phi)*cos(phi)*tpcRPhiRes*tpcRPhiRes ;
  //          cout << endl;
        
  // For no error in R
  float covMat[TRKHITNCOVMATRIX]={sin(phi)*sin(phi)*tpcRPhiRes*tpcRPhiRes,
                                  -cos(phi)*sin(phi)*tpcRPhiRes*tpcRPhiRes,
                                  cos(phi)*cos(phi)*tpcRPhiRes*tpcRPhiRes,
                                  0.,
                                  0.,
                                  float(tpcZRes*tpcZRes)};
                
  trkHit->setCovMatrix(covMat);      
        
  if(pos[0]*pos[0]+pos[1]*pos[1]>0.0){ 
    // 	  push back the SimTHit for this TrackerHit
    _trkhitVec->addElement( trkHit ); 
    ++_nRechits;
  }
  
  delete mergedPoint;
      
}

#ifdef EXPERTCHECKPLOTS
void TPCDigiProcessor::plotHelixHitResidual( MCParticle *mcp, CLHEP::Hep3Vector *thisPoint){

      const double bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
      const double FCT = 2.99792458E-4;
      double charge = mcp->getCharge();
      const double *mom = mcp->getMomentum();
      double pt = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
      double radius = pt / (FCT*bField);
      double tanLambda = mom[2]/pt;
      double omega = charge / radius;

      if(pt>1.0) {

        //FIXME SJA: this is only valid for tracks from the IP and should be done correctly for non prompt tracks
        double Z0 = 0.;      
        double D0  = 0.;
        
        LCVector3D refPoint(0.,0.,0);
        
        SimpleHelix* helix = new SimpleHelix(D0, 
                                             atan2(mom[1],mom[0]), 
                                             omega, 
                                             Z0, 
                                             tanLambda,
                                             refPoint);
        

        // an almost "infinite" cylinder in z        
        LCVector3D startCylinder(0.,0.,-1000000.0);
        LCVector3D endCylinder(0.,0.,1000000.0);
        bool endplane=true;
        
        LCCylinder cylinder(startCylinder,endCylinder,thisPoint->perp(),endplane);
        
        bool pointExists = false;
        
        double pathlength = helix->getIntersectionWithCylinder( cylinder, pointExists);
        
        LCErrorMatrix* errors = new LCErrorMatrix();
        
        if(pointExists){
     
          LCVector3D intersection = helix->getPosition(pathlength, errors); 
     
          double intersectionPhi = atan2(intersection[1],intersection[0]);
          double residualRPhi = ((intersectionPhi-thisPoint->phi())) ;
          _ResidualsRPhiHisto->fill(residualRPhi);
          
        }
        
        delete errors;
        delete helix;

        const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
        const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
        
        int row = padLayout.getRowNumber(padLayout.getNearestPad(thisPoint->perp(),thisPoint->phi()));
        int pad = padLayout.getNearestPad(thisPoint->perp(),thisPoint->phi());
        
        double rHit_diff = thisPoint->perp()  
          - padLayout.getPlaneExtent()[0]  
          - (( row + 0.5 ) 
             * padLayout.getPadHeight(pad)) ;
        
        _radiusCheckHisto->fill(rHit_diff);
        
//      cout << "$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$" << endl;
//      cout << "thisPoint->perp() = " << thisPoint->perp() << endl;
//      cout << "TPC Sensitive rMin = " << padLayout.getPlaneExtent()[0] << endl;
//      cout << "Row number + 0.5 = " <<  row + 0.5 << endl;
//      cout << "Pad Height = " <<  padLayout.getPadHeight(pad) << endl;
//      cout << "Row Height = " <<   padLayout.getRowHeight(row) << endl;
//      cout << "R_hit - TPC Rmin - ((RowIndex + 0.5 )* padheight) = " << rHit_diff << endl;
          
      }
      return;
}
#endif


double TPCDigiProcessor::getPadPhi( CLHEP::Hep3Vector *thisPoint, CLHEP::Hep3Vector* firstPoint, CLHEP::Hep3Vector* middlePoint, CLHEP::Hep3Vector* lastPoint){

  CLHEP::Hep2Vector firstPointRPhi(firstPoint->x(),firstPoint->y());
  CLHEP::Hep2Vector middlePointRPhi(middlePoint->x(),middlePoint->y());
  CLHEP::Hep2Vector lastPointRPhi(lastPoint->x(),lastPoint->y());

  Circle theCircle(&firstPointRPhi, &middlePointRPhi, &lastPointRPhi);
  
  double localPhi = atan2((thisPoint->y() - theCircle.GetCenter()->y()) ,(thisPoint->x() - theCircle.GetCenter()->x())) + (twopi/4.0) ;
  
  if(localPhi>twopi) localPhi=localPhi - twopi;
  if(localPhi<0.0) localPhi=localPhi + twopi;
  if(localPhi>twopi/2.0) localPhi = localPhi - twopi/2.0 ;
  
  double pointPhi = thisPoint->phi();
  
  if(pointPhi>twopi) pointPhi=pointPhi - twopi;
  if(pointPhi<0.0) pointPhi=pointPhi + twopi;
  if(pointPhi>twopi/2.0) pointPhi = pointPhi - twopi/2.0 ;
  
  double padPhi = fabs(pointPhi - localPhi);
  
  return padPhi;

}

double TPCDigiProcessor::getPadTheta(CLHEP::Hep3Vector* firstPoint, CLHEP::Hep3Vector* middlePoint, CLHEP::Hep3Vector* lastPoint){

  // Calculate thetaPad for current hit
  CLHEP::Hep2Vector firstPointRPhi(firstPoint->x(),firstPoint->y()) ;
  CLHEP::Hep2Vector middlePointRPhi(middlePoint->x(),middlePoint->y());
  CLHEP::Hep2Vector lastPointRPhi(lastPoint->x(),lastPoint->y());

  Circle theCircle(&firstPointRPhi, &middlePointRPhi, &lastPointRPhi);
  
  double deltaPhi = firstPoint->deltaPhi(*lastPoint);

  double pathlength = fabs(deltaPhi) *  theCircle.GetRadius();
  
  double padTheta = atan ( pathlength / fabs(lastPoint->z() - firstPoint->z()) ) ;
                        
  double pathlength1 = 2.0 * asin( ( sqrt (
                                           (middlePointRPhi.x() - firstPointRPhi.x()) * (middlePointRPhi.x()-firstPointRPhi.x())
                                           +
                                           (middlePointRPhi.y()-firstPointRPhi.y()) * (middlePointRPhi.y()-firstPointRPhi.y())
                                           ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius()  ;
  
  
  double pathlength2 = 2.0 * asin( ( sqrt (
                                           (lastPointRPhi.x()-middlePointRPhi.x()) * (lastPointRPhi.x()-middlePointRPhi.x())
                                           +
                                           (lastPointRPhi.y()-middlePointRPhi.y()) * (lastPointRPhi.y()-middlePointRPhi.y())
                                           ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius()  ;
  
  
  padTheta = atan ((fabs(pathlength1 + pathlength2)) / (fabs(lastPoint->z() - firstPoint->z())) ) ;

  
  return padTheta;

}
