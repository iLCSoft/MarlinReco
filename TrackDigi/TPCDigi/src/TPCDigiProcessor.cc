/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

#include "TPCDigiProcessor.h"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include <gsl/gsl_randist.h>
#include "marlin/VerbosityLevels.h"


#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>

#include "Circle.h"
#include"constants.h"

//stl exception handler
#include <stdexcept>
#include "constants.h"
#include "voxel.h"

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
//

using namespace lcio ;
using namespace marlin ;
using namespace constants ;
using namespace std ;

#ifdef MARLIN_USE_AIDA
using namespace AIDA ;
#endif


TPCDigiProcessor aTPCDigiProcessor ;

// put this here as I dont't know how to get it into the tpc class or if that is even the write place for it
bool compare_phi( Voxel_tpc * a, Voxel_tpc * b){return ( a->getPhiIndex() < b->getPhiIndex() );}

TPCDigiProcessor::TPCDigiProcessor() : Processor("TPCDigiProcessor") 
{
  
  // modify processor description
  _description = "Produces TPC TrackerHit collection from SimTrackerHit collection, smeared in RPhi and Z" ;
  
  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "CollectionName" , 
                           "Name of the SimTrackerHit collection"  ,
                           _colName ,
                           std::string("STpc01_TPC") ) ;

  registerOutputCollection( LCIO::TRACKERHIT,
                            "TPCTrackerHitsCol" , 
                            "Name of the digitized TrackerHit collection"  ,
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

  registerProcessorParameter( "PixZ" ,
                              "Defines spatial slice in Z"  ,
                              _pixZ ,
                              (float)1.4) ;

  registerProcessorParameter( "PixRP" ,
                              "Defines spatial slice in RP"  ,
                              _pixRP ,
                              (float)1.0) ;




}


void TPCDigiProcessor::init() 
{ 
  
  // From GNU documentation:
  // A replacement for the standard terminate_handler which prints 
  // more information about the terminating exception (if any) on stderr. Call ...
  std::set_terminate (__gnu_cxx::__verbose_terminate_handler);

#ifdef EXPERTCHECKPLOTS
  /// Hook an AIDA implementation -----------------------------------------------

  // First create a pointer to the "IAnalysisFactory" of a specific AIDA
  // implementation. This factory can then be used to produce all other 
  // factories.
  AF = AIDA_createAnalysisFactory();
  
  // Create a ITreeFactory. -----------------------------------------------------
  // A ITree can be used to store AIDA objects in memory or on disk.
  
  TRF = AF->createTreeFactory();
  
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

  TREE = TRF->create("TPCDigi.root",
                     "root",
                     false,
                     true);

  /// Create an IHistogramFactory which is bound to the tree "*TREE". -----------

  /*
   * Create an IHistogramFactory.
   * @param tree The ITree which created histograms will be associated to.
   * @return     The IHistogramFactory.
   */
  // IHistogramFactory * IAnalysisFactory::createHistogramFactory(ITree & tree);

  HF = AF->createHistogramFactory(*TREE);
  
  TREE->mkdir("Histograms");
 
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
  


  phiDiffHisto = HF->createHistogram1D("Histograms/phi_diff",
                                       "Calculated Phi - Track Phi",
                                       201, -0.05, 0.05);

  thetaDiffHisto = HF->createHistogram1D("Histograms/theta_diff",
                                         "Calculated Theta - Track Theta",
                                         201, -0.05, 0.05);
 
  phiRelHisto = HF->createHistogram1D("Histograms/padPhi",
                                      "Phi Relative to the Pad",
                                      201, 0.0, 6.3);

  thetaRelHisto = HF->createHistogram1D("Histograms/padtheta",
                                        "Theta Relative to the pad",
                                        201, 0.0, 6.3);

  rDiffHisto = HF->createHistogram1D("Histograms/rDiff",
                                     "R_rec - R_sim",
                                     201, -1.0, 1.0);

  phiDistHisto = HF->createHistogram1D("Histograms/phiDist",
                                       "phi_rec - Phi_sim",
                                       201, -1.0, 1.0);

  phiPullHisto = HF->createHistogram1D("Histograms/phiPull",
                                       "(Phi_rec - Phi_sim) / Sigma_phi",
                                       201, -10.0, 10.0);


#endif  


  // usually a good idea to
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

  static bool firstEvent = true ;
  int rechits = 0;
  
  // this gets called for every event 
  // usually the working horse ...

  //     int skipToEvent = 0;
  //     if(_nEvt<skipToEvent) {
  //       cout << "skipping event " << _nEvt << endl;
  //      ++_nEvt;
  //       return;
  //     }

  
  if(firstEvent==true) cout << "TPCDigiProcessor called for first event" << endl;

  firstEvent = false ;

  LCCollection* STHcol = 0 ;
  try{
    STHcol = evt->getCollection( _colName ) ;
  }
  catch(DataNotAvailableException &e){
  }
  
  if( STHcol != 0 ){
  
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT )  ;
    
    int n_sim_hits = STHcol->getNumberOfElements()  ;

    streamlog_out(DEBUG) << "number of SimHits = " << n_sim_hits << std::endl;

    // Assume initialy that there is no merging 
    
    //
    const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
    const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

    vector< vector <Voxel_tpc *> > tpcRowHits;

    // set size of row_hits to hold (n_rows) vectors
    tpcRowHits.resize(padLayout.getNRows());

    map< Voxel_tpc *,SimTrackerHit *> tpcHitMap;
  
 
 
    EVENT::MCParticle *mcp(NULL);
    EVENT::MCParticle *previousMCP(NULL);
    EVENT::MCParticle *nextMCP(NULL);

    SimTrackerHit* SimTHit(NULL);
    SimTrackerHit* previousSimTHit(NULL);
    SimTrackerHit* nextSimTHit(NULL);


    for(int i=0; i< n_sim_hits; i++){
      
      SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      mcp = SimTHit->getMCParticle() ; 
      
      nextMCP = NULL;
      if (i<(n_sim_hits-1)) {
        nextSimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i+1 ) ) ;
        nextMCP= nextSimTHit->getMCParticle() ;
      }

      float de_dx;
 
      double padPhi;
      double padTheta = 0.0;

      CLHEP::Hep2Vector *precedingPoint = (NULL);
      CLHEP::Hep2Vector *thisPoint = (NULL);
      CLHEP::Hep2Vector *followingPoint = (NULL);
      CLHEP::Hep2Vector *nMinus2Point = (NULL);
      CLHEP::Hep2Vector *nPlus2Point = (NULL);

      thisPoint = new CLHEP::Hep2Vector(SimTHit->getPosition()[0],SimTHit->getPosition()[1]);

      // Calculate difference in Phi for current hit with that of the pad
      
      if( mcp==previousMCP && mcp==nextMCP) { // still with on the same track

        precedingPoint = new CLHEP::Hep2Vector(previousSimTHit->getPosition()[0],previousSimTHit->getPosition()[1]);
        followingPoint = new CLHEP::Hep2Vector(nextSimTHit->getPosition()[0],nextSimTHit->getPosition()[1]);

        Circle theCircle(precedingPoint, thisPoint, followingPoint);

        //        cout << "the circle has radius = " << theCircle.GetRadius() << endl;
        //        cout << "the circle has center = " << theCircle.GetCenter()->x() << "  " << theCircle.GetCenter()->y() << endl;

        // cout << "point x = " << SimTHit->getPosition()[0] << " point y = " << SimTHit->getPosition()[1] << endl; 

        double localPhi = 
          atan2((thisPoint->y() - theCircle.GetCenter()->y()) ,(thisPoint->x() - theCircle.GetCenter()->x())) + (twopi/4.0) ;
        //          atan2((SimTHit->getPosition()[1] - theCircle.GetCenter()->y()) ,(SimTHit->getPosition()[0] - theCircle.GetCenter()->x())) + (twopi/4.0) ;


        if(localPhi>twopi) localPhi=localPhi - twopi;
        if(localPhi<0.0) localPhi=localPhi + twopi;
        if(localPhi>twopi/2.0) localPhi = localPhi - twopi/2.0 ;

        double pointPhi = thisPoint->phi();

        if(pointPhi>twopi) pointPhi=pointPhi - twopi;
        if(pointPhi<0.0) pointPhi=pointPhi + twopi;
        if(pointPhi>twopi/2.0) pointPhi = pointPhi - twopi/2.0 ;

        padPhi = fabs(pointPhi - localPhi);

#ifdef EXPERTCHECKPLOTS

        const float * mcpMomentum = SimTHit->getMomentum() ;

        CLHEP::Hep3Vector* mom = new CLHEP::Hep3Vector(mcpMomentum[0],mcpMomentum[1],mcpMomentum[2]);

        
        //        cout << "px = " << mcpMomentum[0] << " py = " << mcpMomentum[1] << " pz = " << mcpMomentum[2] << endl;
        
        double trackPhi = mom->phi();
        
        if(trackPhi<0.0) trackPhi=trackPhi+twopi;
        if(trackPhi>twopi) trackPhi=trackPhi-twopi;
        if(trackPhi>twopi/2.0) trackPhi = trackPhi - twopi/2.0 ;

        phiRelHisto->fill(padPhi);
        phiDiffHisto->fill((fabs(localPhi - trackPhi))/trackPhi);

        //        cout << "track Phi = " << trackPhi * (360.0 / twopi) << endl; 
        //        cout << "localPhi = " << localPhi * (360.0 / twopi) << endl; 
        //        cout << "pad Phi from track mom = " << ( pointPhi - trackPhi ) * (360.0 / twopi) << endl; 

#endif        



        // Calculate thetaPad for current hit

        double pathlength1 = 2.0 * asin( ( sqrt (
                                                 (thisPoint->x()-precedingPoint->x()) * (thisPoint->x()-precedingPoint->x())
                                                 +
                                                 (thisPoint->y()-precedingPoint->y()) * (thisPoint->y()-precedingPoint->y())
                                                 ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius()  ;
        

        double pathlength2 = 2.0 * asin( ( sqrt (
                                                 (followingPoint->x()-thisPoint->x()) * (followingPoint->x()-thisPoint->x())
                                                 +
                                                 (followingPoint->y()-thisPoint->y()) * (followingPoint->y()-thisPoint->y())
                                                 ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius()  ;
        

        padTheta = atan ((fabs(pathlength1 + pathlength2)) / (fabs(nextSimTHit->getPosition()[2] - previousSimTHit->getPosition()[2])) ) ;
          

#ifdef EXPERTCHECKPLOTS
        thetaRelHisto->fill(padTheta);
#endif    

#ifdef EXPERTCHECKPLOTS
        thetaDiffHisto->fill( (sin(padTheta) - sin(mom->theta()))/sin(mom->theta()) );
        delete mom;
#endif     

        delete precedingPoint;
        delete followingPoint; 


      }
      
      else if(mcp!=previousMCP && i < (n_sim_hits-2) ) { // first hit with at least two more hits in collection
        // if this is the first hit for this track try to get the next two hits    

        followingPoint = new CLHEP::Hep2Vector(nextSimTHit->getPosition()[0],nextSimTHit->getPosition()[1]);

        SimTrackerHit* nPlus2SimHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i+2 ) ) ;
        EVENT::MCParticle* nPlus2MCP = nPlus2SimHit->getMCParticle() ;

        nPlus2Point = new CLHEP::Hep2Vector(nPlus2SimHit->getPosition()[0], nPlus2SimHit->getPosition()[1]);

        if ( mcp==nextMCP && mcp==nPlus2MCP ){

          Circle theCircle(thisPoint, followingPoint, nPlus2Point);
          
          //        cout << "the circle has radius = " << theCircle.GetRadius() << endl;
          //        cout << "the circle has center = " << theCircle.GetCenter()->x() << "  " << theCircle.GetCenter()->y() << endl;
          
          // cout << "point x = " << SimTHit->getPosition()[0] << " point y = " << SimTHit->getPosition()[1] << endl; 
          
          double localPhi = 
            atan2((thisPoint->y() - theCircle.GetCenter()->y()) ,(thisPoint->x() - theCircle.GetCenter()->x())) + (twopi/4.0) ;
          
          if(localPhi>twopi) localPhi=localPhi - twopi;
          if(localPhi<0.0) localPhi=localPhi + twopi;
          if(localPhi>twopi/2.0) localPhi = localPhi - twopi/2.0 ;

          double pointPhi = thisPoint->phi();

          if(pointPhi<0.0) pointPhi=pointPhi + twopi;          
          if(pointPhi>twopi) pointPhi=pointPhi - twopi;
          if(pointPhi>twopi/2.0) pointPhi = pointPhi - twopi/2.0 ;
          
          padPhi = fabs(pointPhi - localPhi);

#ifdef EXPERTCHECKPLOTS                    
          const float * mcpMomentum = SimTHit->getMomentum() ;
          CLHEP::Hep3Vector *mom = new CLHEP::Hep3Vector(mcpMomentum[0],mcpMomentum[1],mcpMomentum[2]);
          
          //        cout << "px = " << mcpMomentum[0] << " py = " << mcpMomentum[1] << " pz = " << mcpMomentum[2] << endl;
          
          double trackPhi = mom->phi();
          
          if(trackPhi<0.0) trackPhi=trackPhi+twopi;
          if(trackPhi>twopi) trackPhi=trackPhi-twopi;
          if(trackPhi>twopi/2.0) trackPhi = trackPhi - twopi/2.0 ;          


          phiRelHisto->fill(padPhi);
          
          phiDiffHisto->fill((fabs(localPhi - trackPhi))/trackPhi);
#endif          
          
          
          // Calculate thetaPad for current hit
          
          double pathlength1 = 2.0 * asin( ( sqrt (
                                                   (followingPoint->x()-thisPoint->x()) * (followingPoint->x()-thisPoint->x())
                                                   +
                                                   (followingPoint->y()-thisPoint->y()) * (followingPoint->y()-thisPoint->y())
                                                   ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius()  ;
          
          
          double pathlength2 = 2.0 * asin( ( sqrt (
                                                   (followingPoint->x()-nPlus2Point->x()) * (followingPoint->x()-nPlus2Point->x())
                                                   +
                                                   (followingPoint->y()-nPlus2Point->y()) * (followingPoint->y()-nPlus2Point->y())
                                                   ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius()  ;
          
          
          padTheta = atan ((fabs(pathlength1 + pathlength2)) / (fabs(SimTHit->getPosition()[2] - nPlus2SimHit->getPosition()[2])) ) ;
          
  
#ifdef EXPERTCHECKPLOTS
          thetaRelHisto->fill(padTheta);
#endif    
          
#ifdef EXPERTCHECKPLOTS
          thetaDiffHisto->fill( (sin(padTheta) - sin(mom->theta()))/sin(mom->theta()) );

          delete mom;
#endif     
         

 
        }     
        
        else{ 
          // less than three Sim hits for this MC particle 
          // won't be able to fit anything so just set nominal values theta=phi=90 
          padTheta = twopi/4.0 ;
          padPhi = twopi/4.0 ;    
        }

        delete followingPoint;
        delete nPlus2Point;
       
      }
      else if(mcp!=nextMCP && i > 1 ) { // last hit with at least three sim hits in collection
        // if this is the last hit for this track take the last two hits

        precedingPoint = new CLHEP::Hep2Vector(previousSimTHit->getPosition()[0],previousSimTHit->getPosition()[1]);

        SimTrackerHit* nMinus2SimHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i-2 ) ) ;
        EVENT::MCParticle* nMinus2MCP= nMinus2SimHit->getMCParticle() ;

        nMinus2Point = new CLHEP::Hep2Vector(nMinus2SimHit->getPosition()[0], nMinus2SimHit->getPosition()[1]);

        if ( mcp==previousMCP && mcp==nMinus2MCP ){

          Circle theCircle(nMinus2Point, precedingPoint, thisPoint);
          
          //        cout << "the circle has radius = " << theCircle.GetRadius() << endl;
          //        cout << "the circle has center = " << theCircle.GetCenter()->x() << "  " << theCircle.GetCenter()->y() << endl;
          
          // cout << "point x = " << SimTHit->getPosition()[0] << " point y = " << SimTHit->getPosition()[1] << endl; 
          
          double localPhi = 
            atan2((thisPoint->y() - theCircle.GetCenter()->y()) ,(thisPoint->x() - theCircle.GetCenter()->x())) + (twopi/4.0) ;
          
          if(localPhi>twopi) localPhi=localPhi - twopi;
          if(localPhi<0.0) localPhi=localPhi + twopi;
          if(localPhi>twopi/2.0) localPhi = localPhi - twopi/2.0 ;

          double pointPhi = thisPoint->phi();
          

          if(pointPhi<0.0) pointPhi=pointPhi + twopi;
          if(pointPhi>twopi) pointPhi=pointPhi - twopi;
          if(pointPhi>twopi/2.0) pointPhi = pointPhi - twopi/2.0 ;

          padPhi = fabs(pointPhi - localPhi);

#ifdef EXPERTCHECKPLOTS          
          const float * mcpMomentum = SimTHit->getMomentum() ;
          
          CLHEP::Hep3Vector *mom = new CLHEP::Hep3Vector(mcpMomentum[0],mcpMomentum[1],mcpMomentum[2]);
          
          //        cout << "px = " << mcpMomentum[0] << " py = " << mcpMomentum[1] << " pz = " << mcpMomentum[2] << endl;
          
          double trackPhi = mom->phi();
          
          if(trackPhi<0.0) trackPhi=trackPhi+twopi;
          if(trackPhi>twopi) trackPhi=trackPhi-twopi;
          if(trackPhi>twopi/2.0) trackPhi = trackPhi - twopi/2.0 ;          

          
          phiRelHisto->fill(padPhi);
          phiDiffHisto->fill((fabs(localPhi - trackPhi))/trackPhi);
#endif              


          // Calculate thetaPad for current hit
          
          double pathlength1 = 2.0 * asin( ( sqrt (
                                                   (precedingPoint->x()-nMinus2Point->x()) * (precedingPoint->x()-nMinus2Point->x())
                                                   +
                                                   (precedingPoint->y()-nMinus2Point->y()) * (precedingPoint->y()-nMinus2Point->y())
                                                   ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius()  ;
          
          
          double pathlength2 = 2.0 * asin( ( sqrt (
                                                   (thisPoint->x()-precedingPoint->x()) * (thisPoint->x()-precedingPoint->x())
                                                   +
                                                   (thisPoint->y()-precedingPoint->y()) * (thisPoint->y()-precedingPoint->y())
                                                   ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius()  ;
        
          
          padTheta = atan ((fabs(pathlength1 + pathlength2)) / (fabs(nMinus2SimHit->getPosition()[2] - SimTHit->getPosition()[2])) ) ;
          
#ifdef EXPERTCHECKPLOTS
          thetaRelHisto->fill(padTheta);
#endif    
#ifdef EXPERTCHECKPLOTS

          thetaDiffHisto->fill( (sin(padTheta) - sin(mom->theta()))/sin(mom->theta()) );
          //          cout << "Padtheta = " << padTheta << endl;
          //          cout << "Theta from track = " << mom->theta() << endl;
          //          cout << "sin PadTheta  = " <<  sin(padTheta) << endl;
          //          cout << "sin Track Theta  = " <<  sin(mom->theta()) << endl;
          //          
          
          delete mom;

#endif     
          
          
        }
        
        else{ 
          // less than three Sim hits for this MC particle 
          // won't be able to fit anything so just set nominal values theta=phi=PI/2.0 
          padTheta = twopi/4.0 ;
          padPhi = twopi/4.0 ;    
        }    

        delete precedingPoint;
        delete nMinus2Point;

      }

      else {
        // less than three hits Sim hits 
        // won't be able to fit anything so just set nominal values theta=phi=90 
        padTheta = twopi/4.0 ;
        padPhi = twopi/4.0 ;

      }
      








      //       double *pos;
      //       pos = (double*) SimTHit->getPosition();  
      double pos[3] ; // fg: create a copy of the position in order to not modify the sim hit
      pos[0] = SimTHit->getPosition()[0] ;
      pos[1] = SimTHit->getPosition()[1] ; 
      pos[2] = SimTHit->getPosition()[2] ;
      int layerNumber = SimTHit->getCellID();

      if(_rejectCellID0 && (layerNumber<1)) continue;

      de_dx = SimTHit->getdEdx();



      //       cout << "x position for this hit is " << pos[0] << " - " << SimTHit->getPosition()[0] << endl; 
      //       cout << "y position for this hit is " << pos[1] << " - " << SimTHit->getPosition()[1] << endl; 
      //       cout << "z position for this hit is " << pos[2] << " - " << SimTHit->getPosition()[2] << endl; 



      //       cout << "de/dx for this hit is " << de_dx << endl; 
      //       cout << "MCParticle PID for this hit is " << mcp->getPDG() << endl; 
      //      cout << "x =  " << x << endl; 
      
      //  SMEARING
      
      // Calculate Point Resolutions according to Ron's Formula 

      // sigma_{RPhi}^2 = sigma_0^2 + Cd^2/N_{eff} * L_{drift}

      // sigma_0^2 = (50micron)^2 + (900micron*sin(phi))^2
      // Cd^2/N_{eff}} = 25^2/(22/sin(theta)*h/6mm)
      // Cd = 25 ( microns / cm^(1/2) )
      // (this is for B=4T, h is the pad height = pad-row pitch in mm,
      // theta is the polar angle)       

      // sigma_{z}^2 = (400microns)^2 + L_{drift}cm * (80micron/sqrt(cm))^2 

      double aReso =_pointResoRPhi0*_pointResoRPhi0 + (_pointResoPadPhi*_pointResoPadPhi * sin(padPhi)*sin(padPhi)) ;
      double driftLength = gearTPC.getMaxDriftLength() - (fabs(pos[2])-_cathode);
      if (driftLength <0) { 
        std::cout << " TPCDigiProcessor : Warning! driftLength < 0 " << driftLength << " --> Check out your GEAR file!!!!" << std::endl; 
        std::cout << "Setting driftLength to 0.1" << std::endl;
        std::cout << "gearTPC.getMaxDriftLength() = " << gearTPC.getMaxDriftLength() << std::endl; 
        driftLength = 0.10;
      }

      double padheight = padLayout.getPadHeight(padLayout.getNearestPad(thisPoint->r(),thisPoint->phi()));

      double bReso = ((_diffRPhi*_diffRPhi)/_nEff) * sin(padTheta) * (6.0/(padheight));

      double tpcRPhiRes = sqrt(aReso + bReso*(driftLength/10.0)); // driftLength in cm

      double tpcZRes  = sqrt((_pointResoZ0*_pointResoZ0) 
                             + 
                             (_diffZ*_diffZ) *(driftLength/10.0)); // driftLength in cm 

      //       std::cout << "aReso = " << aReso << std::endl;
      //       std::cout << "bReso = " << bReso << std::endl;
      //       std::cout << "_pointResoPadPhi = " <<_pointResoPadPhi << std::endl;
      //       std::cout << "_pointResoRPhi = " <<_pointResoRPhi << std::endl;
      //       std::cout << "_diffRPhi = " << _diffRPhi << std::endl;
      //       std::cout << "tpcRPhiRes = " << tpcRPhiRes << std::endl;
      //       std::cout << "_pointResoZ = " << _pointResoZ << std::endl;

      double randrp = gsl_ran_gaussian(_random,tpcRPhiRes);
      double randz =  gsl_ran_gaussian(_random,tpcZRes);

      // Make sure that the radius is equal to a pad radius
      // Get current hit radius

      double rad = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
      double phi = atan2(pos[1],pos[0]);
      if (phi<0.) phi=phi+twopi;

      phi += randrp/rad;

      //
      //      cout << "x position before smearing " << pos[0] << endl;     
      //      cout << "y position before smearing " << pos[1] << endl;     
      //      cout << "z position before smearing " << pos[2] << endl;     
      //

      pos[0] = rad * cos(phi);
      pos[1] = rad * sin(phi);
      pos[2] = pos[2] + randz;

      //
      //      cout << "x position after smearing " << pos[0] << endl;     
      //      cout << "y position after smearing " << pos[1] << endl;     
      //      cout << "z position after smearing " << pos[2] << endl;     
      //
      
      // At this point hits are mearly smeared now they must be digitised to trackerhits


      // modified to ues GEAR


      int padIndex = padLayout.getNearestPad(rad,phi);

      const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;

      double gearRMin = planeExt[0] ;
      double gearRMax = planeExt[1] ;

      if(fabs(pos[2])>gearTPC.getMaxDriftLength()) pos[2] = (fabs(pos[2])/pos[2])*gearTPC.getMaxDriftLength();

      if(rad>gearRMax){
        pos[0] = gearRMax*cos(phi);
        pos[1] = gearRMax*sin(phi);
        cout << "Radius is greater than TPC RMax " << endl;
        cout << "rad = " << rad << endl; 
        cout << "the tpc OuterRadius = " << gearRMax << endl; 
      }

      if(rad<gearRMin){
        pos[0] = gearRMin*cos(phi);
        pos[1] = gearRMin*sin(phi);
        cout << "Radius is less than TPC RMin" << endl;
        cout << "rad = " << rad << endl; 
        cout << "the tpc InnerRadius = " << gearRMin << endl; 
      }
      
      int iRowHit = padLayout.getRowNumber(padIndex);


      //       cout << "padLayout.getPadWidth(0) = " <<padLayout.getPadWidth(0) << endl;  
      //       cout << "padLayout.getPadWidth(padIndex) = " <<padLayout.getPadWidth(padIndex) << endl;  
      //       cout << "padLayout.getPadHeight(0) = " <<padLayout.getPadHeight(0) << endl;  
      //       cout << "padLayout.getPadHeight(padIndex) = " <<padLayout.getPadHeight(padIndex) << endl;  

      //je: commented out next line as proposed by Kristian Harder
      //gear::Point2D padCoord = padLayout.getPadCenter(padIndex);

      //get phi index of current hit

      int iPhiHit = padLayout.getPadNumber(padIndex);

      int NumberOfTimeSlices =  (int) ((2.0 * gearTPC.getMaxDriftLength()) / _pixZ);

      //get z index of current hit

      int iZHit = (int) ( (float) NumberOfTimeSlices * 
                          ( gearTPC.getMaxDriftLength() + pos[2] ) / ( 2.0 * gearTPC.getMaxDriftLength() ) ) ;


      if(iZHit<0) iZHit=0;
      if(iZHit>NumberOfTimeSlices) iZHit=NumberOfTimeSlices;

      double posRPhi[2];
      
      posRPhi[0] = rad;
      posRPhi[1] = phi;
      
      Voxel_tpc * atpcVoxel = new Voxel_tpc(iRowHit,iPhiHit,iZHit, pos,  posRPhi, de_dx, tpcRPhiRes, tpcZRes);
      
      tpcRowHits.at(iRowHit).push_back(atpcVoxel);
      
      tpcHitMap[atpcVoxel] = SimTHit; 


      //      cout << "a voxel hit "<< i << " has been added to row " << iRowHit << endl;  

      previousMCP = mcp ;
      previousSimTHit = SimTHit;

      delete atpcVoxel; 
      delete thisPoint;

   }    

    streamlog_out(DEBUG) << "finished looping over simhits" << endl;
    // Add background hits here
    
    vector <Voxel_tpc *> row_hits;

    
    //      cout << "get the row hits" << endl;
    
    for (int i = 0; i<padLayout.getNRows(); ++i){
      
      row_hits = tpcRowHits[i];

      //      sort(row_hits.begin(), row_hits.end(), compare_phi );

      //       cout << "row = " << i << "  row_hits.size() = " << row_hits.size() << endl;

      for (unsigned int j = 0; j<row_hits.size(); ++j){
        
        //        cout << "got the row hit j = " << j << "  row_hits.size() = " << row_hits.size() << endl;
        
        for (unsigned int k = j+1; k<row_hits.size(); k++){

          //          cout << "got the row hits k = " << k << endl;
          
          if(row_hits[j]->getPhiIndex() > row_hits[k]->getPhiIndex()){
            //            cout << "phi index of hit ["<<j<<"] = " << row_hits[j]->getPhiIndex() << endl; 
            //            cout << "phi index of hit ["<<k<<"] = " << row_hits[k]->getPhiIndex() << endl; 
          }
          

          
          if(abs(row_hits[j]->getZIndex()-row_hits[k]->getZIndex())<=1){

            if(abs(row_hits[j]->getPhiIndex()-row_hits[k]->getPhiIndex())<=1||abs(row_hits[j]->getPhiIndex()-row_hits[k]->getPhiIndex())==(int)padLayout.getPadsInRow(i).size()){

              /*              
                              cout << "*&^*&*&*&*&*&*&*&**&*&*&*&*&&*&*&**&*&*&" << endl;
                              cout << "double hit candidate found in row: " << i <<  endl;

                              cout << "row_hits[j]->getZ() " << row_hits[j]->getZ() << endl;
                              cout << "row_hits[k]->getZ() " << row_hits[k]->getZ() << endl;
                

                              cout << "row_hits dX^2 + dY^2 " << 
                              ( ( ( row_hits[j]->getX() - row_hits[k]->getX() ) * 
                              ( row_hits[j]->getX() - row_hits[k]->getX() ) ) +
                              ( ( row_hits[j]->getY() - row_hits[k]->getY() ) * 
                              ( row_hits[j]->getY() - row_hits[k]->getY() ) ) )
                              << endl;
                              cout << "padLayout.getPadWidth(0)^2 " << padLayout.getPadWidth(0) * padLayout.getPadWidth(0)
                              << endl;
              */

              if(fabs( row_hits[j]->getZ() - row_hits[k]->getZ() ) < _pixZ ) {

                if((((row_hits[j]->getX()-row_hits[k]->getX())*(row_hits[j]->getX()-row_hits[k]->getX()))
                    +((row_hits[j]->getY()-row_hits[k]->getY())*(row_hits[j]->getY()-row_hits[k]->getY())))
                   //                      <  padLayout.getPadWidth(0) *  padLayout.getPadWidth(0) ){
                   // FIXME: SJA: the function getPadWidth does not return 2.2mm as set in gear_ldc.xml so use hard coded number for now
                   <  2.2 *  2.2 ){


                  //SimTrackerHit* Hit1 = tpcHitMap[row_hits[j]];
                  //SimTrackerHit* Hit2 = tpcHitMap[row_hits[k]];


                  /*
                    
                    cout << "double hit found in row: " << i 
                    << "    MCP1 : " << Hit1->getMCParticle()->id() 
                    << " (" <<  Hit1->getMCParticle()->getPDG() << ")   " 
                    << "Energy : " << Hit1->getMCParticle()->getEnergy() << "     " 
                    << "MCP2 : " << Hit2->getMCParticle()->id() 
                    << " (" <<  Hit1->getMCParticle()->getPDG() << ")   " 
                    << "Energy : " << Hit2->getMCParticle()->getEnergy()
                    << endl;
                  */


                  /*
                    float xd = (float)row_hits[j]->getX();
                    float yd = (float)row_hits[j]->getY();
                    float zd = (float)row_hits[j]->getZ();
	    
                    float xd2 = (float)row_hits[k]->getX();
                    float yd2 = (float)row_hits[k]->getY();
                    float zd2 = (float)row_hits[k]->getZ();
                  */

                  row_hits[j]->setAdjacent(row_hits[k]);
                  row_hits[k]->setAdjacent(row_hits[j]);
                }		  
              }
            }
          }
        }
        
        // FIXME:SJA: 
        // 	At this point the double hits have been identified and at present they are not added to the
        // 	tracker hit collection. What should be done with them will be decided later.
        
        if(row_hits[j]->getNumberOfAdjacent()==0){
          //store hit variables
          TrackerHitImpl* trkHit = new TrackerHitImpl ;
          double pos[3] = {row_hits[j]->getX(),row_hits[j]->getY(),row_hits[j]->getZ()}; 
          trkHit->setPosition(pos);
          trkHit->setdEdx(row_hits[j]->getdEdx());
          trkHit->setType( 500 );
         
          double driftLength = gearTPC.getMaxDriftLength() - (fabs(pos[2])-_cathode);
          if (driftLength <0) { 
            std::cout << " TPCDigiProcessor : Warning! driftLength < 0 "  
                      << driftLength << " --> Check out your GEAR file!!!! "<< std::endl; 
            std::cout << "Setting driftLength to 0.1" << std::endl;
            std::cout << "gearTPC.getMaxDriftLength() = " << gearTPC.getMaxDriftLength() << std::endl; 
            driftLength = 0.10;
          }

          double rSqrd = row_hits[j]->getR()*row_hits[j]->getR();
          double phi = row_hits[j]->getPhi();
          double tpcRPhiRes = row_hits[j]->getRPhiRes();
          double tpcZRes = row_hits[j]->getZRes();

          // For no error in R
          float covMat[TRKHITNCOVMATRIX]={rSqrd*sin(phi)*sin(phi)*tpcRPhiRes*tpcRPhiRes,
                                          -(rSqrd)*cos(phi)*sin(phi)*tpcRPhiRes*tpcRPhiRes,
                                          rSqrd*cos(phi)*cos(phi)*tpcRPhiRes*tpcRPhiRes,
                                          0.,
                                          0.,
                                          float(tpcZRes*tpcZRes)};
          
          
          trkHit->setCovMatrix(covMat);      
          
          if(pos[0]*pos[0]+pos[1]*pos[1]>0.0 * 0.0){
            // 	  push back the SimTHit for this TrackerHit
            trkHit->rawHits().push_back( tpcHitMap[row_hits[j]] );                        
            trkhitVec->addElement( trkHit ); 
            rechits++;
          }

#ifdef EXPERTCHECKPLOTS
          SimTrackerHit* theSimHit = tpcHitMap[row_hits[j]];
          double rSimSqrd = theSimHit->getPosition()[0]*theSimHit->getPosition()[0] + theSimHit->getPosition()[1]*theSimHit->getPosition()[1];
          double phiSim = atan2(theSimHit->getPosition()[1],theSimHit->getPosition()[0]);
          
          double rDiff = (rSqrd - rSimSqrd);
          double phiPull = (phi - phiSim)/(sqrt((covMat[2])/(rSqrd*cos(phi)*cos(phi))));

          rDiffHisto->fill(rDiff);
          phiPullHisto->fill(phiPull);
          phiDistHisto->fill(phi - phiSim);

#endif

        }
      }
    }
    

    streamlog_out(DEBUG) << "number of rec_hits = "  << rechits << endl ;
      streamlog_out(DEBUG) << "finished row hits" << endl;    

        // set the parameters to decode the type information in the collection
        // for the time being this has to be done manually
        // in the future we should provide a more convenient mechanism to 
        // decode this sort of meta information
        StringVec typeNames ;
        IntVec typeValues ;
        typeNames.push_back( LCIO::TPCHIT ) ;
        typeValues.push_back( 1 ) ;
        trkhitVec->parameters().setValues("TrackerHitTypeNames" , typeNames ) ;
        trkhitVec->parameters().setValues("TrackerHitTypeValues" , typeValues ) ;
    
        evt->addCollection( trkhitVec , _TPCTrackerHitsCol ) ;
    
    
        for (int i = 0; i<padLayout.getNRows(); ++i){
      
          vector <Voxel_tpc *> current_row = tpcRowHits[i];
      
          for (unsigned int j = 0; j<current_row.size(); ++j){
        
            delete current_row[j];
          }
        }    
  }
  _nEvt++;
  
}



void TPCDigiProcessor::check( LCEvent * evt ) 
{ 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TPCDigiProcessor::end()
{ 

#ifdef EXPERTCHECKPLOTS
  TREE->commit();
  TREE->cd("/Histograms");
  TREE->ls("..");

  TREE->close();  
  cout << "EXPERTCHECKPLOTS Finished" << endl;
#endif

  gsl_rng_free(_random);
  cout << "TPCDigiProcessor::end()  " << name() 
       << " processed " << _nEvt << " events in " << _nRun << " runs "
       << endl ;
  //  
}
