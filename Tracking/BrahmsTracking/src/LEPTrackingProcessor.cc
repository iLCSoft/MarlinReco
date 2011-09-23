/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
** This file is part of the MarlinReco Project.
** Forming part of the SubPackage: BrahmsTracking.
**
** For the latest version download from Web CVS:
** www.blah.de
**
** $Id: LEPTrackingProcessor.cc,v 1.44 2008-11-13 10:18:17 aplin Exp $
**
** $Log: not supported by cvs2svn $
** Revision 1.43  2008/07/18 11:00:54  aplin
** corrected bug in the clear function of hitbanks
**
** Revision 1.42  2008/07/17 16:53:30  aplin
** added steering parameter AlwaysRunCurlKiller which allows the CurlKiller functionality to be usedfor every event, not just those for whom the F77 tracking code fails
**
** Revision 1.41  2008/07/16 08:39:46  aplin
** made subdetectorHitNumbers for TPCTracks consistent with FullLDCTracking which has implemented ETD and SET hits
**
** Revision 1.40  2008/07/02 11:09:27  aplin
** fixed bug where the wrong TPCHit collection was used to created relations to Tracks
**
** Revision 1.39  2008/07/02 09:01:23  aplin
** Changed default binning for the removal of hits alla CurlKiller
**
** Fixed bug which caused the hit removal to start at the second level of aggression
**
** Revision 1.38  2008/07/01 15:08:31  aplin
** An Error message has been added as an event parameter for the case where the PATREC fails due to too many links.
**
** Only one track collection is now produced, namely "TPCTracks", as LEPTracking now only produces hits from the TPC. FullLDCTracking is now responsible for creating tracks including hits from the silicon detectors.
**
** Revision 1.37  2008/07/01 10:34:39  aplin
** The functionality of CurlKiller has been moved into LEPTracking. It is
** now operated at increasing levels of aggression in order to try to
** save the event from not not having ANY TPC Tracks, due to too many links in
** trslnk.
**
** As a first try, all hits are passed to the PATREC. If this fails three
** more attempts are made where the size of the bins used to find high
** levels of multiplicity in r-phi are increased by factors of 2, 3 and
** finally 4. If this still fails an error message is printed warning
** that NO TPC Tracks will be reconstructed for this event.
**
** The verbosity of LEPTracking is now controlled completely via
** streamlog on the C++ side, and via flags set in tkinit on the F77
** side.
**
** Revision 1.36  2008/06/27 15:32:27  aplin
** corrected conversion of covariance matrix, and made more use of streamlog
**
** Revision 1.35  2008/06/26 09:35:56  aplin
**
** A clean up has been made of the arrays that hold the TE's/Hits for
** track fitting. An effort has been made to ensure that both the arrays
** are big enough to hold the required number of hits, and that the arrays
** are all of the same dimensions between calling routines.
**
** All verbose output, except for the case where the tracking is dropped
** for a complete event, has now been hidden behind debug IF
** Statements. The level of verbosity can be controlled via the variables
** TPCDBG, TKSTDBG, FSFDBG and IDEB  in tkinit.F
**
** Revision 1.34  2008/04/17 14:36:29  aplin
** reduced printoutlevel in f77 and moved cout to debug streamlog_out for C++
**
** Revision 1.33  2008/04/10 21:49:14  aplin
** Corrected covariance matrix definition. Added extra Expert Histos in TPCDigi. Removed VTX and SIT Hits from LEPTracking
**
** Revision 1.32  2008/03/11 15:19:18  engels
** removed spurious dev code (S. Aplin)
**
** Revision 1.31  2008/02/14 15:34:11  aplin
** removed spurious dev code
**
** Revision 1.30  2008/02/14 14:43:30  aplin
** removed spurious dev code
**
** Revision 1.29  2008/02/13 10:23:30  aplin
** fixed mm->cm conversion bug
**
** Revision 1.28  2008/02/11 16:54:51  aplin
** ntroduced smearing which depends on the parameterisation supplied by the LC-TPC group.
** It is parameterised in phi and theta. Not yet in z.
**
** The Covariance matrix corrected so that it is set for XYZ and not r-phi as previously. The corresponding conversion back to r-phi has been set in LEPTracking and FullLDCTracking.
**
** Revision 1.26  2007/09/13 07:39:04  rasp
** Updated version of the LEPTracking processor. For each
** digitized TPC hit spatial resolutions are accessed via
** call to the TrackerHit::getCovMatrix method.
**
** Revision 1.24  2007/09/05 09:47:29  rasp
** Updated version
**
** Revision 1.22  2006/10/17 12:34:19  gaede
** replaced registerProcessorParameter with registerInput/OutputCollection
**
** Revision 1.21  2006/06/28 15:29:04  aplin
** The B-Field is now variable for LEPTracking via the gear xml file. The B-Field is specified in the TPCParameters as follows: <parameter name="tpcBField" type="double"> 4.0  </parameter>
**
** The value is passed internaly to the F77 code via the same function which passes the TPC geometry i.e. gettpcgeom(float* innerrad, float* outerrad, int* npadrows, float* maxdrift, float* tpcpixz, float* ionpoten, float* tpcrpres, float* tpczres, float* tpcbfield). It is set in setmat.F. tpcgeom.F had to be modified as it also uses gettpcgeom, although it does not make use of the B-Field.
**
** Revision 1.22   gaede
** added registerInput/OutputCollection for Marlin -c
** 
** Revision 1.20  2006/05/28 15:22:15  owendt
** changed text for the explanation of a processor parameter
**
** Revision 1.19  2006/04/27 13:07:43  samson
** Fix minor syntax errors to achieve compatibility with gcc4
**
** Revision 1.18  2006/02/09 18:00:41  owendt
** removed cout statements for debugging
**
** Revision 1.17  2006/02/03 15:09:11  owendt
** i) Corrected bug in calculation of weights, relocated brace.
** ii) Weights are now calculated as the percentage of hits that a given MC particle contributes to the reconstructed track's hit collection.
**
** Revision 1.16  2005/12/06 15:26:23  aplin
** corrected erroneous definition of MC Track Relation weight
**
** Revision 1.15  2005/11/03 15:16:14  aplin
** Added the Trackstring creation and the biulding of full Track candiates (TK's) which have passed the Delphi Ambiguity resolver fxambi. The material description of the vtx detector, as for the TPC, is hard coded in setmat. Presently the VTX and SIT resolutions are hard coded in LEPTrackingProcessor. The debug output has been reduced and can be controlled via TKSTDBG etc. in tkinit.F. delsolve contains the delphi ambuguity resolver written in C and is contained in the directory named C. The Tk's are written back into the C++ side in tktrev. The corresponding Tk bank structure analogous to the TE bank structure has been added in tktkbank whilst the access wrapper functions are contained in LEPTracking.
**
** Revision 1.13  2005/08/08 07:09:13  aplin
** Made f77 tracking code use GEAR to define the geomtery of the TPC. LTPDRO now defines the maximum number of rows is used to define the size of arrays, this is limited to 224 due the use of 7 '32 bit' bit registers in trkfnd.F increased, though at present it is not likely that anybody would want more. The number of TPC padrows is defined at run time by NRTPC which should of course not exceed LTPDRO, although this is checked and the programe exits with a verbose error message. A wrapper function gettpcgeom is used to pass the GEAR TPC parameters from C++ to f77. MarlinUtil/include/marlin_tpcgeom.h have MarlinUtil/src/marlin_tpcgeom.cc consequently been removed as they are no longer needed.
**
** Revision 1.12  2005/08/04 12:54:51  aplin
** *** empty log message ***
**
** Revision 1.11  2005/08/03 21:31:09  aplin
** tk*bank structures initialisation move here from BrahmsInitProcessor and BrahmEndProcessor
**
** Revision 1.10  2005/08/03 19:05:24  aplin
** corrected erroneous function declaration of tkmktecpp, by using float * instead of numerous floats and added output collection names as steering parametes
** 
*/
#include "LEPTrackingProcessor.h"
#include <iostream>
#include <string>
#include <stdexcept>
#include "marlin/Exceptions.h"
#include <cmath>


#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>


#include <cfortran.h>


#include"tpchitbank.h"
#include"tkhitbank.h"
#include"tktebank.h"
#include"tktkbank.h"
#include"tkmcbank.h"
//#include"marlin_tpcgeom.h"
#include"constants.h"

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gearimpl/FixedPadSizeDiskLayout.h>
#include <gear/BField.h>
//



PROTOCCALLSFFUN0(INT,TKTREV,tktrev)
#define TKTREV() CCALLSFFUN0(TKTREV,tktrev)

PROTOCCALLSFFUN0(INT,CHECKTPCROWS,checktpcrows)
#define CHECKTPCROWS() CCALLSFFUN0(CHECKTPCROWS,checktpcrows)

  // FIXME:SJA: the namespace should be used explicitly
using namespace lcio ;
using namespace marlin ;
using namespace constants ;
using namespace std ; 

int subdetfirsthitindex(string subdet);

FCALLSCFUN1(INT,subdetfirsthitindex,SUBDETFIRSTHITINDEX,subdetfirsthitindex, STRING)

  int numofsubdethits(string subdet);

FCALLSCFUN1(INT,numofsubdethits,NUMOFSUBDETHITS,numofsubdethits, STRING)

  int writetpccpp(float c, int a, int b); 

FCALLSCFUN3(INT,writetpccpp,WRITETPCCPP,writetpccpp, FLOAT, INT, INT)
      
  float readtpchitscpp(int a, int b); 

FCALLSCFUN2(FLOAT,readtpchitscpp,READTPCHITSCPP,readtpchitscpp, INT, INT)

  int writetkhitcpp(float c, int a, int b); 

FCALLSCFUN3(INT,writetkhitcpp,WRITETKHITCPP,writetkhitcpp, FLOAT, INT, INT)

  int ireadtkhitscpp(int a, int b); 

FCALLSCFUN2(INT,ireadtkhitscpp,IREADTKHITSCPP,ireadtkhitscpp, INT, INT)

  float rreadtkhitscpp(int a, int b); 

FCALLSCFUN2(FLOAT,rreadtkhitscpp,RREADTKHITSCPP,rreadtkhitscpp, INT, INT)


  int writetktkcpp(float c, int a, int b); 

FCALLSCFUN3(INT,writetktkcpp,WRITETKTKCPP,writetktkcpp, FLOAT, INT, INT)

  int writetktecpp(float c, int a, int b); 

FCALLSCFUN3(INT,writetktecpp,WRITETKTECPP,writetktecpp, FLOAT, INT, INT)

  float rreadtktkcpp(int a, int b);

FCALLSCFUN2(FLOAT,rreadtktkcpp,RREADTKTKCPP,rreadtktkcpp, INT, INT)

  float rreadtktecpp(int a, int b); 

FCALLSCFUN2(FLOAT,rreadtktecpp,RREADTKTECPP,rreadtktecpp, INT, INT)

  int ireadtktkcpp(int a, int b); 

FCALLSCFUN2(INT,ireadtktkcpp,IREADTKTKCPP,ireadtktkcpp, INT, INT)

  int ireadtktecpp(int a, int b); 

FCALLSCFUN2(INT,ireadtktecpp,IREADTKTECPP,ireadtktecpp, INT, INT)

  int tkmktkcpp(int modid,int subdetbits,int MesrCode,int tracktype, int numtes,int Charge,int unused,int ndf,float chi2,float L,float xstart, float ystart, float zstart, float xend, float yend, float zend, float cord1,float cord2,float cord3,float theta,float phi,float invp,float* cov);

FCALLSCFUN23(INT,tkmktkcpp,TKMKTKCPP,tkmktkcpp, INT , INT ,INT ,INT ,INT ,INT ,INT ,INT, FLOAT, FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT, FLOAT, FLOAT, FLOATV)
  
  int tkmktecpp(int subid,int submod,int unused,int MesrCode,int PnteTE,int Q,int ndf,float chi2,float L,float cord1,float cord2,float cord3,float theta,float phi,float invp,float dedx,float* cov);

FCALLSCFUN17(INT,tkmktecpp,TKMKTECPP,tkmktecpp, INT , INT ,INT ,INT ,INT ,INT ,INT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOATV)

  int addtetktkcpp(int a , int b);  

FCALLSCFUN2(INT,addtetktkcpp,ADDTETKTKCPP,addtetktkcpp, INT, INT)

  int addhittktecpp(int a, int b); 

FCALLSCFUN2(INT,addhittktecpp,ADDHITTKTECPP,addhittktecpp, INT, INT)

  int readtkitkdatcpp(int a, int b);

FCALLSCFUN2(INT,readtkitkdatcpp,READTKITKDATCPP,readtkitkdatcpp, INT, INT)

  int readtkitedatcpp(int a, int b); 

FCALLSCFUN2(INT,readtkitedatcpp,READTKITEDATCPP,readtkitedatcpp, INT, INT)

  int writetkitkdatcpp(int c, int a, int b);

FCALLSCFUN3(INT,writetkitkdatcpp,WRITETKITKDATCPP,writetkitkdatcpp, INT, INT, INT)

  int writetkitedatcpp(int c, int a, int b);

FCALLSCFUN3(INT,writetkitedatcpp,WRITETKITEDATCPP,writetkitedatcpp, INT, INT, INT)


  // definition of gettpcgeom done here it just sends the geomtertry into to tpcgeom.F
  int gettpcgeom(float* innerrad, float* outerrad, int* npadrows, 
                  float* maxdrift, float* tpcpixz, float* tpcrpres, float* tpczres, float* tpcbfield){

  //  try{

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;
  
  *innerrad = 0.1 * float( planeExt[0] ) ;
  *outerrad = 0.1 *float( planeExt[1] ) ;
  *npadrows = padLayout.getNRows() ;
  *maxdrift = 0.1 * float( gearTPC.getMaxDriftLength() );
  //  *tpcpixz = 0.1 * float(gearTPC.getDoubleVal("tpcPixZ")) ;
  // FIXME:SJA: should this really be multiplied by 0.1 or has it just got caught up in the other mm->cm convertions 
  // *ionpoten = 0.1 * float(gearTPC.getDoubleVal("tpcIonPotential")) ;  
  //  *tpcrpres = 0.1 * float(gearTPC.getDoubleVal("tpcRPhiResConst")) ;  
  //  *tpczres = 0.1 * float(gearTPC.getDoubleVal("tpcZRes")) ;
  //  *tpcbfield = float(gearTPC.getDoubleVal("BField")) ;
  *tpcbfield = float(Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z());

  //  }
//  catch() {return 1} ;

  return 0;
}

FCALLSCFUN8(INT,gettpcgeom,GETTPCGEOM,gettpcgeom, PFLOAT, PFLOAT, PINT, 
            PFLOAT, PFLOAT, PFLOAT, PFLOAT, PFLOAT )




  // end of cfortran.h definitions

// function to generate a unique key for each [r][phi] bin based on a 64bit bitfield
// the most significant 32 bits are used for r 
unsigned long long make_keyNew( unsigned r, unsigned phi ){

  unsigned long long temp = 0xffffffff & r ;

  temp = temp << 32 ;

  return  ( temp )  |  ( 0xffffffff & phi )  ;
} 


// Map to store the enteries for a 2D(r-phi) projection of the Tracker hits 
typedef std::map< unsigned long long  , std::vector<EVENT::TrackerHit*> >  HitMap ; 

  LEPTrackingProcessor aLEPTrackingProcessor ;

LEPTrackingProcessor::LEPTrackingProcessor() : Processor("LEPTrackingProcessor") {
  
  // modify processor description
  _description = "Produces Track collection from TPC TrackerHit collections using LEP tracking algorithms" ;

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::TRACKERHIT,
                           "TPCTrackerHitCollectionName" , 
                           "Name of the TPC TrackerHit collection"  ,
                           _colNameTPC ,
                           std::string("TPCTrackerHits") ) ;

//  registerInputCollection( LCIO::TRACKERHIT,
//                           "VTXTrackerHitCollectionName" , 
//                           "Name of the VTX TrackerHit collection"  ,
//                           _colNameVTX ,
//                           std::string("VTXTrackerHits") ) ;
  
//  registerInputCollection( LCIO::TRACKERHIT,
//                           "SITTrackerHitCollectionName" , 
//                           "Name of the SIT TrackerHit collection"  ,
//                           _colNameSIT ,
//                           std::string("SITTrackerHits") ) ;
  
  registerOutputCollection( LCIO::TRACK,
                            "TPCTrackCollectionName" , 
                            "Name of the TPC Track collection"  ,
                            _colNameTPCTracks ,
                            std::string("TPCTracks") ) ;

//  registerOutputCollection( LCIO::TRACK,
//                            "TrackCollectionName" , 
//                            "Name of the Track collection"  ,
//                            _colNameTracks ,
//                            std::string("Tracks") ) ;
  
  registerOutputCollection( LCIO::LCRELATION,
                            "MCTPCTrackRelCollectionName" , 
                            "Name of the TPC Track MC Relation collection"  ,
                            _colNameMCTPCTracksRel ,
                            std::string("TPCTracksMCP") ) ;
  
//  registerOutputCollection( LCIO::LCRELATION,
//                            "MCTrackRelCollectionName" , 
//                            "Name of the Track MC Relation collection"  ,
//                            _colNameMCTracksRel ,
//                            std::string("TracksMCP") ) ;

  registerOutputCollection( LCIO::TRACKERHIT,
                            "DroppedCollectionName" , 
                            "Name of the TrackerHit collection NOT used in PATREC"  ,
                            _droppedColName ,
                            std::string("DroppedTPCTrackeHits") ) ;
  
  registerOutputCollection( LCIO::TRACKERHIT,
                            "UsedCollectionName" , 
                            "Name of the TrackerHit collection used in PATREC"  ,
                            _usedColName ,
                            std::string("UsedTPCTrackerHits") ) ;
  
  
  registerProcessorParameter( "BinHeight" , 
                              "Bin Height in pad multiples"  ,
                              _binHeight ,
                              int(1) ) ;

  registerProcessorParameter( "BinWidth" , 
                              "Bin Width in pad multiples"  ,
                              _binWidth ,
                              int(3) ) ;

  registerProcessorParameter( "MultiplicityCut" , 
                              "Cut on the number of hits in r-phi bin"  ,
                              _multiplicityCut ,
                              int(8) ) ;


  registerProcessorParameter( "SlicesInZ" , 
                              "Number of slices in Z"  ,
                              _nSlicesInZ ,
                              int(1) ) ;


  registerProcessorParameter( "AlwaysRunCurlKiller" , 
                              "No attempt will be made to run without CurlKiller functionallity"  ,
                              _AlwaysRunCurlKiller,
                              int(0) ) ;


  registerProcessorParameter("Histograms" , 
			     "Hit times histograms" ,
			     _histograms,
			     (int)0);


}


void LEPTrackingProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0; 
  if(_histograms){
    fTPCRRaw    = new TH1F("fTPCRRaw", "Raw TPC Radial profile",334, 0., 2004.0);
    fTPCR       = new TH1F("fTPCR",    "Cut TPC Radial profile",334, 0., 2004.0);
    fTPCZRaw    = new TH1F("fTPCZRaw", "Raw TPC z profile",500, -2500., 2500.);
    fTPCZ       = new TH1F("fTPCZ",    "Cut TPC z profile",500, -2500., 2500.);
    fTPCRZRaw   = new TH2F("fTPCRZRaw","Raw TPC rz profile",5000, -2500., 2500.,334,0.,2004.);
    fTPCRZ      = new TH2F("fTPCRZ",   "Cut TPC rz profile",5000, -2500., 2500.,334,0.,2004.);
    fTPCXYRaw   = new TH2F("fTPCXYRaw","Raw TPC xy profile",1000, -2000., 2000., 1000, -2000., 2000.);
    fTPCXY      = new TH2F("fTPCXY",   "Cut TPC xy profile",1000, -2000., 2500., 1000, -2000., 2000.);
  }

  
}

void LEPTrackingProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void LEPTrackingProcessor::processEvent( LCEvent * evt ) { 

  static bool firstEvent = true ;
  
  //  f77histos->fill1DHist("testhisto",gsl_ran_gaussian(_random,1.0));
  
  // this gets called for every event 
  // usually the working horse ...
  
//   int skipToEvent = 3 ;
//   if(_nEvt<skipToEvent) {
//     ++_nEvt;
//     return;
//   }

  if(firstEvent==true) streamlog_out(MESSAGE) << "LEPTrackingProcessor called for first event" << endl;

  firstEvent = false ;


  // create bank structure
  TkMCBank = new Tk_MC_Bank;
  TPCHitBank = new TPC_Hit_Bank;  
  TkHitBank = new Tk_Hit_Bank;  
  TkTeBank = new Tk_Te_Bank;  
  TkTkBank = new Tk_Tk_Bank;  


  LCCollection* tpcTHcol = 0 ;

  try{
    tpcTHcol = evt->getCollection( _colNameTPC ) ;
  }
  catch(DataNotAvailableException &e){
  }
  

  LCCollection* vtxTHcol = 0 ;
  try{
    vtxTHcol = evt->getCollection( _colNameVTX ) ;
  }
  catch(DataNotAvailableException &e){
  }
  
  LCCollection* sitTHcol = 0 ;
  try{
    sitTHcol = evt->getCollection( _colNameSIT ) ;
  }
  catch(DataNotAvailableException &e){
  }
  
  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

  if (padLayout.getNRows()>CHECKTPCROWS()) 
      {
      std::stringstream errorMsg;
      errorMsg << "\nProcessor: LEPTracking \n" <<
        "The number of TPC padrows specified in the GEAR file is " << padLayout.getNRows() << "\n" 
               << "This is larger than the max number of rows that the code can handle." << "\n"
               << "The maximum number is determined by LTPDRO which is currently set to " 
               << CHECKTPCROWS() << "\n"
               << "LTPDRO must be a multiple of 32, and is defined as N32BITREG*32 \n"
               << "Increase N32BITREG in ./src/f77/include/padrow.inc \n"
               << "gmake clean and recompile"
               << "\n" ;
      throw gear::Exception(errorMsg.str());
    }


  if( tpcTHcol != 0 ){
    
    _tpcTrackVec = new LCCollectionVec( LCIO::TRACK )  ;
    //LCCollectionVec* TrackVec = new LCCollectionVec( LCIO::TRACK )  ;
    
    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit( LCIO::TRBIT_HITS ) ;
    _tpcTrackVec->setFlag( trkFlag.getFlag()  ) ;
    //TrackVec->setFlag( trkFlag.getFlag()  ) ;
    
    //LCCollectionVec* lcRelVec = new LCCollectionVec( LCIO::LCRELATION )  ;
    _tpclcRelVec = new LCCollectionVec( LCIO::LCRELATION )  ;
    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    _tpclcRelVec->setFlag( lcFlag.getFlag()  ) ;

    _droppedCol = NULL;
    _usedCol    = NULL;

    _droppedColForOutput = new LCCollectionVec( LCIO::TRACKERHIT ) ;
    _droppedColForOutput->setSubset() ;     
    _usedColForOutput = new LCCollectionVec( LCIO::TRACKERHIT ) ;
    _usedColForOutput->setSubset() ; 
    evt->addCollection( _usedColForOutput , _usedColName ) ;
    evt->addCollection( _droppedColForOutput , _droppedColName ) ;

    int nTPCHits = tpcTHcol->getNumberOfElements()  ;   
    streamlog_out(WARNING) << "Number of TPCHit before filtering: " << nTPCHits << endl;




    
    _goodHits.reserve(nTPCHits);
    _savedGoodHits.reserve(nTPCHits);

    _goodHits.clear();
    _savedGoodHits.clear();



    float slices = float(_nSlicesInZ);
    if(slices<1)slices = 1.0;
    float maxdrift = (gearTPC.getMaxDriftLength()+10);
    float zSlice = 2*maxdrift/slices;
    float zmin = -maxdrift-zSlice;
    float zmax = -maxdrift;
    for(int ipass=0; ipass<_nSlicesInZ;ipass++){
      zmin += zSlice;
      zmax += zSlice;
      //if(ipass==0)continue;

     if(_usedCol!=NULL){
        delete _droppedCol;
        delete _usedCol;
      }
      _droppedCol = new LCCollectionVec( LCIO::TRACKERHIT ) ;
      _droppedCol->setSubset() ;     
      _usedCol = new LCCollectionVec( LCIO::TRACKERHIT ) ;
      _usedCol->setSubset() ; 
      _goodHits.clear();


      int errTKTREV = 0;
    
      if( _AlwaysRunCurlKiller == 0 ) {
        selectTPCHits(tpcTHcol,_usedCol,zmin,zmax);      
        streamlog_out(WARNING) << "Number of TPCHit passed to PATREC: " << _goodHits.size() << endl;
        FillTPCHitBanks();            
        errTKTREV = TKTREV();       
      }
      
      if( errTKTREV==911 || _AlwaysRunCurlKiller != 0 ){
        
        streamlog_out(WARNING) << endl;
        if( _AlwaysRunCurlKiller == 0 ) streamlog_out(WARNING) << "   LEPTrackingProcessor: TKTREV returns:" << errTKTREV << endl;
        streamlog_out(WARNING) << "   LEPTrackingProcessor: Trying to remove hits alla CurlKiller" << endl;
        streamlog_out(WARNING) << endl;
        
        for(int i=1 ; i<4; ++i){
          

          if(_usedCol!=NULL){
            delete _droppedCol;
            delete _usedCol;
          }
          _droppedCol = new LCCollectionVec( LCIO::TRACKERHIT ) ;
          _droppedCol->setSubset() ;     
          _usedCol = new LCCollectionVec( LCIO::TRACKERHIT ) ;
          _usedCol->setSubset() ; 
          _goodHits.clear();          
          TkMCBank->clear();
          TkTeBank->clear();
          TkTkBank->clear();
          
          selectTPCHits(tpcTHcol, _usedCol, _droppedCol, _binHeight*(i),_binWidth*(i),zmin,zmax);
          
          FillTPCHitBanks();
        
          streamlog_out(WARNING) << "Number of TPCHit after filtering: " << _goodHits.size() << endl;
          
          errTKTREV = TKTREV();  

          if(errTKTREV!=911) {
            break;
          }
        
          else if(i==3){
            streamlog_out(ERROR) << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
            streamlog_out(ERROR) << "   LEPTrackingProcessor: TKTREV returns:" << errTKTREV << endl;
            streamlog_out(ERROR) << "   LEPTrackingProcessor: Removing hits failed to resolve the problem" << endl;          
            streamlog_out(ERROR) << "   LEPTrackingProcessor: This event contains NO TPC TRACKS" << endl;          
            streamlog_out(ERROR) << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
            
            stringstream str;
            str << "_Failure" 
                << "_LEPTracking" ;
          
            stringstream error ;
            error << "Too Many Links" ;
            
            evt->parameters().setValue(str.str(),error.str());

          }
          else {
            streamlog_out(WARNING) << "TKTREV return:" << errTKTREV << " trying to remove more hits using bigger bins" << endl;
          }
    
        }
      }


   
      streamlog_out(WARNING) << "For Event:" << _nEvt << " pass " << ipass << " z = " << zmin << " - " << zmax << " TKTREV return:" << errTKTREV << endl;

      streamlog_out(WARNING) << "number of TE's = " << TkTeBank->size() << endl ;

      streamlog_out(WARNING) << "number of TK's = " << TkTkBank->size() << endl ;


      // copy temporary hit arrays into out put collection
      int n = _usedCol->getNumberOfElements()  ;
      for(int i=0; i< n; i++){
        TrackerHit* THit = dynamic_cast<TrackerHit*>( _usedCol->getElementAt( i ) ) ;
        _usedColForOutput->addElement(THit);
      }
      n = _droppedCol->getNumberOfElements()  ;
      for(int i=0; i< n; i++){
        TrackerHit* THit = dynamic_cast<TrackerHit*>( _droppedCol->getElementAt( i ) ) ;
        _droppedColForOutput->addElement(THit);
      }
      for(unsigned int i=0;i<_goodHits.size();i++)_savedGoodHits.push_back(_goodHits[i]);


      for(int te=0; te<TkTeBank->size();te++){
        if( TkTeBank->getSubdetector_ID(te)==500 ) {
          TrackImpl* tpcTrack = new TrackImpl ; 
          const double ref_r = 10.*TkTeBank->getCoord1_of_ref_point(te);
          const double ref_phi =TkTeBank->getCoord2_of_ref_point(te)/TkTeBank->getCoord1_of_ref_point(te);
          const double ref_z = 10.*TkTeBank->getCoord3_of_ref_point(te);

          //         cout << "ref_r = " << ref_r << endl;
          //         cout << "ref_phi = " << ref_phi << endl;
          
          // transformation from 1/p to 1/R = consb * (1/p) / sin(theta)
          // consb is given by 1/R = (c*B)/(pt*10^9) where B is in T and pt in GeV  
          
          //        const double bField = gearTPC.getDoubleVal("BField") ;
          const double bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
          const double consb = (2.99792458* bField )/(10*1000.) ;     // divide by 1000 m->mm

          // tan lambda and curvature remain unchanged as the track is only extrapolated
          // set negative as 1/p is signed with geometric curvature clockwise negative

          const double omega = ( -consb*TkTeBank->getInvp(te) )/ sin( TkTeBank->getTheta(te) ) ;
          const double tanLambda = tan( (twopi/4.) - TkTeBank->getTheta(te) ) ;

          tpcTrack->setOmega( omega ) ;
          tpcTrack->setTanLambda( tanLambda ) ;      

          // computation of D0 and Z0 taken from fkrtpe.F in Brahms
          
          // xref and yref of ref point 
          const double xref = ref_r*cos(ref_phi) ;
          const double yref = ref_r*sin(ref_phi) ;
          const double zref = ref_z ; 
          const double trkRadius = 1. / omega ;


          ////////////////////////////////
          
          // center of circumference
          const double xc = xref + trkRadius * sin( TkTeBank->getPhi(te) ) ;
          const double yc = yref - trkRadius * cos( TkTeBank->getPhi(te) ) ;
        
          const double xc2 = xc * xc ; 
          const double yc2 = yc * yc ;
        
          const double DCA = ( sqrt( xc2+yc2 ) - fabs(trkRadius) ) ;
          
          double phiOfPCA = atan2 (yc,xc) ;
          
          if (DCA<0.) phiOfPCA = phiOfPCA + twopi/2. ;
          
          if ( phiOfPCA < -twopi/2. ) phiOfPCA = twopi + phiOfPCA ;
          else if ( phiOfPCA > twopi/2. ) phiOfPCA = -twopi + phiOfPCA ;   
          
          const double x0 = fabs( DCA ) * cos( phiOfPCA ) ;
          const double y0 = fabs( DCA ) * sin( phiOfPCA ) ;
          
          double phi = phiOfPCA + twopi/2. - ( fabs(omega)/omega ) * ( fabs(DCA)/DCA ) * (twopi/4.) ;
          
          if (phi<-twopi/2.) phi =  twopi + phi ;
          else if (phi>twopi/2.)  phi = -twopi + phi ;
        
          const double d0 = y0 * cos( phi )  - x0 * sin( phi ) ;
          
          const double alpha = - omega * ( xref - x0 ) * cos( phi ) - omega * ( yref - y0 ) * sin( phi ) ;
          const double beta = 1.0 - omega * ( xref - x0 ) * sin( phi ) + omega * ( yref - y0 ) * cos( phi ) ;
          
          double dphi = atan2( alpha,beta ) ;
          
          double larc =  - dphi/ omega ;
          
          if (larc < 0. ) {
            if ( dphi < 0.0 ) dphi = dphi + twopi ;
            else dphi = dphi - twopi ;
            larc =  - dphi/ omega ;
          }
          
          double z0 = zref - larc * tanLambda ;
          
          float refPoint[3] ;
        
          refPoint[0] = x0 ;
          refPoint[1] = y0 ;
          refPoint[2] = z0 ;
          
          tpcTrack->setPhi( phi ) ;       
          tpcTrack->setD0( d0 ) ;
          tpcTrack->setZ0( z0 ) ;       
          tpcTrack->setReferencePoint( refPoint ) ;
          tpcTrack->setChi2(TkTeBank->getChi2(te)) ;
          tpcTrack->setNdf(TkTeBank->getNdf(te)) ;
          tpcTrack->setdEdx(TkTeBank->getDe_dx(te)) ;
          
          const vector <int> * hits ;
          vector<MCParticle*> mcPointers ;
          vector<int> mcHits ;

          hits = TkTeBank->getHitlist(te) ;

          //std::cout << "the number of the hits on TE = " << hits->size() << std::endl;
          
          for(unsigned int tehit=0; tehit<hits->size();tehit++){

            TrackerHit* trkHitTPC = dynamic_cast<TrackerHit*>( _usedCol->getElementAt( hits->at(tehit) ) ) ;
            //std::cout << hits->at(tehit) << " z = " << trkHitTPC->getPosition()[2] << std::endl;

            tpcTrack->addHit(trkHitTPC) ;
        
            for(unsigned int j=0; j<trkHitTPC->getRawHits().size(); j++){ 
          
              SimTrackerHit * simTrkHitTPC =dynamic_cast<SimTrackerHit*>(trkHitTPC->getRawHits().at(j)) ;
              MCParticle * mcp = dynamic_cast<MCParticle*>(simTrkHitTPC->getMCParticle()) ; 

              //            if(mcp == NULL) cout << "mc particle pointer = null" << endl ; 
              
              bool found = false;
          
              for(unsigned int k=0; k<mcPointers.size();k++)
                {
                  if(mcp==mcPointers[k]){
                    found=true;
                    mcHits[k]++;
                  }
                }
              if(!found){
                mcPointers.push_back(mcp);
                mcHits.push_back(1);
              }
            }
          }
          
          for(unsigned int k=0; k<mcPointers.size();k++){

            MCParticle * mcp = mcPointers[k];

            LCRelationImpl* tpclcRel = new LCRelationImpl;
            tpclcRel->setFrom (tpcTrack);
            tpclcRel->setTo (mcp);
            float weight = (float)(mcHits[k])/(float)(tpcTrack->getTrackerHits().size());
            //float weight = (float)(tpcTrack->getTrackerHits().size())/(float)mcHits[k];


            tpclcRel->setWeight(weight);
        

            _tpclcRelVec->addElement( tpclcRel );
          }
      
          //FIXME:SJA:  Covariance matrix not included yet needs converting for 1/R and TanLambda
      
        
          tpcTrack->subdetectorHitNumbers().resize(12);
          tpcTrack->subdetectorHitNumbers()[0] = int(0);
          tpcTrack->subdetectorHitNumbers()[1] = int(0);
          tpcTrack->subdetectorHitNumbers()[2] = int(0);
          tpcTrack->subdetectorHitNumbers()[3] = int(hits->size());
          tpcTrack->subdetectorHitNumbers()[4] = int(0);
          tpcTrack->subdetectorHitNumbers()[5] = int(0);
          tpcTrack->subdetectorHitNumbers()[6] = int(0);
          tpcTrack->subdetectorHitNumbers()[7] = int(0);
          tpcTrack->subdetectorHitNumbers()[8] = int(0);
          tpcTrack->subdetectorHitNumbers()[9] = int(hits->size());
          tpcTrack->subdetectorHitNumbers()[10] = int(0);
          tpcTrack->subdetectorHitNumbers()[11] = int(0);

          _tpcTrackVec->addElement( tpcTrack );
          
        }
      }
      TkMCBank->clear();
      TkHitBank->clear();
      TPCHitBank->clear();
      TkTeBank->clear();
      TkTkBank->clear();

    }
    delete TkMCBank;
    delete TPCHitBank;
    delete TkHitBank;
    delete TkTeBank;
    delete TkTkBank;
  
    evt->addCollection( _tpcTrackVec , _colNameTPCTracks) ;
    evt->addCollection( _tpclcRelVec , _colNameMCTPCTracksRel) ;
 
    delete _usedCol; _usedCol=NULL;
    delete _droppedCol; _droppedCol=NULL;
  }

  //******************************  

  _nEvt ++ ;
  _goodHits.clear();
  _savedGoodHits.clear();
  
}



void LEPTrackingProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void LEPTrackingProcessor::end(){ 

  if(_histograms){
    TFile *hfile = new TFile("tpcTiming.root","recreate");
    fTPCR->TH1F::Write();
    fTPCZ->TH1F::Write();
    fTPCRZ->TH2F::Write();
    fTPCXY->TH2F::Write();
    fTPCRRaw->TH1F::Write();
    fTPCZRaw->TH1F::Write();
    fTPCRZRaw->TH2F::Write();
    fTPCXYRaw->TH2F::Write();
    hfile->Close();
    delete hfile;
  }


//  std::cout << "LEPTrackingProcessor::end()  " << name() 
//            << " processed " << _nEvt << " events in " << _nRun << " runs "
//            << std::endl ;
//
}

void LEPTrackingProcessor::selectTPCHits(LCCollection* tpcTHcol, LCCollection* _usedCol, float zmin, float zmax){
  if( tpcTHcol != 0 ){
    
    int n_THits = tpcTHcol->getNumberOfElements()  ;
    for(int i=0; i< n_THits; i++){
      TrackerHit* THit = dynamic_cast<TrackerHit*>( tpcTHcol->getElementAt( i ) ) ;

      double *pos;
      pos = (double*) THit->getPosition(); 

      if(pos[2]>zmin && pos[2]<= zmax){
        _usedCol->addElement(THit) ;
        _goodHits.push_back(THit);
      }
    }
  }
}

void LEPTrackingProcessor::selectTPCHits(LCCollection* tpcTHcol, LCCollectionVec* _usedCol, LCCollectionVec* _droppedCol, int nbinHeight, int nbinWidth, float zmin, float zmax){

  if( tpcTHcol != 0 ){
    
    int n_THits = tpcTHcol->getNumberOfElements()  ;
    
    
    const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
    const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
    
    const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;
    
    double gearRMin = planeExt[0] ;
    double gearRMax = planeExt[1] ;
    
    const gear::Vector2D padCoord = padLayout.getPadCenter(1) ;

    double binHeight = padLayout.getRowHeight(1) * (double) nbinHeight ;
    double binWidth = padLayout.getPadWidth(1) * padCoord[0] * (double) nbinWidth ;
        

    // create bined pad layout
    const gear::FixedPadSizeDiskLayout padsAsBins(gearRMin, gearRMax, binHeight, binWidth) ;                                       

    // create hit map
    HitMap hitMap ; 
    

    for(int i=0; i< n_THits; i++){
      
      TrackerHit* THit = dynamic_cast<TrackerHit*>( tpcTHcol->getElementAt( i ) ) ;
      
      double *pos;
      pos = (double*) THit->getPosition(); 
      
      double rad = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);      
      if(_histograms){
        if(fabs(pos[2])>0.1){
          fTPCRRaw->Fill(rad);
          fTPCZRaw->Fill(pos[2]);   
          fTPCRZRaw->Fill(pos[2],rad);
          fTPCXYRaw->Fill(pos[0],pos[1]);
        }
      }


      //get phi of current hit
      float phi = atan2(pos[1],pos[0]);
      if (phi<0.) phi=phi+twopi;     
      
      int padIndex = padsAsBins.getNearestPad(rad,phi);
      unsigned int iRowHit = padsAsBins.getRowNumber(padIndex);      
      unsigned int iPhiHit = padsAsBins.getPadNumber(padIndex);
      
      // enter hit into hitMap
      if(pos[2]>zmin && pos[2]<= zmax){
        hitMap[  make_keyNew( iRowHit, iPhiHit ) ].push_back(  THit ) ;      
      }
    }

    //loop over hitmap and fill both collections of cut and remaining hits

    for( HitMap::iterator it = hitMap.begin() ;it != hitMap.end() ;  ++it ) {
      
      const std::vector<EVENT::TrackerHit*>& v = it->second ;
      
      if(   v.size() >=  (unsigned) _multiplicityCut ) {

        for( unsigned i = 0 ; i < v.size() ; i++){
          _droppedCol->addElement(  v[i] ) ;
        } 

        //                 int color = 0x88ff88 ;
        //                 int layer = 6 ;
//                 int marker = 2 ;
//                 int size = 1 ;
        
        //                MarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer) ;  
        
      }

      if(   v.size() < (unsigned) _multiplicityCut  ) {
        
        for( unsigned i = 0 ; i < v.size() ; i++){
          _usedCol->addElement(  v[i] ) ;

          _goodHits.push_back(v[i]);
          
        }            
        
//                  int color = 0x88ffff ;
//                  int layer = 7 ;
//                  int marker = 2 ;
//                  int size = 1 ;
        
        //                 MarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer) ; 
      }
      
    }   
        
  }

}


void LEPTrackingProcessor::FillTPCHitBanks(){

  TPCHitBank->clear();
  TkHitBank->clear();

  TkHitBank->setFirstHitIndex("TPC");

  for(int i=0; i< _goodHits.size(); ++i){
      
    TrackerHit* trkHitTPC = _goodHits[i];
    
    double *pos;
    float  edep;
    float  time;
    
    //      cellId = 	trkHitTPC->getCellID0();
    pos = (double*) trkHitTPC->getPosition(); 
    edep = trkHitTPC->getEDep();
    time = trkHitTPC->getTime();
    
      // convert to cm needed for BRAHMS(GEANT)
    float x = 0.1*pos[0];
    float y = 0.1*pos[1];
    float z = 0.1*pos[2];
    
    // convert de/dx from GeV (LCIO) to number of electrons 
    const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
    //    double tpcIonisationPotential = gearTPC.getDoubleVal("tpcIonPotential");
    //    edep = edep/tpcIonisationPotential;
    
//       double tpcRPhiResConst = gearTPC.getDoubleVal("tpcRPhiResConst");
//       double tpcRPhiResDiff  = gearTPC.getDoubleVal("tpcRPhiResDiff");
//       double aReso = tpcRPhiResConst*tpcRPhiResConst;
//       double driftLenght = gearTPC.getMaxDriftLength() - fabs(pos[2]);
//       if (driftLenght <0) { 
//         driftLenght = 0.10;
//       }
//       double bReso = tpcRPhiResDiff*tpcRPhiResDiff;
//       double tpcRPhiRes = sqrt(aReso + bReso*driftLenght);
//       double tpcZRes = gearTPC.getDoubleVal("tpcZRes");
    

    // Covariance Matrix in LCIO is defined in XYZ convert to R-Phi-Z
    // For no error in r
    
    
    double rSqrd = pos[0]*pos[0] + pos[1]*pos[1];

    float r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);


    double phi = atan2(pos[1],pos[0]); 
    double tpcRPhiRes = sqrt(trkHitTPC->getCovMatrix()[0] + trkHitTPC->getCovMatrix()[2]);
    double tpcZRes = sqrt(trkHitTPC->getCovMatrix()[5]);
    
      //      cout << "row_hits->getY() = " <<  pos[1] << "  row_hits->getY() = " << pos[0] ; 
//      cout << "  phi = " <<  phi ;  
//      cout << "  tpcRPhiRes = " << tpcRPhiRes ;
//      cout << "  cos(phi)*cos(phi)*tpcRPhiRes*tpcRPhiRes = " << trkHitTPC->getCovMatrix()[2] << endl;

    // convert mm to cm 
    tpcRPhiRes = 0.1 * tpcRPhiRes;
    tpcZRes = 0.1 * tpcZRes;
    
    
    // Brahms resolution code for TPC = 3 REF tkhtpc.F
    int icode = 3;
    int subid = 500;
    
    int mctrack = 0;
    
    if(_histograms){
      fTPCR->Fill(r);
      fTPCZ->Fill(pos[2]);
      fTPCRZ->Fill(pos[2],r);
      fTPCXY->Fill(pos[0],pos[1]);
    }

    TkHitBank->add_hit(x,y,z,edep,subid,mctrack,0,0,icode,tpcRPhiRes,tpcZRes);
    
    TPCHitBank->add_hit(x,y,z,edep,subid,tpcRPhiRes,tpcZRes,mctrack);
    
  }
  
  TkHitBank->setLastHitIndex("TPC"); 

  if(TkHitBank->getNumOfSubDetHits("TPC") > 0) {
    int tpcsubid = TkHitBank->getSubdetectorID(TkHitBank->getFirstHitIndex("TPC")) ;
    streamlog_out(WARNING) << "the first hit for the TPC has id " << tpcsubid << endl ;
    tpcsubid = TkHitBank->getSubdetectorID(TkHitBank->getLastHitIndex("TPC")) ;
    streamlog_out(WARNING) << "the last hit for the TPC has id " << tpcsubid << endl ;
  }
}

