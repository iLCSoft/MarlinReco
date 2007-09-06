/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
** This file is part of the MarlinReco Project.
** Forming part of the SubPackage: BrahmsTracking.
**
** For the latest version download from Web CVS:
** www.blah.de
**
** $Id: LEPTrackingProcessor.cc,v 1.25 2007-09-06 13:24:04 harderk Exp $
**
** $Log: not supported by cvs2svn $
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
#include <gear/BField.h>
//



PROTOCCALLSFFUN0(INT,TKTREV,tktrev)
#define TKTREV() CCALLSFFUN0(TKTREV,tktrev)

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
                  float* maxdrift, float* tpcpixz, float* ionpoten, float* tpcrpres, float* tpczres, float* tpcbfield){

  //  try{

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;
  
  *innerrad = 0.1 * float( planeExt[0] ) ;
  *outerrad = 0.1 *float( planeExt[1] ) ;
  *npadrows = padLayout.getNRows() ;
  *maxdrift = 0.1 * float( gearTPC.getMaxDriftLength() );
  *tpcpixz = 0.1 * float(gearTPC.getDoubleVal("tpcPixZ")) ;
  *ionpoten = 0.1 * float(gearTPC.getDoubleVal("tpcIonPotential")) ;  
  *tpcrpres = 0.1 * float(gearTPC.getDoubleVal("tpcRPhiResConst")) ;  
  *tpczres = 0.1 * float(gearTPC.getDoubleVal("tpcZRes")) ;
  *tpcbfield = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z();

  //  }
//  catch() {return 1} ;

  return 0;
}

FCALLSCFUN9(INT,gettpcgeom,GETTPCGEOM,gettpcgeom, PFLOAT, PFLOAT, PINT, 
            PFLOAT, PFLOAT, PFLOAT, PFLOAT, PFLOAT, PFLOAT )




  // end of cfortran.h definitions

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

  registerInputCollection( LCIO::TRACKERHIT,
                           "VTXTrackerHitCollectionName" , 
                           "Name of the VTX TrackerHit collection"  ,
                           _colNameVTX ,
                           std::string("VTXTrackerHits") ) ;
  
  registerInputCollection( LCIO::TRACKERHIT,
                           "SITTrackerHitCollectionName" , 
                           "Name of the SIT TrackerHit collection"  ,
                           _colNameSIT ,
                           std::string("SITTrackerHits") ) ;
  
  registerOutputCollection( LCIO::TRACK,
                            "TPCTrackCollectionName" , 
                            "Name of the TPC Track collection"  ,
                            _colNameTPCTracks ,
                            std::string("TPCTracks") ) ;

  registerOutputCollection( LCIO::TRACK,
                            "TrackCollectionName" , 
                            "Name of the Track collection"  ,
                            _colNameTracks ,
                            std::string("Tracks") ) ;
  
  registerOutputCollection( LCIO::LCRELATION,
                            "MCTPCTrackRelCollectionName" , 
                            "Name of the TPC Track MC Relation collection"  ,
                            _colNameMCTPCTracksRel ,
                            std::string("TPCTracksMCP") ) ;
  
  registerOutputCollection( LCIO::LCRELATION,
                            "MCTrackRelCollectionName" , 
                            "Name of the Track MC Relation collection"  ,
                            _colNameMCTracksRel ,
                            std::string("TracksMCP") ) ;
}


void LEPTrackingProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void LEPTrackingProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void LEPTrackingProcessor::processEvent( LCEvent * evt ) { 

  static bool firstEvent = true ;
    
  // this gets called for every event 
  // usually the working horse ...
  
//   int skipToEvent = 3 ;
//   if(_nEvt<skipToEvent) {
//     ++_nEvt;
//     return;
//   }

  if(firstEvent==true) cout << "LEPTrackingProcessor called for first event" << endl;

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

  
  if( tpcTHcol != 0 ){
    
    LCCollectionVec* tpcTrackVec = new LCCollectionVec( LCIO::TRACK )  ;
    LCCollectionVec* TrackVec = new LCCollectionVec( LCIO::TRACK )  ;

    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit( LCIO::TRBIT_HITS ) ;
    tpcTrackVec->setFlag( trkFlag.getFlag()  ) ;
    TrackVec->setFlag( trkFlag.getFlag()  ) ;

    LCCollectionVec* lcRelVec = new LCCollectionVec( LCIO::LCRELATION )  ;
    LCCollectionVec* tpclcRelVec = new LCCollectionVec( LCIO::LCRELATION )  ;

    int nTPCHits = tpcTHcol->getNumberOfElements()  ;   
    
    TkHitBank->setFirstHitIndex("TPC"); 
    
    for(int i=0; i< nTPCHits; i++){
      
      TrackerHit* trkHitTPC = dynamic_cast<TrackerHit*>( tpcTHcol->getElementAt( i ) ) ;
      
      double *pos;
      float  de_dx;
      float  time;
  
      //      cellId = 	trkHitTPC->getCellID();
      pos = (double*) trkHitTPC->getPosition(); 
      de_dx = trkHitTPC->getdEdx();
      time = trkHitTPC->getTime();
      
      // convert to cm needed for BRAHMS(GEANT)
      float x = 0.1*pos[0];
      float y = 0.1*pos[1];
      float z = 0.1*pos[2];
      
      // convert de/dx from GeV (LCIO) to number of electrons 
      
      double tpcIonisationPotential = gearTPC.getDoubleVal("tpcIonPotential");
      de_dx = de_dx/tpcIonisationPotential;

      double tpcRPhiResConst = gearTPC.getDoubleVal("tpcRPhiResConst");
      double tpcRPhiResDiff  = gearTPC.getDoubleVal("tpcRPhiResDiff");
      double aReso = tpcRPhiResConst*tpcRPhiResConst;
      double driftLenght = gearTPC.getMaxDriftLength() - fabs(pos[2]);
      if (driftLenght <0) { 
        driftLenght = 0.10;
      }
      double bReso = tpcRPhiResDiff*tpcRPhiResDiff;
      double tpcRPhiRes = sqrt(aReso + bReso*driftLenght);
      double tpcZRes = gearTPC.getDoubleVal("tpcZRes");

      tpcRPhiRes = 0.1 * tpcRPhiRes;
      tpcZRes = 0.1 * tpcZRes;


      // Brahms resolution code for TPC = 3 REF tkhtpc.F
      int icode = 3;
      int subid = 500;

      int mctrack = 0;
      

      TkHitBank->add_hit(x,y,z,de_dx,subid,mctrack,0,0,icode,tpcRPhiRes,tpcZRes);
      


      TPCHitBank->add_hit(x,y,z,de_dx,subid,tpcRPhiRes,tpcZRes,mctrack);

    }

    TkHitBank->setLastHitIndex("TPC"); 
    
    //    cout << "the number of tpc hits sent to brahms = " << TPCHitBank->size() << endl;
    //    CNTPC.ntphits = TPCHitBank->size();
 
    //_____________________________________________________________________

    if( vtxTHcol != 0 ) { 

      int nVTXHits = vtxTHcol->getNumberOfElements()  ;   
      
      TkHitBank->setFirstHitIndex("VTX"); 
      
      for(int i=0; i< nVTXHits; i++){
        
        TrackerHit* trkHitVTX = dynamic_cast<TrackerHit*>( vtxTHcol->getElementAt( i ) ) ;
    
        double *pos;
        float  de_dx;
        float  time;
        
        //      cellId = 	trkHitVTX->getCellID();
        pos = (double*) trkHitVTX->getPosition(); 
        de_dx = trkHitVTX->getdEdx() ;
        time = trkHitVTX->getTime() ;
        
        // convert to cm needed for BRAHMS(GEANT)
        float x = 0.1*pos[0] ;
        float y = 0.1*pos[1] ;
        float z = 0.1*pos[2] ;
        
        
        // Brahms resolution code for VTX = 3 REF tkhtpc.F
        
        int subid = trkHitVTX->getType() ;
        
        // brsimu/brgeom/brtrac/code_f/vxpgeom.F:      VXDPPNT=7.0E-4
        float vtxRes = 0.007 ;
        int resCode = 3 ;
        
        int mctrack = 0 ;
        
        TkHitBank->add_hit(x,y,z,de_dx,subid,mctrack,0,0,resCode,vtxRes,vtxRes) ;
        
      }
      
      TkHitBank->setLastHitIndex("VTX"); 
      
      
    }
    
    if( sitTHcol != 0 ) { 

      int nSITHits = sitTHcol->getNumberOfElements()  ;   
      
      TkHitBank->setFirstHitIndex("SIT"); 
      
      for(int i=0; i< nSITHits; i++){
        
        TrackerHit* trkHitSIT = dynamic_cast<TrackerHit*>( sitTHcol->getElementAt( i ) ) ;
        
        double *pos;
        float  de_dx;
        float  time;
        
        //      cellId = 	trkHitSIT->getCellID();
        pos = (double*) trkHitSIT->getPosition(); 
        de_dx = trkHitSIT->getdEdx() ;
        time = trkHitSIT->getTime() ;
        
        // convert to cm needed for BRAHMS(GEANT)
        float x = 0.1*pos[0] ;
        float y = 0.1*pos[1] ;
        float z = 0.1*pos[2] ;
        
        
        // Brahms resolution code for SIT = 3 REF tkhtpc.F
        
        int subid = trkHitSIT->getType() ;
        
        // brsimu/brgeom/brtrac/code_f/vxpgeom.F:      VXDPPNT=7.0E-4
        float sitRPhiRes = 0.01 ;
        float sitZRes = 0.05 ;
        int resCode = 3 ;
        
        int mctrack = 0 ;
        
        TkHitBank->add_hit(x,y,z,de_dx,subid,mctrack,0,0,resCode,sitRPhiRes,sitZRes) ;
        
      }
      
      TkHitBank->setLastHitIndex("SIT"); 
      
    }
        
      
      //_____________________________________________________________________
      
    if(TkHitBank->getNumOfSubDetHits("TPC") > 0) {
      int tpcsubid = TkHitBank->getSubdetectorID(TkHitBank->getFirstHitIndex("TPC")) ;
      cout << "the first hit for the TPC has id " << tpcsubid << endl ;
    }

    if(TkHitBank->getNumOfSubDetHits("VTX") > 0) {
      int vtxsubid = TkHitBank->getSubdetectorID(TkHitBank->getFirstHitIndex("VTX")) ;
      cout << "the first hit for the vtx has id " << vtxsubid << endl ;
    }
    
    if(TkHitBank->getNumOfSubDetHits("SIT") > 0) {
      int sitsubid = TkHitBank->getSubdetectorID(TkHitBank->getFirstHitIndex("SIT")) ;
      cout << "the first hit for the sit has id " << sitsubid << endl ;
    }
    
    int errTKTREV = TKTREV(); 

    cout << "TKTREV returns:" << errTKTREV << endl;

    
    if(errTKTREV!=0) cout << "have you set the ionisation potential correctly in the gear xml file" << endl;    
    cout << "number of TE's = " << TkTeBank->size() << endl ;

    cout << "number of TK's = " << TkTkBank->size() << endl ;



    for(int te=0; te<TkTeBank->size();te++){


      if( TkTeBank->getSubdetector_ID(te)==500 ) {

        TrackImpl* tpcTrack = new TrackImpl ; 
      
        const double ref_r = 10.*TkTeBank->getCoord1_of_ref_point(te);
        const double ref_phi =TkTeBank->getCoord2_of_ref_point(te)/TkTeBank->getCoord1_of_ref_point(te);
        const double ref_z = 10.*TkTeBank->getCoord3_of_ref_point(te);

        //         cout << "ref_r = " << ref_r << endl;
        //         cout << "ref_phi = " << ref_phi << endl;
      
        //FIXME:SJA: B-field hard coded needs to redeemed
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
        else if ( phiOfPCA > twopi/2. ) phiOfPCA = phiOfPCA = twopi + phiOfPCA ;   
        
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

//         std::cout << "calc value of omega = " << omega;
//         std::cout << " calc value of phi = " << phi;
//         std::cout << " calc value of d0 = " << d0;
//         std::cout << " calc value of z0 = " << z0;
//         std::cout << " calc value of tanLambda = " << tanLambda<< std::endl;

        ////////////////////////////////

        tpcTrack->setIsReferencePointPCA(true) ;
        tpcTrack->setChi2(TkTeBank->getChi2(te)) ;
        tpcTrack->setNdf(TkTeBank->getNdf(te)) ;
        tpcTrack->setdEdx(TkTeBank->getDe_dx(te)) ;

        const vector <int> * hits ;
        vector<MCParticle*> mcPointers ;
        vector<int> mcHits ;

        hits = TkTeBank->getHitlist(te) ;

        //        std::cout << "the number of the hits on TE = " << hits->size() << std::endl;

        for(unsigned int tehit=0; tehit<hits->size();tehit++){

          TrackerHit* trkHitTPC = dynamic_cast<TrackerHit*>( tpcTHcol->getElementAt( hits->at(tehit) ) ) ;

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
        

          tpclcRelVec->addElement( tpclcRel );
        }
      
        //FIXME:SJA:  Covariance matrix not included yet needs converting for 1/R and TanLambda
      
        
        tpcTrack->subdetectorHitNumbers().resize(8);
        tpcTrack->subdetectorHitNumbers()[0] = int(0);
        tpcTrack->subdetectorHitNumbers()[1] = int(0);
        tpcTrack->subdetectorHitNumbers()[2] = int(0);
        tpcTrack->subdetectorHitNumbers()[3] = int(hits->size());
        tpcTrack->subdetectorHitNumbers()[4] = int(0);
        tpcTrack->subdetectorHitNumbers()[5] = int(0);
        tpcTrack->subdetectorHitNumbers()[6] = int(0);
        tpcTrack->subdetectorHitNumbers()[7] = int(hits->size());


        tpcTrackVec->addElement( tpcTrack );

      
        //      cout << "TkTeBank->getSubdetector_ID(te) = " << TkTeBank->getSubdetector_ID(te) << endl;
        //      cout << "TkTeBank->getSubmodule(te) = " << TkTeBank->getSubmodule(te) << endl;
        //      cout << "TkTeBank->getUnused(te) = " << TkTeBank->getUnused(te) << endl;
        //      cout << "TkTeBank->getMeasurement_code(te) = " << TkTeBank->getMeasurement_code(te) << endl;
        //      cout << "TkTeBank->getPointer_to_end_of_TE(te) = " << TkTeBank->getPointer_to_end_of_TE(te) << endl;
        //      cout << "TkTeBank->getNdf(te) = " << TkTeBank->getNdf(te) << endl;
        //      cout << "TkTeBank->getChi2(te) = " << TkTeBank->getChi2(te) << endl;
        //      cout << "TkTeBank->getLength(te) = " << TkTeBank->getLength(te) << endl;
        //      cout << "TkTeBank->getCoord1_of_ref_point(te) = " << TkTeBank->getCoord1_of_ref_point(te) << endl;
        //      cout << "TkTeBank->getCoord2_of_ref_point(te) = " << TkTeBank->getCoord2_of_ref_point(te) << endl;
        //      cout << "TkTeBank->getCoord3_of_ref_point(te) = " << TkTeBank->getCoord3_of_ref_point(te) << endl;
        //      cout << "TkTeBank->getTheta(te) = " << TkTeBank->getTheta(te) << endl;
        //      cout << "TkTeBank->getPhi(te) = " << TkTeBank->getPhi(te) << endl;
        //      cout << "TkTeBank->getInvp(te) = " << TkTeBank->getInvp(te) << endl;
        //      cout << "TkTeBank->getDe_dx(te) = " << TkTeBank->getDe_dx(te) << endl;
        //      cout << "TkTeBank->getCovmatrix1(te) = " << TkTeBank->getCovmatrix1(te) << endl;
        //      cout << "TkTeBank->getCovmatrix2(te) = " << TkTeBank->getCovmatrix2(te) << endl;
        //      cout << "TkTeBank->getCovmatrix3(te) = " << TkTeBank->getCovmatrix3(te) << endl;
        //      cout << "TkTeBank->getCovmatrix4(te) = " << TkTeBank->getCovmatrix4(te) << endl;
        //      cout << "TkTeBank->getCovmatrix5(te) = " << TkTeBank->getCovmatrix5(te) << endl;
        //      cout << "TkTeBank->getCovmatrix6(te) = " << TkTeBank->getCovmatrix6(te) << endl;
        //      cout << "TkTeBank->getCovmatrix7(te) = " << TkTeBank->getCovmatrix7(te) << endl;
        //      cout << "TkTeBank->getCovmatrix8(te) = " << TkTeBank->getCovmatrix8(te) << endl;
        //      cout << "TkTeBank->getCovmatrix9(te) = " << TkTeBank->getCovmatrix9(te) << endl;
        //      cout << "TkTeBank->getCovmatrix0(te) = " << TkTeBank->getCovmatrix10(te) << endl;
        //      cout << "TkTeBank->getCovmatrix1(te) = " << TkTeBank->getCovmatrix11(te) << endl;
        //      cout << "TkTeBank->getCovmatrix2(te) = " << TkTeBank->getCovmatrix12(te) << endl;
        //      cout << "TkTeBank->getCovmatrix3(te) = " << TkTeBank->getCovmatrix13(te) << endl;
        //      cout << "TkTeBank->getCovmatrix4(te) = " << TkTeBank->getCovmatrix14(te) << endl;
        //      cout << "TkTeBank->getCovmatrix5(te) = " << TkTeBank->getCovmatrix15(te) << endl;
        //
      }
    }
    // set the parameters to decode the type information in the collection
    // for the time being this has to be done manually
    // in the future we should provide a more convenient mechanism to 
    // decode this sort of meta information

    //     StringVec typeNames ;
    //     IntVec typeValues ;
    //     typeNames.push_back( LCIO::TRACK ) ;
    //     typeValues.push_back( 1 ) ;
    //     tpcTrackVec->parameters().setValues("TrackTypeNames" , typeNames ) ;
    //     tpcTrackVec->parameters().setValues("TrackTypeValues" , typeValues ) ;

    //    evt->addCollection( tpcTrackVec , _colNameTPCTracks) ;
    //    evt->addCollection( lcRelVec , _colNameMCTracksRel) ;


  //******************************
  // try here for TK's

    
    for(int tk=0; tk<TkTkBank->size();tk++){

      TrackImpl* Track = new TrackImpl ; 
        
      const double ref_r = 10.*TkTkBank->getCoord1_of_ref_point(tk);
      const double ref_phi =TkTkBank->getCoord2_of_ref_point(tk)/TkTkBank->getCoord1_of_ref_point(tk);
      const double ref_z = 10.*TkTkBank->getCoord3_of_ref_point(tk);
      
      //         cout << "ref_r = " << ref_r << endl;
      //         cout << "ref_phi = " << ref_phi << endl;
      
      //FIXME:SJA: B-field hard coded needs to redeemed
      // transformation from 1/p to 1/R = consb * (1/p) / sin(theta)
      // consb is given by 1/R = (c*B)/(pt*10^9) where B is in T and pt in GeV  
      
      //      const double bField = gearTPC.getDoubleVal("BField") ;
      const double bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
      const double consb = (2.99792458* bField )/(10*1000.) ;     // divide by 1000 m->mm
      
      // tan lambda and curvature remain unchanged as the track is only extrapolated
      // set negative as 1/p is signed with geometric curvature clockwise negative
      
      const double omega = ( -consb*TkTkBank->getInvp(tk) )/ sin( TkTkBank->getTheta(tk) ) ;
      const double tanLambda = tan( (twopi/4.) - TkTkBank->getTheta(tk) ) ;
      
      Track->setOmega( omega ) ;
      Track->setTanLambda( tanLambda ) ;      
      
      // computation of D0 and Z0 taken from fkrtpe.F in Brahms
      
      // xref and yref of ref point 
      const double xref = ref_r*cos(ref_phi) ;
      const double yref = ref_r*sin(ref_phi) ;
      const double zref = ref_z ; 
      const double trkRadius = 1. / omega ;
      
      
      ////////////////////////////////
      
      // center of circumference
      const double xc = xref + trkRadius * sin( TkTkBank->getPhi(tk) ) ;
      const double yc = yref - trkRadius * cos( TkTkBank->getPhi(tk) ) ;
      
      const double xc2 = xc * xc ; 
      const double yc2 = yc * yc ;
      
      const double DCA = ( sqrt( xc2+yc2 ) - fabs(trkRadius) ) ;
      
      double phiOfPCA = atan2 (yc,xc) ;
      
      if (DCA<0.) phiOfPCA = phiOfPCA + twopi/2. ;
      
      if ( phiOfPCA < -twopi/2. ) phiOfPCA = twopi + phiOfPCA ;
      else if ( phiOfPCA > twopi/2. ) phiOfPCA = phiOfPCA = twopi + phiOfPCA ;   
      
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
      
      Track->setPhi( phi ) ;       
      Track->setD0( d0 ) ;
      Track->setZ0( z0 ) ;       
      Track->setReferencePoint( refPoint ) ;
      
      //                 std::cout << "calc value of omega = " << omega;
      //                 std::cout << " calc value of phi = " << phi;
      //                 std::cout << " calc value of d0 = " << d0;
      //                 std::cout << " calc value of z0 = " << z0;
      //                 std::cout << " calc value of tanLambda = " << tanLambda<< std::endl;
      //
      ////////////////////////////////
      
      Track->setIsReferencePointPCA(true) ;
      Track->setChi2(TkTkBank->getChi2(tk)) ;
      Track->setNdf(TkTkBank->getNdf(tk)) ;
      //        Track->setdEdx(TkTkBank->getDe_dx(tk)) ;
      
      const vector <int> * hits ;
      const vector <int> * tes ;
      vector<MCParticle*> mcPointers ;
      vector<int> mcHits ;
      
      
      tes = TkTkBank->getTElist(tk) ;
      
      //      std::cout << "the number of TE's in TK " << tk << " = " << tes->size() << std::endl;   
      
      for(unsigned int tkte=0; tkte<tes->size();tkte++){
        
        
        //        std::cout << "the number of the TE in TK = " << tes->at(tkte) << std::endl;
        //        std::cout << "the Subdetector number of the TE in TK " << TkTeBank->getSubdetector_ID(tes->at(tkte)) << std::endl;
        
        hits = TkTeBank->getHitlist(tes->at(tkte)) ;
        
        //        std::cout << "the number of the hits on TE = " << hits->size() << std::endl;
        
        //          for(unsigned int tkhit=0; tkhit<hits->size();tkhit++){
        
        //           std::cout << "the Subdetector number of the Hit is " << TkHitBank->getSubdetectorID(hits->at(tkhit)) << std::endl; 
        //          }
        
        //        }

        for(unsigned int tkhit=0; tkhit<hits->size();tkhit++){
          
          TrackerHit* trkHit ;

          if(TkHitBank->getSubdetectorID(hits->at(tkhit))==500){
            
            int tpchitindex = hits->at(tkhit) - TkHitBank->getFirstHitIndex("TPC") ;
            trkHit = dynamic_cast<TrackerHit*>( tpcTHcol->getElementAt( tpchitindex ) ) ;
          }
          
          if(TkHitBank->getSubdetectorID(hits->at(tkhit))>99
             && TkHitBank->getSubdetectorID(hits->at(tkhit))<106){
           
            int vtxhitindex = hits->at(tkhit) - TkHitBank->getFirstHitIndex("VTX") ;  
            trkHit = dynamic_cast<TrackerHit*>( vtxTHcol->getElementAt( vtxhitindex ) ) ;
          }          

          if(TkHitBank->getSubdetectorID(hits->at(tkhit))>399
             && TkHitBank->getSubdetectorID(hits->at(tkhit))<403){
           
            int sithitindex = hits->at(tkhit) - TkHitBank->getFirstHitIndex("SIT") ;  
            trkHit = dynamic_cast<TrackerHit*>( sitTHcol->getElementAt( sithitindex ) ) ;
          }

          Track->addHit(trkHit) ;
          
          for(unsigned int j=0; j<trkHit->getRawHits().size(); j++){ 
            
            SimTrackerHit * simTrkHit =dynamic_cast<SimTrackerHit*>(trkHit->getRawHits().at(j)) ;
            MCParticle * mcp = dynamic_cast<MCParticle*>(simTrkHit->getMCParticle()) ; 

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
      }

      for(unsigned int k=0; k<mcPointers.size();k++){
        
        MCParticle * mcp = mcPointers[k];
        
        LCRelationImpl* lcRel = new LCRelationImpl;
        lcRel->setFrom (Track);
        lcRel->setTo (mcp);
        float weight = (float)(mcHits[k])/(float)(Track->getTrackerHits().size());
        //float weight = (float)(Track->getTrackerHits().size())/(float)mcHits[k];
        
        // debug
        /*
        std::cout << "TkTkBank->size() : " << TkTkBank->size() << " Track : " << tk 
                  << "  # MCs : " << mcPointers.size() 
                  << "  actual : " << k << "  # TrackerHits : " 
                  << Track->getTrackerHits().size();
        std::cout << "   mcHits[" << k << "] = " << mcHits[k];
        std::cout << "   LEPTR WEIGHT: " 
                  << weight << "  mcp-> " << mcp->getPDG() << " energy = " << mcp->getEnergy() << std::endl;
        */
        
        
        lcRel->setWeight(weight);
        
        
        lcRelVec->addElement( lcRel );
      }
      

        
      //FIXME:SJA:  Covariance matrix not included yet needs converting for 1/R and TanLambda
      
      
      TrackVec->addElement( Track );

    }

    //    // set the parameters to decode the type information in the collection
    //    // for the time being this has to be done manually
    //    // in the future we should provide a more convenient mechanism to 
    //    // decode this sort of meta information
    //
    //    //     StringVec typeNames ;
    //    //     IntVec typeValues ;
    //    //     typeNames.push_back( LCIO::TRACK ) ;
    //    //     typeValues.push_back( 1 ) ;
    //    //     TrackVec->parameters().setValues("TrackTypeNames" , typeNames ) ;
    //    //     TrackVec->parameters().setValues("TrackTypeValues" , typeValues ) ;
    //
    
    evt->addCollection( tpcTrackVec , _colNameTPCTracks) ;
    evt->addCollection( tpclcRelVec , _colNameMCTPCTracksRel) ;
    evt->addCollection( TrackVec , _colNameTracks) ;
    evt->addCollection( lcRelVec , _colNameMCTracksRel) ;
    

  }
  

  //******************************  

  delete TkMCBank;
  delete TPCHitBank;
  delete TkHitBank;
  delete TkTeBank;
  delete TkTkBank;

  _nEvt ++ ;
  
}



void LEPTrackingProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void LEPTrackingProcessor::end(){ 
  
  //   std::cout << "LEPTrackingProcessor::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;

}

