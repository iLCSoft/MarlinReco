/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
** This file is part of the MarlinReco Project.
** Forming part of the SubPackage: BrahmsTracking.
**
** For the latest version download from Web CVS:
** www.blah.de
**
** $Id: LEPTrackingProcessor.cc,v 1.10 2005-08-03 19:05:24 aplin Exp $
**
** $Log: not supported by cvs2svn $ 
*/
#include "LEPTrackingProcessor.h"
#include <iostream>
#include <string>

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
#include <IMPL/TrackImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>

#include <cfortran.h>

#include"tpchitbank.h"
#include"tkhitbank.h"
#include"tktebank.h"
#include"tpc.h"
#include"marlin_tpcgeom.h"
#include"constants.h"

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
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

  int writetktecpp(float c, int a, int b); 

FCALLSCFUN3(INT,writetktecpp,WRITETKTECPP,writetktecpp, FLOAT, INT, INT)

  float rreadtktecpp(int a, int b); 

FCALLSCFUN2(FLOAT,rreadtktecpp,RREADTKTECPP,rreadtktecpp, INT, INT)

  int ireadtktecpp(int a, int b); 

FCALLSCFUN2(INT,ireadtktecpp,IREADTKTECPP,ireadtktecpp, INT, INT)

  int tkmktecpp(int subid,int submod,int unused,int MesrCode,int PnteTE,int Q,int ndf,float chi2,float L,float cord1,float cord2,float cord3,float theta,float phi,float invp,float dedx,float* cov);
  

FCALLSCFUN17(INT,tkmktecpp,TKMKTECPP,tkmktecpp, INT , INT ,INT ,INT ,INT ,INT ,INT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOATV)

  int addhittktecpp(int a, int b); 

FCALLSCFUN2(INT,addhittktecpp,ADDHITTKTECPP,addhittktecpp, INT, INT)

  int readtkitedatcpp(int a, int b); 

FCALLSCFUN2(INT,readtkitedatcpp,READTKITEDATCPP,readtkitedatcpp, INT, INT)

  int writetkitedatcpp(int c, int a, int b);

FCALLSCFUN3(INT,writetkitedatcpp,WRITETKITEDATCPP,writetkitedatcpp, INT, INT, INT)


  // end of cfortran.h definitions

  LEPTrackingProcessor aLEPTrackingProcessor ;

LEPTrackingProcessor::LEPTrackingProcessor() : Processor("LEPTrackingProcessor") {
  
  // modify processor description
  _description = "Produces Track collection from TPC TrackerHit collections using LEP tracking algorithms" ;

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "TPCTrackerHitCollectionName" , 
                              "Name of the TPC TrackerHit collection"  ,
                              _colNameTPC ,
                              std::string("TPCTrackerHits") ) ;

  registerProcessorParameter( "VTXTrackerHitCollectionName" , 
                              "Name of the VTX TrackerHit collection"  ,
                              _colNameVTX ,
                              std::string("VTXTrackerHits") ) ;

  registerProcessorParameter( "TPCTrackCollectionName" , 
                              "Name of the TPC Track collection"  ,
                              _colNameTPCTracks ,
                              std::string("TPCTracks") ) ;

  registerProcessorParameter( "MCTrackRelCollectionName" , 
                              "Name of the TPC Track MC Relation collection"  ,
                              _colNameMCTracksRel ,
                              std::string("MCTracksRel") ) ;
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

  if(firstEvent==true) cout << "LEPTrackingProcessor called for first event" << endl;

  firstEvent = false ;


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
  
  
  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;

  
  if( tpcTHcol != 0 ){
    
    LCCollectionVec* tpcTrackVec = new LCCollectionVec( LCIO::TRACK )  ;
    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit( LCIO::TRBIT_HITS ) ;
    tpcTrackVec->setFlag( trkFlag.getFlag()  ) ;

    LCCollectionVec* lcRelVec = new LCCollectionVec( LCIO::LCRELATION )  ;

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

      double tpcRPhiResMax = (gearTPC.getDoubleVal("tpcRPhiResMax"));
      double tpcRPhiRes = 0.1 * (tpcRPhiResMax - fabs(pos[2])/gearTPC.getMaxDriftLength()*0.01);
      double tpcZRes = 0.1 * (gearTPC.getDoubleVal("tpcZRes"));


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
        float vtxRes = 0.0007 ;
        int resCode = 3 ;
        
        int mctrack = 0 ;
        
        TkHitBank->add_hit(x,y,z,de_dx,subid,mctrack,0,0,resCode,vtxRes,vtxRes) ;
        
      }
      
      TkHitBank->setLastHitIndex("VTX"); 
      
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
    
    int errTKTREV = TKTREV(); 

    cout << "TKTREV returns:" << errTKTREV << endl;

    
    if(errTKTREV!=0) cout << "have you set the ionisation potential correctly in the gear xml file" << endl;    
    cout << "number of TE's = " << TkTeBank->size() << endl ;
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
      
        const double consb = (2.99792458*4.)/(10*1000.) ;     // divide by 1000 m->mm

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

        for(unsigned int tehit=0; tehit<hits->size();tehit++){

          TrackerHit* trkHitTPC = dynamic_cast<TrackerHit*>( tpcTHcol->getElementAt( hits->at(tehit) ) ) ;

          tpcTrack->addHit(trkHitTPC) ;
        
          for(unsigned int j=0; j<trkHitTPC->getRawHits().size(); j++){ 
          
            SimTrackerHit * simTrkHitTPC =dynamic_cast<SimTrackerHit*>(trkHitTPC->getRawHits().at(j)) ;
            MCParticle * mcp = dynamic_cast<MCParticle*>(simTrkHitTPC->getMCParticle()) ; 
            if(mcp == NULL) cout << "mc particle pointer = null" << endl ; 
          
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

          LCRelationImpl* lcRel = new LCRelationImpl;
          lcRel->setFrom (tpcTrack);
          lcRel->setTo (mcp);
          float weight = mcHits[k]/tpcTrack->getTrackerHits().size();
          lcRel->setWeight(weight);
        

          lcRelVec->addElement( lcRel );
        }
      
        //FIXME:SJA:  Covariance matrix not included yet needs converting for 1/R and TanLambda
      

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

    evt->addCollection( tpcTrackVec , _colNameTPCTracks) ;
    evt->addCollection( lcRelVec , _colNameMCTracksRel) ;
    

  }
  



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

