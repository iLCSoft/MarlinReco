#include "FortranProcessor.h"
#include <iostream>

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



PROTOCCALLSFFUN0(INT,TPCRUN,tpcrun)
#define TPCRUN() CCALLSFFUN0(TPCRUN,tpcrun)

  // FIXME:SJA: the namespace should be used explicitly
using namespace lcio ;
using namespace marlin ;
using namespace tpcgeom;
using namespace constants;

int writetpccpp(float c, int a, int b); 

FCALLSCFUN3(INT,writetpccpp,WRITETPCCPP,writetpccpp, FLOAT, INT, INT)

float readtpchitscpp(int a, int b); 

FCALLSCFUN2(FLOAT,readtpchitscpp,READTPCHITSCPP,readtpchitscpp, INT, INT)

int writetkhitcpp(float c, int a, int b); 

FCALLSCFUN3(INT,writetkhitcpp,WRITETKHITCPP,writetkhitcpp, FLOAT, INT, INT)

float readtkhitscpp(int a, int b); 

FCALLSCFUN2(FLOAT,readtkhitscpp,READTKHITSCPP,readtkhitscpp, INT, INT)

int writetktecpp(float c, int a, int b); 

FCALLSCFUN3(INT,writetktecpp,WRITETKTECPP,writetktecpp, FLOAT, INT, INT)

float readtktecpp(int a, int b); 

FCALLSCFUN2(FLOAT,readtktecpp,READTKTECPP,readtktecpp, INT, INT)

int tkmktecpp(int subid,int submod,int unused,int MesrCode,int PnteTE,int Q,int ndf,float chi2,float L,float cord1,float cord2,float cord3,float theta,float phi,float invp,float dedx,float cov[15]);
  

FCALLSCFUN17(INT,tkmktecpp,TKMKTECPP,tkmktecpp, INT , INT ,INT ,INT ,INT ,INT ,INT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOATV)

  int addhittktecpp(int a, int b); 

FCALLSCFUN2(INT,addhittktecpp,ADDHITTKTECPP,addhittktecpp, INT, INT)

FortranProcessor aFortranProcessor ;

FortranProcessor::FortranProcessor() : Processor("FortranProcessor") {
  
  // modify processor description
  _description = "Produces Track collection from TPC TrackerHit collections using LEP tracking algorithms" ;

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "CollectionName" , 
			      "Name of the TrackerHit collection"  ,
			      _colName ,
			      std::string("TPCTrackerHits") ) ;
}


void FortranProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void FortranProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void FortranProcessor::processEvent( LCEvent * evt ) { 

  static bool firstEvent = true ;
    
  // this gets called for every event 
  // usually the working horse ...

  if(firstEvent==true) std::cout << "FortranProcessor called for first event" << std::endl;

  firstEvent = false ;
  
  LCCollection* THcol = evt->getCollection( _colName ) ;

  if( THcol != 0 ){

    LCCollectionVec* TPC_TrackVec = new LCCollectionVec( LCIO::TRACK )  ;
    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit( LCIO::TRBIT_HITS ) ;
    TPC_TrackVec->setFlag( trkFlag.getFlag()  ) ;

    LCCollectionVec* LCRelVec = new LCCollectionVec( LCIO::LCRELATION )  ;

    int n_hits = THcol->getNumberOfElements()  ;   
    
    
    for(int i=0; i< n_hits; i++){
      
      TrackerHit* THit = dynamic_cast<TrackerHit*>( THcol->getElementAt( i ) ) ;
      
      int    cellId;
      double *pos;
      float  de_dx;
      float  time;
  
      //      cellId = 	THit->getCellID();
      pos = (double*) THit->getPosition(); 
      de_dx = THit->getdEdx();
      time = THit->getTime();
      
      // convert to cm needed for BRAHMS(GEANT)
      float x = 0.1*pos[0];
      float y = 0.1*pos[1];
      float z = 0.1*pos[2];
      
      // convert de/dx from GeV (LCIO) to number of electrons 

      de_dx = de_dx/ionisation_potential;

      float Rz   = 0.1*the_tpc->getTpcZRes();
      float Rrphi = 0.1*the_tpc->gettpc_rphi_res(pos[2]);
      float tpc_halfL = 0.1*the_tpc->getHalfLength();

      // Brahms resolution code for TPC = 3 REF tkhtpc.F
      int ICODE = 3;
      int SUBID = 500;

      int mctrack = 0;
      
      TkHitBank->add_hit(x,y,z,de_dx,SUBID,mctrack,0,0,ICODE,Rrphi,Rz);


      TPCHitBank->add_hit(x,y,z,de_dx,SUBID,Rrphi,Rz,mctrack);

    } 
    
    cout << "the number of tpc hits sent to brahms = " << TPCHitBank->size() << endl;
    CNTPC.ntphits = TPCHitBank->size();
 

   int error = TPCRUN();
    
    for(int te=0; te<TkTeBank->size();te++){

      TrackImpl* TPC_Track = new TrackImpl ; 

      double ref_r = 10.*TkTeBank->getCoord1_of_ref_point(te);
      double ref_phi =TkTeBank->getCoord2_of_ref_point(te)/TkTeBank->getCoord1_of_ref_point(te);
      double ref_z = 10.*TkTeBank->getCoord3_of_ref_point(te);

      //FIXME:SJA: B-field hard coded needs to redeemed
      // transformation from 1/p to 1/R = consb * (1/p) / sin(theta)
      // consb is given by 1/R = (c*B)/(pt*10^9) where B is in T and pt in GeV  

      double consb = (2.99792458*4.)/(10*1000.);    // divide by 1000 m->mm

      // computation of D0 and Z0 taken from fkrtpe.F in Brahms
       
      // x0 and y0 of ref point 
      double x0 = ref_r*cos(ref_phi);
      double y0 = ref_r*sin(ref_phi);
      
      // signed track radius
      double trk_radius = sin(TkTeBank->getTheta(te))/(consb*TkTeBank->getInvp(te));

      //      cout << "TkTeBank->getInvp(te) " << TkTeBank->getInvp(te) << endl;
      //      cout << " trk_radius = " << trk_radius << endl;
      
      // center of circumference
      double xc = x0 - trk_radius * sin(TkTeBank->getPhi(te));
      double yc = y0 + trk_radius * cos(TkTeBank->getPhi(te));
      
      // sign of geometric curvature: anti-clockwise == postive
      double geom_curvature = fabs(TkTeBank->getInvp(te))/TkTeBank->getInvp(te);

      //      cout << "geom_curvature = " << geom_curvature << endl;
      
      double xc2 = xc * xc;
      double yc2 = yc * yc;
      
      // Set D0
      TPC_Track->setD0( trk_radius - geom_curvature * sqrt(xc2+yc2));
 
      // Phi at D0
      double phiatD0 = atan2(yc,xc)+(twopi/2.)+geom_curvature*(twopi/4.);
      if (phiatD0<0.) phiatD0 = phiatD0 + twopi;
      if (phiatD0>twopi) phiatD0 = phiatD0 - twopi;

      //      cout << "phi at ref = " << ref_phi << endl;
      //      cout << "phi at D0 = " << phiatD0 << endl;
 
      // difference between phi at ref and D0 
      double dphi = fmod((phiatD0 - TkTeBank->getPhi(te) +  twopi + (twopi/2.) ) , twopi ) - (twopi/2.);
      
      //      cout << "dphi at D0 = " << dphi << endl;

      // signed length of arc
      double larc = trk_radius * dphi;

      // Set Z0
      TPC_Track->setZ0(ref_z+(1/(tan(TkTeBank->getTheta(te)))*larc)); 
      
      //      cout << " D0 = " << trk_radius - geom_curvature * sqrt(xc2+yc2) << endl;
      //      cout << " Z0 = " << ref_z+1/(tan(TkTeBank->getTheta(te)))*larc  << endl;

      //FIXME: SJA: phi is set at PCA not ref this should be resolved
      //      TPC_Track->setPhi(TkTeBank->getPhi(te));
      TPC_Track->setPhi(phiatD0);

      // tan lambda and curvature remain unchanged as the track is only extrapolated
      // set negative as 1/p is signed with geometric curvature clockwise negative
      TPC_Track->setOmega((-consb*TkTeBank->getInvp(te))/sin(TkTeBank->getTheta(te)));
      TPC_Track->setTanLambda(tan((twopi/4.)-TkTeBank->getTheta(te)));

      TPC_Track->setIsReferencePointPCA(false);
      TPC_Track->setChi2(TkTeBank->getChi2(te));
      TPC_Track->setNdf(TkTeBank->getNdf(te));
      TPC_Track->setdEdx(TkTeBank->getDe_dx(te));
     

      const vector <int> * hits;
      vector<MCParticle*> MCPointers;
      vector<int> MChits;

      hits = TkTeBank->getHitlist(te);

      for(int tehit=0; tehit<hits->size();tehit++){
	TrackerHit* THit = dynamic_cast<TrackerHit*>( THcol->getElementAt( hits->at(tehit) ) ) ;

	TPC_Track->addHit(THit);

	for(int j=0; j<THit->getRawHits().size(); j++){ 
	  
	  SimTrackerHit * STHit =dynamic_cast<SimTrackerHit*>(THit->getRawHits().at(j));
	  MCParticle * mcp = dynamic_cast<MCParticle*>(STHit->getMCParticle()); 
	  if(mcp == NULL) cout << "mc particle pointer = null" << endl; 
	  
	  bool found = false;
	  
	  for(int k=0; k<MCPointers.size();k++)
	    {
	      if(mcp==MCPointers[k]){
		found=true;
		MChits[k]++;
	      }
	    }
	  if(!found){
	    MCPointers.push_back(mcp);
	    MChits.push_back(1);
	  }
	}
      }
      
      for(int k=0; k<MCPointers.size();k++){

	MCParticle * mcp = MCPointers[k];
	
	LCRelationImpl* LCRel = new LCRelationImpl;
	LCRel->setFrom (TPC_Track);
	LCRel->setTo (mcp);
	float weight = MChits[k]/TPC_Track->getTrackerHits().size();
	LCRel->setWeight(weight);
	
	LCRelVec->addElement( LCRel );
      }
            
      //FIXME:SJA:  Covariance matrix not included yet needs converting for 1/R and TanLambda
      

      TPC_TrackVec->addElement( TPC_Track );

      
//      std::cout << "TkTeBank->getSubdetector_ID(te) = " << TkTeBank->getSubdetector_ID(te) << std::endl;
//      std::cout << "TkTeBank->getSubmodule(te) = " << TkTeBank->getSubmodule(te) << std::endl;
//      std::cout << "TkTeBank->getUnused(te) = " << TkTeBank->getUnused(te) << std::endl;
//      std::cout << "TkTeBank->getMeasurement_code(te) = " << TkTeBank->getMeasurement_code(te) << std::endl;
//      std::cout << "TkTeBank->getPointer_to_end_of_TE(te) = " << TkTeBank->getPointer_to_end_of_TE(te) << std::endl;
//      std::cout << "TkTeBank->getNdf(te) = " << TkTeBank->getNdf(te) << std::endl;
//      std::cout << "TkTeBank->getChi2(te) = " << TkTeBank->getChi2(te) << std::endl;
//      std::cout << "TkTeBank->getLength(te) = " << TkTeBank->getLength(te) << std::endl;
//      std::cout << "TkTeBank->getCoord1_of_ref_point(te) = " << TkTeBank->getCoord1_of_ref_point(te) << std::endl;
//      std::cout << "TkTeBank->getCoord2_of_ref_point(te) = " << TkTeBank->getCoord2_of_ref_point(te) << std::endl;
//      std::cout << "TkTeBank->getCoord3_of_ref_point(te) = " << TkTeBank->getCoord3_of_ref_point(te) << std::endl;
//      std::cout << "TkTeBank->getTheta(te) = " << TkTeBank->getTheta(te) << std::endl;
//      std::cout << "TkTeBank->getPhi(te) = " << TkTeBank->getPhi(te) << std::endl;
//      std::cout << "TkTeBank->getInvp(te) = " << TkTeBank->getInvp(te) << std::endl;
//      std::cout << "TkTeBank->getDe_dx(te) = " << TkTeBank->getDe_dx(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix1(te) = " << TkTeBank->getCovmatrix1(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix2(te) = " << TkTeBank->getCovmatrix2(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix3(te) = " << TkTeBank->getCovmatrix3(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix4(te) = " << TkTeBank->getCovmatrix4(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix5(te) = " << TkTeBank->getCovmatrix5(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix6(te) = " << TkTeBank->getCovmatrix6(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix7(te) = " << TkTeBank->getCovmatrix7(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix8(te) = " << TkTeBank->getCovmatrix8(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix9(te) = " << TkTeBank->getCovmatrix9(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix0(te) = " << TkTeBank->getCovmatrix10(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix1(te) = " << TkTeBank->getCovmatrix11(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix2(te) = " << TkTeBank->getCovmatrix12(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix3(te) = " << TkTeBank->getCovmatrix13(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix4(te) = " << TkTeBank->getCovmatrix14(te) << std::endl;
//      std::cout << "TkTeBank->getCovmatrix5(te) = " << TkTeBank->getCovmatrix15(te) << std::endl;
//
    }

    // set the parameters to decode the type information in the collection
    // for the time being this has to be done manually
    // in the future we should provide a more convenient mechanism to 
    // decode this sort of meta information

//     StringVec typeNames ;
//     IntVec typeValues ;
//     typeNames.push_back( LCIO::TRACK ) ;
//     typeValues.push_back( 1 ) ;
//     TPC_TrackVec->parameters().setValues("TrackTypeNames" , typeNames ) ;
//     TPC_TrackVec->parameters().setValues("TrackTypeValues" , typeValues ) ;
    
    evt->addCollection( TPC_TrackVec , "TPC_Tracks") ;
    evt->addCollection( LCRelVec , "MC_Track_Relations") ;
    

  }
  
  _nEvt ++ ;
  
}



void FortranProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void FortranProcessor::end(){ 
  
//   std::cout << "FortranProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

