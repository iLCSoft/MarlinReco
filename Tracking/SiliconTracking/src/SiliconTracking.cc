#include "SiliconTracking.h"
#include <iostream>

#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <iostream>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <marlin/Global.h>
#include "ClusterShapes.h"

#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/GearParameters.h>
#include <gear/VXDLayerLayout.h>

using namespace lcio ;
using namespace marlin ;

int min(int a,int b) {
  int result;
  if (a<b) {
    result = a;
  }
  else {
    result = b;
  }    
  return result;
}

int max(int a,int b) {
  int result;
  if (a>b) {
    result = a;
  }
  else {
    result = b;
  }    
  return result;
}

int abs(int x) {
  int result = x;
  if (x<0)
    result = -x;
  return result;
}

float FastTripletCheck(TrackerHitExtended * hitIn,
		       TrackerHitExtended * hitMiddle,
		       TrackerHitExtended * hitOut) {

  float _xHit[3];
  float _yHit[3];
  float _zHit[3];

  _xHit[0] = float(hitIn->getTrackerHit()->getPosition()[0]);
  _yHit[0] = float(hitIn->getTrackerHit()->getPosition()[1]);
  _zHit[0] = float(hitIn->getTrackerHit()->getPosition()[2]);
 
  _xHit[1] = float(hitMiddle->getTrackerHit()->getPosition()[0]);
  _yHit[1] = float(hitMiddle->getTrackerHit()->getPosition()[1]);
  _zHit[1] = float(hitMiddle->getTrackerHit()->getPosition()[2]);
 
  _xHit[2] = float(hitOut->getTrackerHit()->getPosition()[0]);
  _yHit[2] = float(hitOut->getTrackerHit()->getPosition()[1]);
  _zHit[2] = float(hitOut->getTrackerHit()->getPosition()[2]);
 
  float x0  = 0.5*(_xHit[1]+_xHit[0]);
  float y0  = 0.5*(_yHit[1]+_yHit[0]);
  float x0p = 0.5*(_xHit[2]+_xHit[1]);
  float y0p = 0.5*(_yHit[2]+_yHit[1]);
  float ax  = _yHit[1] - _yHit[0];
  float ay  = _xHit[0] - _xHit[1];
  float axp = _yHit[2] - _yHit[1];
  float ayp = _xHit[1] - _xHit[2];
  float det = ax * ayp - axp * ay;
  float time;

  if (det == 0.) {
      time = 500.;
  }
  else {
      gsl_matrix* A = gsl_matrix_alloc(2,2);
      gsl_vector* B = gsl_vector_alloc(2);
      gsl_vector* T = gsl_vector_alloc(2);     
      gsl_matrix_set(A,0,0,ax);
      gsl_matrix_set(A,0,1,-axp);
      gsl_matrix_set(A,1,0,ay);
      gsl_matrix_set(A,1,1,-ayp);
      gsl_vector_set(B,0,x0p-x0);
      gsl_vector_set(B,1,y0p-y0);
      gsl_linalg_HH_solve(A,B,T);
      time = gsl_vector_get(T,0); 
      gsl_matrix_free(A);
      gsl_vector_free(B);
      gsl_vector_free(T);
  }

  float X0 = x0 + ax*time;
  float Y0 = y0 + ay*time;

  float _pi = acos(-1.);

  float phi0 = (float)atan2(_yHit[0]-Y0,_xHit[0]-X0);
  float phi1 = (float)atan2(_yHit[1]-Y0,_xHit[1]-X0);
  float phi2 = (float)atan2(_yHit[2]-Y0,_xHit[2]-X0);


  if ( phi0 > phi1 ) 
      phi1 = phi1 + 2.0*_pi;
  if ( phi0 > phi2 )
      phi2 = phi2 + 2.0*_pi;
  if ( phi1 > phi2 )
      phi2 = phi2 + 2.0*_pi;


  float ZP = _zHit[0] + (phi1-phi0)*(_zHit[2]-_zHit[0])/(phi2-phi0);
  return ZP - _zHit[1];
    
  
}

extern "C" {
  void tfithl_(int & NPT, double * XF,double * YF, float * RF, float *PF, double * WF,float * ZF,
	       float * WZF, int & IOPT, float * VV0,
	       float * EE0, float & CH2PH, float & CH2Z);
}

extern "C" {
  void trkfit_(int & nhits, int * idet, int * itype, 
	       float * x, float * y, float * z, float * phireso, float * zreso,
	       float * ref, int & ierr, float * rfit, float * rfite, float & chi2, int & ndf);
  
}

int TrackFitting(int & nhits, float bField, int * idet, int * itype, // inputs 
		 float chi2PrefitCut, // inputs 
		 float * x, float * y, float * z, float * RPReso, float * ZReso, // inputs 
		 float * param, float * eparam, float & chi2, int & ndf,
		 float & chi2rphi, float & chi2z) { // outputs

  //
  // Subroutine performs track fitting taking into account energy loss and MS 
  // First simple helix fit is performed to define initial track parameters
  // at the PCA to primary IP. If simple helix fit converges and has qood quality 
  // defined by "distMax" proximity criterion, then DELPHI fitting routine trkfit.F is envoked
  // Inputs :
  //          nhits - number of tracker hits used in fit
  //          bField - magnetic field in Tesla
  //          idet - list of the hit detector identifiers 
  //          itype - list of the hit types ( 2 - planar disks, 
  //                                          3 - cyllindrical or laddered detector )
  //          distMax - proximity parameter defining quality of simple
  //                    helix fit. If maximal hit distance to the fitted helix is
  //                    greater than distMax fit is regarded to fail
  //          x,y,z - list of the hit Cartesian coordinates
  //          RPReso - list of the hit resolutions in R-Phi projection
  //          ZReso  - list of the hit Z resolutions
  // Outpus :
  //          param  - list of the track parameters at the PCA
  //          param[0] - Omega (signed curvuture)
  //          param[1] - tan(Lambda) (tangence of the dip angle)
  //          param[2] - Phi angle of the particle momentum @ PCA
  //          param[3] - D0 (IP in the R-Phi projection @ PCA )
  //          param[4] - Z0 (z coordinate displacement @ PCA in the R-Phi projection )
  //          eparam[15] - covariance matrix of the track parameters
  //          chi2 - fit chi2
  //          ndf - number of degrees of freedom
  //          returned integer indicate error flag
  //          error = 0 - no error, fit is successfull
  //          error = -1 - simple helix fit failed
  //  std::cout << "Tracking processor ----------->" << std::endl;
  //  std::cout << "# of hits = " << nhits << std::endl;
  //

  float * xhit = new float[nhits+1];
  float * yhit = new float[nhits+1];
  float * zhit = new float[nhits+1];
  int * iDet = new int[nhits+1];
  int * iTyp = new int[nhits+1];
  float * rphireso = new float[nhits+1];
  float * zreso = new float[nhits+1];
  float * phi = new float[nhits];
  double * xd = new double[nhits];
  double * yd = new double[nhits];
  float * r = new float[nhits];
  double * wfr = new double[nhits];
  float  * wfz = new float[nhits];
  float * ampl = new float[nhits];
  float ref[6];
  float rfit[6];
  float rfite[15];
  float rshape[6];
  
  float dzMin = 1.0e+20;
  float dzMax = -1.0e+20;
  float zBegin = zhit[0];
  float zEnd = zhit[1];
  for (int i=0;i<nhits;++i) {
    xhit[i] = 0.1*x[i];
    yhit[i] = 0.1*y[i];
    zhit[i] = 0.1*z[i];
    ampl[i] = 1.0;
    xd[i] = xhit[i];
    yd[i] = yhit[i];
    phi[i] = atan2(yhit[i],xhit[i]);
    r[i] = sqrt(xhit[i]*xhit[i]+yhit[i]*yhit[i]);
    wfz[i] = 1./(0.01*ZReso[i]*ZReso[i]);
    wfr[i] = 1./(0.01*RPReso[i]*RPReso[i]);
    rphireso[i] = 0.01*RPReso[i]*RPReso[i];
    zreso[i] = 0.01*ZReso[i]*ZReso[i];    
    iDet[i] = idet[i];
    iTyp[i] = itype[i];
    if (fabs(zhit[i])>dzMax) {
      dzMax = fabs(zhit[i]);
      zEnd = zhit[i];
    }
    if (fabs(zhit[i])<dzMin) {
      dzMin = fabs(zhit[i]);
      zBegin = zhit[i];
    }
  }
  int IOPT = 2;
  HelixClass  helix;
  tfithl_(nhits,xd,yd,r,phi,wfr,zhit,
	  wfz,IOPT,param,eparam,chi2rphi,chi2z);

  int ndfSh = 2*nhits-5;  
  float chi2Sh = chi2rphi + chi2z;

  if (chi2Sh/float(ndfSh) > chi2PrefitCut) {
    delete[] xhit;
    delete[] yhit;
    delete[] zhit;
    delete[] rphireso;
    delete[] zreso;
    delete[] phi;
    delete[] xd;
    delete[] yd;
    delete[] r;
    delete[] wfr;
    delete[] wfz;
    delete[] iDet;
    delete[] iTyp;
    delete[] ampl;
    chi2 = chi2Sh/float(ndfSh);
    return -1;
  }

  float tanlambda = param[1];
  float omega = param[0];  
  float phi0 = param[2];
  float z0 = param[4];
  float d0 = param[3];
  float xInitial =  -d0*sin(phi0);
  float yInitial = d0*cos(phi0);
  
  ref[0] = xInitial;
  ref[1] = yInitial;
  ref[2] = z0;
  ref[4] = phi0;
  ref[3] = 0.5*acos(-1.) - atan(tanlambda);
  ref[5] = -100.0*omega/(0.299*bField);
  ref[5] = ref[5]*fabs(sin(ref[3]));

  for (int i=0;i<6;++i)
    rshape[i] = ref[i];

  xhit[nhits] = xInitial;
  yhit[nhits] = yInitial;
  zhit[nhits] = z0;
  rphireso[nhits] = 1.0;
  zreso[nhits] = 1.0;
  iDet[nhits] = 2;
  iTyp[nhits] = 2;
  int nfit = nhits+1 ;

  int ierr;
  
  trkfit_(nfit,iDet,iTyp,xhit,yhit,zhit,rphireso,zreso,ref,ierr,rfit,rfite,chi2,ndf);  


  //  std::cout << "Delphi Fit : chi2 = " << chi2 
  //	    << " ; ndf = " << ndf
  //	    << " ; ierr = " << ierr << std::endl;

  //  std::cout << std::endl;
  if (chi2Sh <0 ) std::cout << "WARNING chi2 = " << chi2Sh << std::endl;
  if (chi2<0 && ierr==0) std::cout << "WARNING chi2 = " << chi2 << std::endl;
  
  if ((ierr != 0) || (chi2 < 0) || (chi2 > chi2Sh)) {
    ierr = 1;
    chi2 = chi2Sh;
    ndf  = ndfSh;
    for (int i=0;i<6;++i)
      rfit[i] = rshape[i];
  }
  float xx = rfit[0]*cos(rfit[1]/rfit[0]);
  float yy = rfit[0]*sin(rfit[1]/rfit[0]);
  float pos[3];
  pos[0] = 10*xx;
  pos[1] = 10*yy;
  pos[2] = 10*rfit[2];
  float mom[3];
  mom[0] = sin(rfit[3])*cos(rfit[4])/fabs(rfit[5]);
  mom[1] = sin(rfit[3])*sin(rfit[4])/fabs(rfit[5]);
  mom[2] = cos(rfit[3])/fabs(rfit[5]);
  float charge = -rfit[5]/fabs(rfit[5]);
  
  helix.Initialize_VP(pos,mom,charge,bField);
  
  param[0] = helix.getOmega();
  param[1] = helix.getTanLambda();
  param[2] = helix.getPhi0();
  param[3] = helix.getD0();
  param[4] = helix.getZ0();
  
  delete[] xhit;
  delete[] yhit;
  delete[] zhit;
  delete[] rphireso;
  delete[] zreso;
  delete[] phi;
  delete[] xd;
  delete[] yd;
  delete[] r;
  delete[] wfr;
  delete[] wfz;
  delete[] iDet;
  delete[] iTyp;
  delete[] ampl;

  return ierr;

}



SiliconTracking aSiliconTracking ;

SiliconTracking::SiliconTracking() : Processor("SiliconTracking") {

  std::vector<int> combinations;

  

  combinations.push_back(6);
  combinations.push_back(4);
  combinations.push_back(3);
  
  combinations.push_back(6);
  combinations.push_back(4);
  combinations.push_back(2);
  
  combinations.push_back(6);
  combinations.push_back(3);
  combinations.push_back(2);
  
  combinations.push_back(5);
  combinations.push_back(4);
  combinations.push_back(3);

  combinations.push_back(5);
  combinations.push_back(4);
  combinations.push_back(2);

  combinations.push_back(5);
  combinations.push_back(3);
  combinations.push_back(2);

  combinations.push_back(4);
  combinations.push_back(3);
  combinations.push_back(2);
  
  combinations.push_back(4);
  combinations.push_back(3);
  combinations.push_back(1);
  
  combinations.push_back(4);
  combinations.push_back(2);
  combinations.push_back(1);
  
  combinations.push_back(3);
  combinations.push_back(2);
  combinations.push_back(1);
  

  registerProcessorParameter("LayerCombinations",
			     "Combinations of Hits in Layers",
			     _Combinations,
			     combinations);

  std::vector<int> combinationsFTD;

  combinationsFTD.push_back(6);
  combinationsFTD.push_back(5);
  combinationsFTD.push_back(4);

  combinationsFTD.push_back(5);
  combinationsFTD.push_back(4);
  combinationsFTD.push_back(3);
  
  combinationsFTD.push_back(5);
  combinationsFTD.push_back(4);
  combinationsFTD.push_back(2);
  
  combinationsFTD.push_back(5);
  combinationsFTD.push_back(4);
  combinationsFTD.push_back(1);

  combinationsFTD.push_back(5);
  combinationsFTD.push_back(3);
  combinationsFTD.push_back(2);
  
  combinationsFTD.push_back(5);
  combinationsFTD.push_back(3);
  combinationsFTD.push_back(1);
  
  combinationsFTD.push_back(5);
  combinationsFTD.push_back(2);
  combinationsFTD.push_back(1);
  
  combinationsFTD.push_back(4);
  combinationsFTD.push_back(3);
  combinationsFTD.push_back(2);
  
  combinationsFTD.push_back(4);
  combinationsFTD.push_back(3);
  combinationsFTD.push_back(1);

  combinationsFTD.push_back(4);
  combinationsFTD.push_back(3);
  combinationsFTD.push_back(0);
  
  combinationsFTD.push_back(4);
  combinationsFTD.push_back(2);
  combinationsFTD.push_back(1);

  combinationsFTD.push_back(4);
  combinationsFTD.push_back(2);
  combinationsFTD.push_back(0);

  combinationsFTD.push_back(4);
  combinationsFTD.push_back(1);
  combinationsFTD.push_back(0);

  combinationsFTD.push_back(3);
  combinationsFTD.push_back(2);
  combinationsFTD.push_back(1);

  combinationsFTD.push_back(3);
  combinationsFTD.push_back(2);
  combinationsFTD.push_back(0);

  combinationsFTD.push_back(3);
  combinationsFTD.push_back(1);
  combinationsFTD.push_back(0);

  combinationsFTD.push_back(2);
  combinationsFTD.push_back(1);
  combinationsFTD.push_back(0);




  registerProcessorParameter("LayerCombinationsFTD",
			     "Combinations of Hits in FTD",
			     _CombinationsFTD,
			     combinationsFTD);


//   registerProcessorParameter("NumberOfLayersVTX",
// 			     "Number of Layers in VTX",
// 			     _nLayersVTX,
// 			     int(5));
  

//   registerProcessorParameter("NumberOfLayersSIT",
// 			     "Number of Layers in SIT",
// 			     _nLayersSIT,
// 			     int(0));
  

//   registerProcessorParameter("NumberOfFTDLayers",
// 			     "Number of FTD Layers",
// 			     _nLayersFTD,
// 			     int(7));
  
  registerProcessorParameter("NDivisionsInPhi",
			     "Number of divisions in Phi",
			     _nDivisionsInPhi,
			     int(60));
  
  registerProcessorParameter("NDivisionsInPhiFTD",
			     "Number of divisions in Phi for FTD",
			     _nPhiFTD,
			     int(40));
  
  registerProcessorParameter("NDivisionsInTheta",
			     "Number of divisions in Theta",
			     _nDivisionsInTheta,
			     int(60));
  
  registerProcessorParameter("VXDHitCollectionName",
			     "VXD Hit Collection Name",
			     _VTXHitCollection,
			     std::string("VTXTrackerHits"));


  registerProcessorParameter("FTDHitCollectionName",
			     "FTD Hit Collection Name",
			     _FTDHitCollection,
			     std::string("FTDTrackerHits"));  


  registerProcessorParameter("SITHitCollectionName",
			     "SIT Hit Collection Name",
			     _SITHitCollection,
			     std::string("SITTrackerHits"));  

   
  registerProcessorParameter("Chi2WRphiTriplet",
			     "Chi2WRphiTriplet",
			     _chi2WRPhiTriplet,
			     float(10.));
   
  registerProcessorParameter("Chi2WRphiQuartet",
			     "Chi2WRphiQuartet",
			     _chi2WRPhiQuartet,
			     float(10.));
   
  registerProcessorParameter("Chi2WRphiSeptet",
			      "Chi2WRphiSeptet",
			     _chi2WRPhiSeptet,
			     float(10.));
   
  registerProcessorParameter("Chi2WZTriplet",
			     "Chi2WZTriplet",
			     _chi2WZTriplet,
			     float(20.));
   
  registerProcessorParameter("Chi2WZQuartet",
			     "Chi2WZQuartet",
			     _chi2WZQuartet,
			     float(20.));
   
  registerProcessorParameter("Chi2WZSeptet",
			     "Chi2WZSeptet",
			     _chi2WZSeptet,
			     float(20.));
   
  registerProcessorParameter("Chi2FitCut",
			     "Chi2 Fit Cut",
			     _chi2FitCut,
			     float(25.0));

  registerProcessorParameter("Chi2PrefitCut",
			     "Chi2 Prefit Cut",
			     _chi2PrefitCut,
			     float(1.0e+4));
  
  // Resolutions used in fit
  //^^^^^^^^^^^^^^^^^^^^^^^^
  registerProcessorParameter("ResolutionRPhiVTX",
			     "ResolutionRPhiVTX",
			     _resolutionRPhiVTX,
			     float(0.004));

  registerProcessorParameter("ResolutionZVTX",
			     "ResolutionZVTX",
			     _resolutionZVTX,
			     float(0.004));

  registerProcessorParameter("ResolutionRPhiFTD",
			     "ResolutionRPhiFTD",
			     _resolutionRPhiFTD,
			     float(0.01));

  registerProcessorParameter("ResolutionZFTD",
			     "ResolutionZFTD",
			     _resolutionZFTD,
			     float(0.10));

  registerProcessorParameter("ResolutionRPhiSIT",
			     "ResolutionRPhiSIT",
			     _resolutionRPhiSIT,
			     float(0.01));

  registerProcessorParameter("ResolutionZSIT",
			     "ResolutionZSIT",
			     _resolutionZSIT,
			     float(0.10));


  // Some Steering Parameters for Patrec
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  registerProcessorParameter("TanLambdaCutForMerging",
			     "TanLambdaCutForMerging",
			     _tanlambdaCutForMerging,
			     float(0.1));
   
  registerProcessorParameter("PhiCutForMerging",
			     "PhiCutForMerging",
			     _phiCutForMerging,
			     float(0.1));

  registerProcessorParameter("AngleCutForMerging",
			     "Angle Cut For Merging",
			     _angleCutForMerging,
			     float(0.2));

     
  registerProcessorParameter("MinDistCutAttach",
			     "MinDistCutAttach",
			     _minDistCutAttach,
			     float(1.0));
   
  registerProcessorParameter("MinLayerToAttach",
			     "MinLayerToAttach",
			     _minimalLayerToAttach,
			     int(-2));
   
  registerProcessorParameter("CutOnD0",
			     "cut on D0 for tracks",
			     _cutOnD0,
			     float(100.0));
   
  registerProcessorParameter("CutOnZ0",
			     "cut on Z0 for tracks",
			     _cutOnZ0,
			     float(100.0));
       
  registerProcessorParameter("CutOnPt",
			     "cut on Pt",
			     _cutOnPt,
			     float(0.1));

  registerProcessorParameter("MinimalHits",
			     "minimal hits",
			     _minimalHits,
			     int(4));


  registerProcessorParameter("FastAttachment",
			     "Fast attachment",
			     _attachFast,
			     int(1));

  registerProcessorParameter("SimpleHelixFit",
			     "Simple Helix Fit ?",
			     _simpleHelixFit,
			     int(1));

  registerProcessorParameter("UseSIT",
			     "Use SIT",
			     _useSIT,
			     int(0));

  
  registerProcessorParameter("FinalRefit",
			     "Final Refit ?",
			     _finalRefit,
			     int(1));


  registerProcessorParameter("CreateMap",
			     "Create Track To MCP Relations",
			     _createMap,
			     int(1));

  //  _nLayers = _nLayersVTX + _nLayersSIT;

}



void SiliconTracking::init() { 

    _nRun = -1 ;
    _nEvt = 0 ;
    printParameters() ;

    const gear::VXDParameters& pVXDDet = Global::GEAR->getVXDParameters();
    const gear::VXDLayerLayout& pVXDLayerLayout = pVXDDet.getVXDLayerLayout();
    _nLayersVTX = int(pVXDLayerLayout.getNLayers());
    
    const gear::GearParameters& pFTDDet = Global::GEAR->getGearParameters("FTD");

    _nLayersFTD = int(pFTDDet.getDoubleVals("FTDZCoordinate").size());

    _zLayerFTD.resize( _nLayersFTD);

    for (int iL=0;iL<_nLayersFTD;++iL) {
      _zLayerFTD[iL] = float(pFTDDet.getDoubleVals("FTDZCoordinate")[iL]);
    }

    const gear::GearParameters& pSITDet = Global::GEAR->getGearParameters("SIT");
    
    _nLayersSIT = int(pSITDet.getDoubleVals("SITLayerRadius").size());

    if (_useSIT == 0)
      _nLayers = _nLayersVTX;
    else 
      _nLayers = _nLayersVTX + _nLayersSIT;
    

}


void SiliconTracking::processRunHeader( LCRunHeader* run) { 

    _nRun++ ;
    _nEvt = 0;

    // Intitialization of some constants and cuts
    PI = acos(double(-1.0));
    TWOPI = double(2.0)*PI;
    PIOVER2 = double(0.5)*PI;
    _dPhi = TWOPI/double(_nDivisionsInPhi);
    _dTheta = double(2.0)/double(_nDivisionsInTheta);
    _bField = Global::parameters->getFloatVal("BField");
    _dPhiFTD = TWOPI/double(_nPhiFTD);
    float cutOnR = _cutOnPt/(0.3*_bField);
    cutOnR = 1000.*cutOnR;
    _cutOnOmega = 1/cutOnR;

} 

void SiliconTracking::processEvent( LCEvent * evt ) { 

//   _nEvt++;
//   if (_nEvt < 211)
//     return;


  // Clearing all working dynamical arrays (vectors)
  _tracks5Hits.clear();
  _tracks4Hits.clear();
  _tracks3Hits.clear();
  _trackImplVec.clear();

  int success = InitialiseVTX( evt );
  int successFTD = InitialiseFTD( evt );
  //  std::cout << "Event = " << _nEvt << " Run = " << _nRun << std::endl;

  if (success == 1) {
    for (int iPhi=0; iPhi<_nDivisionsInPhi; ++iPhi) 
      for (int iTheta=0; iTheta<_nDivisionsInTheta;++iTheta)
	ProcessOneSector(iPhi,iTheta); // Process one VXD sector     
    //    std::cout << "End of Processing VXD sectors" << std::endl;
  }

  if (successFTD == 1) {
    TrackingInFTD(); // Perform tracking in the FTD
    //    std::cout << "End of Processing FTD sectors" << std::endl;
  }

  if (success == 1 || successFTD == 1) {

    Sorting( _tracks5Hits);
    Sorting( _tracks4Hits);
    Sorting( _tracks3Hits);
    //    std::cout << "End of Sorting " << std::endl;

    int nTrk = int(_tracks5Hits.size());
    for (int iTrk=0; iTrk<nTrk;++iTrk) {
      TrackExtended * trackAR = _tracks5Hits[iTrk];
      CreateTrack( trackAR );
    }
    //    std::cout << "End of creating 5 hits tracks " << std::endl;
    

    nTrk = int(_tracks4Hits.size());
    for (int iTrk=0; iTrk<nTrk;++iTrk) {
      TrackExtended * trackAR = _tracks4Hits[iTrk];
      CreateTrack( trackAR );
    }
    //    std::cout << "End of creating 4 hits tracks " << std::endl;

    nTrk = int(_tracks3Hits.size());
    for (int iTrk=0; iTrk<nTrk;++iTrk) {
      TrackExtended * trackAR = _tracks3Hits[iTrk];
      CreateTrack( trackAR );
    }
    //    std::cout << "End of creating 3 hits tracks " << std::endl;

    if (_attachFast == 0) {
      AttachRemainingVTXHitsSlow();
      AttachRemainingFTDHitsSlow();
    }
    else {
      AttachRemainingVTXHitsFast();
      AttachRemainingFTDHitsFast();
    }
    //    std::cout << "End of picking up remaining hits " << std::endl;

    if (_finalRefit > 0) {
      FinalRefit();
    }

    LCCollectionVec * trkCol = new LCCollectionVec(LCIO::TRACK);
    LCCollectionVec * relCol = NULL;
    if (_createMap)
      relCol = new LCCollectionVec(LCIO::LCRELATION);
    nTrk = int(_trackImplVec.size());
    for (int iTrk=0; iTrk<nTrk;++iTrk) {
      TrackExtended * trackAR = _trackImplVec[iTrk];
      TrackerHitExtendedVec hitVec = trackAR->getTrackerHitExtendedVec();
      int nh = int(hitVec.size());
      if (nh >= _minimalHits) {
	TrackImpl * trackImpl = new TrackImpl();
	trackImpl->setOmega(trackAR->getOmega());
	trackImpl->setTanLambda(trackAR->getTanLambda());
	trackImpl->setPhi(trackAR->getPhi());
	trackImpl->setD0(trackAR->getD0());
	trackImpl->setZ0(trackAR->getZ0());
	trackImpl->setChi2(trackAR->getChi2());
	int nHits = int(hitVec.size());
	std::vector<MCParticle*> mcPointers ;
	std::vector<int> mcHits ;
	mcPointers.clear();
	mcHits.clear();
	for (int iHit=0;iHit<nHits;++iHit) {
	  TrackerHit * trkHit = hitVec[iHit]->getTrackerHit();
	  trackImpl->addHit(trkHit);
	  if (_createMap > 0) {
	    int nSH = int(trkHit->getRawHits().size());
	    for (int ish=0;ish<nSH;++ish) {
	      SimTrackerHit * simHit = dynamic_cast<SimTrackerHit*>(trkHit->getRawHits()[ish]);
	      MCParticle * mcp = simHit->getMCParticle();
	      bool found = false;
	      int nMCP = int(mcPointers.size());
	      for (int iMCP=0;iMCP<nMCP;++iMCP) {
		if (mcp == mcPointers[iMCP]) {
		  found = true;
		  mcHits[iMCP]++;
		  break;
		}
	      }
	      if (!found) {
		mcPointers.push_back(mcp);
		mcHits.push_back(1);
	      }
	    }
	  }
	}
	trkCol->addElement(trackImpl);
	if (_createMap > 0) {
	  int nRel = int(mcPointers.size());
	  for (int k=0;k<nRel;++k) {
	    LCRelationImpl* tpclcRel = new LCRelationImpl;
	    MCParticle * mcp = mcPointers[k];
	    tpclcRel->setFrom (trackImpl);
	    tpclcRel->setTo (mcp);
	    float weight = (float)(mcHits[k])/(float)(trackImpl->getTrackerHits().size());
	    tpclcRel->setWeight(weight);
	    relCol->addElement(tpclcRel);
	  }
	}
      }
    }        
    evt->addCollection(trkCol,"SiTracks");     
    if (_createMap>0)
      evt->addCollection(relCol,"SiTracksMCP");
  }
  CleanUp();
  //  std::cout << "Event is done " << std::endl;
  _nEvt++;

}


void SiliconTracking::CleanUp() {
  int nTrk = int(_tracks5Hits.size());
  for (int iTrk=0; iTrk<nTrk;++iTrk) {
    TrackExtended * trackAR = _tracks5Hits[iTrk];
    delete trackAR;
  }
  nTrk = int(_tracks4Hits.size());
  for (int iTrk=0; iTrk<nTrk;++iTrk) {
    TrackExtended * trackAR = _tracks4Hits[iTrk];
    delete trackAR;
  }
  nTrk = int(_tracks3Hits.size());
  for (int iTrk=0; iTrk<nTrk;++iTrk) {
    TrackExtended * trackAR = _tracks3Hits[iTrk];
    delete trackAR;
  }

  for (int il=0;il<_nLayers;++il) {
    for (int ip=0;ip<_nDivisionsInPhi;++ip) {
      for (int it=0;it<_nDivisionsInTheta; ++it) {
	int iCode = il + _nLayers*ip + _nLayers*_nDivisionsInPhi*it;      
	TrackerHitExtendedVec hitVec = _sectors[iCode];
	int nH = int(hitVec.size());
	for (int iH=0; iH<nH; ++iH) {
	  TrackerHitExtended * hit = hitVec[iH];
	  delete hit;
	}
      }
    }
  }

  for (int iS=0;iS<2;++iS) {
    for (int layer=0;layer<_nLayersFTD;++layer) {
      for (int ip=0;ip<_nPhiFTD;++ip) {
	int iCode = iS + 2*layer + 2*_nLayersFTD*ip;
	TrackerHitExtendedVec hitVec = _sectorsFTD[iCode];
	int nH = int(hitVec.size());
	for (int iH=0; iH<nH; ++iH) {
	  TrackerHitExtended * hit = hitVec[iH];
	  delete hit;
	}
      }
    }
  }

}

int SiliconTracking::InitialiseFTD(LCEvent * evt) {
  int success = 1;
  _nTotalFTDHits = 0;
  _sectorsFTD.clear();
  _sectorsFTD.resize(2*_nLayersFTD*_nPhiFTD);
  for (int i=0; i<2*_nLayersFTD*_nPhiFTD;++i) {
    TrackerHitExtendedVec hitVec;
    hitVec.clear();
    _sectorsFTD.push_back(hitVec);    
  }

  // Reading out FTD Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    LCCollection * hitCollection = evt->getCollection(_FTDHitCollection.c_str());
    int nelem = hitCollection->getNumberOfElements();
    _nTotalFTDHits = nelem;
    for (int ielem=0; ielem<nelem; ++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
      hitExt->setResolutionRPhi(_resolutionRPhiFTD);
      hitExt->setResolutionZ(_resolutionZFTD);
      hitExt->setType(int(2));
      hitExt->setDet(int(2));
      double pos[3];
      for (int i=0; i<3; ++i) {
	pos[i] = hit->getPosition()[i];
      }
      double Phi = atan2(pos[1],pos[0]);
      if (Phi < 0.) Phi = Phi + TWOPI;
      int layer = abs(hit->getType()) - 1;
      if (layer < 0 || layer > _nLayersFTD-1) {
	std::cout << "SiliconTracking => fatal error in FTD : layer is outside allowed range : " << layer << std::endl;
	exit(1);
      }
      int iPhi = int(Phi/_dPhiFTD);
      int iSemiSphere = 0;
      if (hit->getType() > 0) 
	iSemiSphere = 1;
      int iCode = iSemiSphere + 2*layer + 2*_nLayersFTD*iPhi;
      _sectorsFTD[iCode].push_back( hitExt );
    }
  }
  catch(DataNotAvailableException &e ) {
    success = 0;
  }
  
  return success;
}

int SiliconTracking::InitialiseVTX(LCEvent * evt) {

  int success = 1;

  _nTotalVTXHits = 0;
  _nTotalSITHits = 0;
  _sectors.clear();
  _sectors.resize(_nLayers*_nDivisionsInPhi*_nDivisionsInTheta);

  for (int i=0; i<_nLayers*_nDivisionsInPhi*_nDivisionsInTheta; ++i) {
    TrackerHitExtendedVec hitVec;
    hitVec.clear();
    _sectors.push_back(hitVec);
  }

  // Reading out VTX Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
  try {
    LCCollection * hitCollection = evt->getCollection(_VTXHitCollection.c_str());
    int nelem = hitCollection->getNumberOfElements();
    _nTotalVTXHits = nelem;
    for (int ielem=0; ielem<nelem; ++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
      hitExt->setResolutionRPhi(_resolutionRPhiVTX);
      hitExt->setResolutionZ(_resolutionZVTX);
      hitExt->setType(int(3));
      hitExt->setDet(int(3));
      double pos[3];
      double radius = 0;
      for (int i=0; i<3; ++i) {
	pos[i] = hit->getPosition()[i];
	radius += pos[i]*pos[i];
      }
      radius = sqrt(radius);
      double cosTheta = pos[2]/radius;
      double Phi = atan2(pos[1],pos[0]);
      if (Phi < 0.) Phi = Phi + TWOPI;
      int layer = (hit->getType() - 100) - 1;
      if (layer < 0 || layer >= _nLayers) {
	std::cout << "SiliconTracking => fatal error in VTX : layer is outside allowed range : " << layer << std::endl;
	exit(1);
      }
      int iPhi = int(Phi/_dPhi);
      int iTheta = int ((cosTheta + double(1.0))/_dTheta);
      int iCode = layer + _nLayers*iPhi + _nLayers*_nDivisionsInPhi*iTheta;      
      _sectors[iCode].push_back( hitExt );
    }
  }
  catch(DataNotAvailableException &e) {
    success = 0;
  }

  if (_useSIT > 0 ) {
  // Reading out SIT Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    try {
      LCCollection *hitCollection = evt->getCollection(_SITHitCollection.c_str());
      int nelem = hitCollection->getNumberOfElements();
      _nTotalSITHits = nelem;
      for (int ielem=0; ielem<nelem; ++ielem) {
	TrackerHit * hit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));
	TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
	hitExt->setResolutionRPhi(_resolutionRPhiSIT);
	hitExt->setResolutionZ(_resolutionZSIT);
	hitExt->setType(int(3));
	hitExt->setDet(int(7));
	double pos[3];
	double radius = 0;
	for (int i=0; i<3; ++i) {
	  pos[i] = hit->getPosition()[i];
	  radius += pos[i]*pos[i];
	}
	radius = sqrt(radius);
	double cosTheta = pos[2]/radius;
	double Phi = atan2(pos[1],pos[0]);
	if (Phi < 0.) Phi = Phi + TWOPI;
	int layer = (hit->getType() - 400) - 1 + _nLayersVTX;
	//	std::cout << "# VTX = " << _nLayersVTX
	//		  << "# SIT = " << _nLayersSIT
	//		  << " Type = " << hit->getType() 
	//		  << " Layer = " << layer << std::endl;
	if (layer < 0 || layer >= _nLayers) {
	  std::cout << "SiliconTracking => fatal error in SIT : layer is outside allowed range : " << layer << std::endl;
	  exit(1);
	}
	int iPhi = int(Phi/_dPhi);
	int iTheta = int ((cosTheta + double(1.0))/_dTheta);
	int iCode = layer + _nLayers*iPhi + _nLayers*_nDivisionsInPhi*iTheta;      
	_sectors[iCode].push_back( hitExt );
      }
    }
    catch(DataNotAvailableException &e) {
      
    }
  }
  return success;

}

void SiliconTracking::check(LCEvent * evt) { };

void SiliconTracking::end() {

}


void SiliconTracking::ProcessOneSector(int iPhi, int iTheta) {
  
  int iPhi_Up    = iPhi + 1;
  int iPhi_Low   = iPhi - 1;
  int iTheta_Up  = iTheta + 1; 
  int iTheta_Low = iTheta - 1;
  if (iTheta_Low < 0) iTheta_Low = 0;
  if (iTheta_Up  >= _nDivisionsInTheta) iTheta_Up = _nDivisionsInTheta-1;

  int nComb = int( _Combinations.size() / 3 ); // number of triplet combinations
  //  std::cout << iPhi << " " << iTheta << " " << _nEvt << std::endl;
  int iNC = 0;
  for (int iComb=0; iComb < nComb; ++iComb) { // loop over triplets
    int nLR[3];
    for (int iS=0; iS<3; ++iS) {
      nLR[iS] = _Combinations[iNC];
      iNC++;
    }    
//     std::cout << iPhi << " " << iTheta << " " << nLR[0] << " " << nLR[1] << " " << nLR[2] << " " 
// 	      << std::endl;
    int iCode = nLR[0] + _nLayers*iPhi +  _nLayers*_nDivisionsInPhi*iTheta;
    TrackerHitExtendedVec hitVecOuter =  _sectors[iCode]; 
    int nHitsOuter = int(hitVecOuter.size());
    if (nHitsOuter > 0) {
      for (int ipMiddle=iPhi_Low; ipMiddle<iPhi_Up+1;ipMiddle++) { // loop over phi in the Middle
	for (int itMiddle=iTheta_Low; itMiddle<iTheta_Up+1;itMiddle++) { // loop over theta in the Middle 
	  int iPhiMiddle = ipMiddle;
	  if (ipMiddle < 0) iPhiMiddle = _nDivisionsInPhi-1;
	  if (ipMiddle >= _nDivisionsInPhi) iPhiMiddle = ipMiddle - _nDivisionsInPhi;
	  iCode = nLR[1] + _nLayers*iPhiMiddle +  _nLayers*_nDivisionsInPhi*itMiddle;
	  TrackerHitExtendedVec hitVecMiddle = _sectors[iCode];
	  int nHitsMiddle = int(hitVecMiddle.size());
	  
	  int iPhiLowInner = iPhi_Low;
	  int iPhiUpInner = iPhi_Up;
	  int iThetaLowInner = iTheta_Low;
	  int iThetaUpInner = iTheta_Up;	
	  
	  if (ipMiddle == iPhi && itMiddle==iTheta) { 
	    iPhiLowInner = iPhi_Low;
	    iPhiUpInner  = iPhi_Up;
	    iThetaLowInner = iTheta_Low;
	    iThetaUpInner = iTheta_Up;
	  }
	  else {
	    int difP = abs(ipMiddle-iPhi);
	    int difT = abs(itMiddle-iTheta);
	    int minP = min(ipMiddle,iPhi);
	    int minT = min(itMiddle,iTheta);
	    int maxP = max(ipMiddle,iPhi);
	    int maxT = max(itMiddle,iTheta);
	    if (difP==1 && difT==1) {
	      iPhiLowInner = minP;
	      iPhiUpInner = maxP;
	      iThetaLowInner = minT;
	      iThetaUpInner = maxT;
	    }
	    if (difP==0) {
	      iPhiLowInner = iPhi_Low;
	      iPhiUpInner  = iPhi_Up;
	      iThetaLowInner = minT;
	      iThetaUpInner = maxT;
	    }
	    if (difT==0) {
	      iPhiLowInner = minP;
	      iPhiUpInner  = maxP;
	      iThetaLowInner = iTheta_Low;
	      iThetaUpInner = iTheta_Up;	    
	    }
	  }		  
	  if (nHitsMiddle > 0) {
	    for (int ipInner=iPhiLowInner; ipInner<iPhiUpInner+1;ipInner++) { // loop over phi in the Inner
	      for (int itInner=iThetaLowInner; itInner<iThetaUpInner+1;itInner++) { // loop over theta in the Inner 
		int iPhiInner = ipInner;
		if (ipInner < 0) iPhiInner = _nDivisionsInPhi-1;
		if (ipInner >= _nDivisionsInPhi) iPhiInner = ipInner - _nDivisionsInPhi;
		iCode = nLR[2] + _nLayers*iPhiInner +  _nLayers*_nDivisionsInPhi*itInner;
		TrackerHitExtendedVec hitVecInner = _sectors[iCode];
		int nHitsInner = int(hitVecInner.size());
 		if (nHitsInner > 0) {
// 		  std::cout << iPhi << " " << ipMiddle << " " << ipInner << "     " 
// 			    << iTheta << " " << itMiddle << " " << itInner << "     " 
// 			    << nLR[0] << " " << nLR[1] << " " << nLR[2] << "     " 
// 			    << nHitsOuter << " " << nHitsMiddle << " " << nHitsInner << "     " 
// 			    << nHitsOuter*nHitsMiddle* nHitsInner << std::endl;
		  for (int iOuter=0; iOuter<nHitsOuter; ++iOuter) { // loop over hits in the outer sector
		    TrackerHitExtended * outerHit = hitVecOuter[iOuter];
		    for (int iMiddle=0;iMiddle<nHitsMiddle;iMiddle++) { // loop over hits in the middle sector
		      TrackerHitExtended * middleHit = hitVecMiddle[iMiddle];
		      for (int iInner=0;iInner<nHitsInner;iInner++) { // loop over hits in the inner sector
			TrackerHitExtended * innerHit = hitVecInner[iInner];
			HelixClass helix;
			TrackExtended * trackAR = TestTriplet(outerHit,middleHit,innerHit,helix);
			if ( trackAR != NULL ) {
			  int nH = BuildTrack(outerHit,middleHit,innerHit,helix,nLR[2],
					      iPhiLowInner,iPhiUpInner,
					      iThetaLowInner,iThetaUpInner,trackAR);
			  if (nH == 3) 
			    _tracks3Hits.push_back(trackAR);
			  if (nH == 4)
			    _tracks4Hits.push_back(trackAR);
			  if (nH >= 5)
			    _tracks5Hits.push_back(trackAR);		
			}	
		      } // endloop over hits in the inner sector
		    } // endloop over hits in the middle sector
		  } // endloop over hits in the outer sector
		} // endif nHitsInner > 0
	      } // endloop over theta in the Inner
	    } // endloop over phi in the Inner	    
	  } // endif nHitsMiddle > 0
	} // endloop over theta in the Middle
      } // endloop over phi in the Middle
    } // endif nHitsOuter > 0
  } // endloop over triplets
}

TrackExtended * SiliconTracking::TestTriplet(TrackerHitExtended * outerHit, 
					    TrackerHitExtended * middleHit,
					    TrackerHitExtended * innerHit,
					    HelixClass & helix) {
  /*
    Methods checks if the triplet of hits satisfies helix hypothesis
   */

  TrackExtended * trackAR = NULL;

  TrackExtendedVec trackOuterVec  = outerHit->getTrackExtendedVec();
  TrackExtendedVec trackMiddleVec = middleHit->getTrackExtendedVec();
  TrackExtendedVec trackInnerVec  = innerHit->getTrackExtendedVec();
  int nTrackOuter  = int (trackOuterVec.size());
  int nTrackMiddle = int (trackMiddleVec.size());
  int nTrackInner  = int (trackInnerVec.size());
  
//   if (nTrackInner > 0 && nTrackMiddle > 0) {
//     for (int iInner=0; iInner<nTrackInner; iInner++) {
//       for (int iMiddle=0; iMiddle<nTrackMiddle; iMiddle++) {
// 	if (trackInnerVec[iInner] == trackMiddleVec[iMiddle]) 
// 	  return trackAR;     
//       }      
//     }    
//   }

//   if (nTrackOuter > 0 && nTrackMiddle > 0) {
//     for (int iOuter=0; iOuter<nTrackOuter; iOuter++) {
//       for (int iMiddle=0; iMiddle<nTrackMiddle; iMiddle++) {
// 	if (trackOuterVec[iOuter] == trackMiddleVec[iMiddle]) 
// 	  return trackAR;     
//       }      
//     }    
//   }

//   if (nTrackOuter > 0 && nTrackInner > 0) {
//     for (int iOuter=0; iOuter<nTrackOuter; iOuter++) {
//       for (int iInner=0; iInner<nTrackInner; iInner++) {
// 	if (trackOuterVec[iOuter] == trackInnerVec[iInner]) 
// 	  return trackAR;     
//       }      
//     }    
//   }

  if (nTrackOuter > 0 && nTrackInner > 0 && nTrackMiddle > 0) {
    for (int iMiddle=0; iMiddle<nTrackMiddle ; iMiddle++) {
      for (int iOuter=0; iOuter<nTrackOuter; iOuter++) {
	for (int iInner=0; iInner<nTrackInner; iInner++) {
	  if (trackOuterVec[iOuter] == trackInnerVec[iInner] && trackInnerVec[iInner] == trackMiddleVec[iMiddle]) 
	    return trackAR;     
	}      
      }    
    }
  }

  float dZ = FastTripletCheck(innerHit, middleHit, outerHit);

  if (fabs(dZ) > _minDistCutAttach)
    return trackAR;
    


  double * xh = new double[3];
  double * yh = new double[3];
  float * zh = new float[3];
  double * wrh = new double[3];
  float * wzh = new float[3];
  float * rh = new float[3];
  float * ph = new float[3];
  float * xh_fl = new float[3];
  float * yh_fl = new float[3];
  float * rphi_reso = new float[3];
  float * z_reso = new float[3];
  int * idet_h = new int[3];
  int * ityp_h = new int[3];
  float par[5];
  float epar[15];

  xh[0] = outerHit->getTrackerHit()->getPosition()[0];
  yh[0] = outerHit->getTrackerHit()->getPosition()[1];
  zh[0] = float(outerHit->getTrackerHit()->getPosition()[2]);
  wrh[0] = double(1.0/(outerHit->getResolutionRPhi()*outerHit->getResolutionRPhi()));
  wzh[0] = 1.0/(outerHit->getResolutionZ()*outerHit->getResolutionZ());
  rphi_reso[0] = outerHit->getResolutionRPhi();
  z_reso[0] = outerHit->getResolutionZ();
  idet_h[0] = outerHit->getDet();
  ityp_h[0] = outerHit->getType();
  
  xh[1] = middleHit->getTrackerHit()->getPosition()[0];
  yh[1] = middleHit->getTrackerHit()->getPosition()[1];
  zh[1] = float(middleHit->getTrackerHit()->getPosition()[2]);
  wrh[1] = double(1.0/(middleHit->getResolutionRPhi()*middleHit->getResolutionRPhi()));
  wzh[1] = 1.0/(middleHit->getResolutionZ()*middleHit->getResolutionZ());
  rphi_reso[1] = middleHit->getResolutionRPhi();
  z_reso[1] = middleHit->getResolutionZ();
  idet_h[1] = middleHit->getDet();
  ityp_h[1] = middleHit->getType();
  

  xh[2] = innerHit->getTrackerHit()->getPosition()[0];
  yh[2] = innerHit->getTrackerHit()->getPosition()[1];
  zh[2] = float(innerHit->getTrackerHit()->getPosition()[2]);
  wrh[2] = double(1.0/(innerHit->getResolutionRPhi()*innerHit->getResolutionRPhi()));
  wzh[2] = 1.0/(innerHit->getResolutionZ()*innerHit->getResolutionZ());
  rphi_reso[2] = innerHit->getResolutionRPhi();
  z_reso[2] = innerHit->getResolutionZ();
  idet_h[2] = innerHit->getDet();
  ityp_h[2] = innerHit->getType();
  
  for (int ih=0; ih<3; ih++) {
      xh_fl[ih] = float(xh[ih]);
      yh_fl[ih] = float(yh[ih]);
      rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
      ph[ih] = atan2(yh[ih],xh[ih]);
      if (ph[ih] < 0.) 
	  ph[ih] = 2.0*acos(-1.0) + ph[ih]; 
  }

  int NPT = 3;
  int IOPT = 2;
  float chi2RPhi;
  float chi2Z;
  float chi2_D;
  int ndf_D;
  int ierr;
  if (_simpleHelixFit > 0) {
    tfithl_(NPT, xh, yh, rh, ph, wrh, zh,
	    wzh, IOPT, par, epar, chi2RPhi, chi2Z);
  }
  else {
    ierr = TrackFitting(NPT,_bField,idet_h,ityp_h,_chi2PrefitCut,
			xh_fl,yh_fl,zh,rphi_reso,z_reso,
			par,epar,chi2_D,ndf_D, chi2RPhi, chi2Z);
  }

  delete[] xh;
  delete[] yh;
  delete[] zh;
  delete[] wrh;
  delete[] wzh;
  delete[] rh;
  delete[] ph;
  delete[] xh_fl;
  delete[] yh_fl;
  delete[] rphi_reso;
  delete[] z_reso;
  delete[] idet_h;
  delete[] ityp_h;
  float omega = par[0];
  float tanlambda = par[1];
  float phi0 = par[2];
  float d0 = par[3];
  float z0 = par[4];

  float Chi2 = chi2_D/float(ndf_D);
  if (_simpleHelixFit > 0) 
    Chi2 = chi2RPhi/_chi2WRPhiTriplet+chi2Z/_chi2WZTriplet;
 
  // Check if track satisfies all conditions

  if ( Chi2 > _chi2FitCut || fabs(d0) > _cutOnD0 || fabs(z0) > _cutOnZ0 )
    return trackAR;

  helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);

  trackAR = new TrackExtended();
  trackAR->addTrackerHitExtended(outerHit);
  trackAR->addTrackerHitExtended(middleHit);
  trackAR->addTrackerHitExtended(innerHit);
  outerHit->addTrackExtended(trackAR);
  middleHit->addTrackExtended(trackAR);
  innerHit->addTrackExtended(trackAR);    
  trackAR->setD0(d0);
  trackAR->setZ0(z0);
  trackAR->setPhi(phi0);
  trackAR->setTanLambda(tanlambda);
  trackAR->setOmega(omega);
  trackAR->setChi2( Chi2 );
  trackAR->setCovMatrix(epar);
  return trackAR;

}

int SiliconTracking::BuildTrack(TrackerHitExtended * outerHit, 
			       TrackerHitExtended * middleHit,
			       TrackerHitExtended * innerHit,
			       HelixClass & helix,
			       int innerLayer,
			       int iPhiLow, int iPhiUp,
			       int iThetaLow, int iThetaUp, 
			       TrackExtended * trackAR) {
  /**
     Method for building up track in the VXD. Method starts from the found triplet and performs
     sequential attachment of hits in other layers 
   */


  for (int layer = innerLayer-1; layer>=0; layer--) { // loop over remaining layers
    float distMin = 1.0e+20;
    TrackerHitExtended * assignedhit = NULL;
    for (int ipInner=iPhiLow; ipInner<iPhiUp+1;ipInner++) { // loop over phi in the Inner region
      for (int itInner=iThetaLow; itInner<iThetaUp+1;itInner++) { // loop over theta in the Inner region 
	int iPhiInner = ipInner;
	if (ipInner < 0) iPhiInner = _nDivisionsInPhi-1;
	if (ipInner >= _nDivisionsInPhi) iPhiInner = ipInner - _nDivisionsInPhi;
	int iCode = layer + _nLayers*iPhiInner +  _nLayers*_nDivisionsInPhi*itInner;
	TrackerHitExtendedVec hitVecInner = _sectors[iCode];
	int nHitsInner = int(hitVecInner.size());
	for (int iInner=0;iInner<nHitsInner;iInner++) { // loop over hits in the Inner sector
	  TrackerHitExtended * currentHit = hitVecInner[iInner];
	  float pos[3];
	  float distance[3];
	  //	  float ref[3];
	  //	  float point[6];
	  for (int i=0; i<3; ++i) {
	    pos[i] = float(currentHit->getTrackerHit()->getPosition()[i]);
	    //	    ref[i] = helix.getReferencePoint()[i];
	  }
	  //	  float radius = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
	  float time = helix.getDistanceToPoint(pos,distance);	  
	  //	  float time = helix.getPointOnCircle(radius,ref,point);
	  if (time < 1.0e+10) {
	    // 	    float distRPhi1 = sqrt((point[0]-pos[0])*(point[0]-pos[0])+
	    // 				   (point[1]-pos[1])*(point[1]-pos[1]));
	    // 	    float distRPhi2 = sqrt((point[3]-pos[0])*(point[3]-pos[0])+
	    // 				   (point[4]-pos[1])*(point[4]-pos[1]));
	    // 	    float distZ1 = fabs(point[2]-pos[2]);
	    // 	    float distZ2 = fabs(point[5]-pos[2]);	    
	    // 	    float dist1 = distRPhi1/_resolutionRPhi + distZ1/_resolutionZ;
	    // 	    float dist2 = distRPhi2/_resolutionRPhi + distZ2/_resolutionZ;
	    // 	    float dist = fmin(dist1,dist2);
	    // 	    _distRPhi = fmin(distRPhi1,distRPhi2);
	    // 	    _distZ = fmin(distZ1,distZ2);
	    if (distance[2] < distMin) {
	      distMin = distance[2];
	      assignedhit = currentHit;
	    }
	  }
	} // endloop over hits in the Inner sector
      } // endloop over theta in the Inner region 
    } // endloop over phi in the Inner region
    if (distMin < _minDistCutAttach) {
      TrackerHitExtendedVec hvec = trackAR->getTrackerHitExtendedVec();
      int  nHits = int(hvec.size());
      double * xh = new double[nHits+1];
      double * yh = new double[nHits+1];
      float * zh = new float[nHits+1];
      double * wrh = new double[nHits+1];
      float * wzh = new float[nHits+1];
      float * rh = new float[nHits+1];
      float * ph = new float[nHits+1];
      float * xh_fl = new float[nHits+1];
      float * yh_fl = new float[nHits+1];
      float * rphi_reso = new float[nHits+1];
      float * z_reso = new float[nHits+1];
      int * ityp_h = new int[nHits+1];
      int * idet_h = new int[nHits+1];
      float par[5];
      float epar[15];
      for (int ih=0;ih<nHits;++ih) {
	TrackerHit * trkHit = hvec[ih]->getTrackerHit();
	xh[ih] = trkHit->getPosition()[0];
	yh[ih] = trkHit->getPosition()[1];
	zh[ih] = float(trkHit->getPosition()[2]);
	wrh[ih] = double(1.0/(hvec[ih]->getResolutionRPhi()*hvec[ih]->getResolutionRPhi()));
	wzh[ih] = 1.0/(hvec[ih]->getResolutionZ()*hvec[ih]->getResolutionZ());
	rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
	ph[ih] = float(atan2(yh[ih],xh[ih]));
	if (ph[ih] < 0.) 
	  ph[ih] = 2.0*acos(-1.0) + ph[ih]; 
	xh_fl[ih] = float(xh[ih]);
	yh_fl[ih] = float(yh[ih]);
	rphi_reso[ih] = hvec[ih]->getResolutionRPhi();
	z_reso[ih] = hvec[ih]->getResolutionZ();
	ityp_h[ih] = hvec[ih]->getType();
	idet_h[ih] = hvec[ih]->getDet();
      }      
      TrackerHit * assignedTrkHit = assignedhit->getTrackerHit();
      xh[nHits] = assignedTrkHit->getPosition()[0];
      yh[nHits] = assignedTrkHit->getPosition()[1];
      xh_fl[nHits] = float(xh[nHits]);
      yh_fl[nHits] = float(yh[nHits]);
      zh[nHits] = float(assignedTrkHit->getPosition()[2]);
      rh[nHits] = float(sqrt(xh[nHits]*xh[nHits]+yh[nHits]*yh[nHits]));
      ph[nHits] = float(atan2(yh[nHits],xh[nHits]));
      if (ph[nHits] < 0.) 
	ph[nHits] = 2.0*acos(-1.0) + ph[nHits]; 
      wrh[nHits] = double(1.0/(assignedhit->getResolutionRPhi()*assignedhit->getResolutionRPhi()));
      wzh[nHits] = 1.0/(assignedhit->getResolutionZ()*assignedhit->getResolutionZ());
      rphi_reso[nHits] = assignedhit->getResolutionRPhi();
      z_reso[nHits] = assignedhit->getResolutionZ();
      ityp_h[nHits] = assignedhit->getType();
      idet_h[nHits] = assignedhit->getDet();
      int NPT = nHits + 1;
      int IOPT = 2;
      float chi2RPhi;
      float chi2Z;
      float chi2_D;
      int ndf_D;
      int ierr;
      if (_simpleHelixFit > 0) {
	tfithl_(NPT, xh, yh, rh, ph, wrh, zh,
		wzh, IOPT, par, epar, chi2RPhi, chi2Z);
      }
      else {
	ierr = TrackFitting(NPT,_bField,idet_h,ityp_h,_chi2PrefitCut, 
			    xh_fl,yh_fl,zh,rphi_reso,z_reso,
			    par,epar,chi2_D,ndf_D,chi2RPhi,chi2Z);
      }

      delete[] xh;
      delete[] yh;
      delete[] zh;
      delete[] wrh;
      delete[] wzh;
      delete[] rh;
      delete[] ph;
      delete[] xh_fl;
      delete[] yh_fl;
      delete[] rphi_reso;
      delete[] z_reso;
      delete[] ityp_h;
      delete[] idet_h;

      bool validCombination = 0;
      float Chi2;
      if (_simpleHelixFit>0) {
	if ((nHits+1) == 4) {
	  Chi2 = chi2RPhi/_chi2WRPhiQuartet+chi2Z/_chi2WZQuartet;
	}	  
	if ((nHits+1) >= 5) {
	  Chi2 = chi2RPhi/_chi2WRPhiSeptet+chi2Z/_chi2WZSeptet;
	}
      }
      else {
	Chi2 = chi2_D/float(ndf_D);
      }

      validCombination = Chi2 < _chi2FitCut;

      if ( validCombination ) {
	trackAR->addTrackerHitExtended(assignedhit);
	assignedhit->addTrackExtended(trackAR);
	float omega = par[0];
	float tanlambda = par[1];
	float phi0 = par[2];
	float d0 = par[3];
	float z0 = par[4];
	helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
	trackAR->setD0(d0);
	trackAR->setZ0(z0);
	trackAR->setPhi(phi0);
	trackAR->setTanLambda(tanlambda);
	trackAR->setOmega(omega);
	trackAR->setChi2( Chi2 );
	trackAR->setCovMatrix(epar);
      }
     
    }
  } // endloop over remaining layers

  TrackerHitExtendedVec hvec = trackAR->getTrackerHitExtendedVec();  
  int nTotalHits = int(hvec.size());
  return nTotalHits;

}


void SiliconTracking::Sorting(TrackExtendedVec & trackVec) {

  int sizeOfVector = int(trackVec.size());
  TrackExtended *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++)
      {
	one = trackVec[j];
	two = trackVec[j+1];
	if( one->getChi2() > two->getChi2() )
	  {
	    Temp = trackVec[j];
	    trackVec[j] = trackVec[j+1];
	    trackVec[j+1] = Temp;
	  }
      }  
  for (int i=0; i<sizeOfVector; ++i) {
    TrackerHitExtendedVec hitVec = trackVec[i]->getTrackerHitExtendedVec();
    int nHits = int(hitVec.size());
    for (int ih=0;ih<nHits;ih++) {
      hitVec[ih]->clearTrackVec();
    }
  }
  
}

void SiliconTracking::CreateTrack(TrackExtended * trackAR ) {

  /**
   Method which creates Track out of TrackExtended objects. Check for possible
   track splitting (separate track segments in VXD and FTD) is done.
   */


  TrackerHitExtendedVec hitVec = trackAR->getTrackerHitExtendedVec();
  int nHits = int(hitVec.size());

  for (int i=0; i<nHits; ++i) {
    TrackExtendedVec trackVec = hitVec[i]->getTrackExtendedVec();
    if (trackVec.size() != 0) 
      return ;
  }

  // First check if the current track is piece of the splitted one
  // look for matching track segment

  int found = 0;
  int nTrk = int(_trackImplVec.size());
  for (int i=0; i<nTrk; ++i) {
    TrackExtended * trackOld = _trackImplVec[i];
    TrackerHitExtendedVec hitVecOld = trackOld->getTrackerHitExtendedVec();

    float phiNew = trackAR->getPhi();
    float phiOld = trackOld->getPhi();
    float thetaNew = 0.5*acos(-1.) - atan(trackAR->getTanLambda());
    float thetaOld = 0.5*acos(-1.) - atan(trackOld->getTanLambda());

    float angle = (cos(phiNew)*cos(phiOld)+sin(phiNew)*sin(phiOld))*sin(thetaNew)*sin(thetaOld)+cos(thetaNew)*cos(thetaOld);
    angle = acos(angle);

    //    if ((fabs(trackAR->getTanLambda()-trackOld->getTanLambda())<_tanlambdaCutForMerging && 
    //	fabs(trackAR->getPhi()-trackOld->getPhi())<_phiCutForMerging)) {
    if (angle < _angleCutForMerging) {
      int nHitsOld = int(hitVecOld.size());
      int nTotHits = nHits + nHitsOld;
      double * xh = new double[nTotHits];
      double * yh = new double[nTotHits];
      float * xh_fl = new float[nTotHits];
      float * yh_fl = new float[nTotHits];
      float * zh = new float[nTotHits];
      double * wrh = new double[nTotHits];
      float * wzh = new float[nTotHits];
      float * rh = new float[nTotHits];
      float * ph = new float[nTotHits];
      float * rphi_reso = new float[nTotHits];
      float * z_reso = new float[nTotHits];
      int * idet_h = new int[nTotHits];
      int * ityp_h = new int[nTotHits];
      float par[5];
      float epar[15];
      for (int ih=0;ih<nHits;++ih) {
	TrackerHit * trkHit = hitVec[ih]->getTrackerHit();
	float rR = hitVec[ih]->getResolutionRPhi();
	float rZ = hitVec[ih]->getResolutionZ();
	if (int(hitVec[ih]->getTrackExtendedVec().size()) != 0)
	  std::cout << "WARNING : HIT POINTS TO TRACK " << std::endl;
	xh[ih] = trkHit->getPosition()[0];
	yh[ih] = trkHit->getPosition()[1];
	xh_fl[ih] = float(xh[ih]);
	yh_fl[ih] = float(yh[ih]);
	zh[ih] = float(trkHit->getPosition()[2]);
	wrh[ih] = double(1.0/(rR*rR));
	wzh[ih] = 1.0/(rZ*rZ);
	rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
	ph[ih] = float(atan2(yh[ih],xh[ih]));
	rphi_reso[ih] = rR;
	z_reso[ih] = rZ;
	ityp_h[ih] =  hitVec[ih]->getType();
	idet_h[ih] =  hitVec[ih]->getDet();
      }      
      for (int ih=0;ih<nHitsOld;++ih) {
	TrackerHit * trkHit = hitVecOld[ih]->getTrackerHit();
	xh[ih+nHits] = trkHit->getPosition()[0];
	yh[ih+nHits] = trkHit->getPosition()[1];
	xh_fl[ih+nHits] = float(xh[ih+nHits]);
	yh_fl[ih+nHits] = float(yh[ih+nHits]);
	zh[ih+nHits] = float(trkHit->getPosition()[2]);
	float rR = hitVecOld[ih]->getResolutionRPhi();
	float rZ = hitVecOld[ih]->getResolutionZ();	
	wrh[ih+nHits] = double(1.0/(rR*rR));
	wzh[ih+nHits] = 1.0/(rZ*rZ);
	rh[ih+nHits] = float(sqrt(xh[ih+nHits]*xh[ih+nHits]+yh[ih+nHits]*yh[ih+nHits]));
	ph[ih+nHits] = float(atan2(yh[ih+nHits],xh[ih+nHits]));
	rphi_reso[ih+nHits] = rR;
	z_reso[ih+nHits] = rZ;
	ityp_h[ih+nHits] = hitVecOld[ih]->getType();
	idet_h[ih+nHits] = hitVecOld[ih]->getDet();
      }
      int NPT = nTotHits;
      int IOPT = 3;
      float chi2RPhi;
      float chi2Z;
      float chi2_D;
      int ndf_D;
      if (_simpleHelixFit > 0) {
	tfithl_(NPT, xh, yh, rh, ph, wrh, zh,
		wzh, IOPT, par, epar, chi2RPhi, chi2Z);
      }
      else {
	int ierr = TrackFitting(NPT,_bField,idet_h,ityp_h,_chi2PrefitCut,
				xh_fl,yh_fl,zh,rphi_reso,z_reso,
				par,epar,chi2_D,ndf_D,chi2RPhi,chi2Z);
      }


      float omega = par[0];
      float tanlambda = par[1];
      float phi0 = par[2];
      float d0 = par[3];
      float z0 = par[4];
     
      float chi2Min = chi2_D/float(ndf_D);
      if (_simpleHelixFit>0)
	  chi2Min = chi2RPhi/_chi2WRPhiSeptet+chi2Z/_chi2WZSeptet;
      float chi2Full = chi2Min;
      float chi2MinRPhi = chi2RPhi;
      float chi2MinZ = chi2Z;
      float chi2DMin = chi2_D;
      int ndfDMin  = ndf_D;
      
      int iBad = -1;
      if (chi2Full < _chi2FitCut) {
	found = 1;
      }
      else {
	float * wzhOld = new float[nTotHits];
	double * wrhOld = new double[nTotHits];
	float * rphi_resoOld = new float[nTotHits];
	float * z_resoOld = new float[nTotHits];
	for (int i=0;i<nTotHits;++i) {
	  wzhOld[i] = wzh[i];
	  wrhOld[i] = wrh[i];
	  rphi_resoOld[i] = rphi_reso[i];
	  z_resoOld[i] = z_reso[i];
	}
	for (int i=0; i<nTotHits; ++i) {
	  for (int j=0;j<nTotHits;++j) {
	    if (i == j) {
	      wrh[j] = 0.0;
	      wzh[j] = 0.0;
	      z_reso[j] = 1.0e+20;
	      rphi_reso[j] = 1.0e+20;
	    } 
	    else {
	      wrh[j] = wrhOld[j];
	      wzh[j] = wzhOld[j];
	      z_reso[j] = z_resoOld[j];
	      rphi_reso[j] = rphi_resoOld[j];
	    }
	  }
	  if (_simpleHelixFit>0){
	    tfithl_(NPT, xh, yh, rh, ph, wrh, zh, 
		    wzh, IOPT, par, epar, chi2RPhi, chi2Z);
	  }
	  else {
	    TrackFitting(NPT,_bField,idet_h,ityp_h, _chi2PrefitCut, 
			 xh_fl,yh_fl,zh,rphi_reso,z_reso,
			 par,epar,chi2_D,ndf_D,chi2RPhi,chi2Z);
	  }
	  float chi2Cur = chi2_D/float(ndf_D);
	  if (_simpleHelixFit>0)
	    chi2Cur = chi2RPhi/_chi2WRPhiSeptet+chi2Z/_chi2WZSeptet;
	  if (chi2Cur < chi2Min) {
	    chi2Min = chi2Cur;
	    chi2MinRPhi = chi2RPhi;
	    chi2MinZ = chi2Z;
	    chi2DMin = chi2_D;
	    ndfDMin = ndf_D;
	    omega = par[0];
	    tanlambda = par[1];
	    phi0 = par[2];
	    d0 = par[3];
	    z0 = par[4];
	    iBad = i;
	  }
	}
	if (chi2Min < _chi2FitCut) {
	  found = 1;
	}
	delete[] wzhOld;
	delete[] wrhOld;
	delete[] rphi_resoOld;
	delete[] z_resoOld;
      }

      // Splitted track is found.
      // Attach hits belonging to the current track segment to  
      // the track already created
      if (found == 1) {
	trackOld->ClearTrackerHitExtendedVec();
	for (int i=0;i<nHits;++i) {
	  TrackerHitExtended * trkHit = hitVec[i];
	  trkHit->clearTrackVec();
	  if (i == iBad) {	    
	  }
	  else {
	    trackOld->addTrackerHitExtended(trkHit);
	    trkHit->addTrackExtended( trackOld );
	  }
	}  
	for (int i=0;i<nHitsOld;++i) {
	  int icur = i+nHits;
	  TrackerHitExtended * trkHit = hitVecOld[i];
	  trkHit->clearTrackVec();
	  if (icur == iBad) {
	  }
	  else {
	    trackOld->addTrackerHitExtended(trkHit);
	    trkHit->addTrackExtended( trackOld );
	  }
	}
	trackOld->setOmega(omega);
	trackOld->setTanLambda(tanlambda);
	trackOld->setPhi(phi0);
	trackOld->setD0(d0);
	trackOld->setZ0(z0);
	trackOld->setChi2(chi2Min );	
      }

      delete[] xh;
      delete[] yh;
      delete[] zh;
      delete[] wrh;
      delete[] wzh;
      delete[] rh;
      delete[] ph;
      delete[] xh_fl;
      delete[] yh_fl;
      delete[] idet_h;
      delete[] ityp_h;
      delete[] rphi_reso;
      delete[] z_reso;

    }
    if (found == 1)
      break;
  }

  // Candidate is a unique track
  // No other segments are found
  if (found == 0 ) {
    _trackImplVec.push_back(trackAR);
    for (int i=0;i<nHits;++i) {
      TrackerHitExtended * hit = hitVec[i];
      hit->addTrackExtended( trackAR );
    }
  }
    

}

void SiliconTracking::AttachRemainingVTXHitsFast() {

  std::vector<TrackerHitExtendedVec> nonAttachedHits;
  nonAttachedHits.resize(_nDivisionsInPhi*_nDivisionsInTheta);
  std::vector<TrackExtendedVec> trackVector;
  trackVector.resize(_nDivisionsInPhi*_nDivisionsInTheta);
  int nTrk = int(_trackImplVec.size());

  for (int iTrk=0;iTrk<nTrk;++iTrk) {
    TrackExtended * track = _trackImplVec[iTrk];
    double Phi = double(track->getPhi());
    if (Phi < 0)
      Phi = Phi + TWOPI;
    float tanlambda = track->getTanLambda();
    double cosTheta = double(tanlambda/sqrt(1+tanlambda*tanlambda));
    int iPhi = int(Phi/_dPhi);
    int iTheta = int ((cosTheta + double(1.0))/_dTheta);
    int iCode = iPhi + _nDivisionsInPhi*iTheta; 
    trackVector[iCode].push_back( track );
  }
    
  for (int il=0;il<_nLayers;++il) {
    for (int ip=0;ip<_nDivisionsInPhi;++ip) {
      for (int it=0;it<_nDivisionsInTheta; ++it) {
	int iCode = il + _nLayers*ip + _nLayers*_nDivisionsInPhi*it;      
	TrackerHitExtendedVec hitVec = _sectors[iCode];
	int nH = int(hitVec.size());
	for (int iH=0; iH<nH; ++iH) {
	  TrackerHitExtended * hitExt = hitVec[iH];
	  TrackExtendedVec trackVec = hitExt->getTrackExtendedVec();
	  if (trackVec.size()==0) {
	    TrackerHit * hit = hitExt->getTrackerHit();
	    double pos[3];
	    double radius = 0;
	    for (int i=0; i<3; ++i) {
	      pos[i] = hit->getPosition()[i];
	      radius += pos[i]*pos[i];
	    }
	    radius = sqrt(radius);
	    double cosTheta = pos[2]/radius;
	    double Phi = atan2(pos[1],pos[0]);
	    if (Phi < 0.) Phi = Phi + TWOPI;
	    int layer = hit->getType() - 1;
	    if (layer < 0) {
	      std::cout << "Fatal error; layer < 0 : " << layer << std::endl;
	      exit(1);
	    }
	    int iPhi = int(Phi/_dPhi);
	    int iTheta = int ((cosTheta + double(1.0))/_dTheta);
	    int iCode = iPhi + _nDivisionsInPhi*iTheta;      
	    nonAttachedHits[iCode].push_back( hitExt );
	  }
	}
      }
    }
  }

  for (int iT=0; iT<_nDivisionsInTheta; ++iT) {
    for (int iP=0; iP<_nDivisionsInPhi; ++iP) {
      int iCode = iP + _nDivisionsInPhi*iT; 
      int nHits = int(nonAttachedHits[iCode].size());
      int iT1 = iT - 1;
      int iT2 = iT + 1; 
      if (iT == 0) {
	iT1 = iT;
	iT2 = iT1 + 1;
      }
      if (iT == _nDivisionsInTheta - 1) {
	iT2 = iT;
	iT1 = iT2 - 1;
      }
      int iPHI[3];
      iPHI[0] = iP - 1;
      iPHI[1] = iP;
      iPHI[2] = iP + 1;
      if (iP == 0) 
	iPHI[0] = _nDivisionsInPhi - 1;
      if (iP == _nDivisionsInPhi - 1 )
	iPHI[2] = 0;

      for (int ihit = 0; ihit<nHits; ++ihit) {
	
	TrackerHitExtended * hit = nonAttachedHits[iCode][ihit];
	TrackExtended * trackToAttach = NULL;
	float minDist = 1.0e+6;

	for (int iTheta = iT1; iTheta <iT2+1; ++iTheta) {
	  for (int indexP=0;indexP<3;++indexP) {
	    int iPhi = iPHI[indexP];	    
	    int iCodeForTrack = iPhi + _nDivisionsInPhi*iTheta;
	    int nTrk = int(trackVector[iCodeForTrack].size());
	    for (int iTrk=0; iTrk<nTrk; ++iTrk) {	  
	      TrackExtended * trackAR = trackVector[iCodeForTrack][iTrk];
	      float phi0 = trackAR->getPhi();
	      float d0 = trackAR->getD0();
	      float z0 = trackAR->getZ0();
	      float omega = trackAR->getOmega();
	      float tanlambda = trackAR->getTanLambda();
	      HelixClass helix;
	      helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
	      int layer = hit->getTrackerHit()->getType() - 1;
	      if (layer > _minimalLayerToAttach) {
		float pos[3];
		for (int i=0; i<3; ++i) 
		  pos[i] = hit->getTrackerHit()->getPosition()[i];      
		float distance[3];
		float time = helix.getDistanceToPoint(pos,distance);
		if (time < 1.0e+10) {
		  if (distance[2] < minDist) {
		    minDist = distance[2];
		    trackToAttach = trackAR;
		  }		  	   
		}    
	      }
	    }
	  }
	}
	if (minDist < _minDistCutAttach && trackToAttach != NULL) {
	  AttachHitToTrack(trackToAttach,hit);
	}      
      }
    }
  }
}

void SiliconTracking::AttachRemainingVTXHitsSlow() {
  TrackerHitExtendedVec nonAttachedHits;
  nonAttachedHits.clear();

  for (int il=0;il<_nLayers;++il) {
    for (int ip=0;ip<_nDivisionsInPhi;++ip) {
      for (int it=0;it<_nDivisionsInTheta; ++it) {
	int iCode = il + _nLayers*ip + _nLayers*_nDivisionsInPhi*it;      
	TrackerHitExtendedVec hitVec = _sectors[iCode];
	int nH = int(hitVec.size());
	for (int iH=0; iH<nH; ++iH) {
	  TrackerHitExtended * hit = hitVec[iH];
	  TrackExtendedVec trackVec = hit->getTrackExtendedVec();
	  if (trackVec.size()==0)
	    nonAttachedHits.push_back( hit );
	}
      }
    }
  }

  int nNotAttached = int(nonAttachedHits.size());

  int nTrk = int(_trackImplVec.size()); 
  for (int iHit=0; iHit<nNotAttached; ++iHit) {
    TrackerHitExtended * hit = nonAttachedHits[iHit];
    int layer = hit->getTrackerHit()->getType() - 1;
    if (layer > _minimalLayerToAttach) {
      float pos[3];
      for (int i=0; i<3; ++i) 
	pos[i] = hit->getTrackerHit()->getPosition()[i];      
      float minDist = 1e+10;
      TrackExtended * trackToAttach = NULL;
      for (int iTrk=0; iTrk<nTrk; ++iTrk) {
	TrackExtended * trackAR = _trackImplVec[iTrk];
	HelixClass helix;
	float phi0 = trackAR->getPhi();
	float d0 = trackAR->getD0();
	float z0 = trackAR->getZ0();
	float omega = trackAR->getOmega();
	float tanlambda = trackAR->getTanLambda();
	helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
	float distance[3];
	float time = helix.getDistanceToPoint(pos,distance);
	if (time < 1.0e+10) {
	  if (distance[2] < minDist) {
	    minDist = distance[2];
	    trackToAttach = trackAR;
	  }
	}
      }
      if (minDist < _minDistCutAttach && trackToAttach != NULL) {
	AttachHitToTrack(trackToAttach,hit);
      }      
    }
  }  
}

void SiliconTracking::AttachRemainingFTDHitsSlow() {
  TrackerHitExtendedVec nonAttachedHits;
  nonAttachedHits.clear();

  for (int iS=0;iS<2;++iS) {
    for (int layer=0;layer<_nLayersFTD;++layer) {
      for (int ip=0;ip<_nPhiFTD;++ip) {
	int iCode = iS + 2*layer + 2*_nLayersFTD*ip;      
	TrackerHitExtendedVec hitVec = _sectorsFTD[iCode];
	int nH = int(hitVec.size());
	for (int iH=0; iH<nH; ++iH) {
	  TrackerHitExtended * hit = hitVec[iH];
	  TrackExtendedVec trackVec = hit->getTrackExtendedVec();
	  if (trackVec.size()==0)
	    nonAttachedHits.push_back( hit );
	}
      }
    }
  }

  int nNotAttached = int(nonAttachedHits.size());

  int nTrk = int(_trackImplVec.size()); 
  for (int iHit=0; iHit<nNotAttached; ++iHit) {
    TrackerHitExtended * hit = nonAttachedHits[iHit];
    int layer = hit->getTrackerHit()->getType();
    float pos[3];
    for (int i=0; i<3; ++i) 
      pos[i] = hit->getTrackerHit()->getPosition()[i];      
    float minDist = 1e+10;
    TrackExtended * trackToAttach = NULL;
    for (int iTrk=0; iTrk<nTrk; ++iTrk) {
      TrackExtended * trackAR = _trackImplVec[iTrk];
      HelixClass helix;
      float phi0 = trackAR->getPhi();
      float d0 = trackAR->getD0();
      float z0 = trackAR->getZ0();
      float omega = trackAR->getOmega();
      float tanlambda = trackAR->getTanLambda();
      if (tanlambda*float(layer) > 0) {
	helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
	float distance[3];
	float time = helix.getDistanceToPoint(pos,distance);
	if (time < 1.0e+10) {
	  if (distance[2] < minDist) {
	    minDist = distance[2];
	    trackToAttach = trackAR;
	  }
	}
      }
    }
    if (minDist < _minDistCutAttach && trackToAttach != NULL) {
      AttachHitToTrack(trackToAttach,hit);
    }      
  }  
}


void SiliconTracking::AttachRemainingFTDHitsFast() {
  int nTrk = _trackImplVec.size();

  for (int iTrk=0; iTrk<nTrk; ++iTrk) {
    TrackExtended * trackAR = _trackImplVec[iTrk];
    HelixClass helix;
    float phi0 = trackAR->getPhi();
    float d0 = trackAR->getD0();
    float z0 = trackAR->getZ0();
    float omega = trackAR->getOmega();
    float tanlambda = trackAR->getTanLambda();
    helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
    int iSemiSphere = 0;
    if (tanlambda > 0) 
      iSemiSphere = 1;
    float ref[3];
    for (int i=0;i<3;++i) 
      ref[i] = helix.getReferencePoint()[i];
    // Start loop over FTD layers
    for (int layer=0;layer<_nLayersFTD;layer++) {
      float ZL = _zLayerFTD[layer];
      if (iSemiSphere == 0)
	ZL = - ZL;
      float point[3];
      helix.getPointInZ(ZL,ref,point);
      float Phi = atan2(point[1],point[0]);
      if (Phi < 0) 
	Phi = Phi + TWOPI;
      int iPhi = int(Phi/_dPhiFTD);
      float distMin = 1e+10;
      TrackerHitExtended * attachedHit = NULL;     
      for (int iP=iPhi-1;iP<=iPhi+1;++iP) {
	int iPP = iP;
	if (iP < 0) 
	  iPP = iP + _nPhiFTD;
	if (iP >= _nPhiFTD)
	  iPP = iP - _nPhiFTD;	
	int iCode = iSemiSphere + 2*layer + 2*_nLayersFTD*iPP;
	int nHits = int(_sectorsFTD[iCode].size());
	for (int iHit=0;iHit<nHits;++iHit) {
	  TrackerHitExtended * hit = _sectorsFTD[iCode][iHit];
	  float pos[3];
	  for (int i=0;i<3;++i) {
	    pos[i] = hit->getTrackerHit()->getPosition()[i];
	  }
	  float distance[3];
	  float time = helix.getDistanceToPoint(pos,distance);
	  if (time < 1.0e+10) {
	    if (distance[2] < distMin) {
	      distMin = distance[2];
	      attachedHit = hit;
	    }
	  }
	}
      }
      if (distMin < _minDistCutAttach && attachedHit != NULL) {
	AttachHitToTrack(trackAR,attachedHit);
      }
    }
  }
}

void SiliconTracking::TrackingInFTD() {
  int nComb = int(_CombinationsFTD.size()) / 3;
  for (int iComb=0;iComb<nComb;++iComb) {
    int nLS[3];
    nLS[0] = _CombinationsFTD[3*iComb];
    nLS[1] = _CombinationsFTD[3*iComb+1];
    nLS[2] = _CombinationsFTD[3*iComb+2];
    for (int iS=0;iS<2;++iS) {
      //      std::cout << "Combinations : " << iS << " " << nLS[0] << " " << nLS[1] << " " << nLS[2] << std::endl;
      //      int iC = iS + 2*nLS[0];
      //      TrackerHitExtendedVec hitVec = _sectorsFTD[iC];
      //      int nO = int(hitVec.size());
      //      iC = iS + 2*nLS[1];
      //      hitVec = _sectorsFTD[iC];
      //      int nM = int(hitVec.size());
      //      iC = iS + 2*nLS[2];
      //      hitVec = _sectorsFTD[iC];
      //      int nI = int(hitVec.size());
      //      std::cout << nO << " " << nM << " " << nI << std::endl;
      for (int ipOuter=0;ipOuter<_nPhiFTD;++ipOuter) {
	int ipMiddleLow = ipOuter - _nPhiFTD / 4 - 1;
	int ipMiddleUp  = ipOuter + _nPhiFTD / 4 + 1;
	int iCodeOuter = iS + 2*nLS[0] + 2*_nLayersFTD*ipOuter;
	TrackerHitExtendedVec hitVecOuter = _sectorsFTD[iCodeOuter];
	int nOuter = int(hitVecOuter.size());
	for (int iOuter=0;iOuter<nOuter;++iOuter) {
	  TrackerHitExtended * hitOuter = hitVecOuter[iOuter];
	  //	  for (int ipMiddle=ipMiddleLow;ipMiddle<=ipMiddleUp;++ipMiddle) {
	  for(int ipMiddle=0;ipMiddle<_nPhiFTD;++ipMiddle) {
	    int ipM = ipMiddle;
	    if (ipM < 0) 
	      ipM = ipMiddle + _nPhiFTD;
	    if (ipM >= _nPhiFTD) 
	      ipM = ipMiddle - _nPhiFTD;
	    int iCodeMiddle = iS + 2*nLS[1] + 2*_nLayersFTD*ipM;
	    TrackerHitExtendedVec hitVecMiddle = _sectorsFTD[iCodeMiddle];
	    int ipInnerLow,ipInnerUp;	    
	    if (ipMiddle < ipOuter) {
	      ipInnerLow = ipMiddleLow;
	      ipInnerUp =  ipOuter + 1;
	    }
	    else {
	      ipInnerLow = ipOuter - 1;
	      ipInnerUp  = ipMiddleUp;
	    }
	    int nMiddle = int(hitVecMiddle.size());
	    for (int iMiddle=0;iMiddle<nMiddle;++iMiddle) {
	      TrackerHitExtended * hitMiddle = hitVecMiddle[iMiddle];
	      //	      for (int ipInner=ipInnerLow;ipInner<=ipInnerUp;++ipInner) {
	      for (int ipInner=0;ipInner<_nPhiFTD;++ipInner) {
		int ipI = ipInner;
		if (ipI < 0)
		  ipI = ipInner + _nPhiFTD;
		if (ipI >= _nPhiFTD) 
		  ipI = ipInner - _nPhiFTD;
		int iCodeInner = iS + 2*nLS[2] + 2*_nLayersFTD*ipI;
		TrackerHitExtendedVec hitVecInner = _sectorsFTD[iCodeInner];
		int nInner = int(hitVecInner.size());
		for (int iInner=0;iInner<nInner;++iInner) {
		  TrackerHitExtended * hitInner = hitVecInner[iInner];
		  HelixClass helix;
		  // 	 std::cout << hitOuter->getTrackerHit()->getType() << " " 
		  // 		   << hitMiddle->getTrackerHit()->getType() << " " 
		  // 		   << hitInner->getTrackerHit()->getType() << std::endl;
		  TrackExtended * trackAR = TestTriplet(hitOuter,hitMiddle,hitInner,helix);
		  if (trackAR != NULL) {
		    //	  std::cout << "FTD triplet found" << std::endl;
		    int nH = BuildTrackFTD(trackAR,nLS,iS);
		    if (nH == 3) 
		      _tracks3Hits.push_back(trackAR);
		    if (nH == 4)
		      _tracks4Hits.push_back(trackAR);
		    if (nH >= 5)
		      _tracks5Hits.push_back(trackAR);
		  }
		}
	      }
	    }
	  }	  
	}
      }
    }
  }
}


int SiliconTracking::BuildTrackFTD(TrackExtended * trackAR, int * nLR, int iS) {
  //  std::cout << "Layers = " << nLR[0] << " " << nLR[1] << " " << nLR[2] << std::endl;
  for (int iL=0;iL<_nLayersFTD;++iL) {
    if (iL != nLR[0] && iL != nLR[1] && iL != nLR[2]) {
      HelixClass helix;
      float d0 = trackAR->getD0();
      float z0 = trackAR->getZ0();
      float phi0 = trackAR->getPhi();
      float tanlambda = trackAR->getTanLambda();
      float omega = trackAR->getOmega();
      helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
      float ref[3];
      for (int i=0;i<3;++i) {
	ref[i] = helix.getReferencePoint()[i];
      }
      float point[3];
      float ZL = _zLayerFTD[iL];
      if (iS == 0) 
	ZL = - ZL;
      helix.getPointInZ(ZL,ref,point);
      float Phi = atan2(point[1],point[0]);
      int iPhi = int(Phi/_dPhiFTD);
      float distMin = 1e+6;
      TrackerHitExtended * attachedHit = NULL;
      for (int ip=0;ip<=_nPhiFTD;++ip) {
	int iP = ip;
	if (iP < 0)
	  iP = ip + _nPhiFTD;
	if (iP >= _nPhiFTD)
	  iP = ip - _nPhiFTD;	
	int iCode = iS + 2*iL + 2*_nLayersFTD*iP;
	TrackerHitExtendedVec hitVec = _sectorsFTD[iCode];
	int nH = int(hitVec.size());
	for (int iH=0; iH<nH; ++iH) {
	  TrackerHitExtended * hit = hitVec[iH];
	  TrackerHit * trkHit = hit->getTrackerHit();
	  float pos[3];
	  for (int i=0;i<3;++i)
	    pos[i] = float(trkHit->getPosition()[i]);
	  float distance[3];
	  float time = helix.getDistanceToPoint(pos,distance);
	  if (time < 1.0e+10) {
	    if (distance[2] < distMin) {
	      distMin = distance[2];
	      attachedHit = hit;
	    }
	  }
	}
      }
      //      std::cout << "Layer = " << iL << "  distMin = " << distMin << std::endl;
      if (distMin < _minDistCutAttach && attachedHit != NULL) {
	AttachHitToTrack( trackAR, attachedHit );
      }
    }
  }
  TrackerHitExtendedVec hitVec = trackAR->getTrackerHitExtendedVec();
  int nH = int (hitVec.size());
  return nH;
}

int SiliconTracking::AttachHitToTrack(TrackExtended * trackAR, TrackerHitExtended * hit) {

  int attached = 0;
  TrackerHitExtendedVec hitVec = trackAR->getTrackerHitExtendedVec();
  int nHits = int(hitVec.size());
  double * xh = new double[nHits+1];
  double * yh = new double[nHits+1];
  float * xh_fl = new float[nHits+1];
  float * yh_fl = new float[nHits+1];
  float * rphi_reso = new float[nHits+1];
  float * z_reso = new float[nHits+1];  
  float * zh = new float[nHits+1];
  double * wrh = new double[nHits+1];
  float * wzh = new float[nHits+1];
  float * rh = new float[nHits+1];
  float * ph = new float[nHits+1];
  int * idet_h = new int[nHits+1];
  int * ityp_h = new int[nHits+1];
  float par[5];
  float epar[15];
  
  for (int i=0; i<nHits; ++i) {
    TrackerHit * trkHit = hitVec[i]->getTrackerHit();
    xh[i] = double(trkHit->getPosition()[0]);
    yh[i] = double(trkHit->getPosition()[1]);
    xh_fl[i] = float(xh[i]);
    yh_fl[i] = float(yh[i]);
    zh[i] = float(trkHit->getPosition()[2]);
    ph[i] = float(atan2(yh[i],xh[i]));
    rh[i] = float(sqrt(xh[i]*xh[i]+yh[i]*yh[i]));
    float rR = hitVec[i]->getResolutionRPhi();
    float rZ = hitVec[i]->getResolutionZ();
    wrh[i] = double(1.0/(rR*rR));
    wzh[i] = 1.0/(rZ*rZ);
    rphi_reso[i] = rR;
    z_reso[i] = rZ;
    ityp_h[i] = hitVec[i]->getType();
    idet_h[i] = hitVec[i]->getDet();
  }
  TrackerHit * trkHit = hit->getTrackerHit();
  xh[nHits] = double(trkHit->getPosition()[0]);
  yh[nHits] = double(trkHit->getPosition()[1]);
  xh_fl[nHits] = float(xh[nHits]);
  yh_fl[nHits] = float(yh[nHits]);  
  zh[nHits] = float(trkHit->getPosition()[2]);
  ph[nHits] = float(atan2(yh[nHits],xh[nHits]));
  rh[nHits] = float(sqrt(xh[nHits]*xh[nHits]+yh[nHits]*yh[nHits]));
  float rR = hit->getResolutionRPhi();
  float rZ = hit->getResolutionZ();
  wrh[nHits] = double(1.0/(rR*rR));
  wzh[nHits] = 1.0/(rZ*rZ);
  rphi_reso[nHits] = rR;
  z_reso[nHits] = rZ;
  ityp_h[nHits] = hit->getType();
  idet_h[nHits] = hit->getDet();
  
  int NPT = nHits + 1;
  int IOPT = 3;
  float chi2RPhi;
  float chi2Z;
  float chi2_D;
  int ndf_D;
  if (_simpleHelixFit>0){
    tfithl_(NPT, xh, yh, rh, ph, wrh, zh,
	    wzh, IOPT, par, epar, chi2RPhi, chi2Z);
  }
  else {
    TrackFitting(NPT,_bField,idet_h,ityp_h, _chi2PrefitCut, 
		 xh_fl,yh_fl,zh,rphi_reso,z_reso,
		 par,epar,chi2_D,ndf_D,chi2RPhi,chi2Z);
  }

  
  float omega = par[0];
  float tanlambda = par[1];
  float phi0 = par[2];
  float d0 = par[3];
  float z0 = par[4];
  float chi2;

  if (_simpleHelixFit>0) {
    if (NPT == 3) {
      chi2 = chi2RPhi/_chi2WRPhiTriplet+chi2Z/_chi2WZTriplet;
    }
    if (NPT == 4) {
      chi2 = chi2RPhi/_chi2WRPhiQuartet+chi2Z/_chi2WZQuartet;
    }
    if (NPT > 4) {
      chi2 = chi2RPhi/_chi2WRPhiSeptet+chi2Z/_chi2WZSeptet;
    }
  }
  else 
      chi2 = chi2_D/float(ndf_D);
  
  if (chi2 < _chi2FitCut) {
    trackAR->addTrackerHitExtended(hit);
    hit->addTrackExtended( trackAR );
    if (_simpleHelixFit>0)
	trackAR->setChi2( chi2RPhi );
    else 
	trackAR->setChi2( chi2_D/float(ndf_D));
    trackAR->setOmega( omega );
    trackAR->setTanLambda( tanlambda );
    trackAR->setD0( d0 );
    trackAR->setZ0( z0 );
    trackAR->setPhi( phi0 );
    attached = 1;
  }	
  delete[] xh;
  delete[] yh;
  delete[] zh;
  delete[] wrh;
  delete[] wzh;
  delete[] rh;
  delete[] ph;
  delete[] xh_fl;
  delete[] yh_fl;
  delete[] rphi_reso;
  delete[] z_reso;
  delete[] ityp_h;
  delete[] idet_h;
  

  return attached;


}

void SiliconTracking::FinalRefit() {

  int nTracks = int(_trackImplVec.size());

  for (int iTrk=0;iTrk<nTracks;++iTrk) {
    TrackExtended * trackAR = _trackImplVec[iTrk];
    TrackerHitExtendedVec hitVec = trackAR->getTrackerHitExtendedVec();
    int nHits = int(hitVec.size());
    double * xh = new double[nHits];
    double * yh = new double[nHits];
    float * xh_fl = new float[nHits];
    float * yh_fl = new float[nHits];
    float * rphi_reso = new float[nHits];
    float * z_reso = new float[nHits];  
    float * zh = new float[nHits];
    double * wrh = new double[nHits];
    float * wzh = new float[nHits];
    float * rh = new float[nHits];
    float * ph = new float[nHits];
    int * idet_h = new int[nHits];
    int * ityp_h = new int[nHits];
    float par[5];
    float epar[15];
  
    for (int i=0; i<nHits; ++i) {
      TrackerHit * trkHit = hitVec[i]->getTrackerHit();
      xh[i] = double(trkHit->getPosition()[0]);
      yh[i] = double(trkHit->getPosition()[1]);
      xh_fl[i] = float(xh[i]);
      yh_fl[i] = float(yh[i]);
      zh[i] = float(trkHit->getPosition()[2]);
      ph[i] = float(atan2(yh[i],xh[i]));
      rh[i] = float(sqrt(xh[i]*xh[i]+yh[i]*yh[i]));
      float rR = hitVec[i]->getResolutionRPhi();
      float rZ = hitVec[i]->getResolutionZ();
      wrh[i] = double(1.0/(rR*rR));
      wzh[i] = 1.0/(rZ*rZ);
      rphi_reso[i] = rR;
      z_reso[i] = rZ;
      ityp_h[i] = hitVec[i]->getType();
      idet_h[i] = hitVec[i]->getDet();
    }

    int NPT = nHits;
    float chi2_D,chi2RPhi,chi2Z;
    int ndf_D;

    TrackFitting(NPT,_bField,idet_h,ityp_h, _chi2PrefitCut, 
		 xh_fl,yh_fl,zh,rphi_reso,z_reso,
		 par,epar,chi2_D,ndf_D,chi2RPhi,chi2Z);
    
    float omega = par[0];
    float tanlambda = par[1];
    float phi0 = par[2];
    float d0 = par[3];
    float z0 = par[4];
    float chi2 = chi2_D/float(ndf_D);
    trackAR->setOmega(omega);
    trackAR->setTanLambda(tanlambda);
    trackAR->setPhi(phi0);
    trackAR->setD0(d0);
    trackAR->setZ0(z0);
    trackAR->setCovMatrix(epar);
    trackAR->setChi2(chi2);

    delete[] xh;
    delete[] yh;
    delete[] xh_fl;
    delete[] yh_fl;
    delete[] rphi_reso;
    delete[] z_reso;
    delete[] zh;
    delete[] wrh;
    delete[] wzh;
    delete[] rh;
    delete[] ph;
    delete[] idet_h;
    delete[] ityp_h;

  }


}
