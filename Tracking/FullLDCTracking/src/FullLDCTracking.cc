#include "FullLDCTracking.h"
#include <iostream>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <iostream>
#include <math.h>
#include <map>
#include <marlin/Global.h>
#include "ClusterShapes.h"
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>

using namespace lcio ;
using namespace marlin ;
using namespace UTIL;

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

int LDCTrackFitting(int & nhits, float bField, int * idet, int * itype, float chi2PrefitCut, // inputs 
		    float * x, float * y, float * z, float * RPReso, float * ZReso, // inputs 
		    float * param, float * eparam, float & chi2, int & ndf,float & chi2Sh) { // outputs

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
    xd[i] = x[i];
    yd[i] = y[i];
    phi[i] = atan2(y[i],x[i]);
    r[i] = sqrt(x[i]*x[i]+y[i]*y[i]);
    wfz[i] = 1./(ZReso[i]*ZReso[i]);
    wfr[i] = 1./(RPReso[i]*RPReso[i]);
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
  
  HelixClass  helix;
  ClusterShapes * shapes = new ClusterShapes(nhits,ampl,xhit,yhit,zhit);
  float parSh[5];
  float dparSh[5];
  //  float chi2Sh;
  float distmax;
  shapes->FitHelix(500, 0, 1, parSh, dparSh, chi2Sh, distmax);
  float x0Sh = parSh[0];
  float y0Sh = parSh[1];
  float r0Sh = parSh[2];
  float bzSh = parSh[3];
  float phi0Sh = parSh[4];
  float signPz = 1;
  delete shapes;

  //  std::cout << "DistMax = " << distmax << std::endl;
  if (zEnd<zBegin)
    signPz = -1;
  //  std::cout << "signPz = " << signPz << std::endl;
  helix.Initialize_BZ(x0Sh, y0Sh, r0Sh, 
		      bzSh, phi0Sh, bField,signPz,
		      zBegin);

  float phi0 = helix.getPhi0();
  if (phi0 < 0) phi0 += 2.0*acos(-1.0);
  float z0 = helix.getZ0();
  float d0 = helix.getD0();
  chi2Sh = 0.0;
  float refp[3];
  refp[0] = -d0*sin(phi0);
  refp[1] = d0*cos(phi0);
  refp[2] = z0;
  for (int i=0;i<nhits;++i) {
    if (iTyp[i] == 2) {
      float point[3];
      helix.getPointInZ(zhit[i],refp,point);
      float dX = xhit[i] - point[0];
      float dY = yhit[i] - point[1];
      float dR2 = dX*dX + dY*dY;
      chi2Sh = chi2Sh + dR2/rphireso[i];
    }
    else {
      float point[6];
      float Radius = sqrt(xhit[i]*xhit[i]+yhit[i]*yhit[i]);
      helix.getPointOnCircle(Radius, refp, point);
      float dX1 = xhit[i] - point[0];
      float dY1 = yhit[i] - point[1];
      float dZ1 = zhit[i] - point[2];
      float dX2 = xhit[i] - point[0];
      float dY2 = yhit[i] - point[1];
      float dZ2 = zhit[i] - point[2];
      float dist1 = dX1*dX1+dY1*dY1+dZ1*dZ1;
      float dist2 = dX2*dX2+dY2*dY2+dZ2*dZ2;
      float dX,dY,dZ;
      if (dist1<dist2) {
	dX = dX1;
	dY = dY1;
	dZ = dZ1;
      }
      else {
	dX = dX2;
	dY = dY2;
	dZ = dZ2;
      } 
      float dRSquare = dX*dX + dY*dY;
      float dZSquare = dZ*dZ;
      chi2Sh = chi2Sh + dRSquare/rphireso[i] + dZSquare/zreso[i];
    }
  }
  int ndfSh = 2*nhits-5;

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
    return -1;
  }


  float tanlambda = helix.getTanLambda();
  float omega = helix.getOmega();
  
  // new definitions
//   omega = 10.*param[0];   
//   tanlambda = param[1];
//   phi0 = param[2];
//   d0  = 0.1*param[3];
//   z0  = 0.1*param[4];

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


  //    std::cout << "Delphi Fit : chi2 = " << chi2 
  //	      << " ; ndf = " << ndf
  //	      << " ; ierr = " << ierr << std::endl;

  //  std::cout << std::endl;

  if (ierr != 0) {
    ierr = 1;
    chi2 = 0.1*chi2Sh;
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




FullLDCTracking aFullLDCTracking ;

FullLDCTracking::FullLDCTracking() : Processor("FullLDCTracking") {  
  _description = "Create rootfile with FTD hit info " ;  

  registerInputCollection( LCIO::TRACKERHIT,
                             "FTDHitCollection",
			     "FTD Hit Collection Name",
			     _FTDTrackerHitCollection,
			     std::string("FTDTrackerHits"));  

  registerInputCollection( LCIO::TRACKERHIT,
                             "VXDHitCollection",
			     "VXD Hit Collection Name",
			     _VTXTrackerHitCollection,
			     std::string("VTXTrackerHits"));  

  registerInputCollection( LCIO::TRACKERHIT,
                             "SITHitCollection",
			     "SIT Hit Collection Name",
			     _SITTrackerHitCollection,
			     std::string("SITTrackerHits"));
  
  registerInputCollection( LCIO::TRACKERHIT,
                             "TPCHitCollection",
			     "TPC Hit Collection Name",
			     _TPCTrackerHitCollection,
			     std::string("TPCTrackerHits"));

  registerInputCollection( LCIO::TRACK,
                             "TPCTracks",
			     "TPC Track Collection",
			     _TPCTrackCollection,
			     std::string("TPCTracks"));

  registerInputCollection( LCIO::TRACK,
                             "SiTracks",
			     "Si Track Collection",
			     _SiTrackCollection,
			     std::string("SiTracks"));

  registerInputCollection( LCIO::LCRELATION,
                             "SiTrackMCPRel",
                             "Name of the Si Track MCP Relation collection",
                             _SiTracksMCP,
                             std::string("SiTracksMCP"));
  
  registerInputCollection( LCIO::LCRELATION,
                            "MCTPCTrackRelCollectionName" ,
                            "Name of the TPC Track MCP Relation collection"  ,
                            _MCTPCTracksRel ,
                            std::string("MCTPCTracksRel") ) ;
  
  registerOutputCollection( LCIO::TRACK,
                             "LDCTrackCollection",
                             "Name of the LDC Track Collection",
                             _LDCTracks,
                             std::string("LDCTracks"));
                                                                                                                                                             
  registerOutputCollection( LCIO::LCRELATION,
                             "LDCTrackMCPRelCollection",
                             "Name of the LDC Track MCP Relation collection",
                             _LDCTracksMCPRel,
                             std::string("LDCTracksMCP"));


  registerProcessorParameter("SimpleHelixFit",
			     "Simple Helix Fit?",
			     _simpleHelixFit,
			     int(0));


  registerProcessorParameter("ResolutionRPhi_VTX",
			     "R-Phi Resolution for VTX",
			     _resolutionRPhi_VTX,
			     float(0.004));

  registerProcessorParameter("ResolutionZ_VTX",
			     "Z Resolution for VTX",
			     _resolutionZ_VTX,
			     float(0.004));


  registerProcessorParameter("ResolutionRPhi_FTD",
			     "R-Phi Resolution for FTD",
			     _resolutionRPhi_FTD,
			     float(0.01));

  registerProcessorParameter("ResolutionZ_FTD",
			     "Z Resolution for FTD",
			     _resolutionZ_FTD,
			     float(0.10));

  registerProcessorParameter("ResolutionRPhi_SIT",
			     "R-Phi Resolution for SIT",
			     _resolutionRPhi_SIT,
			     float(0.01));

  registerProcessorParameter("ResolutionZ_SIT",
			     "Z Resolution for SIT",
			     _resolutionZ_SIT,
			     float(0.01));

  registerProcessorParameter("DeltaPhiForTracks",
			     "Delta Phi for Track Merging",
			     _deltaPhiForTracks,
			     float(0.1));

  registerProcessorParameter("DeltaThetaForTracks",
			     "Delta Theta for Track Merging",
			     _deltaQForTracks,
			     float(0.1));

  registerProcessorParameter("AngleCutForMerging",
			     "Angle Cut For Merging",
			     _angleCutForMerging,
			     float(0.2));

  registerProcessorParameter("Chi2FitCut",
			     "Cut on fit Chi2",
			     _chi2FitCut,
			     float(100.0));

  registerProcessorParameter("Chi2PrefitCut",
			     "Cut on fit Chi2",
			     _chi2PrefitCut,
			     float(1.0e+10));

  registerProcessorParameter("CreateMap",
			     "Create Track to MCP Relations",
			     _createMap,
			     int(1));

}



void FullLDCTracking::init() { 

  printParameters();  
  _nRun = -1 ;
  _nEvt = 0 ;
  _bField = Global::parameters->getFloatVal("BField");
  PI = acos(-1.);
  PIOVER2 = 0.5*PI;
  TWOPI = 2*PI;

}

void FullLDCTracking::processRunHeader( LCRunHeader* run) { 
    _nRun++ ;
    _nEvt = 0;
    std::cout << "New Run : run number = " << _nRun << std::endl;
} 

void FullLDCTracking::processEvent( LCEvent * evt ) { 

  std::cout << "Full LDC Tracking Processor" << std::endl;
  prepareVectors( evt );
  std::cout << "prepareVectors done..." << std::endl;
  MergeTPCandSiTracks();
  std::cout << "Merging done..." << std::endl;
  Sorting(_allCombinedTracks);
  std::cout << "Sorting done..." << std::endl;
  SelectCombinedTracks();
  std::cout << "Selection of combined tracks done..." << std::endl;
  AddNotCombinedTracks( evt );
  std::cout << "Not combined tracks added..." << std::endl;

  LCCollectionVec * colTRK = new LCCollectionVec(LCIO::TRACK);

  LCCollectionVec * colRel = NULL;
  if (_createMap > 0)
    colRel = new LCCollectionVec(LCIO::LCRELATION);

  int nTrkCand = int(_trkImplVec.size());

  for (int iTRK=0;iTRK<nTrkCand;++iTRK) {
    TrackExtended * trkCand = _trkImplVec[iTRK];
    TrackImpl * track = new TrackImpl();
    track->setOmega(trkCand->getOmega());
    track->setTanLambda(trkCand->getTanLambda());
    track->setPhi(trkCand->getPhi());
    track->setZ0(trkCand->getZ0());
    track->setD0(trkCand->getD0());
    track->setChi2(trkCand->getChi2());
    TrackerHitExtendedVec hitVec = trkCand->getTrackerHitExtendedVec();
    int nHits = int(hitVec.size());
    std::vector<MCParticle*> mcPointers ;
    std::vector<int> mcHits ;
    mcPointers.clear();
    mcHits.clear();
    for (int iH=0;iH<nHits;++iH) {
      TrackerHitExtended * hitExt = hitVec[iH];
      TrackerHit * hit = hitExt->getTrackerHit();
      track->addHit(hit);
      if (_createMap > 0) {
	int nSH = int(hit->getRawHits().size());
	for (int ish=0;ish<nSH;++ish) {
	  SimTrackerHit * simHit = dynamic_cast<SimTrackerHit*>(hit->getRawHits()[ish]);
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
    GroupTracks * group = trkCand->getGroupTracks();
    if (group != NULL) {
      TrackExtendedVec trkVec = group->getTrackExtendedVec();
      int nGrTRK = int(trkVec.size());
      for (int iGr=0;iGr<nGrTRK;++iGr) {
	TrackExtended * subTrack = trkVec[iGr];
	track->addTrack(subTrack->getTrack());
      }
    }

    colTRK->addElement(track);
    if (_createMap > 0) {
      int nRel = int(mcPointers.size());
      for (int k=0;k<nRel;++k) {
	LCRelationImpl* tpclcRel = new LCRelationImpl;
	MCParticle * mcp = mcPointers[k];
	tpclcRel->setFrom (track);
	tpclcRel->setTo (mcp);
	float weight = (float)(mcHits[k])/(float)(track->getTrackerHits().size());
	tpclcRel->setWeight(weight);
	colRel->addElement(tpclcRel);
      }
    }

  }

  evt->addCollection(colTRK,_LDCTracks);
  if (_createMap > 0)
    evt->addCollection(colRel,_LDCTracksMCPRel);

  CleanUp();
  _nEvt++;
  //  getchar();

}

void FullLDCTracking::prepareVectors(LCEvent * event ) {


  _allTPCHits.clear();
  _allVTXHits.clear();
  _allFTDHits.clear();
  _allSITHits.clear();
  _allTPCTracks.clear();
  _allSiTracks.clear();
  _allCombinedTracks.clear();
  _trkImplVec.clear();

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  //  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;


  std::map <TrackerHit*,TrackerHitExtended*> mapTrackerHits;

  // Reading TPC hits
  try {
    LCCollection * col = event->getCollection(_TPCTrackerHitCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      double tpcRPhiResMax = gearTPC.getDoubleVal("tpcRPhiResMax");      
      double tpcRPhiRes = tpcRPhiResMax-fabs(hit->getPosition()[2])/gearTPC.getMaxDriftLength()*0.010;
      double tpcZRes = gearTPC.getDoubleVal("tpcZRes");
      hitExt->setResolutionRPhi(float(tpcRPhiRes));
      hitExt->setResolutionZ(float(tpcZRes));
      hitExt->setType(int(3));
      hitExt->setDet(int(1));
      _allTPCHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    std::cout << _TPCTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  };

  // Reading FTD Hits
  try {
    LCCollection * col = event->getCollection(_FTDTrackerHitCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      hitExt->setResolutionRPhi(_resolutionRPhi_FTD);
      hitExt->setResolutionZ(_resolutionZ_FTD);
      hitExt->setType(int(2));
      hitExt->setDet(int(2));
      _allFTDHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    std::cout << _FTDTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading SIT Hits
  try {
    LCCollection * col = event->getCollection(_SITTrackerHitCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      hitExt->setResolutionRPhi(_resolutionRPhi_SIT);
      hitExt->setResolutionZ(_resolutionZ_SIT);
      hitExt->setType(int(3));
      hitExt->setDet(int(7));
      _allSITHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    std::cout << _SITTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading VTX Hits
  try {
    LCCollection * col = event->getCollection(_VTXTrackerHitCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
       TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      hitExt->setResolutionRPhi(_resolutionRPhi_VTX);
      hitExt->setResolutionZ(_resolutionZ_VTX);
      hitExt->setType(int(3));
      hitExt->setDet(int(3));
      _allVTXHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    std::cout << _VTXTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }


  // Reading TPC Tracks
  try {
    LCCollection * col = event->getCollection(_TPCTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    std::cout << "Number of TPC Tracks = " << nelem << std::endl;
    for (int iTrk=0; iTrk<nelem; ++iTrk) {
      Track * tpcTrack = dynamic_cast<Track*>(col->getElementAt(iTrk));
      TrackExtended * trackExt = new TrackExtended( tpcTrack );
           TrackerHitVec hitVec = tpcTrack->getTrackerHits();
      int nHits = int(hitVec.size());
      trackExt->setOmega(tpcTrack->getOmega());
      trackExt->setTanLambda(tpcTrack->getTanLambda());
      trackExt->setPhi(tpcTrack->getPhi());
      trackExt->setD0(tpcTrack->getD0());
      trackExt->setZ0(tpcTrack->getZ0());
      // FIXME Chi2 is defined as the total Chi2 of the fit 
      char strg[200];
      sprintf(strg,"%3i %11.8f  %8.3f  %9.6f  %9.3f  %9.3f",iTrk,
	      tpcTrack->getOmega(),
	      tpcTrack->getTanLambda(),
	      tpcTrack->getPhi(),
	      tpcTrack->getD0(),
	      tpcTrack->getZ0());
      std::cout << strg << std::endl;
      trackExt->setChi2(tpcTrack->getChi2());
      for (int iHit=0;iHit<nHits;++iHit) {
	TrackerHit * hit = hitVec[iHit];
	TrackerHitExtended * hitExt = mapTrackerHits[hit];
	hitExt->setTrackExtended( trackExt );
	trackExt->addTrackerHitExtended( hitExt );	
      }
      _allTPCTracks.push_back( trackExt );      
    }    
  }
  catch ( DataNotAvailableException &e) {
    std::cout << _TPCTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading Si Tracks
  try {
    LCCollection * col = event->getCollection(_SiTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    std::cout << "Number of Si Tracks = " << nelem << std::endl;
    for (int iTrk=0; iTrk<nelem; ++iTrk) {
      Track * siTrack = dynamic_cast<Track*>(col->getElementAt(iTrk));
      TrackExtended * trackExt = new TrackExtended( siTrack );
      TrackerHitVec hitVec = siTrack->getTrackerHits();
      int nHits = int(hitVec.size());
      trackExt->setOmega(siTrack->getOmega());
      trackExt->setTanLambda(siTrack->getTanLambda());
      trackExt->setPhi(siTrack->getPhi());
      trackExt->setD0(siTrack->getD0());
      trackExt->setZ0(siTrack->getZ0());
      char strg[200];
      sprintf(strg,"%3i %11.8f  %8.3f  %9.6f  %9.3f  %9.3f",iTrk,
	      siTrack->getOmega(),
	      siTrack->getTanLambda(),
	      siTrack->getPhi(),
	      siTrack->getD0(),
	      siTrack->getZ0());
      std::cout << strg << std::endl;
      trackExt->setChi2(siTrack->getChi2());
      for (int iHit=0;iHit<nHits;++iHit) {
	TrackerHit * hit = hitVec[iHit];
	TrackerHitExtended * hitExt = mapTrackerHits[hit];
	hitExt->setTrackExtended( trackExt );
	trackExt->addTrackerHitExtended( hitExt );	
      }
      _allSiTracks.push_back( trackExt );      
    }    
  }
  catch ( DataNotAvailableException &e) {
    std::cout << _SiTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }



};

void FullLDCTracking::CleanUp(){

  int nSITHits = int(_allSITHits.size());
  for (int i=0;i<nSITHits;++i) {
    TrackerHitExtended * hitExt = _allSITHits[i];
    delete hitExt;
  }
  _allSITHits.clear();

  int nTPCHits = int(_allTPCHits.size());
  for (int i=0;i<nTPCHits;++i) {
    TrackerHitExtended * hitExt = _allTPCHits[i];
    delete hitExt;
  }
  _allTPCHits.clear();

  int nFTDHits = int(_allFTDHits.size());
  for (int i=0;i<nFTDHits;++i) {
    TrackerHitExtended * hitExt = _allFTDHits[i];
    delete hitExt;
  }
  _allFTDHits.clear();

  int nVTXHits = int(_allVTXHits.size());
  for (int i=0;i<nVTXHits;++i) {
    TrackerHitExtended * hitExt = _allVTXHits[i];
    delete hitExt;
  }
  _allVTXHits.clear();

  int nSiTrk = int(_allSiTracks.size());
  for (int i=0;i<nSiTrk;++i) {
    TrackExtended * trkExt = _allSiTracks[i];
    delete trkExt;
  }
  _allSiTracks.clear();

  int nTPCTrk = int(_allTPCTracks.size());
  for (int i=0;i<nTPCTrk;++i) {
    TrackExtended * trkExt = _allTPCTracks[i];
    delete trkExt;
  }
  _allTPCTracks.clear();

  int nCombTrk = int(_allCombinedTracks.size());
  for (int i=0;i<nCombTrk;++i) {
    TrackExtended * trkExt = _allCombinedTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    delete group;
    delete trkExt;    
  }
  _allCombinedTracks.clear();

//   int nImplTrk = int(_trkImplVec.size());
//   for (int i=0;i<nImplTrk;++i) {
//     TrackExtended * trkImpl = _trkImplVec[i];
//     delete trkImpl;
//   }
  _trkImplVec.clear();

};

void FullLDCTracking::MergeTPCandSiTracks() {

  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks  = int(_allSiTracks.size());

  for (int iTPC=0;iTPC<nTPCTracks;++iTPC) {
    TrackExtended * tpcTrackExt = _allTPCTracks[iTPC];
    Track * tpcTrack = tpcTrackExt->getTrack();
    float phiTPC = tpcTrack->getPhi();
    float tanLambdaTPC = tpcTrack->getTanLambda();
    float qTPC = PIOVER2 - atan(tanLambdaTPC);
    for (int iSi=0;iSi<nSiTracks;++iSi) {
      TrackExtended * siTrackExt = _allSiTracks[iSi];
      Track * siTrack = siTrackExt->getTrack();
      float phiSi = siTrack->getPhi();
      float tanLambdaSi = siTrack->getTanLambda();
      float qSi = PIOVER2 - atan(tanLambdaSi);
      float angle = (cos(phiTPC)*cos(phiSi)+sin(phiTPC)*sin(phiSi))*sin(qTPC)*sin(qSi)+cos(qTPC)*cos(qSi);
      angle = acos(angle);

      if (angle < _angleCutForMerging) {
	std::cout << iTPC << " " << iSi << " Angle = " << angle;
	TrackExtended * combinedTrack = CombineTracks(tpcTrackExt,siTrackExt);
	if (combinedTrack != NULL) {
	  _allCombinedTracks.push_back( combinedTrack );
	}
      }
    }
    std::cout << std::endl;
  }


};

TrackExtended * FullLDCTracking::CombineTracks(TrackExtended * tpcTrack, TrackExtended * siTrack) {

  TrackExtended * OutputTrack = NULL;

  TrackerHitExtendedVec siHitVec = siTrack->getTrackerHitExtendedVec();
  TrackerHitExtendedVec tpcHitVec = tpcTrack->getTrackerHitExtendedVec();

  int nSiHits = int(siHitVec.size());
  int nTPCHits = int(tpcHitVec.size());

  int nHits = nTPCHits + nSiHits;

  double * xh = new double[nHits];
  double * yh = new double[nHits];
  float * zh = new float[nHits];
  double * wrh = new double[nHits];
  float * wzh = new float[nHits];
  float * rh = new float[nHits];
  float * ph = new float[nHits];
  float * xh_fl = new float[nHits];
  float * yh_fl = new float[nHits];
  float * rphi_reso = new float[nHits];
  float * z_reso = new float[nHits];
  int * ityp_h = new int[nHits];
  int * idet_h = new int[nHits];
  float par[5];
  float epar[15];
  for (int ih=0;ih<nSiHits;++ih) {
    TrackerHit * trkHit = siHitVec[ih]->getTrackerHit();
    float rR = siHitVec[ih]->getResolutionRPhi();
    float rZ = siHitVec[ih]->getResolutionZ();
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
    ityp_h[ih] =  siHitVec[ih]->getType();
    idet_h[ih] =  siHitVec[ih]->getDet();
  }      
  for (int ih=0;ih<nTPCHits;++ih) {
    TrackerHit * trkHit = tpcHitVec[ih]->getTrackerHit();
    float rR = tpcHitVec[ih]->getResolutionRPhi();
    float rZ = tpcHitVec[ih]->getResolutionZ();
    xh[ih+nSiHits] = trkHit->getPosition()[0];
    yh[ih+nSiHits] = trkHit->getPosition()[1];
    xh_fl[ih+nSiHits] = float(xh[ih+nSiHits]);
    yh_fl[ih+nSiHits] = float(yh[ih+nSiHits]);
    zh[ih+nSiHits] = float(trkHit->getPosition()[2]);
    wrh[ih+nSiHits] = double(1.0/(rR*rR));
    wzh[ih+nSiHits] = 1.0/(rZ*rZ);
    rh[ih+nSiHits] = float(sqrt(xh[ih+nSiHits]*xh[ih+nSiHits]+yh[ih+nSiHits]*yh[ih+nSiHits]));
    ph[ih+nSiHits] = float(atan2(yh[ih+nSiHits],xh[ih+nSiHits]));
    rphi_reso[ih+nSiHits] = rR;
    z_reso[ih+nSiHits] = rZ;
    ityp_h[ih+nSiHits] =  tpcHitVec[ih]->getType();
    idet_h[ih+nSiHits] =  tpcHitVec[ih]->getDet();
  }      

  int NPT = nHits;
  float chi2_D;
  int ndf_D;
  float chi2Sh;
  par[0] = siTrack->getOmega();
  par[1] = siTrack->getTanLambda();
  par[2] = siTrack->getPhi();
  par[3] = siTrack->getD0();
  par[4] = siTrack->getZ0();
  int ierr = LDCTrackFitting(NPT,_bField,idet_h,ityp_h,_chi2PrefitCut,
			     xh_fl,yh_fl,zh,rphi_reso,z_reso,
			     par,epar,chi2_D,ndf_D,chi2Sh);
  

  if (ierr >= 0) {
    float omega = par[0];
    float tanlambda = par[1];
    float phi0 = par[2];
    float d0 = par[3];
    float z0 = par[4];
    
    float chi2Fit = chi2_D/float(ndf_D);

    std::cout << "  Chi2 = " << chi2Fit << " " << chi2Sh << " " << ierr 
	      << " " << par[0] << " " << par[1] << " " << par[2];

    if (chi2Fit > _chi2FitCut) {
      std::cout << std::endl;
      return OutputTrack;
    }    
    else 
      std::cout << "   *" << std::endl;

    OutputTrack = new TrackExtended();

    GroupTracks * group = new GroupTracks();
    group->addTrackExtended(tpcTrack);
    group->addTrackExtended(siTrack);
    OutputTrack->setGroupTracks(group);
    OutputTrack->setOmega(omega);
    OutputTrack->setTanLambda(tanlambda);
    OutputTrack->setPhi(phi0);
    OutputTrack->setZ0(z0);
    OutputTrack->setD0(d0);
    OutputTrack->setChi2(chi2Fit);
    for (int i=0;i<nTPCHits;++i) {
      TrackerHitExtended * hitExt = tpcHitVec[i];
      OutputTrack->addTrackerHitExtended(hitExt);
    }
    for (int i=0;i<nSiHits;++i) {
      TrackerHitExtended * hitExt = siHitVec[i];
      OutputTrack->addTrackerHitExtended(hitExt);
    }
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

  return OutputTrack;

}

void FullLDCTracking::Sorting(TrackExtendedVec & trackVec) {

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
}

void FullLDCTracking::SelectCombinedTracks() {

  int nCombTrk = int(_allCombinedTracks.size());

  for (int i=0; i<nCombTrk;++i) {
    TrackExtended * trkExt = _allCombinedTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    TrackExtendedVec tracks = group->getTrackExtendedVec();
    int nTracks = int(tracks.size());
    if (nTracks == 2) {
      TrackExtended * firstTrack = tracks[0];
      TrackExtended * secondTrack = tracks[1];
      if ((firstTrack->getGroupTracks() == NULL) &&
	  (secondTrack->getGroupTracks() == NULL) ) {
	firstTrack->setGroupTracks(group);
	secondTrack->setGroupTracks(group);	
	_trkImplVec.push_back(trkExt);	
      }
    }
  }


}

void FullLDCTracking::AddNotCombinedTracks(LCEvent * evt) {  

  LCCollection * colTPCRel = evt->getCollection(_MCTPCTracksRel);
  LCCollection * colSiRel  = evt->getCollection(_SiTracksMCP);
  LCRelationNavigator navTPC(colTPCRel);
  LCRelationNavigator navSi(colSiRel);

  
  

  int nTPCTrk = int(_allTPCTracks.size());
  int nSiTrk = int(_allSiTracks.size());

  for (int i=0;i<nTPCTrk;++i) {
    TrackExtended * trkExtTPC = _allTPCTracks[i];
    GroupTracks * groupTPC = trkExtTPC->getGroupTracks();
    if (groupTPC == NULL) {
      for (int j=0;j<nSiTrk;++j) {
	TrackExtended * trkExtSi = _allSiTracks[j];
	GroupTracks * groupSi = trkExtSi->getGroupTracks();
	if (groupSi == NULL) {
	  float omegaTPC = trkExtTPC->getOmega();
	  float omegaSi = trkExtSi->getOmega();
	  float dOmega = fabs((omegaTPC-omegaSi)/omegaTPC);
	  float phiSi = trkExtSi->getPhi();
	  float tanLambdaSi = trkExtSi->getTanLambda();
	  float qSi = PIOVER2 - atan(tanLambdaSi);
	  float phiTPC = trkExtTPC->getPhi();
	  float tanLambdaTPC = trkExtTPC->getTanLambda();
	  float qTPC = PIOVER2 - atan(tanLambdaTPC);
	  float angle = (cos(phiTPC)*cos(phiSi)+sin(phiTPC)*sin(phiSi))*sin(qTPC)*sin(qSi)+cos(qTPC)*cos(qSi);
	  angle = acos(angle);
	  std::cout << i << " " << j << " " << angle << " " << dOmega << std::endl;	  
	  if (dOmega<0.1) {
	    if (angle<0.07) {
	      TrackExtended * OutputTrack = new TrackExtended();
	      GroupTracks * group = new GroupTracks();
	      group->addTrackExtended(trkExtTPC);
	      group->addTrackExtended(trkExtSi);
	      OutputTrack->setGroupTracks(group);
	      trkExtSi->setGroupTracks(group);
	      trkExtTPC->setGroupTracks(group);	      
	      OutputTrack->setOmega(trkExtTPC->getOmega());
	      OutputTrack->setTanLambda(trkExtSi->getTanLambda());
	      OutputTrack->setPhi(trkExtSi->getPhi());
	      OutputTrack->setZ0(trkExtSi->getZ0());
	      OutputTrack->setD0(trkExtSi->getD0());
	      OutputTrack->setChi2(trkExtTPC->getChi2());
	      TrackerHitExtendedVec tpcHitVec = trkExtTPC->getTrackerHitExtendedVec();
	      TrackerHitExtendedVec siHitVec = trkExtSi->getTrackerHitExtendedVec();	      
	      int nTPCHits = int( tpcHitVec.size());
	      int nSiHits = int( siHitVec.size());	      
	      for (int iH=0;iH<nTPCHits;++iH) {
		TrackerHitExtended * hitExt = tpcHitVec[iH];
		OutputTrack->addTrackerHitExtended(hitExt);
	      }
	      for (int iH=0;iH<nSiHits;++iH) {
		TrackerHitExtended * hitExt = siHitVec[iH];
		OutputTrack->addTrackerHitExtended(hitExt);
	      }	  
	      _trkImplVec.push_back(OutputTrack);
	      break;
	    }	    
	  }
	}
      }
    }
  }  

  for (int i=0;i<nTPCTrk;++i) {
    TrackExtended * trkExt = _allTPCTracks[i];
    Track * track = trkExt->getTrack();
    GroupTracks * group = trkExt->getGroupTracks();
    if (group == NULL) {
      char strg[200];
      TrackerHitVec hitVec = track->getTrackerHits();
      sprintf(strg," TPC %3i   %11.8f  %8.3f  %9.6f  %9.3f  %9.3f  %9.3f",i,
	      track->getOmega(),
	      track->getTanLambda(),
	      cos(track->getPhi()),
	      sin(track->getPhi()),
	      track->getD0(),
	      track->getZ0());      
      std::cout << strg << "   " << navTPC.getRelatedToObjects(track)[0] << " " <<  navTPC.getRelatedToWeights(track)[0] << std::endl;
      _trkImplVec.push_back(trkExt);
    }
  }

  std::cout << std::endl;

  for (int i=0;i<nSiTrk;++i) {
    TrackExtended * trkExt = _allSiTracks[i];
    Track * track = trkExt->getTrack();
    GroupTracks * group = trkExt->getGroupTracks();
    if (group == NULL) {
      char strg[200];
      TrackerHitVec hitVec = track->getTrackerHits();
      sprintf(strg," Si  %3i   %11.8f  %8.3f  %9.6f  %9.3f  %9.3f  %9.3f",i,
	      track->getOmega(),
	      track->getTanLambda(),
	      cos(track->getPhi()),
	      sin(track->getPhi()),
	      track->getD0(),
	      track->getZ0());
      std::cout << strg << "   " << navSi.getRelatedToObjects(track)[0] << " " <<  navSi.getRelatedToWeights(track)[0] << std::endl;

      _trkImplVec.push_back(trkExt);
    }
  }

}

void FullLDCTracking::check(LCEvent * evt) { };
void FullLDCTracking::end() {};
    
