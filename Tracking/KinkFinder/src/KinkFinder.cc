#include "KinkFinder.h"
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include "marlin/Global.h"
#include "EVENT/LCCollection.h"
#include <ClusterShapes.h>
#include "EVENT/Track.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/VertexImpl.h"
#include <math.h>

#include <DD4hep/Detector.h>
#include <DD4hep/DD4hepUnits.h>
#include <DDRec/DetectorData.h>

#include "HelixClass.h"

#include <memory>

using namespace lcio ;
using namespace marlin ;

const float mMuon      = 0.106;    // GeV
//const float cTauMuon   = 657000;   // mm
const float mPion      = 0.140;    // GeV
const float cTauPion   = 7800;     // mm
const float mKaon      = 0.494;    // GeV
const float cTauKaon   = 3714;     // mm
const float mSigma     = 1.189;    // GeV
const float cTauSigma  = 24.0;     // mm
const float mLambda    = 1.115;    // GeV
//const float cTauLambda = 78.0;     // mm
const float mHyperon   = 1.314;    // GeV
const float cTauHyperon= 87.0;     // mm
const float mNeutron   = 0.940;    // GeV
const float mProton    = 0.940;    // GeV

KinkFinder aKinkFinder ;


KinkFinder::KinkFinder() : Processor("KinkFinder") {

  _description = "Kink Finder Processor : finds kinks, prongs and splits" ;
  
  registerInputCollection(LCIO::TRACK,
			  "TrackCollection",
			  "Name of input collection of reconstructed particles",
			  _trackColName,
			  std::string("LDCTracks"));

  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "KinkRecoParticleCollection",
			  "Name of output collection of kink reconstructed particles",
			  _kinkRecoPartColName,
			  std::string("KinkRecoParticles"));

  registerOutputCollection(LCIO::VERTEX,
			  "KinkVertexCollection",
			  "Name of output collection of kink vertices",
			  _kinkVertexColName,
			  std::string("KinkVertices"));

  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "ProngRecoParticleCollection",
			  "Name of output collection of prong reconstructed particles",
			  _prongRecoPartColName,
			  std::string("ProngRecoParticles"));

  registerOutputCollection(LCIO::VERTEX,
			  "ProngVertexCollection",
			  "Name of output collection of prong vertices",
			  _prongVertexColName,
			  std::string("ProngVertices"));


  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "SplitRecoParticleCollection",
			  "Name of output collection of split reconstructed particles",
			  _splitRecoPartColName,
			  std::string("SplitRecoParticles"));

  registerOutputCollection(LCIO::VERTEX,
			  "SplitVertexCollection",
			  "Name of output collection of split vertices",
			  _splitVertexColName,
			  std::string("SplitVertices"));


  registerProcessorParameter("SplitTrackMaxFracDeltaP",
			     "Maximum fractional difference in deltaP for split tracks",
			     _maxSplitTrackFracDeltaP,
			     float(0.02));

  registerProcessorParameter("SplitTrackMaxDeltaP",
			     "Maximum difference in deltaP for split tracks",
			     _maxSplitTrackDeltaP,
			     float(0.2));

  registerProcessorParameter("CutOnRadius",
			     "Cuts on Kink radius",
			     _rKinkCut,
			     float(100.0));

  registerProcessorParameter("MinimumTrackHits",
			     "cuts on number of track hits",
			     _minTrackHits,
			     int(5));

  registerProcessorParameter("MaxDeltaTpcLayers",
			     "Cuts on Kink DeltaRxy for tracks in TPC",
			     _maxDeltaTpcLayers,
			     int(10));


  registerProcessorParameter("PionDecayMassCut",
			     "Cuts on kink mass for pion decay",
			     _pionDecayMassCut,
			     float(0.03));

  registerProcessorParameter("KaonDecayMassCut",
			     "Cuts on kink mass for kaon decay",
			     _kaonDecayMassCut,
			     float(0.075));

  registerProcessorParameter("SigmaDecayMassCut",
			     "Cuts on kink mass for sigma decay",
			     _sigmaDecayMassCut,
			     float(0.15));

  registerProcessorParameter("SigmaTimeCut",
			     "Cuts on lifetime for sigma decay",
			     _sigmaTimeCut,
			     float(6.));

  registerProcessorParameter("HyperonDecayMassCut",
			     "Cuts on kink mass for hyperon decay",
			     _hyperonDecayMassCut,
			     float(0.15));

  registerProcessorParameter("HyperonTimeCut",
			     "Cuts on lifetime for hyperon decay",
			     _hyperonTimeCut,
			     float(6.));

  registerProcessorParameter("KinkProjectionCutTPC",
			     "Cuts on kink distance of closest approach in TPC",
			     _drCutTPC,
			     float(20));

  registerProcessorParameter("TightKinkProjectionCutTPC",
			     "Tight Cut on kink distance of closest approach in TPC",
			     _tightDrCutTPC,
			     float(5));

  registerProcessorParameter("VeryTightKinkProjectionCutTPC",
			     "Very Tight Cut on kink distance of closest approach in TPC",
			     _veryTightDrCutTPC,
			     float(1));

  registerProcessorParameter("KinkProjectionCutSIT",
			     "Cuts on kink distance of closest approach in SIT",
			     _drCutSIT,
			     float(10));

  registerProcessorParameter("LooseProjectionCutSIT",
			     "Cuts on kink distance of closest approach in SIT",
			     _looseDrCutSIT,
			     float(10));

  registerProcessorParameter("MinELambda",
			     "Minimum Lambda Energy",
			     _minELambda,
			     float(2.));


  registerProcessorParameter("DebugPrinting",
			     "Debug level",
			     _debugPrinting,
			     int(0));

}

void KinkFinder::init() {

  // Print algorithm parameters
  printParameters() ;


  dd4hep::Detector& theDet = dd4hep::Detector::getInstance();

  double bfieldV[3] ;
  theDet.field().magneticField( { 0., 0., 0. }  , bfieldV  ) ;
  _bField = bfieldV[2]/dd4hep::tesla ;
  
  dd4hep::DetElement tpcDE = theDet.detector("TPC") ;
  const dd4hep::rec::FixedPadSizeTPCData* tpc = tpcDE.extension<dd4hep::rec::FixedPadSizeTPCData>() ;

  _tpcInnerR = tpc->rMinReadout/dd4hep::mm;
  _tpcOuterR = tpc->rMaxReadout/dd4hep::mm;
  _tpcZmax   = tpc->driftLength/dd4hep::mm;
  _tpcMaxRow = tpc->maxRow;


  dd4hep::DetElement vtxDE = theDet.detector("VXD");
  dd4hep::rec::ZPlanarData* vtx = vtxDE.extension<dd4hep::rec::ZPlanarData>();
  _nLayersVTX=vtx->layers.size(); 

  for(int iL=0;iL<_nLayersVTX;iL++){
    _rVTX.push_back( vtx->layers[iL].distanceSensitive/dd4hep::mm );
  }

  dd4hep::DetElement sitDE = theDet.detector("SIT");
  dd4hep::rec::ZPlanarData* sit = sitDE.extension<dd4hep::rec::ZPlanarData>();
  _nLayersSIT=sit->layers.size(); 

  for(int iL=0;iL<_nLayersSIT;iL++){
    _rSIT.push_back( sit->layers[iL].distanceSensitive/dd4hep::mm );
  }


  _nRun = -1;
  _nEvt = 0;

}


void KinkFinder::processRunHeader( LCRunHeader* ) { 

  _nRun++ ;
  _nEvt = 0;

} 

void KinkFinder::processEvent( LCEvent * evt ) { 

  TrackVec tracks;


  if(_debugPrinting>0)std::cout << " Searching for kinks/prongs/split tracks " << std::endl;


  LCCollection* relationTrackCollection = NULL;
  LCRelationNavigator* trackNavigator=NULL;
  try {
    relationTrackCollection = evt->getCollection("LDCTracksMCP");	
    trackNavigator = new LCRelationNavigator(relationTrackCollection); 
  }
  catch(DataNotAvailableException &e){
  }

  try {
    LCCollection * col = evt->getCollection( _trackColName.c_str() );

    int nelem = col->getNumberOfElements();
    TrackPairVec  trkPairs;
    trkPairs.clear();
    std::map<Track*,int> trackUsed;

    for (int i=0;i<nelem;++i) {
      Track * trk = dynamic_cast<Track*>(col->getElementAt(i));
      tracks.push_back(trk);
      trackUsed[trk] = 0;
    }
  }
  catch(DataNotAvailableException &e) {}

  std::vector<vec3>rin;
  std::vector<vec3>rout;
  std::vector<float>rOuter;
  std::vector<float>rInner;
  std::vector<vec3>zin;
  std::vector<vec3>zout;
  std::vector<float>momentum;
  std::vector<float>momentumZ;
  std::vector<float>momentumS;
  std::vector<float>momentumE;
  std::vector<float>charge;
  std::vector<int>hits;
  std::vector<MCParticle*>mcParticle;
  std::vector<HelixClass*> helixEnd(  tracks.size() );
  std::vector<HelixClass*> helixStart(tracks.size() );
  std::vector< std::vector<twoTrackIntersection_t> > kinkDaughters ( tracks.size() ) ;
  std::vector< std::vector<twoTrackIntersection_t> > prongDaughters( tracks.size() ) ;
  std::vector< std::vector<twoTrackIntersection_t> > splitDaughters( tracks.size() ) ;
  std::vector< float> zAtEnd(  tracks.size() );
  std::vector< float> zAtStart(tracks.size() );

  for(unsigned int  itrack=0;itrack< tracks.size();++itrack){
    TrackerHitVec hitvec = tracks[itrack]->getTrackerHits();
    int nhits = (int)hitvec.size();
    hits.push_back(nhits);
    int outer = 0;
    float router =-99999.;
    int inner = 0;
    float rinner =99999.;
    int max = 0;
    int min = 0;
    float abszmax =-99999.;
    float abszmin =99999.;
    float zmax =-99999.;
    float zmin =99999.;
    for(int ih =0;ih<nhits;++ih){
      float x = (float)hitvec[ih]->getPosition()[0];
      float y = (float)hitvec[ih]->getPosition()[1];
      float z = (float)hitvec[ih]->getPosition()[2];
      float r2 = x*x+y*y;
      float  r = sqrt(r2);
      if(z<zmin)zmin=z;
      if(z>zmax)zmax=z;
      if(fabs(z)<abszmin){
	min=ih;
	abszmin = fabs(z);
      }
      if(fabs(z)>abszmax){
	max=ih;
	abszmax = fabs(z);
      }
      if(r<rinner && r>15.){
	inner=ih;
	rinner = r;
      }
      if(r>router){
	outer=ih;
	router = r;
      }
    }
    vec3 ri;
    ri.x = (float)hitvec[inner]->getPosition()[0];
    ri.y = (float)hitvec[inner]->getPosition()[1];
    ri.z = (float)hitvec[inner]->getPosition()[2];
    rin.push_back(ri);
    vec3 ro;
    ro.x = (float)hitvec[outer]->getPosition()[0];
    ro.y = (float)hitvec[outer]->getPosition()[1];
    ro.z = (float)hitvec[outer]->getPosition()[2];
    rout.push_back(ro);
    vec3 zi;
    zi.x = (float)hitvec[min]->getPosition()[0];
    zi.y = (float)hitvec[min]->getPosition()[1];
    zi.z = (float)hitvec[min]->getPosition()[2];
    vec3 zo;
    zo.x = (float)hitvec[max]->getPosition()[0];
    zo.y = (float)hitvec[max]->getPosition()[1];
    zo.z = (float)hitvec[max]->getPosition()[2];
    zin.push_back(zi);
    zout.push_back(zo);
    rInner.push_back(rinner);
    rOuter.push_back(router);
    float zBegin, zEnd;
    if (fabs(zmin)<fabs(zmax)) {
      zBegin = zmin;
      zEnd   = zmax;
    }else {
      zBegin = zmax;
      zEnd   = zmin;
    }
    float signPz = zEnd - zBegin;		  
    zAtEnd[itrack] = zEnd;
    zAtStart[itrack] = zBegin;

    int _nHitsInFit = 50;
    int nhitsFit;
    if (nhits > _nHitsInFit) {      
      for (int iz = 0 ; iz < nhits-1; iz++)
	for (int jz = 0; jz < nhits-iz-1; jz++){
	  TrackerHit * one = hitvec[jz];
	  TrackerHit * two = hitvec[jz+1];
	  float dist1 = fabs(float(one->getPosition()[2])-zBegin);
	  float dist2 = fabs(float(two->getPosition()[2])-zBegin);
	  if(dist1 > dist2)
	    {
	      TrackerHit * Temp = hitvec[jz];
	      hitvec[jz] = hitvec[jz+1];
	      hitvec[jz+1] = Temp;
	    }
	}
      nhitsFit = _nHitsInFit;
    }else {
      nhitsFit = nhits;
    }

    HelixClass helix;
    float phi0  = tracks[itrack]->getPhi();
    float d0    = tracks[itrack]->getD0();
    float omega = tracks[itrack]->getOmega();
    float z0    = tracks[itrack]->getZ0();
    float tanlambda = tracks[itrack]->getTanLambda();
    helix.Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _bField);  
    float totMomentum =0; 
    float Mom[3]; 
    for (int j=0; j < 3; ++j) {
      Mom[j] = helix.getMomentum()[j];
      totMomentum += Mom[j]*Mom[j];
    }
    totMomentum = sqrt(totMomentum);
    momentum.push_back(totMomentum);
    momentumZ.push_back(Mom[2]);
    charge.push_back(helix.getCharge());

    MCParticle* mcpart = NULL;
    if(trackNavigator!=NULL){
      LCObjectVec objectVec = trackNavigator->getRelatedToObjects(tracks[itrack]);
      if (objectVec.size() > 0)mcpart = dynamic_cast<MCParticle*>(objectVec[0]);
      if (objectVec.size() > 1){
	float dbest = 999999999.;
	for(unsigned im=0;im<objectVec.size();im++ ){
	  MCParticle* mcparti = dynamic_cast<MCParticle*>(objectVec[im]);
	  float p[3];
	  if(mcparti!=NULL){
	    p[0] = mcparti->getMomentum()[0];
	    p[1] = mcparti->getMomentum()[1];
	    p[2] = mcparti->getMomentum()[2];
	    float ptot = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
	    float d = fabs(ptot-totMomentum);
	    if(d<dbest){
	      dbest = d;
	      mcpart = mcparti;
	    }
	  }
	}
      }
    }
    mcParticle.push_back(mcpart);

    float * xh = new float[nhitsFit];
    float * yh = new float[nhitsFit];
    float * zh = new float[nhitsFit];
    float * ah = new float[nhitsFit]; 
    for (int i=0; i<nhitsFit; ++i) {
      ah[i] = 0;
      xh[i] = float(hitvec[i]->getPosition()[0]);
      yh[i] = float(hitvec[i]->getPosition()[1]);
      zh[i] = float(hitvec[i]->getPosition()[2]);
    }      

    ClusterShapes * shapesS = new ClusterShapes(nhitsFit,ah,xh,yh,zh);
    float par[5];
    float dpar[5];
    float chi2;
    float distmax;
    shapesS->FitHelix(500, 0, 1, par, dpar, chi2, distmax);

    delete shapesS;
    float x0 = par[0];
    float y0 = par[1];
    float r0 = par[2];
    float bz = par[3];
    float phiH = par[4]; 
    helixStart[itrack] = new HelixClass();
    helixStart[itrack]->Initialize_BZ(x0, y0, r0, 
			 bz, phiH, _bField,signPz,
			 zBegin);


    float seeds[3];
    float refs[3];
    refs[0]  = helixStart[itrack]->getReferencePoint()[0];
    refs[1]  = helixStart[itrack]->getReferencePoint()[1];
    refs[2]  = helixStart[itrack]->getReferencePoint()[2];

    helixStart[itrack]->getPointInZ(zAtEnd[itrack], refs, seeds);
    for (int i=0; i<nhitsFit; ++i) {
      ah[i] = 0;
      xh[i] = float(hitvec[nhits-1-i]->getPosition()[0]);
      yh[i] = float(hitvec[nhits-1-i]->getPosition()[1]);
      zh[i] = float(hitvec[nhits-1-i]->getPosition()[2]);
    }      
    ClusterShapes * shapesE = new ClusterShapes(nhitsFit,ah,xh,yh,zh);
    shapesE->FitHelix(500, 0, 1, par, dpar, chi2, distmax);
    delete shapesE;
    x0 = par[0];
    y0 = par[1];
    r0 = par[2];
    bz = par[3];
    phiH = par[4]; 
    helixEnd[itrack] = new HelixClass();
    helixEnd[itrack]->Initialize_BZ(x0, y0, r0, 
			 bz, phiH, _bField,signPz,
			 zEnd);

    delete[] xh;
    delete[] yh;
    delete[] zh;
    delete[] ah;

  }

  for(unsigned int i=0;i< tracks.size();++i){ 
//    Track* tracki = tracks[i];
//    float d0i = tracki->getD0();
//    float z0i = tracki->getZ0();
    if(hits[i]>=_minTrackHits){
      float seedi[3];
      float refe[3];
      float z = zAtEnd[i];
      refe[0]  = helixEnd[i]->getReferencePoint()[0];
      refe[1]  = helixEnd[i]->getReferencePoint()[1];
      refe[2]  = helixEnd[i]->getReferencePoint()[2];
      helixEnd[i]->getPointInZ(z, refe, seedi);
      float rendi = sqrt(seedi[0]*seedi[0]+seedi[1]*seedi[1]);
      for(unsigned int  j=0;j< tracks.size();++j){
	if(i!=j && hits[j]>_minTrackHits && momentum[j] < (1.0 + _maxSplitTrackFracDeltaP)*momentum[i]){
//	  Track* trackj = tracks[j];
//        float d0j = trackj->getD0();
//	  float z0j = trackj->getZ0();
	  if(rInner[i]>_rKinkCut || rInner[j]>_rKinkCut){
	    float seedj[3];
	    float refs[3];
	    float deltaz;
	    float ddx;
	    float ddy;
	    float ddz;
	    float deltar;
	    float deltarxy;
	    int ip;
	    int id;
	    float z1;
	    float z2;
	    bool flipped = false;
	    if(fabs(zAtEnd[i]-zAtStart[j])<=fabs(zAtEnd[i]-zAtEnd[j])){
	      // normal kink where both tracks propagate outwards
	      ip = i;
	      id = j;
	      refs[0]  = helixStart[j]->getReferencePoint()[0];
	      refs[1]  = helixStart[j]->getReferencePoint()[1];
	      refs[2]  = helixStart[j]->getReferencePoint()[2];
	      helixStart[j]->getPointInZ(z, refs, seedj);
	      deltaz =  zAtEnd[i] - zAtStart[j];
	      ddx = (zout[i].x-zin[j].x);
	      ddy = (zout[i].y-zin[j].y);
	      ddz = (zout[i].z-zin[j].z);
	      deltar = sqrt(ddx*ddx+ddy*ddy+ddz*ddz);
	      deltarxy = fabs(rInner[j]-rOuter[i]);
	      z1 = zAtEnd[ip];
	      z2 = zAtStart[id];
	    }else{
	      // kink where daughter track propagates inwards
	      // or a V0 !!!
	      // or a backscattered track
	      flipped = true;
	      ip = j;
	      id = i;
	      refs[0]  = helixEnd[j]->getReferencePoint()[0];
	      refs[1]  = helixEnd[j]->getReferencePoint()[1];
	      refs[2]  = helixEnd[j]->getReferencePoint()[2];
	      helixEnd[j]->getPointInZ(z, refs, seedj);	      
	      deltaz =  zAtEnd[i] - zAtEnd[j];
	      ddx = (zout[i].x-zout[j].x);
	      ddy = (zout[i].y-zout[j].y);
	      ddz = (zout[i].z-zout[j].z);
	      deltar = sqrt(ddx*ddx+ddy*ddy+ddz*ddz);
	      deltarxy = fabs(rInner[i]-rOuter[j]);
	      z1 = zAtEnd[ip];
	      z2 = zAtEnd[id];
	    }
	    float dx = seedi[0]-seedj[0];
	    float dy = seedi[1]-seedj[1];
	    float dz = seedi[2]-seedj[2];
	    float dr = sqrt(dx*dx+dy*dy+dz*dz);
	    
	    float xkink = 0;
	    float ykink = 0;
	    float zkink = 0;
	    float rkink = 0;


	    bool mcKink = false;
	    if(mcParticle[j]!=NULL && _debugPrinting>0){
	      EVENT::MCParticleVec parents = mcParticle[j]->getParents();
	      if(parents.size()>0 && mcParticle[i]!=NULL){
		for(unsigned im = 0; im<parents.size();im++){
		  if(parents[im]==mcParticle[i] && momentum[i]>2.0){
		    //		    if(abs(mcParticle[i]->getPDG())==211 &&
		    //  abs(mcParticle[j]->getPDG())==13)mcKink = true;
		    if(abs(mcParticle[i]->getPDG())==321 &&
		       abs(mcParticle[j]->getPDG())==13)mcKink = true;
		    if(abs(mcParticle[i]->getPDG())==321 &&
		       abs(mcParticle[j]->getPDG())==211)mcKink = true;
		    if(abs(mcParticle[i]->getPDG())==3222 &&
		       abs(mcParticle[j]->getPDG())==211)mcKink = true;
		    if(abs(mcParticle[i]->getPDG())==3212 &&
		       abs(mcParticle[j]->getPDG())==211)mcKink = true;
		    if(abs(mcParticle[i]->getPDG())==3112 &&
		       abs(mcParticle[j]->getPDG())==211)mcKink = true;
		  }
		}
	      }
	    }
	  
	  

	    if(dr<100 || mcKink){
	      dr = 999999.;
	      float refsp[3];
	      float refsd[3];
	      refsp[0]  = helixEnd[ip]->getReferencePoint()[0];
	      refsp[1]  = helixEnd[ip]->getReferencePoint()[1];
	      refsp[2]  = helixEnd[ip]->getReferencePoint()[2];
	      refsd[0]  = helixStart[id]->getReferencePoint()[0];
	      refsd[1]  = helixStart[id]->getReferencePoint()[1];
	      refsd[2]  = helixStart[id]->getReferencePoint()[2];
	      float seedp[3];
	      float seedd[3];
	      float zs = z1;
	      float ze = z2;
	      if(z2 < z1){
		zs = z2;
		ze = z1;
	      }
	      float zstep = (ze-zs)/10.;
	      if(zstep<1)zstep=1;
	      
	      for(float zf = zs; zf < ze; zf+= zstep){
		helixEnd[ip]->getPointInZ(zf, refsp, seedp);
		helixStart[id]->getPointInZ(zf, refsd, seedd);
 		float dxf = seedp[0]-seedd[0];
		float dyf = seedp[1]-seedd[1];
		float dzf = seedp[2]-seedd[2];
		float newdr = sqrt(dxf*dxf+dyf*dyf+dzf*dzf);
		if(newdr<dr){
		  dr=newdr;
		  xkink = (seedp[0]+seedd[0])/2;
		  ykink = (seedp[1]+seedd[1])/2;
		  rkink = sqrt(xkink*xkink+ykink*ykink);
		  zkink = zf;
		  if(rkink<rOuter[i])rkink=rOuter[i];
		  if(rkink>rInner[j])rkink=rInner[j];
		}
	      }
	    }
	    
	    bool  ok = true;
	    if(fabs(deltaz)>200)ok=false;
	    if(fabs(deltaz)>100 && dr > 5.0)ok=false;
	    
	    float deltaRxyCut = -100;
	    float drCut   = -100;
//    	    bool goodRadialSep = false;
	    // require kink to be outside VTX detector
	    if(rkink > _rVTX[_nLayersVTX-1]){
	      // for TPC use
	      drCut  = _drCutTPC;
	      deltaRxyCut= (_tpcOuterR-_tpcInnerR)/_tpcMaxRow*_maxDeltaTpcLayers;
	      if(dr<_tightDrCutTPC)deltaRxyCut=deltaRxyCut*1.5;
	      if(dr<_veryTightDrCutTPC)deltaRxyCut=deltaRxyCut/1.5*2.0;

	      // for SIT
	      if(rkink<_tpcInnerR+deltaRxyCut){
		drCut  = _drCutSIT;
		int iSitLayer=0;
		for(int il=0;il<_nLayersSIT;il++){
		  if(rkink>_rSIT[il])iSitLayer=il+1;
		}
		// kink between VTX and SIT
		if(iSitLayer==0)deltaRxyCut = _rSIT[0]-_rVTX[_nLayersVTX-1];
		// kink between SIT and TPC
		if(iSitLayer==_nLayersSIT){
		  deltaRxyCut = _tpcInnerR - _rSIT[_nLayersSIT-1];
		  drCut  = _drCutTPC;
		}
		if(iSitLayer>0 && iSitLayer <_nLayersSIT)deltaRxyCut = _rSIT[iSitLayer] - _rSIT[iSitLayer-1];
		// add some protection for SIT layers - large gap
		int ili = -999;
		int ilo = -999;
		for(int il=0; il<_nLayersSIT;il++){
		  if(fabs(rOuter[i]-_rSIT[il])<10.)ili=il;
		  if(fabs(rInner[j]-_rSIT[il])<10.)ilo=il;
		}
		if(ili>=0&&ilo>=0){
		  float fix = (_rSIT[ilo] - _rSIT[ili]);
		  if(fix<10)fix=10.;
		  if(fix>deltaRxyCut)deltaRxyCut = fix;
		}
		// add some slop to account for geometry
		deltaRxyCut = deltaRxyCut*1.1; 
	      }
	    }

	      
	    // looser cuts for debug
	    //	    std::cout << i << " : " << j << " dr = " << dr << " ( " << drCut << " )    deltaRxy = " << deltarxy << " ( " << deltaRxyCut << " ) " << std::endl; 
	    if( (dr<drCut && deltarxy < deltaRxyCut*2) || mcKink){
	      bool possibleSplit = false;
	      bool split = false;
	      rkink = sqrt(xkink*xkink+ykink*ykink);

	      if( (rkink > _rKinkCut && !flipped) || mcKink){

		float massENu   = this->kinkMass(helixEnd[i],helixStart[j],0., 0.);
		float massMuNu  = this->kinkMass(helixEnd[i],helixStart[j],mMuon,0.);
		float massPiPi  = this->kinkMass(helixEnd[i],helixStart[j],mPion,mPion);
		float massPiN   = this->kinkMass(helixEnd[i],helixStart[j],mPion,mNeutron);
		float massPPi0  = this->kinkMass(helixEnd[i],helixStart[j],mProton,mPion);
		float massPiL   = this->kinkMass(helixEnd[i],helixStart[j],mPion,mLambda);
		float FmassENu   = this->kinkMass(helixEnd[i],helixEnd[j],0., 0.);
		float FmassMuNu  = this->kinkMass(helixEnd[i],helixEnd[j],mMuon,0.);
		float FmassPiPi  = this->kinkMass(helixEnd[i],helixEnd[j],mPion,mPion);
		float FmassPiN   = this->kinkMass(helixEnd[i],helixEnd[j],mPion,mNeutron);
		float FmassPPi0  = this->kinkMass(helixEnd[i],helixEnd[j],mProton,mPion);
		float FmassPiL   = this->kinkMass(helixEnd[i],helixEnd[j],mPion,mLambda);
		float tPion = fabs(zkink-zAtStart[i])*mPion/fabs(momentumZ[i])/cTauPion;
		float tSigma = fabs(zkink-zAtStart[i])*mSigma/fabs(momentumZ[i])/cTauSigma;
		float tKaon  = fabs(zkink-zAtStart[i])*mKaon/fabs(momentumZ[i])/cTauKaon;
		float tHyperon  = fabs(zkink-zAtStart[i])*mHyperon/fabs(momentumZ[i])/cTauHyperon;
		float probPion    = 0;
		float probKaon    = 0;
		float probSigma   = 0;
		float probHyperon = 0;
		if(_pionDecayMassCut>0){
		  float deltaPion = fabs(massMuNu-mPion)/_pionDecayMassCut; 
		  if(deltaPion<1 && tPion > 0.001)probPion = 3.125*deltaPion*deltaPion+tPion;
		}
		if(_kaonDecayMassCut>0){
		  float deltaKaonMuNu = fabs(massMuNu-mKaon)/_kaonDecayMassCut; 
		  float deltaKaonPiPi = fabs(massPiPi-mKaon)/_kaonDecayMassCut;
		  float deltaKaon = deltaKaonMuNu;
		  if(deltaKaonPiPi<deltaKaon)deltaKaon = deltaKaonPiPi;
		  if(deltaKaon<1 && tKaon > 0.005)probKaon = 3.125*deltaKaon*deltaKaon+tKaon;
		}
		if(_sigmaDecayMassCut>0 && momentum[i] > _minELambda){
		  float deltaSigmaPiN  = fabs(massPiN-mSigma)/_sigmaDecayMassCut; 
		  float deltaSigmaPPi0 = fabs(massPPi0-mSigma)/_sigmaDecayMassCut;
		  float deltaSigma     = deltaSigmaPiN;
		  if(deltaSigmaPPi0<deltaSigma)deltaSigma = deltaSigmaPPi0;
		  if(deltaSigma<1 && tSigma < _sigmaTimeCut)probSigma = 3.125*deltaSigma*deltaSigma+tSigma;
		}
		if(_hyperonDecayMassCut>0 && momentum[i] > _minELambda){
		  float deltaHyperon   = fabs(massPiL-mHyperon)/_hyperonDecayMassCut; 
		  if(deltaHyperon<1 && tHyperon < _hyperonTimeCut)probHyperon = 3.125*deltaHyperon*deltaHyperon+tHyperon;
		}
		bool goodKinkMass = false;
		if(probPion>0.0001||probKaon>0.0001||probSigma>0.0001||probHyperon>0.0001)goodKinkMass=true;


		float dpop = 2*fabs(momentum[i]-momentum[j])/(momentum[i]+momentum[j]); 
		float dp   = fabs(momentum[i]-momentum[j]);
		if(dr<drCut && dpop < _maxSplitTrackFracDeltaP && dp < _maxSplitTrackDeltaP){
		  possibleSplit = true;
		}

		if( (deltarxy<deltaRxyCut && dr<drCut) || possibleSplit ){
		  twoTrackIntersection_t kinkij;
		  kinkij.tracki   = i;
		  kinkij.trackj   = j;
		  kinkij.vtx[0]   = xkink;
		  kinkij.vtx[1]   = ykink;
		  kinkij.vtx[2]   = zkink;
		  kinkij.p[0]     = helixStart[i]->getMomentum()[0];
		  kinkij.p[1]     = helixStart[i]->getMomentum()[1];
		  kinkij.p[2]     = helixStart[i]->getMomentum()[2];
		  kinkij.mass     = 0.;
		  kinkij.distance = dr;

		  // split tracks
		  if(deltarxy<2*deltaRxyCut && dr<drCut*2 ){
		    if(possibleSplit && charge[i]*charge[j]>0 ){
		      TrackerHitVec hitveci= tracks[i]->getTrackerHits();
		      TrackerHitVec hitvecj= tracks[j]->getTrackerHits();
		      int nhitsi = (int)hitveci.size();
		      int nhitsj = (int)hitvecj.size();
		      int ntpci = 0;
		      int ntpcj = 0;
		      HelixClass* helixi = helixStart[j];
		      HelixClass* helixj = helixEnd[i]; 
		      float hitxyz[3];
		      float dist[3];
		      float maxdisti=0;
		      float maxdistj=0;
		      int nclosei = 0;
		      int nclosej = 0;
		      float zmini = 99999;
		      float zminj = 99999;
		      float zmaxi = -99999;
		      float zmaxj = -99999;
		      for(int ih =0;ih<nhitsi;++ih){
			float x = (float)hitveci[ih]->getPosition()[0];
			float y = (float)hitveci[ih]->getPosition()[1];
			float zz = (float)hitveci[ih]->getPosition()[2];
			if(fabs(zz)<zmini)zmini=fabs(zz);
			if(fabs(zz)>zmaxi)zmaxi=fabs(zz);
			float r2 = x*x+y*y;
			float  r = sqrt(r2);
			if(r>_tpcInnerR)ntpci++;
			hitxyz[0]=x;
			hitxyz[1]=y;
			hitxyz[2]=zz;
			helixj->getDistanceToPoint(hitxyz, dist);
			if(dist[2]>maxdisti)maxdisti=dist[2];
			if(dist[2]<25.)nclosei++;
		      }
		      for(int ih =0;ih<nhitsj;++ih){
			float x = (float)hitvecj[ih]->getPosition()[0];
			float y = (float)hitvecj[ih]->getPosition()[1];
			float zz = (float)hitvecj[ih]->getPosition()[2];
			if(fabs(zz)<zminj)zminj=fabs(zz);
			if(fabs(zz)>zmaxj)zmaxj=fabs(zz);
			float r2 = x*x+y*y;
			float  r = sqrt(r2);
			if(r>_tpcInnerR)ntpcj++;
			hitxyz[0]=x;
			hitxyz[1]=y;
			hitxyz[2]=zz;
			helixi->getDistanceToPoint(hitxyz, dist);
			if(dist[2]>maxdistj)maxdistj=dist[2];
			if(dist[2]<25.)nclosej++;
		      }
		      float fclosei = (float)nclosei/(float)nhitsi;
		      float fclosej = (float)nclosej/(float)nhitsj;
		      if(_debugPrinting>0){
			std::cout << " CAND SPLIT I : " << nhitsi << " ntpc " << ntpci << " nclose " << nclosei << " max " << maxdisti << " fclose : " << fclosei << std::endl; 
			std::cout << " CAND SPLIT J : " << nhitsj << " ntpc " << ntpcj << " nclose " << nclosej << " max " << maxdistj << " fclose : " << fclosej << std::endl; 
		      }
		      if(maxdistj<50 && maxdisti < 50 && fclosej > 0.95 && fclosei > 0.95 && ntpcj+ntpci < _tpcMaxRow+10.)split = true;
		      splitDaughters[i].push_back(kinkij);
		    }
		  }

		  // kinks
		  if(deltarxy<deltaRxyCut && dr<drCut ){
		    if(charge[i]*charge[j]>0 && goodKinkMass){
		      if(_debugPrinting>0)std::cout << "Found  kink candidate " << i << " " << j  << std::endl;
		      if(probPion>0||probKaon>0||probSigma>0||probHyperon>0)goodKinkMass=true;
		      if(probPion > probKaon  && 
			 probPion > probSigma &&
			 probPion > probHyperon){
			kinkij.pdgCode = 211;
			if(charge[i]<0)kinkij.pdgCode = -211;
			kinkij.mass = mPion;
		      }
		      if(probKaon > probPion  && 
			 probKaon > probSigma &&
			 probKaon > probHyperon){
			kinkij.mass = mKaon;
			kinkij.pdgCode = 321;
			if(charge[i]<0)kinkij.pdgCode = -321;
		      }
		      if(probSigma > probPion  && 
			 probSigma > probKaon &&
			 probSigma > probHyperon){
			kinkij.mass = mSigma;
			kinkij.pdgCode = 3222;
			if(charge[i]<0)kinkij.pdgCode = 3222;
		      }
		      if(probHyperon > probPion  && 
			 probHyperon > probKaon  &&
			 probHyperon > probSigma){
			kinkij.mass = mHyperon;
			kinkij.pdgCode = 3312;
			if(charge[i]>0)kinkij.pdgCode = -3312;
		      }
		      kinkDaughters[i].push_back(kinkij);
		    }
		  }

		  // prongs
		  if(deltarxy<deltaRxyCut && dr<drCut ){
		
		    if(_debugPrinting>0)std::cout << "Found prong candidate " << i << " " << j  << std::endl;
		    prongDaughters[i].push_back(kinkij);
		    kinkij.pdgCode = 211;
		    if(charge[i]<0)kinkij.pdgCode = -211;
		    kinkij.mass = mPion;
		  }

		}
		
		if(_debugPrinting>0){
		  if(mcParticle[j]!=NULL){
		    EVENT::MCParticleVec parents = mcParticle[j]->getParents();
		    if(parents.size()>0 && mcParticle[i]!=NULL){
		      for(unsigned im = 0; im<parents.size();im++){
			if(parents[im]==mcParticle[i])std::cout << " TRUE KINK . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . " << std::endl;
		      }
		    }
		  }
		  


		  std::cout << "   Candidate kink tracks : " << i << "," << j << " p[i] " << momentum[i] << " p[j] " << momentum[j] << std::endl;
		  if(flipped)std::cout << "   Flipped Track " << std::endl;;
		  std::cout << "   MC  : ";
		  if(mcParticle[i]!=NULL)std::cout << mcParticle[i]->getPDG();
		  std::cout << " - "; 
		  if(mcParticle[j]!=NULL)std::cout << mcParticle[j]->getPDG();
		  std::cout << std::endl;	
		  std::cout << "   Mass                  : " << massENu << " " << massMuNu << " " << massPiPi << " " << massPiN << " " << massPPi0 << " " << massPiL << std::endl;  
		  if(flipped)std::cout << "   MassF                 : " << FmassENu << " " << FmassMuNu << " " << FmassPiPi << " " << FmassPiN << " " << FmassPPi0 << FmassPiL << std::endl;  

		  std::cout << "     ProbPion            : " << probPion << " TimePion    : " << tPion  << " Dm = " << massMuNu-mPion << std::endl;
		  std::cout << "     ProbKaon            : " << probKaon << " TimeKaon    : " << tKaon  << 
		    " Dm = " << massMuNu-mKaon << 
		    " Dm = " << massPiPi-mKaon <<  std::endl;
		  std::cout << "     ProbSigma           : " << probSigma << " TimeSigma   : " << tSigma << " Dm = " << massPiN-mSigma << 
		    " Dm = " << massPPi0-mSigma <<  std::endl;
		  std::cout << "     ProbHyperon         : " << probHyperon << " TimeHyperon : " << tHyperon << " Dm = " << massPiL-mHyperon <<  std::endl;
		  std::cout << "   Charge                : " << i << "," << j << " q[i] " << charge[i] << " q[j] " << charge[j] << std::endl;
		  std::cout << "   Rendi                 : " << i << "," << j << "  rend " << rendi << "  dr = " << dr << " (deltaR = " << deltar << ")" << std::endl;
		  std::cout << "   zkink                 : " <<  zkink << " rzykink " << rkink << std::endl;
		  if(deltarxy<deltaRxyCut)std::cout << "   DeltaRxy vs cut       : " << deltarxy << " < " << deltaRxyCut << std::endl;
		  if(deltarxy>deltaRxyCut)std::cout << "   DeltaRxy vs cut       : " << deltarxy << " > " << deltaRxyCut << std::endl;
		  if(dr<drCut)std::cout <<             "   dr vs cut             : " << dr << " < " << drCut << std::endl;
		  if(dr>drCut)std::cout <<             "   dr vs cut             : " << dr << " > " << drCut << std::endl;
		  std::cout << "   zi/zj                 : " <<  zAtStart[i] << " - " << zAtEnd[i] << "     " << zAtStart[j] << " - " << zAtEnd[j] << std::endl;
		  std::cout << "   rinner/router         : " << rInner[i] << " - " << rOuter[i]    << "     " << rInner[j]   << " - " << rOuter[j] << std::endl;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  int countk = 0;
  int countp = 0;
  for(unsigned int  i=0;i< tracks.size();++i){ 


    if(kinkDaughters[i].size()>0 || prongDaughters[i].size()>0){
      if(kinkDaughters[i].size()==1 && prongDaughters[i].size()<2){
	countk++;
      }else{
	countp++;
      }
    }
  }

  auto colKinkRecoPart  = std::make_unique<LCCollectionVec>(LCIO::RECONSTRUCTEDPARTICLE);
  auto colKinkVertex    = std::make_unique<LCCollectionVec>(LCIO::VERTEX);
  auto colProngRecoPart = std::make_unique<LCCollectionVec>(LCIO::RECONSTRUCTEDPARTICLE);
  auto colProngVertex   = std::make_unique<LCCollectionVec>(LCIO::VERTEX);
  auto colSplitRecoPart = std::make_unique<LCCollectionVec>(LCIO::RECONSTRUCTEDPARTICLE);
  auto colSplitVertex   = std::make_unique<LCCollectionVec>(LCIO::VERTEX);

  for(unsigned int i=0;i< tracks.size();++i){ 
    if(kinkDaughters[i].size()>0 || prongDaughters[i].size()>0 || splitDaughters[i].size()>0){
      if(_debugPrinting>0){
	if(kinkDaughters[i].size()>0){
	  std::cout  << " Track " << i << " has " << kinkDaughters[i].size() << " kink daughters : ";
	  for(unsigned j=0;j<kinkDaughters[i].size();j++)std::cout << kinkDaughters[i][j].trackj << " ";
	  std::cout  << std::endl;
	}
	if(prongDaughters[i].size()>0){
	  std::cout  << " Track " << i << " has " << prongDaughters[i].size() << " prong daughters : ";
	  for(unsigned j=0;j<prongDaughters[i].size();j++)std::cout << prongDaughters[i][j].trackj << " ";
	  std::cout  << std::endl;
	}
	if(splitDaughters[i].size()>0){
	  std::cout  << " Track " << i << " has " << splitDaughters[i].size() << " split daughters : ";
	  for(unsigned j=0;j<splitDaughters[i].size();j++)std::cout << splitDaughters[i][j].trackj << " ";
	  std::cout  << std::endl;
	}
      }

      // Split track collection
      if(splitDaughters[i].size()==1){
	if(_debugPrinting>0)std::cout << " Reconstructed SPLIT " << std::endl;
	ReconstructedParticleImpl * part = new ReconstructedParticleImpl();
	VertexImpl * vtx = new VertexImpl();
	float vertex[3];
	float mom[3];
	int code = 211;
	for (int iC=0;iC<3;++iC) {
	  vertex[iC]   = splitDaughters[i][0].vtx[iC];
	  mom[iC] = splitDaughters[i][0].p[iC];
	}
	float distance = splitDaughters[i][0].distance;	
	vtx->setPosition( vertex );
	vtx->addParameter( distance );
	part->setMomentum( mom );
	part->setType( code );
	
	if(_debugPrinting>0){
	  std::cout << "   Code = " << code << "  Distance = " << distance << std::endl;
	  std::cout << "   Vertex = (" 
		    << vertex[0] << "," 
		    << vertex[1] << ","
		    << vertex[2] << ")" << std::endl;
	  std::cout << "   Momentum = ("
		    << mom[0] << ","
		    << mom[1] << ","
		    << mom[2] << ")" << std::endl;
	}
	float mass = splitDaughters[i][0].mass;
	part->setMass( mass );	
	vtx->setAssociatedParticle( part );
	part->setStartVertex( vtx );
	part->addTrack( tracks[splitDaughters[i][0].tracki] );
	part->addTrack( tracks[splitDaughters[i][0].trackj] );
	colSplitRecoPart->addElement( part );
	colSplitVertex->addElement( vtx );
      }

      // Kink collection
      if(splitDaughters[i].size()!=1){
	if(kinkDaughters[i].size()==1 && prongDaughters[i].size()<2){
	  if(_debugPrinting>0)std::cout << " Reconstructed KINKS " << std::endl;
	  ReconstructedParticleImpl * part = new ReconstructedParticleImpl();
	  VertexImpl * vtx = new VertexImpl();
	  float vertex[3];
	  float mom[3];
	  int code = kinkDaughters[i][0].pdgCode;
	  for (int iC=0;iC<3;++iC) {
	    vertex[iC]   = kinkDaughters[i][0].vtx[iC];
	    mom[iC] = kinkDaughters[i][0].p[iC];
	  }
	  float distance = kinkDaughters[i][0].distance;	
	  vtx->setPosition( vertex );
	  vtx->addParameter( distance );
	  part->setMomentum( mom );
	  part->setType( code );
	  if(_debugPrinting>0){
	    std::cout << "   Code = " << code << "  Distance = " << distance << std::endl;
	    std::cout << "   Vertex = (" 
		      << vertex[0] << "," 
		      << vertex[1] << ","
		      << vertex[2] << ")" << std::endl;
	    std::cout << "   Momentum = ("
		      << mom[0] << ","
		      << mom[1] << ","
		      << mom[2] << ")" << std::endl;
	  }
	  float mass = kinkDaughters[i][0].mass;
	  part->setMass( mass );	
	  vtx->setAssociatedParticle( part );
	  part->setStartVertex( vtx );
	  part->addTrack( tracks[kinkDaughters[i][0].tracki] );
	  part->addTrack( tracks[kinkDaughters[i][0].trackj] );
	  colKinkRecoPart->addElement( part );
	  colKinkVertex->addElement( vtx );
	  //   trackUsed[firstTrack] = 1;
	  //   trackUsed[secondTrack] = 1;
	} else if(prongDaughters[i].size() > 0 ) {
	  // Prong track collection
	  if(_debugPrinting>0)std::cout << " Reconstructed PRONG " << std::endl;
	  ReconstructedParticleImpl * part = new ReconstructedParticleImpl();
	  VertexImpl * vtx = new VertexImpl();
	  float vertex[3];
	  float mom[3];
	  int code = 211;
	  for (int iC=0;iC<3;++iC) {
	    vertex[iC]   = prongDaughters[i][0].vtx[iC];
	    mom[iC] = prongDaughters[i][0].p[iC];
	  }
	  float distance = prongDaughters[i][0].distance;	
	  vtx->setPosition( vertex );
	  vtx->addParameter( distance );
	  part->setMomentum( mom );
	  part->setType( code );
	  if(_debugPrinting>0){
	    std::cout << "   Code = " << code << "  Distance = " << distance << std::endl;
	    std::cout << "   Vertex = (" 
		      << vertex[0] << "," 
		      << vertex[1] << ","
		      << vertex[2] << ")" << std::endl;
	    std::cout << "   Momentum = ("
		      << mom[0] << ","
		      << mom[1] << ","
		      << mom[2] << ")" << std::endl;
	  }
	  float mass = prongDaughters[i][0].mass;
	  part->setMass( mass );	
	  vtx->setAssociatedParticle( part );
	  part->setStartVertex( vtx );
	  part->addTrack( tracks[prongDaughters[i][0].tracki] );
	  for(unsigned id=0;id<prongDaughters[i].size();id++)part->addTrack( tracks[prongDaughters[i][id].trackj] );
	  colProngRecoPart->addElement( part );
	  colProngVertex->addElement( vtx );
	}
	//   trackUsed[firstTrack] = 1;
	//   trackUsed[secondTrack] = 1;
      }
    }
   
  }

  evt->addCollection(colKinkRecoPart.release(), _kinkRecoPartColName);
  evt->addCollection(colKinkVertex.release(), _kinkVertexColName);
  evt->addCollection(colProngRecoPart.release(), _prongRecoPartColName);
  evt->addCollection(colProngVertex.release(), _prongVertexColName);
  evt->addCollection(colSplitRecoPart.release(), _splitRecoPartColName);
  evt->addCollection(colSplitVertex.release(), _splitVertexColName);

  for(unsigned int  itrack=0;itrack< tracks.size();++itrack){
    delete helixEnd[itrack];
    delete helixStart[itrack];
  }

  _nEvt++;
  if(_debugPrinting>0)std::cout << " Done " << std::endl;
  if(trackNavigator != NULL)  delete trackNavigator;
}


void KinkFinder::check( LCEvent * ) { }
  
void KinkFinder::end(){ } 

void KinkFinder::Sorting( TrackPairVec & trkPairVec ) {

  int sizeOfVector = int(trkPairVec.size());
  TrackPair *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++)
      {
	one = trkPairVec[j];
	two = trkPairVec[j+1];
	if( one->getDistance() > two->getDistance() )
	  {
	    Temp = trkPairVec[j];
	    trkPairVec[j] = trkPairVec[j+1];
	    trkPairVec[j+1] = Temp;
	  }
      }  

}

float KinkFinder::kinkMass(HelixClass* parent, HelixClass* daughter, float daughterMass, float neutralMass){

  // Calculate the invariant mass for a decaying charged particle
  
  float parentMom[4];
  parentMom[0] = parent->getMomentum()[0];
  parentMom[1] = parent->getMomentum()[1];
  parentMom[2] = parent->getMomentum()[2];
  float daughterMom[4];
  daughterMom[0] = daughter->getMomentum()[0];
  daughterMom[1] = daughter->getMomentum()[1];
  daughterMom[2] = daughter->getMomentum()[2];
  float  daughterE = sqrt(daughterMom[0]*daughterMom[0] + daughterMom[1]*daughterMom[1]+
                          daughterMom[2]*daughterMom[2] + daughterMass*daughterMass);

  float pnu[3];
  pnu[0] = parentMom[0]-daughterMom[0];
  pnu[1] = parentMom[1]-daughterMom[1];
  pnu[2] = parentMom[2]-daughterMom[2];
 
  float Enu = sqrt(pnu[0]*pnu[0]+pnu[1]*pnu[1]+pnu[2]*pnu[2]+neutralMass*neutralMass);
  float mx = sqrt((daughterE+Enu)*(daughterE+Enu)-parentMom[0]*parentMom[0]-
                       parentMom[1]*parentMom[1]-parentMom[2]*parentMom[2]);
  return mx;

}

