#include "VertexTracking.h"
#include <iostream>

#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <iostream>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

using namespace lcio ;
using namespace marlin ;

extern "C" {
  void tfithl_(int & NPT, double * XF,double * YF, float * RF, float *PF, double * WF,float * ZF,
	       float * WZF, int & IOPT, float * VV0,
	       float * EE0, float & CH2PH, float & CH2Z);
}

VertexTracking aVertexTracking ;

VertexTracking::VertexTracking() : Processor("VertexTracking") {

  std::vector<int> combinations;

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

  registerProcessorParameter("NumberOfLayers",
			     "Number of Layers",
			     _nLayers,
			     int(5));
  
  registerProcessorParameter("NDivisionsInPhi",
			     "Number of divisions in Phi",
			     _nDivisionsInPhi,
			     int(60));
  
   registerProcessorParameter("NDivisionsInTheta",
			     "Number of divisions in Theta",
			     _nDivisionsInTheta,
			     int(60));
  
   registerProcessorParameter("HitCollectionName",
			      "Hit Collection Name",
			      _VTXHitCollection,
			      std::string("VTXTrackerHits"));

   registerProcessorParameter("BField",
			      "Magnetic Field",
			      _bField,
			      float(4.0));
   
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
			      float(10.));
   
   registerProcessorParameter("Chi2WZQuartet",
			      "Chi2WZQuartet",
			      _chi2WZQuartet,
			      float(10.));
   
   registerProcessorParameter("Chi2WZSeptet",
			      "Chi2WZSeptet",
			      _chi2WZSeptet,
			      float(10.));
   
   registerProcessorParameter("ResolutionRPhi",
			      "ResolutionRPhi",
			      _resolutionRPhi,
			      float(0.006));

   registerProcessorParameter("ResolutionZ",
			      "ResolutionZ",
			      _resolutionZ,
			      float(0.006));

   registerProcessorParameter("MinimalDistance",
			      "Minimal Distance from Hit to Extrapolated Helix",
			      _distMinCut,
			      float(100.));
   

   registerProcessorParameter("TanLambdaCutForMerging",
			      "TanLambdaCutForMerging",
			      _tanlambdaCutForMerging,
			      float(0.01));
   

   registerProcessorParameter("PhiCutForMerging",
			      "PhiCutForMerging",
			      _phiCutForMerging,
			      float(0.01));
   
  
   registerProcessorParameter("MinDistCutAttach",
			      "MinDistCutAttach",
			      _minDistCutAttach,
			      float(1.0));
   
    registerProcessorParameter("MinLayerToAttach",
			      "MinLayerToAttach",
			      _minimalLayerToAttach,
			      int(0));
   
  
   

}



void VertexTracking::init() { 

    _nRun = -1 ;
    _nEvt = 0 ;
    printParameters() ;


}


void VertexTracking::processRunHeader( LCRunHeader* run) { 

    _nRun++ ;
    _nEvt = 0;

    PI = acos(double(-1.0));
    TWOPI = double(2.0)*PI;
    PIOVER2 = double(0.5)*PI;
    _dPhi = TWOPI/double(_nDivisionsInPhi);
    _dTheta = double(2.0)/double(_nDivisionsInTheta);


} 

void VertexTracking::processEvent( LCEvent * evt ) { 

  _nEvt++;

  //  std::cout << std::endl;
  //  std::cout << "Event number = " << _nEvt << std::endl;

  int success = Initialise( evt );

  if (success == 1) {

  std::cout << "Event = " << _nEvt << std::endl;
    

    for (int iPhi=0; iPhi<_nDivisionsInPhi; ++iPhi) 
      for (int iTheta=0; iTheta<_nDivisionsInTheta;++iTheta)
	ProcessOneSector(iPhi,iTheta);      

    std::cout << "End of Processing sectors" << std::endl;

    Sorting( _tracks5Hits);
    Sorting( _tracks4Hits);
    Sorting( _tracks3Hits);

    std::cout << "End of Sorting " << std::endl;


    int nTrk = int(_tracks5Hits.size());
    for (int iTrk=0; iTrk<nTrk;++iTrk) {
      TrackExtended * trackAR = _tracks5Hits[iTrk];
      CreateVTXTrack( trackAR );
    }

    std::cout << "End of creating 5 hits tracks " << std::endl;
    

    nTrk = int(_tracks4Hits.size());
    for (int iTrk=0; iTrk<nTrk;++iTrk) {
      TrackExtended * trackAR = _tracks4Hits[iTrk];
      CreateVTXTrack( trackAR );
    }

    std::cout << "End of creating 4 hits tracks " << std::endl;

    nTrk = int(_tracks3Hits.size());
    for (int iTrk=0; iTrk<nTrk;++iTrk) {
      TrackExtended * trackAR = _tracks3Hits[iTrk];
      CreateVTXTrack( trackAR );
    }

    std::cout << "End of creating 3 hits tracks " << std::endl;

    AttachRemainingHits();
    
    std::cout << "End of picking up remaining hits " << std::endl;

    LCCollectionVec * trkCol = new LCCollectionVec(LCIO::TRACK);
    nTrk = int(_trackImplVec.size());
    for (int iTrk=0; iTrk<nTrk;++iTrk) {
      TrackExtended * trackAR = _trackImplVec[iTrk];
      TrackerHitExtendedVec hitVec = trackAR->getTrackerHitExtendedVec();
      TrackImpl * trackImpl = new TrackImpl();
      trackImpl->setOmega(trackAR->getOmega());
      trackImpl->setTanLambda(trackAR->getTanLambda());
      trackImpl->setPhi(trackAR->getPhi());
      trackImpl->setD0(trackAR->getD0());
      trackImpl->setZ0(trackAR->getZ0());
      trackImpl->setChi2(trackAR->getChi2());
      int nHits = int(hitVec.size());
      for (int iHit=0;iHit<nHits;++iHit) {
	TrackerHit * trkHit = hitVec[iHit]->getTrackerHit();
	trackImpl->addHit(trkHit);
      }

      trkCol->addElement(trackImpl);
    }    
    
    evt->addCollection(trkCol,"VTXTracks");     

  }

  CleanUp();



}


void VertexTracking::CleanUp() {
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


}

int VertexTracking::Initialise(LCEvent * evt) {

  int success = 1;

  _sectors.clear();
  _sectors.resize(_nLayers*_nDivisionsInPhi*_nDivisionsInTheta);
  _tracks5Hits.clear();
  _tracks4Hits.clear();
  _tracks3Hits.clear();
  _trackImplVec.clear();

  for (int i=0; i<_nLayers*_nDivisionsInPhi*_nDivisionsInTheta; ++i) {
    TrackerHitExtendedVec hitVec;
    hitVec.clear();
    _sectors.push_back(hitVec);
  }
    
  try {
    LCCollection * hitCollection = evt->getCollection(_VTXHitCollection.c_str());
    int nelem = hitCollection->getNumberOfElements();
    //    std::cout << "Number of hits = " << nelem << std::endl;
    for (int ielem=0; ielem<nelem; ++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
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
      int iCode = layer + _nLayers*iPhi + _nLayers*_nDivisionsInPhi*iTheta;      
      _sectors[iCode].push_back( hitExt );
    }
  }
  catch(DataNotAvailableException &e) {

    success = 0;
  }

  return success;

}

void VertexTracking::check(LCEvent * evt) { };

void VertexTracking::end() {

}


void VertexTracking::ProcessOneSector(int iPhi, int iTheta) {
  
  int iPhi_Up    = iPhi + 1;
  int iPhi_Low   = iPhi - 1;
  int iTheta_Up  = iTheta + 1; 
  int iTheta_Low = iTheta - 1;
  if (iTheta_Low < 0) iTheta_Low = 0;
  if (iTheta_Up  >= _nDivisionsInTheta) iTheta_Up = _nDivisionsInTheta-1;

  int nComb = int( _Combinations.size() / 3 ); // number of triplet combinations

  //  std::cout << iPhi << " " << iTheta << std::endl;

  int iNC = 0;
  for (int iComb=0; iComb < nComb; ++iComb) { // loop over triplets
    int nLR[3];
    for (int iS=0; iS<3; ++iS) {
      nLR[iS] = _Combinations[iNC];
      iNC++;
    }
    
    //    std::cout << iComb << " " << nLR[0] << " " << nLR[1] << " " << nLR[2] << std::endl;

    int iCode = nLR[0] + _nLayers*iPhi +  _nLayers*_nDivisionsInPhi*iTheta;
    TrackerHitExtendedVec hitVecOuter =  _sectors[iCode]; 
    int nHitsOuter = int(hitVecOuter.size());
    for (int iOuter=0; iOuter<nHitsOuter; ++iOuter) { // loop over outer hits
      TrackerHitExtended * outerHit = hitVecOuter[iOuter];
      for (int ipMiddle=iPhi_Low; ipMiddle<iPhi_Up+1;ipMiddle++) { // loop over phi in the Middle
	for (int itMiddle=iTheta_Low; itMiddle<iTheta_Up+1;itMiddle++) { // loop over theta in the Middle 
	  int iPhiMiddle = ipMiddle;
	  if (ipMiddle < 0) iPhiMiddle = _nDivisionsInPhi-1;
	  if (ipMiddle >= _nDivisionsInPhi) iPhiMiddle = ipMiddle - _nDivisionsInPhi;
	  iCode = nLR[1] + _nLayers*iPhiMiddle +  _nLayers*_nDivisionsInPhi*itMiddle;
	  TrackerHitExtendedVec hitVecMiddle = _sectors[iCode];
	  int nHitsMiddle = int(hitVecMiddle.size());
	  for (int iMiddle=0;iMiddle<nHitsMiddle;iMiddle++) { // loop over hits in the Middle sector
	    TrackerHitExtended * middleHit = hitVecMiddle[iMiddle];
	    for (int ipInner=iPhi_Low; ipInner<iPhi_Up+1;ipInner++) { // loop over phi in the Inner
	      for (int itInner=iTheta_Low; itInner<iTheta_Up+1;itInner++) { // loop over theta in the Inner 
		int iPhiInner = ipInner;
		if (ipInner < 0) iPhiInner = _nDivisionsInPhi-1;
		if (ipInner >= _nDivisionsInPhi) iPhiInner = ipInner - _nDivisionsInPhi;
		iCode = nLR[2] + _nLayers*iPhiInner +  _nLayers*_nDivisionsInPhi*itInner;
		TrackerHitExtendedVec hitVecInner = _sectors[iCode];
		int nHitsInner = int(hitVecInner.size());
		for (int iInner=0;iInner<nHitsInner;iInner++) { // loop over hits in the Inner sector
		  TrackerHitExtended * innerHit = hitVecInner[iInner];
		  HelixClass helix;
		  //		  std::cout << "Trace step " << std::endl;
		  TrackExtended * trackAR = TestTriplet(outerHit,middleHit,innerHit,helix);
		  if ( trackAR != NULL ) {
		    //		    std::cout << "Track found : " << nLR[0] << " " << nLR[1] << " " << nLR[2] << std::endl;
		    int nH = BuildTrack(outerHit,middleHit,innerHit,helix,nLR[2],iPhi,iTheta,trackAR);
		    if (nH == 3) 
		      _tracks3Hits.push_back(trackAR);
		    if (nH == 4)
		      _tracks4Hits.push_back(trackAR);
		    if (nH >= 5)
		      _tracks5Hits.push_back(trackAR);
		    //		    std::cout << "OK " << nH << std::endl; 
		  }
		} // endloop over hits in the Inner sector
	      } // endloop over theta in the Inner
	    } // endloop over phi in the Inner	    
	  } // endloop over hits in the Middle sector
	} // endloop over theta in the Middle
      } // endloop over phi in the Middle
    } // endloop over outer hits
  } // endloop over triplets

  //  std::cout << "Heya " << std::endl;

}

TrackExtended * VertexTracking::TestTriplet(TrackerHitExtended * outerHit, 
					    TrackerHitExtended * middleHit,
					    TrackerHitExtended * innerHit,
					    HelixClass & helix) {

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


  if (nTrackOuter > 0 && nTrackInner > 0 && nTrackMiddle) {
    for (int iMiddle=0; iMiddle<nTrackMiddle ; iMiddle++) {
      for (int iOuter=0; iOuter<nTrackOuter; iOuter++) {
	for (int iInner=0; iInner<nTrackInner; iInner++) {
	  if (trackOuterVec[iOuter] == trackInnerVec[iInner] && trackInnerVec[iInner] == trackMiddleVec[iMiddle]) 
	    return trackAR;     
	}      
      }    
    }
  }








  double * xh = new double[3];
  double * yh = new double[3];
  float * zh = new float[3];
  double * wrh = new double[3];
  float * wzh = new float[3];
  float * rh = new float[3];
  float * ph = new float[3];
  float par[5];
  float epar[15];

  xh[0] = outerHit->getTrackerHit()->getPosition()[0];
  yh[0] = outerHit->getTrackerHit()->getPosition()[1];
  zh[0] = float(outerHit->getTrackerHit()->getPosition()[2]);
  
  xh[1] = middleHit->getTrackerHit()->getPosition()[0];
  yh[1] = middleHit->getTrackerHit()->getPosition()[1];
  zh[1] = float(middleHit->getTrackerHit()->getPosition()[2]);
  
  xh[2] = innerHit->getTrackerHit()->getPosition()[0];
  yh[2] = innerHit->getTrackerHit()->getPosition()[1];
  zh[2] = float(innerHit->getTrackerHit()->getPosition()[2]);
  
  for (int ih=0; ih<3; ih++) {
    rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
    ph[ih] = atan2(yh[ih],xh[ih]);
    if (ph[ih] < 0.) 
      ph[ih] = 2.0*acos(-1.0) + ph[ih]; 
    wrh[ih] = double(1.0/(_resolutionRPhi*_resolutionRPhi));
    wzh[ih] = 1.0/(_resolutionZ*_resolutionZ);
  }

  int NPT = 3;
  int IOPT = 2;
  float chi2RPhi;
  float chi2Z;
  tfithl_(NPT, xh, yh, rh, ph, wrh, zh,
	  wzh, IOPT, par, epar, chi2RPhi, chi2Z);


  delete[] xh;
  delete[] yh;
  delete[] zh;
  delete[] wrh;
  delete[] wzh;
  delete[] rh;
  delete[] ph;
 

  float Chi2 = chi2RPhi/_chi2WRPhiTriplet+chi2Z/_chi2WZTriplet;
  if ( Chi2 > 2.0 )
    return trackAR;

  float omega = par[0];
  float tanlambda = par[1];
  float phi0 = par[2];
  float d0 = par[3];
  float z0 = par[4];

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
  trackAR->setChi2(chi2Z+chi2RPhi);
  trackAR->setCovMatrix(epar);
  return trackAR;

}

int VertexTracking::BuildTrack(TrackerHitExtended * outerHit, 
			       TrackerHitExtended * middleHit,
			       TrackerHitExtended * innerHit,
			       HelixClass & helix,
			       int innerLayer,
			       int iPhi,
			       int iTheta, 
			       TrackExtended * trackAR) {

  int iPhi_Up    = iPhi + 1;
  int iPhi_Low   = iPhi - 1;
  int iTheta_Up  = iTheta + 1; 
  int iTheta_Low = iTheta - 1;
  if (iTheta_Low < 0) iTheta_Low = 0;
  if (iTheta_Up  >= _nDivisionsInTheta) iTheta_Up = _nDivisionsInTheta-1;

  for (int layer = innerLayer-1; layer>=0; layer--) { // loop over remaining layers
    //    std::cout << "checking layer " << layer << std::endl;
    float distMin = 1.0e+20;
    TrackerHitExtended * assignedhit = NULL;
    for (int ipInner=iPhi_Low; ipInner<iPhi_Up+1;ipInner++) { // loop over phi in the Inner region
      for (int itInner=iTheta_Low; itInner<iTheta_Up+1;itInner++) { // loop over theta in the Inner region 
	int iPhiInner = ipInner;
	if (ipInner < 0) iPhiInner = _nDivisionsInPhi-1;
	if (ipInner >= _nDivisionsInPhi) iPhiInner = ipInner - _nDivisionsInPhi;
	int iCode = layer + _nLayers*iPhiInner +  _nLayers*_nDivisionsInPhi*itInner;
	TrackerHitExtendedVec hitVecInner = _sectors[iCode];
	int nHitsInner = int(hitVecInner.size());
	for (int iInner=0;iInner<nHitsInner;iInner++) { // loop over hits in the Inner sector
	  TrackerHitExtended * currentHit = hitVecInner[iInner];
	  float pos[3];
	  //	  float distance[3];
	  float ref[3];
	  float point[6];
	  for (int i=0; i<3; ++i) {
	    pos[i] = float(currentHit->getTrackerHit()->getPosition()[i]);
	    ref[i] = helix.getReferencePoint()[i];
	  }
	  float radius = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
	  //	  helix.getDistanceToPoint(pos,distance);
	  
	  float time = helix.getPointOnCircle(radius,ref,point);
	  if (time < 1.0e+10) {
	    float distRPhi1 = sqrt((point[0]-pos[0])*(point[0]-pos[0])+
				   (point[1]-pos[1])*(point[1]-pos[1]));
	    float distRPhi2 = sqrt((point[3]-pos[0])*(point[3]-pos[0])+
				   (point[4]-pos[1])*(point[4]-pos[1]));
	    float distZ1 = fabs(point[2]-pos[2]);
	    float distZ2 = fabs(point[5]-pos[2]);	    
	    float dist1 = distRPhi1/_resolutionRPhi + distZ1/_resolutionZ;
	    float dist2 = distRPhi2/_resolutionRPhi + distZ2/_resolutionZ;
	    float dist = fmin(dist1,dist2);
	    _distRPhi = fmin(distRPhi1,distRPhi2);
	    _distZ = fmin(distZ1,distZ2);
	    if (dist < distMin) {
	      distMin = dist;
	      assignedhit = currentHit;
	    }
	  }
	} // endloop over hits in the Inner sector
      } // endloop over theta in the Inner region 
    } // endloop over phi in the Inner region
    if (distMin < _distMinCut) {
      TrackerHitExtendedVec hvec = trackAR->getTrackerHitExtendedVec();
      int  nHits = int(hvec.size());
      double * xh = new double[nHits+1];
      double * yh = new double[nHits+1];
      float * zh = new float[nHits+1];
      double * wrh = new double[nHits+1];
      float * wzh = new float[nHits+1];
      float * rh = new float[nHits+1];
      float * ph = new float[nHits+1];
      float par[5];
      float epar[15];
      for (int ih=0;ih<nHits;++ih) {
	TrackerHit * trkHit = hvec[ih]->getTrackerHit();
	xh[ih] = trkHit->getPosition()[0];
	yh[ih] = trkHit->getPosition()[1];
	zh[ih] = float(trkHit->getPosition()[2]);
	wrh[ih] = double(1.0/(_resolutionRPhi*_resolutionRPhi));
	wzh[ih] = 1.0/(_resolutionZ*_resolutionZ);
	rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
	ph[ih] = float(atan2(yh[ih],xh[ih]));
	if (ph[ih] < 0.) 
	  ph[ih] = 2.0*acos(-1.0) + ph[ih]; 
      }      
      TrackerHit * assignedTrkHit = assignedhit->getTrackerHit();
      xh[nHits] = assignedTrkHit->getPosition()[0];
      yh[nHits] = assignedTrkHit->getPosition()[1];
      zh[nHits] = float(assignedTrkHit->getPosition()[2]);
      rh[nHits] = float(sqrt(xh[nHits]*xh[nHits]+yh[nHits]*yh[nHits]));
      ph[nHits] = float(atan2(yh[nHits],xh[nHits]));
      if (ph[nHits] < 0.) 
	ph[nHits] = 2.0*acos(-1.0) + ph[nHits]; 
      wrh[nHits] = double(1.0/(_resolutionRPhi*_resolutionRPhi));
      wzh[nHits] = 1.0/(_resolutionZ*_resolutionZ);
      int NPT = nHits + 1;
      int IOPT = 2;
      float chi2RPhi;
      float chi2Z;
      tfithl_(NPT, xh, yh, rh, ph, wrh, zh,
	      wzh, IOPT, par, epar, chi2RPhi, chi2Z);

      delete[] xh;
      delete[] yh;
      delete[] zh;
      delete[] wrh;
      delete[] wzh;
      delete[] rh;
      delete[] ph;

      bool validCombination = 0;
      float Chi2;
      if ((nHits+1) == 4) {
	Chi2 = chi2RPhi/_chi2WRPhiQuartet+chi2Z/_chi2WZQuartet;
	validCombination = Chi2 < 2.0;
      }      
      
      if ((nHits+1) >= 5) {
	Chi2 = chi2RPhi/_chi2WRPhiSeptet+chi2Z/_chi2WZSeptet;
	validCombination = Chi2 < 2.0;                 
      }

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
	trackAR->setChi2(chi2Z+chi2RPhi);
	trackAR->setCovMatrix(epar);
      }
     
    }
  } // endloop over remaining layers

  TrackerHitExtendedVec hvec = trackAR->getTrackerHitExtendedVec();  
  int nTotalHits = int(hvec.size());
  return nTotalHits;

}


void VertexTracking::Sorting(TrackExtendedVec & trackVec) {

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

void VertexTracking::CreateVTXTrack(TrackExtended * trackAR ) {

  TrackerHitExtendedVec hitVec = trackAR->getTrackerHitExtendedVec();
  int nHits = int(hitVec.size());

  for (int i=0; i<nHits; ++i) {
    TrackExtendedVec trackVec = hitVec[i]->getTrackExtendedVec();
    if (trackVec.size() != 0) 
      return ;
  }

  // Find splitted tracks

  int found = 0;
  int nTrk = int(_trackImplVec.size());
  for (int i=0; i<nTrk; ++i) {
    TrackExtended * trackOld = _trackImplVec[i];
    TrackerHitExtendedVec hitVecOld = trackOld->getTrackerHitExtendedVec();
    
//     gsl_matrix* A = gsl_matrix_alloc(5,5);
//     gsl_vector* z = gsl_vector_alloc(5);
//     gsl_vector* e = gsl_vector_alloc(5);

//     FloatVec covMatOld;
//     float covMat[15];
//     float XOLD[5];
//     float X[5];

//     XOLD[0] = trackOld->getOmega();
//     XOLD[1] = trackOld->getTanLambda();
//     XOLD[2] = trackOld->getPhi();
//     XOLD[3] = trackOld->getD0();
//     XOLD[4] = trackOld->getZ0();
    
//     X[0] = trackAR->getOmega();
//     X[1] = trackAR->getTanLambda();
//     X[2] = trackAR->getPhi();
//     X[3] = trackAR->getD0();
//     X[4] = trackAR->getZ0();
    
//     covMatOld = trackOld->getCovMatrix();
//     for (int i=0;i<15;++i)
//       covMat[i] = trackAR->getCovMatrix()[i];

//   // initialise matrix and vectors
  
//     for (int i=0;i<5;++i) {
//       float B = 0;      
//       for (int j=0;j<5;++j) {
// 	int code = decode(i,j);
// 	B = B + covMatOld[code]*XOLD[j] + covMat[code]*X[j];
// 	float Element = covMatOld[code]+covMat[code];
// 	gsl_matrix_set(A,i,j,Element);
//       }
//       gsl_vector_set(z,i,B);
//     }

//     gsl_linalg_HH_solve(A,z,e);

//     float Chi2 = 0;
//     for (int i=0;i<5;++i) {
//       for (int j=0;j<5;++j) {
// 	float XI =  gsl_vector_get(e,i);
// 	float XJ =  gsl_vector_get(e,j);
// 	int code = decode(i,j);
// 	Chi2 += (XI-XOLD[i])*covMatOld[code]*(XJ-XOLD[j]);
// 	Chi2 += (XI-X[i])*covMat[code]*(XJ-X[j]);
//       }
//     }

//     float omega = gsl_vector_get(e,0);
//     float tanlambda = gsl_vector_get(e,1);
//     float phi0 = gsl_vector_get(e,2);
//     float d0 = gsl_vector_get(e,3);
//     float z0 = gsl_vector_get(e,4);
    
//     gsl_matrix_free(A);
//     gsl_vector_free(z);
//     gsl_vector_free(e);    

    if ((fabs(trackAR->getTanLambda()-trackOld->getTanLambda())<_tanlambdaCutForMerging && 
	fabs(trackAR->getPhi()-trackOld->getPhi())<_phiCutForMerging)) {
      
      int nHitsOld = int(hitVecOld.size());
      int nTotHits = nHits + nHitsOld;
      double * xh = new double[nTotHits];
      double * yh = new double[nTotHits];
      float * zh = new float[nTotHits];
      double * wrh = new double[nTotHits];
      float * wzh = new float[nTotHits];
      float * rh = new float[nTotHits];
      float * ph = new float[nTotHits];
      float par[5];
      float epar[15];
      for (int ih=0;ih<nHits;++ih) {
	TrackerHit * trkHit = hitVec[ih]->getTrackerHit();
	if (int(hitVec[ih]->getTrackExtendedVec().size()) != 0)
	  std::cout << "WARNING : HIT POINTS TO TRACK " << std::endl;
	xh[ih] = trkHit->getPosition()[0];
	yh[ih] = trkHit->getPosition()[1];
	zh[ih] = float(trkHit->getPosition()[2]);
	wrh[ih] = double(1.0/(_resolutionRPhi*_resolutionRPhi));
	wzh[ih] = 1.0/(_resolutionZ*_resolutionZ);
	rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
	ph[ih] = float(atan2(yh[ih],xh[ih]));
      }      
      for (int ih=0;ih<nHitsOld;++ih) {
	TrackerHit * trkHit = hitVecOld[ih]->getTrackerHit();
	xh[ih+nHits] = trkHit->getPosition()[0];
	yh[ih+nHits] = trkHit->getPosition()[1];
	zh[ih+nHits] = float(trkHit->getPosition()[2]);
	wrh[ih+nHits] = double(1.0/(_resolutionRPhi*_resolutionRPhi));
	wzh[ih+nHits] = 1.0/(_resolutionZ*_resolutionZ);
	rh[ih+nHits] = float(sqrt(xh[ih+nHits]*xh[ih+nHits]+yh[ih+nHits]*yh[ih+nHits]));
	ph[ih+nHits] = float(atan2(yh[ih+nHits],xh[ih+nHits]));
      }
      int NPT = nTotHits;
      int IOPT = 3;
      float chi2RPhi;
      float chi2Z;
      tfithl_(NPT, xh, yh, rh, ph, wrh, zh,
	      wzh, IOPT, par, epar, chi2RPhi, chi2Z);


      float omega = par[0];
      float tanlambda = par[1];
      float phi0 = par[2];
      float d0 = par[3];
      float z0 = par[4];
     
      float chi2Min = chi2RPhi/_chi2WRPhiSeptet+chi2Z/_chi2WZSeptet;
      float chi2Full = chi2Min;
      float chi2MinRPhi = chi2RPhi;
      float chi2MinZ = chi2Z;
      
      int iBad = -1;
      if (chi2Full < 2.0) {
	found = 1;
      }
      else {
	for (int i=0; i<nTotHits; ++i) {
	  for (int j=0;j<nTotHits;++j) {
	    if (i == j) {
	      wrh[j] = 0.0;
	      wzh[j] = 0.0;
	    } 
	    else {
	      wrh[j] = double(1.0/(_resolutionRPhi*_resolutionRPhi));
	      wzh[j] = 1.0/(_resolutionZ*_resolutionZ);
	    }
	  }
	  tfithl_(NPT, xh, yh, rh, ph, wrh, zh,
		  wzh, IOPT, par, epar, chi2RPhi, chi2Z);
	  float chi2Cur = chi2RPhi/_chi2WRPhiSeptet+chi2Z/_chi2WZSeptet;
	  if (chi2Cur < chi2Min) {
	    chi2Min = chi2Cur;
	    chi2MinRPhi = chi2RPhi;
	    chi2MinZ = chi2Z;
	    omega = par[0];
	    tanlambda = par[1];
	    phi0 = par[2];
	    d0 = par[3];
	    z0 = par[4];
	    iBad = i;
	  }
	}
	if (chi2Min < 2.0) {
	  found = 1;
	}
      }

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
	trackOld->setChi2(chi2MinRPhi+chi2MinZ);	
      }

      delete[] xh;
      delete[] yh;
      delete[] zh;
      delete[] wrh;
      delete[] wzh;
      delete[] rh;
      delete[] ph;

    }
    if (found == 1)
      break;

  }
    
  if (found == 0 ) {
    _trackImplVec.push_back(trackAR);
    for (int i=0;i<nHits;++i) {
      TrackerHitExtended * hit = hitVec[i];
      hit->addTrackExtended( trackAR );
    }
  }
    

}

int VertexTracking::decode(int i, int j) {
  
  int imin = i; 
  int imax = j; 
  if (i>j) {
    imin = j;
    imax = i;
  }
  int output = 0;

  if (imin == 0 && imax == 0)
    output = 0;
  else if (imin == 0 && imax == 1)
    output = 1;
  else if (imin == 0 && imax == 2)
    output = 3;
  else if (imin == 0 && imax == 3)
    output = 6;
  else if (imin == 0 && imax == 4)
    output = 10;
  else if (imin == 1 && imax == 1)
    output = 2;
  else if (imin == 1 && imax == 2)
    output = 4;
  else if (imin == 1 && imax == 3)
    output = 7;
  else if (imin == 1 && imax == 4)
    output = 11;
  else if (imin == 2 && imax == 2)
    output = 5;
  else if (imin == 2 && imax == 3)
    output = 8;
  else if (imin == 2 && imax == 4)
    output = 12;
  else if (imin == 3 && imax == 3)
    output = 9;
  else if (imin == 3 && imax == 4)
    output = 13;
  else if (imin == 4 && imax == 4)
    output = 14;

//   if (imin == 0 && imax == 0)
//     output = 0;
//   else if (imin == 0 && imax == 1)
//     output = 1;
//   else if (imin == 0 && imax == 2)
//     output = 2;
//   else if (imin == 0 && imax == 3)
//     output = 3;
//   else if (imin == 0 && imax == 4)
//     output = 4;
//   else if (imin == 1 && imax == 1)
//     output = 5;
//   else if (imin == 1 && imax == 2)
//     output = 6;
//   else if (imin == 1 && imax == 3)
//     output = 7;
//   else if (imin == 1 && imax == 4)
//     output = 8;
//   else if (imin == 2 && imax == 2)
//     output = 9;
//   else if (imin == 2 && imax == 3)
//     output = 10;
//   else if (imin == 2 && imax == 4)
//     output = 11;
//   else if (imin == 3 && imax == 3)
//     output = 12;
//   else if (imin == 3 && imax == 4)
//     output = 13;
//   else if (imin == 4 && imax == 4)
//     output = 14;


  return output;

}


void VertexTracking::AttachRemainingHits() {

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
	helix.getDistanceToPoint(pos,distance);
	if (distance[2] < minDist) {
	  minDist = distance[2];
	  trackToAttach = trackAR;
	}
      }
      if (minDist < _minDistCutAttach && trackToAttach != NULL) {
	TrackerHitExtendedVec hitVec = trackToAttach->getTrackerHitExtendedVec();
	int nHits = int(hitVec.size());
	double * xh = new double[nHits+1];
	double * yh = new double[nHits+1];
	float * zh = new float[nHits+1];
	double * wrh = new double[nHits+1];
	float * wzh = new float[nHits+1];
	float * rh = new float[nHits+1];
	float * ph = new float[nHits+1];
	float par[5];
	float epar[15];
	
	for (int i=0; i<nHits; ++i) {
	  TrackerHit * trkHit = hitVec[i]->getTrackerHit();
	  xh[i] = double(trkHit->getPosition()[0]);
	  yh[i] = double(trkHit->getPosition()[1]);
	  zh[i] = float(trkHit->getPosition()[2]);
	  ph[i] = float(atan2(yh[i],xh[i]));
	  rh[i] = float(sqrt(xh[i]*xh[i]+yh[i]*yh[i]));
	  wrh[i] = double(1.0/(_resolutionRPhi*_resolutionRPhi));
	  wzh[i] = 1.0/(_resolutionZ*_resolutionZ);	  
	}
	TrackerHit * trkHit = hit->getTrackerHit();
	xh[nHits] = double(trkHit->getPosition()[0]);
	yh[nHits] = double(trkHit->getPosition()[1]);
	zh[nHits] = float(trkHit->getPosition()[2]);
	ph[nHits] = float(atan2(yh[nHits],xh[nHits]));
	rh[nHits] = float(sqrt(xh[nHits]*xh[nHits]+yh[nHits]*yh[nHits]));
	wrh[nHits] = double(1.0/(_resolutionRPhi*_resolutionRPhi));
	wzh[nHits] = 1.0/(_resolutionZ*_resolutionZ);	

	int NPT = nHits + 1;
	int IOPT = 3;
	float chi2RPhi;
	float chi2Z;
	tfithl_(NPT, xh, yh, rh, ph, wrh, zh,
		wzh, IOPT, par, epar, chi2RPhi, chi2Z);


	float omega = par[0];
	float tanlambda = par[1];
	float phi0 = par[2];
	float d0 = par[3];
	float z0 = par[4];
	

	float chi2;
	if (NPT == 3) {
	  chi2 = chi2RPhi/_chi2WRPhiTriplet+chi2Z/_chi2WZTriplet;
	}
	if (NPT == 4) {
	  chi2 = chi2RPhi/_chi2WRPhiQuartet+chi2Z/_chi2WZQuartet;
	}
	if (NPT > 4) {
	  chi2 = chi2RPhi/_chi2WRPhiSeptet+chi2Z/_chi2WZSeptet;
	}

	if (chi2 < 2.0) {
	  trackToAttach->addTrackerHitExtended(hit);
	  hit->addTrackExtended( trackToAttach );
	  trackToAttach->setChi2( chi2RPhi + chi2Z );
	  trackToAttach->setOmega( omega );
	  trackToAttach->setTanLambda( tanlambda );
	  trackToAttach->setD0( d0 );
	  trackToAttach->setZ0( z0 );
	  trackToAttach->setPhi( phi0 );
	}	
	delete[] xh;
	delete[] yh;
	delete[] zh;
	delete[] wrh;
	delete[] wzh;
	delete[] rh;
	delete[] ph;
      }      
    }
  }

}
