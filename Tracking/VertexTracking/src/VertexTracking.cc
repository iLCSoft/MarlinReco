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
#include <marlin/Global.h>

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

  std::vector<int> combinationsFTD;

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

  combinationsFTD.push_back(2);
  combinationsFTD.push_back(1);
  combinationsFTD.push_back(0);




  registerProcessorParameter("LayerCombinationsFTD",
			     "Combinations of Hits in FTD",
			     _CombinationsFTD,
			     combinationsFTD);


  registerProcessorParameter("NumberOfLayers",
			     "Number of Layers",
			     _nLayers,
			     int(5));
  
  registerProcessorParameter("NumberOfFTDLayers",
			     "Number of FTD Layers",
			     _nLayersFTD,
			     int(7));
  
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
  
  registerInputCollection( LCIO::TRACKERHIT,
			   "VXDHitCollectionName",
			   "VXD Hit Collection Name",
			   _VTXHitCollection,
			   std::string("VTXTrackerHits"));


  registerInputCollection( LCIO::TRACKERHIT,
			   "FTDHitCollectionName",
			   "FTD Hit Collection Name",
			   _FTDHitCollection,
			   std::string("FTDTrackerHits"));  
  
  registerOutputCollection( LCIO::TRACK,
			    "VTXTrackCollection",
			    "VTX Tracks",
			    _colVTXTracks,
			    std::string("VTXTracks"));  
   


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
   
  registerProcessorParameter("ResolutionRPhi",
			     "ResolutionRPhi",
			     _resolutionRPhi,
			     float(0.005));

  registerProcessorParameter("ResolutionZ",
			     "ResolutionZ",
			     _resolutionZ,
			      float(0.005));

  registerProcessorParameter("ResolutionRPhiFTD",
			     "ResolutionRPhiFTD",
			     _resolutionRPhiFTD,
			      float(0.01));

  registerProcessorParameter("ResolutionZFTD",
			     "ResolutionZFTD",
			     _resolutionZFTD,
			     float(0.10));

  registerProcessorParameter("TanLambdaCutForMerging",
			     "TanLambdaCutForMerging",
			     _tanlambdaCutForMerging,
			     float(0.05));
   

  registerProcessorParameter("PhiCutForMerging",
			     "PhiCutForMerging",
			     _phiCutForMerging,
			     float(0.05));
   
  
  registerProcessorParameter("MinDistCutAttach",
			     "MinDistCutAttach",
			     _minDistCutAttach,
			     float(1.0));
   
  registerProcessorParameter("MinLayerToAttach",
			     "MinLayerToAttach",
			     _minimalLayerToAttach,
			     int(0));
   
  registerProcessorParameter("CutOnD0",
			     "cut on D0 for tracks",
			     _cutOnD0,
			     float(50.0));
   
  registerProcessorParameter("CutOnZ0",
			     "cut on Z0 for tracks",
			     _cutOnZ0,
			     float(50.0));
   
    
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


  std::vector<float> zlayer;
  zlayer.push_back(200);
  zlayer.push_back(320);
  zlayer.push_back(440);
  zlayer.push_back(550);
  zlayer.push_back(800);
  zlayer.push_back(1050);
  zlayer.push_back(1300);
  
  registerProcessorParameter("ZlayerFTD",
			     "Z Layer for FTD",
			     _zLayerFTD,
			     zlayer);


}



void VertexTracking::init() { 

    _nRun = -1 ;
    _nEvt = 0 ;
    printParameters() ;


}


void VertexTracking::processRunHeader( LCRunHeader* run) { 

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

void VertexTracking::processEvent( LCEvent * evt ) { 


  // Clearing all working dynamical arrays (vectors)
  _tracks5Hits.clear();
  _tracks4Hits.clear();
  _tracks3Hits.clear();
  _trackImplVec.clear();

  int success = Initialise( evt );
  int successFTD = InitialiseFTD( evt );
  std::cout << "Event = " << _nEvt << " Run = " << _nRun << std::endl;

  if (success == 1) {
    for (int iPhi=0; iPhi<_nDivisionsInPhi; ++iPhi) 
      for (int iTheta=0; iTheta<_nDivisionsInTheta;++iTheta)
	ProcessOneSector(iPhi,iTheta); // Process one VXD sector     
    std::cout << "End of Processing VXD sectors" << std::endl;
  }

  if (successFTD == 1) {
    TrackingInFTD(); // Perform tracking in the FTD
    std::cout << "End of Processing FTD sectors" << std::endl;
  }

  if (success == 1 || successFTD == 1) {

    Sorting( _tracks5Hits);
    Sorting( _tracks4Hits);
    Sorting( _tracks3Hits);
    std::cout << "End of Sorting " << std::endl;

    int nTrk = int(_tracks5Hits.size());
    for (int iTrk=0; iTrk<nTrk;++iTrk) {
      TrackExtended * trackAR = _tracks5Hits[iTrk];
      CreateTrack( trackAR );
    }
    std::cout << "End of creating 5 hits tracks " << std::endl;
    

    nTrk = int(_tracks4Hits.size());
    for (int iTrk=0; iTrk<nTrk;++iTrk) {
      TrackExtended * trackAR = _tracks4Hits[iTrk];
      CreateTrack( trackAR );
    }
    std::cout << "End of creating 4 hits tracks " << std::endl;

    nTrk = int(_tracks3Hits.size());
    for (int iTrk=0; iTrk<nTrk;++iTrk) {
      TrackExtended * trackAR = _tracks3Hits[iTrk];
      CreateTrack( trackAR );
    }
    std::cout << "End of creating 3 hits tracks " << std::endl;

    if (_attachFast == 0) {
      AttachRemainingVTXHitsSlow();
      AttachRemainingFTDHitsSlow();
    }
    else {
      AttachRemainingVTXHitsFast();
      AttachRemainingFTDHitsFast();
    }
    std::cout << "End of picking up remaining hits " << std::endl;

    LCCollectionVec * trkCol = new LCCollectionVec(LCIO::TRACK);
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
	for (int iHit=0;iHit<nHits;++iHit) {
	  TrackerHit * trkHit = hitVec[iHit]->getTrackerHit();
	  trackImpl->addHit(trkHit);
	}
	trkCol->addElement(trackImpl);
      }
    }        
    evt->addCollection(trkCol,_colVTXTracks );     
  }
  CleanUp();

  _nEvt++;


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

int VertexTracking::InitialiseFTD(LCEvent * evt) {
  int success = 1;
  _sectorsFTD.clear();
  _sectorsFTD.resize(2*_nLayersFTD*_nPhiFTD);
  for (int i=0; i<2*_nLayersFTD*_nPhiFTD;++i) {
    TrackerHitExtendedVec hitVec;
    hitVec.clear();
    _sectorsFTD.push_back(hitVec);    
  }
  try {
    LCCollection * hitCollection = evt->getCollection(_FTDHitCollection.c_str());
    int nelem = hitCollection->getNumberOfElements();
    for (int ielem=0; ielem<nelem; ++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
      hitExt->setResolutionRPhi(_resolutionRPhiFTD);
      hitExt->setResolutionZ(_resolutionZFTD);
      double pos[3];
      for (int i=0; i<3; ++i) {
	pos[i] = hit->getPosition()[i];
      }
      double Phi = atan2(pos[1],pos[0]);
      if (Phi < 0.) Phi = Phi + TWOPI;
      int layer = abs(hit->getType()) - 1;
      if (layer < 0 || layer > _nLayersFTD-1) {
	std::cout << "VertexTracking => fatal error in FTD : layer is outside allowed range : " << layer << std::endl;
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

int VertexTracking::Initialise(LCEvent * evt) {

  int success = 1;

  _sectors.clear();
  _sectors.resize(_nLayers*_nDivisionsInPhi*_nDivisionsInTheta);

  for (int i=0; i<_nLayers*_nDivisionsInPhi*_nDivisionsInTheta; ++i) {
    TrackerHitExtendedVec hitVec;
    hitVec.clear();
    _sectors.push_back(hitVec);
  }
    
  try {
    LCCollection * hitCollection = evt->getCollection(_VTXHitCollection.c_str());
    int nelem = hitCollection->getNumberOfElements();
    for (int ielem=0; ielem<nelem; ++ielem) {
      TrackerHit * hit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
      hitExt->setResolutionRPhi(_resolutionRPhi);
      hitExt->setResolutionZ(_resolutionZ);
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
      if (layer < 0 || layer >= _nLayers) {
	std::cout << "VertexTracking => fatal error in FTD : layer is outside allowed range : " << layer << std::endl;
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
		  TrackExtended * trackAR = TestTriplet(outerHit,middleHit,innerHit,helix);
		  if ( trackAR != NULL ) {
		    int nH = BuildTrack(outerHit,middleHit,innerHit,helix,nLR[2],iPhi,iTheta,trackAR);
		    if (nH == 3) 
		      _tracks3Hits.push_back(trackAR);
		    if (nH == 4)
		      _tracks4Hits.push_back(trackAR);
		    if (nH >= 5)
		      _tracks5Hits.push_back(trackAR);
		  }
		} // endloop over hits in the Inner sector
	      } // endloop over theta in the Inner
	    } // endloop over phi in the Inner	    
	  } // endloop over hits in the Middle sector
	} // endloop over theta in the Middle
      } // endloop over phi in the Middle
    } // endloop over outer hits
  } // endloop over triplets
}

TrackExtended * VertexTracking::TestTriplet(TrackerHitExtended * outerHit, 
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
  wrh[0] = double(1.0/(outerHit->getResolutionRPhi()*outerHit->getResolutionRPhi()));
  wzh[0] = 1.0/(outerHit->getResolutionZ()*outerHit->getResolutionZ());
  
  xh[1] = middleHit->getTrackerHit()->getPosition()[0];
  yh[1] = middleHit->getTrackerHit()->getPosition()[1];
  zh[1] = float(middleHit->getTrackerHit()->getPosition()[2]);
  wrh[1] = double(1.0/(middleHit->getResolutionRPhi()*middleHit->getResolutionRPhi()));
  wzh[1] = 1.0/(middleHit->getResolutionZ()*middleHit->getResolutionZ());
  
  xh[2] = innerHit->getTrackerHit()->getPosition()[0];
  yh[2] = innerHit->getTrackerHit()->getPosition()[1];
  zh[2] = float(innerHit->getTrackerHit()->getPosition()[2]);
  wrh[2] = double(1.0/(innerHit->getResolutionRPhi()*innerHit->getResolutionRPhi()));
  wzh[2] = 1.0/(innerHit->getResolutionZ()*innerHit->getResolutionZ());
  
  for (int ih=0; ih<3; ih++) {
    rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
    ph[ih] = atan2(yh[ih],xh[ih]);
    if (ph[ih] < 0.) 
      ph[ih] = 2.0*acos(-1.0) + ph[ih]; 
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
  float omega = par[0];
  float tanlambda = par[1];
  float phi0 = par[2];
  float d0 = par[3];
  float z0 = par[4];
 
  float Chi2 = chi2RPhi/_chi2WRPhiTriplet+chi2Z/_chi2WZTriplet;
  // Check if track satisfies all conditions
  if ( Chi2 > 2.0 || fabs(d0) > _cutOnD0 || fabs(z0) > _cutOnZ0 || omega > _cutOnOmega)
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
  /**
     Method for building up track in the VXD. Method starts from the found triplet and performs
     sequential attachment of hits in other layers 
   */


  int iPhi_Up    = iPhi + 1;
  int iPhi_Low   = iPhi - 1;
  int iTheta_Up  = iTheta + 1; 
  int iTheta_Low = iTheta - 1;
  if (iTheta_Low < 0) iTheta_Low = 0;
  if (iTheta_Up  >= _nDivisionsInTheta) iTheta_Up = _nDivisionsInTheta-1;

  for (int layer = innerLayer-1; layer>=0; layer--) { // loop over remaining layers
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

void VertexTracking::CreateTrack(TrackExtended * trackAR ) {

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
	float rR = hitVec[ih]->getResolutionRPhi();
	float rZ = hitVec[ih]->getResolutionZ();
	if (int(hitVec[ih]->getTrackExtendedVec().size()) != 0)
	  std::cout << "WARNING : HIT POINTS TO TRACK " << std::endl;
	xh[ih] = trkHit->getPosition()[0];
	yh[ih] = trkHit->getPosition()[1];
	zh[ih] = float(trkHit->getPosition()[2]);
	wrh[ih] = double(1.0/(rR*rR));
	wzh[ih] = 1.0/(rZ*rZ);
	rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
	ph[ih] = float(atan2(yh[ih],xh[ih]));
      }      
      for (int ih=0;ih<nHitsOld;++ih) {
	TrackerHit * trkHit = hitVecOld[ih]->getTrackerHit();
	xh[ih+nHits] = trkHit->getPosition()[0];
	yh[ih+nHits] = trkHit->getPosition()[1];
	zh[ih+nHits] = float(trkHit->getPosition()[2]);
	float rR = hitVecOld[ih]->getResolutionRPhi();
	float rZ = hitVecOld[ih]->getResolutionZ();	
	wrh[ih+nHits] = double(1.0/(rR*rR));
	wzh[ih+nHits] = 1.0/(rZ*rZ);
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
	float * wzhOld = new float[nTotHits];
	double * wrhOld = new double[nTotHits];
	for (int i=0;i<nTotHits;++i) {
	  wzhOld[i] = wzh[i];
	  wrhOld[i] = wrh[i];
	}
	for (int i=0; i<nTotHits; ++i) {
	  for (int j=0;j<nTotHits;++j) {
	    if (i == j) {
	      wrh[j] = 0.0;
	      wzh[j] = 0.0;
	    } 
	    else {
	      wrh[j] = wrhOld[j];
	      wzh[j] = wzhOld[j];
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
	delete[] wzhOld;
	delete[] wrhOld;
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

void VertexTracking::AttachRemainingVTXHitsFast() {

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

void VertexTracking::AttachRemainingVTXHitsSlow() {
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

void VertexTracking::AttachRemainingFTDHitsSlow() {
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


void VertexTracking::AttachRemainingFTDHitsFast() {
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

void VertexTracking::TrackingInFTD() {
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
	  for (int ipMiddle=ipMiddleLow;ipMiddle<=ipMiddleUp;++ipMiddle) {
	    // for(int ipMiddle=0;ipMiddle<_nPhiFTD;++ipMiddle) {
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
	      for (int ipInner=ipInnerLow;ipInner<=ipInnerUp;++ipInner) {
	      //for (int ipInner=0;ipInner<_nPhiFTD;++ipInner) {
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


int VertexTracking::BuildTrackFTD(TrackExtended * trackAR, int * nLR, int iS) {
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

int VertexTracking::AttachHitToTrack(TrackExtended * trackAR, TrackerHitExtended * hit) {

  int attached = 0;
  TrackerHitExtendedVec hitVec = trackAR->getTrackerHitExtendedVec();
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
    float rR = hitVec[i]->getResolutionRPhi();
    float rZ = hitVec[i]->getResolutionZ();
    wrh[i] = double(1.0/(rR*rR));
    wzh[i] = 1.0/(rZ*rZ);	  
  }
  TrackerHit * trkHit = hit->getTrackerHit();
  xh[nHits] = double(trkHit->getPosition()[0]);
  yh[nHits] = double(trkHit->getPosition()[1]);
  zh[nHits] = float(trkHit->getPosition()[2]);
  ph[nHits] = float(atan2(yh[nHits],xh[nHits]));
  rh[nHits] = float(sqrt(xh[nHits]*xh[nHits]+yh[nHits]*yh[nHits]));
  float rR = hit->getResolutionRPhi();
  float rZ = hit->getResolutionZ();
  wrh[nHits] = double(1.0/(rR*rR));
  wzh[nHits] = 1.0/(rZ*rZ);	
  
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
  
  //  std::cout << "Chi2RPhi = " << chi2RPhi << "   Chi2Z = " << chi2Z << std::endl; 

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
    trackAR->addTrackerHitExtended(hit);
    hit->addTrackExtended( trackAR );
    trackAR->setChi2( chi2RPhi + chi2Z );
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

  return attached;


}
