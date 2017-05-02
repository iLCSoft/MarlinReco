#include "CLICPfoSelector.h"
#include <iostream>
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <iostream>
#include <math.h>
#include <map>
#include <marlin/Global.h>
#include <CalorimeterHitType.h>
#include "ClusterShapes.h"
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/VXDParameters.h>
#include <gear/GearParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/BField.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>

#include <limits>

using namespace lcio ;
using namespace marlin ;

const int precision = 2;
const int widthFloat = 7;
const int widthInt = 5;

CLICPfoSelector aCLICPfoSelector ;

CLICPfoSelector::CLICPfoSelector() : Processor("CLICPfoSelector") {  
  _description = "Selects Pfos from full PFO list using timing cuts" ;  

  // Input pfo collections

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "InputPfoCollection",
			  "Input PFO Collection",
			  m_inputPfoCollection,
			  std::string("PandoraPFANewPFOs"));
  
  // Output pfo collection
  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			   "SelectedPfoCollection",
			   "Selected pfo collection name",
			   m_selectedPfoCollection,
			   std::string("SelectedPandoraPFANewPFOs"));

  registerProcessorParameter("Debug",
			     "Activate debugging?",
			     m_debug,
			     int(0));

  registerProcessorParameter("Monitoring",
			     "Monitoring",
			     m_monitoring,
			     int(0));

  registerProcessorParameter("DisplaySelectedPfos",
			     "DisplaySelectedPfos",
			     m_displaySelectedPfos,
			     int(0));

  registerProcessorParameter("DisplayRejectedPfos",
			     "DisplayRejectedPfos",
			     m_displayRejectedPfos,
			     int(0));

  registerProcessorParameter("CorrectHitTimesForTimeOfFlight",
			     "CorrectHitTimesForTimeOfFlight",
                             m_correctHitTimesForTimeOfFlight, 
                             int(0));

  registerProcessorParameter("MonitoringPfoEnergyToDisplay",
			     "MonitoringPfoEnergyToDisplay",
			     m_monitoringPfoEnergyToDisplay,
			     float(1.0));

  registerProcessorParameter("CheckProtonCorrection",
			     "CheckProtonCorrection",
                             m_checkProtonCorrection, 
                             int(0));


  registerProcessorParameter("CheckKaonCorrection",
			     "CheckKaonCorrection",
                             m_checkKaonCorrection, 
                             int(0));

  registerProcessorParameter("KeepKShorts",
			     "KeepKShorts",
                             m_keepKShorts, 
                             int(1));
 
  registerProcessorParameter("UseNeutronTiming",
			     "UseNeutronTiming",
                             m_useNeutronTiming,
			     int(0));
  
  registerProcessorParameter("MinimumEnergyForNeutronTiming",
			     "MinimumEnergyForNeutronTiming",
			     m_minimumEnergyForNeutronTiming,
			     float(1.));

  registerProcessorParameter("ForwardCosThetaForHighEnergyNeutralHadrons",
			     "ForwardCosThetaForHighEnergyNeutralHadrons",
			     m_forwardCosThetaForHighEnergyNeutralHadrons,
			     float(0.95));

  registerProcessorParameter("ForwardHighEnergyNeutralHadronsEnergy",
			     "ForwardHighEnergyNeutralHadronsEnergy",
			     m_forwardHighEnergyNeutralHadronsEnergy,
			     float(10.00));

  registerProcessorParameter("FarForwardCosTheta",
			     "FarForwardCosTheta",
			     m_farForwardCosTheta,
			     float(0.975));
  
  registerProcessorParameter("PtCutForTightTiming",
			     "PtCutForTightTiming",
			     m_ptCutForTightTiming,
			     float(0.75));
  
  registerProcessorParameter("PhotonPtCut",
			     "PhotonPtCut",
			     m_photonPtCut,
			     float(0.0));

  registerProcessorParameter("PhotonPtCutForLooseTiming",
			     "PhotonPtCutForLooseTiming",
			     m_photonPtCutForLooseTiming,
			     float(4.0));
  
  registerProcessorParameter("PhotonLooseTimingCut",
			     "PhotonLooseTimingCut",
			     m_photonLooseTimingCut,
			     float(2.0));
  
  registerProcessorParameter("PhotonTightTimingCut",
			     "PhotonTightTimingCut",
			     m_photonTightTimingCut,
			     float(1.0));

  registerProcessorParameter("ChargedPfoPtCut",
			     "ChargedPfoPtCut",
			     m_chargedPfoPtCut,
			     float(0.0));
 
  registerProcessorParameter("ChargedPfoPtCutForLooseTiming",
			     "ChargedPfoPtCutForLooseTiming",
			     m_chargedPfoPtCutForLooseTiming,
			     float(4.0));

  registerProcessorParameter("ChargedPfoLooseTimingCut",
			     "ChargedPfoLooseTimingCut",
			     m_chargedPfoLooseTimingCut,
			     float(3.0));

  registerProcessorParameter("ChargedPfoTightTimingCut",
			     "ChargedPfoTightTimingCut",
			     m_chargedPfoTightTimingCut,
			     float(1.5)); 
 
  registerProcessorParameter("ChargedPfoNegativeLooseTimingCut",
			     "ChargedPfoNegativeLooseTimingCut",
			     m_chargedPfoNegativeLooseTimingCut,
			     float(-1.0));

  registerProcessorParameter("ChargedPfoNegativeTightTimingCut",
			     "ChargedPfoNegativeTightTimingCut",
			     m_chargedPfoNegativeTightTimingCut,
			     float(-0.5));

  registerProcessorParameter("NeutralHadronPtCut",
			     "NeutralHadronPtCut",
			     m_neutralHadronPtCut,
			     float(0.0));
 
  registerProcessorParameter("NeutralHadronPtCutForLooseTiming",
			     "NeutralHadronPtCutForLooseTiming",
			     m_neutralHadronPtCutForLooseTiming,
			     float(8.0));
 
  registerProcessorParameter("NeutralHadronLooseTimingCut",
			     "NeutralHadronLooseTimingCut",
			     m_neutralHadronLooseTimingCut,
			     float(2.5));
  
  registerProcessorParameter("NeutralHadronTightTimingCut",
			     "NeutralHadronTightTimingCut",
			     m_neutralHadronTightTimingCut,
			     float(1.5));
  
  registerProcessorParameter("NeutralFarForwardLooseTimingCut",
			     "NeutralFarForwardLooseTimingCut",
			     m_neutralFarForwardLooseTimingCut,
			     float(2.0));
  
  registerProcessorParameter("NeutralFarForwardTightTimingCut",  
			     "NeutralFarForwardTightTimingCut",  
			     m_neutralFarForwardTightTimingCut,
			     float(1.0));
 
  registerProcessorParameter("PhotonFarForwardLooseTimingCut",
			     "PhotonFarForwardLooseTimingCut",
			     m_photonFarForwardLooseTimingCut,
			     float(2.0));
  
  registerProcessorParameter("PhotonFarForwardTightTimingCut",  
			     "PhotonFarForwardTightTimingCut",  
			     m_photonFarForwardTightTimingCut,
			     float(1.0));
 
  registerProcessorParameter("HCalBarrelLooseTimingCut",
			     "HCalBarrelLooseTimingCut",
			     m_hCalBarrelLooseTimingCut,
			     float(20.0));
 
  registerProcessorParameter("HCalBarrelTightTimingCut",
			     "HCalBarrelTightTimingCut",
			     m_hCalBarrelTightTimingCut,
			     float(10.0));

  registerProcessorParameter("HCalEndCapTimingFactor",
			     "HCalEndCapTimingFactor",
			     m_hCalEndCapTimingFactor,
			     float(1.0));

  registerProcessorParameter("NeutralHadronBarrelPtCutForLooseTiming",
			     "NeutralHadronBarrelPtCutForLooseTiming",
			     m_neutralHadronBarrelPtCutForLooseTiming,
			     float(3.5));

  registerProcessorParameter("MinECalHitsForTiming",
			     "MinECalHitsForTiming",
			     m_minECalHitsForTiming,
			     int(5));
 
  registerProcessorParameter("MinHCalEndCapHitsForTiming",
			     "MinHCalEndCapHitsForTiming",
			     m_minHCalEndCapHitsForTiming,
			     int(5));

  registerProcessorParameter("UseClusterLessPfos",
			     "UseClusterLessPfos",
			     m_useClusterLessPfos,
			     int(1));
 
  registerProcessorParameter("MinMomentumForClusterLessPfos",
			     "MinMomentumForClusterLessPfos",
			     m_minMomentumForClusterLessPfos,
			     float(0.5));
 

  registerProcessorParameter("MaxMomentumForClusterLessPfos",
			     "MaxMomentumForClusterLessPfos",
			     m_maxMomentumForClusterLessPfos,
			     float(2.0));
 

  registerProcessorParameter("MinPtForClusterLessPfos",
			     "MinPtForClusterLessPfos",
			     m_minPtForClusterLessPfos,
			     float(0.5));

  registerProcessorParameter("ClusterLessPfoTrackTimeCut",
			     "ClusterLessPfoTrackTimeCut",
			     m_clusterLessPfoTrackTimeCut,
			     float(10.0));
  

}



void CLICPfoSelector::init() { 

  printParameters();  
  _nRun = -1 ;
  _nEvt = 0 ;

}

void CLICPfoSelector::processRunHeader( LCRunHeader* ) {

  _nRun++ ;
  _nEvt = 0;
  streamlog_out( MESSAGE ) << std::endl;
  streamlog_out( MESSAGE ) << "CLICPfoSelector ---> new run : run number = " << _nRun << std::endl;

} 

void CLICPfoSelector::processEvent( LCEvent * evt ) { 

  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;

  streamlog_out( DEBUG ) << "CLICPfoSelector -> run = " << _nRun << "  event = " << _nEvt << std::endl;

  LCCollectionVec * colPFO = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  colPFO->setSubset(true);

  // if we want to point back to the hits we need to set the flag

  
  // Reading PFOs

  try {
    LCCollection * col = evt->getCollection(m_inputPfoCollection.c_str());
    int nelem = col->getNumberOfElements();
    PfoUtil::PfoList pfos;

    for (int iPfo=0; iPfo<nelem; ++iPfo)pfos.push_back(dynamic_cast<ReconstructedParticle*>(col->getElementAt(iPfo)));
    std::sort(pfos.begin(),pfos.end(),PfoUtil::PfoSortFunction);

    if (m_monitoring) {
      streamlog_out( MESSAGE ) << std::endl;
      streamlog_out( MESSAGE ) << "Number of Input Pfos = " << nelem << std::endl;
      streamlog_out( MESSAGE ) << "    Type          PDG    E      Pt  cosTheta #trk time  #Clu  time   ecal  hcal  " << std::endl;
    }
//    int nDropped(0);

    float eTotalInput(0.);
    float eTotalOutput(0.);

    for (int iPfo=0; iPfo<nelem; ++iPfo) {
      bool passPfoSelection = true;
      ReconstructedParticle * pPfo = pfos[iPfo];
//      const int id = pPfo->id();
      const int type = pPfo->getType();
//      const bool isCompound = pPfo->isCompound(); 
      float momentum[3];
      for(unsigned int i=0;i<3;i++)momentum[i] = pPfo->getMomentum()[i];
      const float pT_pfo = sqrt(momentum[0]*momentum[0]+momentum[1]*momentum[1]);
      const float p_pfo  = sqrt(pT_pfo*pT_pfo+momentum[2]*momentum[2]);
      const float cosTheta = fabs(momentum[2])/p_pfo;
      const float energy  = pPfo->getEnergy();
      eTotalInput+=energy;
      // float covMatrix[10];
      // for(unsigned int i=0;i<10;i++)covMatrix[i] = pPfo->getCovMatrix()[i];
//      const float mass = pPfo->getMass();
//      const float charge = pPfo->getCharge();
      // float referencePoint[3];
      // for(unsigned int i=0;i<3;i++)referencePoint[i] = pPfo->getReferencePoint()[i];
      //const ParticleIDVec particleIDs = pPfo->getParticleIDs();
//      ParticleID *particleIDUsed = pPfo->getParticleIDUsed();
//      const float goodnessOfPID = pPfo->getGoodnessOfPID();
      //const ReconstructedParticleVec particles = pPfo->getParticles();
      const ClusterVec clusters = pPfo->getClusters();
      const TrackVec   tracks   = pPfo->getTracks();
      //const Vertex startVertex(pPfo->getStartVertex());

      float trackTime = std::numeric_limits<float>::max();
      float clusterTime = 999.;
      float clusterTimeEcal = 999.;
      float clusterTimeHcalEndcap = 999.;
      int   nEcalHits(0);
      int   nHcalEndCapHits(0);
      int   nCaloHits(0);
      float tproton(0.);
      float tkaon(0.);

      for(unsigned int i = 0; i< tracks.size(); i++){
	const Track *track = tracks[i];
	const TrackerHitVec hitVec = track->getTrackerHits();
	const int nHits = int(hitVec.size());
	//	const int nHitsVTX = track->getSubdetectorHitNumbers()[6];
	//const int nHitsFTD = track->getSubdetectorHitNumbers()[7];
	//const int nHitsSIT = track->getSubdetectorHitNumbers()[8];
	//const int nHitsTPC = track->getSubdetectorHitNumbers()[9];
	//const int nHitsSET = track->getSubdetectorHitNumbers()[10];
	//const int nHitsETD = track->getSubdetectorHitNumbers()[11];
	float r2Min = std::numeric_limits<float>::max(); 
	const float d0    = track->getD0();
	const float z0    = track->getZ0();
	const float omega = track->getOmega();
	const float tanL  = track->getTanLambda();
	const float phi0  = track->getPhi();
	
	HelixClass helix;
	helix.Initialize_Canonical(phi0,d0,z0,omega,tanL,_bField);
	const float px = helix.getMomentum()[0];
	const float py = helix.getMomentum()[1];
	const float pz = helix.getMomentum()[2];
	const float pT = sqrt(px*px+py*py);
	const float p  = sqrt(pT*pT+pz*pz);
	float tof;
	const float time = this->TimeAtEcal(track,tof);
	for (int iH=0;iH<nHits;++iH) {
	  const TrackerHit * hit = hitVec[iH];
	  const float hitX = float(hit->getPosition()[0]);
	  const float hitY = float(hit->getPosition()[1]);
	  const float hitR2 = hitX*hitX+hitY*hitY;
	  if (hitR2<r2Min)r2Min = hitR2;
	}
	if(fabs(time)<trackTime){
	  trackTime = time;
	  const float cproton = sqrt((p*p+0.94*0.94)/(p*p+0.14*0.14));
	  const float ckaon = sqrt((p*p+0.49*0.49)/(p*p+0.14*0.14));
	  tproton = (trackTime+tof)*(cproton-1);
	  tkaon   = (trackTime+tof)*(ckaon-1);
	}
      }


      for(unsigned int i = 0; i< clusters.size(); i++){
	float meanTime(999.);
	float meanTimeEcal(999.);
	float meanTimeHcalEndcap(999.);
	int   nEcal(0);
	int   nHcalEnd(0);
	int   nCaloHitsUsed(0);

	const Cluster *cluster = clusters[i];
	this->GetClusterTimes(cluster,meanTime,nCaloHitsUsed,meanTimeEcal,nEcal,meanTimeHcalEndcap,nHcalEnd);

	// correct for track propagation time
	if(!tracks.empty()){
	  meanTime -= trackTime;
	  meanTimeEcal -= trackTime;
	  meanTimeHcalEndcap -= trackTime;
	}

	if(fabs(meanTime)<clusterTime){
	  clusterTime=meanTime;
	  nCaloHits = nCaloHitsUsed;
	}
	if(fabs(meanTimeEcal)<clusterTimeEcal){
	  clusterTimeEcal=meanTimeEcal;
	  nEcalHits = nEcal;
	}
	if(fabs(meanTimeHcalEndcap)<clusterTimeHcalEndcap){
	  clusterTimeHcalEndcap=meanTimeHcalEndcap;
	  nHcalEndCapHits = nHcalEnd;
	}
      }

      // now make selection

      float ptCut(m_neutralHadronPtCut);
      float ptCutForLooseTiming(m_neutralHadronPtCutForLooseTiming);
      float timingCutLow(0.);
      float timingCutHigh(m_neutralHadronLooseTimingCut);
      float hCalBarrelTimingCut(m_hCalBarrelLooseTimingCut);
      if (cosTheta > m_farForwardCosTheta)
	timingCutHigh = m_neutralFarForwardLooseTimingCut;
      
      // Neutral hadron cuts
      if (pT_pfo <= m_ptCutForTightTiming)
        {
	  timingCutHigh = m_neutralHadronTightTimingCut;
	  hCalBarrelTimingCut = m_hCalBarrelTightTimingCut;
	  if (cosTheta > m_farForwardCosTheta)
	    timingCutHigh = m_neutralFarForwardTightTimingCut;
        }
      
      // Photon cuts
      if (type == 22)
        {
	  ptCut = m_photonPtCut;
	  ptCutForLooseTiming = m_photonPtCutForLooseTiming;
	  timingCutHigh = m_photonLooseTimingCut;
	  if (cosTheta > m_farForwardCosTheta)
	    timingCutHigh = m_photonFarForwardLooseTimingCut;

	  if (pT_pfo <= m_ptCutForTightTiming){
	    timingCutHigh = m_photonTightTimingCut;
	    if (cosTheta > m_farForwardCosTheta)
	      timingCutHigh = m_photonFarForwardTightTimingCut;
	  }
        }
      
      // Charged PFO cuts
      if (!tracks.empty())
        {
	  ptCut = m_chargedPfoPtCut;
	  ptCutForLooseTiming = m_chargedPfoPtCutForLooseTiming;
	  timingCutLow = m_chargedPfoNegativeLooseTimingCut;
	  timingCutHigh = m_chargedPfoLooseTimingCut;
	  if (pT_pfo <= m_ptCutForTightTiming)
            {
	      timingCutLow = m_chargedPfoNegativeTightTimingCut;
	      timingCutHigh = m_chargedPfoTightTimingCut;
            }
        }
 
   
      // Reject low pt pfos (default is to set ptcut to zero)
      if (pT_pfo < ptCut)
        {
	  passPfoSelection = false;
	}

      // Reject out of time clusterless tracks
      if (clusters.empty() && fabs(trackTime) > m_clusterLessPfoTrackTimeCut)
        {
	  passPfoSelection = false;
	}


      const float pfoEnergyToDisplay(m_monitoring ? m_monitoringPfoEnergyToDisplay : std::numeric_limits<float>::max());
      


      // Only apply cuts to low pt pfos and very forward neutral hadrons
      bool applyTimingCuts = ( (pT_pfo < ptCutForLooseTiming) || ( (cosTheta > m_forwardCosThetaForHighEnergyNeutralHadrons) && (type == 2112) ) );
      bool useHcalTimingOnly = ( (cosTheta > m_forwardCosThetaForHighEnergyNeutralHadrons) && (type == 2112) && (energy > m_forwardHighEnergyNeutralHadronsEnergy));
    

      if (passPfoSelection && applyTimingCuts){

	// Examine any associated clusters for additional timing information
	bool selectPfo(false);
	// Require any cluster to be "in time" to select pfo
	if (!clusters.empty())
	  {
	    // Make the selection decisions
	    if (!useHcalTimingOnly && ((nEcalHits > m_minECalHitsForTiming) || (nEcalHits >= nCaloHits/2.)) )
	      {
		if ((clusterTimeEcal >= timingCutLow) && (clusterTimeEcal <= timingCutHigh))
		  selectPfo = true;
	      }
	    else if (type == 22)
	      {
		if ((clusterTime >= timingCutLow) && (clusterTime <= timingCutHigh))
		  selectPfo = true;
	      }
	    else if ( (nHcalEndCapHits >= m_minHCalEndCapHitsForTiming) || (nHcalEndCapHits >= nCaloHits/2.) )
	      {
		if ((clusterTimeHcalEndcap >= timingCutLow) && (clusterTimeHcalEndcap <= (m_hCalEndCapTimingFactor * timingCutHigh)))
		  selectPfo = true;
	      }
	    else
	      {
		if ((clusterTime >= timingCutLow) && (clusterTime < hCalBarrelTimingCut))
		  selectPfo = true;
		
		if (tracks.empty() && (pT_pfo > m_neutralHadronBarrelPtCutForLooseTiming))
		  selectPfo = true;
	      }



	    // keep KShorts
	    if(m_keepKShorts && type == 310){
	      if(!selectPfo && m_monitoring && m_displayRejectedPfos)streamlog_out( MESSAGE ) <<   " Recovered KS     : " << energy << std::endl;
	      selectPfo = true;
	    }
	    // check kaon and proton hypotheses
	    if (nEcalHits > m_minECalHitsForTiming)
	    {
	        if(m_checkProtonCorrection && (clusterTimeEcal-tproton >= timingCutLow) && (clusterTimeEcal-tproton<= timingCutHigh))
                {
		    if(!selectPfo && m_monitoring  && m_displayRejectedPfos)streamlog_out( MESSAGE ) << " Recovered proton : " << energy << std::endl;
		    selectPfo = true;
		}
		if(m_checkKaonCorrection &&   (clusterTimeEcal-tkaon >= timingCutLow) && (clusterTimeEcal-tkaon<= timingCutHigh))
                {
		    if(!selectPfo && m_monitoring && m_displayRejectedPfos)streamlog_out( MESSAGE ) << " Recovered kaon   : " << energy << std::endl;
		    selectPfo = true;
		}
	    }
	  }
	else
	  {
	    // No clusters form part of this pfo - no additional timing information
	    if (p_pfo > m_minMomentumForClusterLessPfos &&
                p_pfo < m_maxMomentumForClusterLessPfos &&
                pT_pfo > m_minPtForClusterLessPfos)
	      selectPfo = m_useClusterLessPfos;
	  }	 

	if(!selectPfo)passPfoSelection = false;
      }
      
      if (m_monitoring && (energy > pfoEnergyToDisplay)) {
        if( (passPfoSelection && m_displaySelectedPfos) || (!passPfoSelection && m_displayRejectedPfos)) {
          std::stringstream output;
          output << std::fixed;
          output << std::setprecision(precision);
          if(passPfoSelection) {
            output << " Selected PFO : ";
          } else {
            output << " Rejected PFO : ";
          }
          if(clusters.size()==0)
            FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT_pfo,cosTheta,tracks.size(),trackTime,"-","-","-","-");
          if(tracks.size()==0)
            FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT_pfo,cosTheta,"","-",clusters.size(),clusterTime,clusterTimeEcal,clusterTimeHcalEndcap);
          if(tracks.size()>0&&clusters.size()>0)
            FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT_pfo,cosTheta,tracks.size(),trackTime,clusters.size(),clusterTime,clusterTimeEcal,clusterTimeHcalEndcap);
          streamlog_out( MESSAGE ) << output.str();
        }
      }      
      if(passPfoSelection){
	eTotalOutput+=energy;
	colPFO->addElement(pPfo);
      }else{
	//streamlog_out( MESSAGE ) << " dropped E = " << energy << std::endl;
      }
      
    }
    
    if(m_monitoring){
      streamlog_out( MESSAGE ) << " Total PFO energy in  = " << eTotalInput << " GeV " << std::endl;
      streamlog_out( MESSAGE ) << " Total PFO energy out = " << eTotalOutput << " GeV " << std::endl;
    }
  }
  catch( DataNotAvailableException &e ) {
    streamlog_out( MESSAGE ) << m_inputPfoCollection.c_str() << " collection is unavailable" << std::endl;
  };
  
  
  evt->addCollection(colPFO,m_selectedPfoCollection.c_str());
  


  CleanUp();
  
  _nEvt++;

}



float CLICPfoSelector::TimeAtEcal(const Track* pTrack, float &tof){


  const gear::CalorimeterParameters& pEcalBarrel = Global::GEAR->getEcalBarrelParameters();
  const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
  // determine pseudo-geometry of detector (determined by ECAL barrel)
  // symmetry 0 =cylinder, 1=prorotype, >1 = polygon
  const int symmetry = pEcalBarrel.getSymmetryOrder();
  const float zOfEndCap = (float)pEcalEndcap.getExtent()[2];
  const float phi0 = (float)pEcalBarrel.getPhi0();
  const float rBarrel = (float)pEcalBarrel.getExtent()[0];

  
  HelixClass helix;
  helix.Initialize_Canonical(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(), pTrack->getTanLambda(), _bField);
  
  const EVENT::TrackerHitVec &trackerHitvec(pTrack->getTrackerHits());
  float zMin(std::numeric_limits<float>::max()), zMax(-std::numeric_limits<float>::max());
  
  for (int iz = 0, nTrackHits = trackerHitvec.size(); iz < nTrackHits - 1; ++iz)
    {
      const float hitZ(trackerHitvec[iz]->getPosition()[2]);

      if (hitZ > zMax)
	zMax = hitZ;
      
      if (hitZ < zMin)
	zMin = hitZ;
    }
  
  const int signPz(fabs(zMin) < fabs(zMax) ? 1 : -1);
  
  float referencePoint[3];
  referencePoint[0] = helix.getReferencePoint()[0];
  referencePoint[1] = helix.getReferencePoint()[1];
  referencePoint[2] = helix.getReferencePoint()[2];

  
  // First project to endcap
  float minTime(std::numeric_limits<float>::max());
  //bool isProjectedToEndCap(true);
  
  float bestECalProjection[3];
  minTime =  helix.getPointInZ(static_cast<float>(signPz) * zOfEndCap, referencePoint, bestECalProjection);

  // Then project to barrel surface(s)
  static const float pi(acos(-1.));
  float barrelProjection[3];
  
  // n-sided Polygon
  float twopi_n = 2. * pi / (static_cast<float>(symmetry));
  
  for (int i = 0; i < symmetry; ++i)
    {
      float time(std::numeric_limits<float>::max());
      const float phi(twopi_n * static_cast<float>(i) + phi0);
      
      time = helix.getPointInXY(rBarrel * cos(phi), rBarrel * sin(phi),
			   cos(phi + 0.5 * pi), sin(phi + 0.5 * pi), referencePoint, barrelProjection);
      
      if ((time < minTime))
	{
	    minTime = time;
	    //isProjectedToEndCap = false;
	    bestECalProjection[0] = barrelProjection[0];
	    bestECalProjection[1] = barrelProjection[1];
	    bestECalProjection[2] = barrelProjection[2];

	}
    }
  
  tof    = sqrt( bestECalProjection[0]*bestECalProjection[0]+
		 bestECalProjection[1]*bestECalProjection[1]+
		 bestECalProjection[2]*bestECalProjection[2])/300;
  
  float px = helix.getMomentum()[0];
  float py = helix.getMomentum()[1];
  float pz = helix.getMomentum()[2];
  float E = sqrt(px*px+py*py+pz*pz+0.139*0.139);
  minTime = minTime/300*E-tof;

  return minTime;
  
}

 void CLICPfoSelector::GetClusterTimes(const Cluster* cluster, float &meanTime, int &nCaloHitsUsed, float &meanTimeEcal, int &nEcal, float &meanTimeHcalEndcap, int &nHcalEnd)
{

    // Calculate cluster times

//    const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
    // determine pseudo-geometry of detector (determined by ECAL barrel)
    // symmetry 0 =cylinder, 1=prorotype, >1 = polygon
//    const float zOfEndCap = (float)pEcalEndcap.getExtent()[2];

    float sumTimeEnergy(0.f);
    float sumEnergy(0.f);
    float sumEnergyEcal(0.f);
    float sumTimeEnergyEcal(0.f);
    float sumEnergyHcalEndcap(0.f);
    float sumTimeEnergyHcalEndcap(0.f);
    meanTime     = std::numeric_limits<float>::max();
    meanTimeEcal = std::numeric_limits<float>::max();
    meanTimeHcalEndcap = std::numeric_limits<float>::max();
    nEcal    = 0;
    nHcalEnd = 0;
    nCaloHitsUsed = 0;

    CalorimeterHitVec hits = cluster->getCalorimeterHits();
    std::vector<float> hittimes;
    std::vector<float> tofCorrections;
    std::vector<float> deltaTimes;

    for(unsigned int ihit=0;ihit<hits.size();++ihit){
      // optionally correct hit times for straight line tof (may have already been done in another processor)
      if(m_correctHitTimesForTimeOfFlight){
	const float x = hits[ihit]->getPosition()[0];
	const float y = hits[ihit]->getPosition()[1];
	const float z = hits[ihit]->getPosition()[2];
	const float r = sqrt(x*x+y*y+z*z);
	const float tof = r/300.;
	tofCorrections.push_back(tof);
	hittimes.push_back(hits[ihit]->getTime()-tof);
      }else{
	hittimes.push_back(hits[ihit]->getTime());
      }
    }

    std::sort(hittimes.begin(),hittimes.end());
    
    int iMedian = static_cast<int>(hits.size()/2.);
    float medianTime = hittimes[iMedian];
    //streamlog_out( MESSAGE ) << " Median time : " << medianTime << std::endl;

    for(unsigned int ihit=0;ihit<hits.size();++ihit)deltaTimes.push_back( fabs(hittimes[ihit]-medianTime)); 
    std::sort(deltaTimes.begin(),deltaTimes.end());
    
    unsigned ihit90 = 0;

    if (hits.size() > 1) {
      ihit90 = static_cast<int>((hits.size()*9)/10.);
      if(ihit90>=hits.size()-1)ihit90=hits.size()-2;
    } else {
      ihit90 = 0;
    }
  
    //streamlog_out( MESSAGE ) << " hits " << hits.size() << " hit 90 = " << ihit90 << std::endl;
    float deltaMedian = deltaTimes[ihit90]+0.1;
    //streamlog_out( MESSAGE ) << " deltaCut : " << deltaMedian << std::endl;

    for(unsigned int ihit=0;ihit<hits.size();++ihit){
      CalorimeterHit *hit = hits[ihit];
      float hitTime  = hits[ihit]->getTime();
      if(m_correctHitTimesForTimeOfFlight)hitTime -= tofCorrections[ihit];

      if( (hitTime - medianTime) < deltaMedian){
	sumEnergy += hit->getEnergy();
	sumTimeEnergy += hit->getEnergy()*hitTime;
	nCaloHitsUsed++;
	//streamlog_out( MESSAGE ) << " Using : " << hit->getEnergy() << " : " << hit->getTime() << std::endl;
	CHT ch = hit->getType();
	if(ch.is(CHT::ecal)){
	  nEcal++;
	  sumEnergyEcal += hit->getEnergy();
	  sumTimeEnergyEcal += hit->getEnergy()*hitTime;
	}else{
//	  float z = hit->getPosition()[2]; 
	  if(!ch.is(CHT::barrel)){
	    nHcalEnd++;
	    sumEnergyHcalEndcap += hit->getEnergy();
	    sumTimeEnergyHcalEndcap += hit->getEnergy()*hitTime;
	  }
	}
      }else{
	//streamlog_out( MESSAGE ) << " notus : " << hit->getEnergy() << " : " << hit->getTime() << std::endl;
      }
    }
    
    if (sumEnergy > 0.f)meanTime = sumTimeEnergy/sumEnergy;
    if (sumEnergyEcal > 0.f)meanTimeEcal = sumTimeEnergyEcal/sumEnergyEcal;
    if (sumEnergyHcalEndcap > 0.f)meanTimeHcalEndcap = sumTimeEnergyHcalEndcap/sumEnergyHcalEndcap;

    //streamlog_out( MESSAGE ) << sumEnergy << " " << sumEnergyEcal << " " << nEcal << std::endl;

    return;
}


void CLICPfoSelector::CleanUp(){
 
}

void CLICPfoSelector::check(LCEvent * ) { }

void CLICPfoSelector::end() {

}
