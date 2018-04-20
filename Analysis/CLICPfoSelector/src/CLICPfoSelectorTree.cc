#include "CLICPfoSelectorTree.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include "PfoUtilities.h"

#include "marlin/VerbosityLevels.h"

#include <marlin/AIDAProcessor.h>

using namespace lcio ;
using namespace marlin ;

const int precision = 2;
const int widthFloat = 7;
const int widthInt = 5;

CLICPfoSelectorTree aCLICPfoSelectorTree ;


CLICPfoSelectorTree::CLICPfoSelectorTree() : Processor("CLICPfoSelectorTree") {

  _description = "CLICPfoSelectorTree produces a tree with the variables used in the CLICPfoSelector." ;

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
          "PFOCollectionName" , 
          "Name of the PFO collection"  ,
          _colName_pfo,
          std::string("PandoraPFOs")
  );

}



void CLICPfoSelectorTree::init() { 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  AIDAProcessor::histogramFactory(this);

  m_tree = new TTree("PfoTree","PfoTree");
  m_tree->Branch("type", &type, "type/I");
  m_tree->Branch("p", &p, "p/D");
  m_tree->Branch("px", &px, "px/D");
  m_tree->Branch("py", &py, "py/D");
  m_tree->Branch("pz", &pz, "pz/D");
  m_tree->Branch("pT", &pT, "pT/D");

  m_tree->Branch("costheta", &costheta, "costhetaMC/D");
  m_tree->Branch("energy", &energy, "energy/D");
  m_tree->Branch("mass", &mass, "mass/D");
  m_tree->Branch("charge", &charge, "charge/D");

  m_tree->Branch("clusterTime", &clusterTime, "clusterTime/D");
  m_tree->Branch("clusterTimeEcal", &clusterTimeEcal, "clusterTimeEcal/D");
  m_tree->Branch("clusterTimeHcalEndcap", &clusterTimeHcalEndcap, "clusterTimeHcalEndcap/D");

  m_tree->Branch("nCaloHits", &nCaloHits, "nCaloHits/I");
  m_tree->Branch("nEcalHits", &nEcalHits, "nEcalHits/I");
  m_tree->Branch("nHcalEndCapHits", &nHcalEndCapHits, "nHcalEndCapHits/I");

  m_tree->Branch("eventNumber", &eventNumber, "eventNumber/I");

  m_tree->Branch("nPartPFO", &nPartPFO, "nPartPFO/I");
  m_tree->Branch("nPartMC", &nPartMC, "nPartMC/I");

  streamlog_out(DEBUG) << "   set up ttree  " << std::endl ;
}

void CLICPfoSelectorTree::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void CLICPfoSelectorTree::processEvent( LCEvent * evt ) { 

  streamlog_out( DEBUG ) << "Processing event " << evt->getEventNumber() << " in CLICPfoSelectorTree processor " << std::endl;

  eventNumber=evt->getEventNumber();

  fillTree(evt, _colName_pfo);

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
      << "   in run:  " << evt->getRunNumber() << std::endl ;

  _nEvt ++ ;
}

void CLICPfoSelectorTree::end(){ 
/*
  std::pair<double,double> resoPFO_tot = computeEnReso("energyPFO");
  streamlog_out(DEBUG) << "Total PFO energy resolution with error: " << resoPFO_tot.first << " , " << resoPFO_tot.second << std::endl;

  std::pair<double,double> resoMC_tot = computeEnReso("energyMC");
  streamlog_out(DEBUG) << "Total MC energy resolution with error: " << resoMC_tot.first << " , " << resoMC_tot.second << std::endl;

  const unsigned int nRegionBins(13);
  float pRegionBinEdges[nRegionBins + 1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975, 1.0};

  AIDAProcessor::histogramFactory(this);
  TGraphErrors* resoVsCostheta = new TGraphErrors(nRegionBins);
  resoVsCostheta->SetTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|");
  resoVsCostheta->SetName("enResoVsCosTheta");
  TGraphErrors* scaleVsCostheta = new TGraphErrors(nRegionBins);
  scaleVsCostheta->SetTitle("Mean(E_{TOT,PFO}/E_{TOT,MC}) vs |cos(#theta)|");
  scaleVsCostheta->SetName("enScaleVsCosTheta");

  // filling graph with energy resolution and scale
  std::pair<double,double> resoTemp;
  double scaleTemp;

  for (unsigned int i = 0; i < nRegionBins; ++i){
    std::ostringstream lowBin;
    std::ostringstream uppBin;

    lowBin << pRegionBinEdges[i];
    uppBin << pRegionBinEdges[i+1];

    std::string cosThetaCut = lowBin.str() + " <= costhetaMC && costhetaMC < " + uppBin.str(); 
    resoTemp = computeEnReso("energyPFO", cosThetaCut);

    if(resoTemp.first!=0.0) resoVsCostheta->SetPoint(i,(pRegionBinEdges[i+1]-pRegionBinEdges[i])/2+pRegionBinEdges[i], resoTemp.first);
    if(resoTemp.first!=0.0) resoVsCostheta->SetPointError(i,(pRegionBinEdges[i+1]-pRegionBinEdges[i])/2, resoTemp.second);

    scaleTemp = profiling("energyScale", cosThetaCut);
    scaleVsCostheta->SetPoint(i,(pRegionBinEdges[i+1]-pRegionBinEdges[i])/2+pRegionBinEdges[i], scaleTemp);
    scaleVsCostheta->SetPointError(i,(pRegionBinEdges[i+1]-pRegionBinEdges[i])/2, 0.0);
  }

  resoVsCostheta->Write();
  scaleVsCostheta->Write();
*/
  streamlog_out(DEBUG) << "CLICPfoSelectorTree::end()  " << name() 
    << " processed " << _nEvt << " events in " << _nRun << " runs " << std::endl;

}

void CLICPfoSelectorTree::fillTree(LCEvent * evt, std::string collName){

  streamlog_out(DEBUG) << "CLICPfoSelectorTree::fillTree" << std::endl; 

  LCCollection* col = evt->getCollection( collName ) ;

  streamlog_out( MESSAGE ) << "    Type          PDG    E      Pt  cosTheta #trk time  #Clu  time   ecal  hcal  " << std::endl;

  // this will only be entered if the collection is available
  if( col != NULL){
    int nelem = col->getNumberOfElements();
    streamlog_out(DEBUG) << "Number of PFOs in this event = " << nelem << std::endl;

    // loop on MC/PFO particles
    for(int i=0; i< nelem ; i++){

      ReconstructedParticle* pPfo = static_cast<ReconstructedParticle*>( col->getElementAt( i ) ) ;
      type = pPfo->getType();
      px = pPfo->getMomentum()[0];
      py = pPfo->getMomentum()[1];
      pz = pPfo->getMomentum()[2];
      pT = sqrt(px*px + py*py);
      p  = sqrt(pT*pT + pz*pz);
      costheta = fabs(pz)/p;
      energy   = pPfo->getEnergy();
      mass     = pPfo->getMass();
      charge   = pPfo->getCharge();

      const TrackVec   tracks   = pPfo->getTracks();
      const ClusterVec clusters = pPfo->getClusters();

      //get track time
      float trackTime = std::numeric_limits<float>::max();
      clusterTime = 999.;
      clusterTimeEcal = 999.;
      clusterTimeHcalEndcap = 999.;

//      streamlog_out(DEBUG) << " *** PFO with number of tracks = " << tracks.size() << std::endl;
      for(unsigned int trk = 0; trk < tracks.size(); trk++){
	const Track *track = tracks[trk];
	float tof;
  	const float time = PfoUtil::TimeAtEcal(track,tof);
  	if( fabs(time) < trackTime ){
  	  trackTime = time;
  	}
      }

      //get clusters time
//      streamlog_out(DEBUG) << "PFO with number of clusters = " << clusters.size() << std::endl;
      for(unsigned int clu = 0; clu < clusters.size(); clu++){
	float meanTime(999.);
	float meanTimeEcal(999.);
	float meanTimeHcalEndcap(999.);
	int   nEcal(0);
	int   nHcalEnd(0);
	int   nCaloHitsUsed(0);

	const Cluster *cluster = clusters[clu];
	PfoUtil::GetClusterTimes(cluster,meanTime,nCaloHitsUsed,meanTimeEcal,nEcal,meanTimeHcalEndcap,nHcalEnd,false);

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

        streamlog_out(DEBUG) << "clusterTime: " << clusterTime << std::endl;
      
      }


      std::stringstream output;
      output << std::fixed;
      output << std::setprecision(precision);

     if(clusters.size()==0)
       FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT,costheta,tracks.size(),trackTime,"-","-","-","-");
     if(tracks.size()==0)
       FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT,costheta,"","-",clusters.size(),clusterTime,clusterTimeEcal,clusterTimeHcalEndcap);
     if(tracks.size()>0&&clusters.size()>0)
       FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT,costheta,tracks.size(),trackTime,clusters.size(),clusterTime,clusterTimeEcal,clusterTimeHcalEndcap);

      streamlog_out( MESSAGE ) << output.str();

      m_tree->Fill();

    }

  }

}
