#include "CLICPfoSelectorAnalysis.h"
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

CLICPfoSelectorAnalysis aCLICPfoSelectorAnalysis ;


CLICPfoSelectorAnalysis::CLICPfoSelectorAnalysis() : Processor("CLICPfoSelectorAnalysis") {

  _description = "CLICPfoSelectorAnalysis produces a tree and scatter plots using PFO variables defined in the CLICPfoSelector." ;

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                          "PFOCollectionName" , 
                          "Name of the PFO collection",
                          colNamePFOs,
                          string("PandoraPFOs")
  );

  registerProcessorParameter("TreeName",
                             "Name of output tree",
                             treeName,
                             string("PfoTree") );

  registerProcessorParameter("CosThetaCut",
                             "Cut on the PFO cosTheta to define central/forward region",
                             cutCosTheta,
                             float(0.975));

  registerProcessorParameter("AnalyzePhotons",
                             "Boolean factor to decide if perform the analysis on photons",
                             analyzePhotons,
                             bool(true));

  registerProcessorParameter("AnalyzeChargedPfos",
                             "Boolean factor to decide if perform the analysis on charged PFOs",
                             analyzeChargedPfos,
                             bool(true));

  registerProcessorParameter("AnalyzeNeutralHadrons",
                             "Boolean factor to decide if perform the analysis on neutral hadrons",
                             analyzeNeutralHadrons,
                             bool(true));
}



void CLICPfoSelectorAnalysis::init() { 

  streamlog_out(DEBUG) << "   init called  " << endl ;
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  AIDAProcessor::histogramFactory(this);

  //initializing TTree
  pfo_tree = new TTree(treeName.c_str(), treeName.c_str());
  pfo_tree->Branch("type", &type, "type/I");
  pfo_tree->Branch("p", &p, "p/D");
  pfo_tree->Branch("px", &px, "px/D");
  pfo_tree->Branch("py", &py, "py/D");
  pfo_tree->Branch("pz", &pz, "pz/D");
  pfo_tree->Branch("pT", &pT, "pT/D");

  pfo_tree->Branch("costheta", &costheta, "costhetaMC/D");
  pfo_tree->Branch("energy", &energy, "energy/D");
  pfo_tree->Branch("mass", &mass, "mass/D");
  pfo_tree->Branch("charge", &charge, "charge/D");
  pfo_tree->Branch("nTracks", &nTracks, "nTracks/I");
  pfo_tree->Branch("nClusters", &nClusters, "nClusters/I");

  pfo_tree->Branch("clusterTime", &clusterTime, "clusterTime/D");
  pfo_tree->Branch("clusterTimeEcal", &clusterTimeEcal, "clusterTimeEcal/D");
  pfo_tree->Branch("clusterTimeHcalEndcap", &clusterTimeHcalEndcap, "clusterTimeHcalEndcap/D");

  pfo_tree->Branch("nCaloHits", &nCaloHits, "nCaloHits/I");
  pfo_tree->Branch("nEcalHits", &nEcalHits, "nEcalHits/I");
  pfo_tree->Branch("nHcalEndCapHits", &nHcalEndCapHits, "nHcalEndCapHits/I");

  pfo_tree->Branch("eventNumber", &eventNumber, "eventNumber/I");

  pfo_tree->Branch("nPartPFO", &nPartPFO, "nPartPFO/I");
  pfo_tree->Branch("nPartMC", &nPartMC, "nPartMC/I");

  streamlog_out(DEBUG) << "   set up ttree  " << endl ;

  //initializing TGraphs
  if(analyzePhotons){
    particleCategories.push_back("photons");
  }
  if(analyzeChargedPfos){
    particleCategories.push_back("chargedPfos");
  }
  if(analyzeNeutralHadrons){
    particleCategories.push_back("neutralHadrons");
  }

  for(vector<string>::iterator it = particleCategories.begin(); it != particleCategories.end(); ++it){
    streamlog_out(DEBUG) << "Analysing followig particle category: " << *it << endl;
    timeVsPt[*it] = new TGraph();
    timeVsPt_barrel[*it] = new TGraph();
    timeVsPt_endcap[*it] = new TGraph();
  }

}

void CLICPfoSelectorAnalysis::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void CLICPfoSelectorAnalysis::processEvent( LCEvent * evt ) { 

  streamlog_out( DEBUG ) << "Processing event " << evt->getEventNumber() << " in CLICPfoSelectorAnalysis processor " << endl;

  eventNumber=evt->getEventNumber();

  fillTree(evt, colNamePFOs);

  fillScatterPlots();

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
      << "   in run:  " << evt->getRunNumber() << endl ;

  _nEvt ++ ;
}

void CLICPfoSelectorAnalysis::end(){ 
  streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis::end()  " << name() 
    << " processed " << _nEvt << " events in " << _nRun << " runs " << endl;

  //filling TGraphs for each category
  AIDAProcessor::histogramFactory(this);
  for(vector<string>::iterator it = particleCategories.begin(); it != particleCategories.end(); ++it){
    streamlog_out(DEBUG) << "Analysing followig particle category: " << *it << endl;
    timeVsPt[*it]->Write((*it + "_timeVsPt").c_str());
    timeVsPt_barrel[*it]->Write((*it + "_timeVsPt_central").c_str());
    timeVsPt_endcap[*it]->Write((*it + "_timeVsPt_forward").c_str());
  }
}

void CLICPfoSelectorAnalysis::fillTree(LCEvent * evt, string collName){

  streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis::fillTree" << endl; 

  LCCollection* col = evt->getCollection( collName ) ;

  streamlog_out( MESSAGE ) << "    Type          PDG    E      Pt  cosTheta #trk time  #Clu  time   ecal  hcal  " << endl;

  // this will only be entered if the collection is available
  if( col != NULL){
    int nelem = col->getNumberOfElements();
    streamlog_out(DEBUG) << "Number of PFOs in this event = " << nelem << endl;

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
      nTracks = tracks.size();
      nClusters = clusters.size();

      //get track time
      float trackTime = numeric_limits<float>::max();
      clusterTime = 999.;
      clusterTimeEcal = 999.;
      clusterTimeHcalEndcap = 999.;

//      streamlog_out(DEBUG) << " *** PFO with number of tracks = " << tracks.size() << endl;
      for(unsigned int trk = 0; trk < tracks.size(); trk++){
	const Track *track = tracks[trk];
	float tof;
  	const float time = PfoUtil::TimeAtEcal(track,tof);
  	if( fabs(time) < trackTime ){
  	  trackTime = time;
  	}
      }

      //get clusters time
//      streamlog_out(DEBUG) << "PFO with number of clusters = " << clusters.size() << endl;
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

      }


      stringstream output;
      output << fixed;
      output << setprecision(precision);

     if(clusters.size()==0)
       FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT,costheta,tracks.size(),trackTime,"-","-","-","-");
     if(tracks.size()==0)
       FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT,costheta,"","-",clusters.size(),clusterTime,clusterTimeEcal,clusterTimeHcalEndcap);
     if(tracks.size()>0&&clusters.size()>0)
       FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT,costheta,tracks.size(),trackTime,clusters.size(),clusterTime,clusterTimeEcal,clusterTimeHcalEndcap);

      streamlog_out( MESSAGE ) << output.str();

      pfo_tree->Fill();

    }

  }

}

void CLICPfoSelectorAnalysis::fillScatterPlots(){

  streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis::fillScatterPlots" << endl; 

  Int_t nEntries = pfo_tree->GetEntries();
  streamlog_out(DEBUG) << "Reading TTree with nEntries: " << nEntries << endl;

  for (int ie = 0; ie < nEntries; ie++){
    pfo_tree->GetEntry(ie);
    streamlog_out(DEBUG) << "pfo pT: " << pT << endl;

    //in the case of photons, the time of ECAL clusters are used
    if( analyzePhotons==true && type==22 ){
      if(std::find(particleCategories.begin(), particleCategories.end(), "photons") != particleCategories.end()){
        streamlog_out(DEBUG) << "Filling scatter plot for photon " << endl;
        timeVsPt["photons"]->SetPoint(ie,clusterTime,pT);
        if(costheta < cutCosTheta)
          timeVsPt_barrel["photons"]->SetPoint(ie,clusterTime,pT);
        else
          timeVsPt_endcap["photons"]->SetPoint(ie,clusterTime,pT);
      } else {
        streamlog_out(ERROR) << "Cannot fill scatter plots because TGraph for photons was not created. " << endl;
        exit(0);
      }
    }

    if( analyzeNeutralHadrons==true && type!=22 && charge==0 ){
      if(std::find(particleCategories.begin(), particleCategories.end(), "neutralHadrons") != particleCategories.end()){
        streamlog_out(DEBUG) << "Filling scatter plot for neutralHadrons " << endl;
        timeVsPt["neutralHadrons"]->SetPoint(ie,clusterTime,pT);
        if(costheta < cutCosTheta)
          timeVsPt_barrel["neutralHadrons"]->SetPoint(ie,clusterTime,pT);
        else
          timeVsPt_endcap["neutralHadrons"]->SetPoint(ie,clusterTime,pT);
      } else {
        streamlog_out(ERROR) << "Cannot fill scatter plots because TGraph for neutralHadrons was not created. " << endl;
        exit(0);
      }
    }

    if( analyzeChargedPfos==true && charge!=0 ){
      if(std::find(particleCategories.begin(), particleCategories.end(), "chargedPfos") != particleCategories.end()){
        streamlog_out(DEBUG) << "Filling scatter plot for chargedPfos " << endl;
        timeVsPt["chargedPfos"]->SetPoint(ie,clusterTime,pT);
        if(costheta < cutCosTheta)
          timeVsPt_barrel["chargedPfos"]->SetPoint(ie,clusterTime,pT);
        else
          timeVsPt_endcap["chargedPfos"]->SetPoint(ie,clusterTime,pT);
      } else {
        streamlog_out(ERROR) << "Cannot fill scatter plots because TGraph for chargedPfos was not created. " << endl;
        exit(0);
      }
    }

  }


}
