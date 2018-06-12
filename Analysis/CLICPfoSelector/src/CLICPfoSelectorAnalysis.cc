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

LCRelationNavigator* m_reltrue = 0;
LCRelationNavigator* m_trackreltrue = 0;
LCRelationNavigator* m_clureltrue = 0;

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

  registerProcessorParameter("MinECalHitsForHadrons",
			    			             "Min number of Ecal hits to use clusterTime info from Ecal (for neutral and charged hadrons only)",
			     			             minECalHits,
			     			             int(5));

  registerProcessorParameter("MinHcalEndcapHitsForHadrons",
			    			             "Min number of Hcal Endcap hits to use clusterTime info from Hcal Endcap (for neutral and charged hadrons only)",
			     			             minHcalEndcapHits,
			     			             int(5));

  registerProcessorParameter("ForwardCosThetaForHighEnergyNeutralHadrons",
			                       "ForwardCosThetaForHighEnergyNeutralHadrons",
			                       forwardCosThetaForHighEnergyNeutralHadrons,
			                       float(0.95));

  registerProcessorParameter("ForwardHighEnergyNeutralHadronsEnergy",
			                       "ForwardHighEnergyNeutralHadronsEnergy",
			                       forwardHighEnergyNeutralHadronsEnergy,
			                       float(10.00));
  
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

  registerProcessorParameter("AnalyzeAll",
                             "Boolean factor to decide if perform the analysis on all PFOs",
                             analyzeAll,
                             bool(true));

  registerProcessorParameter("AnalyzeSignal",
                             "Boolean factor to decide if perform the analysis only on PFOs belonging to the signal",
                             analyzeSignal,
                             bool(true));

  registerProcessorParameter("AnalyzeOverlay",
                             "Boolean factor to decide if perform the analysis on only on PFOs belonging to the overlay",
                             analyzeOverlay,
                             bool(true));

  registerInputCollection(LCIO::MCPARTICLE,
                          "MCPhysicsParticleCollectionName",
                          "Name of the MCPhysicsParticle input collection",
                          m_inputPhysicsParticleCollection,
                          std::string("MCPhysicsParticles"));

  registerInputCollection(LCIO::LCRELATION,
                          "RecoMCTruthLink",
                          "Name of the RecoMCTruthLink input collection"  ,
                          m_recoMCTruthLink,
                          std::string("RecoMCTruthLink") ) ;

  registerInputCollection(LCIO::LCRELATION,
                          "SiTracksMCTruthLink",
                          "Name of the SiTracksMCTruthLink input collection"  ,
                          m_SiTracksMCTruthLink,
                          std::string("SiTracksMCTruthLink") ) ;

registerInputCollection(LCIO::LCRELATION,
                        "ClusterMCTruthLink",
                        "Name of the ClusterMCTruthLink input collection"  ,
                        m_ClusterMCTruthLink,
                        std::string("ClusterMCTruthLink") ) ;
}



void CLICPfoSelectorAnalysis::init() { 

  streamlog_out(DEBUG) << "   init called  " << endl ;
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  AIDAProcessor::histogramFactory(this);

  //initializing TTree
  pfo_tree = new TTree(treeName.c_str(), treeName.c_str());
  pfo_tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  pfo_tree->Branch("runNumber", &runNumber, "runNumber/I");
  pfo_tree->Branch("nPartPFO", &nPartPFO, "nPartPFO/I");

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
  pfo_tree->Branch("nCaloHits", &nCaloHits, "nCaloHits/I");
  pfo_tree->Branch("nEcalHits", &nEcalHits, "nEcalHits/I");
  pfo_tree->Branch("nHcalEndCapHits", &nHcalEndCapHits, "nHcalEndCapHits/I");

  pfo_tree->Branch("clusterTime", &clusterTime, "clusterTime/D");
  pfo_tree->Branch("clusterTimeEcal", &clusterTimeEcal, "clusterTimeEcal/D");
  pfo_tree->Branch("clusterTimeHcalEndcap", &clusterTimeHcalEndcap, "clusterTimeHcalEndcap/D");

  pfo_tree->Branch("trk_clu_sameMCPart", &trk_clu_sameMCPart, "trk_clu_sameMCPart/I");
  pfo_tree->Branch("atLeastOneSignal", &atLeastOneSignal, "atLeastOneSignal/I");



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
  if(analyzeAll){
    generationCategories.push_back("all");
  }
  if(analyzeSignal){
    generationCategories.push_back("signal");
  }
  if(analyzeOverlay){
    generationCategories.push_back("overlay");
  }

  std::string graphLabel = "";
  for(auto ipc : particleCategories){
    for(auto igc : generationCategories){
      graphLabel = ipc + "_" + igc;
      streamlog_out(DEBUG) << "Analysing followig particle category: " << graphLabel << endl;
      timeVsPt[graphLabel] = new TGraph();
      timeVsPt_barrel[graphLabel] = new TGraph();
      timeVsPt_endcap[graphLabel] = new TGraph();
    }
  }

}

void CLICPfoSelectorAnalysis::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void CLICPfoSelectorAnalysis::processEvent( LCEvent * evt ) { 

  streamlog_out( DEBUG1 ) << "Processing event " << evt->getEventNumber() << " in CLICPfoSelectorAnalysis processor " << endl;

  eventNumber=evt->getEventNumber();
  runNumber=_nRun;

  // Get the collection of MC physics particles (signal)
  // and store them in a std::set
  LCCollection* physicsParticleCollection = NULL;
  try{
    physicsParticleCollection = evt->getCollection( m_inputPhysicsParticleCollection );
  }
  catch( lcio::DataNotAvailableException e )
  {
    streamlog_out(WARNING) << m_inputPhysicsParticleCollection   << " collection not available" << std::endl;
    physicsParticleCollection = NULL;
  }
  for(int ipart = 0; ipart < physicsParticleCollection->getNumberOfElements(); ipart++){
    MCParticle *signal = static_cast<MCParticle*> ( physicsParticleCollection->getElementAt(ipart) );
    physicsParticles.push_back(signal);
  }
  streamlog_out( DEBUG1 ) << physicsParticles.size() << " MC Particles belong to the signal." << endl;
  //for(auto part : physicsParticles){
  //  streamlog_out( DEBUG1 ) << part->id() << endl;
  //}

  //Get MC Particles associated to the PFOs
  LCCollection* rmclcol = NULL;
  try{
    rmclcol = evt->getCollection( m_recoMCTruthLink );
  }
  catch( lcio::DataNotAvailableException e )
  {
    streamlog_out(WARNING) << m_recoMCTruthLink   << " collection not available" << std::endl;
    rmclcol = NULL;
  }
  if( rmclcol != NULL ){
    m_reltrue = new LCRelationNavigator( rmclcol );
  }

  //Get MC Particles associated to the track of the PFOs
  LCCollection* rtrkclcol = NULL;
  try{
    rtrkclcol = evt->getCollection( m_SiTracksMCTruthLink );
  }
  catch( lcio::DataNotAvailableException e )
  {
    streamlog_out(WARNING) << m_SiTracksMCTruthLink   << " collection not available" << std::endl;
    rtrkclcol = NULL;
  }
  if( rtrkclcol != NULL ){
    m_trackreltrue = new LCRelationNavigator( rtrkclcol );
  }
    
  //Get MC Particles associated to this cluster
  LCCollection* rclulcol = NULL;
  try{
    rclulcol = evt->getCollection( m_ClusterMCTruthLink );
  }
  catch( lcio::DataNotAvailableException e )
  {
    streamlog_out(WARNING) << m_ClusterMCTruthLink   << " collection not available" << std::endl;
    rclulcol = NULL;
  }
  if( rclulcol != NULL ){
    m_clureltrue = new LCRelationNavigator( rclulcol );
  }

  fillTree(evt, colNamePFOs);

  fillScatterPlots();

  streamlog_out( DEBUG1 ) << "   processing event: " << evt->getEventNumber() 
      << "   in run:  " << _nRun << endl ;
  
  physicsParticles.clear();
  _nEvt ++ ;
}

void CLICPfoSelectorAnalysis::end(){ 
  streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis::end()  " << name() 
    << " processed " << _nEvt << " events in " << _nRun << " runs " << endl;

  //filling TGraphs for each category
  AIDAProcessor::histogramFactory(this);
  std::string graphLabel = "";
  for(auto ipc : particleCategories){
    for(auto igc : generationCategories){
      graphLabel = ipc + "_" + igc;
      streamlog_out(DEBUG) << "Analysing followig particle category: " << graphLabel << endl;
      timeVsPt[graphLabel]->Write((graphLabel + "_timeVsPt").c_str());
      timeVsPt_barrel[graphLabel]->Write((graphLabel + "_timeVsPt_central").c_str());
      timeVsPt_endcap[graphLabel]->Write((graphLabel + "_timeVsPt_forward").c_str());
    }
  }

}

void CLICPfoSelectorAnalysis::fillTree(LCEvent * evt, string collName){

  //streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis::fillTree" << endl; 

  LCCollection* col = evt->getCollection( collName ) ;

  streamlog_out( DEBUG1 ) << "    Type          PDG    E      Pt  cosTheta #trk time  #Clu  time   ecal  hcal  " << endl;

  // this will only be entered if the collection is available
  if( col != NULL){
    int nelem = col->getNumberOfElements();
    streamlog_out(DEBUG) << nelem << " PFOs in this event" << endl;

    // loop on PFO particles
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

      //MC linker
      LCObjectVec mcvec;
      mcvec = m_reltrue->getRelatedToObjects( pPfo );
      atLeastOneSignal = 0;
      streamlog_out( DEBUG4 ) << "Pfos Type = " << pPfo->getType() << std::endl;
      streamlog_out( DEBUG4 ) << "Pfos Pt = " << pT << std::endl;
      streamlog_out( DEBUG4 ) << "Pfos Energy = " << pPfo->getEnergy() << std::endl;

      streamlog_out( DEBUG4 ) << "mcvec size = " << mcvec.size() << std::endl;
      if ( mcvec.size() > 0 ) { // if the PFO is associated to MC Particle
        for ( auto imcp : mcvec ) { //  MC Particles loop
          MCParticle* mc_part  = dynamic_cast<MCParticle*>(imcp);
          streamlog_out( DEBUG4 ) << "MC Part id = " << mc_part->id() << std::endl;
          streamlog_out( DEBUG4 ) << "MC Part Type = " << mc_part->getPDG() << std::endl;
          streamlog_out( DEBUG4 ) << "MC Part Pt = " << sqrt(mc_part->getMomentum()[0]*mc_part->getMomentum()[0] + mc_part->getMomentum()[1]*mc_part->getMomentum()[1]) << std::endl;
          streamlog_out( DEBUG4 ) << "MC Part energy = " << mc_part->getEnergy() << std::endl;
          
          if(find(physicsParticles.begin(), physicsParticles.end(), mc_part) != physicsParticles.end() && atLeastOneSignal == 0) {
             streamlog_out( DEBUG4 ) << "MC Part is the first signal" << std::endl;
             atLeastOneSignal = 1;
          } else { 
            streamlog_out( DEBUG4 ) << "MC Part is overlay" << std::endl;
          }

        }
      }


      const TrackVec   tracks   = pPfo->getTracks();
      const ClusterVec clusters = pPfo->getClusters();
      nTracks = tracks.size();
      nClusters = clusters.size();

      //get track time
      float trackTime = numeric_limits<float>::max();
      clusterTime = 999.;
      clusterTimeEcal = 999.;
      clusterTimeHcalEndcap = 999.;

      //streamlog_out(DEBUG1) << " PFO with number of tracks = " << tracks.size() << endl;
      for(unsigned int trk = 0; trk < tracks.size(); trk++){
	const Track *track = tracks[trk];
	float tof;
  	const float time = PfoUtil::TimeAtEcal(track,tof);
  	if( fabs(time) < trackTime ){
  	  trackTime = time;
  	}
      }
      streamlog_out(DEBUG1) << " trackTime = " << trackTime << endl;

      //get clusters time
      //streamlog_out(DEBUG1) << " PFO with number of clusters = " << clusters.size() << endl;
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

      trk_clu_sameMCPart = 0;
      for(unsigned int it = 0; it < tracks.size(); it++){
        Track *track = tracks[it];

        //MC linker to tracks
        LCObjectVec mctrkvec;
        mctrkvec = m_trackreltrue->getRelatedToObjects( track );
        if ( mctrkvec.size() > 0 ) {
          streamlog_out( DEBUG ) << " Track is associated to " << mctrkvec.size() << " MC particles." << std::endl;
          for(unsigned int ic = 0; ic< clusters.size(); ic++){
            Cluster *cluster = clusters[ic];

            //MC linker to clusters
            LCObjectVec mccluvec;
            mccluvec = m_clureltrue->getRelatedToObjects( cluster );
            if ( mccluvec.size() > 0 ) {
              streamlog_out( DEBUG ) << " Cluster is also associated to " << mccluvec.size() << " MC particles." << std::endl;

              for ( auto imctrk : mctrkvec ) { //  loop on MC Particles associated with track x
                streamlog_out( DEBUG ) << " Running on new mc_part_trk " << std::endl;
                if(trk_clu_sameMCPart == 1)
                  break;
                for ( auto imcclu : mccluvec ){
                  streamlog_out( DEBUG ) << " Running on new mc_part_clu " << std::endl;
                  if(trk_clu_sameMCPart == 0 && imcclu->id() == imctrk->id()){
                    streamlog_out( DEBUG ) << " \t SAME MC PARTICLE!" << std::endl;
                    trk_clu_sameMCPart = 1;
                    break;
                  }
                }
              }
            }
          }
        }
      }

      if(trk_clu_sameMCPart == 0)
        streamlog_out( DEBUG ) << " \t NOT same MC particle!" << std::endl;

      std::cout << " " << std::endl;

  
      stringstream output;
      output << fixed;
      output << setprecision(precision);

      if(clusterTime<-20){
      //if((clusterTime>2 || clusterTimeEcal>2) && costheta > 0.95 && tracks.size()==0 && type!=22){
       output << " *** Interesting PFO: ";
      } 
      if(clusters.size()==0)
        FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT,costheta,tracks.size(),trackTime,"-","-","-","-");
      if(tracks.size()==0)
         FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT,costheta,"","-",clusters.size(),clusterTime,clusterTimeEcal,clusterTimeHcalEndcap);
      if(tracks.size()>0&&clusters.size()>0)
         FORMATTED_OUTPUT_TRACK_CLUSTER(output,type,energy,pT,costheta,tracks.size(),trackTime,clusters.size(),clusterTime,clusterTimeEcal,clusterTimeHcalEndcap);
        
      streamlog_out( DEBUG ) << output.str();
      pfo_tree->Fill();

      }

  }

}

void CLICPfoSelectorAnalysis::fillScatterPlots(std::string signalOrBkg){

  streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis::fillScatterPlots" << endl; 
  int tot_charg = 0, tot_charg_over = 0, tot_charg_sig = 0;

  Int_t nEntries = pfo_tree->GetEntries();
  streamlog_out(DEBUG) << "Reading TTree with nEntries: " << nEntries << endl;

  for (int ie = 0; ie < nEntries; ie++){
    pfo_tree->GetEntry(ie);
    double currentClusterTime = clusterTime;
    bool useHcalTimingOnly = ( (costheta > forwardCosThetaForHighEnergyNeutralHadrons) && (type == 2112) && (energy > forwardHighEnergyNeutralHadronsEnergy));

    //in the case of photons, the time of ECAL clusters are used
    if( analyzePhotons==true && type==22 ){
      
      if(std::find(particleCategories.begin(), particleCategories.end(), "photons") != particleCategories.end() &&
         std::find(generationCategories.begin(), generationCategories.end(), "all") != generationCategories.end()){

      	//in the case of photons, the clusterTimeEcal is used
      	currentClusterTime = clusterTimeEcal;
      	streamlog_out(DEBUG2) << "Filling scatter plot for photon with pT and current clusterTime: " << pT << "," << currentClusterTime << endl;

        timeVsPt["photons_all"]->SetPoint(ie,pT,currentClusterTime);
        if(costheta < cutCosTheta)
          timeVsPt_barrel["photons_all"]->SetPoint(ie,pT,currentClusterTime);
        else
          timeVsPt_endcap["photons_all"]->SetPoint(ie,pT,currentClusterTime);

        if(analyzeSignal && atLeastOneSignal &&
           std::find(generationCategories.begin(), generationCategories.end(), "signal") != generationCategories.end()){
          timeVsPt["photons_signal"]->SetPoint(ie,pT,currentClusterTime);
          if(costheta < cutCosTheta)
            timeVsPt_barrel["photons_signal"]->SetPoint(ie,pT,currentClusterTime);
          else
            timeVsPt_endcap["photons_signal"]->SetPoint(ie,pT,currentClusterTime);
        }

        if(analyzeOverlay && !atLeastOneSignal &&
           std::find(generationCategories.begin(), generationCategories.end(), "overlay") != generationCategories.end()){
          timeVsPt["photons_overlay"]->SetPoint(ie,pT,currentClusterTime);
          if(costheta < cutCosTheta)
            timeVsPt_barrel["photons_overlay"]->SetPoint(ie,pT,currentClusterTime);
          else
            timeVsPt_endcap["photons_overlay"]->SetPoint(ie,pT,currentClusterTime);
        }

      } else {
        streamlog_out(ERROR) << "Cannot fill scatter plots because TGraph for photons was not created. " << endl;
        exit(0);
      }
    }

    if( analyzeNeutralHadrons==true && type!=22 && charge==0 ){

      if(std::find(particleCategories.begin(), particleCategories.end(), "neutralHadrons") != particleCategories.end() &&
         std::find(generationCategories.begin(), generationCategories.end(), "all") != generationCategories.end()){
        //streamlog_out( MESSAGE ) << "Has nEcalHits = " << nEcalHits << ", nHcalEndCapHits = " << nHcalEndCapHits << ", nCaloHits/2 = " << nCaloHits/2. << ", "<< std::endl;
        //streamlog_out( MESSAGE ) << "Has clusterTimeEcal = " << clusterTimeEcal << ", clusterTimeHcalEndcap = " << clusterTimeHcalEndcap << std::endl;

        //in the case the nEcalHits is more than expected, the time computed Ecal is used
        if(!useHcalTimingOnly && ( nEcalHits > minECalHits || nEcalHits >= nCaloHits/2.) ){
          currentClusterTime = clusterTimeEcal;
        //in the case the nHcalEndCapHits is more than expected, the time computed Hcal endcap is used
        } else if ( (nHcalEndCapHits >= minHcalEndcapHits) || (nHcalEndCapHits >= nCaloHits/2.) ) {
        	currentClusterTime = clusterTimeHcalEndcap;
        }
        streamlog_out(DEBUG2) << "Filling scatter plot for neutralHadrons with pT and clusterTime: " << pT << "," << currentClusterTime << endl;

        timeVsPt["neutralHadrons_all"]->SetPoint(ie,pT,currentClusterTime);
        if(costheta < cutCosTheta)
          timeVsPt_barrel["neutralHadrons_all"]->SetPoint(ie,pT,currentClusterTime);
        else
          timeVsPt_endcap["neutralHadrons_all"]->SetPoint(ie,pT,currentClusterTime);

        if(analyzeSignal && atLeastOneSignal &&
           std::find(generationCategories.begin(), generationCategories.end(), "signal") != generationCategories.end()){
          timeVsPt["neutralHadrons_signal"]->SetPoint(ie,pT,currentClusterTime);
          if(costheta < cutCosTheta)
            timeVsPt_barrel["neutralHadrons_signal"]->SetPoint(ie,pT,currentClusterTime);
          else
            timeVsPt_endcap["neutralHadrons_signal"]->SetPoint(ie,pT,currentClusterTime);
        }

        if(analyzeOverlay && !atLeastOneSignal &&
           std::find(generationCategories.begin(), generationCategories.end(), "overlay") != generationCategories.end()){
          timeVsPt["neutralHadrons_overlay"]->SetPoint(ie,pT,currentClusterTime);
          if(costheta < cutCosTheta)
            timeVsPt_barrel["neutralHadrons_overlay"]->SetPoint(ie,pT,currentClusterTime);
          else
            timeVsPt_endcap["neutralHadrons_overlay"]->SetPoint(ie,pT,currentClusterTime);
        }

          
      } else {
        streamlog_out(ERROR) << "Cannot fill scatter plots because TGraph for neutralHadrons was not created. " << endl;
        exit(0);
      }
    
    }
    
    if( analyzeChargedPfos==true && charge!=0 ){
      if(std::find(particleCategories.begin(), particleCategories.end(), "chargedPfos") != particleCategories.end() &&
         std::find(generationCategories.begin(), generationCategories.end(), "all") != generationCategories.end()){

        //in the case the nEcalHits is more than expected, the time computed Ecal is used
        if((!useHcalTimingOnly && ( nEcalHits > minECalHits || nEcalHits >= nCaloHits/2.) )){
          currentClusterTime = clusterTimeEcal;
        //in the case the nHcalEndCapHits is more than expected, the time computed Hcal endcap is used
        } else if ( (nHcalEndCapHits >= minHcalEndcapHits) || (nHcalEndCapHits >= nCaloHits/2.) ) {
        	currentClusterTime = clusterTimeHcalEndcap;
        }
        streamlog_out(DEBUG2) << "Filling scatter plot for chargedPfos with pT and clusterTime: " << pT << "," << currentClusterTime << endl;

        timeVsPt["chargedPfos_all"]->SetPoint(ie,pT,currentClusterTime); tot_charg++;
        if(costheta < cutCosTheta)
          timeVsPt_barrel["chargedPfos_all"]->SetPoint(ie,pT,currentClusterTime);
        else
          timeVsPt_endcap["chargedPfos_all"]->SetPoint(ie,pT,currentClusterTime);

        if(analyzeSignal && atLeastOneSignal &&
           std::find(generationCategories.begin(), generationCategories.end(), "signal") != generationCategories.end()){
          timeVsPt["chargedPfos_signal"]->SetPoint(ie,pT,currentClusterTime); tot_charg_sig++;
          if(costheta < cutCosTheta)
            timeVsPt_barrel["chargedPfos_signal"]->SetPoint(ie,pT,currentClusterTime);
          else
            timeVsPt_endcap["chargedPfos_signal"]->SetPoint(ie,pT,currentClusterTime);
        }

        if(analyzeOverlay && !atLeastOneSignal &&
           std::find(generationCategories.begin(), generationCategories.end(), "overlay") != generationCategories.end()){
          timeVsPt["chargedPfos_overlay"]->SetPoint(ie,pT,currentClusterTime); tot_charg_over++;
          if(costheta < cutCosTheta)
            timeVsPt_barrel["chargedPfos_overlay"]->SetPoint(ie,pT,currentClusterTime);
          else
            timeVsPt_endcap["chargedPfos_overlay"]->SetPoint(ie,pT,currentClusterTime);
        }

      } else {
        streamlog_out(ERROR) << "Cannot fill scatter plots because TGraph for chargedPfos was not created. " << endl;
        exit(0);
      }
    }

  }

streamlog_out(DEBUG) << "Tot charged Pfos (all): " << tot_charg << endl;
streamlog_out(DEBUG) << "Tot charged Pfos (sig): " << tot_charg_sig << endl;
streamlog_out(DEBUG) << "Tot charged Pfos (bkg): " << tot_charg_over << endl;

}
