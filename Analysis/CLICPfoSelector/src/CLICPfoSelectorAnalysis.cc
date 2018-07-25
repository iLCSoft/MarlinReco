#include "CLICPfoSelectorAnalysis.h"
#include <iostream>

#include <EVENT/LCCollection.h>

#include "marlin/VerbosityLevels.h"

#include <marlin/AIDAProcessor.h>

using namespace lcio;
using namespace marlin;

const int precision  = 2;
const int widthFloat = 7;
const int widthInt   = 5;

LCRelationNavigator* m_reltrue      = 0;
LCRelationNavigator* m_trackreltrue = 0;
LCRelationNavigator* m_clureltrue   = 0;

CLICPfoSelectorAnalysis aCLICPfoSelectorAnalysis;

CLICPfoSelectorAnalysis::CLICPfoSelectorAnalysis() : Processor("CLICPfoSelectorAnalysis") {
  _description =
      "CLICPfoSelectorAnalysis produces a tree and scatter plots using PFO variables defined in the CLICPfoSelector.";

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "PFOCollectionName", "Name of the PFO collection", colNamePFOs,
                          string("PandoraPFOs"));

  registerProcessorParameter("TreeName", "Name of output tree", treeName, string("PfoTree"));

  registerProcessorParameter("CosThetaCut", "Cut on the PFO cosTheta to define central/forward region", cutCosTheta,
                             float(0.975));

  registerProcessorParameter(
      "MinECalHitsForHadrons",
      "Min number of Ecal hits to use clusterTime info from Ecal (for neutral and charged hadrons only)", minECalHits,
      int(5));

  registerProcessorParameter(
      "MinHcalEndcapHitsForHadrons",
      "Min number of Hcal Endcap hits to use clusterTime info from Hcal Endcap (for neutral and charged hadrons only)",
      minHcalEndcapHits, int(5));

  registerProcessorParameter("ForwardCosThetaForHighEnergyNeutralHadrons", "ForwardCosThetaForHighEnergyNeutralHadrons",
                             forwardCosThetaForHighEnergyNeutralHadrons, float(0.95));

  registerProcessorParameter("ForwardHighEnergyNeutralHadronsEnergy", "ForwardHighEnergyNeutralHadronsEnergy",
                             forwardHighEnergyNeutralHadronsEnergy, float(10.00));

  registerProcessorParameter("AnalyzePhotons", "Boolean factor to decide if perform the analysis on photons", analyzePhotons,
                             bool(true));

  registerProcessorParameter("AnalyzeChargedPfos", "Boolean factor to decide if perform the analysis on charged PFOs",
                             analyzeChargedPfos, bool(true));

  registerProcessorParameter("AnalyzeNeutralHadrons", "Boolean factor to decide if perform the analysis on neutral hadrons",
                             analyzeNeutralHadrons, bool(true));

  registerProcessorParameter("AnalyzeAll", "Boolean factor to decide if perform the analysis on all PFOs", analyzeAll,
                             bool(true));

  registerProcessorParameter("AnalyzeSignal",
                             "Boolean factor to decide if perform the analysis only on PFOs belonging to the signal",
                             analyzeSignal, bool(true));

  registerProcessorParameter("AnalyzeOverlay",
                             "Boolean factor to decide if perform the analysis on only on PFOs belonging to the overlay",
                             analyzeOverlay, bool(true));

  registerProcessorParameter("MinEnergy", "Minimum energy (needed for the energy histos)", en_min, float(0));

  registerProcessorParameter("MaxEnergy", "Maximum energy (needed for the energy histos)", en_max, float(500));

  registerInputCollection(LCIO::MCPARTICLE, "MCPhysicsParticleCollectionName",
                          "Name of the MCPhysicsParticle input collection", m_inputPhysicsParticleCollection,
                          std::string("MCPhysicsParticles"));

  registerInputCollection(LCIO::LCRELATION, "RecoMCTruthLink", "Name of the RecoMCTruthLink input collection",
                          m_recoMCTruthLink, std::string("RecoMCTruthLink"));

  registerInputCollection(LCIO::LCRELATION, "SiTracksMCTruthLink", "Name of the SiTracksMCTruthLink input collection",
                          m_SiTracksMCTruthLink, std::string("SiTracksMCTruthLink"));

  registerInputCollection(LCIO::LCRELATION, "ClusterMCTruthLink", "Name of the ClusterMCTruthLink input collection",
                          m_ClusterMCTruthLink, std::string("ClusterMCTruthLink"));
}

void CLICPfoSelectorAnalysis::init() {
  printParameters();

  _nRun = -1;
  _nEvt = 0;

  AIDAProcessor::histogramFactory(this);

  //initializing TTree
  pfo_tree = new TTree(treeName.c_str(), treeName.c_str());
  pfo_tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  pfo_tree->Branch("runNumber", &runNumber, "runNumber/I");

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

  streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis ---> set up ttree " << std::endl;

  //initializing TGraphs and TH1Fs
  if (analyzePhotons) {
    particleCategories.push_back("photons");
  }
  if (analyzeChargedPfos) {
    particleCategories.push_back("chargedPfos");
  }
  if (analyzeNeutralHadrons) {
    particleCategories.push_back("neutralHadrons");
  }
  if (analyzeAll) {
    generationCategories.push_back("all");
  }
  if (analyzeSignal) {
    generationCategories.push_back("signal");
  }
  if (analyzeOverlay) {
    generationCategories.push_back("overlay");
  }

  h_energy_tot            = new TH1F("h_en_tot_all", "h_en_tot_all", 1000, en_min, en_max);
  h_energy_tot_signal     = new TH1F("h_en_tot_signal", "h_en_tot_signal", 1000, en_min, en_max);
  h_energy_tot_background = new TH1F("h_en_tot_overlay", "h_en_tot_overlay", 1000, en_min, en_max);

  std::string graphLabel = "";
  for (auto ipc : particleCategories) {
    for (auto igc : generationCategories) {
      graphLabel = ipc + "_" + igc;
      streamlog_out(DEBUG) << "Analysing followig particle category: " << graphLabel << endl;
      g_timeVsPt[graphLabel]         = new TGraph();
      g_timeVsPt_central[graphLabel] = new TGraph();
      g_timeVsPt_forward[graphLabel] = new TGraph();
      h_energy[graphLabel] =
          new TH1F((graphLabel + "_energy").c_str(), (graphLabel + "_energy").c_str(), 1000, en_min, en_max);
      h_energy_central[graphLabel] =
          new TH1F((graphLabel + "_energy_central").c_str(), (graphLabel + "_energy_central").c_str(), 1000, en_min, en_max);
      h_energy_forward[graphLabel] =
          new TH1F((graphLabel + "_energy_forward").c_str(), (graphLabel + "_energy_forward").c_str(), 1000, en_min, en_max);
    }
  }
}

void CLICPfoSelectorAnalysis::processRunHeader(LCRunHeader*) {
  _nRun++;
  _nEvt = 0;

  streamlog_out(DEBUG) << std::endl;
  streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis ---> new run : run number = " << _nRun << std::endl;
}

void CLICPfoSelectorAnalysis::processEvent(LCEvent* evt) {
  streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis ---> processing run = " << _nRun << " event = " << _nEvt << std::endl;

  eventNumber = _nEvt;
  runNumber   = _nRun;

  // Get the collection of MC physics particles (signal) and store them in a vector
  LCCollection* physicsParticleCollection = NULL;
  try {
    physicsParticleCollection = evt->getCollection(m_inputPhysicsParticleCollection);
  } catch (lcio::DataNotAvailableException e) {
    streamlog_out(WARNING) << m_inputPhysicsParticleCollection << " collection not available" << std::endl;
    physicsParticleCollection = NULL;
  }
  for (int ipart = 0; ipart < physicsParticleCollection->getNumberOfElements(); ipart++) {
    MCParticle* signal = static_cast<MCParticle*>(physicsParticleCollection->getElementAt(ipart));
    physicsParticles.push_back(signal);
  }
  streamlog_out(DEBUG) << physicsParticles.size() << " MC Particles belong to the signal." << endl;

  //Get MC Particles associated to the PFOs
  LCCollection* rmclcol = NULL;
  try {
    rmclcol = evt->getCollection(m_recoMCTruthLink);
  } catch (lcio::DataNotAvailableException e) {
    streamlog_out(WARNING) << m_recoMCTruthLink << " collection not available" << std::endl;
    rmclcol = NULL;
  }
  if (rmclcol != NULL) {
    m_reltrue = new LCRelationNavigator(rmclcol);
  }

  //Get MC Particles associated to the track of the PFOs
  LCCollection* rtrkclcol = NULL;
  try {
    rtrkclcol = evt->getCollection(m_SiTracksMCTruthLink);
  } catch (lcio::DataNotAvailableException e) {
    streamlog_out(WARNING) << m_SiTracksMCTruthLink << " collection not available" << std::endl;
    rtrkclcol = NULL;
  }
  if (rtrkclcol != NULL) {
    m_trackreltrue = new LCRelationNavigator(rtrkclcol);
  }

  //Get MC Particles associated to this cluster
  LCCollection* rclulcol = NULL;
  try {
    rclulcol = evt->getCollection(m_ClusterMCTruthLink);
  } catch (lcio::DataNotAvailableException e) {
    streamlog_out(WARNING) << m_ClusterMCTruthLink << " collection not available" << std::endl;
    rclulcol = NULL;
  }
  if (rclulcol != NULL) {
    m_clureltrue = new LCRelationNavigator(rclulcol);
  }

  fillTree(evt, colNamePFOs);

  streamlog_out(DEBUG) << " finished event: " << evt->getEventNumber() << "   in run:  " << _nRun << endl;

  physicsParticles.clear();
  _nEvt++;
}

void CLICPfoSelectorAnalysis::end() {
  fillPlots();

  //writing TGraphs for each category
  AIDAProcessor::histogramFactory(this);
  std::string graphLabel = "";
  for (auto ipc : particleCategories) {
    for (auto igc : generationCategories) {
      graphLabel = ipc + "_" + igc;
      g_timeVsPt[graphLabel]->Write((graphLabel + "_timeVsPt").c_str());
      g_timeVsPt_central[graphLabel]->Write((graphLabel + "_timeVsPt_central").c_str());
      g_timeVsPt_forward[graphLabel]->Write((graphLabel + "_timeVsPt_forward").c_str());
    }
  }

  streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis ---> finished " << name() << " processed " << _nEvt << " events in "
                       << _nRun << " runs " << endl;
}

void CLICPfoSelectorAnalysis::fillTree(LCEvent* evt, string collName) {
  streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis ---> filling TTree " << std::endl;

  LCCollection* col = evt->getCollection(collName);

  double en_tot            = 0.0;
  double en_tot_signal     = 0.0;
  double en_tot_background = 0.0;

  std::string graphLabel = "";
  for (auto ipc : particleCategories) {
    for (auto igc : generationCategories) {
      graphLabel                     = ipc + "_" + igc;
      energy_tot[graphLabel]         = 0.0;
      energy_tot_central[graphLabel] = 0.0;
      energy_tot_forward[graphLabel] = 0.0;
    }
  }

  // this will only be entered if the collection is available
  if (col != NULL) {
    int nelem = col->getNumberOfElements();
    streamlog_out(DEBUG) << "Number of Input Pfos = " << nelem << endl;

    // loop on PFO particles
    for (int i = 0; i < nelem; i++) {
      streamlog_out(DEBUG) << " --- " << std::endl;

      ReconstructedParticle* pPfo = static_cast<ReconstructedParticle*>(col->getElementAt(i));
      type                        = pPfo->getType();
      px                          = pPfo->getMomentum()[0];
      py                          = pPfo->getMomentum()[1];
      pz                          = pPfo->getMomentum()[2];
      pT                          = sqrt(px * px + py * py);
      p                           = sqrt(pT * pT + pz * pz);
      costheta                    = fabs(pz) / p;
      energy                      = pPfo->getEnergy();
      mass                        = pPfo->getMass();
      charge                      = pPfo->getCharge();

      //MC linker to check if it belongs to the signal
      LCObjectVec mcvec;
      mcvec            = m_reltrue->getRelatedToObjects(pPfo);
      atLeastOneSignal = 0;

      streamlog_out(DEBUG) << " PFO with number of linked MC particles = " << mcvec.size() << endl;

      // if the PFO is associated to MC Particle, check if the MC Particle belongs to physicsParticles
      if (mcvec.size() > 0) {
        for (auto imcp : mcvec) {
          MCParticle* mc_part = dynamic_cast<MCParticle*>(imcp);

          streamlog_out(DEBUG) << " MC Type PDG    Energy" << endl;
          stringstream output;
          output << fixed;
          output << setprecision(precision);
          if (mc_part != NULL) {
            FORMATTED_OUTPUT_MC(output, mc_part->getPDG(), mc_part->getEnergy());
            streamlog_out(DEBUG) << output.str();
          }

          if (find(physicsParticles.begin(), physicsParticles.end(), mc_part) != physicsParticles.end() &&
              atLeastOneSignal == 0) {
            streamlog_out(DEBUG) << " -> belongs to the PhysicsParticle (skip the others)" << std::endl;
            atLeastOneSignal = 1;
          } else {
            streamlog_out(DEBUG) << "-> does not belong to the PhysicsParticle" << std::endl;
          }
        }
      }

      en_tot += energy;
      if (atLeastOneSignal) {
        en_tot_signal += energy;
      } else {
        en_tot_background += energy;
      }

      const TrackVec   tracks   = pPfo->getTracks();
      const ClusterVec clusters = pPfo->getClusters();
      nTracks                   = tracks.size();
      nClusters                 = clusters.size();

      //get track time
      float trackTime       = numeric_limits<float>::max();
      clusterTime           = 999.;
      clusterTimeEcal       = 999.;
      clusterTimeHcalEndcap = 999.;

      streamlog_out(DEBUG) << " PFO with number of tracks = " << tracks.size() << endl;
      for (unsigned int trk = 0; trk < tracks.size(); trk++) {
        const Track* track = tracks[trk];
        float        tof;
        const float  time = PfoUtil::TimeAtEcal(track, tof);
        if (fabs(time) < trackTime) {
          trackTime = time;
        }
      }

      //get clusters time
      streamlog_out(DEBUG) << " PFO with number of clusters = " << clusters.size() << endl;
      for (unsigned int clu = 0; clu < clusters.size(); clu++) {
        float meanTime(999.);
        float meanTimeEcal(999.);
        float meanTimeHcalEndcap(999.);
        int   nEcal(0);
        int   nHcalEnd(0);
        int   nCaloHitsUsed(0);

        const Cluster* cluster = clusters[clu];
        PfoUtil::GetClusterTimes(cluster, meanTime, nCaloHitsUsed, meanTimeEcal, nEcal, meanTimeHcalEndcap, nHcalEnd, false);

        // correct for track propagation time
        if (!tracks.empty()) {
          meanTime -= trackTime;
          meanTimeEcal -= trackTime;
          meanTimeHcalEndcap -= trackTime;
        }

        if (fabs(meanTime) < clusterTime) {
          clusterTime = meanTime;
          nCaloHits   = nCaloHitsUsed;
        }
        if (fabs(meanTimeEcal) < clusterTimeEcal) {
          clusterTimeEcal = meanTimeEcal;
          nEcalHits       = nEcal;
        }
        if (fabs(meanTimeHcalEndcap) < clusterTimeHcalEndcap) {
          clusterTimeHcalEndcap = meanTimeHcalEndcap;
          nHcalEndCapHits       = nHcalEnd;
        }
      }

      //check if tracks and clusters have at least one MC particle in common
      trk_clu_sameMCPart = 0;
      for (unsigned int it = 0; it < tracks.size(); it++) {
        //MC linker to tracks
        Track*      track = tracks[it];
        LCObjectVec mctrkvec;
        mctrkvec = m_trackreltrue->getRelatedToObjects(track);

        if (mctrkvec.size() > 0) {
          for (unsigned int ic = 0; ic < clusters.size(); ic++) {
            //MC linker to clusters
            Cluster*    cluster = clusters[ic];
            LCObjectVec mccluvec;
            mccluvec = m_clureltrue->getRelatedToObjects(cluster);

            // check if MC Particle associated to the track
            // has the same Id of at least one MC Particle associated to the cluster
            if (mccluvec.size() > 0) {
              for (auto imctrk : mctrkvec) {
                if (trk_clu_sameMCPart == 1)
                  break;
                for (auto imcclu : mccluvec) {
                  if (trk_clu_sameMCPart == 0 && imcclu->id() == imctrk->id()) {
                    streamlog_out(DEBUG) << " -> PFO track and cluster have same MCParticle" << std::endl;
                    trk_clu_sameMCPart = 1;
                    break;
                  }
                }
              }
            }
          }
        }
      }

      if (trk_clu_sameMCPart == 0)
        streamlog_out(DEBUG) << " -> PFO track and cluster have NOT same MCParticle" << std::endl;

      //computing energy for signal and overlay for the particles categories
      std::string h_label;
      if (type == 22 &&
          std::find(particleCategories.begin(), particleCategories.end(), "photons") != particleCategories.end()) {
        h_label = "photons";
      } else if (type != 22 && charge == 0 &&
                 std::find(particleCategories.begin(), particleCategories.end(), "neutralHadrons") !=
                     particleCategories.end()) {
        h_label = "neutralHadrons";
      } else if (charge != 0 &&
                 std::find(particleCategories.begin(), particleCategories.end(), "chargedPfos") !=
                     particleCategories.end()) {
        h_label = "chargedPfos";
      }

      if (std::find(generationCategories.begin(), generationCategories.end(), "all") != generationCategories.end()) {
        energy_tot[h_label + "_all"] += energy;
        if (costheta < cutCosTheta) {
          energy_tot_central[h_label + "_all"] += energy;
        } else {
          energy_tot_forward[h_label + "_all"] += energy;
        }
      }
      if (atLeastOneSignal &&
          std::find(generationCategories.begin(), generationCategories.end(), "signal") != generationCategories.end()) {
        energy_tot[h_label + "_signal"] += energy;
        if (costheta < cutCosTheta) {
          energy_tot_central[h_label + "_signal"] += energy;
        } else {
          energy_tot_forward[h_label + "_signal"] += energy;
        }
      }
      if (!atLeastOneSignal &&
          std::find(generationCategories.begin(), generationCategories.end(), "overlay") != generationCategories.end()) {
        energy_tot[h_label + "_overlay"] += energy;
        if (costheta < cutCosTheta) {
          energy_tot_central[h_label + "_overlay"] += energy;
        } else {
          energy_tot_forward[h_label + "_overlay"] += energy;
        }
      }

      streamlog_out(DEBUG) << " Type          PDG    E      Pt  cosTheta #trk time  #Clu  time   ecal  hcal  " << endl;

      stringstream output;
      output << fixed;
      output << setprecision(precision);

      if (clusters.size() == 0)
        FORMATTED_OUTPUT_TRACK_CLUSTER_full(output, type, energy, pT, costheta, tracks.size(), trackTime, "-", "-", "-",
                                            "-");
      if (tracks.size() == 0)
        FORMATTED_OUTPUT_TRACK_CLUSTER_full(output, type, energy, pT, costheta, "", "-", clusters.size(), clusterTime,
                                            clusterTimeEcal, clusterTimeHcalEndcap);
      if (tracks.size() > 0 && clusters.size() > 0)
        FORMATTED_OUTPUT_TRACK_CLUSTER_full(output, type, energy, pT, costheta, tracks.size(), trackTime, clusters.size(),
                                            clusterTime, clusterTimeEcal, clusterTimeHcalEndcap);

      streamlog_out(DEBUG) << output.str();

      pfo_tree->Fill();

      streamlog_out(DEBUG) << " " << std::endl;
    }
  }

  h_energy_tot->Fill(en_tot);
  h_energy_tot_signal->Fill(en_tot_signal);
  h_energy_tot_background->Fill(en_tot_background);
  streamlog_out(DEBUG) << " Energy: TOTAL = " << en_tot << ", SIGNAL = " << en_tot_signal
                       << ", OVERLAY = " << en_tot_background << std::endl;

  for (auto ipc : particleCategories) {
    for (auto igc : generationCategories) {
      graphLabel = ipc + "_" + igc;
      h_energy[graphLabel]->Fill(energy_tot[graphLabel]);
      h_energy_central[graphLabel]->Fill(energy_tot_central[graphLabel]);
      h_energy_forward[graphLabel]->Fill(energy_tot_forward[graphLabel]);
      streamlog_out(DEBUG) << " Energy (" << graphLabel << " - all    ): " << energy_tot[graphLabel] << std::endl;
      streamlog_out(DEBUG) << " Energy (" << graphLabel << " - central): " << energy_tot_central[graphLabel] << std::endl;
      streamlog_out(DEBUG) << " Energy (" << graphLabel << " - forward): " << energy_tot_forward[graphLabel] << std::endl;
    }
  }
}

void CLICPfoSelectorAnalysis::fillPlots() {
  Int_t nEntries = pfo_tree->GetEntries();
  streamlog_out(DEBUG) << "CLICPfoSelectorAnalysis ---> filling plots with TTree (nEntries = " << nEntries << ")"
                       << std::endl;

  for (int ie = 0; ie < nEntries; ie++) {
    pfo_tree->GetEntry(ie);
    double currentClusterTime = clusterTime;
    bool   useHcalTimingOnly  = ((costheta > forwardCosThetaForHighEnergyNeutralHadrons) && (type == 2112) &&
                              (energy > forwardHighEnergyNeutralHadronsEnergy));

    //in the case of photons, the time of ECAL clusters are used
    if (analyzePhotons == true && type == 22) {
      if (std::find(particleCategories.begin(), particleCategories.end(), "photons") != particleCategories.end() &&
          std::find(generationCategories.begin(), generationCategories.end(), "all") != generationCategories.end()) {
        currentClusterTime = clusterTimeEcal;
        //streamlog_out( DEBUG ) << "Filling scatter plot for photon with pT and current clusterTime: " << pT << "," << currentClusterTime << endl;

        g_timeVsPt["photons_all"]->SetPoint(ie, pT, currentClusterTime);
        if (costheta < cutCosTheta) {
          g_timeVsPt_central["photons_all"]->SetPoint(ie, pT, currentClusterTime);
        } else {
          g_timeVsPt_forward["photons_all"]->SetPoint(ie, pT, currentClusterTime);
        }

        if (analyzeSignal && atLeastOneSignal &&
            std::find(generationCategories.begin(), generationCategories.end(), "signal") != generationCategories.end()) {
          g_timeVsPt["photons_signal"]->SetPoint(ie, pT, currentClusterTime);
          if (costheta < cutCosTheta) {
            g_timeVsPt_central["photons_signal"]->SetPoint(ie, pT, currentClusterTime);
          } else {
            g_timeVsPt_forward["photons_signal"]->SetPoint(ie, pT, currentClusterTime);
          }
        }

        if (analyzeOverlay && !atLeastOneSignal &&
            std::find(generationCategories.begin(), generationCategories.end(), "overlay") != generationCategories.end()) {
          g_timeVsPt["photons_overlay"]->SetPoint(ie, pT, currentClusterTime);
          if (costheta < cutCosTheta) {
            g_timeVsPt_central["photons_overlay"]->SetPoint(ie, pT, currentClusterTime);
          } else {
            g_timeVsPt_forward["photons_overlay"]->SetPoint(ie, pT, currentClusterTime);
          }
        }

      } else {
        streamlog_out(ERROR) << "Cannot fill scatter plots because TGraph for photons was not created. " << endl;
        exit(0);
      }
    }

    if (analyzeNeutralHadrons == true && type != 22 && charge == 0) {
      if (std::find(particleCategories.begin(), particleCategories.end(), "neutralHadrons") != particleCategories.end() &&
          std::find(generationCategories.begin(), generationCategories.end(), "all") != generationCategories.end()) {
        //in the case the nEcalHits is more than expected, the time computed Ecal is used
        if (!useHcalTimingOnly && (nEcalHits > minECalHits || nEcalHits >= nCaloHits / 2.)) {
          currentClusterTime = clusterTimeEcal;
          //in the case the nHcalEndCapHits is more than expected, the time computed Hcal endcap is used
        } else if ((nHcalEndCapHits >= minHcalEndcapHits) || (nHcalEndCapHits >= nCaloHits / 2.)) {
          currentClusterTime = clusterTimeHcalEndcap;
        }
        //streamlog_out( DEBUG ) << "Filling scatter plot for neutralHadrons with pT and clusterTime: " << pT << "," << currentClusterTime << endl;

        g_timeVsPt["neutralHadrons_all"]->SetPoint(ie, pT, currentClusterTime);
        if (costheta < cutCosTheta)
          g_timeVsPt_central["neutralHadrons_all"]->SetPoint(ie, pT, currentClusterTime);
        else
          g_timeVsPt_forward["neutralHadrons_all"]->SetPoint(ie, pT, currentClusterTime);

        if (analyzeSignal && atLeastOneSignal &&
            std::find(generationCategories.begin(), generationCategories.end(), "signal") != generationCategories.end()) {
          g_timeVsPt["neutralHadrons_signal"]->SetPoint(ie, pT, currentClusterTime);
          if (costheta < cutCosTheta)
            g_timeVsPt_central["neutralHadrons_signal"]->SetPoint(ie, pT, currentClusterTime);
          else
            g_timeVsPt_forward["neutralHadrons_signal"]->SetPoint(ie, pT, currentClusterTime);
        }

        if (analyzeOverlay && !atLeastOneSignal &&
            std::find(generationCategories.begin(), generationCategories.end(), "overlay") != generationCategories.end()) {
          g_timeVsPt["neutralHadrons_overlay"]->SetPoint(ie, pT, currentClusterTime);
          if (costheta < cutCosTheta)
            g_timeVsPt_central["neutralHadrons_overlay"]->SetPoint(ie, pT, currentClusterTime);
          else
            g_timeVsPt_forward["neutralHadrons_overlay"]->SetPoint(ie, pT, currentClusterTime);
        }

      } else {
        streamlog_out(ERROR) << "Cannot fill scatter plots because TGraph for neutralHadrons was not created. " << endl;
        exit(0);
      }
    }

    if (analyzeChargedPfos == true && charge != 0) {
      if (std::find(particleCategories.begin(), particleCategories.end(), "chargedPfos") != particleCategories.end() &&
          std::find(generationCategories.begin(), generationCategories.end(), "all") != generationCategories.end()) {
        //in the case the nEcalHits is more than expected, the time computed Ecal is used
        if ((!useHcalTimingOnly && (nEcalHits > minECalHits || nEcalHits >= nCaloHits / 2.))) {
          currentClusterTime = clusterTimeEcal;
          //in the case the nHcalEndCapHits is more than expected, the time computed Hcal endcap is used
        } else if ((nHcalEndCapHits >= minHcalEndcapHits) || (nHcalEndCapHits >= nCaloHits / 2.)) {
          currentClusterTime = clusterTimeHcalEndcap;
        }
        //streamlog_out( DEBUG ) << "Filling scatter plot for chargedPfos with pT and clusterTime: " << pT << "," << currentClusterTime << endl;

        g_timeVsPt["chargedPfos_all"]->SetPoint(ie, pT, currentClusterTime);
        if (costheta < cutCosTheta)
          g_timeVsPt_central["chargedPfos_all"]->SetPoint(ie, pT, currentClusterTime);
        else
          g_timeVsPt_forward["chargedPfos_all"]->SetPoint(ie, pT, currentClusterTime);

        if (analyzeSignal && atLeastOneSignal &&
            std::find(generationCategories.begin(), generationCategories.end(), "signal") != generationCategories.end()) {
          g_timeVsPt["chargedPfos_signal"]->SetPoint(ie, pT, currentClusterTime);
          if (costheta < cutCosTheta)
            g_timeVsPt_central["chargedPfos_signal"]->SetPoint(ie, pT, currentClusterTime);
          else
            g_timeVsPt_forward["chargedPfos_signal"]->SetPoint(ie, pT, currentClusterTime);
        }

        if (analyzeOverlay && !atLeastOneSignal &&
            std::find(generationCategories.begin(), generationCategories.end(), "overlay") != generationCategories.end()) {
          g_timeVsPt["chargedPfos_overlay"]->SetPoint(ie, pT, currentClusterTime);
          if (costheta < cutCosTheta)
            g_timeVsPt_central["chargedPfos_overlay"]->SetPoint(ie, pT, currentClusterTime);
          else
            g_timeVsPt_forward["chargedPfos_overlay"]->SetPoint(ie, pT, currentClusterTime);
        }

      } else {
        streamlog_out(ERROR) << "Cannot fill scatter plots because TGraph for chargedPfos was not created. " << endl;
        exit(0);
      }
    }
  }
}
