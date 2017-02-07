/*
 * AnalyseSidEdxProcessor.cc
 *
 *  Created on: Dec 15, 2016
 *      Author: S. Lukic, Vinca Belgrade
 */




#include "AnalyseSidEdxProcessor.h"

#include <EVENT/LCCollection.h>
#include <IMPL/MCParticleImpl.h>
#include <IMPL/TrackImpl.h>
//#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include <UTIL/LCRelationNavigator.h>

// ----- include for verbosity dependent logging ---------
#include "marlin/VerbosityLevels.h"

//#include "marlin/Global.h"
//#include <MarlinTrk/IMarlinTrack.h>
//#include <MarlinTrk/IMarlinTrkSystem.h>
//#include "MarlinTrk/MarlinTrkUtils.h"

#include <TVector3.h>
#include <TROOT.h>



AnalyseSidEdxProcessor aAnalyseSidEdxProcessor ;


AnalyseSidEdxProcessor::AnalyseSidEdxProcessor() : Processor("AnalyseSidEdxProcessor"),
    m_rootFileName(""), m_trackColName(""), m_linkColName(""),
    m_trkHitCollNames(),
    rootfile(NULL), tree(NULL),
    nTracks(0), lastRunHeaderProcessed(-1)
    {

  // modify processor description
  _description = "AnalyseSidEdxProcessor reads Si tracker dEdx data from the .slcio file"
      " and stores them in a root file for analysis." ;


  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter("RootFilename" ,
                             "Name of output root file",
                             m_rootFileName ,
                             std::string("SiTracker_dEdx.root") ) ;

  registerInputCollection( LCIO::TRACK,
      "TrackCollectionName" ,
      "Name of the input Track collection"  ,
      m_trackColName ,
      std::string("SiTracks")
  );

  registerInputCollection( LCIO::LCRELATION,
      "LinkCollectionName" ,
      "Name of the LCRelation collection"  ,
            m_linkColName ,
      std::string("TrackMCTruthLink")
  );

  StringVec defaultTrkHitCollections;
  defaultTrkHitCollections.push_back(std::string("ITrackerHits"));
  defaultTrkHitCollections.push_back(std::string("ITrackerEndcapHits"));
  defaultTrkHitCollections.push_back(std::string("OTrackerHits"));
  defaultTrkHitCollections.push_back(std::string("OTrackerEndcapHits"));
  defaultTrkHitCollections.push_back(std::string("VXDTrackerHits"));
  defaultTrkHitCollections.push_back(std::string("VXDEndcapTrackerHits"));

  registerProcessorParameter("TrkHitCollections" ,
                             "Tracker hit collections that will be analysed",
                             m_trkHitCollNames ,
                             defaultTrkHitCollections ) ;

}



void AnalyseSidEdxProcessor::init() {

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  gROOT->ProcessLine("#include <vector>");

  rootfile = new TFile(m_rootFileName.c_str(), "RECREATE");
  tree = new TTree ("SiTracks", "SiTracks");
  tree->Branch("nTracks", &nTracks, "nTracks/I");
  tree->Branch("pMC", &pMC);
  tree->Branch("thetaMC", &thetaMC);
  tree->Branch("eTrack", &eTrack);
  tree->Branch("dEdxTrack", &dEdxTrack);
  tree->Branch("dEdxError", &dEdxError);
  tree->Branch("eEvt", &eEvt);
  tree->Branch("d0", &d0);
  tree->Branch("m", &m);
  tree->Branch("nTrkHits", &nTrkHits);
  tree->Branch("nTrkRelatedParticles", &nTrkRelatedParticles);
  tree->Branch("zHit", &zHit);
  tree->Branch("xHit", &xHit);
  tree->Branch("yHit", &yHit);
  tree->Branch("eHit", &eHit);
  tree->Branch("typeHit", &typeHit);
  tree->Branch("zTrackHit", &zTrackHit);
  tree->Branch("xTrackHit", &xTrackHit);
  tree->Branch("yTrackHit", &yTrackHit);
  tree->Branch("eTrackHit", &eTrackHit);
  tree->Branch("typeTrackHit", &typeTrackHit);

  lastRunHeaderProcessed = -1;

  streamlog_out(DEBUG) << "   init done  " << std::endl ;

}


void AnalyseSidEdxProcessor::processRunHeader( LCRunHeader* run) {

  lastRunHeaderProcessed = run->getRunNumber();

}



void AnalyseSidEdxProcessor::processEvent( LCEvent * evt ) {

  if(evt->getEventNumber()%1000 == 0) {
    streamlog_out(MESSAGE) << "   processing event: " << evt->getEventNumber()
        << "   in run:  " << evt->getRunNumber() << std::endl ;
  }


  /******************************************/
  /***  Collections of all tracker hits   ***/
  /******************************************/

  zHit.clear(); xHit.clear(); yHit.clear(); eHit.clear(); typeHit.clear();
  double eThisEvt = 0.;

  for(unsigned icoll=0; icoll<m_trkHitCollNames.size(); icoll++) {

    LCCollection* hits = NULL;
    try {
      hits = evt->getCollection(m_trkHitCollNames.at(icoll));
    }
    catch(EVENT::DataNotAvailableException &dataex) {
      streamlog_out(DEBUG) << "Collection " << m_trkHitCollNames.at(icoll) << " not found in event #" << evt->getEventNumber() << ".\n";
      hits = NULL;
      continue;
    }

    int nHits = hits->getNumberOfElements();
    for(int ihit=0; ihit<nHits; ihit++) {
      TrackerHit *hit = dynamic_cast<TrackerHit*>(hits->getElementAt(ihit));
      TVector3 hitpos(hit->getPosition());

      xHit.push_back(hitpos.X());
      yHit.push_back(hitpos.Y());
      zHit.push_back(hitpos.Z());
      eHit.push_back(hit->getEDep());
      typeHit.push_back(hit->getType());
      eThisEvt += hit->getEDep();
    } // Loop over hits in collection

  } // Loop over hit collections



  /************************************/
  /***    Get track collections     ***/
  /************************************/

  LCCollection* tracks = NULL;
  try {
    tracks = evt->getCollection(m_trackColName);
  }
  catch(EVENT::DataNotAvailableException &dataex) {
    streamlog_out(MESSAGE) << "Collection " << m_trackColName << " not found. Skipping event #" << evt->getEventNumber() << ".\n";
    tracks = NULL;
    return;
  }

  LCRelationNavigator *track2mcNav = NULL;
  try {
    track2mcNav = new LCRelationNavigator(evt->getCollection( m_linkColName ));
  }
  catch(EVENT::DataNotAvailableException &dataex) {
    streamlog_out(MESSAGE) << "Collection " << m_linkColName << " not found. Skipping event #" << evt->getEventNumber() << ".\n";
    track2mcNav = NULL;
    return;
  }

  nTracks = tracks->getNumberOfElements()  ;

  /***  Reset variables ***/
  pMC.clear(); thetaMC.clear();
  eTrack.clear(); dEdxTrack.clear(); dEdxError.clear(); eEvt.clear();
  nTrkHits.clear();
  nTrkRelatedParticles.clear();
  zTrackHit.clear(); xTrackHit.clear(); yTrackHit.clear(); eTrackHit.clear(); typeTrackHit.clear();
  d0.clear(); m.clear();

  for (int i = 0; i < nTracks; i++)
  {

    TrackImpl * track = dynamic_cast<TrackImpl*>( tracks->getElementAt(i) );

    /*** Find the associated MC particle ***/

    const LCObjectVec& mcParticles = track2mcNav->getRelatedToObjects(track);
    streamlog_out(DEBUG) << " Number of MCParticles connected to this track " << mcParticles.size() << std::endl;
    const FloatVec& mcpWeights = track2mcNav->getRelatedToWeights(track);
    streamlog_out(DEBUG) << " Number of weights of MCParticle connections " << mcpWeights.size() << std::endl;

    if(mcParticles.size() != mcpWeights.size()) {
      streamlog_out(ERROR) << "Different sizes of MCParticle and weight vectors! Aborting.\n";
      exit(0);
    }

    nTrkRelatedParticles.push_back(mcpWeights.size());

    double maxw = 0.;
    double d_0, mass;
    d_0 = -1.;
    mass = -1.;
    MCParticleImpl *mcp = NULL;
    for (unsigned int iw = 0; iw < mcpWeights.size(); iw++)
    {
      if (mcpWeights[iw] > maxw) {
        mcp = dynamic_cast<MCParticleImpl*>(mcParticles[iw]);
        maxw = mcpWeights[iw];
        d_0 = (TVector3(mcp->getVertex())).Mag();
        mass = mcp->getMass();
      }
    }
    TVector3 p3(mcp->getMomentum());

    pMC.push_back(float(p3.Mag()));
    thetaMC.push_back(float(p3.Theta()));
    d0.push_back(d_0);/**/
    m.push_back(mass);/**/
    dEdxTrack.push_back(track->getdEdx());
    dEdxError.push_back(track->getdEdxError());

    /*** Individual hits belonging to this track ***/

    EVENT::TrackerHitVec trackhits = track->getTrackerHits();
    nTrkHits.push_back(trackhits.size());

    float eDepHit = 0.;

    for(unsigned int ihit = 0; ihit < trackhits.size(); ihit++) {

      // Hit position
      TVector3 hitpos(trackhits[ihit]->getPosition());

      zTrackHit.push_back(hitpos.z()); xTrackHit.push_back(hitpos.x()); yTrackHit.push_back(hitpos.y());
      eTrackHit.push_back(trackhits[ihit]->getEDep());
      typeTrackHit.push_back(double(trackhits[ihit]->getType()));

      /****************
      // Repeatedly store track-related variables with each hit
      // for easier analysis in ROOT
      pMC.push_back(float(p3.Mag()));
      thetaMC.push_back(float(p3.Theta()));
      d0.push_back(d_0);
      m.push_back(mass);

      nTrkHits.push_back(trackhits.size());
      /**********************/

      eDepHit += trackhits[ihit]->getEDep();
    }

    eTrack.push_back(eDepHit);
    // Repeatedly store total energy in event for easier analysis in ROOT
    eEvt.push_back(eThisEvt);
  }






  tree->Fill();

}





void AnalyseSidEdxProcessor::check( LCEvent * evt ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void AnalyseSidEdxProcessor::end(){

  //   std::cout << "AnalyseSidEdxProcessor::end()  " << name()
  //     << " processed " << _nEvt << " events in " << _nRun << " runs "
  //     << std::endl ;
  rootfile->Write();

}

