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
    m_trackColName(""), m_linkColName(""),
    m_rootFileName(""), rootfile(NULL), tree(NULL),
    nTracks(0), lastRunHeaderProcessed(-1)
    {

  // modify processor description
  _description = "AnalyseSidEdxProcessor reads Si tracker dEdx data from the .slcio file"
      " and stores them in a root file for analysis." ;


  // register steering parameters: name, description, class-variable, default value

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


  registerProcessorParameter("RootFilename" ,
                             "Name of output root file",
                             m_rootFileName ,
                             std::string("SiTracker_dEdx.root") ) ;


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
  tree->Branch("Edep", &eDep);
  tree->Branch("dEdx", &dEdx);
  tree->Branch("d0", &d0);
  tree->Branch("nTrkHits", &nTrkHits);
  tree->Branch("nTrkRelatedParticles", &nTrkRelatedParticles);
  tree->Branch("zHit", &zHit);
  tree->Branch("xHit", &xHit);
  tree->Branch("yHit", &yHit);
  tree->Branch("eHit", &eHit);
  tree->Branch("typeHit", &typeHit);

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


  /************************************/
  /***       Get collections        ***/
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

  tracks->getFlag();

  nTracks = tracks->getNumberOfElements()  ;

  for (int i = 0; i < nTracks; i++)
  {
    /***       Reset variables ***/
    pMC.clear();
    thetaMC.clear();
    eDep.clear();
    dEdx.clear();
    nTrkHits.clear();
    nTrkRelatedParticles.clear();
    zHit.clear(); xHit.clear(); yHit.clear(); eHit.clear(); typeHit.clear();
    d0.clear();

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
    double d_0;
    MCParticleImpl *mcp = NULL;
    for (unsigned int iw = 0; iw < mcpWeights.size(); iw++)
    {
      if (mcpWeights[iw] > maxw) {
        mcp = dynamic_cast<MCParticleImpl*>(mcParticles[iw]);
        maxw = mcpWeights[iw];
        d_0 = *mcp->getVertex();
      }
    }
    TVector3 p3(mcp->getMomentum());

/*    pMC.push_back(float(p3.Mag()));
    thetaMC.push_back(float(p3.Theta()));
    d0.push_back(d_0);

    dEdx.push_back(track->getdEdx());
*/

    /*** Individual hits ***/

    EVENT::TrackerHitVec trackhits = track->getTrackerHits();
//    nTrkHits.push_back(trackhits.size());

    float eDepHit = 0.;

    for(unsigned int ihit = 0; ihit < trackhits.size(); ihit++) {

      // Hit position
      TVector3 hitpos(trackhits[ihit]->getPosition());

      zHit.push_back(hitpos.z()); xHit.push_back(hitpos.x()); yHit.push_back(hitpos.y());
      eHit.push_back(trackhits[ihit]->getEDep());
      typeHit.push_back(double(trackhits[ihit]->getType()));

      /****************/
      pMC.push_back(float(p3.Mag()));
      thetaMC.push_back(float(p3.Theta()));
      d0.push_back(d_0);

      dEdx.push_back(track->getdEdx());
      nTrkHits.push_back(trackhits.size());
      /**********************/

      eDepHit += trackhits[ihit]->getEDep();
    }

    eDep.push_back(eDepHit);
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

