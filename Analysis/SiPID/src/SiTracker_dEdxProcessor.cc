#include "SiTracker_dEdxProcessor.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackImpl.h>
//#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include <UTIL/LCRelationNavigator.h>
#include <IMPL/TrackStateImpl.h>
#include <IMPL/TrackerHitImpl.h>

// ----- include for verbosity dependent logging ---------
#include "marlin/VerbosityLevels.h"

#include <DD4hep/LCDD.h>
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DDGear.h"
#include "DDRec/DetectorData.h"
#include <DD4hep/DetType.h>
#include <DDSurfaces/Vector3D.h>
//#include "DDRec/SurfaceHelper.h"

#include "marlin/Global.h"
#include <MarlinTrk/IMarlinTrack.h>
#include <MarlinTrk/IMarlinTrkSystem.h>
#include "MarlinTrk/MarlinTrkUtils.h"

#include <TVector3.h>
#include <TROOT.h>

#include <climits>



SiTracker_dEdxProcessor aSiTracker_dEdxProcessor ;


SiTracker_dEdxProcessor::SiTracker_dEdxProcessor() : Processor("SiTracker_dEdxProcessor"),
    m_trackColName(""), m_trkHitCollNames(), m_elementMask(0),
    surfMap(NULL), trkSystem(NULL), _bField(0),
    collFinder(NULL),
    lastRunHeaderProcessed(-1)
    {

  // modify processor description
  _description = "SiTracker_dEdxProcessor calculates dE/dx for planar silicon trackers" ;


  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::TRACK,
      "TrackCollectionName" ,
      "Name of the input Track collection"  ,
      m_trackColName ,
      std::string("SiTracks")
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

  int elementMask = 0;
  for (int ibit=0; ibit<sizeof(int)*CHAR_BIT; ibit++) {
    elementMask += 1 << ibit;
  }

  registerProcessorParameter("ElementMask" ,
                             "Bit mask which tracker detector elements to use (respecting the order in LCDD)",
                             m_elementMask ,
                             elementMask ) ;

  registerProcessorParameter("CheatSensorThicknesses" ,
                             "Shall we use the sensitive thicknesses from parameters?",
                             m_cheatSensorThicknesses ,
                             false ) ;

  FloatVec sensThicknessCheatVals;
  for (unsigned i=0; i<defaultTrkHitCollections.size(); i++) {
    sensThicknessCheatVals.push_back(-1);
  }

  registerProcessorParameter("SensorThicknessCheatValues" ,
                             "Sensor thicknesses to use instead of automatic values from DD4hep (if CheatSensorThicknesses==true).",
                             m_sensThicknessCheatVals ,
                             sensThicknessCheatVals ) ;

}



void SiTracker_dEdxProcessor::init() {

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  Geometry::LCDD& lcdd = Geometry::LCDD::getInstance();

  DDRec::SurfaceManager& surfMan = *lcdd.extension< DDRec::SurfaceManager >() ;
  surfMap = surfMan.map( "tracker" ) ;

  if(!m_cheatSensorThicknesses) {
    m_sensThicknessCheatVals.clear();
    for(unsigned isub=0; isub<m_trkHitCollNames.size(); isub++) {
      // This is how we tell the collection finder that we do not want to cheat any thickness values
      m_sensThicknessCheatVals.push_back(-1.);
    }
  }

  collFinder = new LayerFinder(m_trkHitCollNames, lcdd, m_sensThicknessCheatVals, m_elementMask);
 // exit(0);

  const double pos[3]={0,0,0};
  double bFieldVec[3]={0,0,0};
  lcdd.field().magneticField(pos,bFieldVec); // get the magnetic field vector from DD4hep
  _bField = bFieldVec[2]/dd4hep::tesla; // z component at (0,0,0)

  //trksystem for marlin track

  trkSystem =  MarlinTrk::Factory::createMarlinTrkSystem( "DDKalTest" , marlin::Global::GEAR , "" ) ;

  if( trkSystem == 0 ) throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("DDKalTest") );

  trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS,      true );
      //_MSOn ) ;
  trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx,       false) ;
  trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing,  true) ;
  trkSystem->init() ;

  gROOT->ProcessLine("#include <vector>");

  lastRunHeaderProcessed = -1;

  streamlog_out(DEBUG) << "   init done  " << std::endl ;

}


void SiTracker_dEdxProcessor::processRunHeader( LCRunHeader* run) {

  lastRunHeaderProcessed = run->getRunNumber();

} 



void SiTracker_dEdxProcessor::processEvent( LCEvent * evt ) {

/*  if (_lastRunHeaderProcessed < evt->getRunNumber()) {
//    streamlog_out(ERROR) << "Run header has not been processed for this run! Exiting.\n";
//    exit(0);
  }*/

  if(evt->getEventNumber()%100 == 0) {
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

  tracks->getFlag();

  /*** Collection finder for hit collections ***/
  if (collFinder->ReadCollections(evt) == 0) {
    streamlog_out(WARNING) << "None of the requested collections found in event #" << evt->getEventNumber() << ". Skipping event.\n";
    return;
  }

  collFinder->ReportKnownDetectors();


  int nTracks = tracks->getNumberOfElements()  ;

  for (int i = 0; i < nTracks; i++)
  {
    TrackImpl * track = dynamic_cast<TrackImpl*>( tracks->getElementAt(i) );

    /*** Analyse hits and get dE/dx from each ***/

    EVENT::TrackerHitVec trackhits = track->getTrackerHits();

    MarlinTrk::IMarlinTrack* marlin_trk = trkSystem->createTrack();

    for( EVENT::TrackerHitVec::iterator it = trackhits.begin() ; it != trackhits.end() ; ++it ){
      marlin_trk->addHit(*it);
    }//end loop on hits
    TrackStateImpl trackState( *(track->getTrackState(TrackState::AtFirstHit)) );

    marlin_trk->initialise( trackState, _bField, MarlinTrk::IMarlinTrack::forward ) ;

    float eDepHit = 0.;
    float traversedThickness = 0.;
    float dEdxSum = 0.;
    unsigned iRegHits = 0;

    for(unsigned int ihit = 0; ihit < trackhits.size(); ihit++) {

      // Tangent to the track at hit position
      DDSurfaces::Vector3D hitpos(trackhits[ihit]->getPosition());

      IMPL::TrackStateImpl ts;
      double chi2 = 0.;
      int ndf = 0;
      marlin_trk->propagate(hitpos, ts, chi2, ndf);

      DDSurfaces::Vector3D rp(ts.getReferencePoint());

      float tanLambda = ts.getTanLambda();
      float sinTheta = 1. / sqrt(1.+pow(tanLambda,2));
      float phi = ts.getPhi();
      float trackDirX = cos(phi)*sinTheta;
      float trackDirY = sin(phi)*sinTheta;
      float trackDirZ = tanLambda*sinTheta;
      DDSurfaces::Vector3D trackDir(trackDirX, trackDirY, trackDirZ);

      // Normal to the surface of hit
      unsigned long cellid = trackhits[ihit]->getCellID0();
      DDRec::SurfaceMap::const_iterator surface = surfMap->find(cellid);
      if (surface == surfMap->end()) {
        streamlog_out(ERROR) << "Cannot find the surface corresponding to track hit ID " << cellid
             << " in event " << evt->getEventNumber() << "!\n";
        exit(0);
      }
      DDSurfaces::Vector3D surfaceNormal = surface->second->normal();
      streamlog_out(DEBUG) << "Found surface corresponding to track hit ID " << cellid << ".\n";

      double norm = sqrt(trackDir.dot(trackDir)*surfaceNormal.dot(surfaceNormal));
      if (norm < FLT_MIN) continue;
      double cosAngle = fabs(trackDir.dot(surfaceNormal)) / norm ;

      int detTypeFlag = 0;
      double thickness = collFinder->SensitiveThickness(dynamic_cast<TrackerHitPlane*>(trackhits[ihit]), detTypeFlag);
      if (thickness < 0.) {
        streamlog_out(DEBUG) << "Could not find hit collection corresponding to hit CellID " << cellid
                             << ", hit ID " << trackhits.at(ihit)->id() << " .\n";
        streamlog_out(DEBUG) << "Event #" << evt->getEventNumber() << ".\n";
        continue;
      }
      if (thickness < 0.00001) {
        streamlog_out(ERROR) << "GetThickness returned zero!\n";
        exit(0);
      }

      traversedThickness += thickness / cosAngle;
      eDepHit += trackhits[ihit]->getEDep();
      dEdxSum += trackhits[ihit]->getEDep() / (thickness / cosAngle);
      // At present there is no method to set dEdx for a EVENT::TrackerHit, IMPL::TrackerHitImpl or IMPL::TrackerHitPlaneImpl

      // I am not sure whether the following is the intended use of the hit "type".
      // The hit "type" is being overwritten here, but it was apparently not used before.
      ((IMPL::TrackerHitImpl*)(trackhits[ihit]))->setType(detTypeFlag);

      iRegHits++;
    }
    if(iRegHits == 0) continue;

    double track_dEdx = eDepHit/traversedThickness;
//    double track_dEdx = dEdxSum/iRegHits; // Unweighted average dE/dx
    track->setdEdx(track_dEdx);
  }

}





void SiTracker_dEdxProcessor::check( LCEvent * evt ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}
/**/

void SiTracker_dEdxProcessor::end(){

  //   std::cout << "SiTracker_dEdxProcessor::end()  " << name()
  //     << " processed " << _nEvt << " events in " << _nRun << " runs "
  //     << std::endl ;
}

