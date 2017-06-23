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

#include <DD4hep/Detector.h>
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


dEdxPoint::dEdxPoint(const double _dE, const double _dx) :
    dE(_dE), dx(_dx), dEdx(_dE/_dx) {}

dEdxPoint::dEdxPoint(const dEdxPoint& orig) :
    dE(orig.Get_dE()), dx(orig.Get_dx()), dEdx(orig.Get_dE()/orig.Get_dx()) {}

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

  /* Type of estimator for dEdx
   * Available estimators: mean, median, truncMean, harmonic, harmonic-2, weighted-harmonic, weighted-harmonic-2
   */
  registerProcessorParameter("dEdxEstimator" ,
                             "Type of estimator for dEdx.",
                             m_dEdxEstimator ,
                             std::string("mean") ) ;

}



void SiTracker_dEdxProcessor::init() {

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  if(m_dEdxEstimator.compare("mean") == 0) {
    dEdxEval = &SiTracker_dEdxProcessor::dEdxMean;
  }
  else if(m_dEdxEstimator.compare("median") == 0) {
    dEdxEval = &SiTracker_dEdxProcessor::dEdxMedian;
  }
  else if (m_dEdxEstimator.compare("truncMean") == 0) {
    dEdxEval = &SiTracker_dEdxProcessor::dEdxTruncMean;
  }
  else if (m_dEdxEstimator.compare("harmonic") == 0) {
    dEdxEval = &SiTracker_dEdxProcessor::dEdxHarmonic;
  }
  else if (m_dEdxEstimator.compare("harmonic-2") == 0) {
    dEdxEval = &SiTracker_dEdxProcessor::dEdxHarmonic2;
  }
  else if (m_dEdxEstimator.compare("weighted-harmonic") == 0) {
    dEdxEval = &SiTracker_dEdxProcessor::dEdxWgtHarmonic;
  }
  else if (m_dEdxEstimator.compare("weighted-harmonic-2") == 0) {
    dEdxEval = &SiTracker_dEdxProcessor::dEdxWgtHarmonic2;
  }
  else {
    streamlog_out(ERROR) << "Unknown dE/dx evaluation method " << m_dEdxEstimator << ". Exiting.\n";
    exit(0);
  }

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension< dd4hep::rec::SurfaceManager >() ;
  surfMap = surfMan.map( "tracker" ) ;

  if(!m_cheatSensorThicknesses) {
    m_sensThicknessCheatVals.clear();
    for(unsigned isub=0; isub<m_trkHitCollNames.size(); isub++) {
      // This is how we tell the collection finder that we do not want to cheat any thickness values
      m_sensThicknessCheatVals.push_back(-1.);
    }
  }

  collFinder = new LayerFinder(m_trkHitCollNames, theDetector, m_sensThicknessCheatVals, m_elementMask);
 // exit(0);

  const double pos[3]={0,0,0};
  double bFieldVec[3]={0,0,0};
  theDetector.field().magneticField(pos,bFieldVec); // get the magnetic field vector from DD4hep
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

    unsigned iRegHits = 0;

    dEdxVec dEdxHitVec;

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
      dd4hep::rec::SurfaceMap::const_iterator surface = surfMap->find(cellid);
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
        streamlog_out(ERROR) << "Could not find hit collection corresponding to hit CellID " << cellid
                             << ", hit ID " << trackhits.at(ihit)->id() << " .\n";
        streamlog_out(ERROR) << "Event #" << evt->getEventNumber() << ".\n";
        exit(0);
      }
      if (thickness < 0.00001) {
        streamlog_out(ERROR) << "GetThickness returned zero!\n";
        exit(0);
      }

      double effThickness = thickness / cosAngle;

      dEdxHitVec.push_back( dEdxPoint(trackhits[ihit]->getEDep(), effThickness) );
      // At present there is no method to set dEdx for a EVENT::TrackerHit, IMPL::TrackerHitImpl or IMPL::TrackerHitPlaneImpl

      // I am not sure whether the following is the intended use of the hit "type".
      // The hit type value is being overwritten here, but it was apparently not used before.
      ((IMPL::TrackerHitImpl*)(trackhits[ihit]))->setType(detTypeFlag);

      iRegHits++;
    }
    if(iRegHits == 0) continue;

    double dEdx, dEdxError;
    dEdx = dEdxEval(dEdxHitVec, dEdxError);
    track->setdEdx(dEdx);
    track->setdEdxError(dEdxError);
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



/*************************************************************
 *
 *   Evaluation methods for dE/dx
 *
 ************************************************************/

double SiTracker_dEdxProcessor::truncFractionUp = 0.3;
double SiTracker_dEdxProcessor::truncFractionLo = 0.1;

// Weighted mean
// Weight of a measurement is the material thickness traversed in the hit
double SiTracker_dEdxProcessor::dEdxMean(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) return 0;

  double eDepSum = 0.;
  double thickness = 0.;
  double mu2dEdx = 0.;
  for (unsigned i=0; i<n; i++) {
    eDepSum += hitVec.at(i).Get_dE();
    thickness =+ hitVec.at(i).Get_dx();
    mu2dEdx += pow(hitVec.at(i).Get_dE(), 2) / hitVec.at(i).Get_dx();
  }

  double track_dEdx = eDepSum / thickness;
  dEdxError = sqrt((mu2dEdx - pow(track_dEdx, 2))/n);
  return track_dEdx;
}


// Median
double SiTracker_dEdxProcessor::dEdxMedian(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) return 0;

  sort(hitVec.begin(), hitVec.end(), dEdxOrder);
  double median=0.;
  if (n%2 ==1) {
    median = hitVec.at(n/2).Get_dEdx();
  }
  else {
    median = (hitVec.at(n/2-1).Get_dEdx() + hitVec.at(n/2).Get_dEdx()) / 2;
  }

  double eDepSum = 0.;
  double thickness = 0.;
  double mu2dEdx = 0.;
  for (unsigned i=0; i<n; i++) {
    eDepSum += hitVec.at(i).Get_dE();
    thickness =+ hitVec.at(i).Get_dx();
    mu2dEdx += pow(hitVec.at(i).Get_dE(), 2) / hitVec.at(i).Get_dx();
  }

  double mean = eDepSum / thickness;
  // Temporarily using sigma/sqrt(n) for the error of the median
  // TODO: Implement a more accurate estimate of the error of the median
  // using a correction factor tbd from toy MC with Landau distribution
  dEdxError = sqrt((mu2dEdx - pow(mean, 2))/n);
  return median;
}


// Weighted truncated mean
// Weight of a measurement is the material thickness traversed in the hit
double SiTracker_dEdxProcessor::dEdxTruncMean(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) return 0;

  sort(hitVec.begin(), hitVec.end(), dEdxOrder);
  const unsigned iStart = unsigned(n*truncFractionLo + 0.5);
  const unsigned iEnd = unsigned(n*(1-truncFractionUp) + 0.5);

  double eDepSum = 0.;
  double thickness = 0.;
  double mu2dEdx = 0.;
  for (unsigned i=iStart; i<iEnd; i++) {
    eDepSum += hitVec.at(i).Get_dE();
    thickness =+ hitVec.at(i).Get_dx();
    mu2dEdx += pow(hitVec.at(i).Get_dE(), 2) / hitVec.at(i).Get_dx();
  }

  double track_dEdx = eDepSum / thickness;
  dEdxError = sqrt((mu2dEdx - pow(track_dEdx, 2))/(iEnd-iStart));
  return track_dEdx;
}


// Simple harmonic mean
double SiTracker_dEdxProcessor::dEdxHarmonic(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) return 0;

  // Calculation of the first and the second moment of
  // 1 / (dE/dx)
  double mu1sum = 0.;
  double mu2sum = 0.;
  for (unsigned i=0; i<n; i++) {
    double inverse = 1/hitVec.at(i).Get_dEdx();
    mu1sum += inverse;
    mu2sum += pow( inverse, 2 );
  }

  double mu2 = mu2sum / n;
  double mu1 = mu1sum / n;
  double sigma = sqrt( mu2 - pow(mu1, 2) );

  double dEdx = 1 / mu1;
  dEdxError = sigma * pow(dEdx, 2) / sqrt(n) ;

  return dEdx;
}


// Simple harmonic-squared mean
double SiTracker_dEdxProcessor::dEdxHarmonic2(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) return 0;

  // Calculation of the first and the second moment of
  // 1 / (dE/dx)^2
  double mu1sum = 0.;
  double mu2sum = 0.;
  for (unsigned i=0; i<n; i++) {
    double sqinverse = pow(hitVec.at(i).Get_dEdx(), -2);
    mu1sum += sqinverse;
    mu2sum += pow( sqinverse, 2 );
  }

  double mu2 = mu2sum / n;
  double mu1 = mu1sum / n;
  double sigma = sqrt( mu2 - pow(mu1, 2) );

  double dEdx = 1 / sqrt(mu1);
  dEdxError = sigma * pow(dEdx, 3) / 2 / sqrt(n) ;

  return dEdx;
}


// Weighted harmonic mean
// Weight of a measurement is the material thickness traversed in the hit
double SiTracker_dEdxProcessor::dEdxWgtHarmonic(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) return 0;

  // Calculation of the first and the second moment of
  // 1 / (dE/dx)
  double mu1sum = 0.;
  double mu2sum = 0.;
  double wgtsum = 0.;
  for (unsigned i=0; i<n; i++) {
    double inverse = 1/hitVec.at(i).Get_dEdx();
    double wgt = hitVec.at(i).Get_dx();
    mu1sum += wgt*inverse;
    mu2sum += wgt*pow( inverse, 2 );
    wgtsum += wgt;
  }

  double mu2 = mu2sum / wgtsum;
  double mu1 = mu1sum / wgtsum;
  double sigma = sqrt( mu2 - pow(mu1, 2) );

  double dEdx = 1 / mu1;
  dEdxError = sigma * pow(dEdx, 2) / sqrt(n) ;

  return dEdx;
}


// Weighted harmonic-squared mean
// Weight of a measurement is the material thickness traversed in the hit
double SiTracker_dEdxProcessor::dEdxWgtHarmonic2(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) return 0;

  // Calculation of the first and the second moment of
  // 1 / (dE/dx)^2
  double mu1sum = 0.;
  double mu2sum = 0.;
  double wgtsum = 0.;
  for (unsigned i=0; i<n; i++) {
    double sqinverse = pow(hitVec.at(i).Get_dEdx(), -2);
    double wgt = hitVec.at(i).Get_dx();
    mu1sum += wgt*sqinverse;
    mu2sum += wgt*pow( sqinverse, 2 );
    wgtsum += wgt;
  }

  double mu2 = mu2sum / wgtsum;
  double mu1 = mu1sum / wgtsum;
  double sigma = sqrt( mu2 - pow(mu1, 2) );

  double dEdx = 1 / sqrt(mu1);
  dEdxError = sigma * pow(dEdx, 3) / 2 / sqrt(n) ;

  return dEdx;
}

