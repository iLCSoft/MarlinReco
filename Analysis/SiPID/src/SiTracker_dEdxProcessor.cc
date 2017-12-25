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
#include <DDRec/Vector3D.h>

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

/**/

SiTracker_dEdxProcessor & SiTracker_dEdxProcessor::operator = (const SiTracker_dEdxProcessor &orig) {

  if (this == &orig) return *this;

  this->m_trackCollName = orig.m_trackCollName;
  this->m_trkHitCollNames = orig.m_trkHitCollNames;
  this->m_elementMask = orig.m_elementMask;
  this->surfMap = orig.surfMap;
  this->trkSystem = orig.trkSystem;
  this->_bField = orig._bField;
  this->layerFinder = orig.layerFinder;
  this->lastRunHeaderProcessed = orig.lastRunHeaderProcessed;
  this->dEdxEval = orig.dEdxEval;

  return *this;
}

SiTracker_dEdxProcessor::SiTracker_dEdxProcessor(const SiTracker_dEdxProcessor& orig) :
    Processor("SiTracker_dEdxProcessor"),
    dEdxEval(orig.dEdxEval),
    m_trackCollName(orig.m_trackCollName), m_trkHitCollNames(orig.m_trkHitCollNames),
    m_elementMask(orig.m_elementMask), surfMap(orig.surfMap), trkSystem(orig.trkSystem),
    _bField(orig._bField), layerFinder(orig.layerFinder),
    lastRunHeaderProcessed(orig.lastRunHeaderProcessed),
    timers(orig.timers),
    lastTP(std::chrono::high_resolution_clock::now()),
    newTP(std::chrono::high_resolution_clock::now())
{}

SiTracker_dEdxProcessor::SiTracker_dEdxProcessor() : Processor("SiTracker_dEdxProcessor"),
    m_trackCollName(""), m_trkHitCollNames(), m_elementMask(0),
    surfMap(NULL), trkSystem(NULL), _bField(0),
    layerFinder(NULL),
    lastRunHeaderProcessed(-1),
    timers(),
    lastTP(std::chrono::high_resolution_clock::now()),
    newTP(std::chrono::high_resolution_clock::now())
    {

  // modify processor description
  _description = "SiTracker_dEdxProcessor calculates dE/dx for planar silicon trackers" ;


  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::TRACK,
      "TrackCollectionName" ,
      "Name of the input Track collection"  ,
      m_trackCollName ,
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
  for (unsigned ibit=0; ibit<sizeof(int)*CHAR_BIT; ibit++) {
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
                             std::string("median") ) ;

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

  layerFinder = new LayerFinder(m_trkHitCollNames, theDetector, m_sensThicknessCheatVals, m_elementMask);
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

  for (unsigned i=0; i<nTimers; i++) {
    timers.push_back(std::chrono::duration<double>(std::chrono::duration_values<double>::zero()));
  }

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

  lastTP = std::chrono::high_resolution_clock::now();

  LCCollection* tracks = NULL;
  try {
    tracks = evt->getCollection(m_trackCollName);
  }
  catch(EVENT::DataNotAvailableException &dataex) {
    streamlog_out(MESSAGE) << "Collection " << m_trackCollName << " not found. Skipping event #" << evt->getEventNumber() << ".\n";
    tracks = NULL;
    return;
  }

  tracks->getFlag();


  addTime(0);


  /*** Read collections to find a valid decoder for CELLID encoding ***/
  if (layerFinder->ReadCollections(evt) != 0) {
    streamlog_out(WARNING) << "None of the requested collections found in event #" << evt->getEventNumber() << ". Skipping event.\n";
    return;
  }

  layerFinder->ReportKnownDetectors();

  addTime(1);

  int nTracks = tracks->getNumberOfElements()  ;

  for (int i = 0; i < nTracks; i++)
  {
    streamlog_out(DEBUG5) << "Processing track #" << i << ".\n";
    TrackImpl * track = dynamic_cast<TrackImpl*>( tracks->getElementAt(i) );

    /*** Analyse hits and get dE/dx from each ***/

    EVENT::TrackerHitVec trackhits = track->getTrackerHits();

    MarlinTrk::IMarlinTrack* marlin_trk = trkSystem->createTrack();
    for( EVENT::TrackerHitVec::iterator it = trackhits.begin() ; it != trackhits.end() ; ++it ){
      marlin_trk->addHit(*it);
    }//end loop on hits

    const TrackStateImpl *trackState = dynamic_cast<const TrackStateImpl*>(track->getTrackState(TrackState::AtFirstHit));
    if (!trackState) {
      streamlog_out(WARNING) << "Cannot get track state for track #" << i
                             << " in event " << evt->getEventNumber() << std::endl;
      streamlog_out(WARNING) << "Skipping track.\n";
      continue;
    }

    marlin_trk->initialise( *trackState, _bField, MarlinTrk::IMarlinTrack::forward ) ;

    addTime(2);

    dEdxVec dEdxHitVec;

    for(unsigned int ihit = 0; ihit < trackhits.size(); ihit++) {

      // Tangent to the track at hit position
      dd4hep::rec::Vector3D hitpos(trackhits[ihit]->getPosition());

      IMPL::TrackStateImpl ts;
      double chi2 = 0.;
      int ndf = 0;
      marlin_trk->extrapolate(hitpos, ts, chi2, ndf);

      addTime(3);

      dd4hep::rec::Vector3D rp(ts.getReferencePoint());

      float tanLambda = ts.getTanLambda();
      float sinTheta = 1. / sqrt(1.+pow(tanLambda,2));
      float phi = ts.getPhi();
      float trackDirX = cos(phi)*sinTheta;
      float trackDirY = sin(phi)*sinTheta;
      float trackDirZ = tanLambda*sinTheta;
      dd4hep::rec::Vector3D trackDir(trackDirX, trackDirY, trackDirZ);

      addTime(4);

      // Normal to the surface of hit
      unsigned long cellid = trackhits[ihit]->getCellID0();
      dd4hep::rec::SurfaceMap::const_iterator surface = surfMap->find(cellid);
      if (surface == surfMap->end()) {
        streamlog_out(ERROR) << "Cannot find the surface corresponding to track hit ID " << cellid
             << " in event " << evt->getEventNumber() << "!\n";
        exit(0);
      }
      dd4hep::rec::Vector3D surfaceNormal = surface->second->normal();

      addTime(5);

      double norm = sqrt(trackDir.dot(trackDir)*surfaceNormal.dot(surfaceNormal));
      if (norm < FLT_MIN) continue;
      double cosAngle = fabs(trackDir.dot(surfaceNormal)) / norm ;

      double thickness = layerFinder->SensitiveThickness(dynamic_cast<TrackerHitPlane*>(trackhits[ihit]));
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

      addTime(6);

    }


    if(dEdxHitVec.size() == 0) continue;

    double dEdx, dEdxError;
    dEdx = dEdxEval(dEdxHitVec, dEdxError);
    // Todo: This is read-only if track is read from existing lcio file!
    // Is there a way to process tracks that are read from the input file?
    track->setdEdx(dEdx);
    track->setdEdxError(dEdxError);

    addTime(7);
  }

}





void SiTracker_dEdxProcessor::check( LCEvent *  /*evt*/ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}
/**/

void SiTracker_dEdxProcessor::end(){

  for (unsigned i=0; i<timers.size(); i++) {
    streamlog_out(MESSAGE) << "Total time in timer #" << i << ": " << timers.at(i).count() << " s\n";
  }
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

// Weighted truncated mean with arbitrary truncation
// Weight of a measurement is the material thickness traversed in the hit
double SiTracker_dEdxProcessor::dEdxGeneralTruncMean(dEdxVec hitVec, double &dEdxError,
    const double truncLo, const double truncHi) {

  const unsigned n = hitVec.size();
  const unsigned iStart = static_cast<unsigned>(floor(n*truncLo + 0.5));
  const unsigned iEnd = static_cast<unsigned>(floor(n*(1-truncHi) + 0.5));

  if(iEnd-iStart == 0) {
    dEdxError = 0;
    return 0;
  }
  if(iEnd-iStart == 1) {
    dEdxError = hitVec.at(iStart).Get_dEdx() ;
    return hitVec.at(iStart).Get_dEdx() ;
  }

  sort(hitVec.begin(), hitVec.end(), dEdxOrder);

  double eDepSum = 0.;
  double thickness = 0.;
  double mu2dEdx = 0.;
  for (unsigned i=iStart; i<iEnd; i++) {
    eDepSum += hitVec.at(i).Get_dE();
    thickness += hitVec.at(i).Get_dx();
    mu2dEdx += pow(hitVec.at(i).Get_dE(), 2) / hitVec.at(i).Get_dx();
  }

  mu2dEdx /= thickness;

  double track_dEdx = eDepSum / thickness;
  dEdxError = sqrt((mu2dEdx - pow(track_dEdx, 2))/(iEnd-iStart));
  return track_dEdx;
}


// Weighted mean
// Weight of a measurement is the material thickness traversed in the hit
double SiTracker_dEdxProcessor::dEdxMean(dEdxVec hitVec, double &dEdxError) {

  return dEdxGeneralTruncMean(hitVec, dEdxError, 0., 0.);
}


// Median
double SiTracker_dEdxProcessor::dEdxMedian(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) {
    dEdxError = 0;
    return 0;
  }
  if(n == 1) {
    dEdxError = hitVec.at(0).Get_dEdx();
    return hitVec.at(0).Get_dEdx();
  }

  sort(hitVec.begin(), hitVec.end(), dEdxOrder);
  double median=0.;
  if (n%2 ==1) {
    median = hitVec.at(n/2).Get_dEdx();
  }
  else {
    median = (hitVec.at(n/2-1).Get_dEdx() + hitVec.at(n/2).Get_dEdx()) / 2;
  }

  // Substituting error of the mean for dEdxError here
  // instead of bootstrapping the error of the median
  dEdxMean(hitVec, dEdxError);
  return median;
}


// Weighted truncated mean with standard truncation
// Weight of a measurement is the material thickness traversed in the hit
double SiTracker_dEdxProcessor::dEdxTruncMean(dEdxVec hitVec, double &dEdxError) {

  return dEdxGeneralTruncMean(hitVec, dEdxError, truncFractionLo, truncFractionUp);
}


// Simple harmonic mean
double SiTracker_dEdxProcessor::dEdxHarmonic(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) {
    dEdxError = 0;
    return 0;
  }
  if(n == 1) {
    dEdxError = hitVec.at(0).Get_dEdx();
    return hitVec.at(0).Get_dEdx();
  }

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
  double sigma = sqrt( (mu2 - pow(mu1, 2)) / n );

  double dEdx = 1 / mu1;
  dEdxError = sigma / pow(mu1, 2) ;

  return dEdx;
}


// Simple harmonic-squared mean
double SiTracker_dEdxProcessor::dEdxHarmonic2(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) {
    dEdxError = 0;
    return 0;
  }
  if(n == 1) {
    dEdxError = hitVec.at(0).Get_dEdx();
    return hitVec.at(0).Get_dEdx();
  }

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
  double sigma = sqrt( (mu2 - pow(mu1, 2)) / n );

  double dEdx = 1 / sqrt(mu1);
  dEdxError = sigma * pow(dEdx, 3) / 2 ;

  return dEdx;
}


// Weighted harmonic mean
// Weight of a measurement is the material thickness traversed in the hit
double SiTracker_dEdxProcessor::dEdxWgtHarmonic(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) {
    dEdxError = 0;
    return 0;
  }
  if(n == 1) {
    dEdxError = hitVec.at(0).Get_dEdx();
    return hitVec.at(0).Get_dEdx();
  }

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
  double sigma = sqrt( (mu2 - pow(mu1, 2)) / n );

  double dEdx = 1 / mu1;
  dEdxError = sigma * pow(dEdx, 2) ;

  return dEdx;
}


// Weighted harmonic-squared mean
// Weight of a measurement is the material thickness traversed in the hit
double SiTracker_dEdxProcessor::dEdxWgtHarmonic2(dEdxVec hitVec, double &dEdxError) {

  const unsigned n = hitVec.size();
  if(n == 0) {
    dEdxError = 0;
    return 0;
  }
  if(n == 1) {
    dEdxError = hitVec.at(0).Get_dEdx();
    return hitVec.at(0).Get_dEdx();
  }

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
  double sigma = sqrt( (mu2 - pow(mu1, 2)) / n );

  double dEdx = 1 / sqrt(mu1);
  dEdxError = sigma * pow(dEdx, 3) / 2 ;

  return dEdx;
}

