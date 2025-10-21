#include "TOFEstimators.h"
#include "TOFUtils.h"

#include "EVENT/LCCollection.h"
#include "EVENT/TrackerHitPlane.h"

#include "UTIL/PIDHandler.h"

#include "CLHEP/Random/Randomize.h"
#include "EVENT/SimTrackerHit.h"
#include "MarlinTrk/Factory.h"
#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"
#include "marlin/VerbosityLevels.h"
#include "marlinutil/GeometryUtil.h"

using namespace TOFUtils;
using CLHEP::RandGauss;
using dd4hep::rec::Vector3D;
using EVENT::CalorimeterHit;
using EVENT::Cluster;
using EVENT::LCCollection;
using EVENT::LCObject;
using EVENT::ReconstructedParticle;
using EVENT::SimTrackerHit;
using EVENT::Track;
using EVENT::TrackerHit;
using EVENT::TrackerHitPlane;
using EVENT::TrackState;
using std::string;
using std::vector;

TOFEstimators aTOFEstimators;

TOFEstimators::TOFEstimators() : marlin::Processor("TOFEstimators") {
  _description = "TOFEstimators processor computes time-of-flight of the chosen ReconstructedParticle to the \
                    specified end point (SET hit or Ecal surface). To be used for a further particle ID";

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "ReconstructedParticleCollection",
                          "Name of the ReconstructedParticle collection", _pfoCollectionName,
                          std::string("PandoraPFOs"));

  registerProcessorParameter("ExtrapolateToEcal",
                             "If true, track is extrapolated to the Ecal surface for track length calculation, \
                            time of flight estimated using Ecal hits. If false, track length calculated to the last tracker hit.\
                            Time of flight estimated using SET hit if exists.",
                             _extrapolateToEcal, bool(true));

  registerProcessorParameter("TofMethod", "name of the algorithm which estimates time of flight\
                            to the Ecal surface based on Ecal hits time information.\
                            Available options are: closest, frankAvg, frankFit.\
                            In case of _extrapolateToEcal==false is ignored",
                             _tofMethod, std::string("closest"));

  registerProcessorParameter("TimeResolution", "Time resolution of individual SET strips or Ecal hits in ps",
                             _timeResolution, double(0.));

  registerProcessorParameter("MaxEcalLayer", "Time of flight is calculated using Ecal hits only up to MaxLayer",
                             _maxEcalLayer, int(10));
}

void TOFEstimators::init() {

  if (_tofMethod != "closest" && _tofMethod != "frankAvg" && _tofMethod != "frankFit") {
    throw EVENT::Exception("Invalid steering parameter for TofMethod is passed: " + _tofMethod +
                           "\n Available options are: closest, frankAvg, frankFit");
  }

  marlin::Global::EVENTSEEDER->registerProcessor(this);

  _outputParNames = {"timeOfFlight"};
  _bField = MarlinUtil::getBzAtOrigin();
  // internally we use time resolution in nanoseconds
  _timeResolution = _timeResolution / 1000.;
}

void TOFEstimators::processEvent(EVENT::LCEvent* evt) {
  RandGauss::setTheSeed(marlin::Global::EVENTSEEDER->getSeed(this));
  ++_nEvent;
  streamlog_out(DEBUG9) << std::endl << "==========Event========== " << _nEvent << std::endl;

  LCCollection* pfos = evt->getCollection(_pfoCollectionName);
  PIDHandler pidHandler(pfos);
  int algoID = pidHandler.addAlgorithm(name(), _outputParNames);

  for (int i = 0; i < pfos->getNumberOfElements(); ++i) {
    streamlog_out(DEBUG9) << std::endl << "Starting to analyze " << i + 1 << " PFO" << std::endl;
    ReconstructedParticle* pfo = static_cast<ReconstructedParticle*>(pfos->getElementAt(i));

    int nClusters = pfo->getClusters().size();
    int nTracks = pfo->getTracks().size();

    if (nClusters != 1 || nTracks != 1) {
      // Analyze only simple pfos. Otherwise write dummy zeros
      vector<float> results{0.};
      pidHandler.setParticleID(pfo, 0, 0, 0., algoID, results);
      continue;
    }
    Track* track = pfo->getTracks()[0];
    Cluster* cluster = pfo->getClusters()[0];

    double timeOfFlight = 0.;
    if (_extrapolateToEcal) {
      if (_tofMethod == "closest") {
        timeOfFlight = getTofClosest(cluster, track, _timeResolution);
      } else if (_tofMethod == "frankAvg") {
        vector<CalorimeterHit*> frankHits = selectFrankEcalHits(cluster, track, _maxEcalLayer, _bField);
        timeOfFlight = getTofFrankAvg(frankHits, track, _timeResolution);
      } else if (_tofMethod == "frankFit") {
        vector<CalorimeterHit*> frankHits = selectFrankEcalHits(cluster, track, _maxEcalLayer, _bField);
        timeOfFlight = getTofFrankFit(frankHits, track, _timeResolution);
      }
    } else {
      // define tof as an average time between two SET strips
      // if no SET hits found, tof alreasy is 0, just skip
      TrackerHit* hitSET = getSETHit(track);

      if (hitSET != nullptr) {
        const std::vector<LCObject*>& rawObjects = hitSET->getRawHits();
        if (rawObjects.empty()) {
          streamlog_out(WARNING) << "Found no raw SET strip hits, but space point is built!? Writing TOF as 0."
                                 << std::endl;
        } else if (rawObjects.size() == 1) {
          streamlog_out(ERROR)
              << "Found only one SET strip hit, but space point is built!? Writing TOF from a single strip."
              << std::endl;
          auto* rawHit = dynamic_cast<TrackerHitPlane*>(rawObjects[0]);
          if (rawHit != nullptr)
            timeOfFlight = RandGauss::shoot(rawHit->getTime(), _timeResolution);
        } else {
          if (rawObjects.size() > 2)
            streamlog_out(WARNING)
                << "Found more than two SET strip hits! Writing TOF as an average of the first two elements."
                << std::endl;
          auto* rawHitFront = dynamic_cast<TrackerHitPlane*>(rawObjects[0]);
          auto* rawHitBack = dynamic_cast<TrackerHitPlane*>(rawObjects[1]);
          if (rawHitFront != nullptr && rawHitBack != nullptr) {
            double timeFront = RandGauss::shoot(rawHitFront->getTime(), _timeResolution);
            double timeBack = RandGauss::shoot(rawHitBack->getTime(), _timeResolution);
            timeOfFlight = (timeFront + timeBack) / 2.;
          }
        }
      }
    }
    vector<float> results{float(timeOfFlight)};
    pidHandler.setParticleID(pfo, 0, 0, 0., algoID, results);
    streamlog_out(DEBUG9) << "Final results for the " << i + 1 << " PFO" << std::endl;
    streamlog_out(DEBUG9) << "time-of-flight: " << float(timeOfFlight) << " ns" << std::endl;
    streamlog_out(DEBUG9) << std::endl << std::endl;
  }
}
