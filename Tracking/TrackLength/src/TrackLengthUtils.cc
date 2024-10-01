#include "TrackLengthUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "HelixClass.h"
#include "IMPL/TrackImpl.h"
#include "MarlinTrk/MarlinTrkUtils.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "UTIL/ILDConf.h"
#include "marlin/VerbosityLevels.h"
#include "marlinutil/CalorimeterHitType.h"

using CLHEP::RandGauss;
using dd4hep::rec::Vector3D;
using EVENT::CalorimeterHit;
using EVENT::Track;
using EVENT::TrackerHit;
using EVENT::TrackState;
using IMPL::TrackImpl;
using IMPL::TrackStateImpl;
using MarlinTrk::IMarlinTrack;
using MarlinTrk::IMarlinTrkSystem;
using std::numeric_limits;
using std::pair;
using std::vector;
using UTIL::ILDDetID;

bool TrackLengthUtils::sortByRho(EVENT::TrackerHit* a, EVENT::TrackerHit* b) {
  Vector3D posA(a->getPosition()), posB(b->getPosition());
  return posA.rho() < posB.rho();
}

IMPL::TrackStateImpl TrackLengthUtils::getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrk, EVENT::TrackerHit* hit) {
  TrackStateImpl ts;
  double chi2Dummy;
  int ndfDummy;
  marlinTrk->getTrackState(hit, ts, chi2Dummy, ndfDummy);
  return ts;
}

double TrackLengthUtils::getHelixLength(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2) {
  double tanL = ts1.getTanLambda();
  double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
  double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();
  return std::abs((z2 - z1) / tanL) * std::sqrt(1. + tanL * tanL);
}


double TrackLengthUtils::getHelixLengthOption1(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
    // Suboptimal and not used. Use getHelixLength()
    double omega = std::abs( ts1.getOmega() );
    double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
    double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();
    double dz = std::abs(z2 - z1);
    double dPhi = std::abs( ts2.getPhi() - ts1.getPhi() );
    // We are never sure whether the track indeed curled for more than pi
    // or it was just a crossing of the singularity point (-pi/+pi) in the phi coordinate system.
    // We always assume it was the crossing of the (-pi/+pi) and correct for this.
    // Thus this formula only applicable for small distances between the track states or for the non-curly tracks (dPhi < pi).
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
    return std::sqrt(dPhi*dPhi/(omega*omega)+dz*dz);
};

double TrackLengthUtils::getHelixLengthOption2(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
    // Suboptimal and not used. Use getHelixLength()
    double omega = std::abs( ts1.getOmega() );
    double tanL = std::abs( ts1.getTanLambda() );
    double dPhi = std::abs( ts2.getPhi() - ts1.getPhi() );
    // We are never sure whether the track indeed curled for more than pi
    // or it was just a crossing of the singularity point (-pi/+pi) in the phi coordinate system.
    // We always assume it was the crossing of the (-pi/+pi) and correct for this.
    // Thus this formula only applicable for small distances between the track states or for the non-curly tracks (dPhi < pi).
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
    return dPhi/omega*std::sqrt(1+tanL*tanL);
};

std::vector<EVENT::Track*> TrackLengthUtils::getSubTracks(EVENT::Track* track){
    vector<Track*> subTracks;
    // OPTIMIZE: Remove geometry dependency.
    // The geometry dependency comes from the ILD tracking.
    // The way subtracks are stored inside the Track object is special to ILD tracking.
    // The geometry dependency can be easily removed, but at the cost of a few vertex hits sometimes missing with ILD tracking...
    // Ideally this function should be just: return track->getTracks(), if the subtracks stored consistently in ILD tracking...

    // add track itself, which contains VXD+FTD+SIT+TPC hits of the first curl.
    subTracks.push_back(track);

  int nSubTracks = track->getTracks().size();
  if (nSubTracks <= 1)
    return subTracks;

  UTIL::BitField64 encoder(UTIL::LCTrackerCellID::encoding_string());
  auto isTPCHit = [&encoder](TrackerHit* hit) -> bool {
    encoder.setValue(hit->getCellID0());
    int subdet = encoder[UTIL::LCTrackerCellID::subdet()];
    return subdet == UTIL::ILDDetID::TPC;
  };

  int indexOfFirstTPCCurl = 0;
  for (int i = 0; i < nSubTracks; ++i) {
    Track* subTrack = track->getTracks()[i];
    auto hits = subTrack->getTrackerHits();
    if (std::find_if(hits.begin(), hits.end(), isTPCHit) != hits.end()) {
      indexOfFirstTPCCurl = i;
      break;
    }
  }

  for (int j = indexOfFirstTPCCurl + 1; j < nSubTracks; ++j)
    subTracks.push_back(track->getTracks()[j]);
  return subTracks;
}

std::vector<IMPL::TrackStateImpl> TrackLengthUtils::getTrackStates(EVENT::ReconstructedParticle* pfo,
                                                                   MarlinTrk::IMarlinTrkSystem* trkSystem,
                                                                   double bField) {
  // Refit the track and extract track state at every tracker hit along the track
  vector<TrackStateImpl> trackStates;
  if (pfo->getTracks().empty())
    return trackStates;
  vector<Track*> subTracks = getSubTracks(pfo->getTracks()[0]);

  TrackImpl lastGoodRefittedTrack;

  streamlog_out(DEBUG8) << "PFOs track has " << subTracks.size() << " subTracks." << std::endl;
  for (size_t i = 0; i < subTracks.size(); ++i) {
    vector<TrackerHit*> hits = subTracks[i]->getTrackerHits();

    streamlog_out(DEBUG8) << "Subtrack " << i + 1 << " has " << hits.size() << " hits." << std::endl;
    std::sort(hits.begin(), hits.end(), sortByRho);

    // setup initial dummy covariance matrix
    vector<float> covMatrix(15);
    // initialize variances
    covMatrix[0] = 1e+06;   // sigma_d0^2
    covMatrix[2] = 100.;    // sigma_phi0^2
    covMatrix[5] = 0.00001; // sigma_omega^2
    covMatrix[9] = 1e+06;   // sigma_z0^2
    covMatrix[14] = 100.;   // sigma_tanl^2
    double maxChi2PerHit = 100.;
    std::unique_ptr<IMarlinTrack> marlinTrk(trkSystem->createTrack());
    TrackImpl refittedTrack;

    // Need to initialize trackState at last hit
    TrackStateImpl preFit = *(subTracks[i]->getTrackState(TrackState::AtLastHit));
    preFit.setCovMatrix(covMatrix);
    int errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::backward,
                                                       &preFit, bField, maxChi2PerHit);
    // if fit fails, try also fit forward
    if (errorFit != 0) {
      streamlog_out(DEBUG8) << "Fit backward fails! Trying to fit forward for " << i + 1 << " subTrack in this PFO!"
                            << std::endl;
      marlinTrk.reset(trkSystem->createTrack());
      errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::forward,
                                                     &preFit, bField, maxChi2PerHit);
    }
    if (errorFit != 0) {
      streamlog_out(WARNING) << "Fit fails in both directions. Skipping " << i + 1 << " subTrack in this PFO!"
                             << std::endl;
      continue;
    }
    lastGoodRefittedTrack = refittedTrack;

    // here hits are sorted by rho=(x^2+y^2) in the fit direction. forward - increasing rho, backward - decreasing rho
    vector<pair<TrackerHit*, double>> hitsInFit;
    marlinTrk->getHitsInFit(hitsInFit);

    // Find which way to loop over the array of hits. We need to loop in the direction of the track.
    bool loopForward = true;
    double zFirst = std::abs(hitsInFit.front().first->getPosition()[2]);
    double zLast = std::abs(hitsInFit.back().first->getPosition()[2]);

    // OPTIMIZE: 10 mm is just a round number. With very small z difference it is more robust to use rho, to be sure z
    // difference is not caused by tpc Z resolution or multiple scattering
    if (std::abs(zLast - zFirst) > 10.) {
      if (zLast < zFirst)
        loopForward = false;
      streamlog_out(DEBUG8) << "Using Z to define loop direction over subTrack hits." << std::endl;
      streamlog_out(DEBUG8) << "subTrack " << i + 1 << " zFirst: " << hitsInFit.front().first->getPosition()[2]
                            << " zLast: " << hitsInFit.back().first->getPosition()[2]
                            << " loop forward: " << loopForward << std::endl;
    } else {
      double rhoFirst =
          std::hypot(hitsInFit.front().first->getPosition()[0], hitsInFit.front().first->getPosition()[1]);
      double rhoLast = std::hypot(hitsInFit.back().first->getPosition()[0], hitsInFit.back().first->getPosition()[1]);
      if (rhoLast < rhoFirst)
        loopForward = false;
      streamlog_out(DEBUG8)
          << "Track is very perpendicular (dz < 10 mm). Using rho to define loop direction over subTrack hits."
          << std::endl;
      streamlog_out(DEBUG8) << "subTrack " << i + 1 << " zFirst: " << hitsInFit.front().first->getPosition()[2]
                            << " zLast: " << hitsInFit.back().first->getPosition()[2] << std::endl;
      streamlog_out(DEBUG8) << "subTrack " << i + 1 << " rhoFirst: " << rhoFirst << " rhoLast: " << rhoLast
                            << " loop forward: " << loopForward << std::endl;
    }

    int nHitsInFit = hitsInFit.size();
    // if first successfully fitted subTrack add IP track state
    if (trackStates.empty())
      trackStates.push_back(*(static_cast<const TrackStateImpl*>(refittedTrack.getTrackState(TrackState::AtIP))));

    // NOTE: although we use z to understand subTrack's direction, subTrack's hits are still sorted by rho
    if (loopForward) {
      for (int j = 0; j < nHitsInFit; ++j) {
        TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
        trackStates.push_back(ts);
      }
    } else {
      for (int j = nHitsInFit - 1; j >= 0; --j) {
        TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
        trackStates.push_back(ts);
      }
    }
  }

    const TrackStateImpl* tsCaloBugged = static_cast<const TrackStateImpl*> (lastGoodRefittedTrack.getTrackState(TrackState::AtCalorimeter) );
    if ( pfo->getClusters().size() > 0 && tsCaloBugged != nullptr ){
        // d0 and z0 of the track state at calo MUST be 0. This is a bug! Fix manually here for old produced MC samples.
        IMPL::TrackStateImpl tsCalo = *(dynamic_cast<const IMPL::TrackStateImpl*> (tsCaloBugged));
        tsCalo.setD0(0.);
        tsCalo.setZ0(0.);
        trackStates.push_back( tsCalo );
    }
    return trackStates;
}

double TrackLengthUtils::getHelixNRevolutions(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2) {
  double omega = std::abs(ts1.getOmega());
  double tanL = std::abs(ts1.getTanLambda());
  double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
  double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();
  double dz = std::abs(z2 - z1);
  // helix length projected on xy
  double circHelix = dz / tanL;
  double circFull = 2 * M_PI / omega;
  return circHelix / circFull;
}
