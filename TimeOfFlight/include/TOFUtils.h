#ifndef TOFUtils_h
#define TOFUtils_h 1

/******************************************************
 * Utility functions that are used by the TOFEstimators processor
 *
 * @author F. gaede, DESY, 2018
 * @author B. Dudar, DESY, 2021
 ******************************************************
 */

#include <vector>
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "DDRec/Vector3D.h"

namespace TOFUtils{

    bool sortByRho(EVENT::TrackerHit* a, EVENT::TrackerHit* b);

    IMPL::TrackStateImpl getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrk, EVENT::TrackerHit* hit);


    // gets momentum of the helix from parameters of the trackState
    dd4hep::rec::Vector3D getHelixMomAtTrackState(const EVENT::TrackState& ts, double bField);


    double getHelixArcLength(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    double getHelixLengthAlongZ(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    double getHelixNRevolutions(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    double getTPCOuterR();

    EVENT::TrackerHit* getSETHit(EVENT::Track* track, double tpcOuterR);

    std::vector<EVENT::CalorimeterHit*> selectFrankEcalHits( EVENT::Cluster* cluster, EVENT::Track* track, int maxEcalLayer, double bField );

    // This is not guarantied to work 100% of times,
    // but I didn't find a better way to collect subTracks...
    // 1) Always add main track as a first subtrack
    // 2) If nTPCHits == nHits of SubTrack0 +-1 assume subTrack0 stores TPC hits
    // Deviation on 1 hit may happen because of the SET and only God knows what other reasons..
    // This would mean that VXD+SIT subTrack is not stored in the track!
    // So we skip initial subTrack0 with TPC hits which we have anyhow added with global Track
    // and loop over remaining subTracks
    // 3) If 2) is false and nTPCHits == nHits of SubTrack1 +-1 assume subTrack1 stores TPC hits
    // This would mean that subTrack0 stores VXD+SIT hits
    // So we skip both VXD+SIT and 1st TPC subTracks which we have anyhow added with global Track
    // and loop over remaining subTracks
    std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track);

    // return vector of trackStates for the IP, every hit in the whole track and ECAL if extrapolate to the ECAL option is chosen
    // track state are sorted along the helix
    std::vector<IMPL::TrackStateImpl> getTrackStatesPerHit(std::vector<EVENT::Track*> tracks, MarlinTrk::IMarlinTrkSystem* trkSystem, bool extrapolateToEcal, double bField);

    double getTofClosest( EVENT::Cluster* cluster, EVENT::Track* track, double timeResolution);

    double getTofFrankAvg( std::vector<EVENT::CalorimeterHit*> selectedHits, EVENT::Track* track, double timeResolution);

    double getTofFrankFit( std::vector<EVENT::CalorimeterHit*> selectedHits, EVENT::Track* track, double timeResolution);


}



#endif
