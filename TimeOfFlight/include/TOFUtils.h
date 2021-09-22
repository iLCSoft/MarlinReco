#ifndef TOFUtils_h
#define TOFUtils_h 1

/**
 * Utility functions that are used by the TOFEstimators processor.
 *
 * \author F. Gaede, DESY, 2018
 * \author B. Dudar, DESY, 2021
 */

#include <vector>
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "DDRec/Vector3D.h"

namespace TOFUtils{

    /** Comparator function for tracker hits.

    Returns `true` if \f$ \rho_{a} < \rho_{b} \f$.
    With \f$ \rho \f$ being a radius of the tracker hit projected on the \f$ xy \f$ plane: \f$ \rho = \sqrt{x^2 + y^2} \f$.

    Primarily used to sort vector of tracker hits by \f$ \rho \f$ for the Kalman Filter fit.
    */
    bool sortByRho(EVENT::TrackerHit* a, EVENT::TrackerHit* b);


    /** Extracts track state at the given tracker hit used by Kalman Filter in the refit.
    */
    IMPL::TrackStateImpl getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrk, EVENT::TrackerHit* hit);


    /** Returns momentum vector from the given track state.
    */
    dd4hep::rec::Vector3D getHelixMomAtTrackState(const EVENT::TrackState& ts, double bField);


    /** Estimate helix arc length between two track states.

    Works only for arcs with \Delta \varphi < \pi
    Here is the formula it uses

    */
    double getHelixArcLength(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);


    /** Estimate helix arc length between two track states without \varphi.

    Relies on tanL. Less precise than the getHelixMomAtTrackState.
    Used to calculate track length with \Delta \varphi > \pi

    Here is the formula it uses

    */
    double getHelixLengthAlongZ(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);


    /** Get number of helix revolutions between two track states assuming constant momentum.
    */
    double getHelixNRevolutions(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    /** Return TPC outer radius from DD4hep detector geometry.
    */
    double getTPCOuterR();

    /** Return SET hit.
    SET hit is assumed to be a tracker hit with a radius larger than TPC outer radius
    */
    EVENT::TrackerHit* getSETHit(EVENT::Track* track, double tpcOuterR);

    /** TEST.
    SET hit is assumed to be a tracker hit with a radius larger than TPC outer radius
    */
    std::vector<EVENT::CalorimeterHit*> selectFrankEcalHits( EVENT::Cluster* cluster, EVENT::Track* track, int maxEcalLayer, double bField );


    /** TEST.
    SET hit is assumed to be a tracker hit with a radius larger than TPC outer radius

    This is not guarantied to work 100% of times,
    but I didn't find a better way to collect subTracks...
    1) Always add main track as a first subtrack
    2) If nTPCHits == nHits of SubTrack0 +-1 assume subTrack0 stores TPC hits
    Deviation on 1 hit may happen because of the SET and only God knows what other reasons..
    This would mean that VXD+SIT subTrack is not stored in the track!
    So we skip initial subTrack0 with TPC hits which we have anyhow added with global Track
    and loop over remaining subTracks
    3) If 2) is false and nTPCHits == nHits of SubTrack1 +-1 assume subTrack1 stores TPC hits
    This would mean that subTrack0 stores VXD+SIT hits
    So we skip both VXD+SIT and 1st TPC subTracks which we have anyhow added with global Track
    and loop over remaining subTracks

    */
    std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track);

    /** test.
     return vector of trackStates for the IP, every hit in the whole track and ECAL if extrapolate to the ECAL option is chosen
    track state are sorted along the helix
    */
    std::vector<IMPL::TrackStateImpl> getTrackStatesPerHit(std::vector<EVENT::Track*> tracks, MarlinTrk::IMarlinTrkSystem* trkSystem, bool extrapolateToEcal, double bField);

    /** test.
    */
    double getTofClosest( EVENT::Cluster* cluster, EVENT::Track* track, double timeResolution);

    /** test.
    */
    double getTofFrankAvg( std::vector<EVENT::CalorimeterHit*> selectedHits, EVENT::Track* track, double timeResolution);

    /** test.
    */
    double getTofFrankFit( std::vector<EVENT::CalorimeterHit*> selectedHits, EVENT::Track* track, double timeResolution);


}



#endif
