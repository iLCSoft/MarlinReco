#ifndef TrackLengthUtils_h
#define TrackLengthUtils_h 1

#include <vector>
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "DDRec/Vector3D.h"

/**
 * Utility functions that are used by the TrackLengthProcessor.
 *
 * \author B. Dudar, DESY, 2022
*/
namespace TrackLengthUtils{

    /** Comparator function by radius for tracker hits.
    Returns `true` if \f$ \rho_{a} < \rho_{b} \f$.
    With \f$ \rho \f$ being a radius of the tracker hit projected on the \f$ xy \f$ plane: \f$ \rho = \sqrt{x^2 + y^2} \f$.

    Primarily used to sort vector of tracker hits by \f$ \rho \f$ for the Kalman Filter fit.
    */
    bool sortByRho(EVENT::TrackerHit* a, EVENT::TrackerHit* b);


    /** Get track state at tracker hit.
    Returns track state at tracker hit from the Kalman Filter fit.

    Tracker hit has to be used in the fit by the Kalman Filter.
    */
    IMPL::TrackStateImpl getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrk, EVENT::TrackerHit* hit);


    /** Get track length.
    Returns the track length between two track states estimated by the helix length formula:

    \f$ \ell = \frac{\left |z_{i+1} - z_{i}\right |}{| \tan{\lambda}| } \sqrt{1 + \tan{\lambda}^{2} } \f$
    */
    double getHelixLength(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);


    /** Get all subtracks of the Track.
    Returns a vector of the subTracks of the main track that is passed as an argument.

    The main purpose of this function is to capture all hits of the track.
    This requires iteration over the subTracks as the main Track stores hits only
    of the first half turn due to the our fit procedure.

    The first subTrack in the returned vector is always main Track itself, which
    containes VXD, SIT and TPC tracker hits for the first half turn.

    Then additional subTracks are added which contain hits for additional track revolutions if such exist.

    If nTPCHits \f$ \pm  1\f$ = nHits of SubTrack0 then assume subTrack0 stores TPC hits.
    This would mean that VXD and SIT subTrack is not stored!
    So we need skip *only first* subtrack: initial subTrack0 with TPC hits which we have anyhow added with the main Track
    Else if nTPCHits \f$ \pm  1\f$ = nHits of SubTrack1 assume subTrack1 stores TPC hits
    This would mean that subTrack0 stores VXD and SIT hits.
    So we need to *skip both* subtracks VXD+SIT and 1st TPC half-turn subTracks which we have anyhow added with the main Track

    Note: We consider deviations for \f$\pm 1\f$ hit may happen because of the
    SET and only God knows what other reasons...

    Note: This function is not guarantied to properly work 100% of times,
    but I didn't find a better way to collect subTracks.
    */
    std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track);

    /** Get list of track states.
    Returns a vector of track states at the IP, track state for every tracker hit
    inside all provided subTracks as tracks argument and the track state at the ECal surface if extrapolateToEcal argument is set to true.
    */
    std::vector<IMPL::TrackStateImpl> getTrackStatesPerHit(std::vector<EVENT::Track*> tracks, MarlinTrk::IMarlinTrkSystem* trkSystem, double bField);
}



#endif
