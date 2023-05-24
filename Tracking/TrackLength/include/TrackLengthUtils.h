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

    The main purpose of this function is to capture all hits of the track, not only the first curl which is stored with the main track.
    */
    std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track);

    /** Get list of track states.
    Returns a vector of track states at the IP, track state for every tracker hit of the pfo
    and the track state at the ECal surface.
    */
    std::vector<IMPL::TrackStateImpl> getTrackStates(EVENT::ReconstructedParticle* pfo, MarlinTrk::IMarlinTrkSystem* trkSystem, double bField);
}



#endif
