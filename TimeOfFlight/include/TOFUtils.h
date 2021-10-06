#ifndef TOFUtils_h
#define TOFUtils_h 1

#include <vector>
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "DDRec/Vector3D.h"

/**
 * Utility functions that are used by the TOFEstimators processor.
 *
 * \author F. Gaede, DESY, 2018
 * \author B. Dudar, DESY, 2021
*/
namespace TOFUtils{

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


    /** Get track momentum at the track state.
    Returns momentum Vector3D from the given track state.
    */
    dd4hep::rec::Vector3D getHelixMomAtTrackState(const EVENT::TrackState& ts, double bField);


    /** Get track length.
    Returns the track length between two track states estimated by the helix length formula:

    \f$ \ell = \sqrt{\left( \frac{\varphi_{i+1} - \varphi_{i}}{\Omega}\right)^{2} + \left( z_{i+1} - z_{i} \right)^{2} } \f$

    Note: The formula above works only for the arcs with \f$ \Delta \varphi < \pi \f$.
    */
    double getHelixArcLength(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);



    /** Get track length.
    Returns the track length between two track states estimated by the helix length formula:

    \f$ \ell = \frac{\left |z_{i+1} - z_{i}\right |}{\tan{\lambda}} \sqrt{1 + \tan{\lambda}^{2} } \f$

    Note: The formula above works for any \f$ \Delta \varphi \f$.

    However it is less precise than getHelixArcLength() due to the less precise \f$ \tan{\lambda} \f$.
    Also helix formula implies constant momentum assumption which may show higher discrepancy for long tracks with low momentum.
    */
    double getHelixLengthAlongZ(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);


    /** Get number of helix revolutions.
    Returns number of helix revolutions between two track states.

    The calculation is done with:
    \f$  N_{\mathrm{turns}} = \frac{\left |z_{i+1} - z_{i}\right |}{\tan{\lambda}} \bigg / (2 \pi \frac{1}{|\Omega|}) \f$
    */
    double getHelixNRevolutions(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    /** Returns TPC outer radius from the DD4hep detector geometry.
    */
    double getTPCOuterR();

    /** Returns SET hit.
    Checks the \f$ \rho = \sqrt{x^{2}+y^{2}} \f$ of the tracker hit.
    If \f$ \rho > R_{\mathrm{TPC, outer}} \f$  then it is the SET hit.
    */
    EVENT::TrackerHit* getSETHit(EVENT::Track* track, double tpcOuterR);

    /** Get a subset of the cluster calorimeter hits.
    Returns the subset of the cluster calorimeter hits which we call "Frank"
    by historical reasons as these hits were used to calculate time-of-flight as
    a default method for the IDR production.

    "Frank" hits are the closest ECal hits to the linearly extrapolated
    track line inside the ECal in each of the first maxEcalLayer ECal layers.
    */
    std::vector<EVENT::CalorimeterHit*> selectFrankEcalHits( EVENT::Cluster* cluster, EVENT::Track* track, int maxEcalLayer, double bField );


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
    std::vector<IMPL::TrackStateImpl> getTrackStatesPerHit(std::vector<EVENT::Track*> tracks, MarlinTrk::IMarlinTrkSystem* trkSystem, bool extrapolateToEcal, double bField);

    /** Get the time-of-flight using the closest ECal hit.
    Returns time measured by the closest ECal hit to the extrapolated track position at the ECal surface.

    \f$ \mathrm{TOF} = t_{\mathrm{closest}} - \frac{\left| \vec{r}_{\mathrm{track}} - \vec{r}_{\mathrm{closest}} \right|}{c} \f$

    If no ECal hits are found returns `0.0`.
    */
    double getTofClosest( EVENT::Cluster* cluster, EVENT::Track* track, double timeResolution);

    /** Get the time-of-flight using average of the Frank ECal hits.
    Returns the time-of-flight as an average of the Frank hits time corrected for the distance to
    the extrapolated track position at the ECal surface assuming speed of flight is the speed of light.

    \f$ \mathrm{TOF} = \frac{1}{\mathrm{MaxEcalLayer}}\sum_{i}^{\mathrm{MaxEcalLayer}} \left( t_{i} - \frac{\left|\vec{r}_{\mathrm{track}} - \vec{r}_{i} \right|}{c} \right) \f$

    If no ECal hits are found within selectedHits, then returns `0.0`.
    */
    double getTofFrankAvg( std::vector<EVENT::CalorimeterHit*> selectedHits, EVENT::Track* track, double timeResolution);

    /** Get the time-of-flight using fit of the Frank ECal hits.
    Returns the  time-of-flight as an extrapolated time at the extrapolated track
    position at the ECal surface by using a linear fit of the Frank hits time as a
    function of the distance to the extrapolated track position at the ECal surface.

    \f$ \mathrm{TOF}=f(|\vec{r}_{\mathrm{track}} - \vec{r}_{\mathrm{hit}} |=0) \f$

    If no ECal hits are found within selectedHits, then returns `0.0`.
    If *only one* ECal hit is found, which is not enough to perform the linear fit, then returns the same as getTofClosest() and getTofFrankAvg().
    */
    double getTofFrankFit( std::vector<EVENT::CalorimeterHit*> selectedHits, EVENT::Track* track, double timeResolution);

    /** Print debug info.
    Prints a current usage of the virtual memory (VM) and resident set size (RSS).
    */
    void debugPrint();

}



#endif
