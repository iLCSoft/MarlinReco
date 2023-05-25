#ifndef TOFUtils_h
#define TOFUtils_h 1

#include <vector>
#include "EVENT/ReconstructedParticle.h"
#include "DDRec/Vector3D.h"

/**
 * Utility functions that are used by the TOFEstimators processor.
 *
 * \author F. Gaede, DESY, 2018
 * \author B. Dudar, DESY, 2022
*/
namespace TOFUtils{

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

}

#endif
