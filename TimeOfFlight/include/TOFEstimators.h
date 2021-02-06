#ifndef TOFEstimators_h
#define TOFEstimators_h 1

#include "marlin/Processor.h"
#include "lcio.h"

//std stuff
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <functional>
#include <memory>
using std::string, std::vector, std::endl;

#include <gsl/gsl_rng.h>

#include <EVENT/Track.h>
#include <EVENT/Cluster.h>

#include "DDRec/Vector3D.h"
using dd4hep::rec::Vector3D;

using namespace lcio ;
using namespace marlin ;

#include "TH2F.h"
// class TH2F ;

/** Compute estimators for the time of flight from the CalorimeterHits in Clusters.
 *  The estimators are stored in a PID object with the name of the processor.
 *  The following parameters are stored:
 *   TOFFirstHit         -  first hit ( closest to calo entry point)
 *   TOFClosestHits      -  closest hit in every layer (<lMax)
 *   TOFClosestHitsError -  error of above
 *   TOFFlightLength     -  trajectory length to reference point
 *   TOFLastTrkHit       -  (unsmeared) time of last tracker hit (in SET)
 *   TOFLastTrkHitFlightLength -  trajectory length to last trk hit
 *
 *  processor parameters:
 *  MaxLayerNumber  - restrict the hits to the first MaxLayerNumbers
 *  ReconstructedParticleCollection - input collection
 *  TimeResolution  - assumed single hit resolution in ps
 *
 *  Only Ecal hits are considered.
 *  For charged particles the TrackState at the calorimeter is used as reference point and for the
 *  extrapolation into the calorimeter. For neutral particles the position of the first hit (the one closest to the IP)
 *  is used and a straight flight path from the IP is taken for the extrapolation.
 *  The hit time is corrected for the flight time from the reference assuming speed of light.
 *
 *  See ../scripts/tofestimators.xml for example steering file.
 *
 * @author F.Gaede, DESY, April 2018
 * @author B.Dudar, DESY, February 2021
 */

class TOFEstimators : public Processor {
    public:

    virtual Processor*  newProcessor() { return new TOFEstimators ; }

    TOFEstimators(const TOFEstimators&) = delete;
    TOFEstimators& operator=(const TOFEstimators&) = delete;

    TOFEstimators() ;

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    virtual void init() ;

    /** Called for every event - the working horse.
    */
    virtual void processEvent(LCEvent* evt) ;

    virtual void check(LCEvent* evt) ;

    /** Called after data processing for clean up.
    */
    virtual void end() ;

    Vector3D getMomAtCalo(const Track*);
    double getFlightLength(const Track*);
    double getTOFClosest(const Track*, const Cluster*);
    double getTOFFastest(const Track*, const Cluster*);
    double getTOFCylFit(const Track*, const Cluster*);
    double getTOFClosestFit(const Track*, const Cluster*);

    protected:

    /** Input collection name with ReconstructedParticles (PFOs).
    */
    string _colNamePFO{};

    string _procVersion{};
    int _maxLayerNum{};
    float _resolution{};
    vector<string> _TOFNames{};

    gsl_rng* _rng = nullptr;
    vector<TH2F*> _h{};

    std::default_random_engine _generator{};
    std::normal_distribution<double> _smearing{};
};

#endif
