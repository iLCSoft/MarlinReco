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
#include <limits>

#include <gsl/gsl_rng.h>

#include <EVENT/Track.h>
#include <EVENT/Cluster.h>

#include "DDRec/Vector3D.h"

using namespace lcio ;
using namespace marlin ;

#include "TH2F.h"

/**
 * @section DESCRIPTION
 * Compute time of flight (ToF) from the CalorimeterHits in the ECAL Cluster.
 * The ToFs are stored in a PID object with the name of the processor attached
 * to the PandoraPFOs.
 * Steering parameters:
 *  MaxLayerNumber  - Estimate ToF using only MaxLayerNumber first layers in the ECAL
 *  ReconstructedParticleCollection - input collection (usually PandoraPFOs)
 *  TimeResolution  - assumed single Calorimeter hit time resolution in ps
 *  CylRadius - radius within which hits are accepted for the fit. Relevant only
                for new TOFCylFit algorithm. Default 5 mm
 *  ProcessorVersion - controls output of the processor between idr
 *                (interim design report) and dev (current development) versions
 * Output for ProcessorVersion = "idr"
 *  The following parameters are stored with the IDR ProcessorVersion:
 *   TOFFirstHit         -  ToF of the CalorimeterHit closest to track entry
 *                          point to the ECAL corrected with the distance
 *                          between cell center and entry point assuming
 *                          speed of light propagation
 *   TOFClosestHits      -  ToF estimated as the average of corrected
 *                          CalorimeterHit times from hits that are closest to
 *                          the linear prolongation of the trajectory line
 *                          of the track in MaxLayerNumber first layers
 *   TOFClosestHitsError -  error of the TOFClosestHits
 *   TOFFlightLength     -  Track length of the helix between IP and ECAL entrance
 *                          assuming IP at (0, 0, 0) and ECAL entry point is
 *                          obtained from Kalman filter by extrapolation.
 *                          Omega/lambda track parameters are taken from IP hit
 *                          and assumed to be constant along the track.
 *                          Results are reasonable for delta phi < 2pi
 *                          meaning not more than one curl.
 *   TOFLastTrkHit       -  (unsmeared) time of last tracker hit (in the SET)
 *   TOFLastTrkHitFlightLength - track length similar to the TOFFlightLength
 *                               but to the SET hit
 * Notes:
 *  For charged particles the TrackState at the calorimeter is used as reference point and for the
 *  extrapolation into the calorimeter. For neutral particles the position of the first hit (the one closest to the IP)
 *  is used and a straight flight path from the IP is taken for the extrapolation.
 *  The hit time is corrected for the flight time from the reference assuming speed of light.
 *
 * Output for ProcessorVersion = "dev"
 *  The following parameters are stored with the IDR ProcessorVersion:
 *   TOFClosest         -  ToF of the CalorimeterHit closest to track entry
 *                          point to the ECAL corrected with the distance
 *                          between cell center and entry point assuming
 *                          speed of light propagation
 *   TOFFastest         -  ToF of the CalorimeterHit arrived fastest corrected
 *                          with the distance between cell center and entry
 *                          point assuming speed of light propagation
 *   TOFCylFit          -  ToF estimated from extrapolation from a fit of the
 *                            hit times vs distance to to the track entry point
 *                            to the ECAL at distance = 0
 *   TOFClosestFit    -  ToF estimated with a fit as TOFCylFit but from hits
 *                            that are closest to the linear prolongation of the
 *                            trajectory line of the track in MaxLayerNumber first layers
 *   FlightLength       -  Track length of the helix between IP and ECAL entrance
 *                          assuming IP at (0, 0, 0) and ECAL entry point is
 *                          obtained from Kalman filter by extrapolation.
 *                          Omega/lambda track parameters are taken from the
 *                            closest track hit to the ECAL and assumed to be
 *                            constant along the track.
 *                          Results are reasonable for delta phi < 2pi
 *                          meaning not more than one curl.
 *   MomAtCalo - momentum of the track estimated from the track curvature near
 *                 ECAL entry point.
 * Notes:
 * FlightLength shows less bias in mass than TOFFlightLength
 * MomAtCalo shows less bias in mass than taking momentum from IP as in IDR
 * Every method in both idr and dev considers only Ecal hits.
 * Methods work only for charged PFOs with ONE associated track and ONE associated cluster
 * Therefore neutral PFOs and PFOs with multiple tracks/clusters are ignored
 *
 *  See ../scripts/tofestimators.xml for example steering file.
 *
 * @file
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

    /** Called after data processing to plot some histograms to debug.
    */
    virtual void check(LCEvent* evt) ;

    /** Called after data processing for clean up.
    */
    virtual void end() ;

    /** Called in dev version to write MomAtCalo parameter.
    Can be used for more precise mass calculation
    */
    dd4hep::rec::Vector3D getMomAtCalo(const Track*);

    /** Track length of the helix between IP and ECAL entrance
    * assuming IP at (0, 0, 0) and ECAL entry point is
    * obtained from the Kalman filter by extrapolation.
    * Omega/lambda track parameters are taken from the
    * closest track hit to the ECAL and assumed to be
    * constant along the track.
    * Length is calculated by the formula that is reasonable only for delta phi < 2pi
    * meaning only tracks with less than one curl.
    */
    double getFlightLength(const Track*);

    /**
    *   TOFClosest         -  ToF of the CalorimeterHit closest to track entry
    *                          point to the ECAL corrected with the distance
    *                          between cell center and entry point assuming
    *                          speed of light propagation
    */
    double getTOFClosest(const Track*, const Cluster*);

    /**
    *   TOFFastest         -  ToF of the CalorimeterHit arrived fastest corrected
    *                          with the distance between cell center and entry
    *                          point assuming speed of light propagation
    */
    double getTOFFastest(const Track*, const Cluster*);

    /**
    *   TOFCylFit          -  ToF estimated from extrapolation from a fit of the
    *                            hit times vs distance to to the track entry point
    *                            to the ECAL at distance = 0
    */
    double getTOFCylFit(const Track*, const Cluster*);

    /**
    *   TOFClosestFit    -  ToF estimated with a fit as TOFCylFit but from hits
    *                            that are closest to the linear prolongation of the
    *                            trajectory line of the track in MaxLayerNumber first layers
    */
    double getTOFClosestFit(const Track*, const Cluster*);

    protected:

    /** Input collection name with ReconstructedParticles (PFOs).
        used for ReconstructedParticleCollection steering parameter. usually PandoraPFOs
    */
    std::string _colNamePFO{};

    /** Input version idr (default) or dev (current development).
    */
    std::string _procVersion{};

    /** Select ECAL hits only from _maxLayerNum first layers of ECAL for ToF calculations
    */
    int _maxLayerNum{};


    /** Select hits for tofCylFit method from the cylinder limited by this radius
    */
    double _cylRadiusCut{};

    /** Single ECAL hit time resolution in ps
    */
    float _resolution{};

    /** vector of names of output PID object parameters which are written as output
    */
    std::vector<std::string> _TOFNames{};

    gsl_rng* _rng = nullptr;
    std::vector<TH2F*> _h{};

    std::mt19937 _generator{};
    std::normal_distribution<double> _smearing{};

    /** B field from the DD4HEP geometry processor for MomAtCalo calculations
    */
    double _bField[3];
};

#endif
