#ifndef TOFEstimators_h
#define TOFEstimators_h 1

/**
\author F. Gaede, DESY, April 2018
\author B. Dudar, DESY, September 2021
*/

#include <string>
#include <vector>
#include "marlin/Processor.h"
#include "MarlinTrk/IMarlinTrkSystem.h"

class TOFEstimators : public marlin::Processor {
    public:
        /**
        Copy constructor is removed to avoid W-effc++ warnings.
        */
        TOFEstimators(const TOFEstimators&) = delete;

        /**
        Copy operator is removed to avoid W-effc++ warnings.
        */
        TOFEstimators& operator=(const TOFEstimators&) = delete;

        /**
        Method required by the Marlin to register processor in the global scope.
        */
        marlin::Processor* newProcessor() { return new TOFEstimators; }

        /**
        Default constructor.
        Defines steering parameters from the xml steering file.
        */
        TOFEstimators();

        /** Called at the begin of the job before anything is read.
        Extracts geometry details and initializes Kalman Filter System.
        */
        void init();

        /** Called for every event - the working horse.
        Calculates momentum, track length and time of flight and
        writes them into PIDHandler of the input collection object
        */
        void processEvent(EVENT::LCEvent* evt);

    private:
        /** Steering parameter: name of the input collection of ReconstructedParticles to analyse.
        Usually PandoraPFOs
        */
        std::string _pfoCollectionName{};

        /**  Steering parameter: an option to indicate whether to do the calculations.
        to the last tracker hit or extrapolate to the ECal surface.
        */
        bool _extrapolateToEcal{};

        /**  Steering parameter: name of the method to calculate time of flight based on the ECal hits.
        If _extrapolateToEcal == false then ignored
        */
        std::string _tofMethod{};

        /** Steering parameter: Assumed time resolution of the detector elements in ps.
        In case _extrapolateToEcal == false, then smearing applied for both SET strips.
        In case _extrapolateToEcal == true, then smearing applied for every ECal hit.
        */
        double _timeResolution{};

        /** Steering parameter: Number of inner ECal layers to consider for TOF methods.
        In case _extrapolateToEcal == true and (_tofMethod == "frankAvg" or _tofMethod == "frankFit").
        This parameter indicates how many ECal layers these tof Methods will consider to estimate TOF.
        */
        int _maxEcalLayer{};


        /** Number of the current event.
        */
        int _nEvent{};

        /** Names of the output parameters written to the PIDHandler.
        These parameters are: momentumHM, trackLength and timeOfFlight
        */
        std::vector<std::string> _outputParNames{};


        /** Kalman Filter System for refitting.
        Release notes of MarlinTrk v02-00:
        USERS SHOULD NO LONGER DELETE THE IMarlinTrkSystem POINTER IN THEIR CODE
        */
        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;

        /** B_{z} at the origin in Tesla.
        */
        double _bField{};

        /** TPC outer radius in mm.
        */
        double _tpcOuterR{};
};

#endif
