#ifndef TOFEstimators_h
#define TOFEstimators_h 1

/**
Marlin processor than calculates output parameters.
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
        Copy constructor .
        We remove it to avoid W-effc++ warnings.
        Copying objects with pointer members is a bad idea.
        */
        TOFEstimators(const TOFEstimators&) = delete;

        /**
        Copy assignment operator.
        We remove it to avoid W-effc++ warnings.
        Copying objects with pointer members is a bad idea.
        */
        TOFEstimators& operator=(const TOFEstimators&) = delete;

        /**
        Method required by the Marlin to register processor in the global scope.
        */
        marlin::Processor* newProcessor() { return new TOFEstimators; }

        /**
        Default constructor.
        Registers steering parameters from the xml steering file.
        */
        TOFEstimators();

        /** Called at the begin of the job before anything is read.
        Extracts geometry details and initializes Kalman Filter System.
        */
        void init();

        /** Called for every event - the working horse.
        Calculates momentum, track length and time of flight and
        writes them into PIDHandler of the input collection object.
        */
        void processEvent(EVENT::LCEvent* evt);

    private:
        /** Stores ReconstructedParticleCollection steering parameter.
        */
        std::string _pfoCollectionName{};

        /** Stores ExtrapolateToEcal steering parameter.
        */
        bool _extrapolateToEcal{};

        /** Stores TofMethod steering parameter.
        */
        std::string _tofMethod{};

        /** Stores TimeResolution steering parameter.
        */
        double _timeResolution{};

        /** Stores MaxEcalLayer steering parameter.
        */
        int _maxEcalLayer{};

        /** Stores current event number.
        */
        int _nEvent{};

        /** Stores names of the output parameters.
        These are "momentumHM", "trackLength" and "timeOfFlight".
        */
        std::vector<std::string> _outputParNames{};


        /** Kalman Filter System.
        \note Release notes of MarlinTrk v02-00:
        \note users should no longer delete the IMarlinTrkSystem pointer in their code
        */
        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;

        /** Stores z component of the magnetic field at the origin in Tesla.
        */
        double _bField{};

        /** Stores outer TPC radius in mm.
        */
        double _tpcOuterR{};
};

#endif
