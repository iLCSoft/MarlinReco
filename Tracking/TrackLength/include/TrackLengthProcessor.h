#ifndef TrackLengthProcessor_h
#define TrackLengthProcessor_h 1

#include "MarlinTrk/IMarlinTrkSystem.h"
#include "marlin/Processor.h"
#include <string>
#include <vector>

/**
Marlin processor that calculates harmonic mean momentum and track length of the track.
\author B. Dudar, DESY, 2022
*/
class TrackLengthProcessor : public marlin::Processor {
public:
  /**
  Copy constructor.
  */
  TrackLengthProcessor(const TrackLengthProcessor&) = delete;

  /**
  Copy assignment operator.
  */
  TrackLengthProcessor& operator=(const TrackLengthProcessor&) = delete;

  /**
  Method required by the Marlin to register processor in the global scope.
  */
  marlin::Processor* newProcessor() { return new TrackLengthProcessor; }

  /**
  Default constructor.
  Registers steering parameters from the xml steering file.
  */
  TrackLengthProcessor();

  /** Called at the begin of the job before anything is read.
  Extracts geometry details and initializes Kalman Filter System.
  */
  void init();

  /** Called for every event - the working horse.
  Calculates momentum and track length and writes them into PIDHandler of the input collection object.
  */
  void processEvent(EVENT::LCEvent* evt);

private:
  /** Stores ReconstructedParticleCollection steering parameter.
   */
  std::string _pfoCollectionName{};

  /** Stores current event number.
   */
  int _nEvent{};

  /** Stores names of the output parameters.
  These are "trackLengthToSET", "trackLengthToEcal", "momentumHMToSET", "momentumHMToEcal".
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
};

#endif
