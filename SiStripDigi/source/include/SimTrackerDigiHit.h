#ifndef SIMTRACKERDIGIHIT_H
#define SIMTRACKERDIGIHIT_H 1

// Include CLHEP header files
#include <CLHEP/Vector/ThreeVector.h>

// Include LCIO header files
#include <EVENT/SimTrackerHit.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <lcio.h>

namespace sistrip {

//! Digitization hit inheritad from LCIO SimTrackerHitImpl, which naturally
//! extends basic features of SimTrackerHitImpl. It defines so-called
//! preStep position (step initial position) and posStep position (step
//! final position).
//!
//! @author Z. Drasal, Charles University Prague
//!

class SimTrackerDigiHit : public IMPL::SimTrackerHitImpl {

public:
  //! Constructor
  SimTrackerDigiHit();

  //! Destructor
  ~SimTrackerDigiHit();

  // SET METHODS

  //! Set preStep position of a hit - DEPRECATED method, necessary for
  //! backwards compatibility
  void setPosition(double pos[3]);

  //! Set preStep position of a hit
  void setPrePosition(double pos[3], float momentum[3], float pathLength);

  //! Set preStep position of a hit
  void setPrePosition(double preX, double preY, double preZ);

  //! Set preStep position Three vector
  void set3PrePosition(const CLHEP::Hep3Vector& prePosition);

  //! Set posStep position of a hit (parameters: preStep position, momentum
  //! at this position and total path length), necessary for backwards
  //! compatibility
  void setPosPosition(double pos[3], float momentum[3], float pathLength);

  //! Set posStep position of a hit
  void setPosPosition(double posX, double posY, double posZ);

  //! Set posStep Three vector
  void set3PosPosition(const CLHEP::Hep3Vector& posPosition);

  //! Set particle momentum at preStep position
  void setMomentum(float p[3]);

  //! Set particle momentum at preStep position
  void setMomentum(float pX, float pY, float pZ);

  //! Set particle Three momentum
  void set3Momentum(const CLHEP::Hep3Vector& momentum);

  //! Set layer ID
  inline void setLayerID(short int iLayer) { _iLayer = iLayer; }

  //! Set ladder ID
  inline void setLadderID(short int iLadder) { _iLadder = iLadder; }

  //! Set sensor ID
  inline void setSensorID(short int iSensor) { _iSensor = iSensor; }

  //! Set original SimTrackerHit
  inline void setSimTrackerHit(EVENT::SimTrackerHit* simHit) { _simHit = simHit; }

  // GET METHODS

  //! Get preStep position Three vector
  inline CLHEP::Hep3Vector get3PrePosition() const { return _prePosition; }

  //! Get preStep position X
  inline double getPreX() const { return _prePosition.getX(); }

  //! Get preStep position Y
  inline double getPreY() const { return _prePosition.getY(); }

  //! Get preStep position Z
  inline double getPreZ() const { return _prePosition.getZ(); }

  //! Get posStep position Three vector
  inline CLHEP::Hep3Vector get3PosPosition() const { return _posPosition; }

  //! Get posStep position X
  inline double getPosX() const { return _posPosition.getX(); }

  //! Get posStep position Y
  inline double getPosY() const { return _posPosition.getY(); }

  //! Get posStep position Z
  inline double getPosZ() const { return _posPosition.getZ(); }

  //! Get step
  inline CLHEP::Hep3Vector get3Step() const { return (_posPosition - _prePosition); }

  //! Get step size
  double getStepSize() const;

  //! Get momentum Three vector
  inline CLHEP::Hep3Vector get3Momentum() const { return _momentum; }

  //! Get momentum X
  inline float getPx() const { return _momentum.getX(); }

  //! Get momentum X
  inline float getPy() const { return _momentum.getY(); }

  //! Get momentum Z
  inline float getPz() const { return _momentum.getZ(); }

  //! Get layer ID
  inline short int getLayerID() const { return _iLayer; }

  //! Get ladder ID
  inline short int getLadderID() const { return _iLadder; }

  //! Get sensor ID
  inline short int getSensorID() const { return _iSensor; }

  //! Get original SimTrackerHit
  inline EVENT::SimTrackerHit* getSimTrackerHit() const { return _simHit; }

protected:
  CLHEP::Hep3Vector _prePosition; //!< PreStep position in cm.
  CLHEP::Hep3Vector _posPosition; //!< PosStep position in cm.
  CLHEP::Hep3Vector _momentum;    //!< Momentum in GeV

  short int _iLayer;  //!< ID number of a layer
  short int _iLadder; //!< ID number of a ladder
  short int _iSensor; //!< ID number of a sensor

  EVENT::SimTrackerHit* _simHit; //!< Original SimTrackerHit

}; // Class

} // namespace sistrip

#endif // SIMTRACKERDIGIHIT_H
