#include "DigiCluster.h"

// Namespaces
using namespace CLHEP;
using namespace lcio;

namespace sistrip {

//
// Constructor
//
DigiCluster::DigiCluster() {
  _charge = 0;
  _time = 0.;
  _driftTime = 0.;
  _nCarriers = 0;

  _iLayer = 0;
  _iLadder = 0;
  _iSensor = 0;
  _iCell = 0;

  _position.setX(0.);
  _position.setY(0.);
  _position.setZ(0.);

  _velocity.setX(0.);
  _velocity.setY(0.);
  _velocity.setZ(0.);
}

//
// Destructor
//
DigiCluster::~DigiCluster() {
  //	std::cout << "Deleting DigiCluster" << std::endl;
}

//
// Set cluster position Three vector
//
void DigiCluster::set3Position(const Hep3Vector& position) {
  _position.setX(position.getX());
  _position.setY(position.getY());
  _position.setZ(position.getZ());
}

//
// Set cluster velocity Three vector
//
void DigiCluster::set3Velocity(const Hep3Vector& velocity) {
  _velocity.setX(velocity.getX());
  _velocity.setY(velocity.getY());
  _velocity.setZ(velocity.getZ());
}

} // namespace sistrip
