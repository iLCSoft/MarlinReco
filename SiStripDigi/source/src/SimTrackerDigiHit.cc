#include "SimTrackerDigiHit.h"

// Namespaces
using namespace CLHEP;
using namespace lcio;

namespace sistrip {

//
// Constructor
//
SimTrackerDigiHit::SimTrackerDigiHit()
{
   _cellID0    = 0.;
   _EDep      = 0.;
   _time       = 0.;
   _particle   = 0;
   _pathLength = 0.;

   _pos[0]   = 0.;
   _pos[1]   = 0.;
   _pos[2]   = 0.;

   _p[0]     = 0.;
   _p[1]     = 0.;
   _p[2]     = 0.;

   _prePosition.setX(0);
   _prePosition.setY(0);
   _prePosition.setZ(0);

   _posPosition.setX(0);
   _posPosition.setY(0);
   _posPosition.setZ(0);
}

//
// Destructor
//
SimTrackerDigiHit::~SimTrackerDigiHit()
{
//	std::cout << "Deleting SimTrackerDigiHit" << std::endl;
}


//
// Set preStep position of a hit - DEPRECATED method
//
void SimTrackerDigiHit::setPosition( double pos[3])
{
   setPrePosition(pos[0], pos[1], pos[2]);
}

//
// Set preStep position of a hit
//
void SimTrackerDigiHit::setPrePosition( double pos[3], float momentum[3], float pathLength)
{
   Hep3Vector direction(momentum[0], momentum[1], momentum[2]);
   direction = direction.unit();
   direction *= pathLength/2.;

   _prePosition.setX( pos[0] - direction.getX() );
   _prePosition.setY( pos[1] - direction.getY() );
   _prePosition.setZ( pos[2] - direction.getZ() );
}

//
// Set preStep position of a hit
//
void SimTrackerDigiHit::setPrePosition( double preX, double preY, double preZ)
{
   _pos[0] = preX;
   _pos[1] = preY;
   _pos[2] = preZ;

   _prePosition.setX(preX);
   _prePosition.setY(preY);
   _prePosition.setZ(preZ);

}

//
// Set preStep Three vector
//
void SimTrackerDigiHit::set3PrePosition( const Hep3Vector & prePosition)
{
   _pos[0] = prePosition.getX();
   _pos[1] = prePosition.getY();
   _pos[2] = prePosition.getZ();

   _prePosition = prePosition;
}

//
// Set posStep position of a hit
//
void SimTrackerDigiHit::setPosPosition( double pos[3], float momentum[3], float pathLength)
{
   Hep3Vector direction(momentum[0], momentum[1], momentum[2]);
   direction = direction.unit();
   direction *= pathLength/2.;

   _posPosition.setX( pos[0] + direction.getX() );
   _posPosition.setY( pos[1] + direction.getY() );
   _posPosition.setZ( pos[2] + direction.getZ() );
}

//
// Set posStep position of a hit
//
void SimTrackerDigiHit::setPosPosition( double posX, double posY, double posZ)
{
   _posPosition.setX(posX);
   _posPosition.setY(posY);
   _posPosition.setZ(posZ);
}

//
// Set posStep position of a hit
//
void SimTrackerDigiHit::set3PosPosition( const Hep3Vector & posPosition)
{
   _posPosition = posPosition;
}

//
// Set particle momentum at preStep position
//
void SimTrackerDigiHit::setMomentum( float p[3])
{
   _p[0] = p[0];
   _p[1] = p[1];
   _p[2] = p[2];

   _momentum.setX(p[0]);
   _momentum.setY(p[1]);
   _momentum.setZ(p[2]);
}

//
// Set particle momentum at preStep position
//
void SimTrackerDigiHit::setMomentum( float pX, float pY, float pZ)
{
   _p[0] = pX;
   _p[1] = pY;
   _p[2] = pZ;

   _momentum.setX(pX);
   _momentum.setY(pY);
   _momentum.setZ(pZ);
}

//
// Set particle Three momentum
//
void SimTrackerDigiHit::set3Momentum( const Hep3Vector & momentum)
{
   _p[0] = momentum.getX();
   _p[1] = momentum.getY();
   _p[2] = momentum.getZ();

   _momentum = momentum;
}

//
// Get step size
//
double SimTrackerDigiHit::getStepSize() const
{
   Hep3Vector step = _posPosition - _prePosition;

   return step.mag();
}

} // Namespace


