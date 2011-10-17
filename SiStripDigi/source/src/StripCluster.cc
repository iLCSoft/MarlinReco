#include "StripCluster.h"
#include "Colours.h"
#include <stdlib.h>

// Include Marlin
#include <streamlog/streamlog.h>
#include <cstdlib>

// Namespaces
using namespace CLHEP;
using namespace lcio;

namespace sistrip {

//
// Destructor
//
StripCluster::~StripCluster()
{
//	std::cout << "Deleting StripCluster" << std::endl;
}

//
// Set cluster position Three vector
//
void StripCluster::set3Position( const Hep3Vector & position)
{
   _position.setX( position.getX());
   _position.setY( position.getY());
   _position.setZ( position.getZ());
}

//
// Set cluster - position sigma Three vector
//
void StripCluster::set3PosSigma( const Hep3Vector & posSigma)
{
   _posSigma.setX( posSigma.getX());
   _posSigma.setY( posSigma.getY());
   _posSigma.setZ( posSigma.getZ());
}

//
// Update MC truth information - SimTrackerHit
//
void StripCluster::updateSimHitMap(SimTrackerHitMap simHitMap)
{
   for (SimTrackerHitMap::const_iterator iterSHM=simHitMap.begin(); iterSHM!=simHitMap.end(); iterSHM++) {

      EVENT::SimTrackerHit * simHit = iterSHM->first;
      float                  weight = iterSHM->second;

      if (_simHitMap.find(simHit) != _simHitMap.end()) _simHitMap[simHit] += weight;
      else                                             _simHitMap[simHit]  = weight;
   }
}

//
// Get MC truth information - sum of weights
//
float StripCluster::getSimHitWeightSum()
{
   float weightSum = 0;

   for (SimTrackerHitMap::const_iterator iterSHM=_simHitMap.begin(); iterSHM!=_simHitMap.end(); iterSHM++) {

      weightSum += iterSHM->second;
   }

   return weightSum;
}

} // Namespace;

