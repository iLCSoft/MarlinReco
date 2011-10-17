#include "Signal.h"
#include "Colours.h"
#include <cstdlib>
#include <iomanip>

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Namespaces
using namespace lcio;
using namespace marlin;

namespace sistrip {

//
// Update MC truth information - SimTrackerHit
//
void Signal::updateSimHitMap(EVENT::SimTrackerHit * simHit, float weight)
{
	if(_simHitMap.find(simHit) != _simHitMap.end())
	{
		_simHitMap[simHit] += weight;
	}
	else
	{
		_simHitMap[simHit]  = weight;
	}
}

void Signal::updateSimHitMap(SimTrackerHitMap simHitMap)
{
	for(SimTrackerHitMap::const_iterator iterSHM=simHitMap.begin(); iterSHM!=simHitMap.end(); ++iterSHM) 
	{
		EVENT::SimTrackerHit * simHit = iterSHM->first;
		float                  weight = iterSHM->second;
		
		if(_simHitMap.find(simHit) != _simHitMap.end()) _simHitMap[simHit] += weight;
		else                                             _simHitMap[simHit]  = weight;
   }
}

//
// Get MC truth information - sum of weights
//
float Signal::getSimHitWeightSum()
{
   float weightSum = 0;

   for (SimTrackerHitMap::const_iterator iterSHM=_simHitMap.begin(); iterSHM!=_simHitMap.end(); ++iterSHM) {

      weightSum += iterSHM->second;
   }

   return weightSum;
}

} // Namespace

