#ifndef SimDigitalToEnergy_h
#define SimDigitalToEnergy_h

#include <marlin/Processor.h>

#include "RealisticCaloReco.h"

/**
\addtogroup CaloDigi CaloDigi
@{
*/

class SimDigitalToEnergy : public RealisticCaloReco
{
	public :
		virtual Processor* newProcessor() { return new SimDigitalToEnergy ; }
		SimDigitalToEnergy() ;

	protected :
		virtual void init();
		virtual float reconstructEnergy(const CalorimeterHit* hit) ;

		std::vector<float> _energyCoefficients{};
} ;
/** @} */

#endif //SimDigitalToEnergy_h
