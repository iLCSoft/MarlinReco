#include "SimDigitalToEnergy.h"

#include <cassert>

#include <EVENT/CalorimeterHit.h>

using namespace lcio ;
using namespace marlin ;

SimDigitalToEnergy aSimDigitalToEnergy ;

SimDigitalToEnergy::SimDigitalToEnergy()
	: RealisticCaloReco::Processor("SimDigitalToEnergy")
{
	_description = "This processor transforms the threshold value of digitized SDHCAL CalorimeterHits into energy" ;

	std::vector<float> energyCoefficients = {0.4f} ;
	registerProcessorParameter("EnergyCalibration" ,
							   "Threshold to Energy correspondace" ,
							   _energyCoefficients,
							   energyCoefficients) ;
}

void SimDigitalToEnergy::init()
{
	//to avoid crash
	_calibrCoeff.push_back(0) ;
	_calLayers.push_back(0) ;

	RealisticCaloReco::init() ;
	assert( _energyCoefficients.size() > 0 ) ;
}

float SimDigitalToEnergy::reconstructEnergy(const CalorimeterHit* hit)
{
	float threshold = hit->getEnergy() - 1.0f ;

	if ( threshold < 0 )
		return 0.f ;

	unsigned int iCoeff = static_cast<unsigned int>( threshold ) ;
	if ( iCoeff >= _energyCoefficients.size() )
		return _energyCoefficients.at( _energyCoefficients.size()-1 ) ;
	else
		return _energyCoefficients.at( iCoeff ) ;
}

