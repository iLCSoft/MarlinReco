#ifndef SimDigitalLinkToParticles_h
#define SimDigitalLinkToParticles_h

#include <marlin/Processor.h>
#include <IMPL/LCCollectionVec.h>

//used for standalone SDHCAL Simulation

class SimDigitalLinkToParticles : public marlin::Processor
{
	public :
		virtual marlin::Processor* newProcessor() { return new SimDigitalLinkToParticles ; }
		SimDigitalLinkToParticles() ;

		virtual void processEvent( LCEvent * evt ) ;

	protected :
		virtual void init() ;

		LCCollectionVec* processCollection(LCCollection* inputCol , LCCollection* inputRelCol) ;

		std::vector<std::string> _inputCollections{}; // input CalorimeterHit collection
		std::vector<std::string> _inputRelCollections{}; // input CalorimeterHit to SimCalorimeterHit relation collection
		std::vector<std::string> _outputRelCollections{}; // output CalorimeterHit to MCParticle relation collection
} ;

#endif //SimDigitalLinkToParticles_h
