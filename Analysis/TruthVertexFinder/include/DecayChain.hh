#include <stdlib.h>
#include <cmath>
#include <vector>
#include <string>
#include <EVENT/MCParticle.h>
#include "ConstantStorage.hh"
#ifndef _DecayChain_hh
#define _DecayChain_hh
namespace TTbarAnalysis 
{
	class DecayChain 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			DecayChain (const std::vector<EVENT::MCParticle *> * particles, std::string name, int pdg);
			DecayChain (std::string name, int pdg);
			virtual ~DecayChain () {};
		//
		//	Methods
		//
			void Add(EVENT::MCParticle * particle);
			EVENT::MCParticle * Get(int i) ;
			int GetSize() const;
			int GetParentPDG() const;
			std::string GetName() const;
			const std::vector< EVENT::MCParticle * > & GetAll() const;
			void Merge(DecayChain & other);
			EVENT::MCParticle * Find(PDGTYPE type) const;
			//oid Print();
		private:
		//
		//	Data
		//
			std::vector< EVENT::MCParticle * > myParticles;
			std::string myName;
			int myPDG;
		//
		//	Private methods
		//
	};
} /* DecayChain */
#endif

