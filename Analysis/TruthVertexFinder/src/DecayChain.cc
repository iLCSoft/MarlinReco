#include "DecayChain.hh"
using std::vector;
using std::string;
using EVENT::MCParticle;
namespace TTbarAnalysis
{
	DecayChain ::DecayChain (const vector<MCParticle *> * particles, string name, int pdg)
	{
		myParticles = * particles;
		myPDG = pdg;
		myName = name;
	}
	DecayChain::DecayChain (string name, int pdg)
	{
		myName = name;
		myPDG = pdg;
	}
	
	void DecayChain::Add(MCParticle * particle)
	{
		myParticles.push_back(particle);
	}
	MCParticle * DecayChain::Get(int i) 
	{
		if (i > -1 && i < GetSize()) 
		{
			return myParticles[i];
		}
		return NULL; 
	}
	int DecayChain::GetSize() const
	{
		return myParticles.size();
	}
	const std::vector< MCParticle * > & DecayChain::GetAll() const
	{
		return myParticles;
	}
	void DecayChain::Merge(DecayChain & other)
	{
		for (int i = 0; i < other.GetSize(); i++) 
		{
			Add(other.Get(i));
		}
	}
	std::string DecayChain::GetName() const
	{
		return myName;
	}
	int DecayChain::GetParentPDG() const
	{
		return myPDG;
	}
	EVENT::MCParticle * DecayChain::Find(PDGTYPE type) const
	{
		vector<int> pdg =  ConstantStorage::GET_PDG(type);
		for (int i = 0; i < GetSize(); i++) 
		{
			for (int j = 0; j < pdg.size(); j++) 
			{
				if (abs(myParticles[i]->getPDG()) == pdg[i]) 
				{
					return myParticles[i];
				}
			}
		}
	}
	
	
}

