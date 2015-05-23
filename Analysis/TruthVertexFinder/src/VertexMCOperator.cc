#include "VertexMCOperator.hh"
using EVENT::Vertex;
using std::vector;
using std::string;
using EVENT::MCParticle;
using EVENT::ReconstructedParticle;
using IMPL::VertexImpl;
using IMPL::ReconstructedParticleImpl;
using EVENT::LCCollection;
using UTIL::LCRelationNavigator;
namespace TTbarAnalysis 
{
	VertexMCOperator:: VertexMCOperator(LCCollection * rel)
	{
		myRelCollection = rel;
	}
	vector< Vertex * > * VertexMCOperator::Construct(DecayChain * chain)
	{
		if (!chain || chain->GetSize() == 0) 
		{
			return NULL;
		}
		
		vector< Vertex * > * result = new vector< Vertex * >();
		const double * ip = chain->Get(0)->getVertex();
		for (unsigned int i = 1; i < chain->GetSize(); i++) // <<==============================
		{
			result->push_back(construct(chain->Get(i), ip, chain->GetParentPDG(), i+1));
		}
		for (unsigned int i = 0; i < result->size(); i++) 
		{
			addParticle(result->at(i), chain->Get(i));
		}
		return result;
	}

	Vertex * VertexMCOperator::construct(EVENT::MCParticle * particle, const double * ip, int pdg, int number)
	{
		//VertexImpl * result = new VertexImpl();
		MyVertex * result = new MyVertex();
		
		const double * initial;
			initial = particle->getVertex();
		float distance = MathOperator::getDistance(ip, initial);
		result->setPrimary(false);
		result->setAlgorithmType("VertexMCOperator");
		
		result->setPosition(initial[0], initial[1], initial[2]);
		result->addParameter (distance);
		result->addParameter (pdg);
		result->addParameter (number);

		return result;
	}
	void VertexMCOperator::AddProngs(Vertex * vertex, vector< MCParticle * > & particles, bool usingRelation)
	{
		if (!vertex || particles.size() == 0) 
		{
			streamlog_out(DEBUG) << "ERRORMC: argument is null!\n";
			return;
		}
		ReconstructedParticle * reco = vertex->getAssociatedParticle();
		for (unsigned int i = 0; i < particles.size(); i++) 
		{
			ReconstructedParticle * prong =(usingRelation)? getReco(particles[i]) : translate(particles[i]);
			if (prong) 
			{
				reco->addParticle(prong);
			}
			else 
			{
				streamlog_out(DEBUG) << "ERROR: Corrupted vertex!\n";

			}
		}
		MyVertex * myvertex = static_cast< MyVertex * >(vertex);
		myvertex->__SetMCParticles(particles);
		//streamlog_out(DEBUG) << "Added " << myvertex->__GetMCParticles().size() << " particles!\n";

	}
	void VertexMCOperator::addParticle(Vertex * vertex, MCParticle * particle)
	{
		if (!vertex || !particle) 
		{
			streamlog_out(DEBUG) << "ERRORMC: argument is null!\n";
			return;
		}
		ReconstructedParticle * reco = translate(particle);
		VertexImpl * ivertex = static_cast<VertexImpl*>(vertex);
		ivertex->setAssociatedParticle(reco);
	}
	ReconstructedParticle * VertexMCOperator::getReco(EVENT::MCParticle * particle)
	{
		LCRelationNavigator navigator(myRelCollection);
		int nvtx = navigator.getRelatedFromObjects(particle).size();
		streamlog_out(DEBUG) << "Particles: " << nvtx <<'\n';
		ReconstructedParticle * reco = NULL; // check!!!
		int winner = -1;
		float weight = 0.0;
		const vector< float > weights = navigator.getRelatedFromWeights(particle);
		if (nvtx > 0) 
		{
			for (int i = 0; i < nvtx; i++) 
			{
				if (weights[i] > weight) 
				{
					winner = i;
					weight = weights[i];
				}
			}
			reco = dynamic_cast<ReconstructedParticle*>(navigator.getRelatedFromObjects(particle)[winner]);
		}

		return reco;
	}
	ReconstructedParticle * VertexMCOperator::translate(EVENT::MCParticle * particle)
	{
		ReconstructedParticleImpl * reco = new ReconstructedParticleImpl();
		reco->setType(particle->getPDG());
		reco->setMass(particle->getMass());
		reco->setCharge(particle->getCharge());
		reco->setMomentum(particle->getMomentum());
		reco->setEnergy(particle->getEnergy());
		return reco;
	}
} /* TTbarAnalysis */
