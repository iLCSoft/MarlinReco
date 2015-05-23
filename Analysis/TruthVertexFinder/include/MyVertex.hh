#include <EVENT/Vertex.h>
#include <IMPL/VertexImpl.h>
#include <EVENT/MCParticle.h>
#include <vector>

#ifndef MyVertex_h
#define MyVertex_h 1
namespace TTbarAnalysis 
{
	class MyVertex : public IMPL::VertexImpl
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			MyVertex ();
			virtual ~MyVertex () {};
		//
		//	Methods
		//
			void __SetMCParticles(const std::vector< EVENT::MCParticle * > particles);
			std::vector< EVENT::MCParticle * > & __GetMCParticles();
		private:
		//
		//	Data
		//
			std::vector< EVENT::MCParticle * > myMCParticles;
		//
		//	Private methods
		//
	};
}
#endif
