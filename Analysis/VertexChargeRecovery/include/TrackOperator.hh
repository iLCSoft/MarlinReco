#include <iostream>
#include <vector>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include "MathOperator.hh"
#ifndef TrackOperator_hh
#define TrackOperator_hh 1
namespace TTbarAnalysis
{
	struct GConfig
	{
		const double * Momentum;
		const double * C;
		const float * IP;
		const double m;
		const double p;
		const double * v;
	};
	class TrackOperator 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			TrackOperator ();
			virtual ~TrackOperator () {};
			float GetOffset(EVENT::ReconstructedParticle * particle);
			void PrintTrack(EVENT::Track * track);
			EVENT::ReconstructedParticle * ReconstructParticle(EVENT::Track *  track);
			float GetDistanceBtw(const EVENT::Vertex * ip, const EVENT::Vertex * sec, const EVENT::ReconstructedParticle * particle);
			float GetError(EVENT::ReconstructedParticle * particle);
			float GetOffsetError(EVENT::ReconstructedParticle * particle, double * trackPosition, const EVENT::Vertex * ip, double offset);
			
			double * GetStartPoint(const EVENT::ReconstructedParticle * particle);
			float GetDprime(const EVENT::ReconstructedParticle * particle1, const EVENT::ReconstructedParticle * particle2, double * primaryPosition);
			void test();
			void test(EVENT::ReconstructedParticle * particle);
		//
		//	Methods
		//
		
		private:
		//
		//	Data
		//
			
		//
		//	Private methods
		//
			const std::vector<float> getErrorPoint(const EVENT::ReconstructedParticle * particle) const;
			float getError(GConfig & conf, const std::vector< float > pcovMatrix, const std::vector< float > ccovMatrix, const std::vector< float > ipcovMatrix);
			void printConf(GConfig & conf);
			double doffsetdC(GConfig & conf, int i);
			double doffsetdp(GConfig & conf, int i);

	};
}
#endif
