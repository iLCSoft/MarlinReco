#include <iostream>
#include <vector>
#include <EVENT/ReconstructedParticle.h>
#ifndef MyReconstructedParticle_hh
#define MyReconstructedParticle_hh 1
namespace TTbarAnalysis
{
	class MyReconstructedParticle 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			MyReconstructedParticle (EVENT::ReconstructedParticle * particle);
			virtual ~MyReconstructedParticle () {};
		//
		//	Methods
		//
			EVENT::ReconstructedParticle * Get();
			void SetAngle(float value);
			void SetOffset(float value);
			void SetSecOffset(float value);
			void SetAccuracy(float value);
			void SetObservable(float value);
			void SetCostheta(float value);

			float GetAngle();
			float GetOffset();
			float GetSecOffset();
			float GetAccuracy();
			float GetObservable();
			float GetCostheta();
			//void SetAngle(float angle);
		private:
		//
		//	Data
		//
			EVENT::ReconstructedParticle * myParticle;
			float myAngle;
			float myOffset;
			float mySecOffset;
			float myAccuracy;
			float myObservable;
			float myCos;
		//
		//	Private methods
		//
	};
}
#endif
