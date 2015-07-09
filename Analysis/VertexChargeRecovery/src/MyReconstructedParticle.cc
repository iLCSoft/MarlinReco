#include "MyReconstructedParticle.hh"
using std::vector;
using std::string;
using EVENT::ReconstructedParticle;
namespace TTbarAnalysis
{
	MyReconstructedParticle:: MyReconstructedParticle(ReconstructedParticle * particle)
	{
		myParticle = particle;
		float value = 0.0;
		myAngle = value;
		myOffset = value;
		mySecOffset = value;
		myAccuracy = value;
		myObservable = value;
	}

	ReconstructedParticle * MyReconstructedParticle::Get()
	{
		return myParticle;
	}
	void MyReconstructedParticle::SetAngle(float value)
	{
		myAngle = value;
	}
	void MyReconstructedParticle::SetOffset(float value)
	{
		myOffset = value;
	}
	void MyReconstructedParticle::SetSecOffset(float value)
	{
		mySecOffset = value;
	}
	void MyReconstructedParticle::SetAccuracy(float value)
	{
		myAccuracy = value;
	}
	void MyReconstructedParticle::SetObservable(float value)
	{
		myObservable = value;
	}
	void MyReconstructedParticle::SetCostheta(float value)
	{
		myCos = value;
	}

	float MyReconstructedParticle::GetCostheta()
	{
		return myCos;
	}
	float MyReconstructedParticle::GetAngle()
	{
		return myAngle;
	}
	float MyReconstructedParticle::GetOffset()
	{
		return myOffset;
	}
	float MyReconstructedParticle::GetSecOffset()
	{
		return mySecOffset;
	}
	float MyReconstructedParticle::GetAccuracy()
	{
		return myAccuracy;
	}
	float MyReconstructedParticle::GetObservable()
	{
		return myObservable;
	}
	//void SetAngle(float angle);
	
}
