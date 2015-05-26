#include "ConstantStorage.hh"
using std::vector;
namespace TTbarAnalysis 
{
	vector<int> ConstantStorage::myCharmedMesonsPDGs;
	vector<int> ConstantStorage::myBottomMesonsPDGs;
	vector<int> ConstantStorage::myStrangeMesonsPDGs;
	vector<int> ConstantStorage::myBottomBaryonsPDGs;
	vector<int> ConstantStorage::myBottomHadronsPDGs;
	vector<int> ConstantStorage::myCharmedBaryonsPDGs;
	vector<int> ConstantStorage::myCharmedHadronsPDGs;
	vector<int> ConstantStorage::myTrackableParticlesPDGs;
	vector<int> ConstantStorage::myNonTrackableParticlesPDGs;

	ConstantStorage::_init ConstantStorage::_initializer;

	const std::vector<int> & ConstantStorage::CHARMED_MESONS_PGD()
	{
		return myCharmedMesonsPDGs;
	}
	const std::vector<int> & ConstantStorage::STRANGE_MESONS_PGD()
	{
		return myStrangeMesonsPDGs;
	}
	const std::vector<int> & ConstantStorage::BOTTOM_MESONS_PGD()
	{
		return myBottomMesonsPDGs;
	}
	const std::vector<int> & ConstantStorage::BOTTOM_BARYONS_PGD()
	{
		return myBottomBaryonsPDGs;
	}
	const std::vector<int> & ConstantStorage::BOTTOM_HADRONS_PGD()
	{
		return myBottomHadronsPDGs;
	}
	const std::vector<int> & ConstantStorage::CHARMED_BARYONS_PGD()
	{
		return myCharmedBaryonsPDGs;
	}
	const std::vector<int> & ConstantStorage::CHARMED_HADRONS_PGD()
	{
		return myCharmedHadronsPDGs;
	}
	const std::vector<int> & ConstantStorage::TRACKABLE_PARTICLES_PGD()
	{
		return myTrackableParticlesPDGs;
	}
	const std::vector<int> & ConstantStorage::NONTRACKABLE_PARTICLES_PGD()
	{
		return myNonTrackableParticlesPDGs;
	}
	const std::vector<int> & ConstantStorage::GET_PDG(PDGTYPE type)
	{
		switch(type)
		{
			case BOTTOM_MESONS:
			return BOTTOM_MESONS_PGD();
			case CHARMED_MESONS:
			return CHARMED_MESONS_PGD();
			case STRANGE_MESONS:
			return STRANGE_MESONS_PGD();
			case BOTTOM_BARYONS:
			return BOTTOM_BARYONS_PGD();
			case BOTTOM_HADRONS:
			return BOTTOM_HADRONS_PGD();
			case CHARMED_HADRONS:
			return CHARMED_HADRONS_PGD();
			case CHARMED_BARYONS:
			return CHARMED_BARYONS_PGD();
			case TRACKABLE_PARTICLES:
			return TRACKABLE_PARTICLES_PGD();
			case NONTRACKABLE_PARTICLES:
			return NONTRACKABLE_PARTICLES_PGD();
			default:
			return *(new vector<int>());
		}
	}

} /* TTbarAnalysis */
