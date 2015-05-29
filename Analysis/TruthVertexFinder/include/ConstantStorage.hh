#include <iostream>
#include <vector>
#ifndef _TTbarAnalysis_h
#define _TTbarAnalysis_h 1
namespace TTbarAnalysis 
{
	enum PDGTYPE
	{
		EXCEPTIONAL_PDGTYPE = 0,
		BOTTOM_MESONS = 500,
		CHARMED_MESONS = 400,
		STRANGE_MESONS = 300,
		BOTTOM_BARYONS = 5000,
		BOTTOM_HADRONS = 5500,
		CHARMED_BARYONS = 4000,
		CHARMED_HADRONS = 4400,
		TAU_LEPTON = 15,
		TRACKABLE_PARTICLES = 1,
		NONTRACKABLE_PARTICLES = 1000000,     //Should be a large number
		_max_MESONS = NONTRACKABLE_PARTICLES  //Because of this line
	};
	class ConstantStorage 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			ConstantStorage ( ){}
			virtual ~ConstantStorage (){}
		//
		//	Methods
		//
		static const std::vector<int> & CHARMED_MESONS_PGD();
		static const std::vector<int> & STRANGE_MESONS_PGD();
		static const std::vector<int> & BOTTOM_MESONS_PGD();
		static const std::vector<int> & BOTTOM_BARYONS_PGD();
		static const std::vector<int> & BOTTOM_HADRONS_PGD();
		static const std::vector<int> & CHARMED_BARYONS_PGD();
		static const std::vector<int> & CHARMED_HADRONS_PGD();
		static const std::vector<int> & TAU_LEPTON_PGD();
		static const std::vector<int> & TRACKABLE_PARTICLES_PGD();
		static const std::vector<int> & NONTRACKABLE_PARTICLES_PGD();
		
		static const std::vector<int> & GET_PDG(PDGTYPE type);
		//static constructor pattern
		static class _init
		{
			public:
			_init()
			{
				///////////////CHARMEDMESONS/////////////////////
				myCharmedMesonsPDGs.push_back(421);
				myCharmedMesonsPDGs.push_back(411);
				//*
				myCharmedMesonsPDGs.push_back(433);
				myCharmedMesonsPDGs.push_back(413);
				myCharmedMesonsPDGs.push_back(431);
				myCharmedMesonsPDGs.push_back(10411);
				myCharmedMesonsPDGs.push_back(10421);
				myCharmedMesonsPDGs.push_back(423);
				myCharmedMesonsPDGs.push_back(10413);
				myCharmedMesonsPDGs.push_back(10423);
				myCharmedMesonsPDGs.push_back(20413);
				myCharmedMesonsPDGs.push_back(20423);
				myCharmedMesonsPDGs.push_back(415);
				myCharmedMesonsPDGs.push_back(425);
				myCharmedMesonsPDGs.push_back(10431);
				myCharmedMesonsPDGs.push_back(10433);
				myCharmedMesonsPDGs.push_back(20433);
				myCharmedMesonsPDGs.push_back(435);
				//*/
				///////////////BOTTOMMESONS//////////////////////
				myBottomMesonsPDGs.push_back(511);
				myBottomMesonsPDGs.push_back(521);
				///*
				myBottomMesonsPDGs.push_back(513);
				myBottomMesonsPDGs.push_back(523);
				myBottomMesonsPDGs.push_back(515);
				myBottomMesonsPDGs.push_back(525);
				myBottomMesonsPDGs.push_back(531);
				myBottomMesonsPDGs.push_back(533);
				myBottomMesonsPDGs.push_back(535);
				myBottomMesonsPDGs.push_back(541);
				myBottomMesonsPDGs.push_back(543);
				myBottomMesonsPDGs.push_back(545);
				myBottomMesonsPDGs.push_back(10511);
				myBottomMesonsPDGs.push_back(10521);
				myBottomMesonsPDGs.push_back(10523);
				myBottomMesonsPDGs.push_back(10513);
				myBottomMesonsPDGs.push_back(20513);
				myBottomMesonsPDGs.push_back(20523);
				myBottomMesonsPDGs.push_back(10531);
				myBottomMesonsPDGs.push_back(10533);
				myBottomMesonsPDGs.push_back(20533);
				myBottomMesonsPDGs.push_back(10541);
				myBottomMesonsPDGs.push_back(10543);
				myBottomMesonsPDGs.push_back(20543);
				///*/
				///////////////STRANGEMESONS//////////////////
				myStrangeMesonsPDGs.push_back(130);
				myStrangeMesonsPDGs.push_back(310);
				myStrangeMesonsPDGs.push_back(311);
				myStrangeMesonsPDGs.push_back(321);
				myStrangeMesonsPDGs.push_back(10311);
				myStrangeMesonsPDGs.push_back(10321);
				myStrangeMesonsPDGs.push_back(100311);
				myStrangeMesonsPDGs.push_back(100321);
				myStrangeMesonsPDGs.push_back(200311);
				myStrangeMesonsPDGs.push_back(200321);
				myStrangeMesonsPDGs.push_back(9000311);
				myStrangeMesonsPDGs.push_back(9000321);
				myStrangeMesonsPDGs.push_back(313);
				myStrangeMesonsPDGs.push_back(323);
				myStrangeMesonsPDGs.push_back(10313);
				myStrangeMesonsPDGs.push_back(10323);
				myStrangeMesonsPDGs.push_back(20313);
				myStrangeMesonsPDGs.push_back(20323);
				myStrangeMesonsPDGs.push_back(100313);
				myStrangeMesonsPDGs.push_back(10323);
				myStrangeMesonsPDGs.push_back(9000313);
				myStrangeMesonsPDGs.push_back(9000323);
				myStrangeMesonsPDGs.push_back(30313);
				myStrangeMesonsPDGs.push_back(30323);
				myStrangeMesonsPDGs.push_back(315);
				myStrangeMesonsPDGs.push_back(325);
				myStrangeMesonsPDGs.push_back(9000315);
				myStrangeMesonsPDGs.push_back(9000325);
				myStrangeMesonsPDGs.push_back(10315);
				myStrangeMesonsPDGs.push_back(10325);
				myStrangeMesonsPDGs.push_back(20315);
				myStrangeMesonsPDGs.push_back(20325);
				myStrangeMesonsPDGs.push_back(100315);
				myStrangeMesonsPDGs.push_back(100325);
				myStrangeMesonsPDGs.push_back(9010315);
				myStrangeMesonsPDGs.push_back(9010325);
				myStrangeMesonsPDGs.push_back(317);
				myStrangeMesonsPDGs.push_back(327);
				myStrangeMesonsPDGs.push_back(9010317);
				myStrangeMesonsPDGs.push_back(9010327);
				myStrangeMesonsPDGs.push_back(319);
				myStrangeMesonsPDGs.push_back(329);
				myStrangeMesonsPDGs.push_back(9000319);
				myStrangeMesonsPDGs.push_back(9000329);
				///////////////CHARMEDBARYONS//////////////////////
				myCharmedBaryonsPDGs.push_back(4122);
				myCharmedBaryonsPDGs.push_back(4222);
				myCharmedBaryonsPDGs.push_back(4212);
				myCharmedBaryonsPDGs.push_back(4112);
				myCharmedBaryonsPDGs.push_back(4224);
				myCharmedBaryonsPDGs.push_back(4214);
				myCharmedBaryonsPDGs.push_back(4114);
				myCharmedBaryonsPDGs.push_back(4232);
				myCharmedBaryonsPDGs.push_back(4132);
				myCharmedBaryonsPDGs.push_back(4322);
				myCharmedBaryonsPDGs.push_back(4312);
				myCharmedBaryonsPDGs.push_back(4324);
				myCharmedBaryonsPDGs.push_back(4314);
				myCharmedBaryonsPDGs.push_back(4332);
				myCharmedBaryonsPDGs.push_back(4334);
				myCharmedBaryonsPDGs.push_back(4412);
				myCharmedBaryonsPDGs.push_back(4422);
				myCharmedBaryonsPDGs.push_back(4414);
				myCharmedBaryonsPDGs.push_back(4424);
				myCharmedBaryonsPDGs.push_back(4432);
				myCharmedBaryonsPDGs.push_back(4434);
				myCharmedBaryonsPDGs.push_back(4444);
				///////////////CHARMEDHADRONS//////////////////////
				myCharmedHadronsPDGs.reserve( myCharmedBaryonsPDGs.size() + myCharmedMesonsPDGs.size() );
				myCharmedHadronsPDGs.insert(myCharmedHadronsPDGs.end(), myCharmedMesonsPDGs.begin(), myCharmedMesonsPDGs.end());
				myCharmedHadronsPDGs.insert(myCharmedHadronsPDGs.end(), myCharmedBaryonsPDGs.begin(), myCharmedBaryonsPDGs.end());
				///////////////BOTTOMBARYONS//////////////////////
				myBottomBaryonsPDGs.push_back(5122);
				myBottomBaryonsPDGs.push_back(5112);
				myBottomBaryonsPDGs.push_back(5212);
				myBottomBaryonsPDGs.push_back(5222);
				myBottomBaryonsPDGs.push_back(5114);
				myBottomBaryonsPDGs.push_back(5214);
				myBottomBaryonsPDGs.push_back(5224);
				myBottomBaryonsPDGs.push_back(5132);
				myBottomBaryonsPDGs.push_back(5232);
				myBottomBaryonsPDGs.push_back(5312);
				myBottomBaryonsPDGs.push_back(5322);
				myBottomBaryonsPDGs.push_back(5314);
				myBottomBaryonsPDGs.push_back(5324);
				myBottomBaryonsPDGs.push_back(5332);
				myBottomBaryonsPDGs.push_back(5334);
				myBottomBaryonsPDGs.push_back(5142);
				myBottomBaryonsPDGs.push_back(5242);
				myBottomBaryonsPDGs.push_back(5412);
				myBottomBaryonsPDGs.push_back(5422);
				myBottomBaryonsPDGs.push_back(5414);
				myBottomBaryonsPDGs.push_back(5424);
				myBottomBaryonsPDGs.push_back(5342);
				myBottomBaryonsPDGs.push_back(5432);
				myBottomBaryonsPDGs.push_back(5434);
				myBottomBaryonsPDGs.push_back(5442);
				myBottomBaryonsPDGs.push_back(5444);
				myBottomBaryonsPDGs.push_back(5512);
				myBottomBaryonsPDGs.push_back(5522);
				myBottomBaryonsPDGs.push_back(5514);
				myBottomBaryonsPDGs.push_back(5524);
				myBottomBaryonsPDGs.push_back(5532);
				myBottomBaryonsPDGs.push_back(5534);
				myBottomBaryonsPDGs.push_back(5542);
				myBottomBaryonsPDGs.push_back(5544);
				///////////////BOTTOMHADRONS//////////////////////
				myBottomHadronsPDGs.reserve( myBottomBaryonsPDGs.size() + myBottomMesonsPDGs.size() );
				myBottomHadronsPDGs.insert(myBottomHadronsPDGs.end(), myBottomMesonsPDGs.begin(), myBottomMesonsPDGs.end());
				myBottomHadronsPDGs.insert(myBottomHadronsPDGs.end(), myBottomBaryonsPDGs.begin(), myBottomBaryonsPDGs.end());
				///////////////NONTRACKABLE//////////////////////
				myTauLeptonPDGs.push_back(15);
				///////////////TRACKABLE//////////////////////
				myTrackableParticlesPDGs.push_back(11);
				myTrackableParticlesPDGs.push_back(13);
				myTrackableParticlesPDGs.push_back(211);
				myTrackableParticlesPDGs.push_back(321);
				myTrackableParticlesPDGs.push_back(2212);
				///////////////NONTRACKABLE//////////////////////
				myNonTrackableParticlesPDGs.push_back(111);
				myNonTrackableParticlesPDGs.push_back(22);
				myNonTrackableParticlesPDGs.push_back(2112);
				myNonTrackableParticlesPDGs.push_back(3122);
				myNonTrackableParticlesPDGs.push_back(130);
				myNonTrackableParticlesPDGs.push_back(310);
			}
		} _initializer;
		
		private:
		//
		//	Data
		//
			static std::vector<int> myCharmedMesonsPDGs; 
			static std::vector<int> myBottomMesonsPDGs; 
			static std::vector<int> myStrangeMesonsPDGs; 
			static std::vector<int> myBottomBaryonsPDGs; 
			static std::vector<int> myBottomHadronsPDGs; 
			static std::vector<int> myCharmedBaryonsPDGs; 
			static std::vector<int> myCharmedHadronsPDGs; 
			static std::vector<int> myTrackableParticlesPDGs; 
			static std::vector<int> myNonTrackableParticlesPDGs; 
			static std::vector<int> myTauLeptonPDGs; 
			static std::vector<int> myEmptyPDGs; 
		//
		//	Private methods
		//


	};
} /* TTbarAnalysis */
#endif
