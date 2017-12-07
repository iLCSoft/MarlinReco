/**
 * @brief Marlin processor for finding isolated leptons.
 * @author Ryo Yonamine <yonamine@post.kek.jp>
 * @author Tomohiko Tanabe <tomohiko@icepp.s.u-tokyo.ac.jp>
 * @version $Id:$
 *
 * Given a list of ReconstructedParticle, identify isolated leptons
 * based on the track cone energy, lepton identification,
 * and the track impact parameters (optional).
 */
#ifndef IsolatedLeptonFinderProcessor_h
#define IsolatedLeptonFinderProcessor_h 1

#include <algorithm>
#include <string>
#include <map>

#include <marlin/Processor.h>
#include <lcio.h>

#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include "IMPL/ReconstructedParticleImpl.h"
#include <UTIL/LCRelationNavigator.h>

using namespace lcio ;
using namespace marlin ;

class IsolatedLeptonFinderProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new IsolatedLeptonFinderProcessor ; }

		IsolatedLeptonFinderProcessor() ;
		IsolatedLeptonFinderProcessor(const IsolatedLeptonFinderProcessor &);
		IsolatedLeptonFinderProcessor & operator = (const IsolatedLeptonFinderProcessor &);

		virtual void init() ;
		virtual void processEvent( LCEvent * evt ) ;
		virtual void end() ;

	protected:

		/** Returns true if pfo is a lepton */
		bool IsGoodLepton( ReconstructedParticle* pfo ) ;

		/** Returns true if pfo is an isolated lepton */
		bool IsIsolatedLepton( ReconstructedParticle* pfo , bool omitDressed) ;

		/** Returns true if isolated, as defined by the cone energy */
		bool IsIsolatedRectangular( ReconstructedParticle* pfo , bool omitDressed) ;
		bool IsIsolatedPolynomial( ReconstructedParticle* pfo , bool omitDressed) ;
		bool IsIsolatedJet( ReconstructedParticle* pfo ) ;

		/** Returns true if charged */
		bool IsCharged( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes muon or electron ID cuts */
		bool IsLepton( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes electron ID cuts */
		bool IsElectron( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes muon ID cuts */
		bool IsMuon( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes photon ID cuts */
		bool IsPhoton( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes impact parameter cuts */
		bool PassesImpactParameterCuts( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes impact parameter significance cuts */
		bool PassesImpactParameterSignificanceCuts( ReconstructedParticle* pfo ) ;

		/** Helper function to order PFOS by energy */
		bool isMoreEnergetic (int i, int j) {
			ReconstructedParticle* pfo_i = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(i) );
			ReconstructedParticle* pfo_j = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(j) );
			return (pfo_i->getEnergy()>pfo_j->getEnergy());
		}

		/** Adds photons around lepton to four vector */
		void dressLepton( ReconstructedParticleImpl* pfo, int PFO_idx ) ;

		/** Calculates the cone energy */
		float getConeEnergy( ReconstructedParticle* pfo, bool omitDressed) ;

		/** [0]:Ecal energy, [1]:Hcal energy */
		void getCalEnergy( ReconstructedParticle* pfo , float* cale) ;

		/** Replace missing copy constructor by hand */
		ReconstructedParticleImpl* CopyReconstructedParticle ( ReconstructedParticle* pfo ) ;

		/** Input collection */
		std::string _inputPFOsCollection{};

		/** Output collection (all input with isolated leptons removed) */
		std::string _outputPFOsRemovedIsoLepCollection{};

		/** Output collection of isolated leptons */
		std::string _outputIsoLepCollection{};

		/** Output collection (all input with dressed isolated leptons removed) */
		std::string _outputPFOsRemovedDressedIsoLepCollection{};

		/** Output collection of dressed isolated leptons */
		std::string _outputDressedIsoLepCollection{};

		LCCollection* _pfoCol=nullptr;
		float _cosConeAngle = 0;

		/** If set to true, uses PID cuts */
		bool _usePID = false;
		float _electronMinEnergyDepositByMomentum = 0;
		float _electronMaxEnergyDepositByMomentum = 0;
		float _electronMinEcalToHcalFraction = 0;
		float _electronMaxEcalToHcalFraction = 0;
		float _muonMinEnergyDepositByMomentum = 0;
		float _muonMaxEnergyDepositByMomentum = 0;
		float _muonMinEcalToHcalFraction = 0;
		float _muonMaxEcalToHcalFraction = 0;

		/** If set to true, uses impact parameter cuts */
		bool _useImpactParameter = false;
		float _minD0 = 0;
		float _maxD0 = 0;
		float _minZ0 = 0;
		float _maxZ0 = 0;
		float _minR0 = 0;
		float _maxR0 = 0;

		/** If set to true, uses impact parameter significance cuts */
		bool _useImpactParameterSignificance = false;
		float _minD0Sig = 0;
		float _maxD0Sig = 0;
		float _minZ0Sig = 0;
		float _maxZ0Sig = 0;
		float _minR0Sig = 0;
		float _maxR0Sig = 0;

		/** If set to true, uses rectangular cuts for isolation */
		bool _useRectangularIsolation = false;
		float _isoMinTrackEnergy = 0;
		float _isoMaxTrackEnergy = 0;
		float _isoMinConeEnergy = 0;
		float _isoMaxConeEnergy = 0;

		/** If set to true, uses polynomial cuts for isolation */
		bool _usePolynomialIsolation = false;
		float _isoPolynomialA = 0;
		float _isoPolynomialB = 0;
		float _isoPolynomialC = 0;

		/** If set to true, uses jet-based isolation (LAL algorithm) */
		bool _useJetIsolation = false;
		std::string _jetCollectionName {};
		std::map<ReconstructedParticle*,ReconstructedParticle*> _rpJetMap {};
		float _jetIsoVetoMinXt = 0;
		float _jetIsoVetoMaxXt = 0;
		float _jetIsoVetoMinZ = 0;
		float _jetIsoVetoMaxZ = 0;


		/** If set to true, uses lepton dressing */
		std::string _whichLeptons = "UNDRESSED";
		bool _mergeCloseElectrons = false;
		float _dressPhotonConeAngle = 0;
		float _mergeLeptonConeAngle = 0;
		std::vector<int> _dressedPFOs {};

		/** If set to true, uses Pandora particle IDs */
		bool _usePandoraIDs = false;
} ;

#endif

