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

#include <string>

#include <marlin/Processor.h>
#include <lcio.h>

#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>

using namespace lcio ;
using namespace marlin ;

class IsolatedLeptonFinderProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new IsolatedLeptonFinderProcessor ; }

		IsolatedLeptonFinderProcessor() ;

		virtual void init() ;
		virtual void processRunHeader( LCRunHeader* run ) ;
		virtual void processEvent( LCEvent * evt ) ; 
		virtual void check( LCEvent * evt ) ; 
		virtual void end() ;

	protected:

		/** Returns true if pfo is an isolated lepton */
		bool IsIsolatedLepton( ReconstructedParticle* pfo ) ;

		/** Returns true if isolated, as defined by the cone energy */
		bool IsIsolated( ReconstructedParticle* pfo ) ;

		/** Returns true if charged */
		bool IsCharged( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes lepton ID cuts */
		bool IsLepton( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes impact parameter cuts */
		bool PassesImpactParameterCuts( ReconstructedParticle* pfo ) ; 

		/** Calculates the cone energy */
		float getConeEnergy( ReconstructedParticle* pfo ) ;

		/** [0]:Ecal energy, [1]:Hcal energy */
		void getCalEnergy( ReconstructedParticle* pfo , float* cale) ;

		/** Input collection */
		std::string _inputPFOsCollection;

		/** Output collection (all input with isolated leptons removed) */
		std::string _outputPFOsRemovedIsoLepCollection;

		/** Output collection of isolated leptons */
		std::string _outputIsoLepCollection;

		LCCollection* _pfoCol;
		float _cosConeAngle;

		/** If set to true, uses PID cuts */
		bool _usePID;
		float _electronMinEnergyDepositByMomentum;
		float _electronMaxEnergyDepositByMomentum;
		float _electronMinEcalToHcalFraction;
		float _electronMaxEcalToHcalFraction;
		float _muonMinEnergyDepositByMomentum;
		float _muonMaxEnergyDepositByMomentum;
		float _muonMinEcalToHcalFraction;
		float _muonMaxEcalToHcalFraction;

		/** If set to true, uses impact parameter cuts */
		bool _useImpactParameter;
		float _minD0;
		float _maxD0;
		float _minZ0;
		float _maxZ0;
		float _minR0;
		float _maxR0;

		float _isoMinTrackEnergy;
		float _isoMaxTrackEnergy;
		float _isoMinConeEnergy;
		float _isoMaxConeEnergy;
		/** If set to true, uses polynomial cuts for isolation */
		bool _usePolynomialIsolation;
		float _isoPolynomialA;
		float _isoPolynomialB;
		float _isoPolynomialC;
} ;

#endif

