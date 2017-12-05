#include "IsolatedLeptonFinderProcessor.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>

// ----- include for verbosity dependend logging ---------
#include <marlin/VerbosityLevels.h>

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

using namespace lcio ;
using namespace marlin ;

IsolatedLeptonFinderProcessor aIsolatedLeptonFinderProcessor ;

IsolatedLeptonFinderProcessor::IsolatedLeptonFinderProcessor()
	: Processor("IsolatedLeptonFinderProcessor") {

		// Processor description
		_description = "Isolated Lepton Finder Processor" ;

		// register steering parameters: name, description, class-variable, default value
		registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"InputCollection" ,
				"Input collection of ReconstructedParticles",
				_inputPFOsCollection,
				std::string("PandoraPFOs"));

		registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"OutputCollectionWithoutIsolatedLepton",
				"Copy of input collection but without the isolated leptons",
				_outputPFOsRemovedIsoLepCollection,
				std::string("PandoraPFOsWithoutIsoLep") );

		registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"OutputCollectionIsolatedLeptons",
				"Output collection of isolated leptons",
				_outputIsoLepCollection,
				std::string("Isolep") );

		registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"OutputCollectionWithoutDressedIsolatedLepton",
				"Copy of input collection but without the dressed isolated leptons",
				_outputPFOsRemovedDressedIsoLepCollection,
				std::string("PandoraPFOsWithoutDressedIsoLep") );

		registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"OutputCollectionDressedIsolatedLeptons",
				"Output collection of dressed isolated leptons",
				_outputDressedIsoLepCollection,
				std::string("DressedIsolep") );

		registerProcessorParameter( "CosConeAngle",
				"Cosine of the half-angle of the cone used in isolation criteria",
				_cosConeAngle,
				float(0.98));

		registerProcessorParameter( "UsePID",
				"Use primitive particle ID based on calorimeter energy deposits",
				_usePID,
				bool(true));

		registerProcessorParameter( "ElectronMinEnergyDepositByMomentum",
				"Electron ID: Minimum energy deposit divided by momentum",
				_electronMinEnergyDepositByMomentum,
				float(0.7));

		registerProcessorParameter( "ElectronMaxEnergyDepositByMomentum",
				"Electron ID: Maximum energy deposit divided by momentum",
				_electronMaxEnergyDepositByMomentum,
				float(1.4));

		registerProcessorParameter( "ElectronMinEcalToHcalFraction",
				"Electron ID: Minimum Ecal deposit divided by sum of Ecal and Hcal deposits",
				_electronMinEcalToHcalFraction,
				float(0.9));

		registerProcessorParameter( "ElectronMaxEcalToHcalFraction",
				"Electron ID: Maximum Ecal deposit divided by sum of Ecal and Hcal deposits",
				_electronMaxEcalToHcalFraction,
				float(1.0));

		registerProcessorParameter( "MuonMinEnergyDepositByMomentum",
				"Muon ID: Minimum energy deposit divided by momentum",
				_muonMinEnergyDepositByMomentum,
				float(0.0));

		registerProcessorParameter( "MuonMaxEnergyDepositByMomentum",
				"Muon ID: Maximum energy deposit divided by momentum",
				_muonMaxEnergyDepositByMomentum,
				float(0.3));

		registerProcessorParameter( "MuonMinEcalToHcalFraction",
				"Muon ID: Minimum Ecal deposit divided by sum of Ecal and Hcal deposits",
				_muonMinEcalToHcalFraction,
				float(0.0));

		registerProcessorParameter( "MuonMaxEcalToHcalFraction",
				"Muon ID: Maximum Ecal deposit divided by sum of Ecal and Hcal deposits",
				_muonMaxEcalToHcalFraction,
				float(0.4));

		registerProcessorParameter( "UseImpactParameter",
				"Use impact parameter cuts for consistency with primary/secondary track",
				_useImpactParameter,
				bool(true));

		registerProcessorParameter( "ImpactParameterMinD0",
				"Minimum d0 impact parameter",
				_minD0,
				float(0.0));

		registerProcessorParameter( "ImpactParameterMaxD0",
				"Maximum d0 impact parameter",
				_maxD0,
				float(1e20));

		registerProcessorParameter( "ImpactParameterMinZ0",
				"Minimum z0 impact parameter",
				_minZ0,
				float(0.0));

		registerProcessorParameter( "ImpactParameterMaxZ0",
				"Maximum z0 impact parameter",
				_maxZ0,
				float(1e20));

		registerProcessorParameter( "ImpactParameterMin3D",
				"Minimum impact parameter in 3D",
				_minR0,
				float(0.0));

		registerProcessorParameter( "ImpactParameterMax3D",
				"Maximum impact parameter in 3D",
				_maxR0,
				float(0.01));

		registerProcessorParameter( "UseImpactParameterSignificance",
				"Use impact parameter significance cuts for consistency with primary/secondary track",
				_useImpactParameterSignificance,
				bool(false));

		registerProcessorParameter( "ImpactParameterMinD0Significance",
				"Minimum d0 impact parameter significance",
				_minD0Sig,
				float(0.0));

		registerProcessorParameter( "ImpactParameterMaxD0Significance",
				"Maximum d0 impact parameter significance",
				_maxD0Sig,
				float(1e20));

		registerProcessorParameter( "ImpactParameterMinZ0Significance",
				"Minimum z0 impact parameter significance",
				_minZ0Sig,
				float(0.0));

		registerProcessorParameter( "ImpactParameterMaxZ0Significance",
				"Maximum z0 impact parameter significance",
				_maxZ0Sig,
				float(1e20));

		registerProcessorParameter( "ImpactParameterMin3DSignificance",
				"Minimum impact parameter significance in 3D",
				_minR0Sig,
				float(0.0));

		registerProcessorParameter( "ImpactParameterMax3DSignificance",
				"Maximum impact parameter significance in 3D",
				_maxR0Sig,
				float(1e20));

		registerProcessorParameter( "UseRectangularIsolation",
				"Use rectangular cuts on track and cone energy",
				_useRectangularIsolation,
				bool(true));

		registerProcessorParameter( "IsolationMinimumTrackEnergy",
				"Minimum track energy for isolation requirement",
				_isoMinTrackEnergy,
				float(15));

		registerProcessorParameter( "IsolationMaximumTrackEnergy",
				"Maximum track energy for isolation requirement",
				_isoMaxTrackEnergy,
				float(1e20));

		registerProcessorParameter( "IsolationMinimumConeEnergy",
				"Minimum cone energy for isolation requirement",
				_isoMinConeEnergy,
				float(0));

		registerProcessorParameter( "IsolationMaximumConeEnergy",
				"Maximum cone energy for isolation requirement",
				_isoMaxConeEnergy,
				float(1e20));

		registerProcessorParameter( "UsePolynomialIsolation",
				"Use polynomial cuts on track and cone energy",
				_usePolynomialIsolation,
				bool(true));

		registerProcessorParameter( "IsolationPolynomialCutA",
				"Polynomial cut (A) on track energy and cone energy: Econe^2 < A*Etrk^2+B*Etrk+C",
				_isoPolynomialA,
				float(0));

		registerProcessorParameter( "IsolationPolynomialCutB",
				"Polynomial cut (B) on track energy and cone energy: Econe^2 < A*Etrk^2+B*Etrk+C",
				_isoPolynomialB,
				float(20));

		registerProcessorParameter( "IsolationPolynomialCutC",
				"Polynomial cut (C) on track energy and cone energy: Econe^2 < A*Etrk^2+B*Etrk+C",
				_isoPolynomialC,
				float(-300));

		registerProcessorParameter( "UseJetIsolation",
				"Use jet-based isolation",
				_useJetIsolation,
				bool(false));

		registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"JetCollection" ,
				"Input collection of jets for isolation",
				_jetCollectionName,
				std::string("JetsForIsolation"));

		registerProcessorParameter( "JetIsolationVetoMinimumXt",
				"Minimum Xt in jet-based isolation",
				_jetIsoVetoMinXt,
				float(0.));

		registerProcessorParameter( "JetIsolationVetoMaximumXt",
				"Maximum Xt in jet-based isolation",
				_jetIsoVetoMaxXt,
				float(0.25));

		registerProcessorParameter( "JetIsolationVetoMinimumZ",
				"Mininum Z in jet-based isolation",
				_jetIsoVetoMinZ,
				float(0.));

		registerProcessorParameter( "JetIsolationVetoMaximumZ",
				"Maximum Z in jet-based isolation",
				_jetIsoVetoMaxZ,
				float(0.6));

		registerProcessorParameter( "UseLeptonDressing",
				"Use lepton dressing",
				_useLeptonDressing,
				bool(true));

		registerProcessorParameter( "DressPhotonCosConeAngle",
				"Cosine of the half-angle of the cone used for lepton dressing with photons",
				_dressPhotonCosConeAngle,
				float(0.999848));

		registerProcessorParameter( "MergeLeptonCosConeAngle",
				"Cosine of the half-angle of the cone used for lepton merging",
				_mergeLeptonCosConeAngle,
				float(0.999048));
	}


void IsolatedLeptonFinderProcessor::init() {
	streamlog_out(DEBUG) << "   init called  " << std::endl ;
	printParameters() ;
}

void IsolatedLeptonFinderProcessor::processEvent( LCEvent * evt ) {

	streamlog_out(MESSAGE) <<std::endl;
	streamlog_out(MESSAGE) << "processing event: " << evt->getEventNumber() << "   in run:  " << evt->getRunNumber() << std::endl ;

	_dressedPFOs.clear();
	_rpJetMap.clear();

	_pfoCol = evt->getCollection( _inputPFOsCollection ) ;

	// Output PFOs removed isolated leptons
	LCCollectionVec* otPFOsRemovedIsoLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
	otPFOsRemovedIsoLepCol->setSubset(true) ;

	// Output PFOs of isolated leptons
	LCCollectionVec* otIsoLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
	otIsoLepCol->setSubset(true);

	// Output PFOs removed dressed isolated leptons
	LCCollectionVec* otPFOsRemovedDressedIsoLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
	otPFOsRemovedDressedIsoLepCol->setSubset(true) ;

	// Output PFOs of dressed isolated leptons
	LCCollectionVec* otDressedIsoLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
	otDressedIsoLepCol->setSubset(true);

	// Prepare jet/recoparticle map for jet-based isolation
	if (_useJetIsolation) {
		LCCollection *colJet = evt->getCollection(_jetCollectionName);
		int njet = colJet->getNumberOfElements();
		for (int i=0; i<njet; ++i) {
			ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( colJet->getElementAt(i) );
			for (ReconstructedParticleVec::const_iterator iter = jet->getParticles().begin();
					iter != jet->getParticles().end(); ++iter) {
				_rpJetMap.insert( std::make_pair( *iter, jet ) );
			}
		}
	}


	// Undressed leptons
	std::vector<int> goodLeptonIndices;
	int npfo = _pfoCol->getNumberOfElements();
	for (int i = 0; i < npfo; i++ ) {
		ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(i) );

		if ( !IsGoodLepton( pfo )){
			otPFOsRemovedIsoLepCol->addElement( pfo );
			continue;
		}

		if  ( IsIsolatedLepton( pfo, false ) ){
			streamlog_out(DEBUG) << "ISOLATION undressed "<<pfo->getType()<<": SUCCESS !" << std::endl;
			otIsoLepCol->addElement( pfo );
		}
		else{
			streamlog_out(DEBUG) << "ISOLATION undressed "<<pfo->getType()<<": FAILED" << std::endl;
			otPFOsRemovedIsoLepCol->addElement( pfo );
		}

		// remember lepton indices for dressing
		// goodLeptonIndices.push_back(i);
		if (fabs(pfo->getType()) == 11 || fabs(pfo->getType()) == 13) goodLeptonIndices.push_back(i);
	}


	// Dressed leptons
	// order by energy
	for (unsigned int i = 1; i < goodLeptonIndices.size(); ++i)
	{
		if (isMoreEnergetic(i, i-1)) {
			std::swap(goodLeptonIndices.at(i), goodLeptonIndices.at(i-1));
			i = 1;
		}
	}
	// dress them
	for (unsigned int i = 0; i < goodLeptonIndices.size(); ++i)
	{
		// copy this in an ugly fashion to be modofiable - a versatile copy constructor would be much better!
		ReconstructedParticle* pfo_tmp = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(goodLeptonIndices.at(i) ));
		ReconstructedParticleImpl* pfo = new ReconstructedParticleImpl();
		pfo->setMomentum(pfo_tmp->getMomentum());
		pfo->setEnergy(pfo_tmp->getEnergy());
		pfo->setType(pfo_tmp->getType());
		pfo->setCovMatrix(pfo_tmp->getCovMatrix());
		pfo->setMass(pfo_tmp->getMass());
		pfo->setCharge(pfo_tmp->getCharge());
		pfo->setParticleIDUsed(pfo_tmp->getParticleIDUsed());
		pfo->setGoodnessOfPID(pfo_tmp->getGoodnessOfPID());
		pfo->setStartVertex(pfo_tmp->getStartVertex());

		// test how close they are to the other leptons

		streamlog_out(DEBUG) << "Processing lep "<<i<<" with type "<< pfo->getType()<<" with E = "<<pfo->getEnergy()<<std::endl;
		// for (unsigned int j = i+1; j < goodLeptonIndices.size(); ++j)
		// {
		// 	ReconstructedParticle* pfo_other = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(goodLeptonIndices.at(j) ));
		// 	TVector3 P_this( pfo->getMomentum() );
		// 	TVector3 P_other( pfo_other->getMomentum() );
		// 	float cosTheta = P_this.Dot( P_other )/(P_this.Mag()*P_other.Mag());
		// 	streamlog_out(DEBUG) << "Lep "<<i<<"("<<goodLeptonIndices.at(i)<<") - "<<j<<"("<<goodLeptonIndices.at(j)<<"): "<<cosTheta<<std::endl;
		// }

		// don't reprocess merged leptons
		if (std::find(_dressedPFOs.begin(), _dressedPFOs.end(), i) != _dressedPFOs.end()){
			continue;
		}

		dressLepton(pfo, goodLeptonIndices.at(i));
  		// printf("dressedMomentum: %.2f -> %.2f, %.2f -> %.2f ,%.2f -> %.2f ,%.2f -> %.2f\n", pfo_tmp->getMomentum()[0], pfo->getMomentum()[0], pfo_tmp->getMomentum()[1], pfo->getMomentum()[1], pfo_tmp->getMomentum()[2], pfo->getMomentum()[2], pfo->getEnergy(), pfo_tmp->getEnergy());

		if  ( IsIsolatedLepton( pfo, true ) ){
			streamlog_out(DEBUG) << "ISOLATION dressed "<<pfo->getType()<<": SUCCESS !" << std::endl;
			otDressedIsoLepCol->addElement( pfo );
		}
		else{
			streamlog_out(DEBUG) << "ISOLATION dressed "<<pfo->getType()<<": FAILED" << std::endl;
			otPFOsRemovedDressedIsoLepCol->addElement( pfo );
		}
	}

	// PFO loop for filling remaining PFOs
	for (int i = 0; i < npfo; i++ ) {

		// don't add leptons again
		if (std::find(goodLeptonIndices.begin(), goodLeptonIndices.end(), i) != goodLeptonIndices.end()){
			continue;
		}

		// don't add dressed PFOs
		if (std::find(_dressedPFOs.begin(), _dressedPFOs.end(), i) != _dressedPFOs.end()){
			continue;
		}

		ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(i) );
		otPFOsRemovedDressedIsoLepCol->addElement( pfo );
	}


	streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber()
		<< "   in run:  " << evt->getRunNumber()
		<< std::endl ;


	// calculate total energy
	double etot_iso(0);
	for (int i = 0; i < otPFOsRemovedIsoLepCol->getNumberOfElements(); ++i)
	{
		ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( otPFOsRemovedIsoLepCol->getElementAt(i) );
		etot_iso += pfo->getEnergy();
	}
	for (int i = 0; i < otIsoLepCol->getNumberOfElements(); ++i)
	{
		ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( otIsoLepCol->getElementAt(i) );
		etot_iso += pfo->getEnergy();
	}
	double etot_dressediso(0);
	for (int i = 0; i < otPFOsRemovedDressedIsoLepCol->getNumberOfElements(); ++i)
	{
		ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( otPFOsRemovedDressedIsoLepCol->getElementAt(i) );
		etot_dressediso += pfo->getEnergy();
	}
	for (int i = 0; i < otDressedIsoLepCol->getNumberOfElements(); ++i)
	{
		ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( otDressedIsoLepCol->getElementAt(i) );
		etot_dressediso += pfo->getEnergy();
	}

	streamlog_out(DEBUG) << "npfo:                     " <<npfo<<std::endl;
	streamlog_out(DEBUG) << "nLepton:                  " <<goodLeptonIndices.size()<<std::endl;
	streamlog_out(DEBUG) << "nDressed:                 " <<_dressedPFOs.size()<<std::endl;
	streamlog_out(DEBUG) << "Energy:                   " <<etot_iso<<" and "<<etot_dressediso<<std::endl;
	streamlog_out(DEBUG) << "Sizes removed collection: " <<otPFOsRemovedIsoLepCol->getNumberOfElements()<<" and "<<otPFOsRemovedDressedIsoLepCol->getNumberOfElements()<<std::endl;
	streamlog_out(DEBUG) << "Sizes lepton collection:  " <<otIsoLepCol->getNumberOfElements()<<" and "<<otDressedIsoLepCol->getNumberOfElements()<<std::endl;


	// Add PFOs to new collections
	evt->addCollection( otPFOsRemovedIsoLepCol, _outputPFOsRemovedIsoLepCollection.c_str() );
	evt->addCollection( otIsoLepCol, _outputIsoLepCollection.c_str() );
	evt->addCollection( otPFOsRemovedDressedIsoLepCol, _outputPFOsRemovedDressedIsoLepCollection.c_str() );
	evt->addCollection( otDressedIsoLepCol, _outputDressedIsoLepCollection.c_str() );
}
void IsolatedLeptonFinderProcessor::dressLepton( ReconstructedParticleImpl* pfo, int PFO_idx ) {
	TVector3 P_lep( pfo->getMomentum() );
	int npfo = _pfoCol->getNumberOfElements();
	for ( int i = 0; i < npfo; i++ ) {
		ReconstructedParticle* pfo_dress = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(i) );

		// only add photons and electrons
		bool isPhoton = (pfo_dress->getType() == 22);
		bool isElectron = (fabs(pfo_dress->getType()) == 11);
		if (!isPhoton && !isElectron) continue;

		// don't add itself to itself
		if ( i == PFO_idx ) continue;

		TVector3 Pdress( pfo_dress->getMomentum() );
		float cosTheta = P_lep.Dot( Pdress )/(P_lep.Mag()*Pdress.Mag());

		if ( (isPhoton && cosTheta >= _dressPhotonCosConeAngle) ||
			 (isElectron && cosTheta >= _mergeLeptonCosConeAngle) ){
			if (std::find(_dressedPFOs.begin(), _dressedPFOs.end(), i) != _dressedPFOs.end()){
				if (pfo_dress->getType() == 22) streamlog_out(DEBUG) << "WARNING: photon "<<i<<" with cosTheta "<<cosTheta<<" already close to another lepton!"<<std::endl;
				if (fabs(pfo_dress->getType()) == 11) streamlog_out(DEBUG) << "WARNING: lepton "<<i<<" with cosTheta "<<cosTheta<<" already close to another lepton!"<<std::endl;
				// printf(" -- this lep: %.2f, %.2f ,%.2f ,%.2f\n", pfo->getMomentum()[0], pfo->getMomentum()[1], pfo->getMomentum()[2], pfo->getEnergy());
				continue;
			}
			if (pfo_dress->getType() == 22) streamlog_out(DEBUG) << "MESSAGE: dressing photon "<<i<<" with cosTheta "<<cosTheta<<std::endl;
			if (fabs(pfo_dress->getType()) == 11) streamlog_out(DEBUG) << "MESSAGE: dressing lepton "<<i<<" with cosTheta "<<cosTheta<<std::endl;
			_dressedPFOs.push_back(i);
			double dressedMomentum[3] = {pfo->getMomentum()[0] + pfo_dress->getMomentum()[0],
								  		 pfo->getMomentum()[1] + pfo_dress->getMomentum()[1],
								  		 pfo->getMomentum()[2] + pfo_dress->getMomentum()[2]};
			double dressedE = pfo->getEnergy() + pfo_dress->getEnergy();
			pfo->setMomentum(dressedMomentum);
			pfo->setEnergy(dressedE);
		}
	}
}
void IsolatedLeptonFinderProcessor::end() {
}

bool IsolatedLeptonFinderProcessor::IsCharged( ReconstructedParticle* pfo ) {
	if ( pfo->getCharge() == 0 ) return false;
	return true;
}

bool IsolatedLeptonFinderProcessor::IsLepton( ReconstructedParticle* pfo ) {

	float CalE[2];
	getCalEnergy( pfo , CalE );
	double ecale  = CalE[0];
	double hcale  = CalE[1];
	double p      = TVector3( pfo->getMomentum() ).Mag();
	double calByP = p>0 ? (ecale + hcale)/p : 0;
	double calSum = ecale+hcale;
	double ecalFrac = calSum>0 ? ecale / calSum : 0;

	// electron
	if ( calByP >= _electronMinEnergyDepositByMomentum
			&& calByP <= _electronMaxEnergyDepositByMomentum
			&& ecalFrac >= _electronMinEcalToHcalFraction
			&& ecalFrac <= _electronMaxEcalToHcalFraction )
		return true;

	// muon
	if ( calByP >= _muonMinEnergyDepositByMomentum
			&& calByP <= _muonMaxEnergyDepositByMomentum
			&& ecalFrac >= _muonMinEcalToHcalFraction
			&& ecalFrac <= _muonMaxEcalToHcalFraction )
		return true;

	return false;
}

bool IsolatedLeptonFinderProcessor::IsGoodLepton( ReconstructedParticle* pfo ) {

	if ( !IsCharged(pfo) )
		return false;

	if ( _usePID && !IsLepton(pfo) )
		return false;

	if ( _useImpactParameter && !PassesImpactParameterCuts(pfo) )
		return false ;

	if ( _useImpactParameterSignificance && !PassesImpactParameterSignificanceCuts(pfo) )
		return false ;

	return true;
}

bool IsolatedLeptonFinderProcessor::IsIsolatedLepton( ReconstructedParticle* pfo, bool omitDressed ) {

	if ( _useRectangularIsolation && !IsIsolatedRectangular(pfo, omitDressed) )
		return false;

	if ( _usePolynomialIsolation && !IsIsolatedPolynomial(pfo, omitDressed) )
		return false;

	if ( _useJetIsolation && !IsIsolatedJet(pfo) )
		return false;

	return true;
}

bool IsolatedLeptonFinderProcessor::IsIsolatedRectangular( ReconstructedParticle* pfo, bool omitDressed ) {
	float E     = pfo->getEnergy() ;
	float coneE = getConeEnergy( pfo, omitDressed );

	if (E < _isoMinTrackEnergy) return false;
	if (E > _isoMaxTrackEnergy) return false;
	if (coneE < _isoMinConeEnergy) return false;
	if (coneE > _isoMaxConeEnergy) return false;

	return true;
}

bool IsolatedLeptonFinderProcessor::IsIsolatedPolynomial( ReconstructedParticle* pfo, bool omitDressed ) {
	float E     = pfo->getEnergy() ;
	float coneE = getConeEnergy( pfo, omitDressed );

	if ( coneE*coneE <= _isoPolynomialA*E*E + _isoPolynomialB*E + _isoPolynomialC )
		return true ;
	return false;
}

bool IsolatedLeptonFinderProcessor::IsIsolatedJet( ReconstructedParticle* pfo ) {
	// jet-based isolated lepton (LAL algorithm)

	if ( _rpJetMap.find( pfo ) == _rpJetMap.end() ) {
		// this is often the case when jet finding fails e.g. due to too few particles in event
		return false;
	}

	ReconstructedParticle* jet = _rpJetMap[pfo];
	TVector3 vec1( pfo->getMomentum() );
	TVector3 jetmom( jet->getMomentum() );
	TLorentzVector jetmom4( jet->getMomentum(), jet->getEnergy() );

	float jetxt = vec1.Pt( jetmom )/jetmom4.M();
	float jetz = pfo->getEnergy()/jet->getEnergy();

	if (jetxt >= _jetIsoVetoMinXt && jetxt < _jetIsoVetoMaxXt
			&& jetz >= _jetIsoVetoMinZ && jetz < _jetIsoVetoMaxZ) {
		//printf("xt=%f z=%f (not pass)\n",jetxt,jetz);
		return false;
	}

	//printf("xt=%f z=%f (PASS)\n",jetxt,jetz);
	return true;
}

bool IsolatedLeptonFinderProcessor::PassesImpactParameterCuts( ReconstructedParticle* pfo ) {
	const EVENT::TrackVec & trkvec = pfo->getTracks();

	if (trkvec.size()==0) return false;

	// TODO: more sophisticated pfo/track matching
	float d0 = fabs(trkvec[0]->getD0());
	float z0 = fabs(trkvec[0]->getZ0());
	float r0 = sqrt( d0*d0 + z0*z0 );

	if ( d0 < _minD0 ) return false;
	if ( d0 > _maxD0 ) return false;
	if ( z0 < _minZ0 ) return false;
	if ( z0 > _maxZ0 ) return false;
	if ( r0 < _minR0 ) return false;
	if ( r0 > _maxR0 ) return false;

	return true;
}

bool IsolatedLeptonFinderProcessor::PassesImpactParameterSignificanceCuts( ReconstructedParticle* pfo ) {
	const EVENT::TrackVec & trkvec = pfo->getTracks();

	if (trkvec.size()==0) return false;

	// TODO: more sophisticated pfo/track matching
	float d0 = fabs(trkvec[0]->getD0());
	float z0 = fabs(trkvec[0]->getZ0());
	float d0err = sqrt(trkvec[0]->getCovMatrix()[0]);
	float z0err = sqrt(trkvec[0]->getCovMatrix()[9]);

	float d0sig = d0err != 0 ? d0/d0err : 0;
	float z0sig = z0err != 0 ? z0/z0err : 0;
	float r0sig = sqrt( d0sig*d0sig + z0sig*z0sig );

	if ( d0sig < _minD0Sig ) return false;
	if ( d0sig > _maxD0Sig ) return false;
	if ( z0sig < _minZ0Sig ) return false;
	if ( z0sig > _maxZ0Sig ) return false;
	if ( r0sig < _minR0Sig ) return false;
	if ( r0sig > _maxR0Sig ) return false;

	return true;
}

float IsolatedLeptonFinderProcessor::getConeEnergy( ReconstructedParticle* pfo, bool omitDressed ) {
	float coneE = 0;

	TVector3 P( pfo->getMomentum() );
	int npfo = _pfoCol->getNumberOfElements();
	for ( int i = 0; i < npfo; i++ ) {
		ReconstructedParticle* pfo_i = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(i) );

		// don't add itself to the cone energy
		if ( pfo == pfo_i ) continue;

		// don't add dressed PFOs to the cone energy
		if (omitDressed && std::find(_dressedPFOs.begin(), _dressedPFOs.end(), i) != _dressedPFOs.end()){
			continue;
		}

		TVector3 P_i( pfo_i->getMomentum() );
		float cosTheta = P.Dot( P_i )/(P.Mag()*P_i.Mag());
		if ( cosTheta >= _cosConeAngle )
			coneE += pfo_i->getEnergy();
	}

	return coneE;
}

void IsolatedLeptonFinderProcessor::getCalEnergy( ReconstructedParticle* pfo , float* cale) {
	float ecal = 0;
	float hcal = 0;
	std::vector<lcio::Cluster*> clusters = pfo->getClusters();
	for ( std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();
			iCluster!=clusters.end();
			++iCluster) {
		ecal += (*iCluster)->getSubdetectorEnergies()[0];
		hcal += (*iCluster)->getSubdetectorEnergies()[1];
	}
	cale[0] = ecal;
	cale[1] = hcal;
}


