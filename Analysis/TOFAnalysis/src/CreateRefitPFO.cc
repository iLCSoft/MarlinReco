#include "CreateRefitPFO.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/ParticleIDImpl.h"
#include "EVENT/MCParticle.h"
#include "GeometryUtil.h"
#include "marlin/VerbosityLevels.h"
#include "TMatrixD.h"
using namespace lcio;

CreateRefitPFO aCreateRefitPFO;

CreateRefitPFO::CreateRefitPFO():marlin::Processor("CreateRefitPFO"){
    _description = "Create a copy of PandoraPFO collection, where for each Kaon/Proton particle the tracks are substituted from MarlinTrkTracksKaon and MarlinTrkTracksProton colelctions respectively. This simulates scenario as all pi/k/p tracks are fitted with their true mass assuming perfect PID.";
}

void CreateRefitPFO::init(){
	_bField = MarlinUtil::getBzAtOrigin();
}

void CreateRefitPFO::processEvent( EVENT::LCEvent* event ){
	LCCollection* pfos = event->getCollection("PandoraPFOs");
    LCCollection* tracksCol = event->getCollection("MarlinTrkTracks");
	LCCollection* tracksKaonRefitCol = event->getCollection("MarlinTrkTracksKaon");
	LCCollection* tracksProtonRefitCol = event->getCollection("MarlinTrkTracksProton");
    LCRelationNavigator trackToMcNav( event->getCollection("MarlinTrkTracksMCTruthLink") );

	LCCollectionVec* outputPfosCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

	for(int i=0; i<pfos->getNumberOfElements(); ++i){
		ReconstructedParticle* pfo = static_cast<ReconstructedParticle*>( pfos->getElementAt(i) );
        const std::vector<Track*>& tracks = pfo->getTracks();
        int nTracks = tracks.size();

        //Substitute ALL Kaon/Proton tracks inside this PFO
        std::vector<Track*> refittedTracks;
        std::vector<int> tracksPdg;
        for(int j=0; j<nTracks; ++j){
            Track* track = tracks[j];
            int pdg = getTrackPDG(track, trackToMcNav);
            int idx = getTrackIndex(tracksCol, track);
            if ( std::abs(pdg) == 321 ) track = static_cast<Track*> ( tracksKaonRefitCol->getElementAt( idx ) );
            else if( std::abs(pdg) == 2212 ) track = static_cast<Track*> ( tracksProtonRefitCol->getElementAt( idx ) );
            tracksPdg.push_back(pdg);
            refittedTracks.push_back(track);
        }

        // create new PFO substituting tracks and 4mom and cov matrix in simple cases. Everything else just copy from the old pfo
        ReconstructedParticleImpl* outputPFO = new ReconstructedParticleImpl();
        for(size_t j=0; j<pfo->getTracks().size(); ++j) outputPFO->addTrack(refittedTracks[j]);

        if (nTracks == 1 && ( std::abs(tracksPdg[0]) == 321 || std::abs(tracksPdg[0]) == 2212 ) ){
            double mass = 0.13957039;
            if ( std::abs(tracksPdg[0]) == 321 ) mass = 0.493677;
            else if ( std::abs(tracksPdg[0]) == 2212 ) mass = 0.938272088;
        
            //do calculations
            std::vector<float> covMatrix = updateChargedPFOCovMat(refittedTracks[0], mass);
            TLorentzVector fourMomentum = getTrackFourMomentum(refittedTracks[0], mass);
            const double momentum[3] = {fourMomentum.Px(), fourMomentum.Py(), fourMomentum.Pz()};
        
            //update these parameters, copy everything else
            outputPFO->setType(tracksPdg[0]);
            outputPFO->setMomentum(momentum);
            outputPFO->setEnergy(fourMomentum.E());
            outputPFO->setCovMatrix(covMatrix);
            outputPFO->setMass(fourMomentum.M());
        }
        else if ( nTracks == 2 && pfo->getCharge() == 0. ){
            double mass1 = 0.13957039;
            if ( std::abs(tracksPdg[0]) == 321 ) mass1 = 0.493677;
            else if ( std::abs(tracksPdg[0]) == 2212 ) mass1 = 0.938272088;
            double mass2 = 0.13957039;
            if ( std::abs(tracksPdg[1]) == 321 ) mass2 = 0.493677;
            else if ( std::abs(tracksPdg[1]) == 2212 ) mass2 = 0.938272088;
        
            TLorentzVector fourMomentum1 = getTrackFourMomentum(refittedTracks[0], mass1);
            TLorentzVector fourMomentum2 = getTrackFourMomentum(refittedTracks[1], mass2);
            TLorentzVector fourMomentum = fourMomentum1 + fourMomentum2;
            const double momentum[3] = {fourMomentum.Px(), fourMomentum.Py(), fourMomentum.Pz()};
        
            std::vector<float> covMatrix1 = updateChargedPFOCovMat(refittedTracks[0], mass1);
            std::vector<float> covMatrix2 = updateChargedPFOCovMat(refittedTracks[1], mass2);
            std::vector<float> covMatrix;
            for(size_t j=0; j<covMatrix1.size(); ++j) covMatrix.push_back( covMatrix1[j] + covMatrix2[j] );
        
            //update these parameters, copy everything else
            outputPFO->setType(pfo->getType());
            outputPFO->setMomentum(momentum);
            outputPFO->setEnergy(fourMomentum.E());
            outputPFO->setCovMatrix(covMatrix);
            outputPFO->setMass(fourMomentum.M());
        }
        else{
            outputPFO->setType(pfo->getType());
            outputPFO->setMomentum(pfo->getMomentum());
            outputPFO->setEnergy(pfo->getEnergy());
            outputPFO->setCovMatrix(pfo->getCovMatrix());
            outputPFO->setMass(pfo->getMass());
        }
        outputPFO->setCharge(pfo->getCharge());
        outputPFO->setReferencePoint(pfo->getReferencePoint());
        outputPFO->setParticleIDUsed(pfo->getParticleIDUsed());
        outputPFO->setGoodnessOfPID(pfo->getGoodnessOfPID());
        outputPFO->setStartVertex(pfo->getStartVertex());
        for(size_t j=0; j<pfo->getParticles().size(); ++j) outputPFO->addParticle(pfo->getParticles()[j]);
        for(size_t j=0; j<pfo->getClusters().size(); ++j) outputPFO->addCluster(pfo->getClusters()[j]);
        for(size_t j=0; j<pfo->getParticleIDs().size(); ++j){
            ParticleIDImpl* inPID = static_cast<ParticleIDImpl*>( pfo->getParticleIDs()[j] );
            ParticleIDImpl* outPID = new ParticleIDImpl;
            outPID->setType(inPID->getType());
            outPID->setPDG(inPID->getPDG());
            outPID->setLikelihood(inPID->getLikelihood());
            outPID->setAlgorithmType(inPID->getAlgorithmType()) ;
            for(size_t k=0; k<inPID->getParameters().size(); ++k) outPID->addParameter( inPID->getParameters()[k] );
            outputPFO->addParticleID( outPID );
        }
        outputPfosCol->addElement( outputPFO );
    }
    event->addCollection( outputPfosCol , "updatedPandoraPFOs" );
}

int CreateRefitPFO::getTrackIndex(EVENT::LCCollection* trackCollection, EVENT::Track* selectedTrack){
	for (int i=0; i<trackCollection->getNumberOfElements(); ++i){
		Track* track = static_cast<Track*>( trackCollection->getElementAt(i) );
		if ( track == selectedTrack ) return i;
    }
    return -1;
}

int CreateRefitPFO::getTrackPDG(EVENT::Track* track, UTIL::LCRelationNavigator& nav){
    const std::vector<LCObject*>& particles = nav.getRelatedToObjects(track);
    const std::vector<float>&  weights = nav.getRelatedToWeights(track);
    if( particles.size() == 0 ) return 0;

    int i = std::max_element( weights.begin(), weights.end() ) - weights.begin();
	MCParticle* mc = static_cast<MCParticle*> (particles[i]);
    return mc->getPDG();
}

TLorentzVector CreateRefitPFO::getTrackFourMomentum(EVENT::Track* track, double mass){
    if ( track->getOmega() == 0. ) return TLorentzVector();
	double pt = 0.299792458e-3 * _bField / std::abs( track->getOmega() );
	double px = pt*std::cos( track->getPhi() );
	double py = pt*std::sin( track->getPhi() );
	double pz = pt*track->getTanLambda();
	double energy = std::sqrt(mass*mass + px*px + py*py + pz*pz);
	return TLorentzVector(px, py, pz, energy);
}

std::vector<float> CreateRefitPFO::updateChargedPFOCovMat(EVENT::Track* track, double mass){
    //	Obtain covariance matrix on (px,py,pz,E) from the covariance matrix on helix parameters.
    //			Dpx/DTanLambda		Dpy/DTanLambda		Dpz/DTanLambda		DE/DTanLambda
    //	 J =	Dpx/DPhi		      Dpy/DPhi		      Dpz/DPhi		      DE/DPhi
    //			Dpx/DOmega		     Dpy/DOmega		     Dpz/DOmega		     DE/DOmega
    //
    //			0			         0			     Pt	          Pz.Pt/E
    //	 J =	-Py			        Px			     0			    0
    //			-Px/Omega		-Py/Omega		-Pz/Omega		-P2/E.Omega
    //
    //	Order in the covariance matrix on helix parameters:
    //			TanLambda.TanLambda	       TanLambda.phi		   TanLambda.Omega
    //	Cov =	phi.TanLambda		         phi.phi			      phi.Omega
    //			Omega.TanLambda		        Omega.phi		         Omega.Omega
	const int rows = 3; // n rows jacobian
	const int columns = 4; // n columns jacobian

	double omega = track->getOmega();
    if (omega == 0.0) return std::vector<float>(10, 0); 

    double pt = 0.299792458e-3 * _bField / std::abs( track->getOmega() );
	double px = pt*std::cos( track->getPhi() );
	double py = pt*std::sin( track->getPhi() );
    double pz = pt*track->getTanLambda();
    double p = std::sqrt(pt*pt + pz*pz);

	double energy = std::sqrt(mass*mass + p*p);

    double jacobianByRows[rows*columns] = {
        0, 0, pt, pz*pt/energy,
		-py, px, 0,	0,
		-px/omega, -py/omega, -pz/omega, -p*p/(energy*omega)
    };
    // Create Jacobian ("F" fill by columns, "C" fill by rows)
    TMatrixD jacobian(rows, columns, jacobianByRows, "C");

    // we ignore d0, z0 parameters and using only phi, omega, tanL
    std::vector<float> flatTrackCovMat = track->getCovMatrix();
	double trackCovMatByRows[rows*rows] = {
        flatTrackCovMat[14], flatTrackCovMat[11], flatTrackCovMat[12],
		flatTrackCovMat[11], flatTrackCovMat[2], flatTrackCovMat[4],
		flatTrackCovMat[12], flatTrackCovMat[4], flatTrackCovMat[5]
	};
	TMatrixD trackCovMat(rows, rows, trackCovMatByRows, "C");

    //4momentum cov matrix and its diagonal 10 elements flattened into array
	TMatrixD covMat4Mom(4, 4);
    covMat4Mom.Mult( TMatrixD(jacobian, TMatrixD::kTransposeMult, trackCovMat), jacobian );

	std::vector<float> flatCovMat4Mom;
	flatCovMat4Mom.push_back( covMat4Mom(0,0) ); // x-x
	flatCovMat4Mom.push_back( covMat4Mom(1,0) ); // y-x
	flatCovMat4Mom.push_back( covMat4Mom(1,1) ); // y-y
	flatCovMat4Mom.push_back( covMat4Mom(2,0) ); // z-x
	flatCovMat4Mom.push_back( covMat4Mom(2,1) ); // z-y
	flatCovMat4Mom.push_back( covMat4Mom(2,2) ); // z-z
	flatCovMat4Mom.push_back( covMat4Mom(3,0) ); // e-x
	flatCovMat4Mom.push_back( covMat4Mom(3,1) ); // e-y
	flatCovMat4Mom.push_back( covMat4Mom(3,2) ); // e-z
	flatCovMat4Mom.push_back( covMat4Mom(3,3) ); // e-e

	return flatCovMat4Mom;
}
