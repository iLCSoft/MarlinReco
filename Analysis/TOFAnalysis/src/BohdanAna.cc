#include "BohdanAna.h"
#include "BohdanDrawing.h"
#include "TOF.h"
#include "TrackLength.h"

#include "marlin/Global.h"
#include "marlinutil/GeometryUtil.h"
#include "marlinutil/CalorimeterHitType.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/TrackTools.h"
#include "MarlinTrk/Factory.h"
#include "CLHEP/Random/Randomize.h"

using dd4hep::rec::Vector3D;
using ROOT::Math::PxPyPzEVector;

BohdanAna theBohdanAna;

BohdanAna::BohdanAna() : marlin::Processor("BohdanAna"), EventDisplayer(this){
    _description = "Main analysis of the track length and time-of-flight and momentum methods for time-of-flight pID";

    registerProcessorParameter( "produce_csv_output",
                                "Produce csv output file for Konrad or not",
                                _produce_csv_output,
                                false);

    registerProcessorParameter( "dst_mode",
                                "Write only DST level output",
                                _dst_mode,
                                false);

    // Rerunning LCFIPlus vertexing on statistics of DST files is very slow and is a bottleneck.
    registerProcessorParameter( "produce_refit_output",
                                "Store info about refitted tracks and vertices",
                                _produce_refit_output,
                                true);

}


void BohdanAna::init(){
    _bField = MarlinUtil::getBzAtOrigin();
    _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
    _trkSystem->init();

    initialiseTTree();

    std::vector<std::string> fileNames ;
    Global::parameters->getStringVals("LCIOInputFiles" , fileNames ) ;
    if( fileNames.empty() ) streamlog_out(ERROR)<<" No slcio files provided!"<<std::endl;
    for(auto& name : fileNames){
        if ( _dst_mode && name.find("rec") != std::string::npos  ){
            streamlog_out(WARNING)<<" Found \"rec\" in a filename while using DST mode! use DST files with full statistics!"<<std::endl;
            break;
        }
        else if ( not _dst_mode && name.find("dst") != std::string::npos ){
            streamlog_out(ERROR)<<" Found \"dst\" in a filename while using REC mode! Use DST mode or a REC file. Exiting to prevent a crash."<<std::endl;
            exit(1);
        }
    }

    if(_produce_csv_output){
        if (_dst_mode) streamlog_out(ERROR)<<" I cannot produce Konrad CSV file with hit information while running in DST mode!"<<std::endl;
        _csv_output_file = std::ofstream("output.csv");
        _csv_output_file<<"PFO #,"
                    "PDG,"
                    "trk length (mm),"
                    "trk p (GeV),"
                    "trk pT (GeV),"
                    "trk px (GeV),"
                    "trk py (GeV),"
                    "trk pz (GeV),"
                    "trk Ecal x (mm),"
                    "trk Ecal y (mm),"
                    "trk Ecal z (mm),"
                    "true TOF (ns),"
                    "Hit type,"
                    "Hit Layout,"
                    "true hit time (ns),"
                    "hit time 50ps (ns),"
                    "hit time 100ps (ns),"
                    "hit energy (GeV),"
                    "hit layer,"
                    "x hit pos (mm),"
                    "y hit pos (mm),"
                    "z hit pos (mm)\n";
    }
}

void BohdanAna::fillMCTrueInfo(EVENT::MCParticle* mc, EVENT::ReconstructedParticle* pfo){
    _pdg = mc->getPDG();
    if ( not mc->getParents().empty() ){
        _imidiateParentPDG = mc->getParents()[0]->getPDG();
        MCParticle* firstStableParent = getFirstStableParent(mc);
        if ( firstStableParent != nullptr ) _firstStableParentPDG = firstStableParent->getPDG();
    }
    _mcVtx = Vector3D( mc->getVertex() );
    _timeTrue = mc->getTime();

    _generatorStatus = mc->getGeneratorStatus();
    _isSimulated = mc->isCreatedInSimulation();
    _isBackscatter = mc->isBackscatter();
    _vertexIsNotEndpointOfParent = mc->vertexIsNotEndpointOfParent();
    _isDecayedInTracker = mc->isDecayedInTracker();
    _isDecayedInCalorimeter = mc->isDecayedInCalorimeter();
    _hasLeftDetector = mc->hasLeftDetector();
    _isStopped = mc->isStopped();
    _isOverlay = mc->isOverlay();

    unsigned int quarkTypeDecay = getQuarkTypeDecay(mc);
    _isBottomQuarkDecay = quarkTypeDecay == 5;
    _isCharmQuarkDecay = quarkTypeDecay == 4;
    _isHadronisationDecay = not _isBottomQuarkDecay && not _isCharmQuarkDecay;
    _isV0DecayTrue = std::abs(_firstStableParentPDG) == 310  || std::abs(_firstStableParentPDG) == 3122;

    _isReconstructed = pfo != nullptr;
    if ( _isReconstructed ){
        // Have a WELL-DEFINED track/cluster. I do not store reco infromation of weird PFOs with 2 tracks/showers...
        _hasTrack = pfo->getTracks().size() == 1;
        _hasShower = pfo->getClusters().size() == 1;
    }
    _omegaTrue = getOmegaTrue(mc);
    _tanLambdaTrue = getTanLTrue(mc);
    _d0True = getD0True(mc);
    _z0True = getZ0True(mc);
    _phiTrue = getPhiTrue(mc);
    _mcMom = Vector3D( mc->getMomentum() );


    streamlog_out(DEBUG5)<<"Checking MC particle information"<<std::endl;
    streamlog_out(DEBUG5)<<"_pdg: "<<_pdg<<std::endl;
    streamlog_out(DEBUG5)<<"_imidiateParentPDG: "<<_imidiateParentPDG<<std::endl;
    streamlog_out(DEBUG5)<<"_firstStableParentPDG: "<<_firstStableParentPDG<<std::endl;
    streamlog_out(DEBUG5)<<"_vtxPos: ("<<_mcVtx[0]<<", "<<_mcVtx[1]<<", "<<_mcVtx[2]<<")"<<std::endl;
    streamlog_out(DEBUG5)<<"_isOverlay: "<<_isOverlay<<std::endl;
    streamlog_out(DEBUG5)<<"_isSimulated: "<<_isSimulated<<std::endl;
    streamlog_out(DEBUG5)<<"_isBottomQuarkDecay: "<<_isBottomQuarkDecay<<std::endl;
    streamlog_out(DEBUG5)<<"_isCharmQuarkDecay: "<<_isCharmQuarkDecay<<std::endl;
    streamlog_out(DEBUG5)<<"_isHadronisationDecay: "<<_isHadronisationDecay<<std::endl;
    streamlog_out(DEBUG5)<<"_isV0DecayTrue: "<<_isV0DecayTrue<<std::endl;
    streamlog_out(DEBUG5)<<"_isReconstructed: "<<_isReconstructed<<std::endl;
    streamlog_out(DEBUG5)<<"_hasTrack: "<<_hasTrack<<std::endl;
    streamlog_out(DEBUG5)<<"_hasShower: "<<_hasShower<<std::endl;

    streamlog_out(DEBUG5)<<"_d0True: "<<_d0True<<std::endl;
    streamlog_out(DEBUG5)<<"_z0True: "<<_z0True<<std::endl;
    streamlog_out(DEBUG5)<<"_omegaTrue: "<<_omegaTrue<<std::endl;
    streamlog_out(DEBUG5)<<"_tanLambdaTrue: "<<_tanLambdaTrue<<std::endl;
    streamlog_out(DEBUG5)<<"_timeTrue: "<<_timeTrue<<std::endl;
    streamlog_out(DEBUG5)<<"_vtxMom: ("<<_mcMom[0]<<", "<<_mcMom[1]<<", "<<_mcMom[2]<<")"<<std::endl;

}


void BohdanAna::fillTrueVertexInfo(EVENT::MCParticle* mc, const std::vector<VertexData>& trueVertices){
    for(auto& vtx : trueVertices) {
        auto it = std::find(vtx.mcs.begin(), vtx.mcs.end(), mc);
        if (it == vtx.mcs.end()) continue;
        if (vtx.isPrimary) _isInTruePrimaryVertex = true;
        else _isInTrueSecondaryVertex = true;
        _isInTrueV0Vertex = vtx.isV0;
        _trueVertexPos = vtx.pos;
        _nTracksAtTrueVertex = vtx.mcs.size();

        PxPyPzEVector fourMomentum{};
        for (auto* vtxMc : vtx.mcs){
            fourMomentum += PxPyPzEVector( vtxMc->getMomentum()[0], vtxMc->getMomentum()[1], vtxMc->getMomentum()[2], vtxMc->getEnergy() );
        }
        _invMassOfTrueVertex = fourMomentum.M();
        auto compEnergy = [](MCParticle* a, MCParticle* b){ return a->getEnergy() < b->getEnergy(); };
        _minEnergyOfTrueVertex = (*std::min_element(vtx.mcs.begin(), vtx.mcs.end(), compEnergy ))->getEnergy();

        // These are used for V0 check only if there are exactly two tracks
        if (_nTracksAtTrueVertex == 2){
            _oppositeChargeOfTrueVertex = vtx.mcs[0]->getCharge()*vtx.mcs[1]->getCharge() < 0.f ;
            Vector3D mom1( vtx.mcs[0]->getMomentum() );
            Vector3D mom2( vtx.mcs[1]->getMomentum() );
            // Make sure fillMCTrueInfo() run first to fill _mcVtx, _ipTrue!!! It is very bad coding, sorry...
            _cosThetaOfTrueVertex = (_mcVtx - _ipTrue).unit().dot((mom1+mom2).unit());
            // Assuming IP is at 0,0,0
            _cosThetaToIpOfTrueVertex = (_mcVtx).unit().dot((mom1+mom2).unit());
        }
        return;
    }
}


void BohdanAna::fillRecoVertexInfo(EVENT::LCEvent* evt, EVENT::MCParticle* mc, const UTIL::LCRelationNavigator& pfo2mc){
    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCCollection* primVtxCol = evt->getCollection("PrimaryVertex");
    LCCollection* secondaryVtxCol = evt->getCollection("BuildUpVertex");
    LCCollection* secondaryV0VtxCol = evt->getCollection("BuildUpVertex_V0");

    EVENT::Vertex* matchedVertex = nullptr;
    for(int i=0; i<primVtxCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (primVtxCol->getElementAt(i));
        auto vtxPfos = vertex->getAssociatedParticle()->getParticles();
        auto it = std::find_if(vtxPfos.begin(), vtxPfos.end(), [&](EVENT::ReconstructedParticle* pfo) {return getMC(pfo, pfo2mc) == mc; } );
        if ( it == vtxPfos.end() ) continue;
        _isInRecoPrimaryVertex = true;
        matchedVertex = vertex;
        break;
    }

    // NOTE: Ideally MC particle must be only in secondary or primary vertex. However, due to the bug cleaning up tracks from primary vertex it can appear in both. In this case I write the properties of the secondary vertex as the most interesting one, but still mark recoInPrimary as true.
    for(int i=0; i<secondaryVtxCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (secondaryVtxCol->getElementAt(i));
        auto vtxPfos = vertex->getAssociatedParticle()->getParticles();
        auto it = std::find_if(vtxPfos.begin(), vtxPfos.end(), [&](EVENT::ReconstructedParticle* pfo) {return getMC(pfo, pfo2mc) == mc; } );
        if ( it == vtxPfos.end() ) continue;
        _isInRecoSecondaryVertex = true;
        matchedVertex = vertex;
        break;
    }

    for(int i=0; i<secondaryV0VtxCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (secondaryV0VtxCol->getElementAt(i));
        auto vtxPfos = vertex->getAssociatedParticle()->getParticles();
        auto it = std::find_if(vtxPfos.begin(), vtxPfos.end(), [&](EVENT::ReconstructedParticle* pfo) {return getMC(pfo, pfo2mc) == mc; } );
        if ( it == vtxPfos.end() ) continue;
        _isInRecoSecondaryVertex = true;
        _isV0DecayReco = true;
        matchedVertex = vertex;
        break;
    }

    if ( matchedVertex != nullptr){
        auto prongs = matchedVertex->getAssociatedParticle()->getParticles();

        _recoVertexPos = Vector3D( matchedVertex->getPosition() );
        float sigmaX = std::sqrt( matchedVertex->getCovMatrix()[0] );
        float sigmaY = std::sqrt( matchedVertex->getCovMatrix()[2] );
        float sigmaZ = std::sqrt( matchedVertex->getCovMatrix()[5] );
        _recoVertexPosErr = Vector3D( sigmaX, sigmaY, sigmaZ );
        _nTracksAtRecoVertex = prongs.size();
        PxPyPzEVector fourMomentum{};
        for (auto* pfo : prongs){
            fourMomentum += PxPyPzEVector( pfo->getMomentum()[0], pfo->getMomentum()[1], pfo->getMomentum()[2], pfo->getEnergy() );
        }
        _invMassOfRecoVertex = fourMomentum.M();
        auto compEnergy = [](ReconstructedParticle* a, ReconstructedParticle* b){ return a->getEnergy() < b->getEnergy(); };
        _minEnergyOfRecoVertex = (*std::min_element(prongs.begin(), prongs.end(), compEnergy ))->getEnergy();
        _chi2OfRecoVertex = matchedVertex->getChi2();

        // These are used for V0 check only if there are exactly two tracks
        if (_nTracksAtRecoVertex == 2){
            _oppositeChargeOfRecoVertex = prongs[0]->getCharge()*prongs[1]->getCharge() < 0.f ;
            Vector3D mom1( prongs[0]->getMomentum() );
            Vector3D mom2( prongs[1]->getMomentum() );
            // Make sure _recoVertexPos and _primVertexReco are filled first!!! It is very bad coding, sorry...
            _cosThetaOfRecoVertex = (_recoVertexPos - _primVertexReco).unit().dot((mom1+mom2).unit());
            // Assuming IP is at 0,0,0
            _cosThetaToIpOfRecoVertex = (_recoVertexPos).unit().dot((mom1+mom2).unit());
        }
    }

    if ( not _produce_refit_output ) return;


    // Do literaly the same for the refitted vertex collections. Should be refactored better in functions...
    LCCollection* updatedPfos = evt->getCollection("updatedPandoraPFOs");
    LCCollection* primVtxRefitCol = evt->getCollection("PrimaryVertex_refit");
    LCCollection* secondaryVtxRefitCol = evt->getCollection("BuildUpVertex_refit");
    LCCollection* secondaryV0VtxRefitCol = evt->getCollection("BuildUpVertex_V0_refit");
    EVENT::Vertex* matchedRefittedVertex = nullptr;

    for(int i=0; i<primVtxRefitCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (primVtxRefitCol->getElementAt(i));
        auto refittedPfos = vertex->getAssociatedParticle()->getParticles();
        std::vector<ReconstructedParticle*> vtxPfos;
        // NOTE: I need to find corresponding PFOs from non-refitted collection as only they have links to the MC particles.
        for (auto* refittedPfo : refittedPfos){
            auto* pfoObject = getMatchingElement(updatedPfos, refittedPfo, pfos);
            if( pfoObject == nullptr ) continue;
            auto* pfo = static_cast<ReconstructedParticle*> ( pfoObject );
            vtxPfos.push_back(pfo);
        }
        auto it = std::find_if(vtxPfos.begin(), vtxPfos.end(), [&](EVENT::ReconstructedParticle* pfo) {return getMC(pfo, pfo2mc) == mc; } );
        if ( it == vtxPfos.end() ) continue;
        _isInRecoPrimaryRefitVertex = true;
        matchedRefittedVertex = vertex;
        break;
    }

    for(int i=0; i<secondaryVtxRefitCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (secondaryVtxRefitCol->getElementAt(i));
        auto refittedPfos = vertex->getAssociatedParticle()->getParticles();
        std::vector<ReconstructedParticle*> vtxPfos;
        // NOTE: I need to find corresponding PFOs from non-refitted collection as only they have links to the MC particles.
        for (auto* refittedPfo : refittedPfos){
            auto* pfoObject = getMatchingElement(updatedPfos, refittedPfo, pfos);
            if( pfoObject == nullptr ) continue;
            auto* pfo = static_cast<ReconstructedParticle*> ( pfoObject );
            vtxPfos.push_back(pfo);
        }
        auto it = std::find_if(vtxPfos.begin(), vtxPfos.end(), [&](EVENT::ReconstructedParticle* pfo) {return getMC(pfo, pfo2mc) == mc; } );
        if ( it == vtxPfos.end() ) continue;
        _isInRecoSecondaryRefitVertex = true;
        matchedRefittedVertex = vertex;
        break;
    }

    for(int i=0; i<secondaryV0VtxRefitCol->getNumberOfElements(); ++i){
        auto vertex = static_cast<Vertex*> (secondaryV0VtxRefitCol->getElementAt(i));
        auto refittedPfos = vertex->getAssociatedParticle()->getParticles();
        std::vector<ReconstructedParticle*> vtxPfos;
        // NOTE: I need to find corresponding PFOs from non-refitted collection as only they have links to the MC particles.
        for (auto* refittedPfo : refittedPfos){
            auto* pfoObject = getMatchingElement(updatedPfos, refittedPfo, pfos);
            if( pfoObject == nullptr ) continue;
            auto* pfo = static_cast<ReconstructedParticle*> ( pfoObject );
            vtxPfos.push_back(pfo);
        }
        auto it = std::find_if(vtxPfos.begin(), vtxPfos.end(), [&](EVENT::ReconstructedParticle* pfo) {return getMC(pfo, pfo2mc) == mc; } );
        if ( it == vtxPfos.end() ) continue;
        _isInRecoSecondaryRefitVertex = true;
        _isV0DecayRefitReco = true;
        matchedRefittedVertex = vertex;
        break;
    }

    if ( matchedRefittedVertex != nullptr ){
        auto refittedProngs = matchedRefittedVertex->getAssociatedParticle()->getParticles();

        _recoRefitVertexPos = Vector3D( matchedRefittedVertex->getPosition() );
        float refittedSigmaX = std::sqrt( matchedRefittedVertex->getCovMatrix()[0] );
        float refittedSigmaY = std::sqrt( matchedRefittedVertex->getCovMatrix()[2] );
        float refittedSigmaZ = std::sqrt( matchedRefittedVertex->getCovMatrix()[5] );
        _recoRefitVertexPosErr = Vector3D( refittedSigmaX, refittedSigmaY, refittedSigmaZ );
        _nTracksAtRecoRefitVertex = refittedProngs.size();
        PxPyPzEVector refittedFourMomentum{};
        for (auto* pfo : refittedProngs){
            refittedFourMomentum += PxPyPzEVector( pfo->getMomentum()[0], pfo->getMomentum()[1], pfo->getMomentum()[2], pfo->getEnergy() );
        }
        _invMassOfRecoRefitVertex = refittedFourMomentum.M();
        auto compEnergy = [](ReconstructedParticle* a, ReconstructedParticle* b){ return a->getEnergy() < b->getEnergy(); };
        _minEnergyOfRecoRefitVertex = (*std::min_element(refittedProngs.begin(), refittedProngs.end(), compEnergy ))->getEnergy();
        _chi2OfRecoRefitVertex = matchedRefittedVertex->getChi2();

        // These are used for V0 check only if there are exactly two tracks
        if (_nTracksAtRecoRefitVertex == 2){
            _oppositeChargeOfRecoRefitVertex = refittedProngs[0]->getCharge()*refittedProngs[1]->getCharge() < 0.f ;
            Vector3D mom1( refittedProngs[0]->getMomentum() );
            Vector3D mom2( refittedProngs[1]->getMomentum() );
            // Make sure _recoRefitVertexPos and _primRefitVertexReco are filled first!!! It is very bad coding, sorry...
            _cosThetaOfRecoRefitVertex = (_recoRefitVertexPos - _primRefitVertexReco).unit().dot((mom1+mom2).unit());
            // Assuming IP is at 0,0,0
            _cosThetaToIpOfRecoRefitVertex = (_recoRefitVertexPos).unit().dot((mom1+mom2).unit());
        }
    }
}



void BohdanAna::fillTrackStates(EVENT::ReconstructedParticle* pfo, EVENT::ReconstructedParticle* refittedPfo){
    if (pfo == nullptr || pfo->getTracks().size() != 1) return;
    auto track = pfo->getTracks()[0];

    auto tsIP = track->getTrackState( TrackState::AtIP );
    _omegaIP = tsIP->getOmega();
    _tanLambdaIP = tsIP->getTanLambda();
    _d0IP = tsIP->getD0();
    _z0IP = tsIP->getZ0();
    _phiIP = tsIP->getPhi();
    _omegaErrIP = std::sqrt( tsIP->getCovMatrix()[5] );
    _tanLambdaErrIP = std::sqrt( tsIP->getCovMatrix()[14] );
    _d0ErrIP = std::sqrt( tsIP->getCovMatrix()[0] );
    _z0ErrIP = std::sqrt( tsIP->getCovMatrix()[9] );
    _phiErrIP = std::sqrt( tsIP->getCovMatrix()[2] );
    // _recoIpPos is always (0, 0, 0) so I do not store it...
    auto momIP = UTIL::getTrackMomentum(tsIP, _bField);
    _recoIpMom = Vector3D( momIP[0], momIP[1], momIP[2] ); // should be identical to pfo->getMomentum()

    auto tsECAL = getTrackStateAtCalorimeter( track );
    _omegaECAL = tsECAL->getOmega();
    _tanLambdaECAL = tsECAL->getTanLambda();
    _d0ECAL = tsECAL->getD0(); // NOTE: must be 0 by definition at ECAL
    _z0ECAL = tsECAL->getZ0(); // NOTE: must be 0 by definition at ECAL
    _phiECAL = tsECAL->getPhi();
    _omegaErrECAL = std::sqrt( tsECAL->getCovMatrix()[5] );
    _tanLambdaErrECAL = std::sqrt( tsECAL->getCovMatrix()[14] );
    _d0ErrECAL = std::sqrt( tsECAL->getCovMatrix()[0] );
    _z0ErrECAL = std::sqrt( tsECAL->getCovMatrix()[9] );
    _phiErrECAL = std::sqrt( tsECAL->getCovMatrix()[2] );
    _recoCaloPos = Vector3D( tsECAL->getReferencePoint() );
    auto momECAL = UTIL::getTrackMomentum(tsECAL, _bField);
    _recoCaloMom = Vector3D( momECAL[0], momECAL[1], momECAL[2] );

    if (not _produce_refit_output) return;

    if (refittedPfo == nullptr || refittedPfo->getTracks().size() != 1 ||  refittedPfo->getTracks()[0]->getOmega() == 0.f || refittedPfo->getTracks()[0]->getTanLambda() == 0.f) return;

    auto refittedTrack = refittedPfo->getTracks()[0];

    auto refittedTsIP = refittedTrack->getTrackState( TrackState::AtIP );
    _refittedOmegaIP = refittedTsIP->getOmega();
    _refittedTanLambdaIP = refittedTsIP->getTanLambda();
    _refittedD0IP = refittedTsIP->getD0();
    _refittedZ0IP = refittedTsIP->getZ0();
    _refittedPhiIP = refittedTsIP->getPhi();
    _refittedOmegaErrIP = std::sqrt( refittedTsIP->getCovMatrix()[5] );
    _refittedTanLambdaErrIP = std::sqrt( refittedTsIP->getCovMatrix()[14] );
    _refittedD0ErrIP = std::sqrt( refittedTsIP->getCovMatrix()[0] );
    _refittedZ0ErrIP = std::sqrt( refittedTsIP->getCovMatrix()[9] );
    _refittedPhiErrIP = std::sqrt( refittedTsIP->getCovMatrix()[2] );
    // _recoIpPos is always (0, 0, 0)
    auto refittedMomIP = UTIL::getTrackMomentum(refittedTsIP, _bField);
    _refittedRecoIpMom = Vector3D( refittedMomIP[0], refittedMomIP[1], refittedMomIP[2] ); // should be identical to refittedPfo->getMomentum();

    auto refittedTsECAL = getTrackStateAtCalorimeter( refittedTrack );
    _refittedOmegaECAL = refittedTsECAL->getOmega();
    _refittedTanLambdaECAL = refittedTsECAL->getTanLambda();
    _refittedD0ECAL = refittedTsECAL->getD0(); // NOTE: must be 0 by definition at ECAL
    _refittedZ0ECAL = refittedTsECAL->getZ0(); // NOTE: must be 0 by definition at ECAL
    _refittedPhiECAL = refittedTsECAL->getPhi();
    _refittedOmegaErrECAL = std::sqrt( refittedTsECAL->getCovMatrix()[5] );
    _refittedTanLambdaErrECAL = std::sqrt( refittedTsECAL->getCovMatrix()[14] );
    _refittedD0ErrECAL = std::sqrt( refittedTsECAL->getCovMatrix()[0] );
    _refittedZ0ErrECAL = std::sqrt( refittedTsECAL->getCovMatrix()[9] );
    _refittedPhiErrECAL = std::sqrt( refittedTsECAL->getCovMatrix()[2] );
    _refittedRecoCaloPos = Vector3D( refittedTsECAL->getReferencePoint() );
    auto refittedMomECAL = UTIL::getTrackMomentum(refittedTsECAL, _bField);
    _refittedRecoCaloMom = Vector3D( refittedMomECAL[0], refittedMomECAL[1], refittedMomECAL[2] );
}


void BohdanAna::fillTrackLengthInfo(EVENT::MCParticle* mc, EVENT::ReconstructedParticle* pfo, const UTIL::LCRelationNavigator& navToSimTrackerHits){
    if (pfo == nullptr || pfo->getTracks().size() != 1) return;
    auto track = pfo->getTracks()[0];
    _trackLength_IDR = getTrackLengthIDR(track);
    _trackLength_SHA_phiLambda_IP = getTrackLengthSHA(track, TrackState::AtIP, TrackLengthOption::phiLambda);
    _trackLength_SHA_phiZed_IP = getTrackLengthSHA(track, TrackState::AtIP, TrackLengthOption::phiZed);
    _trackLength_SHA_zedLambda_IP = getTrackLengthSHA(track, TrackState::AtIP, TrackLengthOption::zedLambda);
    _trackLength_SHA_phiLambda_ECAL = getTrackLengthSHA(track, TrackState::AtCalorimeter, TrackLengthOption::phiLambda);
    _trackLength_SHA_phiZed_ECAL = getTrackLengthSHA(track, TrackState::AtCalorimeter, TrackLengthOption::phiZed);
    _trackLength_SHA_zedLambda_ECAL = getTrackLengthSHA(track, TrackState::AtCalorimeter, TrackLengthOption::zedLambda);

    std::vector<HitState> trackHitStates = getTrackStates(pfo, _bField, _trkSystem, navToSimTrackerHits);
    std::vector<IMPL::TrackStateImpl> trackStates;
    for(auto& hitState: trackHitStates){
        trackStates.push_back(hitState.ts);

        if ( hitState.simHit != nullptr && hitState.simHit->getMCParticle() != mc ) _cleanTrack = false;
    }

    std::tie(_trackLength_IKF_phiLambda, _harmonicMom_IKF_phiLambda) = getTrackLengthIKF(trackStates, _bField, TrackLengthOption::phiLambda);
    std::tie(_trackLength_IKF_phiZed, _harmonicMom_IKF_phiZed) = getTrackLengthIKF(trackStates, _bField, TrackLengthOption::phiZed);
    std::tie(_trackLength_IKF_zedLambda, _harmonicMom_IKF_zedLambda) = getTrackLengthIKF(trackStates, _bField, TrackLengthOption::zedLambda);

    auto it = std::find_if(trackHitStates.begin(), trackHitStates.end(), [](const HitState& hitState){return isSETHit(hitState.hit);});
    bool foundSETHit = it != trackHitStates.end();
    if (foundSETHit){
        int idx = it - trackHitStates.begin();
        std::vector<IMPL::TrackStateImpl> trackStatesToSET = std::vector<IMPL::TrackStateImpl>(trackStates.begin(), trackStates.begin() + idx + 1 );
        std::tie(_trackLengthToSET_IKF_phiLambda, _harmonicMomToSET_IKF_phiLambda) = getTrackLengthIKF(trackStatesToSET, _bField, TrackLengthOption::phiLambda);
        std::tie(_trackLengthToSET_IKF_phiZed, _harmonicMomToSET_IKF_phiZed) = getTrackLengthIKF(trackStatesToSET, _bField, TrackLengthOption::phiZed);
        std::tie(_trackLengthToSET_IKF_zedLambda, _harmonicMomToSET_IKF_zedLambda) = getTrackLengthIKF(trackStatesToSET, _bField, TrackLengthOption::zedLambda);
    }

}


void BohdanAna::fillTOFInfo(EVENT::MCParticle* mc, EVENT::ReconstructedParticle* pfo, const UTIL::LCRelationNavigator& navToSimCalorimeterHits){
    if (pfo == nullptr || pfo->getTracks().size() != 1 || pfo->getClusters().size() != 1) return;

    auto tsCalo = getTrackStateAtCalorimeter( pfo->getTracks()[0] );
    Vector3D trackPosAtCalo( tsCalo->getReferencePoint() );
    auto closestHit = getClosestHit( pfo->getClusters()[0], trackPosAtCalo);

    _typeClosest = getHitCaloType(closestHit);
    _caloIDClosest = getHitCaloID(closestHit);
    _layoutClosest = getHitCaloLayout(closestHit);
    _layerClosest = getHitCaloLayer(closestHit);
    _cleanClosestHit = getHitEarliestMC(closestHit, navToSimCalorimeterHits) == mc;

    //Do not fill TOF info when closest hit is in LumiCal
    if ( not isEcalHit(closestHit) ){
        return;
    }

    auto mom = UTIL::getTrackMomentum(tsCalo, _bField);
    Vector3D trackMomAtCalo( mom[0], mom[1], mom[2] );

    auto selectedFrankHits = selectFrankEcalHits( pfo->getClusters()[0] , trackPosAtCalo, trackMomAtCalo, 10);

    for (size_t j = 0; j < _resolutions.size(); j++){
        float res = _resolutions[j]/1000.; // in ns
        _tofClosest.at(j) = getHitTof(closestHit, trackPosAtCalo, res);
        _tofAverage.at(j) = getTofFrankAvg(selectedFrankHits, trackPosAtCalo, res);
        _tofFit.at(j) = getTofFrankFit(selectedFrankHits, trackPosAtCalo, res);
        std::tie(_tofSETFront.at(j), _tofSETBack.at(j)) = getTofSET( pfo->getTracks()[0], res);
    }

    for (auto* hit:  pfo->getClusters()[0]->getCalorimeterHits()){
        //Count only ECAL hits. No LumiCal, BeamCal, HCAL, Yoke hits are recorded!
        if (not isEcalHit(hit) ) continue;
        _xHit.push_back(hit->getPosition()[0]);
        _yHit.push_back(hit->getPosition()[1]);
        _zHit.push_back(hit->getPosition()[2]);
        _tHit.push_back(hit->getTime());
        _layerHit.push_back( CHT( hit->getType() ).layer() );
        _energyHit.push_back( hit->getEnergy() );
    }
    _nHits = _tHit.size();
}


void BohdanAna::fillCsvForKonrad(EVENT::ReconstructedParticle* pfo, int pdg, double trueTOF, double trackLength){
    if (pfo == nullptr || pfo->getTracks().size() != 1 || pfo->getClusters().size() != 1) return;

    auto tsCalo = getTrackStateAtCalorimeter( pfo->getTracks()[0] );
    Vector3D trackPosAtCalo( tsCalo->getReferencePoint() );
    auto mom = UTIL::getTrackMomentum(tsCalo, _bField);
    Vector3D trackMomAtCalo( mom[0], mom[1], mom[2] );

    _global_pfo_number++;
    for (auto* hit : pfo->getClusters()[0]->getCalorimeterHits()){
        //Count only ECAL hits. No LumiCal, BeamCal, HCAL, Yoke hits are recorded!
        if ( not isEcalHit(hit) ) continue;

        auto hitCaloID = getHitCaloID(hit);
        auto hitLayout = getHitCaloLayout(hit);
        auto hitLayer = getHitCaloLayer(hit);

        std::stringstream ss;
        ss<<_global_pfo_number<<", ";
        ss<<pdg<<", ";
        ss<<std::scientific<<std::setprecision(5)<<trackLength<<", ";

        ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.r()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.trans()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.x()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.y()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackMomAtCalo.z()<<", ";

        ss<<std::scientific<<std::setprecision(4)<<trackPosAtCalo.x()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackPosAtCalo.y()<<", ";
        ss<<std::scientific<<std::setprecision(4)<<trackPosAtCalo.z()<<", ";
        ss<<std::scientific<<std::setprecision(6)<<trueTOF<<", ";
        ss<<hitCaloID<<", "; // type of hit, i.e. ECal, LumiCal etc. hit
        ss<<hitLayout<<", ";
        ss<<std::scientific<<std::setprecision(6)<<hit->getTime()<<", ";
        ss<<std::scientific<<std::setprecision(6)<<CLHEP::RandGauss::shoot(hit->getTime(), 0.05)<<", ";
        ss<<std::scientific<<std::setprecision(6)<<CLHEP::RandGauss::shoot(hit->getTime(), 0.1)<<", ";
        ss<<std::scientific<<std::setprecision(6)<<hit->getEnergy()<<", ";
        ss<<hitLayer<<", ";
        ss<<std::scientific<<std::setprecision(6)<<hit->getPosition()[0]<<", ";
        ss<<std::scientific<<std::setprecision(6)<<hit->getPosition()[1]<<", ";
        ss<<std::scientific<<std::setprecision(6)<<hit->getPosition()[2]<<"\n";

        _csv_output_file << ss.str();
    }
}

void BohdanAna::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(MESSAGE)<<"==================== Event: "<<_nEvent<<std::endl;

    LCCollection* mcs = getMCParticleCollection(evt);
    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCRelationNavigator mc2pfo ( evt->getCollection("MCTruthRecoLink") );
    LCRelationNavigator pfo2mc ( evt->getCollection("RecoMCTruthLink") );
    LCCollection* primVtxCol = evt->getCollection("PrimaryVertex");
    LCCollection* updatedPfos = nullptr;
    LCCollection* primVtxRefitCol = nullptr;
    if ( _produce_refit_output ){
        updatedPfos = evt->getCollection("updatedPandoraPFOs");
        primVtxRefitCol = evt->getCollection("PrimaryVertex_refit");
    }
    std::unique_ptr<LCRelationNavigator> navToSimTrackerHits = nullptr;
    std::unique_ptr<LCRelationNavigator> navToSimCalorimeterHits = nullptr;
    if ( not _dst_mode){
        navToSimTrackerHits = std::make_unique<LCRelationNavigator> ( evt->getCollection("TrackerHitsRelations") );
        navToSimCalorimeterHits = std::make_unique<LCRelationNavigator> ( evt->getCollection("CalorimeterHitsRelations") );
    }

    std::vector<VertexData> trueVertices = getReconstructableTrueVertices(evt);
    _crossSection = evt->getParameters().getFloatVal("crossSection");
    _quarksToPythia = getQuarksToPythia(evt);
    _isHiggsProcess = checkIsHiggsProcess(evt);
    _higgsDaughters = getHiggsDaughters(evt);
    _primVertexTrue = Vector3D();
    _primVertexReco = Vector3D();
    _primRefitVertexReco = Vector3D();

    _ipTrue = Vector3D( static_cast< EVENT::MCParticle* >(mcs->getElementAt(0))->getVertex() );

    auto primVertexIt = std::find_if(trueVertices.begin(), trueVertices.end(), [](const VertexData& vtx){return vtx.isPrimary;});
    if ( primVertexIt != trueVertices.end() ) _primVertexTrue = (*primVertexIt).pos;

    if ( primVtxCol->getNumberOfElements() > 0 ){
        auto vtx = static_cast<Vertex*> (primVtxCol->getElementAt(0));
        _primVertexReco = Vector3D( vtx->getPosition() );
    }

    if (_produce_refit_output){
        if ( primVtxRefitCol->getNumberOfElements() > 0 ){
            auto vtx = static_cast<Vertex*> (primVtxRefitCol->getElementAt(0));
            _primRefitVertexReco = Vector3D( vtx->getPosition() );
        }
    }

    _nPionsGenerated = 0;
    _nKaonsGenerated = 0;
    _nProtonsGenerated = 0;
    _nPionsTrack = 0;
    _nKaonsTrack = 0;
    _nProtonsTrack = 0;
    _nPionsShower = 0;
    _nKaonsShower = 0;
    _nProtonsShower = 0;

    //Loop over MC particles
    for (int i=0; i<mcs->getNumberOfElements(); ++i){
        // Storing only charged hadrons - pi/K/p
        MCParticle* mc = static_cast <MCParticle*> ( mcs->getElementAt(i) );
        if (mc == nullptr) continue;
        bool isPion = std::abs( mc->getPDG() ) == 211;
        bool isKaon = std::abs( mc->getPDG() ) == 321;
        bool isProton = std::abs( mc->getPDG() ) == 2212;
        bool isHadron = isPion || isKaon || isProton;
        if ( not isHadron ) continue;
        if (isPion) _nPionsGenerated++;
        if (isKaon) _nKaonsGenerated++;
        if (isProton) _nProtonsGenerated++;

        resetVariables();

        //IMPORTANT: Run exactly in the order below there are interdependencies in fill functions...
        auto pfo = getRelatedReconstructedParticle(mc, mc2pfo, pfo2mc);
        fillMCTrueInfo( mc, pfo );
        fillTrueVertexInfo( mc, trueVertices );

        // Work only with simple reconstructed particles (1 track and 1 shower). Ignore rest ( ~ O( 0.1%) )
        if ( pfo == nullptr || pfo->getTracks().size() != 1){
            _tree->Fill();
            continue;
        }
        if (isPion) _nPionsTrack++;
        if (isKaon) _nKaonsTrack++;
        if (isProton) _nProtonsTrack++;

        _dEdx = pfo->getTracks()[0]->getdEdx();

        fillRecoVertexInfo(evt, mc, pfo2mc);

        ReconstructedParticle* refittedPFO = nullptr;
        if ( _produce_refit_output ) refittedPFO = static_cast<ReconstructedParticle* > ( getMatchingElement(pfos, pfo, updatedPfos) );

        fillTrackStates(pfo, refittedPFO);

        if (not _dst_mode) fillTrackLengthInfo(mc, pfo, *navToSimTrackerHits);

        if ( pfo->getClusters().size() != 1){
            _tree->Fill();
            continue;
        }
        if (isPion) _nPionsShower++;
        if (isKaon) _nKaonsShower++;
        if (isProton) _nProtonsShower++;

        if (not _dst_mode){
            fillTOFInfo(mc, pfo, *navToSimCalorimeterHits);
            if(_produce_csv_output && _tofClosest[0] > 6. ) fillCsvForKonrad( pfo, _pdg, _tofClosest[0], _trackLength_IKF_zedLambda);
        }

        // drawDisplay(this, evt, displayPFO, pfo);
        _tree->Fill();
    }
    _treeEvents->Fill();


}

void BohdanAna::end(){
    _file->Write();
    if (_produce_csv_output) _csv_output_file.close();
    // _application.Run(true);
}

void BohdanAna::initialiseTTree(){
    _file.reset( new TFile("results.root", "RECREATE") );
    _tree.reset( new TTree("treename", "treename") );
    _treeEvents.reset( new TTree("treeEvents", "treeEvents") );

    // Event information
    _tree->Branch("crossSection", &_crossSection);
    _tree->Branch("quarksToPythia", &_quarksToPythia);
    _tree->Branch("isHiggsProcess", &_isHiggsProcess);
    _tree->Branch("higgsDaughters", &_higgsDaughters);
    _tree->Branch("ipTrueX", &(_ipTrue[0]) );
    _tree->Branch("ipTrueY", &(_ipTrue[1]) );
    _tree->Branch("ipTrueZ", &(_ipTrue[2]) );
    _tree->Branch("primVertexTrueX", &(_primVertexTrue[0]) );
    _tree->Branch("primVertexTrueY", &(_primVertexTrue[1]) );
    _tree->Branch("primVertexTrueZ", &(_primVertexTrue[2]) );
    _tree->Branch("primVertexRecoX", &(_primVertexReco[0]) );
    _tree->Branch("primVertexRecoY", &(_primVertexReco[1]) );
    _tree->Branch("primVertexRecoZ", &(_primVertexReco[2]) );
    if ( _produce_refit_output ){
        _tree->Branch("primRefitVertexRecoX", &(_primRefitVertexReco[0]) );
        _tree->Branch("primRefitVertexRecoY", &(_primRefitVertexReco[1]) );
        _tree->Branch("primRefitVertexRecoZ", &(_primRefitVertexReco[2]) );
    }

    // True infromation of MCParticle
    _tree->Branch("pdg", &_pdg);
    _tree->Branch("imidiateParentPDG", &_imidiateParentPDG);
    _tree->Branch("firstStableParentPDG", &_firstStableParentPDG);
    _tree->Branch("mcVtxX", &(_mcVtx[0]) );
    _tree->Branch("mcVtxY", &(_mcVtx[1]) );
    _tree->Branch("mcVtxZ", &(_mcVtx[2]) );
    _tree->Branch("timeTrue", &_timeTrue);
    _tree->Branch("generatorStatus", &_generatorStatus);
    _tree->Branch("isSimulated", &_isSimulated);
    _tree->Branch("isBackscatter", &_isBackscatter);
    _tree->Branch("vertexIsNotEndpointOfParent", &_vertexIsNotEndpointOfParent);
    _tree->Branch("isDecayedInTracker", &_isDecayedInTracker);
    _tree->Branch("isDecayedInCalorimeter", &_isDecayedInCalorimeter);
    _tree->Branch("hasLeftDetector", &_hasLeftDetector);
    _tree->Branch("isStopped", &_isStopped);
    _tree->Branch("isOverlay", &_isOverlay);
    _tree->Branch("isBottomQuarkDecay", &_isBottomQuarkDecay);
    _tree->Branch("isCharmQuarkDecay", &_isCharmQuarkDecay);
    _tree->Branch("isHadronisationDecay", &_isHadronisationDecay);
    _tree->Branch("isV0DecayTrue", &_isV0DecayTrue);
    _tree->Branch("isReconstructed", &_isReconstructed);
    _tree->Branch("hasTrack", &_hasTrack);
    _tree->Branch("hasShower", &_hasShower);

    // TRUE VERTEX
    _tree->Branch("isInTruePrimaryVertex", &_isInTruePrimaryVertex);
    _tree->Branch("isInTrueSecondaryVertex", &_isInTrueSecondaryVertex);
    _tree->Branch("isInTrueV0Vertex", &_isInTrueV0Vertex);
    _tree->Branch("trueVertexPosX", &(_trueVertexPos[0]) );
    _tree->Branch("trueVertexPosY", &(_trueVertexPos[1]) );
    _tree->Branch("trueVertexPosZ", &(_trueVertexPos[2]) );
    _tree->Branch("nTracksAtTrueVertex", &_nTracksAtTrueVertex);
    _tree->Branch("invMassOfTrueVertex", &_invMassOfTrueVertex);
    _tree->Branch("minEnergyOfTrueVertex", &_minEnergyOfTrueVertex);
    _tree->Branch("oppositeChargeOfTrueVertex", &_oppositeChargeOfTrueVertex);
    _tree->Branch("cosThetaOfTrueVertex", &_cosThetaOfTrueVertex);
    _tree->Branch("cosThetaToIpOfTrueVertex", &_cosThetaToIpOfTrueVertex);

    // RECO VERTEX
    _tree->Branch("isInRecoPrimaryVertex", &_isInRecoPrimaryVertex);
    _tree->Branch("isInRecoSecondaryVertex", &_isInRecoSecondaryVertex);
    _tree->Branch("isV0DecayReco", &_isV0DecayReco);
    _tree->Branch("recoVertexPosX", &(_recoVertexPos[0]) );
    _tree->Branch("recoVertexPosY", &(_recoVertexPos[1]) );
    _tree->Branch("recoVertexPosZ", &(_recoVertexPos[2]) );
    _tree->Branch("recoVertexPosErrX", &(_recoVertexPosErr[0]) );
    _tree->Branch("recoVertexPosErrY", &(_recoVertexPosErr[1]) );
    _tree->Branch("recoVertexPosErrZ", &(_recoVertexPosErr[2]) );
    _tree->Branch("nTracksAtRecoVertex", &_nTracksAtRecoVertex);
    _tree->Branch("invMassOfRecoVertex", &_invMassOfRecoVertex);
    _tree->Branch("minEnergyOfRecoVertex", &_minEnergyOfRecoVertex);
    _tree->Branch("oppositeChargeOfRecoVertex", &_oppositeChargeOfRecoVertex);
    _tree->Branch("cosThetaOfRecoVertex", &_cosThetaOfRecoVertex);
    _tree->Branch("cosThetaToIpOfRecoVertex", &_cosThetaToIpOfRecoVertex);

    // RECO REFITTED VERTEX
    if ( _produce_refit_output ){
        _tree->Branch("isInRecoPrimaryRefitVertex", &_isInRecoPrimaryRefitVertex);
        _tree->Branch("isInRecoSecondaryRefitVertex", &_isInRecoSecondaryRefitVertex);
        _tree->Branch("isV0DecayRefitReco", &_isV0DecayRefitReco);
        _tree->Branch("recoRefitVertexPosX", &(_recoRefitVertexPos[0]) );
        _tree->Branch("recoRefitVertexPosY", &(_recoRefitVertexPos[1]) );
        _tree->Branch("recoRefitVertexPosZ", &(_recoRefitVertexPos[2]) );
        _tree->Branch("recoRefitVertexPosErrX", &(_recoRefitVertexPosErr[0]) );
        _tree->Branch("recoRefitVertexPosErrY", &(_recoRefitVertexPosErr[1]) );
        _tree->Branch("recoRefitVertexPosErrZ", &(_recoRefitVertexPosErr[2]) );
        _tree->Branch("nTracksAtRecoRefitVertex", &_nTracksAtRecoRefitVertex);
        _tree->Branch("invMassOfRecoRefitVertex", &_invMassOfRecoRefitVertex);
        _tree->Branch("minEnergyOfRecoRefitVertex", &_minEnergyOfRecoRefitVertex);
        _tree->Branch("oppositeChargeOfRecoRefitVertex", &_oppositeChargeOfRecoRefitVertex);
        _tree->Branch("cosThetaOfRecoRefitVertex", &_cosThetaOfRecoRefitVertex);
        _tree->Branch("cosThetaToIpOfRecoRefitVertex", &_cosThetaToIpOfRecoRefitVertex);
    }

    // TRACK
    _tree->Branch("dEdx", &_dEdx);

    // TRUE TRACK STATE AT IP
    _tree->Branch("omegaTrue", &_omegaTrue);
    _tree->Branch("tanLambdaTrue", &_tanLambdaTrue);
    _tree->Branch("d0True", &_d0True);
    _tree->Branch("z0True", &_z0True);
    _tree->Branch("phiTrue", &_phiTrue);
    _tree->Branch("mcPx", &(_mcMom[0]) );
    _tree->Branch("mcPy", &(_mcMom[1]) );
    _tree->Branch("mcPz", &(_mcMom[2]) );

    // TRACK STATE AT IP
    _tree->Branch("omegaIP", &_omegaIP);
    _tree->Branch("tanLambdaIP", &_tanLambdaIP);
    _tree->Branch("d0IP", &_d0IP);
    _tree->Branch("z0IP", &_z0IP);
    _tree->Branch("phiIP", &_phiIP);
    _tree->Branch("omegaErrIP", &_omegaErrIP);
    _tree->Branch("tanLambdaErrIP", &_tanLambdaErrIP);
    _tree->Branch("d0ErrIP", &_d0ErrIP);
    _tree->Branch("z0ErrIP", &_z0ErrIP);
    _tree->Branch("phiErrIP", &_phiErrIP);
    _tree->Branch("recoIpPx", &(_recoIpMom[0]) );
    _tree->Branch("recoIpPy", &(_recoIpMom[1]) );
    _tree->Branch("recoIpPz", &(_recoIpMom[2]) );

    //TRACK STATE AT ECAL
    _tree->Branch("omegaECAL", &_omegaECAL);
    _tree->Branch("tanLambdaECAL", &_tanLambdaECAL);
    _tree->Branch("d0ECAL", &_d0ECAL);
    _tree->Branch("z0ECAL", &_z0ECAL);
    _tree->Branch("phiECAL", &_phiECAL);
    _tree->Branch("omegaErrECAL", &_omegaErrECAL);
    _tree->Branch("tanLambdaErrECAL", &_tanLambdaErrECAL);
    _tree->Branch("d0ErrECAL", &_d0ErrECAL);
    _tree->Branch("z0ErrECAL", &_z0ErrECAL);
    _tree->Branch("phiErrECAL", &_phiErrECAL);
    _tree->Branch("recoCaloX", &(_recoCaloPos[0]) );
    _tree->Branch("recoCaloY", &(_recoCaloPos[1]) );
    _tree->Branch("recoCaloZ", &(_recoCaloPos[2]) );
    _tree->Branch("recoCaloPx", &(_recoCaloMom[0]) );
    _tree->Branch("recoCaloPy", &(_recoCaloMom[1]) );
    _tree->Branch("recoCaloPz", &(_recoCaloMom[2]) );

    if ( _produce_refit_output ){
        //TRACK STATE AT IP REFITTED
        _tree->Branch("refittedOmegaIP", &_refittedOmegaIP);
        _tree->Branch("refittedTanLambdaIP", &_refittedTanLambdaIP);
        _tree->Branch("refittedD0IP", &_refittedD0IP);
        _tree->Branch("refittedZ0IP", &_refittedZ0IP);
        _tree->Branch("refittedPhiIP", &_refittedPhiIP);
        _tree->Branch("refittedOmegaErrIP", &_refittedOmegaErrIP);
        _tree->Branch("refittedTanLambdaErrIP", &_refittedTanLambdaErrIP);
        _tree->Branch("refittedD0ErrIP", &_refittedD0ErrIP);
        _tree->Branch("refittedZ0ErrIP", &_refittedZ0ErrIP);
        _tree->Branch("refittedPhiErrIP", &_refittedPhiErrIP);
        _tree->Branch("refittedRecoIpPx", &(_refittedRecoIpMom[0]) );
        _tree->Branch("refittedRecoIpPy", &(_refittedRecoIpMom[1]) );
        _tree->Branch("refittedRecoIpPz", &(_refittedRecoIpMom[2]) );

        //TRACK STATE AT ECAL REFITTED
        _tree->Branch("refittedOmegaECAL", &_refittedOmegaECAL);
        _tree->Branch("refittedTanLambdaECAL", &_refittedTanLambdaECAL);
        _tree->Branch("refittedD0ECAL", &_refittedD0ECAL);
        _tree->Branch("refittedZ0ECAL", &_refittedZ0ECAL);
        _tree->Branch("refittedPhiECAL", &_refittedPhiECAL);
        _tree->Branch("refittedOmegaErrECAL", &_refittedOmegaErrECAL);
        _tree->Branch("refittedTanLambdaErrECAL", &_refittedTanLambdaErrECAL);
        _tree->Branch("refittedD0ErrECAL", &_refittedD0ErrECAL);
        _tree->Branch("refittedZ0ErrECAL", &_refittedZ0ErrECAL);
        _tree->Branch("refittedPhiErrECAL", &_refittedPhiErrECAL);
        _tree->Branch("refittedRecoCaloX", &(_refittedRecoCaloPos[0]) );
        _tree->Branch("refittedRecoCaloY", &(_refittedRecoCaloPos[1]) );
        _tree->Branch("refittedRecoCaloZ", &(_refittedRecoCaloPos[2]) );
        _tree->Branch("refittedRecoCaloPx", &(_refittedRecoCaloMom[0]) );
        _tree->Branch("refittedRecoCaloPy", &(_refittedRecoCaloMom[1]) );
        _tree->Branch("refittedRecoCaloPz", &(_refittedRecoCaloMom[2]) );
    }

    // TRACK LENGTH and HARMONIC MEAN MOMENTUM ESTIMATORS
    if (not _dst_mode){
        _tree->Branch("trackLength_IDR", &_trackLength_IDR);
        _tree->Branch("trackLengthToEcal_SHA_phiLambda_IP", &_trackLength_SHA_phiLambda_IP);
        _tree->Branch("trackLengthToEcal_SHA_phiZed_IP", &_trackLength_SHA_phiZed_IP);
        _tree->Branch("trackLengthToEcal_SHA_zedLambda_IP", &_trackLength_SHA_zedLambda_IP);
        _tree->Branch("trackLengthToEcal_SHA_phiLambda_ECAL", &_trackLength_SHA_phiLambda_ECAL);
        _tree->Branch("trackLengthToEcal_SHA_phiZed_ECAL", &_trackLength_SHA_phiZed_ECAL);
        _tree->Branch("trackLengthToEcal_SHA_zedLambda_ECAL", &_trackLength_SHA_zedLambda_ECAL);

        _tree->Branch("trackLengthToEcal_IKF_phiLambda", &_trackLength_IKF_phiLambda);
        _tree->Branch("trackLengthToEcal_IKF_phiZed", &_trackLength_IKF_phiZed);
        _tree->Branch("trackLengthToEcal_IKF_zedLambda", &_trackLength_IKF_zedLambda);
        _tree->Branch("harmonicMomToEcal_IKF_phiLambda", &_harmonicMom_IKF_phiLambda);
        _tree->Branch("harmonicMomToEcal_IKF_phiZed", &_harmonicMom_IKF_phiZed);
        _tree->Branch("harmonicMomToEcal_IKF_zedLambda", &_harmonicMom_IKF_zedLambda);

        _tree->Branch("trackLengthToSET_IKF_phiLambda", &_trackLengthToSET_IKF_phiLambda);
        _tree->Branch("trackLengthToSET_IKF_phiZed", &_trackLengthToSET_IKF_phiZed);
        _tree->Branch("trackLengthToSET_IKF_zedLambda", &_trackLengthToSET_IKF_zedLambda);
        _tree->Branch("harmonicMomToSET_IKF_phiLambda", &_harmonicMomToSET_IKF_phiLambda);
        _tree->Branch("harmonicMomToSET_IKF_phiZed", &_harmonicMomToSET_IKF_phiZed);
        _tree->Branch("harmonicMomToSET_IKF_zedLambda", &_harmonicMomToSET_IKF_zedLambda);

        _tree->Branch("cleanTrack", &_cleanTrack);
        //tofs

        _tree->Branch("typeClosest", &_typeClosest);
        _tree->Branch("caloIDClosest", &_caloIDClosest);
        _tree->Branch("layoutClosest", &_layoutClosest);
        _tree->Branch("layerClosest", &_layerClosest);
        _tree->Branch("cleanClosestHit", &_cleanClosestHit);
        for (size_t i = 0; i < _resolutions.size(); i++){
            int res = int(_resolutions[i]);
            _tree->Branch(( "tofClosest"+std::to_string(res) ).c_str(), &( _tofClosest[i]) );
            _tree->Branch(( "tofAverage"+std::to_string(res) ).c_str(), &( _tofAverage[i]) );
            _tree->Branch(( "tofSETFront"+std::to_string(res) ).c_str(), &( _tofSETFront[i]) );
            _tree->Branch(( "tofSETBack"+std::to_string(res) ).c_str(), &( _tofSETBack[i]) );
            _tree->Branch(( "tofFit"+std::to_string(res) ).c_str(), &( _tofFit[i]) );
        }

        _tree->Branch("nHits", &_nHits);
        _tree->Branch("xHit", &_xHit);
        _tree->Branch("yHit", &_yHit);
        _tree->Branch("zHit", &_zHit);
        _tree->Branch("tHit", &_tHit);
        _tree->Branch("layerHit", &_layerHit);
        _tree->Branch("energyHit", &_energyHit);
    }

    _treeEvents->Branch("quarksToPythia", &_quarksToPythia);
    _treeEvents->Branch("isHiggsProcess", &_isHiggsProcess);
    _treeEvents->Branch("higgsDaughters", &_higgsDaughters);
    _treeEvents->Branch("nPionsGenerated", &_nPionsGenerated);
    _treeEvents->Branch("nKaonsGenerated", &_nKaonsGenerated);
    _treeEvents->Branch("nProtonsGenerated", &_nProtonsGenerated);
    _treeEvents->Branch("nPionsTrack", &_nPionsTrack);
    _treeEvents->Branch("nKaonsTrack", &_nKaonsTrack);
    _treeEvents->Branch("nProtonsTrack", &_nProtonsTrack);
    _treeEvents->Branch("nPionsShower", &_nPionsShower);
    _treeEvents->Branch("nKaonsShower", &_nKaonsShower);
    _treeEvents->Branch("nProtonsShower", &_nProtonsShower);


}

void BohdanAna::resetVariables(){
    // True infromation of MCParticle
    _pdg = 0;
    _imidiateParentPDG = 0;
    _firstStableParentPDG = 0;
    _mcVtx = Vector3D();
    _timeTrue = 0.f;
    _generatorStatus = -1;
    _isSimulated = false;
    _isBackscatter = false;
    _vertexIsNotEndpointOfParent = false;
    _isDecayedInTracker = false;
    _isDecayedInCalorimeter = false;
    _hasLeftDetector = false;
    _isStopped = false;
    _isOverlay = false;
    _isBottomQuarkDecay = false;
    _isCharmQuarkDecay = false;
    _isHadronisationDecay = false;
    _isV0DecayTrue = false;
    _isReconstructed = false;
    _hasTrack = false;
    _hasShower = false;

    // TRUE VERTEX
    _isInTruePrimaryVertex = false;
    _isInTrueSecondaryVertex = false;
    _isInTrueV0Vertex = false;
    _trueVertexPos = Vector3D();
    _nTracksAtTrueVertex = 0;
    _invMassOfTrueVertex = 0.f;
    _minEnergyOfTrueVertex = 0.f;
    _oppositeChargeOfTrueVertex = false;
    _cosThetaOfTrueVertex = 0.f;
    _cosThetaToIpOfTrueVertex = 0.f;

    // RECO VERTEX
    _isInRecoPrimaryVertex = false;
    _isInRecoSecondaryVertex = false;
    _isV0DecayReco = false;
    _recoVertexPos = Vector3D();
    _recoVertexPosErr = Vector3D();
    _nTracksAtRecoVertex = 0;
    _invMassOfRecoVertex = 0.f;
    _minEnergyOfRecoVertex = 0.f;
    _oppositeChargeOfRecoVertex = false;
    _chi2OfRecoVertex = 0.f;
    _cosThetaOfRecoVertex = 0.f;
    _cosThetaToIpOfRecoVertex = 0.f;

    // RECO REFITTED VERTEX
    _isInRecoPrimaryRefitVertex = false;
    _isInRecoSecondaryRefitVertex = false;
    _isV0DecayRefitReco = false;
    _recoRefitVertexPos = Vector3D();
    _recoRefitVertexPosErr = Vector3D();
    _nTracksAtRecoRefitVertex = 0;
    _invMassOfRecoRefitVertex = 0.f;
    _minEnergyOfRecoRefitVertex = 0.f;
    _oppositeChargeOfRecoRefitVertex = false;
    _chi2OfRecoRefitVertex = 0.f;
    _cosThetaOfRecoRefitVertex = 0.f;
    _cosThetaToIpOfRecoRefitVertex = 0.f;

    //TRACK
    _dEdx = 0.f;

    // TRUE TRACK STATE AT IP
    _omegaTrue = 0.f;
    _tanLambdaTrue = 0.f;
    _d0True = 0.f;
    _z0True = 0.f;
    _phiTrue = 0.f;
    _mcMom = Vector3D();

    // TRACK STATE AT IP
    _omegaIP = 0.f;
    _tanLambdaIP = 0.f;
    _d0IP = 0.f;
    _z0IP = 0.f;
    _phiIP = 0.f;
    _omegaErrIP = 0.f;
    _tanLambdaErrIP = 0.f;
    _d0ErrIP = 0.f;
    _z0ErrIP = 0.f;
    _phiErrIP = 0.f;
    // _recoIpPos is (0,0,0) no point of storing it...
    _recoIpMom = Vector3D();

    //TRACK STATE AT ECAL
    _omegaECAL = 0.f;
    _tanLambdaECAL = 0.f;
    _d0ECAL = 0.f;
    _z0ECAL = 0.f;
    _phiECAL = 0.f;
    _omegaErrECAL = 0.f;
    _tanLambdaErrECAL = 0.f;
    _d0ErrECAL = 0.f;
    _z0ErrECAL = 0.f;
    _phiErrECAL = 0.f;
    _recoCaloPos = Vector3D();
    _recoCaloMom = Vector3D();

    //TRACK STATE AT IP REFITTED
    _refittedOmegaIP = 0.f;
    _refittedTanLambdaIP = 0.f;
    _refittedD0IP = 0.f;
    _refittedZ0IP = 0.f;
    _refittedPhiIP = 0.f;
    _refittedOmegaErrIP = 0.f;
    _refittedTanLambdaErrIP = 0.f;
    _refittedD0ErrIP = 0.f;
    _refittedZ0ErrIP = 0.f;
    _refittedPhiErrIP = 0.f;
    _refittedRecoIpMom = Vector3D();

    //TRACK STATE AT ECAL REFITTED
    _refittedOmegaECAL = 0.f;
    _refittedTanLambdaECAL = 0.f;
    _refittedD0ECAL = 0.f;
    _refittedZ0ECAL = 0.f;
    _refittedPhiECAL = 0.f;
    _refittedOmegaErrECAL = 0.f;
    _refittedTanLambdaErrECAL = 0.f;
    _refittedD0ErrECAL = 0.f;
    _refittedZ0ErrECAL = 0.f;
    _refittedPhiErrECAL = 0.f;
    _refittedRecoCaloPos = Vector3D();
    _refittedRecoCaloMom = Vector3D();

    // TRACK LENGTH and HARMONIC MEAN MOMENTUM ESTIMATORS
    _trackLength_IDR = 0.f;
    _trackLength_SHA_phiLambda_IP = 0.f;
    _trackLength_SHA_phiZed_IP = 0.f;
    _trackLength_SHA_zedLambda_IP = 0.f;
    _trackLength_SHA_phiLambda_ECAL = 0.f;
    _trackLength_SHA_phiZed_ECAL = 0.f;
    _trackLength_SHA_zedLambda_ECAL = 0.f;

    _trackLength_IKF_phiLambda = 0.f;
    _trackLength_IKF_phiZed = 0.f;
    _trackLength_IKF_zedLambda = 0.f;
    _harmonicMom_IKF_phiLambda = 0.f;
    _harmonicMom_IKF_phiZed = 0.f;
    _harmonicMom_IKF_zedLambda = 0.f;

    _trackLengthToSET_IKF_phiLambda = 0.f;
    _trackLengthToSET_IKF_phiZed = 0.f;
    _trackLengthToSET_IKF_zedLambda = 0.f;
    _harmonicMomToSET_IKF_phiLambda = 0.f;
    _harmonicMomToSET_IKF_phiZed = 0.f;
    _harmonicMomToSET_IKF_zedLambda = 0.f;

    _cleanTrack = true;

    //TOF RECONSTRUCTED FOR A FEW TOF RESOLUTIONS
    _typeClosest = -1;
    _caloIDClosest = -1;
    _layoutClosest = -1;
    _layerClosest = -1;
    _cleanClosestHit = false;
    _tofClosest.fill(0.f);
    _tofAverage.fill(0.f);
    _tofSETFront.fill(0.f);
    _tofSETBack.fill(0.f);
    _tofFit.fill(0.f);

    // ECAL HITS FOR FURTHER TOF RECONSTRUCTION
    _nHits = 0;
    _xHit.clear();
    _yHit.clear();
    _zHit.clear();
    _tHit.clear();
    _layerHit.clear();
    _energyHit.clear();
}