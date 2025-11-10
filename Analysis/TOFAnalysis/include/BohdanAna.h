#ifndef BohdanAna_h
#define BohdanAna_h 1

#include "marlin/Processor.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "EventDisplayer.h"
#include "TrackLength.h"

#include "BohdanUtils.h"

class BohdanAna : public marlin::Processor, EventDisplayer {
    friend class EventDisplayer;
    public:
        BohdanAna(const BohdanAna&) = delete;
        BohdanAna& operator=(const BohdanAna&) = delete;

        marlin::Processor* newProcessor() { return new BohdanAna; }

        BohdanAna();
        void init();
        void processEvent(EVENT::LCEvent* evt);
        void end();

        void fillMCTrueInfo(EVENT::MCParticle* mc, EVENT::ReconstructedParticle* pfo);

        void fillTrueVertexInfo(EVENT::MCParticle* mc, const std::vector<VertexData>& trueVertices);

        void fillRecoVertexInfo(EVENT::LCEvent* evt, EVENT::MCParticle* mc, const UTIL::LCRelationNavigator& pfo2mc);

        void fillTrackStates(EVENT::ReconstructedParticle* pfo, EVENT::ReconstructedParticle* refittedPfo);

        void fillTrackLengthInfo(EVENT::MCParticle* mc, EVENT::ReconstructedParticle* pfo, const UTIL::LCRelationNavigator& navToSimTrackerHits);

        void fillTOFInfo(EVENT::MCParticle* mc, EVENT::ReconstructedParticle* pfo, const UTIL::LCRelationNavigator& navToSimCalorimeterHits);

        void fillCsvForKonrad(EVENT::ReconstructedParticle* pfo, int pdg, double trueTOF, double trackLength);

        void initialiseTTree();
        void resetVariables();
    private:
        int _nEvent{};
        float _bField{};

        std::unique_ptr<TFile> _file{};
        std::unique_ptr<TTree> _tree{};
        std::unique_ptr<TTree> _treeEvents{};

        bool _dst_mode{};
        bool _produce_refit_output{};
        bool _produce_csv_output{};
        std::ofstream _csv_output_file{};
        int _global_pfo_number{};

        // Event information
        float _crossSection{};
        int _nPionsGenerated{};
        int _nKaonsGenerated{};
        int _nProtonsGenerated{};
        int _nPionsTrack{};
        int _nKaonsTrack{};
        int _nProtonsTrack{};
        int _nPionsShower{};
        int _nKaonsShower{};
        int _nProtonsShower{};
        int _quarksToPythia{};
        bool _isHiggsProcess{};
        int _higgsDaughters{};
        dd4hep::rec::Vector3D _ipTrue{};
        dd4hep::rec::Vector3D _primVertexTrue{};
        dd4hep::rec::Vector3D _primVertexReco{};
        dd4hep::rec::Vector3D _primRefitVertexReco{};

        // True infromation of MCParticle
        int _pdg{};
        int _imidiateParentPDG{};
        int _firstStableParentPDG{};
        dd4hep::rec::Vector3D _mcVtx{};
        float _timeTrue{};
        int _generatorStatus{};
        bool _isSimulated{};
        bool _isBackscatter{};
        bool _vertexIsNotEndpointOfParent{};
        bool _isDecayedInTracker{};
        bool _isDecayedInCalorimeter{};
        bool _hasLeftDetector{};
        bool _isStopped{};
        bool _isOverlay{};
        bool _isBottomQuarkDecay{};
        bool _isCharmQuarkDecay{};
        bool _isHadronisationDecay{};
        bool _isV0DecayTrue{};
        bool _isReconstructed{};
        bool _hasTrack{};
        bool _hasShower{};
        // TRUE TRACK STATE AT IP
        float _omegaTrue{};
        float _tanLambdaTrue{};
        float _d0True{};
        float _z0True{};
        float _phiTrue{};
        dd4hep::rec::Vector3D _mcMom{};

        // TRUE VERTEX (NOTE: TRUE here means TRUE ON MC LEVEL BUT RECONSTRUCTABLE IN PRINCIPLE. E.g. has two reco tracks!)
        bool _isInTruePrimaryVertex{};
        bool _isInTrueSecondaryVertex{};
        bool _isInTrueV0Vertex{};
        dd4hep::rec::Vector3D _trueVertexPos{};
        unsigned int _nTracksAtTrueVertex{};
        float _invMassOfTrueVertex{};
        float _minEnergyOfTrueVertex{};
        bool _oppositeChargeOfTrueVertex{};
        float _cosThetaOfTrueVertex{};
        float _cosThetaToIpOfTrueVertex{};

        // RECO VERTEX
        bool _isInRecoPrimaryVertex{};
        bool _isInRecoSecondaryVertex{};
        bool _isV0DecayReco{};
        dd4hep::rec::Vector3D _recoVertexPos{};
        dd4hep::rec::Vector3D _recoVertexPosErr{};
        unsigned int _nTracksAtRecoVertex{};
        float _invMassOfRecoVertex{};
        float _minEnergyOfRecoVertex{};
        bool _oppositeChargeOfRecoVertex{};
        float _chi2OfRecoVertex{};
        float _cosThetaOfRecoVertex{};
        float _cosThetaToIpOfRecoVertex{};

        // RECO REFITTED VERTEX
        bool _isInRecoPrimaryRefitVertex{};
        bool _isInRecoSecondaryRefitVertex{};
        bool _isV0DecayRefitReco{};
        dd4hep::rec::Vector3D _recoRefitVertexPos{};
        dd4hep::rec::Vector3D _recoRefitVertexPosErr{};
        unsigned int _nTracksAtRecoRefitVertex{};
        float _invMassOfRecoRefitVertex{};
        float _minEnergyOfRecoRefitVertex{};
        bool _oppositeChargeOfRecoRefitVertex{};
        float _chi2OfRecoRefitVertex{};
        float _cosThetaOfRecoRefitVertex{};
        float _cosThetaToIpOfRecoRefitVertex{};

        // TRACK.
        float _dEdx{};

        // TRACK STATE AT IP
        float _omegaIP{};
        float _tanLambdaIP{};
        float _d0IP{};
        float _z0IP{};
        float _phiIP{};
        float _omegaErrIP{};
        float _tanLambdaErrIP{};
        float _d0ErrIP{};
        float _z0ErrIP{};
        float _phiErrIP{};
        // _recoIpPos is (0,0,0) no point of storing it...
        dd4hep::rec::Vector3D _recoIpMom{};

        //TRACK STATE AT ECAL
        float _omegaECAL{};
        float _tanLambdaECAL{};
        float _d0ECAL{};
        float _z0ECAL{};
        float _phiECAL{};
        float _omegaErrECAL{};
        float _tanLambdaErrECAL{};
        float _d0ErrECAL{};
        float _z0ErrECAL{};
        float _phiErrECAL{};
        dd4hep::rec::Vector3D _recoCaloPos{};
        dd4hep::rec::Vector3D _recoCaloMom{};

        //TRACK STATE AT IP REFITTED
        float _refittedOmegaIP{};
        float _refittedTanLambdaIP{};
        float _refittedD0IP{};
        float _refittedZ0IP{};
        float _refittedPhiIP{};
        float _refittedOmegaErrIP{};
        float _refittedTanLambdaErrIP{};
        float _refittedD0ErrIP{};
        float _refittedZ0ErrIP{};
        float _refittedPhiErrIP{};
        dd4hep::rec::Vector3D _refittedRecoIpMom{};

        //TRACK STATE AT ECAL REFITTED
        float _refittedOmegaECAL{};
        float _refittedTanLambdaECAL{};
        float _refittedD0ECAL{};
        float _refittedZ0ECAL{};
        float _refittedPhiECAL{};
        float _refittedOmegaErrECAL{};
        float _refittedTanLambdaErrECAL{};
        float _refittedD0ErrECAL{};
        float _refittedZ0ErrECAL{};
        float _refittedPhiErrECAL{};
        dd4hep::rec::Vector3D _refittedRecoCaloPos{};
        dd4hep::rec::Vector3D _refittedRecoCaloMom{};

        // TRACK LENGTH and HARMONIC MEAN MOMENTUM ESTIMATORS
        float _trackLength_IDR{};
        float _trackLength_SHA_phiLambda_IP{};
        float _trackLength_SHA_phiZed_IP{};
        float _trackLength_SHA_zedLambda_IP{};
        float _trackLength_SHA_phiLambda_ECAL{};
        float _trackLength_SHA_phiZed_ECAL{};
        float _trackLength_SHA_zedLambda_ECAL{};

        float _trackLength_IKF_phiLambda{};
        float _trackLength_IKF_phiZed{};
        float _trackLength_IKF_zedLambda{};
        float _harmonicMom_IKF_phiLambda{};
        float _harmonicMom_IKF_phiZed{};
        float _harmonicMom_IKF_zedLambda{};

        float _trackLengthToSET_IKF_phiLambda{};
        float _trackLengthToSET_IKF_phiZed{};
        float _trackLengthToSET_IKF_zedLambda{};
        float _harmonicMomToSET_IKF_phiLambda{};
        float _harmonicMomToSET_IKF_phiZed{};
        float _harmonicMomToSET_IKF_zedLambda{};

        bool _cleanTrack = true;

        //TOF RECONSTRUCTED FOR A FEW TOF RESOLUTIONS
        int _typeClosest = -1;
        int _caloIDClosest = -1;
        int _layoutClosest = -1;
        int _layerClosest = -1;
        bool _cleanClosestHit = false;
        std::vector<float> _resolutions = {0, 1, 5, 10, 17, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300};
        std::array<float, 16> _tofClosest{};
        std::array<float, 16> _tofAverage{};
        std::array<float, 16> _tofSETFront{};
        std::array<float, 16> _tofSETBack{};
        std::array<float, 16> _tofFit{};

        // ECAL HITS FOR FURTHER TOF RECONSTRUCTION
        int _nHits{};
        std::vector<float> _xHit{};
        std::vector<float> _yHit{};
        std::vector<float> _zHit{};
        std::vector<float> _tHit{};
        std::vector<int> _layerHit{};
        std::vector<float> _energyHit{};

        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;
        TApplication _application = TApplication("app", 0, nullptr);
};

#endif
