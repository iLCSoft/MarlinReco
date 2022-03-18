#include "TrackLengthProcessor.h"
#include "TOFUtils.h"

#include <chrono>

#include "EVENT/LCCollection.h"
#include "UTIL/PIDHandler.h"

#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"
#include "marlin/VerbosityLevels.h"
#include "marlinutil/GeometryUtil.h"
#include "MarlinTrk/Factory.h"
#include "EVENT/SimTrackerHit.h"
#include "UTIL/LCRelationNavigator.h"
#include "CLHEP/Random/Randomize.h"

using namespace TOFUtils;
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::ReconstructedParticle;
using EVENT::TrackerHit;
using EVENT::Track;
using EVENT::SimTrackerHit;
using EVENT::Cluster;
using EVENT::CalorimeterHit;
using EVENT::TrackState;
using EVENT::LCObject;
using UTIL::LCRelationNavigator;
using CLHEP::RandGauss;
using dd4hep::rec::Vector3D;

TrackLengthProcessor aTrackLengthProcessor ;


TrackLengthProcessor::TrackLengthProcessor() : marlin::Processor("TrackLengthProcessor") {
    _description = "TrackLengthProcessor computes track length and harmonic mean of the momentum square up to the SET hit and Ecal surface.\
                    It should give a relatively good approximation of the track length of any particle reaching the Ecal to use in the time-of-flight particle identification";

    registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                            "ReconstructedParticleCollection",
                            "Name of the ReconstructedParticle collection",
                            _pfoCollectionName,
                            std::string("PandoraPFOs") );
}


void TrackLengthProcessor::init(){
    marlin::Global::EVENTSEEDER->registerProcessor(this);

    _outputParNames = {"trackLengthToEcal", "trackLengthToSET", "momentumHMToEcal", "momentumHMToSET"};
    _bField = MarlinUtil::getBzAtOrigin();

    _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
    _trkSystem->init();
}


void TrackLengthProcessor::processEvent(EVENT::LCEvent * evt){
    RandGauss::setTheSeed( marlin::Global::EVENTSEEDER->getSeed(this) );
    ++_nEvent;
    streamlog_out(DEBUG9)<<std::endl<<"==========Event========== "<<_nEvent<<std::endl;

    LCCollection* pfos = evt->getCollection(_pfoCollectionName);

    PIDHandler pidHandler( pfos );
    int algoID = pidHandler.addAlgorithm( name(), _outputParNames );


    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        streamlog_out(DEBUG9)<<std::endl<<"Starting to analyze "<<i+1<<" PFO"<<std::endl;
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );

        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();

        if( nClusters != 1 || nTracks != 1){
            // Analyze only simple pfos. Otherwise write dummy zeros
            vector<float> results{0., 0., 0.};
            pidHandler.setParticleID(pfo , 0, 0, 0., algoID, results);
            continue;
        }
        Track* track = pfo->getTracks()[0];
        Cluster* cluster = pfo->getClusters()[0];

        vector<Track*> subTracks = getSubTracks(track);
        bool extrapolateToEcal = false;
        vector<TrackStateImpl> trackStates = getTrackStatesPerHit(subTracks, _trkSystem, _bField);

        double trackLengthToSET = 0.;
        double harmonicMomToSET = 0.;
        double trackLengthToEcal = 0.;
        double harmonicMomToEcal = 0.;
        int nTrackStates = trackStates.size();
        //exclude last track state at the ECal
        for( int j=1; j < nTrackStates - 1; ++j ){
            //we check which track length formula to use
            double nTurns = getHelixNRevolutions( trackStates[j-1], trackStates[j] );
            double arcLength;
            // we cannot calculate arc length for more than pi revolution using delta phi. Use formula with only z
            if ( nTurns <= 0.5 ) arcLength = getHelixArcLength( trackStates[j-1], trackStates[j] );
            else arcLength = getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );

            Vector3D mom = getHelixMomAtTrackState( trackStates[j-1], _bField );
            trackLengthToSET += arcLength;
            harmonicMomToSET += arcLength/mom.r2();
        }
        harmonicMomToSET = std::sqrt(trackLengthToSET/harmonicMomToSET);
        
        trackLengthToEcal = trackLengthToSET;
        harmonicMomToEcal = harmonicMomToSET;

        //calculate 1 more step from SET to the Ecal
        double nTurns = getHelixNRevolutions( trackStates[nTrackStates - 2], trackStates[nTrackStates - 1] );
        double arcLength;
        if ( nTurns <= 0.5 ) arcLength = getHelixArcLength( trackStates[nTrackStates - 2], trackStates[nTrackStates - 1] );
        else arcLength = getHelixLengthAlongZ( trackStates[nTrackStates - 2], trackStates[nTrackStates - 1] );
        Vector3D mom = getHelixMomAtTrackState( trackStates[nTrackStates - 2], _bField );
        trackLengthToEcal += arcLength;
        harmonicMomToEcal += arcLength/mom.r2();


        vector<float> results{float(trackLengthToSET), float(trackLengthToEcal), float(harmonicMomToSET), float(harmonicMomToEcal)};
        pidHandler.setParticleID(pfo , 0, 0, 0., algoID, results);
        streamlog_out(DEBUG9)<<"Final results for the "<<i+1<<" PFO"<<std::endl;
        streamlog_out(DEBUG9)<<"Track length to the SET: "<< float(trackLengthToSET)<<" mm"<<std::endl;
        streamlog_out(DEBUG9)<<"Track length to the ECal: "<< float(trackLengthToEcal)<<" mm"<<std::endl;
        streamlog_out(DEBUG9)<<"Harmonic mean momentum to the SET: "<< float(harmonicMomToSET)<<" Gev"<<std::endl;
        streamlog_out(DEBUG9)<<"Harmonic mean momentum to the Ecal: "<< float(harmonicMomToEcal)<<" Gev"<<std::endl;
        streamlog_out(DEBUG9)<<std::endl<<std::endl;
    }
}
