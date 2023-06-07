#include "TrackLengthProcessor.h"
#include "TrackLengthUtils.h"

#include "EVENT/LCCollection.h"
#include "UTIL/PIDHandler.h"

#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"
#include "marlin/VerbosityLevels.h"
#include "marlinutil/GeometryUtil.h"
#include "MarlinTrk/Factory.h"
#include "EVENT/SimTrackerHit.h"
#include "UTIL/TrackTools.h"

using namespace TrackLengthUtils;
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::ReconstructedParticle;
using EVENT::TrackerHit;
using EVENT::Track;
using EVENT::SimTrackerHit;
using EVENT::TrackState;
using EVENT::LCObject;
using dd4hep::rec::Vector3D;

TrackLengthProcessor aTrackLengthProcessor ;


TrackLengthProcessor::TrackLengthProcessor() : marlin::Processor("TrackLengthProcessor") {
    _description = "TrackLengthProcessor computes track length and harmonic mean of the momentum square up to the last tracker hit and Ecal surface.\
                    It should give a relatively good approximation of the track length of any particle reaching the Ecal to use in the time-of-flight particle identification";

    registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                            "ReconstructedParticleCollection",
                            "Name of the ReconstructedParticle collection",
                            _pfoCollectionName,
                            std::string("PandoraPFOs") );
}


void TrackLengthProcessor::init(){
    marlin::Global::EVENTSEEDER->registerProcessor(this);

    //NOTE: although the name includes SET, it is actually last tracker hit and not necessarily SET.
    _outputParNames = {"trackLengthToSET", "trackLengthToEcal", "momentumHMToSET", "momentumHMToEcal"};
    _bField = MarlinUtil::getBzAtOrigin();

    _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
    _trkSystem->init();
}


void TrackLengthProcessor::processEvent(EVENT::LCEvent * evt){
    ++_nEvent;
    streamlog_out(DEBUG9)<<std::endl<<"==========Event========== "<<_nEvent<<std::endl;

    LCCollection* pfos = evt->getCollection(_pfoCollectionName);

    PIDHandler pidHandler( pfos );
    int algoID = pidHandler.addAlgorithm( name(), _outputParNames );


    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        streamlog_out(DEBUG9)<<std::endl<<"Starting to analyze "<<i+1<<" PFO"<<std::endl;
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );

        if( pfo->getTracks().empty() ){
            // Calculate track length only for particles with at least one track.
            // Note: If the pfo has more than one track attached (e.g. kink, v0), only the first track will be used. The results may be not accurate in these cases. 
            streamlog_out(DEBUG9)<<"PFO has no tracks. Writing zeros."<<std::endl;
            vector<float> results{0., 0., 0., 0.};
            pidHandler.setParticleID(pfo , 0, 0, 0., algoID, results);
            continue;
        }
        vector<TrackStateImpl> trackStates = getTrackStates(pfo, _trkSystem, _bField);

        double trackLengthToSET = 0.;
        double harmonicMomToSET = 0.;
        double trackLengthToEcal = 0.;
        double harmonicMomToEcal = 0.;
        int nTrackStates = trackStates.size();
        streamlog_out(DEBUG9)<<"PFO has "<<nTrackStates<<" track states to calculate the length"<<std::endl;
        if (nTrackStates <= 1){
            streamlog_out(DEBUG9)<<" Not enough track states to calculate the track length. Writing zeros."<<std::endl;
            vector<float> results{0., 0., 0., 0.};
            pidHandler.setParticleID(pfo , 0, 0, 0., algoID, results);
            continue;
        }

        //exclude last track state at the ECal
        for( int j=1; j < nTrackStates-1; ++j ){
            double arcLength = getHelixLength( trackStates[j-1], trackStates[j] );

            std::array<double, 3> momArr = UTIL::getTrackMomentum( &(trackStates[j-1]), _bField);
            Vector3D mom(momArr[0], momArr[1], momArr[2]);
            trackLengthToSET += arcLength;
            harmonicMomToSET += arcLength/mom.r2();
        }
        
        //now calculate to the Ecal one more step
        double arcLength = getHelixLength( trackStates[nTrackStates - 2], trackStates[nTrackStates - 1] );
        std::array<double, 3> momArr = UTIL::getTrackMomentum( &(trackStates[nTrackStates - 2]), _bField );
        Vector3D mom(momArr[0], momArr[1], momArr[2]);
        trackLengthToEcal = trackLengthToSET + arcLength;
        harmonicMomToEcal = harmonicMomToSET + arcLength/mom.r2();

        harmonicMomToSET = std::sqrt(trackLengthToSET/harmonicMomToSET);
        harmonicMomToEcal = std::sqrt(trackLengthToEcal/harmonicMomToEcal);

        vector<float> results{float(trackLengthToSET), float(trackLengthToEcal), float(harmonicMomToSET), float(harmonicMomToEcal)};
        pidHandler.setParticleID(pfo , 0, 0, 0., algoID, results);
        streamlog_out(DEBUG9)<<"Final results for the "<<i+1<<" PFO"<<std::endl;
        streamlog_out(DEBUG9)<<"Track length to the last tracker hit: "<< float(trackLengthToSET)<<" mm"<<std::endl;
        streamlog_out(DEBUG9)<<"Track length to the ECal: "<< float(trackLengthToEcal)<<" mm"<<std::endl;
        streamlog_out(DEBUG9)<<"Harmonic mean momentum to the last tracker hit: "<< float(harmonicMomToSET)<<" GeV"<<std::endl;
        streamlog_out(DEBUG9)<<"Harmonic mean momentum to the Ecal: "<< float(harmonicMomToEcal)<<" GeV"<<std::endl;
        streamlog_out(DEBUG9)<<std::endl<<std::endl;
    }
}
