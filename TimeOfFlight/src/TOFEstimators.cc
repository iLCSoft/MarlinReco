#include "TOFEstimators.h"
#include "TOFUtils.h"

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

TOFEstimators aTOFEstimators ;


TOFEstimators::TOFEstimators() : marlin::Processor("TOFEstimators") {
    _description = "TOFEstimators processor computes momentum harmonic mean, track length \
                    and time-of-flight of the chosen ReconstructedParticle to the \
                    specified end point (SET hit or Ecal surface). To be used for a further particle ID";

    registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                            "ReconstructedParticleCollection",
                            "Name of the ReconstructedParticle collection",
                            _pfoCollectionName,
                            std::string("PandoraPFOs") );

    registerProcessorParameter("ExtrapolateToEcal",
                            "If true, track is extrapolated to the Ecal surface for track length calculation, \
                            time of flight estimated using Ecal hits. If false, track length calculated to the last tracker hit.\
                            Time of flight estimated using SET hit if exists.",
                            _extrapolateToEcal,
                            bool(true));

    registerProcessorParameter("TofMethod",
                            "name of the algorithm which estimates time of flight\
                            to the Ecal surface based on Ecal hits time information.\
                            Available options are: closest, frankAvg, frankFit.\
                            In case of _extrapolateToEcal==false is ignored",
                            _tofMethod,
                            std::string("closest") );

    registerProcessorParameter("TimeResolution",
                            "Time resolution of individual SET strips or Ecal hits in ps",
                            _timeResolution,
                            double(0.));

    registerProcessorParameter("MaxEcalLayer",
                            "Time of flight is calculated using Ecal hits only up to MaxLayer",
                            _maxEcalLayer,
                            int(10) );
}


void TOFEstimators::init(){
    marlin::Global::EVENTSEEDER->registerProcessor(this);

    _outputParNames = {"momentumHM", "trackLength", "timeOfFlight"};
    _bField = MarlinUtil::getBzAtOrigin();
    _tpcOuterR = getTPCOuterR();

    _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
    _trkSystem->init();
}


void TOFEstimators::processEvent(LCEvent * evt){
    RandGauss::setTheSeed( marlin::Global::EVENTSEEDER->getSeed(this) );
    ++_nEvent;
    streamlog_out(MESSAGE)<<"******Event****** "<<_nEvent<<std::endl;


    LCCollection* pfos = evt->getCollection(_pfoCollectionName);
    LCRelationNavigator navigatorSET = LCRelationNavigator( evt->getCollection("SETSpacePointRelations") );

    PIDHandler pidHandler( pfos );
    int algoID = pidHandler.addAlgorithm( name(), _outputParNames );


    for (int i=0; i<pfos->getNumberOfElements(); ++i){
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

        ///////////////////////////////////////////////////////////////
        // This part calculates track length and momentum harmonic mean
        ///////////////////////////////////////////////////////////////
        vector<Track*> subTracks = getSubTracks(track);
        vector<TrackStateImpl> trackStates = getTrackStatesPerHit(subTracks, _trkSystem, _extrapolateToEcal, _bField);

        double trackLength = 0.;
        double harmonicMom = 0.;
        int nTrackStates = trackStates.size();
        for( int j=1; j < nTrackStates; ++j ){
            if (j == nTrackStates - 1){
                double nTurns = getHelixNRevolutions( trackStates[j-1], trackStates[j] );
                if( nTurns > 0.5 ){
                    // we cannot calculate arc length for more than pi revolution. Use formula with only z
                    double arcLength = getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );
                    Vector3D mom = getHelixMomAtTrackState( trackStates[j-1], _bField );
                    trackLength += arcLength;
                    harmonicMom += arcLength/mom.r2();
                    continue;
                }
            }

            double arcLength = getHelixArcLength( trackStates[j-1], trackStates[j] );
            Vector3D mom = getHelixMomAtTrackState( trackStates[j-1], _bField );
            trackLength += arcLength;
            harmonicMom += arcLength/mom.r2();
        }
        harmonicMom = std::sqrt(trackLength/harmonicMom);


        //////////////////////////////////////
        // This part calculates Time of flight
        //////////////////////////////////////
        double timeOfFlight = 0.;
        if( _extrapolateToEcal ){
            // if we extrapolate to the calorimeter check TOF method steering parameter
            if (_tofMethod == "closest"){
                timeOfFlight = getTofClosest(cluster, track, _timeResolution/1000.);
            }
            else if (_tofMethod == "frankAvg"){
                vector<CalorimeterHit*> frankHits = selectFrankEcalHits(cluster, track, _maxEcalLayer, _bField);
                timeOfFlight = getTofFrankAvg(frankHits, track, _timeResolution/1000.);
            }
            else if (_tofMethod == "frankFit"){
                vector<CalorimeterHit*> frankHits = selectFrankEcalHits(cluster, track, _maxEcalLayer, _bField);
                timeOfFlight = getTofFrankFit(frankHits, track, _timeResolution/1000.);
            }
            else{
                //wrong parameter silly!
            }
        }
        else{
            //define tof as average time between two SET strips
            //if no SET hits found, tof=0
            TrackerHit* hitSET = getSETHit(track, _tpcOuterR);
            if ( hitSET != nullptr ){
                const vector<LCObject*>& simHitsSET = navigatorSET.getRelatedToObjects( hitSET );
                if ( simHitsSET.size() == 2 ){
                    //it should be 2 always.. but just to be safe
                    SimTrackerHit* simHitSETFront = static_cast <SimTrackerHit*>( simHitsSET[0] );
                    SimTrackerHit* simHitSETBack = static_cast <SimTrackerHit*>( simHitsSET[1] );
                    double timeFront = RandGauss::shoot( simHitSETFront->getTime(), _timeResolution/1000. );
                    double timeBack = RandGauss::shoot( simHitSETBack->getTime(), _timeResolution/1000. );
                        timeOfFlight = (timeFront + timeBack)/2.;
                }
            }
        }

        vector<float> results{float(harmonicMom), float(trackLength), float(timeOfFlight)};
        pidHandler.setParticleID(pfo , 0, 0, 0., algoID, results);
    }
}
