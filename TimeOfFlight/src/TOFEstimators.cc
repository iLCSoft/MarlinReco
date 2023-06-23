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
    _description = "TOFEstimators processor computes time-of-flight of the chosen ReconstructedParticle to the \
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

    if(_tofMethod != "closest" && _tofMethod != "frankAvg" && _tofMethod != "frankFit"){
        throw EVENT::Exception( "Invalid steering parameter for TofMethod is passed: " + _tofMethod + "\n Available options are: closest, frankAvg, frankFit" );
    }

    marlin::Global::EVENTSEEDER->registerProcessor(this);

    _outputParNames = {"timeOfFlight"};
    _bField = MarlinUtil::getBzAtOrigin();
    // internally we use time resolution in nanoseconds
    _timeResolution = _timeResolution/1000.;
}


void TOFEstimators::processEvent(EVENT::LCEvent * evt){
    RandGauss::setTheSeed( marlin::Global::EVENTSEEDER->getSeed(this) );
    ++_nEvent;
    streamlog_out(DEBUG9)<<std::endl<<"==========Event========== "<<_nEvent<<std::endl;

    LCCollection* pfos = evt->getCollection(_pfoCollectionName);
    LCCollection* setRelations = evt->getCollection("SETSpacePointRelations");
    
    LCRelationNavigator navigatorSET = LCRelationNavigator( setRelations );

    PIDHandler pidHandler( pfos );
    int algoID = pidHandler.addAlgorithm( name(), _outputParNames );


    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        streamlog_out(DEBUG9)<<std::endl<<"Starting to analyze "<<i+1<<" PFO"<<std::endl;
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );

        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();

        if( nClusters != 1 || nTracks != 1){
            // Analyze only simple pfos. Otherwise write dummy zeros
            vector<float> results{0.};
            pidHandler.setParticleID(pfo , 0, 0, 0., algoID, results);
            continue;
        }
        Track* track = pfo->getTracks()[0];
        Cluster* cluster = pfo->getClusters()[0];

        double timeOfFlight = 0.;
        if( _extrapolateToEcal ){
            if (_tofMethod == "closest"){
                timeOfFlight = getTofClosest(cluster, track, _timeResolution);
            }
            else if (_tofMethod == "frankAvg"){
                vector<CalorimeterHit*> frankHits = selectFrankEcalHits(cluster, track, _maxEcalLayer, _bField);
                timeOfFlight = getTofFrankAvg(frankHits, track, _timeResolution);
            }
            else if (_tofMethod == "frankFit"){
                vector<CalorimeterHit*> frankHits = selectFrankEcalHits(cluster, track, _maxEcalLayer, _bField);
                timeOfFlight = getTofFrankFit(frankHits, track, _timeResolution);
            }
        }
        else{
            //define tof as an average time between two SET strips
            //if no SET hits found, tof alreasy is 0, just skip
            TrackerHit* hitSET = getSETHit(track);
            if ( hitSET != nullptr ){
                const vector<LCObject*>& simHitsSET = navigatorSET.getRelatedToObjects( hitSET );
                if ( simHitsSET.size() >= 2 ){
                    //It must be always 2, but just in case...
                    if (simHitsSET.size() > 2) streamlog_out(WARNING)<<"Found more than two SET strip hits! Writing TOF as an average of the first two elements in the array."<<std::endl;

                    SimTrackerHit* simHitSETFront = static_cast <SimTrackerHit*>( simHitsSET[0] );
                    SimTrackerHit* simHitSETBack = static_cast <SimTrackerHit*>( simHitsSET[1] );
                    double timeFront = RandGauss::shoot(simHitSETFront->getTime(), _timeResolution);
                    double timeBack = RandGauss::shoot(simHitSETBack->getTime(), _timeResolution);
                        timeOfFlight = (timeFront + timeBack)/2.;
                }
                else if (simHitsSET.size() == 1){
                    streamlog_out(WARNING)<<"Found only one SET strip hit! Writing TOF from a single strip."<<std::endl;
                    SimTrackerHit* simHitSET = static_cast <SimTrackerHit*>(simHitsSET[0]);
                    timeOfFlight = RandGauss::shoot(simHitSET->getTime(), _timeResolution);
                }
                else{
                    // this happens very rarily (0.1%). When >1 simHits associated with a single strip none simHits are written by the DDSpacePointBuilder.
                    streamlog_out(WARNING)<<"Found NO simHits associated with the found SET hit! Writing TOF as 0."<<std::endl;
                }
            }
        }
        vector<float> results{float(timeOfFlight)};
        pidHandler.setParticleID(pfo , 0, 0, 0., algoID, results);
        streamlog_out(DEBUG9)<<"Final results for the "<<i+1<<" PFO"<<std::endl;
        streamlog_out(DEBUG9)<<"time-of-flight: "<< float(timeOfFlight)<<" ns"<<std::endl;
        streamlog_out(DEBUG9)<<std::endl<<std::endl;
    }
}
