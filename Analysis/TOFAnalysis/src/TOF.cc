#include "TOF.h"
#include "BohdanUtils.h"
#include "marlinutil/CalorimeterHitType.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/TrackState.h"
#include "IMPL/TrackStateImpl.h"
#include "UTIL/ILDConf.h"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "TGraph.h"
#include "TF1.h"
#include "marlin/VerbosityLevels.h"
#include <limits>

using namespace EVENT;
using dd4hep::rec::Vector3D;
using CLHEP::RandGauss;



std::vector<EVENT::CalorimeterHit*> selectFrankEcalHits( EVENT::Cluster* cluster, const Vector3D& posAtEcal, const Vector3D& momAtEcal, int maxEcalLayer ){
    std::vector<CalorimeterHit*> selectedHits(maxEcalLayer, nullptr);
    std::vector<float> minDistances(maxEcalLayer, std::numeric_limits<float>::max());

    for ( auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        int layer = hitType.layer();
        if ( (!isECALHit) || ( layer >= maxEcalLayer) ) continue;

        Vector3D hitPos( hit->getPosition() );
        float dToLine = (hitPos - posAtEcal).cross(momAtEcal.unit()).r();
        if ( dToLine < minDistances[layer] ){
            minDistances[layer] = dToLine;
            selectedHits[layer] = hit;
        }
    }
    selectedHits.erase( std::remove_if( selectedHits.begin(), selectedHits.end(), [](CalorimeterHit* h) { return h == nullptr; } ), selectedHits.end() );

    return selectedHits;
}


EVENT::CalorimeterHit* getClosestHit( EVENT::Cluster* cluster, const dd4hep::rec::Vector3D& posAtEcal){
    auto hits = cluster->getCalorimeterHits();
    if ( hits.empty() ) return nullptr;
    auto sortByDistanceToTrack = [&posAtEcal] (const EVENT::CalorimeterHit* lhs, const EVENT::CalorimeterHit* rhs){
        return (Vector3D(lhs->getPosition()) - posAtEcal).r() < (Vector3D(rhs->getPosition()) - posAtEcal).r();
    };
    return *std::min_element( hits.begin(), hits.end(), sortByDistanceToTrack);
}

float getHitTof( EVENT::CalorimeterHit* hit, const dd4hep::rec::Vector3D& posAtEcal, float timeResolution){
    if (hit == nullptr) return -1.;
    Vector3D hitPos( hit->getPosition() );
    return RandGauss::shoot(hit->getTime(), timeResolution) - (hitPos - posAtEcal).r()/CLHEP::c_light;
}


EVENT::MCParticle* getHitEarliestMC( EVENT::CalorimeterHit* hit, const UTIL::LCRelationNavigator& navToSimCalorimeterHits ){
    // I merge all Calorimeter hit relation collections in the steering file. ENSURE this happens!
    // Otherwise I need to check every hit relation collection separately, which makes this code x10 longer.
    // In case collection doesn't exist, merging is still happens (I think..) with a warning, which is good.
    if (navToSimCalorimeterHits.getRelatedToObjects(hit).empty()) return nullptr;
    
    // There should be really only one sim hit for calo rec hit, but we still do it thourougly.
    const std::vector<float>& weights = navToSimCalorimeterHits.getRelatedToWeights(hit);
    if ( weights.empty() ) return nullptr;
    int max_i = std::max_element(weights.begin(), weights.end()) - weights.begin();
    EVENT::SimCalorimeterHit* simHit = static_cast<EVENT::SimCalorimeterHit*> (navToSimCalorimeterHits.getRelatedToObjects(hit)[max_i]);

    EVENT::MCParticle* mc = nullptr;
    float earliestTime = std::numeric_limits<float>::max();
    for(int i=0; i < simHit->getNMCContributions(); i++){
        if (simHit->getTimeCont(i) < earliestTime){
            earliestTime = simHit->getTimeCont(i);
            mc = simHit->getParticleCont(i);
        }
    }
    return mc;
}

float getTofFrankAvg( const std::vector<EVENT::CalorimeterHit*>& selectedHits, const Vector3D& posAtEcal, float timeResolution){
    float tof = 0.;
    if ( selectedHits.empty() ) return tof;

    for ( auto hit : selectedHits ){
        Vector3D hitPos( hit->getPosition() );
        float dToTrack = (hitPos - posAtEcal).r();
        tof += RandGauss::shoot(hit->getTime(), timeResolution) - dToTrack/CLHEP::c_light;
    }
    return tof/selectedHits.size();
}


float getTofFrankFit( const std::vector<EVENT::CalorimeterHit*>& selectedHits, const Vector3D& posAtEcal, float timeResolution){
    float tof = 0.;
    if ( selectedHits.empty() ) return tof;
    else if ( selectedHits.size() == 1 ){
        //we can't fit 1 point, but lets return something reasonable
        Vector3D hitPos( selectedHits[0]->getPosition() );
        float dToTrack = (hitPos - posAtEcal).r();
        return RandGauss::shoot(selectedHits[0]->getTime(), timeResolution) - dToTrack/CLHEP::c_light;
    }

    std::vector <float> x, y;
    for ( auto hit : selectedHits ){
        Vector3D hitPos( hit->getPosition() );
        float dToTrack = (hitPos - posAtEcal).r();
        x.push_back(dToTrack);
        float time = RandGauss::shoot(hit->getTime(), timeResolution);
        y.push_back(time);
    }

    TGraph gr(x.size(), x.data(), y.data());
    gr.Fit("pol1", "Q");
    return gr.GetFunction("pol1")->GetParameter(0);
}


std::tuple<float, float> getTofSET(EVENT::Track* track, float timeResolution){
    EVENT::TrackerHit* setHit = getSETHit(track);
    if (setHit == nullptr) return std::tuple(0.f, 0.f);
    auto stripObjects = setHit->getRawHits();

    if ( stripObjects.empty() ) return std::tuple(0.f, 0.f);
    else if (stripObjects.size() == 1){
        streamlog_out(WARNING)<<"Found only one SET strip hit, how is this possible!? Writing TOF from a single strip."<<std::endl;
        auto strip = static_cast<EVENT::TrackerHitPlane*> (stripObjects[0]);
        float time = RandGauss::shoot(strip->getTime(), timeResolution);
        return std::tuple(time, time);
    }
    if (stripObjects.size() > 2) streamlog_out(WARNING)<<"Found more than two SET strip hits, how is this possible!?"<<std::endl;
    auto stripFront = static_cast<EVENT::TrackerHitPlane*> (stripObjects[0]);
    auto stripBack = static_cast<EVENT::TrackerHitPlane*> (stripObjects[1]);
    float timeFront = RandGauss::shoot(stripFront->getTime(), timeResolution);
    float timeBack = RandGauss::shoot(stripBack->getTime(), timeResolution);
    return std::tuple(timeFront, timeBack);
}


float getTofPhotonTrue(EVENT::MCParticle* mc){
    // ignore non photons and those who doesn't have a straight path to the calorimeter.
    if (mc->getPDG() != 22 || mc->isDecayedInTracker()) return 0.;

    Vector3D startPos( mc->getVertex() );
    //find intersection points between photon momentum line and ECAL surface planes
    Vector3D finishPos = getPhotonAtCalorimeter(mc);

    //NOTE: photon might not be created at event_time = 0
    return mc->getTime() + (finishPos - startPos).r()/CLHEP::c_light; // in ns
}
