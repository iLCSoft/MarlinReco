#include "TOFUtils.h"

#include <cmath>
#include <algorithm>
#include <limits>

#include "marlin/VerbosityLevels.h"
#include "marlinutil/CalorimeterHitType.h"
#include "HelixClass.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "TGraph.h"
#include "TF1.h"

using std::vector;
using std::numeric_limits;
using EVENT::TrackerHit;
using EVENT::Track;
using EVENT::Cluster;
using EVENT::CalorimeterHit;
using EVENT::TrackState;
using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::Vector3D;
using dd4hep::rec::FixedPadSizeTPCData;
using CLHEP::RandGauss;


dd4hep::rec::Vector3D TOFUtils::getHelixMomAtTrackState(const EVENT::TrackState& ts, double bField){
    double phi = ts.getPhi();
    double d0 = ts.getD0();
    double z0 = ts.getZ0();
    double omega = ts.getOmega();
    double tanL = ts.getTanLambda();

    HelixClass helix;
    helix.Initialize_Canonical(phi, d0, z0, omega, tanL, bField);
    return helix.getMomentum();
}


double TOFUtils::getTPCOuterR(){
    auto& detector = Detector::getInstance();
    DetElement tpcDet = detector.detector("TPC");
    FixedPadSizeTPCData* tpc = tpcDet.extension <FixedPadSizeTPCData>();
    return tpc->rMaxReadout/dd4hep::mm;
}


EVENT::TrackerHit* TOFUtils::getSETHit(EVENT::Track* track, double tpcOuterR){
    vector<TrackerHit*> hits = track->getTrackerHits();
    TrackerHit* lastHit = hits.back();
    Vector3D pos (lastHit->getPosition());

    if ( pos.rho() > tpcOuterR ) return lastHit;
    return nullptr;
}


std::vector<EVENT::CalorimeterHit*> TOFUtils::selectFrankEcalHits( EVENT::Cluster* cluster, EVENT::Track* track, int maxEcalLayer, double bField ){
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D trackPosAtEcal ( tsEcal->getReferencePoint() );
    Vector3D trackMomAtEcal = TOFUtils::getHelixMomAtTrackState(*tsEcal, bField);

    vector<CalorimeterHit*> selectedHits(maxEcalLayer, nullptr);
    vector<double> minDistances(maxEcalLayer, numeric_limits<double>::max());

    for ( auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        int layer = hitType.layer();
        if ( (!isECALHit) || ( layer >= maxEcalLayer) ) continue;

        Vector3D hitPos( hit->getPosition() );
        double dToLine = (hitPos - trackPosAtEcal).cross(trackMomAtEcal.unit()).r();
        if ( dToLine < minDistances[layer] ){
            minDistances[layer] = dToLine;
            selectedHits[layer] = hit;
        }
    }
    selectedHits.erase( std::remove_if( selectedHits.begin(), selectedHits.end(), [](CalorimeterHit* h) { return h == nullptr; } ), selectedHits.end() );

    return selectedHits;
}


double TOFUtils::getTofClosest( EVENT::Cluster* cluster, EVENT::Track* track, double timeResolution){
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D trackPosAtEcal ( tsEcal->getReferencePoint() );

    double hitTime = numeric_limits<double>::max();
    double closestDistance = numeric_limits<double>::max();
    for( auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if (! isECALHit) continue;

        Vector3D hitPos( hit->getPosition() );
        double dToTrack = (hitPos - trackPosAtEcal).r();
        if( dToTrack < closestDistance ){
            closestDistance = dToTrack;
            hitTime = hit->getTime();
        }
    }

    if ( hitTime == numeric_limits<double>::max() ) return 0.;
    return RandGauss::shoot(hitTime, timeResolution) - closestDistance/CLHEP::c_light;
}


double TOFUtils::getTofFrankAvg( std::vector<EVENT::CalorimeterHit*> selectedHits, EVENT::Track* track, double timeResolution){
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D trackPosAtEcal ( tsEcal->getReferencePoint() );

    int nHits = selectedHits.size();
    if (nHits == 0) return 0.;

    double tof = 0.;
    for ( auto hit : selectedHits ){
        Vector3D hitPos( hit->getPosition() );
        double dToTrack = (hitPos - trackPosAtEcal).r();
        tof += RandGauss::shoot(hit->getTime(), timeResolution) - dToTrack/CLHEP::c_light;
    }
    return tof/nHits;
}


double TOFUtils::getTofFrankFit( std::vector<EVENT::CalorimeterHit*> selectedHits, EVENT::Track* track, double timeResolution){
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D trackPosAtEcal ( tsEcal->getReferencePoint() );

    int nHits = selectedHits.size();
    if (nHits == 0) return 0.;
    else if (nHits == 1){
        //we can't fit 1 point, but lets return something reasonable
        Vector3D hitPos( selectedHits[0]->getPosition() );
        double dToTrack = (hitPos - trackPosAtEcal).r();
        return RandGauss::shoot(selectedHits[0]->getTime(), timeResolution) - dToTrack/CLHEP::c_light;
    }

    vector <double> x, y;
    for ( auto hit : selectedHits ){
        Vector3D hitPos( hit->getPosition() );
        double dToTrack = (hitPos - trackPosAtEcal).r();
        x.push_back(dToTrack);
        double time = RandGauss::shoot(hit->getTime(), timeResolution);
        y.push_back(time);
    }

    TGraph gr(nHits, x.data(), y.data());
    gr.Fit("pol1", "Q");
    return gr.GetFunction("pol1")->GetParameter(0);
}
