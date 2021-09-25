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
#include "TGraphErrors.h"
#include "TF1.h"
#include "MarlinTrk/MarlinTrkUtils.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/ILDConf.h"
#include "IMPL/TrackImpl.h"

using std::vector;
using std::numeric_limits;
using std::pair;
using EVENT::TrackerHit;
using EVENT::Track;
using EVENT::Cluster;
using EVENT::CalorimeterHit;
using EVENT::TrackState;
using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::Vector3D;
using dd4hep::rec::FixedPadSizeTPCData;
using MarlinTrk::IMarlinTrack;
using MarlinTrk::IMarlinTrkSystem;
using UTIL::LCRelationNavigator;
using UTIL::ILDDetID;
using IMPL::TrackImpl;
using IMPL::TrackStateImpl;
using CLHEP::RandGauss;

bool TOFUtils::sortByRho(EVENT::TrackerHit* a, EVENT::TrackerHit* b){
    Vector3D posA( a->getPosition() );
    Vector3D posB( b->getPosition() );
    return posA.rho() < posB.rho();
}


IMPL::TrackStateImpl TOFUtils::getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrk, EVENT::TrackerHit* hit){
    TrackStateImpl ts;
    double chi2Dummy;
    int ndfDummy;
    marlinTrk->getTrackState(hit, ts, chi2Dummy, ndfDummy);
    return ts;
}


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


double TOFUtils::getHelixArcLength(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
    double omega = ts1.getOmega();
    double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
    double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();
    double dPhi = std::abs( ts2.getPhi() - ts1.getPhi() );
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;

    return std::sqrt( std::pow(dPhi/omega, 2) + std::pow(z2-z1, 2) );
}


double TOFUtils::getHelixLengthAlongZ(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
    double tanL = ts1.getTanLambda();
    double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
    double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();

    return std::abs( (z2-z1)/tanL ) * std::sqrt( 1.+tanL*tanL );
}

double TOFUtils::getHelixNRevolutions(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
    double omega = ts1.getOmega();
    double tanL = ts1.getTanLambda();
    double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
    double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();

    // helix length projected on xy
    double circHelix = std::abs( (z2-z1)/tanL );
    double circFull = 2*M_PI/std::abs(omega);

    return circHelix/circFull;
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
    vector<CalorimeterHit*> selectedHits;

    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D trackPosAtEcal ( tsEcal->getReferencePoint() );
    Vector3D trackMomAtEcal = TOFUtils::getHelixMomAtTrackState(*tsEcal, bField);

    for (int l=0; l<maxEcalLayer; ++l){
        // find the closest hit to the linearly extrapolated track inside ecal for each layer
        // OPTIMIZE: this is probably not the most efficient way to loop over ALL hits "maxEcalLayer" times
        CalorimeterHit* selectedHit = nullptr;
        double closestDistanceToLine = numeric_limits<double>::max();
        for ( auto hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isECALHit = ( hitType.caloID() == CHT::ecal );
            if ( (!isECALHit) || ( int( hitType.layer() ) != l) ) continue;

            Vector3D hitPos( hit->getPosition() );
            double dToLine = (hitPos - trackPosAtEcal).cross(trackMomAtEcal.unit()).r();
            if (dToLine < closestDistanceToLine){
                closestDistanceToLine = dToLine;
                selectedHit = hit;
            }
        }
        if ( selectedHit != nullptr ) selectedHits.push_back(selectedHit);
    }
    return selectedHits;
}


std::vector<EVENT::Track*> TOFUtils::getSubTracks(EVENT::Track* track){
    vector<Track*> subTracks;
    subTracks.push_back(track);

    int nSubTracks = track->getTracks().size();
    if (nSubTracks <= 1) return subTracks;

    int nTPCHits = track->getSubdetectorHitNumbers()[(ILDDetID::TPC)*2-1];
    int nSubTrack0Hits = track->getTracks()[0]->getTrackerHits().size();
    int nSubTrack1Hits = track->getTracks()[1]->getTrackerHits().size();

    //OPTIMIZE: this is not reliable, but seems no other way is possible.
    //Read documentation in the header file for details.
    int startIdx;
    if( std::abs(nTPCHits - nSubTrack0Hits) <= 1  ) startIdx = 1;
    else if ( std::abs(nTPCHits - nSubTrack1Hits) <= 1 ) startIdx = 2;
    else{
        // this shouldn't happen in princinple at all...
        streamlog_out(WARNING)<<"Can't understand which subTrack is responsible for the first TPC hits! Skip adding subTracks."<<std::endl;
        return subTracks;
    }
    for(int j=startIdx; j < nSubTracks; ++j) subTracks.push_back( track->getTracks()[j] );
    return subTracks;
}


std::vector<IMPL::TrackStateImpl> TOFUtils::getTrackStatesPerHit(std::vector<EVENT::Track*> tracks, MarlinTrk::IMarlinTrkSystem* trkSystem, bool extrapolateToEcal, double bField){
    vector<TrackStateImpl> trackStatesPerHit;
    int nTracks = tracks.size();
    for(int i=0; i<nTracks; ++i){
        Track* track = tracks[i];
        vector <TrackerHit*> hits = track->getTrackerHits();
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        vector<float> covMatrix(15);
        // initialize variances
        covMatrix[0]  = 1e+06; //sigma_d0^2
        covMatrix[2]  = 100.; //sigma_phi0^2
        covMatrix[5]  = 0.00001; //sigma_omega^2
        covMatrix[9]  = 1e+06; //sigma_z0^2
        covMatrix[14] = 100.; //sigma_tanl^2
        double maxChi2PerHit = 100.;
        std::unique_ptr<IMarlinTrack> marlinTrk( trkSystem->createTrack() );
        TrackImpl refittedTrack;

        //Need to initialize trackState at last hit
        TrackStateImpl preFit = *track->getTrackState(TrackState::AtLastHit);
        preFit.setCovMatrix( covMatrix );
        int errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::backward, &preFit, bField, maxChi2PerHit);
        if (errorFit != 0) continue;

        vector< pair<TrackerHit*, double> > hitsInFit;
        marlinTrk->getHitsInFit(hitsInFit);

        int nHitsInFit = hitsInFit.size();
        // if first subTrack
        if (i == 0){
            trackStatesPerHit.push_back(*(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));

            //add hits in increasing rho for the FIRST subTrack!!!!!
            for( int j=nHitsInFit-1; j>=0; --j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStatesPerHit.push_back(ts);
            }
        }
        else{
            // check which hit is closer to the last hit of previous fit.
            // and iterate starting from the closest
            Vector3D innerHit ( hitsInFit.back().first->getPosition() );
            Vector3D outerHit ( hitsInFit.front().first->getPosition() );
            Vector3D prevHit ( trackStatesPerHit.back().getReferencePoint() );

            if ( (innerHit - prevHit).r() < (outerHit - prevHit).r() ){
                for( int j=nHitsInFit-1; j>=0; --j ){
                    //iterate in increasing rho
                    TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                    trackStatesPerHit.push_back(ts);
                }
            }
            else{
                for( int j=0; j<nHitsInFit; ++j ){
                    //iterate in decreasing rho
                    TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                    trackStatesPerHit.push_back(ts);
                }
            }
        }
        //if last subTrack
        if (i == nTracks - 1){
            // SET hit is not present in hitsInFit as it is composite hit from strips
            // Add ts at the SET hit manualy which fitter returns with reffited track
            // If LastHit != SET hit, then we duplicate previous track state at last TPC hit
            // isn't pretty, but shouldn't affect the track length
            trackStatesPerHit.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );
            if (extrapolateToEcal) trackStatesPerHit.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtCalorimeter) ) ) );
        }
    }
    // one can maybe use hits of refittedTrack, but they include also hits that had failed in the fit
    // code would look cleaner, but using hits that are failed in fit probably would have worse performance..
    // needs to be checked.
    return trackStatesPerHit;
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

    vector <double> x, y, xErr, yErr;
    for ( auto hit : selectedHits ){
        Vector3D hitPos( hit->getPosition() );
        double dToTrack = (hitPos - trackPosAtEcal).r();
        x.push_back(dToTrack);
        double time = RandGauss::shoot(hit->getTime(), timeResolution);
        y.push_back(time);
        xErr.push_back(0.);
        yErr.push_back(0.3);
        //OPTIMIZE: setting this to 0 is not good for the fit.. So I put random 300ps...
        //Changing this doesn't seem to affect results, although one may want to check it more carefully
    }

    TGraphErrors gr(nHits, &x[0], &y[0], &xErr[0], &yErr[0]);
    gr.Fit("pol1", "Q");
    return gr.GetFunction("pol1")->GetParameter(0);
}


void TOFUtils::debugPrint(){

    auto parseLine = [](char* line){
        // This assumes that a digit will be found and the line ends in " Kb".
        int i = strlen(line);
        const char* p = line;
        while (*p <'0' || *p > '9') p++;
        line[i-3] = '\0';
        i = atoi(p);
        return i;
    };

    auto getValue = [&parseLine](std::string memType="VmSize:"){
        //Note: this value is in KB!
        int memTypeIdx;
        if(memType == "VmSize:") memTypeIdx = 7;
        if(memType == "VmRSS:") memTypeIdx = 6;

        FILE* file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != NULL){
            if (strncmp(line, memType.c_str(), memTypeIdx) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    };

    streamlog_out(DEBUG7)<<"Using VM (MB): "<<getValue("VmSize:")/1000.<<std::endl;
    streamlog_out(DEBUG7)<<"Using RAM (MB): "<<getValue("VmRSS:")/1000.<<std::endl;
}
