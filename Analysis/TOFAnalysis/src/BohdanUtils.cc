#include "BohdanUtils.h"
#include "marlinutil/GeometryUtil.h"
#include "marlinutil/CalorimeterHitType.h"
#include "marlin/VerbosityLevels.h"

#include <cstring>
#include <iomanip>

using dd4hep::rec::LayeredCalorimeterData;
using dd4hep::DetType;
using dd4hep::rec::Vector3D;

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getVirtualMemoryUsage(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int getPhysicalMemoryUsage(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

void printMemoryUsage(){
    int vm = getVirtualMemoryUsage();
    int rm = getPhysicalMemoryUsage();
    streamlog_out(MESSAGE)<<"VM usage: "<<vm/1000.<<"    PM usage: "<<rm/1000.<<"  MB"<<std::endl;
}

std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track){
    std::vector<EVENT::Track*> subTracks;
    // add track itself, which contains VXD+FTD+SIT+TPC hits of the first curl.
    subTracks.push_back(track);

    int nSubTracks = track->getTracks().size();
    if (nSubTracks <= 1) return subTracks;

    UTIL::BitField64 encoder( UTIL::LCTrackerCellID::encoding_string() ) ;
    auto isTPCHit = [&encoder](EVENT::TrackerHit* hit) -> bool {
        encoder.setValue( hit->getCellID0() ) ;
        int subdet = encoder[ UTIL::LCTrackerCellID::subdet() ];
        return subdet == UTIL::ILDDetID::TPC;
    };

    int indexOfFirstTPCCurl = 0;
    for(int i = 0; i < nSubTracks; ++i){
        EVENT::Track* subTrack = track->getTracks()[i];
        //DST files do not store subTracks. Return imidiately.
        if (subTrack == nullptr) return subTracks;
        auto hits = subTrack->getTrackerHits();
        if ( std::find_if(hits.begin(), hits.end(), isTPCHit) != hits.end() ){
            indexOfFirstTPCCurl = i;
            break;
        }
    }

    for(int j=indexOfFirstTPCCurl+1; j < nSubTracks; ++j) subTracks.push_back( track->getTracks()[j] );
    return subTracks;
}

float getParameterFromPID(EVENT::ReconstructedParticle* pfo, UTIL::PIDHandler& pidHandler, std::string algorithmName, std::string parameterName){
    int algorithmID = pidHandler.getAlgorithmID(algorithmName);
    const EVENT::ParticleID& pfoPID = pidHandler.getParticleID(pfo, algorithmID);
    const std::vector<float>& parameters = pfoPID.getParameters();
    int parIdx = pidHandler.getParameterIndex(algorithmID, parameterName);
    return parameters[parIdx]; 
}


EVENT::TrackerHit* getSETHit(EVENT::Track* track){
    std::vector<EVENT::TrackerHit*> hits = track->getTrackerHits();
    UTIL::BitField64 encoder( UTIL::LCTrackerCellID::encoding_string() ) ;
    auto isSETHit = [&encoder](EVENT::TrackerHit* hit) -> bool {
        encoder.setValue( hit->getCellID0() ) ;
        int subdet = encoder[ UTIL::LCTrackerCellID::subdet() ];
        return subdet == UTIL::ILDDetID::SET;
    };
    auto it = std::find_if(hits.begin(), hits.end(), isSETHit);
    if ( it == hits.end() ) return nullptr;
    return *it;
}

const EVENT::TrackState* getTrackStateAtCalorimeter(EVENT::Track* track){
    return getSubTracks(track).back()->getTrackState( EVENT::TrackState::AtCalorimeter );
}

IMPL::TrackStateImpl getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrack, EVENT::TrackerHit* hit){
    IMPL::TrackStateImpl ts;
    double chi2Dummy;
    int ndfDummy;
    marlinTrack->getTrackState(hit, ts, chi2Dummy, ndfDummy);
    return ts;
}

float getTrackWeight(float encodedWeight){
    return float( int(encodedWeight) % 10000 ) / 1000.f;
}

float getClusterWeight(float encodedWeight){
    return float( int(encodedWeight) / 10000 ) / 1000.f;
}

EVENT::ReconstructedParticle* getRelatedReconstructedParticle(EVENT::MCParticle* mc, const UTIL::LCRelationNavigator& mc2pfo, const UTIL::LCRelationNavigator& pfo2mc){
    const std::vector<EVENT::LCObject*>& objects = mc2pfo.getRelatedToObjects(mc);
    const std::vector<float>& encodedWeights = mc2pfo.getRelatedToWeights(mc);
    if ( objects.empty() ) return nullptr; // has NO linked recosntructed particle

    //Decode track and cluster weights
    std::vector<float> trackWeights{};
    std::vector<float> clusterWeights{};
    for (auto encodedWeight : encodedWeights ){
        trackWeights.push_back( getTrackWeight(encodedWeight) );
        clusterWeights.push_back( getClusterWeight(encodedWeight) );
    }

    //Find PFO with the highest track weight. If no track, find highest cluster weight
    auto maxTrackWeightIt = std::max_element( trackWeights.begin(), trackWeights.end() );
    auto maxClusterWeightIt = std::max_element( clusterWeights.begin(), clusterWeights.end() );
    int maxIdx = std::distance(trackWeights.begin(), maxTrackWeightIt);
    if ( *maxTrackWeightIt == 0.f ) maxIdx = std::distance(clusterWeights.begin(), maxClusterWeightIt);

    EVENT::ReconstructedParticle* pfo = static_cast<EVENT::ReconstructedParticle*> (objects[maxIdx]);
    EVENT::MCParticle* mcTest = getMC(pfo, pfo2mc);
    // If MCParticle contributes, but it is not the largest contribution, consider it lost.
    if ( mc != mcTest ) return nullptr;
    return pfo;
}

EVENT::MCParticle* getMC(EVENT::ReconstructedParticle* pfo, const UTIL::LCRelationNavigator& pfo2mc){
    const std::vector<EVENT::LCObject*>& objects = pfo2mc.getRelatedToObjects(pfo);
    const std::vector<float>& encodedWeights = pfo2mc.getRelatedToWeights(pfo);
    if ( objects.empty() ) return nullptr; // has NO linked MCParticles

    //Decode track and cluster weights
    std::vector<float> trackWeights{};
    std::vector<float> clusterWeights{};
    for (auto encodedWeight : encodedWeights ){
        trackWeights.push_back( getTrackWeight(encodedWeight) );
        clusterWeights.push_back( getClusterWeight(encodedWeight) );
    }

    //Find MC with the highest track weight. If no track, find highest cluster weight
    auto maxTrackWeightIt = std::max_element( trackWeights.begin(), trackWeights.end() );
    auto maxClusterWeightIt = std::max_element( clusterWeights.begin(), clusterWeights.end() );
    int maxIdx = std::distance(trackWeights.begin(), maxTrackWeightIt);
    if ( *maxTrackWeightIt == 0.f ) maxIdx = std::distance(clusterWeights.begin(), maxClusterWeightIt);
    return static_cast<EVENT::MCParticle*> (objects[maxIdx]);
}


float getECALBarelRMin(){
    float cm2mm = 10.;
    auto ecalBarrelData = MarlinUtil::getLayeredCalorimeterData(( DetType::CALORIMETER | DetType::ELECTROMAGNETIC | DetType::BARREL),
                                                                ( DetType::AUXILIARY | DetType::FORWARD ) );

    return ecalBarrelData->extent[0]*cm2mm; // rmin in mm
}

float getECALEndcapZMin(){
    float cm2mm = 10.;
    auto ecalEndcapData = MarlinUtil::getLayeredCalorimeterData( ( DetType::CALORIMETER | DetType::ELECTROMAGNETIC | DetType::ENDCAP),
                                                                 ( DetType::AUXILIARY | DetType::FORWARD ) );
    return ecalEndcapData->extent[2]*cm2mm; // zmin in mm
}

dd4hep::rec::Vector3D getPhotonAtCalorimeter(EVENT::MCParticle* mc){
    // find intersection point between photon momentum line and ECAL surface planes
    // https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection

    // startPos - starting position of the photon (l0)
    // mom - momentum of the photon (l)
    Vector3D startPos( mc->getVertex() );
    Vector3D mom( mc->getMomentum() );

    // rMin - minimal radial distance to the barrel ECAL surface to deduce normal vector n and point at the surface p0.
    // zMin - minimal longitudinal distance to the endcap ECAL surface to deduce normal vector n and point at the surface p0.
    float rMin = getECALBarelRMin();
    float zMin = getECALEndcapZMin();


    // ENDCAP plane parameters:
    int direction = (mom.z() > 0) ? 1 : -1;
    Vector3D p0Endcap(0, 0, direction*zMin);
    Vector3D nEndcap = p0Endcap.unit();

    // BARREL plane parameters:
    Vector3D p0Barrel, nBarrel;
    int nSides = 8;
    float step = M_PI/nSides;
    float phi = mom.phi();
    // phi is in the range of singularity point [pi, -pi]. Check this first.

    if( phi < (- M_PI + step) || phi > (M_PI - step) ){
        p0Barrel = Vector3D(rMin, M_PI, M_PI/2., Vector3D::spherical);
        nBarrel = p0Barrel.unit();
    }
    else{
        float ecalPhi = -M_PI + 2*step;
        for ( int i=0; i < nSides-1; ++i ){
            if ( ecalPhi-step <= phi && phi < ecalPhi+step ){
                p0Barrel = Vector3D(rMin, ecalPhi, M_PI/2., Vector3D::spherical);
                nBarrel = p0Barrel.unit();
                break;
            }
            else ecalPhi += 2*step;
        }
    }

    //find intersection point, but don't divide by zero
    if ( mom.z() == 0 ){
        float d = (p0Barrel - startPos).dot(nBarrel)/(mom.dot(nBarrel));
        return startPos + d*mom;
    }
    else if( mom.rho() == 0 ){
        float d = (p0Endcap - startPos).dot(nEndcap)/(mom.dot(nEndcap));
        return startPos + d*mom;
    }
    //choose closest intersection point to the 0,0,0
    float dBarrel = (p0Barrel - startPos).dot(nBarrel)/(mom.dot(nBarrel));
    Vector3D intersectionBarrel = startPos + dBarrel*mom;
    float dEndcap = (p0Endcap - startPos).dot(nEndcap)/(mom.dot(nEndcap));
    Vector3D intersectionEndcap = startPos + dEndcap*mom;
    if ( intersectionBarrel.r() <= intersectionEndcap.r() ) return intersectionBarrel;
    return intersectionEndcap;
};

EVENT::SimTrackerHit* getSimTrackerHit(EVENT::TrackerHit* hit, const UTIL::LCRelationNavigator& navToSimTrackerHits){
    // I merge all tracker hit relation collections in the steering file. ENSURE this happens!
    // Otherwise I need to check every possible tracker hit relation collection, which makes this code x10 longer.
    // In case collection doesn't exist, merging is still happens (I think..) with a warning, which is good.
    if (navToSimTrackerHits.getRelatedToObjects(hit).empty()) return nullptr;
    
    const std::vector<float>& weights = navToSimTrackerHits.getRelatedToWeights(hit);
    int max_i = std::max_element(weights.begin(), weights.end()) - weights.begin();
    EVENT::SimTrackerHit* simHit = static_cast<EVENT::SimTrackerHit*> (navToSimTrackerHits.getRelatedToObjects(hit)[max_i]);
    return simHit;
}

bool isSETHit(const EVENT::TrackerHit* hit){
    if (hit == nullptr) return false;
    UTIL::BitField64 encoder( UTIL::LCTrackerCellID::encoding_string() ) ;
    encoder.setValue( hit->getCellID0() ) ;
    int subdet = encoder[ UTIL::LCTrackerCellID::subdet() ];
    return subdet == UTIL::ILDDetID::SET;
}

float getPhiTrue(EVENT::MCParticle* mc){
    Vector3D mom( mc->getMomentum() );
    return std::atan2( mom.y(), mom.x() );
}

float getD0True(EVENT::MCParticle* mc){
    //Return d0 of the particle given the MC truth infromation of the vertex and momentum
    //Copied from somewhere I hope it has no bugs...
    double bField = MarlinUtil::getBzAtOrigin();
    double ct = 2.99792458E-4;
    Vector3D pos( mc->getVertex() );
    Vector3D mom( mc->getMomentum() );
    double q = mc->getCharge();
    if (std::abs(q) != 1) return -1;
    double radius = mom.rho() / (bField*ct);
    double momPhi = std::atan2( mom.y(), mom.x() );
    double xC = pos.x() + radius*std::cos(momPhi - q * M_PI/2.);
    double yC = pos.y() + radius*std::sin(momPhi - q * M_PI/2.);

    return q*(radius - std::hypot(xC,yC));
}

float getZ0True(EVENT::MCParticle* mc){
    //Return z0 of the particle given the MC truth infromation of the vertex and momentum
    //Copied from somewhere I hope it has no bugs...
    double bField = MarlinUtil::getBzAtOrigin();
    double ct = 2.99792458E-4;
    Vector3D pos( mc->getVertex() );
    Vector3D mom( mc->getMomentum() );

    double q = mc->getCharge();
    if (std::abs(q) != 1) return -1;
    double radius = mom.rho() / (bField*ct);
    double tanL = mom.z()/mom.rho();
    double momPhi = std::atan2( mom.y(), mom.x() );
    double xC = pos.x() + radius*std::cos(momPhi - q * M_PI/2.);
    double yC = pos.y() + radius*std::sin(momPhi - q * M_PI/2.);

    double phiRefPoint = std::atan2( pos.y()-yC, pos.x()-xC );
    double phiAtPCA = std::atan2( -yC, -xC );
    double deltaPhi = phiRefPoint - phiAtPCA;    
    double xCircles = ( -pos.z()*q/(radius*tanL) - deltaPhi) / (2.*M_PI);
    int nCircles, n1, n2;
    if (xCircles >= 0.){
        n1 = int(xCircles);
        n2 = n1+1;
    }
    else{
        n1 = int(xCircles)-1;
        n2 = n1+1;
    }
    if ( std::abs(n1-xCircles) < std::abs(n2-xCircles) ) nCircles = n1;
    else nCircles = n2;

    return pos.z() + q*radius*tanL*(deltaPhi + 2.*M_PI*nCircles);
}

float getOmegaTrue(EVENT::MCParticle* mc){
    double bField = MarlinUtil::getBzAtOrigin();
    double ct = 2.99792458E-4;
    Vector3D mom( mc->getMomentum() );
    return ct*bField / mom.rho();    
}

float getTanLTrue(EVENT::MCParticle* mc){
    Vector3D mom( mc->getMomentum() );
    return mom.z()/mom.rho();
}

EVENT::MCParticle* getFirstStableParent(EVENT::MCParticle* mc){
    //Return first parent up the chain, which is not short-lived resonance but somewhat stable
    // pi/k/p are always assumed to have exactly ONE parent! If not overlay.
    if ( mc->getParents().size() != 1 ) return nullptr;
    EVENT::MCParticle* parent = mc->getParents()[0];
    Vector3D parentStartPoint( parent->getVertex() );
    Vector3D parentEndPoint( parent->getEndpoint() );
    bool isMCInternal =  80 < parent->getPDG() && parent->getPDG() < 101;
    bool isNotShortResonance = (parentEndPoint - parentStartPoint).r() > 0.003;
    if ( isNotShortResonance || isMCInternal ) return parent;
    return getFirstStableParent(parent);
}

bool isBottomPDG(unsigned int pdg){
    std::vector<unsigned int> bottomPDGs = { 511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523, 515, 525, 531, 10531, 533, 10533, 20533, 535, 541, 10541, 543, 10543, 20543, 545, 551, 10551, 100551, 110551, 200551, 210551, 553, 10553, 20553, 30553, 100553, 110553, 120553, 130553, 200553, 210553, 220553, 300553, 9000553, 9010553, 555, 10555, 20555, 100555, 110555, 120555, 200555, 557, 100557, 5122, 5112, 5212, 5222, 5114, 5214, 5224, 5132, 5232, 5312, 5322, 5314, 5324, 5332, 5334, 5142, 5242, 5412, 5422, 5414, 5424, 5342, 5432, 5434, 5442, 5444, 5512, 5522, 5514, 5524, 5532, 5534, 5542, 5544, 5554 };
    return std::find(bottomPDGs.begin(), bottomPDGs.end(), pdg) != bottomPDGs.end();
}

bool isCharmPDG(unsigned int pdg){
    std::vector<unsigned int> charmPDGs = { 411, 421, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 415, 425, 431, 10431, 433, 10433, 20433, 435, 441, 10441, 100441, 443, 10443, 20443, 100443, 30443, 9000443, 9010443, 9020443, 445, 100445, 4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 4324, 4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444 };
    return std::find(charmPDGs.begin(), charmPDGs.end(), pdg) != charmPDGs.end();
}

bool isStrangePDG(unsigned int pdg){
    std::vector<unsigned int> strangePDGs = { 3122, 3222, 3212, 3112, 3224, 3214, 3114, 3322, 3312, 3324, 3314, 3334, 130, 310, 311, 321, 9000311, 9000321, 10311, 10321, 100311, 100321, 9010311, 9010321, 9020311, 9020321, 313, 323, 10313, 10323, 20313, 20323, 100313, 100323, 9000313, 9000323, 30313, 30323, 315, 325, 9000315, 9000325, 10315, 10325, 20315, 20325, 9010315, 9010325, 9020315, 9020325, 317, 327, 9010317, 9010327, 319, 329, 9000319, 9000329 };
    return std::find(strangePDGs.begin(), strangePDGs.end(), pdg) != strangePDGs.end();
}

bool isMCInternalPDG(unsigned int pdg){
    return (80 < pdg && pdg < 101) || (900 < pdg && pdg < 931) || (997 < pdg && pdg < 999) || (1900 < pdg && pdg < 1931) || (2900 < pdg && pdg < 2931) || (3900 < pdg && pdg < 3931);
}

unsigned int getQuarkTypeDecay(EVENT::MCParticle* mc){
    //guess quark type decay. 0 - overlay, 2 - light, 3 - strange, 4 - charm, 5 - bottom
    std::vector<unsigned int> pdgs = {};

    //get PDGs of all the particles up the chain until MC internal or no parents is met.
    EVENT::MCParticle* parent = mc;
    pdgs.push_back( std::abs( parent->getPDG() ) );
    if ( parent->getParents().size() != 1 ) return 0; // probably overlay
    while ( true ){
        parent = parent->getParents()[0];
        if ( isMCInternalPDG( std::abs( parent->getPDG() ) ) ) break;
        pdgs.push_back( std::abs( parent->getPDG() ) );
        if ( parent->getParents().size() != 1 ) return 0; // probably overlay
    }

    if ( std::find_if(pdgs.begin(), pdgs.end(), [](unsigned int pdg){return isBottomPDG(pdg);}) != pdgs.end() ) return 5;
    else if ( std::find_if(pdgs.begin(), pdgs.end(), [](unsigned int pdg){return isCharmPDG(pdg);}) != pdgs.end() ) return 4;
    return 1;
}

std::vector<VertexData> getReconstructableTrueVertices(EVENT::LCEvent * evt){
    EVENT::LCCollection* mcs = getMCParticleCollection(evt);
    UTIL::LCRelationNavigator mc2pfo ( evt->getCollection("MCTruthRecoLink") );
    UTIL::LCRelationNavigator pfo2mc ( evt->getCollection("RecoMCTruthLink") );

    std::vector<VertexData> verticies;
    //Find all MCParticles with reconstructed track.
    for (int i=0; i<mcs->getNumberOfElements(); ++i){
        EVENT::MCParticle* mc = static_cast <EVENT::MCParticle*> ( mcs->getElementAt(i) );
        if (mc == nullptr) continue;
        EVENT::ReconstructedParticle* pfo = getRelatedReconstructedParticle(mc, mc2pfo, pfo2mc);
        if( pfo == nullptr || pfo->getTracks().size() == 0) continue;

        Vector3D pos( mc->getVertex() );
        auto it = std::find_if(verticies.begin(), verticies.end(), [&pos](VertexData vtx){return vtx.pos == pos;} );
        if ( it != verticies.end() ){
            (*it).mcs.push_back(mc);
        }
        else{
            VertexData vtxData;
            vtxData.pos = Vector3D( mc->getVertex() );
            vtxData.mcs.push_back(mc);
            verticies.push_back(vtxData);
        }
    }

    // DEBUG OUTPUT
    streamlog_out(DEBUG2)<<"List of raw true vertices:"<<std::endl;
    for( auto& vtxData : verticies ){
        streamlog_out(DEBUG2)<<" Vertex at ("<<vtxData.pos.r()<<")  MCs:  ";
        for(auto mc : vtxData.mcs) streamlog_out(DEBUG2)<<mc<<"    ";
        streamlog_out(DEBUG2)<<std::endl;
    }

    // Merge all vertices within 3 um and combine their MC particles. And recalculate vertex position as a weighted average
    mergeCloseVerticies(verticies, 0.003);

    // DEBUG OUTPUT
    streamlog_out(DEBUG2)<<"List of raw true vertices:"<<std::endl;
    for( auto& vtxData : verticies ){
        streamlog_out(DEBUG2)<<" Vertex at ("<<vtxData.pos.r()<<")  MCs:  ";
        for(auto* mc : vtxData.mcs) streamlog_out(DEBUG2)<<mc<<"    ";
        streamlog_out(DEBUG2)<<std::endl;
    }


    // Remove any vertices with only single track as such verticies are unreconstructable
    verticies.erase(std::remove_if(verticies.begin(), verticies.end(), [](const VertexData& vtx) {return vtx.mcs.size() < 2;}), verticies.end());

    // DEBUG OUTPUT
    streamlog_out(DEBUG2)<<"List of raw true vertices:"<<std::endl;
    for( auto& vtxData : verticies ){
        streamlog_out(DEBUG2)<<" Vertex at ("<<vtxData.pos.r()<<")  MCs:  ";
        for(auto* mc : vtxData.mcs) streamlog_out(DEBUG2)<<mc<<"    ";
        streamlog_out(DEBUG2)<<std::endl;
    }

    Vector3D primVertexPosTrue( static_cast< EVENT::MCParticle* >(mcs->getElementAt(0))->getVertex() );
    auto compDistanceToTrue = [&](const VertexData& a, const VertexData& b){return (a.pos - primVertexPosTrue).r() < (b.pos - primVertexPosTrue).r();};

    if ( verticies.empty() ) return verticies; // no TRUE primary vertex may hapen. This happens once in ~200k events.
    (*std::min_element(verticies.begin(), verticies.end(), compDistanceToTrue)).isPrimary = true;

    for (auto& vtx : verticies){
        if ( vtx.mcs.size() != 2 ) continue;
        if ( vtx.mcs[0]->getCharge() == vtx.mcs[1]->getCharge() ) continue;
        auto* parent1 = getFirstStableParent( vtx.mcs[0] );
        auto* parent2 = getFirstStableParent( vtx.mcs[1] );
        if ( parent1 == nullptr || parent2 == nullptr || (parent1 != parent2) ) continue;

        unsigned int parentPDG = std::abs( parent1->getPDG() );
        vtx.isV0 = parentPDG == 310 || parentPDG == 3122;
    }
    return verticies;
}

void mergeCloseVerticies(std::vector<VertexData>& verticies, double mergeDistance=0.003){
    for(auto iter1=verticies.begin(); iter1!=verticies.end();){
        bool isMerged = false;
        for(auto iter2=std::next(iter1); iter2!=verticies.end();){
            // streamlog_out(DEBUG1)<<"Beginning merge loop Iter1: "<<std::distance(vtx2mcs.begin(), iter1)<<"/"<<vtx2mcs.size()<<endl;
            // streamlog_out(DEBUG1)<<"Beginning merge loop Iter2: "<<std::distance(vtx2mcs.begin(), iter2)<<"/"<<vtx2mcs.size()<<endl;
            Vector3D pos1 = (*iter1).pos;
            Vector3D pos2 = (*iter2).pos;
            double distance = (pos1-pos2).r();
            //Assuming vertex detector resolution 3um, merge vertices within 3um and recalculate vertex position as a weighted average.
            //NOTE: I merge them in arbitrary order. Practially it doesn't matter, but this may cause a difference in some cases.
            if (distance < mergeDistance){
                isMerged = true;
                streamlog_out(DEBUG1)<<" Vertices are being merged! Resetting the loop iterators!"<<std::endl;
                auto mcs1 = (*iter1).mcs;
                int nMCs1 = mcs1.size();
                auto mcs2 = (*iter2).mcs;
                int nMCs2 = mcs2.size();

                //removing two vertices that are being merged. erase-remove idiom: https://stackoverflow.com/questions/347441/how-can-you-erase-elements-from-a-vector-while-iterating
                verticies.erase(std::remove_if(verticies.begin(), verticies.end(), [&](const VertexData& vtx) {return vtx.mcs == mcs1 || vtx.mcs == mcs2;}), verticies.end());

                VertexData vtxData;
                vtxData.pos = (1./(nMCs1+nMCs2))*(nMCs1*pos1 + nMCs2*pos2);
                vtxData.mcs = mcs1;
                vtxData.mcs.insert( vtxData.mcs.end(), mcs2.begin(), mcs2.end() );
                verticies.push_back(vtxData);
                break;
            }
            else{
                ++iter2;
            }
        }
        if (isMerged) iter1 = verticies.begin();
        else ++iter1;
    }

}


unsigned int getQuarksToPythia(EVENT::LCEvent * evt){
    // Return int, where each quark flavour goes as a separate input digit.
    // 55 - bb, 44 - cc, 33 - ss, ..., 5423 - bcus, etc... Signs are ignored
    std::vector<unsigned int> pdgs;
    EVENT::LCCollection* mcs = getMCParticleCollection(evt);
    for (int i=0; i<mcs->getNumberOfElements(); ++i){
        EVENT::MCParticle* mc = static_cast <EVENT::MCParticle*> ( mcs->getElementAt(i) );
        if ( mc == nullptr || mc->getPDG() != 94 ) continue;
        auto parents = mc->getParents();
        for (auto parent:parents){
            unsigned int quarkPDG = std::abs( parent->getPDG() );
            if ( not ( 0 < quarkPDG && quarkPDG < 10 ) && not (quarkPDG ==  21) ){
                streamlog_out(DEBUG8)<<"Ignoring parents of pdg=94 (pythia input) in the event "<<evt->getEventNumber()<<" with abs(pdg) : "<<quarkPDG<<std::endl;
                continue;
            }
            if ( quarkPDG ==  21 ) pdgs.push_back( 9 );
            else pdgs.push_back( quarkPDG );
        }
    }
    unsigned int encodedPDGs = 0;
    for(size_t i=0; i < pdgs.size(); i++) encodedPDGs += pdgs[i] * int(std::pow(10, i));
    return encodedPDGs;
}

bool checkIsHiggsProcess(EVENT::LCEvent * evt){
    EVENT::LCCollection* mcs = getMCParticleCollection(evt);
    for (int i=0; i<mcs->getNumberOfElements(); ++i){
        auto mc = static_cast <EVENT::MCParticle*> ( mcs->getElementAt(i) );
        if ( std::abs( mc->getPDG() ) == 25 ) return true;
    }
    return false;
}

unsigned int getHiggsDaughters(EVENT::LCEvent * evt){
    std::vector<unsigned int> pdgs;
    EVENT::LCCollection* mcs = getMCParticleCollection(evt);
    for (int i=0; i<mcs->getNumberOfElements(); ++i){
        EVENT::MCParticle* mc = static_cast <EVENT::MCParticle*> ( mcs->getElementAt(i) );
        //IMPORTANT: This assumes only SINGLE Higgs in the event
        if ( std::abs( mc->getPDG() ) == 25 ){
            auto daughters = mc->getDaughters();
            for (auto* daughter : daughters){
                pdgs.push_back( std::abs( daughter->getPDG() ) );
            }
            break;
        }
    }
    unsigned int encodedPDGs = 0;
    for(size_t i=0; i < pdgs.size(); i++) encodedPDGs += pdgs[i] * int(std::pow(10, 2*i));
    return encodedPDGs;
}

int getHitCaloType( EVENT::CalorimeterHit* hit ){
    if (hit == nullptr) return -1;
    return CHT( hit->getType() ).caloType();
}
int getHitCaloID( EVENT::CalorimeterHit* hit ){
    if (hit == nullptr) return -1;
    return CHT( hit->getType() ).caloID();
}
int getHitCaloLayout( EVENT::CalorimeterHit* hit ){
    if (hit == nullptr) return -1;
    return CHT( hit->getType() ).layout();
}

int getHitCaloLayer( EVENT::CalorimeterHit* hit ){
    if (hit == nullptr) return -1;
    return CHT( hit->getType() ).layer();
}

bool isEcalHit(EVENT::CalorimeterHit* hit){
    //NOTE: This does NOT include LumiCal/BeamCal/LHCAL
    //It is not the same as hitType.caloType() == CHT::em
    return ( getHitCaloID(hit) == CHT::ecal);
}

EVENT::LCObject* getMatchingElement(EVENT::LCCollection* col1, EVENT::LCObject* requestedObject, EVENT::LCCollection* col2){
    //Return LCObject from col2 that matches the index of object in col1.
    for (int i=0; i<col1->getNumberOfElements(); ++i){
        if ( requestedObject == col1->getElementAt(i) && i < col2->getNumberOfElements() ) return col2->getElementAt(i);
    }
    //IMPROVE: print collection names. I hope they will be stored in the corresponding in the future. Now they are accesible only through LCEvent
    streamlog_out(WARNING)<<"Could not find a matching object in the collection "<<col2<<" for the requested object "<<requestedObject<<" in the collection "<<col1<<std::endl;
    return nullptr;
}

EVENT::LCCollection* getMCParticleCollection(EVENT::LCEvent* evt){
    // Return collection of MCParticles. Which is "MCParticle" for REC files and "MCParticlesSkimmed" for DST files.
    const std::vector<std::string>* names = evt->getCollectionNames();
    if ( std::find(names->begin(), names->end(), "MCParticle") != names->end() ) return evt->getCollection("MCParticle");
    else if ( std::find(names->begin(), names->end(), "MCParticlesSkimmed") != names->end() ) return evt->getCollection("MCParticlesSkimmed");
    else return nullptr;    
}
