#ifndef BohdanUtils_h
#define BohdanUtils_h 1

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "EVENT/TrackState.h"
#include "EVENT/MCParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/PIDHandler.h"
#include "UTIL/ILDConf.h"
#include "UTIL/Operators.h" // for debuging <<
#include "DDRec/Vector3D.h"
#include "EVENT/SimTrackerHit.h"
#include <string>
#include <vector>

struct VertexData{
    dd4hep::rec::Vector3D pos;
    std::vector<EVENT::MCParticle* > mcs{};
    bool isPrimary = false;
    bool isV0 = false;
};

int parseLine(char* line);

int getVirtualMemoryUsage();

int getPhysicalMemoryUsage();

void printMemoryUsage();

std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track);

float getParameterFromPID(EVENT::ReconstructedParticle* pfo, UTIL::PIDHandler& pidHandler, std::string algorithmName, std::string parameterName);

EVENT::TrackerHit* getSETHit(EVENT::Track* track);

const EVENT::TrackState* getTrackStateAtCalorimeter(EVENT::Track* track);

IMPL::TrackStateImpl getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrack, EVENT::TrackerHit* hit);

float getTrackWeight(float encodedWeight);

float getClusterWeight(float encodedWeight);

EVENT::ReconstructedParticle* getRelatedReconstructedParticle(EVENT::MCParticle* mc, const UTIL::LCRelationNavigator& mc2pfo, const UTIL::LCRelationNavigator& pfo2mc);

EVENT::MCParticle* getMC(EVENT::ReconstructedParticle* pfo, const UTIL::LCRelationNavigator& pfo2mc);

float getECALBarelRMin();

float getECALEndcapZMin();

dd4hep::rec::Vector3D getPhotonAtCalorimeter(EVENT::MCParticle* mc);

EVENT::SimTrackerHit* getSimTrackerHit(EVENT::TrackerHit* hit, const UTIL::LCRelationNavigator& navToSimTrackerHits);

unsigned long interpolateHexColor(unsigned long startColor, unsigned long endColor, float ratio);

bool isSETHit(const EVENT::TrackerHit* hit);

float getPhiTrue(EVENT::MCParticle* mc);
float getD0True(EVENT::MCParticle* mc);
float getZ0True(EVENT::MCParticle* mc);
float getOmegaTrue(EVENT::MCParticle* mc);
float getTanLTrue(EVENT::MCParticle* mc);
EVENT::MCParticle* getFirstStableParent(EVENT::MCParticle* mc);

bool isBottomPDG(unsigned int pdg);
bool isCharmPDG(unsigned int pdg);
bool isStrangePDG(unsigned int pdg);
bool isMCInternalPDG(unsigned int pdg);
unsigned int getQuarkTypeDecay(EVENT::MCParticle* mc);
unsigned int getQuarksToPythia(EVENT::LCEvent * evt);
bool checkIsHiggsProcess(EVENT::LCEvent * evt);
unsigned int getHiggsDaughters(EVENT::LCEvent * evt);

int getHitCaloType( EVENT::CalorimeterHit* hit );
int getHitCaloID( EVENT::CalorimeterHit* hit );
int getHitCaloLayout( EVENT::CalorimeterHit* hit );
int getHitCaloLayer( EVENT::CalorimeterHit* hit );
bool isEcalHit(EVENT::CalorimeterHit* hit);

std::vector<VertexData> getReconstructableTrueVertices(EVENT::LCEvent* evt);

void mergeCloseVerticies(std::vector<VertexData>& verticies, double mergeDistance);

EVENT::LCObject* getMatchingElement(EVENT::LCCollection* col1, EVENT::LCObject* requestedObject, EVENT::LCCollection* col2);

EVENT::LCCollection* getMCParticleCollection(EVENT::LCEvent* evt);

#endif