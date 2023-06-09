#include <lcio.h>
#include <EVENT/Track.h>
#include <EVENT/Cluster.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Vertex.h>
#include <UTIL/LCRelationNavigator.h>
#include <marlinutil/WeightedPoints3D.h>
#include <TVector3.h>

using namespace lcio;

LCObject *getRelated(LCObject *p, LCRelationNavigator const &nav, int &_cw, int &_tw);

WeightedPoints3D getWeightedPoints3D(const Cluster * clu, const LCCollection *PandoraClusters = nullptr);

MCParticle *getMCParent(const MCParticle *mcp);

MCParticle *getMCParticle(ReconstructedParticle *p, LCRelationNavigator const &RecoMCTruthNavigtor, LCRelationNavigator const &MCTruthRecoNavigator);