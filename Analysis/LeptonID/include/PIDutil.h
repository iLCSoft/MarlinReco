#include <EVENT/Cluster.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include <TVector3.h>
#include <UTIL/LCRelationNavigator.h>
#include <lcio.h>
#include <marlinutil/WeightedPoints3D.h>

using namespace lcio;

std::tuple<LCObject*, float, float> getRelated(LCObject* p, LCRelationNavigator const& nav);

MCParticle* getMCParent(const MCParticle* mcp);

MCParticle* getMCParticle(ReconstructedParticle* p, LCRelationNavigator const& RecoMCTruthNavigtor,
                          LCRelationNavigator const& MCTruthRecoNavigator, float RecoMCTruthClusterWeightCut,
                          float RecoMCTruthTrackWeightCut, float MCTruthRecoClusterWeightCut,
                          float MCTruthRecoTrackWeightCut);

WeightedPoints3D getWeightedPoints3D(const Cluster* clu, const LCCollection* PandoraClusters);