#include <lcio.h>
#include <EVENT/Track.h>
#include <EVENT/Cluster.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Vertex.h>
#include <UTIL/LCRelationNavigator.h>
#include <marlinutil/WeightedPoints3D.h>
#include <TVector3.h>

using namespace lcio;

float getCorrectdEdx(float dEdx, Track *track);

LCObject *getRelated(LCObject *p, LCRelationNavigator const &nav, int &_cw, int &_tw);

WeightedPoints3D getWeightedPoints3D(const Cluster * clu, const LCCollection *PandoraClusters = nullptr);

MCParticle *getMCParent(const MCParticle *mcp);

std::vector<float> getPosMomFromTrackState(const TrackState *trackstate, double B);

double intersectTrackLine(Track *track, TVector3 momentumOfLine,
                          std::vector<double> pointOnLine,
                          Vertex *primaryVertex, TVector3 &PCAatTrack,
                          TVector3 &PCAatLine);

const MCParticle *checkSLD(const MCParticle *mcp, const MCParticle *parent);

VertexVec getJetVertices(const ReconstructedParticle *jet);

Vertex *getClosestVtx(const Vertex *primVtx, VertexVec jetVertices);

MCParticle *getMCParticle(ReconstructedParticle *p, LCRelationNavigator const &RecoMCTruthNavigtor, LCRelationNavigator const &MCTruthRecoNavigator);