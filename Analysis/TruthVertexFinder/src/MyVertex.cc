#include "MyVertex.hh"
using EVENT::MCParticle;
using IMPL::VertexImpl;
using std::vector;
namespace TTbarAnalysis {
MyVertex::MyVertex() {}

void MyVertex::__SetMCParticles(const vector<MCParticle*> particles) { myMCParticles = particles; }
vector<MCParticle*>& MyVertex::__GetMCParticles() { return myMCParticles; }
} // namespace TTbarAnalysis
