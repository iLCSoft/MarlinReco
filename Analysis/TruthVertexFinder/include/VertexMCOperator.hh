#include "marlin/VerbosityLevels.h"
#include <EVENT/MCParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include "ConstantStorage.hh"
#include "DecayChain.hh"
#include "MathOperator.hh"
#include "MyVertex.hh"
#include <EVENT/LCCollection.h>
#include <IMPL/VertexImpl.h>
#include <UTIL/LCRelationNavigator.h>
#ifndef _VertexMCOperator_hh
#define _VertexMCOperator_hh
namespace TTbarAnalysis {
class VertexMCOperator {
public:
  //
  //	Constants
  //

  //
  //	Constructors
  //
  VertexMCOperator(EVENT::LCCollection* rel);
  virtual ~VertexMCOperator() {};
  VertexMCOperator(const VertexMCOperator&) = delete;
  VertexMCOperator& operator=(const VertexMCOperator&) = delete;
  //
  //	Methods
  //
  std::vector<EVENT::Vertex*>* Construct(DecayChain* chain);
  void AddProngs(EVENT::Vertex* vertex, std::vector<EVENT::MCParticle*>& particles, bool usingRelation = false);

private:
  //
  //	Data
  //
  EVENT::LCCollection* myRelCollection{};
  /* data */
  //
  //	Private methods
  //
  EVENT::Vertex* construct(EVENT::MCParticle* particle, const double* ip, int pdg, int number);
  void addParticle(EVENT::Vertex* vertex, EVENT::MCParticle* particle);
  EVENT::ReconstructedParticle* translate(EVENT::MCParticle* particle);
  EVENT::ReconstructedParticle* getReco(EVENT::MCParticle* particle);
};
} // namespace TTbarAnalysis
#endif
