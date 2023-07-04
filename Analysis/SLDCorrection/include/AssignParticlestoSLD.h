#ifndef AssignParticlestoSLD_h_1
#define AssignParticlestoSLD_h_1


#include "lcio.h"
#include "streamlog/streamlog.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "DDMarlinCED.h"
#include <math.h>
#include "SLDCorrectionTypes.h"

void assignParticlesToSemiLeptonicDecay( PFOVector &assignedParticles , PFOVector &availableParticles , const double &invariantMass , const TVector3 &direction , TLorentzVector &visibleFourMomentum );
//	associates particles (availableParticles: un-sorted, sorts wrt the association cone axis: direction) to the vertex of semi-leptonic decay until invariant mass of associated particles is less than invariantMass. associated particles are assignedParticles

void assignVerticesToSemiLeptonicDecay( PFOVector &assignedParticles , VertexVector &availableVertices , const double &invariantMass , const TVector3 &direction , const Vtx &startVertex , TLorentzVector &visibleFourMomentum );
//	associates particles from vertices (availableParticles: un-sorted, sorts wrt the association cone axis: direction) to the vertex of semi-leptonic decay until invariant mass of associated particles is less than invariantMass. associated particles are assignedParticles

void sortParticles( PFOVector &sortedParticles , PFOVector &unSortedParticles , const TVector3 &direction );
//	sort particles (unSortedParticles) wrt a given direction and stores in sortedParticles. the particle with lowest angle of its momentum to the direction is the first

void sortVertices( VertexVector &sortedVertices , VertexVector &availableVertices , const TVector3 &direction , const Vtx &startVertex );
//	sort vertices (availableVertices) wrt a given direction and stores in sortedParticles. the vertex with lowest angle (direction from startVertex to availableVertices[i]) to the direction is the first

bool isParticleInVertex( const RecoParticle &particle , const Vtx &vertex );
//	checks whether a particle (particle) is associated to vertex or not, true: particle is associated to vertex, false: particle is not associated

PFOVector getParticlesWithAloneTracks( const RecoParticle &particle , const RecoParticle &jet , const Vtx &primaryVertex , const VertexVector &vertexVector );
//	returns a vector of reconstructed particles with at leat one track that are not associated to the primary vertex or buil-up vertices

Vtx getParticleVertex( const RecoParticle &particle , const VertexVector &vertexVector , bool &foundParticleInVertex );
//	returns the vertex (from vertexVector) that particle is associated to. if such a vertex is found foundParticleInVertex is modified to true

RecoParticle getJetAssignedToParticle( const RecoParticle &particle , const PFOVector &jetVector , bool &foundParticleInJet );
//	returns the jet (from jetVector) that particle is associated to. if such a jet is found foundParticleInJet is modified to true

VertexVector getVerticesInJet( const RecoParticle &assignedJet , const VertexVector &buildUpVertexVector );
//	returns a vector of vertices (seb set of buildUpVertexVector) of which particles are associated to assignedJet

RecoParticle getLeadingChargedParticle( const RecoParticle &Jet );
//	returns the leading charged particle in Jet

RecoParticle getLeadingNeutralParticle( const RecoParticle &Jet );
//	returns the leading neutral particle in Jet

#endif
