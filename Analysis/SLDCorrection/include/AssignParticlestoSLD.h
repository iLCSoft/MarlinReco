#ifndef AssignParticlestoSLD_h_1
#define AssignParticlestoSLD_h_1


#include "lcio.h"
#include "streamlog/streamlog.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "DDMarlinCED.h"
#include <math.h>
#include "SLDCorrectionTypes.h"

//	assignParticlesToSemiLeptonicDecay associates particles (availableParticles: un-sorted, sorts wrt the association cone axis: direction) to the vertex of semi-leptonic decay until invariant mass of associated particles is less than invariantMass. associated particles are assignedParticles
void assignParticlesToSemiLeptonicDecay( PFOVector &assignedParticles , PFOVector &availableParticles , const double &invariantMass , const TVector3 &direction , TLorentzVector &visibleFourMomentum );

//	assignVerticesToSemiLeptonicDecay associates particles from vertices (availableParticles: un-sorted, sorts wrt the association cone axis: direction) to the vertex of semi-leptonic decay until invariant mass of associated particles is less than invariantMass. associated particles are assignedParticles
void assignVerticesToSemiLeptonicDecay( PFOVector &assignedParticles , VertexVector &availableVertices , const double &invariantMass , const TVector3 &direction , const Vtx &startVertex , TLorentzVector &visibleFourMomentum );

//	sortParticles sort particles (unSortedParticles) wrt a given direction and stores in sortedParticles. the particle with lowest angle of its momentum to the direction is the first
void sortParticles( PFOVector &sortedParticles , PFOVector &unSortedParticles , const TVector3 &direction );

//	sortVertices sort vertices (availableVertices) wrt a given direction and stores in sortedParticles. the vertex with lowest angle (direction from startVertex to availableVertices[i]) to the direction is the first
void sortVertices( VertexVector &sortedVertices , VertexVector &availableVertices , const TVector3 &direction , const Vtx &startVertex );

//	isParticleInVertex checks whether a particle (particle) is associated to vertex or not, true: particle is associated to vertex, false: particle is not associated
bool isParticleInVertex( const RecoParticle &particle , const Vtx &vertex );

//	getParticlesWithAloneTracks returns a vector of reconstructed particles with at leat one track that are not associated to the primary vertex or buil-up vertices
PFOVector getParticlesWithAloneTracks( const RecoParticle &particle , const RecoParticle &jet , const Vtx &primaryVertex , const VertexVector &vertexVector );

//	getParticleVertex returns the vertex (from vertexVector) that particle is associated to. if such a vertex is found foundParticleInVertex is modified to true
Vtx getParticleVertex( const RecoParticle &particle , const VertexVector &vertexVector , bool &foundParticleInVertex );

//	getJetAssignedToParticle returns the jet (from jetVector) that particle is associated to. if such a jet is found foundParticleInJet is modified to true
RecoParticle getJetAssignedToParticle( const RecoParticle &particle , const PFOVector &jetVector , bool &foundParticleInJet );

//	getVerticesInJet returns a vector of vertices (seb set of buildUpVertexVector) of which particles are associated to assignedJet
VertexVector getVerticesInJet( const RecoParticle &assignedJet , const VertexVector &buildUpVertexVector );

//	getLeadingChargedParticle returns the leading charged particle in Jet
RecoParticle getLeadingChargedParticle( const RecoParticle &Jet );

//	getLeadingNeutralParticle returns the leading neutral particle in Jet
RecoParticle getLeadingNeutralParticle( const RecoParticle &Jet );

#endif
