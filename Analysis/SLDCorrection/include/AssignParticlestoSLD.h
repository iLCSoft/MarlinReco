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

void assignVerticesToSemiLeptonicDecay( PFOVector &assignedParticles , VertexVector &availableVertices , const double &invariantMass , const TVector3 &direction , const Vtx &startVertex , TLorentzVector &visibleFourMomentum );

void sortParticles( PFOVector &sortedParticles , PFOVector &unSortedParticles , const TVector3 &direction );

void sortVertices( VertexVector &sortedVertices , VertexVector &availableVertices , const TVector3 &direction , const Vtx &startVertex );

bool isParticleInVertex( const RecoParticle &particle , const Vtx &vertex );

PFOVector getParticlesWithAloneTracks( const RecoParticle &particle , const RecoParticle &jet , const Vtx &primaryVertex , const VertexVector &vertexVector );

Vtx getParticleVertex( const RecoParticle &particle , const VertexVector &vertexVector , bool &foundParticleInVertex );

RecoParticle getJetAssignedToParticle( const RecoParticle &particle , const PFOVector &jetVector , bool &foundParticleInJet );

VertexVector getVerticesInJet( const RecoParticle &assignedJet , const VertexVector &buildUpVertexVector );

RecoParticle getLeadingChargedParticle( const RecoParticle &Jet );

RecoParticle getLeadingNeutralParticle( const RecoParticle &Jet );

#endif
