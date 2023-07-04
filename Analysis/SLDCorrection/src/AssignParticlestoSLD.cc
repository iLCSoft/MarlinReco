#include "AssignParticlestoSLD.h"

using namespace lcio;
using namespace marlin;

void assignParticlesToSemiLeptonicDecay( PFOVector &assignedParticles , PFOVector &availableParticles , const double &invariantMass , const TVector3 &direction , TLorentzVector &visibleFourMomentum )
{
	PFOVector sortedParticles{};
	sortParticles( sortedParticles , availableParticles , direction );
	for ( unsigned int i_par = 0 ; i_par < sortedParticles.size() ; ++i_par )
	{
		RecoParticle particle = sortedParticles[ i_par ];
		TLorentzVector particleFourMomentum = TLorentzVector( particle->getMomentum()[ 0 ] , particle->getMomentum()[ 1 ] , particle->getMomentum()[ 2 ] , particle->getEnergy() );
		if ( ( visibleFourMomentum + particleFourMomentum ).M() <= invariantMass )
		{
			visibleFourMomentum += particleFourMomentum;
			assignedParticles.push_back( particle );
			streamlog_out(DEBUG1) << "		added one particle to the semi-leptonic decay" << std::endl;
			streamlog_out(DEBUG1) << *particle << std::endl;
		}
		else
		{
			return;
		}
	}
}

void assignVerticesToSemiLeptonicDecay( PFOVector &assignedParticles , VertexVector &availableVertices , const double &invariantMass , const TVector3 &direction , const Vtx &startVertex , TLorentzVector &visibleFourMomentum )
{
	VertexVector sortedVertices{};
	streamlog_out(DEBUG8) << "			Sorting vertices with respect to direction	:	(	" << direction.X() << "	, " << direction.Y() << "	, " << direction.Z() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	Association vertices to the vertex of the semi-leptonic decay	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	streamlog_out(DEBUG8) << "			Initial fourMomentum (Px,Py,Pz,M):	(	" << visibleFourMomentum.Px() << "	, " << visibleFourMomentum.Py() << "	, " << visibleFourMomentum.Pz() << "	, " << visibleFourMomentum.E() << "	, " << visibleFourMomentum.M() << "	)" << std::endl;
	sortVertices( sortedVertices , availableVertices , direction , startVertex );
	PFOVector availableParticles{};
	for ( unsigned int i_vtx = 0 ; i_vtx < sortedVertices.size() ; ++i_vtx )
	{
		TLorentzVector vertexFourMomentum = TLorentzVector( 0.0 , 0.0 , 0.0 , 0.0 );
		for ( unsigned int i_par = 0 ; i_par < ( sortedVertices[ i_vtx ]->getAssociatedParticle() )->getParticles().size() ; ++i_par )
		{
			vertexFourMomentum += TLorentzVector( ( sortedVertices[ i_vtx ]->getAssociatedParticle() )->getParticles()[ i_par ]->getMomentum() , ( sortedVertices[ i_vtx ]->getAssociatedParticle() )->getParticles()[ i_par ]->getEnergy() );
		}
		if ( ( visibleFourMomentum + vertexFourMomentum ).M() <= invariantMass )
		{
			visibleFourMomentum += vertexFourMomentum;
			streamlog_out(DEBUG8) << "		added one vertex to the semi-leptonic decay" << std::endl;
			streamlog_out(DEBUG8) << "			New fourMomentum (Px,Py,Pz,M):	(	" << visibleFourMomentum.Px() << "	, " << visibleFourMomentum.Py() << "	, " << visibleFourMomentum.Pz() << "	, " << visibleFourMomentum.E() << "	, " << visibleFourMomentum.M() << "	)" << std::endl;
			for ( unsigned int i_par = 0 ; i_par < ( sortedVertices[ i_vtx ]->getAssociatedParticle() )->getParticles().size() ; ++i_par )
			{
				assignedParticles.push_back( ( sortedVertices[ i_vtx ]->getAssociatedParticle() )->getParticles()[ i_par ] );
				streamlog_out(DEBUG1) << *( sortedVertices[ i_vtx ]->getAssociatedParticle() )->getParticles()[ i_par ] << std::endl;
			}
		}
		else
		{
			return;
		}
	}
}

void sortParticles( PFOVector &sortedParticles , PFOVector &unSortedParticles , const TVector3 &direction )
{
	double closestAngleCos = -1.0;
	RecoParticle closestParticle = NULL;
	streamlog_out(DEBUG0) << "	" << sortedParticles.size() << " of " << unSortedParticles.size() + sortedParticles.size() << " particles have been sorted wrt ( " << direction.X() << " , " << direction.Y() << " , " << direction.Z() <<" )" << std::endl;
	for ( unsigned int i_par = 0 ; i_par < unSortedParticles.size() ; ++i_par )
	{
		RecoParticle testParticle = unSortedParticles[ i_par ];
		TVector3 momentum = TVector3( testParticle->getMomentum() );
		momentum.SetMag( 1.0 );
		if ( momentum.Dot( direction ) >= closestAngleCos )
		{
			closestAngleCos = momentum.Dot( direction );
			closestParticle = testParticle;

		}
	}
	PFOVector remainingUnSortedParticles{};
	for ( unsigned int i_par = 0 ; i_par < unSortedParticles.size() ; ++i_par )
	{
		RecoParticle testParticle = unSortedParticles[ i_par ];
		if ( testParticle == closestParticle )
		{
			sortedParticles.push_back( testParticle );
		}
		else
		{
			remainingUnSortedParticles.push_back( testParticle );
		}
	}
	streamlog_out(DEBUG0) << "	" << sortedParticles.size() << " particles sorted and " << remainingUnSortedParticles.size() << " particles have not been sorted wrt ( " << direction.X() << " , " << direction.Y() << " , " << direction.Z() <<" )" << std::endl;
	if ( remainingUnSortedParticles.size() != 0 ) sortParticles( sortedParticles , remainingUnSortedParticles , direction );
}

void sortVertices( VertexVector &sortedVertices , VertexVector &unSortedVertices , const TVector3 &direction , const Vtx &startVertex )
{
	double closestAngleCos = -1.0;
	Vtx closestVertex = NULL;
	streamlog_out(DEBUG0) << "	" << sortedVertices.size() << " of " << unSortedVertices.size() + sortedVertices.size() << " vertices have been sorted wrt ( " << direction.X() << " , " << direction.Y() << " , " << direction.Z() <<" )" << std::endl;
	for ( unsigned int i_vtx = 0 ; i_vtx < unSortedVertices.size() ; ++i_vtx )
	{
		TVector3 vertexPosition = TVector3( unSortedVertices[ i_vtx ]->getPosition()[ 0 ] - startVertex->getPosition()[ 0 ] , unSortedVertices[ i_vtx ]->getPosition()[ 1 ] - startVertex->getPosition()[ 1 ] , unSortedVertices[ i_vtx ]->getPosition()[ 2 ] - startVertex->getPosition()[ 2 ] );
		vertexPosition.SetMag( 1.0 );
		if ( vertexPosition.Dot( direction ) >= closestAngleCos )
		{
			closestAngleCos = vertexPosition.Dot( direction );
			closestVertex = unSortedVertices[ i_vtx ];
		}
	}
	VertexVector remainingUnSortedVertices{};
	for ( unsigned int i_vtx = 0 ; i_vtx < unSortedVertices.size() ; ++i_vtx )
	{
		if ( unSortedVertices[ i_vtx ] == closestVertex )
		{
			sortedVertices.push_back( unSortedVertices[ i_vtx ] );
			TVector3 vertexPosition = TVector3( unSortedVertices[ i_vtx ]->getPosition()[ 0 ] - startVertex->getPosition()[ 0 ] , unSortedVertices[ i_vtx ]->getPosition()[ 1 ] - startVertex->getPosition()[ 1 ] , unSortedVertices[ i_vtx ]->getPosition()[ 2 ] - startVertex->getPosition()[ 2 ] );
			streamlog_out(DEBUG8) << "			Sorted One vertex at (x,y,z):	(	" << vertexPosition.X() << "	, " << vertexPosition.Y() << "	, " << vertexPosition.Z() << "	) with Cos(Theta) = " << closestAngleCos << std::endl;
		}
		else
		{
			remainingUnSortedVertices.push_back( unSortedVertices[ i_vtx ] );
		}
	}
	streamlog_out(DEBUG0) << "	" << sortedVertices.size() << " vertices sorted and " << remainingUnSortedVertices.size() << " vertices have not been sorted wrt ( " << direction.X() << " , " << direction.Y() << " , " << direction.Z() <<" )" << std::endl;
	if ( remainingUnSortedVertices.size() != 0 ) sortVertices( sortedVertices , remainingUnSortedVertices , direction , startVertex );
}

bool isParticleInVertex( const RecoParticle &particle , const Vtx &vertex )
{
	bool particleIsInVertex = false;
	streamlog_out(DEBUG0) << "" << std::endl;
	streamlog_out(DEBUG1) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "------------ Looking for particle (" << particle->id() << ") in vertex -------------" << std::endl;
	streamlog_out(DEBUG1) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG0) << "	Looking for particle:" << std::endl;
	streamlog_out(DEBUG0) << *particle << std::endl;
	streamlog_out(DEBUG0) << "	in vertex:" << std::endl;
	streamlog_out(DEBUG0) << *vertex << std::endl;
	RecoParticle associatedParticle = vertex->getAssociatedParticle();
	int nPar = ( associatedParticle->getParticles() ).size();
	streamlog_out(DEBUG0) << "	Associated Particle of the Vertex: " << std::endl;
	streamlog_out(DEBUG0) << *associatedParticle << std::endl;
	for ( int i_par = 0 ; i_par < nPar ; ++i_par )
	{
		RecoParticle testParticle = associatedParticle->getParticles()[ i_par ];
		streamlog_out(DEBUG0) << "	Checking particle " << i_par << " (" << testParticle->id() << ") in vertex" << std::endl;
		if ( particle == testParticle )
		{
			particleIsInVertex = true;
			streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<  Found particle in vertex  >>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			streamlog_out(DEBUG2) << *testParticle << std::endl;
			streamlog_out(DEBUG2) << "" << std::endl;
		}
	}
	return particleIsInVertex;
}

PFOVector getParticlesWithAloneTracks( const RecoParticle &particle , const RecoParticle &jet , const Vtx &primaryVertex , const VertexVector &vertexVector )
{
	PFOVector particlesWithAloneTracks{};
	streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG0) << "-------- Looking for alone tracks in jet " << jet << " in jets ----------" << std::endl;
	streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG0) << *jet << std::endl;
	int nParticles = ( jet->getParticles() ).size();
	bool particleIsInPrimaryVertex = false;
	bool particleIsInAVertex = false;
	for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
	{
		RecoParticle testParticle = jet->getParticles()[ i_particle ];
		if ( ( testParticle->getTracks() ).size() == 0 ) continue;
		if ( testParticle == particle ) continue;
		streamlog_out(DEBUG0) << "	Checking particle " << i_particle << " (" << testParticle << ") in jet" << std::endl;
		streamlog_out(DEBUG0) << "	Looking for particle in primary vertex " << std::endl;
		streamlog_out(DEBUG0) << *primaryVertex << std::endl;
		particleIsInPrimaryVertex = isParticleInVertex( testParticle , primaryVertex );
		for ( unsigned int i_vtx = 0 ; i_vtx < vertexVector.size() ; ++i_vtx )
		{
			if ( !particleIsInAVertex )
			{
				Vertex *vertex = vertexVector[ i_vtx ];
				streamlog_out(DEBUG0) << "	Looking for particle in build up vertex[ " << i_vtx << " ]: " << vertex << std::endl;
				streamlog_out(DEBUG0) << *vertex << std::endl;
				particleIsInAVertex = isParticleInVertex( testParticle , vertex );
			}
		}
		if ( !particleIsInPrimaryVertex && !particleIsInAVertex ) particlesWithAloneTracks.push_back( testParticle );
	}
	return particlesWithAloneTracks;
}

Vtx getParticleVertex( const RecoParticle &particle , const VertexVector &vertexVector , bool &foundParticleInVertex )
{
	foundParticleInVertex = false;
	Vtx particleVertex = NULL;
	streamlog_out(DEBUG0) << "" << std::endl;
	streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG0) << "----------- Looking for particle (" << particle << ") in " << vertexVector.size() << " vertices----------" << std::endl;
	streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG0) << *particle << std::endl;
	for ( unsigned int i_vtx = 0 ; i_vtx < vertexVector.size() ; ++i_vtx )
	{
		Vertex *vertex = vertexVector[ i_vtx ];
		streamlog_out(DEBUG0) << "	Looking for particle in vertex " << i_vtx << std::endl;
		streamlog_out(DEBUG0) << *vertex << std::endl;
		foundParticleInVertex = isParticleInVertex( particle , vertex );
		if ( foundParticleInVertex ) particleVertex = vertex;
	}
	return particleVertex;
}

RecoParticle getJetAssignedToParticle( const RecoParticle &particle , const PFOVector &jetVector , bool &foundParticleInJet )
{
	foundParticleInJet = false;
	RecoParticle assignedJet = NULL;
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "------------- Looking for particle (" << particle->id() << ") in jets --------------" << std::endl;
	streamlog_out(DEBUG1) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << *particle << std::endl;
	for ( unsigned int i_jet = 0 ; i_jet < jetVector.size() ; ++i_jet )
	{
		RecoParticle jet = jetVector[ i_jet ];
		streamlog_out(DEBUG0) << "	Looking for particle in jet " << i_jet << std::endl;
		streamlog_out(DEBUG0) << *jet << std::endl;
		int nParticles = ( jet->getParticles() ).size();
		for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
		{
			RecoParticle testParticle = jet->getParticles()[ i_particle ];
			streamlog_out(DEBUG0) << "	Checking particle " << i_particle << " (" << testParticle->id() << ") in jet" << std::endl;
			if ( testParticle == particle )
			{
				streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<  Found particle in jet  >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG2) << *testParticle << std::endl;
				streamlog_out(DEBUG2) << "" << std::endl;
				assignedJet = jet;
				foundParticleInJet = true;
			}
		}
	}
	return assignedJet;
}

VertexVector getVerticesInJet( const RecoParticle &jet , const VertexVector &vertexVector )
{
	VertexVector jetVertices;
	for ( unsigned int i_vtx = 0 ; i_vtx < vertexVector.size() ; ++i_vtx )
	{
		Vtx vertex = vertexVector[ i_vtx ];
		streamlog_out(DEBUG0) << "	Checking Vertex ( " << i_vtx << " ) " << vertex->id() << std::endl;
		streamlog_out(DEBUG0) << *vertex << std::endl;
		RecoParticle associatedParticle = vertex->getAssociatedParticle();
		streamlog_out(DEBUG0) << "	with associated particle " << associatedParticle->id() << std::endl;
		streamlog_out(DEBUG0) << *associatedParticle << std::endl;
		for ( unsigned int i_par = 0 ; i_par < associatedParticle->getParticles().size() ; ++i_par )
		{
			RecoParticle vertexParticle = associatedParticle->getParticles()[ i_par ];
			int nParticles = ( jet->getParticles() ).size();
			for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
			{
				RecoParticle jetParticle = jet->getParticles()[ i_particle ];
				streamlog_out(DEBUG0) << "	Checking particle " << i_particle << " (" << jetParticle->id() << ") in jet" << std::endl;
				if ( vertexParticle == jetParticle )
				{
					streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<  Found Vertex in Jet  >>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					streamlog_out(DEBUG2) << "---------------------------     Vertex     ---------------------------" << std::endl;
					streamlog_out(DEBUG2) << *vertex << std::endl;
					streamlog_out(DEBUG2) << "---------------     Associated Particle Of Vertex     ----------------" << std::endl;
					streamlog_out(DEBUG2) << *associatedParticle << std::endl;
					streamlog_out(DEBUG2) << "----------------------------     Jet     -----------------------------" << std::endl;
					streamlog_out(DEBUG2) << *jet << std::endl;
					streamlog_out(DEBUG2) << "" << std::endl;
					bool vertexWasAdded = false;
					for ( unsigned int i_v = 0 ; i_v < jetVertices.size() ; ++i_v )
					{
						if ( jetVertices[ i_v ] == vertex ) vertexWasAdded = true;
					}
					if ( !vertexWasAdded ) jetVertices.push_back( vertex );
				}
			}
		}
	}
	return jetVertices;
}

RecoParticle getLeadingChargedParticle( const RecoParticle &Jet )
{
	RecoParticle leadingParticle = NULL;
	double leadingEnergy = 0.0;
	for ( unsigned int i_par = 0 ; i_par < Jet->getParticles().size() ; ++i_par )
	{
		RecoParticle testParticle = Jet->getParticles()[ i_par ];
		if ( testParticle->getTracks().size() >= 1 && testParticle->getEnergy() > leadingEnergy )
		{
			leadingEnergy = testParticle->getEnergy();
			leadingParticle = testParticle;
		}
	}
	return leadingParticle;
}

RecoParticle getLeadingNeutralParticle( const RecoParticle &Jet )
{
	RecoParticle leadingParticle = NULL;
	double leadingEnergy = 0.0;
	for ( unsigned int i_par = 0 ; i_par < Jet->getParticles().size() ; ++i_par )
	{
		RecoParticle testParticle = Jet->getParticles()[ i_par ];
		if ( testParticle->getTracks().size() == 0 && testParticle->getEnergy() > leadingEnergy )
		{
			leadingEnergy = testParticle->getEnergy();
			leadingParticle = testParticle;
		}
	}
	return leadingParticle;
}
