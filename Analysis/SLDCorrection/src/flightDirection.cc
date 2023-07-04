#include "flightDirection.h"
using namespace lcio;
using namespace marlin;

void getTrueFlightDirection( EVENT::MCParticle *SLDLepton , TVector3 &trueFlightDirection , std::vector<double> &trueStartVertex , std::vector<double> &trueSLDVertex )
{
	const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
	trueFlightDirection = TVector3( MotherHadron->getMomentumAtEndpoint()[ 0 ] , MotherHadron->getMomentumAtEndpoint()[ 1 ] , MotherHadron->getMomentumAtEndpoint()[ 2 ] );
	trueFlightDirection.SetMag( 1.0 );
	streamlog_out(DEBUG1) << "	True Flight Direction ( nx , ny , nz ):	(	" << trueFlightDirection.X() << "	,	" << trueFlightDirection.Y() << "	,	" << trueFlightDirection.Z() << "	)" << std::endl;
	trueStartVertex.clear();
	trueStartVertex.push_back( MotherHadron->getVertex()[ 0 ] );
	trueStartVertex.push_back( MotherHadron->getVertex()[ 1 ] );
	trueStartVertex.push_back( MotherHadron->getVertex()[ 2 ] );
	streamlog_out(DEBUG1) << "	True Start Vertex ( x , y , z ):	(	" << trueStartVertex[ 0 ] << "	,	" << trueStartVertex[ 1 ] << "	,	" << trueStartVertex[ 2 ] << "	)" << std::endl;
	trueSLDVertex.clear();
	trueSLDVertex.push_back( MotherHadron->getEndpoint()[ 0 ] );
	trueSLDVertex.push_back( MotherHadron->getEndpoint()[ 1 ] );
	trueSLDVertex.push_back( MotherHadron->getEndpoint()[ 2 ] );
	streamlog_out(DEBUG1) << "	True End Vertex ( x , y , z ):		(	" << trueSLDVertex[ 0 ] << "	,	" << trueSLDVertex[ 1 ] << "	,	" << trueSLDVertex[ 2 ] << "	)" << std::endl;
	return;
}

int getRecoFlightDirection( EVENT::ReconstructedParticle *linkedRecoLepton , TVector3 &recoFlightDirection ,
							double &hadronFlightLength , EVENT::Vertex *primaryVertex , EVENT::Vertex *startVertex ,
							vertex &SLDVertex , recoParticle &SLDVertexRP , EVENT::ReconstructedParticle *assignedJet ,
							std::vector<EVENT::Vertex*> verticesInJet , std::vector<EVENT::ReconstructedParticle*> PFOswithAloneTracks ,
							float &helicesDistance , int vertexingScenario , TVector3 &daughterHadronFlightDirection ,
							double &daughterHadronFlightDistance , doubleVector& sldVertexPosition ,
							TVector3 &PCAatLepton , TVector3 &PCAatOtherParticle , double &lepton3DImpactParameter , double &OtherParticle3DImpactParameter )
{
	int SLDStatus = -999;
	bool recoLeptonIsInVertex = false;
	TVector3 startVertexPosition( startVertex->getPosition()[ 0 ] , startVertex->getPosition()[ 1 ] , startVertex->getPosition()[ 2 ] );
	EVENT::Vertex* recoLeptonVertex = getParticleVertex( linkedRecoLepton , verticesInJet , recoLeptonIsInVertex );
	recoFlightDirection = TVector3( 0.0 , 0.0 , 0.0 );
	hadronFlightLength = 0.0;
	helicesDistance = 0.0;
	daughterHadronFlightDistance = 0.0;
	daughterHadronFlightDirection = TVector3( 0.0 , 0.0 , 0.0 );
	sldVertexPosition.clear();
	if ( recoLeptonIsInVertex )
	{
		SLDStatus = 4;
		streamlog_out(DEBUG1) << "	(" << SLDStatus << ") Lepton from semi-leptonic decay found in a BuildUp Vertex, BuildUp Vertex is used as vertex of semi-leptonic decay!" << std::endl;
		SLDVertex = recoLeptonVertex;
		SLDVertexRP = recoLeptonVertex->getAssociatedParticle();
		recoFlightDirection = TVector3( recoLeptonVertex->getPosition()[ 0 ] - startVertex->getPosition()[ 0 ] , recoLeptonVertex->getPosition()[ 1 ] - startVertex->getPosition()[ 1 ] , recoLeptonVertex->getPosition()[ 2 ] - startVertex->getPosition()[ 2 ] );
		hadronFlightLength = recoFlightDirection.Mag();
		recoFlightDirection.SetMag( 1.0 );
		ReconstructedParticle* recoLeptonVertexRP = recoLeptonVertex->getAssociatedParticle();
		for ( unsigned int i_par = 0 ; i_par < recoLeptonVertexRP->getParticles().size() ; ++i_par )
		{
			if ( recoLeptonVertexRP->getParticles()[ i_par ] != linkedRecoLepton ) daughterHadronFlightDirection += TVector3( ( recoLeptonVertexRP->getParticles()[ i_par ] )->getMomentum() );
		}
		sldVertexPosition.push_back( recoLeptonVertex->getPosition()[ 0 ] );
		sldVertexPosition.push_back( recoLeptonVertex->getPosition()[ 1 ] );
		sldVertexPosition.push_back( recoLeptonVertex->getPosition()[ 2 ] );
		daughterHadronFlightDistance = daughterHadronFlightDirection.Mag();
		TVector3 PCAatTrackInIP( 0.0 , 0.0 , 0.0 );
		lepton3DImpactParameter = get3DImpactParameter( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , PCAatTrackInIP );
		TVector3 PCAatLineInIP( 0.0 , 0.0 , 0.0 );
		OtherParticle3DImpactParameter = get3DImpactParameter( daughterHadronFlightDirection , sldVertexPosition , primaryVertex , PCAatLineInIP );
		streamlog_out(DEBUG1) << "	3D-Impact parameter of Lepton at Interaction Point = " << lepton3DImpactParameter << std::endl;
		streamlog_out(DEBUG1) << "	3D-Impact parameter of direction of secondary Vertex at Interaction Point = " << OtherParticle3DImpactParameter << std::endl;
		daughterHadronFlightDirection.SetMag( 1.0 );
	}
	else if ( verticesInJet.size() != 0 )
	{
		streamlog_out(DEBUG1) << "	Lepton from semi-leptonic decay not found in a BuildUp Vertex, Investigating BuildUp vertices in jet" << std::endl;
		EVENT::Vertex* thirdVertex = NULL;
		//float minDistanceToPrimaryVertex = 1000000.0;
		TVector3 jetAxis( assignedJet->getMomentum() ); jetAxis.SetMag( 1.0 );
		float cosAngleFlightDirectionwJetAxis = -1.0;
		for ( unsigned int i_vtx = 0 ; i_vtx < verticesInJet.size() ; ++i_vtx )
		{
			thirdVertex = verticesInJet[ i_vtx ];

			Track* leptonTrack = linkedRecoLepton->getTracks()[ 0 ];
			TVector3 momentumOfLine( ( thirdVertex->getAssociatedParticle() )->getMomentum() );
			TVector3 primaryVertexPosition( primaryVertex->getPosition()[ 0 ] , primaryVertex->getPosition()[ 1 ] , primaryVertex->getPosition()[ 2 ] );
			TVector3 tertiaryVertexPosition( thirdVertex->getPosition()[ 0 ] , thirdVertex->getPosition()[ 1 ] , thirdVertex->getPosition()[ 2 ] );
			std::vector<double> pointOnLine{}; pointOnLine.clear();
			pointOnLine.push_back( thirdVertex->getPosition()[ 0 ] );
			pointOnLine.push_back( thirdVertex->getPosition()[ 1 ] );
			pointOnLine.push_back( thirdVertex->getPosition()[ 2 ] );
			helicesDistance = intersectTrackLine( leptonTrack , primaryVertex , momentumOfLine , pointOnLine , PCAatLepton , PCAatOtherParticle );
			TVector3 testRecoFlightDirection = PCAatLepton - startVertexPosition; testRecoFlightDirection.SetMag( 1.0 );
			if ( testRecoFlightDirection.Dot( jetAxis ) > cosAngleFlightDirectionwJetAxis )
			{
				cosAngleFlightDirectionwJetAxis = testRecoFlightDirection.Dot( jetAxis );
				recoFlightDirection = PCAatLepton - startVertexPosition;
				hadronFlightLength = recoFlightDirection.Mag();
				daughterHadronFlightDirection = TVector3( thirdVertex->getPosition()[ 0 ] - PCAatLepton.X() , thirdVertex->getPosition()[ 1 ] - PCAatLepton.Y() , thirdVertex->getPosition()[ 2 ] - PCAatLepton.Z() );
				daughterHadronFlightDistance = daughterHadronFlightDirection.Mag();
				daughterHadronFlightDirection.SetMag( 1.0 );
				double sldVertexPositionX = PCAatLepton.X();
				double sldVertexPositionY = PCAatLepton.Y();
				double sldVertexPositionZ = PCAatLepton.Z();
				recoFlightDirection.SetMag( 1.0 );
				SLDVertex = thirdVertex;
				SLDVertexRP = thirdVertex->getAssociatedParticle();
				sldVertexPosition.clear();
				sldVertexPosition.push_back( sldVertexPositionX );
				sldVertexPosition.push_back( sldVertexPositionY );
				sldVertexPosition.push_back( sldVertexPositionZ );
				TVector3 PCAatTrackInIP( 0.0 , 0.0 , 0.0 );
				lepton3DImpactParameter = get3DImpactParameter( leptonTrack , primaryVertex , PCAatTrackInIP );
				TVector3 PCAatLineInIP( 0.0 , 0.0 , 0.0 );
				OtherParticle3DImpactParameter = get3DImpactParameter( momentumOfLine , pointOnLine , primaryVertex , PCAatLineInIP );

				streamlog_out(DEBUG1) << "	(" << SLDStatus << ") There is One BuildUp Vertex in jet" << std::endl;
				streamlog_out(DEBUG1) << "		Intersection point of Lepton and other BuildUp Vertices in jet is used as vertex of semi-leptonic decay!" << std::endl;
				SLDStatus = 5;
				streamlog_out(DEBUG1) << "	3D-Impact parameter of Lepton at Interaction Point = " << lepton3DImpactParameter << std::endl;
				streamlog_out(DEBUG1) << "	3D-Impact parameter of direction of secondary Vertex at Interaction Point = " << OtherParticle3DImpactParameter << std::endl;
			}
			//float distanceToPrimaryVertex = std::sqrt( pow( testVertex->getPosition()[ 0 ] - startVertex->getPosition()[ 0 ] , 2 ) + pow( testVertex->getPosition()[ 1 ] - startVertex->getPosition()[ 1 ] , 2 ) + pow( testVertex->getPosition()[ 2 ] - startVertex->getPosition()[ 2 ] , 2 ) );
			//if ( distanceToPrimaryVertex < minDistanceToPrimaryVertex )
			//{
			//	minDistanceToPrimaryVertex = distanceToPrimaryVertex;
			//	thirdVertex = testVertex;
			//}
		}
		if ( thirdVertex != NULL )
		{
		}
	}
	else
	{
		streamlog_out(DEBUG1) << "	There is NO BuildUp Vertex in jet" << std::endl;
		streamlog_out(DEBUG1) << "	Other scenarios are investigated, but Nu-correction will not be performed" << std::endl;
		//TVector3 PCAatLepton;
		//TVector3 PCAatOtherParticle;
		std::vector<double> chargedParticlePosition( 3 , 0.0 );
		TVector3 chargedParticleMomentum( 0.0 , 0.0 , 0.0 );
		ReconstructedParticle* leadingChargedParticle = getLeadingChargedParticle( assignedJet );
		if ( leadingChargedParticle != NULL )
		{
			chargedParticleMomentum = TVector3( leadingChargedParticle->getMomentum() );
			chargedParticlePosition[ 0 ] = leadingChargedParticle->getReferencePoint()[ 0 ];
			chargedParticlePosition[ 1 ] = leadingChargedParticle->getReferencePoint()[ 1 ];
			chargedParticlePosition[ 2 ] = leadingChargedParticle->getReferencePoint()[ 2 ];
		}
		std::vector<double> neutralParticlePosition( 3 , 0.0 );
		TVector3 neutralParticleMomentum( 0.0 , 0.0 , 0.0 );
		ReconstructedParticle* leadingNeutralParticle = getLeadingNeutralParticle( assignedJet );
		if ( leadingNeutralParticle != NULL )
		{
			neutralParticleMomentum = TVector3( leadingNeutralParticle->getMomentum() );
			neutralParticlePosition[ 0 ] = leadingNeutralParticle->getReferencePoint()[ 0 ];
			neutralParticlePosition[ 1 ] = leadingNeutralParticle->getReferencePoint()[ 1 ];
			neutralParticlePosition[ 1 ] = leadingNeutralParticle->getReferencePoint()[ 2 ];
		}
		if ( vertexingScenario == 2 )
		{
			if ( leadingChargedParticle != NULL )
			{
				SLDStatus = 7;
				if ( leadingChargedParticle->getTracks().size() == 1 )
				{
					helicesDistance = intersectTrackTrack( linkedRecoLepton->getTracks()[ 0 ] , leadingChargedParticle->getTracks()[ 0 ] , PCAatLepton , PCAatOtherParticle );
					TVector3 PCAatTrackInIP( 0.0 , 0.0 , 0.0 );
					lepton3DImpactParameter = get3DImpactParameter( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , PCAatTrackInIP );
					TVector3 PCAatLineInIP( 0.0 , 0.0 , 0.0 );
					OtherParticle3DImpactParameter = get3DImpactParameter( leadingChargedParticle->getTracks()[ 0 ] , primaryVertex , PCAatLineInIP );
					streamlog_out(DEBUG1) << "	3D-Impact parameter of Lepton at Interaction Point = " << lepton3DImpactParameter << std::endl;
					streamlog_out(DEBUG1) << "	3D-Impact parameter of direction of secondary Vertex at Interaction Point = " << OtherParticle3DImpactParameter << std::endl;
				}
				else
				{
					helicesDistance = intersectTrackLine( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , chargedParticleMomentum , chargedParticlePosition , PCAatLepton , PCAatOtherParticle );
					TVector3 PCAatTrackInIP( 0.0 , 0.0 , 0.0 );
					lepton3DImpactParameter = get3DImpactParameter( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , PCAatTrackInIP );
					TVector3 PCAatLineInIP( 0.0 , 0.0 , 0.0 );
					OtherParticle3DImpactParameter = get3DImpactParameter( chargedParticleMomentum , chargedParticlePosition , primaryVertex , PCAatLineInIP );
					streamlog_out(DEBUG1) << "	3D-Impact parameter of Lepton at Interaction Point = " << lepton3DImpactParameter << std::endl;
					streamlog_out(DEBUG1) << "	3D-Impact parameter of direction of secondary Vertex at Interaction Point = " << OtherParticle3DImpactParameter << std::endl;
				}
				recoFlightDirection = PCAatLepton - startVertexPosition;
				hadronFlightLength = recoFlightDirection.Mag();
				sldVertexPosition.push_back( PCAatLepton.X() );
				sldVertexPosition.push_back( PCAatLepton.Y() );
				sldVertexPosition.push_back( PCAatLepton.Z() );
			}
			else
			{
				vertexingScenario = 1;
			}
		}
		else if ( vertexingScenario == 3 )
		{
			if ( leadingNeutralParticle != NULL )
			{
				SLDStatus = 7;
				helicesDistance = intersectTrackLine( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , neutralParticleMomentum , neutralParticlePosition , PCAatLepton , PCAatOtherParticle );
				recoFlightDirection = PCAatLepton - startVertexPosition;
				hadronFlightLength = recoFlightDirection.Mag();
				sldVertexPosition.push_back( PCAatLepton.X() );
				sldVertexPosition.push_back( PCAatLepton.Y() );
				sldVertexPosition.push_back( PCAatLepton.Z() );
				TVector3 PCAatTrackInIP( 0.0 , 0.0 , 0.0 );
				lepton3DImpactParameter = get3DImpactParameter( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , PCAatTrackInIP );
				TVector3 PCAatLineInIP( 0.0 , 0.0 , 0.0 );
				OtherParticle3DImpactParameter = get3DImpactParameter( neutralParticleMomentum , neutralParticlePosition , primaryVertex , PCAatLineInIP );
				streamlog_out(DEBUG1) << "	3D-Impact parameter of Lepton at Interaction Point = " << lepton3DImpactParameter << std::endl;
				streamlog_out(DEBUG1) << "	3D-Impact parameter of direction of secondary Vertex at Interaction Point = " << OtherParticle3DImpactParameter << std::endl;
			}
			else
			{
				vertexingScenario = 1;
			}
		}
		else if ( vertexingScenario == 4  )
		{
			SLDStatus = 6;
			helicesDistance = 10000.0;
			double temphelicesDistance = 0.0;
			for ( unsigned int i_par = 0 ; i_par < PFOswithAloneTracks.size() ; ++i_par )
			{
				TVector3 tempPCAatLepton;
				TVector3 tempPCAatOtherParticle;
				ReconstructedParticle* alonePFO = PFOswithAloneTracks[ i_par ];
				TVector3 alonePFOMomentum( alonePFO->getMomentum() );
				std::vector<double> alonePFOPosition{};
				alonePFOPosition.push_back( alonePFO->getReferencePoint()[ 0 ] );
				alonePFOPosition.push_back( alonePFO->getReferencePoint()[ 1 ] );
				alonePFOPosition.push_back( alonePFO->getReferencePoint()[ 2 ] );
				if ( alonePFO->getTracks().size() == 1 )
				{
					temphelicesDistance = intersectTrackTrack( linkedRecoLepton->getTracks()[ 0 ] , alonePFO->getTracks()[ 0 ] , tempPCAatLepton , tempPCAatOtherParticle );
					TVector3 PCAatTrackInIP( 0.0 , 0.0 , 0.0 );
					lepton3DImpactParameter = get3DImpactParameter( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , PCAatTrackInIP );
					TVector3 PCAatLineInIP( 0.0 , 0.0 , 0.0 );
					OtherParticle3DImpactParameter = get3DImpactParameter( alonePFO->getTracks()[ 0 ] , primaryVertex , PCAatLineInIP );
					streamlog_out(DEBUG1) << "	3D-Impact parameter of Lepton at Interaction Point = " << lepton3DImpactParameter << std::endl;
					streamlog_out(DEBUG1) << "	3D-Impact parameter of direction of secondary Vertex at Interaction Point = " << OtherParticle3DImpactParameter << std::endl;
				}
				else
				{
					temphelicesDistance = intersectTrackLine( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , alonePFOMomentum , alonePFOPosition , tempPCAatLepton , tempPCAatOtherParticle );
					TVector3 PCAatTrackInIP( 0.0 , 0.0 , 0.0 );
					lepton3DImpactParameter = get3DImpactParameter( linkedRecoLepton->getTracks()[ 0 ] , primaryVertex , PCAatTrackInIP );
					TVector3 PCAatLineInIP( 0.0 , 0.0 , 0.0 );
					OtherParticle3DImpactParameter = get3DImpactParameter( alonePFOMomentum , alonePFOPosition , primaryVertex , PCAatLineInIP );
					streamlog_out(DEBUG1) << "	3D-Impact parameter of Lepton at Interaction Point = " << lepton3DImpactParameter << std::endl;
					streamlog_out(DEBUG1) << "	3D-Impact parameter of direction of secondary Vertex at Interaction Point = " << OtherParticle3DImpactParameter << std::endl;
				}
				if ( temphelicesDistance < helicesDistance )
				{
					helicesDistance = temphelicesDistance;
					PCAatLepton = tempPCAatLepton;
					PCAatOtherParticle = tempPCAatOtherParticle;
				}
			}
			recoFlightDirection = PCAatLepton - startVertexPosition;
			hadronFlightLength = recoFlightDirection.Mag();
			sldVertexPosition.push_back( PCAatLepton.X() );
			sldVertexPosition.push_back( PCAatLepton.Y() );
			sldVertexPosition.push_back( PCAatLepton.Z() );
		}
		if ( vertexingScenario == 1 )
		{
			recoFlightDirection = TVector3( assignedJet->getMomentum() );
			hadronFlightLength = -1.0;
			SLDStatus = 7;
		}
		daughterHadronFlightDirection = recoFlightDirection;
		daughterHadronFlightDistance = daughterHadronFlightDirection.Mag();
		recoFlightDirection.SetMag( 1.0 );
		daughterHadronFlightDirection.SetMag( 1.0 );
	}
	if ( SLDVertex != NULL )
	{
		streamlog_out(DEBUG2) << "Lepton is in this vertex:" << std::endl;
		streamlog_out(DEBUG2) << *SLDVertex << std::endl;
	}
	recoFlightDirection.SetMag( 1.0 );
	return SLDStatus;
}

double intersectTrackLine( 	EVENT::Track *track , EVENT::Vertex* primaryVertex ,
				TVector3 momentumOfLine , std::vector<double> pointOnLine ,
				TVector3 &PCAatTrack , TVector3 &PCAatLine )
{
	streamlog_out(DEBUG1) << "<<<<<<<<<<<<<<<<<<<<   Intersecting a Track and a Line   >>>>>>>>>>>>>>>>>>>>" << std::endl;
	streamlog_out(DEBUG1) << "	TRACK:" << std::endl;
	streamlog_out(DEBUG1) << *track << std::endl;
	double trackD0		= track->getD0();
	double trackZ0		= track->getZ0();
	double trackPhi		= track->getPhi();
	double trackOmega	= track->getOmega();
	double trackTanLambda	= track->getTanLambda();
	double trackcharge	= ( track->getOmega() > 0.0 ?  1.0 : -1.0 );
	double xReference	= track->getReferencePoint()[ 0 ];
	double yReference	= track->getReferencePoint()[ 1 ];
	double zReference	= track->getReferencePoint()[ 2 ];
	double xCenter		= xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi );
	double yCenter		= yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi );
	double zCenter		= zReference + trackZ0;

	while ( trackPhi < 0.0 )
	{
		trackPhi += 2.0 * 3.14159265359;
	}
	while ( trackPhi >= 2.0 * 3.14159265359 )
	{
		trackPhi -= 2.0 * 3.14159265359;
	}

	double m_Bfield = MarlinUtil::getBzAtOrigin();
	double trackPt = m_Bfield * 3.0e-4 / std::fabs( track->getOmega() );
	double trackPx = trackPt * std::cos( track->getPhi() ) ;
	double trackPy = trackPt * std::sin( track->getPhi() ) ;
	double trackPz = trackPt * track->getTanLambda() ;
	double trackXs = track->getReferencePoint()[ 0 ] - track->getD0() * std::sin( track->getPhi() );
	double trackYs = track->getReferencePoint()[ 1 ] + track->getD0() * std::cos( track->getPhi() );
	double trackZs = track->getReferencePoint()[ 2 ] + track->getZ0();
	streamlog_out(DEBUG1) << "	Charge of Track :	" << trackcharge << std::endl;
	streamlog_out(DEBUG1) << "	Center of Track ( Xc , Yc , Zc ) :	(	" << xCenter << "	,	" << yCenter << "	,	" << zCenter << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	StartPoint of Track ( Xc , Yc , Zc ) :	(	" << trackXs << "	,	" << trackYs << "	,	" << trackZs << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Momentum of Track ( Px , Py , Pz ) :	(	" << trackPx << "	,	" << trackPy << "	,	" << trackPz << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Track Reference Point (x,y,z) :		(	" << xReference << "	,	" << yReference << "	,	" << zReference << "	)" << std::endl;

////////////////////////////////////////////////////////////////////////////////
//	double X = helixFreeParameter;
//
//	double xTrack		= xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) - std::sin( X ) / trackOmega;
//	double yTrack		= yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) + std::cos( X ) / trackOmega;
//	double zTrack		= zReference + trackZ0 - ( X - trackPhi ) * trackTanLambda / std::fabs( trackOmega );
////////////////////////////////////////////////////////////////////////////////

//	double vMin = -1.0 * ( momentumOfLine.Px() * ( pointOnLine[ 0 ] - primaryVertex->getPosition()[ 0 ] ) + momentumOfLine.Py() * ( pointOnLine[ 1 ] - primaryVertex->getPosition()[ 1 ] ) + momentumOfLine.Pz() * ( pointOnLine[ 2 ] - primaryVertex->getPosition()[ 2 ] ) ) / momentumOfLine.Mag2();
	double vMin = -1.0 * ( ( pointOnLine[ 0 ] - primaryVertex->getPosition()[ 0 ] ) / momentumOfLine.X() + ( pointOnLine[ 1 ] - primaryVertex->getPosition()[ 1 ] ) / momentumOfLine.Y() + ( pointOnLine[ 2 ] - primaryVertex->getPosition()[ 2 ] ) / momentumOfLine.Z() );
	double vMax = 0.0;

	streamlog_out(DEBUG1) << "	DownStreamVertex Position:  (	" << pointOnLine[ 0 ] << "	,	" << pointOnLine[ 1 ] << "	,	" << pointOnLine[ 2 ] << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	DownStreamVertex Momentum:  (	" << momentumOfLine.Px() << "	,	" << momentumOfLine.Py() << "	,	" << momentumOfLine.Pz() << "	)" << std::endl;


////////////////////////////////////////////////////////////////////////////////
//	double y = lineFreeParameter;
//
//	double xDownStream	= point[ 0 ] + Momentum.Px() * y;
//	double yDownStream	= point[ 1 ] + Momentum.Py() * y;
//	double zDownStream	= point[ 2 ] + Momentum.Pz() * y;
////////////////////////////////////////////////////////////////////////////////

	double xMin = trackPhi - 3.14159265359 / 40.0;
	double xMax = trackPhi + 3.14159265359 / 40.0;
	double yMin = vMin;
	double yMax = vMax;

	DDMarlinCED::drawHelix( m_Bfield , trackcharge , trackXs , trackYs , trackZs , trackPx , trackPy , trackPz , 2 , 1 , 0xff0000 , 0.0 , 1500.0 , 2000.0 , 0 );
	double t = 20.0;
	double endPointX = pointOnLine[ 0 ] + t * momentumOfLine.Px();
	double endPointY = pointOnLine[ 1 ] + t * momentumOfLine.Py();
	double endPointZ = pointOnLine[ 2 ] + t * momentumOfLine.Pz();
	ced_line( endPointX , endPointY , endPointZ , pointOnLine[ 0 ] , pointOnLine[ 1 ] , pointOnLine[ 2 ] , 2 , 1 , 0x0061ff );
/*
	double dphi = ( xMax - xMin ) / 100.0;
	double dv = ( yMax - yMin ) / 10.0;
	for ( double phi = xMin ; phi <= xMax ; phi += dphi )
	{
		double trackx = xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) - std::sin( phi ) / trackOmega;
		double tracky = yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) + std::cos( phi ) / trackOmega;
		double trackz = zReference + trackZ0 - ( phi - trackPhi ) * trackTanLambda / trackOmega;
		ced_hit( trackx , tracky , trackz , 2 , 3 , 0xf44336 );
	}
	for ( double v = yMin ; v <= yMax ; v += dv )
	{
		double linex = pointOnLine[ 0 ] + v * momentumOfLine.Px();
		double liney = pointOnLine[ 1 ] + v * momentumOfLine.Py();
		double linez = pointOnLine[ 2 ] + v * momentumOfLine.Pz();
		ced_hit( linex , liney , linez , 2 , 3 , 0x0075df );

	}
*/


	TF2 *distance = new TF2( "distance" , "std::sqrt( std::pow( ( [ 0 ] + [ 1 ] * std::sin( x ) ) - ( [ 7 ] + [ 8 ] * y ) , 2 ) + std::pow( ( [ 2 ] + [ 3 ] * std::cos( x ) ) - ( [ 9 ] + [ 10 ] * y ) , 2 ) + std::pow( ( [ 4 ] + [ 5 ] * ( x - [ 6 ] ) ) - ( [ 11 ] + [ 12 ] * y ) , 2 ) )" , xMin , xMax , yMin , yMax );

	distance->SetParameter( 0 , xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) );
	streamlog_out(DEBUG1) << "	parameter [ 0 ] =  " << xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) << std::endl;
	distance->SetParameter( 1 , -1.0 / trackOmega );
	streamlog_out(DEBUG1) << "	parameter [ 1 ] =  " << -1.0 / trackOmega << std::endl;
	distance->SetParameter( 2 , yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) );
	streamlog_out(DEBUG1) << "	parameter [ 2 ] =  " << yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) << std::endl;
	distance->SetParameter( 3 , 1.0 / trackOmega );
	streamlog_out(DEBUG1) << "	parameter [ 3 ] =  " << 1.0 / trackOmega << std::endl;
	distance->SetParameter( 4 , zReference + trackZ0 );
	streamlog_out(DEBUG1) << "	parameter [ 4 ] =  " << zReference + trackZ0 << std::endl;
	distance->SetParameter( 5 , -1.0 * trackTanLambda / trackOmega );
	streamlog_out(DEBUG1) << "	parameter [ 5 ] =  " << -1.0 * trackTanLambda / ( trackOmega ) << std::endl;
	distance->SetParameter( 6 , trackPhi );
	streamlog_out(DEBUG1) << "	parameter [ 3 ] =  " << trackPhi << std::endl;
	distance->SetParameter( 7 , pointOnLine[ 0 ] );
	streamlog_out(DEBUG1) << "	parameter [ 7 ] =  " << pointOnLine[ 0 ] << std::endl;
	distance->SetParameter( 8 , momentumOfLine.Px() );
	streamlog_out(DEBUG1) << "	parameter [ 8 ] =  " << momentumOfLine.Px() << std::endl;
	distance->SetParameter( 9 , pointOnLine[ 1 ] );
	streamlog_out(DEBUG1) << "	parameter [ 9 ] =  " << pointOnLine[ 1 ] << std::endl;
	distance->SetParameter( 10 , momentumOfLine.Py() );
	streamlog_out(DEBUG1) << "	parameter [ 10 ] =  " << momentumOfLine.Py() << std::endl;
	distance->SetParameter( 11 , pointOnLine[ 2 ] );
	streamlog_out(DEBUG1) << "	parameter [ 11 ] =  " << pointOnLine[ 2 ] << std::endl;
	distance->SetParameter( 12 , momentumOfLine.Pz() );
	streamlog_out(DEBUG1) << "	parameter [ 12 ] = " << momentumOfLine.Pz() << std::endl;

	distance->SetRange( xMin , yMin , xMax , yMax );
//	distance->SetNpx( 1000 );
//	distance->SetNpy( 1000 );
	double minPhi = trackPhi;
	double minV = 0.0;
	double minDistance = distance->GetMinimumXY( minPhi , minV );

	streamlog_out(DEBUG1) << "	Phi at min Distance = " << minPhi << std::endl;
	streamlog_out(DEBUG1) << "	V at min Distance = " << minV << std::endl;

	double xTrackPCA	= xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) - std::sin( minPhi ) / trackOmega;
	double yTrackPCA	= yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) + std::cos( minPhi ) / trackOmega;
	double zTrackPCA	= zReference + trackZ0 - ( minPhi - trackPhi ) * trackTanLambda / trackOmega;
	streamlog_out(DEBUG1) << "	Lepton PCA (x,y,z) = 	( " << xTrackPCA << "	,	" << yTrackPCA << "	,	" << zTrackPCA << "	)" << std::endl;
	ced_hit( xTrackPCA , yTrackPCA , zTrackPCA , 5 , 5 , 0xff0000 );
	PCAatTrack = TVector3( xTrackPCA , yTrackPCA , zTrackPCA );

	double xLine	= pointOnLine[ 0 ] + momentumOfLine.Px() * minV;
	double yLine	= pointOnLine[ 1 ] + momentumOfLine.Py() * minV;
	double zLine	= pointOnLine[ 2 ] + momentumOfLine.Pz() * minV;
	streamlog_out(DEBUG1) << "	DS PCA (x,y,z) = 	( " << xLine << "	,	" << yLine << "	,	" << zLine << "	)" << std::endl;
	ced_line( xLine , yLine , zLine , pointOnLine[ 0 ] , pointOnLine[ 1 ] , pointOnLine[ 2 ] , 5 , 1 , 0xf66bef );
	ced_hit( xLine , yLine , zLine , 5 , 5 , 0xf66bef );
	PCAatLine = TVector3( xLine , yLine , zLine );
	streamlog_out(DEBUG1) << "	Distance of Track and Line = " << minDistance << std::endl;
	return minDistance;
}

double intersectTrackTrack(	EVENT::Track *track1 , EVENT::Track *track2 ,
				TVector3 &PCAatTrack1 ,
				TVector3 &PCAatTrack2 )
{
	double m_Bfield = MarlinUtil::getBzAtOrigin();
	streamlog_out(DEBUG1) << "<<<<<<<<<<<<<<<<<<<<<<<<   Intersecting Two Tracks   >>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	streamlog_out(DEBUG1) << "	TRACK1:" << std::endl;
	streamlog_out(DEBUG1) << *track1 << std::endl;
	double track1D0		= track1->getD0();
	double track1Z0		= track1->getZ0();
	double track1Phi	= track1->getPhi();
	double track1Omega	= track1->getOmega();
	double track1TanLambda	= track1->getTanLambda();
	double track1charge	= ( track1->getOmega() > 0.0 ?  1.0 : -1.0 );
	double xReference1	= track1->getReferencePoint()[ 0 ];
	double yReference1	= track1->getReferencePoint()[ 1 ];
	double zReference1	= track1->getReferencePoint()[ 2 ];
	double xCenter1		= xReference1 + ( 1.0 / track1Omega - track1D0 ) * std::sin( track1Phi );
	double yCenter1		= yReference1 - ( 1.0 / track1Omega - track1D0 ) * std::cos( track1Phi );
	double zCenter1		= zReference1 + track1Z0;

	while ( track1Phi < 0.0 )
	{
		track1Phi += 2.0 * 3.14159265359;
	}
	while ( track1Phi >= 2.0 * 3.14159265359 )
	{
		track1Phi -= 2.0 * 3.14159265359;
	}

	double track1Pt = m_Bfield * 3.0e-4 / std::fabs( track1->getOmega() );
	double track1Px = track1Pt * std::cos( track1->getPhi() ) ;
	double track1Py = track1Pt * std::sin( track1->getPhi() ) ;
	double track1Pz = track1Pt * track1->getTanLambda() ;
	double track1Xs = track1->getReferencePoint()[ 0 ] - track1->getD0() * std::sin( track1->getPhi() );
	double track1Ys = track1->getReferencePoint()[ 1 ] + track1->getD0() * std::cos( track1->getPhi() );
	double track1Zs = track1->getReferencePoint()[ 2 ] + track1->getZ0();
	streamlog_out(DEBUG1) << "	Charge of Track1 :	" << track1charge << std::endl;
	streamlog_out(DEBUG1) << "	Center of Track1 ( Xc , Yc , Zc ) :	(	" << xCenter1 << "	,	" << yCenter1 << "	,	" << zCenter1 << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	StartPoint of Track1 ( Xc , Yc , Zc ) :	(	" << track1Xs << "	,	" << track1Ys << "	,	" << track1Zs << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Momentum of Track1 ( Px , Py , Pz ) :	(	" << track1Px << "	,	" << track1Py << "	,	" << track1Pz << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Track1 Reference Point (x,y,z) :		(	" << xReference1 << "	,	" << yReference1 << "	,	" << zReference1 << "	)" << std::endl;
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "-----------------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "	TRACK2:" << std::endl;
	streamlog_out(DEBUG1) << *track2 << std::endl;
	double track2D0		= track2->getD0();
	double track2Z0		= track2->getZ0();
	double track2Phi	= track2->getPhi();
	double track2Omega	= track2->getOmega();
	double track2TanLambda	= track2->getTanLambda();
	double track2charge	= ( track2->getOmega() > 0.0 ?  1.0 : -1.0 );
	double xReference2	= track2->getReferencePoint()[ 0 ];
	double yReference2	= track2->getReferencePoint()[ 1 ];
	double zReference2	= track2->getReferencePoint()[ 2 ];
	double xCenter2		= xReference2 + ( 1.0 / track2Omega - track2D0 ) * std::sin( track2Phi );
	double yCenter2		= yReference2 - ( 1.0 / track2Omega - track2D0 ) * std::cos( track2Phi );
	double zCenter2		= zReference2 + track2Z0;

	while ( track2Phi < 0.0 )
	{
		track2Phi += 2.0 * 3.14159265359;
	}
	while ( track2Phi >= 2.0 * 3.14159265359 )
	{
		track2Phi -= 2.0 * 3.14159265359;
	}

	double track2Pt = m_Bfield * 3.0e-4 / std::fabs( track2->getOmega() );
	double track2Px = track2Pt * std::cos( track2->getPhi() ) ;
	double track2Py = track2Pt * std::sin( track2->getPhi() ) ;
	double track2Pz = track2Pt * track2->getTanLambda() ;
	double track2Xs = track2->getReferencePoint()[ 0 ] - track2->getD0() * std::sin( track2->getPhi() );
	double track2Ys = track2->getReferencePoint()[ 1 ] + track2->getD0() * std::cos( track2->getPhi() );
	double track2Zs = track2->getReferencePoint()[ 2 ] + track2->getZ0();
	streamlog_out(DEBUG1) << "	Charge of Track2 :	" << track2charge << std::endl;
	streamlog_out(DEBUG1) << "	Center of Track2 ( Xc , Yc , Zc ) :	(	" << xCenter2 << "	,	" << yCenter2 << "	,	" << zCenter2 << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	StartPoint of Track2 ( Xc , Yc , Zc ) :	(	" << track2Xs << "	,	" << track2Ys << "	,	" << track2Zs << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Momentum of Track2 ( Px , Py , Pz ) :	(	" << track2Px << "	,	" << track2Py << "	,	" << track2Pz << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Track2 Reference Point (x,y,z) :		(	" << xReference1 << "	,	" << yReference1 << "	,	" << zReference1 << "	)" << std::endl;
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "-----------------------------------------------------------------------------" << std::endl;

	////////////////////////////////////////////////////////////////////////////////
	//	double X = helixFreeParameter;
	//
	//	double xTrack		= xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) - std::sin( X ) / trackOmega;
	//	double yTrack		= yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) + std::cos( X ) / trackOmega;
	//	double zTrack		= zReference + trackZ0 - ( X - trackPhi ) * trackTanLambda / trackOmega;
	////////////////////////////////////////////////////////////////////////////////

//	double trackLengthToScan2D = 1000.0;
//	double deltaPhi1 = trackLengthToScan2D * track1Omega;
//	double deltaPhi2 = trackLengthToScan2D * track2Omega;
	double xMin = track1Phi - 3.14159265359 / 100.0;
	double xMax = track1Phi + 3.14159265359 / 100.0;
	double yMin = track2Phi - 3.14159265359 / 100.0;
	double yMax = track2Phi + 3.14159265359 / 100.0;


	DDMarlinCED::drawHelix( m_Bfield , track1charge , track1Xs , track1Ys , track1Zs , track1Px , track1Py , track1Pz , 2 , 1 , 0xff0000 , 0.0 , 1500.0 , 2000.0 , 0 );
	DDMarlinCED::drawHelix( m_Bfield , track2charge , track2Xs , track2Ys , track2Zs , track2Px , track2Py , track2Pz , 2 , 1 , 0x0061ff , 0.0 , 1500.0 , 2000.0 , 0 );

	double dphi1 = ( xMax - xMin ) / 100.0;
	double dphi2 = ( yMax - yMin ) / 100.0;
	for ( double phi1 = xMin ; phi1 <= xMax ; phi1 += dphi1 )
	{
		double track1x = xReference1 + ( 1.0 / track1Omega - track1D0 ) * std::sin( track1Phi ) - std::sin( phi1 ) / track1Omega;
		double track1y = yReference1 - ( 1.0 / track1Omega - track1D0 ) * std::cos( track1Phi ) + std::cos( phi1 ) / track1Omega;
		double track1z = zReference1 + track1Z0 - ( phi1 - track1Phi ) * track1TanLambda / track1Omega;
		ced_hit( track1x , track1y , track1z , 2 , 5 , 0xff0000 );
	}
	for ( double phi2 = yMin ; phi2 <= yMax ; phi2 += dphi2 )
	{
		double track2x = xReference2 + ( 1.0 / track2Omega - track2D0 ) * std::sin( track2Phi ) - std::sin( phi2 ) / track2Omega;
		double track2y = yReference2 - ( 1.0 / track2Omega - track2D0 ) * std::cos( track2Phi ) + std::cos( phi2 ) / track2Omega;
		double track2z = zReference2 + track2Z0 - ( phi2 - track2Phi ) * track2TanLambda / track2Omega;
		ced_hit( track2x , track2y , track2z , 2 , 5 , 0x0061ff );

	}

//	dist = std::sqrt( std::pow( ( [ 0 ] + [ 1 ] * std::sin( x ) ) - ( [ 7 ] + [ 8 ] * std::sin( x ) ) , 2 ) + std::pow( ( [ 2 ] + [ 3 ] * std::cos( x ) ) - ( [ 9 ] + [ 10 ] * std::cos( x ) ) , 2 ) + std::pow( ( [ 4 ] + [ 5 ] * x ) - ( [ 11 ] + [ 12 ] * y ) , 2 ) )

	TF2 *distance = new TF2( "distance" , "std::sqrt( std::pow( ( [ 0 ] + [ 1 ] * std::sin( x ) ) - ( [ 6 ] + [ 7 ] * std::sin( y ) ) , 2 ) + std::pow( ( [ 2 ] + [ 3 ] * std::cos( x ) ) - ( [ 8 ] + [ 9 ] * std::cos( y ) ) , 2 ) + std::pow( ( [ 4 ] + [ 5 ] * x ) - ( [ 10 ] + [ 11 ] * y ) , 2 ) )" , xMin , xMax , yMin , yMax );

	distance->SetParameter( 0 , xReference1 + ( 1.0 / track1Omega - track1D0 ) * std::sin( track1Phi ) );
	streamlog_out(DEBUG1) << "	parameter [ 0 ] =  " << xReference1 + ( 1.0 / track1Omega - track1D0 ) * std::sin( track1Phi ) << std::endl;
	distance->SetParameter( 1 , -1.0 / track1Omega );
	streamlog_out(DEBUG1) << "	parameter [ 1 ] =  " << -1.0 / track1Omega << std::endl;
	distance->SetParameter( 2 , yReference1 - ( 1.0 / track1Omega - track1D0 ) * std::cos( track1Phi ) );
	streamlog_out(DEBUG1) << "	parameter [ 2 ] =  " << yReference1 - ( 1.0 / track1Omega - track1D0 ) * std::cos( track1Phi ) << std::endl;
	distance->SetParameter( 3 , 1.0 / track1Omega );
	streamlog_out(DEBUG1) << "	parameter [ 3 ] =  " << 1.0 / track1Omega << std::endl;
	distance->SetParameter( 4 , zReference1 + track1Z0 + track1Phi * track1TanLambda / track1Omega );
	streamlog_out(DEBUG1) << "	parameter [ 4 ] =  " << zReference1 + track1Z0 + track1Phi * track1TanLambda / track1Omega;
	distance->SetParameter( 5 , -1.0 * track1TanLambda / track1Omega );
	streamlog_out(DEBUG1) << "	parameter [ 5 ] =  " << -1.0 * track1TanLambda / ( track1Omega ) << std::endl;
	distance->SetParameter( 6 , xReference2 + ( 1.0 / track2Omega - track2D0 ) * std::sin( track2Phi ) );
	streamlog_out(DEBUG1) << "	parameter [ 6 ] =  " << xReference2 + ( 1.0 / track2Omega - track2D0 ) * std::sin( track2Phi ) << std::endl;
	distance->SetParameter( 7 , -1.0 / track2Omega );
	streamlog_out(DEBUG1) << "	parameter [ 7 ] =  " << -1.0 / track2Omega << std::endl;
	distance->SetParameter( 8 , yReference2 - ( 1.0 / track2Omega - track2D0 ) * std::cos( track2Phi ) );
	streamlog_out(DEBUG1) << "	parameter [ 8 ] =  " << yReference1 - ( 1.0 / track2Omega - track2D0 ) * std::cos( track2Phi ) << std::endl;
	distance->SetParameter( 9 , 1.0 / track2Omega );
	streamlog_out(DEBUG1) << "	parameter [ 9 ] =  " << 1.0 / track2Omega << std::endl;
	distance->SetParameter( 10 , zReference2 + track2Z0 + track2Phi * track2TanLambda / track2Omega );
	streamlog_out(DEBUG1) << "	parameter [ 10 ] =  " << zReference2 + track2Z0 + track2Phi * track2TanLambda / track2Omega;
	distance->SetParameter( 11 , -1.0 * track2TanLambda / track2Omega );
	streamlog_out(DEBUG1) << "	parameter [ 11 ] =  " << -1.0 * track2TanLambda / ( track2Omega ) << std::endl;

	distance->SetRange( xMin , yMin , xMax , yMax );
//	distance->SetNpx( 1000 );
//	distance->SetNpy( 1000 );
	double minPhi1 = track1Phi;
	double minPhi2 = track2Phi;
	double minDistance = distance->GetMinimumXY( minPhi1 , minPhi2 );

	streamlog_out(DEBUG1) << "	Phi at Track 1 = " << minPhi1 << std::endl;
	streamlog_out(DEBUG1) << "	Phi at Track 2 = " << minPhi2 << std::endl;

	double xTrack1PCA	= xReference1 + ( 1.0 / track1Omega - track1D0 ) * std::sin( track1Phi ) - std::sin( minPhi1 ) / track1Omega;
	double yTrack1PCA	= yReference1 - ( 1.0 / track1Omega - track1D0 ) * std::cos( track1Phi ) + std::cos( minPhi1 ) / track1Omega;
	double zTrack1PCA	= zReference1 + track1Z0 - ( minPhi1 - track1Phi ) * track1TanLambda / track1Omega;
	streamlog_out(DEBUG1) << "	Track 1 PCA (x,y,z) = 	( " << xTrack1PCA << "	,	" << yTrack1PCA << "	,	" << zTrack1PCA << "	)" << std::endl;
	ced_hit( xTrack1PCA , yTrack1PCA , zTrack1PCA , 5 , 5 , 0xff0000 );
	PCAatTrack1 = TVector3( xTrack1PCA , yTrack1PCA , zTrack1PCA );

	double xTrack2PCA	= xReference2 + ( 1.0 / track2Omega - track2D0 ) * std::sin( track2Phi ) - std::sin( minPhi2 ) / track2Omega;
	double yTrack2PCA	= yReference2 - ( 1.0 / track2Omega - track2D0 ) * std::cos( track2Phi ) + std::cos( minPhi2 ) / track2Omega;
	double zTrack2PCA	= zReference2 + track2Z0 - ( minPhi2 - track2Phi ) * track2TanLambda / track2Omega;
	streamlog_out(DEBUG1) << "	Track 2 PCA (x,y,z) = 	( " << xTrack2PCA << "	,	" << yTrack2PCA << "	,	" << zTrack2PCA << "	)" << std::endl;
	ced_hit( xTrack2PCA , yTrack2PCA , zTrack2PCA , 5 , 5 , 0x0061ff );
	PCAatTrack2 = TVector3( xTrack2PCA , yTrack2PCA , zTrack2PCA );


	streamlog_out(DEBUG1) << "	Distance of Two Tracks = " << minDistance << std::endl;
	return minDistance;

}

double get3DImpactParameter( EVENT::Track *track , EVENT::Vertex* primaryVertex , TVector3 &PCAatTrack )
{
	streamlog_out(DEBUG1) << "<<<<<<<<<<<<<<<<<<<<   Get 3D Impact Parameter of Track   >>>>>>>>>>>>>>>>>>>>" << std::endl;
	streamlog_out(DEBUG1) << "	TRACK:" << std::endl;
	streamlog_out(DEBUG1) << *track << std::endl;
	double trackD0		= track->getD0();
	double trackZ0		= track->getZ0();
	double trackPhi		= track->getPhi();
	double trackOmega	= track->getOmega();
	double trackTanLambda	= track->getTanLambda();
	double trackcharge	= ( track->getOmega() > 0.0 ?  1.0 : -1.0 );
	double xReference	= track->getReferencePoint()[ 0 ];
	double yReference	= track->getReferencePoint()[ 1 ];
	double zReference	= track->getReferencePoint()[ 2 ];
	double xCenter		= xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi );
	double yCenter		= yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi );
	double zCenter		= zReference + trackZ0;

	while ( trackPhi < 0.0 )
	{
		trackPhi += 2.0 * 3.14159265359;
	}
	while ( trackPhi >= 2.0 * 3.14159265359 )
	{
		trackPhi -= 2.0 * 3.14159265359;
	}

	double m_Bfield = MarlinUtil::getBzAtOrigin();
	double trackPt = m_Bfield * 3.0e-4 / std::fabs( track->getOmega() );
	double trackPx = trackPt * std::cos( track->getPhi() ) ;
	double trackPy = trackPt * std::sin( track->getPhi() ) ;
	double trackPz = trackPt * track->getTanLambda() ;
	double trackXs = track->getReferencePoint()[ 0 ] - track->getD0() * std::sin( track->getPhi() );
	double trackYs = track->getReferencePoint()[ 1 ] + track->getD0() * std::cos( track->getPhi() );
	double trackZs = track->getReferencePoint()[ 2 ] + track->getZ0();
	streamlog_out(DEBUG1) << "	Charge of Track :	" << trackcharge << std::endl;
	streamlog_out(DEBUG1) << "	Center of Track ( Xc , Yc , Zc ) :	(	" << xCenter << "	,	" << yCenter << "	,	" << zCenter << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	StartPoint of Track ( Xc , Yc , Zc ) :	(	" << trackXs << "	,	" << trackYs << "	,	" << trackZs << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Momentum of Track ( Px , Py , Pz ) :	(	" << trackPx << "	,	" << trackPy << "	,	" << trackPz << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Track Reference Point (x,y,z) :		(	" << xReference << "	,	" << yReference << "	,	" << zReference << "	)" << std::endl;

////////////////////////////////////////////////////////////////////////////////
//	double X = helixFreeParameter;
//
//	double xTrack		= xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) - std::sin( X ) / trackOmega;
//	double yTrack		= yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) + std::cos( X ) / trackOmega;
//	double zTrack		= zReference + trackZ0 - ( X - trackPhi ) * trackTanLambda / std::fabs( trackOmega );
////////////////////////////////////////////////////////////////////////////////

	double xMin = trackPhi - 3.14159265359 / 40.0;
	double xMax = trackPhi + 3.14159265359 / 40.0;

	DDMarlinCED::drawHelix( m_Bfield , trackcharge , trackXs , trackYs , trackZs , trackPx , trackPy , trackPz , 2 , 1 , 0xff0000 , 0.0 , 1500.0 , 2000.0 , 0 );

/*
	double dphi = ( xMax - xMin ) / 100.0;
	for ( double phi = xMin ; phi <= xMax ; phi += dphi )
	{
		double trackx = xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) - std::sin( phi ) / trackOmega;
		double tracky = yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) + std::cos( phi ) / trackOmega;
		double trackz = zReference + trackZ0 - ( phi - trackPhi ) * trackTanLambda / trackOmega;
		ced_hit( trackx , tracky , trackz , 1 , 1 , 0xf44336 );
	}
*/
	TF1 *distance = new TF1( "distance" , "std::sqrt( std::pow( ( [ 0 ] + [ 1 ] * std::sin( x ) ) - [ 7 ] , 2 ) + std::pow( ( [ 2 ] + [ 3 ] * std::cos( x ) ) - [ 8 ] , 2 ) + std::pow( ( [ 4 ] + [ 5 ] * ( x - [ 6 ] ) ) - [ 9 ] , 2 ) )" , xMin , xMax );

	distance->SetParameter( 0 , xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) );
	streamlog_out(DEBUG1) << "	parameter [ 0 ] =  " << xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) << std::endl;
	distance->SetParameter( 1 , -1.0 / trackOmega );
	streamlog_out(DEBUG1) << "	parameter [ 1 ] =  " << -1.0 / trackOmega << std::endl;
	distance->SetParameter( 2 , yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) );
	streamlog_out(DEBUG1) << "	parameter [ 2 ] =  " << yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) << std::endl;
	distance->SetParameter( 3 , 1.0 / trackOmega );
	streamlog_out(DEBUG1) << "	parameter [ 3 ] =  " << 1.0 / trackOmega << std::endl;
	distance->SetParameter( 4 , zReference + trackZ0 );
	streamlog_out(DEBUG1) << "	parameter [ 4 ] =  " << zReference + trackZ0 << std::endl;
	distance->SetParameter( 5 , -1.0 * trackTanLambda / trackOmega );
	streamlog_out(DEBUG1) << "	parameter [ 5 ] =  " << -1.0 * trackTanLambda / ( trackOmega ) << std::endl;
	distance->SetParameter( 6 , trackPhi );
	streamlog_out(DEBUG1) << "	parameter [ 6 ] =  " << trackPhi << std::endl;
	distance->SetParameter( 7 , primaryVertex->getPosition()[ 0 ] );
	streamlog_out(DEBUG1) << "	parameter [ 7 ] =  " << primaryVertex->getPosition()[ 0 ] << std::endl;
	distance->SetParameter( 8 , primaryVertex->getPosition()[ 1 ] );
	streamlog_out(DEBUG1) << "	parameter [ 8 ] =  " << primaryVertex->getPosition()[ 1 ] << std::endl;
	distance->SetParameter( 9 , primaryVertex->getPosition()[ 2 ] );
	streamlog_out(DEBUG1) << "	parameter [ 9 ] =  " << primaryVertex->getPosition()[ 2 ] << std::endl;

	distance->SetRange( xMin , xMax );
//	distance->SetNpx( 1000 );
//	distance->SetNpy( 1000 );
	double minPhi = distance->GetMinimumX( xMin , xMax );
	double minDistance = distance->GetMinimum( xMin , xMax );

	streamlog_out(DEBUG1) << "	Phi at min Distance = " << minPhi << std::endl;

	double xTrackPCA	= xReference + ( 1.0 / trackOmega - trackD0 ) * std::sin( trackPhi ) - std::sin( minPhi ) / trackOmega;
	double yTrackPCA	= yReference - ( 1.0 / trackOmega - trackD0 ) * std::cos( trackPhi ) + std::cos( minPhi ) / trackOmega;
	double zTrackPCA	= zReference + trackZ0 - ( minPhi - trackPhi ) * trackTanLambda / trackOmega;
	streamlog_out(DEBUG1) << "	Lepton PCA (x,y,z) = 	( " << xTrackPCA << "	,	" << yTrackPCA << "	,	" << zTrackPCA << "	)" << std::endl;
	ced_hit( xTrackPCA , yTrackPCA , zTrackPCA , 10 , 1 , 0xff0000 );
	PCAatTrack = TVector3( xTrackPCA , yTrackPCA , zTrackPCA );

	streamlog_out(DEBUG1) << "	Distance of Track and Primary Vertex (3D Impact Parameter) = " << minDistance << std::endl;
	return minDistance;
}

double get3DImpactParameter( TVector3 momentumOfLine , std::vector<double> pointOnLine , EVENT::Vertex* primaryVertex , TVector3 &PCAatLine )
{
	streamlog_out(DEBUG1) << "<<<<<<<<<<<<<<<<<<<<   Get 3D Impact Parameter of Line   >>>>>>>>>>>>>>>>>>>>" << std::endl;
	double vMin = -1.0 * ( ( pointOnLine[ 0 ] - primaryVertex->getPosition()[ 0 ] ) / momentumOfLine.X() + ( pointOnLine[ 1 ] - primaryVertex->getPosition()[ 1 ] ) / momentumOfLine.Y() + ( pointOnLine[ 2 ] - primaryVertex->getPosition()[ 2 ] ) / momentumOfLine.Z() );
//	double vMin = -1.0 * ( momentumOfLine.Px() * ( pointOnLine[ 0 ] - primaryVertex->getPosition()[ 0 ] ) + momentumOfLine.Py() * ( pointOnLine[ 1 ] - primaryVertex->getPosition()[ 1 ] ) + momentumOfLine.Pz() * ( pointOnLine[ 2 ] - primaryVertex->getPosition()[ 2 ] ) ) / momentumOfLine.Mag2();
	double vMax = 0.0;

	streamlog_out(DEBUG1) << "	Position of Reference Point on Line:  (	" << pointOnLine[ 0 ] << "	,	" << pointOnLine[ 1 ] << "	,	" << pointOnLine[ 2 ] << "	)" << std::endl;
	streamlog_out(DEBUG1) << "	Direction of Line :  (	" << momentumOfLine.Px() << "	,	" << momentumOfLine.Py() << "	,	" << momentumOfLine.Pz() << "	)" << std::endl;


////////////////////////////////////////////////////////////////////////////////
//	double y = lineFreeParameter;
//
//	double xDownStream	= point[ 0 ] + Momentum.Px() * y;
//	double yDownStream	= point[ 1 ] + Momentum.Py() * y;
//	double zDownStream	= point[ 2 ] + Momentum.Pz() * y;
////////////////////////////////////////////////////////////////////////////////

	double yMin = vMin;
	double yMax = vMax;

	double t = 20.0;
	double endPointX = pointOnLine[ 0 ] + t * momentumOfLine.Px();
	double endPointY = pointOnLine[ 1 ] + t * momentumOfLine.Py();
	double endPointZ = pointOnLine[ 2 ] + t * momentumOfLine.Pz();
	ced_line( endPointX , endPointY , endPointZ , pointOnLine[ 0 ] , pointOnLine[ 1 ] , pointOnLine[ 2 ] , 2 , 1 , 0x0061ff );

/*
	double dv = ( yMax - yMin ) / 10.0;
	for ( double v = yMin ; v <= yMax ; v += dv )
	{
		double linex = pointOnLine[ 0 ] + v * momentumOfLine.Px();
		double liney = pointOnLine[ 1 ] + v * momentumOfLine.Py();
		double linez = pointOnLine[ 2 ] + v * momentumOfLine.Pz();
		ced_hit( linex , liney , linez , 1 , 1 , 0x0075df );
	}
*/


	TF1 *distance = new TF1( "distance" , "std::sqrt( std::pow( ( [ 0 ] + [ 1 ] * x ) - [ 6 ]  , 2 ) + std::pow( ( [ 2 ] + [ 3 ] * x ) - [ 7 ]  , 2 ) + std::pow( ( [ 4 ] + [ 5 ] * x ) - [ 8 ]  , 2 ) )" , yMin , yMax );
	distance->SetParameter( 0 , pointOnLine[ 0 ] );
	streamlog_out(DEBUG1) << "	parameter [ 0 ] =  " << pointOnLine[ 0 ] << std::endl;
	distance->SetParameter( 1 , momentumOfLine.Px() );
	streamlog_out(DEBUG1) << "	parameter [ 1 ] =  " << momentumOfLine.Px() << std::endl;
	distance->SetParameter( 2 , pointOnLine[ 1 ] );
	streamlog_out(DEBUG1) << "	parameter [ 2 ] =  " << pointOnLine[ 1 ] << std::endl;
	distance->SetParameter( 3 , momentumOfLine.Py() );
	streamlog_out(DEBUG1) << "	parameter [ 3 ] =  " << momentumOfLine.Py() << std::endl;
	distance->SetParameter( 4 , pointOnLine[ 2 ] );
	streamlog_out(DEBUG1) << "	parameter [ 4 ] =  " << pointOnLine[ 2 ] << std::endl;
	distance->SetParameter( 5 , momentumOfLine.Pz() );
	streamlog_out(DEBUG1) << "	parameter [ 5 ] = " << momentumOfLine.Pz() << std::endl;
	distance->SetParameter( 6 , primaryVertex->getPosition()[ 0 ] );
	streamlog_out(DEBUG1) << "	parameter [ 6 ] =  " << primaryVertex->getPosition()[ 0 ] << std::endl;
	distance->SetParameter( 7 , primaryVertex->getPosition()[ 1 ] );
	streamlog_out(DEBUG1) << "	parameter [ 7 ] =  " << primaryVertex->getPosition()[ 1 ] << std::endl;
	distance->SetParameter( 8 , primaryVertex->getPosition()[ 2 ] );
	streamlog_out(DEBUG1) << "	parameter [ 8 ] =  " << primaryVertex->getPosition()[ 2 ] << std::endl;

	distance->SetRange( yMin , yMax );
//	distance->SetNpx( 1000 );
//	distance->SetNpy( 1000 );
	double minV = distance->GetMinimumX( yMin , yMax );
	double minDistance = distance->GetMinimum( yMin , yMax );

	streamlog_out(DEBUG1) << "	V at min Distance = " << minV << std::endl;

	double xLine	= pointOnLine[ 0 ] + momentumOfLine.Px() * minV;
	double yLine	= pointOnLine[ 1 ] + momentumOfLine.Py() * minV;
	double zLine	= pointOnLine[ 2 ] + momentumOfLine.Pz() * minV;
	streamlog_out(DEBUG1) << "	DS PCA (x,y,z) = 	( " << xLine << "	,	" << yLine << "	,	" << zLine << "	)" << std::endl;
	ced_hit( xLine , yLine , zLine , 10 , 1 , 0x0061ff );
	PCAatLine = TVector3( xLine , yLine , zLine );
	streamlog_out(DEBUG1) << "	Distance of Track and Line = " << minDistance << std::endl;
	return minDistance;
}



void drawMCParticles( EVENT::MCParticle *MotherHadron , EVENT::MCParticle *mcParticle )
{
	double m_Bfield = MarlinUtil::getBzAtOrigin();
	if ( MotherHadron == mcParticle )
	{
		ced_line( MotherHadron->getEndpoint()[ 0 ] , MotherHadron->getEndpoint()[ 1 ] , MotherHadron->getEndpoint()[ 2 ] , MotherHadron->getVertex()[ 0 ] , MotherHadron->getVertex()[ 1 ] , MotherHadron->getVertex()[ 2 ] , 1 , 3 , 0x7b00ff ); //DRAW MOTHER HADRON IN PURPLE
	}
	else
	{
		ced_line( mcParticle->getEndpoint()[ 0 ] , mcParticle->getEndpoint()[ 1 ] , mcParticle->getEndpoint()[ 2 ] , mcParticle->getVertex()[ 0 ] , mcParticle->getVertex()[ 1 ] , mcParticle->getVertex()[ 2 ] , 1 , 1 , 0x00E0FF ); //DRAW UNSTABLE PARTICLES IN CYAN
	}
	streamlog_out(DEBUG0) << " An unstable decay product: " << std::endl;
	streamlog_out(DEBUG0) << *mcParticle << std::endl;
	for ( unsigned int i_daughter = 0 ; i_daughter < mcParticle->getDaughters().size() ; ++i_daughter )
	{
		MCParticle *duaughter = mcParticle->getDaughters()[ i_daughter ];
		if ( duaughter->getGeneratorStatus() == 1 )
		{
			if ( std::fabs( duaughter->getPDG() ) == 12 || std::fabs( duaughter->getPDG() ) == 14 || std::fabs( duaughter->getPDG() ) == 16 ) //DRAW Neutrinos in GRAY
			{
				ced_line( duaughter->getEndpoint()[ 0 ] , duaughter->getEndpoint()[ 1 ] , duaughter->getEndpoint()[ 2 ] , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , 2 , 2 , 0x949494 );
				streamlog_out(DEBUG0) << " True Neutrino: " << std::endl;
				streamlog_out(DEBUG0) << *duaughter << std::endl;
			}
			else if ( std::fabs( duaughter->getPDG() ) == 11 || std::fabs( duaughter->getPDG() ) == 13 || std::fabs( duaughter->getPDG() ) == 15 ) // DRAW Leptons
			{
				for ( unsigned int i_td = 0 ; i_td < MotherHadron->getDaughters().size() ; ++i_td )
				{
					MCParticle *testDuaughter = MotherHadron->getDaughters()[ i_td ];
					if ( testDuaughter == duaughter ) // DRAW Leptons from SLD in RED
					{
						DDMarlinCED::drawHelix( m_Bfield , +1.0 , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , 100.0 * duaughter->getMomentum()[ 0 ] , 100.0 * duaughter->getMomentum()[ 1 ] , 100.0 * duaughter->getMomentum()[ 2 ] , 1 , 1 , 0xfe1100 , 0.0 , 2100.0 , 3000.0 , 0 );
						streamlog_out(DEBUG0) << " True Lepton: " << std::endl;
						streamlog_out(DEBUG0) << *duaughter << std::endl;
					}
					else // DRAW Other Leptons in ORANG
					{
						DDMarlinCED::drawHelix( m_Bfield , +1.0 , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , 100.0 * duaughter->getMomentum()[ 0 ] , 100.0 * duaughter->getMomentum()[ 1 ] , 100.0 * duaughter->getMomentum()[ 2 ] , 1 , 1 , 0xf6b26b , 0.0 , 2100.0 , 3000.0 , 0 );
						streamlog_out(DEBUG0) << " True Lepton: " << std::endl;
						streamlog_out(DEBUG0) << *duaughter << std::endl;
					}
				}
			}
			else if ( std::fabs( duaughter->getCharge() ) <=0.1 ) //DRAW Neutral Particles in GREEN
			{
				ced_line( duaughter->getEndpoint()[ 0 ] , duaughter->getEndpoint()[ 1 ] , duaughter->getEndpoint()[ 2 ] , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , 1 , 1 , 0x01be4b );
				streamlog_out(DEBUG0) << " True Neutral Particle: " << std::endl;
				streamlog_out(DEBUG0) << *duaughter << std::endl;
			}
			else //DRAW Charged Particles in BLUE
			{
				DDMarlinCED::drawHelix( m_Bfield , duaughter->getCharge() , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , duaughter->getMomentum()[ 0 ] , duaughter->getMomentum()[ 1 ] , duaughter->getMomentum()[ 2 ] , 1 , 1 , 0x5a6ffa , 0.0 , 2100.0 , 3000.0 , 0 );
				streamlog_out(DEBUG0) << " True Charged Particle: " << std::endl;
				streamlog_out(DEBUG0) << *duaughter << std::endl;
			}
		}
		else if ( duaughter->getGeneratorStatus() == 2 )
		{
			drawMCParticles( MotherHadron , duaughter );
		}
		else //DRAW Unstable Particles in GRAY (Neutral) / BLACK (Charged)
		{
			if ( std::fabs( duaughter->getCharge() ) <=0.1 )
			{
				ced_line( duaughter->getEndpoint()[ 0 ] , duaughter->getEndpoint()[ 1 ] , duaughter->getEndpoint()[ 2 ] , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , 1 , 1 , 0xbcbcbc );
				streamlog_out(DEBUG0) << " Other Unstable Neutral Particles (genStatus != 1 && genStatus != 2 ): " << std::endl;
				streamlog_out(DEBUG0) << *duaughter << std::endl;
			}
			else
			{
				DDMarlinCED::drawHelix( m_Bfield , duaughter->getCharge() , duaughter->getVertex()[ 0 ] , duaughter->getVertex()[ 1 ] , duaughter->getVertex()[ 2 ] , duaughter->getMomentum()[ 0 ] , duaughter->getMomentum()[ 1 ] , duaughter->getMomentum()[ 2 ] , 1 , 1 , 0x000000 , 0.0 , 2100.0 , 3000.0 , 0 );
				streamlog_out(DEBUG0) << " Other Unstable Charged Particles (genStatus != 1 && genStatus != 2 ): " << std::endl;
				streamlog_out(DEBUG0) << *duaughter << std::endl;
			}
		}
	}
}

void drawReconstructedParticle( EVENT::ReconstructedParticle *reconstructedParticle , EVENT::Vertex *primaryVertex , int colorCharged , int colorNeutral )
{
	double m_Bfield = MarlinUtil::getBzAtOrigin();
	int nTracks = ( reconstructedParticle->getTracks() ).size();
	if ( nTracks > 0 )
	{
		for ( int i_trk = 0 ; i_trk < nTracks ; ++i_trk )
		{
			Track *track = reconstructedParticle->getTracks()[ i_trk ];
			double trackCharge = ( track->getOmega() > 0.0 ?  1.0 : -1.0 );
			double trackPt = m_Bfield * 3.0e-4 / std::fabs( track->getOmega() );
			double trackPx = trackPt * std::cos( track->getPhi() ) ;
			double trackPy = trackPt * std::sin( track->getPhi() ) ;
			double trackPz = trackPt * track->getTanLambda() ;
			double trackXs = track->getReferencePoint()[ 0 ] - track->getD0() * std::sin( track->getPhi() );
			double trackYs = track->getReferencePoint()[ 1 ] + track->getD0() * std::cos( track->getPhi() );
			double trackZs = track->getReferencePoint()[ 2 ] + track->getZ0();
			DDMarlinCED::drawHelix( m_Bfield , trackCharge , trackXs , trackYs , trackZs , trackPx , trackPy , trackPz , 2 , 2 , colorCharged , 0.0 , 1500.0 , 2000.0 , 0 );
		}
	}
	else
	{
		double startPointX = primaryVertex->getPosition()[ 0 ];
		double startPointY = primaryVertex->getPosition()[ 1 ];
		double startPointZ = primaryVertex->getPosition()[ 2 ];
		double endPointX = ( reconstructedParticle->getClusters()[ 0 ]->getPosition()[ 0 ] - primaryVertex->getPosition()[ 0 ] ) * 2.0 / 3.0;
		double endPointY = ( reconstructedParticle->getClusters()[ 0 ]->getPosition()[ 1 ] - primaryVertex->getPosition()[ 1 ] ) * 2.0 / 3.0;
		double endPointZ = ( reconstructedParticle->getClusters()[ 0 ]->getPosition()[ 2 ] - primaryVertex->getPosition()[ 2 ] ) * 2.0 / 3.0;
		ced_line( endPointX , endPointY , endPointZ , startPointX , startPointY , startPointZ , 2 , 2 , colorNeutral );
	}
}
