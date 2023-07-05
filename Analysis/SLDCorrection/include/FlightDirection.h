#ifndef FlightDirection_h_1
#define FlightDirection_h_1

#include "lcio.h"
#include "streamlog/streamlog.h"
#include "DDMarlinCED.h"
#include "GeometryUtil.h"

#include "TVector3.h"
#include "TF2.h"

#include "SLDCorrection.h"
#include "AssignParticlestoSLD.h"
#include "SLDCorrectionTypes.h"

//	getTrueFlightDirection gives the true flight direction of parent hadron: trueFlightDirection and the position of start and end vertex of the hadron decays semi-leptonically: trueStartVertex & trueSLDVertex
void getTrueFlightDirection( const MCP &SLDLepton , TVector3 &trueFlightDirection , DoubleVector &trueStartVertex , DoubleVector &trueSLDVertex );

//	getRecoFlightDirection reconstructs the flight direction of parent hadron and returnes the status of the charged lepton: 
//	1- The charged lepton is not reconstructed
//	2- The charged lepton is not associated to a jet
//	3- The charged lepton is associated to the primary vertex
//	4- The charged lepton is associated to a secondary vertex: the flight direction is the direction from the position of the primary vertex to the position of the secondary vertex
//	5- The charged lepton is not associated to a secondary vertex, but there is another vertex in the jet: The intersection point of the charged lepton and the direction of momentum of the vertex is chosen as the secondary vertex. The direction from the primary vertex to the secondary vertex is the flight direction of parent hadron
//	6- There is no vertex in the jet, the charged lepton is not associated to a vertex. the secondary vertex is found by intersecting the track of charged lepton and other single tracks in the jet
//	7- There is no vertex/single track in the jet, the charged lepton is not associated to a vertex. the secondary vertex is found by intersecting the track of charged lepton and leading neutral particle in the jet
int getRecoFlightDirection( const RecoParticle &linkedRecoLepton , TVector3 &recoFlightDirection , double &hadronFlightLength , const Vtx &primaryVertex , const Vtx &startVertex , Vtx &SLDVertices , RecoParticle &SLDVerticesRP , const RecoParticle &assignedJet , const VertexVector &verticesInJet , const PFOVector &PFOswithAloneTracks , float &helicesDistance , int vertexingScenario , TVector3 &daughterHadronFlightDirection , double &daughterHadronFlightDistance , DoubleVector& sldVertexPosition , TVector3 &PCAatLepton , TVector3 &PCAatOtherParticle , double &lepton3DImpactParameter , double &OtherParticle3DImpactParameter );

//	intersectTrackLine finds the minimum distance between a track (track) and a straight line (passing from pointOnLine with the direction of momentumOfLine). position of intersection point on the track (PCAatTrack) and the line (PCAatLine) are found
double intersectTrackLine( const Trk &track , const Vtx &primaryVertex , const TVector3 &momentumOfLine , const DoubleVector &pointOnLine , TVector3 &PCAatTrack , TVector3 &PCAatLine );

//	intersectTrackTrack finds the minimum distance between two tracks (track1 & track2). position of intersection point on each track (PCAatTrack1 & PCAatTrack2) are found
double intersectTrackTrack( const Trk &track1 , const Trk &track2 , TVector3 &PCAatTrack1 , TVector3 &PCAatTrack2 );

//	get3DImpactParameter gives the 3D distance of a track (track) from the primary vertex. Point of the closest approach on the track (PCAatTrack) is found
double get3DImpactParameter( const Trk &track , const Vtx &primaryVertex , TVector3 &PCAatTrack );

//	get3DImpactParameter gives the 3D distance of a line passing from a point (pointOnLine) with direction (momentumOfLine) from the primary vertex. Point of the closest approach on the line (PCAatLine) is found
double get3DImpactParameter( const TVector3 &momentumOfLine , const DoubleVector &pointOnLine , const Vtx &primaryVertex , TVector3 &PCAatLine );

//	drawMCParticles draws a MCParticle in event display
void drawMCParticles( const MCP &motherHadron , const MCP &mcParticle );

//	drawReconstructedParticle draws a Reconstructed Particle in event display
void drawReconstructedParticle( const RecoParticle &RecoParticle , const Vtx &primaryVertex , int colorCharged , int colorNeutral );

#endif
