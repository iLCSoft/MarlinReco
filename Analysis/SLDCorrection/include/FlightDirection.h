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

void getTrueFlightDirection( const MCP &SLDLepton , TVector3 &trueFlightDirection , DoubleVector &trueStartVertex , DoubleVector &trueSLDVertex );

int getRecoFlightDirection( const RecoParticle &linkedRecoLepton , TVector3 &recoFlightDirection , double &hadronFlightLength , const Vtx &primaryVertex , const Vtx &startVertex , Vtx &SLDVertices , RecoParticle &SLDVerticesRP , const RecoParticle &assignedJet , const VertexVector &verticesInJet , const PFOVector &PFOswithAloneTracks , float &helicesDistance , int vertexingScenario , TVector3 &daughterHadronFlightDirection , double &daughterHadronFlightDistance , DoubleVector& sldVertexPosition , TVector3 &PCAatLepton , TVector3 &PCAatOtherParticle , double &lepton3DImpactParameter , double &OtherParticle3DImpactParameter );

double intersectTrackLine( const Trk &track , const Vtx &primaryVertex , const TVector3 &momentumOfLine , const DoubleVector &pointOnLine , TVector3 &PCAatTrack , TVector3 &PCAatLine );

double intersectTrackTrack( const Trk &track1 , const Trk &track2 , TVector3 &PCAatTrack1 , TVector3 &PCAatTrack2 );

double get3DImpactParameter( const Trk &track , const Vtx &primaryVertex , TVector3 &PCAatTrack );

double get3DImpactParameter( const TVector3 &momentumOfLine , const DoubleVector &pointOnLine , const Vtx &primaryVertex , TVector3 &PCAatLine );

void drawMCParticles( const MCP &motherHadron , const MCP &mcParticle );

void drawReconstructedParticle( const RecoParticle &RecoParticle , const Vtx &primaryVertex , int colorCharged , int colorNeutral );

#endif
