#ifndef flightDirection_h_1
#define flightDirection_h_1


#include "lcio.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Vertex.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/VertexImpl.h>
#include "marlinutil/HelixClass.h"
#include "UTIL/LCRelationNavigator.h"
#include "streamlog/streamlog.h"
#include "DDMarlinCED.h"
#include "GeometryUtil.h"

#include "TVector3.h"
#include "TF2.h"

#include "SLDCorrection.h"
#include "AssignParticlestoSLD.h"
//#include "linkedPFO.h"

typedef std::vector<EVENT::MCParticle*>				mcpVector;
typedef std::vector<EVENT::ReconstructedParticle*>	pfoVector;
typedef std::vector<EVENT::Vertex*>					vtxVector;
typedef std::vector<float>							floatVector;
typedef std::vector<double>							doubleVector;
typedef EVENT::Vertex* 								vertex;
typedef EVENT::ReconstructedParticle*				recoParticle;

void getTrueFlightDirection( EVENT::MCParticle *SLDLepton , TVector3 &trueFlightDirection , std::vector<double> &trueStartVertex , std::vector<double> &trueSLDVertex );

int getRecoFlightDirection( EVENT::ReconstructedParticle *linkedRecoLepton , TVector3 &recoFlightDirection , double & hadronFlightLength , EVENT::Vertex *primaryVertex , EVENT::Vertex *startVertex , vertex &SLDVertices , recoParticle &SLDVerticesRP , EVENT::ReconstructedParticle *assignedJet , vtxVector verticesInJet , pfoVector PFOswithAloneTracks , float &helicesDistance , int vertexingScenario , TVector3 &daughterHadronFlightDirection , double &daughterHadronFlightDistance , doubleVector& sldVertexPosition , TVector3 &PCAatLepton , TVector3 &PCAatOtherParticle , double &lepton3DImpactParameter , double &OtherParticle3DImpactParameter );

double intersectTrackLine( EVENT::Track *track , EVENT::Vertex* primaryVertex , TVector3 momentumOfLine , std::vector<double> pointOnLine , TVector3 &PCAatTrack , TVector3 &PCAatLine );

double intersectTrackTrack( EVENT::Track *track1 , EVENT::Track *track2 , TVector3 &PCAatTrack1 , TVector3 &PCAatTrack2 );

double get3DImpactParameter( EVENT::Track *track , EVENT::Vertex* primaryVertex , TVector3 &PCAatTrack );

double get3DImpactParameter( TVector3 momentumOfLine , std::vector<double> pointOnLine , EVENT::Vertex* primaryVertex , TVector3 &PCAatLine );

void drawMCParticles( EVENT::MCParticle *MotherHadron , EVENT::MCParticle *mcParticle );

void drawReconstructedParticle( EVENT::ReconstructedParticle *recoParticle , EVENT::Vertex *primaryVertex , int colorCharged , int colorNeutral );

#endif
