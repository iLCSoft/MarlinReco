#ifndef SLDCorrectionTypes_h_1
#define SLDCorrectionTypes_h_1

#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Vertex.h>
#include <GeometryUtil.h>

typedef std::vector<int> IntVector;
typedef std::vector<float> FloatVector;
typedef std::vector<double> DoubleVector;
typedef EVENT::MCParticle* MCP;
typedef EVENT::ReconstructedParticle* RecoParticle;
typedef EVENT::Vertex* Vtx;
typedef EVENT::Track* Trk;
typedef std::vector<EVENT::MCParticle*> MCPVector;
typedef std::vector<EVENT::ReconstructedParticle*> PFOVector;
typedef std::vector<EVENT::Vertex*> VertexVector;
typedef std::vector<std::vector<EVENT::ReconstructedParticle*>> PFOVectorVector;

#endif
