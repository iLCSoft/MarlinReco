#ifndef BohdanDrawing_h
#define BohdanDrawing_h 1

#include "BohdanUtils.h"
#include "TrackLength.h"

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/MCParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "DDRec/Vector3D.h"
#include "TStyle.h"

#include <vector>

TStyle* getMyStyle();
void displayPFO(EVENT::ReconstructedParticle* pfo);
void displayFTDSimHits(EVENT::LCEvent* evt);
void plotECALTimes(EVENT::Cluster* cluster, dd4hep::rec::Vector3D posAtEcal, dd4hep::rec::Vector3D momAtEcal, EVENT::MCParticle* mc);
void plotTrackParams(const std::vector<HitState>& trackStates, EVENT::MCParticle* mc, float bField);
void displayTOFExplanation(std::vector<EVENT::CalorimeterHit*> allHits, std::vector<EVENT::CalorimeterHit*> selectedHits, double x, double y, double z, double px, double py, double pz);


#endif