#include "TrackToRecoParticleConverter.h"

// ----- include for verbosity dependent logging ---------
// #include "marlin/VerbosityLevels.h"
// #include "marlin/StringParameters.h"
// #define SLM streamlog_out(MESSAGE)

#include <gear/BField.h>
#include <marlin/Global.h>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <cmath>

using namespace std;
using namespace lcio;
using namespace marlin;

TrackToRecoParticleConverter aTrackToRecoParticleConverter;

TrackToRecoParticleConverter::TrackToRecoParticleConverter() : Processor("TrackToRecoParticleConverter") {

  // modify processor description
  _description = "Processor to convert track collection into ReconstructedParticle collection";

  // input collections
  registerInputCollection(LCIO::TRACK, "InputTrackCollection", "Input track collection to convert",
                          _inputTrackCollectionName, std::string("InputTracks"));

  // output collection
  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "OutputPFOCollection", "Ouptut collection",
                           _outputPFOCollectionName, std::string("ConvertedTrackAsPFO"));
}

TrackToRecoParticleConverter::~TrackToRecoParticleConverter() {}

void TrackToRecoParticleConverter::init() {
  // usually a good idea to
  printParameters();
}

void TrackToRecoParticleConverter::processRunHeader(LCRunHeader* /*run*/) {}

void TrackToRecoParticleConverter::processEvent(LCEvent* evt) {
  LCCollection* trkCol = evt->getCollection(_inputTrackCollectionName);
  LCCollectionVec* outCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  const double c = 2.99792458e8; // m*s^-1
  const double B = 3.5;          // Tesla
  const double mm2m = 1e-3;
  const double eV2GeV = 1e-9;
  const double eB = B * c * mm2m * eV2GeV;
  const double pi_mass = 0.13957018;
  const double pi_mass_sq = pi_mass * pi_mass;

  unsigned int ntrk = trkCol->getNumberOfElements();
  for (unsigned int i = 0; i < ntrk; i++) {
    Track* trk = dynamic_cast<Track*>(trkCol->getElementAt(i));
    ReconstructedParticleImpl* rp = new ReconstructedParticleImpl();
    rp->addTrack(trk);

    double om = trk->getOmega();
    double td = trk->getTanLambda();
    double cd = 1. / sqrt(1 + td * td);
    double sd = td * cd;
    double cph = cos(trk->getPhi());
    double sph = sin(trk->getPhi());
    double pT = eB / fabs(om);
    double p = pT / cd;

    double px = p * cd * cph;
    double py = p * cd * sph;
    double pz = p * sd;
    double e = sqrt(pi_mass_sq + px * px + py * py + pz * pz);
    double mom[3] = {px, py, pz};
    rp->setMomentum(mom);
    rp->setEnergy(e);

    // FIXME verify charge calculation
    double bz = marlin::Global::GEAR->getBField().at(gear::Vector3D(0., 0., 0.)).z();
    int signbz = (bz > 0 ? 1 : -1);
    int signom = (om > 0 ? 1 : -1);

    double chrg = signbz * signom;
    rp->setCharge(chrg);

    outCol->addElement(rp);
  }

  evt->addCollection(outCol, _outputPFOCollectionName.c_str());
}

void TrackToRecoParticleConverter::check(LCEvent* /*evt*/) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void TrackToRecoParticleConverter::end() {}
