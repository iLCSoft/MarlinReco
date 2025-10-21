#include "AngularCorrection_dEdxProcessor.hh"

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

using namespace lcio;
using namespace marlin;

AngularCorrection_dEdxProcessor aAngularCorrection_dEdxProcessor;

AngularCorrection_dEdxProcessor::AngularCorrection_dEdxProcessor() : Processor("AngularCorrection_dEdxProcessor") {

  // Processor description
  _description = "Correct_AngularCorrection_dEdxProcessor: Makes a hard angular-based correction of dEdx for all the "
                 "Tracks in the event. ATTENTION: this processor rewrites the MarlinTrk Collection and it is to be "
                 "used only for simulations produced with ILCSoft from v02 to v02-02-01";

  registerInputCollection(LCIO::TRACK, "LDCTrackCollection", "LDC track collection name", _LDCTrackCollection,
                          std::string("MarlinTrkTracks"));

  std::vector<float> _newpar = {0.970205, 0.0007506, 4.41781e-8, 5.8222e-8};

  // ***************************************
  // Fit of NormLamdaFullAll_1 (dEdxAnalyser) using single particle samples recosntructed with v02-02-01
  // 2021/04
  // (including the default previous angular correction)
  // Minimizer is Linear
  // Chi2                      =      807.398
  // NDf                       =           27
  // p0                        =     0.970205   +/-   0.000127468
  // p1                        =   0.00075065   +/-   1.55853e-05
  // p2                        =  4.41781e-08   +/-   5.09662e-07
  // p3                        =   5.8222e-08   +/-   4.71913e-09
  //
  registerProcessorParameter("AngularCorrectionParameters",
                             "parameter for new angular correction dedx= uncorrected_dedx  / f, with f= pol3(lambda)",
                             _par, _newpar);
}

void AngularCorrection_dEdxProcessor::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;
  // it's usually a good idea to
  printParameters();
}

void AngularCorrection_dEdxProcessor::processRunHeader(LCRunHeader*) {}

void AngularCorrection_dEdxProcessor::modifyEvent(LCEvent* evt) {

  // fill values
  _LDCCol = evt->getCollection(_LDCTrackCollection);
  int nTrkCand = _LDCCol->getNumberOfElements();

  for (int iTRK = 0; iTRK < nTrkCand; ++iTRK) {

    TrackImpl* trkCand = (TrackImpl*)_LDCCol->getElementAt(iTRK);

    float dedx = trkCand->getdEdx();
    float dedx_error = trkCand->getdEdxError();
    float trklambda = trkCand->getTanLambda();

    float lambda = fabs(atan(trklambda) * 180. / M_PI);
    double f3 = 1 / (_par[0] + _par[1] * lambda + _par[2] * pow(lambda, 2) + _par[3] * pow(lambda, 3));

    double new_dedx = dedx * f3;

    streamlog_out(DEBUG) << "Original dEdx: " << dedx << " Error: " << dedx_error << std::endl;
    streamlog_out(DEBUG) << "NeW dEdx: " << new_dedx << " Error: " << dedx_error << std::endl;

    // fill values
    trkCand->setdEdx(new_dedx);
    trkCand->setdEdxError(dedx_error);
  }
}

void AngularCorrection_dEdxProcessor::check(LCEvent*) {}

void AngularCorrection_dEdxProcessor::end() {}
