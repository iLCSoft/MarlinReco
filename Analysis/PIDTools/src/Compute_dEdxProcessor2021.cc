#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

#include <marlin/AIDAProcessor.h>
#include <marlin/Global.h>
#include <marlin/Processor.h>

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/TrackImpl.h>
#include <lcio.h>

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitField64.h>

// MarlinUtil
#include "HelixClass.h"
#include "LCGeometryTypes.h"
#include "SimpleHelix.h"

#include "TCanvas.h"
#include "TImage.h"
#include "TStyle.h"

#include "Compute_dEdxProcessor2021.hh"

Compute_dEdxProcessor2021 aCompute_dEdxProcessor2021;

Compute_dEdxProcessor2021::Compute_dEdxProcessor2021() : Processor("Compute_dEdxProcessor2021") {

  // Processor description
  _description = "dE/dx calculation using truncation method, angular correction applied after smearing";

  registerInputCollection(LCIO::TRACK, "LDCTrackCollection", "LDC track collection name", _LDCTrackCollection,
                          std::string("MarlinTrkTracks"));

  registerProcessorParameter(
      "Write_dEdx", "If set, the calculated dEdx value will be attached to its corresponding track (default: true).",
      _writedEdx, bool(true));

  registerProcessorParameter("EnergyLossErrorTPC", "Fractional error of dEdx in the TPC", _energyLossErrorTPC,
                             float(0.054));

  registerProcessorParameter("LowerTruncationFraction", "Lower truncation fraction, default: 0.08 (8%)", _lowerTrunFrac,
                             float(0.08));

  registerProcessorParameter("UpperTruncationFraction", "Upper truncation fraction, default: 0.3 (30%)", _upperTrunFrac,
                             float(0.3));

  registerProcessorParameter("isSmearing", "Flag for dEdx Smearing", _isSmearing, bool(0));

  registerProcessorParameter("smearingFactor", "Smearing factor for dEdx smearing", _smearingFactor, float(0.035));

  std::vector<float> errexp = {-0.34, -0.45};
  registerProcessorParameter("dEdxErrorScalingExponents",
                             "scaling exponents of the dE/dx error for path length and number of hits, respectively",
                             _errexp, errexp);

  registerProcessorParameter("dxStrategy",
                             "ID of the track length calculation: 1 - hit-to-hit distance; 2 - hit-to-hit path length "
                             "of projected hits; 3 - path over hit row",
                             _dxStrategy, int(1));

  registerProcessorParameter("StrategyCompHist",
                             "If set, an AIDA root file will be generated with Bethe-Bloch histograms (dE/dx over "
                             "momentum) for each dx strategy (default: false).",
                             _StratCompHist, bool(false));

  registerProcessorParameter("StrategyCompHistWeight",
                             "If set, the Bethe-Bloch histograms will have a sqrt(nHits) weighting (default: false).",
                             _StratCompHistWeight, bool(false));

  registerProcessorParameter("StrategyCompHistFiles",
                             "File names of the generated dx strategy comparison histograms (if chosen). The "
                             "respective strategy number and '.png' is added.",
                             _StratCompHistFiles, std::string("dEdx_Histo_Strategy"));

  registerProcessorParameter("doAngularCorr_dEdx",
                             "If set, the calculated dEdx value is corrected as function of a pol3(lambda)",
                             _angularcorrdEdx, bool(false));

  std::vector<float> _newpar = {0.948185, 0.000520502, 3.33231e-5, -1.51052e-7};
  registerProcessorParameter("AngularCorrectionParameters",
                             "parameter for new angular correction dedx= uncorrected_dedx  / f, with f= pol3(lambda)",
                             _par, _newpar);
}

void Compute_dEdxProcessor2021::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;

  // it's usually a good idea to
  printParameters();

  // get TPC radii, pad height and B field
  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
  dd4hep::DetElement tpcDE = lcdd.detector("TPC");
  const dd4hep::rec::FixedPadSizeTPCData* tpc = tpcDE.extension<dd4hep::rec::FixedPadSizeTPCData>();
  _TPC_inner = tpc->rMinReadout / dd4hep::mm;
  _TPC_outer = tpc->rMaxReadout / dd4hep::mm;
  _TPC_padHeight = tpc->padHeight / dd4hep::mm;
  double bfieldV[3];
  lcdd.field().magneticField({0., 0., 0.}, bfieldV);
  _bField = bfieldV[2] / dd4hep::tesla;

  engine = new std::default_random_engine();

  streamlog_out(MESSAGE) << "TPC inner radius = " << _TPC_inner << " mm   TPC outer radius = " << _TPC_outer << " mm"
                         << std::endl;
  streamlog_out(MESSAGE) << "TPC pad height = " << _TPC_padHeight << " mm   B field = " << _bField << " T" << std::endl;
  streamlog_out(DEBUG) << "dd4hep::mm = " << dd4hep::mm << "   dd4hep::cm = " << dd4hep::cm << std::endl;

  if (_StratCompHist) {
    // prepare the histogram bins - this could also be implemented as processor parameters
    const int nBinsX = 1000, nBinsY = 300;
    double minBinX = -1, maxBinX = 2; // log10(min/max momentum / GeV)
    double minBinY = 0, maxBinY = 1e-6;
    double histbinsX[nBinsX + 1], histbinsY[nBinsY + 1];
    for (int i = 0; i < nBinsX + 1; i++)
      histbinsX[i] = pow(10, minBinX + (maxBinX - minBinX) * i / (float)(nBinsX));
    for (int i = 0; i < nBinsY + 1; i++)
      histbinsY[i] = minBinY + (maxBinY - minBinY) * i / (float)(nBinsY);

    // create the histograms
    AIDAProcessor::tree(this);
    _BBHist_Strategy1 = new TH2D("BBHist_Strategy1", "Bethe-Bloch curve for dx strategy 1: hit-to-hit distance", nBinsX,
                                 histbinsX, nBinsY, histbinsY);
    _BBHist_Strategy2 =
        new TH2D("BBHist_Strategy2", "Bethe-Bloch curve for dx strategy 2: hit-to-hit path length of projected hits",
                 nBinsX, histbinsX, nBinsY, histbinsY);
    _BBHist_Strategy3 = new TH2D("BBHist_Strategy3", "Bethe-Bloch curve for dx strategy 3: path over hit row", nBinsX,
                                 histbinsX, nBinsY, histbinsY);
  }
}

void Compute_dEdxProcessor2021::processRunHeader(LCRunHeader*) {}

void Compute_dEdxProcessor2021::processEvent(LCEvent* evt) {
  _LDCCol = evt->getCollection(_LDCTrackCollection);
  int nTrkCand = _LDCCol->getNumberOfElements();

  for (int iTRK = 0; iTRK < nTrkCand; ++iTRK) {

    TrackImpl* trkCand = (TrackImpl*)_LDCCol->getElementAt(iTRK);
    TrackerHitVec trkHits = trkCand->getTrackerHits();

    // calculate dEdx here when the subtrack is the TPC track!
    // silicon dE/dx can also be calculated. is it necessary??
    float dedx = 0.0;
    float unum = 0.0;

    // get tpc hits only
    TrackerHitVec tpcHitVec;
    for (unsigned int sdet = 0; sdet < trkHits.size(); sdet++) {
      if (trkHits[sdet]->getType() == 0) { // it is very suspicious condition...
        // check whether this hit is in TPC
        float x = trkHits[sdet]->getPosition()[0];
        float y = trkHits[sdet]->getPosition()[1];
        if (sqrt(x * x + y * y) >= _TPC_inner && sqrt(x * x + y * y) <= _TPC_outer)
          tpcHitVec.push_back(trkHits[sdet]); // add hit
        else
          streamlog_out(DEBUG) << "hit not in TPC: " << x << " " << y << " " << trkHits[sdet]->getPosition()[2] << " "
                               << sqrt(x * x + y * y) << std::endl;
      }
    }

    std::pair<double, double> dummy = CalculateEnergyLoss(tpcHitVec, trkCand);
    dedx = dummy.first;
    unum = dummy.second * _energyLossErrorTPC;

    if (_angularcorrdEdx == true) {
      float trklambda = trkCand->getTanLambda();
      dedx = get_Corrected_dEdx(dedx, trklambda);
    }

    // fill values
    if (_writedEdx) {
      trkCand->setdEdx(dedx);
      trkCand->setdEdxError(unum);
    }
  }
}

void Compute_dEdxProcessor2021::check(LCEvent*) {}

void Compute_dEdxProcessor2021::end() {

  if (_StratCompHist) {
    TCanvas* can = new TCanvas;
    TImage* img = TImage::Create();
    gStyle->SetPalette(kBird);
    can->SetGrid();
    can->SetLogx();
    _BBHist_Strategy1->SetXTitle("momentum / GeV");
    _BBHist_Strategy1->SetYTitle("dE/dx / (GeV/mm)");
    _BBHist_Strategy2->SetXTitle("momentum / GeV");
    _BBHist_Strategy2->SetYTitle("dE/dx / (GeV/mm)");
    _BBHist_Strategy3->SetXTitle("momentum / GeV");
    _BBHist_Strategy3->SetYTitle("dE/dx / (GeV/mm)");

    std::string StratCompHistFiles2 = _StratCompHistFiles;
    _BBHist_Strategy1->Draw("colz");
    img->FromPad(can);
    _StratCompHistFiles.append("1.png");
    img->WriteImage(_StratCompHistFiles.c_str());
    _BBHist_Strategy2->Draw("colz");
    img->FromPad(can);
    _StratCompHistFiles.replace(_StratCompHistFiles.end() - 5, _StratCompHistFiles.end() - 4, "2");
    img->WriteImage(_StratCompHistFiles.c_str());
    _BBHist_Strategy3->Draw("colz");
    img->FromPad(can);
    _StratCompHistFiles.replace(_StratCompHistFiles.end() - 5, _StratCompHistFiles.end() - 4, "3");
    img->WriteImage(_StratCompHistFiles.c_str());

    delete can;
    delete img;
  }

  delete engine;
}

std::pair<double, double> Compute_dEdxProcessor2021::CalculateEnergyLoss(TrackerHitVec& hitVec, Track* trk) {
  int tH = hitVec.size();

  // fill the track into a SimpleHelix and also into a HelixClass to make use of their specific geometric utilities
  const float* refPoint = trk->getReferencePoint();
  LCVector3D refPointV3D =
      LCVector3D(refPoint[0], refPoint[1], refPoint[2]); // MarlinUtil::LCVector3D is CLHEP::Hep3Vector
  SimpleHelix trkHelix =
      SimpleHelix(trk->getD0(), trk->getPhi(), trk->getOmega(), trk->getZ0(), trk->getTanLambda(), refPointV3D);
  HelixClass trkHelixC;
  trkHelixC.Initialize_Canonical(trk->getPhi(), trk->getD0(), trk->getZ0(), trk->getOmega(), trk->getTanLambda(),
                                 _bField);

  // set up some variables
  std::vector<double> dEdxVec1, dEdxVec2, dEdxVec3; // one vector per dx strategy
  double dx = 0, dy = 0, dz = 0, length = 0, totlength = 0;
  double dxy = 0, r = 0, theta = 0, l = 0;
  double path = 0, totpath = 0, cpath = 0;
  double over = 0, totover = 0, cover = 0;
  double hitX = 0, hitY = 0, hitZ = 0;
  double hit0X = 0, hit0Y = 0, hit0Z = 0;

  for (int iH = 0; iH < tH; iH++) {
    TrackerHit* hit = hitVec[iH];
    hitX = hit->getPosition()[0];
    hitY = hit->getPosition()[1];
    hitZ = hit->getPosition()[2];

    if (iH > 0) {
      if (_dxStrategy == 1 || _dxStrategy == 2 || _StratCompHist) {
        TrackerHit* hit0 = hitVec[iH - 1];
        hit0X = hit0->getPosition()[0];
        hit0Y = hit0->getPosition()[1];
        hit0Z = hit0->getPosition()[2];

        if (_dxStrategy == 1 || _StratCompHist) { // strategy 1: dx = distance between 2 hits
          dx = hitX - hit0X;
          dy = hitY - hit0Y;
          dz = hitZ - hit0Z;

          // cal. length of curvature
          dxy = sqrt(dx * dx + dy * dy);
          r = 1.0 / trk->getOmega();
          theta = 2.0 * asin(0.5 * dxy / r);
          l = r * theta;
          // total helix length
          length = sqrt(l * l + dz * dz);

          // sum of flight length
          totlength += length;
          dEdxVec1.push_back(hit->getEDep() / length);
        }

        if (_dxStrategy == 2 || _StratCompHist) { // strategy 2: dx = path difference on the track via SimpleHelix
          double path1 = trkHelix.getPathAt(LCVector3D(hitX, hitY, hitZ));
          double path2 = trkHelix.getPathAt(LCVector3D(hit0X, hit0Y, hit0Z));
          path = fabs(path2 - path1);
          totpath += path;
          dEdxVec2.push_back(hit->getEDep() / path);
          cpath++;
        }
      }
    }

    if (_dxStrategy == 3 ||
        _StratCompHist) { // strategy 3: dx = path of track over the hit's row via SimpleHelix and HelixClass
      // since this strategy does not require a distance between two hits, it can be used for every single hit
      double hitRadius = sqrt(hitX * hitX + hitY * hitY);

      // create cylinders of the row's inner and outer edge, check if the track crosses both, then get the path
      // difference of these crossing points
      double cyl1_radius =
          hitRadius - _TPC_padHeight * .5; // assuming the radial hit position is always in the center of the row
      double cyl2_radius = hitRadius + _TPC_padHeight * .5;
      float trkRefP[3];
      trkRefP[0] = (trkHelixC.getReferencePoint())[0];
      trkRefP[1] = (trkHelixC.getReferencePoint())[1];
      trkRefP[2] = (trkHelixC.getReferencePoint())[2];
      float cross1[6], cross2[6];
      float time1 = trkHelixC.getPointOnCircle(cyl1_radius, trkRefP, cross1);
      float time2 = trkHelixC.getPointOnCircle(cyl2_radius, trkRefP, cross2);

      bool rejected = false;
      if (fabs(time1) < 1e+10 &&
          fabs(time2) < 1e+10) { // getPointOnCircle returns -1e+20 when the helix doesn't cross the cylinder
        if (sqrt(pow(cross1[0] - cross1[3], 2) + pow(cross1[1] - cross1[4], 2) + pow(cross1[2] - cross1[5], 2)) <
            2 * _TPC_padHeight)
          rejected = true; // curler protection
        if (sqrt(pow(cross2[0] - cross2[3], 2) + pow(cross2[1] - cross2[4], 2) + pow(cross2[2] - cross2[5], 2)) <
            2 * _TPC_padHeight)
          rejected = true;
      } else
        rejected = true;

      if (!rejected) {
        float* cross_1 = cross1; // check which crossing point is the correct one by comparing to the hit position
        double dist1 = sqrt(pow(cross1[0] - hitX, 2) + pow(cross1[1] - hitY, 2) + pow(cross1[2] - hitZ, 2));
        double dist2 = sqrt(pow(cross1[3] - hitX, 2) + pow(cross1[4] - hitY, 2) + pow(cross1[5] - hitZ, 2));
        if (dist2 < dist1)
          cross_1 = &cross1[3];
        float* cross_2 = cross2;
        dist1 = sqrt(pow(cross2[0] - hitX, 2) + pow(cross2[1] - hitY, 2) + pow(cross2[2] - hitZ, 2));
        dist2 = sqrt(pow(cross2[3] - hitX, 2) + pow(cross2[4] - hitY, 2) + pow(cross2[5] - hitZ, 2));
        if (dist2 < dist1)
          cross_2 = &cross2[3];

        // again correct for curvature
        dxy = sqrt(pow(cross_1[0] - cross_2[0], 2) + pow(cross_1[1] - cross_2[1], 2));
        r = 1. / trk->getOmega();
        l = r * 2. * asin(0.5 * dxy / r);
        over = sqrt(l * l + pow(cross_1[2] - cross_2[2], 2));

        totover += over;
        cover++;
        dEdxVec3.push_back(hit->getEDep() / over);
      }
    }

    streamlog_out(DEBUG3) << "1: " << length << "  2: " << path << "  3: " << over << "  EDep: " << hit->getEDep()
                          << std::endl;
    if (length / over > 1.5 || over / length > 1.5)
      streamlog_out(DEBUG2) << "  " << hitX << " " << hitY << " " << hitZ << " " << sqrt(hitX * hitX + hitY * hitY)
                            << std::endl;
  }
  streamlog_out(DEBUG4) << "1: " << totlength << "  2: " << totpath << " of " << cpath << "  3: " << totover << " of "
                        << cover << std::endl;

  int nTruncate = 0; // number of remaining track hits after truncation
  double normdedx = 0, totway = 0;
  // double trkcos = sqrt(trk->getTanLambda()*trk->getTanLambda()/(1.0+trk->getTanLambda()*trk->getTanLambda()));

  // if comparison plots should be made, the dE/dx for all three strategies is calculated, filled in histograms,
  // and the one of the chosen strategy is also transferred to be returned
  if (_StratCompHist) {
    double dedx1 = 0, dedx2 = 0, dedx3 = 0;
    int nTruncate1 = 0, nTruncate2 = 0, nTruncate3 = 0;

    std::sort(dEdxVec1.begin(), dEdxVec1.end());
    std::sort(dEdxVec2.begin(), dEdxVec2.end());
    std::sort(dEdxVec3.begin(), dEdxVec3.end());

    for (int idE = dEdxVec1.size() * _lowerTrunFrac; idE < dEdxVec1.size() * (1 - _upperTrunFrac); idE++) {
      dedx1 += dEdxVec1[idE];
      nTruncate1++;
    }
    for (int idE = dEdxVec2.size() * _lowerTrunFrac; idE < dEdxVec2.size() * (1 - _upperTrunFrac); idE++) {
      dedx2 += dEdxVec2[idE];
      nTruncate2++;
    }
    for (int idE = dEdxVec3.size() * _lowerTrunFrac; idE < dEdxVec3.size() * (1 - _upperTrunFrac); idE++) {
      dedx3 += dEdxVec3[idE];
      nTruncate3++;
    }

    // cal. truncated mean
    dedx1 = dedx1 / (float)(nTruncate1);
    dedx2 = dedx2 / (float)(nTruncate2);
    dedx3 = dedx3 / (float)(nTruncate3);
    streamlog_out(DEBUG4) << "  check dedx: " << dedx1 << " " << dedx2 << " " << dedx3 << std::endl;

    double normdedx1 = dedx1; // getNormalization(dedx1, (float)nTruncate1, trkcos);
    double normdedx2 = dedx2; // getNormalization(dedx2, (float)nTruncate2, trkcos);
    double normdedx3 = dedx3; // getNormalization(dedx3, (float)nTruncate3, trkcos);
    if (_isSmearing) {
      normdedx1 += getSmearing(normdedx1);
      normdedx2 += getSmearing(normdedx2);
      normdedx3 += getSmearing(normdedx3);
    }

    double lmom = sqrt(pow(trkHelixC.getMomentum()[0], 2) + pow(trkHelixC.getMomentum()[1], 2) +
                       pow(trkHelixC.getMomentum()[2], 2));

    if (_StratCompHistWeight) {
      _BBHist_Strategy1->Fill(lmom, normdedx1, sqrt(nTruncate1));
      _BBHist_Strategy2->Fill(lmom, normdedx2, sqrt(nTruncate2));
      _BBHist_Strategy3->Fill(lmom, normdedx3, sqrt(nTruncate3));
    } else {
      _BBHist_Strategy1->Fill(lmom, normdedx1);
      _BBHist_Strategy2->Fill(lmom, normdedx2);
      _BBHist_Strategy3->Fill(lmom, normdedx3);
    }

    switch (_dxStrategy) {
    case 1:
      normdedx = normdedx1;
      nTruncate = nTruncate1;
      totway = totlength;
      break;
    case 2:
      normdedx = normdedx2;
      nTruncate = nTruncate2;
      totway = totpath;
      break;
    case 3:
      normdedx = normdedx3;
      nTruncate = nTruncate3;
      totway = totover;
      break;
    }
  } else // to if (_StratCompHist)
  {
    double dedx = 0;
    std::vector<double> dEdxVec;
    switch (_dxStrategy) {
    case 1:
      dEdxVec = dEdxVec1;
      totway = totlength;
      break;
    case 2:
      dEdxVec = dEdxVec2;
      totway = totpath;
      break;
    case 3:
      dEdxVec = dEdxVec3;
      totway = totover;
      break;
    }

    std::sort(dEdxVec.begin(), dEdxVec.end());
    for (int idE = dEdxVec.size() * _lowerTrunFrac; idE < dEdxVec.size() * (1 - _upperTrunFrac); idE++) {
      dedx += dEdxVec[idE];
      nTruncate++;
    }
    dedx = dedx / (float)(nTruncate);
    normdedx = dedx; // getNormalization(dedx, (float)nTruncate, trkcos);
    if (_isSmearing)
      normdedx += getSmearing(normdedx);
  }

  std::pair<double, double> ret;
  if (std::isnan(normdedx)) {
    ret.first = 0;
    ret.second = 0;
    streamlog_out(DEBUG4) << "NAN!" << std::endl;
  } else {
    ret.first = normdedx;
    ret.second = normdedx * std::pow(0.001 * totway, _errexp[0]) *
                 std::pow(nTruncate, _errexp[1]); // todo: what is gas pressure in TPC?
  }

  return ret;
}

// create smearing factor
double Compute_dEdxProcessor2021::getSmearing(double dEdx) {
  double z = sqrt(-2.0 * std::log(dist(*engine))) * std::sin(2.0 * M_PI * dist(*engine)); // 3.141592

  // std::cout << dist(*engine) << " " << z << std::endl;

  return 0.0 + dEdx * _smearingFactor * z;
}

double Compute_dEdxProcessor2021::get_Corrected_dEdx(float dedx, float trklambda) {

  // // DEFAULT CORRECTION APPLIED FOR THE MC2020 production obtained using versions v02-02-01 and previous
  // //cal. hit dep.
  // double f1 = 1 + std::exp(-nHit/_ncorrpar);
  // //cal. polar angle dep.
  // // double c=1.0/sqrt(1.0-trkcos*trkcos);
  // // double f2=1.0/(1.0-0.08887*std::log(c));

  // //cal. polar angle dep.   20160702
  // //double c = std::acos(trkcos);
  // //if(c>3.141592/2.0) c= 3.141592-c;
  // //double f2 = 1.0/std::pow(c, 0.0703);

  // //change polar angle dep. 20170901
  // double f2 = _acorrpar[0] / (_acorrpar[0] + _acorrpar[1] * trkcos * trkcos);
  // return dedx/(f1*f2);

  // ****************************************
  // NEW CORRECTION
  // Parametrization by A. Irles, U. Einhaus  2021/04
  // Using single particle mc2020 samples v02-02-01
  // Analyzed with dEdxAnalyser, fitting histogram NormLambdaFullAll_1
  float lambda = fabs(atan(trklambda) * 180. / M_PI);
  double f3 = 1 / (_par[0] + _par[1] * lambda + _par[2] * pow(lambda, 2) + _par[3] * pow(lambda, 3));

  return dedx * f3;
}
