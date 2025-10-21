#include "photonCorrectionProcessor.h"
#include "lcio.h"
#include <iostream>
#include <marlin/Global.h>

#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"
#include "DDRec/DetectorData.h"

#ifdef MARLIN_USE_AIDA
#include <AIDA/AIDA.h>
#include <marlin/AIDAProcessor.h>
#endif

using namespace lcio;
using namespace marlin;
using std::cout;
using std::endl;

photonCorrectionProcessor aphotonCorrectionProcessor;

photonCorrectionProcessor::photonCorrectionProcessor() : Processor("photonCorrectionProcessor") {
  // processor description
  _description = "photonCorrectionProcessor applies an energy correction to photon-like PFOs";

  registerProcessorParameter("inputCollection", "name of input PFO collection", _inputCollection,
                             std::string("PandoraPFOs"));

  registerProcessorParameter("modifyPFOenergies", "apply the corrected energies to the PFOs", _modifyPFOenergies, true);
  registerProcessorParameter("modifyPFOdirection", "apply the corrected direction to the PFOs", _modifyPFOdirections,
                             true);

  std::vector<float> linearise;
  linearise.push_back(9.870e-01);
  linearise.push_back(1.426e-02);
  registerProcessorParameter("energyCor_Linearise", "parameters to linearise overall energy response",
                             _energyCorr_linearise, linearise);

  std::vector<float> eCorr_barPhi;
  eCorr_barPhi.push_back(0.412249);
  eCorr_barPhi.push_back(0.0142289);
  eCorr_barPhi.push_back(-0.0933687);
  eCorr_barPhi.push_back(0.01345);
  eCorr_barPhi.push_back(0.0408156);
  registerProcessorParameter("energyCorr_barrelPhi", "paramters to correct energy response vs. phi in barrel",
                             _energyCorr_barrelPhi, eCorr_barPhi);

  std::vector<float> eCorr_costh;
  eCorr_costh.push_back(-0.0900);
  eCorr_costh.push_back(0.);
  eCorr_costh.push_back(0.235);
  eCorr_costh.push_back(0.007256);
  eCorr_costh.push_back(-0.0369648);
  eCorr_costh.push_back(0.);
  eCorr_costh.push_back(0.588);
  eCorr_costh.push_back(0.0121604);
  eCorr_costh.push_back(-0.0422968);
  eCorr_costh.push_back(0.774);
  eCorr_costh.push_back(0.009);
  eCorr_costh.push_back(1.002);
  registerProcessorParameter("energyCorr_costh", "paramters to correct energy response vs. cos(theta)",
                             _energyCorr_costheta, eCorr_costh);

  std::vector<float> eCorr_endcap;
  eCorr_endcap.push_back(-0.025);
  eCorr_endcap.push_back(855.);
  eCorr_endcap.push_back(23.);
  eCorr_endcap.push_back(-0.07);
  eCorr_endcap.push_back(1489.);
  eCorr_endcap.push_back(18.);
  registerProcessorParameter("energyCorr_endcap", "paramters to correct energy response vs. endcap cracks",
                             _energyCorr_endcap, eCorr_endcap);

  std::vector<float> phiCorr_barrel;
  phiCorr_barrel.push_back(2.36517e-05);
  phiCorr_barrel.push_back(1.32090e-04);
  phiCorr_barrel.push_back(-3.86883e+00);
  phiCorr_barrel.push_back(-1.67809e-01);
  phiCorr_barrel.push_back(2.28614e-05);
  phiCorr_barrel.push_back(6.03495e-05);
  phiCorr_barrel.push_back(0.419);
  phiCorr_barrel.push_back(0.00728);
  phiCorr_barrel.push_back(0.025);
  phiCorr_barrel.push_back(0.00);
  phiCorr_barrel.push_back(2.86667e-05);
  phiCorr_barrel.push_back(2.49371e-05);
  phiCorr_barrel.push_back(-7.71684e-06);
  phiCorr_barrel.push_back(-1.48118e-05);
  phiCorr_barrel.push_back(-5.63786e-06);
  phiCorr_barrel.push_back(-9.38376e-06);
  phiCorr_barrel.push_back(-4.96296e-06);
  phiCorr_barrel.push_back(2.91262e-06);
  registerProcessorParameter("phiCorr_barrel", "paramters to correct phi bias in barrel", _phiCorr_barrel,
                             phiCorr_barrel);

  std::vector<float> thetaCorr_barrel;
  thetaCorr_barrel.push_back(-0.000166568);
  thetaCorr_barrel.push_back(-7.119e-05);
  thetaCorr_barrel.push_back(0.000223618);
  thetaCorr_barrel.push_back(-3.95915e-05);
  registerProcessorParameter("thetaCorr_barrel", "paramters to correct theta bias in barrel", _thetaCorr_barrel,
                             thetaCorr_barrel);

  std::vector<float> thetaCorr_endcap;
  thetaCorr_endcap.push_back(0.000129478);
  thetaCorr_endcap.push_back(-3.73863e-05);
  thetaCorr_endcap.push_back(-0.000847783);
  thetaCorr_endcap.push_back(0.000153646);
  thetaCorr_endcap.push_back(0.000806605);
  thetaCorr_endcap.push_back(-0.000132608);
  registerProcessorParameter("thetaCorr_endcap", "paramters to correct theta bias in endcap", _thetaCorr_endcap,
                             thetaCorr_endcap);

  registerProcessorParameter("validationPlots", "produce validation plots", _validationPlots, false);
  registerProcessorParameter("nominalEnergy", "nominal photon energy (for validation plots)", _nominalEnergy,
                             float(200));
}

void photonCorrectionProcessor::init() {
  streamlog_out(MESSAGE) << "hello from photonCorrectionProcessor::init" << endl;

  printParameters();

  // first get some geometry information
  dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
  // endcap ecal
  unsigned int includeFlag =
      (dd4hep::DetType::CALORIMETER | dd4hep::DetType::ENDCAP | dd4hep::DetType::ELECTROMAGNETIC);
  unsigned int excludeFlag = (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD);
  const std::vector<dd4hep::DetElement>& endcapEcalDetectors =
      dd4hep::DetectorSelector(mainDetector).detectors(includeFlag, excludeFlag);
  if (endcapEcalDetectors.size() == 1) {
    _assumed_boxsize =
        endcapEcalDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>()->extent[0] / dd4hep::mm; // r-min
    _assumed_endZ =
        endcapEcalDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>()->extent[2] / dd4hep::mm; // z-min
  } else {
    streamlog_out(ERROR) << "did not find exactly one endcap ECAL! found " << endcapEcalDetectors.size()
                         << "; refusing to continue!" << endl;
    assert(0);
  }

  // barrel ecal
  includeFlag = (dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL | dd4hep::DetType::ELECTROMAGNETIC);
  excludeFlag = (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD);
  const std::vector<dd4hep::DetElement>& barrelEcalDetectors =
      dd4hep::DetectorSelector(mainDetector).detectors(includeFlag, excludeFlag);
  if (barrelEcalDetectors.size() == 1) {
    float barrelLength =
        barrelEcalDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>()->extent[3] / dd4hep::mm; // z-max
    float barrelInRad =
        barrelEcalDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>()->extent[0] / dd4hep::mm; // r-min
    _barrelendcap_limit_costh = barrelLength / sqrt(pow(barrelLength, 2) + pow(barrelInRad, 2));
  } else {
    streamlog_out(ERROR) << "did not find exactly one barrel ECAL! found " << barrelEcalDetectors.size()
                         << "; refusing to continue!" << endl;
    assert(0);
  }

  streamlog_out(MESSAGE) << "*** ECAL endcap dimensions: minZ, minR, barrel/endcap transition " << _assumed_endZ
                         << " , " << _assumed_boxsize << " , " << _barrelendcap_limit_costh << endl;

  _photonCorrector = new photonCorrector();

  _photonCorrector->set_barrelendcap_limit(_barrelendcap_limit_costh);
  _photonCorrector->set_assumed_boxsize(_assumed_boxsize);
  _photonCorrector->set_assumed_endZ(_assumed_endZ);

  _photonCorrector->set_energyCorr_linearise(_energyCorr_linearise);
  _photonCorrector->set_energyCorr_barrelPhi(_energyCorr_barrelPhi);
  _photonCorrector->set_energyCorr_costheta(_energyCorr_costheta);
  _photonCorrector->set_energyCorr_endcap(_energyCorr_endcap);
  _photonCorrector->set_phiCorr_barrel(_phiCorr_barrel);
  _photonCorrector->set_thetaCorr_barrel(_thetaCorr_barrel);
  _photonCorrector->set_thetaCorr_endcap(_thetaCorr_endcap);

  return;
}

void photonCorrectionProcessor::processRunHeader(LCRunHeader* /* run */) {
  streamlog_out(MESSAGE) << "hello from photonCorrectionProcessor::processRunHeader" << endl;
}

void photonCorrectionProcessor::processEvent(LCEvent* evt) {
  // adjust energy of all type==22 objects in the input collection

  if (!_modifyPFOdirections && !_modifyPFOenergies) {
    streamlog_out(DEBUG) << "not asked to modify any PFO properties: so doing nothing" << std::endl;
    return;
  }

  // for validation plots: sum the pfo energies before and after correction
  float totPfoEn[2] = {0};
  float totPfoGammaEn[2] = {0};
  float totomom[3] = {0};

  try {
    LCCollectionVec* pfocol = dynamic_cast<LCCollectionVec*>(evt->getCollection(_inputCollection));
    for (auto lcobj : (*pfocol)) {
      ReconstructedParticleImpl* pfo = dynamic_cast<ReconstructedParticleImpl*>(lcobj);

      if (_validationPlots) {
        for (int j = 0; j < 3; j++) {
          totomom[j] += pfo->getMomentum()[j];
        }
      }

      if (pfo->getType() == 22) {

        float origEn = pfo->getEnergy();
        float origMom(0);
        for (int i = 0; i < 3; i++) {
          origMom += pow(pfo->getMomentum()[i], 2);
        }
        origMom = sqrt(origMom);
        float origTheta = acos(pfo->getMomentum()[2] / origMom);
        float origPhi = atan2(pfo->getMomentum()[1], pfo->getMomentum()[0]);

        float correctedEnergy = _photonCorrector->photonEnergyCorrection(pfo);

        float correctedMom = origMom * correctedEnergy /
                             origEn; // just scale. in principle should be same as energy for massless photon PFO

        float correctedTheta, correctedPhi;
        _photonCorrector->photonDirectionCorrection(pfo, correctedTheta, correctedPhi);

        // these are the properties we will set for the PFO
        float newEnergy = _modifyPFOenergies ? correctedEnergy : origEn;
        float newMom = _modifyPFOenergies ? correctedMom : origMom;
        float newTheta = _modifyPFOdirections ? correctedTheta : origTheta;
        float newPhi = _modifyPFOdirections ? correctedPhi : origPhi;

        float newMomentum[3];
        newMomentum[0] = newMom * sin(newTheta) * cos(newPhi);
        newMomentum[1] = newMom * sin(newTheta) * sin(newPhi);
        newMomentum[2] = newMom * cos(newTheta);

        streamlog_out(DEBUG) << "updating energy/momentum from " << origEn << " / (" << pfo->getMomentum()[0] << " "
                             << pfo->getMomentum()[1] << " " << pfo->getMomentum()[2] << ") to " << newEnergy << " / ("
                             << newMomentum[0] << " " << newMomentum[1] << " " << newMomentum[2] << ")" << std::endl;

        pfo->setEnergy(newEnergy);
        pfo->setMomentum(newMomentum);

        if (_validationPlots) {
          totPfoEn[0] += pfo->getEnergy();
          totPfoGammaEn[0] += pfo->getEnergy();
          totPfoEn[1] += correctedEnergy;
          totPfoGammaEn[1] += correctedEnergy;
        }

      } else {
        totPfoEn[0] += pfo->getEnergy();
        totPfoEn[1] += pfo->getEnergy();
      }
    }
  } catch (DataNotAvailableException& e) {
  };

#ifdef MARLIN_USE_AIDA

  if (_validationPlots) {
    static AIDA::IHistogram2D* h_sumPfoE_orig;
    static AIDA::IHistogram2D* h_sumPfoE_corr;
    static AIDA::IHistogram2D* h_sumGamPfoE_orig;
    static AIDA::IHistogram2D* h_sumGamPfoE_corr;
    if (isFirstEvent()) {
      h_sumPfoE_orig = AIDAProcessor::histogramFactory(this)->createHistogram2D("sumPfoE_orig", "energy [GeV]", 200, -1,
                                                                                1, 100, 0., 1.5 * _nominalEnergy);
      h_sumPfoE_corr = AIDAProcessor::histogramFactory(this)->createHistogram2D("sumPfoE_corr", "energy [GeV]", 200, -1,
                                                                                1, 100, 0., 1.5 * _nominalEnergy);
      h_sumGamPfoE_orig = AIDAProcessor::histogramFactory(this)->createHistogram2D(
          "sumGamPfoE_orig", "energy [GeV]", 200, -1, 1, 100, 0., 1.5 * _nominalEnergy);
      h_sumGamPfoE_corr = AIDAProcessor::histogramFactory(this)->createHistogram2D(
          "sumGamPfoE_corr", "energy [GeV]", 200, -1, 1, 100, 0., 1.5 * _nominalEnergy);
    }

    float costh = totomom[2] / sqrt(pow(totomom[0], 2) + pow(totomom[1], 2) + pow(totomom[2], 2));

    h_sumPfoE_orig->fill(costh, totPfoEn[0]);
    h_sumPfoE_corr->fill(costh, totPfoEn[1]);
    h_sumGamPfoE_orig->fill(costh, totPfoGammaEn[0]);
    h_sumGamPfoE_corr->fill(costh, totPfoGammaEn[1]);
  }

#endif

  return;
}

void photonCorrectionProcessor::check(LCEvent* /* evt */) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void photonCorrectionProcessor::end() {
  streamlog_out(MESSAGE) << "photonCorrectionProcessor::end()  " << name() << std::endl;
}
