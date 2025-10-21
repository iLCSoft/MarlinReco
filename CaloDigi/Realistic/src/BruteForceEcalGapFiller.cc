#include "BruteForceEcalGapFiller.h"

#include <IMPL/CalorimeterHitImpl.h>

#include <CalorimeterHitType.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/LCCollection.h>
#include <UTIL/CellIDDecoder.h>

#include "DD4hep/DetType.h"

#include <iostream>
using std::cout;
using std::endl;

#include <cassert>
#include <math.h>

#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"

BruteForceEcalGapFiller aBruteForceEcalGapFiller;

BruteForceEcalGapFiller::BruteForceEcalGapFiller() : Processor("BruteForceEcalGapFiller") {

  _description = "makes a collection of ECAL gap hits";

  // input collection of calohits
  std::string inputCollection;
  registerInputCollection(LCIO::CALORIMETERHIT, "inputHitCollection", "input simcalhit Collection Name",
                          _inputHitCollection, inputCollection);

  std::string outputCollection;
  registerOutputCollection(LCIO::CALORIMETERHIT, "outputHitCollection", "output calorimeterhit Collection Name",
                           _outputHitCollection, outputCollection);

  registerProcessorParameter("CellIDLayerString", "name of the part of the cellID that holds the layer",
                             _cellIDLayerString, std::string("layer"));

  registerProcessorParameter("CellIDModuleString", "name of the part of the cellID that holds the module",
                             _cellIDModuleString, std::string("module"));

  registerProcessorParameter("CellIDStaveString", "name of the part of the cellID that holds the stave",
                             _cellIDStaveString, std::string("stave"));

  registerProcessorParameter(
      "expectedInterModuleDistance",
      "size of gap across module boundaries (from edge to edge of cells, in mm ; accuracy < cell size)",
      _interModuleDist, float(7.));

  registerProcessorParameter("interModuleNonlinearFactor",
                             "nonlin factor f: E_corr = interModuleCorrectionFactor*(1/f)*log(1 + f*E_calc)",
                             _interModuleNonlinearFactor, float(1.));

  registerProcessorParameter("intraModuleNonlinearFactor",
                             "nonlin factor f: E_corr = intraModuleCorrectionFactor*(1/f)*log(1 + f*E_calc)",
                             _intraModuleNonlinearFactor, float(1.));

  registerProcessorParameter("interModuleCorrectionFactor",
                             "factor applied to calculated energy of inter-module gap hits", _interModuleFactor,
                             float(0.35));

  registerProcessorParameter("intraModuleCorrectionFactor",
                             "factor applied to calculated energy of intra-module gap hits", _intraModuleFactor,
                             float(1.0));

  registerProcessorParameter("applyInterModuleCorrection", "apply correction for gaps between modules?",
                             _applyInterModuleCor, bool(true));
}

void BruteForceEcalGapFiller::init() {

  printParameters();

  _flag.setBit(LCIO::CHBIT_LONG);
  _flag.setBit(LCIO::RCHBIT_TIME); // store timing on output hits.

  _currentLayout = -99;
}

void BruteForceEcalGapFiller::processEvent(LCEvent* evt) {

  streamlog_out(DEBUG1) << "looking for collection: " << _inputHitCollection << endl;
  try {
    LCCollection* col = evt->getCollection(_inputHitCollection.c_str());

    int numElements = col->getNumberOfElements();
    streamlog_out(DEBUG1) << _inputHitCollection << " number of elements = " << numElements << endl;

    if (numElements > 0) {

      CellIDDecoder<CalorimeterHit> idDecoder(col);

      for (int il = 0; il < MAXLAYER; il++)
        for (int is = 0; is < MAXSTAVE; is++)
          for (int im = 0; im < MAXMODULE; im++)
            hitsByLayerModuleStave[il][is][im].clear();

      // loop over input hits
      for (int j(0); j < numElements; ++j) {
        CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(j));

        if (j == 0)
          getGeometryData(hit->getType()); // update geom info for first hit in collection (assumes that single
                                           // collection doesn't mix endcap and barrel hits)

        int layer = idDecoder(hit)[_cellIDLayerString];
        int stave = idDecoder(hit)[_cellIDStaveString];
        int module = idDecoder(hit)[_cellIDModuleString];
        assert(layer >= 0 && layer < MAXLAYER);
        assert(stave >= 0 && stave < MAXSTAVE);
        assert(module >= 0 && module < MAXMODULE);
        hitsByLayerModuleStave[layer][stave][module].push_back(hit);
      }

      // now make the gap hits

      // create new collection: hits
      std::string encodingString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      LCCollectionVec* newcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
      newcol->parameters().setValue(LCIO::CellIDEncoding, encodingString);
      newcol->setFlag(_flag.getFlag());

      addIntraModuleGapHits(newcol); // gaps within a module
      if (_applyInterModuleCor)
        addInterModuleGapHits(newcol); // gaps between modules

      evt->addCollection(newcol, _outputHitCollection);
    }

  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG1) << "could not find input collection " << _inputHitCollection << std::endl;
  }
}

void BruteForceEcalGapFiller::getGeometryData(const int ihitType) {
  // get information about geometry
  // calorimeter hit type used to decide if it's in barrel or endcap

  CHT calHitType(ihitType);
  if (!calHitType.is(CHT::ecal)) {
    streamlog_out(ERROR) << "this is not an ECAL hit!" << endl;
    assert(0);
  }

  int iLayout(-99);
  if (calHitType.is(CHT::barrel))
    iLayout = ECALBARREL;
  else if (calHitType.is(CHT::endcap))
    iLayout = ECALENDCAP;
  else {
    streamlog_out(ERROR) << "this hit is in neither barrel not endcap" << endl;
    assert(0);
  }

  if (iLayout != _currentLayout) { // layout (barrel/endcap region) has changed, get appropriate geom data
    _currentLayout = iLayout;

    unsigned int includeFlag(0);
    unsigned int excludeFlag(0);

    if (calHitType.is(CHT::barrel)) {
      includeFlag = (dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL);
      excludeFlag = (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD);
    } else if (calHitType.is(CHT::endcap)) {
      includeFlag = (dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP);
      excludeFlag = (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD);
    }

    dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
    const std::vector<dd4hep::DetElement>& theDetectors =
        dd4hep::DetectorSelector(lcdd).detectors(includeFlag, excludeFlag);
    streamlog_out(DEBUG2) << " getExtension :  includeFlag: " << dd4hep::DetType(includeFlag)
                          << " excludeFlag: " << dd4hep::DetType(excludeFlag) << "  found : " << theDetectors.size()
                          << "  - first det: " << theDetectors.at(0).name() << std::endl;

    if (theDetectors.size() != 1) {
      std::stringstream es;
      streamlog_out(ERROR) << " getExtension: selection is not unique (or empty)  includeFlag: "
                           << dd4hep::DetType(includeFlag) << " excludeFlag: " << dd4hep::DetType(excludeFlag)
                           << " --- found detectors : ";
      for (unsigned i = 0, N = theDetectors.size(); i < N; ++i) {
        streamlog_out(ERROR) << theDetectors.at(i).name() << ", ";
      }
      assert(0);
    }
    _caloGeomData = theDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

    if (!_caloGeomData) {
      streamlog_out(WARNING) << "could not get calorimeter geometry information!" << endl;
      assert(0);
    }
  }

  return;
}

void BruteForceEcalGapFiller::addIntraModuleGapHits(LCCollectionVec* newcol) {
  // look for gaps within modules
  // i.e. between wafers, between towers

  streamlog_out(DEBUG1) << " starting addIntraModuleGapHits" << endl;

  const float verySmallDist = 0.01; // don't consider differences below this distance to be a gap
  const float slop = 0.01;          // flexibility, as ratio

  for (int il = 0; il < MAXLAYER; il++) {

    // we have to get the cell sizes here
    float cellsizeA, cellsizeB;
    if (_currentLayout == ECALBARREL) {
      cellsizeA = _caloGeomData->layers[il].cellSize0; // phi dir, in CM!!
      cellsizeB = _caloGeomData->layers[il].cellSize1; // z dir
    } else { // endcap - we should take care of the rotation with stave number if thse are different...
      cellsizeA = _caloGeomData->layers[il].cellSize0; // x dir
      cellsizeB = _caloGeomData->layers[il].cellSize1; // y dir
    }
    cellsizeA *= 10;
    cellsizeB *= 10; // convert to mm

    streamlog_out(DEBUG1) << "cell sizes in layer " << il << " = " << cellsizeA << " " << cellsizeB << endl;

    for (int is = 0; is < MAXSTAVE; is++) {
      for (int im = 0; im < MAXMODULE; im++) {
        std::vector<CalorimeterHit*> theseHits = hitsByLayerModuleStave[il][is][im];
        if (theseHits.size() > 1) {
          bool gap(false);
          float enFrac(0);
          for (size_t ih = 0; ih < theseHits.size() - 1; ih++) {
            for (size_t jh = ih + 1; jh < theseHits.size(); jh++) {
              float dist1d[3];
              for (int i = 0; i < 3; i++)
                dist1d[i] = fabs(theseHits[ih]->getPosition()[i] - theseHits[jh]->getPosition()[i]);
              float distXY = sqrt(pow(dist1d[0], 2) + pow(dist1d[1], 2));

              gap = false;
              if (_currentLayout == ECALBARREL) {
                if (dist1d[2] < verySmallDist &&        // same z coord
                    distXY > (1. + slop) * cellsizeA && // bigger than one cell period, smaller than two
                    distXY < (2. - slop) * cellsizeA) {
                  gap = true;
                  enFrac = (distXY - cellsizeA) / cellsizeA;
                } else if (distXY < verySmallDist && // same x-y coord
                           dist1d[2] > (1. + slop) * cellsizeB && dist1d[2] < (2. - slop) * cellsizeB) {
                  gap = true;
                  enFrac = (dist1d[2] - cellsizeB) / cellsizeB;
                }

              } else { // endcap
                if (dist1d[1] < verySmallDist && dist1d[0] > (1. + slop) * cellsizeA &&
                    dist1d[0] <
                        (2. - slop) * cellsizeA) { // be careful, if different size in x,y may have to worry about stave
                  gap = true;
                  enFrac = (dist1d[0] - cellsizeA) / cellsizeA;
                } else if (dist1d[0] < verySmallDist && dist1d[1] > (1. + slop) * cellsizeB &&
                           dist1d[1] <
                               (2. - slop) *
                                   cellsizeB) { // be careful, if different size in x,y may have to worry about stave
                  gap = true;
                  enFrac = (dist1d[1] - cellsizeB) / cellsizeB;
                }
              }
              if (gap) {

                streamlog_out(DEBUG1) << " GOT A GAP " << endl;

                float position[3] = {0.};
                for (int k = 0; k < 3; k++)
                  position[k] = 0.5 * (theseHits[ih]->getPosition()[k] + theseHits[jh]->getPosition()[k]);
                float extraEnergy = enFrac * (theseHits[ih]->getEnergy() + theseHits[jh]->getEnergy()) / 2.;
                float mintime = std::min(theseHits[ih]->getTime(), theseHits[jh]->getTime());
                CHT::CaloType cht_type = CHT::em;
                CHT::CaloID cht_id = CHT::ecal;
                CHT::Layout cht_lay = (_currentLayout == ECALBARREL)   ? CHT::barrel
                                      : (_currentLayout == ECALENDCAP) ? CHT::endcap
                                                                       : CHT::any;
                CalorimeterHitImpl* newGapHit = new CalorimeterHitImpl();
                newGapHit->setEnergy(_intraModuleFactor * std::log(1 + _intraModuleNonlinearFactor * extraEnergy) /
                                     _intraModuleNonlinearFactor);
                newGapHit->setPosition(position);
                newGapHit->setTime(mintime);
                newGapHit->setType(CHT(cht_type, cht_id, cht_lay, il));
                newcol->addElement(newGapHit);
              } // if gap
            } // jh
          } // ih
        } // >1 hit
      } // im
    } // is
  } // ilayer

  streamlog_out(DEBUG1) << " done addIntraModuleGapHits " << newcol->getNumberOfElements() << endl;
}

void BruteForceEcalGapFiller::addInterModuleGapHits(LCCollectionVec* newcol) {
  // look for gaps between modules
  //  compare hits in same stave, same layer

  streamlog_out(DEBUG1) << " starting addInterModuleGapHits" << endl;

  const float verySmallDist = 0.01; // don't consider differences below this distance to be a gap
  //  const float expectedGap = 12.3; // expected gap between cell centres across module boundary (mm)

  for (int il = 0; il < MAXLAYER; il++) {

    // we have to get the cell sizes here
    float cellsizeA, cellsizeB;
    if (_currentLayout == ECALBARREL) {
      cellsizeA = _caloGeomData->layers[il].cellSize0; // phi dir in CM
      cellsizeB = _caloGeomData->layers[il].cellSize1; // z dir
    } else { // endcap - we should take care of the rotation with stave number if thse are different...
      cellsizeA = _caloGeomData->layers[il].cellSize0; // x dir
      cellsizeB = _caloGeomData->layers[il].cellSize1; // y dir
    }
    cellsizeA *= 10; // to mm
    cellsizeB *= 10; // to mm

    for (int is = 0; is < MAXSTAVE; is++) {

      for (int im = 0; im < MAXMODULE; im++) {
        std::vector<CalorimeterHit*> theseHits = hitsByLayerModuleStave[il][is][im];

        if (theseHits.size() == 0)
          continue;

        // look in next module
        if (im + 1 >= 0 && im + 1 < MAXMODULE) {
          std::vector<CalorimeterHit*> nextHits = hitsByLayerModuleStave[il][is][im + 1];
          if (nextHits.size() == 0)
            continue;

          bool gap(false);
          float enFrac(0);

          for (size_t ih = 0; ih < theseHits.size(); ih++) {

            for (size_t jh = 0; jh < nextHits.size(); jh++) {

              float dist1d[3];
              for (int i = 0; i < 3; i++)
                dist1d[i] = fabs(theseHits[ih]->getPosition()[i] - nextHits[jh]->getPosition()[i]);
              float distXY = sqrt(pow(dist1d[0], 2) + pow(dist1d[1], 2));
              gap = false;

              if (_currentLayout == ECALBARREL) { // intermodule gaps only along z
                if (distXY < verySmallDist &&     // same phi coord
                    dist1d[2] < _interModuleDist +
                                    cellsizeB * 1.9) { // _interModuleDist is expected distance between sensor edges
                  gap = true;
                  enFrac = dist1d[2] / cellsizeB;
                }

              } else {                           // endcap
                if (dist1d[1] < verySmallDist && // same y
                    dist1d[0] <
                        _interModuleDist +
                            1.9 * cellsizeA) { // be careful, if different size in x,y may have to worrk about stave
                  gap = true;
                  enFrac = dist1d[0] / cellsizeA;
                } else if (dist1d[0] < verySmallDist &&                      // same x
                           dist1d[1] < _interModuleDist + 1.9 * cellsizeB) { // be careful, if different size in x,y may
                                                                             // have to worrk about stave
                  gap = true;
                  enFrac = dist1d[1] / cellsizeB;
                }
              }
              if (gap) {

                streamlog_out(DEBUG1) << " addInterModuleGapHits: found gap " << dist1d[0] << " " << dist1d[1] << " "
                                      << dist1d[2] << endl;

                float position[3] = {0.};
                for (int k = 0; k < 3; k++)
                  position[k] = 0.5 * (theseHits[ih]->getPosition()[k] + nextHits[jh]->getPosition()[k]);
                float extraEnergy = enFrac * (theseHits[ih]->getEnergy() + nextHits[jh]->getEnergy()) / 2.;
                float mintime = std::min(theseHits[ih]->getTime(), nextHits[jh]->getTime());
                CHT::CaloType cht_type = CHT::em;
                CHT::CaloID cht_id = CHT::ecal;
                CHT::Layout cht_lay = (_currentLayout == ECALBARREL)   ? CHT::barrel
                                      : (_currentLayout == ECALENDCAP) ? CHT::endcap
                                                                       : CHT::any;
                CalorimeterHitImpl* newGapHit = new CalorimeterHitImpl();
                newGapHit->setEnergy(_interModuleFactor * std::log(1 + _interModuleNonlinearFactor * extraEnergy) /
                                     _interModuleNonlinearFactor);
                newGapHit->setPosition(position);
                newGapHit->setTime(mintime);
                newGapHit->setType(CHT(cht_type, cht_id, cht_lay, il));
                newcol->addElement(newGapHit);
              } // if gap
            } // jh
          } // ih
        } // >1 hit
      } // im
    } // is
  } // ilayer

  streamlog_out(DEBUG1) << " done addInterModuleGapHits " << newcol->getNumberOfElements() << endl;

  return;
}

void BruteForceEcalGapFiller::end() {}
