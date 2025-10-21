/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/*
** This file is part of the MarlinReco Project.
** Forming part of the SubPackage: SatoruJetFinder.
**
** For the latest version download from Web CVS:
** http://www-zeuthen.desy.de/lc-cgi-bin/cvsweb.cgi/marlinreco/?cvsroot=MarlinReco
**
** $Id: SatoruJetFinderProcessor.cc,v 1.8 2008-05-20 08:48:06 aplin Exp $
**
**
*/

#include "SatoruJetFinderProcessor.h"
#include <iostream>

#include <algorithm>
#include <math.h>

#include "EVENT/LCCollection.h"
#include "EVENT/LCIO.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"

using namespace std;

namespace marlin {

SatoruJetFinderProcessor aSatoruJetFinderProcessor;

SatoruJetFinderProcessor::SatoruJetFinderProcessor() : Processor("SatoruJetFinderProcessor") {
  _description = "A multi algorithm jet finder";

  // general steering parameters

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "InputCollection", "Collection of reconstructed particles",
                          _inputCollection, std::string("Unset"));

  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "OutputCollection", "Name of collection with the found jets",
                           _outputCollection, std::string("Unset"));

  registerOptionalParameter("Debug", "Set debug level", _debug, 0);

  registerOptionalParameter("SuccessTag", "Name of parameter added to event in case of successful jet finding",
                            _successTag, std::string("JetsFound"));

  // mode steering parameters

  registerProcessorParameter("Mode",
                             "Select predefined algorithms for jet finding"
                             "(or \"manual\")",
                             _jetFindingMode, std::string("manual"));

  registerOptionalParameter("NJetRequested",
                            "Force everything to N jets"
                            "(if supported by current mode)",
                            _nJetRequested, 4);

  registerOptionalParameter("YCut",
                            "YCut for jet finding algorithm"
                            "(if supported by current mode)",
                            _yCutParam, (float)0.);

  registerOptionalParameter("RCone",
                            "Half cone opening angle for cone jet finding "
                            "algorithm with variable number of jet",
                            _rConeParam, (float)0.7);
  registerOptionalParameter("EpsCone",
                            "Jet energycut for cone jet finding algorithm "
                            "with variable number of jets",
                            _epsConeParam, (float)7);

  // steering parameters for manual mode
  registerOptionalParameter("GlobalMode",
                            "mode for manual alogorithm selection, "
                            "see documentation for details",
                            _globalMode, std::string("0A"));

  registerOptionalParameter("Threshold", "Threshold, if mode is \"manual\"", _threshold, (float)0.);

  registerOptionalParameter("PrimaryJetFindingMode", "Primary jet finding mode, if mode is \"manual\"",
                            _primaryJetFindingMode, 0);

  registerOptionalParameter("SecondJetFindingMode", "Secong jet finding mode, if mode is \"manual\"",
                            _secondJetFindingMode, 0);

  registerOptionalParameter("MergingMode", "Merging mode, if mode is \"manual\"", _mergingMode, 0);

  registerOptionalParameter("MergingThreshold", "Merging threshold, if mode is \"manual\"", _mergingThreshold,
                            (float)0.);
}

void SatoruJetFinderProcessor::init() {

  streamlog_out(DEBUG4) << "SatoruJetFinderProcessor::init()  " << name() << std::endl
                        << "  parameters: " << std::endl
                        << *parameters();

  // write success tag only if parameter "SuccessTag" has been set
  _writeTag = parameters()->isParameterSet("SuccessTag");

  transform(_jetFindingMode.begin(), _jetFindingMode.end(), _jetFindingMode.begin(), (int (*)(int))std::tolower);
  streamlog_out(DEBUG4) << "jet finding mode:" << _jetFindingMode << endl;

  // Fixme! catch wrong mode names
  if (_jetFindingMode != "manual") {
    if (_jetFindingMode == "durhamnjet") {
      _globalMode = "0B";
      _primaryJetFindingMode = 5;
      streamlog_out(DEBUG4) << "GlobalMode " << _globalMode << " NJetRequested " << _nJetRequested << endl;
    } else if (_jetFindingMode == "durhamycut") {
      _globalMode = "0A";
      _primaryJetFindingMode = 5;
      _yCut[0] = _yCutParam;
    } else if (_jetFindingMode == "coneblanka") {
      _nJetRequested = 0;
      _threshold = 0.7;
      _primaryJetFindingMode = 12;
      _yCut[0] = 0.2;
      _yCut[1] = 0.7;
      _globalMode = "0A";
    } else if (_jetFindingMode == "satoru") {
      _threshold = 0.01;
      _nJetRequested = parameters()->getIntVal("NJetRequested");
      _primaryJetFindingMode = 5;
      _yCut[0] = 0.0;
      _mergingMode = 2;
      _mergingThreshold = 0.;
      _secondJetFindingMode = 4;
      _globalMode = "2C";
    }
  } else {
    // fixme: check whether values are set?
    _yCut[0] = _yCutParam;

    // different meaning for _yCut[] in cone mode
    if (_primaryJetFindingMode == MD::CONE) {
      _yCut[0] = _rConeParam;
      _yCut[1] = _epsConeParam;
    }
  }
}

void SatoruJetFinderProcessor::processRunHeader(LCRunHeader* run) {
  //    std::cout << "SatoruJetFinderProcessor::processRun()  " << name()
  //        << " in run " << run->getRunNumber()
  //        << std::endl ;
}

void SatoruJetFinderProcessor::processEvent(LCEvent* evt) {
  LCCollectionVec* JetsCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  // just a simple example
  // first write the name of all collection included ...

  // cout << " start" << endl;
  goSatoru(evt, JetsCol);
  // cout << " after" << endl;
}

void SatoruJetFinderProcessor::end() {}

/* *********************************************************************** */

void SatoruJetFinderProcessor::goSatoru(LCEvent* evt, LCCollection* JetsCol) {
  putPartons(evt);
  //    WritePartons();
  callSatoru(evt);
  getJets(evt, JetsCol);
  // GetPointerFromPartonToJet();
}

void SatoruJetFinderProcessor::putPartons(LCEvent* evt) {
  LCCollection* enflowcol = evt->getCollection(_inputCollection);
  int nenflow = enflowcol->getNumberOfElements();
  _partonsWorkArray.NumberOfPartons = nenflow;
  for (int ienflow = 0; ienflow < nenflow; ienflow++) {
    ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt(ienflow));
    for (int i = 0; i < 3; i++) {
      _partonsWorkArray.Momentum[ienflow * 4 + i] = (enflow->getMomentum())[i];
    }
    _partonsWorkArray.Momentum[ienflow * 4 + 3] = enflow->getEnergy();
  }
}

void SatoruJetFinderProcessor::writePartons() {
  for (int iparton = 0; iparton < _partonsWorkArray.NumberOfPartons; iparton++) {
    streamlog_out(DEBUG4) << "Px,Py,Pz,E: " << _partonsWorkArray.Momentum[iparton * 4 + 0] << ", "
                          << _partonsWorkArray.Momentum[iparton * 4 + 1] << ", "
                          << _partonsWorkArray.Momentum[iparton * 4 + 2] << ", "
                          << _partonsWorkArray.Momentum[iparton * 4 + 3] << endl;
  }
}

void SatoruJetFinderProcessor::callSatoru(LCEvent* evt) {
  int DimensionOfInputArray = 4;
  int DimensionOfOutputArray = 4;
  int MaximalNumberOfJets = 20;
  //    float YMinus,YPlus;
  int IError;
  int GlobalModusLength = _globalMode.length();

  _YMinus = -9999;
  _YPlus = -9999;
  syjkrn_(_globalMode.c_str(), _nJetRequested, _threshold, _primaryJetFindingMode, _yCut, _mergingMode,
          _mergingThreshold, _secondJetFindingMode, _partonsWorkArray.NumberOfPartons, DimensionOfInputArray,
          _partonsWorkArray.Momentum, MaximalNumberOfJets, _jetsWorkArray.NumberOfJets,
          _partonsWorkArray.PointerParticleToJets, DimensionOfOutputArray, _jetsWorkArray.Momentum, _YMinus, _YPlus,
          IError, GlobalModusLength);
}

void SatoruJetFinderProcessor::getJets(LCEvent* evt, LCCollection* JetsCol) {
  LCCollection* enflowcol = evt->getCollection(_inputCollection);
  streamlog_out(DEBUG4) << " " << endl;
  streamlog_out(DEBUG4) << " number of jets found: " << _jetsWorkArray.NumberOfJets << endl;
  for (int ijets = 0; ijets < _jetsWorkArray.NumberOfJets; ijets++) {
    ReconstructedParticleImpl* Jets = new ReconstructedParticleImpl;
    float momentum[3], energy;
    momentum[0] = _jetsWorkArray.Momentum[4 * ijets + 0];
    momentum[1] = _jetsWorkArray.Momentum[4 * ijets + 1];
    momentum[2] = _jetsWorkArray.Momentum[4 * ijets + 2];
    energy = _jetsWorkArray.Momentum[4 * ijets + 3];
    Jets->setMomentum(momentum);
    Jets->setEnergy(energy);
    // JL June 20, 2016: add jet mass!
    double mass = energy * energy - momentum[0] * momentum[0] - momentum[1] * momentum[1] - momentum[2] * momentum[2];
    mass = (mass > 0) ? sqrt(mass) : 0;
    Jets->setMass(mass);
    // end JL
    for (int iobj = 0; iobj < _partonsWorkArray.NumberOfPartons; iobj++) {
      if (_partonsWorkArray.PointerParticleToJets[iobj] == (ijets + 1)) {
        ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt(iobj));
        Jets->addParticle(enflow);
      }
    }
    JetsCol->addElement(Jets);
  }

  JetsCol->parameters().setValue("YMinus", _YMinus);
  JetsCol->parameters().setValue("YPlus", _YPlus);

  evt->addCollection(JetsCol, _outputCollection);
  if (_writeTag == true) {
    /// \todo fixme: check for requested number of Jets...
    if (JetsCol->getNumberOfElements() > 0)
      evt->parameters().setValue(_successTag, 1);
    else
      evt->parameters().setValue(_successTag, 0);
  }
}

void SatoruJetFinderProcessor::writeJets() {}
void SatoruJetFinderProcessor::writeParameters() {}

} // namespace marlin
