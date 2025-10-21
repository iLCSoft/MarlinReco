#include "MCTruthJetEnergy.h"

// #include <iostream>
// #include <vector>
#include <set>
#include <sstream>
#include <string>

#include "UTIL/LCRelationNavigator.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <gearimpl/Vector3D.h>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#endif

using namespace lcio;
using namespace marlin;

MCTruthJetEnergy aMCTruthJetEnergy;

MCTruthJetEnergy::MCTruthJetEnergy() : Processor("MCTruthJetEnergy") {

  // modify processor description
  _description = "MCTruthJetEnergy calculates the MC truth jet energy of the MCParticles "
                 "associated to the reconstructed particles in the jet";

  // register steering parameters: name, description, class-variable, default value

  StringVec jCols;
  jCols.push_back("FTSelected_2Jets");
  jCols.push_back("FTSelected_3Jets");
  jCols.push_back("FTSelected_4Jets");
  jCols.push_back("FTSelected_5Jets");
  jCols.push_back("FTSelected_6Jets");

  registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE, "JetCollectionNames", "Names of the Jet collections",
                           _jetcolNames, jCols);

  registerInputCollection(LCIO::LCRELATION, "MCTruthRelationName",
                          "Name of the ReconstructedParticle-MCParticle  relation collection", _relName,
                          std::string("RecoMCTruthLink"));
}

void MCTruthJetEnergy::init() {

  // usually a good idea to
  printParameters();

  _nRun = 0;
  _nEvt = 0;

#ifdef MARLIN_USE_AIDA

  // create histograms for jet energy resolution

  int nCol = _jetcolNames.size();

  _jetEnergyHists.resize(nCol);

  for (int i = 0; i < nCol; i++) {

    std::stringstream str;
    str << _jetcolNames[i] << " _E ";

    _jetEnergyHists[i] = AIDAProcessor::histogramFactory(this)->createHistogram1D(str.str().c_str(),
                                                                                  "jet energy / Gev  ] ", 50, 0., 500.);

    //      createHistogram1D( str.str().c_str() , "jet energy resolution [ sigma_E/sqrt(E/Gev) ] ", 100 ,  ) ;
  }

  _jetEnergyTruthReco = AIDAProcessor::histogramFactory(this)->createHistogram2D(
      "EJetTruthVsRec ", " jet energy Reconstructed vs MCTruth", 500, 0., 300., 500, 0., 300.);

#endif
}

void MCTruthJetEnergy::processRunHeader(LCRunHeader* /*run*/) { _nRun++; }

void MCTruthJetEnergy::processEvent(LCEvent* evt) {

  int nCol = _jetcolNames.size();

  for (int i = 0; i < nCol; ++i) { // loop over jet collections

    LCCollection* jetcol = 0;

    try {

      jetcol = evt->getCollection(_jetcolNames[i]);
    }

    catch (lcio::DataNotAvailableException&) {

      streamlog_out(WARNING) << " collection " << _jetcolNames[i] << "  not found ! " << std::endl;

      continue; // try next collection name
    }

    LCCollection* relcol = evt->getCollection(_relName);

    int nJETS = jetcol->getNumberOfElements();

    streamlog_out(DEBUG1) << " found " << nJETS << " jets in event " << evt->getEventNumber() << "  in run "
                          << evt->getRunNumber() << std::endl;

    LCRelationNavigator rel(relcol);

    // loop over jets in this collection

    FloatVec jetEnergies;

    for (int j = 0; j < nJETS; j++) {

      ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(jetcol->getElementAt(j));

      streamlog_out(DEBUG0) << "   jet energy = " << jet->getEnergy() << std::endl;

      int nPart = jet->getParticles().size();

      // save all MCParticles used in this jet in a set (to avoid double counting)
      std::set<MCParticle*> mcpset;
      std::set<MCParticle*>::iterator mcpsetiter;

      for (int k = 0; k < nPart; k++) {

        const LCObjectVec& mcpvec = rel.getRelatedToObjects(jet->getParticles()[k]);

        if (mcpvec.size() > 0 && mcpvec[0] != 0) {

          MCParticle* mcp = dynamic_cast<MCParticle*>(mcpvec[0]);

          if (mcp != 0)
            mcpset.insert(mcp);
        } else {

          ReconstructedParticle* p = dynamic_cast<ReconstructedParticle*>(jet->getParticles()[k]);

          gear::Vector3D v(p->getMomentum()[0], p->getMomentum()[1], p->getMomentum()[2]);

          streamlog_out(DEBUG) << " found broken MCParticle link for particle [E:" << p->getEnergy() << "]"
                               << " theta: " << (v.theta() * 180. / 3.141592) << " charge: " << p->getCharge()
                               << " nDaught: " << p->getParticles().size() << std::endl;
        }
      }

      double mcJetEnergy = 0;

      for (mcpsetiter = mcpset.begin(); mcpsetiter != mcpset.end(); mcpsetiter++) {

        mcJetEnergy += (*mcpsetiter)->getEnergy();
      }

      streamlog_out(DEBUG0) << "MC particle jet energy:" << mcJetEnergy << std::endl;

      jetEnergies.push_back(mcJetEnergy);
    }

    // save the mc truth jet energies in the collection (same order as particles/jets )

    jetcol->parameters().setValues("MCTruthJetEnergies", jetEnergies);

  } // end loop over jet collections

  if (_nEvt % 100 == 0) {

    streamlog_out(MESSAGE2) << "processed event/run   " << evt->getEventNumber() << " / " << evt->getRunNumber()
                            << std::endl;
  }

  _nEvt++;
}

void MCTruthJetEnergy::check(LCEvent* evt) {

#ifdef MARLIN_USE_AIDA

  //     hMCPEnergy =
  //       AIDAProcessor::histogramFactory(this)->
  //       createCloud1D( "hMCPEnergy", "energy of the MCParticles", 100 ) ;

  int nCol = _jetcolNames.size();

  for (int i = 0; i < nCol; ++i) { // loop over jet collections

    LCCollection* jetcol = 0;

    try {

      jetcol = evt->getCollection(_jetcolNames[i]);
    }

    catch (lcio::DataNotAvailableException&) {

      streamlog_out(WARNING) << " collection " << _jetcolNames[i] << "  not found ! " << std::endl;

      continue; // try next collection name
    }

    //    AIDA::IHistogram1D* _eHist = _jetEnergyHists[ i ] ;

    FloatVec jetEnergies;
    jetcol->getParameters().getFloatVals("MCTruthJetEnergies", jetEnergies);

    int nJETS = jetcol->getNumberOfElements();

    if (jetEnergies.size() != unsigned(nJETS)) {

      streamlog_out(WARNING) << " wrong number of jet energy values in 'MCTruthJetEnergies' : " << jetEnergies.size()
                             << " with " << nJETS << " jets in the collection ! " << std::endl;

      continue;
    }

    for (int j = 0; j < nJETS; j++) {

      ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(jetcol->getElementAt(j));

      _jetEnergyHists[i]->fill(jet->getEnergy());

      _jetEnergyTruthReco->fill(jetEnergies[j], jet->getEnergy());
    }
  }
#endif
}

void MCTruthJetEnergy::end() {}
