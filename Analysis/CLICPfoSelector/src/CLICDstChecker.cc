#include "CLICDstChecker.h"
#include <CalorimeterHitType.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCObject.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/TrackImpl.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCTOOLS.h>
#include <gear/BField.h>
#include <gear/CalorimeterParameters.h>
#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/TPCParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/VXDParameters.h>
#include <marlin/Global.h>
#include <math.h>
#include <iostream>
#include <iostream>
#include <map>
#include "ClusterShapes.h"
#include "PfoUtilities.h"

#include <limits>

using namespace lcio;
using namespace marlin;

const int   precision        = 2;
const int   widthFloat       = 7;
const int   widthInt         = 5;
const int   widthSmallInt    = 3;
const float radiansToDegrees = 57.2957;

CLICDstChecker aCLICDstChecker;

CLICDstChecker::CLICDstChecker() : Processor("CLICDstChecker") {
  _description = "Prints out information in DST in a useful format";

  registerProcessorParameter("Debug", "Activate debugging?", m_debug, int(0));

  registerProcessorParameter("Monitoring", "Monitoring", m_monitoring, int(1));

  registerProcessorParameter("ShowBackground", "Show background pfo information", m_showBackground, int(1));

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "InputPfoCollection", "Input PFO Collection", m_inputPfoCollection,
                          std::string("PandoraPFANewPFOs"));

  registerProcessorParameter("PfoToMcRelationCollection", "Input PFO to MC particle relation Collection",
                             m_inputPfoToMcRelationCollection, std::string("RecoMCTruthLink"));

  registerInputCollection(LCIO::MCPARTICLE, "McParticleCollection", "Input MC particle collection",
                          m_inputMcParticleCollection, std::string("MCParticlesSkimmed"));
}

void CLICDstChecker::init() {
  printParameters();
  m_nRun = -1;
  m_nEvt = 0;
}

void CLICDstChecker::processRunHeader(LCRunHeader*) {
  m_nRun++;
  m_nEvt = 0;
  streamlog_out(MESSAGE) << std::endl;
  streamlog_out(MESSAGE) << "CLICDstChecker ---> new run : run number = " << m_nRun << std::endl;
}

void CLICDstChecker::processEvent(LCEvent* evt) {
  if (m_debug >= 1) {
    streamlog_out(MESSAGE) << std::endl;
    streamlog_out(MESSAGE) << "CLICDstChecker -> run = " << m_nRun << "  event = " << m_nEvt << std::endl;
    streamlog_out(MESSAGE) << std::endl;
  }

  std::vector<ReconstructedParticle*> physicsPfos;
  std::vector<MCParticle*>            physicsMatchedMcParticle;
  std::vector<ReconstructedParticle*> backgroundPfos;

  try {
    LCCollection* col   = evt->getCollection(m_inputPfoCollection.c_str());
    const int     nelem = col->getNumberOfElements();
    for (int iPfo = 0; iPfo < nelem; ++iPfo)
      m_pfoVector.push_back(static_cast<ReconstructedParticle*>(col->getElementAt(iPfo)));
    std::sort(m_pfoVector.begin(), m_pfoVector.end(), PfoUtil::PfoSortFunction);
  } catch (DataNotAvailableException& e) {
    streamlog_out(MESSAGE) << m_inputPfoCollection.c_str() << " collection is unavailable" << std::endl;
  };

  try {
    LCCollection* col   = evt->getCollection(m_inputMcParticleCollection.c_str());
    int           nelem = col->getNumberOfElements();
    for (int iMc = 0; iMc < nelem; ++iMc) {
      MCParticle* pMc = static_cast<MCParticle*>(col->getElementAt(iMc));
      // insert all existing mc pointers into set *** MUST CHECK THIS BEFORE USING RELATION ***
      m_mcSet.insert(pMc);
      const float px(pMc->getMomentum()[0]);
      const float py(pMc->getMomentum()[1]);
      const float pz(pMc->getMomentum()[2]);
      const float p          = sqrt(px * px + py * py + pz * pz);
      const float cosThetaMc = fabs(pz) / p;
      float       thetaMc    = 0.;
      if (cosThetaMc < 1.0)
        thetaMc = acos(cosThetaMc) * radiansToDegrees;
      if ((pMc->getParents()).size() == 0)
        streamlog_out(MESSAGE) << " Primary MC particle : " << pMc->getPDG() << " E = " << pMc->getEnergy()
                               << " cosT = " << cosThetaMc << " theta = " << thetaMc << std::endl;
    }
  } catch (DataNotAvailableException& e) {
    streamlog_out(MESSAGE) << m_inputMcParticleCollection.c_str() << " collection is unavailable" << std::endl;
  };

  m_pfoToMcNavigator = NULL;
  try {
    m_pfoToMcNavigator = new LCRelationNavigator(evt->getCollection(m_inputPfoToMcRelationCollection.c_str()));
  } catch (DataNotAvailableException& e) {
    streamlog_out(MESSAGE) << m_inputPfoToMcRelationCollection.c_str() << " collection is unavailable" << std::endl;
  }

  for (unsigned int iPfo = 0; iPfo < m_pfoVector.size(); ++iPfo) {
    ReconstructedParticle* pPfo      = m_pfoVector[iPfo];
    LCObjectVec            objectVec = m_pfoToMcNavigator->getRelatedToObjects(pPfo);
    if (objectVec.size() > 0) {
      for (unsigned int imc = 0; imc < objectVec.size(); imc++) {
        MCParticle* pMC = static_cast<MCParticle*>(objectVec[imc]);
        // since only saving skimmed set of MC particles, check pMC points to an existing object
        if (m_mcSet.count(pMC) != 0) {
          physicsPfos.push_back(pPfo);
          physicsMatchedMcParticle.push_back(pMC);
        } else {
          backgroundPfos.push_back(pPfo);
        }
      }
    }
  }

  float eSignal     = 0;
  float eBackground = 0;
  float eMc         = 0;

  for (unsigned int ipass = 0; ipass < 2; ipass++) {
    unsigned int nPfos;
    (ipass == 0) ? nPfos = physicsPfos.size() : nPfos = backgroundPfos.size();
    for (unsigned int iPfo = 0; iPfo < nPfos; ++iPfo) {
      if (iPfo == 0) {
        if (m_monitoring == 1 && ipass == 0) {
          streamlog_out(MESSAGE) << " *********************** Physics    Pfos ************************" << std::endl;
          streamlog_out(MESSAGE) << "   pfo   E       pT   cost  nt nc  mc   Emc " << std::endl;
        } else {
          if (m_monitoring == 1 && m_showBackground == 1 && ipass == 1)
            streamlog_out(MESSAGE) << " ********************** Background Pfos ************************" << std::endl;
          streamlog_out(MESSAGE) << "   pfo   E       pT   cost  nt nc " << std::endl;
        }
      }

      ReconstructedParticle* pPfo;
      (ipass == 0) ? pPfo = physicsPfos[iPfo] : pPfo = backgroundPfos[iPfo];
      MCParticle* pMc                                = NULL;
      if (ipass == 0) {
        pMc = physicsMatchedMcParticle[iPfo];
        eMc += pMc->getEnergy();
      }
      //      const int id = pPfo->id();
      const int type = pPfo->getType();
      //      const bool isCompound = pPfo->isCompound();
      float momentum[3];
      for (unsigned int i  = 0; i < 3; i++)
        momentum[i]        = pPfo->getMomentum()[i];
      const float pT       = sqrt(momentum[0] * momentum[0] + momentum[1] * momentum[1]);
      const float p        = sqrt(pT * pT + momentum[2] * momentum[2]);
      const float cosTheta = fabs(momentum[2]) / p;
      const float energy   = pPfo->getEnergy();
      //eTotalInput+=energy;
      //      const float mass = pPfo->getMass();
      //      const float charge = pPfo->getCharge();
      const ParticleIDVec particleIDs = pPfo->getParticleIDs();
      //      ParticleID *particleIDUsed = pPfo->getParticleIDUsed();
      //      const float goodnessOfPID = pPfo->getGoodnessOfPID();
      const ReconstructedParticleVec particles = pPfo->getParticles();
      const ClusterVec               clusters  = pPfo->getClusters();
      const TrackVec                 tracks    = pPfo->getTracks();

      //      for(unsigned int i = 0; i< tracks.size(); i++){
      //	Track *track = tracks[i];
      //	const int nHitsVTX = track->getSubdetectorHitNumbers()[6];
      //	const int nHitsFTD = track->getSubdetectorHitNumbers()[7];
      //	const int nHitsSIT = track->getSubdetectorHitNumbers()[8];
      //	const int nHitsTPC = track->getSubdetectorHitNumbers()[9];
      //	const int nHitsSET = track->getSubdetectorHitNumbers()[10];
      //	const int nHitsETD = track->getSubdetectorHitNumbers()[11];
      //	const float d0    = track->getD0();
      //	const float z0    = track->getZ0();
      //	const float omega = track->getOmega();
      //	const float tanL  = track->getTanLambda();
      //	const float phi0  = track->getPhi();
      //      }

      //      for(unsigned int i = 0; i< clusters.size(); i++){
      //	const Cluster *cluster = clusters[i];
      //      }

      (ipass == 0) ? eSignal += energy : eBackground += energy;

      if (m_monitoring && (ipass == 0 || (ipass == 1 && m_showBackground))) {
        std::stringstream output;
        output << std::fixed;
        output << std::setprecision(precision);
        if (clusters.size() == 0)
          FORMATTED_OUTPUT_TRACK(output, type, energy, pT, cosTheta, tracks.size(), 0);
        if (tracks.size() == 0)
          FORMATTED_OUTPUT_CLUSTER(output, type, energy, pT, cosTheta, 0, clusters.size());

        if (tracks.size() > 0 && clusters.size() > 0)
          FORMATTED_OUTPUT_TRACK_CLUSTER(output, type, energy, pT, cosTheta, tracks.size(), clusters.size());
        if (pMc != NULL)
          FORMATTED_OUTPUT_MC(output, pMc->getPDG(), pMc->getEnergy());

        streamlog_out(MESSAGE) << output.str();
        streamlog_out(MESSAGE) << std::endl;
      }
    }
  }

  if (m_monitoring) {
    streamlog_out(MESSAGE) << " ESignal = " << eSignal << "( " << eMc << " )" << std::endl;
    streamlog_out(MESSAGE) << " EBack   = " << eBackground << std::endl;
  }

  CleanUp();

  m_nEvt++;
}

void CLICDstChecker::CleanUp() {
  m_pfoVector.clear();
  m_mcSet.clear();
  delete m_pfoToMcNavigator;
  m_pfoToMcNavigator = NULL;
}

void CLICDstChecker::check(LCEvent*) {}

void CLICDstChecker::end() {}
