#include "Fox.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#ifdef MARLIN_USE_AIDA
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogramFactory.h>
#include <marlin/AIDAProcessor.h>
// #include <AIDA/IHistogram1D.h>
#endif
#include "IMPL/CalorimeterHitImpl.h"
#include "IMPL/ClusterImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCRelationImpl.h"
#include "IMPL/MCParticleImpl.h"
#include "IMPL/ParticleIDImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/SimCalorimeterHitImpl.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackerHitImpl.h"
#include "MCTree.h"

#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>

using namespace lcio;
using namespace marlin;
using namespace std;

vector<double> legendre_recursive(const double& x, const int& n, const vector<int> moments);

Fox aFox;

Fox::Fox() : Processor("Fox") {

  // modify processor description
  _description = "Fox calculates Fox-Wolfram moments";

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "NameOfReconstructedParticlesCollection",
                          "Name of the ReconstructedParticle collection", _colName, std::string("RecoParticles"));

  vector<int> momentsToCalculate;

  registerProcessorParameter("CalculateFoxWolframMoments",
                             "Numbers of the moments that are to be calculate 0-th is calculate by default",
                             _momentsToCalculate, momentsToCalculate);
}

void Fox::init() {

  // usually a good idea to
  printParameters();

  _nRun = 0;
  _nEvt = 0;
}

void Fox::processRunHeader(LCRunHeader* /*run*/) { _nRun++; }

void Fox::processEvent(LCEvent* evt) {

  // this gets called for every event
  // usually the working horse ...

  vector<ReconstructedParticle*> rpart;

  // fill histogram from LCIO data :
  try {
    LCCollection* col = evt->getCollection(_colName);

    vector<int> moments;
    moments.push_back(0);
    // to be shure that 0-th moment is calculated
    for (unsigned int kk = 0; kk < _momentsToCalculate.size(); kk++)
      moments.push_back(_momentsToCalculate[kk]);

    sort(moments.begin(), moments.end());
    unsigned int largest_moment;
    largest_moment = *max_element(moments.begin(), moments.end());

    vector<double> outputMoments;
    for (unsigned int jj = 0; jj < moments.size(); jj++)
      outputMoments.push_back(0.0);

    if (col != 0) {
      unsigned int nRecP = col->getNumberOfElements();

      for (unsigned int kk = 0; kk < nRecP; kk++) {
        rpart.push_back(dynamic_cast<ReconstructedParticle*>(col->getElementAt(kk)));
      }
      double* pi;
      double* pj;
      double eVisible = 0.0;
      double pimod;
      double pjmod;

      for (unsigned int i = 0; i < nRecP; i++) {
        pi = const_cast<double*>((rpart[i])->getMomentum());
        pimod = sqrt(pi[0] * pi[0] + pi[1] * pi[1] + pi[2] * pi[2]);
        for (unsigned int j = 0; j < nRecP; j++) {
          double cosThetaIJ = 0.0;
          pj = const_cast<double*>((rpart[j])->getMomentum());
          pjmod = sqrt(pj[0] * pj[0] + pj[1] * pj[1] + pj[2] * pj[2]);
          eVisible += (rpart[j])->getEnergy();
          if ((pimod != 0.0) && (pjmod != 0))
            cosThetaIJ = (pi[0] * pj[0] + pi[1] * pj[1] + pi[2] * pj[2]) / (pimod * pjmod);

          vector<double> legendreOutput(legendre_recursive(cosThetaIJ, largest_moment, moments));
          for (unsigned int kk = 0; kk < legendreOutput.size(); kk++) {
            outputMoments[kk] += pimod * pjmod * legendreOutput[kk];
          }
        }
      }

      for (unsigned int kk = 0; kk < outputMoments.size(); kk++) {
        outputMoments[kk] = outputMoments[kk] / (eVisible * eVisible);
      }

      for (unsigned int kk = 0; kk < outputMoments.size(); kk++) {
        ostringstream dtes;
        dtes << moments[kk];
        string nth = dtes.str();
        string base_name = "FoxWolfram_moment(";
        base_name = base_name + nth + ")";
        col->parameters().setValue(base_name.c_str(), (float)outputMoments[kk]);
      }

    } // col
  } catch (DataNotAvailableException& e) {
  }
  //  getchar();

  _nEvt++;
}

void Fox::check(LCEvent* /*evt*/) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void Fox::end() {

  //   std::cout << "MyProcessor::end()  " << name()
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
}

vector<double> legendre_recursive(const double& x, const int& n, const vector<int> moments) {

  vector<double> lr, rr;
  lr.push_back(1.0);
  lr.push_back(x);
  for (int i = 2; i < n + 1; i++)
    lr.push_back(((1.0 / i) * ((2.0 * i - 1.0) * x * lr[i - 1] - (i - 1.0) * lr[i - 2])));
  for (unsigned int i = 0; i < moments.size(); i++)
    rr.push_back(lr[moments[i]]);
  return rr;
}
