#ifndef STANDARD_INCLUDES_H
#define STANDARD_INCLUDES_H 1

#include "lcio.h"
#include <fstream>
#include <iostream>
#include <set>

// #include <math>
#include <algorithm> // for std::find
#include <iterator>  // for std::begin, std::end
#include <set>
#include <string>
#include <vector>

#include "CalorimeterHitType.h"

// ROOT
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TTree.h"

// LCIO
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Vertex.h>

//----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

#endif
