/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "anaPix.h"
#include "FPCCDData.h"
#include "FPCCDPixelHit.h"

#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <gear/VXDLayerLayout.h>
#include <gear/VXDParameters.h>

#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <vector>

using namespace lcio;
using namespace marlin;
using namespace std;

static int g_event = 0;
static double g_x = 0.;
static double g_y = 0.;
static double g_z = 0.;
static int g_layer = 0;
static int g_ladder = 0;
static double g_zeta = 0.;
static double g_xi = 0.;
static double g_zeta2 = 0.;
static double g_xi2 = 0.;
static double g_dE = 0.;
static double g_digi_dE = 0.;
static int g_quality = 0;
static int g_nparents = 0;

anaPix aanaPix;

// =====================================================================
anaPix::anaPix() : Processor("anaPix") {

  // modify processor description
  _description = "anaPix icteats TTree from FPCCDPixelHits";

  // register steering parameters: name, description, class-variable, default value
  FloatVec PixelSizeVec; // This is a temporary variable. Each ladder's PixelSize are different.
  for (int i = 0; i < 6; i++) {
    PixelSizeVec.push_back(0.005);
  }

  registerProcessorParameter("Each_FPCCD_pixelSize(mm)", "Each ladder's Pixel size of FPCCD (unit:mm) (default:0.005)",
                             _pixelSizeVec, PixelSizeVec);

  registerProcessorParameter("FPCCD_PixelSize", "Pixel size of FPCCD (unit:mm) (default: 0.005)", _pixelSize,
                             float(0.0050));
  /*
   registerProcessorParameter( "PointResolutionRPhi" ,
                               "R-Phi Resolution in VTX"  ,
                               _pointResoRPhi ,
                                float(0.001440)) ; // 5/sqrt(12)

   registerProcessorParameter( "PointResolutionZ" ,
                               "Z Resolution in VTX" ,
                               _pointResoZ ,
                               float(0.001440));   // 5/sqrt(12)
 */

  registerProcessorParameter("OutputROOTfileName", "Name of output ROOT file", _rootFileName,
                             std::string("./pixel.root"));

  // Input collections
  registerInputCollection(LCIO::LCGENERICOBJECT, "VTXPixelHitCollection", "Name of the VTX PixelHit collection",
                          _colNameVTX, std::string("VTXPixelHits"));
}

// =====================================================================
void anaPix::init() {

  // usually a good idea to
  printParameters();

  _nRun = 0;
  _nEvt = 0;

  InitGeometry();

  outroot = new TFile(_rootFileName.c_str(), "RECREATE");

  // hTreePix = new TTree("hTreePix","");
  hTreePix = new TTree("t", "");
  hTreePix->Branch("event", &g_event, "event/I");
  hTreePix->Branch("x", &g_x, "x/D");
  hTreePix->Branch("y", &g_y, "y/D");
  hTreePix->Branch("z", &g_z, "z/D");
  hTreePix->Branch("layer", &g_layer, "layer/I");
  hTreePix->Branch("ladder", &g_ladder, "ladder/I");
  hTreePix->Branch("dE", &g_dE, "dE/D");
  hTreePix->Branch("digi_dE", &g_digi_dE, "digi_dE/D");
  hTreePix->Branch("quality", &g_quality, "quality/I");
  hTreePix->Branch("multiplicity", &g_nparents, "multiplicity/I");
  hTreePix->Branch("xi", &g_xi, "xi/D");
  hTreePix->Branch("zeta", &g_zeta, "zeta/D");
  hTreePix->Branch("xi2", &g_xi2, "xi2/D");
  hTreePix->Branch("zeta2", &g_zeta2, "zeta2/D");
  /*
   hTreeLocalPix = new TTree("hTreeLocalPix","");
   hTreeLocalPix->Branch("event", &g_event, "event/I");
   hTreeLocalPix->Branch("zeta", &g_zeta,"zeta/D");
   hTreeLocalPix->Branch("xi", &g_xi, "xi/D");
   hTreeLocalPix->Branch("layer", &g_layer, "layer/I");
   hTreeLocalPix->Branch("ladder", &g_ladder, "ladder/I");
   hTreeLocalPix->Branch("dE", &g_dE, "dE/D");
   hTreeLocalPix->Branch("quality", &g_quality, "quality/I");
  */
}

// =====================================================================
void anaPix::InitGeometry() {
  // Save frequently used parameters.

  const gear::VXDParameters& gearVXD = Global::GEAR->getVXDParameters();
  const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout();

  _nLayer = layerVXD.getNLayers();
  _geodata.resize(_nLayer);
  _maxLadder = 0;

  for (int ly = 0; ly < _nLayer; ly++) {
    _geodata[ly].nladder = layerVXD.getNLadders(ly); // Number of ladders in this layer
    if (_maxLadder < _geodata[ly].nladder) {
      _maxLadder = _geodata[ly].nladder;
    }
    _geodata[ly].rmin = layerVXD.getSensitiveDistance(ly); // Distance of sensitive area from IP
    _geodata[ly].dphi = (2 * M_PI) / (double)_geodata[ly].nladder;
    _geodata[ly].phi0 = layerVXD.getPhi0(ly); // phi offset.
    _geodata[ly].sthick = layerVXD.getSensitiveThickness(ly);
    _geodata[ly].sximin = -layerVXD.getSensitiveOffset(ly) - layerVXD.getSensitiveWidth(ly) / 2.0;
    _geodata[ly].sximax = -layerVXD.getSensitiveOffset(ly) + layerVXD.getSensitiveWidth(ly) / 2.0;
    _geodata[ly].hlength = layerVXD.getSensitiveLength(ly);
    _geodata[ly].cosphi.resize(_geodata[ly].nladder);
    _geodata[ly].sinphi.resize(_geodata[ly].nladder);
    for (int ld = 0; ld < _geodata[ly].nladder; ld++) {
      double phi = _geodata[ly].phi0 + _geodata[ly].dphi * ld;
      _geodata[ly].cosphi[ld] = cos(phi);
      _geodata[ly].sinphi[ld] = sin(phi);
    }
  }
}

// =====================================================================
void anaPix::processRunHeader(LCRunHeader* /*run*/) { _nRun++; }

// =====================================================================
void anaPix::processEvent(LCEvent* evt) {
  LCCollection* pHitCol = 0;
  try {
    pHitCol = evt->getCollection(_colNameVTX);
  } catch (DataNotAvailableException& e) {
    if (_debug == 1) {
      std::cout << "Collection " << _colNameVTX.c_str() << " is unavailable in event " << _nEvt << std::endl;
    }
  }
  if (_debug == 1) {
    std::cout << " Collection =" << _colNameVTX << " nevt=" << _nEvt << std::endl;
    std::cout << " number of elements is " << pHitCol->getNumberOfElements() << std::endl;
  }
  if (pHitCol != 0) {
    FPCCDData theData(_nLayer, _maxLadder); // prepare object to make pixelhits
    int nhit = theData.unpackPixelHits(*pHitCol);
    if (_debug == 1) {
      theData.dump();
    }
    if (nhit > 0) { // Output Trackhit, if there are pixel hits
      g_event = _nEvt;
      fillTTree(theData);
    }
  } // End of process when VXD has hits
  _nEvt++;
}

// ====================================================================
void anaPix::fillTTree(FPCCDData& pHitCol) {

  // Convert pixel hits to TrackerHits

  for (int layer = 0; layer < _nLayer; layer++) {
    int nladder = _geodata[layer].nladder; // Number of ladders in this layer
    double rmin = _geodata[layer].rmin;    // Distance of sensitive area from IP
    double sthick = _geodata[layer].sthick;
    double sximin = _geodata[layer].sximin;
    for (int ladder = 0; ladder < nladder; ladder++) {
      PixelHitMap_t::iterator it = pHitCol.itBegin(layer, ladder);
      while (it != pHitCol.itEnd(layer, ladder)) {
        FPCCDPixelHit* aHit = (*it).second;
        int xiID = aHit->getXiID();
        int zetaID = aHit->getZetaID();
        int quality = aHit->getQuality();
        int nparents = aHit->getNMCParticles();
        double Edep = aHit->getEdep();
        double eta = sthick * 0.5;
        _pixelSize = _pixelSizeVec[layer];
        double xi = (xiID + 0.5) * _pixelSize;
        double newpos[3];
        g_zeta = (zetaID + 0.5) * _pixelSize;
        g_xi = xi;
        g_layer = layer;
        g_ladder = ladder;
        g_dE = Edep;
        // mori added//////(cf. FPCCDClustering)////////////////////
        double electronsPerKeV = 276.;
        int electronsPerStep = 25;
        int nbitsForEdep = 7;
        int nEle = static_cast<int>(g_dE * 1e+06 * electronsPerKeV);
        int nStep = 1 << nbitsForEdep;
        int count = nEle / electronsPerStep;
        if (count > nStep)
          count = nStep;
        g_digi_dE = (count * electronsPerStep) * 1e-6 / electronsPerKeV;
        /////////////////////////////////////////////////////////////

        g_quality = quality;
        g_nparents = nparents;

        newpos[0] = (xi + sximin) * _geodata[layer].sinphi[ladder] + (rmin + eta) * _geodata[layer].cosphi[ladder];
        newpos[1] = -(xi + sximin) * _geodata[layer].cosphi[ladder] + (rmin + eta) * _geodata[layer].sinphi[ladder];
        newpos[2] = (zetaID + 0.5) * _pixelSize - _geodata[layer].hlength;

        it++;

        g_x = newpos[0];
        g_y = newpos[1];
        g_z = newpos[2];

        // mori added (cf. FPCCDClustering)//////////////////////////////
        double lpos[3]; // xi,eta,zeta
        double gpos[3];
        gpos[0] = g_x;
        gpos[1] = g_y;
        gpos[2] = g_z;
        int tlayer = layer;
        int tladder = ladder;
        double sinphi = _geodata[tlayer].sinphi[tladder];
        double cosphi = _geodata[tlayer].cosphi[tladder];
        lpos[0] = gpos[0] * sinphi - gpos[1] * cosphi + gpos[2] * 0;
        lpos[1] = gpos[0] * cosphi + gpos[1] * sinphi + gpos[2] * 0;
        lpos[2] = gpos[0] * 0 + gpos[1] * 0 + gpos[2] * 1;
        g_xi2 = lpos[0];
        g_zeta2 = lpos[2];
        ////////////////////////////////////////////////////////////////////

        hTreePix->Fill();
        // hTreeLocalPix->Fill();
      }
    }
  }
}

// =====================================================================
void anaPix::check(LCEvent* /*evt*/) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

// =====================================================================
void anaPix::end() {
  streamlog_out(MESSAGE) << " end()  " << name() << " processed " << _nEvt << " events in " << _nRun << " runs "
                         << std::endl;
  outroot->Write();
}
