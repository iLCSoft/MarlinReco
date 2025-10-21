// History:
//- improved clustering procedure, two different algorithms implemented: COG & analog head-tail, Z. Drasal Apr 2010
//- only 1 relation TrackerPulse <--> SimTrackerHit read in, Z. Drasal May 2010
//- output relation TrackerHit <--> SimTrackerHit not used, SimTrackerHits (without weights) accessible via
// TrackerHit->getRawHits() method, Z. Drasal May 2010
//- corrected charge units, Z. Drasal May 2010
//- TrackerHit -> rawHits() contains 1 SimTrackerHit, which contributed with highest weight, Z. Drasal Oct 2010

#include "SiStripClus.h"
#include "Colours.h"
#include "PhysicalConstants.h"

#include "SiStripGeomBuilder.h"

#include <iomanip>
#include <math.h>
#include <stdlib.h>

// Include LCIO header files
#include "UTIL/LCTrackerConf.h"
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
#include <UTIL/LCRelationNavigator.h>
#include <lcio.h>

// Include Marlin
#include <streamlog/streamlog.h>

// Namespaces
using namespace CLHEP;

namespace sistrip {

//
// Instantiate this object
//
SiStripClus anSiStripClus;

//
// Constructor
//
SiStripClus::SiStripClus() : Processor("SiStripClus") {
  // Processor description
  _description =
      "SiStripClus: Marlin processor intended for cluster finding - using digitized data worked out with SiStripDigi";

  // Define compulsory parameters
  registerProcessorParameter("Subdetector", "Subdetector to be digitised", _subdetector, std::string("FTD"));

  registerProcessorParameter("CMSnoise", "Common mode subtracted noise, set in electrons", _CMSnoise, float(1000));

  registerProcessorParameter("FloatingStripsRPhi", "Is every even strip floating in R-Phi?", _floatStripsRPhi,
                             bool(false));

  registerProcessorParameter("FloatingStripsZ", "Is every even strip floating in Z?", _floatStripsZ, bool(false));

  std::vector<float> resSVDFirstInRPhi;
  resSVDFirstInRPhi.push_back(6.5); // 20 degrees
  resSVDFirstInRPhi.push_back(6.0); // 30 degrees
  resSVDFirstInRPhi.push_back(6.0); // 40  degrees

  resSVDFirstInRPhi.push_back(7.0);
  resSVDFirstInRPhi.push_back(7.0);
  resSVDFirstInRPhi.push_back(6.0); // 50 , 60 , 70  degrees

  resSVDFirstInRPhi.push_back(6.5);
  resSVDFirstInRPhi.push_back(6.5);
  resSVDFirstInRPhi.push_back(6.5); // 80 , 90 , 100 degrees

  resSVDFirstInRPhi.push_back(6.5);
  resSVDFirstInRPhi.push_back(6.5);
  resSVDFirstInRPhi.push_back(6.5); // 110, 120, 130 degrees

  resSVDFirstInRPhi.push_back(6.5);
  resSVDFirstInRPhi.push_back(6.0); // 140, 150 degrees

  registerProcessorParameter("ResolutionOfSVDFirstInRPhi", "Resolution of first SVD layer in R-Phi (in um)",
                             _resSVDFirstInRPhi, resSVDFirstInRPhi);

  std::vector<float> resSVDOtherInRPhi;
  resSVDOtherInRPhi.push_back(7.5);
  resSVDOtherInRPhi.push_back(9.5);
  resSVDOtherInRPhi.push_back(10.5); // 20 , 30 , 40  degrees

  resSVDOtherInRPhi.push_back(11.0);
  resSVDOtherInRPhi.push_back(11.0);
  resSVDOtherInRPhi.push_back(11.0); // 50 , 60 , 70  degrees

  resSVDOtherInRPhi.push_back(11.0);
  resSVDOtherInRPhi.push_back(11.0);
  resSVDOtherInRPhi.push_back(11.0); // 80 , 90 , 100 degrees

  resSVDOtherInRPhi.push_back(11.0);
  resSVDOtherInRPhi.push_back(11.0);
  resSVDOtherInRPhi.push_back(11.0); // 110, 120, 130 degrees

  resSVDOtherInRPhi.push_back(11.0);
  resSVDOtherInRPhi.push_back(10.0); // 140, 150 degrees

  registerProcessorParameter("ResolutionOfSVDOtherInRPhi", "Resolution of other SVD layers in R-Phi (in um) ",
                             _resSVDOtherInRPhi, resSVDOtherInRPhi);

  std::vector<float> resSVDFirstInZ;
  resSVDFirstInZ.push_back(20.0);
  resSVDFirstInZ.push_back(16.0);
  resSVDFirstInZ.push_back(14.0); // 20 , 30 , 40  degrees

  resSVDFirstInZ.push_back(14.0);
  resSVDFirstInZ.push_back(13.0);
  resSVDFirstInZ.push_back(13.0); // 50 , 60 , 70  degrees

  resSVDFirstInZ.push_back(28.0);
  resSVDFirstInZ.push_back(28.0);
  resSVDFirstInZ.push_back(28.0); // 80 , 90 , 100 degrees

  resSVDFirstInZ.push_back(14.0);
  resSVDFirstInZ.push_back(13.0);
  resSVDFirstInZ.push_back(13.0); // 110, 120, 130 degrees

  resSVDFirstInZ.push_back(14.0);
  resSVDFirstInZ.push_back(16.0); // 140, 150 degrees

  registerProcessorParameter("ResolutionOfSVDFirstInZ", "Resolution of first SVD layer along Z axis (in um)",
                             _resSVDFirstInZ, resSVDFirstInZ);

  std::vector<float> resSVDOtherInZ;
  resSVDOtherInZ.push_back(16.0);
  resSVDOtherInZ.push_back(16.0);
  resSVDOtherInZ.push_back(15.0); // 20 , 30 , 40  degrees

  resSVDOtherInZ.push_back(14.0);
  resSVDOtherInZ.push_back(17.0);
  resSVDOtherInZ.push_back(34.0); // 50 , 60 , 70  degrees

  resSVDOtherInZ.push_back(50.0);
  resSVDOtherInZ.push_back(50.0);
  resSVDOtherInZ.push_back(50.0); // 80 , 90 , 100 degrees

  resSVDOtherInZ.push_back(34.0);
  resSVDOtherInZ.push_back(14.0);
  resSVDOtherInZ.push_back(14.0); // 110, 120, 130 degrees

  resSVDOtherInZ.push_back(16.0);
  resSVDOtherInZ.push_back(16.0); // 140, 150 degrees

  registerProcessorParameter("ResolutionOfSVDOtherInZ", "Resolution of other SVD layers along Z axis (in um) ",
                             _resSVDOtherInZ, resSVDOtherInZ);

  registerProcessorParameter("S/NSeedStrips", "Signal to noise ratio cut for seed strips", _SNseed, float(5));

  registerProcessorParameter("S/NAdjacentStrips", "Signal to noise ratio cut for adjacent strips", _SNadjacent,
                             float(3));

  registerProcessorParameter("S/NCluster", "Signal to noise ratio cut for total cluster", _SNtotal, float(8));

  //   registerProcessorParameter( "TanOfAvgELorentzShift",
  //                             "Tangent of electrons' average Lorentz shift (0.20 for 273 K, 0.175 for 300 K)",
  //                               _TanOfAvgELorentzShift,
  //                               float(0.175) );

  registerProcessorParameter("TanOfAvgHLorentzShift",
                             "Tangent of holes' average Lorentz shift (0.05 for 273 K, 0.039 for 300 K)",
                             _TanOfAvgHLorentzShift, float(0.039));

  registerProcessorParameter("InputCollectionName", "Name of TrackerPulse input collection", _inColName,
                             std::string("FTDDigits"));

  registerProcessorParameter("OutputCollectionName", "Name of TrackerHitPlane output collection", _outColName,
                             std::string("FTDTrackerHits"));

  registerProcessorParameter("RelCollectionNamePlsToSim",
                             "Name of input relation collection - TrackerPulse to SimTrackerHit (if nonzero, required)",
                             _relColNamePlsToSim, std::string("FTDDigitsToSimHitsRel"));

  registerProcessorParameter("pitchFront", "Pitch of the front sensor family in the center of the sensors, set in um",
                             _pitchFront, double(50.0));

  registerProcessorParameter("pitchRear", "Pitch of the rear sensor family in the center of the sensors, set in um",
                             _pitchRear, double(50.0));
}

//
// Method called at the beginning of data processing
//
void SiStripClus::init() {
  // Set variables in appropriate physical units
  _CMSnoise *= e;

  for (unsigned int i = 0; i < _resSVDFirstInRPhi.size(); ++i) {
    _resSVDFirstInRPhi[i] *= um;
    _resSVDOtherInRPhi[i] *= um;
    _resSVDFirstInZ[i] *= um;
    _resSVDOtherInZ[i] *= um;
  }

  _pitchFront *= um;
  _pitchRear *= um;

  // Get geometry parameters from Gear xml file
  _geometry = SiStripGeomBuilder::Build("FTD", _pitchFront, _pitchRear);
  // FIXME: Provisional para compilar
  _geometry->initGearParams();

  // Print set parameters
  printProcessorParams();

  //
  // ROOT variables
  //
#ifdef ROOT_OUTPUT
  _rootFile = new TFile("BelleII_SVD_Hits.root", "recreate");

  _rootFile->cd("");

  // Declare Tree
  _rootTree = new TTree("Hits", "Hit info");

  _rootTree->Branch("Layer", &_rootLayerID, "Layer/I");
  _rootTree->Branch("Ladder", &_rootLadderID, "Ladder/I");
  _rootTree->Branch("Sensor", &_rootSensorID, "Sensor/I");
  _rootTree->Branch("RPhiSim", &_rootSimRPhi, "SimRPhi/D");
  _rootTree->Branch("RPhiRec", &_rootRecRPhi, "RecRPhi/D");
  _rootTree->Branch("ResRPhi", &_rootResRPhi, "ResRPhi/D");
  _rootTree->Branch("ResR", &_rootResR, "ResR/D");
  _rootTree->Branch("ResModule", &_rootResModule, "Resmodule/D");
  _rootTree->Branch("RPhiClsSize", &_rootClsSizeRPhi, "ClsSizeRPhi/D");
  _rootTree->Branch("RPhiMCPDG", &_rootMCPDGRPhi, "MCPDGRPhi/I");
  _rootTree->Branch("ZSim", &_rootSimZ, "SimZ/D");
  _rootTree->Branch("ZRec", &_rootRecZ, "RecZ/D");
  _rootTree->Branch("ZClsSize", &_rootClsSizeZ, "ClsSizeZ/D");
  _rootTree->Branch("ZMCPDG", &_rootMCPDGZ, "MCPDGZ/I");
  _rootTree->Branch("EvtNum", &_rootEvtNum, "EvtNum/I");
#endif
}

//
// Method called for each run
//
void SiStripClus::processRunHeader(LCRunHeader* run) {
  // Print run number
  //   streamlog_out(MESSAGE3) << DGREEN
  //                           << " Processing run: "
  //                           << ENDCOLOR
  //                           << (run->getRunNumber()+1)
  //                           << std::endl << std::endl;

  _nRun++;
}

//
// Method called for each event
//
void SiStripClus::processEvent(LCEvent* event) {
  // Get names of all collections saved in LCIO file
  ConstStringVec* strVec = event->getCollectionNames();

  // Initialize relation navigators
  _navigatorPls = NULL;

  // Initialize
#ifdef ROOT_OUTPUT
  _rootEvtNum = event->getEventNumber();
#endif

  //
  // Get TrackerPulse collections and relation to MCParticles collection
  LCCollection* colOfTrkPulses = 0;

  // Print header info about collections (for 1. event)
  if (event->getEventNumber() == 0) {
    streamlog_out(MESSAGE3) << " " << DUNDERL << DBLUE << "LCCollection(s) processed:" << ENDCOLOR << std::endl
                            << std::endl;
  }

  // Go through all collection names
  for (ConstStringVec::const_iterator colName = strVec->begin(); colName != strVec->end(); ++colName) {
    // Tracker pulses
    if (_inColName == (*colName)) {
      // Collection must be of TrackerPulse type
      if ((event->getCollection(*colName))->getTypeName() == LCIO::TRACKERPULSE) {
        // Save collection
        colOfTrkPulses = event->getCollection(*colName);
        // Print info
        if (event->getEventNumber() == 0) {
          streamlog_out(MESSAGE3) << "  " << *colName << std::endl;
        }
      }
      // Collection NOT of TrackerPulse type
      else {
        streamlog_out(ERROR) << "Required collection: " << _inColName << " found, but NOT of TRACKERPULSE type!!!"
                             << std::endl;
        exit(0);
      }
    }

    // Relation TrkPulses to SimTrackerHits
    if ((!_relColNamePlsToSim.empty()) && (_relColNamePlsToSim == (*colName))) {
      // Collection must be of relation type
      if ((event->getCollection(*colName))->getTypeName() == LCIO::LCRELATION) {
        // Save collection into relation navigator
        _navigatorPls = new LCRelationNavigator(event->getCollection(*colName));
        // Print info
        if (event->getEventNumber() == 0) {
          streamlog_out(MESSAGE3) << "  " << *colName << std::endl;
        }
      }
      // Collection NOT of relation type
      else {
        streamlog_out(ERROR) << "Required collection: " << _relColNamePlsToSim
                             << " found, but NOT of LCRELATION type!!!" << std::endl;
        exit(0);
      }
    }
  } // For - collection names

  // Input collection NOT found
  if (colOfTrkPulses == 0) {
    streamlog_out(WARNING) << "Required collection: " << _inColName << " not found!!!" << std::endl;
  }

  if ((!_relColNamePlsToSim.empty()) && (_navigatorPls == NULL)) {
    streamlog_out(WARNING) << "Required collection: " << _relColNamePlsToSim << " not found!!!" << std::endl;
  }

  //
  // Process TrackerPulse collections
  if (colOfTrkPulses != 0) {
    CellIDDecoder<TrackerPulse> cellIDDec(colOfTrkPulses);

    // Initialize variables
    TrackerPulseImpl* pulse = 0;

    // Sensor map with map of hit strips (corresponding to the given sensor)
    SensorStripMap sensorMap;

    // Get number of elements in each collection
    int nPulses = colOfTrkPulses->getNumberOfElements();

    // Process pulses
    for (int i = 0; i < nPulses; ++i) {
      // Copy the content of collection TrackerPulse to a pulse
      pulse = dynamic_cast<TrackerPulseImpl*>(colOfTrkPulses->getElementAt(i));

      // Update the sensor strip map with the pulse
      updateMap(pulse, sensorMap);
    }

    //
    // Find clusters
    ClsVec clsVec = findClus(sensorMap);

    // Releasing memory and clearing
    releaseMap(sensorMap);

    //
    // Calculate real + ghost hits from clusters - in global ref. system +
    // create relations to MCParticles
    IMPL::LCCollectionVec* colOfTrkHits = new IMPL::LCCollectionVec(LCIO::TRACKERHITPLANE);

    // Calculate the hits
    calcHits(clsVec, colOfTrkHits);

    //
    // Save the collection (vector) of hits + relation to MC
    event->addCollection(colOfTrkHits, _outColName);

  } // If (colOfTrkPulses != 0)
  // Release memory
  if (_navigatorPls != NULL) {
    delete _navigatorPls;
    _navigatorPls = NULL;
  }

  _nEvent++;
}

//
// Method called after each event to check the data processed
//
void SiStripClus::check(LCEvent* event) {}

//
// Method called after all data processing
//
void SiStripClus::end() {
  // Release memory
  delete _geometry;
  _geometry = 0;

#ifdef ROOT_OUTPUT

  // Close file
  _rootFile->cd("");
  _rootFile->Write();
  _rootFile->Close();

#endif

  // Print message
  streamlog_out(MESSAGE3) << std::endl << " " << DGREEN << "Processor successfully finished!" << ENDCOLOR << std::endl;
}

// MAIN CLUSTER METHOD

//
// Method searching for clusters
//
typedef std::pair<int, StripCluster*> StripClusterPair;
typedef std::map<int, std::map<StripType, std::vector<StripClusterPair>>> SensorStripClusterMap;

ClsVec SiStripClus::findClus(SensorStripMap& sensorMap) {
  ClsVec clsVec;

  // Cluster vectors -
  SensorStripClusterMap clsvectFrontRear;

  std::vector<StripType> stv;
  stv.push_back(STRIPFRONT);
  stv.push_back(STRIPREAR);
  // As the sensorID was lost (see SiStripClus::updateMap method)
  // assigning it
  std::map<StripType, int> stSensorIDmap;
  stSensorIDmap[STRIPFRONT] = 1;
  stSensorIDmap[STRIPREAR] = 3;

  // Bunch of strips forming cluster
  StripChargeMap clsStrips;

  //
  // Search complete sensor map - find seeds & their neghbouring strips
  for (SensorStripMap::iterator iterSMap = sensorMap.begin(); iterSMap != sensorMap.end(); ++iterSMap) {
    // Save layer ID , ...
    const int cellID = iterSMap->first;
    std::map<std::string, int> bfmap = _geometry->decodeCellID(cellID);
    const int layerID = bfmap["layer"];
    const int ladderID = bfmap["module"];

    // Storing the two types of clusters
    for (std::vector<StripType>::iterator itST = stv.begin(); itST != stv.end(); ++itST) {
      StripType STRIPTYPE = *itST;
      for (StripChargeMap::iterator iterChMap = iterSMap->second[STRIPTYPE].begin();
           iterChMap != iterSMap->second[STRIPTYPE].end(); iterChMap++) {
        const int sensorID = stSensorIDmap[STRIPTYPE];

        // Begin algorithm
        // Zero: Save new candidate for seed strip + MC true info
        const int seedStrip = iterChMap->first;
        const double seedCharge = iterChMap->second->getCharge();

        const SimTrackerHitMap& seedSimHitMap = iterChMap->second->getSimHitMap();

        if (seedCharge < (_SNseed * _CMSnoise)) {
          continue;
        }

        // First: New cluster and its seed strip has been found
        clsStrips[seedStrip] = new Signal(seedCharge, 0.);
        clsStrips[seedStrip]->updateSimHitMap(seedSimHitMap);

        // Set this charge as zero to avoid double counting
        iterChMap->second->setCharge(0.);

        // Map ordered from lower to higher
        // Continue searching - find left and right neighbours
        // Taking account also the floating strips case...
        StripChargeMap::iterator seedIt = iterSMap->second[STRIPTYPE].find(seedStrip);
        // Second: search for left neighbours
        //-- Left neighbours
        StripChargeMap::reverse_iterator lschMap(seedIt);
        clsStrips = storeHitsAdjacents<StripChargeMap::reverse_iterator>(
            clsStrips, lschMap, iterSMap->second[STRIPTYPE].rend(), iterSMap->second[STRIPTYPE]);

        // Third: search for rigth neighbours
        //-- Right neighbours
        StripChargeMap::iterator rschMap(++seedIt);
        clsStrips = storeHitsAdjacents<StripChargeMap::iterator>(clsStrips, rschMap, iterSMap->second[STRIPTYPE].end(),
                                                                 iterSMap->second[STRIPTYPE]);

        // Fourth: Calculate mean position of a new cluster
        SimTrackerHitMap clsSimHitMap;

        // Cluster: position, charge & size
        double clsCharge = 0.0;

        double xLeftSignal = 0.0;
        double qLeftSignal = 0.0;

        double xRightSignal = 0.0;
        double qRightSignal = 0.0;

        double qIntermSignal = 0.0;

        int stripID = 0;
        for (StripChargeMap::iterator iterChMap2 = clsStrips.begin(); iterChMap2 != clsStrips.end(); iterChMap2++) {
          // Current strip ID, posZ & charge
          stripID = iterChMap2->first;
          const double stripPosYatz0 = _geometry->getStripPosY(layerID, sensorID, stripID, 0.0);

          const double stripCharge = iterChMap2->second->getCharge();

          // Update info about MC particles which contributed
          const SimTrackerHitMap& simHitMap = iterChMap2->second->getSimHitMap();

          if (simHitMap.size() != 0) {
            for (SimTrackerHitMap::const_iterator iterSHM = simHitMap.begin(); iterSHM != simHitMap.end(); ++iterSHM) {
              EVENT::SimTrackerHit* simHit = dynamic_cast<EVENT::SimTrackerHit*>(iterSHM->first);
              float weight = iterSHM->second;
              if (clsSimHitMap.find(simHit) != clsSimHitMap.end()) {
                clsSimHitMap[simHit] += weight;
              } else {
                clsSimHitMap[simHit] = weight;
              }
            }
          }

          // Get last iterator

          // Get leftmost signal
          if (iterChMap2 == clsStrips.begin()) {
            xLeftSignal = stripPosYatz0;
            qLeftSignal = stripCharge;
          }

          // Get rightmost signal
          else if (iterChMap2 == --(clsStrips.end())) {
            xRightSignal = stripPosYatz0;
            qRightSignal = stripCharge;
          }
          // Get intermediate signal
          else {
            qIntermSignal += stripCharge;
          }

          // Update total charge
          clsCharge += stripCharge;
        }

        // Number of strips being part of cluster
        int clsSize = clsStrips.size();
        // Get average intermediate signal
        if (clsSize > 2) {
          qIntermSignal /= (clsSize - 2.0);
        } else {
          qIntermSignal = 0.;
        }

        // Building the cluster
        double geomPitchAtz0 = _geometry->getSensorPitch(layerID, sensorID, 0.0);
        double readOutPitch = 0.;
        if (_floatStripsZ) // FIXME---> CAMBIAR POR mapa _floatStrips[type]
        {
          readOutPitch = 2.0 * geomPitchAtz0;
        } else {
          readOutPitch = geomPitchAtz0;
        }

        double clsPosYatz0 = 0.0;
        // Analog head-tail algorithm
        // qIntermSignal >= 1/eps*qRighSignal || 1/eps*qLeftSignal
        //    --> if not don't use head-tail (delta electron ...) - use
        //        epsilon ~ 1.5
        if ((clsSize > 2) && (qLeftSignal < 1.5 * qIntermSignal) && (qRightSignal < 1.5 * qIntermSignal)) {
          clsPosYatz0 =
              (xRightSignal + xLeftSignal) / 2. + (qRightSignal - qLeftSignal) / 2. / qIntermSignal * readOutPitch;
        } else // COG algorithm
        {
          clsPosYatz0 = (xRightSignal * qRightSignal + xLeftSignal * qLeftSignal) / (qRightSignal + qLeftSignal);
        }

        //
        // Fifth: Correct mean position to Lorentz shift
        // As theta is Pi/2.0 --> this correction is 0: FIXME-> TO BE REMOVED
        clsPosYatz0 +=
            _TanOfAvgHLorentzShift * _geometry->getSensorThick(layerID) / 2. * cos(_geometry->getLadderTheta(layerID));

        //
        // Sixth: Save information about new cluster
        if (clsCharge >= (_SNtotal * _CMSnoise)) {
          StripCluster* pCluster = new StripCluster(
              layerID, ladderID, sensorID, Hep3Vector(_geometry->getSensorThick(layerID) / 2., clsPosYatz0, 0.0),
              Hep3Vector(_geometry->getSensorThick(layerID) / 2., 0., 0.), clsCharge, clsSize);
          pCluster->updateSimHitMap(clsSimHitMap);

          clsvectFrontRear[cellID][STRIPTYPE].push_back(std::pair<int, StripCluster*>(stripID, pCluster));
        }

        // Release memory
        for (StripChargeMap::iterator iterChMap2 = clsStrips.begin(); iterChMap2 != clsStrips.end(); iterChMap2++) {
          delete iterChMap2->second; // Signal pointers
        }

        // Clear content
        clsStrips.clear();
      }
    } // For cluster strip type (front-rear)
  } // For sensor map

  // Stores the hit (using the STRIPFRONT as init)
  for (SensorStripClusterMap::iterator it = clsvectFrontRear.begin(); it != clsvectFrontRear.end(); ++it) {
    const int cellID = it->first;
    std::map<std::string, int> bfmap = _geometry->decodeCellID(cellID);
    const int layerID = bfmap["layer"];
    const int ladderID = bfmap["module"];

    // the front sensors
    for (std::vector<StripClusterPair>::iterator itFront = it->second[STRIPFRONT].begin();
         itFront != it->second[STRIPFRONT].end(); ++itFront) {
      const int stripIDFront = itFront->first;
      // Extracting the points in the edges
      const double yatz0Front = _geometry->getStripPosY(layerID, 1, stripIDFront, 0.0);
      const double yatzLFront = _geometry->getStripPosY(layerID, 1, stripIDFront, _geometry->getSensorLength(layerID));

      // Checking with the Rear sensors
      for (std::vector<StripClusterPair>::iterator itRear = it->second[STRIPREAR].begin();
           itRear != it->second[STRIPREAR].end(); ++itRear) {
        const int stripIDRear = itRear->first;
        // Extracting the points in the edges and put them in a common
        // reference system (Front petal local ref., remember the y=0
        // of one system is SensorWidthMax distance from the other)
        const double yatz0Rear =
            _geometry->getSensorWidthMax(layerID) - _geometry->getStripPosY(layerID, 3, stripIDRear, 0.0);
        const double yatzLRear = _geometry->getSensorWidthMax(layerID) -
                                 _geometry->getStripPosY(layerID, 3, stripIDRear, _geometry->getSensorLength(layerID));

        // If not crossing  continue
        if ((yatz0Front - yatz0Rear) * (yatzLFront - yatzLRear) > 0.) {
          continue;
        }

        // There are intersection, so store the Hit
        // Cluster in Front
        StripCluster* pclusterFront = itFront->second;
        // Cluster in Rear
        StripCluster* pclusterRear = itRear->second;

        // Get the intersection point (note that the returning value
        // is placed in the FRONT sensor and in the x_local = thick/2
        CLHEP::Hep3Vector position = _geometry->getCrossLinePoint(layerID, ladderID, stripIDFront, stripIDRear);

        // FIXME: Sigma calculation!??
        const double sigmaY = sqrt(pclusterFront->getPosY() * pclusterFront->getPosY() +
                                   pclusterRear->getPosY() * pclusterRear->getPosY());
        // The Z indetermination: putting z in the middle of the disk
        // (i.e. middle position between front and rear sensors)
        // FIXME---> Hay que ver que es esto ???!
        const double sigmaZ = _geometry->getSensorThick(layerID);
        Hep3Vector posSigma(_geometry->getSensorThick(layerID) / 2., sigmaY, sigmaZ);

        double totalCharge = (pclusterFront->getCharge() + pclusterRear->getCharge()) / 2.0;

        // Always postioned as it was in the front
        StripCluster* pCluster3D = new StripCluster(layerID, ladderID, 1, position, posSigma, totalCharge, 0);

        // Update MC true info for 3D
        SimTrackerHitMap cls3DSimHitMap = pclusterFront->getSimHitMap();

        if (cls3DSimHitMap.size() != 0) {
          SimTrackerHitMap clsRearSimHitMap = pclusterRear->getSimHitMap();
          for (SimTrackerHitMap::iterator iterSHM = clsRearSimHitMap.begin(); iterSHM != clsRearSimHitMap.end();
               iterSHM++) {
            EVENT::SimTrackerHit* simHit = dynamic_cast<EVENT::SimTrackerHit*>(iterSHM->first);
            float weight = iterSHM->second;
            if (cls3DSimHitMap.find(simHit) != cls3DSimHitMap.end()) {
              cls3DSimHitMap[simHit] += weight;
            } else {
              cls3DSimHitMap[simHit] = weight;
            }
          }
          pCluster3D->updateSimHitMap(cls3DSimHitMap);
        }
        //
        // if defined ROOT-OUTPUT, save info
#ifdef ROOT_OUTPUT
        // Set layer, ladder, sensor ID
        _rootLayerID = layerID;
        _rootLadderID = ladderID;
        _rootSensorID = 0;

        // Set reconstructed position (local frame)
        //_rootRecRPhi     = pCluster3D->getPosY() / mm;
        //_rootRecZ        = pCluster3D->getPosZ() / mm;
        // Set cluster sizes
        _rootClsSizeRPhi = pclusterFront->getSize();
        _rootClsSizeZ = pclusterRear->getSize();

        // Set simulated position in Front & MC particle

        // Find SimTrackerHit with highest weight
        EVENT::SimTrackerHit* simHit = 0;
        float weight = 0;

        const SimTrackerHitMap& simHitMapFront = pclusterFront->getSimHitMap();
        for (SimTrackerHitMap::const_iterator iterSHM = simHitMapFront.begin(); iterSHM != simHitMapFront.end();
             ++iterSHM) {
          // Find contribution with highest weight
          if ((iterSHM->first != 0) && ((iterSHM->second) > weight)) {
            simHit = iterSHM->first;
            weight = iterSHM->second;
          }
        }

        // SimHit global position
        Hep3Vector simPosGlob;
        Hep3Vector simPosLoc;

        simPosGlob.setX(simHit->getPosition()[0] * mm);
        simPosGlob.setY(simHit->getPosition()[1] * mm);
        simPosGlob.setZ(simHit->getPosition()[2] * mm);
        simPosLoc = _geometry->transformPointToLocal(layerID, ladderID, 1, simPosGlob);
        _rootMCPDGRPhi = simHit->getMCParticle()->getPDG();
        _rootSimRPhi = sqrt(simPosGlob.getY() * simPosGlob.getY() + simPosGlob.getX() * simPosGlob.getX()) / mm;
        _rootSimZ = simPosGlob.getZ() / mm;

        // Residuals: ref. frame difined in the simHit-- rPhi, r
        const double theta = atan2(simPosGlob.getY(), simPosGlob.getX());
        //-- getting the reconstructed hit to the global ref.
        CLHEP::Hep3Vector recPoint = _geometry->transformPointToGlobal(layerID, ladderID, 1, position);
        _rootRecRPhi = sqrt(recPoint.getX() * recPoint.getX() + recPoint.getY() * recPoint.getY()) / mm;
        _rootRecZ = recPoint.getZ() / mm;

        recPoint -= simPosGlob;
        recPoint.rotateZ(theta); // FIXME: Esto es correcto??
        // Extracting the residuals: In the local frame
        // Sim position - rec position
        _rootResRPhi = (simPosLoc.getY() - position.getY()) / mm;
        _rootResR = (simPosLoc.getZ() - position.getZ()) / mm;
        _rootResModule = sqrt(_rootResRPhi * _rootResRPhi + _rootResR * _rootResR) / mm;

        // FIXME: Que hago con esta doble informacion???
        /*			      // Set simulated position in Rear & MC particle

                                      // Find SimTrackerHit with highest weight
                                      const SimTrackerHitMap & simHitMapRear =
                                              pclusterRear->getSimHitMap();

                                      for(SimTrackerHitMap::const_iterator iterSHM=
                                                      simHitMapRear.begin();
                                                      iterSHM!=simHitMapRear.end(); ++iterSHM)
                                      {
                                              // Find contribution with highest weight
                                              if( (iterSHM->first!=0) && ((iterSHM->second)>weight) )
                                              {
                                                      simHit = iterSHM->first;
                                                      weight = iterSHM->second;
                                              }
                                      }

                                      // SimHit global position
                                      simPosGlob.setX(simHit->getPosition()[0]*mm);
                                      simPosGlob.setY(simHit->getPosition()[1]*mm);
                                      simPosGlob.setZ(simHit->getPosition()[2]*mm);
                                      simPosLoc = _geometry->transformPointToLocal(layerID, ladderID,
                                                      3, simPosGlob);
                                      _rootMCPDGZ = simHit->getMCParticle()->getPDG();
                                      _rootSimZ   = simPosLoc.getZ() / mm;*/

        // Fill the tree
        _rootFile->cd("");
        _rootTree->Fill();

#endif

        clsVec.push_back(pCluster3D);
      } // for rear sensors
    } // for front sensors
  } // For cluster strip vector pairs

  /*	for(SensorStripClusterMap::iterator it = clsvectFrontRear.begin();
                          it != clsvectFrontRear.end(); ++it)
          {
                  std::cout << "========================================== " << std::endl;
                  std::cout << "CELLID: " << it->first << std::endl;
                  for(std::map<StripType,std::vector<StripClusterPair> >::iterator itM
                                  = it->second.begin(); itM != it->second.end(); ++itM)
                  {
                          std::cout << "  Strip Type: " << itM->first << std::endl;
                          for(std::vector<StripClusterPair>::iterator itP = itM->second.begin();
                                          itP != itM->second.end(); ++itP)
                          {
                                  std::cout << "    Cluster Info " << std::endl;
                                  std::cout << "       - Position[mm]:"
                                          << itP->second->get3Position()/mm << std::endl;
                                  std::cout << "       - Charge      :"
                                          << itP->second->getCharge() << std::endl;
                          }
                  }
          }
          */

  // Release memory
  for (SensorStripClusterMap::iterator it = clsvectFrontRear.begin(); it != clsvectFrontRear.end(); ++it) {
    for (std::map<StripType, std::vector<StripClusterPair>>::iterator itMap = it->second.begin();
         itMap != it->second.end(); ++itMap) {
      for (std::vector<StripClusterPair>::iterator itSCl = itMap->second.begin(); itSCl != itMap->second.end();
           ++itSCl) {
        if (itSCl->second != 0) {
          delete itSCl->second;
        }
      }
    }
  }

  return clsVec;
}

//-- Right neighbours
/*			adjCharge = 0.0;
                        goNextStrip = true;
                        while( goNextStrip && LschMap != iterSMap->second[STRIPTYPE].end() )
                        {
                              adjCharge = LschMap->getCharge();
                              const SimTrackerHitMap & adjSimHitMap = LschMap->getSimHitMap();
std::cout << "  Strip: " << LschMap->first<< " -- carga:" << adjCharge << "threshold:" <<
        _SNadjacent*_CMSnoise<< " " ;
                              // Charge higher than threshold set
                              if( adjCharge >= (_SNadjacent*_CMSnoise) )
                              {
                                      clsStrips[leftStrip] = 	new Signal(adjCharge, 0.);
                                      clsStrips[leftStrip]->updateSimHitMap(adjSimHitMap);
                                      // Set this charge as zero to avoid
                                      // double counting
                                      iterSMap->second[STRIPZ][leftStrip]->setCharge(0.);
                                      // And go on to the next strip
                                      ++LschMap;
std::cout << " ---- OK!!" << std::endl;
                              }
                              else  // Charge lower - stop searching
                              {
                                      goNextStrip = false;
                              }
                        }*/

//
// Cluster vector iterators
//---ClsVec::iterator iterClsVec, iterClsVec2;

// Initiate variables - layerID, ladderID, sensorID, stripID, seed strip, seed charge
/*short int layerID   = 0;
short int ladderID  = 0;
short int sensorID  = 0;*/

// int      seedStrip  = 0;
// double   seedCharge = 0;

// Bunch of strips forming cluster
// StripChargeMap clsStrips;

//
// Search complete sensor map - find seeds & their neighbouring strips
/*	for(SensorStripMap::iterator iterSMap=sensorMap.begin();
                        iterSMap!=sensorMap.end(); iterSMap++)
        {
                // Save layer ID , ...
                int cellID = iterSMap->first;

                std::map<std::string,int> bfmap = _geometry->decodeCellID(cellID);
                layerID = bfmap["layer"];
                ladderID= bfmap["module"];
                sensorID= bfmap["sensor"]; // NOTA QUE YA NO TIENE SENTIDO
                // SI FRONT --> sensorID = 1
                // SI REAR  --> sensorID = 3

                //
                // Clusters in Z
                for(StripChargeMap::iterator iterChMap=iterSMap->second[STRIPZ].begin();
                                iterChMap!=iterSMap->second[STRIPZ].end(); iterChMap++)
                {
                        sensorID = 3;
                        //
                        // Zero: Save new candidate for seed strip + MC true info
                        int seedStrip  = iterChMap->first;
                        double seedCharge = iterChMap->second->getCharge();

                        const SimTrackerHitMap & seedSimHitMap =
                                iterChMap->second->getSimHitMap();


                        if( seedCharge < (_SNseed*_CMSnoise) )
                        {
                                continue;
                        }
                        //
                        // First: New cluster and its seed strip has been found

                        clsStrips[seedStrip] = new Signal(seedCharge, 0.);
                        clsStrips[seedStrip]->updateSimHitMap(seedSimHitMap);

                        // Set this charge as zero to avoid double counting
                        iterChMap->second->setCharge(0.);

                        // Continue searching - find left and right neighbours
                        // If floating strips exist --> go to every even strip
                        int leftStrip  = seedStrip -1;
                        int rightStrip = seedStrip +1;

                        if(_floatStripsZ)
                        {
                                leftStrip--;
                                rightStrip--;
                        }

                        double adjCharge = 0;

                        bool goLeft      = true;
                        bool goRight     = true;

                        // Control if left strip and right strip are in range
                        if(leftStrip  <  0)
                        {
                                goLeft  = false;
                        }
//			if(rightStrip >= _geometry->getSensorNStripsInZ(layerID))
                        if(rightStrip >= _geometry->getSensorNStrips(layerID,sensorID)) //FIXME: >= o solo > ???
                        {
                                goRight = false;
                        }

                        // ---> Abajo codigo repetido
                        //      goleft, leftStrip , siguiente strip -
                        //      gorigth, rightStrip, siguiente strip +
                        // ---> Puede hacerse tambien aprovechando que el mapa ordena
                        //      de menor a mayor (o al reves no me acuerdo..)
                        //      podemos pues ir de izq a derecha y de derecha a izq usando ++ -- en el iter y parando en
rend() end()
                        // Second: search for left neighbours
                        while( goLeft && ((iterSMap->second[STRIPZ].find(leftStrip))
                                                !=iterSMap->second[STRIPZ].end()) )
                        {
                                adjCharge = iterSMap->second[STRIPZ][leftStrip]->getCharge();

                                const SimTrackerHitMap & adjSimHitMap =
iterSMap->second[STRIPZ][leftStrip]->getSimHitMap();
                                // Charge higher than threshold set
                                if( adjCharge >= (_SNadjacent*_CMSnoise) )
                                {
                                        clsStrips[leftStrip] = 	new Signal(adjCharge, 0.);
                                        clsStrips[leftStrip]->updateSimHitMap(adjSimHitMap);
                                        // Set this charge as zero to avoid
                                        // double counting
                                        iterSMap->second[STRIPZ][leftStrip]->setCharge(0.);
                                        // Go on to the next strip
                                        if(_floatStripsZ)
                                        {
                                                leftStrip -= 2;
                                        }
                                        else
                                        {
                                                leftStrip -= 1;
                                        }
                                }
                                // Charge lower - stop searching
                                else
                                {
                                        goLeft = false;
                                }
                        }

                        //
                        // Third: search for right neighbours
                        while(goRight && ((iterSMap->second[STRIPZ].find(rightStrip))
                                                !=iterSMap->second[STRIPZ].end()) )
                        {
                                adjCharge = iterSMap->second[STRIPZ][rightStrip]->getCharge();
                                const SimTrackerHitMap & adjSimHitMap =
iterSMap->second[STRIPZ][rightStrip]->getSimHitMap();
                                // Charge higher than threshold set
                                if( adjCharge >= (_SNadjacent*_CMSnoise) )
                                {
                                        clsStrips[rightStrip] = new Signal(adjCharge, 0.);
                                        clsStrips[rightStrip]->updateSimHitMap(adjSimHitMap);
                                        // Set this charge as zero to avoid double counting
                                        iterSMap->second[STRIPZ][rightStrip]->setCharge(0.);

                                        // Go on to the next strip
                                        if(_floatStripsZ)
                                        {
                                                rightStrip += 2;
                                        }
                                        else
                                        {
                                                rightStrip += 1;
                                        }
                                }
                                // Charge lower - stop searching
                                else
                                {
                                        goRight = false;
                                }
                        }

                        //
                        // Fourth: Calculate mean position of a new cluster
                        SimTrackerHitMap clsZSimHitMap;
                        SimTrackerHitMap::const_iterator iterSHM;

                        // Cluster: position, charge & size
                        short int clsSizeZ   = clsStrips.size();
                        double    clsPosZ    = 0.;
                        double    clsChargeZ = 0.;

                        double xLeftSignal   = 0 ;
                        double qLeftSignal   = 0 ;

                        double xRightSignal  = 0 ;
                        double qRightSignal  = 0 ;

                        double qIntermSignal = 0 ;

                        for(StripChargeMap::iterator iterChMap2=clsStrips.begin();
                                        iterChMap2!=clsStrips.end(); iterChMap2++)
                        {
                                // Current strip ID, posZ & charge
                                int stripID        = iterChMap2->first;
                                double  stripPosZ   = _geometry->getStripPosY(layerID, sensorID,stripID,0.0);
                                double stripCharge = iterChMap2->second->getCharge();

                                // Update info about MC particles which contributed
                                const SimTrackerHitMap & simHitMap = iterChMap2->second->getSimHitMap();

                                if(simHitMap.size() != 0)
                                {
                                        for(iterSHM=simHitMap.begin();
                                                        iterSHM!=simHitMap.end(); iterSHM++)
                                        {
                                                EVENT::SimTrackerHit * simHit = dynamic_cast<EVENT::SimTrackerHit
*>(iterSHM->first); float                  weight = iterSHM->second;

                                                if(clsZSimHitMap.find(simHit)!=clsZSimHitMap.end())
                                                {
                                                        clsZSimHitMap[simHit] += weight;
                                                }
                                                else
                                                {
                                                        clsZSimHitMap[simHit]  = weight;
                                                }
                                        }
                                }

                                // Get last iterator
                                StripChargeMap::iterator iterChMapLast = clsStrips.end();
                                -- iterChMapLast;

                                // Get leftmost signal
                                if (iterChMap2 == clsStrips.begin())
                                {
                                        xLeftSignal = stripPosZ;
                                        qLeftSignal = stripCharge;
                                }

                                // Get rightmost signal
                                else if (iterChMap2 == iterChMapLast)
                                {
                                        xRightSignal = stripPosZ;
                                        qRightSignal = stripCharge;
                                }
                                // Get intermediate signal
                                else
                                {
                                        qIntermSignal += stripCharge;
                                }

                                // Update total charge
                                clsChargeZ += stripCharge;

                        }

                        // Get average intermediate signal
                        if (clsSizeZ > 2)
                        {
                                qIntermSignal /= (clsSizeZ - 2);
                        }
                        else
                        {
                                qIntermSignal  = 0.;
                        }
                        // Analog head-tail algorithm
        //		double geomPitchInZ    = _geometry->getSensorPitchInZ(layerID);
                        double geomPitchInZ    = _geometry->getSensorPitch(layerID,sensorID,0);
                        double readOutPitchInZ = 0.;
                        if (_floatStripsZ)
                        {
                                readOutPitchInZ = 2*geomPitchInZ;
                        }
                        else
                        {
                                readOutPitchInZ =   geomPitchInZ;
                        }

                        // qIntermSignal >= 1/eps*qRighSignal || 1/eps*qLeftSignal --> if not don't use head-tail (delta
electron ...) - use epsilon ~ 1.5 if ( (clsSizeZ>2) && (qLeftSignal<1.5*qIntermSignal) &&
(qRightSignal<1.5*qIntermSignal) )
                        {
                                clsPosZ = (xRightSignal + xLeftSignal)/2. + (qRightSignal -
qLeftSignal)/2./qIntermSignal * readOutPitchInZ;
                        }
                        // COG algorithm
                        else
                        {
                                clsPosZ = (xRightSignal*qRightSignal + xLeftSignal*qLeftSignal)/(qRightSignal +
qLeftSignal);
                        }

                        //
                        // Fifth: Correct cluster position to Lorentz shift (Bz is expected non-zero only --> no
correction) clsPosZ += 0.;

                        //
                        // Sixth: Save information about new cluster (cluster in Z)
                        if (clsChargeZ >= (_SNtotal*_CMSnoise))
                        {
                                StripCluster * pClusterZ = new StripCluster( layerID, ladderID, sensorID,
Hep3Vector(_geometry->getSensorThick(layerID)/2., 0., clsPosZ), Hep3Vector(_geometry->getSensorThick(layerID)/2., 0.,
0.), clsChargeZ, clsSizeZ ); pClusterZ->updateSimHitMap(clsZSimHitMap);

                                clsVecInZ.push_back(pClusterZ);
                        }

                        // Release memory
                        for (StripChargeMap::iterator iterChMap2=clsStrips.begin(); iterChMap2!=clsStrips.end();
iterChMap2++)
                        {
                                delete iterChMap2->second;
                        }
                        // Clear content
                        clsStrips.clear();

                } // For clusters in Z

                //
                // Clusters in R-Phi + 3D clusters utilizing Z (use Z position of "clusters in Z" and create 3D clusters
from them) for(StripChargeMap::iterator iterChMap=iterSMap->second[STRIPRPHI].begin();
iterChMap!=iterSMap->second[STRIPRPHI].end(); iterChMap++) {

         sensorID = 1;

         //
         // Zero: Save new candidate for seed strip + MC true info
         int seedStrip  = iterChMap->first;
         double seedCharge = iterChMap->second->getCharge();

         const SimTrackerHitMap & seedSimHitMap = iterChMap->second->getSimHitMap();

         //
         // First: New cluster and its seed strip has been found
         if ( seedCharge >= (_SNseed*_CMSnoise) ) {

            clsStrips[seedStrip] = new Signal(seedCharge, 0.);
            clsStrips[seedStrip]->updateSimHitMap(seedSimHitMap);

            // Set this charge as zero to avoid double counting
            iterChMap->second->setCharge(0.);

            // Continue searching - find left and right neighbours
            bool goLeft      = true;
            bool goRight     = true;

            // If floating strips exist --> go to every even strip
            int leftStrip  = -1;
            int rightStrip = -1;

            if (_floatStripsRPhi) {

               leftStrip    = seedStrip - 2;
               rightStrip   = seedStrip + 2;
            }
            else {

               leftStrip    = seedStrip - 1;
               rightStrip   = seedStrip + 1;
            }

            double adjCharge = 0;

            // Control if left strip and right strip are in range
            if (leftStrip  <  0)                                          goLeft  = false;
//            if (rightStrip >= _geometry->getSensorNStripsInRPhi(layerID)) goRight = false;
            if (rightStrip >= _geometry->getSensorNStrips(layerID,sensorID)) goRight = false;

            //
            // Second: Search for left neighbours
            while ( goLeft && ((iterSMap->second[STRIPRPHI].find(leftStrip))!=iterSMap->second[STRIPRPHI].end()) ) {

               adjCharge   = iterSMap->second[STRIPRPHI][leftStrip]->getCharge();

               const SimTrackerHitMap & adjSimHitMap = iterSMap->second[STRIPRPHI][leftStrip]->getSimHitMap();

               // Charge higher than threshold set
               if ( adjCharge >= (_SNadjacent*_CMSnoise) ) {

                  clsStrips[leftStrip] = new Signal(adjCharge, 0.);
                  clsStrips[leftStrip]->updateSimHitMap(adjSimHitMap);

                  // Set this charge as zero to avoid double counting
                  iterSMap->second[STRIPRPHI][leftStrip]->setCharge(0.);

                  // Go on to the next strip
                  if (_floatStripsRPhi) leftStrip -= 2;
                  else                  leftStrip -= 1;
               }
               // Charge lower - stop searching
               else goLeft = false;
            }

            //
            // Third: Search for right neighbours
            while (goRight && ((iterSMap->second[STRIPRPHI].find(rightStrip))!=iterSMap->second[STRIPRPHI].end()) ) {

               adjCharge = iterSMap->second[STRIPRPHI][rightStrip]->getCharge();

               const SimTrackerHitMap & adjSimHitMap = iterSMap->second[STRIPRPHI][rightStrip]->getSimHitMap();

               // Charge higher than threshold set
               if ( adjCharge >= (_SNadjacent*_CMSnoise) ) {

                  clsStrips[rightStrip] = new Signal(adjCharge, 0.);
                  clsStrips[rightStrip]->updateSimHitMap(adjSimHitMap);

                  // Set this charge as zero to avoid double counting
                  iterSMap->second[STRIPRPHI][rightStrip]->setCharge(0.);

                  // Go on to the next strip
                  if (_floatStripsRPhi) rightStrip += 2;
                  else                  rightStrip += 1;
               }
               // Charge lower - stop searching
               else goRight = false;
            }

            //
            // Fourth: Go through all Z clusters and calculate mean position in R-Phi and save them together as a 3D
cluster for (iterClsVec=clsVecInZ.begin(); iterClsVec!=clsVecInZ.end(); iterClsVec++) {

               // Get and save info about cluster in Z
               StripCluster * pClusterZ    = *iterClsVec;

               const SimTrackerHitMap & clsZSimHitMap = pClusterZ->getSimHitMap();

               // Cluster: position, charge & size
               short int clsSizeRPhi   = clsStrips.size();
               double    clsPosRPhi    = 0.;
               double    clsChargeRPhi = 0.;

               double xLeftSignal   = 0 ;
               double qLeftSignal   = 0 ;

               double xRightSignal  = 0 ;
               double qRightSignal  = 0 ;

               double qIntermSignal = 0 ;

               SimTrackerHitMap clsRPhiSimHitMap;
               SimTrackerHitMap::const_iterator iterSHM;

               // Find clusters in R-Phi
               for (StripChargeMap::iterator iterChMap2=clsStrips.begin(); iterChMap2!=clsStrips.end(); iterChMap2++) {

                  // Current strip ID, posRPhi & charge
                  int stripID         = iterChMap2->first;
                  //double stripPosRPhi = _geometry->getStripPosInRPhi(layerID, stripID, pClusterZ->getPosZ());
                  double stripPosRPhi = _geometry->getStripPosY(layerID, sensorID, stripID,pClusterZ->getPosZ());
                  double stripCharge  = iterChMap2->second->getCharge();

                  // Update info about MC particles which contributed
                  const SimTrackerHitMap & simHitMap = iterChMap2->second->getSimHitMap();

                  if (simHitMap.size() != 0) {

                     for (iterSHM=simHitMap.begin(); iterSHM!=simHitMap.end(); iterSHM++) {

                        EVENT::SimTrackerHit * simHit = dynamic_cast<EVENT::SimTrackerHit *>(iterSHM->first);
                        float                  weight = iterSHM->second;

                        if (clsRPhiSimHitMap.find(simHit)!=clsRPhiSimHitMap.end()) clsRPhiSimHitMap[simHit] += weight;
                        else                                                       clsRPhiSimHitMap[simHit]  = weight;
                     }
                  }

                  // Get last iterator
                  StripChargeMap::iterator iterChMapLast = clsStrips.end();
                   -- iterChMapLast;

                   // Get leftmost signal
                   if (iterChMap2 == clsStrips.begin()) {

                      xLeftSignal = stripPosRPhi;
                      qLeftSignal = stripCharge;
                   }

                   // Get rightmost signal
                   else if (iterChMap2 == iterChMapLast) {

                      xRightSignal = stripPosRPhi;
                      qRightSignal = stripCharge;
                   }
                   // Get intermediate signal
                   else qIntermSignal += stripCharge;

                   // Update total charge
                   clsChargeRPhi += stripCharge;

               }

               // Get average intermediate signal
               if (clsSizeRPhi > 2) qIntermSignal /= (clsSizeRPhi - 2);
               else                 qIntermSignal  = 0.;

               // Analog head-tail algorithm
               //double geomPitchInRPhi    = _geometry->getSensorPitchInRPhi(layerID, pClusterZ->getPosZ());
               double geomPitchInRPhi    = _geometry->getSensorPitch(layerID,sensorID, pClusterZ->getPosZ());
               double readOutPitchInRPhi = 0.;
               if (_floatStripsRPhi) readOutPitchInRPhi = 2*geomPitchInRPhi;
               else                  readOutPitchInRPhi =   geomPitchInRPhi;

               // qIntermSignal >= 1/eps*qRighSignal || 1/eps*qLeftSignal --> if not don't use head-tail (delta electron
...) - use epsilon ~ 1.5 if ( (clsSizeRPhi>2) && (qLeftSignal<1.5*qIntermSignal) && (qRightSignal<1.5*qIntermSignal) ) {

                    clsPosRPhi = (xRightSignal + xLeftSignal)/2. + (qRightSignal - qLeftSignal)/2./qIntermSignal *
readOutPitchInRPhi;
               }
               // COG algorithm
               else clsPosRPhi = (xRightSignal*qRightSignal + xLeftSignal*qLeftSignal)/(qRightSignal + qLeftSignal);


               // Fifth: Correct mean position to Lorentz shift (Bz is expected non-zero only; correct to theta angle of
each sensor) clsPosRPhi += _TanOfAvgELorentzShift * _geometry->getSensorThick(layerID)/2. *
cos(_geometry->getLadderTheta(layerID));

               //
               // Sixth: Save information about new cluster (cluster in RPhi), if cluster charge higher than threshold
set if (clsChargeRPhi >= (_SNtotal*_CMSnoise)) {

                  StripCluster * pClusterRPhi = new StripCluster( layerID, ladderID, sensorID,
Hep3Vector(_geometry->getSensorThick(layerID)/2., clsPosRPhi, 0.), Hep3Vector(_geometry->getSensorThick(layerID)/2., 0.,
0.), clsChargeRPhi, clsSizeRPhi );

                  pClusterRPhi->updateSimHitMap(clsRPhiSimHitMap);

                  clsVecInRPhi.push_back(pClusterRPhi);

                  // Create final 3D cluster
                  Hep3Vector position( _geometry->getSensorThick(layerID)/2., pClusterRPhi->getPosY()     ,
pClusterZ->getPosZ()); Hep3Vector posSigma( _geometry->getSensorThick(layerID)/2., pClusterRPhi->getPosSigmaY(),
pClusterZ->getPosSigmaZ());

                  double     totalCharge = (pClusterRPhi->getCharge() + pClusterZ->getCharge() )/2.;

                  StripCluster * pCluster3D = new StripCluster( layerID, ladderID, sensorID, position, posSigma,
totalCharge, 0);

                  SimTrackerHitMap cls3DSimHitMap = clsRPhiSimHitMap;

                  // Update MC true info for 3D
                  if (cls3DSimHitMap.size() != 0) {

                     for (iterSHM=clsZSimHitMap.begin(); iterSHM!=clsZSimHitMap.end(); iterSHM++) {

                        EVENT::SimTrackerHit * simHit = dynamic_cast<EVENT::SimTrackerHit *>(iterSHM->first);
                        float                  weight = iterSHM->second;

                        if (cls3DSimHitMap.find(simHit)!=cls3DSimHitMap.end()) cls3DSimHitMap[simHit] += weight;
                        else                                                   cls3DSimHitMap[simHit]  = weight;
                     }

                     pCluster3D->updateSimHitMap(cls3DSimHitMap);
                  }

                  //
                  // if defined ROOT-OUTPUT, save info
#ifdef ROOT_OUTPUT

                  // Set layer, ladder, sensor ID
                  _rootLayerID     = layerID;
                  _rootLadderID    = ladderID;
                  _rootSensorID    = sensorID;

                  // Set reconstructed positions
                  _rootRecRPhi     = pClusterRPhi->getPosY() / mm;
                  _rootRecZ        = pClusterZ->getPosZ()    / mm;

                  // Set cluster sizes
                  _rootClsSizeRPhi = pClusterRPhi->getSize();
                  _rootClsSizeZ    = pClusterZ->getSize();

                  // Set simulated position in RPhi & MC particle

                  // Find SimTrackerHit with highest weight
                  EVENT::SimTrackerHit * simHit = 0;
                  float                  weight = 0;

                  const SimTrackerHitMap & simHitMapRPhi = pClusterRPhi->getSimHitMap();

                  for (iterSHM=simHitMapRPhi.begin(); iterSHM!=simHitMapRPhi.end(); iterSHM++) {

                     // Find contribution with highest weight
                     if ( (iterSHM->first!=0) && ((iterSHM->second)>weight) ) {

                        simHit = iterSHM->first;
                        weight = iterSHM->second;
                     }
                  }

                  // SimHit global position
                  Hep3Vector simPosGlob;
                  Hep3Vector simPosLoc;

                  simPosGlob.setX(simHit->getPosition()[0]*mm);
                  simPosGlob.setY(simHit->getPosition()[1]*mm);
                  simPosGlob.setZ(simHit->getPosition()[2]*mm);
                  simPosLoc = _geometry->transformPointToLocal(layerID, ladderID, sensorID, simPosGlob);

                  _rootMCPDGRPhi = simHit->getMCParticle()->getPDG();
                  _rootSimRPhi   = simPosLoc.getY() / mm;

                  // Set simulated position in Z & MC particle

                  // Find SimTrackerHit with highest weight
                  const SimTrackerHitMap & simHitMapZ = pClusterZ->getSimHitMap();

                  for (iterSHM=simHitMapZ.begin(); iterSHM!=simHitMapZ.end(); iterSHM++) {

                     // Find contribution with highest weight
                     if ( (iterSHM->first!=0) && ((iterSHM->second)>weight) ) {

                        simHit = iterSHM->first;
                        weight = iterSHM->second;
                     }
                  }

                  // SimHit global position
                  simPosGlob.setX(simHit->getPosition()[0]*mm);
                  simPosGlob.setY(simHit->getPosition()[1]*mm);
                  simPosGlob.setZ(simHit->getPosition()[2]*mm);
                  simPosLoc = _geometry->transformPointToLocal(layerID, ladderID, sensorID, simPosGlob);

                  _rootMCPDGZ = simHit->getMCParticle()->getPDG();
                  _rootSimZ   = simPosLoc.getZ() / mm;

                  // Fill the tree
                  _rootFile->cd("");
                  _rootTree->Fill();

#endif

                  clsVec.push_back(pCluster3D);
               }

            } // Go through all Z clusters

            // Release memory
            for (StripChargeMap::iterator iterChMap2=clsStrips.begin(); iterChMap2!=clsStrips.end(); iterChMap2++)
delete iterChMap2->second;

            // Clear content
            clsStrips.clear();

         } // If found new cluster

      } // For clusters in R-Phi

      //
      // Release memory
      for (iterClsVec=clsVecInRPhi.begin(); iterClsVec!=clsVecInRPhi.end(); iterClsVec++) delete *iterClsVec;
      for (iterClsVec=clsVecInZ.begin();    iterClsVec!=clsVecInZ.end();    iterClsVec++) delete *iterClsVec;

      // Clear
      clsVecInRPhi.clear();
      clsVecInZ.clear();

   } // For sensor map*/

// OTHER METHODS

//
// Calculated and stored the hits adjacents
// Template for the reverse_iterator (leftStrips) and iterator (rightStrips)
//

template <class It>
StripChargeMap& SiStripClus::storeHitsAdjacents(StripChargeMap& clsStrips, It schMap, const It& endIt,
                                                StripChargeMap& currentMap) {
  // Map ordered from lower to higher
  double adjCharge = 0.0;
  bool goNextStrip = true;
  while (goNextStrip && schMap != endIt) {
    adjCharge = schMap->second->getCharge();
    const SimTrackerHitMap& adjSimHitMap = schMap->second->getSimHitMap();
    // Charge higher than threshold set
    if (adjCharge >= (_SNadjacent * _CMSnoise)) {
      clsStrips[schMap->first] = new Signal(adjCharge, 0.);
      clsStrips[schMap->first]->updateSimHitMap(adjSimHitMap);
      // Set this charge as zero to avoid
      // double counting
      currentMap[schMap->first]->setCharge(0.);
      // And go on to the next strip
      ++schMap;
    } else // Charge lower - stop searching
    {
      goNextStrip = false;
    }
  }

  return clsStrips;
}

//
// Method calculating hits from given clusters
//
// FIXME: Returning colOfTrkHits, clsVec reference?? --> the clsVec.clear
void SiStripClus::calcHits(ClsVec& clsVec, IMPL::LCCollectionVec* colOfTrkHits) {
  // FIXME: Why passing clsVec by reference?
  //        return the LCCol
  // Cluster - position, covariance matrix (3x3 = 6 parameters = lower triangle matrix)
  ClsVec::iterator iterClsVec;

  short int layerID;
  short int ladderID;
  short int sensorID;

  Hep3Vector position;
  Hep3Vector posSigma;

  double totalCharge;

  // Print sensor info
  streamlog_out(MESSAGE2) << std::endl
                          << "   Total number of reconstructed hits: " << clsVec.size() << " hit(s)" << std::endl;

  // Set collection flag - cellID 1 will be stored,
  // and LCrelations
  LCFlagImpl flag1(0);
  flag1.setBit(LCIO::RTHPBIT_ID1);
  colOfTrkHits->setFlag(flag1.getFlag());

  // CODIFICATION --->FIXME: METHOD IN GEAR (Centralizing...)
  CellIDEncoder<TrackerHitPlaneImpl> cellEnc(LCTrackerCellID::encoding_string() + ",stripFront:11,stripRear:11",
                                             colOfTrkHits);

  // Go through all clusters
  for (iterClsVec = clsVec.begin(); iterClsVec != clsVec.end(); iterClsVec++) {
    // Get cluster info
    StripCluster* pCluster = *iterClsVec;

    layerID = pCluster->getLayerID();
    ladderID = pCluster->getLadderID();
    sensorID = pCluster->getSensorID();

    position = pCluster->get3Position();
    posSigma = pCluster->get3PosSigma();

    totalCharge = pCluster->getCharge();

    CLHEP::Hep3Vector localPos = position;

    // Transform hit and covariance to the global ref. system
    position = _geometry->transformPointToGlobal(layerID, ladderID, sensorID, position);
    // Reasigning the z-position to be placed in the middle between the
    // front and rear sensor
    const int sign = abs(_geometry->getLayerRealID(layerID)) / _geometry->getLayerRealID(layerID);
    const double offsetZ = _geometry->getLadderThick(layerID) / 2. + _geometry->getSensorThick(layerID) / 2.;
    position.setZ(position.getZ() + sign * offsetZ);

    // Save into LCIO native variables and in appropriate units
    double posLCIO[3] = {position.getX() / mm, position.getY() / mm, position.getZ() / mm};

    // Calculate resolution
    float* covLCIO = calcResolution(layerID, position.theta() / pi * 180., localPos.getZ());

    // Create new Tracker hit
    TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl();

    // Set hit type: SVD type is 201 + layer ID ??
    trkHit->setType(layerID + 201);
    trkHit->setPosition(posLCIO);
    // U and V error, not using covMatrix anymore
    // trkHit->setCovMatrix(covLCIO);
    trkHit->setdU(covLCIO[2]);
    trkHit->setdV(covLCIO[5]);
    trkHit->setEDep(totalCharge / e); // in electrons
    // FIXME: Still missing EDepError
    trkHit->setTime(pCluster->getTime()); // FIXME: NOT IMPLEMENTED

    // The unit vector of the local reference frame of the front sensor
    // u: rphi direction, v: radial direction
    CLHEP::Hep3Vector Uvector = _geometry->transformVecToGlobal(layerID, ladderID, 1, CLHEP::Hep3Vector(0.0, 1.0, 0.0));
    float Uf[2] = {Uvector.getX(), Uvector.getY()};
    CLHEP::Hep3Vector Vvector = _geometry->transformVecToGlobal(layerID, ladderID, 1, CLHEP::Hep3Vector(0.0, 0.0, 1.0));
    float Vf[2] = {Vvector.getX(), Vvector.getY()};
    trkHit->setU(Uf);
    trkHit->setV(Vf);

    // Set simTrkHits which contributed & find hit with highest weight
    const SimTrackerHitMap& simHitMap = pCluster->getSimHitMap();
    // float                    weightSum = pCluster->getSimHitWeightSum();
    SimTrackerHit* simTrkHit = 0;
    float weight = 0;

    for (SimTrackerHitMap::const_iterator iterSHM = simHitMap.begin(); iterSHM != simHitMap.end(); iterSHM++) {
      // Don't save "noise hits", i.e. zero pointers
      // Find contribution with highest weight
      if ((iterSHM->first != 0) && ((iterSHM->second) > weight)) {
        simTrkHit = iterSHM->first;
        weight = iterSHM->second;
      }
    }

    int realLayer = _geometry->getLayerRealID(layerID);
    // Codification FIXME: Call FTD method?
    cellEnc["subdet"] = ILDDetID::FTD;
    cellEnc["side"] = abs(realLayer) / realLayer;
    cellEnc["layer"] = abs(realLayer);
    cellEnc["module"] = ladderID + 1;
    cellEnc["sensor"] = 1;
    cellEnc["stripFront"] = pCluster->getStripFront();
    cellEnc["stripRear"] = pCluster->getStripRear();

    cellEnc.setCellID(trkHit);

    // Save only hit with highest weight
    trkHit->rawHits().push_back(simTrkHit);

    // Save the hit to the collection
    colOfTrkHits->addElement(trkHit);

    // Print infor
    printHitInfo(pCluster);

    // Release memory
    delete pCluster;

  } // For

  streamlog_out(MESSAGE2) << std::endl;

  // Clear content
  clsVec.clear();
}

//
// Method calculating hit resolution, i.e. covariance matrix
// The posZ is in local frame coordinate system
float* SiStripClus::calcResolution(const int& layerID, const double& hitTheta, const double& posZ) {
  static float covMatrix[6] = {0., 0., 0., 0., 0., 0.};
  // Assuming resolution as pitch/2
  // we need the zposlocal in the local reference system
  const double halfPitch = _geometry->getSensorPitch(layerID, 1, posZ);

  const double sigmax = halfPitch / (2. * cos(hitTheta));
  const double sigmay = halfPitch * sin(hitTheta) / 2.;

  // Set covariance matrix in appropriate units*/
  // Assuming no correlations between the front and rear sensor
  covMatrix[2] = sigmax / mm * sigmax / mm;          // Local frame Y (u-vector)
  covMatrix[5] = sigmay / mm * sigmay / mm;          // Local frame Z (v-vector)
  covMatrix[0] = (_geometry->getLadderThick(layerID) // Local frame x
                  + _geometry->getSensorThick(layerID)) /
                 2.0;

  return covMatrix;
  // Define array of theta angle
  /*static float theta[14] = {20., 30., 40., 50., 60., 70., 80., 90.,
          100., 110., 120., 130., 140., 150.};

  // Find correct bin for theta (-1 for < 20deg; 14 for > 150deg)
  int iTheta = -1;
  for(int i=0; i<14; i++)
  {
          if (hitTheta>=theta[i])
          {
                  iTheta++;
          }
  }

  // Define tracker hit resolution
  float resInZ    = 0.;
  float resInRPhi = 0.;

  //
  // Calculate resolution - interpolate between given theta angles

  // Theta is less than 20 degrees
  if (iTheta == -1)
  {
          // SVD - First layer
          if (layerID == 2)
          {
                  resInZ    = _resSVDFirstInZ[0];
                  resInRPhi = _resSVDFirstInRPhi[0];
          }

          // SVD - Other layers
          else
          {
                  resInZ    = _resSVDOtherInZ[0];
                  resInRPhi = _resSVDOtherInRPhi[0];
          }
  }
  // Theta is higher than 150 degrees
  else if (iTheta == 14)
  {
          // SVD - First layer
          if (layerID == 2)
          {
                  resInZ    = _resSVDFirstInZ[13];
                  resInRPhi = _resSVDFirstInRPhi[13];
          }

          // SVD - Other layers
          else
          {
                  resInZ    = _resSVDOtherInZ[13];
                  resInRPhi = _resSVDOtherInRPhi[13];
          }
  }

  // Interpolate between 20 and 150 degrees
  else
  {
          // SVD - First layer
          if (layerID == 2)
          {
                  resInZ    = _resSVDFirstInZ[iTheta]    +
                          (_resSVDFirstInZ[iTheta+1]    - _resSVDFirstInZ[iTheta]   )
                          /(theta[iTheta+1] - theta[iTheta])*(hitTheta - theta[iTheta]);
                  resInRPhi = _resSVDFirstInRPhi[iTheta] +
                          (_resSVDFirstInRPhi[iTheta+1] - _resSVDFirstInRPhi[iTheta])
                          /(theta[iTheta+1] - theta[iTheta])*(hitTheta - theta[iTheta]);
          }

          // SVD - Other layers
          else
          {
                  resInZ    = _resSVDOtherInZ[iTheta]    +
                          (_resSVDOtherInZ[iTheta+1]    - _resSVDOtherInZ[iTheta]   )
                          /(theta[iTheta+1] - theta[iTheta])*(hitTheta - theta[iTheta]);
                  resInRPhi = _resSVDOtherInRPhi[iTheta] +
                          (_resSVDOtherInRPhi[iTheta+1] - _resSVDOtherInRPhi[iTheta])
                          /(theta[iTheta+1] - theta[iTheta])*(hitTheta - theta[iTheta]);
          }
  }

  // Set covariance matrix in appropriate units
  covMatrix[2] = resInRPhi/mm * resInRPhi/mm;
  covMatrix[5] = resInZ/mm * resInZ/mm;

  return covMatrix;*/
}

//
// Method to save the signal
//
// FTD stuff --> Need to change the sensorMap in order to store in the
//               same StripCharge map the two opposites sensors of the same disk-petal
//               as in FTD the Hits are built using two Single side sensors of the same
//               disk-petal
//
// FIXME: Returns the SensorStripMap
void SiStripClus::updateMap(TrackerPulseImpl* pulse, SensorStripMap& sensorMap) {
  // CellID0 encodes layerID, ladderID and sensorID;
  // cellID1 encodes strip
  int cellID0 = pulse->getCellID0();
  // the new cellID0 is as the old one but extracted the sensor
  // information (cellID["sensor"] = 0
  cellID0 = _geometry->cellID0withSensor0(cellID0);

  int cellID1 = pulse->getCellID1();

  double charge = pulse->getCharge() * fC;

  // Decode stripType and stripID
  std::pair<StripType, int> tpIdpair = _geometry->decodeStripID(cellID1);

  const StripType stripType = tpIdpair.first;
  const int stripID = tpIdpair.second;

  // Controlling some errors
  if (stripType != STRIPFRONT && stripType != STRIPREAR) {
    streamlog_out(ERROR) << "cellID1 - problem to identify if strips in FRONT or REAR!!!" << std::endl;
    exit(0);
  }

  Signal* signal = 0;

  // Find if sensor already saved in map
  SensorStripMap::iterator iterSMap = sensorMap.find(cellID0);
  // If not, create map  for the sensor (cellID0)
  if (iterSMap == sensorMap.end()) {
    // Create map
    sensorMap[cellID0] = new StripChargeMap[2];
    // And update the iterator
    iterSMap = sensorMap.find(cellID0);
  }

  // Find if strip already saved in map
  StripChargeMap::iterator iterChMap = iterSMap->second[stripType].find(stripID);
  // IF not, create entry for strip, otherwise update
  if (iterChMap == iterSMap->second[stripType].end()) {
    signal = new Signal(charge, 0.0);
    iterSMap->second[stripType][stripID] = signal;
    // And update the iterator
    iterChMap = iterSMap->second[stripType].find(stripID);
  } else {
    iterChMap->second->updateCharge(charge);
  }

  // Get MCParticles which contributed and update MCParticles
  if (_navigatorPls != NULL) {
    LCObjectVec lcObjVec = _navigatorPls->getRelatedToObjects(pulse);
    FloatVec floatVec = _navigatorPls->getRelatedToWeights(pulse);

    LCObjectVec::iterator iterLCObjVec = lcObjVec.begin();
    FloatVec::iterator iterFloatVec = floatVec.begin();

    for (iterLCObjVec = lcObjVec.begin(); iterLCObjVec != lcObjVec.end(); iterLCObjVec++, iterFloatVec++) {
      EVENT::SimTrackerHit* simHit = dynamic_cast<EVENT::SimTrackerHit*>(*iterLCObjVec);
      float weight = *iterFloatVec;

      iterChMap->second->updateSimHitMap(simHit, weight);
    }
  }
  // FIXME CONTROL DE ERRORES
  // return true;
}

void SiStripClus::releaseMap(SensorStripMap& sensorMap) {
  // Release memory
  for (SensorStripMap::iterator iterSMap = sensorMap.begin(); iterSMap != sensorMap.end(); iterSMap++) {
    // Array contents
    //  Strips in Front
    for (StripChargeMap::iterator iterChMap = iterSMap->second[STRIPFRONT].begin();
         iterChMap != iterSMap->second[STRIPFRONT].end(); iterChMap++) {
      delete iterChMap->second;
    }
    // Strips in Rear
    for (StripChargeMap::iterator iterChMap = iterSMap->second[STRIPREAR].begin();
         iterChMap != iterSMap->second[STRIPREAR].end(); iterChMap++) {
      delete iterChMap->second;
    }

    // Release memory for array
    delete[] iterSMap->second;
  }

  // Clearing
  sensorMap.clear();
}

// PRINT METHODS

//
// Method printing processor parameters
//
void SiStripClus::printProcessorParams() const {
  streamlog_out(MESSAGE3) << std::endl
                          << " " << DUNDERL << DBLUE << "SiStripClus parameters:" << ENDCOLOR << " " << std::endl
                          << std::endl;

  streamlog_out(MESSAGE3) << std::setiosflags(std::ios::fixed | std::ios::internal) << std::setprecision(2)
                          << "  CMS noise [fC]:              " << std::setw(4) << _CMSnoise / fC << std::endl
                          << "                               " << std::endl
                          << std::setprecision(1) << "  S/N cut for seed strips:     " << std::setw(3) << _SNseed
                          << std::endl
                          << "  S/N cut for adjacent strips: " << std::setw(3) << _SNadjacent << std::endl
                          << "  S/N cut for total charge:    " << std::setw(3) << _SNtotal << std::endl
                          << std::resetiosflags(std::ios::showpos) << std::setprecision(0) << std::endl;

  if (_floatStripsRPhi)
    streamlog_out(MESSAGE3) << "  Read-out pitch in R-Phi is 2x geom. pitch." << std::endl;
  if (_floatStripsZ)
    streamlog_out(MESSAGE3) << "  Read-out pitch in Z     is 2x geom. pitch." << std::endl << std::endl;
}

//
// Method printing hit info
//
void SiStripClus::printHitInfo(const StripCluster* pCluster) const {
  short int layerID = pCluster->getLayerID();
  short int ladderID = pCluster->getLadderID();
  short int sensorID = pCluster->getSensorID();

  Hep3Vector position = pCluster->get3Position();

  // Print sensor info
  streamlog_out(MESSAGE2) << "    Layer " << _geometry->getLayerRealID(layerID) << " "
                          << "Ladder " << ladderID << " "
                          << "Sensor " << sensorID << std::endl;
  streamlog_out(MESSAGE2) << "     Hit in local ref. system: " << std::fixed << std::setprecision(3)
                          << std::setiosflags(std::ios::showpos) << "PosX [mm]: " << position.getX() / mm << " "
                          << "PosY [mm]: " << position.getY() / mm << " "
                          << "PosZ [mm]: " << position.getZ() / mm << std::setprecision(0)
                          << std::resetiosflags(std::ios::internal | std::ios::showpos) << std::endl;

  position = _geometry->transformPointToGlobal(layerID, ladderID, sensorID, position);

  streamlog_out(MESSAGE2) << "     Hit in global ref. system:  " << std::fixed << std::setprecision(3)
                          << std::setiosflags(std::ios::showpos) << "PosX [mm]: " << position.getX() / mm << " "
                          << "PosY [mm]: " << position.getY() / mm << " "
                          << "PosZ [mm]: " << position.getZ() / mm << std::setprecision(0)
                          << std::resetiosflags(std::ios::internal | std::ios::showpos) << std::endl;
}

} // namespace sistrip
