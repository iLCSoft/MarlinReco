// History:
//- added possibility to digitize forward slanted sensors with non-uniform pitch in RPhi (along Z-axiz), Z. Drasal Jul
// 2009
//- implemented new parameter for cut on time to emulate the integration time, K. Prothmann Dec. 2009
//- added possibility to digitize floating strips (in RPhi or Z) too, Z. Drasal Feb 2010
//- added detailed landau fluctuations in Si based on Geant 4 G4UniversalFluctuation, Z. Drasal Apr 2010
//- created only 1 relation TrackerPuls <--> SimTrackerHit, removed relations TrackerPulse <--> MCParticle, MCParticle
//<--> SimTrackerHit, Z. Drasal May 2010
//- corrected charge units, Z. Drasal May 2010

// FIXME: ?? No guarda los clusters con carga negativa??

#include "SiStripDigi.h"

#include "SiStripGeomBuilder.h"

#include <cstdlib>
#include <iomanip>
#include <math.h>
#include <time.h>

// Include CLHEP header files
#include <CLHEP/Random/Randomize.h>
#include <CLHEP/Vector/ThreeVector.h>

// Include Digi header files
#include "Colours.h"
#include "DigiCluster.h"
#include "PhysicalConstants.h"
#include "RombIntSolver.h"

// Include LCIO header files
#include "UTIL/LCTrackerConf.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <lcio.h>

// Include Marlin
#include <streamlog/streamlog.h>

// Namespaces
using namespace CLHEP;
using namespace lcio;
using namespace marlin;

namespace sistrip {

//
// Instantiate this object
//
SiStripDigi anSiStripDigi;

// enum for variable map
enum { ELOSSDIGI, ELOSSG4 };

// PUBLIC METHODS

//
// Constructor
//
SiStripDigi::SiStripDigi() : Processor("SiStripDigi") {
  // Processor description
  _description = "SiStripDigi: Marlin processor intended for detailed digitization of silicon strip sensors";

  //
  // Processor parameters

  // Define compulsory parameters
  // NEW -----
  registerProcessorParameter("Subdetector", "Subdetector to be digitised", _subdetector, std::string("FTD"));

  registerProcessorParameter("DeplVoltage", "Depletion voltage of a silicon strip detector, set in volts", _Vdepl,
                             float(60));

  registerProcessorParameter("BiasVoltage", "Bias voltage set on a silicon strip detector, set in volts", _Vbias,
                             float(150));

  registerProcessorParameter("Temperature", "Temperature measured on a sensor, set in Kelvins", _temp, float(300));

  registerProcessorParameter("LandauFluct", "Use internal Landau fluctuations (instead of Geant4)?", _landauFluct,
                             bool(true));

  registerProcessorParameter("LandauBetaGammaCut", "Below this beta*gamma factor internal Landau fluctuations not used",
                             _landauBetaGammaCut, float(0.7));

  registerProcessorParameter("ProductionThreshOnDeltaRays",
                             "Production threshold cut on delta electrons in keV  (for Landau fluct.) - use the same "
                             "as in Geant4 (80keV ~ 0.05 mm)",
                             _prodThreshOnDeltaRays, double(80));

  registerProcessorParameter("ElectronicEffects", "Apply electronic effects?", _electronicEffects, bool(true));

  registerProcessorParameter("InterStripCapacitance", "Interstrip capacitance, set in pF", _capInterStrip, float(6.));

  registerProcessorParameter("BackplaneCapacitance", "Backplane capacitance, set in pF", _capBackPlane, float(0.));

  registerProcessorParameter("CouplingCapacitance", "Coupling capacitance, set in pF", _capCoupl, float(120.));

  registerProcessorParameter("ElectronicsNoise", "Noise added by the electronics, set in electrons", _elNoise,
                             float(1000));

  registerProcessorParameter("FloatingStripsRPhi", "Is every even strip floating in R-Phi?", _floatStripsRPhi,
                             bool(false));

  registerProcessorParameter("FloatingStripsZ", "Is every even strip floating in Z?", _floatStripsZ, bool(false));

  registerProcessorParameter("AbsoluteSpacePrecision", "Absolute digitization space precision, set in microns",
                             _epsSpace, float(20));

  registerProcessorParameter("RelativeAnglePrecision",
                             "Relative digitization precision when calculating tan(Lorentz angle)", _epsAngle,
                             float(0.01));

  registerProcessorParameter("RelativeDriftTimePrecision",
                             "Relative digitization precision when calculating cluster drift time", _epsTime,
                             float(0.01));

  registerProcessorParameter("InputCollectionName", "Name of SimTrackerHit input collection", _inColName,
                             std::string("FTDCollection"));

  registerProcessorParameter("OutputCollectionName", "Name of TrackerPulse output collection", _outColName,
                             std::string("FTDDigits"));

  registerProcessorParameter("RelCollectionNameTrkPlsToSimHit",
                             "Name of relation collection - TrkPulses to SimTrackerHit (if nonzero, created)",
                             _relColNamePlsToSim, std::string("FTDDigitsToSimHitsRel"));

  registerProcessorParameter("IntegrationWindow", "Use integration window?", _integrationWindow, bool(false));

  registerProcessorParameter("StartIntegration",
                             "Only Simulated hits after the StartIntegration time in ns will be digitized",
                             _startIntegration, double(-10.0));

  registerProcessorParameter("StopIntegration",
                             "Only Simulated hits before the StopIntegration time in ns will be digitized",
                             _stopIntegration, double(10.0));

  registerProcessorParameter("pitchFront", "Pitch of the front sensor family in the center of the sensors, set in um",
                             _pitchFront, double(50.0));

  registerProcessorParameter("pitchRear", "Pitch of the rear sensor family in the center of the sensors, set in um",
                             _pitchRear, double(50.0));
}

//
// Method called at the beginning of data processing
//
void SiStripDigi::init() {
  // Initialize variables
  _nRun = 0;
  _nEvent = 0;

  _currentLayerID = 0;
  _currentLadderID = 0;
  _currentSensorID = 0;

  _sensorThick = -1.;

  // Set variables in appropriate physical units
  _Vdepl *= V;
  _Vbias *= V;

  _temp *= K;
  _elNoise *= e;
  _epsSpace *= um;

  _prodThreshOnDeltaRays *= keV;

  _startIntegration *= ns;
  _stopIntegration *= ns;

  _pitchFront *= um;
  _pitchRear *= um;

  // Get geometry parameters from Gear xml file
  _geometry = SiStripGeomBuilder::Build(_subdetector, _pitchFront, _pitchRear);
  _geometry->initGearParams();
  //	_geometry -> printGearParams();  FIXME

  // Initialize random generator (engine, mean, sigma)
  _genGauss = new RandGauss(new RandEngine(SEED), 0., (double)_elNoise);

  // Print set parameters
  printProcessorParams();

  // CPU time start
  _timeCPU = clock() * us;

  // Internal Landau fluctuations
  if (_landauFluct) {
    _fluctuate = new SiEnergyFluct(_prodThreshOnDeltaRays);
  }

  //
  // ROOT variables
  //
#ifdef ROOT_OUTPUT_LAND
  _rootFile = new TFile("Digi_FTD.root", "RECREATE");
  _rootFile->cd("");

  _tree = new TTree("Digi", "Digi info");

  std::map<int, std::string> names;

  // Variables to store
  _variables[ELOSSDIGI] = 0.0;
  names[ELOSSDIGI] = "ELossDigi";

  _variables[ELOSSG4] = 0.0;
  names[ELOSSG4] = "ELossG4";

  for (std::map<int, double>::iterator it = _variables.begin(); it != _variables.end(); ++it) {
    _tree->Branch(names[it->first].c_str(), &(it->second), (names[it->first] + "/D").c_str());
  }
#endif
}

//
// Method called for each run
//
void SiStripDigi::processRunHeader(LCRunHeader* run) {
  // Print run number
  streamlog_out(MESSAGE3) << DGREEN << " Processing run: " << ENDCOLOR << (run->getRunNumber()) << std::endl
                          << std::endl;

  // Add info about this process to LCIO run header
  run->parameters().setValue("[ProcessorName]", name());
  run->parameters().setValue("[" + name() + "_gearType]", _geometry->getGearType());
  //   run->parameters().setValue ( "["+name()+"_inColName]" , _inColName              );
  //   run->parameters().setValue ( "["+name()+"_outColName]", _outColName             );

  _nRun++;
}

//
// Method called for each event
//
void SiStripDigi::processEvent(LCEvent* event) {
  // Local error count variable
  static int errors_occured = 0;

  // Get names of all collections saved in LCIO file
  ConstStringVec* strVec = event->getCollectionNames();

  // Initialize number of found muon hits versus all hits - efficiency map, resp.
  // dep. energy in each event
#ifdef ROOT_OUTPUT_LAND
  for (std::map<int, double>::iterator it = _variables.begin(); it != _variables.end(); ++it) {
    it->second = 0.0;
  }
#endif

  //
  // Get SimTrackerHit collection
  LCCollection* colOfSimHits = 0;

  // Print header info about collections (for 1. event)  // FIXME: To be deleted
  if (event->getEventNumber() == 0) {
    streamlog_out(MESSAGE3) << " " << DUNDERL << DBLUE << "LCCollection processed:" << ENDCOLOR << std::endl
                            << std::endl;
  }

  // Go through all collection names and save the one that is defined as an input
  // collection
  for (ConstStringVec::const_iterator colName = strVec->begin(); colName != strVec->end(); ++colName) {
    if (_inColName == (*colName)) {
      // Collection must be of SimTrackerHit type
      if ((event->getCollection(*colName))->getTypeName() == LCIO::SIMTRACKERHIT) {
        // Save collection
        colOfSimHits = event->getCollection(*colName);
        // Print info
        if (event->getEventNumber() == 0) {
          streamlog_out(MESSAGE3) << "  " << *colName << std::endl;
        }
      }
      // Collection NOT of SimTrackerHit type
      else {
        streamlog_out(ERROR) << "Required collection: " << _inColName << " found, but NOT of SIMTRACKERHIT type!!!"
                             << std::endl;
        exit(-1);
      }
    }
  } // For - collection names

  // Input collection NOT found
  if (colOfSimHits == 0) {
    streamlog_out(WARNING) << "Required collection: " << _inColName << " not found!!!" << std::endl;
  }

  // Print info
  if (event->getEventNumber() == 0) {
    streamlog_out(MESSAGE3) << std::endl;
  }
  //
  // Print event number
  //   if ( ((event->getEventNumber()+1)%100 == 0) || ((event->getEventNumber()+1) == 1))
  if ((event->getEventNumber() + 1) % 100 == 0) {
    streamlog_out(MESSAGE3) << DYELLOW << "  Processing event: " << ENDCOLOR << std::setw(5)
                            << (event->getEventNumber() + 1) << std::endl
                            << std::endl;
  }

  //
  // Process SimTrackerHit collection
  if (colOfSimHits != 0) {
    // Sensor map of strips with total integrated charge
    SensorStripMap sensorMap;

    // Set collection decoder
    CellIDDecoder<SimTrackerHit> cellIDDec(colOfSimHits);
    // Guarda la string encoder...??

    // Get number of SimTrackerHits in the collection
    int nHits = colOfSimHits->getNumberOfElements();

    // Process hits
    for (int i = 0; i < nHits; ++i) {
      //
      // Copy the content of current simhit to a DigiHit
      SimTrackerHit* simHit = dynamic_cast<SimTrackerHitImpl*>(colOfSimHits->getElementAt(i));

#ifndef ROOT_OUTPUT_LAND
      // Cut on simHit creation time --> simulate integration time of a sensor (if option switched on))
      if ((simHit != 0) && (_integrationWindow)) {
        if (simHit->getTime() * ns < _startIntegration || simHit->getTime() * ns > _stopIntegration) {
          errors_occured++;
          continue;
        }
      }
#endif
      // Hit phys. info
      double hitPos[3] = {*(simHit->getPosition()) * mm, *(simHit->getPosition() + 1) * mm,
                          *(simHit->getPosition() + 2) * mm};
      float hitMom[3] = {*(simHit->getMomentum()) * GeV, *(simHit->getMomentum() + 1) * GeV,
                         *(simHit->getMomentum() + 2) * GeV};
      float hitdEdx = simHit->getEDep() * GeV;
      float hitTime = simHit->getTime() * ns;
      float hitPath = simHit->getPathLength() * mm;

      // Hit geom. info
      // Recall: codification in Mokka/G4:
      //         Layer: 1,...7
      //         Side:  -1  +1
      //         Petal: 1...16
      //         Sensor: 1 2 3 4
      //---- The notation used in this code:
      //     Layer:  0,...., 6 (Positive Z)
      //             7,....,13 (Negative Z)
      //     Petal:  0,....,15
      //     Sensor: 1 2 3 4
      std::map<std::string, int> cellIDmap = _geometry->decodeCellID(cellIDDec(simHit)); // decodeCellID0
      const int hitLayerID = cellIDmap["layer"];
      const int hitLadderID = cellIDmap["module"];
      const int hitSensorID = cellIDmap["sensor"];
      // FIXME: The pixel disks digitization must be done
      //        by some other package...
      if (abs(_geometry->getLayerRealID(hitLayerID)) < 3) {
        continue;
      }
      int hitCellID = cellIDDec(simHit).lowWord();

      // MC information
      MCParticle* mcPart = simHit->getMCParticle();

      //
      // Set information of currently digitized hit
      SimTrackerDigiHit* simDigiHit = new SimTrackerDigiHit();

      simDigiHit->setPrePosition(hitPos, hitMom, hitPath);
      simDigiHit->setPosPosition(hitPos, hitMom, hitPath);
      simDigiHit->setMomentum(hitMom);
      simDigiHit->setEDep(hitdEdx);
      simDigiHit->setTime(hitTime);
      simDigiHit->setPathLength(hitPath);
      simDigiHit->setCellID0(hitCellID);
      simDigiHit->setLayerID(hitLayerID);
      simDigiHit->setLadderID(hitLadderID);
      simDigiHit->setSensorID(hitSensorID);
      simDigiHit->setMCParticle(mcPart);
      simDigiHit->setSimTrackerHit(dynamic_cast<EVENT::SimTrackerHit*>(simHit));

      // Set current layer --> ladder --> sensor ID
      _currentLayerID = hitLayerID;
      _currentLadderID = hitLadderID;
      _currentSensorID = hitSensorID;

      // Set actual sensor geometry
      _sensorThick = _geometry->getSensorThick(_currentLayerID);
      // Set magnetic field
      _magField = _geometry->get3MagField();

      // Transform hit to local ref. system of a sensor
      transformSimHit(simDigiHit);

      // Transform magnetic field to local ref. system of a sensor
      _magField = _geometry->transformVecToLocal(_currentLayerID, _currentLadderID, _currentSensorID, _magField);

      // Digitize the given hit and get sensor map of strips with total
      // integrated charge and time when a particle crossed the sensor
      digitize(simDigiHit, sensorMap);

      // Unset actual sensor parameters
      _currentLayerID = 0;
      _currentLadderID = 0;
      _currentSensorID = 0;

      _magField.set(0., 0., 0.);

      // Clear memory
      delete simDigiHit;
      simDigiHit = 0;

    } // For - process hits

    // Add electronics effects
    if (_electronicEffects) {
      // Calculate crosstalk and add this effect
      calcCrossTalk(sensorMap);
      // Generate noise and add this effect
      genNoise(sensorMap);
    }
    // Print final info
    printStripsInfo("all effects included", sensorMap);

    // Strips Types to loop overt there
    std::vector<StripType> stripTypesvect;
    stripTypesvect.push_back(STRIPFRONT);
    stripTypesvect.push_back(STRIPREAR);

    /*#ifdef ROOT_OUTPUT_LAND
                    for(SensorStripMap::iterator itSM = sensorMap.begin(); itSM != sensorMap.end();
                                    ++itSM)
                    {
                            const int cellID = itSM->first;
                            std::map<std::string,int> bfM = _geometry->decodeCellID(cellID);
                            const int diskID=bfM["layer"];
                            const int petalID=bfM["module"];
                            const int sensorID=bfM["sensor"];
                            std::cout << "-|| cellID: " << cellID
                                    << " (layer:" << diskID << ", petal:"<< petalID
                                    << ", sensor:"<< sensorID << ")" << std::endl;
                            for(std::vector<StripType>::iterator stT = stripTypesvect.begin();
                                            stT != stripTypesvect.end(); ++stT)
                            {
                                    for(StripChargeMap::iterator itSC = itSM->second[*stT].begin();
                                                    itSC != itSM->second[*stT].end(); ++itSC)
                                    {
                                            const double ylocal
    =_geometry->getStripPosY(diskID,sensorID,itSC->first,0.0); const CLHEP::Hep3Vector PosGlobal0 =
    _geometry->transformPointToGlobal(diskID,petalID,sensorID,CLHEP::Hep3Vector(0.00,ylocal,0.0)); const
    CLHEP::Hep3Vector unitV = _geometry->getStripUnitVector(diskID,sensorID,itSC->first); const double alpha =
    fabs(acos(unitV.getZ())); const double ylocalL = ylocal+tan(alpha)*_geometry->getSensorLength(diskID);

                                            const CLHEP::Hep3Vector PosGlobal1 =
    _geometry->transformPointToGlobal(diskID,petalID,sensorID,CLHEP::Hep3Vector(0.00,ylocalL,_geometry->getSensorLength(diskID)));
                                            std::cout << " |-|- stripID: " <<  itSC->first
                                                    << " -- Position (at zlocal=0):" << PosGlobal0
                                                    << "  (at z=L:" << PosGlobal1 << ")" << std::endl;
                                            std::cout << "   |--- Total integrated charge: " <<
    itSC->second->getCharge()/fC
                                                    << " fC (" << itSC->second->getCharge() << ")"  << std::endl;
                                            std::cout << "   |--- Time                   : " <<  itSC->second->getTime()
                                                    << std::endl;
                                            std::cout << "   |-|- Simulated hits " <<  std::endl;
                                            SimTrackerHitMap hm = itSC->second->getSimHitMap();
                                            for(SimTrackerHitMap::iterator it = hm.begin(); it != hm.end();
                                                            ++it)
                                            {
                                                    std::cout << "     |-- Position Global[mm]  : "
                                                            << "(" << (it->first->getPosition()[0]*mm)/mm
                                                            << "," << (it->first->getPosition()[1]*mm)/mm
                                                            << "," << (it->first->getPosition()[2]*mm)/mm << ")" <<
    std::endl; std::cout << "     |-- Energy deposited[keV]: "
                                                            << it->first->getEDep()/keV << "  (w=" << it->second
                                                            << ")" << std::endl;
                                                    std::cout << "     |-- Path Length[um]      : "
                                                            << it->first->getPathLength()/um << std::endl;
                                                    std::cout << "     |-- CellID0              : "
                                                            << it->first->getCellID0() << std::endl;
                                                    std::cout << "     |---------------------------" << std::endl;
                                            }
                                    }
                            }
                            std::cout << std::endl;
                    }
    #endif*/

    //
    // Create new collection (TrackerPulses + relations) from obtained results
    IMPL::LCCollectionVec* colOfTrkPulses = new IMPL::LCCollectionVec(LCIO::TRACKERPULSE);

    IMPL::LCCollectionVec* colOfRelPlsToSim = NULL;

    if (!_relColNamePlsToSim.empty()) {
      colOfRelPlsToSim = new IMPL::LCCollectionVec(LCIO::LCRELATION);
    }

    // Set collection flag - cellID 1 will be stored, relations will contain weights
    LCFlagImpl flag1(0), flag2(0);

    flag1.setBit(LCIO::TRAWBIT_ID1);
    flag2.setBit(LCIO::LCREL_WEIGHTED);

    colOfTrkPulses->setFlag(flag1.getFlag());
    if (!_relColNamePlsToSim.empty()) {
      colOfRelPlsToSim->setFlag(flag2.getFlag());
    }

    // CODIFICATION --->FIXME: METHOD IN GEAR (Centralizing...)
    CellIDEncoder<TrackerPulseImpl> cellEnc(LCTrackerCellID::encoding_string() + ",stripType:2,stripID:11",
                                            colOfTrkPulses);

    // TrackerPulses
    SensorStripMap::iterator iterSMap;
    StripChargeMap::iterator iterChMap;

    for (iterSMap = sensorMap.begin(); iterSMap != sensorMap.end(); ++iterSMap) {
      int cellID = iterSMap->first;

      for (std::vector<StripType>::iterator stripTypeIt = stripTypesvect.begin(); stripTypeIt != stripTypesvect.end();
           ++stripTypeIt) {
        StripType stripType = *stripTypeIt;
        for (iterChMap = iterSMap->second[stripType].begin(); iterChMap != iterSMap->second[stripType].end();
             ++iterChMap) {
          int stripID = iterChMap->first;

          // Create new Tracker pulse (time in ns, charge in fC), if charge
          // positive, otherwise will be cut-out in clustering anyway
          if (iterChMap->second->getCharge() > 0. * fC) {
            TrackerPulseImpl* trkPulse = new TrackerPulseImpl();

            // Storing the cellID0 and cellID1
            _geometry->updateCanonicalCellID(cellID, stripType, stripID, &cellEnc);
            cellEnc.setCellID(trkPulse);

            trkPulse->setTime(iterChMap->second->getTime() / ns);
            trkPulse->setCharge(iterChMap->second->getCharge() / fC);
            trkPulse->setTrackerData(0);

            // Create relation from Tracker pulse to SimTrackerHit
            if (!_relColNamePlsToSim.empty()) {
              const SimTrackerHitMap& simHitMap = iterChMap->second->getSimHitMap();
              float weightSum = iterChMap->second->getSimHitWeightSum();

              for (SimTrackerHitMap::const_iterator iterSHM = simHitMap.begin(); iterSHM != simHitMap.end();
                   ++iterSHM) {
                // Create LC relation
                LCRelationImpl* relation = new LCRelationImpl;
                SimTrackerHit* simHit = iterSHM->first;
                float weight = 0.;

                // Set from TrkPulse to MCParticle
                relation->setFrom(trkPulse);
                relation->setTo(simHit);
                // Set weight
                if (weightSum != 0.) {
                  weight = float(iterSHM->second) / weightSum;
                } else {
                  weight = 0.;
                }
                relation->setWeight(weight);

                // Add
                if (weight != 0.) {
                  colOfRelPlsToSim->addElement(relation);
                } else {
                  delete relation;
                }
              }
            }
            // Save the pulse to the collection
            colOfTrkPulses->addElement(trkPulse);
          }
          // FIXME: que pinta esto aqui
          /*			     short int layerID  = 0;
                                       short int ladderID = 0;
                                       short int sensorID = 0;
                                       _geometry->decodeCellID(layerID, ladderID, sensorID, cellID);*/

          // Release memory
          delete iterChMap->second;
        }

      } // For for strips type
      // Release memory
      delete[] iterSMap->second;
    } // For */
    //
    // Save the collection (vector) of pulses + relations to MC
    event->addCollection(colOfTrkPulses, _outColName);

    if (!_relColNamePlsToSim.empty()) {
      event->addCollection(colOfRelPlsToSim, _relColNamePlsToSim);
    }

  } // If (colOfSimHits != 0)
  if (errors_occured != 0) {
    streamlog_out(DEBUG4) << "SimTrackerHit outside " << errors_occured
                          << " times outside of integration time of sensor" << std::endl;
  }
  errors_occured = 0;

#ifdef ROOT_OUTPUT_LAND
  _tree->Fill();
#endif
  _nEvent++;
}

//
// Method called after each event to check the data processed
//
void SiStripDigi::check(LCEvent* event) {}

//
// Method called after all data processing
//
void SiStripDigi::end() {
  // Release memory
  delete _geometry;
  _geometry = 0;

  // CPU time end
  _timeCPU = clock() * us - _timeCPU;

  // Print message
  streamlog_out(MESSAGE3) << std::endl
                          << " "
                          << "Time per event: " << std::setiosflags(std::ios::fixed | std::ios::internal)
                          << std::setprecision(3) << _timeCPU / _nEvent / ms << " ms" << std::endl
                          << std::setprecision(3) << std::endl
                          << DGREEN << " "
                          << "Processor succesfully finished!" << ENDCOLOR << std::endl;

  // Clean memory
  if (_landauFluct)
    delete _fluctuate;

#ifdef ROOT_OUTPUT_LAND

  // Close file
  _rootFile->cd("");
  _rootFile->Write();
  _rootFile->Close();

  delete _rootFile;

#endif
}

// MAIN DIGI METHOD

//
// The main method digitizing given hit (input parameter: a pointer to
// SimTrackerDigiHit, output parameter: sensor map of strips with total
// integrated charge and time when particle crossed the sensor)
//
void SiStripDigi::digitize(const SimTrackerDigiHit* simDigiHit, SensorStripMap& sensorMap) {
  //
  // Calculate electron, resp. hole, clusters from obtained hits & provide Landau fluctuations
  // DigiClusterVec eClusterVec;
  DigiClusterVec hClusterVec;

  // FIXME: Return the vector, not passing per reference
  calcClusters(simDigiHit, hClusterVec);

  //
  // Go through all e, resp. h, clusters, perform drift, diffusion and Lorentz shift
  DigiClusterVec::iterator iterVec;

  // Cluster - total diffusion - sigma
  double diffSigma = 0.;

  // Cluster - Lorentz shift
  Hep3Vector shiftLorentz(0., 0., 0.);

  // Cluster - total drift time
  double driftTime = 0.;

  // Cluster - charge collected by a strip
  double chargeCollect = 0.;

  // Get the type of strip (sensors 1,2 or sensors 3,4)
  int stripType = -1;
  if (_currentSensorID == 1 || _currentSensorID == 2) {
    stripType = STRIPFRONT;
  } else if (_currentSensorID == 3 || _currentSensorID == 4) {
    stripType = STRIPREAR;
  } else {
    streamlog_out(ERROR) << "SiStripDigi::digitize - "
                         << " Incoherent sensor ID! ID = " << _currentSensorID << " do not exist. " << std::endl;
    exit(-1);
  }

  // Hole clusters
  for (iterVec = hClusterVec.begin(); iterVec != hClusterVec.end(); ++iterVec) {
    DigiCluster* cluster = (*iterVec);

    // Calculate charge collected by a strip (using Shockley-Ramo theorem
    // dQh = |qClusterH * x / d|)
    // chargeCollect = cluster->getNCarriers() *
    //        (cluster->getPosX())/_sensorThick * fC; // Presume that Z strips
    // collect only holes
    chargeCollect = cluster->getNCarriers() * e;

    // Calculate drift time
    driftTime = getHoleDriftTime(cluster->getPosX());
    cluster->setDriftTime(driftTime);

    // Calculate mean diffusivity and then diffusion sigma
    diffSigma = getHoleDiffusivity(cluster->getPosX() / 2.);
    diffSigma *= 2 * driftTime;
    diffSigma = sqrt(diffSigma);

    // Calculate Lorentz shift + evaluate cluster new position, i.e. X = 0.
    shiftLorentz = getHoleLorentzShift(cluster->getPosX());
    cluster->set3Position(shiftLorentz + Hep3Vector(0., cluster->getPosY(), cluster->getPosZ()));

    //
    // Add diffusion effect and save results as signals (IN Z)
    //
    int iMinStrip = _geometry->getStripID(cluster->getLayerID(), cluster->getSensorID(),
                                          (cluster->getPosY() - 3. * diffSigma), cluster->getPosZ());
    int iMaxStrip = _geometry->getStripID(cluster->getLayerID(), cluster->getSensorID(),
                                          (cluster->getPosY() + 3. * diffSigma), cluster->getPosZ());
    double sensorPitch = _geometry->getSensorPitch(_currentLayerID, _currentSensorID, cluster->getPosZ());

    // The calculus for the strips must be done in the (rotated by the
    // stereo angle) local frame
    const double tanPhi = tan(_geometry->getLayerHalfPhi(_currentLayerID));
    // const double stAngle = _geometry->getStereoAngle(_currentLayerID,_currentSensorID);
    const CLHEP::Hep3Vector pointRot = _geometry->transformPointToRotatedLocal(
        _currentLayerID, _currentSensorID, CLHEP::Hep3Vector(0., cluster->getPosY(), cluster->getPosZ()));
    const double zRot = pointRot.getZ();
    const double yorigenRot = (_geometry->getSensorLength(_currentLayerID) - zRot) * tanPhi;

    //  Gauss distr. - primitive function: from A to B
    double meanRot = pointRot.getY();

    //  Gauss distr. - primitive function: from A to B
    // double mean       = cluster->getPosY(); //cluster->getPosZ();

    double sigmaSqrt2 = diffSigma * sqrt(2.);
    double primAtA = 0.5 * (1. + erf(((yorigenRot + iMinStrip * sensorPitch) - meanRot) / sigmaSqrt2));
    double primAtB = 0.;

    //  Sensor map of strips with total integrated charge and time when particle
    //  crossed the detector
    SensorStripMap::iterator iterSMap;
    StripChargeMap::iterator iterChMap;

    //  Strip signal
    Signal* signal = 0;

    //  Calculate signal at each strip and save
    for (int i = iMinStrip; i <= iMaxStrip; ++i) {
      // Charge collected by a strip i
      double charge = 0.;

      // Gauss distr. - prim. function at B
      primAtB = 0.5 * (1. + erf(((yorigenRot + (i + 1) * sensorPitch) - meanRot) / sigmaSqrt2));

      // Integration result
      charge = (primAtB - primAtA) * chargeCollect;

      // New integration starting point
      primAtA = primAtB;

      // Find sensor with given ID
      iterSMap = sensorMap.find(cluster->getCellID());

      // Sensor has already collected some charge
      if (iterSMap != sensorMap.end()) {
        // Find strip i in Z
        iterChMap = iterSMap->second[stripType].find(i);

        // Strip i in Z has already collected some charge
        if (iterChMap != iterSMap->second[stripType].end()) {
          // Get signal at strip i
          signal = iterChMap->second;

          // Update signal information (add current charge to
          // the total)
          signal->updateCharge(charge);

          // Update MC truth information
          signal->updateSimHitMap(cluster->getSimTrackerHit(), charge);

        }
        // Strip i in Z hasn't collected any charge yet
        else {
          // Create signal, i.e. total charge + time when
          // particle crossed the detector, all corresponding
          // to strip i
          signal = new Signal(charge, cluster->getTime());

          // Save MC truth information
          signal->updateSimHitMap(cluster->getSimTrackerHit(), charge);

          // Save information about strip i in Z
          iterSMap->second[stripType][i] = signal;
        }
      }
      // Sensor has not collected any charge yet
      else {
        // Create map of strips with total collected charge
        StripChargeMap* stripMap = new StripChargeMap[2];

        // Create signal, i.e. total charge + time when particle crossed the
        // detector, all corresponding to strip i
        signal = new Signal(charge, cluster->getTime());

        // Save MC truth information
        signal->updateSimHitMap(cluster->getSimTrackerHit(), charge);

        // Save information about strip i in Z
        stripMap[stripType][i] = signal;
        sensorMap[cluster->getCellID()] = stripMap;
      }
    } // Calculate signal at each strip
    // Release memory
    delete cluster;
    cluster = 0;
  } // For holes

  hClusterVec.clear();
}

// OTHER METHODS

//
// Method that calculates e-h clusters from simulated hits (parameters: a pointer
// to SimTrackerDigiHit and two vectors of pointers to electron, resp. hole
// clusters (pairs) )
//
void SiStripDigi::calcClusters(const SimTrackerDigiHit* hit, DigiClusterVec& hClusterVec)
// DigiClusterVec & eClusterVec, DigiClusterVec & hClusterVec)
{
  // DigiCluster * eCluster  = 0;
  DigiCluster* hCluster = 0;

  // Calculate the number of clusters (each cluster sits in the middle of a substep,
  // where the substep increases precision of a step and corresponds to space precision)
  int nClusters = 0;
  int nSubSteps = (hit->getStepSize()) / _epsSpace;

  if ((((hit->getStepSize()) - nSubSteps * _epsSpace) >= (0.5 * _epsSpace)) || ((hit->getStepSize()) < _epsSpace)) {
    nClusters = nSubSteps + 1;
  } else {
    nClusters = nSubSteps;
  }

  // Calculate cluster step
  Hep3Vector clusterStep = hit->get3Step() / (nClusters);

  // Calculate variables relevant for Landau fluctuations - info about MC particle
  MCParticle* mcPart = dynamic_cast<MCParticle*>(hit->getMCParticle());
  double partMomMag = double(hit->get3Momentum().mag());
  double partMass = 0.;
  double betaGamma = 0.;
  if (mcPart != 0) {
    partMass = hit->getMCParticle()->getMass() * GeV;
  }
  if (mcPart != 0 && partMass != 0.) {
    betaGamma = partMomMag / partMass;
  }

  double clusterStepSize = clusterStep.mag();

#ifdef ROOT_OUTPUT_LAND

  // Update info
  _variables[ELOSSG4] += hit->getEDep() / keV;
  //_rootDepEG4 += hit->getEDep();
#endif
  // Decide if particle is low-energy --> Geant 4 Landau fluctuations used in
  // low energy regime
  bool lowEnergyPart = false;

  if (betaGamma < _landauBetaGammaCut) {
    lowEnergyPart = true;
  }

  // Create elec-clusters and hole-clusters
  for (int i = 0; i < nClusters; ++i) {
    // Landau fluctuated deposited energy
    double eLoss = 0.;

    // High energy regime --> use internal Landau fluctuation (if defined)
    if (_landauFluct && !lowEnergyPart && mcPart != 0) {
      eLoss = _fluctuate->SampleFluctuations(mcPart, clusterStepSize);
    }
    // Low energy regime --> use Geant4 info and distribute it uniformly
    else {
      eLoss = hit->getEDep() / nClusters;
    }

#ifdef ROOT_OUTPUT_LAND
    // Update info
    _variables[ELOSSDIGI] += eLoss / keV;
    //_rootDepEDigi += eLoss;
#endif

    // Electron clusters
    /*eCluster = new DigiCluster(-1, hit->getTime(), eLoss,
                    hit->get3PrePosition() + (i+0.5)*clusterStep,
                    hit->getLayerID()   , hit->getLadderID(),
                    hit->getSensorID()  , hit->getCellID0(),
                    hit->getMCParticle(), hit->getSimTrackerHit());
    eClusterVec.push_back(eCluster);*/

    // Hole clusters
    hCluster = new DigiCluster(+1, hit->getTime(), eLoss, hit->get3PrePosition() + (i + 0.5) * clusterStep,
                               hit->getLayerID(), hit->getLadderID(), hit->getSensorID(), hit->getCellID0(),
                               hit->getMCParticle(), hit->getSimTrackerHit());
    hClusterVec.push_back(hCluster);
  }
  // Print detailed info about clusters
  // printClustersInfo("Electron clusters in local ref. system", eClusterVec);
  printClustersInfo("Hole clusters in local ref. system", hClusterVec);
  streamlog_out(MESSAGE2) << std::endl;
}

//
// Method that calculates crosstalk effect, i.e. total charge redistribution
// (input parameter: sensor map of strips with total integrated charge)
//
void SiStripDigi::calcCrossTalk(SensorStripMap& sensorMap) {
  // Calculate how much charge is collected by adjacent strips
  // If read-out pitch = geom. pitch
  static float capRatio = _capInterStrip / (_capInterStrip + _capBackPlane + _capCoupl);
  // If read-out pitch = 2x geom. pitch
  static float capFloatRatio = _capInterStrip / 2. / (_capInterStrip / 2. + _capBackPlane + _capCoupl);

  // Iterators
  SensorStripMap::iterator iterSMap;
  StripChargeMap::iterator iterChMap;

  short int layerID = 0;
  short int ladderID = 0;
  short int sensorID = 0;

  Signal* signal = 0;

  // Strips types
  std::vector<StripType> stvec;
  stvec.push_back(STRIPFRONT);
  stvec.push_back(STRIPREAR);

  // Go through all sensors
  for (iterSMap = sensorMap.begin(); iterSMap != sensorMap.end(); ++iterSMap) {
    // Find corresponding layerID
    std::map<std::string, int> bfMap = _geometry->decodeCellID(iterSMap->first);
    layerID = bfMap["layer"];
    ladderID = bfMap["module"];
    sensorID = bfMap["sensor"];
    //_geometry->decodeCellID(layerID, ladderID, sensorID, iterSMap->first);

    // Define strip map with recalculated signals
    StripChargeMap* stripMap = new StripChargeMap[2];

    // For all the strips types
    for (std::vector<StripType>::iterator itT = stvec.begin(); itT != stvec.end(); ++itT) {
      const StripType STRIP = *itT;
      // Read map of strips in R-Phi and save new results to recalculated strip map
      for (iterChMap = (iterSMap->second[STRIP]).begin(); iterChMap != (iterSMap->second[STRIP]).end(); ++iterChMap) {
        const int iStrip = iterChMap->first;
        int iStripLeft = iStrip - 1;
        int iStripRight = iStrip + 1;
        double Kf = 0.;

        // Calculate charge redistribution
        signal = iterChMap->second;

        const double chargeTotal = signal->getCharge();

        SimTrackerHitMap simHitMapLeft;
        SimTrackerHitMap simHitMapCentr;
        SimTrackerHitMap simHitMapRight;

        // Floating strip - capRatio = 0.5
        if ((_floatStripsRPhi) && (iStrip % 2 == 1)) {
          Kf = 0.5;
        }
        // Read-out strip, adjacent floating - capRatio = capFloatRatio
        else if ((_floatStripsRPhi) && (iStrip % 2 == 0)) {
          Kf = capFloatRatio;

          iStripLeft--;
          iStripRight--;
        }
        // All strips are read-out - capRatio = catRatio
        else {
          Kf = capRatio;
        }

        // Calculating charge redistribution
        const double chargeLeft = chargeTotal * Kf;
        const double chargeRight = chargeTotal * Kf;
        const double chargeCentr = chargeTotal - chargeLeft - chargeRight;

        // Reasigning weights to hits
        for (SimTrackerHitMap::const_iterator iterSHM = signal->getSimHitMap().begin();
             iterSHM != signal->getSimHitMap().end(); ++iterSHM) {
          simHitMapLeft[iterSHM->first] = iterSHM->second * Kf;
          simHitMapCentr[iterSHM->first] = iterSHM->second * chargeCentr / chargeTotal;
          simHitMapRight[iterSHM->first] = iterSHM->second * Kf;
        }

        // Time
        const double time = signal->getTime();
        // Left neighbour (crosstalk cannot go to non-existing strips, i.e. stripID >=0 && stripID < NSTRIPS
        if ((iStripLeft) >= 0) {
          if (stripMap[STRIP].find(iStripLeft) != stripMap[STRIP].end()) {
            stripMap[STRIP][iStripLeft]->updateCharge(chargeLeft);
            stripMap[STRIP][iStripLeft]->updateSimHitMap(simHitMapLeft);
          } else {
            stripMap[STRIP][iStripLeft] = new Signal(chargeLeft, time);
            stripMap[STRIP][iStripLeft]->updateSimHitMap(simHitMapLeft);
          }
        }
        // Central strip
        if (chargeCentr != 0) {
          if (stripMap[STRIP].find(iStrip) != stripMap[STRIP].end()) {
            stripMap[STRIP][iStrip]->updateCharge(chargeCentr);
            stripMap[STRIP][iStrip]->updateSimHitMap(simHitMapCentr);
          } else {
            stripMap[STRIP][iStrip] = new Signal(chargeCentr, time);
            stripMap[STRIP][iStrip]->updateSimHitMap(simHitMapCentr);
          }
        }

        // Right neighbour (crosstalk cannot go to non-existing strips, i.e. stripID >=0 && stripID < NSTRIPS
        if ((iStripRight) < _geometry->getSensorNStrips(layerID, sensorID)) {
          if (stripMap[STRIP].find(iStripRight) != stripMap[STRIP].end()) {
            stripMap[STRIP][iStripRight]->updateCharge(chargeRight);
            stripMap[STRIP][iStripRight]->updateSimHitMap(simHitMapRight);
          } else {
            stripMap[STRIP][iStripRight] = new Signal(chargeRight, time);
            stripMap[STRIP][iStripRight]->updateSimHitMap(simHitMapRight);
          }
        }

        // Release memory
        delete iterChMap->second;
      }
    }

    // Delete previous strip map
    delete[] iterSMap->second;

    // Save new results
    iterSMap->second = stripMap;
  } // For
}

//
// Method generating random noise using Gaussian distribution (input parameter:
// sensor map of strips with total integrated charge)
//
void SiStripDigi::genNoise(SensorStripMap& sensorMap) {
  // Add noise, only if set nonzero!
  if (_elNoise < 1e-4) {
    return;
  }

  // Strips type
  std::vector<StripType> stvec;
  stvec.push_back(STRIPFRONT);
  stvec.push_back(STRIPREAR);

  for (SensorStripMap::iterator iterSMap = sensorMap.begin(); iterSMap != sensorMap.end(); ++iterSMap) {
    for (std::vector<StripType>::iterator itType = stvec.begin(); itType != stvec.end(); ++itType) {
      for (StripChargeMap::iterator iterChMap = iterSMap->second[(*itType)].begin();
           iterChMap != iterSMap->second[(*itType)].end(); ++iterChMap) {
        const double elNoise = _genGauss->fire();

        iterChMap->second->updateCharge(elNoise);
      }
    }
  }
}

//
// Method transforming given SimTrackerHit into local ref. system of each sensor, where
// the system is positioned such as x, y and z coordinates are always positive;
// positive X corresponds to the side collecting holes and X = 0 corresponds to
// the strip collection electrons (parameter: a vector of pointers to
// SimTrackerDigiHits).
//
void SiStripDigi::transformSimHit(SimTrackerDigiHit* simDigiHit) {
  // Print sensor info
  streamlog_out(MESSAGE2) << " * Layer " << _geometry->getLayerRealID(_currentLayerID) << " "
                          << "Ladder " << _currentLadderID << " "
                          << "Sensor " << _currentSensorID << " " << std::fixed << std::setprecision(2) << std::endl
                          << "    Local B field: " << _magField / T << " T" << std::setprecision(0) << std::endl;

  _geometry->printSensorParams(_currentLayerID);

  // Print detailed info about hits in global ref. system
  printHitInfo("Hit in global ref. system", simDigiHit);

  // Hit global preStep, resp. posStep, positiona and momentum
  Hep3Vector globPrePosition = simDigiHit->get3PrePosition();
  Hep3Vector globPosPosition = simDigiHit->get3PosPosition();
  Hep3Vector globMomentum = simDigiHit->get3Momentum();
  // Hit local preStep, resp. posStep, position and momentum
  streamlog_out(DEBUG4) << "-- PRE-STEP Point" << std::endl;
  Hep3Vector locPrePosition =
      _geometry->transformPointToLocal(_currentLayerID, _currentLadderID, _currentSensorID, globPrePosition);
  streamlog_out(DEBUG4) << "-- POST-STEP Point" << std::endl;
  Hep3Vector locPosPosition =
      _geometry->transformPointToLocal(_currentLayerID, _currentLadderID, _currentSensorID, globPosPosition);

  Hep3Vector locMomentum =
      _geometry->transformVecToLocal(_currentLayerID, _currentLadderID, _currentSensorID, globMomentum);

  //-- Checking coherence
  // Avoid preStep rounding errors, checking inside sensitive
  const bool isPreIn =
      _geometry->isPointInsideSensor(_currentLayerID, _currentLadderID, _currentSensorID, locPrePosition);
  if (!isPreIn) {
    streamlog_out(ERROR) << "SiStripDigi::transformSimHit: "
                         << "Pre-step Position out of Sensor!\n"
                         << std::endl;
    exit(-1);
  }
  // Avoid posStep rounding errors, checking inside sensitive
  const bool isPosIn =
      _geometry->isPointInsideSensor(_currentLayerID, _currentLadderID, _currentSensorID, locPosPosition);
  if (!isPosIn) {
    streamlog_out(ERROR) << "SiStripDigi::transformSimHit: "
                         << "Post-step Position out of Sensor!\n"
                         << std::endl;
    exit(-1);
  }

  // Save hit local preStep, resp. posStep, position and momentum
  simDigiHit->set3PrePosition(locPrePosition);
  simDigiHit->set3PosPosition(locPosPosition);
  simDigiHit->set3Momentum(locMomentum);

  // Print detailed info about hit in local ref. system
  printHitInfo("Hit in local ref. system", simDigiHit);
}

// GET METHODS

//
// Method returning electron diffusivity (parameters: position in cm)
//
double SiStripDigi::getElecDiffusivity(double pos) {
  // Electron mobility
  double mobility = getElecMobility(pos);

  return ((k * _temp) / ePlus * mobility);
}

//
// Method returning hole diffusivity (parameters: position in cm)
//
double SiStripDigi::getHoleDiffusivity(double pos) {
  // Hole mobility
  double mobility = getHoleMobility(pos);

  return ((k * _temp) / ePlus * mobility);
}

//
// Method returning electron drift time (parameters: position in cm)
//
double SiStripDigi::getElecDriftTime(double pos) {
  // Set pointer to velocity function and create instance of IntSolver
  double (SiStripDigi::*pfce)(double) = &SiStripDigi::getElecInvVelocity;
  RombIntSolver<SiStripDigi> elecIntSolver(pfce, this, _epsTime);

  if (pos > _sensorThick) {
    streamlog_out(ERROR) << "Problem to calculate total drift time. "
                         << "Electrons at position: " << pos << " are out of range!" << std::endl;
    exit(-1);
  }

  // Result (Be carefull about rounding errors)
  if ((_sensorThick - pos) <= ROUNDEPS * um) {
    return 0.;
  } else {
    return (elecIntSolver.Integrate(pos, _sensorThick));
  }
}

//
// Method returning hole drift time (parameters: position X in cm)
//
double SiStripDigi::getHoleDriftTime(double pos) {
  // Set pointer to velocity function and create instance of IntSolver
  double (SiStripDigi::*pfce)(double) = &SiStripDigi::getHoleInvVelocity;
  RombIntSolver<SiStripDigi> holeIntSolver(pfce, this, _epsTime);

  if (pos < -ROUNDEPS * um) {
    streamlog_out(ERROR) << "Problem to calculate total drift time."
                         << " Holes at position: " << std::setprecision(5) << pos << " are out of range!" << std::endl;
    exit(-1);
  }

  // Result (Be carefull about rounding errors)
  if ((pos) <= ROUNDEPS * um) {
    return 0.;
  } else {
    return (holeIntSolver.Integrate(pos, 0.));
  }
}

//
// Method returning electric intensity in V/cm (parameters: position X in cm)
//
double SiStripDigi::getEField(double pos) {
  // Si wafer thickness in cm (Added the rounding error)
  if (pos < -ROUNDEPS * um || pos > _sensorThick + ROUNDEPS * um) {
    streamlog_out(ERROR) << std::setprecision(3) << "Electric field at required position[um]: " << pos / um
                         << " not defined!" << std::endl;
    exit(-1);
  }

  // Return electric intensity
  return (-(_Vbias + _Vdepl) / _sensorThick + 2 * pos / (_sensorThick * _sensorThick) * _Vdepl);
}

//
// Method returning electron mobility in cm^2/V.s (parameters: position X in cm)
//
double SiStripDigi::getElecMobility(double pos) {
  // Electron parameters - maximum velocity, critical intenzity, beta factor
  static double vmElec = 1.53 * pow(_temp, -0.87) * 1.E9 * cm / s;
  static double EcElec = 1.01 * pow(_temp, +1.55) * V / cm;
  static double betaElec = 2.57 * pow(_temp, +0.66) * 1.E-2;

  // Absolute value of electric intenzity
  double E = getEField(pos);
  if (E < 0.)
    E = -E;

  // Return mobility
  return (vmElec / EcElec * 1. / (pow(1. + pow((E / EcElec), betaElec), (1. / betaElec))));
}

//
// Method returning hole mobility in cm^2/V.s (parameters: position X in cm)
//
double SiStripDigi::getHoleMobility(double pos) {
  // Hole parameters - maximum velocity, critical intenzity, beta factor
  static double vmHole = 1.62 * pow(_temp, -0.52) * 1.E8 * cm / s;
  static double EcHole = 1.24 * pow(_temp, +1.68) * V / cm;
  static double betaHole = 0.46 * pow(_temp, +0.17);

  // Absolute value of electric intenzity
  double E = getEField(pos);
  if (E < 0.) {
    E = -E;
  }

  // Return mobility
  return (vmHole / EcHole * 1. / (pow(1. + pow((E / EcHole), betaHole), (1. / betaHole))));
}

//
// Method that calculates Lorentz angle for electrons (parameters: position X
// in cm)
//
Hep3Vector SiStripDigi::getElecLorentzShift(double pos) {
  // Set pointer to mobility function and create static instance to IntSolver
  double (SiStripDigi::*pfce)(double) = &SiStripDigi::getElecMobility;
  RombIntSolver<SiStripDigi> elecIntSolver(pfce, this, _epsAngle);

  // Hall scattering factor for electrons
  static float rElec = 1.13 + 0.0008 * (_temp - 273);

  // Si wafer thickness in cm
  if (pos > _sensorThick) {
    streamlog_out(ERROR) << "Problem to calculate Lorentz angle. Electrons at position: " << pos << " are out of range!"
                         << std::endl;
    exit(-1);
  }

  // Final Lorentz shift for electrons (Be careful about rounding errors)
  Hep3Vector shiftLorentz(0., 0., 0.);

  if ((_sensorThick - pos) >= ROUNDEPS * um) {
    double integral = elecIntSolver.Integrate(pos, _sensorThick);

    shiftLorentz.setY(integral * rElec * _magField.getZ() * -1.); //????
    shiftLorentz.setZ(integral * rElec * _magField.getY());       //????

    // std::cout << "Mag field: " << _magField/T << " Lorentz shift: " << shiftLorentz/(_sensorThick-pos) << " " <<
    // (_sensorThick-pos)/um << " " << shiftLorentz << std::endl;
  }

  return (shiftLorentz);
}

//
// Method that calculates Lorentz angle for holes (parameters: position X in cm)
//
Hep3Vector SiStripDigi::getHoleLorentzShift(double pos) {
  // Set pointer to mobility function and create static instance to IntSolver
  double (SiStripDigi::*pfce)(double) = &SiStripDigi::getHoleMobility;
  RombIntSolver<SiStripDigi> holeIntSolver(pfce, this);

  // Hall scattering factor for holes
  static float rHole = 0.72 - 0.0005 * (_temp - 273);

  // Si wafer thickness in cm
  if (pos < -ROUNDEPS * um) {
    streamlog_out(ERROR) << "Problem to calculate Lorentz angle."
                         << " Holes at position: " << pos << " are out of range!" << std::endl;
    exit(-1);
  }

  // Final Lorentz shift for holes (Be careful about rounding errors)
  Hep3Vector shiftLorentz(0., 0., 0.);
  if (pos >= ROUNDEPS * um) {
    double integral = holeIntSolver.Integrate(0., pos);

    shiftLorentz.setY(integral * rHole * _magField.getZ());       //????
    shiftLorentz.setZ(integral * rHole * _magField.getY() * -1.); //????
  }

  return (shiftLorentz);
}

//
// Method that returns actual electron velocity (parameters: position in cm)
//
double SiStripDigi::getElecVelocity(double pos) { return (-1. * getElecMobility(pos) * getEField(pos)); }

//
// Method that returns actual electron inverse velocity (parameters: position in cm)
//
double SiStripDigi::getElecInvVelocity(double pos) { return (-1. / (getElecMobility(pos) * getEField(pos))); }

//
// Method that returns actual hole velocity (parameters: position in cm)
//
double SiStripDigi::getHoleVelocity(double pos) { return (getHoleMobility(pos) * getEField(pos)); }

//
// Method that returns actual hole inverse velocity (parameters: position in cm)
//
double SiStripDigi::getHoleInvVelocity(double pos) { return (1. / (getHoleMobility(pos) * getEField(pos))); }

// PRINT METHODS

//
// Method printing cluster info
//
void SiStripDigi::printClusterInfo(const DigiCluster& cluster) const {
  streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal | std::ios::showpos)
                          << std::setprecision(2) << "    Cluster:"
                          << " Q [fC]: " << cluster.getNCarriers() * cluster.getCharge() / fC << std::setprecision(1)
                          << " Pos X [um]: " << std::setw(6) << cluster.getPosX() / um << std::setprecision(3)
                          << " Pos Y [mm]: " << std::setw(6) << cluster.getPosY() / mm
                          << " Pos Z [mm]: " << std::setw(6) << cluster.getPosZ() / mm
                          << std::resetiosflags(std::ios::internal | std::ios::showpos) << std::setprecision(0)
                          << std::endl;
}

//
// Method printing clusters info (clusters in std::vector)
//
void SiStripDigi::printClustersInfo(std::string info, const DigiClusterVec& clusterVec) const {
  DigiClusterVec::const_iterator iterClusters;

  streamlog_out(MESSAGE1) << "   " << info << ": " << clusterVec.size() << " clusters" << std::endl;
  for (iterClusters = clusterVec.begin(); iterClusters != clusterVec.end(); ++iterClusters) {
    printClusterInfo(**iterClusters);
  }
}

//
// Method printing hit info
//
void SiStripDigi::printHitInfo(std::string info, const SimTrackerDigiHit* hit) const {
  streamlog_out(MESSAGE2) << "   " << info << ":" << std::endl;
  streamlog_out(MESSAGE2) << std::fixed << std::setprecision(1) << "    Hit:"
                          << " DepE [keV]: " << hit->getEDep() / keV << std::setprecision(3)
                          << std::setiosflags(std::ios::showpos)
                          << " Pos X [mm]: " << (hit->getPreX() + hit->getPosX()) / 2. / mm
                          << " Pos Y [mm]: " << (hit->getPreY() + hit->getPosY()) / 2. / mm
                          << " Pos Z [mm]: " << (hit->getPreZ() + hit->getPosZ()) / 2. / mm << " Path [mm]: "
                          << hit->getPathLength() / mm
                          //<< " NDaughters: "   << hit->getMCParticle()->getDaughters().size()
                          //<< " VertexZ: "      << *(hit->getMCParticle()->getVertex()+2)
                          << std::resetiosflags(std::ios::showpos) << std::setprecision(0) << std::endl;
}

//
// Method printing processor parameters
//
void SiStripDigi::printProcessorParams() const {
  streamlog_out(MESSAGE3) << std::endl
                          << " " << DUNDERL << DBLUE << "SiStripDigi parameters:" << ENDCOLOR << " " << std::endl
                          << std::endl;

  streamlog_out(MESSAGE3) << std::setiosflags(std::ios::fixed | std::ios::internal) << std::setprecision(2)
                          << "  Sensor depletion voltage [V]:          " << std::setw(6) << _Vdepl / V << std::endl
                          << "  Sensor bias voltage [V]:               " << std::setw(6) << _Vbias / V << std::endl
                          << "  Temperature set on sensor [K]:         " << std::setw(6) << _temp / K << std::endl;
  if (_electronicEffects)
    streamlog_out(MESSAGE3) << "  Mutual interstrip capacitance [pF]:    " << std::setw(6) << _capInterStrip
                            << std::endl
                            << "  Strip-to-backplane capacitance [pF]:   " << std::setw(6) << _capBackPlane << std::endl
                            << "  AC coupling - capacitance [pF]:        " << std::setw(6) << _capCoupl << std::endl
                            << "  Electronics noise - ENC [fC]:          " << std::setw(6) << _elNoise / fC << std::endl
                            << std::endl;
  if (_floatStripsRPhi)
    streamlog_out(MESSAGE3) << "  Read-out R-Phi pitch is 2x geom. pitch." << std::endl;
  if (_floatStripsZ)
    streamlog_out(MESSAGE3) << "  Read-out Z     pitch is 2x geom. pitch." << std::endl;
  if (_landauFluct) {
    streamlog_out(MESSAGE3) << "                                         " << std::endl
                            << "  Internal Landau fluctuations?:         " << "   yes" << std::endl
                            << "  Prod. threhold on second. e [keV]:     " << std::setw(6)
                            << _prodThreshOnDeltaRays / keV << std::endl
                            << "  Beta*Gamma limit for Landau fluct:     " << std::setw(6) << _landauBetaGammaCut
                            << std::endl;
  } else {
    streamlog_out(MESSAGE3) << "                                         " << std::endl
                            << "  Internal Landau fluctuations?:         " << "    no" << std::endl;
  }
  streamlog_out(MESSAGE3) << "                                         " << std::endl
                          << "  Digi precision in space [um]:          " << std::setw(6) << _epsSpace / um << std::endl
                          << "  Digi rel. precision in Lorentz angle:  " << std::setw(6) << _epsAngle << std::endl
                          << "  Digi rel. precision in Drift time:     " << std::setw(6) << _epsTime << std::endl
                          << std::resetiosflags(std::ios::showpos) << std::setprecision(0) << std::endl;
}

//
// Method printing info about signals at each strip
//
void SiStripDigi::printStripsInfo(std::string info, const SensorStripMap& sensorMap) const {
  // Sensor map of strips with total integrated charge
  SensorStripMap::const_iterator iterSMap;
  StripChargeMap::const_iterator iterChMap;

  int cellID = 0;

  short int layerID = 0;
  short int ladderID = 0;
  short int sensorID = 0;

  int stripID = 0;

  streamlog_out(MESSAGE1) << "  Digi results - " << info << ":" << std::endl;

  for (iterSMap = sensorMap.begin(); iterSMap != sensorMap.end(); ++iterSMap) {

    cellID = iterSMap->first;
    std::map<std::string, int> bfMap = _geometry->decodeCellID(cellID);
    layerID = bfMap["layer"];
    ladderID = bfMap["module"];
    sensorID = bfMap["sensor"];
    //_geometry->decodeCellID(layerID, ladderID, sensorID, cellID);
    layerID = _geometry->getLayerRealID(layerID);

    streamlog_out(MESSAGE1) << std::endl
                            << "   Layer: " << layerID << " Ladder: " << ladderID << " Sensor: " << sensorID
                            << std::endl;

    // Strips in R-Phi
    for (iterChMap = iterSMap->second[STRIPRPHI].begin(); iterChMap != iterSMap->second[STRIPRPHI].end(); ++iterChMap) {

      stripID = iterChMap->first;
      streamlog_out(MESSAGE1) << "    Strip number in R-Phi: " << stripID
                              << std::setiosflags(std::ios::fixed | std::ios::internal) << std::setprecision(2)
                              << " Total charge [fC]: "
                              << iterChMap->second->getCharge() / fC
                              //<< " Generation time [ns]: " << iterChMap->second->getTime()/ns
                              << std::setprecision(0) << std::endl;
    }

    // Strips in Z
    for (iterChMap = iterSMap->second[STRIPZ].begin(); iterChMap != iterSMap->second[STRIPZ].end(); ++iterChMap) {

      stripID = iterChMap->first;
      streamlog_out(MESSAGE1) << "    Strip number in Z: " << stripID
                              << std::setiosflags(std::ios::fixed | std::ios::internal) << std::setprecision(2)
                              << " Total charge [fC]: "
                              << iterChMap->second->getCharge() / fC
                              //<< " Generation time [ns]: " << iterChMap->second->getTime()/ns
                              << std::setprecision(0) << std::endl;
    }
    streamlog_out(MESSAGE1) << std::endl;

  } // For
}

} // namespace sistrip
