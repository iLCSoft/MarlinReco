/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "VTXDigitizer.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCFlagImpl.h>
#include "random.h"
#include "gsl/gsl_sf_erf.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h" 
#include "CLHEP/Random/RandFlat.h"


// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>

namespace CLHEP{}    // declare namespace CLHEP for backward compatibility
using namespace CLHEP ;

using namespace lcio ;
using namespace marlin ;
using namespace std ;

typedef std::vector<SimTrackerHit*> SimTrackerHitVec;

VTXDigitizer aVTXDigitizer ;


VTXDigitizer::VTXDigitizer() : Processor("VTXDigitizer") {
  
  // modify processor description
  _description = "VTXDigitizer should create VTX TrackerHits from SimTrackerHits" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "CollectionName" , 
                           "Name of the SimTrackerHit collection"  ,
                           _colName ,
                           std::string("vxd01_VXD") ) ;

  registerOutputCollection( LCIO::TRACKERHIT,
                            "OutputCollectionName" , 
                            "Name of the output TrackerHit collection"  ,
                            _outputCollectionName ,
                            std::string("VTXTrackerHits") ) ;
  
  registerOutputCollection( LCIO::LCRELATION,
                            "RelationColName" , 
                            "Name of the output VTX trackerhit relation collection"  ,
                            _colVTXRelation ,
                            std::string("VTXRelation") ) ;

  registerProcessorParameter("TanLorentz",
                             "Tangent of Lorentz Angle",
                             _tanLorentzAngle,
                             (double)0.8);

  registerProcessorParameter("CutOnDeltaRays",
                             "Cut on delta-ray energy (MeV)",
                             _cutOnDeltaRays,
                             (double)0.030);

  registerProcessorParameter("Diffusion",
                             "Diffusion coefficient (in mm) for layer thickness",
                             _diffusionCoefficient,
                             (double)0.002);

  

  registerProcessorParameter("PixelSizeX",
                             "Pixel Size X",
                             _pixelSizeX,
                             (double)0.025);


  registerProcessorParameter("PixelSizeY",
                             "Pixel Size Y",
                             _pixelSizeY,
                             (double)0.025);

  registerProcessorParameter("Debug",
                             "Debug option",
                             _debug,
                             int(0));


   
  registerProcessorParameter("ElectronsPerKeV",
                             "Electrons per keV",
                             _electronsPerKeV,
                             (double)270.3);




 
  std::vector<float> bkgdHitsInLayer;
  bkgdHitsInLayer.push_back(34400.);
  bkgdHitsInLayer.push_back(23900.);
  bkgdHitsInLayer.push_back(9600.);
  bkgdHitsInLayer.push_back(5500.);
  bkgdHitsInLayer.push_back(3100.);    
  registerProcessorParameter("BackgroundHitsPerLayer",
                             "Background Hits per Layer",
                             _bkgdHitsInLayer,
                             bkgdHitsInLayer);

  registerProcessorParameter("SegmentLength",
                             "Segment Length",
                             _segmentLength,
                             double(0.005));

  registerProcessorParameter("WidthOfCluster",
                             "Width of cluster",
                             _widthOfCluster,
                             double(3.0));


  registerProcessorParameter("Threshold",
                             "Cell Threshold in electrons",
                             _threshold,
                             200.);

  registerProcessorParameter("PoissonSmearing",
                             "Apply Poisson smearing of electrons collected on pixels",
                             _PoissonSmearing,
                             1);

  registerProcessorParameter("ElectronicEffects",
                             "Apply Electronic Effects",
                             _electronicEffects,
                             int(1));

  registerProcessorParameter("ElectronicNoise",
                             "electronic noise in electrons",
                             _electronicNoise,
                             100.);

  registerProcessorParameter("AdditionalCollections",
                             "Additional Collections to store hit position in the local ladder frame",
                             _produceFullPattern,
                             int(0));

  registerProcessorParameter("UseMCPMomentum",
                             "Use Particle Momentum",
                             _useMCPMomentum,
                             int(1));

  registerProcessorParameter("EnergyLoss",
                             "Energy Loss keV/mm",
                             _energyLoss,
                             double(280.0));

  registerProcessorParameter("RemoveDRayPixels",
                             "Remove D-Ray Pixels",
                             _removeDrays,
                             int(1));

  registerProcessorParameter("GenerateBackground",
                             "Generate Background",
                             _generateBackground,
                             int(0));

}




void VTXDigitizer::init() { 

  // usually a good idea to
  printParameters();
  // internal parameters
  PI = (double)acos((double)(-1.0));
  TWOPI = (double)(2.0)*PI;	
  PI2 = 0.5*PI;
  _nRun = 0 ;
  _nEvt = 0 ;
  _totEntries = 0;
  _fluctuate = new MyG4UniversalFluctuationForSi();
}

void VTXDigitizer::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;

  //------Get the geometry from the gear file-----//

  const gear::VXDParameters& gearVXD = Global::GEAR->getVXDParameters() ;
  const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout(); 
 
 //number of layers
  _numberOfLayers  = layerVXD.getNLayers();

  //the number of ladders within layers
  _laddersInLayer.resize(_numberOfLayers);

 //Azimuthal offset of the whole structure of each layer.
  _layerHalfPhi.resize(_numberOfLayers);
 
  //layer thickness and half-thickness
  _layerHalfThickness.resize(_numberOfLayers);
  _layerThickness.resize(_numberOfLayers);

  //Distance from the middle of the Si sensor to IP in each layer.
  _layerRadius.resize(_numberOfLayers);

 //The length of the Si sensor in each layer
  _layerLadderLength.resize(_numberOfLayers);

  //The width and half-width of the Si sensor in each layer
  _layerLadderWidth.resize(_numberOfLayers);
  _layerLadderHalfWidth.resize(_numberOfLayers);

 //The offset of the sensitive area in ladder in each layer
 _layerActiveSiOffset.resize(_numberOfLayers);

 // The gaps in z between two subladders within each layer
  const gear::GearParameters& gearVXDInfra = Global::GEAR->getGearParameters("VXDInfra") ;
  const std::vector<double> laddergaps = gearVXDInfra.getDoubleVals("LadderGaps");
  std::cout<<"laddergaps size "<<laddergaps.size()<<std::endl;
  _layerLadderGap.resize(laddergaps.size());

  //ladder offset in phi
 _layerPhiOffset.resize(_numberOfLayers);

  for(int layer = 0; layer < _numberOfLayers; layer++)
    {
      _laddersInLayer[layer] = layerVXD.getNLadders(layer);  
      _layerHalfPhi[layer] = layerVXD.getPhi0(layer); 

      //should half phi be replaced by this from Alexei's original code?
      _layerHalfPhi[layer] = PI/((double)_laddersInLayer[layer]);

      _layerThickness[layer] =layerVXD.getSensitiveThickness(layer);
      _layerHalfThickness[layer] = 0.5*_layerThickness[layer];      
      _layerRadius[layer] =layerVXD.getSensitiveDistance(layer) + 0.5 * _layerThickness[layer];
      _layerLadderLength[layer] = 2*layerVXD.getSensitiveLength(layer);
      _layerLadderWidth[layer] = layerVXD.getSensitiveWidth(layer);
      _layerLadderHalfWidth[layer] = _layerLadderWidth[layer]/2.;
      _layerActiveSiOffset[layer] = - (layerVXD.getSensitiveOffset(layer));
      _layerLadderGap[layer] = laddergaps[layer];
      _layerPhiOffset[layer] = layerVXD.getPhi0(layer);
    }
 
 
 
  cout<<" _numberOfLayers "<<_numberOfLayers<<endl;
  cout<<" _pixelSizeX "<<_pixelSizeX<<endl;
  cout<<" _pixelSizeY "<<_pixelSizeY<<endl;
  cout<<" _electronsPerKeV "<<_electronsPerKeV<<endl;
  cout<<" _segmentDepth "<<_segmentDepth<<endl;
  cout<<" _currentTotalCharge "<<_currentTotalCharge<<endl;
  for (int i=0; i<_numberOfLayers; ++i) 
    {
      cout<<"layer "<<i<<endl;
      cout<<" _laddersInLayer "<<_laddersInLayer[i]<<endl;
      cout<<" _layerRadius "<<_layerRadius[i]<<endl;
      cout<<" _layerLadderLength "<<_layerLadderLength[i]<<endl;
      cout<<" _layerLadderHalfWidth "<<_layerLadderHalfWidth[i]<<endl;
      cout<<" _layerPhiOffset "<<_layerPhiOffset[i]<<endl;
      cout<<" _layerActiveSiOffset "<< _layerActiveSiOffset[i]<<endl;
      cout<<" _layerHalfPhi "<<_layerHalfPhi[i]<<endl;
      cout<<" _layerLadderGap "<<_layerLadderGap[i]<<endl;
      cout<<" _bkgdHitsInLayer "<<_bkgdHitsInLayer[i]<<endl;
      cout<<" _layerLadderWidth "<<_layerLadderWidth[i]<<endl;
      cout<<" _layerThickness "<<_layerThickness[i]<<endl;
      
      cout<<" _layerHalfThickness "<<_layerHalfThickness[i]<<endl;


    }


  SCALING = 25000.;
} 

void VTXDigitizer::processEvent( LCEvent * evt ) { 

  try{
    LCCollection * STHcol = evt->getCollection( _colName ) ;
 //   LCFlagImpl flag;
//    flag.setBit(LCIO::THBIT_MOMENTUM);
//    STHcol->setFlag(flag.getFlag());
    LCCollectionVec * THcol = new LCCollectionVec(LCIO::TRACKERHIT);
    LCCollectionVec * RelCol = new LCCollectionVec(LCIO::LCRELATION);
    LCCollectionVec * STHLocCol = NULL;
    LCCollectionVec * THLocCol = NULL;
    LCCollectionVec * RelLocCol = NULL;
    if (_produceFullPattern == 1) {
      STHLocCol = new LCCollectionVec(LCIO::SIMTRACKERHIT);
      THLocCol = new LCCollectionVec(LCIO::TRACKERHIT);
      RelLocCol = new LCCollectionVec(LCIO::LCRELATION);
    }




    int nSTH = STHcol->getNumberOfElements();
    // Loop over sim tracker hits;
    for (int i=0; i<nSTH; ++i) {
      //  std::cout << "Hit number : " << i << std::endl;
      SimTrackerHit * simTrkHit = 
        dynamic_cast<SimTrackerHit*>(STHcol->getElementAt(i));
 //     _momX = simTrkHit->getMomentum()[0];
 //     _momY = simTrkHit->getMomentum()[1];
 //     _momZ = simTrkHit->getMomentum()[2];
 //     _eDep = simTrkHit->getdEdx();
 //     float totPartMomentum = sqrt(_momX*_momX+_momY*_momY+_momZ*_momZ);

      // do digitization for one SimTrackerHit 
      // Produce ionisation points along track
      // std::cout << "Beginning of Cycle " << std::endl;
      ProduceIonisationPoints( simTrkHit );      
      // std::cout << "End of ProduceIonisationPoints( ) " << std::endl;
      // std::cout << "Current layer = " << _currentLayer << std::endl;
      // std::cout << "Total particle momentum = " << totPartMomentum << std::endl;
      bool accept = _currentLayer >= 0 && _currentLayer < int(_layerRadius.size());
 //     if (_removeDrays == 1)
 //       accept = accept && totPartMomentum > 10;
      if ( accept ) {
        // Produce signal points on collection plane
        ProduceSignalPoints( );
        // std::cout << "End of ProduceSignalPoints( ) " << std::endl;
        SimTrackerHitImplVec simTrkHitVec;
        // Produce fired pixels 
        ProduceHits( simTrkHitVec );

   //      for (int iHit=0; iHit<int(simTrkHitVec.size()); ++iHit) {
//           SimTrackerHit * hit = simTrkHitVec[iHit];
//           //cout<<"hit 1"<<iHit <<" "<<hit->getdEdx()<<endl;
//         }
        // std::cout << "E


        // std::cout << "End of ProduceHits( ) " << std::endl;
        // Apply Poisson Smearing to deposited charges
        if (_PoissonSmearing != 0) PoissonSmearer( simTrkHitVec );
        if (_electronicEffects != 0) GainSmearer( simTrkHitVec );
        //std::cout << "Amplitude 6 = " << _ampl << std::endl;
        TrackerHitImpl * recoHit = ReconstructTrackerHit( simTrkHitVec );
        //std::cout << "Amplitude 7 = " << _ampl << std::endl;
        // std::cout << "End of ReconstructTrackerHit( ) " << std::endl;
        if (_produceFullPattern == 1 && recoHit !=0 ) {
          SimTrackerHitImpl * sth = new SimTrackerHitImpl();
          sth->setCellID(simTrkHit->getCellID());
          sth->setdEdx(simTrkHit->getdEdx());
          double currentPosition[3];
          for (int j=0; j<3; ++j)
            currentPosition[j] = 0.5*(_currentExitPoint[j] + _currentEntryPoint[j]);
          sth->setPosition(currentPosition);
          TrackerHitImpl * th = new TrackerHitImpl();
          th->setdEdx(recoHit->getdEdx());
          double xp[3];
          for (int j=0; j<3; ++j)
            xp[j] = recoHit->getPosition()[j];
          th->setPosition(xp);          
          STHLocCol->addElement(sth); 
          THLocCol->addElement(th);
          LCRelationImpl * rel = new LCRelationImpl(th,sth,float(1.0));
          RelLocCol->addElement(rel);
        }
        if (recoHit != NULL) {
          TrackerHitToLab( recoHit );
        }
        if (_debug != 0) {
          if (recoHit == NULL) 
           std::cout << "Number of pixels above threshold = 0 " << std::endl;
          else
            PrintInfo( simTrkHit, recoHit );
        }

        if ( recoHit != NULL) {


          recoHit->rawHits().push_back(simTrkHit);


          float pointResoRPhi=0.004;
          float pointResoZ=0.004;
          float covMat[TRKHITNCOVMATRIX]={0.,0.,pointResoRPhi*pointResoRPhi,0.,0.,pointResoZ*pointResoZ};
          recoHit->setCovMatrix(covMat);      


          recoHit->setType(100+simTrkHit->getCellID());
          THcol->addElement( recoHit );
          //          std::cout << "Hit is added to collection " << _nEvt << std::endl;
          LCRelationImpl * rel = new LCRelationImpl(recoHit,simTrkHit,float(1.0));
          RelCol->addElement(rel);
        }
// Clean Up        
        for (int i=0; i < int(simTrkHitVec.size()); ++i) {
          SimTrackerHit * hit = simTrkHitVec[i];
          delete hit;
        }     
      } 
    }
    if (_generateBackground == 1) 
      generateBackground( THcol );

    evt->addCollection(THcol,_outputCollectionName.c_str());
    evt->addCollection(RelCol, _colVTXRelation );
    if (_produceFullPattern == 1) {
      evt->addCollection(STHLocCol,"VTXLocalSimTrackerHits");
      evt->addCollection(THLocCol,"VTXLocalTrackerHits");
      evt->addCollection(RelLocCol,"VTXLocalRelation");
    }
  }
  catch(DataNotAvailableException &e){}

    
  _nEvt ++ ;
  
  if (_nEvt % 100 == 0)
   std::cout << "Processed " << _nEvt << "events " << std::endl;

}



void VTXDigitizer::check( LCEvent * evt ) { 

}


void VTXDigitizer::end(){ 
  
  std::cout << "VTXDigitizer::end()  " << name() 
             << " processed " << _nEvt << " events in " << _nRun << " runs "
             << std::endl ;

  delete _fluctuate;

}

void VTXDigitizer::FindLocalPosition(SimTrackerHit * hit,
                                         double * localPosition,
                                         double * localDirection) { 
  /** Function calculates local coordinates of the sim hit 
   * in the given ladder and local momentum of particle. 
   * Also returns module number and ladder 
   * number.
   * Local coordinate system within the ladder 
   * is defined as following :  <br> 
   *    - x axis lies in the ladder plane and orthogonal to the beam axis <br>
   *    - y axis is perpendicular to the ladder plane <br>
   *    - z axis lies in the ladder plane and parallel to the beam axis <br>
   * 
   *    Encoding of modules: <br>
   *    ======================  <br>
   *    - 0 = left endcap <br>
   *    - 1 = left ladder in the barrel <br>
   *    - 2 = right ladder in the barrel <br>
   *    - 3 = right endcap <br>
   * 
   */

 
  double xLab[3] = { hit->getPosition()[0],
                     hit->getPosition()[1],
                     hit->getPosition()[2]};

  int layer = -1;

// FIXME : Assume for the moment only barrel in VDX

// Find layer number by coordinates 

  double RXY = sqrt(xLab[0]*xLab[0]+
                    xLab[1]*xLab[1]);      
/*
  for (int i=0; i<_numberOfLayers; ++i) {
    double xmin = _layerRadius[i] - _layerThickness;
    double xmax;
    if (i<_numberOfLayers-1)
      xmax = _layerRadius[i+1] - _layerThickness;
    else
      xmax = 2.0*_layerRadius[i];
    // LADDERED STRUCTURE IMPLIES NUMBER OF LADDERS > 2 !!!!!
    if (_laddersInLayer[i] > 2) {
      xmax = xmax/(double)cos(_layerHalfPhi[i]);
    }
    if (RXY > xmin && RXY < xmax) {
      layer = i;
      break;
    }
  }
*/


// layer is encoded in CellID; 
   layer = hit->getCellID() - 1;
                                          
  _currentLayer = layer;

// Check boundary of layers
  if (layer < 0 || layer > _numberOfLayers) 
    return;

// Find module by z coordinate
// -z : module 1
// +z : module 2 
  int module;
  if (xLab[2] < 0.0 ) {
    module = 1;
  }
  else {
    module = 2;
  }
  _currentModule = module;

  double Momentum[3];
  if (hit->getMCParticle()) {
    for (int j=0; j<3; ++j)
      Momentum[j] = hit->getMCParticle()->getMomentum()[j];
  }
  else {
    for (int j=0; j<3; ++j) {
      Momentum[j] = hit->getMomentum()[j];

    }
  }

  _currentParticleMass = 0;
  if (hit->getMCParticle())
    _currentParticleMass   = hit->getMCParticle()->getMass();
  if (_currentParticleMass < 0.510e-3) 
    _currentParticleMass = 0.510e-3;  
  _currentParticleMomentum = 0.0;
  for (int i=0; i<3; ++i)
    _currentParticleMomentum += Momentum[i]*Momentum[i];
  _currentParticleMomentum = sqrt(_currentParticleMomentum);
  


  double PXY = sqrt(Momentum[0]*Momentum[0]+
                    Momentum[1]*Momentum[1]);

  double PhiInLab = (double)atan2(xLab[1],xLab[0]);
  if (PhiInLab < 0.0) PhiInLab += TWOPI;
  double PhiInLabMom = atan2(Momentum[1],Momentum[0]);
  if (PhiInLabMom < 0.0) PhiInLabMom += TWOPI;
  double Radius = _layerRadius[layer];
  //  << " Radius = " << Radius << std::endl;

  double Phi0 = _layerPhiOffset[layer];

  int nLadders = _laddersInLayer[layer];

  double dPhi = 2.0*_layerHalfPhi[layer];
  // And now compute local coordinates
  // and local momentum of particle

  double PhiLadder=0;
  double PhiInLocal=0;
  //cout<<"nLadders "<<nLadders<<" "<<dPhi<<" "<<Phi0<<" "<<endl;
  
  if (nLadders > 2) { // laddered structure
    //std::cout<<"laddered structure "<<std::endl;
    int iLadder=0;
    for (int ic=0; ic<nLadders; ++ic) {
      //      PhiLadder = - PI2 + double(ic)*dPhi + Phi0;
      PhiLadder = double(ic)*dPhi + Phi0;
      PhiInLocal = PhiInLab - PhiLadder;
      //cout<<"Phi "<<PhiLadder<<" "<<PhiInLocal<<" "<<PhiInLab<<" "<<_layerThickness[layer]<<" "<<Radius<<endl;
      if (RXY*cos(PhiInLocal)-Radius > -_layerThickness[layer] && 
          RXY*cos(PhiInLocal)-Radius < _layerThickness[layer]) {
        iLadder = ic;
        break;
      }
      //cout<<"phi ladder "<<PhiLadder<<endl;
    }
    double PhiLocalMom = PhiInLabMom - PhiLadder;
    localPosition[0] = RXY*sin(PhiInLocal);
    localPosition[1] = xLab[2];
    localPosition[2] = RXY*cos(PhiInLocal)-Radius;
    localDirection[0]=PXY*sin(PhiLocalMom);
    localDirection[1]=Momentum[2];
    localDirection[2]=PXY*cos(PhiLocalMom);
    _currentPhi = PhiLadder;
    //cout<<"local direction "<<localDirection[0]<<" "<<localDirection[1]<<" "<<localDirection[2]<<endl;
    //cout<<"phi local mom "<<PhiLocalMom<<" "<<PhiInLabMom<<" "<<PhiLadder<<endl;
  }  
  else { // cyllindrical structure
    //std::cout<<"cyllindrical structure "<<std::endl;
    localPosition[0]=0.0;
    localPosition[1]=xLab[2];
    localPosition[2]=RXY-Radius;
    double PhiLocalMom = PhiInLabMom - PhiInLab;
    localDirection[0]=PXY*sin(PhiLocalMom);
    localDirection[1]=Momentum[2];
    localDirection[2]=PXY*cos(PhiLocalMom); 
    _currentPhi = PhiInLab;
  }

 
}

void VTXDigitizer::TransformToLab(double * xLoc, double * xLab) {
  /** Function transforms local coordinates in the ladder
   * into global coordinates
   */

  int layer = _currentLayer;
  int nLadders = _laddersInLayer[layer];
  double Phi = _currentPhi;
  double Radius = _layerRadius[layer];

  if (nLadders > 2 ) { // laddered structure
    double baseLine = Radius + xLoc[2];
    double PhiInLab = Phi + atan2(xLoc[0],baseLine);
    double RXY = sqrt(baseLine*baseLine+xLoc[0]*xLoc[0]);
    xLab[2] = xLoc[1];
    xLab[0] = RXY*cos(PhiInLab);
    xLab[1] = RXY*sin(PhiInLab);
  }
  else { // cyllindrical structure
    double baseLine = Radius + xLoc[2];
    double PhiInLab = Phi + xLoc[0]/baseLine;
    xLab[0] = baseLine*cos(PhiInLab);
    xLab[1] = baseLine*sin(PhiInLab);
    xLab[2] = xLoc[1];    
  } 

}

void VTXDigitizer::ProduceIonisationPoints( SimTrackerHit * hit) {
  /** Produces ionisation points along track segment within active Silicon layer.
   */
  //std::cout << "Amplitude 1 = " << _ampl << std::endl;
  // Position of hit in the Lab frame
  double pos[3];
  double dir[3];
  double entry[3];
  double exit[3];

  // Find local position in the ladder frame
  FindLocalPosition( hit, pos, dir);
  //  std::cout << "End of FindLocalPosition" << std::endl;

  // check layer boundaries
  if (_currentLayer < 0 || _currentLayer > _numberOfLayers) 
    return;

  
  // find entry and exit points of track  
  // x = pos[0] + dir[0]*time
  // y = pos[1] + dir[1]*time
  // z = pos[2] + dir[2]*time  

  entry[2] = -_layerHalfThickness[_currentLayer]; 
  exit[2] = _layerHalfThickness[_currentLayer];

  for (int i=0; i<2; ++i) {
    entry[i]=pos[i]+dir[i]*(entry[2]-pos[2])/dir[2];
    exit[i]=pos[i]+dir[i]*(exit[2]-pos[2])/dir[2];
  }

  for (int i=0; i<3; ++i) {
    _currentLocalPosition[i] = pos[i];
    _currentEntryPoint[i] = entry[i];
    _currentExitPoint[i] = exit[i];
  }

 
  double tanx = dir[0]/dir[2];
  double tany = dir[1]/dir[2];  
  double trackLength = min(1.0e+3,_layerThickness[_currentLayer]*sqrt(1.0+tanx*tanx+tany*tany));
  //std::cout << "HERE WE ARE " << trackLength << " " << _segmentLength << " "<<_layerThickness[_currentLayer]<<" "<<tanx<<" "<<tany<<std::endl;

  double dEmean = 1e-6*_energyLoss * trackLength;  
//   _ampl = _fluctuate->SampleFluctuations(double(1000.*_currentParticleMomentum),
//                                               double(1000.*_currentParticleMass),
//                                               _cutOnDeltaRays,segmentLength,
//                                               double(1000.*dEmean))/1000.;


  //  dEmean = hit->getdEdx()/((double)_numberOfSegments);
  _numberOfSegments = int(trackLength/_segmentLength) + 1;
  //std::cout << "number of segments = " << _numberOfSegments << std::endl;
  dEmean = dEmean/((double)_numberOfSegments);
  _ionisationPoints.resize(_numberOfSegments);

  _eSum = 0.0;

  double segmentLength = trackLength/((double)_numberOfSegments);
  _segmentDepth = _layerThickness[_currentLayer]/((double)_numberOfSegments);

  for (int i=0; i<_numberOfSegments; ++i) {
    double z = -_layerHalfThickness[_currentLayer] + ((double)(i)+0.5)*_segmentDepth;
    double x = pos[0]+dir[0]*(z-pos[2])/dir[2];
    double y = pos[1]+dir[1]*(z-pos[2])/dir[2];
    IonisationPoint ipoint;
    double de = _fluctuate->SampleFluctuations(double(1000.*_currentParticleMomentum),
                                              double(1000.*_currentParticleMass),
                                              _cutOnDeltaRays,segmentLength,
                                              double(1000.*dEmean))/1000.;
     //std::cout << "segment " << i << " dE = " << de << std::endl;
    _eSum = _eSum + de;
    ipoint.eloss = de;
    ipoint.x = x;
    ipoint.y = y;
    ipoint.z = z;
    _ionisationPoints[i] = ipoint;
  }

  //std::cout << "Amplitude = " << _ampl << std::endl;

}


void VTXDigitizer::ProduceSignalPoints() {
  /** Produces signal points on the collection plane.
   */

  double TanLorentzX = 0;
  double TanLorentzY = 0;

  if (_currentModule == 1 || _currentModule == 2) {
    TanLorentzX = _tanLorentzAngle;
  }

  double inverseCosLorentzX = sqrt(1.0+TanLorentzX*TanLorentzX);
  double inverseCosLorentzY = sqrt(1.0+TanLorentzY*TanLorentzY);

  _signalPoints.resize(_numberOfSegments);

  // run over ionisation points
  for (int i=0; i<_numberOfSegments; ++i) {
    IonisationPoint ipoint = _ionisationPoints[i];
    double z = ipoint.z;
    double x = ipoint.x;
    double y = ipoint.y;
    double de = ipoint.eloss;
    double DistanceToPlane = _layerHalfThickness[_currentLayer] - z;
    double xOnPlane = x + TanLorentzX*DistanceToPlane;
    double yOnPlane = y + TanLorentzY*DistanceToPlane;
    double DriftLength = DistanceToPlane*sqrt(1.0+TanLorentzX*TanLorentzX+TanLorentzY*TanLorentzY);
    double SigmaDiff = sqrt(DriftLength/_layerThickness[_currentLayer])*_diffusionCoefficient;
    double SigmaX = SigmaDiff*inverseCosLorentzX;
    double SigmaY = SigmaDiff*inverseCosLorentzY;
    double charge = 1.0e+6*de*_electronsPerKeV;
    SignalPoint  spoint;
    spoint.x = xOnPlane;
    spoint.y = yOnPlane;
    spoint.sigmaX = SigmaX;
    spoint.sigmaY = SigmaY;
    spoint.charge = charge;
    _signalPoints[i] = spoint;
  }


}




void VTXDigitizer::ProduceHits( SimTrackerHitImplVec & vectorOfHits) {
  /** Simulation of fired pixels. Each fired pixel is considered 
   * as SimTrackerHit 
   */

  vectorOfHits.clear();

  _currentTotalCharge = 0.0;

  //cout<<"width "<<_widthOfCluster<<" "<<_numberOfSegments<<endl;

  for (int i=0; i<_numberOfSegments; ++i) {
    SignalPoint spoint = _signalPoints[i];
    double xCentre = spoint.x;
    double yCentre = spoint.y;
    double sigmaX = spoint.sigmaX;
    double sigmaY = spoint.sigmaY;
    double xLo = spoint.x - _widthOfCluster*spoint.sigmaX;
    double xUp = spoint.x + _widthOfCluster*spoint.sigmaX;
    double yLo = spoint.y - _widthOfCluster*spoint.sigmaY;
    double yUp = spoint.y + _widthOfCluster*spoint.sigmaY;
    
    

    _currentTotalCharge += spoint.charge;

    //cout<<"spoint "<<xCentre<<" "<<yCentre<<" "<<sigmaX<<" "<<sigmaY<<" "<<xLo<<" "<<xUp<<" "<<yLo<<" "<<yUp<<endl;
    //    cout<<"charge "<<_currentTotalCharge<<endl;
    int ixLo, ixUp, iyLo, iyUp;

    TransformXYToCellID(xLo,yLo,ixLo,iyLo);
    TransformXYToCellID(xUp,yUp,ixUp,iyUp);

     //   std::cout << i << std::endl;
//        std::cout << xLo << " " << xUp << std::endl;    
//        std::cout << yLo << " " << yUp << std::endl;
//        std::cout << ixLo << " " << ixUp << std::endl;
//        std::cout << iyLo << " " << iyUp << std::endl;

    // Loop over all fired pads 
    // and calculate deposited charges
    for (int ix = ixLo; ix<ixUp+1; ++ix) {
      if (ix >= 0) {
        for (int iy = iyLo; iy<iyUp+1; ++iy) {
          if (iy >=0) {
            double xCurrent,yCurrent;
            TransformCellIDToXY(ix,iy,xCurrent,yCurrent);
            gsl_sf_result result;
            int status = gsl_sf_erf_Q_e(float((xCurrent - 0.5*_pixelSizeX - xCentre)/sigmaX), &result);
            float LowerBound = 1 - result.val;
            status = gsl_sf_erf_Q_e(float((xCurrent + 0.5*_pixelSizeX - xCentre)/sigmaX), &result);
            float UpperBound = 1 - result.val;
            float integralX = UpperBound - LowerBound;
            status = gsl_sf_erf_Q_e(float((yCurrent - 0.5*_pixelSizeY - yCentre)/sigmaY), &result);
            LowerBound = 1 - result.val;
            status = gsl_sf_erf_Q_e(float((yCurrent + 0.5*_pixelSizeY - yCentre)/sigmaY), &result);
            UpperBound = 1 - result.val;
            float integralY = UpperBound - LowerBound;
            float totCharge = float(spoint.charge)*integralX*integralY;
            int iexist = 0;
            int cellID = 100000*ix + iy;
            SimTrackerHitImpl * existingHit = 0;
            for (int iHits=0; iHits<int(vectorOfHits.size()); ++iHits) {
              existingHit = vectorOfHits[iHits];
              int cellid = existingHit->getCellID();
              if (cellid == cellID) {
                iexist = 1;
                break;
              }
            }
            if (iexist == 1) {
              float de = existingHit->getdEdx();
              de += totCharge;
              existingHit->setdEdx( de );
            }
            else {
              SimTrackerHitImpl * hit = new SimTrackerHitImpl();
              double pos[3] = {xCurrent, yCurrent, _layerHalfThickness[_currentLayer]};
              hit->setPosition( pos );
              hit->setCellID( cellID );
              hit->setdEdx( totCharge );
              vectorOfHits.push_back( hit );
            }
          }
        }
      }
    }


  }

}


void VTXDigitizer::TransformXYToCellID(double x, double y, 
                                           int & ix, 
                                           int & iy) {
  /**
   * Function calculates position in pixel matrix based on the 
   * local coordinates of point in the ladder.
   */

  int layer = _currentLayer;
  int nladders = _laddersInLayer[layer];
  double ladderGap = _layerLadderGap[layer];
  double Phi = _currentPhi;
  double ladderLength = _layerLadderLength[layer];

  double yInLadder = 0.0;
  if (y < 0.0) {
    yInLadder = y + ladderLength;
  }
  else {
    yInLadder = y - ladderGap;
  }

  if (yInLadder < 0.0) {
    iy = -1;
    //    std::cout << "warning " << std::endl;
  }
  else {
    iy = int(yInLadder/_pixelSizeY);
  }

  double xInLadder = 0.0;
  if (nladders > 2) { // laddered structure
    xInLadder = x + _layerLadderHalfWidth[layer] + _layerActiveSiOffset[layer];
  }
  else { // cyllindrical structure
    xInLadder = x + (_layerRadius[layer]+_layerHalfThickness[layer])*Phi;
  }

  if (xInLadder < 0.0) {
    ix = -1;
  }
  else {
    ix = int(xInLadder/_pixelSizeX);
  }
  

}

void VTXDigitizer::PositionWithinCell(double x, double y, 
                                          int & ix, int & iy, 
                                          double & xCell, double & yCell) {

  int layer = _currentLayer;
  int nladders = _laddersInLayer[layer];
  double ladderGap = _layerLadderGap[layer];
  double Phi = _currentPhi;

  double yInLadder = 0.0;
  if (y < 0.0) {
    yInLadder = -y - ladderGap;
  }
  else {
    yInLadder = y - ladderGap;
  }

  if (yInLadder < 0.0) {
    iy = -1;
  }
  else {
    iy = int(yInLadder/_pixelSizeY);
  }
  yCell = yInLadder - iy*_pixelSizeY - 0.5*_pixelSizeY;


  double xInLadder = 0.0;
  if (nladders > 2) { // laddered structure
    xInLadder = x + _layerLadderHalfWidth[layer] + _layerActiveSiOffset[layer];
  }
  else { // cyllindrical structure
    xInLadder = x + (_layerRadius[layer]+_layerHalfThickness[layer])*Phi;
  }

  if (xInLadder < 0.0) {
    ix = -1;
  }
  else {
    ix = int(xInLadder/_pixelSizeX);
  }
  xCell = xInLadder - ix*_pixelSizeX - 0.5*_pixelSizeX;


}

void VTXDigitizer::TransformCellIDToXY(int ix, int iy,
                                           double & x, double & y) {
  
  /**
     Function calculates position in the local frame 
     based on the index of pixel in the ladder.
  */

  int layer = _currentLayer;
  int module = _currentModule;
  int nladders = _laddersInLayer[layer];
  double ladderGap = _layerLadderGap[layer];
  double ladderLength = _layerLadderLength[layer];
  double Phi = _currentPhi;

  y = (0.5+double(iy))*_pixelSizeY;

  if (module == 1) 
    y = -ladderLength + y;
  else
    y = ladderGap + y;

  x = (0.5+double(ix))*_pixelSizeX;
  if (nladders > 2) { // laddered structure
    x = x - _layerLadderHalfWidth[layer] - _layerActiveSiOffset[layer];
  }  
  else { // cyllindrical structure
    x = x - (_layerRadius[layer]+_layerHalfThickness[layer])*Phi;
  }
  

}


void VTXDigitizer::PoissonSmearer( SimTrackerHitImplVec & simTrkVec ) {
/**
 * Function that fluctuates charge (in units of electrons)
 * deposited on the fired pixels according to the Poisson
 * distribution...
 */

  for (int ihit=0; ihit<int(simTrkVec.size()); ++ihit) {
    SimTrackerHitImpl * hit = simTrkVec[ihit];
    double charge = hit->getdEdx();
    double rng;
    if (charge > 1000.) { // assume Gaussian 
      double sigma = sqrt(charge);
      rng = double(RandGauss::shoot(charge,sigma));
      hit->setdEdx(rng); 
    }
    else { // assume Poisson
      rng = double(RandPoisson::shoot(charge));
    }
    hit->setdEdx(float(rng));
  }  
}

void VTXDigitizer::GainSmearer( SimTrackerHitImplVec & simTrkVec ) {
/**
 * Simulation of electronic noise.
 */

  int nPixels = int( simTrkVec.size() );
  //  std::cout << "Gain smearer applied" << std::endl;

  for (int i=0;i<nPixels;++i) {
    double Noise = RandGauss::shoot(0.,_electronicNoise);
    SimTrackerHitImpl * hit = simTrkVec[i];
    double charge = hit->getdEdx() + Noise ;
    hit->setdEdx( charge );
  }

}


TrackerHitImpl * VTXDigitizer::ReconstructTrackerHit( SimTrackerHitImplVec & simTrkVec ) {
  /**
   * Emulates reconstruction of Tracker Hit 
   * Tracker hit position is reconstructed as center-of-gravity 
   * of cluster of fired cells. Position is corrected for Lorentz shift.
   */

  // new Tracker Hit
  double pos[3] = {0,0,0};
  double charge = 0;
  int nPixels = 0;
  int ixmin =  1000000;
  int ixmax = -1000000;
  int iymin =  1000000;
  int iymax = -1000000;
  int ixSeed = 0;
  int iySeed = 0;
  _amplMax = 0.0;
//cout<<"size "<<simTrkVec.size()<<endl;
//cout<<"threshold "<<_threshold<<endl;
  for (int iHit=0; iHit<int(simTrkVec.size()); ++iHit) {
    SimTrackerHit * hit = simTrkVec[iHit];
    //cout<<"hit "<<iHit <<" "<<hit->getdEdx()<<endl;
    if (hit->getdEdx() > _threshold) {
      if (nPixels < 100)
        _amplC[nPixels] = hit->getdEdx();
      nPixels++;
      charge += hit->getdEdx();
      int cellID = hit->getCellID();
      int ix = cellID / 100000 ;
      int iy = cellID - 100000 * ix;      
      if (hit->getdEdx() > _amplMax) {
        _amplMax = hit->getdEdx();
        ixSeed = ix;
        iySeed = iy;
      }
      
      if (ix > ixmax)
        ixmax = ix;
      if (ix < ixmin)
        ixmin = ix;
      if (iy > iymax)
        iymax = iy;
      if (iy < iymin)
        iymin = iy;
      for (int j=0; j<2; ++j)
        pos[j] += hit->getdEdx()*hit->getPosition()[j];
    }
  }

  //cout<<"charge "<<charge<<endl;
  if (charge > 0.) {
    for (int j=0; j<2; ++j)
      pos[j] /= charge;
  }

  _nCells = nPixels;
  _ampl = charge;
  _nCoveredX = ixmax - ixmin + 1;
  _nCoveredY = iymax - iymin + 1;

  //cout<<"ampl "<<_ampl<<endl;

  double tanXLorentz = _tanLorentzAngle;
  double tanYLorentz = 0;

  pos[0] = pos[0] - _layerHalfThickness[_currentLayer]*tanXLorentz;
  pos[1] = pos[1] - _layerHalfThickness[_currentLayer]*tanYLorentz;

  //  cout<<"pos "<<pos[0]<<" "<<pos[1]<<endl;


  _clusterWidthX = 0.;
  _clusterWidthY = 0.;
  _ampl33 = 0.;
  _ampl55 = 0.;
  _ampl77 = 0.;
  _ncell33 = 0;
  _ncell55 = 0;
  _ncell77 = 0;

  if (charge > 0. && nPixels > 0) {
    TrackerHitImpl * recoHit = new TrackerHitImpl();
    recoHit->setdEdx( charge );
    for (int iY=0;iY<20;++iY) {
      _amplY[iY] = 0.0;
      _amplX[iY] = 0.0;
    }
    for (int iHit=0; iHit<int(simTrkVec.size()); ++iHit) {
      SimTrackerHit * hit = simTrkVec[iHit];
      if (hit->getdEdx() > _threshold) {
        float deltaX = hit->getPosition()[0]-pos[0];
        _clusterWidthX += deltaX*deltaX*hit->getdEdx();
        deltaX = hit->getPosition()[1]-pos[1];
        _clusterWidthY += deltaX*deltaX*hit->getdEdx();
        int cellID = hit->getCellID();
        int ix = cellID / 100000 ;
        int iy = cellID - 100000 * ix;
        if ((iy - iymin) < 20)
          _amplY[iy-iymin] = _amplY[iy-iymin] + hit->getdEdx();
        if ((ix - ixmin) < 20) 
          _amplX[ix-ixmin] = _amplX[ix-ixmin] + hit->getdEdx();        
        bool frame = abs(ix-ixSeed) < 2;
        frame = frame && abs(iy-iySeed) < 2;
        if (frame) {
          _ncell33++;
          _ampl33 += hit->getdEdx();
        } 
        frame = abs(ix-ixSeed) < 3;
        frame = frame && abs(iy-iySeed) < 3;
        if (frame) {
          _ncell55++;
          _ampl55 += hit->getdEdx();
        }
        frame = abs(ix-ixSeed) < 4;
        frame = frame && abs(iy-iySeed) < 4;
        if (frame) {
          _ncell77++;
          _ampl77 += hit->getdEdx();
        }
      }
    }
    _clusterWidthX = sqrt( _clusterWidthX / charge);
    _clusterWidthY = sqrt( _clusterWidthY / charge);
    double aXCentre = 0;
    double aYCentre = 0;
    for (int i=ixmin+1;i<ixmax;++i) {
      aXCentre += _amplX[i-ixmin];
    }
    for (int i=iymin+1;i<iymax;++i) {
      aYCentre += _amplY[i-iymin];
    }
    aXCentre = aXCentre/max(1,ixmax-ixmin-1);
    aYCentre = aYCentre/max(1,iymax-iymin-1);
    double _xRecoS = 0;
    double _yRecoS = 0;
    double aTot = 0;
    for (int i=ixmin;i<ixmax+1;++i) {
      double xx,yy;
      aTot += _amplX[i-ixmin];      
       TransformCellIDToXY(i,2,xx,yy);      
       if (i != ixmin && i != ixmax) {
         _xRecoS += xx*aXCentre;
       }
       else {
         _xRecoS += xx*_amplX[i-ixmin];
       }
    }
    _xRecoS = _xRecoS / aTot;
    _xRecoS = _xRecoS - _layerHalfThickness[_currentLayer]*tanXLorentz;
    aTot = 0;
    for (int i=iymin;i<iymax+1;++i) {
      double xx,yy;
      TransformCellIDToXY(i,i,xx,yy);      
      aTot += _amplY[i-iymin];
      if (i != iymin && i != iymax) {
        _yRecoS += yy*aYCentre;
      }
      else {
        _yRecoS += yy*_amplY[i-iymin];
      }
    }
    _yRecoS = _yRecoS / aTot;
    _xLocalRecoCOG = pos[0];
    _yLocalRecoCOG = pos[1];
    _xLocalRecoEdge = _xRecoS;
    _yLocalRecoEdge = _yRecoS;
    pos[0] = _xRecoS;
    pos[1] = _yRecoS;
    _xLocalSim = _currentLocalPosition[0];
    _yLocalSim = _currentLocalPosition[1];
    recoHit->setPosition( pos );
    return recoHit;
  }
  else 
    return NULL;
  

}

void VTXDigitizer::TrackerHitToLab( TrackerHitImpl * recoHit) {

  double pos[3];
  for (int i=0; i<3; ++i) 
    pos[i] = recoHit->getPosition()[i];

  double xLab[3];
  TransformToLab( pos, xLab);
  
  recoHit->setPosition( xLab );


} 

void VTXDigitizer::PrintInfo( SimTrackerHit * simHit, TrackerHitImpl * recoHit) {

  std::cout << std::endl;

  std::cout << "Simulated hit position = " 
            << " " << simHit->getPosition()[0]
            << " " << simHit->getPosition()[1]
            << " " << simHit->getPosition()[2] << std::endl;
    
  std::cout << "Reconstructed hit position = " 
            << " " << recoHit->getPosition()[0]
            << " " << recoHit->getPosition()[1]
            << " " << recoHit->getPosition()[2] << std::endl;
    
  
  std::cout << std::endl;
  std::cout << "Type Q to exit display mode" << std::endl;   
  char q = getchar();
  if (q=='q' || q=='Q')
    _debug = 0;
}


void VTXDigitizer::generateBackground(LCCollectionVec * col) {

  for (int ilayer=0;ilayer<_numberOfLayers;++ilayer) {
    double mean = _bkgdHitsInLayer[ilayer];
    int nHits = int(RandPoisson::shoot(mean));
    for (int ihit=0;ihit<nHits;++ihit) {
      double pos[3];
      pos[0] = RandFlat::shoot(_layerLadderWidth[ilayer]);
      pos[1] = _layerLadderGap[ilayer] + RandFlat::shoot(_layerLadderLength[ilayer]);
      pos[2] = 0.0;
      if (RandFlat::shoot(double(1.)) > 0.5) {
        pos[1] = -pos[1];
      }
      pos[0] = pos[0] -  _layerLadderHalfWidth[ilayer] - _layerActiveSiOffset[ilayer];
      _currentLayer = ilayer;
      double xLadders = _laddersInLayer[ilayer];
      double Phi0 = _layerPhiOffset[ilayer];
      double dPhi = 2.0*_layerHalfPhi[ilayer];
      int nPhi = int(RandFlat::shoot(xLadders));
      _currentPhi = - PI2 + double(nPhi)*dPhi + Phi0;
      double xLab[3];
      TransformToLab(pos,xLab);
      TrackerHitImpl * trkHit = new TrackerHitImpl();
      trkHit->setPosition( xLab );
      trkHit->setdEdx(1000.);
      trkHit->setType(ilayer+1);
      col->addElement( trkHit );
    }
    
  }

}
