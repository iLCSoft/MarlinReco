/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "CCDDigitizer.h"
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


CCDDigitizer aCCDDigitizer ;


CCDDigitizer::CCDDigitizer() : Processor("CCDDigitizer") {
  
  // modify processor description
  _description = "CCDDigitizer should create VTX TrackerHits from SimTrackerHits" ;
  

  // register steering parameters: name, description, class-variable, default value


  // name of input SimTrackerHit collection
  // (default parameter value : "vxd01_VXD", taken from Mokka)
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


  //  cut on the energy of delta-electrons (in MeV)
  // (default parameter value : 0.03)
  registerProcessorParameter("CutOnDeltaRays",
                             "Cut on delta-ray energy (MeV)",
                             _cutOnDeltaRays,
                             (double)0.030); 


  // pixel size along direction perpendicular to beam axis (in mm)
  registerProcessorParameter("PixelSizeX",
                             "Pixel Size X",
                             _pixelSizeX,
                             (double)0.020);
  

 // pixel size along beam axis (in mm) <br>
  registerProcessorParameter("PixelSizeY",
                             "Pixel Size Y",
                             _pixelSizeY,
                             (double)0.020);//0.025
 


  registerProcessorParameter("Debug",
                             "Debug option",
                             _debug,
                             int(0));


  //  number of electrons produced per MeV of deposited energy
  // (default parameter value : 270.3)
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


  // segment length along track path which is used to subdivide track into segments (in mm).
  // The number of track subsegments is calculated as int(TrackLengthWithinActiveLayer/SegmentLength)+1
  //acoording to first tests  the default value should be set smaller or at least equal to 0.003 
  registerProcessorParameter("SegmentLength",
                             "Segment Length",
                             _segmentLength,
                             double(0.003));


 // flag to switch on smearing of signal (signal noise)
  registerProcessorParameter("PoissonSmearing",
                             "Apply Poisson smearing of electrons collected on pixels",
                             _PoissonSmearing,
                             1);

  // flag to switch on gaussian smearing of signal (electronic noise)
  registerProcessorParameter("ElectronicEffects",
                             "Apply Electronic Effects",
                             _electronicEffects,
                             int(1));


  registerProcessorParameter("ElectronicNoise",
                             "electronic noise in electrons",
                             _electronicNoise,
                             double(100));

 registerProcessorParameter("Saturation",
                             "maximum number of electrons, which can be stroed in a pixel",
                             _saturation,
                             double(10000));


 //  registerProcessorParameter("AdditionalCollections",
//                              "Additional Collections to store hit position in the local ladder frame",
//                              _produceFullPattern,
//                              int(0));

  //mean energy loss of MIP per tracklength
  //(default parameter value 280)
  registerProcessorParameter("EnergyLoss",
                             "Energy Loss keV/mm",
                             _energyLoss,
                             double(280.0));


 // flag to switch on additional background signals
  registerProcessorParameter("GenerateBackground",
                             "Generate Background",
                             _generateBackground,
                             int(0));


 registerProcessorParameter("depletedDepth",
                             "Thickness of depleted zone",
                             depdep,
                             double(0.01));


 registerProcessorParameter("undepletedDepth",
                             "Thickness of undepleted zone",
                            undep,
                            double(0.01)); // 0.03744 layerthickness



 registerProcessorParameter("BField",
                             "BField in Tesla",
                             _bfield,
                             (double)4);

 //Voltage applied on the n layer is used to calculate electric field in depleted zone,this parameter is used only when parameter Electric Field is set to be negative
 registerProcessorParameter("BiasVolt",
                            "BiasVoltage in V",
                            _biasvolt,
                            (double)10);

 //electric field in depleted zone 
 //(default parameter value: 10000)
 registerProcessorParameter("ElectricField",
                            "Electric Field in V/cm",
                            _efield,
                            (double)10000);//when using input biasvoltage, set efield negative


 registerProcessorParameter("mu",
                             "electrons mobility in depleted zoneat low elcectric fields, in  cm^2*V^-1*S^-1",
                             _mu,
                            (double)1400);//value from sinevs code:1340;the paper mentioned in void settanlorentzangle suggested 1417


 registerProcessorParameter("Temperatur",
                             "Temperatur in Kelvin",
                             _T,
                            (double)250);//220-270 according to Konstantin


//default parameter value from sinevs code: 34; value in "Vertex Detectors:The State of the Art and Future Prospect",by Damerell,Dec 95: 34.6 cm^2/s at room temperature
 registerProcessorParameter("diffusioncoefficient",
                             "diffcoeff in depleted zone in cm^2*S^-1",
                             _difcoef,
                            (double)34);


// flag to choose reconstruction method;
// method 0: centre of gravity finder for all pixels above threshold
// method 1: centre of gravity finder for all pixels within certain distance to the pixel with highest amplitude 
  registerProcessorParameter("reconstructmethod",
                             "reconstruction method",
                             _recmethod,
                             (int)1);//

 //threshold on charge deposited on one pixel (in electons),only needed reconstruction method 0
  registerProcessorParameter("Threshold",
                             "Cell Threshold in electrons",
                             _threshold,
                             200.);


  //maximum distance between pixel within the considered frame and the pixel with highest amplitude
  //only needed reconstruction method 1
  registerProcessorParameter("framesize",
                             "size of cluster in reconstruction",
                             _framesize,
                             (int)2);
}






void CCDDigitizer::init() { 

  // usually a good idea to
  printParameters();
  // internal parameters
  PI = (double)acos((double)(-1.0));
  TWOPI = (double)(2.0)*PI;	
  PI2 = 0.5*PI;
  _nRun = 0 ;
  _nEvt = 0 ;
  _fluctuate = new MyG4UniversalFluctuationForSi();

  sigmacoefficient=0.6;//coefficient for diffusion in undepleted layer; Sinevs default parameter value 0.6  
  epitaxdep=depdep+undep;//thickness of epitaxel layer 
  midpixx=(maxpixx-1)/2;
  midpixy=(maxpixy-1)/2;
  stepx= _pixelSizeX / Numstepx;//distance between points at which amplitude of diffusion is calculated within one pixel in x direction
  stepy= _pixelSizeY / Numstepy;//distance between points  at which amplitude of diffusion is calculated within one pixel in y direction  
  xobsoffset= stepx/2;
  yobsoffset= stepy/2;//offset for positions, so that positions are arranged symmetrically in the pixel
  maxnionpoint=10000;//maximum number of ionisationpoionts per track

  if(_efield<0) _efield=_biasvolt/(depdep*0.1);
  settanlorentzangle(_bfield,_efield,_mu,_T);
  // settanlorentzangleb(_bfield,_efield,_mu,_T);//lorentzangle according to the paper Christian found

 

#ifdef CCD_diagnostics
  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
 
  histdist=  pHistogramFactory->createHistogram1D( "dist","distance (mm) between rawhits and reconstructed hits",100, 0, 0.05 );  
   histcluster= pHistogramFactory->createHistogram1D("clustersize","number of fired pixels in cluster",20, 0, 20 );  
  histclustxy= pHistogramFactory->createHistogram2D("clustxy","width and heigth of cluster",30,0,15 , 30,0,15 );
   histcharge= pHistogramFactory->createHistogram1D( "clustcharge","number of electrons per hit",120,0,12000);
  histdistxy= pHistogramFactory->createHistogram2D("distxy","distance(mm) between rawhits and reconstructed hit (x,y- direction in local coordinates of ladder)",100,-0.05,0.05 , 100,-0.05,0.05 );
  histNionpoint= pHistogramFactory->createHistogram1D( "numberofionpoints","number of ionisationpoints per hit",30,0,30);
  // histzcoord= pHistogramFactory->createHistogram1D( "locraw-z-coordinate","locraw-z-coordinate",100,-0.01,+0.01);

 histenergy= pHistogramFactory->createHistogram1D( "energyloss","energy loss per hit",100,0,10);

 histsignal= pHistogramFactory->createHistogram1D( "signal","electronnumber per pixel",100,0,500);
 histsignalframe= pHistogramFactory->createHistogram1D( "signalframe","electronnumber per pixel",100,0,1000);

 histenergycentre= pHistogramFactory->createHistogram1D( "signalcentre","electronnumber per pixel",100,0,5000);
#endif


 //table
 // sigmin=0;//minimum sigma for diffusion,is sigma smaller than sigmin, the total charge is stored in the central pixel
 // sigstep=0.005;
   // settable();
 //table//
  
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

}

void CCDDigitizer::processRunHeader( LCRunHeader*  /*run*/) { 
  _nRun++ ;


} 

void CCDDigitizer::processEvent( LCEvent * evt ) { 

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


      // do digitization for one SimTrackerHit 
      // Produce ionisation points along track
     

      ProduceIonisationPoints( simTrkHit );   
   
      // std::cout << "End of ProduceIonisationPoints( ) " << std::endl;
      // std::cout << "Current layer = " << _currentLayer << std::endl;
      // std::cout << "Total particle momentum = " << totPartMomentum << std::endl;

      bool accept = _currentLayer >= 0 && _currentLayer < int(_layerRadius.size()); 

      if ( accept ) {
       
        SimTrackerHitImplVec simTrkHitVec;
        // Produce fired pixels,ionpoint diffuses to plane 
        ProduceHits( simTrkHitVec );
      
   //      for (int iHit=0; iHit<int(simTrkHitVec.size()); ++iHit) {
//           SimTrackerHit * hit = simTrkHitVec[iHit];
//           //cout<<"hit 1"<<iHit <<" "<<hit->getEDep()<<endl;
//         }
       


        // std::cout << "End of ProduceHits( ) " << std::endl;
        // Apply Poisson Smearing to deposited charges
        if (_PoissonSmearing != 0) PoissonSmearer( simTrkHitVec );
        if (_electronicEffects != 0) GainSmearer( simTrkHitVec );
        //std::cout << "Amplitude 6 = " << _ampl << std::endl;

        TrackerHitImpl * recoHit = ReconstructTrackerHit( simTrkHitVec ); 
              
        // std::cout << "End of ReconstructTrackerHit( ) " << std::endl;
       //  if (_produceFullPattern == 1 && recoHit !=0 ) {
//           SimTrackerHitImpl * sth = new SimTrackerHitImpl();
//           sth->setCellID(simTrkHit->getCellID0());
//           sth->setEDep(simTrkHit->getEDep());
//           double currentPosition[3];
//           for (int j=0; j<3; ++j)
//             currentPosition[j] = 0.5*(_currentExitPoint[j] + _currentEntryPoint[j]);
//           sth->setPosition(currentPosition);
//           TrackerHitImpl * th = new TrackerHitImpl();
//           th->setEDep(recoHit->getEDep());
//           double xp[3];
//           for (int j=0; j<3; ++j)
//             xp[j] = recoHit->getPosition()[j];
//           th->setPosition(xp);          
//           STHLocCol->addElement(sth); 
//           THLocCol->addElement(th);
//           LCRelationImpl * rel = new LCRelationImpl(th,sth,float(1.0));
//           RelLocCol->addElement(rel);
//         }

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


          recoHit->setType(100+simTrkHit->getCellID0());
          THcol->addElement( recoHit );
          //          std::cout << "Hit is added to collection " << _nEvt << std::endl;
          LCRelationImpl * rel = new LCRelationImpl(recoHit,simTrkHit,float(1.0));
          RelCol->addElement(rel);
        }
// Clean Up        
        for (int j=0; j < int(simTrkHitVec.size()); ++j) {
          SimTrackerHit * hit = simTrkHitVec[j];
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



void CCDDigitizer::check( LCEvent *  /*evt*/ ) { 

}


void CCDDigitizer::end(){ 

  delete _fluctuate;
  
  std::cout << "CCDDigitizer::end()  " << name() 
             << " processed " << _nEvt << " events in " << _nRun << " runs "
             << std::endl ;
}

void CCDDigitizer::FindLocalPosition(SimTrackerHit * hit,
                                         double * localPosition,
                                         double * localDirection) { 
  /** Function calculates local coordinates of the sim hit 
   * in the given ladder and local momentum of particle. 
   * Also returns module number and ladder 
   * number.
   * Local coordinate system within the ladder 
   * is defined as following :  <br> 
   *    - x (0) axis lies in the ladder plane and orthogonal to the beam axis <br>
   *    - y (1) axis lies in the ladder plane parallel to beam
   *    - z (2) axis is perpendicular to the ladder plane <br>
  

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


// layer is encoded in CellID; 
   layer = hit->getCellID0() - 1;
                                          
  _currentLayer = layer;

// Check boundary of layers
  if (layer < 0 || layer > _numberOfLayers-1) 
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
  
  if(_currentParticleMomentum==0.){
    for (int i=0; i<3; ++i){
      Momentum[i] = 0.001;
      _currentParticleMomentum += Momentum[i]*Momentum[i];
    }
    _currentParticleMomentum = sqrt(_currentParticleMomentum);

  }

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
  

  for (int ic=0; ic<nLadders; ++ic) {
    PhiLadder = double(ic)*dPhi + Phi0;
    PhiInLocal = PhiInLab - PhiLadder;
    //cout<<"Phi "<<PhiLadder<<" "<<PhiInLocal<<" "<<PhiInLab<<" "<<_layerThickness[layer]<<" "<<Radius<<endl;
    if (RXY*cos(PhiInLocal)-Radius > -_layerThickness[layer] && 
        RXY*cos(PhiInLocal)-Radius < _layerThickness[layer]) {
      break;
    }
    //cout<<"phi ladder "<<PhiLadder<<endl;
  }
  double PhiLocalMom = PhiInLabMom - PhiLadder;
  localPosition[0] = RXY*sin(PhiInLocal);
  localPosition[1] = xLab[2];
  localPosition[2] = RXY*cos(PhiInLocal)-Radius;
  //geant sometimes (about 1 from 10 times) produces a z coordinate, which is not in the exact centre of the active layer, but which is still in the active layer. 

  localDirection[0]=PXY*sin(PhiLocalMom);
  localDirection[1]=Momentum[2];
  localDirection[2]=PXY*cos(PhiLocalMom);
  _currentPhi = PhiLadder;

//   if(_debug){
//      cout<<"x: "<<localPosition[0]<<"  ";
//     cout<<"y: "<<localPosition[1]<<"  ";
//     cout<<"z: "<<localPosition[2]<<"  "<<endl;

#ifdef CCD_diagnostics
  for(int i=0;i<3;i++){
    dirraw[i]=localDirection[i];
    posraw[i]=localPosition[i];
  }
#endif
  
}  


 


void CCDDigitizer::TransformToLab(double * xLoc, double * xLab) {
  /** Function transforms local coordinates in the ladder
   * into global coordinates
   */

  int layer = _currentLayer;
  double Phi = _currentPhi;
  double Radius = _layerRadius[layer];

 
    double baseLine = Radius + xLoc[2];
    double PhiInLab = Phi + atan2(xLoc[0],baseLine);
    double RXY = sqrt(baseLine*baseLine+xLoc[0]*xLoc[0]);
    xLab[2] = xLoc[1];
    xLab[0] = RXY*cos(PhiInLab);
    xLab[1] = RXY*sin(PhiInLab);  

}

void CCDDigitizer::ProduceIonisationPoints( SimTrackerHit * hit) {
  /** Produces ionisation points along track segment within active Silicon layer.
   */
  //std::cout << "Amplitude 1 = " << _ampl << std::endl;
  // Position of hit in the Lab frame
  double pos[3];
  double dir[3];
  // double dir[3];


  // Find local position in the ladder frame
  FindLocalPosition( hit, pos, dir);
  //  std::cout << "End of FindLocalPosition" << std::endl;

  // check layer boundaries
  if (_currentLayer < 0 || _currentLayer > _numberOfLayers) 
    return;

  

 
  double tanx = dir[0]/dir[2];
  double tany = dir[1]/dir[2];  
  double trackLength = min(10.,epitaxdep*sqrt(1.0+tanx*tanx+tany*tany));
  // double trackLength = min(1.0e+3,_layerThickness[_currentLayer]*sqrt(1.0+tanx*tanx+tany*tany));//code without bulk

  double dEmean = 1e-6*_energyLoss * trackLength;  


  //  dEmean = hit->getEDep()/((double)_numberOfSegments);


  _numberOfSegments = int(trackLength/_segmentLength) + 1;
 
  dEmean = dEmean/((double)_numberOfSegments);
// if(_numberOfSegments>1000) cout <<"check1000,ionpointnumber: "<< _numberOfSegments<<endl;
// if(_numberOfSegments>500) cout <<"check500,ionpointnumber: "<< _numberOfSegments<<endl;
// if(_numberOfSegments>100) cout <<"check100,ionpointnumber: "<< _numberOfSegments<<endl;
  if(_numberOfSegments>maxnionpoint){   
    _numberOfSegments=maxnionpoint;
     dEmean = hit->getEDep()/((double)_numberOfSegments)*epitaxdep/_layerThickness[_currentLayer];
    }


#ifdef CCD_diagnostics
  histNionpoint->fill(_numberOfSegments,1);
  Nionpoint=_numberOfSegments;
#endif



  _ionisationPoints.resize(_numberOfSegments);

  double segmentLength = trackLength/((double)_numberOfSegments);
  _segmentDepth =epitaxdep/((double)_numberOfSegments);
  //_segmentDepth =_layerThickness[_currentLayer]/((double)_numberOfSegments);//code without bulk

#ifdef CCD_diagnostics
  double energy=0;
#endif



  for (int i=0; i<_numberOfSegments; ++i) {
    double z = (_layerHalfThickness[_currentLayer]-epitaxdep) + ((double)(i)+0.5)*_segmentDepth;
    // double z =  - _layerHalfThickness[_currentLayer] + ((double)(i)+0.5)*_segmentDepth;//code without bulk
    double x = pos[0]+dir[0]*(z-pos[2])/dir[2];
    double y = pos[1]+dir[1]*(z-pos[2])/dir[2];
    IonisationPoint ipoint;
    double de = _fluctuate->SampleFluctuations(double(1000.*_currentParticleMomentum),
                                              double(1000.*_currentParticleMass),
                                              _cutOnDeltaRays,segmentLength,
                                              double(1000.*dEmean))/1000.;
    // if(_debug) std::cout << "segment " << i << " dE = " << de << std::endl;
#ifdef CCD_diagnostics
    energy+=1e+6*de;
#endif

    ipoint.eloss = de;
    ipoint.x = x;
    ipoint.y = y;
    ipoint.z = z;
    _ionisationPoints[i] = ipoint;
  }


#ifdef CCD_diagnostics
  energy=energy/(trackLength)*epitaxdep;
  histenergy->fill(energy,1);
#endif

  //std::cout << "Amplitude = " << _ampl << std::endl;

}


void CCDDigitizer::ProduceHits( SimTrackerHitImplVec & vectorOfHits) {
  // Produces signal points on the collection plane.
 
  double TanLorentzX = TanLorentzAngle;
  double TanLorentzY = 0;
 
  //loop over all ionisisationpoints of the hit
  for(int i=0; i<_numberOfSegments; i++) {
    IonisationPoint ipoint = _ionisationPoints[i];
    double z = ipoint.z;
    double x = ipoint.x;
    double y = ipoint.y;
    double energy=ipoint.eloss;
    // if(_debug)cout<<"x: "<<x<<" y: "<<y<< "z: "<<z<<"  energy  "<< energy<<endl;
  
   

    //distance to plane,where charge is accumulated
    double distancetoplane;
    distancetoplane = _layerHalfThickness[_currentLayer] - z;
      
    
    double xdif, ydif;//distance between position of hit and edge of pixel
    int xcell, ycell;//number of rows and columns of pixel,in which the ionpoint is located

    //check in which area ionpoint is located
   
  //  //  //if the ionpoint is in the bulk, no charge is deposited in detector 
//      if(distancetoplane>epitaxdep){
//        continue;
//      }


    //in undepleted zone: reflected and direct part contributes to charge in pixel, both parts are weighted by a factor depending on distance to the undepleted zone
    // the diffusion in the depleted zone is neglected for that part
    if(distancetoplane<=epitaxdep && distancetoplane>depdep){

      x+= TanLorentzX*depdep;//magnetic effects in undepleted zone are neglected
      y+= TanLorentzY*depdep;
     
      double sigmadirect= sigmacoefficient*(distancetoplane-depdep);
      double sigmareflect=sigmacoefficient*(2*undep-(distancetoplane-depdep));
      double weight=(distancetoplane-depdep)/undep;      
      //assuming a linear relation between distance to plane and probability that charge reaches reflecting plane first
      //reflectionprocesses of higher order are neglected
      
      
     //  double ** spxl;
//       spxl= new double* [maxpixx];
//       for(int i=0;i<maxpixx;i++) spxl[i]=new double [maxpixy];
      double spxl[maxpixx][maxpixy];
      //array for summation of direct and reflected parts
      
      TransformXYToCellID(x, y, xcell, ycell, xdif, ydif);
      
       diffusion(xdif, ydif,sigmadirect);
       //diffusiontable(xdif, ydif,sigmadirect);
      
      for(int j=0;j<maxpixx;j++){
        for(int k=0;k<maxpixy;k++){
          
          spxl[j][k]= (1-weight) * pxl[j][k];
        }
      }
      
      diffusion(xdif, ydif,sigmareflect);
      //diffusiontable(xdif, ydif,sigmareflect);
      for(int j=0;j<maxpixx;j++){
        for(int k=0;k<maxpixy;k++){
          pxl[j][k]= spxl[j][k]+ weight * pxl[j][k]; 
        }
      }      
      //  delete spxl;
    }
    
    //in undepleted zone the width of the diffusion depends on the traveltime of the charge
    if(distancetoplane<=depdep){
      
      x+= TanLorentzX*distancetoplane;
      y+= TanLorentzY*distancetoplane;
      
      double trt=(distancetoplane*0.1)/(_mu*_efield);// trt in s,convert distance into cm
      
      // double sigma=1.375*sqrt(_difcoef*trt)*10;// convert in mm with a factor 10, according to gear data;Sinevs formular
      double sigma=sqrt(2*_difcoef*trt)*10;// convert in mm with a factor 10, according to gear data;formular taken from "Vertex Detectors:The State of the Art and Future Prospects",by Damerell,Dec 95
            
      TransformXYToCellID(x, y, xcell, ycell, xdif, ydif);
      
      diffusion(xdif, ydif,sigma); 
      //diffusiontable(xdif, ydif,sigma);
    }
    
    
    //the values of the charge distribution are stored in pxl 
    // now add the values of pxl to the list of fired pixel
   
    double ladderwidth =  _layerLadderWidth[_currentLayer];
    double Numladderpixx=(int)(ladderwidth/_pixelSizeX);
    double ladderlength = _layerLadderLength[_currentLayer];
    double Numladderpixy=(int)(ladderlength/_pixelSizeY); 
     

    for (int j = 0; j<maxpixx; j++) {
      int ix=j+xcell-midpixx;
      // if(_debug)cout<<ix<<","<<endl;

      if (ix >= 0 && ix<=Numladderpixx) {//test, whether pixel exists

        for (int k = 0; k<maxpixy; k++) {
          int iy=k+ycell-midpixy;
          if (iy >=0 && iy<=Numladderpixy) {


            double charge=(1e+6*energy*_electronsPerKeV) *pxl[j][k]; 
            // if(_debug)cout<<"charge   " <<i<<" "<<k<<" "<<charge<<endl;
            int iexist = 0;

            //  double xCurrent,yCurrent;
            //  TransformCellIDToXY(ix,iy,xCurrent,yCurrent);

            int currentcellid = 100000*ix + iy;
            SimTrackerHitImpl * existingHit = 0;

            for (int iHits=0; iHits<int(vectorOfHits.size()); ++iHits) {
              existingHit = vectorOfHits[iHits];
              int cellid = existingHit->getCellID0();
              if (cellid == currentcellid) {
                iexist = 1;
                break;
              }
            }
            if (iexist == 1) {
              float edep = existingHit->getEDep();
              edep += charge;
              existingHit->setEDep( edep );
            }
            else {
              SimTrackerHitImpl * hit = new SimTrackerHitImpl();

              // double pos[3] = {xCurrent, yCurrent, 0};
              // hit->setPosition( pos );

              hit->setCellID0( currentcellid );
              hit->setEDep( charge );
              vectorOfHits.push_back( hit );
            }
          }
        }
      }
    }
  }
}


void CCDDigitizer::TransformXYToCellID(double x, double y, 
                                           int & ix, 
                                           int & iy,double & xdif, double & ydif) {
  /**
   * Function calculates position in pixel matrix based on the 
   * local coordinates of point in the ladder. Also calculates
   * the position within this pixel.
   */

  int layer = _currentLayer;
  double ladderGap = _layerLadderGap[layer];
  double ladderLength = _layerLadderLength[layer];

  double yInLadder = 0.0;
  if (y < 0.0) {
    yInLadder = y + ladderLength;
  }
  else {
    yInLadder = y - ladderGap;
  }


  iy = int(yInLadder/_pixelSizeY);//position in pixelmatrix
  ydif= yInLadder-(((double)iy)*_pixelSizeY); // position within the pixel
  
  double xInLadder = 0.0;
  
  xInLadder = x + _layerLadderHalfWidth[layer] + _layerActiveSiOffset[layer];
  
  ix = int(xInLadder/_pixelSizeX);
  xdif= xInLadder-(((double)ix)*_pixelSizeX);  

}

void CCDDigitizer::TransformCellIDToXY(int ix, int iy,
                                           double & x, double & y) {
  
  /**
     Function calculates position in the local frame 
     based on the index of pixel in the ladder.
  */

  int layer = _currentLayer;
  int module = _currentModule;
  double ladderGap = _layerLadderGap[layer];
  double ladderLength = _layerLadderLength[layer];

  y = (0.5+double(iy))*_pixelSizeY;

  if (module == 1) 
    y = -ladderLength + y;
  else
    y = ladderGap + y;

  x = (0.5+double(ix))*_pixelSizeX;
 
  x = x - _layerLadderHalfWidth[layer] - _layerActiveSiOffset[layer];  

}


void CCDDigitizer::PoissonSmearer( SimTrackerHitImplVec & simTrkVec ) {
/**
 * Function that fluctuates charge (in units of electrons)
 * deposited on the fired pixels according to the Poisson
 * distribution...
 */

  for (int ihit=0; ihit<int(simTrkVec.size()); ++ihit) {
    SimTrackerHitImpl * hit = simTrkVec[ihit];
    double charge = hit->getEDep();
    double rng;
    if (charge > 1000.) { // assume Gaussian 
      double sigma = sqrt(charge);
      rng = double(RandGauss::shoot(charge,sigma));
      hit->setEDep(rng); 
    }
    else { // assume Poisson
      rng = double(RandPoisson::shoot(charge));
    }
    hit->setEDep(float(rng));
  }  
}

void CCDDigitizer::GainSmearer( SimTrackerHitImplVec & simTrkVec ) {
/**
 * Simulation of electronic noise.
 */

  int nPixels = int( simTrkVec.size() );
  //  std::cout << "Gain smearer applied" << std::endl;

  for (int i=0;i<nPixels;++i) {
    double Noise = RandGauss::shoot(0.,_electronicNoise);
    SimTrackerHitImpl * hit = simTrkVec[i];
    double charge = hit->getEDep() + Noise ;
    if (charge> _saturation) charge = _saturation;
    hit->setEDep( charge );
  }

}


TrackerHitImpl * CCDDigitizer::ReconstructTrackerHit( SimTrackerHitImplVec & simTrkVec ) {
   /**
   * Emulates reconstruction of Tracker Hit 
   * Position is corrected for Lorentz shift.
   */
 
  // new Tracker Hit
  double pos[3] = {0,0,0};
  double charge = 0;
  

  //generic centre of gravity finder for all pixels above threshold
  if(_recmethod==0){
    for (int iHit=0; iHit<int(simTrkVec.size()); ++iHit) {
      SimTrackerHit * hit = simTrkVec[iHit];
      
      //  if(_debug) cout<<"hit "<<iHit <<" "<<hit->getEDep()<<endl;
      // if(_debug) cout<< "edep:  "<< hit->getEDep()<<endl;
      //istsignal->fill(hit->getEDep(),1);
      
      if (hit->getEDep() > _threshold) {
        charge+= hit->getEDep();
        
        
        int cellID = hit->getCellID0();
        int ix = cellID / 100000 ;
        int iy = cellID - 100000 * ix;                 
        double xCurrent,yCurrent;
        TransformCellIDToXY(ix,iy,xCurrent,yCurrent);
        pos[0]+=xCurrent* hit->getEDep();
        pos[1]+=yCurrent* hit->getEDep();
        
        // for (int j=0; j<2; ++j){
        // pos[j] += hit->getEDep()*hit->getPosition()[j];
      }
      
    }
  }  
  
  //looks for the pixel with highest amplitude and computes centre of gravity within a certain grid around this pixel
  if(_recmethod==1){
    
    int xcentre=-10000;
    int ycentre=-10000;
    double emax=0;
    
    for (int iHit=0; iHit<int(simTrkVec.size()); ++iHit) {
      SimTrackerHit * hit = simTrkVec[iHit];
                     
      if (hit->getEDep() > emax) {
        emax=hit->getEDep();        
        int cellID = hit->getCellID0();
        xcentre = cellID / 100000 ;
        ycentre = cellID - 100000 * xcentre;
      }
    }
#ifdef CCD_diagnostics
    histenergycentre->fill(emax,1);//number of electrons in pixel with highest amplitude
#endif
    
    for (int iHit=0; iHit<int(simTrkVec.size()); ++iHit) {
      SimTrackerHit * hit = simTrkVec[iHit];
      
      int cellID = hit->getCellID0();
      int ix = cellID / 100000 ;
      int iy = cellID - 100000 * ix;
      
      bool inframe = (abs(ix-xcentre) < _framesize) && (abs(iy-ycentre) < _framesize);
      if (inframe) {
                       
        charge+= hit->getEDep();
        double xCurrent,yCurrent;
        TransformCellIDToXY(ix,iy,xCurrent,yCurrent);
        pos[0]+=xCurrent* hit->getEDep();
        pos[1]+=yCurrent* hit->getEDep();

#ifdef CCD_diagnostics
        double a=hit->getEDep();
        histsignalframe->fill(a,1);
#endif
      }
      
    }
    
  }
  
  
  if(charge>0){
    TrackerHitImpl * recoHit = new TrackerHitImpl();
    recoHit->setEDep( charge );
    for (int j=0; j<2; ++j)
      pos[j]/=charge;
    
    
    
    
    //correction due to magnetic effects (shift of depdep weighted by undep,shift of depdep/2 is weighted by depdep)
    double tanXLorentz = TanLorentzAngle;
    double meanlorentzdepth= ((undep*depdep)+depdep*(depdep/2))/(undep+depdep);
    
    pos[0] = pos[0] - meanlorentzdepth * tanXLorentz;


    //  cout<<"xuncor: "<<pos[0]<<"  ";
//         cout<<"yuncor: "<<pos[1]<<"  ";
//         cout<<"zuncor: "<<pos[2]<<"  "<<endl;
  
      
    // correction of z coordinate,because of the bulk:
    double bulkthickness= _layerThickness[_currentLayer] - epitaxdep;
    pos[2] += bulkthickness/2;
    
                
#ifdef CCD_diagnostics
    //shift along the track back to direction initial interaction plane to make results comparable      
    double shift[3];
    shift[2]=-bulkthickness/2;

    for (int i=0; i<2; ++i) {      
      if (dirraw[2]!=0)
      shift[i]=(dirraw[i]/dirraw[2])*shift[2];
    }
    
    for (int i=0;i<3;i++){
      pos[i]+=shift[i];
    }
    
    //  cout<<"xcor: "<<pos[0]<<"  ";
//       cout<<"ycor: "<<pos[1]<<"  ";
//       cout<<"zcor: "<<pos[2]<<"  "<<endl;
    
    
    int ixmin =  1000000;
    int ixmax = -1000000;
    int iymin =  1000000;
    int iymax = -1000000;
    int nPixels=0;     
   
       
    for (int iHit=0; iHit<int(simTrkVec.size()); ++iHit) {
      SimTrackerHit * hit = simTrkVec[iHit];  
      if (hit->getEDep() > _threshold) {
        nPixels++;
  
        int cellID = hit->getCellID0();
        int ix = cellID / 100000 ;
        int iy = cellID - 100000 * ix     ;
        
        if (ix > ixmax)
          ixmax = ix;
        if (ix < ixmin)
          ixmin = ix;
        if (iy > iymax)
          iymax = iy;
        if (iy < iymin)
          iymin = iy;
               
      }
    }

    int xwidth=ixmax-ixmin+1;
    int ywidth=iymax-iymin+1;
    
   //   cout<<"xmin "<< ixmin<<" xmax "<<ixmax<<endl;
//      cout <<"xwidth: "<<xwidth<<" "<<"ywidth "<<ywidth<<endl;
    
    
    
//       cout<<"npixels:" << nPixels<<endl;
//        cout<<"charge: " <<charge<<endl;
    //   cout<<"posraw: " <<posraw[0]<<" ,"<<posraw[1]<<", "<<posraw[2]<<endl
    //       <<"posreco: "<<pos[0]<<" ,"<<pos[1]<<", "<<pos[2]<<endl
    //       <<"numberofionpoints: "<<Nionpoint<<endl;
    


    // cout<<"clustxy:  "<<xwidth<<", "<<ywidth<<endl;
    histclustxy->fill(xwidth,ywidth,1);//maximum distance between two pixels above threshold
    // cout<<"firedpixelnumber:  "<<nPixels<<endl;
    histcluster->fill(nPixels,1);//number of pixels above threshold
    // cout<<"charge:  "<<charge<<endl;
    histcharge->fill(charge,1);//collected electrons per hit
    double dist=sqrt(pow(posraw[0]-pos[0],2)
                     +pow(posraw[1]-pos[1],2)
                     +pow(posraw[2]-pos[2],2));
    //  cout<<"dist:  "<<dist<<endl;
    histdist->fill(dist,1);
    
    double distx=posraw[0]-pos[0];
    double disty=posraw[1]-pos[1];
    //cout<<"histdistxy:  "<<distx<<" ,"<<disty<<endl;
    histdistxy->fill(distx,disty,1);
    //   histzcoord->fill(posraw[2],1);



   //  if(dist>0.1){
//       cout<<" dist>0.1"<<endl;

//       cout<<"dist "<<dist<<endl<<"ionpoints "<<Nionpoint<<endl;
//       cout<<"modul "<<_currentModule<<endl;
//         cout<<"clustxy:  "<<xwidth<<", "<<ywidth<<endl;
//       cout<<"xcor: "<<pos[0]<<"  ";
//       cout<<"ycor: "<<pos[1]<<"  ";
//       cout<<"zcor: "<<pos[2]<<"  "<<endl;
//       cout<<"xuncor: "<<pos[0]<<"  ";
//       cout<<"yuncor: "<<pos[1]<<"  ";
//       cout<<"zuncor: "<<pos[2]<<"  "<<endl; 
//     cout<<"posraw: " <<posraw[0]<<" ,"<<posraw[1]<<", "<<posraw[2]<<endl;
    // char q = getchar();
    // }
    
#endif
    
        
    recoHit->setPosition( pos );
    return recoHit;
    } 
  
  else 
    return NULL;
    
}


void CCDDigitizer::TrackerHitToLab( TrackerHitImpl * recoHit) {

  double pos[3];
  for (int i=0; i<3; ++i) 
    pos[i] = recoHit->getPosition()[i];

  double xLab[3];
  TransformToLab( pos, xLab);
  
  recoHit->setPosition( xLab );


} 

void CCDDigitizer::PrintInfo( SimTrackerHit * simHit, TrackerHitImpl * recoHit) {

  std::cout << std::endl;

  std::cout << "Simulated hit position = " 
            << " " << simHit->getPosition()[0]
            << " " << simHit->getPosition()[1]
            << " " << simHit->getPosition()[2] << std::endl;
    
  std::cout << "Reconstructed hit position = " 
            << " " << recoHit->getPosition()[0]
            << " " << recoHit->getPosition()[1]
            << " " << recoHit->getPosition()[2]
            << " total number of electrons " << recoHit->getEDep() <<std::endl;
    
  
  std::cout << std::endl;
  std::cout << "Type Q to exit display mode" << std::endl;   
  char q = getchar();
  if (q=='q' || q=='Q')
    _debug = 0;
}


void CCDDigitizer::generateBackground(LCCollectionVec * col) {

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
      trkHit->setEDep(1000.);
      trkHit->setType(ilayer+1);
      col->addElement( trkHit );
    }
    
  }

}

void CCDDigitizer::diffusion(double xdif,double ydif, double sigma){

  //function computes the part of charge, which diffuses in each pixel of a maxpixx*maxpixy- array
  //it computes a dampfactor(dependent on the distance between ionpoint and observatiopoint (in the xy plane) and the functionparameter sigma) (sinev) for Numstepx*Numstepy observationpoints in each pixel and summarizes them (the distribution is approximated by discrete points); after that values are normed 
  
  double xobs0,yobs0,xobs,yobs,sum,damp,dist,distsquare;

  double xhit=xdif+(midpixx * _pixelSizeX);
  double yhit=ydif+(midpixy * _pixelSizeY);
  double norm=0;
  double twosigmasquare=2*sigma*sigma;  
  
  for(int xpix=0;xpix<maxpixx;xpix++){
    xobs0=xpix* _pixelSizeX + xobsoffset;
    
    for(int ypix=0;ypix<maxpixy;ypix++){
      yobs0=ypix*_pixelSizeY+ yobsoffset;
      sum=0;
      
      for(int i=0;i<Numstepx;i++){
        xobs=xobs0+i*stepx;
        
        for(int k=0;k<Numstepx;k++){
          yobs=yobs0+k*stepy;
          
          distsquare=(xobs-xhit)*(xobs-xhit)+(yobs-yhit)*(yobs-yhit);
          dist=sqrt(distsquare);
          damp=exp(-distsquare/twosigmasquare)/dist;
          sum+=damp;
        }
        
      } 
      
      pxl[xpix][ypix]=sum;      
      norm+=sum;
      // cout<<"norm: "<<norm<<endl;
      
    }
    
  }
  for (int i=0;i<maxpixx;i++){
    for(int k=0;k<maxpixy;k++){
      pxl[i][k]/=norm;
      // if(_debug)cout <<"pxl "<<i<<k<<" "<<pxl[i][k]<<endl;
    }
    
  }      
  
}


    
void CCDDigitizer::settanlorentzangleb  (double B,double E,double mu,double T){


/* sets lorentzangle for depleted zone

   taken from sinevs code:
   Lorentz angle is calculated according to parametrisation given by
   T.Lari (INFN, ATLAS Pixel Collaboration) in February of 2003 talk
   Strictly speaking this parametrization is valid only for the
   silicon with cariers mobility cited in the above talk
   ( 1920. cm**2V**-1S**-1) 



   B - b field in Tesla
   E - electric field in V/cm 
   mu - initial carriers mobility in cm**2 x V**-1 x S**-1
   T - temperature (in Kelvin)
*/

  double exp2=1.55;
  double exp3=0.66;
  double Ec0 = 6030.;
  double r = 1.13 + 0.0008 * (T - 273.);
  double beta = 1.04 * pow(T/273.0 ,exp3);
 
  double Ec = Ec0 * pow(T/273.0 ,exp2);
  double mud = mu /  pow((1.+pow(E/Ec,beta)),1./beta);
  TanLorentzAngle = 0.0001 * r * mud *B;
 //  cout<<"beta"<<beta<<endl;
//   cout<<"Ec= "<<Ec<<endl;
  cout<<"lorentzangle"<<TanLorentzAngle<<endl;
 
}  

void CCDDigitizer::settanlorentzangle  (double B,double E,double mu,double T){


/* sets lorentzangle for depleted zone

 //according to a paper Kristian Harder found


   B - b field in Tesla
   E - electric field in V/cm 
   mu - initial carriers mobility(at low electric field) in cm**2 x V**-1 x S**-1
   T - temperature (in Kelvin)
*/

  double exp2=0.87;
  double exp3=0.66;//equal to sinevs value
  double r = 1.15;//similar to sinev value
  double beta = 1.109 * pow(T/300.0 ,exp3);//similar to sinevs formular
  // double mulow=1417*pow((T/300),-2.2);
 
  double vsat =1.07e+7 * pow((T/300),exp2);
  double mud = mu /  pow((1.+pow(E*mu/vsat,beta)),1./beta);
  TanLorentzAngle = 0.0001 * r * mud *B;
  // double Ec=vsat/mu;
 //  cout<<"Ec= vsat/mu= "<<Ec<<endl;
//   cout<<"beta"<<beta<<endl;
  cout<<"lorentzangle "<<TanLorentzAngle<<endl;

}  

/*

I started to create a lookuptable to decrease computing time. this part is not debugged. it should reduce computing time to 70% of the time needed without. 
void CCDDigitizer::settable(){

  //creates a lookuptable for diffusionprocess

  int maxx= ((midpixx+2)* Numstepx)+1;
  int maxy= ((midpixy+2)* Numstepy)+1;

  for(int isigma=0; isigma<numsigstep;isigma++){
    
    double norm=0;
    double sigma= sigmin+sigstep*isigma;
    double twosigmasquare=2*sigma*sigma; 
    
    double ** htable;
    htable= new double* [maxx];
    for(int i=0;i<maxx;i++) htable[i]=new double [maxy];

    for(int i=0;i<maxx;i++){
      for(int k=0;k<maxy;k++){

        double distsquare=pow(( i* stepx)-xobsoffset,2)+pow(( k* stepy-yobsoffset),2);
        double dist=sqrt(distsquare);
        htable[i][k]=exp(-distsquare/twosigmasquare)/dist;
      }
    }
  
    for(int i=0; i<numhitstepx ;i++){
      for(int k=0 ;i<numhitstepy ;k++){
        for(int ipixx=0;ipixx<maxpixx;ipixx++){
          for(int ipixy=0;ipixy<maxpixy;ipixy++){          
           
            double sum=0;
           
            cout<<"t " <<i +( (ipixx-midpixx)*Numstepx)<< "y "<<k + (ipixy-midpixy) * Numstepy <<endl;

            for(int x=0;x<Numstepx;x++){
              for(int y=0;y<Numstepy;y++){
                //  cout<<"t " <<i +( (ipixx-midpixx)*Numstepx) +x<< "y "<<k + (ipixy-midpixy) * Numstepy + y<<endl;
                sum+=htable[abs(i +( (ipixx-midpixx)*Numstepx) +x)][abs(k +( (ipixy-midpixy) * Numstepy) + y)]; 
                
              }                  
            }
            table[isigma][i][k][ipixx][ipixy]=sum;
            norm+=sum;        
          }
      
        }
        for (int ipixx=0;ipixx<maxpixx;ipixx++){
          for(int ipixy=0;ipixy<maxpixy;ipixy++){
            table[isigma][i][k][ipixx][ipixy]/=norm;
          }
        }
      }
    }  
    delete htable;
  }
}

void CCDDigitizer::diffusiontable(double xdif, double ydif,double sigma){


  //the diffusionprocess using the lookuptable

  if(sigma<sigmin){
    for (int i=0;i<maxpixx;i++){
      for(int k=0;k<maxpixy;k++){
        pxl[i][k]=0;
      }
    }  
    pxl[midpixx][midpixy]=1;
  }
  else{

    bool mirrowx=0;
    bool mirrowy=0;
    double pix [maxpixx][maxpixy];
    //double ** pix;
    //pix= new double* [maxpixx];
    //for(int i=0;i<maxpixx;i++) pix[i]=new double [maxpixy];
    
    if(xdif> 0.5*_pixelSizeX){
      mirrowx=1;
      xdif=_pixelSizeX-xdif;
    }
    if(ydif> 0.5*_pixelSizeY){
      mirrowy=1;
      ydif=_pixelSizeY-ydif;
    }
    
    int xhit=(int)(xdif/stepx);
    double xweight=(xdif-xhit)/stepx;
    int yhit=(int)(ydif/stepy);
    double yweight=(xdif-xhit)/stepx;

    int isigma=(int)((sigma-sigmin)/sigstep);
    double sigweight= sigma-(isigma*sigstep);
    //   if (sigma<sigmin){
    //     isigma=0;
    //   }
 
    if (sigma>(sigmin+(numsigstep*sigstep))){
      isigma=numsigstep;
      sigweight=0;
    }

    for(int x=0;x<maxpixx;x++){
      for(int y=0;y<maxpixy;y++){    

        double x_y  =table[isigma][xhit]  [yhit]  [x][y];
        double xp_y =table[isigma][xhit+1][yhit]  [x][y];
        double x_yp =table[isigma][xhit]  [yhit+1][x][y];
        double xp_yp=table[isigma][xhit+1][yhit+1][x][y];
        double xsum_y =(1-xweight) *x_y  + xweight * xp_y;
        double xsum_yp=(1-xweight) *x_yp + xweight * xp_yp;
        double xsum_ysum_d=(1-yweight)*xsum_y+ yweight * xsum_yp;
        
        x_y   =table[isigma+1][xhit]  [yhit]  [x][y];
        xp_y  =table[isigma+1][xhit+1][yhit]  [x][y];
        x_yp  =table[isigma+1][xhit]  [yhit+1][x][y];
        xp_yp =table[isigma+1][xhit+1][yhit+1][x][y];
        xsum_y =(1-xweight)*x_y+ xweight * xp_y;
        xsum_yp=(1-xweight)*x_yp+ xweight * xp_yp;
        double xsum_ysum_u=(1-yweight)*xsum_y+ yweight * xsum_yp;
        
        pix[x][y]= (1-sigweight)*xsum_ysum_d+ sigweight * xsum_ysum_u;        
      }
    }
        
    for(int x=0;x<maxpixx;x++){
      double ix=x;
      if(mirrowx) ix=maxpixx-1-x;
      
      for(int y=0;x<maxpixy;y++){
        double iy=y;
        if(mirrowy) iy=maxpixy-1-y;        
        pxl[x][y]=pix[x][y];
      }
    } 
        // delete pix;                        
  }
}
*/
