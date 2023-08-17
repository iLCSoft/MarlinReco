// Calorimeter digitiser for the LDC ECAL and HCAL 
// For other detectors/models SimpleCaloDigi should be used
#include "LDCCaloDigi.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;
using namespace lcio ;
using namespace marlin ;

// protect agains rounding errors 
// will not find caps smaller than this
const float slop = 0.25; // (mm)
const float pi = acos(-1.0);
const float twopi = 2.0*pi;

LDCCaloDigi aLDCCaloDigi ;


LDCCaloDigi::LDCCaloDigi() : Processor("LDCCaloDigi") {

  _description = "Performs simple digitization of sim calo hits..." ;
  
  std::vector<std::string> ecalCollections;

  ecalCollections.push_back(std::string("EcalBarrelCollection"));
  ecalCollections.push_back(std::string("EcalEndcapCollection"));
  ecalCollections.push_back(std::string("EcalRingCollection"));

  registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			    "ECALCollections" , 
			    "ECAL Collection Names" ,
			    _ecalCollections ,
			    ecalCollections);
  
  std::vector<std::string> hcalCollections;

  hcalCollections.push_back(std::string("HcalBarrelRegCollection"));
  hcalCollections.push_back(std::string("HcalEndcapRingsCollection"));
  hcalCollections.push_back(std::string("HcalEndcapsCollection"));

  registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			    "HCALCollections" , 
			    "HCAL Collection Names" , 
			    _hcalCollections , 
			    hcalCollections);
  
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "ECALOutputCollection" , 
			    "ECAL Collection of real Hits" , 
			    _outputEcalCollection , 
			    std::string("ECAL")) ; 
  
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "HCALOutputCollection" , 
			    "HCAL Collection of real Hits" , 
			    _outputHcalCollection , 
			    std::string("HCAL")) ; 
  
  
  registerOutputCollection( LCIO::LCRELATION, 
			    "RelationOutputCollection" , 
			    "CaloHit Relation Collection" , 
			    _outputRelCollection , 
			    std::string("RelationCaloHit")) ; 
  
  registerProcessorParameter("ECALThreshold" , 
			     "Threshold for ECAL Hits in GeV" ,
			     _thresholdEcal,
			     (float)5.0e-5);

  registerProcessorParameter("HCALThreshold" , 
			     "Threshold for HCAL Hits in GeV" ,
			     _thresholdHcal,
			     (float)2.5e-4);


  std::vector<int> ecalLayers;
  ecalLayers.push_back(20);
  ecalLayers.push_back(100);


  registerProcessorParameter("ECALLayers" , 
			     "Index of ECal Layers" ,
			     _ecalLayers,
			     ecalLayers);

  

  std::vector<int> hcalLayers;
  hcalLayers.push_back(100);

  registerProcessorParameter("HCALLayers" , 
			     "Index of HCal Layers" ,
			     _hcalLayers,
			     hcalLayers);


  std::vector<float> calibrEcal;
  calibrEcal.push_back(40.91);
  calibrEcal.push_back(81.81);


  registerProcessorParameter("CalibrECAL" , 
			     "Calibration coefficients for ECAL" ,
			     _calibrCoeffEcal,
			     calibrEcal);
  

  std::vector<float> calibrHcal;
  calibrHcal.push_back(34.8);

  registerProcessorParameter("CalibrHCAL" , 
			     "Calibration coefficients for HCAL" ,
			     _calibrCoeffHcal,
			     calibrHcal);


  registerProcessorParameter("IfDigitalEcal" ,
			     "Digital Ecal" , 
			     _digitalEcal , 
			     0);


  registerProcessorParameter("IfDigitalHcal" ,
			     "Digital Hcal" , 
			     _digitalHcal , 
			     0);

  registerProcessorParameter("ECALGapCorrection" , 
			     "Correct for ECAL gaps" ,
			     _ecalGapCorrection,
			     (int)1);

  registerProcessorParameter("ECAlEndcapCorrectionFactor" , 
			     "Energy correction for endcap" ,
			     _ecalEndcapCorrectionFactor,
			     (float)1.025);

  registerProcessorParameter("ECALGapCorrectionFactor" , 
			     "Factor applied to gap correction" ,
			     _ecalGapCorrectionFactor,
			     (float)1.0);

  registerProcessorParameter("ECALModuleGapCorrectionFactor" , 
			     "Factor applied to module gap correction" ,
			     _ecalModuleGapCorrectionFactor,
			     (float)0.5);


}

void LDCCaloDigi::init() {

  _nRun = -1;
  _nEvt = 0;

  //fg: need to set default encoding in for reading old files...
  // dudarboh: does nothing anymore, default encoding spedified in the constructor
  // CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;

}


void LDCCaloDigi::processRunHeader( LCRunHeader*  /*run*/) { 

  _nRun++ ;
  _nEvt = 0;

  // Calorimeter geometry from GEAR
  const gear::CalorimeterParameters& pEcalBarrel = Global::GEAR->getEcalBarrelParameters();
  const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
  //  const gear::CalorimeterParameters& pHcalBarrel = Global::GEAR->getHcalBarrelParameters();
  //  const gear::CalorimeterParameters& pHcalEndcap = Global::GEAR->getHcalEndcapParameters();
  const gear::LayerLayout& ecalBarrelLayout = pEcalBarrel.getLayerLayout();
  const gear::LayerLayout& ecalEndcapLayout = pEcalEndcap.getLayerLayout();
  // const gear::LayerLayout& hcalBarrelLayout = pHcalBarrel.getLayerLayout();
  // const gear::LayerLayout& hcalEndcapLayout = pHcalEndcap.getLayerLayout();

  // determine geometry of ECAL
  int symmetry = pEcalBarrel.getSymmetryOrder();
  _zOfEcalEndcap = (float)pEcalEndcap.getExtent()[2];

  // Determine ECAL polygon angles
  // Store radial vectors perpendicular to stave layers in _ecalBarrelStaveDir 
  // ASSUMES Mokka Stave numbering 0 = top, then numbering increases anti-clockwise
  if(symmetry>1){
    float nFoldSymmetry = static_cast<float>(symmetry);
    float phi0 = pEcalBarrel.getPhi0();
    for(int i=0;i<symmetry;++i){
      float phi  = phi0 + i*twopi/nFoldSymmetry;
      _barrelStaveDir[i][0] = cos(phi);
      _barrelStaveDir[i][1] = sin(phi);
    }
  }  

  for(int i=0;i<ecalBarrelLayout.getNLayers();++i){
    _barrelPixelSizeT[i] = ecalBarrelLayout.getCellSize0(i);
    _barrelPixelSizeZ[i] = ecalBarrelLayout.getCellSize1(i);
   }

  for(int i=0;i<ecalEndcapLayout.getNLayers();++i){
    _endcapPixelSizeX[i] = ecalEndcapLayout.getCellSize0(i);
    _endcapPixelSizeY[i] = ecalEndcapLayout.getCellSize1(i);
  }

} 

void LDCCaloDigi::processEvent( LCEvent * evt ) { 
    
  // create the output collections
  LCCollectionVec *ecalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  LCCollectionVec *hcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);

  // copy the flags from the input collection
  LCFlagImpl flag;
  flag.setBit(LCIO::CHBIT_LONG);
  ecalcol->setFlag(flag.getFlag());
  hcalcol->setFlag(flag.getFlag());

  // if making gap corrections clear the vectors holding pointers to calhits
  if(_ecalGapCorrection!=0){
    for(int is=0;is<MAX_STAVES;is++){
      for(int il=0;il<MAX_LAYERS;il++){	
	_calHitsByStaveLayer[is][il].clear();
	_calHitsByStaveLayerModule[is][il].clear();
      }
    }
  }

  // 
  // * Reading Collections of ECAL Simulated Hits * 
  // 
  string initString;
  for (unsigned int i(0); i < _ecalCollections.size(); ++i) {
    try{
      LCCollection * col = evt->getCollection( _ecalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      int numElements = col->getNumberOfElements();
      CellIDDecoder<SimCalorimeterHit> idDecoder( col );
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	float energy = hit->getEnergy();
	// apply threshold cut
	if (energy > _thresholdEcal) {
	  CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	  int cellid = hit->getCellID0();
	  int cellid1 = hit->getCellID1();
	  float calibr_coeff(1.);
	  int layer = idDecoder(hit)["K-1"];
	  int stave = idDecoder(hit)["S-1"];
	  int module= idDecoder(hit)["M"];
	  // save hits by module/stave/layer if required later
	  if(_ecalGapCorrection!=0){
	    _calHitsByStaveLayer[stave][layer].push_back(calhit);
	    _calHitsByStaveLayerModule[stave][layer].push_back(module);
	  }

	  // retrieve calibration constants
	  for (unsigned int k(0); k < _ecalLayers.size(); ++k) {
	    int min,max;
	    if (k == 0){ 
	      min = 0;		      
	    }else{ 
	      min = _ecalLayers[k-1];
	    } 
	    max = _ecalLayers[k];
	    if (layer >= min && layer < max) {
	      calibr_coeff = _calibrCoeffEcal[k];
	      break;
	    }
	  } 
	  // apply calibration
	  if (_digitalEcal) {
	    calhit->setEnergy(calibr_coeff); 
	  }
	  else {
	    // if in endcap apply additional factor to calibration to account for
	    // the difference in response due to the orientation of B wrt absorber
	    if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)energy=energy*_ecalEndcapCorrectionFactor;
	    calhit->setEnergy(calibr_coeff*energy);
	  }
	  // set other ECAL quanties
	  calhit->setPosition(hit->getPosition());
	  calhit->setType((int)0);
	  calhit->setRawHit(hit);
	  calhit->setCellID0(cellid);
	  calhit->setCellID1(cellid1);
	  ecalcol->addElement(calhit);
	  // make relation between hit and sim hit
	  LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.);
	  relcol->addElement( rel );
	}

      }
    }
    catch(DataNotAvailableException &e){ 

    }
  }

  // if requested apply gap corrections in ECAL ? 
  if(_ecalGapCorrection!=0)this->fillECALGaps();

  // add ECAL collection to event
  ecalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
  evt->addCollection(ecalcol,_outputEcalCollection.c_str());
  
  //
  // * Reading HCAL Collections of Simulated Hits * 
  //

  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    try{
      LCCollection * col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      int numElements = col->getNumberOfElements();
      CellIDDecoder<SimCalorimeterHit> idDecoder(col);
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	float energy = hit->getEnergy();


	if (energy > _thresholdHcal) {

	  CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	  int cellid = hit->getCellID0();
	  int cellid1 = hit->getCellID1();
	  float calibr_coeff(1.);
	  int layer =idDecoder(hit)["K-1"]; 
	  for (unsigned int k(0); k < _hcalLayers.size(); ++k) {
	    int min,max;
	    if (k == 0) 
	      min = 0;
	    else 
	      min = _hcalLayers[k-1];
	    max = _hcalLayers[k];
	    if (layer >= min && layer < max) {
	      calibr_coeff = _calibrCoeffHcal[k];
	      break;
	    }
	  } 
	  calhit->setCellID0(cellid);		  
	  calhit->setCellID1(cellid1);
	  if (_digitalHcal) {
	    calhit->setEnergy(calibr_coeff); 
	  }
	  else {
	    calhit->setEnergy(calibr_coeff*energy);
	  }
	  calhit->setPosition(hit->getPosition());
	  calhit->setType(int(1));
	  calhit->setRawHit(hit);
	  hcalcol->addElement(calhit);
	  LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.0);
	  relcol->addElement( rel );
	}

      }
    }
    catch(DataNotAvailableException &e){ 
    }
  }

  // add HCAL collection to event
  hcalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
  evt->addCollection(hcalcol,_outputHcalCollection.c_str());

  // add relation collection for ECAL/HCAL to event
  evt->addCollection(relcol,_outputRelCollection.c_str());

  _nEvt++;

}


void LDCCaloDigi::check( LCEvent *  /*evt*/ ) { }
  
void LDCCaloDigi::end(){ } 


void LDCCaloDigi::fillECALGaps( ) { 

  // Loop over hits in the Barrel
  // For each layer calculated differences in hit positions
  // Look for gaps based on expected separation of adjacent hits
  // loop over staves and layers

  for (int is=0; is < MAX_STAVES; ++is) {
    for (int il=0; il < MAX_LAYERS; ++il) {
      if(_calHitsByStaveLayer[is][il].size()>1){
	// compare all pairs of hits just once (j>i)

	for (unsigned int i=0;i<_calHitsByStaveLayer[is][il].size()-1;++i){
	  CalorimeterHitImpl* hiti = _calHitsByStaveLayer[is][il][i]; 
	  int modulei = _calHitsByStaveLayerModule[is][il][i];
	  float xi = hiti->getPosition()[0];
	  float yi = hiti->getPosition()[1];
	  float zi = hiti->getPosition()[2];

	  for (unsigned int j=i+1;j<_calHitsByStaveLayer[is][il].size();++j){
	    CalorimeterHitImpl* hitj = _calHitsByStaveLayer[is][il][j]; 
	    int modulej = _calHitsByStaveLayerModule[is][il][j];
	    float xj = hitj->getPosition()[0];
	    float yj = hitj->getPosition()[1];
	    float zj = hitj->getPosition()[2];
	    float dz = fabs(zi-zj);
	    // *** BARREL CORRECTION ***
	    if( fabs(zi)<_zOfEcalEndcap && fabs(zj)<_zOfEcalEndcap){
	      // account for stave directions using normals
	    // calculate difference in hit postions in z and along stave
	      float dx = xi-xj;
	      float dy = yi-yj;
	      float dt = fabs(dx*_barrelStaveDir[is][0] + dy*_barrelStaveDir[is][1]);
	      // flags for evidence for gaps
	      bool zgap = false;   // in z direction
	      bool tgap = false;   // along stave 
	      bool ztgap = false;  // in both z and along stave 
	      bool mgap = false;   // gaps between ECAL modules
	      
	      // criteria gaps in the z and t direction
	      float zminm = 1.0*_barrelPixelSizeZ[il]-slop;
	      float zmin = 1.0*_barrelPixelSizeZ[il]+slop;
	      float zmax = 2.0*_barrelPixelSizeZ[il]-slop;
	      float tminm = 1.0*_barrelPixelSizeT[il]-slop;
	      float tmin = 1.0*_barrelPixelSizeT[il]+slop;
	      float tmax = 2.0*_barrelPixelSizeT[il]-slop;
	      
	      // criteria for gaps
	      // WOULD BE BETTER TO USE GEAR TO CHECK GAPS ARE OF EXPECTED SIZE
	      if( dz > zmin  && dz < zmax && dt < tminm )zgap = true;
	      if( dz < zminm && dt > tmin && dt < tmax )tgap = true;
	      if( dz > zmin && dz < zmax && dt > tmin && dt < tmax )ztgap=true;

	      if(modulei!=modulej){
		if( dz > zmin && dz < 3.0*_barrelPixelSizeZ[il]-slop && dt < tmin)mgap = true;
	      }
 



	      // found a gap now apply a correction based on area of gap/area of pixel
	      if(zgap||tgap||ztgap||mgap){
		float ecor = 1.;
		float f = _ecalGapCorrectionFactor; // fudge
		if(mgap)f = _ecalModuleGapCorrectionFactor;
		if(zgap||mgap)ecor = 1.+f*(dz - _barrelPixelSizeZ[il])/2./_barrelPixelSizeZ[il];
		if(tgap)ecor = 1.+f*(dt - _barrelPixelSizeT[il])/2./_barrelPixelSizeT[il];
		if(ztgap)ecor= 1.+f*(dt - _barrelPixelSizeT[il])*(dz - _barrelPixelSizeZ[il])/4./_barrelPixelSizeT[il]/_barrelPixelSizeZ[il];     
		float ei = hiti->getEnergy()*ecor;
		float ej = hitj->getEnergy()*ecor;
		hiti->setEnergy(ei);
		hitj->setEnergy(ej);
	      }
	      
	    // *** ENDCAP CORRECTION ***
	    }else if(fabs(zi)>_zOfEcalEndcap && fabs(zj)>_zOfEcalEndcap&&dz<100){
	      float dx = fabs(xi-xj);
	      float dy = fabs(yi-yj);
	      bool xgap = false;
	      bool ygap = false;
	      bool xygap = false;
	      // criteria gaps in the z and t direction
	      float xmin = 1.0*_endcapPixelSizeX[il]+slop;
	      float xminm = 1.0*_endcapPixelSizeX[il]-slop;
	      float xmax = 2.0*_endcapPixelSizeX[il]-slop;
	      float ymin = 1.0*_endcapPixelSizeY[il]+slop;
	      float yminm = 1.0*_endcapPixelSizeY[il]-slop;
	      float ymax = 2.0*_endcapPixelSizeY[il]-slop;
	      // look for gaps
	      if(dx > xmin && dx < xmax && dy < yminm )xgap = true;
	      if(dx < xminm && dy > ymin && dy < ymax )ygap = true;
	      if(dx > xmin && dx < xmax && dy > ymin && dy < ymax )xygap=true;
	    
	      if(xgap||ygap||xygap){
		// found a gap make correction
		float ecor = 1.;
		float f = _ecalGapCorrectionFactor; // fudge
		if(xgap)ecor = 1.+f*(dx - _endcapPixelSizeX[il])/2./_endcapPixelSizeX[il];
		if(ygap)ecor = 1.+f*(dy - _endcapPixelSizeY[il])/2./_endcapPixelSizeY[il];
		if(xygap)ecor= 1.+f*(dx - _endcapPixelSizeX[il])*(dy - _endcapPixelSizeY[il])/4./_endcapPixelSizeX[il]/_endcapPixelSizeY[il];     
		hiti->setEnergy( hiti->getEnergy()*ecor );
		hitj->setEnergy( hitj->getEnergy()*ecor );
	      }
	    }
	  }
	}
      }
    }
  }

  return;

}
