// Calorimeter digitiser for the IDC ECAL and HCAL 
// For other detectors/models SimpleCaloDigi should be used
#include "ILDCaloDigi.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/LayerLayout.h>
#include <gear/CalorimeterParameters.h>
#include <gear/GearParameters.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace lcio ;
using namespace marlin ;

// protect agains rounding errors 
// will not find caps smaller than this
const float slop = 0.25; // (mm)
const float pi = acos(-1.0);
const float twopi = 2.0*pi;

ILDCaloDigi aILDCaloDigi ;



// helper struct for string comparision
struct XToLower{
  int operator() ( int ch ) {
    return std::tolower ( ch );
  }  
}; 

ILDCaloDigi::ILDCaloDigi() : Processor("ILDCaloDigi") {

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
  
  _outputEcalCollections.push_back(std::string("ECALBarrel"));
  _outputEcalCollections.push_back(std::string("ECALEndcap"));
  _outputEcalCollections.push_back(std::string("ECALOther"));
  _outputHcalCollections.push_back(std::string("HCALBarrel"));
  _outputHcalCollections.push_back(std::string("HCALEndcap"));
  _outputHcalCollections.push_back(std::string("HCALOther"));

  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "ECALOutputCollection0" , 
			    "ECAL Collection of real Hits" , 
			    _outputEcalCollections[0], 
			    std::string("ECALBarrel") ); 
  

  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "ECALOutputCollection1" , 
			    "ECAL Collection of real Hits" , 
			     _outputEcalCollections[1], 
			    std::string("ECALEndcap") ); 
  
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "ECALOutputCollection2" , 
			    "ECAL Collection of real Hits" , 
			    _outputEcalCollections[2], 
			    std::string("ECALOther") ) ; 
  
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "HCALOutputCollection0" , 
			    "HCAL Collection of real Hits" , 
			    _outputHcalCollections[0], 
			    std::string("HCALBarrel")  ); 
  
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "HCALOutputCollection1" , 
			    "HCAL Collection of real Hits" , 
			    _outputHcalCollections[1], 
			    std::string("HCALEndcap") ); 
  
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "HCALOutputCollection2" , 
			    "HCAL Collection of real Hits" , 
			    _outputHcalCollections[2], 
			    std::string("HCALOther") ) ; 
  
 registerOutputCollection( LCIO::LCRELATION, 
			    "RelationOutputCollection" , 
			    "CaloHit Relation Collection" , 
			    _outputRelCollection , 
			    std::string("RelationCaloHit")) ; 
  
  registerProcessorParameter("ECALThreshold" , 
			     "Threshold for ECAL Hits in GeV" ,
			     _thresholdEcal,
			     (float)5.0e-5);



  std::vector<float> hcalThresholds;
  hcalThresholds.push_back(0.00004);
  registerProcessorParameter("HCALThreshold" , 
			     "Threshold for HCAL Hits in GeV" ,
			     _thresholdHcal,
			     hcalThresholds);


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
  

  std::vector<float> calibrHcalBarrel;
  calibrHcalBarrel.push_back(0.);
  registerProcessorParameter("CalibrHCALBarrel" , 
			     "Calibration coefficients for Barrel HCAL" ,
			     _calibrCoeffHcalBarrel,
			     calibrHcalBarrel);


  std::vector<float> calibrHcalEndCap;
  calibrHcalEndCap.push_back(0.);
  registerProcessorParameter("CalibrHCALEndcap", 
			     "Calibration coefficients for EndCap HCAL" ,
			     _calibrCoeffHcalEndCap,
			     calibrHcalEndCap);

  std::vector<float> calibrHcalOther;
  calibrHcalOther.push_back(0.);
  registerProcessorParameter("CalibrHCALOther" , 
			     "Calibration coefficients for Other (Ring) HCAL" ,
			     _calibrCoeffHcalOther,
			     calibrHcalOther);

  registerProcessorParameter("IfDigitalEcal" ,
			     "Digital Ecal" , 
			     _digitalEcal , 
			     0);

  registerProcessorParameter("MapsEcalCorrection" ,
			     "Ecal correction for theta dependency of calibration for MAPS" , 
			     _mapsEcalCorrection , 
			     0);


  registerProcessorParameter("IfDigitalHcal" ,
			     "Digital Hcal" , 
			     _digitalHcal , 
			     0);

  registerProcessorParameter("ECALGapCorrection" , 
			     "Correct for ECAL gaps" ,
			     _ecalGapCorrection,
			     (int)1);

  registerProcessorParameter("ECALEndcapCorrectionFactor" , 
			     "Energy correction for ECAL endcap" ,
			     _ecalEndcapCorrectionFactor,
			     (float)1.025);

  registerProcessorParameter("HCALEndcapCorrectionFactor" , 
			     "Energy correction for HCAL endcap" ,
			     _hcalEndcapCorrectionFactor,
			     (float)1.025);

  registerProcessorParameter("ECALGapCorrectionFactor" , 
			     "Factor applied to gap correction" ,
			     _ecalGapCorrectionFactor,
			     (float)1.0);

  registerProcessorParameter("ECALModuleGapCorrectionFactor" , 
			     "Factor applied to module gap correction" ,
			     _ecalModuleGapCorrectionFactor,
			     (float)0.5);


  registerProcessorParameter("Histograms" , 
			     "Hit times histograms" ,
			     _histograms,
			     (int)0);

  registerProcessorParameter("UseEcalTiming" , 
			     "Use ECAL hit times" ,
			     _useEcalTiming               ,
			     (int)0);

  registerProcessorParameter("ECALCorrectTimesForPropagation" , 
			     "Correct ECAL hit times for propagation: radial distance/c" ,
			     _ecalCorrectTimesForPropagation,
			     (int)0);

  registerProcessorParameter("ECALTimeWindowMin" , 
			     "ECAL Time Window minimum time in ns" ,
			     _ecalTimeWindowMin,
			     (float)-10.0);

  registerProcessorParameter("ECALTimeWindowMax" , 
			     "ECAL Time Window maximum time in ns" ,
			     _ecalTimeWindowMax,
			     (float)+100.0);

  registerProcessorParameter("ECALDeltaTimeHitResolution" , 
			     "ECAL Minimum Delta Time in ns for resolving two hits" ,
			     _ecalDeltaTimeHitResolution,
			     (float)+10.0);

  registerProcessorParameter("ECALTimeResolution" , 
			     "ECAL Time Resolution used to smear hit times" ,
			     _ecalTimeResolution,
			     (float)10.);



  registerProcessorParameter("UseHcalTiming" , 
			     "Use HCAL hit times" ,
			     _useHcalTiming               ,
			     (int)1);

  registerProcessorParameter("HCALCorrectTimesForPropagation" , 
			     "Correct HCAL hit times for propagation: radial distance/c" ,
			     _hcalCorrectTimesForPropagation,
			     (int)0);


  registerProcessorParameter("HCALTimeWindowMin" , 
			     "HCAL Time Window minimum time in ns" ,
			     _hcalTimeWindowMin,
			     (float)-10.0);

  registerProcessorParameter("HCALTimeWindowMax" , 
			     "HCAL Time Window maximum time in ns" ,
			     _hcalTimeWindowMax,
			     (float)+100.0);

  registerProcessorParameter("HCALDeltaTimeHitResolution" , 
			     "HCAL Minimum Delta Time in ns for resolving two hits" ,
			     _hcalDeltaTimeHitResolution,
			     (float)+10.0);

  registerProcessorParameter("HCALTimeResolution" , 
			     "HCAL Time Resolution used to smear hit times" ,
			     _hcalTimeResolution,
			     (float)10.);


}

void ILDCaloDigi::init() {

  _nRun = -1;
  _nEvt = 0;

  if(_histograms){
    fEcal    = new TH1F("fEcal", "Ecal time profile",1000, 0., 1000.0);
    fHcal    = new TH1F("fHcal", "Hcal time profile",1000, 0., 1000.0);
    
    fEcalC    = new TH1F("fEcalC", "Ecal time profile cor",1000, 0., 1000.0);
    fHcalC    = new TH1F("fHcalC", "Hcal time profile cor",1000, 0., 1000.0);
    
    fEcalC1   = new TH1F("fEcalC1", "Ecal time profile cor",100, 0., 1000.0);
    fHcalC1   = new TH1F("fHcalC1", "Hcal time profile cor",100, 0., 1000.0);
    
    fEcalC2   = new TH1F("fEcalC2", "Ecal time profile cor",10, 0., 1000.0);
    fHcalC2   = new TH1F("fHcalC2", "Hcal time profile cor",10, 0., 1000.0);

    fHcalCvsE = new TH2F("fHcalCvsE", "Hcal time profile cor",100, 0., 500.0,100,0.,10.);
  }

  //fg: need to set default encoding in for reading old files...
  CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;

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


void ILDCaloDigi::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
  _nEvt = 0;

} 

void ILDCaloDigi::processEvent( LCEvent * evt ) { 

    
  // create the output collections
  LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);

  // copy the flags from the input collection
  LCFlagImpl flag;
  flag.setBit(LCIO::CHBIT_LONG);

  // 
  // * Reading Collections of ECAL Simulated Hits * 
  // 
  for (unsigned int i(0); i < _ecalCollections.size(); ++i) {

    std::string colName =  _ecalCollections[i] ;
    std::transform( colName.begin() , colName.end() , colName.begin(), XToLower() ) ;
    
    //fg: need to establish the subdetetcor part here 
    //    use collection name as cellID does not seem to have that information
    CHT::Layout caloLayout = CHT::any ;
    if( colName == "barrel" )
      caloLayout = CHT::barrel ;
    else 
      if( colName == "endcap" )
	caloLayout = CHT::endcap ;
      else
	if( colName == "plug" )
	  caloLayout = CHT::plug ;
    
    try{
      LCCollection * col = evt->getCollection( _ecalCollections[i].c_str() ) ;
      string initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);

      CellIDDecoder<SimCalorimeterHit> idDecoder( col );

      // create new collection
      LCCollectionVec *ecalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
      ecalcol->setFlag(flag.getFlag());

      // if making gap corrections clear the vectors holding pointers to calhits
      if(_ecalGapCorrection!=0){
	for(int is=0;is<MAX_STAVES;is++){
	  for(int il=0;il<MAX_LAYERS;il++){	
	    _calHitsByStaveLayer[is][il].clear();
	    _calHitsByStaveLayerModule[is][il].clear();
	  }
	}
      }

      int numElements = col->getNumberOfElements();
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	float energy = hit->getEnergy();


	// apply threshold cut
	if (energy > _thresholdEcal) {
	  int cellid = hit->getCellID0();
	  int cellid1 = hit->getCellID1();
	  int layer = idDecoder(hit)["K-1"];
	  int stave = idDecoder(hit)["S-1"];
	  int module= idDecoder(hit)["M"];
	  // save hits by module/stave/layer if required later
	  float calibr_coeff(1.);
	  float x = hit->getPosition()[0];
	  float y = hit->getPosition()[1];
	  float z = hit->getPosition()[2];
	  float r = sqrt(x*x+y*y+z*z);
	  float cost = fabs(z)/r;
	  if(_digitalEcal){
	    calibr_coeff = this->digitalEcalCalibCoeff(layer);
	    if(_mapsEcalCorrection){ 
	      if(caloLayout == CHT::barrel){
		float correction = 1.1387 - 0.068*cost - 0.191*cost*cost;
		calibr_coeff/=correction;
	      }else{
		float correction = 0.592 + 0.590*cost;
		calibr_coeff/=correction;
	      }
	    }
	  }else{
	    calibr_coeff = this->analogueEcalCalibCoeff(layer);
	  }
	  if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff *= _ecalEndcapCorrectionFactor;

	  //	  float energyCal = energy*calibr_coeff;

	  if(_useEcalTiming){
	    float dt = r/300-0.1;
	    const unsigned int n = hit->getNMCContributions();
	    
	    std::vector<bool> used;
	    for(unsigned int i =0; i<n;i++)used.push_back(fal se);
	    
	    int count = 0;
	    float eCellInTime = 0.;
	    float eCellOutput = 0.;
	    
	    for(unsigned int i =0; i<n;i++){
	      float timei   = hit->getTimeCont(i);
	      float energyi = hit->getEnergyCont(i);
	      float deltat = 0;
	      if(_ecalCorrectTimesForPropagation)deltat=dt;
	      if(timei-deltat > _ecalTimeWindowMin && timei-deltat < _ecalTimeWindowMax){
		float ecor = energyi*calibr_coeff;
		eCellInTime+=ecor;
	      }
	      
	      if(!used[i]){
		// merge with other hits?
		used[i] = true;
		for(unsigned int j =i+1; j<n;j++){
		  if(!used[j]){
		    float timej   = hit->getTimeCont(j);
		    float energyj = hit->getEnergyCont(j);
		    float deltat = fabs(timei-timej);
		    //		    std::cout << " ECAL  deltat : " << deltat << std::endl;
		    if(deltat<_ecalDeltaTimeHitResolution){
		      if(energyj>energyi)timei=timej;
		      energyi+=energyj;
		      used[j] = true;
		    }
		  }
		}
		
		if(_digitalEcal){
		  calibr_coeff = this->digitalEcalCalibCoeff(layer);
		  if(_mapsEcalCorrection){
		    if(caloLayout == CHT::barrel){
		      float correction = 1.1387 - 0.068*cost - 0.191*cost*cost;
		      calibr_coeff/=correction;
		    }else{
		      float correction = 0.592 + 0.590*cost;
		      calibr_coeff/=correction;
		    }
		  }
		}else{
		  calibr_coeff = this->analogueEcalCalibCoeff(layer);
		}
		if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff *= _ecalEndcapCorrectionFactor;

		if(_histograms){
		  fEcal->Fill(timei,energyi*calibr_coeff);
		  fEcalC->Fill(timei-dt,energyi*calibr_coeff);
		  fEcalC1->Fill(timei-dt,energyi*calibr_coeff);
		  fEcalC2->Fill(timei-dt,energyi*calibr_coeff);
		}
		if (energyi > _thresholdEcal) {
		  float timeCor=0;
		  if(_ecalCorrectTimesForPropagation)timeCor=dt;
		  timei = timei - timeCor;
		  if(timei > _ecalTimeWindowMin && timei < _ecalTimeWindowMax){
		    count++;
		    CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
		    if(_ecalGapCorrection!=0){
		      _calHitsByStaveLayer[stave][layer].push_back(calhit);
		      _calHitsByStaveLayerModule[stave][layer].push_back(module);
		    }
		    calhit->setCellID0(cellid);		  
		    calhit->setCellID1(cellid1);
		    if(_digitalEcal){
		      calhit->setEnergy(calibr_coeff);
		    }else{
		      calhit->setEnergy(calibr_coeff*energyi);
		    }
		    eCellOutput+= energyi*calibr_coeff;
		    calhit->setTime(timei);
		    calhit->setPosition(hit->getPosition());
		    calhit->setType( CHT( CHT::had, CHT::hcal , caloLayout ,  layer ) );
		    calhit->setRawHit(hit);
		    ecalcol->addElement(calhit);
		    LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.0);
		    relcol->addElement( rel );
		  }
		}
	      }
	    }
	  }else{
	    CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	    if(_ecalGapCorrection!=0){
	      _calHitsByStaveLayer[stave][layer].push_back(calhit);
	      _calHitsByStaveLayerModule[stave][layer].push_back(module);
	    }
	    float energyi = hit->getEnergy();
	    calhit->setCellID0(cellid);		  
	    calhit->setCellID1(cellid1);
	    if(_digitalEcal){
	      calhit->setEnergy(calibr_coeff);
	    }else{
	      calhit->setEnergy(calibr_coeff*energyi);
	    }
	    calhit->setTime(0);
	    calhit->setPosition(hit->getPosition());
	    calhit->setType( CHT( CHT::had, CHT::hcal , caloLayout ,  layer ) );
	    calhit->setRawHit(hit);
	    ecalcol->addElement(calhit);
	    LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.0);
	    relcol->addElement( rel );
	  }
	  

	  //	  std::cout << hit->getTimeCont(0) << " count = " << count <<  " E ECAL = " << energyCal << " - " << eCellInTime << " - " << eCellOutput << std::endl;
	}
      }
      // if requested apply gap corrections in ECAL ? 
      if(_ecalGapCorrection!=0)this->fillECALGaps();
      // add ECAL collection to event
      ecalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
      evt->addCollection(ecalcol,_outputEcalCollections[i].c_str());
    }
    catch(DataNotAvailableException &e){ 
    }
  }

  
  //
  // * Reading HCAL Collections of Simulated Hits * 
  //

  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {

    std::string colName =  _hcalCollections[i] ;
    std::transform( colName.begin() , colName.end() , colName.begin(), XToLower() ) ;

    //fg: need to establish the subdetetcor part here 
    //    use collection name as cellID does not seem to have that information
    CHT::Layout caloLayout = CHT::any ;
    if( colName.find("barrel")!=std::string::npos )
      caloLayout = CHT::barrel ;
    else 
      if( colName.find("endcap")!=std::string::npos )
	caloLayout = CHT::endcap ;
      else
	if(colName.find("ring")!=std::string::npos )
	  caloLayout = CHT::plug ;
    
    try{
      LCCollection * col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      string initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      int numElements = col->getNumberOfElements();
      CellIDDecoder<SimCalorimeterHit> idDecoder(col);
      LCCollectionVec *hcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
      hcalcol->setFlag(flag.getFlag());
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	float energy = hit->getEnergy();
	//std::cout << " Hit energy " << energy << std::endl;


       	if (energy > _thresholdHcal[0]) {
	  int cellid = hit->getCellID0();
	  int cellid1 = hit->getCellID1();
	  float calibr_coeff(1.);
	  int layer =idDecoder(hit)["K-1"];
	  // NOTE : for a digital HCAL this does not allow for varying layer thickness
	  // with depth - would need a simple mod to make response proportional to layer thickness
	  if(_digitalHcal){
	    calibr_coeff = this->digitalHcalCalibCoeff(caloLayout,energy);
	  }else{
	    calibr_coeff = this->analogueHcalCalibCoeff(caloLayout,layer);
	  }
	  if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff*=_hcalEndcapCorrectionFactor;

	  //float energyCal = energy*calibr_coeff;
	  float x = hit->getPosition()[0];
	  float y = hit->getPosition()[1];
	  float z = hit->getPosition()[2];


	  if(_useHcalTiming){
	    float r = sqrt(x*x+y*y+z*z);
	    float dt = r/300-0.1;
	    const unsigned int n = hit->getNMCContributions();
	    
	    std::vector<bool> used;
	    for(unsigned int i =0; i<n;i++)used.push_back(fal se);
	   
	    int count = 0;
	    float eCellInTime = 0.;
	    float eCellOutput = 0.;

	    for(unsigned int i =0; i<n;i++){
	      float timei   = hit->getTimeCont(i);
	      float energyi = hit->getEnergyCont(i);
	      float deltat = 0;
	      if(_hcalCorrectTimesForPropagation)deltat=dt;
	      if(timei-deltat > _hcalTimeWindowMin && timei-deltat < _hcalTimeWindowMax)eCellInTime+=energyi*calibr_coeff;
	      
	      //std::cout << i << " " << timei << " energy = " <<  energyi*calibr_coeff*1000 << std::endl;
	      if(!used[i]){
		// merge with other hits?
		used[i] = true;
		for(unsigned int j =i+1; j<n;j++){
		  if(!used[j]){
		    float timej   = hit->getTimeCont(j);
		    float energyj = hit->getEnergyCont(j);
		    float deltat = fabs(timei-timej);
		    //		    std::cout << " HCAL  deltat : " << deltat << std::endl;
		    if(deltat<_hcalDeltaTimeHitResolution){
		      if(energyj>energyi)timei=timej;
		      //std::cout << timei << " - " << timej << std::endl;
		      //std::cout << energyi << " - " << energyj << std::endl;
		      energyi+=energyj;
		      used[j] = true;
		      //std::cout << timei << " " << energyi << std::endl;
		    }
		  }
		}
		
		if(_digitalHcal){
		  calibr_coeff = this->digitalHcalCalibCoeff(caloLayout,energyi);
		}else{
		  calibr_coeff = this->analogueHcalCalibCoeff(caloLayout,layer);
		}
		if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff*=_hcalEndcapCorrectionFactor;

		if(_histograms){
		  fHcal->Fill(timei,energyi*calibr_coeff);
		  fHcalC->Fill(timei-dt,energyi*calibr_coeff);
		  fHcalC1->Fill(timei-dt,energyi*calibr_coeff);
		  fHcalC2->Fill(timei-dt,energyi*calibr_coeff);
		  fHcalCvsE->Fill(timei-dt,energyi*calibr_coeff);
		}
		if (energyi > _thresholdHcal[0]) {
		  float timeCor=0;
		  if(_hcalCorrectTimesForPropagation)timeCor=dt;
		  timei = timei - timeCor;
		  if(timei > _hcalTimeWindowMin && timei < _hcalTimeWindowMax){
		    count++;
		    CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
		    calhit->setCellID0(cellid);		  
		    calhit->setCellID1(cellid1);
		    if(_digitalHcal){
		      calhit->setEnergy(calibr_coeff);
		      eCellOutput+= calibr_coeff;
		    }else{
		      calhit->setEnergy(calibr_coeff*energyi);
		      eCellOutput+= energyi*calibr_coeff;
		    }
		    calhit->setTime(timei);
		    calhit->setPosition(hit->getPosition());
		    calhit->setType( CHT( CHT::had, CHT::hcal , caloLayout ,  layer ) );
		    calhit->setRawHit(hit);
		    hcalcol->addElement(calhit);
		    LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.0);
		    relcol->addElement( rel );
		  }
		}
	      }
	    }
	  }else{
	    CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	    calhit->setCellID0(cellid);		  
	    calhit->setCellID1(cellid1);
	    float energyi = hit->getEnergy();
	    if(_digitalHcal){
	      calhit->setEnergy(calibr_coeff);
	    }else{
	      calhit->setEnergy(calibr_coeff*energyi);
	    }
	    //eCellOutput+= energyi*calibr_coeff;
	    calhit->setTime(0);
	    calhit->setPosition(hit->getPosition());
	    calhit->setType( CHT( CHT::had, CHT::hcal , caloLayout ,  layer ) );
	    calhit->setRawHit(hit);
	    hcalcol->addElement(calhit);
	    LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.0);
	    relcol->addElement( rel );
	  }
	  
	  // std::cout << hit->getTimeCont(0) << " count = " << count <<  " EHCAL = " << energyCal << " - " << eCellInTime << " - " << eCellOutput << std::endl;
	}
      }
      // add HCAL collection to event
      hcalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
      evt->addCollection(hcalcol,_outputHcalCollections[i].c_str());
    }
    catch(DataNotAvailableException &e){ 
    }
  }

  // add relation collection for ECAL/HCAL to event
  evt->addCollection(relcol,_outputRelCollection.c_str());

  _nEvt++;

}


void ILDCaloDigi::check( LCEvent * evt ) { }
  
void ILDCaloDigi::end(){ 

  if(_histograms){
    TFile *hfile = new TFile("calTiming.root","recreate");
    fEcal->TH1F::Write();
    fHcal->TH1F::Write();
    fEcalC->TH1F::Write();
    fHcalC->TH1F::Write();
    fEcalC1->TH1F::Write();
    fHcalC1->TH1F::Write();
    fEcalC2->TH1F::Write();
    fHcalC2->TH1F::Write();
    fHcalCvsE->TH2F::Write();
    hfile->Close();
    delete hfile;
  }

} 


void ILDCaloDigi::fillECALGaps( ) { 

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

float ILDCaloDigi::digitalHcalCalibCoeff(CHT::Layout caloLayout, float energy ) {

  float calib_coeff = 0;
  unsigned int ilevel = 0;
  for(unsigned int ithresh=1;ithresh<_thresholdHcal.size();ithresh++){
    // Assume!!!  hit energies are stored as floats, i.e. 1, 2 or 3
    if(energy>_thresholdHcal[ithresh])ilevel=ithresh;   // ilevel = 0 , 1, 2
  }

  switch(caloLayout){
  case CHT::barrel:
    if(ilevel>_calibrCoeffHcalBarrel.size()-1){
      streamlog_out(ERROR)  << " Semi-digital level " << ilevel  << " greater than number of HCAL Calibration Constants (" <<_calibrCoeffHcalBarrel.size() << ")" << std::endl;
    }else{
      calib_coeff = _calibrCoeffHcalBarrel[ilevel];
    }
    break;
  case CHT::endcap:
    if(ilevel>_calibrCoeffHcalEndCap.size()-1){
      streamlog_out(ERROR)  << " Semi-digital level " << ilevel  << " greater than number of HCAL Calibration Constants (" <<_calibrCoeffHcalEndCap.size() << ")" << std::endl;
    }else{
      calib_coeff = _calibrCoeffHcalEndCap[ilevel];
    }
    break;
  case CHT::plug:
    if(ilevel>_calibrCoeffHcalOther.size()-1){
      streamlog_out(ERROR)  << " Semi-digital level " << ilevel  << " greater than number of HCAL Calibration Constants (" <<_calibrCoeffHcalOther.size() << ")" << std::endl;
    }else{
      calib_coeff = _calibrCoeffHcalOther[ilevel];
    }
    break;
  default:
    streamlog_out(ERROR)  << " Unknown HCAL Hit Type " << std::endl;
    break;
  }





  return calib_coeff;
}


float ILDCaloDigi::analogueHcalCalibCoeff(CHT::Layout caloLayout, int layer ) {

  float calib_coeff = 0;

  for (unsigned int k(0); k < _hcalLayers.size(); ++k) {
    int min,max;
    if (k == 0) 
      min = 0;
    else 
      min = _hcalLayers[k-1];

    max = _hcalLayers[k];
    if (layer >= min && layer < max) {
      switch(caloLayout){
      case CHT::barrel:
	calib_coeff = _calibrCoeffHcalBarrel[k];
	break;
      case CHT::endcap:
	calib_coeff = _calibrCoeffHcalEndCap[k];
	break;
      case CHT::plug:
	calib_coeff = _calibrCoeffHcalOther[k];
	break;
      default:
	streamlog_out(ERROR)  << " Unknown HCAL Hit Type " << std::endl;
	break;
      }
    }
  } 
  
  return calib_coeff;
}

float ILDCaloDigi::digitalEcalCalibCoeff(int layer ) {

  float calib_coeff = 0;
  
  for (unsigned int k(0); k < _ecalLayers.size(); ++k) {
    int min,max;
    if (k == 0) 
      min = 0;
    else 
      min = _ecalLayers[k-1];

    max = _ecalLayers[k];
    if (layer >= min && layer < max) {
      calib_coeff = _calibrCoeffEcal[k];
      break;
    }
  } 

  return calib_coeff;

}

float ILDCaloDigi::analogueEcalCalibCoeff(int layer ) {

  float calib_coeff = 0;
  
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
      calib_coeff = _calibrCoeffEcal[k];
      break;
    }
  } 
  return calib_coeff;
  
}
