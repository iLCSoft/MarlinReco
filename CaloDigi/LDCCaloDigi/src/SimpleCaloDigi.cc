#include "SimpleCaloDigi.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <iostream>

#define SHIFT_M 0
#define SHIFT_S 3
#define SHIFT_I 6
#define SHIFT_J 15
#define SHIFT_K 24
#define SHIFT_2 30
#define SHIFT_1 31

#define MASK_M (unsigned int) 0x00000007
#define MASK_S (unsigned int) 0x00000038
#define MASK_I (unsigned int) 0x00007FC0
#define MASK_J (unsigned int) 0x00FF8000
#define MASK_K (unsigned int) 0x3F000000
#define MASK_2 (unsigned int) 0x40000000
#define MASK_1 (unsigned int) 0x80000000


using namespace lcio ;
using namespace marlin ;


SimpleCaloDigi aSimpleCaloDigi ;


SimpleCaloDigi::SimpleCaloDigi() : Processor("SimpleCaloDigi") {

  _description = "Performs simple digitization of sim calo hits..." ;
  
  std::vector<std::string> ecalCollections;

  ecalCollections.push_back(std::string("ecal02_EcalBarrel"));
  ecalCollections.push_back(std::string("ecal02_EcalEndcap"));

  registerProcessorParameter( "ECALCollections" , 
			      "ECAL Collection Names" ,
			      _ecalCollections ,
			       ecalCollections);

  std::vector<std::string> hcalCollections;

  hcalCollections.push_back(std::string("hcalFeRPC1_HcalBarrelEnd"));
  hcalCollections.push_back(std::string("hcalFeRPC1_HcalBarrelReg"));
  hcalCollections.push_back(std::string("hcalFeRPC1_HcalEndCaps"));

  registerProcessorParameter("HCALCollections" , 
			     "HCAL Collection Names" , 
			     _hcalCollections , 
			     hcalCollections);

  registerProcessorParameter("ECALOutputCollection" , 
			     "ECAL Collection of real Hits" , 
			     _outputEcalCollection , 
			     std::string("ECAL")) ; 

  registerProcessorParameter("HCALOutputCollection" , 
			     "HCAL Collection of real Hits" , 
			     _outputHcalCollection , 
			     std::string("HCAL")) ; 

  
  registerProcessorParameter("RelationOutputCollection" , 
			     "CaloHit Relation Collection" , 
			     _outputRelCollection , 
			     std::string("RelationCaloHit")) ; 

  registerProcessorParameter("ECALThreshold" , 
			     "Threshold for ECAL Hits in GeV" ,
			     _thresholdEcal,
			     (float)1.0e-4);

  registerProcessorParameter("HCALThreshold" , 
			     "Threshold for HCAL Hits in GeV" ,
			     _thresholdHcal,
			     (float)0.25e-6);

  std::vector<int> ecalLayers;
  ecalLayers.push_back(30);
  ecalLayers.push_back(40);


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
  calibrEcal.push_back(1.);
  calibrEcal.push_back(1.);


  registerProcessorParameter("CalibrECAL" , 
			     "Calibration coefficients for ECAL" ,
			     _calibrCoeffEcal,
			     calibrEcal);
  

  std::vector<float> calibrHcal;
  calibrHcal.push_back(1.);

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


}

void SimpleCaloDigi::init() {

    _nRun = -1;
    _nEvt = 0;

}


void SimpleCaloDigi::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 

void SimpleCaloDigi::processEvent( LCEvent * evt ) { 
    

  LCCollectionVec *ecalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  LCCollectionVec *hcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);

  LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);

  LCFlagImpl flag;

  flag.setBit(LCIO::CHBIT_LONG);

  ecalcol->setFlag(flag.getFlag());
  hcalcol->setFlag(flag.getFlag());


// 
// * Reading Collections of ECAL Simulated Hits * 
// 

  for (unsigned int i(0); i < _ecalCollections.size(); ++i) {
      try{
	  LCCollection * col = evt->getCollection( _ecalCollections[i].c_str() ) ;
	  int numElements = col->getNumberOfElements();
	  for (int j(0); j < numElements; ++j) {
	      SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	      float energy = hit->getEnergy();

	      if (energy > _thresholdEcal) {
		  CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
		  int cellid = hit->getCellID0();
		  float calibr_coeff(1.);
		  int layer = cellid >> 24;
		  for (unsigned int k(0); k < _ecalLayers.size(); ++k) {
		      int min,max;
		      if (k == 0) 
			  min = 0;		      
		      else 
			  min = _ecalLayers[k-1];		      
		      max = _ecalLayers[k];
		      if (layer >= min && layer < max) {
			  calibr_coeff = _calibrCoeffEcal[k];
			  break;
		      }
		  } 
		  calhit->setCellID0(cellid);
		  if (_digitalEcal) {
		      calhit->setEnergy(calibr_coeff); 
		  }
		  else {
		      calhit->setEnergy(calibr_coeff*energy);
		  }
		  calhit->setPosition(hit->getPosition());
		  calhit->setType((int)0);
		  ecalcol->addElement(calhit);
		  LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.);
		  relcol->addElement( rel );
	      }

	  }
      }
      catch(DataNotAvailableException &e){ 
      }
  }

  evt->addCollection(ecalcol,_outputEcalCollection.c_str());


//
// * Reading HCAL Collections of Simulated Hits * 
//

  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
      try{
	  LCCollection * col = evt->getCollection( _hcalCollections[i].c_str() ) ;
	  int numElements = col->getNumberOfElements();
	  for (int j(0); j < numElements; ++j) {
	      SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	      float energy = hit->getEnergy();


	      if (energy > _thresholdHcal) {
		  CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
		  int cellid = hit->getCellID0();
		  float calibr_coeff(1.);
		  int layer = cellid >> 24 ;
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
		  if (_digitalHcal) {
		      calhit->setEnergy(calibr_coeff); 
		  }
		  else {
		      calhit->setEnergy(calibr_coeff*energy);
		  }
		  calhit->setPosition(hit->getPosition());
		  calhit->setType(int(1));
		  hcalcol->addElement(calhit);
		  LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.0);
		  relcol->addElement( rel );
	      }

	  }
      }
      catch(DataNotAvailableException &e){ 
      }
  }

  evt->addCollection(hcalcol,_outputHcalCollection.c_str());
  evt->addCollection(relcol,_outputRelCollection.c_str());

  _nEvt++;

}


void SimpleCaloDigi::check( LCEvent * evt ) { }
  
void SimpleCaloDigi::end(){ } 
