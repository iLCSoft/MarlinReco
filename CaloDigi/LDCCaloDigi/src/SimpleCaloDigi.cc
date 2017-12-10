#include "SimpleCaloDigi.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
// #include <string>
// #include <algorithm>

#include "CalorimeterHitType.h"


using namespace std;
using namespace lcio ;
using namespace marlin ;


SimpleCaloDigi aSimpleCaloDigi ;



SimpleCaloDigi::SimpleCaloDigi() : Processor("SimpleCaloDigi") {

  _description = "Performs simple digitization of sim calo hits..." ;
  
  std::vector<std::string> ecalCollections;

  ecalCollections.push_back(std::string("ecal02_EcalBarrel"));
  ecalCollections.push_back(std::string("ecal02_EcalEndcap"));

  registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			    "ECALCollections" , 
			    "ECAL Collection Names" ,
			    _ecalCollections ,
			    ecalCollections);
  
  std::vector<std::string> hcalCollections;

  hcalCollections.push_back(std::string("hcalFeScintillator_HcalBarrelEnd"));
  hcalCollections.push_back(std::string("hcalFeScintillator_HcalBarrelReg"));
  hcalCollections.push_back(std::string("hcalFeScintillator_HcalEndCaps"));

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
			     (float)1.0e-4);

  registerProcessorParameter("HCALThreshold" , 
			     "Threshold for HCAL Hits in GeV" ,
			     _thresholdHcal,
			     (float)4.0e-4);


  std::vector<int> ecalLayers;
  ecalLayers.push_back(30);
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
  calibrEcal.push_back(31.3);
  calibrEcal.push_back(83.0);


  registerProcessorParameter("CalibrECAL" , 
			     "Calibration coefficients for ECAL" ,
			     _calibrCoeffEcal,
			     calibrEcal);
  

  std::vector<float> calibrHcal;
  calibrHcal.push_back(27.3);

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

  //fg: need to set default encoding in for reading old files...
  CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;

}


void SimpleCaloDigi::processRunHeader( LCRunHeader*  /*run*/) { 
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
  string initString;
  for (unsigned int i(0); i < _ecalCollections.size(); ++i) {

    std::string colName =  _ecalCollections[i] ;
    
    //fg: need to establish the subdetetcor part here 
    //    use collection name as cellID does not seem to have that information
    CHT::Layout caloLayout = layoutFromString( colName ) ;

    try{
      LCCollection * col = evt->getCollection( _ecalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      int numElements = col->getNumberOfElements();
      CellIDDecoder<SimCalorimeterHit> idDecoder( col );
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	float energy = hit->getEnergy();
	
	if (energy > _thresholdEcal) {
	  CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	  int cellid = hit->getCellID0();
	  int cellid1 = hit->getCellID1();
	  float calibr_coeff(1.);
	  int layer = idDecoder(hit)["K-1"];
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
	  calhit->setCellID1(cellid1);
	  if (_digitalEcal) {
	    calhit->setEnergy(calibr_coeff); 
	  }
	  else {
	    calhit->setEnergy(calibr_coeff*energy);
	  }
	  calhit->setPosition(hit->getPosition());

	  calhit->setType( CHT( CHT::em, CHT::ecal, caloLayout ,  layer ) );
	  
	  calhit->setRawHit(hit);
	  ecalcol->addElement(calhit);
	  LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.);
	  relcol->addElement( rel );
	}

      }
    }
    catch(DataNotAvailableException &e){ 
    }
  }
  ecalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
  evt->addCollection(ecalcol,_outputEcalCollection.c_str());


  //
  // * Reading HCAL Collections of Simulated Hits * 
  //

  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {

    std::string colName =  _hcalCollections[i] ;

    //fg: need to establish the subdetetcor part here 
    //    use collection name as cellID does not seem to have that information
    CHT::Layout caloLayout = layoutFromString( colName ) ; 


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

	  calhit->setType( CHT( CHT::had, CHT::hcal , caloLayout ,  layer ) );

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
  hcalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
  evt->addCollection(hcalcol,_outputHcalCollection.c_str());
  evt->addCollection(relcol,_outputRelCollection.c_str());

  _nEvt++;

}


void SimpleCaloDigi::check( LCEvent *  /*evt*/ ) { }
  
void SimpleCaloDigi::end(){ } 
