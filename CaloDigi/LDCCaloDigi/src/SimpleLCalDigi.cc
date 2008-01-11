#include "SimpleLCalDigi.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
#include <string>
#include <math.h>

using namespace std;
using namespace lcio ;
using namespace marlin ;


SimpleLCalDigi aSimpleLCalDigi ;


SimpleLCalDigi::SimpleLCalDigi() : Processor("SimpleLCalDigi") {

  _description = "Performs simple digitization of sim lcal hits..." ;
  
  std::vector<std::string> lcalCollections;

  lcalCollections.push_back(std::string("SLcal01_LumiCal"));
  
  registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			    "LCALCollections" , 
			    "LCal Collection Names" ,
			    _lcalCollections ,
			    lcalCollections);
  
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "LCALOutputCollection" , 
			    "LCal Collection of real Hits" , 
			    _outputLCalCollection , 
			    std::string("LCAL")) ; 
  
  registerOutputCollection( LCIO::LCRELATION, 
			    "RelationOutputCollection" , 
			    "CaloHit Relation Collection" , 
			    _outputRelCollection , 
			    std::string("RelationLCalHit")) ; 
  
  registerProcessorParameter("LCalThreshold" , 
			     "Threshold for LCal Hits in GeV" ,
			     _thresholdLCal,
			     (float)0.0);

  registerProcessorParameter("CalibrLCAL" , 
			     "Calibration coefficients for LCAL" ,
			     _calibrCoeffLCal,
			     (float)31.);
}

void SimpleLCalDigi::init() {

  _nRun = -1;
  _nEvt = 0;

  //fg: need to set default encoding in for reading old files...
  CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;

}


void SimpleLCalDigi::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 

void SimpleLCalDigi::processEvent( LCEvent * evt ) { 
    

  LCCollectionVec *lcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);

  LCFlagImpl flag;

  flag.setBit(LCIO::CHBIT_LONG);

  lcalcol->setFlag(flag.getFlag());

  // 
  // * Reading Collections of LCAL Simulated Hits * 
  // 
  string initString;
  for (unsigned int i(0); i < _lcalCollections.size(); ++i) {
    try{
      LCCollection * col = evt->getCollection( _lcalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      int numElements = col->getNumberOfElements();
      CellIDDecoder<SimCalorimeterHit> idDecoder( col );
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	float energy = hit->getEnergy();

	if (energy > _thresholdLCal) {
	  CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	  int cellid = hit->getCellID0();
	  int cellid1 = hit->getCellID1();
	  float calibr_coeff(1.);
	  calibr_coeff = _calibrCoeffLCal;
	  calhit->setCellID0(cellid);
	  calhit->setCellID1(cellid1);
	  calhit->setEnergy(calibr_coeff*energy);
	  float pos[3];

	  pos[0] = hit->getPosition()[0];
	  pos[1] = hit->getPosition()[1];
	  pos[2] = hit->getPosition()[2];
	  calhit->setPosition(pos);
	  calhit->setType((int)0);
	  calhit->setRawHit(hit);
	  lcalcol->addElement(calhit);
	  LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.);
	  relcol->addElement( rel );
	}

      }
    }
    catch(DataNotAvailableException &e){ 
    }
  }
  lcalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
  evt->addCollection(lcalcol,_outputLCalCollection.c_str());


  _nEvt++;

}


void SimpleLCalDigi::check( LCEvent * evt ) { }
  
void SimpleLCalDigi::end(){ } 
