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

#include "CalorimeterHitType.h"


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

  registerProcessorParameter("CellIDLayerString" ,
			     "name of the part of the cellID that holds the layer" , 
			     _cellIDLayerString , 
			     std::string("K")
			     );

}

void SimpleLCalDigi::init() {

  _nRun = -1;
  _nEvt = 0;


  //fg set to 'default' value of Mokka/Geometry/CGA/src/Encoder32Fcal.cc
  CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("I:10,J:10,K:10,S-1:2") ; 

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

	  calhit->setType( CHT( CHT::em, CHT::lcal, CHT::endcap ,  idDecoder(hit)[ _cellIDLayerString ] ) );

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
