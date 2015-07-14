#include "SimpleLHCalDigi.h"
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


SimpleLHCalDigi aSimpleLHCalDigi ;


SimpleLHCalDigi::SimpleLHCalDigi() : Processor("SimpleLHCalDigi") {

  _description = "Performs simple digitization of sim lhcal hits..." ;
  
  std::vector<std::string> lhcalCollections;

  lhcalCollections.push_back(std::string("LHcalCollection"));
  
  registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			    "LHCALCollections" , 
			    "LHCal Collection Names" ,
			    _lhcalCollections ,
			    lhcalCollections);
    
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "LHCALOutputCollection" , 
			    "LHCal Collection of real Hits" , 
			    _outputLHCalCollection , 
			    std::string("LHCAL")) ; 
  
  registerOutputCollection( LCIO::LCRELATION, 
			    "RelationOutputCollection" , 
			    "CaloHit Relation Collection" , 
			    _outputRelCollection , 
			    std::string("RelationLHCalHit")) ; 
  
  registerProcessorParameter("LHCalThreshold" , 
			     "Threshold for LHCal Hits in GeV" ,
			     _thresholdLHCal,
			     (float)0.0);

  registerProcessorParameter("CalibrLHCAL" , 
			     "Calibration coefficients for LHCAL" ,
			     _calibrCoeffLHCal,
			     (float)31.);

  registerProcessorParameter("CellIDLayerString" ,
			     "name of the part of the cellID that holds the layer" , 
			     _cellIDLayerString , 
			     std::string("K-1")
			     );

  registerProcessorParameter("CaloType" ,
			     "type of calorimeter: em = 0, had = 1, muon = 2" , 
			     _caloType , 
			     int(0)
			     );
//   registerProcessorParameter("CaloType" ,
// 			     "type of calorimeter: em, had, muon" , 
// 			     _caloType , 
// 			     std::string("had")
// 			     );

  registerProcessorParameter("CaloID" ,
			     "ID of calorimeter: lcal, lhcal, bcal", 
			     _caloID , 
			     std::string("lhcal")
			     );

  registerProcessorParameter("CaloLayout" ,
			     "subdetector layout: barrel, endcap, plug, ring", 
			     _caloLayout , 
			     std::string("endcap")
			     );

  registerProcessorParameter("DefaultEncoding" ,
			     "string defining cell encoding" , 
			     _defaultEncoding , 
			     std::string("M:3,S-1:3,I:9,J:9,K-1:6")
			     );

}

void SimpleLHCalDigi::init() {

  _nRun = -1;
  _nEvt = 0;

  //fg: need to set default encoding in for reading old files...
  //CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;
  CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding(_defaultEncoding.c_str()) ;

}


void SimpleLHCalDigi::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 

void SimpleLHCalDigi::processEvent( LCEvent * evt ) { 
    

  LCCollectionVec *lcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);

  LCFlagImpl flag;

  flag.setBit(LCIO::CHBIT_LONG);

  lcalcol->setFlag(flag.getFlag());

  // 
  // * Reading Collections of LHCAL Simulated Hits * 
  // 
  string initString;
  for (unsigned int i(0); i < _lhcalCollections.size(); ++i) {
    try{
      LCCollection * col = evt->getCollection( _lhcalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      int numElements = col->getNumberOfElements();
      CellIDDecoder<SimCalorimeterHit> idDecoder( col );
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	float energy = hit->getEnergy();

	if (energy > _thresholdLHCal) {
	  CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	  int cellid = hit->getCellID0();
	  int cellid1 = hit->getCellID1();
	  float calibr_coeff(1.);
	  calibr_coeff = _calibrCoeffLHCal;
	  calhit->setCellID0(cellid);
	  calhit->setCellID1(cellid1);
	  calhit->setEnergy(calibr_coeff*energy);
	  float pos[3];

	  pos[0] = hit->getPosition()[0];
	  pos[1] = hit->getPosition()[1];
	  pos[2] = hit->getPosition()[2];
	  calhit->setPosition(pos);

	  //calhit->setType( CHT( CHT::had, _CHType, CHT::endcap ,  idDecoder(hit)[ _cellIDLayerString ] ) );
	  // jl: change to following line as soon as MarlinUtil is updated
          // calhit->setType( caloTypeFromString(_caloType), caloIDFromString(_caloID), layoutFromString(_caloLayout.c_str()),  idDecoder(hit)[ _cellIDLayerString ] ) );

          calhit->setType( CHT( CHT::CaloType(_caloType), caloIDFromString(_caloID), layoutFromString(_caloLayout.c_str()),  idDecoder(hit)[ _cellIDLayerString ] ) );

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
  evt->addCollection(lcalcol,_outputLHCalCollection.c_str());
  evt->addCollection(relcol,_outputRelCollection.c_str());


  _nEvt++;

}


void SimpleLHCalDigi::check( LCEvent * evt ) { }
  
void SimpleLHCalDigi::end(){ } 
