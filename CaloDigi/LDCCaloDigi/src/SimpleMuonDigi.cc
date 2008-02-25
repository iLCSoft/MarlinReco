#include "SimpleMuonDigi.h"
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

using namespace std;
using namespace lcio ;
using namespace marlin ;


SimpleMuonDigi aSimpleMuonDigi ;


SimpleMuonDigi::SimpleMuonDigi() : Processor("SimpleMuonDigi") {

  _description = "Performs simple digitization of sim muon hits..." ;
  
  std::vector<std::string> muonCollections;

  muonCollections.push_back(std::string("yoke03_MuonBarrel"));
  muonCollections.push_back(std::string("yoke03_MuonEndCap"));
  muonCollections.push_back(std::string("yoke03_MuonPlug"));
    

  registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			    "MUONCollections" , 
			    "Muon Collection Names" ,
			    _muonCollections ,
			    muonCollections);
  
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "MUONOutputCollection" , 
			    "Muon Collection of real Hits" , 
			    _outputMuonCollection , 
			    std::string("MUON")) ; 
  
  registerOutputCollection( LCIO::LCRELATION, 
			    "RelationOutputCollection" , 
			    "CaloHit Relation Collection" , 
			    _outputRelCollection , 
			    std::string("RelationMuonHit")) ; 
  
  registerProcessorParameter("MuonThreshold" , 
			     "Threshold for Muon Hits in GeV" ,
			     _thresholdMuon,
			     (float)0.0);

  registerProcessorParameter("CalibrMUON" , 
			     "Calibration coefficients for MUON" ,
			     _calibrCoeffMuon,
			     (float)31.);
}

void SimpleMuonDigi::init() {

  _nRun = -1;
  _nEvt = 0;

  //fg: need to set default encoding in for reading old files...
  CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;

}


void SimpleMuonDigi::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 

void SimpleMuonDigi::processEvent( LCEvent * evt ) { 
    

  LCCollectionVec *muoncol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);

  LCFlagImpl flag;

  flag.setBit(LCIO::CHBIT_LONG);

  muoncol->setFlag(flag.getFlag());

  // 
  // * Reading Collections of MUON Simulated Hits * 
  // 
  string initString;
  for (unsigned int i(0); i < _muonCollections.size(); ++i) {
    try{
      LCCollection * col = evt->getCollection( _muonCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      int numElements = col->getNumberOfElements();
      CellIDDecoder<SimCalorimeterHit> idDecoder( col );
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	float energy = hit->getEnergy();

	if (energy > _thresholdMuon) {
	  CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	  int cellid = hit->getCellID0();
	  int cellid1 = hit->getCellID1();
	  float calibr_coeff(1.);
	  //int layer = idDecoder(hit)["K-1"];
	  calibr_coeff = _calibrCoeffMuon;
	  calhit->setCellID0(cellid);
	  calhit->setCellID1(cellid1);
	  calhit->setEnergy(calibr_coeff*energy);
	  calhit->setPosition(hit->getPosition());
	  calhit->setType((int)0);
	  calhit->setRawHit(hit);
	  muoncol->addElement(calhit);
	  LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.);
	  relcol->addElement( rel );
	}

      }
    }
    catch(DataNotAvailableException &e){ 
    }
  }
  muoncol->parameters().setValue(LCIO::CellIDEncoding,initString);
  evt->addCollection(muoncol,_outputMuonCollection.c_str());


  _nEvt++;

}


void SimpleMuonDigi::check( LCEvent * evt ) { }
  
void SimpleMuonDigi::end(){ } 
