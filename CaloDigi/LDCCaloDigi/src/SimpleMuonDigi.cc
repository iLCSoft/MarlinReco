#include "SimpleMuonDigi.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// #include <algorithm>
// #include <string>
#include <cctype> 
#include <cstdlib>  // abs

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>


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
			     (float)0.025);

  registerProcessorParameter("CalibrMUON" , 
			     "Calibration coefficients for MUON" ,
			     _calibrCoeffMuon,
			     (float)120000.);

  registerProcessorParameter("MaxHitEnergyMUON", 
			     "maximum hit energy for a MUON hit" ,
			     _maxHitEnergyMuon,
			     (float)2.0);

  IntVec keepBarrelLayersVec, keepEndcapLayersVec;

  registerProcessorParameter("KeepBarrelLayersVec" , 
			     "Vector of Barrel layers to be kept. Layers start at 1!",
			     _layersToKeepBarrelVec,
			     keepBarrelLayersVec);

  registerProcessorParameter("KeepEndcapLayersVec" , 
			     "Vector of Endcap layers to be kept. Layers start at 1!",
			     _layersToKeepEndcapVec,
			     keepEndcapLayersVec);


  registerProcessorParameter("CellIDLayerString" ,
			     "name of the part of the cellID that holds the layer" , 
			     _cellIDLayerString , 
			     std::string("K-1")
			     );
  
}

void SimpleMuonDigi::init() {

  _nRun = -1;
  _nEvt = 0;

  //fg: need to set default encoding in for reading old files...
  CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;


  //Get the number of Layers in the Endcap
  int layersEndcap=0, layersBarrel=0;

  try{
    layersBarrel =  Global::GEAR->getYokeBarrelParameters().getLayerLayout().getNLayers();
  }catch( gear::UnknownParameterException& e ){
    streamlog_out(WARNING) << "  oops - no Yoke Barrel available " << std::endl ;
  }
  try{
    layersEndcap =  Global::GEAR->getYokeEndcapParameters().getLayerLayout().getNLayers();
  }catch( gear::UnknownParameterException& e ){
    streamlog_out(WARNING) << "  oops - no Yoke Endcap available " << std::endl ;
  }

  //If the vectors are empty, we are keeping everything 
  if(_layersToKeepBarrelVec.size() > 0) {
    //layers start at 0
    for(int i = 0; i < layersBarrel; ++i) {
      _useLayersBarrelVec.push_back(false);
      for(IntVec::iterator iter = _layersToKeepBarrelVec.begin(); iter < _layersToKeepBarrelVec.end(); ++iter) {
	if (i == *iter-1){
	  _useLayersBarrelVec[i]=true; break;
	}
      }
    }
  }

  if(_layersToKeepEndcapVec.size() > 0) {
    //layers start at 0
    for(int i = 0; i < layersEndcap; ++i) {
      _useLayersEndcapVec.push_back(false);
      for(IntVec::iterator iter = _layersToKeepEndcapVec.begin(); iter < _layersToKeepEndcapVec.end(); ++iter) {
	if (i == *iter-1){
	  _useLayersEndcapVec[i]=true; break;
	}
      }
    }
  }



}


void SimpleMuonDigi::processRunHeader( LCRunHeader*  /*run*/) { 
  _nRun++ ;
  _nEvt = 0;
} 

void SimpleMuonDigi::processEvent( LCEvent * evt ) { 
    

  streamlog_out( DEBUG ) << " process event : " << evt->getEventNumber() 
			 << " - run  " << evt->getRunNumber() << std::endl ;


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

    std::string colName =  _muonCollections[i] ;
    
    //fg: need to establish the yoke subdetetcor part here 
    //    use collection name as cellID does not seem to have that information
    CHT::Layout caloLayout = layoutFromString( colName ) ; 

    try{
      LCCollection * col = evt->getCollection( _muonCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      int numElements = col->getNumberOfElements();
      CellIDDecoder<SimCalorimeterHit> idDecoder( col );
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	float energy = hit->getEnergy();
	int cellid = hit->getCellID0();
	int cellid1 = hit->getCellID1();
	//Get The LayerNumber 
	unsigned int layer = abs( idDecoder(hit)[ _cellIDLayerString ] ) ;
	//Check if we want to use this layer, else go to the next hit
	if( !useLayer(caloLayout, layer) ) continue;
	float calibr_coeff(1.);
	calibr_coeff = _calibrCoeffMuon;
	float hitEnergy = calibr_coeff*energy;
	if(hitEnergy>_maxHitEnergyMuon)hitEnergy=_maxHitEnergyMuon;
	if (hitEnergy > _thresholdMuon) {
	  CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	  calhit->setCellID0(cellid);
	  calhit->setCellID1(cellid1);
	  calhit->setEnergy(hitEnergy);
	  calhit->setPosition(hit->getPosition());
	  calhit->setType( CHT( CHT::muon, CHT::yoke, caloLayout ,  idDecoder(hit)[ _cellIDLayerString ] ) );
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
  evt->addCollection(relcol,_outputRelCollection.c_str());


  _nEvt++;

}


void SimpleMuonDigi::check( LCEvent *  /*evt*/ ) { }
  
void SimpleMuonDigi::end(){ } 

bool SimpleMuonDigi::useLayer(CHT::Layout caloLayout,  unsigned int layer) {
  switch (caloLayout){
  case CHT::barrel:
    if(layer > _useLayersBarrelVec.size() || _useLayersBarrelVec.size() == 0) return true;
    return _useLayersBarrelVec[layer]; //break not needed, because of return
  case CHT::endcap:
    if(layer > _useLayersEndcapVec.size() || _useLayersEndcapVec.size() == 0) return true;
    return _useLayersEndcapVec[layer]; //break not needed, because of return
  //For all other cases, always keep the hit
  default:
    return true;
  }
}//useLayer
