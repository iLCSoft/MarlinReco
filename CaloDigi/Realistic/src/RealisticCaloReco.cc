#include "RealisticCaloReco.h"

#include <marlin/Global.h>

#include <EVENT/LCCollection.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCRelationImpl.h>

#include <UTIL/LCRelationNavigator.h> 

#include "CalorimeterHitType.h"

#include <EVENT/LCParameters.h>
#include <EVENT/SimCalorimeterHit.h>

#include <iostream>
#include <string>
#include <algorithm>
#include <assert.h>
#include <cmath>

using namespace std;
using namespace lcio ;
using namespace marlin ;

RealisticCaloReco::RealisticCaloReco() : Processor("RealisticCaloReco") {

  _description = "Performs simple reconstruction of calo hits..." ;

  std::vector < std::string > inputHitCollections;
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "inputHitCollections",
			    "input hit collection names",
			    _inputHitCollections,
			    inputHitCollections);

  std::vector < std::string > inputRelCollections;
  registerInputCollections( LCIO::LCRELATION,
			    "inputRelationCollections",
			    "input relation collection names (digi<->sim), one per inputHitCollection",
			    _inputRelCollections,
			    inputRelCollections);

  // output collection names
  std::vector < std::string > outputHitCollections;
  registerProcessorParameter( "outputHitCollections",
			      "output hit collection names",
			      _outputHitCollections,
			      outputHitCollections);

  std::vector < std::string > outputRelCollections;
  registerProcessorParameter( "outputRelationCollections",
			      "output hit collection names",
			      _outputRelCollections,
			      outputRelCollections);

  std::vector<int> calLayers;
  registerProcessorParameter("calibration_layergroups" ,
                             "grouping of calo layers" ,
                             _calLayers,
                             calLayers);

  std::vector<float> calibrCoeff;
  registerProcessorParameter("calibration_factorsMipGev" ,
                             "Calibration coefficients (MIP->shower GeV) of layers groups" ,
                             _calibrCoeff,
                             calibrCoeff);

  registerProcessorParameter("CellIDLayerString" ,
                             "name of the part of the cellID that holds the layer" , 
                             _cellIDLayerString , 
                             std::string("K-1")
                             );


}

void RealisticCaloReco::init() {
  printParameters();
  _idDecoder=NULL;

  // if no output collection names specified, set some default based on the input collection names
  if ( _outputHitCollections.size()==0 ) {
    for (size_t i=0; i<_inputHitCollections.size(); i++) {
      _outputHitCollections.push_back( _inputHitCollections[i] + "Reco" );
    }
  }
  if ( _outputRelCollections.size()==0 ) {
    for (size_t i=0; i<_inputHitCollections.size(); i++) {
      _outputRelCollections.push_back( _inputHitCollections[i] + "DigiRelation" );
    }
  }

  // should be one input relation collection per input hit collection
  assert( _inputRelCollections.size() == _inputHitCollections.size() );

  // check that number of input and output collections names are the same
  assert ( _outputHitCollections.size() == _inputHitCollections.size() );
  assert ( _outputRelCollections.size() == _inputHitCollections.size() );

  assert ( _calibrCoeff.size()>0 );
  assert ( _calibrCoeff.size() == _calLayers.size() );

}


void RealisticCaloReco::processRunHeader( LCRunHeader* run) {
}

void RealisticCaloReco::processEvent( LCEvent * evt ) {

  _flag.setBit(LCIO::CHBIT_LONG);
  _flag.setBit(LCIO::RCHBIT_TIME); //store timing on output hits.

  // * Reading Collections of digitised calorimeter Hits *

  for (unsigned int i(0); i < _inputHitCollections.size(); ++i) {
    std::string colName =  _inputHitCollections[i] ;
    std::string relName =  _inputRelCollections[i] ;
    streamlog_out ( DEBUG ) << "looking for hit, relation collection: " << colName << " " << relName << endl;

    try{
      LCCollection * col = evt->getCollection( colName.c_str() ) ;
      string initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);

      LCCollection * inrelcol = evt->getCollection( relName.c_str() ) ;
      LCRelationNavigator navi(inrelcol);

      if ( _idDecoder ) delete _idDecoder;
      _idDecoder = new CellIDDecoder<CalorimeterHit> ( col );

      // create new collection
      LCCollectionVec *newcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
      newcol->setFlag(_flag.getFlag());

      // relation between digitised and reconstructed hits
      LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);

      int numElements = col->getNumberOfElements();
      streamlog_out ( DEBUG ) << colName << " number of elements = " << numElements << endl;

      for (int j(0); j < numElements; ++j) {
        CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;

	int cellid0 = hit->getCellID0();
	int cellid1 = hit->getCellID1();

	CalorimeterHitImpl * calhit = new CalorimeterHitImpl(); // make new hit

	float energy = reconstructEnergy( hit ); // overloaded method, technology dependent

	calhit->setCellID0(cellid0);
	calhit->setCellID1(cellid1);
	calhit->setEnergy(energy);
	calhit->setRawHit( hit->getRawHit() );
	calhit->setTime( hit->getTime() );
	calhit->setPosition( hit->getPosition() );
	calhit->setType( hit->getType() );

	newcol->addElement( calhit );

	// get the simcalohit corresponding to this digitised hit
	if ( navi.getRelatedFromObjects( hit ) .size() > 0 ) {
	  SimCalorimeterHit* simhit = (SimCalorimeterHit*) navi.getRelatedFromObjects(hit)[0]; // assume the first one (should be only one)
	  // make a relation, add to collection - now keep relations between sim and reco hits
	  relcol->addElement( new LCRelationImpl(simhit,calhit,1.0) );
	} else {
	  streamlog_out ( WARNING ) << "could not find relation to sim calo hit!" << endl;
	}

      }

      // add collection to event
      newcol->parameters().setValue(LCIO::CellIDEncoding,initString);
      
      evt->addCollection(newcol,_outputHitCollections[i].c_str());
      evt->addCollection(relcol, _outputRelCollections[i].c_str());
    }
    catch(DataNotAvailableException &e){
      streamlog_out(DEBUG) << "could not find input ECAL collection " << colName << std::endl;
    }
  }
  return;
}

float RealisticCaloReco::getLayerCalib( int ilayer ) {
  float calib_coeff = 0;
  // retrieve calibration constants
  for (unsigned int k(0); k < _calLayers.size(); ++k) {
    int min,max;
    if (k == 0){
      min = 0;
    }else{
      min = _calLayers[k-1];
    }
    max = _calLayers[k];
    if (ilayer >= min && ilayer < max) {
      calib_coeff = _calibrCoeff[k];
      break;
    }
  }
  assert( calib_coeff>0 );
  return calib_coeff;
}

void RealisticCaloReco::check( LCEvent * evt ) { }

void RealisticCaloReco::end(){ }

