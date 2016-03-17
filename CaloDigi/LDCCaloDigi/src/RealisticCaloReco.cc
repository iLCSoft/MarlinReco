#include "RealisticCaloReco.h"

#include <marlin/Global.h>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>

#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCRelationImpl.h>

#include "CalorimeterHitType.h"

#include <EVENT/LCParameters.h>
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

  // output collection names
  std::vector < std::string > outputHitCollections;
  registerProcessorParameter( "outputHitCollections",
			      "output hit collection names",
			      _outputHitCollections,
			      outputHitCollections);

  std::vector < std::string > outputRelCollections;
  registerProcessorParameter( "outputRelCollections",
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

  registerProcessorParameter("gap_correction" ,
                             "account for gaps" ,
                             _gapCorrection,
                             (int)0);

  registerProcessorParameter("CellIDLayerString" ,
                             "name of the part of the cellID that holds the layer" , 
                             _cellIDLayerString , 
                             std::string("K-1")
                             );

  registerProcessorParameter("CellIDModuleString" ,
                             "name of the part of the cellID that holds the module" , 
                             _cellIDModuleString , 
                             std::string("M")
                             );

  registerProcessorParameter("CellIDStaveString" ,
                             "name of the part of the cellID that holds the stave" , 
                             _cellIDStaveString , 
                             std::string("S-1")
                             );

  registerProcessorParameter("CellIDWaferString" ,
                             "name of the part of the cellID that holds the wafer" , 
                             _cellIDWaferString , 
                             std::string("wafer")
                             );

  registerProcessorParameter("CellIDTowerString" ,
                             "name of the part of the cellID that holds the tower" , 
                             _cellIDTowerString , 
                             std::string("tower")
                             );

  registerProcessorParameter("CellIDIndexIString" ,
                             "name of the part of the cellID that holds the index I" , 
                             _cellIDIndexIString , 
                             std::string("I")
                             );

  registerProcessorParameter("CellIDIndexJString" ,
                             "name of the part of the cellID that holds the index J" , 
                             _cellIDIndexJString , 
                             std::string("J")
                             );

}

void RealisticCaloReco::init() {
  printParameters();
  _countWarnings=0;
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

  // check that number of input and output collections names are the same
  assert ( _outputHitCollections.size() == _inputHitCollections.size() );
  assert ( _outputRelCollections.size() == _inputHitCollections.size() );

  assert ( _calibrCoeff.size()>0 );
  assert ( _calibrCoeff.size() == _calLayers.size() );

}


void RealisticCaloReco::processRunHeader( LCRunHeader* run) {
}

void RealisticCaloReco::processEvent( LCEvent * evt ) {

  // create the output relation collection -> don't really need this, since each calohit has a link to rawcalohit
  //  LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);

  _flag.setBit(LCIO::CHBIT_LONG);
  _flag.setBit(LCIO::RCHBIT_TIME); //store timing on output hits.

  resetGaps();

  // * Reading Collections of digitised calorimeter Hits *

  for (unsigned int i(0); i < _inputHitCollections.size(); ++i) {
    std::string colName =  _inputHitCollections[i] ;
    streamlog_out ( DEBUG ) << "looking for collection: " << colName << endl;
    try{
      LCCollection * col = evt->getCollection( colName.c_str() ) ;
      string initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);

      if ( _idDecoder ) delete _idDecoder;
      _idDecoder = new CellIDDecoder<CalorimeterHit> ( col );


      _cht_layout = layoutFromString(colName);
      _cht_caloid = caloIDFromString(colName);
      _cht_calotype = caloTypeFromString(colName);

      if ( _cht_layout==CHT::any     ) streamlog_out ( WARNING ) << "could not determine CHT::layout for " << colName << std::endl;
      if ( _cht_caloid==CHT::unknown ) streamlog_out ( WARNING ) << "could not determine CHT::caloID for " << colName << std::endl;

      getGeometryInformation();

      // create new collection
      LCCollectionVec *newcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
      newcol->setFlag(_flag.getFlag());

      // relation between digitised and reconstructed hits
      LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);

      // if making gap corrections clear the vectors holding pointers to calhits
      if(_gapCorrection>0){
	for(int is=0;is<MAX_STAVES;is++)
          for(int im=0;im<MAX_MODULES;im++)
	    for(int il=0;il<MAX_LAYERS;il++)
	      _calHitsByStaveModuleLayer[is][im][il].clear();
      }

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

	if (_gapCorrection>0) prepareForGaps(calhit);

	// output hit relations to input hits
	relcol->addElement( new LCRelationImpl(hit,calhit,1.0) );

      }

      // if requested apply gap correction
      if(_gapCorrection>0) {
	fillGaps();
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

