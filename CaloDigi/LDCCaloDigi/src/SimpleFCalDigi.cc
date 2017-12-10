#include "SimpleFCalDigi.h"
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
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/CalorimeterParameters.h>
#include <gear/GearParameters.h>
#include <gear/LayerLayout.h>
#include "CalorimeterHitType.h"

using namespace std;
using namespace lcio ;
using namespace marlin ;


SimpleFCalDigi aSimpleFCalDigi ;


SimpleFCalDigi::SimpleFCalDigi() : Processor("SimpleFCalDigi") {

  _description = "Performs simple digitization of SimCalorimeter hits in forward calorimeters ..." ;
  
  std::vector<std::string> fcalCollections;

  fcalCollections.push_back(std::string("LcalCollection"));
  
  registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			    "FCALCollections" , 
			    "Fcal Collection Names" ,
			    _fcalCollections ,
			    fcalCollections);
    
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "FCALOutputCollection" , 
			    "Fcal Collection of real Hits" , 
			    _outputFcalCollection , 
			    std::string("FCAL")) ; 
  
  registerOutputCollection( LCIO::LCRELATION, 
			    "RelationOutputCollection" , 
			    "CaloHit Relation Collection" , 
			    _outputRelCollection , 
			    std::string("RelationFcalHit")) ; 
  
  registerProcessorParameter("FcalThreshold" , 
			     "Threshold for Fcal Hits in GeV" ,
			     _thresholdFcal,
			     (float)0.0);

  registerProcessorParameter("CalibrFCAL" , 
			     "Calibration coefficients for FCAL" ,
			     _calibrCoeffFcal,
			     (float)31.);

  registerProcessorParameter("CellIDLayerString" ,
			     "name of the part of the cellID that holds the layer" , 
			     _cellIDLayerString , 
			     std::string("K-1")
			     );

  registerProcessorParameter("FixLCalHits" ,
			     "Fix the hit positions in LCal using the cellID (for DBD simulated samples)" , 
			     _fixLCalHits , 
			     bool(false)
			     );
                                                          
  registerProcessorParameter("CaloType" ,
                             "type of calorimeter: em, had, muon" , 
                             _caloType , 
                             std::string("had")
                             );

  registerProcessorParameter("CaloID" ,
			     "ID of calorimeter: lcal, fcal, bcal", 
			     _caloID , 
			     std::string("fcal")
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

void SimpleFCalDigi::init() {

  _nRun = -1;
  _nEvt = 0;

  //fg: need to set default encoding in for reading old files...
  //CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;
  CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding(_defaultEncoding.c_str()) ;
  if ( ! _caloID.compare("lcal")  && // true if it is false ... 
                 _fixLCalHits          ) {  
            // parametrs for fixing wrong cellID to xyz coding in LCal Mokka

    const gear::CalorimeterParameters& lcalparam = Global::GEAR->getLcalParameters();

    xing_angle =   lcalparam.getDoubleVal("beam_crossing_angle");       // 1.400000000e+01 ;  beam_crossing_angle (mrad)
    zMin       =   lcalparam.getExtent()[2] ;                           // 2.506900000e+03 ;  inner_z
    dZ         =   lcalparam.getLayerLayout().getThickness(0) ;         // 4.290000000e+00 ;  thickness
    rMin       =   lcalparam.getExtent()[0] ;                           // 8.402209443e+01 ;  inner_r
    cellDimR   =   lcalparam.getLayerLayout().getCellSize0(0);          // 1.718404774e+00 ;  cellSize0
    cellDimPhi =   lcalparam.getLayerLayout().getCellSize1(0);          // 1.308996939e-01 ;  cellSize1
    WThickness =   lcalparam.getLayerLayout().getAbsorberThickness(0);  // 3.500000000e+00 ;  absorberThickness
  }

}


void SimpleFCalDigi::processRunHeader( LCRunHeader*  /*run*/) { 
  _nRun++ ;
  _nEvt = 0;
} 

void SimpleFCalDigi::processEvent( LCEvent * evt ) { 
    

  LCCollectionVec *lcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);

  LCFlagImpl flag;

  flag.setBit(LCIO::CHBIT_LONG);

  lcalcol->setFlag(flag.getFlag());

  // 
  // * Reading Collections of FCAL Simulated Hits * 
  // 
  string initString;
  for (unsigned int i(0); i < _fcalCollections.size(); ++i) {
    try{
      LCCollection * col = evt->getCollection( _fcalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      int numElements = col->getNumberOfElements();
      CellIDDecoder<SimCalorimeterHit> idDecoder( col );
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	float energy = hit->getEnergy();

	if (energy > _thresholdFcal) {
	  CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	  int cellid = hit->getCellID0();
	  int cellid1 = hit->getCellID1();
	  float calibr_coeff(1.);
	  calibr_coeff = _calibrCoeffFcal;
	  calhit->setCellID0(cellid);
	  calhit->setCellID1(cellid1);
	  calhit->setEnergy(calibr_coeff*energy);
	  float pos[3];
	  if ( ! _caloID.compare("lcal")  && // true if it is false ... 
                 _fixLCalHits          ) {  
            // fix for wrong cellID to xyz coding in LCal Mokka

            int i = idDecoder(hit)[ "I" ] ;
            int j = idDecoder(hit)[ "J" ] ;
            int k = idDecoder(hit)[ "K" ] ;
            int s = idDecoder(hit)[ "S-1" ] ;

	    float z_from_cell = ( zMin + (k-1) * dZ + (dZ- WThickness)/2.0 ) - 0.11 ; // 0.11 is a guess on the sensorthikness
            float r_from_cell =  i*cellDimR + cellDimR/2.0 + rMin; 
            int   oddeven = ((k+1)%2) ;
            float phi_from_cell = j*cellDimPhi + cellDimPhi/2.0 - oddeven*cellDimPhi/2.0;

            float angle= -(xing_angle*1.0e-3 / 2.0)* (2 * ( s - 0.5 ) ) - ( s - 1 )*3.14159  ; 
            float rotated_z_from_cell=    z_from_cell*cos( angle) + r_from_cell*cos(phi_from_cell)*sin( angle);
	    float rotated_x_from_cell=  - z_from_cell*sin( angle) + r_from_cell*cos(phi_from_cell)*cos( angle);

            pos[0] = rotated_x_from_cell ; 
            pos[1] =  r_from_cell*sin(phi_from_cell) ; 
            pos[2]=rotated_z_from_cell;

	    streamlog_out( DEBUG3 ) << "i,j,k,s : " << i << " " << j << " " << k << " " << s << " " << oddeven << std::endl;
            streamlog_out( DEBUG3 ) << "xyzr  input                         : " << hit->getPosition()[0] << " " << hit->getPosition()[1]  << " " << hit->getPosition()[2] << " " << 
	      sqrt(hit->getPosition()[0]*hit->getPosition()[0]+hit->getPosition()[1]*hit->getPosition()[1]) << std::endl;
	    streamlog_out( DEBUG3 ) << "xyzr with r and z fr cellID rotated : " << pos[0] << " " << pos[1] << " " << pos[2] << " " << 
              sqrt(pos[0]*pos[0]+pos[1]*pos[1]) << std::endl;
	    streamlog_out( DEBUG2 ) << " hit and diffs "  << hit->getPosition()[0] << " " << hit->getPosition()[1]  << " " << hit->getPosition()[2] << " " 
                                    <<  hit->getPosition()[0]-pos[0] << " " << hit->getPosition()[1]-pos[1]  << " " << hit->getPosition()[2]-pos[2] <<std::endl;
           streamlog_out( DEBUG3 ) << std::endl;
            

          } else {
  	    pos[0] = hit->getPosition()[0];
	    pos[1] = hit->getPosition()[1];
	    pos[2] = hit->getPosition()[2];
          }
	  calhit->setPosition(pos);

          calhit->setType( CHT ( caloTypeFromString(_caloType), caloIDFromString(_caloID), layoutFromString(_caloLayout.c_str()),  idDecoder(hit)[ _cellIDLayerString ] ) );

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
  evt->addCollection(lcalcol,_outputFcalCollection.c_str());
  evt->addCollection(relcol,_outputRelCollection.c_str());


  _nEvt++;

}


void SimpleFCalDigi::check( LCEvent *  /*evt*/ ) { }
  
void SimpleFCalDigi::end(){ } 
