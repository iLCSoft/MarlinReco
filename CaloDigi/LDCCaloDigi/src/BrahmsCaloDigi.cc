#include "BrahmsCaloDigi.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <iostream>

using namespace lcio ;
using namespace marlin ;

BrahmsCaloDigi aBrahmsCaloDigi ;


BrahmsCaloDigi::BrahmsCaloDigi() : Processor("BrahmsCaloDigi") {

  _description = "BRAHMS Digitization..." ;
  
  std::vector<std::string> ecalCollections;

  ecalCollections.push_back(std::string("ecal"));

  registerProcessorParameter( "ECALCollections" , 
			      "ECAL Collection Names" ,
			      _ecalCollections ,
			       ecalCollections);

  std::vector<std::string> hcalCollections;

  hcalCollections.push_back(std::string("hcal"));

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
			     "Collection of relations" , 
			     _outputRelCollection , 
			     std::string("RelationCaloHit"));


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

  registerProcessorParameter("R1BarrelHcal",
			     "Inner Radius of HCAL Barrel",
			     _r1BarrelHcal,
			     (float)1908.);

  registerProcessorParameter("R2BarrelHcal",
			     "Outer Radius of HCAL Barrel",
			     _r2BarrelHcal,
			     (float)3000.);

  registerProcessorParameter("ModuleWallThick",
			     "Thickness of Module Wall",
			     _tWall,
			     (float)2.);

  registerProcessorParameter("TileSize",
			     "Transverse Tile Size",
			     _tileSize,
			     30);

  registerProcessorParameter("halfEndcapHole", 
			     "Half Size Endcap Hole",
			     _halfEndcapHole,
			     (float)300.);

  registerProcessorParameter("LayerThicknessHcal",
			     "Thickness ECAL Layer",
			     _layerThickHcal,
			     (float)26.5);

  registerProcessorParameter("zInLayer",
			     "z c-o-g in Layer",
			     _zInLayer, 
			     (float)20.0);

  registerProcessorParameter("nLayerBarrel1",
			     "Layers in 1 Barrel Section",
			     _nLayerBarrel1,
			     32);

  registerProcessorParameter("nLayerBarrel",
			     "Layers in Barrel Section",
			     _nLayerBarrel,
			     40);


  registerProcessorParameter("zBarrelHcal",
			     "z coordinate of Hcal Barrel",
			     _zBarrelHcal,
			     (float)2658.5);

  registerProcessorParameter("zEndcapHcal",
			     "z coordinate of Hcal Endcap",
			     _zEndcapHcal,
			     (float)2826.0);

}

void BrahmsCaloDigi::init() {
    _nRun = -1;
    _const_pi = acos(-1.0);
    _const_twopi = 2.0*_const_pi;
    _const_pi4 = _const_pi/4.;
    _const_pi8 = _const_pi/8.;
    _y0BarrelModule = _r1BarrelHcal * tan(_const_pi8) - 0.4;
    _x1BarrelModule = _layerThickHcal * _nLayerBarrel1;
    _y1BarrelModule = (_r1BarrelHcal + _x1BarrelModule) * tan(_const_pi8) - 0.4;
    _x2BarrelModule = _r2BarrelHcal - _r1BarrelHcal - 5.0 ;

    std::cout << std::endl;
    std::cout << "            HCAL Geometry " << std::endl;
    std::cout << "=======================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Z Barrel HCAL  = " << _zBarrelHcal << std::endl;
    std::cout << "Z Endcap HCAL  = " << _zEndcapHcal << std::endl;    
    std::cout << "R inner HCAL   = " << _r1BarrelHcal << std::endl;
    std::cout << "R outer HCAL   = " << _r2BarrelHcal << std::endl;
    std::cout << "Hole Half size = " << _halfEndcapHole << std::endl; 
    std::cout << std::endl;
    std::cout << "y0BarrelModule = " << _y0BarrelModule << std::endl;
    std::cout << "x1BarrelModule = " << _x1BarrelModule << std::endl;
    std::cout << "y1BarrelModule = " << _y1BarrelModule << std::endl;
    std::cout << "x2BarrelModule = " << _x2BarrelModule << std::endl;
    std::cout << std::endl;
    std::cout << "Tile Size      = " << _tileSize << std::endl; 
    
    if (_tileSize % 10 != 0) {
	std::cout << "Warning : tile size not multiple of 10 " << std::endl;
	std::cout << "Quitting program " << std::endl;
	exit(1);

    }

    std::cout << std::endl;
    std::cout << std::endl;
    


}


void BrahmsCaloDigi::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 

void BrahmsCaloDigi::processEvent( LCEvent * evt ) { 
    

    _HitVector.clear();
    
    LCCollectionVec *ecalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
    LCCollectionVec *hcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);

    LCCollectionVec *relcol = new LCCollectionVec(LCIO::LCRELATION);

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
		    LCRelationImpl * rel = new LCRelationImpl(calhit,hit,1.0);
		    ecalcol->addElement(calhit);
		    relcol->addElement(rel);
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

    int totSimulated(0);

    for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
      try{
	  LCCollection * col = evt->getCollection( _hcalCollections[i].c_str() ) ;
	  int numElements = col->getNumberOfElements();
	  for (int j(0); j < numElements; ++j) {
	      totSimulated++;
	      SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	      float energy = hit->getEnergy();

	      float xhit =  hit->getPosition()[0];
	      float yhit =  hit->getPosition()[1];
	      float zhit =  hit->getPosition()[2];
	      
	      int kcol = 0x00ffff;

	      float xnew;
	      float ynew;
	      float znew;

	      int Module;
	      int Stave;
	      int SubModule;
	      int icell;
	      int jcell;
	      int layer;

	      getCell(xhit,yhit,zhit,Module,Stave,SubModule,icell,jcell,layer);
	      int cellid_x = EncodeCellID(Module,Stave,SubModule,icell,jcell,layer);

	      int kExist = -1;
	      for (int k(0); k < (int)_HitVector.size(); ++k) {
		  int cellid_e = _HitVector[k]->getCellID();
		  if (cellid_x == cellid_e) {
		      kExist = k;
		      break;
		  }		  
	      }

	      if (kExist < 0) {
		  DigiHitExtended * newhit = new DigiHitExtended();
		  newhit->addSimHit(hit);
		  newhit->setAmpl(energy);
		  newhit->setCellID(cellid_x);
		  float pos[3];
		  getCoordinates(Module,Stave,SubModule,icell,jcell,layer,xnew,ynew,znew);
		  pos[0] = xnew;
		  pos[1] = ynew;
		  pos[2] = znew;
		  newhit->setPosition(pos);
		  _HitVector.push_back(newhit);
	      }
	      else {
		  float xx =  _HitVector[kExist]->getPosition()[0];
		  float yy =  _HitVector[kExist]->getPosition()[1];
		  float zz =  _HitVector[kExist]->getPosition()[2];
		  float distance = sqrt((xx-xhit)*(xx-xhit)+(yy-yhit)*(yy-yhit)+(zz-zhit)*(zz-zhit));
		  
		  float htsize = 0.5*(float)_tileSize;
		  htsize *= sqrt(2.0); 

		  if ( distance > htsize) 
		      std::cout << "Warning : found existing hit sith distance to cell center : " << distance << std::endl; 
 
		  float hitenergy = _HitVector[kExist]->getAmpl();
		  hitenergy += energy;
		  _HitVector[kExist]->setAmpl( hitenergy );
		  _HitVector[kExist]->addSimHit( hit ) ;
	      }



	  }
      }
      catch(DataNotAvailableException &e){ 
      }
  }


    for (int i(0); i < (int)_HitVector.size(); ++i) {

	DigiHitExtended * hit = _HitVector[i];
	float energy = hit->getAmpl();
      
	float xnew = hit->getPosition()[0];
	float ynew = hit->getPosition()[1];
	float znew = hit->getPosition()[2];
      
	int kcol = 0xff0000;


	if (energy > _thresholdHcal) {
	    CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	    int cellid = hit->getCellID();
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
	    SimCalorimeterHitVec hitvec = hit->getSimHitVector();
	    for (unsigned int ish(0); ish < hitvec.size(); ++ish) {
		SimCalorimeterHit * shit = hitvec[ish];
		LCRelationImpl * rel = new LCRelationImpl(calhit,shit,1.);
		relcol->addElement(rel);
	    }
	    hcalcol->addElement(calhit);
	}
      
    }


    std::cout << "Event number = " << _nEvt << std::endl; 
    std::cout << "Total number of simulated hits : " << totSimulated << std::endl;
    std::cout << "Total number of digitised hits : " << _HitVector.size() << std::endl;
    
    
    evt->addCollection(hcalcol,_outputHcalCollection.c_str());
    evt->addCollection(relcol,_outputRelCollection.c_str());
    
    getchar();

    CleanUp();

    _nEvt++;

}


void BrahmsCaloDigi::getCell(float xhit, float yhit, float zhit, int & Module, int & Stave, int & SubModule, int & icell, int & jcell, int & layer) {

    if (zhit < - _zBarrelHcal - 0.1 ) {
	Module = 0;
    }
    else if (zhit >= -_zBarrelHcal && zhit < 0.) {
	Module = 1;
    }	
    else if (zhit >= 0. &&  zhit < _zBarrelHcal + 0.1) { 
	Module = 2;
    }
    else {
	Module = 3;
    }


    if (Module == 0 || Module == 3) {
	getCellEndcap(xhit, yhit, zhit, Stave, icell, jcell, layer);
	SubModule = 0;
    }
    else {
	getCellBarrel(xhit, yhit, zhit, Stave, SubModule, icell, jcell, layer);
    }

}

void BrahmsCaloDigi::getCellEndcap(float xhit, float yhit, float zhit, int & Stave, int & icell, int & jcell, int & layer) {

    if (xhit > _halfEndcapHole && yhit < _halfEndcapHole) {
	Stave = 0;
	icell = (int)(( _halfEndcapHole - yhit )/(float)_tileSize);
	jcell = (int)(( xhit - _halfEndcapHole )/(float)_tileSize);
    }
    else if (xhit > -_halfEndcapHole  && yhit > _halfEndcapHole) {
	Stave = 1;
	icell = (int)(( xhit + _halfEndcapHole )/(float)_tileSize);
	jcell = (int)(( yhit - _halfEndcapHole )/(float)_tileSize);
    }
    else if (xhit < -_halfEndcapHole  && yhit > -_halfEndcapHole) {
	Stave = 2;
	icell = (int)(( yhit + _halfEndcapHole )/(float)_tileSize);
	jcell = (int)(( - _halfEndcapHole - xhit )/(float)_tileSize);
    }
    else if (xhit < _halfEndcapHole  && yhit < -_halfEndcapHole) {
	Stave = 3;
	icell = (int)(( _halfEndcapHole - xhit )/(float)_tileSize);
	jcell = (int)((- _halfEndcapHole - yhit )/(float)_tileSize);
    }
    
    if (zhit < 0.) {
	layer = (int)((-_zEndcapHcal - zhit)/_layerThickHcal);
    }
    else {
	layer = (int)((zhit - _zEndcapHcal)/_layerThickHcal);
    }

}

void BrahmsCaloDigi::getCellBarrel(float xhit, float yhit, float zhit, int & Stave, int & SubModule, int & icell, int & jcell, int & layer) {

    jcell = (int)((fabs(zhit) - _tWall)/(float)_tileSize);
   
    float phi = atan2(yhit,xhit) + _const_pi8;

    if (phi < 0.) 
	phi = phi + _const_twopi;

    Stave = (int)(phi/_const_pi4);

    float radius = sqrt(xhit*xhit + yhit*yhit);

    float phiInStave = phi - Stave * _const_pi4 - _const_pi8;

    float xInStave = radius * cos(phiInStave) - _r1BarrelHcal;
    float yInStave = radius * sin(phiInStave);

    if (yInStave < 0.) {
	SubModule = 0;
    }
    else {
	SubModule = 1;
    }

    layer = (int)(xInStave /_layerThickHcal); 
    
    yInStave = fabs(yInStave);

    icell = (int)((yInStave - _tWall)/(float)_tileSize);

}


void BrahmsCaloDigi::getCoordinates(int Module, int Stave, int SubModule, int icell, int jcell, int layer, float & xnew, float & ynew, float & znew) {

    if (Module == 0 | Module == 3) {
	getCoordinatesEndcap(Module,Stave,icell,jcell,layer,xnew,ynew,znew);
    }
    else {
	getCoordinatesBarrel(Module,Stave,SubModule,icell,jcell,layer,xnew,ynew,znew);
    }


}

void BrahmsCaloDigi::getCoordinatesEndcap(int Module, int Stave, int icell, int jcell, int layer, float & xnew, float & ynew, float & znew) {

    if (Module == 0) {
	znew = -_zEndcapHcal - layer * _layerThickHcal - _zInLayer;
    }
    else if (Module == 3) {
	znew = _zEndcapHcal + layer * _layerThickHcal + _zInLayer;
    }

    if (Stave == 0) {
	xnew = _halfEndcapHole + 0.5*(float)_tileSize + jcell * _tileSize;
	ynew = _halfEndcapHole - 0.5*(float)_tileSize - icell * _tileSize;
    }
    else if (Stave == 1) {
	xnew = -_halfEndcapHole + 0.5*(float)_tileSize + icell * _tileSize;
	ynew =  _halfEndcapHole + 0.5*(float)_tileSize + jcell * _tileSize;
    }
    else if (Stave == 2) {
	xnew = -_halfEndcapHole - 0.5*(float)_tileSize - jcell * _tileSize;
	ynew = -_halfEndcapHole + 0.5*(float)_tileSize + icell * _tileSize;
    }
    else if (Stave == 3) {
	xnew = _halfEndcapHole - 0.5*(float)_tileSize - icell * _tileSize;
	ynew = -_halfEndcapHole - 0.5*(float)_tileSize - jcell * _tileSize;
    }



}

void BrahmsCaloDigi::getCoordinatesBarrel(int Module, int Stave, int SubModule, int icell, int jcell, int layer, float & xnew, float & ynew, float & znew) {

    if (Module == 1) {
	znew = -_tWall - 0.5*(float)_tileSize - jcell * float(_tileSize);
    }
    else {
	znew = +_tWall + 0.5*(float)_tileSize + jcell * float(_tileSize);
    }

    float xInStave = layer * _layerThickHcal + _zInLayer;
    float yInStave = _tWall + 0.5*(float)_tileSize + icell * float(_tileSize);

    float yEdge;

    if (layer < _nLayerBarrel1) {
	float stepy = _y1BarrelModule - _y0BarrelModule;
	yEdge = _y0BarrelModule + xInStave * stepy / _x1BarrelModule;
    }
    else {
	float stepx = (xInStave - _x1BarrelModule)/(_x2BarrelModule - _x1BarrelModule);
	yEdge = _y1BarrelModule - _y1BarrelModule * stepx;
    }

    if (yInStave > yEdge)
	yInStave = yEdge - 5.;

    xInStave = xInStave + _r1BarrelHcal;

    if (SubModule == 0) 
	yInStave = - yInStave;

    float radius = sqrt(xInStave*xInStave+yInStave*yInStave);
    float phi = atan2(yInStave,xInStave);

    phi = phi + Stave * _const_pi4;

    xnew = radius * cos(phi);
    ynew = radius * sin(phi);

}


int BrahmsCaloDigi::EncodeCellID(int Module, int Stave, int SubModule, int icell, int jcell, int layer) {

    int Code = layer << 24;
    Code = Code | Module << 22;
    Code = Code | Stave << 19;
    Code = Code | SubModule << 18;
    Code = Code | icell << 9;
    Code = Code | jcell ;

    return Code;

}

void BrahmsCaloDigi::CleanUp() {

    for (unsigned int i(0); i < _HitVector.size(); ++i) {
	DigiHitExtended * hit =_HitVector[i];
	delete hit;
    }

    _HitVector.clear();

}

void BrahmsCaloDigi::check( LCEvent * evt ) { }
  
void BrahmsCaloDigi::end(){ } 
