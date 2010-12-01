#include "MokkaCaloDigi.h"
#include <iostream>
#include <cmath>
#include <vector>
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif
#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <lcio.h>
#include <iterator>
#define SHIFT_M 0
#define SHIFT_S 3
#define SHIFT_I 6
#define SHIFT_J 15
#define SHIFT_K 24
#define SHIFT_2 30
#define SHIFT_1 31

#define MASK_M (unsigned int) 0x00000007
#define MASK_S (unsigned int) 0x00000038
#define MASK_I (unsigned int) 0x00007FC0
#define MASK_J (unsigned int) 0x00FF8000
#define MASK_K (unsigned int) 0x3F000000
#define MASK_2 (unsigned int) 0x40000000
#define MASK_1 (unsigned int) 0x80000000

using namespace lcio ;
using namespace marlin ;
using namespace std;
using namespace IMPL;

#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>

MokkaCaloDigi aMokkaCaloDigi ;

MokkaCaloDigi::MokkaCaloDigi() : Processor("MokkaCaloDigi") {
  
  // modify processor description
  _description = "Mokka digitizer..." ;
  
  // register steering parameters: name, description, class-variable, default value
 
  std::vector<std::string> ecalCollections;
  ecalCollections.push_back(std::string("SEcal01_EcalEndcap"));
  ecalCollections.push_back(std::string("SEcal01_EcalBarrel"));

  registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			    "ECALCollections" , 
			    "ECAL Collection Names" ,
			    _ecalCollections ,
			    ecalCollections);
  
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("SHcal01_HcalBarrelEnd"));
  hcalCollections.push_back(std::string("SHcal01_HcalBarrelReg"));
  hcalCollections.push_back(std::string("SHcal01_HcalEndCaps"));
  
  registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			    "HCALCollections" , 
			    "HCAL Collection Names" , 
			    _hcalCollections , 
			    hcalCollections);
  
  
  registerOutputCollection(  LCIO::CALORIMETERHIT,
			     "NewECALCollName" , 
			     "name for the new collection of ECAL hits"  ,
			     _newCollNameECAL ,
			     std::string("ECAL")) ;
  
  
  registerOutputCollection(  LCIO::CALORIMETERHIT,
			     "NewHCALCollName" , 
			    "name for the new collection of HCAL hits"  ,
			    _newCollNameHCAL ,
			    std::string("HCAL"));
  

  registerOutputCollection(  LCIO::LCRELATION,
			     "RelationCollName" ,
			     "name for collection of relations between CalorimeterHits and SimCalorimeterHits", 
			     _relationCollName , 
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


void MokkaCaloDigi::init() { 

   printParameters() ;

   const gear::CalorimeterParameters& pHcalBarrel = Global::GEAR->getHcalBarrelParameters();
   const gear::CalorimeterParameters& pHcalEndcap = Global::GEAR->getHcalEndcapParameters();

  _nRun = -1 ;
  _nEvt = 0 ;
  int end_module_type = pHcalBarrel.getIntVal("Hcal_barrel_end_module_type");
  float tpc_ecal_hcalbarrel_halfz = float(pHcalBarrel.getDoubleVal("TPC_Ecal_Hcal_barrel_halfZ"));
  _lateralPlateThickness = float(pHcalBarrel.getDoubleVal("Hcal_lateral_structure_thickness"));
  _modulesGap = float(pHcalBarrel.getDoubleVal("Hcal_modules_gap"));
  float staves_gap = float(pHcalBarrel.getDoubleVal("Hcal_stave_gaps")); 
  _innerHcalRadius = float(pHcalBarrel.getExtent()[0]);
  float top_end_dim_z;
  _numberOfHcalLayers = int(pHcalBarrel.getLayerLayout().getNLayers());
  _hcalLayerThickness = float(pHcalBarrel.getLayerLayout().getThickness(0));
  _hcalAbsorberThickness = float(pHcalBarrel.getLayerLayout().getAbsorberThickness(0));
  _hcalSensitiveThickness = _hcalLayerThickness - _hcalAbsorberThickness;
  float back_plate_thickness = float(pHcalBarrel.getDoubleVal("Hcal_back_plate_thickness"));
  if(end_module_type == 1 ) {
    _regularBarrelModuleLength = 2*tpc_ecal_hcalbarrel_halfz/5.-1.;
    top_end_dim_z = 1180.0;
  }
  else {
      // the 140 and the proportions comes from Tesla TDR
    float total_z_size = 140 +
      2 * tpc_ecal_hcalbarrel_halfz;
    _regularBarrelModuleLength = total_z_size / 5. * 1080./1120.;
    top_end_dim_z = 
      (total_z_size - 3 * _regularBarrelModuleLength)/2;      
  }
  
  _regularBarrelChamberLength = _regularBarrelModuleLength - 2*_lateralPlateThickness;
  _virtualCellSizeX = float(pHcalBarrel.getDoubleVal("Hcal_virtual_cell_size"));
  _virtualCellSizeZ = _regularBarrelChamberLength/int(_regularBarrelChamberLength/_virtualCellSizeX);  
  float Pi = acos(-1.0);
  float total_dim_y = _numberOfHcalLayers*_hcalLayerThickness+back_plate_thickness;
  float module_radius = _innerHcalRadius + total_dim_y ;
  float y_dim2_for_x = module_radius - module_radius*cos(Pi/8.);
  float y_dim1_for_x = total_dim_y - y_dim2_for_x;  
  float bottom_dim_x = 2.*_innerHcalRadius*tan(Pi/8.)-staves_gap;
  float middle_dim_x = bottom_dim_x + 2*y_dim1_for_x*tan(Pi/8.);
  float y_dim1_for_z = 134.8;
  float y_dim2_for_z = (top_end_dim_z-_regularBarrelModuleLength)*tan(Pi*39.28/180.);
  _cellScaleX = int((float(pHcalBarrel.getLayerLayout().getCellSize0(0)+0.01))/_virtualCellSizeX);
  _cellScaleZ = int((float(pHcalBarrel.getLayerLayout().getCellSize1(0)+0.01))/_virtualCellSizeX);
  _newCellSizeX  = _virtualCellSizeX * _cellScaleX;
  _newCellSizeZ  = _virtualCellSizeZ * _cellScaleZ;
  _endBarrelChamberLength = new float[_numberOfHcalLayers];
  _barrelLateralWidth = new float[_numberOfHcalLayers];
  _barrelOffsetMaxX = new float[_numberOfHcalLayers];
  _endBarrelOffsetMaxZ = new float[_numberOfHcalLayers];

  for (int i = 0; i < _numberOfHcalLayers; i++) {
    float x_width = 0;
    float z_width = 0;
    int id = i + 1;
    float z = id*_hcalLayerThickness;
    if (z < y_dim1_for_x) {
      float x = bottom_dim_x + 2*(id*_hcalAbsorberThickness+(id-1)*_hcalSensitiveThickness)*tan(Pi/8.);
      x_width = int(x/_virtualCellSizeX)*_virtualCellSizeX;
    }
    else {
      float x = middle_dim_x -2.*(id*_hcalAbsorberThickness+(id-1)*_hcalSensitiveThickness - y_dim1_for_x)/tan(Pi/8.)
	- 2*_hcalSensitiveThickness/tan(Pi/8.) -2*back_plate_thickness;
      x_width = int(x/_virtualCellSizeX)*_virtualCellSizeX;
    }

    if (z < y_dim1_for_z) {
      z_width = _regularBarrelChamberLength;
    }

    if (z >= y_dim1_for_z && z < (y_dim1_for_z+y_dim2_for_z)) {
      float z = 
	_regularBarrelModuleLength-2*_lateralPlateThickness+(id*_hcalAbsorberThickness+(id-1)*_hcalSensitiveThickness-y_dim1_for_z)*(top_end_dim_z-_regularBarrelModuleLength)/y_dim2_for_z;
      z_width = int(z/_virtualCellSizeZ)*_virtualCellSizeZ;
    }

    if (z >= (y_dim1_for_z+y_dim2_for_z) ) {
      float z = top_end_dim_z - 2*_lateralPlateThickness;
      z_width = int(z/_virtualCellSizeZ)*_virtualCellSizeZ;
    }


    _endBarrelChamberLength[i] = z_width;
    _barrelLateralWidth[i] = x_width;
    int imax = int(x_width/_newCellSizeX);
    float offset = imax * _newCellSizeX;
    _barrelOffsetMaxX[i] = offset + 0.5*(x_width-offset);
    imax = int(z_width/_newCellSizeZ);
    offset = imax * _newCellSizeZ;
    _endBarrelOffsetMaxZ[i] = offset + 0.5*(z_width-offset);

  }

  int imax = int(_regularBarrelChamberLength/_newCellSizeZ);
  float offset = imax * _newCellSizeZ;
  _regularBarrelOffsetMaxZ = offset + 0.5*(_regularBarrelChamberLength-offset);

  _nModules = 7;
  _nStaves = 8;
  _deltaPhi = 2.*Pi/8.;

  float half_endcap_hole = float(pHcalEndcap.getExtent()[0]);
  _startIEndcap = int ((half_endcap_hole+0.00001) / _virtualCellSizeX);
  _startJEndcap = int ((half_endcap_hole+0.00001) / _virtualCellSizeZ);
  _startXEndcap = _startIEndcap * _virtualCellSizeX ;
  _startZEndcap = _startJEndcap * _virtualCellSizeZ ;

}

void MokkaCaloDigi::processRunHeader( LCRunHeader* run) { 

  _nRun++  ;
} 

void MokkaCaloDigi::processEvent( LCEvent * evt ) { 
   
  // set flag to store more information in the output file
  LCFlagImpl flag;
  flag.setBit(LCIO::CHBIT_LONG);

  //
  // * Reading Collections of HCAL Simulated Hits * 
  //
  _calorimeterHitVec.clear();
  int numberOfZones = _nModules*_nStaves*_numberOfHcalLayers;
  _calorimeterHitVec.resize(numberOfZones);
  _relationCollection = new LCCollectionVec(LCIO::LCRELATION);

  float simEnergy = 0;

  for (unsigned int i=0; i < _hcalCollections.size(); ++i) {
    try {
      LCCollection * col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      int numElements = col->getNumberOfElements();
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	simEnergy += hit->getEnergy();
	int cellid = hit->getCellID0();
	int Module=(cellid & MASK_M) >> SHIFT_M; // reed module number on it depends further calculation
	int Stave =(cellid & MASK_S) >> SHIFT_S; // stave	
	int Layer =(cellid & MASK_K) >> SHIFT_K; // layer  
	MyHit * newHit = NULL;
	if (Module == 0 || Module == 6) {
	  newHit = ProcessHitInEndcap( hit );
	}
	else {
	  newHit = ProcessHitInBarrel( hit );
	}
	if (newHit != NULL) {
	  newHit->hit->setType(2);
	  int sector = Layer + _numberOfHcalLayers*Stave + _numberOfHcalLayers*_nStaves*Module;
	  _calorimeterHitVec[sector].push_back( newHit );
	}

      }
      
    }
    catch(DataNotAvailableException &e){ }
  }

  float digitizedEnergy = 0.;
  LCCollectionVec * hcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  hcalcol->setFlag(flag.getFlag());
  for (int i=0; i<numberOfZones; ++i) {
    int nHits = int(_calorimeterHitVec[i].size());
    for (int iHit=0; iHit<nHits; ++iHit) {
      MyHit * myh = _calorimeterHitVec[i][iHit];
      CalorimeterHitImpl * calhit = myh->hit;
      
      float energy = calhit->getEnergy();
      digitizedEnergy += energy;
      if (energy > _thresholdHcal) {
	int Cellid = calhit->getCellID0();
	float calibr_coeff(1.);
	int layer = Cellid >> 24;
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
	if (_digitalHcal) {
	  calhit->setEnergy(calibr_coeff); 
	}
	else {
	  calhit->setEnergy(calibr_coeff*energy);
	}
	int nSimHit = myh->simHits.size();
	float x = calhit->getPosition()[0];
	float y = calhit->getPosition()[1];
	float z = calhit->getPosition()[2];

	for (int iS = 0; iS < nSimHit; ++iS) {
	  float weight = 1.0;
	  SimCalorimeterHit * simHit = myh->simHits[iS];
	  LCRelationImpl * rel = new LCRelationImpl(calhit, simHit , weight);
	  _relationCollection->addElement(rel);
	  float x1 =  simHit->getPosition()[0];
	  float y1 =  simHit->getPosition()[1];
	  float z1 =  simHit->getPosition()[2];
	  float dist = sqrt((x-x1)*(x-x1)+
			    (y-y1)*(y-y1)+
			    (z-z1)*(z-z1));
	  float cutX = 0.5*(_newCellSizeX-_virtualCellSizeX) ;
	  float cutZ = 0.5*(_newCellSizeZ-_virtualCellSizeZ) ;

	  float cut = sqrt(cutX*cutX+cutZ*cutZ) + 0.01;
	  if (dist > cut) {
	    std::cout << "WARNING ==> " << std::endl;
	    std::cout << "Distance = " << dist << " > " << cut << std::endl;
	    std::cout << x << " " << y << " " << z << std::endl;
	    std::cout << x1 << " " << y1 << " " << z1 << std::endl;
	  }
	}
	hcalcol->addElement( calhit );
      }
      else{
	// std::cout << " >>>>>>>> deleting spurious Calohit ! " << std::endl ;
	//fg: fix memory leak - if energy below threshold the hit is never added to the collection,
	//    thus we have to delete it !
	delete calhit ;
      }

      delete myh;
    }
  }

 //  std::cout << "Initial energy = " << simEnergy 
// 	    << " compared to digi energy =" << digitizedEnergy << std::endl;

  // 
  // * Reading Collections of ECAL Simulated Hits * 
  // 

  LCCollectionVec * ecalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  ecalcol->setFlag(flag.getFlag());

  for (unsigned int i(0); i < _ecalCollections.size(); ++i) {
      try{
	  LCCollection * col = evt->getCollection( _ecalCollections[i].c_str() ) ;
	  int numElements = col->getNumberOfElements();

	  for (int j(0); j < numElements; ++j) {
	      SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	     
                    
	      float energy = hit->getEnergy();
	      if (energy > _thresholdEcal) {
		CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
		int Cellid = hit->getCellID0();
		float calibr_coeff(1.);
		int layer = Cellid >> 24;
		int type = 0;
		for (unsigned int k(0); k < _ecalLayers.size(); ++k) {
		  int min,max;
		  if (k == 0) 
		    min = 0;		      
		  else 
		    min = _ecalLayers[k-1];		      
		  max = _ecalLayers[k];
		  if (layer >= min && layer < max) {
		    calibr_coeff = _calibrCoeffEcal[k];
		    type = k;
		    break;
		  }
		} 
		calhit->setCellID0(Cellid);
		if (_digitalEcal) {
		  calhit->setEnergy(calibr_coeff); 
		}
		else {
		  calhit->setEnergy(calibr_coeff*energy);
		}
		calhit->setType(type);
		calhit->setPosition(hit->getPosition());
		ecalcol->addElement(calhit);
		float weight = 1.0;
		LCRelationImpl * rel = new LCRelationImpl(calhit,hit,weight);
		_relationCollection->addElement( rel );
	      }
	  }
      }
      catch(DataNotAvailableException &e){ 
      }
  }

  evt->addCollection(ecalcol, _newCollNameECAL.c_str());
  evt->addCollection(hcalcol, _newCollNameHCAL.c_str());
  evt->addCollection(_relationCollection,_relationCollName.c_str());
  _nEvt ++ ;

} // end of event processor 

MyHit * MokkaCaloDigi::ProcessHitInBarrel( SimCalorimeterHit * hit ) {

  MyHit * newMyHit = NULL;

  int cellid = hit->getCellID0();
  float pos[3];
  for (int i=0;i<3;++i)
    pos[i] = hit->getPosition()[i];
  int Module=(cellid & MASK_M) >> SHIFT_M; // reed module number on it depends further calculation
  int Stave=(cellid & MASK_S) >> SHIFT_S; // stave	
  int Layer=(cellid & MASK_K) >> SHIFT_K; // layer      
  //int J=(cellid & MASK_J) >> SHIFT_J; // J 
  int I=(cellid & MASK_I) >> SHIFT_I; // I

  float zBegin = 0.;
  float chamberLength = 0.;
  float offsetMaxZ;

  // calculation of the lower z coordinate of the sensitive part of barrel module
  if (Module == 1) {
    zBegin = -1.5*_regularBarrelModuleLength - 2*_modulesGap - _lateralPlateThickness - _endBarrelChamberLength[Layer];
  }
  if (Module == 2) {
    zBegin = - 1.5*_regularBarrelModuleLength - _modulesGap + _lateralPlateThickness;
  }
  if (Module == 3) {
    zBegin = -0.5*_regularBarrelChamberLength;
  }
  if (Module == 4) {
    zBegin = 0.5*_regularBarrelModuleLength + _modulesGap + _lateralPlateThickness;
  }
  if (Module == 5) {
    zBegin = 1.5*_regularBarrelModuleLength + 2*_modulesGap + _lateralPlateThickness;
  }
  if (Module == 1 || Module == 5) {
    chamberLength = _endBarrelChamberLength[Layer];
    offsetMaxZ = _endBarrelOffsetMaxZ[Layer];
  }
  else {
    chamberLength = _regularBarrelChamberLength;
    offsetMaxZ = _regularBarrelOffsetMaxZ;
  }


  float xBegin = -0.5*_barrelLateralWidth[Layer];
  int Inew = I / _cellScaleX;
  int Jnew = int((pos[2] - zBegin)/ _newCellSizeZ);
  float offsetI = (Inew+0.5) * _newCellSizeX ;
  float offsetJ = (Jnew+0.5) * _newCellSizeZ ;

  int cellidNew=((Module<<SHIFT_M)&MASK_M)
    |    ((Stave<<SHIFT_S)&MASK_S)
    |    ((Jnew<<SHIFT_J)&MASK_J)
    |    ((Layer<<SHIFT_K)&MASK_K)
    |    ((Inew<<SHIFT_I)&MASK_I);

  int sector = Layer + _numberOfHcalLayers*Stave + _numberOfHcalLayers*_nStaves*Module;
  
  int iExist = 0;

  int SIZE = int(_calorimeterHitVec[sector].size());
  for (int i=0; i<SIZE; ++i) {
    MyHit * myh = _calorimeterHitVec[sector][i];
    CalorimeterHitImpl * hitImpl = myh->hit;
    int cellidImpl = hitImpl->getCellID0();
    if (cellidNew == cellidImpl) {
      iExist = 1;
      float energy = hitImpl->getEnergy() + hit->getEnergy();
      hitImpl->setEnergy( energy );
      myh->simHits.push_back( hit );
      break;
    }
  }

  if (iExist == 0) {
    CalorimeterHitImpl * newHit = new CalorimeterHitImpl();
    newHit->setCellID0(cellidNew);
    newHit->setEnergy(hit->getEnergy());
    float newPos[3];
    newPos[1] = pos[1];
    if (Stave > 0) {
      float Radius = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
      float Phi = atan2(pos[1],pos[0]);    
      Phi = Phi - Stave*_deltaPhi;
      newPos[1] = Radius*sin(Phi);      
    }
    newPos[0] = xBegin + offsetI ;
    newPos[2] = zBegin + offsetJ ;
    if (Stave > 0) {
      float Radius = sqrt(newPos[0]*newPos[0]+newPos[1]*newPos[1]);
      float Phi = atan2(newPos[1],newPos[0]);
      Phi = Phi + Stave*_deltaPhi;
      newPos[0] = Radius*cos(Phi);
      newPos[1] = Radius*sin(Phi);
    }
    newHit->setPosition( newPos );
    newMyHit = new MyHit();
    newMyHit->hit = newHit;
    newMyHit->simHits.push_back( hit );

//    float dist = sqrt((pos[0]-newPos[0])*(pos[0]-newPos[0])+
//		      (pos[1]-newPos[1])*(pos[1]-newPos[1])+
//		      (pos[2]-newPos[2])*(pos[2]-newPos[2]));
//     if (dist > 0.01) {
//       std::cout << "Layer width = " << _barrelLateralWidth[Layer] << std::endl;
//       std::cout << I << " " << Inew << std::endl;
//       float Radius = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
//       float Phi = atan2(pos[1],pos[0]);    
//       Phi = Phi - Stave*_deltaPhi;
//       float xpos = Radius*cos(Phi);
//       float xpos1 = xBegin + offsetI ;
//       std::cout << xpos << " " << xpos1 << std::endl;
//       std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
//       std::cout << newPos[0] << " " << newPos[1] << " " << newPos[2] << " " << std::endl;
//     }

  }


  return newMyHit ;

}

MyHit * MokkaCaloDigi::ProcessHitInEndcap(SimCalorimeterHit * hit) {

  MyHit * newMyHit = NULL;

  int cellid = hit->getCellID0();
  float pos[3];
  for (int i=0;i<3;++i)
    pos[i] = hit->getPosition()[i];
  int Module=(cellid & MASK_M) >> SHIFT_M; // reed module number on it depends further calculation
  int Stave=(cellid & MASK_S) >> SHIFT_S; // stave	
  int Layer=(cellid & MASK_K) >> SHIFT_K; // layer      
  int J=(cellid & MASK_J) >> SHIFT_J; // J 
  int I=(cellid & MASK_I) >> SHIFT_I; // I
  
  int Inew = I / _cellScaleX;
  int Jnew = J / _cellScaleZ;

  int cellidNew=((Module<<SHIFT_M)&MASK_M)
    |    ((Stave<<SHIFT_S)&MASK_S)
    |    ((Jnew<<SHIFT_J)&MASK_J)
    |    ((Layer<<SHIFT_K)&MASK_K)
    |    ((Inew<<SHIFT_I)&MASK_I);

  int sector = Layer + _numberOfHcalLayers*Stave + _numberOfHcalLayers*_nStaves*Module;
  
  int iExist = 0;

  int SIZE = int(_calorimeterHitVec[sector].size());
  for (int i=0; i<SIZE; ++i) {
    MyHit * myh = _calorimeterHitVec[sector][i];
    CalorimeterHitImpl * hitImpl = myh->hit;
    int cellidImpl = hitImpl->getCellID0();
    if (cellidNew == cellidImpl) {
      iExist = 1;
      float energy = hitImpl->getEnergy() + hit->getEnergy();
      hitImpl->setEnergy( energy );
      myh->simHits.push_back( hit );
      break;
    }
  }

  if (iExist == 0) {
    CalorimeterHitImpl * newHit = new CalorimeterHitImpl();
    newHit->setCellID0(cellidNew);
    newHit->setEnergy(hit->getEnergy());
    float newPos[3];
    float offsetX = (Inew + 0.5)*_newCellSizeX; 
    float offsetZ = (Jnew + 0.5)*_newCellSizeZ;
    if (Module == 0) {
      if (Stave == 0) {
	newPos[0] = - offsetZ;
	newPos[1] = offsetX;
      }
      else if (Stave == 1) {
	newPos[0] = -offsetX;
	newPos[1] = -offsetZ;
      }
      else if (Stave == 2) {
	newPos[0] = offsetZ;
	newPos[1] = -offsetX;
      }
      else if (Stave == 3) {
	newPos[0] = offsetX;
	newPos[1] = offsetZ;
      }
    }
    else {
      if (Stave == 0) {
	newPos[0] = offsetZ;
	newPos[1] = offsetX;
      }
      else if (Stave == 1) {
	newPos[0] = offsetX;
	newPos[1] = -offsetZ;
      }
      else if (Stave == 2) {
	newPos[0] = -offsetZ;
	newPos[1] = -offsetX;
      }
      else if (Stave == 3) {
	newPos[0] = -offsetX;
	newPos[1] = offsetZ;
      }
    }
    newPos[2] = pos[2] ;    
    newHit->setPosition( newPos );
    newMyHit = new MyHit();
    newMyHit->hit = newHit;
    newMyHit->simHits.push_back( hit );
  }
  return newMyHit;
}


void MokkaCaloDigi::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void MokkaCaloDigi::end(){ 
  
  std::cout << "MokkaCaloDigi::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

}

