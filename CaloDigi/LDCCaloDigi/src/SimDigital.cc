#include "SimDigital.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/BitField64.h>

#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "CalorimeterHitType.h"


#include <TROOT.h>
#include <TMath.h>
#include "TTree.h"
#include "TH1F.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;



SimDigital aSimDigital ;

SimDigital::SimDigital () : Processor("SimDigital"), _Ranm(0) {
   _description = "the transfer Between Energy and Induced Charge for SDHCAL" ;
 
   // the realtion should be with ECAL

   std::vector<std::string> ecalCollections;
   ecalCollections.push_back(std::string("EcalBarrelCollection"));
   ecalCollections.push_back(std::string("EcalEndcapCollection"));
   ecalCollections.push_back(std::string("EcalRingCollection"));
   ecalCollections.push_back(std::string("EcalBarrelPreShowerCollection"));
   ecalCollections.push_back(std::string("EcalEndcapPreShowerCollection"));
   ecalCollections.push_back(std::string("EcalRingPreShowerCollection"));
   registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			    "ECALCollections" , 
			    "ECAL Collection Names" ,
			    _ecalCollections ,
			    ecalCollections);
   
   ///////////////////////////////////////////////////////////
   std::vector<std::string> hcalCollections;
   // CaloHitCollections.push_back(std::string("Dummy"));
  
   hcalCollections.push_back(std::string("HcalBarrelCollection"));
   hcalCollections.push_back(std::string("HcalEndCapRingsCollection"));
   hcalCollections.push_back(std::string("HcalEndCapsCollection"));

   registerInputCollections( LCIO::SIMCALORIMETERHIT,
                            "HCALCollections" ,
                            "Sim Calorimeter Hit Collections" ,
                            _hcalCollections ,
                             hcalCollections);
   /////////ECAL OUTPUT

   _outputEcalCollections.push_back(std::string("ECALBarrel"));
   _outputEcalCollections.push_back(std::string("ECALEndcap"));
   _outputEcalCollections.push_back(std::string("ECALOther"));
   _outputEcalCollections.push_back(std::string("ECALBarrelPreShower"));
   _outputEcalCollections.push_back(std::string("ECALEndcapPreShower"));
   _outputEcalCollections.push_back(std::string("ECALOtherPreShower"));

   ///////////////////////////////////////////////////////////


   
   _outputHcalCollections.push_back(std::string("HCALBarrel"));
   _outputHcalCollections.push_back(std::string("HCALEndcap"));
   _outputHcalCollections.push_back(std::string("HCALOther"));

//    registerProcessorParameter("SDHCALThreshold" , 
// 			     "Threshold for SDHCAL Hits in GeV" ,
// 			     _thresholdSDHcal,
// 			     (float)2.0e-10);
   std::vector<float> hcalThresholds;
   hcalThresholds.push_back(0.00004);
   registerProcessorParameter("HCALThreshold" ,
                             "Threshold for HCAL Hits in pC" ,
			      _thresholdHcal,
			      hcalThresholds);
   
   std::vector<float> calibrHcal;
   calibrHcal.push_back(34.8);
   //////////////////ECAL//////////////////
registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "ECALOutputCollection0" , 
			    "ECAL Collection of real Hits" , 
			    _outputEcalCollections[0], 
			    std::string("ECALBarrel") ); 
  

  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "ECALOutputCollection1" , 
			    "ECAL Collection of real Hits" , 
			     _outputEcalCollections[1], 
			    std::string("ECALEndcap") ); 
  
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "ECALOutputCollection2" , 
			    "ECAL Collection of real Hits" , 
			    _outputEcalCollections[2], 
			    std::string("ECALOther") ) ; 

  registerOutputCollection( LCIO::CALORIMETERHIT,
                            "ECALOutputCollection3" ,
                            "ECAL Collection of real Hits" ,
                            _outputEcalCollections[3],
                            std::string("ECALBarrelPreShower") );


  registerOutputCollection( LCIO::CALORIMETERHIT,
                            "ECALOutputCollection4" ,
                            "ECAL Collection of real Hits" ,
                             _outputEcalCollections[4],
                            std::string("ECALEndcapPreShower") );

  registerOutputCollection( LCIO::CALORIMETERHIT,
                            "ECALOutputCollection5" ,
                            "ECAL Collection of real Hits" ,
                            _outputEcalCollections[5],
                            std::string("ECALOtherPreShower") ) ;
  
  registerProcessorParameter("ECALThreshold" , 
			     "Threshold for ECAL Hits in GeV" ,
			     _thresholdEcal,
			     (float)5.0e-5);

  std::vector<int> ecalLayers;
  ecalLayers.push_back(20);
  ecalLayers.push_back(100);


  registerProcessorParameter("ECALLayers" , 
			     "Index of ECal Layers" ,
			     _ecalLayers,
			     ecalLayers);

   std::vector<float> calibrEcal;
  calibrEcal.push_back(40.91);
  calibrEcal.push_back(81.81);


  registerProcessorParameter("CalibrECAL" , 
			     "Calibration coefficients for ECAL" ,
			     _calibrCoeffEcal,
			     calibrEcal);

   registerProcessorParameter("IfDigitalEcal" ,
			     "Digital Ecal" , 
			     _digitalEcal , 
			     0);

 registerProcessorParameter("ECALGapCorrection" , 
			     "Correct for ECAL gaps" ,
			     _ecalGapCorrection,
			     (int)1);

  registerProcessorParameter("ECAlEndcapCorrectionFactor" , 
			     "Energy correction for endcap" ,
			     _ecalEndcapCorrectionFactor,
			     (float)1.025);

  registerProcessorParameter("ECALGapCorrectionFactor" , 
			     "Factor applied to gap correction" ,
			     _ecalGapCorrectionFactor,
			     (float)1.0);

  registerProcessorParameter("ECALModuleGapCorrectionFactor" , 
			     "Factor applied to module gap correction" ,
			     _ecalModuleGapCorrectionFactor,
			     (float)0.5);
//////////////////////////////////////////////////////////
//   std::vector<int> hcalLayers;
//   hcalLayers.push_back(100);
//   registerProcessorParameter("HCALLayers" , 
// 			     "Index of HCal Layers" ,
// 			     _hcalLayers,
// 			     hcalLayers);


   registerProcessorParameter("CalibrHCAL" , 
			     "Calibration coefficients for HCAL" ,
			      _calibrCoeffHcal,
			      calibrHcal);
//    registerProcessorParameter("IfDigitalHcal" ,
// 			      "Digital Hcal" , 
// 			      _digitalHcal , 
// 			      1);

   registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "HCALOutputCollection0" , 
			    "HCAL Collection of real Hits" , 
			    _outputHcalCollections[0], 
			    std::string("HCALBarrel")  ); 
  
   registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "HCALOutputCollection1" , 
			    "HCAL Collection of real Hits" , 
			    _outputHcalCollections[1], 
			    std::string("HCALEndcap") );

   registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "HCALOutputCollection2" , 
			    "HCAL Collection of real Hits" , 
			    _outputHcalCollections[2], 
			    std::string("HCALOther") ) ; 

   registerOutputCollection( LCIO::LCRELATION, 
			    "RelationOutputCollection" , 
			    "CaloHit Relation Collection" , 
			     _outputRelCollection , 
			    std::string("RelationCaloHit")) ;
   bool defaultPrinting=false;
   _printSimDigital=false;
   registerOptionalParameter("printSimDigital",
			     "Print or not SimDigital",
			     _printSimDigital,
			     defaultPrinting);
   
   registerProcessorParameter( "CellEdgeDistance" ,
                               "Distance (in mm) from the border of the cell that implies hit on neighbouring cell"  ,
                               _edgedistance,
                               (float)2.0 ) ;

   _doThresholds=true;
   registerOptionalParameter("doThresholds",
			     "Replace analog hit energy by value given in CalibrHCAL according to thresholds given in HCALThreshold",
			     _doThresholds,
			     true);

   _relcol=0;
   _currentHCALCollectionCaloLayout=CHT::any;
}
 //define a Polya Function 

  Double_t Polya(double* x, double* params)
  {
    // Polya function.
    //0 = Polya - average charge
    // 1 = Polya - free parameter
    //  2 = Polya - amplitude
    // Parameters:
    // 0 = average charge (gain)
    // 1 = free parameter determining the shape of the distribution
    Double_t q, t;
    if(*x <= 0)
      return 0;
    q = *x / params[0];
    t = 1.0 + params[1];
    return params[2]*(TMath::Power(t, t) 
		      / TMath::Gamma(t)
		      * TMath::Power(q, params[1])
		      * TMath::Exp(-t*q));
  } 
  // for integral to get charge distribution 
  Double_t ChaDis(Double_t *x, Double_t *par)
  {

    float a=3.66*TMath::Pi();
    float xo=x[0]-par[1];
    float yo=x[1]-par[2];
    float r=sqrt(xo*xo+yo*yo);
    float c=(par[0]*par[0])/a;
    return c/(TMath::CosH(r*par[0]));
    
  }
//

void SimDigital::init(){
  const float pi = acos(-1.0);
  const float twopi = 2.0*pi;

  //a=1.05,b=-2.1,c=1.05,d=-1.05;//for 1 neighbor pad fired
  //previousLayer=1 ,previousmodule=1, previousstave=1,previousIy=1,previousJz=1;
   //define a polya function to get induced charge
   _QPolya = new TF1("QPolya",Polya,0.0,10.0,3);
   _QPolya->SetParameters(1.6,16.3,1);
   //define the function for charge distribution
   //c2d=new TF2("c2d",ChaDis,-10,10,-10,10,3);
   //c2d->SetParNames("width","x'value","y'value");
  

   //fg: need to set default encoding in for reading old files...
   CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;

  // Calorimeter geometry from GEAR
   const gear::CalorimeterParameters& pEcalBarrel = Global::GEAR->getEcalBarrelParameters();
   const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
   //  const gear::CalorimeterParameters& pHcalBarrel = Global::GEAR->getHcalBarrelParameters();
  //  const gear::CalorimeterParameters& pHcalEndcap = Global::GEAR->getHcalEndcapParameters();
   const gear::LayerLayout& ecalBarrelLayout = pEcalBarrel.getLayerLayout();
   const gear::LayerLayout& ecalEndcapLayout = pEcalEndcap.getLayerLayout();
   // const gear::LayerLayout& hcalBarrelLayout = pHcalBarrel.getLayerLayout();
   // const gear::LayerLayout& hcalEndcapLayout = pHcalEndcap.getLayerLayout();


   // determine geometry of ECAL
   int symmetry = pEcalBarrel.getSymmetryOrder();
   _zOfEcalEndcap = (float)pEcalEndcap.getExtent()[2];

   // Determine ECAL polygon angles
   // Store radial vectors perpendicular to stave layers in _ecalBarrelStaveDir 
   // ASSUMES Mokka Stave numbering 0 = top, then numbering increases anti-clockwise
   if(symmetry>1){
     float nFoldSymmetry = static_cast<float>(symmetry);
     float phi0 = pEcalBarrel.getPhi0();
     for(int i=0;i<symmetry;++i){
       float phi  = phi0 + i*twopi/nFoldSymmetry;
       _barrelStaveDir[i][0] = cos(phi);
       _barrelStaveDir[i][1] = sin(phi);

       //	cout<<i<<"\t"<<symmetry<<endl;
     }
   }  


   for(int i=0;i<ecalBarrelLayout.getNLayers();++i){
     _barrelPixelSizeT[i] = ecalBarrelLayout.getCellSize0(i);
     _barrelPixelSizeZ[i] = ecalBarrelLayout.getCellSize1(i);
   }

   for(int i=0;i<ecalEndcapLayout.getNLayers();++i){
     _endcapPixelSizeX[i] = ecalEndcapLayout.getCellSize0(i);
     _endcapPixelSizeY[i] = ecalEndcapLayout.getCellSize1(i);
   }


   //assure SDHCAL _thresholdHcal are in increasing order
   std::sort(_thresholdHcal.begin(),_thresholdHcal.end());
}
 

void SimDigital::processRunHeader( LCRunHeader* run) {}

void SimDigital::processECAL(LCEvent* evt,LCFlagImpl& flag)
{
  // 
  // * Reading Collections of ECAL Simulated Hits * 
  // 
  for (unsigned int i(0); i < _ecalCollections.size(); ++i) {

    std::string colName =  _ecalCollections[i] ;

    //fg: need to establish the subdetetcor part here 
    //    use collection name as cellID does not seem to have that information
    CHT::Layout caloLayout = layoutFromString( colName ) ;

    try{
      LCCollection * col = evt->getCollection( _ecalCollections[i].c_str() ) ;
      string initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);

      CellIDDecoder<SimCalorimeterHit> idDecoder( col );

      // create new collection
      LCCollectionVec *ecalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
      ecalcol->setFlag(flag.getFlag());

      // if making gap corrections clear the vectors holding pointers to calhits
      if(_ecalGapCorrection!=0){
	for(int is=0;is<MAX_STAVES;is++){
	  for(int il=0;il<MAX_LAYERS;il++){	
	    _calHitsByStaveLayer[is][il].clear();
	    _calHitsByStaveLayerModule[is][il].clear();
	  }
	}
      }

      int numElements = col->getNumberOfElements();
      for (int j(0); j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	float energy = hit->getEnergy();
	// apply threshold cut
	if (energy > _thresholdEcal) {
	  CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
	  int cellid = hit->getCellID0();
	  int cellid1 = hit->getCellID1();
	  float calibr_coeff(1.);
	  int layer = idDecoder(hit)["K-1"];
	  int stave = idDecoder(hit)["S-1"];
	  int module= idDecoder(hit)["M"];
	  // save hits by module/stave/layer if required later
	  if(_ecalGapCorrection!=0){
	    _calHitsByStaveLayer[stave][layer].push_back(calhit);
	    _calHitsByStaveLayerModule[stave][layer].push_back(module);
	  }

	  // retrieve calibration constants
	  for (unsigned int k(0); k < _ecalLayers.size(); ++k) {
	    int min,max;
	    if (k == 0){ 
	      min = 0;		      
	    }else{ 
	      min = _ecalLayers[k-1];
	    } 
	    max = _ecalLayers[k];
	    if (layer >= min && layer < max) {
	      calibr_coeff = _calibrCoeffEcal[k];
	      break;
	    }
	  } 
	  // apply calibration
	  if (_digitalEcal) {
	    calhit->setEnergy(calibr_coeff); 
	  }
	  else {
	    // if in endcap apply additional factor to calibration to account for
	    // the difference in response due to the orientation of B wrt absorber
	    if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)energy=energy*_ecalEndcapCorrectionFactor;
	    calhit->setEnergy(calibr_coeff*energy);
	  }
	  // set other ECAL quanties
	  calhit->setPosition(hit->getPosition());

	  calhit->setType( CHT( CHT::em, CHT::ecal, caloLayout ,  layer ) );

	  calhit->setRawHit(hit);
	  calhit->setCellID0(cellid);
	  calhit->setCellID1(cellid1);
	  ecalcol->addElement(calhit);
	  // make relation between hit and sim hit
	  LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.);
	  _relcol->addElement( rel );
	}
      }
      // if requested apply gap corrections in ECAL ? 
      if(_ecalGapCorrection!=0)this->fillECALGaps();
      // add ECAL collection to event
      ecalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
      evt->addCollection(ecalcol,_outputEcalCollections[i].c_str());      
    }
    catch(DataNotAvailableException &e){ 
    }
  
  }
}

void SimDigital::processHCAL(LCEvent* evt, LCFlagImpl& flag)
{
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    std::string colName =  _hcalCollections[i] ;    
    _currentHCALCollectionCaloLayout = layoutFromString( colName ) ;
    try{
      LCCollection * col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      LCCollectionVec *hcalcol = processHCALCollection(col,flag);
      evt->addCollection(hcalcol,_outputHcalCollections[i].c_str());
    }
    catch(DataNotAvailableException &e){  
    }     
  }
}



SimDigital::multiplicityChargeSplitter::multiplicityChargeSplitter(float cell_Size, float edge_Size)
{
  cellSize=cell_Size;
  edgeSize=edge_Size;
  chargeMap[LeftRight_LowerUpper_Coordinates(Central,Central)]=0;
  
}

void SimDigital::multiplicityChargeSplitter::addCharge(float charge, float posLeftRight, float posLowerUpper)
{
  //It is assumed that within a cell coordinate 0,0 is the Lower Left Corner
  std::set<LeftRight_LowerUpper_Coordinates> ModifyCharge;
  ModifyCharge.insert(LeftRight_LowerUpper_Coordinates(Central,Central));
  DIRECTION LeftRight=Central,LowerUpper=Central;
  if (edgeSize<=cellSize)
    {
      if (posLeftRight<edgeSize) LeftRight=Left;
      if (posLeftRight>cellSize-edgeSize) LeftRight=Right;
      if (posLowerUpper<edgeSize) LowerUpper=Lower;
      if (posLowerUpper>cellSize-edgeSize) LowerUpper=Upper;
      ModifyCharge.insert(LeftRight_LowerUpper_Coordinates(Central,LowerUpper));
      ModifyCharge.insert(LeftRight_LowerUpper_Coordinates(LeftRight,Central));
      ModifyCharge.insert(LeftRight_LowerUpper_Coordinates(LeftRight,LowerUpper));
    }
  else
    {
      int Ncells=int(edgeSize/cellSize);
      for (int iLR=-Ncells; iLR<=Ncells; iLR++)
	for (int iLU=-Ncells; iLU<=Ncells; iLU++)
	  ModifyCharge.insert(LeftRight_LowerUpper_Coordinates(iLR,iLU));
    }
  float SplitCharge=charge/ModifyCharge.size();
  for (std::set<LeftRight_LowerUpper_Coordinates>::iterator it=ModifyCharge.begin(); 
       it != ModifyCharge.end(); it++)
    chargeMap[*it]+=SplitCharge;
}


const gear::LayerLayout& SimDigital::getGearLayout()
{
  if (_currentHCALCollectionCaloLayout == CHT::endcap) return Global::GEAR->getHcalEndcapParameters().getLayerLayout();
  if (_currentHCALCollectionCaloLayout == CHT::ring) return Global::GEAR->getHcalRingParameters().getLayerLayout();
  return Global::GEAR->getHcalBarrelParameters().getLayerLayout();
}

void SimDigital::createPotentialOutputHits(std::map<int,hitMemory>& myHitMap, LCCollection *col, CellIDEncoder<CalorimeterHitImpl>& encodid )
{
  const gear::LayerLayout& hcalLayout=getGearLayout();

  CellIDDecoder<SimCalorimeterHit> idDecoder(col);

  int numElements = col->getNumberOfElements();
  for (int j(0); j < numElements; ++j) {
    SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
    //for moudle
    int trueLayer=-999;int stave=-999;int module=-999;int Iy=-999;int Jz=-999;
    trueLayer = idDecoder( hit )["K-1"]; // + 1;
    stave     = idDecoder( hit )["S-1"];
    module    = idDecoder( hit )["M"];
    Iy        = idDecoder( hit )["I"];
    Jz        = idDecoder( hit )["J"];
	
    const float * pos=hit->getPosition();
    int nmcparticles= hit->getNMCContributions();
    float cellSize=hcalLayout.getCellSize0(trueLayer);
    multiplicityChargeSplitter mult(cellSize,_edgedistance); 
    for(int imcp=0;imcp<nmcparticles;imcp++){
      float Minduced=_QPolya->GetRandom();
      if(Minduced<0.4)streamlog_out(DEBUG) <<" "<<Minduced<<std::endl;
      
      //SimCalorimeterHit::getStepPostion(int i);
      //const float * steppos = hit->getStepPostion(j);
      //do stuff to get random draw
      float x= _Ranm.Uniform(cellSize);
      float y= _Ranm.Uniform(cellSize);
   
      mult.addCharge(Minduced,x,y);
    }
    
    
    for (std::map<multiplicityChargeSplitter::LeftRight_LowerUpper_Coordinates,float>::iterator it=mult.chargeMap.begin();
	 it != mult.chargeMap.end(); it++)
      {
	if (it->second > 0)
	  {
	    float posB[3];
	    posB[0]=pos[0]+it->first.first*cellSize;
	    posB[1]=pos[1]+it->first.second*cellSize;
	    posB[2]=pos[2];
	    CalorimeterHitImpl *tmp=new CalorimeterHitImpl();
	    int RealIy=Iy+it->first.first;
	    int RealJz=Jz+it->first.second;
	    if(RealIy<0||RealIy>330)RealIy=0; //FIXME the 330 value should depend on the cellSize and on the Layer
	    if(RealJz<0||RealJz>330)RealJz=0; //FIXME the 330 value should depend on the cellSize and on the Layer
	    encodid["M"]=module;
	    encodid["S-1"]=stave;
	    encodid["I"]=RealIy;
	    encodid["J"]=RealJz;
	    encodid["K-1"]=trueLayer;
	    encodid.setCellID( tmp );
	    hitMemory &calhitMem=myHitMap[tmp->getCellID0()];	
	    if (calhitMem.ahit==0){
	      calhitMem.ahit=tmp;
	      calhitMem.ahit->setCellID1(hit->getCellID1());
	      calhitMem.ahit->setPosition(posB);
	      calhitMem.ahit->setType( CHT( CHT::had, CHT::hcal , _currentHCALCollectionCaloLayout ,  trueLayer ) );
	      calhitMem.ahit->setEnergy(0); 
	    }
	    else
	      delete tmp;
	    calhitMem.rawHit=j; //for (int j(0); j < numElements; ++j)   //FIXME : rawhit is the last one, would be better if the one contibuting most to the hit energy.
	    calhitMem.ahit->setEnergy(calhitMem.ahit->getEnergy()+it->second);
	    calhitMem.relatedHits.insert(j); //for (int j(0); j < numElements; ++j)
	  }	
	else
	  {
	    streamlog_out(ERROR) << "BUG in charge splitter, got a non positive charge : " << it->second << std::endl;
	  }
      }
	
  } // end of for (int j(0); j < numElements; ++j)  //loop on elements in collection

}


void SimDigital::removeHitsBelowThreshold(std::map<int,hitMemory>& myHitMap, float threshold)
{
  std::map<int,hitMemory>::iterator itr;
  ThresholdIsBelow t(threshold);
  do {
    itr=find_if(myHitMap.begin(), myHitMap.end(), t);
    if (itr != myHitMap.end())
      {
	delete itr->second.ahit;
	myHitMap.erase(itr);
      }
  } while (itr != myHitMap.end());
}


void SimDigital::applyThresholds(std::map<int,hitMemory>& myHitMap)
{
  for (std::map<int,hitMemory>::iterator it=myHitMap.begin(); it!=myHitMap.end(); it++) {
    hitMemory & currentHitMem=it->second;
    float Tcharge=currentHitMem.ahit->getEnergy();
    float calibr_coeff(1.);
    unsigned int ilevel=0;
    for(unsigned int ithresh=1;ithresh<_thresholdHcal.size();ithresh++){ 
      if(Tcharge>_thresholdHcal[ithresh])ilevel=ithresh;   // ilevel = 0 , 1, 2
    }
    if(ilevel>_calibrCoeffHcal.size()-1){
      streamlog_out(ERROR)  << " Semi-digital level " << ilevel  << " greater than number of HCAL Calibration Constants (" <<_calibrCoeffHcal.size() << ")" << std::endl;
      ilevel=_calibrCoeffHcal.size()-1;
    }else{
      calibr_coeff = _calibrCoeffHcal[ilevel];
    }
    currentHitMem.ahit->setEnergy(calibr_coeff);
  }
}

LCCollectionVec * SimDigital::processHCALCollection(LCCollection * col, LCFlagImpl& flag)
{


  string initString = col->getParameters().getStringVal(LCIO::CellIDEncoding); 
  LCCollectionVec *hcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  hcalcol->setFlag(flag.getFlag());
  CellIDEncoder<CalorimeterHitImpl> encodid(initString,hcalcol);//for adding new code id
  std::map<int,hitMemory> myHitMap;

  createPotentialOutputHits(myHitMap,col,encodid);
  removeHitsBelowThreshold(myHitMap,_thresholdHcal[0]);
  if (_doThresholds) applyThresholds(myHitMap); 
	  
  //Store element to output collection
  for (std::map<int,hitMemory>::iterator it=myHitMap.begin(); it!=myHitMap.end(); it++) {
    hitMemory & currentHitMem=it->second;
    if (currentHitMem.rawHit != -1) 
      {
	// streamlog_out(DEBUG)  << " rawHit= " << currentHitMem.rawHit << std::endl;
	SimCalorimeterHit * hitraw = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( currentHitMem.rawHit ) ) ;
	currentHitMem.ahit->setRawHit(hitraw);
      }
    hcalcol->addElement(currentHitMem.ahit);
    for (std::set<int>::iterator itset=currentHitMem.relatedHits.begin(); itset != currentHitMem.relatedHits.end(); itset++)
      {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( *itset ) ) ;
	LCRelationImpl *rel = new LCRelationImpl(currentHitMem.ahit,hit,1.0);
	_relcol->addElement( rel );
      }
  } //end of loop on myHitMap 
      
  // add HCAL collection to event
  hcalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
  return hcalcol;
}


void SimDigital::processEvent( LCEvent * evt ) {

  _counters["|ALL"]++;

   // create the output collections
  _relcol = new LCCollectionVec(LCIO::LCRELATION); 
 
  /////////////////for ECAL---------------------------------------------------
  // copy the flags from the input collection
  //GG : it should be checked why we put the flag like this.
  LCFlagImpl flag;
  flag.setBit(LCIO::CHBIT_LONG);
  flag.setBit(LCIO::RCHBIT_ENERGY_ERROR);    //open the energy error flag to store the MC Truth (for easy comparison == not a eligent way!!)
  processECAL(evt,flag);
  processHCAL(evt,flag);

 
  evt->addCollection(_relcol,_outputRelCollection.c_str());
}

void SimDigital::check( LCEvent * evt ) { }

void SimDigital::end(){
 
}


void SimDigital::fillECALGaps( ) { 
  const float slop = 0.25; // (mm)

  // Loop over hits in the Barrel
  // For each layer calculated differences in hit positions
  // Look for gaps based on expected separation of adjacent hits
  // loop over staves and layers

  for (int is=0; is < MAX_STAVES; ++is) {
    for (int il=0; il < MAX_LAYERS; ++il) {
      if(_calHitsByStaveLayer[is][il].size()>1){
	// compare all pairs of hits just once (j>i)

	for (unsigned int i=0;i<_calHitsByStaveLayer[is][il].size()-1;++i){
	  CalorimeterHitImpl* hiti = _calHitsByStaveLayer[is][il][i]; 
	  int modulei = _calHitsByStaveLayerModule[is][il][i];
	  float xi = hiti->getPosition()[0];
	  float yi = hiti->getPosition()[1];
	  float zi = hiti->getPosition()[2];

	  for (unsigned int j=i+1;j<_calHitsByStaveLayer[is][il].size();++j){
	    CalorimeterHitImpl* hitj = _calHitsByStaveLayer[is][il][j]; 
	    int modulej = _calHitsByStaveLayerModule[is][il][j];
	    float xj = hitj->getPosition()[0];
	    float yj = hitj->getPosition()[1];
	    float zj = hitj->getPosition()[2];
	    float dz = fabs(zi-zj);
	    // *** BARREL CORRECTION ***
	    if( fabs(zi)<_zOfEcalEndcap && fabs(zj)<_zOfEcalEndcap){
	      // account for stave directions using normals
	    // calculate difference in hit postions in z and along stave
	      float dx = xi-xj;
	      float dy = yi-yj;
	      float dt = fabs(dx*_barrelStaveDir[is][0] + dy*_barrelStaveDir[is][1]);
	      // flags for evidence for gaps
	      bool zgap = false;   // in z direction
	      bool tgap = false;   // along stave 
	      bool ztgap = false;  // in both z and along stave 
	      bool mgap = false;   // gaps between ECAL modules
	      
	      // criteria gaps in the z and t direction
	      float zminm = 1.0*_barrelPixelSizeZ[il]-slop;
	      float zmin = 1.0*_barrelPixelSizeZ[il]+slop;
	      float zmax = 2.0*_barrelPixelSizeZ[il]-slop;
	      float tminm = 1.0*_barrelPixelSizeT[il]-slop;
	      float tmin = 1.0*_barrelPixelSizeT[il]+slop;
	      float tmax = 2.0*_barrelPixelSizeT[il]-slop;
	      
	      // criteria for gaps
	      // WOULD BE BETTER TO USE GEAR TO CHECK GAPS ARE OF EXPECTED SIZE
	      if( dz > zmin  && dz < zmax && dt < tminm )zgap = true;
	      if( dz < zminm && dt > tmin && dt < tmax )tgap = true;
	      if( dz > zmin && dz < zmax && dt > tmin && dt < tmax )ztgap=true;

	      if(modulei!=modulej){
		if( dz > zmin && dz < 3.0*_barrelPixelSizeZ[il]-slop && dt < tmin)mgap = true;
	      }
 



	      // found a gap now apply a correction based on area of gap/area of pixel
	      if(zgap||tgap||ztgap||mgap){
		float ecor = 1.;
		float f = _ecalGapCorrectionFactor; // fudge
		if(mgap)f = _ecalModuleGapCorrectionFactor;
		if(zgap||mgap)ecor = 1.+f*(dz - _barrelPixelSizeZ[il])/2./_barrelPixelSizeZ[il];
		if(tgap)ecor = 1.+f*(dt - _barrelPixelSizeT[il])/2./_barrelPixelSizeT[il];
		if(ztgap)ecor= 1.+f*(dt - _barrelPixelSizeT[il])*(dz - _barrelPixelSizeZ[il])/4./_barrelPixelSizeT[il]/_barrelPixelSizeZ[il];     
		float ei = hiti->getEnergy()*ecor;
		float ej = hitj->getEnergy()*ecor;
		hiti->setEnergy(ei);
		hitj->setEnergy(ej);
	      }
	      
	    // *** ENDCAP CORRECTION ***
	    }else if(fabs(zi)>_zOfEcalEndcap && fabs(zj)>_zOfEcalEndcap&&dz<100){
	      float dx = fabs(xi-xj);
	      float dy = fabs(yi-yj);
	      bool xgap = false;
	      bool ygap = false;
	      bool xygap = false;
	      // criteria gaps in the z and t direction
	      float xmin = 1.0*_endcapPixelSizeX[il]+slop;
	      float xminm = 1.0*_endcapPixelSizeX[il]-slop;
	      float xmax = 2.0*_endcapPixelSizeX[il]-slop;
	      float ymin = 1.0*_endcapPixelSizeY[il]+slop;
	      float yminm = 1.0*_endcapPixelSizeY[il]-slop;
	      float ymax = 2.0*_endcapPixelSizeY[il]-slop;
	      // look for gaps
	      if(dx > xmin && dx < xmax && dy < yminm )xgap = true;
	      if(dx < xminm && dy > ymin && dy < ymax )ygap = true;
	      if(dx > xmin && dx < xmax && dy > ymin && dy < ymax )xygap=true;
	    
	      if(xgap||ygap||xygap){
		// found a gap make correction
		float ecor = 1.;
		float f = _ecalGapCorrectionFactor; // fudge
		if(xgap)ecor = 1.+f*(dx - _endcapPixelSizeX[il])/2./_endcapPixelSizeX[il];
		if(ygap)ecor = 1.+f*(dy - _endcapPixelSizeY[il])/2./_endcapPixelSizeY[il];
		if(xygap)ecor= 1.+f*(dx - _endcapPixelSizeX[il])*(dy - _endcapPixelSizeY[il])/4./_endcapPixelSizeX[il]/_endcapPixelSizeY[il];     
		hiti->setEnergy( hiti->getEnergy()*ecor );
		hitj->setEnergy( hitj->getEnergy()*ecor );
	      }
	    }
	  }
	}
      }
    }
  }

  return;


  
}



