#include "SimDigital.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <marlin/Global.h>
#include <marlin/Exceptions.h>
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
#include <fstream>
#include <time.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "CalorimeterHitType.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <marlin/AIDAProcessor.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>



#include <TROOT.h>
#include <TMath.h>
#include "TTree.h"
#include "TH1F.h"
#include "TRandom.h"

#include <DDSegmentation/BitField64.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/DetectorSelector.h>

using namespace lcio ;
using namespace marlin ;
using namespace std;

//FIXME to be removed when not more needed
//#define SDHCAL_MARLINUTIL_BUGFIX 1
//#define SDHCAL_MOKKA_BUGFIX 1

std::string SimDigitalGeomCellId::_encodingStrings[ENCODINGTYPES][ENCODINGSTRINGLENGTH] = 
{ 
	// The encoding string for lcgeo: barrel and endcap ring of hcal ( change "y" to "z" for hcal endcap )
	{ "layer", "stave", "module", "tower", "x", "y" },

    // The encoding string for Mokka
	{ "K-1", "S-1", "M", "", "I", "J" }
};

std::string SimDigitalGeomCellId::_hcalOption;


SimDigital aSimDigital ;


SimDigital::SimDigital () : Processor("SimDigital"), _chargeSplitterUniform(), _chargeSplitterFunction(), _chargeSplitterErfFunction(),_theChosenSplitter(NULL), _effMap(NULL) {
   _description = "the transfer Between Energy and Induced Charge for SDHCAL" ;
 
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

   
   _outputHcalCollections.push_back(std::string("HCALBarrel"));
   _outputHcalCollections.push_back(std::string("HCALEndcap"));
   _outputHcalCollections.push_back(std::string("HCALOther"));

//    registerProcessorParameter("SDHCALThreshold" , 
// 			     "Threshold for SDHCAL Hits in GeV" ,
// 			     _thresholdSDHcal,
// 			     (float)2.0e-10);
   std::vector<float> hcalThresholds;
   hcalThresholds.push_back(0.01);
   registerProcessorParameter("HCALThreshold" ,
                              "Threshold for HCAL Hits in pC" ,
			      _thresholdHcal,
			      hcalThresholds);
   
   std::vector<float> calibrHcal;
   calibrHcal.push_back(1.0);
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
   

  registerProcessorParameter("EffMapOption",
			     "Step efficiency correction method : should be Constant or PrototypeMap",
			     _effMapOption,
			     std::string("Constant")) ;

  registerProcessorParameter("EffMapConstantValue",
			    "Value of the constant term for efficiency correction if EffMapOption==Constant",
			    _constEffMapValue,
			    (float)0.97);

  registerProcessorParameter("EffMapPrototypeFileName",
			    "File name where prototype efficiency corresction is stored if EffMapOption==PrototypeMap",
			    _effMapFileName,
			    std::string("map.txt"));

  //  multiplicityChargeSplitterUniform parameters
  registerProcessorParameter( "CellEdgeDistance" ,
			      "Distance (in mm) from the border of the cell that implies hit on neighbouring cell for 'Uniform' charge splitter : used if ChargeSplitterOption==Uniform"  ,
			      _chargeSplitterUniform._edgeSize,
			      (float)2.0 ) ;

   
  //  multiplicityChargeSplitterFuntion parameters
  registerProcessorParameter( "functionFormula" ,
			      "ROOT TF2 formula for function describing the induced charge distribution, x and y are coordinates from the step position in I and J direction : used if ChargeSplitterOption==Function",
			      _chargeSplitterFunction._function_description,
			      std::string("1/cosh([0]*sqrt(x*x+y*y))") ); 
   
  registerProcessorParameter( "functionRange",
			      "maximal distance (in mm) at which a step can induce charge using the 2D function defined with functionFormula or when using ChargeSplitterOption==Erf",
			      _chargeSplitterFunction._range,
			      (float) 30 );

  std::vector<float> chadisParameter;
  chadisParameter.push_back(3.1415/0.24);
  registerProcessorParameter( "functionParameters",
			      "parameter values for the function defined with functionFormula",
			      _chargeSplitterFunction._functionParameters,
			      chadisParameter );


  registerProcessorParameter( "RPC_PadSeparation",
			      "distance in mm between two RPC pads : used if ChargeSplitterOption==Function or Erf", 
			      _chargeSplitterFunction._RPC_PadSeparation,
                             (float) 0.0 );

  std::vector<float> erfWidth;
  erfWidth.push_back(2);
  registerProcessorParameter( "erfWidth",
			      "Width values for the different Erf functions",
			      _chargeSplitterErfFunction._erfWidth,
			      erfWidth );

  std::vector<float> erfWeigth;
  erfWeigth.push_back(1);
  registerProcessorParameter( "erfWeigth",
			      "Weigth for the different Erf functions",
			      _chargeSplitterErfFunction._erfWeigth,
			      erfWeigth );
   
  _chargeSplitterOption="Erf";
  registerProcessorParameter( "ChargeSplitterOption",
			      "Define the charge splitter method. Possible option : Uniform, Function, Erf",
			      _chargeSplitterOption,
			      std::string("Erf") );

  registerProcessorParameter( "ChargeSplitterRandomSeed",
			      "The seed of the srand in the multiplicityChargeSplitterFunction",
			      _chargeSplitterRandomSeed,
			      int(1) );

  registerProcessorParameter( "PolyaRandomSeed",
			      "The seed of the polya function",
			      _polyaRandomSeed,
			      int(1) );

  registerProcessorParameter( "CellIDEncodingStringType",
			      "The type of the encoding, LCGEO or MOKKA",
			      _encodingType,
			      std::string("LCGEO"));

  registerProcessorParameter( "HCALOption",
			      "The HCAL mechanical options, TESLA or VIDEAU, which is effective when gear is absent",
			      _hcalOption,
			      std::string("VIDEAU"));

  _doThresholds=true;
  registerOptionalParameter("doThresholds",
			    "Replace analog hit energy by value given in CalibrHCAL according to thresholds given in HCALThreshold",
			    _doThresholds,
			    true);

  registerProcessorParameter( "PolyaAverageCharge" ,
			      "Parameter for the Polya distribution used to simulate the induced charge distribution : mean of the distribution", 
			      _polyaAverageCharge, 
			      (double) 1.6 );

  registerProcessorParameter( "PolyaWidthParameter" ,
			      "Parameter for the Polya distribution used to simulate the induced charge distribution : related to the distribution width ",
			      _polyaFunctionWidthParameter, 
			      (double) 16.3 );


  registerProcessorParameter( "StepCellCenterMaxDistanceLayerDirection",
			      "Maximum distance (mm) between the Geant4 step position and the cell center, in the RPC width direction, to keep a step for digitization",
			      _absZstepFilter,
			      (float) 0.0005 );

  registerProcessorParameter( "StepsMinDistanceRPCplaneDirection",
			      "Minimum distance (mm) between 2 Geant4 steps, in the RPC plane, to keep the 2 steps",
			      _minXYdistanceBetweenStep,
			      (float) 0.5 );

  registerProcessorParameter( "KeepAtLeastOneStep",
			      "if true, ensure that each hit will keep at least one step for digitisation independatly of filtering conditions (StepCellCenterMaxDistanceLayerDirection)",
			      _keepAtLeastOneStep,
			      true );

  _relcol=0;
  _debugTupleStepFilter=0;
  _tupleStepFilter=0;

}
//define a Polya Function 

Double_t Polya(double* x, double* params)
{
  // Polya function.
  //0 = Polya - average charge
  // 1 = Polya - free parameter related to the width
  Double_t q, t;
  if(*x <= 0)
    return 0;
  q = *x / params[0];
  t = 1.0 + params[1];
  return TMath::Power(t*q, params[1])* TMath::Exp(-t*q);
} 
  
void SimDigital::init(){
  std::cout << "SimDigital: init" << std::endl;

  SimDigitalGeomCellId::setEncodingType(_encodingType);
  SimDigitalGeomCellId::setHcalOption(_hcalOption);

  //a=1.05,b=-2.1,c=1.05,d=-1.05;//for 1 neighbor pad fired
  //previousLayer=1 ,previousmodule=1, previousstave=1,previousIy=1,previousJz=1;
  //define a polya function to get induced charge
  _QPolya = new TF1("QPolya",Polya,0.0,30.0,2);
  _QPolya->SetNpx(200);
  _QPolya->SetParameters(_polyaAverageCharge,_polyaFunctionWidthParameter);
  //define the function for charge distribution
  //c2d=new TF2("c2d",ChaDis,-10,10,-10,10,3);
  //c2d->SetParNames("width","x'value","y'value");
  
  _chargeSplitterFunction.init();
  _chargeSplitterErfFunction._range=_chargeSplitterFunction._range;
  _chargeSplitterErfFunction._RPC_PadSeparation=_chargeSplitterFunction._RPC_PadSeparation;
  _chargeSplitterErfFunction.init();
  
  //fg: need to set default encoding in for reading old files...
// CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;

  if(_chargeSplitterOption=="Erf") _theChosenSplitter=&_chargeSplitterErfFunction;
  else if(_chargeSplitterOption=="Function") _theChosenSplitter=&_chargeSplitterFunction;
  else if(_chargeSplitterOption=="Uniform") _theChosenSplitter=&_chargeSplitterUniform;
  else throw ParseException( _chargeSplitterOption + std::string(" option for charge splitting is not available ") );

  if (_effMap!=NULL) delete _effMap;
  if (_effMapOption=="Constant") _effMap=new effMapConstant(_constEffMapValue);
  else if (_effMapOption=="PrototypeMap") _effMap=new effMapProtoByAsic(_effMapFileName);
  else throw ParseException( _effMapOption + std::string(" option for efficiency correction is not available") );

  streamlog_out( MESSAGE ) << "_effMapOption = " << _effMapOption << std::endl;
  if (_effMapOption=="Constant")   
    streamlog_out( MESSAGE ) << "_constEffMapValue = " << _constEffMapValue << std::endl;
  else if (_effMapOption=="PrototypeMap") 
    streamlog_out( MESSAGE ) << "_effMapFileName = " << _effMapFileName << std::endl;

  //assure SDHCAL _thresholdHcal are in increasing order
  std::sort(_thresholdHcal.begin(),_thresholdHcal.end());

  //book tuples
  _debugTupleStepFilter  = AIDAProcessor::tupleFactory( this )->create("SimDigitalStepDebug",
								       "SimDigital_StepDebug",
								       "int filterlevel, float deltaI,deltaJ,deltaLayer,minIJdist"); 
  streamlog_out(DEBUG) << "Tuple for step debug has been initialized to " << _debugTupleStepFilter << std::endl;
  streamlog_out(DEBUG) << "it has " << _debugTupleStepFilter->columns() << " columns" <<std::endl;


  _tupleStepFilter   = AIDAProcessor::tupleFactory( this )->create("SimDigitalStepStat",
								   "SimDigital_StepStat",
								   "int allsteps, absZfiltered, IJdistancefiltered");
  streamlog_out(DEBUG) << "Tuple for step stat has been initialized to " << _tupleStepFilter << std::endl;
  streamlog_out(DEBUG) << "it has " << _tupleStepFilter->columns() << " columns" <<std::endl;


  _tupleCollection  = AIDAProcessor::tupleFactory( this )->create("CollectionStat",
								  "Collection_statistics",
								  "int NsimHit, NrecoHit, N1, N2, N3"); 
  streamlog_out(DEBUG) << "Tuple for collection stat has been initialized to " << _tupleCollection << std::endl;
  streamlog_out(DEBUG) << "it has " << _tupleCollection->columns() << " columns" <<std::endl;
}
 

void SimDigital::processRunHeader( LCRunHeader* run) {}


void SimDigital::processHCAL(LCEvent* evt, LCFlagImpl& flag)
{
  depositedEnergyInRPC=0.;
//  cout<< "_hcalCollections.size() = "<< _hcalCollections.size() << endl;
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    try{
      std::string colName =  _hcalCollections[i] ;    
      cout<< "colName[i] = "<< colName <<" "<< i << endl;
      //CHT::Layout layout1 = layoutFromString( colName ); 
      LCCollection * col = evt->getCollection( colName.c_str() ) ;
      //CHT::Layout layout2 = layoutFromString( colName ); 
      _counters["NSim"]+=col->getNumberOfElements();
      CHT::Layout layout = layoutFromString( colName ); 
      LCCollectionVec *hcalcol = processHCALCollection(col,layout,flag);
	  std::cout << " ------ " << hcalcol << std::endl;
//      cout<< " CHT::any,barrel,encap, ring " << CHT::any<<" "<<CHT::barrel<<" "<<CHT::endcap<<" "<< CHT::ring<< endl;
      _counters["NReco"]+=hcalcol->getNumberOfElements();
      evt->addCollection(hcalcol,_outputHcalCollections[i].c_str());
    }
    catch(DataNotAvailableException &e){  
    }     
  }
  evt->parameters().setValue("totalVisibleEnergy",depositedEnergyInRPC);
}

//multiplicityChargeSplitterBase& SimDigital::getSplitter()
//{
//  if(_chargeSplitterOption=="Erf") return _chargeSplitterErfFunction;
//  else if(_chargeSplitterOption=="Uniform") return _chargeSplitterUniform;
//  else if(_chargeSplitterOption=="Function") return _chargeSplitterFunction;
//  else throw ParseException( _chargeSplitterOption + std::string(" option for charge splitting is not available ") );
//}

multiplicityChargeSplitterBase::multiplicityChargeSplitterBase() : _chargeMap()
{}


multiplicityChargeSplitterUniform::multiplicityChargeSplitterUniform()
  : multiplicityChargeSplitterBase()
{
  _edgeSize=1;
}


void multiplicityChargeSplitterUniform::addCharge(float charge, float pos_I, float pos_J)
{
  //It is assumed that within a cell coordinate 0,0 is the cell center
  std::set<I_J_Coordinates> ModifyCharge;
  ModifyCharge.insert(I_J_Coordinates(0,0));
  int I=0; int J=0;
  if (_edgeSize<=_cellSize)
    {
      if (pos_I<-_cellSize/2+_edgeSize) I=-1;
      if (pos_I>_cellSize/2-_edgeSize) I=+1;
      if (pos_J<-_cellSize/2+_edgeSize) J=-1;
      if (pos_J>_cellSize/2-_edgeSize) J=+1;
      ModifyCharge.insert(I_J_Coordinates(0,J));
      ModifyCharge.insert(I_J_Coordinates(I,0));
      ModifyCharge.insert(I_J_Coordinates(I,J));
    }
  else
    {
      int Ncells=int(_edgeSize/_cellSize);
      for (int iLR=-Ncells; iLR<=Ncells; iLR++)
	for (int iLU=-Ncells; iLU<=Ncells; iLU++)
	  ModifyCharge.insert(I_J_Coordinates(iLR,iLU));
    }
  float SplitCharge=charge/ModifyCharge.size();
  for (std::set<I_J_Coordinates>::iterator it=ModifyCharge.begin(); 
       it != ModifyCharge.end(); it++)
    _chargeMap[*it]+=SplitCharge;
}


multiplicityChargeSplitterFunction::multiplicityChargeSplitterFunction()
  : multiplicityChargeSplitterBase(),
    _f2(NULL),_range(1),_function_description(""),_functionParameters(),
    _normalisation(1),_RPC_PadSeparation(0)
{}

multiplicityChargeSplitterFunction::~multiplicityChargeSplitterFunction()
{
  if (_f2) delete _f2;
}

void multiplicityChargeSplitterFunction::init()
{
  if (_f2) delete _f2;
  _f2=new TF2("chadis",_function_description.c_str(),-_range,+_range,-_range,+_range);
  for (unsigned int i=0; i<_functionParameters.size(); i++)
    _f2->SetParameter(i,_functionParameters[i]);
  _normalisation=_f2->Integral(-_range,+_range,-_range,+_range);

  double a[2],b[2];
  a[0]=-_range; b[0]=_range;
  a[1]=-_range; b[1]=_range;
  double relerr(0); int nfnevl(0), ifail(1);
  double precision=1e-5;
  int maxpts=1000;
  float _normalisation2=0;
  while (ifail !=0 && maxpts<=20000)
    {
      _normalisation2=_f2->IntegralMultiple(2,a,b,17,maxpts,precision,relerr,nfnevl,ifail);
      if (ifail!=0)  streamlog_out(DEBUG) << "ifail= " << ifail << " maxpts= " << maxpts << "  relerr=" << relerr << std::endl;
      maxpts*=10;
    }
  streamlog_out( DEBUG ) << "Charge splitter normalisation factor: " << _normalisation << std::endl;
  streamlog_out( DEBUG ) << "range : " << _range << " ; padseparation : " << _RPC_PadSeparation << std::endl;


}

void multiplicityChargeSplitterFunction:: addCharge(float charge, float pos_I, float pos_J)
{
  streamlog_out( DEBUG ) << "Charge splitter cellSize: " << _cellSize << std::endl;
  if (_RPC_PadSeparation>_cellSize) return;
  int icell=int(_range/_cellSize);
  for (int I=-icell; I<=icell; I++)
    {
      float distI=(I-0.5)*_cellSize-pos_I+_RPC_PadSeparation/2;
      float distII=(I+0.5)*_cellSize-pos_I-_RPC_PadSeparation/2;
      if (distI<-_range) distI=-_range;
      if (distII>_range) distII=_range;
      for (int J=-icell; J<=icell; J++)
	{
	  float distJ=(J-0.5)*_cellSize-pos_J+_RPC_PadSeparation/2;
	  float distJJ=(J+0.5)*_cellSize-pos_J-_RPC_PadSeparation/2;
	  if (distJ<-_range) distJ=-_range;
	  if (distJJ>_range) distJJ=_range;
	  double a[2],b[2];
	  a[0]=distI; b[0]=distII;
	  a[1]=distJ; b[1]=distJJ;
	  double relerr(0); int nfnevl(0), ifail(1);
	  double precision=1e-5;
	  double integralResult(0);
	  int maxpts=1000;
	  while (ifail !=0 && maxpts<=20000)
	    {
	      integralResult=_f2->IntegralMultiple(2,a,b,17,maxpts,precision,relerr,nfnevl,ifail);
	      if (ifail!=0)  streamlog_out(DEBUG) << "ifail= " << ifail << " maxpts= " << maxpts << "  relerr=" << relerr << std::endl;
	      maxpts*=10;
	    }
	  if(integralResult>=0) _chargeMap[I_J_Coordinates(I,J)]+=charge*integralResult/_normalisation;
	  else _chargeMap[I_J_Coordinates(I,J)]+=charge*_f2->Integral(distI,distII,distJ,distJJ)/_normalisation;
	}
    }
}


multiplicityChargeSplitterErfFunction::multiplicityChargeSplitterErfFunction()
  : multiplicityChargeSplitterBase(),
    _range(1),_erfWidth(),_erfWeigth(),
    _normalisation(0),_RPC_PadSeparation(0)
{}

multiplicityChargeSplitterErfFunction::~multiplicityChargeSplitterErfFunction()
{

}

void multiplicityChargeSplitterErfFunction::init()
{
  if(_erfWidth.size()!=_erfWeigth.size()) throw ParseException( " Different size for erfWidth erfWeigth " );
  for(unsigned int i=0; i<_erfWidth.size(); i++){
    streamlog_out( DEBUG ) << "Erf function parameters " << i+1 << " : " << _erfWidth[i] << ", " << _erfWeigth[i] << std::endl;
    _normalisation+=_erfWeigth[i]*_erfWidth[i]*_erfWidth[i]*TMath::Pi()*TMath::Erf(_range/_erfWidth[i])*TMath::Erf(_range/_erfWidth[i]);
  }
  streamlog_out( DEBUG ) << "Charge splitter normalisation factor: " << _normalisation << std::endl;
  streamlog_out( DEBUG ) << "range : " << _range << " ; padseparation : " << _RPC_PadSeparation << std::endl;

}

void multiplicityChargeSplitterErfFunction:: addCharge(float charge, float pos_I, float pos_J)
{
  if (_RPC_PadSeparation>_cellSize) return;
  int icell=int(_range/_cellSize);
  double chargeTotCheck=0;
  for (int I=-icell; I<=icell; I++)
    {
      float distI=(I-0.5)*_cellSize-pos_I+_RPC_PadSeparation/2;
      float distII=(I+0.5)*_cellSize-pos_I-_RPC_PadSeparation/2;
      if (distI<-_range) distI=-_range;
      if (distII>_range) distII=_range;
      for (int J=-icell; J<=icell; J++)
	{
	  float distJ=(J-0.5)*_cellSize-pos_J+_RPC_PadSeparation/2;
	  float distJJ=(J+0.5)*_cellSize-pos_J-_RPC_PadSeparation/2;
	  if (distJ<-_range) distJ=-_range;
	  if (distJJ>_range) distJJ=_range;
	  double a[2],b[2];
	  a[0]=distI; b[0]=distII;
	  a[1]=distJ; b[1]=distJJ;
	  
	  double integralResult=0;
	  for(unsigned int n=0; n<_erfWidth.size(); n++){
	    integralResult+=fabs(TMath::Erf(distII/_erfWidth[n])-TMath::Erf(distI/_erfWidth[n]))*
	      fabs(TMath::Erf(distJJ/_erfWidth[n])-TMath::Erf(distJ/_erfWidth[n]))*
	      _erfWeigth[n]*TMath::Pi()*_erfWidth[n]*_erfWidth[n]/4;
	  }
	  _chargeMap[I_J_Coordinates(I,J)]+=charge*integralResult/_normalisation ;
	  if(_chargeMap[I_J_Coordinates(I,J)]<0) 
	    streamlog_out( MESSAGE ) << "!!!!!!!!!!Negative Charge!!!!!!!!!!" << std::endl
				     << " X " << pos_J << " " << distJ << " " << distJJ << std::endl
				     << " Y " << pos_I << " " << distI << " " << distII << std::endl;
	  chargeTotCheck+=_chargeMap[I_J_Coordinates(I,J)];
	}
    }
  streamlog_out( DEBUG ) << " Charge = " << charge
			 << " ; total splitted charge = " << chargeTotCheck 
			 << std::endl;
}


effMapProtoByAsic::effMapProtoByAsic(std::string fileName)
{
  std::ifstream in;
  in.open(fileName.c_str());
  if(!in.is_open())
    {streamlog_out(ERROR) << fileName << "\t NO SUCH FILE IN CURRENT DIRECTORY" << std::endl;}
  else {streamlog_out(DEBUG) << "MAP FILE IN \t " << fileName << std::endl;}
  int asickey,nevent;
  double efficiency,multiplicity,efficiencyError;
  while(1){
    if(!in.good()) break;
    in >> asickey >> nevent >> efficiency >> efficiencyError >> multiplicity;
    _effMap[asickey]=efficiency;
  } 
}



AIDA::ITuple* SimDigitalGeomCellId::_tupleHit = NULL;
AIDA::ITuple* SimDigitalGeomCellId::_tupleStep = NULL;


void SimDigitalGeomCellId::bookTuples(const marlin::Processor* proc)
{
  _tupleHit  = AIDAProcessor::tupleFactory( proc )->create("SimDigitalGeom",
							   "SimDigital_Debug",
							   "int detector,chtlayout,module,tower,stave,layer,I,J, float x,y,z, normalx,normaly,normalz, Ix,Iy,Iz,Jx,Jy,Jz"); 
  streamlog_out(DEBUG) << "Tuple for Hit has been initialized to " << _tupleHit << std::endl;
  streamlog_out(DEBUG)<< "it has " << _tupleHit->columns() << " columns" <<std::endl;
  
  _tupleStep = AIDAProcessor::tupleFactory( proc )->create("SimDigitalStep",
							   "SimDigital_DebugStep",
							   "int detector,chtlayout,hitcellid,nstep, float hitx,hity,hitz,stepx,stepy,stepz,deltaI,deltaJ,deltaLayer");
  streamlog_out(DEBUG) << "Tuple for Step has been initialized to " << _tupleStep << std::endl;
  streamlog_out(DEBUG) << "it has " << _tupleStep->columns() << " columns" <<std::endl;
}


//    _normal_I_J_setter(new SimDigitalGeomRPCFrame_TESLA_BARREL(*this)),
//    _layerLayout( & Global::GEAR->getHcalEndcapParameters().getLayerLayout() ),
//    _normal_I_J_setter(new SimDigitalGeomRPCFrame_VIDEAU_ENDCAP(*this)),
SimDigitalGeomCellId::SimDigitalGeomCellId(LCCollection *inputCol, LCCollectionVec *outputCol) 
  : _hitPosition(NULL),
    _decoder(inputCol), _encoder(inputCol->getParameters().getStringVal(LCIO::CellIDEncoding),outputCol),
    _layerLayout(NULL), caloData(NULL),
    _normal_I_J_setter(NULL),
    _currentHCALCollectionCaloLayout(CHT::any),
    _normal(),_Iaxis(),_Jaxis(),
	_useGear(false),
	_isInEndcap(false)
{
    _trueLayer=_stave=_module=_Iy=_Jz=-999;
    outputCol->parameters().setValue(LCIO::CellIDEncoding,inputCol->getParameters().getStringVal(LCIO::CellIDEncoding));

	_cellIDEncodingString = inputCol->getParameters().getStringVal(LCIO::CellIDEncoding);

	if(_cellIDEncodingString.find("z:") != std::string::npos ) _isInEndcap = true;

	std::cout << "_isInEndcap: " << _isInEndcap << std::endl;

	std::string gearFile = Global::parameters->getStringVal("GearXMLFile") ;
	if(gearFile.size() > 0) _useGear = true;

    if(_useGear) 
    {
        _geom=VIDEAU;

        try
        { 
			std::cout << "we will use gear!" << std::endl;
            cout << "gear: " << Global::GEAR << endl;
            Global::GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_order"); //it is VIDEAU geometry if it is OK
        }
        catch (gear::Exception &)
        {
            _geom=TESLA;
        }

    }
    else
    {
		std::cout << "we will use lcgeo!" << std::endl;
        if(_hcalOption == std::string("TESLA"))  _geom = TESLA;
        if(_hcalOption == std::string("VIDEAU")) _geom = VIDEAU;
    }

    streamlog_out( MESSAGE )<< "!!!!!!(Videau=0, TESLA=1) Geometry is _geom= : "<<_geom << std::endl;
}


SimDigitalGeomCellId::~SimDigitalGeomCellId()
{
	if(_normal_I_J_setter != NULL) delete _normal_I_J_setter;
}

void SimDigitalGeomRPCFrame_TESLA_BARREL::setRPCFrame()
{
  normal().set(1,0,0);
  Iaxis().set(0,1,0);
  Jaxis().set(0,0,1);
  double angle=(stave()/2)*(-45)*CLHEP::degree;
  normal().rotateZ(angle);
  Iaxis().rotateZ(angle);
}

void SimDigitalGeomRPCFrame_VIDEAU_BARREL::setRPCFrame()
{
  normal().set(0,1,0);
  Iaxis().set(1,0,0);
  Jaxis().set(0,0,1);
  double angle=(stave())*(45)*CLHEP::degree;
  normal().rotateZ(angle);
  Iaxis().rotateZ(angle);
}

void SimDigitalGeomRPCFrame_TESLA_ENDCAP::setRPCFrame()
{
  if (module()==6) 
    {
      normal().set(0,0,1);
      Iaxis().set(1,0,0);
      Jaxis().set(0,1,0);
      Iaxis().rotateZ(stave()*(-90)*CLHEP::degree);
      Jaxis().rotateZ(stave()*(-90)*CLHEP::degree);
    }
  else if (module()==4) 
    {
      normal().set(0,0,-1);
      Iaxis().set(-1,0,0);
      Jaxis().set(0,1,0);
      Iaxis().rotateZ(stave()*90*CLHEP::degree);
      Jaxis().rotateZ(stave()*90*CLHEP::degree);
    }
  else
    {
      streamlog_out(ERROR) << "ERROR ; TESLA detector : unknown module for endcap " << module() << std::endl;
    }      
}

//valid also for all endcap rings (VIDEAU and TESLA)
void SimDigitalGeomRPCFrame_VIDEAU_ENDCAP::setRPCFrame()
{
  if (module()==6)
    {
//      cout<< "!!!!!!!!!! VIDEAU_ENDCAP : module=6 "<< std::endl;
      normal().set(0,0,1);
      Iaxis().set(1,0,0);
      Jaxis().set(0,1,0);
      Iaxis().rotateZ(stave()*(90)*CLHEP::degree);
      Jaxis().rotateZ(stave()*(90)*CLHEP::degree);
    }
  else if (module()==0)
    {
//      cout<< "!!!!!!!!!! VIDEAU_ENDCAP : module 0"<< std::endl;
      normal().set(0,0,-1);
      Iaxis().set(-1,0,0);
      Jaxis().set(0,1,0);
      Iaxis().rotateZ(stave()*(-90)*CLHEP::degree);
      Jaxis().rotateZ(stave()*(-90)*CLHEP::degree);
    }
  else
    {
      streamlog_out(ERROR) << "ERROR : unknown module for endcap or endcapring " << module() << std::endl;
    }
}



std::vector<StepAndCharge> SimDigitalGeomCellId::decode(SimCalorimeterHit *hit)
{
  if(_isInEndcap) 
  {
     _encodingStrings[_encodingType][5] = "z";
  }
  else
  {
     _encodingStrings[_encodingType][5] = "y";
  }

  _trueLayer = _decoder( hit )[_encodingStrings[_encodingType][0]] - 1; // -1;
  _stave     = _decoder( hit )[_encodingStrings[_encodingType][1]];     // +1
  _module    = _decoder( hit )[_encodingStrings[_encodingType][2]];

  if(_encodingStrings[_encodingType][3].size() != 0)
  _tower     = _decoder( hit )[_encodingStrings[_encodingType][3]];

  _Iy        = _decoder( hit )[_encodingStrings[_encodingType][4]];
  _Jz        = _decoder( hit )[_encodingStrings[_encodingType][5]];



//
// _slice     = _decoder( hit )["slice"];
  _hitPosition = hit->getPosition();
  if(abs(_Iy)<1 && abs(_Iy!=0.0)) streamlog_out(MESSAGE) << "_Iy, _Jz:"<<_Iy <<" "<<_Jz<< std::endl;
//if(_module==0||_module==6) streamlog_out(MESSAGE)<<"tower "<<_tower<<" layer "<<_trueLayer<<" stave "<<_stave<<" module "<<_module<<std::endl;
//<<" Iy " << _Iy <<"  Jz "<<_Jz<<" hitPosition "<<_hitPosition<<std::endl
//<<" _hitPosition[0] "<<_hitPosition[0]<<" _hitPosition[1] "<<_hitPosition[1]<<" _hitPosition[2] "<<_hitPosition[2]<<std::endl;

  if(_useGear)
  {
	  _normal_I_J_setter->setRPCFrame();
  }
  else
  { 
	  // this part is for lcgeo

	  dd4hep::Detector & ild = dd4hep::Detector::getInstance();
 	  dd4hep::rec::CellIDPositionConverter idposConv( ild )  ;


 	  dd4hep::BitField64 idDecoder( _cellIDEncodingString ) ;
 	  
 	  dd4hep::long64 id0 = hit->getCellID0() ;
 	  dd4hep::long64 id1 = hit->getCellID1() ;
 	  
 	  idDecoder.setValue( id0 , id1 ) ;
 	  
 	  dd4hep::long64 id = idDecoder.getValue() ;
 	  dd4hep::Position pos_0 = idposConv.position( id ) ;   

 	  const float* hitPos = hit->getPosition();

#if 1 
 	  std::cout << "hit pos: " << hitPos[0] << " " << hitPos[1] << " " << hitPos[2] << std::endl;
 	  std::cout << "cell pos: " << pos_0.X() << " " << pos_0.Y() << " " << pos_0.Z() << std::endl;

 	  std::cout << "layer: "    << idDecoder[_encodingStrings[_encodingType][0]] 
 	            << ", stave: "  << idDecoder[_encodingStrings[_encodingType][1]] 
 	    		<< ", module: " << idDecoder[_encodingStrings[_encodingType][2]]
 	            << ", tower: "  << idDecoder[_encodingStrings[_encodingType][3]] 
 	    		<< ", x: "      << idDecoder[_encodingStrings[_encodingType][4]]
 	    		<< ", y: "      << idDecoder[_encodingStrings[_encodingType][5]] << std::endl;
#endif

 	  double const epsilon = 1.e-3;

 	  ///// for direction x
 	  dd4hep::Position pos_i_plus_1;
 	  dd4hep::Position dir_i;

 	  std::string xEncoding = _encodingStrings[_encodingType][4];
 	  std::string yEncoding = _encodingStrings[_encodingType][5];

 	  idDecoder[ xEncoding ]   =  idDecoder[ xEncoding ]  + 1 ;
 	  pos_i_plus_1 = idposConv.position(   idDecoder.getValue()   ) ;

 	  dir_i =  pos_i_plus_1  - pos_0  ; 

 	  if(dir_i.R() < epsilon)
 	  { 
 	      idDecoder[ xEncoding ]   =  idDecoder[ xEncoding ] - 2 ;
 	      pos_i_plus_1 = idposConv.position(   idDecoder.getValue()   ) ;
 	      dir_i =  - ( pos_i_plus_1  - pos_0 )  ; 
 	  }

 	  // reset
 	  idDecoder.setValue( id0 , id1 ) ;	

 	  ////// for direction y
 	  dd4hep::Position pos_j_plus_1;
 	  dd4hep::Position dir_j;
 	  
 	  idDecoder[ yEncoding ]   =  idDecoder[ yEncoding ]  + 1 ;
 	  pos_j_plus_1 = idposConv.position(   idDecoder.getValue()   ) ;
 	  
 	  dir_j =  pos_j_plus_1  - pos_0  ; 

 	  if(dir_j.R() < epsilon)
 	  { 
 	      idDecoder[ yEncoding ]   =  idDecoder[ yEncoding ] - 2 ;
 	      pos_j_plus_1 = idposConv.position(   idDecoder.getValue()   ) ;
 	      dir_j =  - ( pos_j_plus_1  - pos_0 )  ; 
 	  }
 	  
 	  dd4hep::Position dir_layer = dir_i.Cross( dir_j ) ;

 	  dir_layer = - dir_layer.Unit();

 	  std::cout << "layer dir: " << dir_layer.X() << " " << dir_layer.Y() << " " << dir_layer.Z() << std::endl;

	  dir_i = dir_i.Unit();
	  dir_j = dir_j.Unit();

      _normal.set(dir_layer.X(), dir_layer.Y(), dir_layer.Z());
      _Iaxis.set(dir_i.X(), dir_i.Y(), dir_i.Z());
      _Jaxis.set(dir_j.X(), dir_j.Y(), dir_j.Z());
  }

  std::vector<StepAndCharge> stepsInIJZcoord;
  LCVector3D hitpos;
  if (NULL != _hitPosition) hitpos.set(_hitPosition[0],_hitPosition[1],_hitPosition[2]);
  for (int imcp=0; imcp<hit->getNMCContributions(); imcp++)
    {
      LCVector3D stepposvec;
      const float* steppos=hit->getStepPosition(imcp);
      if (NULL != steppos) stepposvec.set(steppos[0],steppos[1],steppos[2]);
      if (stepposvec.mag2() != 0)
	{
	  stepposvec-=hitpos;
	  stepsInIJZcoord.push_back( StepAndCharge(LCVector3D(stepposvec*_Iaxis,stepposvec*_Jaxis,stepposvec*_normal)) );
	}
      else
	streamlog_out(WARNING) << "DIGITISATION : STEP POSITION IS (0,0,0)" << std::endl;
    }

  //if no steps have been found, then put one step at the center of the cell :
  if (stepsInIJZcoord.size()==0) stepsInIJZcoord.push_back(StepAndCharge(LCVector3D(0,0,0)));
  

  //these tuples are for debugging geometry aspects
  if (_tupleHit!=NULL) 
    {
      _tupleHit->fill(TH_DETECTOR,int(_geom));
      _tupleHit->fill(TH_CHTLAYOUT,int(_currentHCALCollectionCaloLayout));
      _tupleHit->fill(TH_MODULE,_module);
      _tupleHit->fill(TH_TOWER,_tower);
      _tupleHit->fill(TH_STAVE,_stave);
      _tupleHit->fill(TH_LAYER,_trueLayer);
      _tupleHit->fill(TH_I,_Iy);
      _tupleHit->fill(TH_J,_Jz);
      if (_hitPosition != NULL)
	{
	  _tupleHit->fill(TH_X,_hitPosition[0]); //x
	  _tupleHit->fill(TH_Y,_hitPosition[1]); //y
	  _tupleHit->fill(TH_Z,_hitPosition[2]); //z
	}
      else
	{
	  float notset=-88888;
	  _tupleHit->fill(TH_X,notset);
	  _tupleHit->fill(TH_Y,notset);
          _tupleHit->fill(TH_Z,notset);
	}
      for (int i=0; i<3; i++)
	{
	  _tupleHit->fill(TH_NORMALX+i,_normal[i]);
	  _tupleHit->fill(TH_IX+i,_Iaxis[i]);
	  _tupleHit->fill(TH_JX+i,_Jaxis[i]);
	}
      _tupleHit->addRow();
    }
  if (_tupleStep!=NULL)
    {
      int nsteps=hit->getNMCContributions();
      float notset=-88888;
      for (int imcp=0; imcp<nsteps; imcp++)
	{
	  _tupleStep->fill(TS_DETECTOR,int(_geom));
	  _tupleStep->fill(TS_CHTLAYOUT,int(_currentHCALCollectionCaloLayout));
	  _tupleStep->fill(TS_HITCELLID,hit->getCellID0());
	  _tupleStep->fill(TS_NSTEP,nsteps);
	  const float *steppos=hit->getStepPosition(imcp);
	  for (int i=0; i<3; i++)
	    {
	      if (_hitPosition != NULL ) _tupleStep->fill(TS_HITX+i,_hitPosition[i]);
	      else _tupleStep->fill(TS_HITX+i,notset);
	      if (steppos != NULL) {_tupleStep->fill(TS_STEPX+i,steppos[i]);}
	      else {_tupleStep->fill(TS_STEPX+i,notset);}
	      {
		if (imcp<(int)stepsInIJZcoord.size())
		  _tupleStep->fill(TS_DELTAI+i,stepsInIJZcoord[imcp].step[i]);
		else
		  _tupleStep->fill(TS_DELTAI+i,notset);
	      }
	    }
	  _tupleStep->addRow();
	}
    }

  return stepsInIJZcoord;
}

void SimDigitalGeomCellId::encode(CalorimeterHitImpl *hit, int delta_I, int delta_J)
{
  _encoder[_encodingStrings[_encodingType][2]]=_module;

  if(_encodingStrings[_encodingType][3].size() != 0)
  _encoder[_encodingStrings[_encodingType][3]]=_tower;

  _encoder[_encodingStrings[_encodingType][1]]=_stave;

  int RealIy=_Iy+delta_I;
//  streamlog_out( MESSAGE ) << "RealIy, _Iy, delta_I" << std::endl
//       <<RealIy <<" "<<_Iy <<" " <<delta_I<<std::endl;
//  if (RealIy<0||RealIy>330)RealIy=0; //FIXME the 330 value should depend on the cellSize and on the Layer
//  if (RealIy<0)RealIy=RealIy+176; //FIXME the 330 value should depend on the cellSize and on the Layer
  if (abs(RealIy)>330)RealIy=0; //FIXME the 330 value should depend on the cellSize and on the Layer
  _encoder[_encodingStrings[_encodingType][4]]=RealIy;
  int RealJz=_Jz+delta_J;
//  streamlog_out( MESSAGE ) << "RealIy, _Iy, delta_I" <<RealIy<<" "<<_Iy<<" "<<delta_I<< std::endl;
//  streamlog_out( MESSAGE ) << "RealJz, _Jz, delta_J" <<RealJz<<" "<<_Jz<<" "<<delta_J<< std::endl;
//       <<RealJz <<" "<<_Jz <<" " <<delta_J<<std::endl;
//  if(RealJz<0||RealJz>235)RealJz=0; //FIXME the 330 value should depend on the cellSize and on the Layer
//  if (RealJz<0)RealJz=RealJz+47; //FIXME the 330 value should depend on the cellSize and on the Layer
  if (abs(RealJz)>330)RealJz=0; //FIXME the 330 value should depend on the cellSize and on the Layer
//  if (abs(RealIy)>330||abs(RealJz)>330)streamlog_out( MESSAGE ) << "RealIy, RealJz" << std::endl;
  if (abs(RealIy)>330||abs(RealJz)>330)streamlog_out( MESSAGE ) << "RealIy, RealJz" << std::endl
                                      <<RealIy <<"   "<<RealJz <<std::endl;
  _encoder[_encodingStrings[_encodingType][0]]=_trueLayer   ;
// Depending on the segmentation type:   Barrel,EndcapRing - CartesianGridXY;  EndCaps- CartesianGridXZ !!!

  if(!_isInEndcap)
  {
	  _encoder[_encodingStrings[_encodingType][5]]=RealJz;
  }
  else
  {
	  _encoder["z"]=RealJz;
  }

//  _encoder["z"]=RealJz;
  _encoder.setCellID( hit );
  hit->setType( CHT( CHT::had, CHT::hcal , _currentHCALCollectionCaloLayout,  _trueLayer ) );
//  streamlog_out( MESSAGE )   <<"getCellSize() "<<getCellSize()<<" layer " <<_trueLayer<< std::endl;
// streamlog_out( MESSAGE ) << "!!!!!!!!!! before  posB set TK!!!!!!!!!!" << std::endl
//       <<"delta_I " <<delta_I<< "_Iaxis.z() " <<_Iaxis.z()<<" delta_J "<<delta_J<<" _Jaxis.z() "<<_Jaxis.z()<<std::endl
//       <<"_Iaxis.x() "<<_Iaxis.x()<<" _Iaxis.y() "<<_Iaxis.y()<< " _Iaxis.z() "<<_Iaxis.z()<<std::endl
//       <<"_Jaxis.x() "<<_Jaxis.x()<<" _Jaxis.y() "<<_Jaxis.y()<< " _Jaxis.z() "<< _Jaxis.z()<<std::endl;
  float posB[3];
/// posB[1]=_hitPosition[1]+1.2*(delta_I*_Iaxis.y()+delta_J*_Jaxis.y());
  posB[0]=_hitPosition[0]+getCellSize()*(delta_I*_Iaxis.x()+delta_J*_Jaxis.x());
  posB[1]=_hitPosition[1]+getCellSize()*(delta_I*_Iaxis.y()+delta_J*_Jaxis.y());
  posB[2]=_hitPosition[2]+getCellSize()*(delta_I*_Iaxis.z()+delta_J*_Jaxis.z());
  hit->setPosition(posB);
//  cout<<"_hitPosition[0]= " <<_hitPosition[0]<<" RealIy= "<< RealIy <<" RealJz= " << RealJz<<endl; 
//  cout<<"posB[0]= " <<posB[0]<<" posB[1]= "<< posB[1]  <<" posB[2]= " << posB[2]<<endl; 
 // streamlog_out( MESSAGE ) << "!!! 1.2*(delta_I*_Iaxis.x()+delta_J*_Jaxis.x())!!!" << 1.2*(delta_I*_Iaxis.x()+delta_J*_Jaxis.x())<< std::endl;
//  streamlog_out( MESSAGE ) << "!!! getCellSize()*(delta_I*_Iaxis.y()+delta_J*_Jaxis.y())!!!" << getCellSize()*(delta_I*_Iaxis.y()+delta_J*_Jaxis.y())<< std::endl;
//  streamlog_out( MESSAGE ) << "!!! getCellSize()*(delta_I*_Iaxis.z()+delta_J*_Jaxis.z())!!!" << getCellSize()*(delta_I*_Iaxis.z()+delta_J*_Jaxis.z())<< std::endl;
//  streamlog_out( MESSAGE ) << "!!!!!!!!!! after posB set TK!!!!!!!!!!" << std::endl
//         		   << "  hit-posB[0] " << _hitPosition[0]-posB[0] << " hit-posB[1] " << _hitPosition[1]-posB[1] << "hit-posB[2] " << _hitPosition[2]-posB[2] << std::endl;
//          		   << " _hitPosition[0] " << _hitPosition[0] << "  _hitPosition[1] " << _hitPosition[1] << " _hitPosition[2] " << _hitPosition[2] << std::endl
}



void SimDigitalGeomCellId::setLayerLayout(CHT::Layout layout)
{
  if(_normal_I_J_setter != NULL) delete _normal_I_J_setter;

  _currentHCALCollectionCaloLayout=layout;

  if(_useGear)
  {
	  if (_currentHCALCollectionCaloLayout == CHT::endcap) 
  	  {
		  _layerLayout= & Global::GEAR->getHcalEndcapParameters().getLayerLayout();

  	      if (_geom==TESLA)
			  _normal_I_J_setter= new SimDigitalGeomRPCFrame_TESLA_ENDCAP(*this);
  	      else
			  _normal_I_J_setter= new SimDigitalGeomRPCFrame_VIDEAU_ENDCAP(*this);
  	  } 
	  else if (_currentHCALCollectionCaloLayout == CHT::ring) 
  	  {
		  _layerLayout= & Global::GEAR->getHcalRingParameters().getLayerLayout();
		  
		  // no TESLA ring ?
		  _normal_I_J_setter= new SimDigitalGeomRPCFrame_VIDEAU_ENDCAP(*this);
  	  }  
  	  else 
  	  {
		  _layerLayout= & Global::GEAR->getHcalBarrelParameters().getLayerLayout();

  	      if (_geom==TESLA)
			  _normal_I_J_setter= new SimDigitalGeomRPCFrame_TESLA_BARREL(*this);
  	      else
			  _normal_I_J_setter= new SimDigitalGeomRPCFrame_VIDEAU_BARREL(*this);
  	  }
  }
  else
  {
	  dd4hep::Detector & ild = dd4hep::Detector::getInstance();

	  if(_currentHCALCollectionCaloLayout == CHT::barrel) 
	  { 
		  const std::vector< dd4hep::DetElement>& det = dd4hep::DetectorSelector(ild).detectors(
 	            (dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL),
 	    	    (dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD) ) ;

	      caloData = det.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();
		  cout << "det size: " << det.size() << ", type: " << det.at(0).type() << endl;
	  }

	  if(_currentHCALCollectionCaloLayout == CHT::ring) 
	  { 
		  const std::vector< dd4hep::DetElement>& det = dd4hep::DetectorSelector(ild).detectors(
 	            (dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::AUXILIARY),
 	    	    dd4hep::DetType::FORWARD) ;

	      caloData = det.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

		  cout << "det size: " << det.size() << ", type: " << det.at(0).type() << endl;
	  }

	  if(_currentHCALCollectionCaloLayout == CHT::endcap) 
	  { 
		  const std::vector< dd4hep::DetElement>& det = dd4hep::DetectorSelector(ild).detectors(
 	            (dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP),
 	    	    (dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD) ) ;

	      caloData = det.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();
		  cout << "det size: " << det.size() << ", type: " << det.at(0).type() << endl;
	  }
  }
}

int SimDigitalGeomCellId::_encodingType = 0;

const string EncodingType[ENCODINGTYPES] = {"LCGEO", "MOKKA"};

void SimDigitalGeomCellId::setEncodingType(std::string type)
{
	for(int iType = 0; iType < ENCODINGTYPES; ++iType)
	{
		if(type == EncodingType[iType]) 
		{
			_encodingType = iType;
			break;
		}
	}

	std::cout << "the encoding type: " << _encodingType << std::endl;
}

void SimDigitalGeomCellId::setHcalOption(std::string hcalOption)
{
	_hcalOption = hcalOption;
}

float SimDigitalGeomCellId::getCellSize() 
{
	float cellSize = 0.;

    //cout << "get cell size start -----------------------" << endl;

	if(_useGear && _layerLayout) 
	{
		cellSize = _layerLayout->getCellSize0(_trueLayer);
	}
	else if( caloData )
	{
		const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& hcalBarrelLayers = caloData->layers;
		const double CM2MM = 10.;

		cellSize = hcalBarrelLayers[_trueLayer].cellSize0 * CM2MM;
	}

    cout << "cellSize: " << cellSize << endl;

	return cellSize; 
}


void SimDigital::remove_adjacent_step(std::vector<StepAndCharge>& vec)
{
  if (vec.size()==0) return;
  std::vector<StepAndCharge>::iterator first=vec.begin();
  std::vector<StepAndCharge>::iterator lasttobekept=vec.end();
  lasttobekept--;
  while (int(first-lasttobekept)<0)
    {
      std::vector<StepAndCharge>::iterator second=first;
      second++;
      while (int(second-lasttobekept) < 0)
	{
	  if ( ((*first).step-(*second).step).perp() > _minXYdistanceBetweenStep ) // do nothing
	    second++;
	  else //second is too close of first : second should be removed so put it at the end 
	    {
	      std::iter_swap(second,lasttobekept);
	      lasttobekept--;
	    }
	}
      if ( ((*first).step-(*lasttobekept).step).perp() <= _minXYdistanceBetweenStep ) lasttobekept--; 
      first++;
    }
  std::vector<StepAndCharge>::iterator firstToremove=lasttobekept;
  firstToremove++;
  if (_keepAtLeastOneStep && firstToremove==vec.begin()) firstToremove++;
  vec.erase(firstToremove,vec.end());
}


void SimDigital::fillTupleStep(std::vector<StepAndCharge>& vec,int level)
{
  _tupleStepFilter->fill(level,int(vec.size()));
  for (std::vector<StepAndCharge>::iterator it=vec.begin(); it != vec.end(); it++)
    {
      _debugTupleStepFilter->fill(0,level);
      _debugTupleStepFilter->fill(1,it->step.x());
      _debugTupleStepFilter->fill(2,it->step.y());
      _debugTupleStepFilter->fill(3,it->step.z());
      float minDist=20000;
      for (std::vector<StepAndCharge>::iterator itB=vec.begin(); itB != vec.end(); itB++)
	{
	  if (itB==it) continue;
	  float dist=((it->step)-(itB->step)).perp();
	  if (dist<minDist) minDist=dist;
	}
      _debugTupleStepFilter->fill(4,minDist);
      _debugTupleStepFilter->addRow();
    }
}

void SimDigital::createPotentialOutputHits(std::map<int,hitMemory>& myHitMap, LCCollection *col, SimDigitalGeomCellId& aGeomCellId )
{
  int numElements = col->getNumberOfElements();
  for (int j(0); j < numElements; ++j) 
    {
      SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
      depositedEnergyInRPC+=hit->getEnergy()/1e6;
      std::vector<StepAndCharge> steps=aGeomCellId.decode(hit);
      fillTupleStep(steps,0);


      float cellSize=aGeomCellId.getCellSize();
      multiplicityChargeSplitterBase& mult=getSplitter();
      mult.newHit(cellSize);
     

      std::vector<StepAndCharge>::iterator remPos=std::remove_if(steps.begin(), steps.end(), absZGreaterThan(_absZstepFilter) );
      if (steps.size() > 0 &&_keepAtLeastOneStep && remPos==steps.begin()) remPos++;
      steps.erase(remPos,steps.end());
      fillTupleStep(steps,1);

      float effCorr=_effMap->getEfficiency(aGeomCellId.I(), aGeomCellId.J(), aGeomCellId.K(), aGeomCellId.stave(), aGeomCellId.module());
      steps.erase(std::remove_if(steps.begin(), steps.end(),randomGreater(effCorr)), steps.end());
      fillTupleStep(steps,2);

      for (std::vector<StepAndCharge>::iterator itstep=steps.begin(); itstep != steps.end(); itstep++){
	itstep->charge=_QPolya->GetRandom();
	if(itstep->charge<0.4)streamlog_out(DEBUG) <<" "<<itstep->charge<<std::endl;  
	streamlog_out( DEBUG ) << "step at : " << itstep->step
			       << "\t with a charge of : " << itstep->charge
			       << std::endl;
      } //loop on itstep
      
      std::sort(steps.begin(), steps.end(), sortStepWithCharge );
      streamlog_out( DEBUG ) << "sim hit at " << hit << std::endl;
      if (streamlog::out.write< DEBUG >() )
	{
	  for(std::vector<StepAndCharge>::iterator it=steps.begin(); it!=steps.end(); ++it){
	    streamlog_out( DEBUG ) << "step at : " << (*it).step
				   << "\t with a charge of : " << (*it).charge 
				   << std::endl;
	  }
	}
      
      

      remove_adjacent_step(steps);
      fillTupleStep(steps,3);
      _tupleStepFilter->addRow();
      
      for (std::vector<StepAndCharge>::iterator itstep=steps.begin(); itstep != steps.end(); itstep++){
	mult.addCharge((*itstep).charge,(*itstep).step.x(),(*itstep).step.y());
      } //loop on itstep
	  
    
    
      for (std::map<multiplicityChargeSplitterBase::I_J_Coordinates,float>::const_iterator it=mult.chargeMap().begin();
	   it != mult.chargeMap().end(); it++)
	{
	  if (it->second >= 0)
	    {
	      CalorimeterHitImpl *tmp=new CalorimeterHitImpl();
	      aGeomCellId.encode(tmp,it->first.first,it->first.second);
	      
	      hitMemory &calhitMem=myHitMap[tmp->getCellID0()];	
	      if (calhitMem.ahit==0){
		calhitMem.ahit=tmp;
		calhitMem.ahit->setCellID1(hit->getCellID1());
		calhitMem.ahit->setEnergy(0); 
	      }
	      else
		delete tmp;
	      if (calhitMem.maxEnergydueToHit < it->second)
		{
		  calhitMem.rawHit=j; //for (int j(0); j < numElements; ++j)
		  calhitMem.maxEnergydueToHit=it->second;
		}
	      calhitMem.ahit->setEnergy(calhitMem.ahit->getEnergy()+it->second);
	      calhitMem.relatedHits.insert(j); //for (int j(0); j < numElements; ++j)
	    }	
	  else
	    {
	      streamlog_out(ERROR) << "BUG in charge splitter, got a non positive charge : " << it->second << std::endl;
	    }
	} //loop on added hits for this hit
	
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
    if (ilevel==0) _counters["N1"]++;
    if (ilevel==1) _counters["N2"]++;
    if (ilevel==2) _counters["N3"]++;
    currentHitMem.ahit->setEnergy(calibr_coeff);
  }
}

LCCollectionVec * SimDigital::processHCALCollection(LCCollection * col, CHT::Layout layout, LCFlagImpl& flag)
{
  LCCollectionVec *hcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  hcalcol->setFlag(flag.getFlag());
  std::map<int,hitMemory> myHitMap;

//  cout<<"LCCollectionVec * SimDigital::processHCALCollection: layout= "<< layout<< endl;
  
  SimDigitalGeomCellId g(col,hcalcol);
  g.setLayerLayout(layout);
  createPotentialOutputHits(myHitMap,col, g );
  removeHitsBelowThreshold(myHitMap,_thresholdHcal[0]);
  if (_doThresholds) applyThresholds(myHitMap); 
	  
  //Store element to output collection
  for (std::map<int,hitMemory>::iterator it=myHitMap.begin(); it!=myHitMap.end(); it++) {
    hitMemory & currentHitMem=it->second;
    if (currentHitMem.rawHit != -1) 
      {
	streamlog_out(DEBUG)  << " rawHit= " << currentHitMem.rawHit << std::endl;
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
  return hcalcol;
}

void SimDigital::processEvent( LCEvent * evt ) {

  if( isFirstEvent() ) {
    srand(_chargeSplitterRandomSeed);
    gRandom->SetSeed(_polyaRandomSeed);

    SimDigitalGeomCellId::bookTuples(this);
  }

  _counters["|ALL"]++;
  _counters["NSim"]=0;
  _counters["NReco"]=0;
  _counters["N1"]=0;
  _counters["N2"]=0;
  _counters["N3"]=0;

  // create the output collections
  _relcol = new LCCollectionVec(LCIO::LCRELATION); 
 
  /////////////////for ECAL---------------------------------------------------
  // copy the flags from the input collection
  //GG : it should be checked why we put the flag like this.
  LCFlagImpl flag;
  flag.setBit(LCIO::CHBIT_LONG);
  flag.setBit(LCIO::RCHBIT_ENERGY_ERROR);    //open the energy error flag to store the MC Truth (for easy comparison == not a eligent way!!)

  processHCAL(evt,flag);

 
  evt->addCollection(_relcol,_outputRelCollection.c_str());

  _tupleCollection->fill(0,_counters["NSim"]);
  _tupleCollection->fill(1,_counters["NReco"]);
  _tupleCollection->fill(2,_counters["N1"]);
  _tupleCollection->fill(3,_counters["N2"]);
  _tupleCollection->fill(4,_counters["N3"]);
  _tupleCollection->addRow();

  streamlog_out(MESSAGE) << "have processed " << _counters["|ALL"] << " events" << std::endl;
}

void SimDigital::check( LCEvent * evt ) { }

void SimDigital::end(){
 
}
