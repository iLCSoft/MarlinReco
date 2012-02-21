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
#include "CLHEP/Units/SystemOfUnits.h"

#include <marlin/AIDAProcessor.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>



#include <TROOT.h>
#include <TMath.h>
#include "TTree.h"
#include "TH1F.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;

//FIXME to be removed when not more needed
//#define SDHCAL_MARLINUTIL_BUGFIX 1
//#define SDHCAL_MOKKA_BUGFIX 1

SimDigital aSimDigital ;

SimDigital::SimDigital () : Processor("SimDigital"), _chargeSplitterUniform(), _chargeSplitterFunction() {
   _description = "the transfer Between Energy and Induced Charge for SDHCAL" ;
 
#ifdef SIMDIGITAL_WITHECAL
   registerECALparameters(); 
 #endif 
  
   
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
   hcalThresholds.push_back(0.00004);
   registerProcessorParameter("HCALThreshold" ,
                             "Threshold for HCAL Hits in pC" ,
			      _thresholdHcal,
			      hcalThresholds);
   
   std::vector<float> calibrHcal;
   calibrHcal.push_back(34.8);
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
                               "Distance (in mm) from the border of the cell that implies hit on neighbouring cell for 'Uniform' charge splitter : used if UseFunctionForChargeSplitting is false"  ,
                               _chargeSplitterUniform._edgeSize,
                               (float)2.0 ) ;

   
   registerProcessorParameter( "functionFormula" ,
			       "ROOT TF2 formula for function describing the induced charge distribution, x and y are coordinates from the step position in I and J direction : : used if UseFunctionForChargeSplitting is true",
			       _chargeSplitterFunction._function_description,
			       std::string("1/cosh([0]*sqrt(x*x+y*y))") ); 
   
   registerProcessorParameter( "functionRange",
			       "maximal distance (in mm) at which a step can induce charge using the 2D function defined with functionFormula",
			       _chargeSplitterFunction._range,
			       (float) 30 );

   std::vector<float> chadisParameter;
   chadisParameter.push_back(0.24);
   registerProcessorParameter( "functionParameters",
			       "parameter values for the function defined with functionFormula",
			       _chargeSplitterFunction._functionParameters,
			       chadisParameter );

   registerProcessorParameter( "RPC_PadSeparation",
			       "distance in mm between two RPC pads : used if UseFunctionForChargeSplitting is true ", 
			       _chargeSplitterFunction._RPC_PadSeparation,
			       (float) 0.16 );
   

   _splitChargeWithFunction=true;
   registerProcessorParameter( "UseFunctionForChargeSplitting",
			       "If true, use function to split the charge between pads. If false split uniformely",
			       _splitChargeWithFunction,
			       true );

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

   registerProcessorParameter( "PolyaAmplitude",
			       "Parameter for the Polya distribution used to simulate the induced charge distribution : amplitude or normalisation of the distribution", 
			       _polyaAmplitude,
			       (double) 1 );

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

  //a=1.05,b=-2.1,c=1.05,d=-1.05;//for 1 neighbor pad fired
  //previousLayer=1 ,previousmodule=1, previousstave=1,previousIy=1,previousJz=1;
   //define a polya function to get induced charge
   _QPolya = new TF1("QPolya",Polya,0.0,10.0,3);
   _QPolya->SetParameters(_polyaAverageCharge,_polyaFunctionWidthParameter,_polyaAmplitude);
   //define the function for charge distribution
   //c2d=new TF2("c2d",ChaDis,-10,10,-10,10,3);
   //c2d->SetParNames("width","x'value","y'value");

   _chargeSplitterFunction.init();

   //fg: need to set default encoding in for reading old files...
   CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;

#ifdef SIMDIGITAL_WITHECAL
   setECALgeometry();
#endif 

   //assure SDHCAL _thresholdHcal are in increasing order
   std::sort(_thresholdHcal.begin(),_thresholdHcal.end());

   //book tuples
   _debugTupleStepFilter  = AIDAProcessor::tupleFactory( this )->create("SimDigitalStepDebug",
									"SimDigital_StepDebug",
									"int filterlevel, float deltaI,deltaJ,deltaLayer,minIJdist"); 
   std::cout << "Tuple for step debug has been initialized to " << _debugTupleStepFilter << std::endl;
   std::cout << "it has " << _debugTupleStepFilter->columns() << " columns" <<std::endl;


   _tupleStepFilter   = AIDAProcessor::tupleFactory( this )->create("SimDigitalStepStat",
								    "SimDigital_StepStat",
								    "int allsteps, absZfiltered, IJdistancefiltered");
   std::cout << "Tuple for step stat has been initialized to " << _tupleStepFilter << std::endl;
   std::cout << "it has " << _tupleStepFilter->columns() << " columns" <<std::endl;


}
 

void SimDigital::processRunHeader( LCRunHeader* run) {}


void SimDigital::processHCAL(LCEvent* evt, LCFlagImpl& flag)
{
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    try{
      std::string colName =  _hcalCollections[i] ;    
      LCCollection * col = evt->getCollection( colName.c_str() ) ;
      CHT::Layout layout = layoutFromString( colName ); 
#ifdef SDHCAL_MARLINUTIL_BUGFIX 
      //FIXME should correct calorimeterHitType.cc in MarlinUtil instead
      if (colName=="HcalEndCapRingsCollection") layout=CHT::ring;
#endif
      LCCollectionVec *hcalcol = processHCALCollection(col,layout,flag);
      evt->addCollection(hcalcol,_outputHcalCollections[i].c_str());
    }
    catch(DataNotAvailableException &e){  
    }     
  }
}



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
}

void multiplicityChargeSplitterFunction:: addCharge(float charge, float pos_I, float pos_J)
{
  if (_RPC_PadSeparation<_cellSize) return;
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
	  _chargeMap[I_J_Coordinates(I,J)]+=charge*_f2->Integral(distI,distII,distJ,distJJ)/_normalisation;
	}
    }
}


AIDA::ITuple* SimDigitalGeomCellId::_tupleHit = NULL;
AIDA::ITuple* SimDigitalGeomCellId::_tupleStep = NULL;

void SimDigitalGeomCellId::bookTuples(const marlin::Processor* proc)
{
  _tupleHit  = AIDAProcessor::tupleFactory( proc )->create("SimDigitalGeom",
							   "SimDigital_Debug",
							   "int detector,chtlayout,module,stave,layer,I,J, float x,y,z, normalx,normaly,normalz, Ix,Iy,Iz,Jx,Jy,Jz"); 
  std::cout << "Tuple for Hit has been initialized to " << _tupleHit << std::endl;
  std::cout << "it has " << _tupleHit->columns() << " columns" <<std::endl;
  
  _tupleStep = AIDAProcessor::tupleFactory( proc )->create("SimDigitalStep",
							   "SimDigital_DebugStep",
							   "int detector,chtlayout,hitcellid,nstep, float hitx,hity,hitz,stepx,stepy,stepz,deltaI,deltaJ,deltaLayer");
  std::cout << "Tuple for Step has been initialized to " << _tupleStep << std::endl;
  std::cout << "it has " << _tupleStep->columns() << " columns" <<std::endl;

}


SimDigitalGeomCellId::SimDigitalGeomCellId(LCCollection *inputCol, LCCollectionVec *outputCol) 
  : _hitPosition(NULL),
    _decoder(inputCol), _encoder(inputCol->getParameters().getStringVal(LCIO::CellIDEncoding),outputCol),
    _layerLayout( & Global::GEAR->getHcalBarrelParameters().getLayerLayout() ),
    _normal_I_J_setter(new SimDigitalGeomRPCFrame_TESLA_BARREL(*this)),
    _currentHCALCollectionCaloLayout(CHT::any),
    _normal(),_Iaxis(),_Jaxis()
{
  _trueLayer=_stave=_module=_Iy=_Jz=-999;
  outputCol->parameters().setValue(LCIO::CellIDEncoding,inputCol->getParameters().getStringVal(LCIO::CellIDEncoding));
  _geom=TESLA;
  try
    { 
      Global::GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_order"); //it is TESLA geometry if it is OK
    }
  catch (gear::Exception &)
    {
      _geom=VIDEAU;
    }
      
}

SimDigitalGeomCellId::~SimDigitalGeomCellId()
{
  delete _normal_I_J_setter;
}

void SimDigitalGeomRPCFrame_TESLA_BARREL::setRPCFrame()
{
  normal().set(0,1,0);
  Iaxis().set(1,0,0);
  Jaxis().set(0,0,1);
  double angle=(stave()/2)*(-45)*CLHEP::degree;
  Iaxis().rotateZ(angle);
  normal().rotateZ(angle);
}
void SimDigitalGeomRPCFrame_VIDEAU_BARREL::setRPCFrame()
{
  normal().set(0,1,0);
  Iaxis().set(-1,0,0);
  Jaxis().set(0,0,1);
  double angle=stave()*45*CLHEP::degree;
  Iaxis().rotateZ(angle);
  normal().rotateZ(angle);
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
      std::cout << "ERROR ; TESLA detector : unknown module for endcap " << module() << std::endl;
    }      
}

//valid also for all endcap rings (VIDEAU and TESLA)
void SimDigitalGeomRPCFrame_VIDEAU_ENDCAP::setRPCFrame()
{
  if (module()==6)
    {
      normal().set(0,0,1);
      Iaxis().set(0,1,0);
      Jaxis().set(1,0,0);
      Iaxis().rotateZ(stave()*(-90)*CLHEP::degree);
      Jaxis().rotateZ(stave()*(-90)*CLHEP::degree);
    }
  else if (module()==0)
    {
      normal().set(0,0,-1);
      Iaxis().set(0,1,0);
      Jaxis().set(-1,0,0);
      Iaxis().rotateZ(stave()*90*CLHEP::degree);
      Jaxis().rotateZ(stave()*90*CLHEP::degree);
    }
  else
    {
      std::cout << "ERROR : unknown module for endcap or endcapring " << module() << std::endl;
    }
}



std::vector<LCVector3D> SimDigitalGeomCellId::decode(SimCalorimeterHit *hit)
{
  _trueLayer = _decoder( hit )["K-1"]; // + 1;
  _stave     = _decoder( hit )["S-1"];
  _module    = _decoder( hit )["M"];
  _Iy        = _decoder( hit )["I"];
  _Jz        = _decoder( hit )["J"];
  _hitPosition = hit->getPosition();
  

  _normal_I_J_setter->setRPCFrame();

  std::vector<LCVector3D> stepsInIJZcoord;
  LCVector3D hitpos;
  if (NULL != _hitPosition) hitpos.set(_hitPosition[0],_hitPosition[1],_hitPosition[2]);
  for (int imcp=0; imcp<hit->getNMCContributions(); imcp++)
    {
      LCVector3D stepposvec;
      const float* steppos=hit->getStepPosition(imcp);
      if (NULL != steppos) stepposvec.set(steppos[0],steppos[1],steppos[2]);
#ifdef SDHCAL_MOKKA_BUGFIX 
      //FIXME : Mokka bug
      if (_geom==VIDEAU && _currentHCALCollectionCaloLayout>=2 && hitpos.z()>0) stepposvec.set(-stepposvec.x(),stepposvec.y(),-stepposvec.z());
#endif
      if (stepposvec.mag2() != 0)
	{
	  stepposvec-=hitpos;
	  stepsInIJZcoord.push_back( LCVector3D(stepposvec*_Iaxis,stepposvec*_Jaxis,stepposvec*_normal) );
	}
      else
	streamlog_out(WARNING) << "DIGITISATION : STEP POSITION IS (0,0,0)" << std::endl;
    }

  //if no steps have been found, then put one step at the center of the cell :
  if (stepsInIJZcoord.size()==0) stepsInIJZcoord.push_back(LCVector3D(0,0,0));
  

  //these tuples are for debugging geometry aspects
  if (_tupleHit!=NULL) 
    {
      _tupleHit->fill(TH_DETECTOR,int(_geom));
      _tupleHit->fill(TH_CHTLAYOUT,int(_currentHCALCollectionCaloLayout));
      _tupleHit->fill(TH_MODULE,_module);
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
		  _tupleStep->fill(TS_DELTAI+i,stepsInIJZcoord[imcp][i]);
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
  _encoder["M"]=_module;
  _encoder["S-1"]=_stave;
  int RealIy=_Iy+delta_I;
  if (RealIy<0||RealIy>330)RealIy=0; //FIXME the 330 value should depend on the cellSize and on the Layer
  _encoder["I"]=RealIy;
  int RealJz=_Jz+delta_J;
  if(RealJz<0||RealJz>330)RealJz=0; //FIXME the 330 value should depend on the cellSize and on the Layer
  _encoder["J"]=RealJz;
  _encoder["K-1"]=_trueLayer;
  _encoder.setCellID( hit );
  hit->setType( CHT( CHT::had, CHT::hcal , _currentHCALCollectionCaloLayout,  _trueLayer ) );

  float posB[3];
  posB[0]=_hitPosition[0]+getCellSize()*(delta_I*_Iaxis.x()+delta_J*_Jaxis.x());
  posB[1]=_hitPosition[1]+getCellSize()*(delta_I*_Iaxis.y()+delta_J*_Jaxis.y());
  posB[2]=_hitPosition[2]+getCellSize()*(delta_I*_Iaxis.z()+delta_J*_Jaxis.z());
  hit->setPosition(posB);
}


void SimDigitalGeomCellId::setLayerLayout(CHT::Layout layout)
{
  delete _normal_I_J_setter;
  _currentHCALCollectionCaloLayout=layout;
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

float SimDigitalGeomCellId::getCellSize() {return _layerLayout->getCellSize0(_trueLayer); }



void SimDigital::remove_adjacent_step(std::vector<LCVector3D>& vec)
{
  if (vec.size()==0) return;
  std::vector<LCVector3D>::iterator first=vec.begin();
  std::vector<LCVector3D>::iterator lasttobekept=vec.end();
  lasttobekept--;
  while (int(first-lasttobekept)<0)
    {
      std::vector<LCVector3D>::iterator second=first;
      second++;
      while (int(second-lasttobekept) < 0)
	{
	  if ( ((*first)-(*second)).perp() > _minXYdistanceBetweenStep ) // do nothing
	    second++;
	  else //second is too close of first : second should be removed so put it at the end 
	    {
	      std::iter_swap(second,lasttobekept);
	      lasttobekept--;
	    }
	}
      if ( ((*first)-(*lasttobekept)).perp() <= _minXYdistanceBetweenStep ) lasttobekept--; 
      first++;
    }
  std::vector<LCVector3D>::iterator firstToremove=lasttobekept;
  firstToremove++;
  if (_keepAtLeastOneStep && firstToremove==vec.begin()) firstToremove++;
  vec.erase(firstToremove,vec.end());
}


void SimDigital::fillTupleStep(std::vector<LCVector3D>& vec,int level)
{
  _tupleStepFilter->fill(level,int(vec.size()));
  for (std::vector<LCVector3D>::iterator it=vec.begin(); it != vec.end(); it++)
    {
      _debugTupleStepFilter->fill(0,level);
      _debugTupleStepFilter->fill(1,it->x());
      _debugTupleStepFilter->fill(2,it->y());
      _debugTupleStepFilter->fill(3,it->z());
      float minDist=20000;
      for (std::vector<LCVector3D>::iterator itB=vec.begin(); itB != vec.end(); itB++)
	{
	  if (itB==it) continue;
	  float dist=((*it)-(*itB)).perp();
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
      std::vector<LCVector3D> steps=aGeomCellId.decode(hit);
      fillTupleStep(steps,0);

      float cellSize=aGeomCellId.getCellSize();
      multiplicityChargeSplitterBase& mult=getSplitter();
      mult.newHit(cellSize);
      
      std::vector<LCVector3D>::iterator remPos=std::remove_if(steps.begin(), steps.end(), absZGreaterThan(_absZstepFilter) );
      if (steps.size() > 0 &&_keepAtLeastOneStep && remPos==steps.begin()) remPos++;
      steps.erase(remPos,steps.end());
      fillTupleStep(steps,1);
      remove_adjacent_step(steps);
      fillTupleStep(steps,2);
      _tupleStepFilter->addRow();
      

      for (std::vector<LCVector3D>::iterator itstep=steps.begin(); itstep != steps.end(); itstep++)
	{      	  
	  float Minduced=_QPolya->GetRandom();
	  if(Minduced<0.4)streamlog_out(DEBUG) <<" "<<Minduced<<std::endl;   
	  mult.addCharge(Minduced,itstep->x(),itstep->y());
	} //loop on itstep
    
    
      for (std::map<multiplicityChargeSplitterBase::I_J_Coordinates,float>::const_iterator it=mult.chargeMap().begin();
	   it != mult.chargeMap().end(); it++)
	{
	  if (it->second > 0)
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
    currentHitMem.ahit->setEnergy(calibr_coeff);
  }
}

LCCollectionVec * SimDigital::processHCALCollection(LCCollection * col, CHT::Layout layout, LCFlagImpl& flag)
{
  LCCollectionVec *hcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  hcalcol->setFlag(flag.getFlag());
  std::map<int,hitMemory> myHitMap;

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
  return hcalcol;
}


void SimDigital::processEvent( LCEvent * evt ) {

  if( isFirstEvent() ) {
    SimDigitalGeomCellId::bookTuples(this);
      }

  _counters["|ALL"]++;

   // create the output collections
  _relcol = new LCCollectionVec(LCIO::LCRELATION); 
 
  /////////////////for ECAL---------------------------------------------------
  // copy the flags from the input collection
  //GG : it should be checked why we put the flag like this.
  LCFlagImpl flag;
  flag.setBit(LCIO::CHBIT_LONG);
  flag.setBit(LCIO::RCHBIT_ENERGY_ERROR);    //open the energy error flag to store the MC Truth (for easy comparison == not a eligent way!!)
 #ifdef SIMDIGITAL_WITHECAL
 processECAL(evt,flag);
#endif 
  processHCAL(evt,flag);

 
  evt->addCollection(_relcol,_outputRelCollection.c_str());
}

void SimDigital::check( LCEvent * evt ) { }

void SimDigital::end(){
 
}

#ifdef SIMDIGITAL_WITHECAL


//ECAL specific code, stolen from NewLDCCaloDigi at some time in the past

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



void SimDigital::registerECALparameters()
{
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
  
  /////////ECAL OUTPUT
  
  _outputEcalCollections.push_back(std::string("ECALBarrel"));
  _outputEcalCollections.push_back(std::string("ECALEndcap"));
  _outputEcalCollections.push_back(std::string("ECALOther"));
  _outputEcalCollections.push_back(std::string("ECALBarrelPreShower"));
  _outputEcalCollections.push_back(std::string("ECALEndcapPreShower"));
  _outputEcalCollections.push_back(std::string("ECALOtherPreShower"));

  ///////////////////////////////////////////////////////////
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


}

void SimDigital::setECALgeometry()
{
  const float pi = acos(-1.0);
  const float twopi = 2.0*pi;

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

}

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

#endif 
