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

#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

#include <iostream>
#include <string>
#include <algorithm>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "CalorimeterHitType.h"

#include <marlin/AIDAProcessor.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>

#include <DDSegmentation/BitField64.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/DetectorSelector.h>

using namespace lcio ;
using namespace marlin ;
using namespace std;


SimDigital aSimDigital ;


SimDigital::SimDigital()
	: Processor("SimDigital") ,
	  chargeSpreaderParameters()
{
	_description = "This processor creates SDHCAL digitized CalorimeterHits from SDHCAL SimCalorimeterHits" ;


	std::vector<std::string> hcalCollections = { "HcalBarrelCollection" , "HcalEndCapRingsCollection" , "HcalEndCapsCollection" } ;
	registerInputCollections( LCIO::SIMCALORIMETERHIT,
							  "inputHitCollections" ,
							  "Sim Calorimeter Hit Collections" ,
							  _inputCollections ,
							  hcalCollections) ;


	std::vector<std::string> outputHcalCollections = { "HCALBarrelDigi" , "HCALEndcapDigi" , "HCALOtherDigi" } ;
	registerProcessorParameter( "outputHitCollections",
								"output hit collection names",
								_outputCollections,
								outputHcalCollections) ;

	std::vector<std::string> outputRelCollections = {} ;
	registerProcessorParameter( "outputRelationCollections",
								"output hit relation Collection Names" ,
								_outputRelCollections ,
								outputRelCollections ) ;


	std::vector<float> hcalThresholds = {0.1f} ;
	registerProcessorParameter("HCALThreshold" ,
							   "Threshold for HCAL Hits in pC" ,
							   _thresholdHcal,
							   hcalThresholds) ;

	registerProcessorParameter( "HCALCellSize" ,
								"Cell size (mm) of HCAL, if it is equal or less than zero then the value is taken from dd4hep" ,
								_cellSize ,
								0.0f ) ;


	registerProcessorParameter("EffMapOption" ,
							   "Step efficiency correction method : should be Uniform" ,
							   efficiencyOption ,
							   std::string("Uniform") ) ;

	registerProcessorParameter("EffMapConstantValue",
							   "Value of the constant term for efficiency correction if EffMapOption==Uniform",
							   _constEffMapValue,
							   0.97f) ;

	registerProcessorParameter("EffMapFile" ,
							   "Efficiency map file" ,
							   effMapFile ,
							   std::string("")) ;


	//charge spreader parameters
	registerProcessorParameter("SpreaderMapFile",
							   "Charge spreader map file",
							   spreaderMapFile ,
							   std::string("")) ;

	registerProcessorParameter( "functionRange" ,
								"maximal distance (in mm) at which a step can induce charge using the 2D function defined with functionFormula or when using ChargeSplitterOption==Erf",
								chargeSpreaderParameters.range ,
								30.0f ) ;


	registerProcessorParameter( "RPC_PadSeparation",
								"distance in mm between two RPC pads : used if ChargeSplitterOption==Function or Erf",
								chargeSpreaderParameters.padSeparation ,
								0.0f ) ;

	std::vector<float> erfWidth = {2} ;
	registerProcessorParameter( "erfWidth",
								"Width values for the different Erf functions",
								chargeSpreaderParameters.erfWidth ,
								erfWidth ) ;

	std::vector<float> erfWeigth = {1} ;
	registerProcessorParameter( "erfWeigth",
								"Weigth for the different Erf functions",
								chargeSpreaderParameters.erfWeigth ,
								erfWeigth ) ;


	registerProcessorParameter( "ChargeSplitterd",
								"d parameter for exact splitter",
								chargeSpreaderParameters.d,
								1.0f ) ;

	registerProcessorParameter( "ChargeSplitterOption",
								"Define the charge splitter method. Possible option : Erf , Exact",
								chargeSpreaderOption,
								std::string("Erf") ) ;




	registerProcessorParameter( "CellIDEncodingStringType",
								"The type of the encoding, LCGEO or PROTO",
								_encodingType,
								std::string("LCGEO")) ;



	registerProcessorParameter("doThresholds",
							   "Replace analog hit energy by value given in CalibrHCAL according to thresholds given in HCALThreshold",
							   _doThresholds,
							   true) ;


	registerProcessorParameter( "PolyaOption" ,
								"Uniform polya or different polya per Asic",
								polyaOption,
								std::string("Uniform") ) ;

	registerProcessorParameter( "PolyaMapFile" ,
								"Polya map file",
								polyaMapFile,
								std::string("") ) ;

	registerProcessorParameter( "PolyaRandomSeed",
								"The seed of the polya function",
								_polyaRandomSeed ,
								1 ) ;

	registerProcessorParameter( "PolyaAverageCharge" ,
								"Parameter for the Polya distribution used to simulate the induced charge distribution : mean of the distribution",
								polyaQbar ,
								1.6f ) ;

	registerProcessorParameter( "PolyaWidthParameter" ,
								"Parameter for the Polya distribution used to simulate the induced charge distribution : related to the distribution width ",
								polyaTheta ,
								16.3f ) ;


	registerProcessorParameter( "GasGapWidth" ,
								"Width of the RPC gas gap",
								_gasGapWidth ,
								1.2f ) ;



	registerProcessorParameter( "LinkSteps" ,
								"Parameter for angle correction",
								_linkSteps ,
								false ) ;


	registerProcessorParameter( "AngleCorrectionPower" ,
								"Parameter for angle correction",
								_angleCorrPow ,
								0.4f ) ;

	registerProcessorParameter( "TimeCut" ,
								"Time cut",
								timeCut ,
								std::numeric_limits<double>::max() ) ;

	registerProcessorParameter( "StepLengthCut" ,
								"Step length cut",
								stepLengthCut ,
								-1.0 ) ;





	registerProcessorParameter( "StepCellCenterMaxDistanceLayerDirection",
								"Maximum distance (mm) between the Geant4 step position and the cell center, in the RPC width direction, to keep a step for digitization",
								_absZstepFilter,
								0.0005f ) ;

	registerProcessorParameter( "StepsMinDistanceRPCplaneDirection",
								"Minimum distance (mm) between 2 Geant4 steps, in the RPC plane, to keep the 2 steps",
								_minXYdistanceBetweenStep,
								0.5f ) ;

	registerProcessorParameter( "KeepAtLeastOneStep",
								"if true, ensure that each hit will keep at least one step for digitisation independatly of filtering conditions (StepCellCenterMaxDistanceLayerDirection)",
								_keepAtLeastOneStep,
								true ) ;
}

void SimDigital::init()
{
	printParameters() ;
	// check that number of input and output collections names are the same
	assert ( _outputCollections.size() == _inputCollections.size() ) ;
	assert ( _outputRelCollections.size() == _inputCollections.size() ) ;
	assert ( _encodingType == std::string("LCGEO") || _encodingType == std::string("PROTO") ) ;

	//init charge inducer
	if ( polyaOption == std::string("Uniform") )
		chargeInducer = new UniformPolya(polyaQbar , polyaTheta) ;
	else if ( polyaOption == std::string("PerAsic") )
		chargeInducer = new AsicPolya(polyaQbar , polyaTheta , polyaMapFile) ;
	else
		throw ParseException( chargeSpreaderOption + std::string(" option for charge inducing is not available ") ) ;

	chargeInducer->setSeed(static_cast<unsigned int>(_polyaRandomSeed) ) ;
	srand( static_cast<unsigned int>(_polyaRandomSeed) ) ;


	//init charge spreader
	if (chargeSpreaderOption == "Erf")
		chargeSpreader = new GaussianSpreader ;
	else if (chargeSpreaderOption == "Exact")
		chargeSpreader = new ExactSpreader ;
	else if (chargeSpreaderOption == "ExactPerAsic")
		chargeSpreader = new ExactSpreaderPerAsic(spreaderMapFile) ;
	else
		throw ParseException( chargeSpreaderOption + std::string(" option for charge splitting is not available ") ) ;

	chargeSpreader->setParameters( chargeSpreaderParameters ) ;
	chargeSpreader->init() ;


	//init efficiency manager
	if (efficiencyOption == "Uniform")
		efficiency = new UniformEfficiency(_constEffMapValue) ;
	else if (efficiencyOption == "PerAsic")
		efficiency = new AsicEfficiency(effMapFile , _constEffMapValue) ;
	else
		throw ParseException( efficiencyOption + std::string(" option for efficiency correction is not available") ) ;



	//assure SDHCAL _thresholdHcal are in increasing order
	std::sort( _thresholdHcal.begin() , _thresholdHcal.end() ) ;

	//book tuples
	_debugTupleStepFilter  = AIDAProcessor::tupleFactory( this )->create("SimDigitalStepDebug",
																		 "SimDigital_StepDebug",
																		 "int filterlevel, float stepTime,deltaI,deltaJ,deltaLayer,minIJdist,length,charge");
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

	_histoCellCharge = AIDAProcessor::histogramFactory( this )->createHistogram1D("CellCharge","CellCharge",10000,0,100) ;

}

void SimDigital::removeAdjacentStep(std::vector<StepAndCharge>& vec)
{
	if ( vec.size() == 0 )
		return ;
	std::vector<StepAndCharge>::iterator first = vec.begin() ;
	std::vector<StepAndCharge>::iterator lasttobekept = vec.end() ;
	lasttobekept-- ;

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
		if ( ((*first).step-(*lasttobekept).step).perp() <= _minXYdistanceBetweenStep )
			lasttobekept--;
		first++;
	}
	std::vector<StepAndCharge>::iterator firstToremove=lasttobekept;
	firstToremove++;
	if (_keepAtLeastOneStep && firstToremove==vec.begin())
		firstToremove++;
	vec.erase(firstToremove,vec.end());
}


void SimDigital::fillTupleStep(const std::vector<StepAndCharge>& vec,int level)
{
	_tupleStepFilter->fill(level,int(vec.size()));
	for (std::vector<StepAndCharge>::const_iterator it=vec.begin(); it != vec.end(); it++)
	{
		_debugTupleStepFilter->fill(0,level);
		_debugTupleStepFilter->fill(1,it->time) ;
		_debugTupleStepFilter->fill(2,it->step.x()) ;
		_debugTupleStepFilter->fill(3,it->step.y()) ;
		_debugTupleStepFilter->fill(4,it->step.z()) ;
		float minDist=20000;
		for (std::vector<StepAndCharge>::const_iterator itB=vec.begin(); itB != vec.end(); itB++)
		{
			if (itB == it)
				continue ;
			float dist = ( (it->step)-(itB->step) ).perp() ;
			if (dist < minDist)
				minDist=dist ;
		}
		_debugTupleStepFilter->fill(5,minDist);
		_debugTupleStepFilter->fill(6,it->stepLength) ;
		_debugTupleStepFilter->fill(7,it->charge) ;
		_debugTupleStepFilter->addRow() ;
	}
}


SimDigital::cellIDHitMap SimDigital::createPotentialOutputHits(LCCollection* col, SimDigitalGeomCellId* aGeomCellId)
{
	cellIDHitMap myHitMap ;

	int numElements = col->getNumberOfElements() ;

	for (int j = 0 ; j < numElements ; ++j )
	{
		SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;

		std::vector<StepAndCharge> steps ;

		steps = aGeomCellId->decode(hit , _linkSteps) ;

		fillTupleStep(steps,0) ;

		float cellSize = aGeomCellId->getCellSize() ;
		chargeSpreader->newHit(cellSize) ;


		auto timeGreaterThan = [&](const StepAndCharge& v) -> bool { return std::fabs( v.time ) > timeCut ; } ;
		std::vector<StepAndCharge>::iterator remPos = std::remove_if(steps.begin() , steps.end() , timeGreaterThan ) ;
		steps.erase( remPos , steps.end() ) ;
		fillTupleStep(steps,1) ;


		auto stepSmallerThan = [&](const StepAndCharge& v) -> bool { return v.stepLength < stepLengthCut ; } ;
		remPos = std::remove_if( steps.begin() , steps.end() , stepSmallerThan ) ;
		steps.erase( remPos , steps.end() ) ;


		fillTupleStep(steps,2) ;

		auto absZGreaterThan = [&](const StepAndCharge& v) -> bool { return std::abs( v.step.z() ) > _absZstepFilter ; } ;
		remPos = std::remove_if(steps.begin() , steps.end() , absZGreaterThan ) ;

		if (steps.size() > 0 &&_keepAtLeastOneStep && remPos == steps.begin() )
			remPos++ ;
		steps.erase( remPos , steps.end() ) ;

		float eff = efficiency->getEfficiency(aGeomCellId) ;
		auto randomGreater = [eff](const StepAndCharge&) -> bool { return static_cast<double>(rand())/RAND_MAX > eff ; } ;
		steps.erase( std::remove_if(steps.begin() , steps.end() , randomGreater ) , steps.end() ) ;
		fillTupleStep(steps,3) ;


		float invGasGapWidth = 1.f/_gasGapWidth ;
		for ( auto& itstep : steps )
		{
			float angleCorr = 1 ;
			if ( itstep.stepLength*invGasGapWidth > 1 )
				angleCorr = std::pow( itstep.stepLength*invGasGapWidth , _angleCorrPow ) ;

			itstep.charge = chargeInducer->getCharge(aGeomCellId)*angleCorr ;

			streamlog_out( DEBUG ) << "step at : " << itstep.step << "\t with a charge of : " << itstep.charge << std::endl ;
		}


		auto sortStepWithCharge = [](const StepAndCharge& s1 , const StepAndCharge& s2) -> bool { return s1.charge > s2.charge ; } ;
		std::sort(steps.begin(), steps.end(), sortStepWithCharge ) ;

		streamlog_out( DEBUG ) << "sim hit at " << hit << std::endl ;
		if (streamlog::out.write< DEBUG >() )
		{
			for(std::vector<StepAndCharge>::iterator it=steps.begin(); it!=steps.end(); ++it)
				streamlog_out( DEBUG ) << "step at : " << (*it).step << "\t with a charge of : " << (*it).charge << std::endl;

		}

		removeAdjacentStep(steps) ;
		fillTupleStep(steps,4) ;
		_tupleStepFilter->addRow() ;

		float time = std::numeric_limits<float>::max() ;
		for ( const auto& step : steps )
			time = std::min(time , step.time) ;


		for ( const StepAndCharge& itstep : steps )
			chargeSpreader->addCharge( itstep.charge , static_cast<float>(itstep.step.x()) , static_cast<float>(itstep.step.y()) , aGeomCellId ) ;


		for ( const std::pair<ChargeSpreader::I_J_Coordinates,float>& it : chargeSpreader->getChargeMap() )
		{
			if (it.second >= 0)
			{
				std::unique_ptr<CalorimeterHitImpl> tmp = aGeomCellId->encode(it.first.first , it.first.second) ;

				if (tmp == nullptr)
					continue ;

				dd4hep::long64 index = tmp->getCellID1() ;
				index = index << 32 ;
				index += tmp->getCellID0() ;

				if ( myHitMap.find(index) == myHitMap.end() ) //create hit
				{
					hitMemory& toto = myHitMap[index] ;
					toto.ahit = std::move(tmp) ;
					toto.ahit->setEnergy(0) ;
					toto.ahit->setTime(time) ;
				}

				hitMemory& calhitMem = myHitMap.at(index) ;

				if (calhitMem.maxEnergydueToHit < it.second)
				{
					calhitMem.rawHit = j ;
					calhitMem.maxEnergydueToHit = it.second ;
				}

				calhitMem.ahit->setEnergy( calhitMem.ahit->getEnergy() + it.second ) ;
				calhitMem.relatedHits.insert(j) ;
			}
			else
			{
				streamlog_out(ERROR) << "BUG in charge splitter, got a non positive charge : " << it.second << std::endl ;
			}
		} //loop on added hits for this hit

	} // end of for (int j(0); j < numElements; ++j)  //loop on elements in collection

	for( const auto& it : myHitMap )
		_hitCharge.push_back( it.second.ahit->getEnergy() ) ;

	return myHitMap ;
}


void SimDigital::removeHitsBelowThreshold(cellIDHitMap& myHitMap, float threshold)
{
	for ( auto it = myHitMap.cbegin() ; it != myHitMap.cend() ; )
	{
		if ( (it->second).ahit->getEnergy() < threshold )
			it = myHitMap.erase(it) ;
		else
			++it ;
	}
}


void SimDigital::applyThresholds(cellIDHitMap& myHitMap)
{
	for (cellIDHitMap::iterator it = myHitMap.begin() ; it != myHitMap.end() ; it++)
	{
		hitMemory& currentHitMem = it->second ;
		float hitCharge = currentHitMem.ahit->getEnergy() ;

		unsigned int iThr = 0 ;
		for ( unsigned int i = 0 ; i < _thresholdHcal.size() ; ++i )
		{
			if ( hitCharge >= _thresholdHcal.at(i) )
				iThr = i ;
		}

		if (iThr == 0)
			_counters["N1"]++ ;
		if (iThr == 1)
			_counters["N2"]++ ;
		if (iThr == 2)
			_counters["N3"]++ ;

		currentHitMem.ahit->setEnergy( static_cast<float>( iThr+1 ) ) ;
	}
}

void SimDigital::processCollection(LCCollection* inputCol , LCCollectionVec*& outputCol , LCCollectionVec*& outputRelCol , CHT::Layout layout, LCFlagImpl& flag)
{
	outputCol = new LCCollectionVec(LCIO::CALORIMETERHIT) ;
	outputRelCol = new LCCollectionVec(LCIO::LCRELATION) ;

	outputCol->setFlag(flag.getFlag()) ;

	SimDigitalGeomCellId* geomCellId = nullptr ;

	if ( _encodingType == std::string("LCGEO") )
		geomCellId = new SimDigitalGeomCellIdLCGEO(inputCol,outputCol) ;
	else if ( _encodingType == std::string("PROTO") )
		geomCellId = new SimDigitalGeomCellIdPROTO(inputCol,outputCol) ;

	geomCellId->setCellSize(_cellSize) ;

	geomCellId->setLayerLayout(layout) ;
	cellIDHitMap myHitMap = createPotentialOutputHits(inputCol , geomCellId) ;

	removeHitsBelowThreshold(myHitMap , _thresholdHcal.at(0) ) ;

	if (_doThresholds)
		applyThresholds(myHitMap) ;

	//Store element to output collection
	for (cellIDHitMap::iterator it = myHitMap.begin() ; it != myHitMap.end() ; it++)
	{
		hitMemory& currentHitMem = it->second ;
		if (currentHitMem.rawHit != -1)
		{
			streamlog_out(DEBUG) << " rawHit= " << currentHitMem.rawHit << std::endl ;
			SimCalorimeterHit* hitraw = dynamic_cast<SimCalorimeterHit*>( inputCol->getElementAt( currentHitMem.rawHit ) ) ;
			currentHitMem.ahit->setRawHit(hitraw) ;
		}

		auto caloHit = currentHitMem.ahit.release() ;
		outputCol->addElement(caloHit) ;

		//put only one relation with the SimCalorimeterHit which contributes most
		SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>( inputCol->getElementAt( currentHitMem.rawHit ) ) ;
		LCRelationImpl* rel = new LCRelationImpl(hit , caloHit , 1.0) ;
		outputRelCol->addElement( rel ) ;

	} //end of loop on myHitMap

	delete geomCellId ;
}

void SimDigital::processEvent( LCEvent* evt )
{
	if( isFirstEvent() )
		SimDigitalGeomCellId::bookTuples(this) ;

	_counters["|ALL"]++;
	_counters["NSim"]=0;
	_counters["NReco"]=0;
	_counters["N1"]=0;
	_counters["N2"]=0;
	_counters["N3"]=0;

	geneMap.clear() ;
	_hitCharge.clear() ;



	/////////////////for ECAL---------------------------------------------------
	// copy the flags from the input collection
	//GG : it should be checked why we put the flag like this.
	LCFlagImpl flag ;
	flag.setBit(LCIO::CHBIT_LONG) ;
	flag.setBit(LCIO::RCHBIT_TIME);

	for (unsigned int i(0) ; i < _inputCollections.size() ; ++i)
	{
		try
		{
			std::string inputColName = _inputCollections.at(i) ;
			std::string outputColName = _outputCollections.at(i) ;
			std::string outputRelColName = _outputRelCollections.at(i) ;

			LCCollection* inputCol = evt->getCollection( inputColName.c_str() ) ;
			_counters["NSim"] += inputCol->getNumberOfElements() ;
			CHT::Layout layout = layoutFromString( inputColName ) ;

			LCCollectionVec* outputCol = nullptr ;
			LCCollectionVec* outputRelCol = nullptr ;

			processCollection(inputCol , outputCol , outputRelCol , layout , flag) ;

			_counters["NReco"] += outputCol->getNumberOfElements() ;

			evt->addCollection(outputCol , outputColName.c_str()) ;
			evt->addCollection(outputRelCol , outputRelColName.c_str()) ;
		}
		catch(DataNotAvailableException& )
		{
		}
	}

	_tupleCollection->fill(0,_counters["NSim"]);
	_tupleCollection->fill(1,_counters["NReco"]);
	_tupleCollection->fill(2,_counters["N1"]);
	_tupleCollection->fill(3,_counters["N2"]);
	_tupleCollection->fill(4,_counters["N3"]);
	_tupleCollection->addRow();

	for( const auto& it : _hitCharge)
		_histoCellCharge->fill(it) ;

	streamlog_out(MESSAGE) << "have processed " << _counters["|ALL"] << " events" << std::endl;
}

