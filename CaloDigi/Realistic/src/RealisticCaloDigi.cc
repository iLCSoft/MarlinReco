// Calorimeter digitiser
#include "RealisticCaloDigi.h"

#include <marlin/Global.h>
#include <marlin/Exceptions.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <CalorimeterHitType.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

#include <iostream>
#include <string>
#include <assert.h>
#include <cmath>
#include <set>

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

#define RELATIONFROMTYPESTR "FromType"
#define RELATIONTOTYPESTR "ToType"

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace std::placeholders;

RealisticCaloDigi::RealisticCaloDigi( ) : Processor( "RealisticCaloDigi" ) {

  _description = "Performs digitization of sim calo hits. Virtual class." ;

  // input collections of simcalohits
  std::vector<std::string> inputCollections;
  inputCollections.push_back("SimCalorimeterHits");
  
  registerInputCollections( LCIO::SIMCALORIMETERHIT,
                            "inputHitCollections" ,
                            "input simcalhit Collection Names" ,
                            _inputCollections ,
                            inputCollections);

  // since we don't know number of output collections a priori, 
  // we can't use registerOutputCollection, but have to use registerProcessorParameter
  std::vector<std::string>  outputCollections;
  registerProcessorParameter( "outputHitCollections",
			      "output calorimeterhit Collection Names" ,
			      _outputCollections,
			      outputCollections );

  // the collections of relations between sim and digitised hits
  std::vector<std::string>  outputRelCollections;
  registerProcessorParameter( "outputRelationCollections",
			      "output hit relatiob Collection Names" ,
			      _outputRelCollections,
			      outputRelCollections );

  // energy threshold
  registerProcessorParameter("threshold" ,
                             "Threshold for Hit" ,
                             _threshold_value,
                             (float) 0.5 );

  registerProcessorParameter("thresholdUnit" ,
                             "Unit for threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants" ,
                             _threshold_unit,
                             std::string("MIP"));

  // timing
  registerProcessorParameter("timingCut" ,
                             "Use hit times" ,
                             _time_apply               ,
                             (int)0);

  registerProcessorParameter("timingCorrectForPropagation" ,
                             "Correct hit times for propagation: radial distance/c" ,
                             _time_correctForPropagation,
                             (int)0);

  registerProcessorParameter("timingWindowMin" ,
                             "Time Window minimum time in ns" ,
                             _time_windowMin,
                             (float)-10.0);

  registerProcessorParameter("timingWindowMax" ,
                             "Time Window maximum time in ns" ,
                             _time_windowMax,
                             (float)+100.0);

  registerProcessorParameter("integrationMethod" ,
                             "Energy integration and time calculation method. Options: Standard, ROC" ,
                             _integration_method,
                             std::string("Standard"));
                             
  registerProcessorParameter("fastShaper" ,
                             "Fast shaper value. Unit in ns" ,
                             _fast_shaper,
                             (float)0.f);
                             
  registerProcessorParameter("slowShaper" ,
                             "Slow shaper value. Unit in ns" ,
                             _slow_shaper,
                             (float)0.f);
                             
  registerProcessorParameter("timingResolution" ,
                             "Time resolution to apply (gaussian smearing). Unit in ns" ,
                             _time_resol,
                             (float)0);
                             

  // additional digi effects

  registerProcessorParameter("calibration_mip" ,
                             "average G4 deposited energy by MIP for calibration" ,
                             _calib_mip,
                             (float)1.0e-4);

  registerProcessorParameter("miscalibration_uncorrel" ,
                             "uncorrelated random gaussian miscalibration (as a fraction: 1.0 = 100%) " ,
                             _misCalib_uncorrel,
                             (float)0.0);

  registerProcessorParameter("miscalibration_uncorrel_memorise" ,
                             "store oncorrelated miscalbrations in memory? (i.e. use same miscalibrations in each event. WARNING: can take a lot of memory if used...) " ,
                             _misCalib_uncorrel_keep,
                             (bool)false);

  registerProcessorParameter("miscalibration_correl" ,
                             "correlated random gaussian miscalibration (as a fraction: 1.0 = 100%) " ,
                             _misCalib_correl,
                             (float)0.0);

  registerProcessorParameter("deadCell_fraction" ,
                             "random dead cell fraction (as a fraction: 0->1) " ,
                             _deadCell_fraction,
                             (float)0.);

  registerProcessorParameter("deadCell_memorise" ,
                             "store dead cells in memory? (i.e. use same dead cells in each event. WARNING: can take a lot of memory if used...) " ,
                             _deadCell_keep,
                             (bool)false);

  // simple model of electronics properties
  registerProcessorParameter("elec_noise_mip",
                             "typical electronics noise (in MIP units)",
                             _elec_noiseMip,
                             float (0));

  registerProcessorParameter("elec_range_mip",
                             "maximum of dynamic range of electronics (in MIPs)",
                             _elec_rangeMip,
                             float (2500) );

  // code for layer info for cellID decoder
  registerProcessorParameter("CellIDLayerString" ,
                             "name of the part of the cellID that holds the layer" , 
                             _cellIDLayerString , 
                             std::string("K-1")
                             );

}

void RealisticCaloDigi::init() {

  printParameters();

  // if no output collection names specified, set some default based on the input collection names
  if ( _outputCollections.size()==0 ) {
    for (size_t i=0; i<_inputCollections.size(); i++) {
      _outputCollections.push_back( _inputCollections[i] + "Digi" );
    }
  }
  if ( _outputRelCollections.size()==0 ) {
    for (size_t i=0; i<_inputCollections.size(); i++) {
      _outputRelCollections.push_back( _inputCollections[i] + "DigiRelation" );
    }
  }

  // check that number of input and output collections names are the same
  assert ( _outputCollections.size()    == _inputCollections.size() );
  assert ( _outputRelCollections.size() == _inputCollections.size() );

  // unit in which threshold is specified
  if (_threshold_unit.compare("MIP") == 0){
    _threshold_iunit=MIP;
  } else if (_threshold_unit.compare("GeV") == 0){
    _threshold_iunit=GEVDEP;
  } else if (_threshold_unit.compare("px") == 0){
    _threshold_iunit=NPE;
  } else {
    streamlog_out(ERROR) << "could not identify threshold unit. Please use \"GeV\", \"MIP\" or \"px\"! Aborting." << std::endl;
    marlin::StopProcessingException(this);
  }
  
  // convert the threshold to the approriate units (i.e. MIP for silicon, NPE for scint)
  _threshold_value = convertEnergy( _threshold_value, _threshold_iunit );

  // deal with timing calculations  
  std::map<std::string, integr_function> integrations = {
    {"Standard", std::bind(&RealisticCaloDigi::StandardIntegration, this, _1)},
    {"ROC", std::bind(&RealisticCaloDigi::ROCIntegration, this, _1)}
  };  
  auto findIter = integrations.find( _integration_method ) ;
  if(integrations.end() == findIter) {
    streamlog_out(ERROR) << "Could not guess timing calculation method!" << std::endl;
    streamlog_out(ERROR) << "Available are: Standard, ROC. Provided: " << _integration_method << std::endl;
    streamlog_out(ERROR) << "Aborting..." << std::endl;
    marlin::StopProcessingException(this);
  }
  _integr_function = findIter->second;
  
  // check if parameters are correctly set for the ROC integration
  if("ROC" == _integration_method) {
    if(!parameterSet("fastShaper") || !parameterSet("slowShaper")) {
      streamlog_out(ERROR) << "Fast/slow shaper parameter(s) not set. Required for ROC integration!" << std::endl;
      streamlog_out(ERROR) << "Aborting..." << std::endl;
      marlin::StopProcessingException(this);
    }
  }

  _flag.setBit(LCIO::CHBIT_LONG);
  _flag.setBit(LCIO::RCHBIT_TIME); //store timing on output hits.

  _flag_rel.setBit(LCIO::LCREL_WEIGHTED); // for the hit relations


  return;
}


void RealisticCaloDigi::processRunHeader( LCRunHeader*  /*run*/) {
}

void RealisticCaloDigi::processEvent( LCEvent * evt ) {

  // create the output collections

  // decide on this event's correlated miscalibration
  if ( _misCalib_correl>0 ) _event_correl_miscalib = CLHEP::RandGauss::shoot( 1.0, _misCalib_correl );

  //
  // * Reading Collections of Simulated Hits *
  //

  for (unsigned int i(0); i < _inputCollections.size(); ++i) {
    std::string colName =  _inputCollections[i] ;
    streamlog_out ( DEBUG1 ) << "looking for collection: " << colName << endl;
    try{
      LCCollection * col = evt->getCollection( colName.c_str() ) ;
      string initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      CHT::CaloType cht_type = caloTypeFromString(colName);
      CHT::CaloID   cht_id   = caloIDFromString(colName);
      CHT::Layout   cht_lay  = layoutFromString(colName);

      CellIDDecoder<SimCalorimeterHit> idDecoder( col );

      int numElements = col->getNumberOfElements();
      streamlog_out ( DEBUG1 ) << colName << " number of elements = " << numElements << endl;

      if ( numElements==0 ) continue;

      // create new collection: hits
      LCCollectionVec *newcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
      newcol->setFlag(_flag.getFlag());

      // hit relations to simhits [calo -> sim]
      LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);
      relcol->setFlag(_flag_rel.getFlag());
      relcol->parameters().setValue( RELATIONFROMTYPESTR , LCIO::CALORIMETERHIT ) ;
      relcol->parameters().setValue( RELATIONTOTYPESTR   , LCIO::SIMCALORIMETERHIT ) ;

      // loop over input hits
      for (int j=0; j < numElements; ++j) {
        SimCalorimeterHit * simhit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;

        // deal with energy integration and timing aspects
        auto integrationResult = Integrate(simhit);
        if( ! integrationResult.has_value() ) {
          continue;
        }
        float time      = integrationResult.value().first;
        float energyDep = integrationResult.value().second;
        // apply extra energy digitisation onto the energy
        float energyDig = EnergyDigi(energyDep, simhit->getCellID0() , simhit->getCellID1() );

	      streamlog_out ( DEBUG0 ) << " hit " << j << " time: " << time << " eDep: " << energyDep << " eDigi: " << energyDig << " " << _threshold_value << endl;

        if (energyDig > _threshold_value) { // write out this hit
          CalorimeterHitImpl* newhit = new CalorimeterHitImpl();
          newhit->setCellID0( simhit->getCellID0() );
          newhit->setCellID1( simhit->getCellID1() );
          newhit->setTime( time );
          newhit->setPosition( simhit->getPosition() );
    	    newhit->setEnergy( energyDig );
    	    int layer = idDecoder(simhit)[_cellIDLayerString];
    	    newhit->setType( CHT( cht_type, cht_id, cht_lay, layer ) );
          newhit->setRawHit( simhit );
          newcol->addElement( newhit ); // add hit to output collection

	        streamlog_out ( DEBUG1 ) << "orig/new hit energy: " << simhit->getEnergy() << " " << newhit->getEnergy() << endl;

          LCRelationImpl *rel = new LCRelationImpl(newhit,simhit,1.0);
          relcol->addElement( rel );

        } // theshold
      } // input hits

      // add collection to event
      newcol->parameters().setValue(LCIO::CellIDEncoding,initString);
      evt->addCollection(newcol,_outputCollections[i].c_str());

      // add relation collection to event
      evt->addCollection(relcol,_outputRelCollections[i].c_str());

    } catch(DataNotAvailableException &e){
      streamlog_out(DEBUG1) << "could not find input collection " << colName << std::endl;
    }
  }


  return;
}

//------------------------------------------------------------------------------

void RealisticCaloDigi::check( LCEvent *  /*evt*/ ) {
  
}

//------------------------------------------------------------------------------

void RealisticCaloDigi::end(){
}

//------------------------------------------------------------------------------

RealisticCaloDigi::integr_res_opt RealisticCaloDigi::Integrate( const SimCalorimeterHit * hit ) const {
  return _integr_function(hit);
}

//------------------------------------------------------------------------------

float RealisticCaloDigi::EnergyDigi(float energy, int id0, int id1) {
  // some extra digi effects
  // controlled by _applyDigi = 0 (none), 1 (apply)
  // input parameters: hit energy ( in any unit: effects are all relative )
  //                   id0,1 - cell IDs (used to memorise miscalibrations/dead cells between events if requested)
  // returns energy ( in units determined by the overloaded digitiseDetectorEnergy )

  float e_out(energy);
  e_out = digitiseDetectorEnergy(energy); // this is an overloaded method, provides energy in technology-dependent units

  // the following make only relative changes to the energy

  // random miscalib, uncorrelated in cells
  if (_misCalib_uncorrel>0) {
    float miscal(0);
    if ( _misCalib_uncorrel_keep ) { // memorise the miscalibrations from event-to-event
      std::pair <int, int> id(id0, id1);
      if ( _cell_miscalibs.find(id)!=_cell_miscalibs.end() ) { // this cell was previously seen, and a miscalib stored
        miscal = _cell_miscalibs[id];
      } else {                                                 // we haven't seen this one yet, get a new miscalib for it
        miscal = CLHEP::RandGauss::shoot( 1.0, _misCalib_uncorrel );
        _cell_miscalibs[id]=miscal;
      }
    } else { // get a new calibration for each hit
      miscal = CLHEP::RandGauss::shoot( 1.0, _misCalib_uncorrel );
    }
    e_out*=miscal;
  }

  // random miscalib, correlated across cells in one event
  if (_misCalib_correl>0) e_out*=_event_correl_miscalib;

  float oneMipInMyUnits = convertEnergy( 1.0, MIP );
  // limited electronics dynamic range
  if ( _elec_rangeMip > 0 ) e_out = std::min ( e_out, _elec_rangeMip*oneMipInMyUnits );
  // add electronics noise
  if ( _elec_noiseMip > 0 ) e_out += CLHEP::RandGauss::shoot(0, _elec_noiseMip*oneMipInMyUnits );

  // random cell kill
  if (_deadCell_fraction>0){
    if (_deadCell_keep == true){ // memorise dead cells from event to event
      std::pair <int, int> id(id0, id1);
      if (_cell_dead.find(id)!=_cell_dead.end() ) { // this cell was previously seen
        if (_cell_dead[id] == true){
          e_out = 0;
        }
      } else { // we haven't seen this one yet, get a miscalib for it
        bool thisDead = (CLHEP::RandFlat::shoot(0.0, 1.0) < _deadCell_fraction);
        _cell_dead[id] = thisDead;
        if (thisDead == true){
          e_out = 0;
        }
      }
    } else { // decide on cell-by-cell basis
      if (CLHEP::RandFlat::shoot(0.0,1.0)<_deadCell_fraction ) e_out=0;
    }
  }
  return e_out;
}

//------------------------------------------------------------------------------

RealisticCaloDigi::integr_res_opt RealisticCaloDigi::StandardIntegration( const SimCalorimeterHit * hit ) const {
  // apply timing cuts on simhit contributions
  // outputs a (time,energy) pair
  float timeCorrection(0);
  if ( _time_correctForPropagation ) { // time of flight from IP to this point
    float r(0);
    for (int i=0; i<3; i++) 
      r+=pow(hit->getPosition()[i],2); 
    timeCorrection = sqrt(r)/CLHEP::c_light; // [speed of light in mm/ns]
  }
  // this is Oskar's simple (and probably the most correct) method for treatment of timing
  //  - collect energy in some predefined time window around collision time (possibly corrected for TOF)
  //  - assign time of earliest contribution to hit
  float energySum = 0;
  float earliestTime=std::numeric_limits<float>::max();
  for(int i = 0; i<hit->getNMCContributions(); i++){ // loop over all contributions
    float timei   = hit->getTimeCont(i); //absolute hit timing of current subhit
    float energyi = hit->getEnergyCont(i); //energy of current subhit
    float relativetime = timei - timeCorrection; // wrt time of flight
    if (relativetime>_time_windowMin && relativetime<_time_windowMax){
      energySum += energyi;
      if (relativetime<earliestTime){
	       earliestTime = relativetime; //use earliest hit time for simpletimingcut
      }
    }
  }      
  if(earliestTime > _time_windowMin && earliestTime < _time_windowMax){ //accept this hit
    return integr_res{SmearTime(earliestTime), energySum};
  }
  return integr_res{0, hit->getEnergy()};
}

//------------------------------------------------------------------------------

RealisticCaloDigi::integr_res_opt RealisticCaloDigi::ROCIntegration( const SimCalorimeterHit * hit ) const {
  const unsigned int ncontrib = hit->getNMCContributions() ;
  // Sort MC contribution by time
  std::vector<EVENT::MCParticle*> mcconts{ncontrib, nullptr};
  for(unsigned int i = 0; i<ncontrib; i++){
    mcconts[i]=hit->getParticleCont(i);
  }
  std::sort(mcconts.begin(), mcconts.end(), [](auto lhs, auto rhs){
    return (lhs->getTime() < rhs->getTime()); 
  });
  // Accumulate energy until threshold is reached.
  // The first MC contriubtion after the threshold has been reached sets the hit time 
  bool passThreshold = false;
  float epar=0.f, hitTime=0.f;
  unsigned int thresholdIndex=0;
  // First determine the hit time (hitTime) and the initial hit index 
  // at which we need to start the integration (thresholdIndex)
  for(unsigned int i=0; i<ncontrib ; ++i) {
    const auto timei = mcconts[i]->getTime();
    thresholdIndex = i ;
    epar = 0.f;
    for(unsigned int j=i; j<ncontrib ; ++j) {
      const auto timej = mcconts[j]->getTime();
      if( (timej-timei) < _fast_shaper) {
        epar += mcconts[j]->getEnergy();
      }
      else {
        break;
      }
      if( epar > _threshold_value ) {
        hitTime = timej;
        passThreshold = true ;
        break;
      }
    }
    if(passThreshold) {
      break;
    }
  }
  // If we've found a hit above the threshold, accumulate the energy until
  // until the maximum time given by the slow shaper
  if(passThreshold) {
    const float thresholdTime = mcconts[thresholdIndex]->getTime();
    float energySum = 0.f ; 
    for(unsigned int i=thresholdIndex ; i<ncontrib ; ++i) {
      if( mcconts[i]->getTime() < thresholdTime + _slow_shaper) {
        energySum += mcconts[i]->getEnergy();
      }
    }
    hitTime = SmearTime(hitTime);
    return integr_res{thresholdTime, energySum};
  }
  // else no hit dude !
  else {
    return std::nullopt;
  }
}

//------------------------------------------------------------------------------

float RealisticCaloDigi::SmearTime(float time) const{
  return _time_resol>0.f ? time + CLHEP::RandGauss::shoot(0, _time_resol ) : time;
}

