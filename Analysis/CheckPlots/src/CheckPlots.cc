#include "CheckPlots.h"


using namespace lcio ;
using namespace marlin ;




CheckPlots aCheckPlots ;



CheckPlots::CheckPlots() : Processor("CheckPlots") {
  
  _description = "produces Check Plots" ;



    
  registerProcessorParameter("FillMC","fill clouds for MC particles",
			     _fillMC,
			     (int)1);
  
  registerProcessorParameter("FillMCSim","fill clouds for MC particles created during simulation",
			     _fillMCSim,
			     (int)0);

  
  registerProcessorParameter("FillSimCalo","fill clouds for SimCalorimeter hits",
			     _fillSimCalo,
			     (int)1);
  registerProcessorParameter("SimECut","energy cut for filling clouds for SimCalorimeter hits",
			     _simECut,
			     (float)0.0001);
  
  
  registerProcessorParameter("FillCalo","fill clouds for Calorimeter hits",
			     _fillCalo,
			     (int)1);
  
  registerProcessorParameter("ECut","energy cut for filling clouds for Calorimeter hits",
			     _ECut,
			     (float)0.0001);

  
}




void CheckPlots::init() { 

  // printParameters();

  _nRun = -1;
  _nEvt = 0;

  
}



void CheckPlots::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
  _nEvt = 0;

} 




void CheckPlots::processEvent( LCEvent * evt ) { 

  
  static bool firstEvent = true;
  
  if(firstEvent) std::cout << "Check Plot processor called for first event" << std::endl;


  if (_fillMC) fillMCCheckPlots(evt);

  //  fillSimTrackerCheckPlots(evt);
  if (_fillSimCalo) fillSimCaloCheckPlots(evt);
  
  //  fillTrackerCheckPlots(evt);
  if (_fillCalo) fillCaloCheckPlots(evt);
  

  _nEvt ++;
  firstEvent = false;
 
}




void CheckPlots::check( LCEvent * evt ) { 
  
}




void CheckPlots::end(){ 


}










void CheckPlots::fillMCCheckPlots(LCEvent * evt){

  
  #ifdef MARLIN_USE_AIDA
  
  
  try {

    // MCPs with generator status != 1
    // accumulated numbers
    //  static AIDA::ICloud1D* cMCPIDAccSim;
    
    // MCPs with generator status == 1
    // accumulated numbers
    //  static AIDA::ICloud1D* cMCPIDAccGen;
    
    
	
    // MCPs with generator status != 1
    // numbers per event
    static AIDA::ICloud1D* cMCNumberSim;
    static AIDA::ICloud1D* cMCEnergySumSim;
    
    static AIDA::ICloud1D* cMCNumberElectronsSim;
    static AIDA::ICloud1D* cMCNumberPositronsSim;
    static AIDA::ICloud1D* cMCNumberMuonsSim;
    static AIDA::ICloud1D* cMCNumberMuonBarsSim;
    static AIDA::ICloud1D* cMCNumberTausSim;
    static AIDA::ICloud1D* cMCNumberTauBarsSim;
    static AIDA::ICloud1D* cMCNumberNusSim;
    
    static AIDA::ICloud1D* cMCNumberPiPlusSim;
    static AIDA::ICloud1D* cMCNumberPiMinusSim;
    static AIDA::ICloud1D* cMCNumberKPlusSim;
    static AIDA::ICloud1D* cMCNumberKMinusSim;
    static AIDA::ICloud1D* cMCNumberProtonsSim;
    static AIDA::ICloud1D* cMCNumberProtonBarsSim;
    static AIDA::ICloud1D* cMCNumberPi0Sim;
    static AIDA::ICloud1D* cMCNumberK0lSim;
    static AIDA::ICloud1D* cMCNumberK0sSim;
    static AIDA::ICloud1D* cMCNumberNeutronsSim;
    static AIDA::ICloud1D* cMCNumberGammasSim;
    static AIDA::ICloud1D* cMCNumberLambdasSim;
    
    static AIDA::ICloud1D* cMCNumberRemainingSim;
    
    
    
    // MCPs with generator status == 1
    // numbers per event
    static AIDA::ICloud1D* cMCNumberGen;
    static AIDA::ICloud1D* cMCEnergySumGen;
    
    static AIDA::ICloud1D* cMCNumberElectronsGen;
    static AIDA::ICloud1D* cMCNumberPositronsGen;
    static AIDA::ICloud1D* cMCNumberMuonsGen;
    static AIDA::ICloud1D* cMCNumberMuonBarsGen;
    static AIDA::ICloud1D* cMCNumberTausGen;
    static AIDA::ICloud1D* cMCNumberTauBarsGen;
    static AIDA::ICloud1D* cMCNumberNusGen;
    
    static AIDA::ICloud1D* cMCNumberPiPlusGen;
    static AIDA::ICloud1D* cMCNumberPiMinusGen;
    static AIDA::ICloud1D* cMCNumberKPlusGen;
    static AIDA::ICloud1D* cMCNumberKMinusGen;
    static AIDA::ICloud1D* cMCNumberProtonsGen;
    static AIDA::ICloud1D* cMCNumberProtonBarsGen;
    static AIDA::ICloud1D* cMCNumberPi0Gen;
    static AIDA::ICloud1D* cMCNumberK0lGen;
    static AIDA::ICloud1D* cMCNumberK0sGen;
    static AIDA::ICloud1D* cMCNumberNeutronsGen;
    static AIDA::ICloud1D* cMCNumberGammasGen;
    static AIDA::ICloud1D* cMCNumberLambdasGen;
    
    static AIDA::ICloud1D* cMCNumberRemainingGen;
    
    
    
    // MCPs with generator status != 1
    // numbers per single particle
    static AIDA::ICloud1D* cMCEnergySim;
    static AIDA::ICloud1D* cMCEnergyElectronsSim;
    static AIDA::ICloud1D* cMCEnergyPositronsSim;
    static AIDA::ICloud1D* cMCEnergyMuonsSim;
    static AIDA::ICloud1D* cMCEnergyMuonBarsSim;
    static AIDA::ICloud1D* cMCEnergyTausSim;
    static AIDA::ICloud1D* cMCEnergyTauBarsSim;
    static AIDA::ICloud1D* cMCEnergyNusSim;
    
    static AIDA::ICloud1D* cMCEnergyPiPlusSim;
    static AIDA::ICloud1D* cMCEnergyPiMinusSim;
    static AIDA::ICloud1D* cMCEnergyKPlusSim;
    static AIDA::ICloud1D* cMCEnergyKMinusSim;
    static AIDA::ICloud1D* cMCEnergyProtonsSim;
    static AIDA::ICloud1D* cMCEnergyProtonBarsSim;
    static AIDA::ICloud1D* cMCEnergyPi0Sim;
    static AIDA::ICloud1D* cMCEnergyK0lSim;
    static AIDA::ICloud1D* cMCEnergyK0sSim;
    static AIDA::ICloud1D* cMCEnergyNeutronsSim;
    static AIDA::ICloud1D* cMCEnergyGammasSim;
    static AIDA::ICloud1D* cMCEnergyLambdasSim;
    
    static AIDA::ICloud1D* cMCEnergyRemainingSim;



    // MCPs with generator status == 1
    // numbers per single particle
    static AIDA::ICloud1D* cMCEnergyGen;
    static AIDA::ICloud1D* cMCEnergyElectronsGen;
    static AIDA::ICloud1D* cMCEnergyPositronsGen;
    static AIDA::ICloud1D* cMCEnergyMuonsGen;
    static AIDA::ICloud1D* cMCEnergyMuonBarsGen;
    static AIDA::ICloud1D* cMCEnergyTausGen;
    static AIDA::ICloud1D* cMCEnergyTauBarsGen;
    static AIDA::ICloud1D* cMCEnergyNusGen;
    
    static AIDA::ICloud1D* cMCEnergyPiPlusGen;
    static AIDA::ICloud1D* cMCEnergyPiMinusGen;
    static AIDA::ICloud1D* cMCEnergyKPlusGen;
    static AIDA::ICloud1D* cMCEnergyKMinusGen;
    static AIDA::ICloud1D* cMCEnergyProtonsGen;
    static AIDA::ICloud1D* cMCEnergyProtonBarsGen;
    static AIDA::ICloud1D* cMCEnergyPi0Gen;
    static AIDA::ICloud1D* cMCEnergyK0lGen;
    static AIDA::ICloud1D* cMCEnergyK0sGen;
    static AIDA::ICloud1D* cMCEnergyNeutronsGen;
    static AIDA::ICloud1D* cMCEnergyGammasGen;
    static AIDA::ICloud1D* cMCEnergyLambdasGen;
    
    static AIDA::ICloud1D* cMCEnergyRemainingGen;




	
    if( isFirstEvent() ) {
      
      if (_fillMCSim) {
	    
	// accumulated numbers
	// cMCPIDAccSim             = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCPIDAccSim", "MC particle IDs accumulated over all events (generator status != 1)", -1 );
	
      }
	  
      
      // cMCPIDAccGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCPIDAccGen", "MC particle IDs accumulated over all events (generator status == 1)", -1 );
      
      
      
      if (_fillMCSim) {

	// numbers per event
	cMCNumberSim           = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberSim", "number of all MC particles per event (generator status != 1)", -1 );
	cMCEnergySumSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergySumSim", "MC particle energy sum per event (generator status != 1)", -1 );
	
	cMCNumberElectronsSim  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberElectronsSim", "number of the electrons per event (generator status != 1)", -1 );
	cMCNumberPositronsSim  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPositronsSim", "number of the positrons per event (generator status != 1)", -1 );
	cMCNumberMuonsSim      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberMuonsSim", "number of the muons per event (generator status != 1)", -1 );
	cMCNumberMuonBarsSim   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberMuonBarsSim", "number of the anti-muons per event (generator status != 1)", -1 );
	cMCNumberTausSim       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberTausSim", "number of the taus per event (generator status != 1)", -1 );
	cMCNumberTauBarsSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberTauBarsSim", "number of the anti-taus per event (generator status != 1)", -1 );
	cMCNumberNusSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberNusSim", "number of neutrinos per event (generator status != 1)", -1 );
	
	cMCNumberPiPlusSim     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPiPlusSim", "number of Pi+ per event (generator status != 1)", -1 );
	cMCNumberPiMinusSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPiMinusSim", "number of Pi- per event (generator status != 1)", -1 );
	cMCNumberKPlusSim      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberKPlusSim", "number of K+ per event (generator status != 1)", -1 );
	cMCNumberKMinusSim     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberKMinusSim", "number of K- per event (generator status != 1)", -1 );
	cMCNumberProtonsSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberProtonsSim", "number of protons per event (generator status != 1)", -1 );
	cMCNumberProtonBarsSim = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberProtonBarsSim", "number of anti-protons per event (generator status != 1)", -1 );
	cMCNumberPi0Sim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPi0Sim", "number of Pi0 per event (generator status != 1)", -1 );
	cMCNumberK0lSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberK0lSim", "number of K0l per event (generator status != 1)", -1 );
	cMCNumberK0sSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberK0sSim", "number of K0s per event (generator status != 1)", -1 );
	cMCNumberNeutronsSim   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberNeutronsSim", "number of neutrons per event (generator status != 1)", -1 );
	cMCNumberGammasSim     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberGammasSim", "number of anti-neutrons per event (generator status != 1)", -1 );
	cMCNumberLambdasSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberLambdasSim", "number of lambdas per event (generator status != 1)", -1 );
	
	cMCNumberRemainingSim  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberRemainingSim", "number of the remaining particles per event (generator status != 1)", -1 );
	
      }
      
      
      cMCNumberGen             = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberGen", "number of all MC particles per event (generator status == 1)", -1 );
      cMCEnergySumGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergySumGen", "MC particle energy sum per event (generator status == 1)", -1 );
      
      cMCNumberElectronsGen    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberElectronsGen", "number of the electrons per event (generator status == 1)", -1 );
      cMCNumberPositronsGen    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPositronsGen", "number of the positrons per event (generator status == 1)", -1 );
      cMCNumberMuonsGen        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberMuonsGen", "number of the muons per event (generator status == 1)", -1 );
      cMCNumberMuonBarsGen     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberMuonBarsGen", "number of the anti-muons per event (generator status == 1)", -1 );
      cMCNumberTausGen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberTausGen", "number of the taus per event (generator status == 1)", -1 );
      cMCNumberTauBarsGen      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberTauBarsGen", "number of the anti-taus per event (generator status == 1)", -1 );
      cMCNumberNusGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberNusGen", "number of neutrinos per event (generator status == 1)", -1 );
      
      cMCNumberPiPlusGen       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPiPlusGen", "number of Pi+ per event (generator status == 1)", -1 );
      cMCNumberPiMinusGen      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPiMinusGen", "number of Pi- per event (generator status == 1)", -1 );
      cMCNumberKPlusGen        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberKPlusGen", "number of K+ per event (generator status == 1)", -1 );
      cMCNumberKMinusGen       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberKMinusGen", "number of K- per event (generator status == 1)", -1 );
      cMCNumberProtonsGen      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberProtonsGen", "number of protons per event (generator status == 1)", -1 );
      cMCNumberProtonBarsGen   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberProtonBarsGen", "number of anti-protons per event (generator status == 1)", -1 );
      cMCNumberPi0Gen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPi0Gen", "number of Pi0 per event (generator status == 1)", -1 );
      cMCNumberK0lGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberK0lGen", "number of K0l per event (generator status == 1)", -1 );
      cMCNumberK0sGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberK0sGen", "number of K0s per event (generator status == 1)", -1 );
      cMCNumberNeutronsGen     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberNeutronsGen", "number of neutrons per event (generator status == 1)", -1 );
      cMCNumberGammasGen       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberGammasGen", "number of anti-neutrons per event (generator status == 1)", -1 );
      cMCNumberLambdasGen      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberLambdasGen", "number of lambdas per event (generator status == 1)", -1 );
      
      cMCNumberRemainingGen    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberRemainingGen", "number of the remaining particles per event (generator status == 1)", -1 );
      
      
      
      
      if (_fillMCSim) {
	
	// numbers per single particle
	cMCEnergySim           = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergySim", "energy spectrum of all MC particles in one event (generator status != 1)", -1 );
	cMCEnergyElectronsSim  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyElectronsSim", "energy spectrum of the electrons in one event (generator status != 1)", -1 );
	cMCEnergyPositronsSim  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPositronsSim", "energy spectrum of the positrons in one event (generator status != 1)", -1 );
	cMCEnergyMuonsSim      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyMuonsSim", "energy spectrum of the muons in one event (generator status != 1)", -1 );
	cMCEnergyMuonBarsSim   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyMuonBarsSim", "energy spectrum of the anti-muons in one event (generator status != 1)", -1 );
	cMCEnergyTausSim       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyTausSim", "energy spectrum of the taus in one event (generator status != 1)", -1 );
	cMCEnergyTauBarsSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyTauBarsSim", "energy spectrum of the anti-taus in one event (generator status != 1)", -1 );
	cMCEnergyNusSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyNusSim", "energy spectrum of the neutrinos in one event (generator status != 1)", -1 );
	    
	cMCEnergyPiPlusSim     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPiPlusSim", "energy spectrum of Pi+ in one event (generator status != 1)", -1 );
	cMCEnergyPiMinusSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPiMinusSim", "energy spectrum of Pi- in one event (generator status != 1)", -1 );
	cMCEnergyKPlusSim      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyKPlusSim", "energy spectrum of K+ in one event (generator status != 1)", -1 );
	cMCEnergyKMinusSim     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyKMinusSim", "energy spectrum of K- in one event (generator status != 1)", -1 );
	cMCEnergyProtonsSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyProtonsSim", "energy spectrum of protons in one event (generator status != 1)", -1 );
	cMCEnergyProtonBarsSim = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyProtonBarsSim", "energy spectrum of anti-protons in one event (generator status != 1)", -1 );
	cMCEnergyPi0Sim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPi0Sim", "energy spectrum of Pi0 in one event (generator status != 1)", -1 );
	cMCEnergyK0lSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyK0lSim", "energy spectrum of K0l in one event (generator status != 1)", -1 );
	cMCEnergyK0sSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyK0sSim", "energy spectrum of K0s in one event (generator status != 1)", -1 );
	cMCEnergyNeutronsSim   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyNeutronsSim", "energy spectrum of neutrons in one event (generator status != 1)", -1 );
	cMCEnergyGammasSim     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyGammasSim", "energy spectrum of anti-neutrons in one event (generator status != 1)", -1 );
	cMCEnergyLambdasSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyLambdasSim", "energy spectrum of lambdas in one event (generator status != 1)", -1 );
	
	cMCEnergyRemainingSim  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyRemainingSim", "energy spectrum of the remaining particles in one event (generator status != 1)", -1 );
	
      }
      
      
      cMCEnergyGen             = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyGen", "energy spectrum of all MC particles in one event (generator status == 1)", -1 );
      cMCEnergyElectronsGen    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyElectronsGen", "energy spectrum of the electrons in one event (generator status == 1)", -1 );
      cMCEnergyPositronsGen    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPositronsGen", "energy spectrum of the positrons in one event (generator status == 1)", -1 );
      cMCEnergyMuonsGen        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyMuonsGen", "energy spectrum of the muons in one event (generator status == 1)", -1 );
      cMCEnergyMuonBarsGen     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyMuonBarsGen", "energy spectrum of the anti-muons in one event (generator status == 1)", -1 );
      cMCEnergyTausGen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyTausGen", "energy spectrum of the taus in one event (generator status == 1)", -1 );
      cMCEnergyTauBarsGen      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyTauBarsGen", "energy spectrum of the anti-taus in one event (generator status == 1)", -1 );
      cMCEnergyNusGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyNusGen", "energy spectrum of the neutrinos in one event (generator status == 1)", -1 );
      
      cMCEnergyPiPlusGen       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPiPlusGen", "energy spectrum of Pi+ in one event (generator status == 1)", -1 );
      cMCEnergyPiMinusGen      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPiMinusGen", "energy spectrum of Pi- in one event (generator status == 1)", -1 );
      cMCEnergyKPlusGen        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyKPlusGen", "energy spectrum of K+ in one event (generator status == 1)", -1 );
      cMCEnergyKMinusGen       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyKMinusGen", "energy spectrum of K- in one event (generator status == 1)", -1 );
      cMCEnergyProtonsGen      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyProtonsGen", "energy spectrum of protons in one event (generator status == 1)", -1 );
      cMCEnergyProtonBarsGen   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyProtonBarsGen", "energy spectrum of anti-protons in one event (generator status == 1)", -1 );
      cMCEnergyPi0Gen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPi0Gen", "energy spectrum of Pi0 in one event (generator status == 1)", -1 );
      cMCEnergyK0lGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyK0lGen", "energy spectrum of K0l in one event (generator status == 1)", -1 );
      cMCEnergyK0sGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyK0sGen", "energy spectrum of K0s in one event (generator status == 1)", -1 );
      cMCEnergyNeutronsGen     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyNeutronsGen", "energy spectrum of neutrons in one event (generator status == 1)", -1 );
      cMCEnergyGammasGen       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyGammasGen", "energy spectrum of anti-neutrons in one event (generator status == 1)", -1 );
      cMCEnergyLambdasGen      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyLambdasGen", "energy spectrum of lambdas in one event (generator status == 1)", -1 );
      
      cMCEnergyRemainingGen    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyRemainingGen", "energy spectrum of the remaining particles in one event (generator status == 1)", -1 );
      
      
    }
    
    
    
    const std::vector< std::string >* strVec = evt->getCollectionNames() ;
    std::vector< std::string >::const_iterator name ;
    
    for( name = strVec->begin() ; name != strVec->end() ; name++) {
      
      LCCollection* col = evt->getCollection( *name ) ;

      if ( col->getTypeName() == LCIO::MCPARTICLE ) {
	
	if( col != 0 ){
	  
	  unsigned int NMCSim             = 0;
	  unsigned int NMCElectronsSim    = 0;
	  unsigned int NMCPositronsSim    = 0;
	  unsigned int NMCMuonsSim        = 0;
	  unsigned int NMCMuonBarsSim     = 0;
	  unsigned int NMCTausSim         = 0;
	  unsigned int NMCTauBarsSim      = 0;
	  unsigned int NMCNusSim          = 0;
	  
	  unsigned int NMCPiPlusSim       = 0;
	  unsigned int NMCPiMinusSim      = 0;
	  unsigned int NMCKPlusSim        = 0;
	  unsigned int NMCKMinusSim       = 0;
	  unsigned int NMCProtonsSim      = 0;
	  unsigned int NMCProtonBarsSim   = 0;
	  unsigned int NMCPi0Sim          = 0;
	  unsigned int NMCK0lSim          = 0;
	  unsigned int NMCK0sSim          = 0;
	  unsigned int NMCNeutronsSim     = 0;
	  unsigned int NMCGammasSim       = 0;
	  unsigned int NMCLambdasSim      = 0;
	  
	  unsigned int NMCRemainingSim    = 0;
	  
	  double energySim                = 0.0;

	  
	  unsigned int NMCGen             = 0;
	  unsigned int NMCElectronsGen    = 0;
	  unsigned int NMCPositronsGen    = 0;
	  unsigned int NMCMuonsGen        = 0;
	  unsigned int NMCMuonBarsGen     = 0;
	  unsigned int NMCTausGen         = 0;
	  unsigned int NMCTauBarsGen      = 0;
	  unsigned int NMCNusGen          = 0;
	  
	  unsigned int NMCPiPlusGen       = 0;
	  unsigned int NMCPiMinusGen      = 0;
	  unsigned int NMCKPlusGen        = 0;
	  unsigned int NMCKMinusGen       = 0;
	  unsigned int NMCProtonsGen      = 0;
	  unsigned int NMCProtonBarsGen   = 0;
	  unsigned int NMCPi0Gen          = 0;
	  unsigned int NMCK0lGen          = 0;
	  unsigned int NMCK0sGen          = 0;
	  unsigned int NMCNeutronsGen     = 0;
	  unsigned int NMCGammasGen       = 0;
	  unsigned int NMCLambdasGen      = 0;
	  
	  unsigned int NMCRemainingGen    = 0;
	  
	  double energyGen                = 0.0;

	  
	  
	  int nMCP = col->getNumberOfElements();
	  
	  for(int i = 0; i < nMCP ; ++i){
	    
	    MCParticle* mcp = dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;
	    
	    if (mcp->getGeneratorStatus() != 1 ) {
	      
	      if (_fillMCSim) {
		
		++NMCSim;
		energySim += mcp->getEnergy();
		
		cMCEnergySim->fill(mcp->getEnergy());
		
		switch (mcp->getPDG()) {
		  
		case 11  : {
		  ++NMCElectronsSim;
		  cMCEnergyElectronsSim->fill(mcp->getEnergy());
		}
		case -11 : {
		  ++NMCPositronsSim;
		  cMCEnergyPositronsSim->fill(mcp->getEnergy());
		}
		case 13  : {
		  ++NMCMuonsSim;
		  cMCEnergyMuonsSim->fill(mcp->getEnergy());
		}
		case -13 : {
		  ++NMCMuonBarsSim; 
		  cMCEnergyMuonBarsSim->fill(mcp->getEnergy());
		}
		case 15  : {
		  ++NMCTausSim;
		  cMCEnergyTausSim->fill(mcp->getEnergy());
		}
		case -15 : {
		  ++NMCTauBarsSim;
		  cMCEnergyTauBarsSim->fill(mcp->getEnergy());
		}
		case  12 :
		case -12 :
		case  14 :
		case -14 :
		case  16 :
		case -16 : {
		  ++NMCNusSim;
		  cMCEnergyNusSim->fill(mcp->getEnergy());
		}
		  
		case 211  : {
		  ++NMCPiPlusSim;
		  cMCEnergyPiPlusSim->fill(mcp->getEnergy());
		}
		case -211 : {
		  ++NMCPiMinusSim;
		  cMCEnergyPiMinusSim->fill(mcp->getEnergy());
		}
		case 321  : {
		  ++NMCKPlusSim;
		  cMCEnergyKPlusSim->fill(mcp->getEnergy());
		}
		case -321 : {
		  ++NMCKMinusSim; 
		  cMCEnergyKMinusSim->fill(mcp->getEnergy());
		}
		case 2212  : {
		  ++NMCProtonsSim;
		  cMCEnergyProtonsSim->fill(mcp->getEnergy());
		}
		case -2212 : {
		  ++NMCProtonBarsSim;
		  cMCEnergyProtonBarsSim->fill(mcp->getEnergy());
		}
		case 111  : {
		  ++NMCPi0Sim;
		  cMCEnergyPi0Sim->fill(mcp->getEnergy());
		}
		case 130 : {
		  ++NMCK0lSim;
		  cMCEnergyK0lSim->fill(mcp->getEnergy());
		}
		case 310  : {
		  ++NMCK0sSim;
		  cMCEnergyK0sSim->fill(mcp->getEnergy());
		}
		case 2112 : {
		  ++NMCNeutronsSim; 
		  cMCEnergyNeutronsSim->fill(mcp->getEnergy());
		}
		case 22  : {
		  ++NMCGammasSim;
		  cMCEnergyGammasSim->fill(mcp->getEnergy());
		}
		case 3122 : {
		  ++NMCLambdasSim;
		  cMCEnergyLambdasSim->fill(mcp->getEnergy());
		}
		  
		  
		default  : {
		  ++NMCRemainingSim;
		  cMCEnergyRemainingSim->fill(mcp->getEnergy());
		}
		  
		}
		
	      }
	      
	    }
	    else {
	      
	      ++NMCGen;
	      energyGen += mcp->getEnergy();
	      
	      cMCEnergyGen->fill(mcp->getEnergy());
	      
	      switch (mcp->getPDG()) {
		
	      case 11  : {
		++NMCElectronsGen;
		cMCEnergyElectronsGen->fill(mcp->getEnergy());
	      }
	      case -11 : {
		++NMCPositronsGen;
		cMCEnergyPositronsGen->fill(mcp->getEnergy());
	      }
	      case 13  : {
		++NMCMuonsGen;
		cMCEnergyMuonsGen->fill(mcp->getEnergy());
	      }
	      case -13 : {
		++NMCMuonBarsGen;
		cMCEnergyMuonBarsGen->fill(mcp->getEnergy());
	      }
	      case 15  : {
		++NMCTausGen;
		cMCEnergyTausGen->fill(mcp->getEnergy());
	      }
	      case -15 : {
		++NMCTauBarsGen;
		cMCEnergyTauBarsGen->fill(mcp->getEnergy());
	      }
	      case  12 :
	      case -12 :
	      case  14 :
	      case -14 :
	      case  16 :
	      case -16 : {
		++NMCNusGen;
		cMCEnergyNusGen->fill(mcp->getEnergy());
	      }
		
	      case 211  : {
		++NMCPiPlusGen;
		cMCEnergyPiPlusGen->fill(mcp->getEnergy());
	      }
	      case -211 : {
		++NMCPiMinusGen;
		cMCEnergyPiMinusGen->fill(mcp->getEnergy());
	      }
	      case 321  : {
		++NMCKPlusGen;
		cMCEnergyKPlusGen->fill(mcp->getEnergy());
	      }
	      case -321 : {
		++NMCKMinusGen;
		cMCEnergyKMinusGen->fill(mcp->getEnergy());
	      }
	      case 2212  : {
		++NMCProtonsGen;
		cMCEnergyProtonsGen->fill(mcp->getEnergy());
	      }
	      case -2212 : {
		++NMCProtonBarsGen;
		cMCEnergyProtonBarsGen->fill(mcp->getEnergy());
	      }
	      case 111  : {
		++NMCPi0Gen;
		cMCEnergyPi0Gen->fill(mcp->getEnergy());
	      }
	      case 130 : {
		++NMCK0lGen;
		cMCEnergyK0lGen->fill(mcp->getEnergy());
	      }
	      case 310  : {
		++NMCK0sGen;
		cMCEnergyK0sGen->fill(mcp->getEnergy());
	      }
	      case 2112 : {
		++NMCNeutronsGen;
		cMCEnergyNeutronsGen->fill(mcp->getEnergy());
	      }
	      case 22 : {
		++NMCGammasGen;
		cMCEnergyGammasGen->fill(mcp->getEnergy());
	      }
	      case 3122 : {
		++NMCLambdasGen;
		cMCEnergyLambdasGen->fill(mcp->getEnergy());
	      }
		
	      default  : {
		++NMCRemainingGen;
		cMCEnergyRemainingGen->fill(mcp->getEnergy());
	      }
		
	      }
	      
	    }
	    
	  }
	  
	  
	  if (_fillMCSim) {
	    
	    cMCNumberSim             -> fill(NMCSim);
	    cMCNumberElectronsSim    -> fill(NMCElectronsSim);
	    cMCNumberPositronsSim    -> fill(NMCPositronsSim);
	    cMCNumberMuonsSim        -> fill(NMCMuonsSim);   
	    cMCNumberMuonBarsSim     -> fill(NMCMuonBarsSim); 
	    cMCNumberTausSim         -> fill(NMCTausSim);    
	    cMCNumberTauBarsSim      -> fill(NMCTauBarsSim);  
	    cMCNumberNusSim          -> fill(NMCNusSim);
	    
	    cMCNumberPiPlusSim       -> fill(NMCPiPlusSim);
	    cMCNumberPiMinusSim      -> fill(NMCPiMinusSim);
	    cMCNumberKPlusSim        -> fill(NMCKPlusSim);   
	    cMCNumberKMinusSim       -> fill(NMCKMinusSim); 
	    cMCNumberProtonsSim      -> fill(NMCProtonsSim);    
	    cMCNumberProtonBarsSim   -> fill(NMCProtonBarsSim);  
	    cMCNumberPi0Sim          -> fill(NMCPi0Sim);
	    cMCNumberK0lSim          -> fill(NMCK0lSim);
	    cMCNumberK0sSim          -> fill(NMCK0sSim);   
	    cMCNumberNeutronsSim     -> fill(NMCNeutronsSim); 
	    cMCNumberGammasSim       -> fill(NMCGammasSim);    
	    cMCNumberLambdasSim      -> fill(NMCLambdasSim);  
	    
	    cMCNumberRemainingSim    -> fill(NMCRemainingSim);       
	    
	    cMCEnergySumSim          -> fill(energySim);
	    
	  }

	  
	  cMCNumberGen           -> fill(NMCGen);
	  cMCNumberElectronsGen  -> fill(NMCElectronsGen);
	  cMCNumberPositronsGen  -> fill(NMCPositronsGen);
	  cMCNumberMuonsGen      -> fill(NMCMuonsGen);   
	  cMCNumberMuonBarsGen   -> fill(NMCMuonBarsGen); 
	  cMCNumberTausGen       -> fill(NMCTausGen);    
	  cMCNumberTauBarsGen    -> fill(NMCTauBarsGen);
	  cMCNumberNusGen        -> fill(NMCNusGen);  
	  
	  cMCNumberPiPlusGen     -> fill(NMCPiPlusGen);
	  cMCNumberPiMinusGen    -> fill(NMCPiMinusGen);
	  cMCNumberKPlusGen      -> fill(NMCKPlusGen);   
	  cMCNumberKMinusGen     -> fill(NMCKMinusGen); 
	  cMCNumberProtonsGen    -> fill(NMCProtonsGen);    
	  cMCNumberProtonBarsGen -> fill(NMCProtonBarsGen);
	  cMCNumberPi0Gen        -> fill(NMCPi0Gen);
	  cMCNumberK0lGen        -> fill(NMCK0lGen);
	  cMCNumberK0sGen        -> fill(NMCK0sGen);   
	  cMCNumberNeutronsGen   -> fill(NMCNeutronsGen); 
	  cMCNumberGammasGen     -> fill(NMCGammasGen);    
	  cMCNumberLambdasGen    -> fill(NMCLambdasGen);
	  
	  cMCNumberRemainingGen  -> fill(NMCRemainingGen);       
	  
	  cMCEnergySumGen        -> fill(energyGen);
	  
	}
		
      }

    }
      
  }
  catch(DataNotAvailableException &e){
    std::cout << "MC particle collection not available in Check Plot processor" << std::endl ;
  };
  
  #endif
  
}








void CheckPlots::fillSimCaloCheckPlots(LCEvent * evt) {


  #ifdef MARLIN_USE_AIDA


  try {

    // numbers per event
    static AIDA::ICloud1D* cNumberSimCaloHits;
    static AIDA::ICloud1D* cEnergySimCaloHitsSum;
    
    // numbers per single particle
    static AIDA::ICloud1D* cEnergySimCaloHits;
    
    
    
    if( isFirstEvent() ) {
      
      cNumberSimCaloHits    = AIDAProcessor::histogramFactory(this)->createCloud1D( "NumberSimCaloHits", "number of SimCaloHits per event", -1 );
      cEnergySimCaloHitsSum = AIDAProcessor::histogramFactory(this)->createCloud1D( "EnergySimCaloHitsSum", "energy sum of the SimCaloHits per event", -1 );
      
      cEnergySimCaloHits    = AIDAProcessor::histogramFactory(this)->createCloud1D( "EnergySimCaloHits", "energy spectrum of the SimCaloHits per event", -1 );
      
    }
    

    int numberSimHits = 0;
    float energySimSum = 0.0;


    const std::vector< std::string >* strVec = evt->getCollectionNames() ;
    std::vector< std::string >::const_iterator name ;
    
    for( name = strVec->begin() ; name != strVec->end() ; name++) {
      
      LCCollection* col = evt->getCollection( *name ) ;

      if ( col->getTypeName() == LCIO::SIMCALORIMETERHIT ) {

	if( col != 0 ){
	 	  
	  int nSimHits = col->getNumberOfElements();
	  
	  numberSimHits += nSimHits;

	  for(int i = 0; i < nSimHits; ++i){
	    
	    SimCalorimeterHit* SimCaloHit = dynamic_cast<SimCalorimeterHit*>(col->getElementAt(i));
	    
	    float simEnergy = SimCaloHit->getEnergy();
	    
	    if (simEnergy > _simECut) {
	      cEnergySimCaloHits->fill(simEnergy);
	      energySimSum += simEnergy;
	    }
	    
	  }
	  
	}
	
      }
      
    }

    cNumberSimCaloHits->fill(numberSimHits);
    cEnergySimCaloHitsSum->fill(energySimSum);

  }
  catch(DataNotAvailableException &e){
      std::cout << "SimCalorimeterHit collection not available in Check Plot processor" << std::endl ;
  };
  

  #endif

}




void CheckPlots::fillCaloCheckPlots(LCEvent * evt) {


  #ifdef MARLIN_USE_AIDA


  try {

    // numbers per event
    static AIDA::ICloud1D* cNumberCaloHits;
    static AIDA::ICloud1D* cEnergyCaloHitsSum;
    
    // numbers per single particle
    static AIDA::ICloud1D* cEnergyCaloHits;
    
    
    
    if( isFirstEvent() ) {
      
      cNumberCaloHits    = AIDAProcessor::histogramFactory(this)->createCloud1D( "NumberCaloHits", "number of CaloHits per event", -1 );
      cEnergyCaloHitsSum = AIDAProcessor::histogramFactory(this)->createCloud1D( "EnergyCaloHitsSum", "energy sum of the CaloHits per event", -1 );
      
      cEnergyCaloHits    = AIDAProcessor::histogramFactory(this)->createCloud1D( "EnergyCaloHits", "energy spectrum of the CaloHits per event", -1 );
      
    }
    

    int numberHits = 0;
    float energySum = 0.0;


    const std::vector< std::string >* strVec = evt->getCollectionNames() ;
    std::vector< std::string >::const_iterator name ;
    
    for( name = strVec->begin() ; name != strVec->end() ; name++) {
      
      LCCollection* col = evt->getCollection( *name ) ;

      if ( col->getTypeName() == LCIO::CALORIMETERHIT ) {

	if( col != 0 ){
	 	  
	  int nHits = col->getNumberOfElements();
	  
	  numberHits += nHits;

	  for(int i = 0; i < nHits; ++i){
	    
	    CalorimeterHit* CaloHit = dynamic_cast<CalorimeterHit*>(col->getElementAt(i));
	    
	    float Energy = CaloHit->getEnergy();
	    
	    if (Energy > _ECut) {
	      cEnergyCaloHits->fill(Energy);
	      energySum += Energy;
	    }
	    
	  }
	  
	}
	
      }
      
    }

    cNumberCaloHits->fill(numberHits);
    cEnergyCaloHitsSum->fill(energySum);

  }
  catch(DataNotAvailableException &e){
      std::cout << "CalorimeterHit collection not available in Check Plot processor" << std::endl ;
  };
  

  #endif

}
