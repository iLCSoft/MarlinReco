#include "CheckPlots.h"


using namespace lcio ;
using namespace marlin ;




CheckPlots aCheckPlots ;



CheckPlots::CheckPlots() : Processor("CheckPlots") {

  _nMC = 0;
  _nMCCh = 0;
  _nMCN = 0;
  _nReco = 0;
  _nRecoCh = 0;
  _nRecoN = 0;
  
  _energyMC = 0.0;
  _energyMCCh = 0.0;
  _energyMCN = 0.0;
  _energyReco = 0.0;
  _energyRecoCh = 0.0;
  _energyRecoN = 0.0;
  

  _description = "produces Check Plots" ;


    
  registerProcessorParameter("FillMCGen","fill clouds for MC particles generated at IP",
			     _fillMCGen,
			     (int)1);
  
  registerProcessorParameter("FillMCSim","fill clouds for MC particles created during simulation",
			     _fillMCSim,
			     (int)0);

  
  registerProcessorParameter("FillSimCaloHit","fill clouds for SimCalorimeter hits",
			     _fillSimCaloHit,
			     (int)1);
  registerProcessorParameter("SimECut","energy cut for filling clouds for SimCalorimeter hits",
			     _simECut,
			     (float)0.0001);
  
  
  registerProcessorParameter("FillCaloHit","fill clouds for Calorimeter hits",
			     _fillCaloHit,
			     (int)1);
  
  registerProcessorParameter("ECut","energy cut for filling clouds for Calorimeter hits",
			     _ECut,
			     (float)0.0001);

  registerProcessorParameter("ThetaCut","Polar angle cut, given in rad. Particles with a smaller polar angle will be treated as 'lost in the beam pipe'",
			     _thetaCut,
			     (float)0.1);

    
  registerProcessorParameter("FillTracks","fill clouds for tracks",
			     _fillTracks,
			     (int)1);

  registerProcessorParameter("ColNameTracks" ,
			     "name of the Track collection" ,
			     _colNameTracks ,
			     std::string("Tracks") ); 

  registerProcessorParameter( "ColNameRelationTrackToMCP" , 
			      "name of the LC Relation collection between Tracks and MC particles"  ,
			      _colNameRelationTrackToMCP,
			      std::string("TrueTrackToMCP") );
  

  registerProcessorParameter("FillReconstructedParticles","fill clouds for ReconstructedParticles",
			     _fillReconstructedParticles,
			     (int)1);

  registerProcessorParameter("ColNameReconstructedParticles" ,
			     "name of the ReconstructedParticles collection" ,
			     _colNameReconstructedParticles ,
			     std::string("ReconstructedParticles") ); 


  registerProcessorParameter("Fill","fill clouds for comparison of MC tree and reconstructed particles",
			     _fillComparisonMCReco,
			     (int)1);

}




void CheckPlots::init() { 

  // printParameters();
 
  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z();

  createClouds();

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


  if (_fillMCGen) fillMCGenCheckPlots(evt);
  if (_fillMCSim) fillMCSimCheckPlots(evt);

  //  fillSimTrackerHitCheckPlots(evt);
  if (_fillSimCaloHit) fillSimCaloHitCheckPlots(evt);
  
  //  fillTrackerHitCheckPlots(evt);
  if (_fillCaloHit) fillCaloHitCheckPlots(evt);
  
  if (_fillTracks) fillTrackCheckPlots(evt);

  if (_fillReconstructedParticles) fillReconstructedParticlesCheckPlots(evt);

  if (_fillComparisonMCReco) fillComparisonMCRecoPlots();
  

  _nEvt ++;
  firstEvent = false;
 
}




void CheckPlots::check( LCEvent * evt ) { 
  
}




void CheckPlots::end(){ 


}




void CheckPlots::createClouds() {

  int nMax = 100 ;

  if (_fillMCSim) {

    // numbers per event
    _cMCNumberSim           = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberSim", "number of all MC particles per event (generator status != 1)", nMax );
    _cMCEnergySumSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergySumSim", "MC particle energy sum per event (generator status != 1)", nMax );
    
    _cMCNumberElectronsSim  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberElectronsSim", "number of the e+/- per event (generator status != 1)", nMax );
    _cMCNumberMuonsSim      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberMuonsSim", "number of the mu+/- per event (generator status != 1)", nMax );
    _cMCNumberTausSim       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberTausSim", "number of the tau+/- per event (generator status != 1)", nMax );
    _cMCNumberNusSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberNusSim", "number of neutrinos per event (generator status != 1)", nMax );
    
    _cMCNumberPiChSim       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPiChSim", "number of Pi+/- per event (generator status != 1)", nMax );
    _cMCNumberKChSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberKChSim", "number of K+/- per event (generator status != 1)", nMax );
    _cMCNumberProtonsSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberProtonsSim", "number of protons per event (generator status != 1)", nMax );
    _cMCNumberPi0Sim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPi0Sim", "number of Pi0 per event (generator status != 1)", nMax );
    _cMCNumberK0lSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberK0lSim", "number of K0l per event (generator status != 1)", nMax );
    _cMCNumberK0sSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberK0sSim", "number of K0s per event (generator status != 1)", nMax );
    _cMCNumberNeutronsSim   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberNeutronsSim", "number of neutrons per event (generator status != 1)", nMax );
    _cMCNumberGammasSim     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberGammasSim", "number of gammas per event (generator status != 1)", nMax );
    _cMCNumberLambda0sSim   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberLambda0sSim", "number of Lambda0s per event (generator status != 1)", nMax );
    _cMCNumberSigma0sSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberSigma0sSim", "number of Sigma0s per event (generator status != 1)", nMax );
    _cMCNumberXi0sSim       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberXi0sSim", "number of Xi0s per event (generator status != 1)", nMax );
    
    _cMCNumberRemainingSim  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberRemainingSim", "number of the remaining particles per event (generator status != 1)", nMax );
    
  }
  

  if (_fillMCGen) {

    _cMCNumberGen             = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberGen", "number of all MC particles per event (generator status == 1)", nMax );
    _cMCEnergySumGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergySumGen", "MC particle energy sum per event (generator status == 1)", nMax );

    _cMCNumberHChGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberHChGen", "number of charged particles/hadrons (generator status == 1)", nMax );
    _cMCNumberH0Gen           = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberH0Gen", "number of neutral particles/hadrons (generator status == 1)", nMax );
    _cMCNumberGGen            = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberGGen", "number of gammas and pi0s (generator status == 1)", nMax );
    _cMCFractionHChGen        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCFractionHChGen", "fraction of charged particles/hadrons (generator status == 1)", nMax );
    _cMCFractionH0Gen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCFractionH0Gen", "fraction of neutral particles/hadrons (generator status == 1)", nMax );
    _cMCFractionGGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCFractionGGen", "fraction of gammas and pi0s (generator status == 1)", nMax );

    _cMCNumberElectronsGen    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberElectronsGen", "number of the e+/- per event (generator status == 1)", nMax );
    _cMCNumberMuonsGen        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberMuonsGen", "number of the mu+/- per event (generator status == 1)", nMax );
    _cMCNumberTausGen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberTausGen", "number of the tau+/- per event (generator status == 1)", nMax );
    _cMCNumberNusGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberNusGen", "number of neutrinos per event (generator status == 1)", nMax );
    
    _cMCNumberPiChGen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPiChGen", "number of Pi+/- per event (generator status == 1)", nMax );
    _cMCNumberKChGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberKChGen", "number of K+/- per event (generator status == 1)", nMax );
    _cMCNumberProtonsGen      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberProtonsGen", "number of protons per event (generator status == 1)", nMax );
    _cMCNumberPi0Gen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberPi0Gen", "number of Pi0 per event (generator status == 1)", nMax );
    _cMCNumberK0lGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberK0lGen", "number of K0l per event (generator status == 1)", nMax );
    _cMCNumberK0sGen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberK0sGen", "number of K0s per event (generator status == 1)", nMax );
    _cMCNumberNeutronsGen     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberNeutronsGen", "number of neutrons per event (generator status == 1)", nMax );
    _cMCNumberGammasGen       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberGammasGen", "number of gammas per event (generator status == 1)", nMax );
    _cMCNumberLambda0sGen     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberLambda0sGen", "number of Lambda0s per event (generator status == 1)", nMax );
    _cMCNumberSigma0sGen      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberSigma0sSim", "number of Sigma0s per event (generator status == 1)", nMax );
    _cMCNumberXi0sGen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberXi0sSim", "number of Xi0s per event (generator status == 1)", nMax );


    _cMCNumberLostInBeamPipe  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberLostInBeamPipe", "number of particles lost in beam pipe (generator status == 1)", nMax );
    _cMCNumberRemainingGen    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCNumberRemainingGen", "number of the remaining particles per event (generator status == 1)", nMax );
      
  }  
    
      
  if (_fillMCSim) {
	
    // numbers per single particle
    _cMCEnergySim           = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergySim", "energy spectrum of all MC particles in one event (generator status != 1)", nMax );
    _cMCEnergyElectronsSim  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyElectronsSim", "energy spectrum of the e+/- in one event (generator status != 1)", nMax );
    _cMCEnergyMuonsSim      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyMuonsSim", "energy spectrum of the mu+/- in one event (generator status != 1)", nMax );
    _cMCEnergyTausSim       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyTausSim", "energy spectrum of the tau+/- in one event (generator status != 1)", nMax );
    _cMCEnergyNusSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyNusSim", "energy spectrum of the neutrinos in one event (generator status != 1)", nMax );
    
    _cMCEnergyPiChSim       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPiChSim", "energy spectrum of Pi+/- in one event (generator status != 1)", nMax );
    _cMCEnergyKChSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyKChSim", "energy spectrum of K+/- in one event (generator status != 1)", nMax );
    _cMCEnergyProtonsSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyProtonsSim", "energy spectrum of protons in one event (generator status != 1)", nMax );
    _cMCEnergyPi0Sim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPi0Sim", "energy spectrum of Pi0 in one event (generator status != 1)", nMax );
    _cMCEnergyK0lSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyK0lSim", "energy spectrum of K0l in one event (generator status != 1)", nMax );
    _cMCEnergyK0sSim        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyK0sSim", "energy spectrum of K0s in one event (generator status != 1)", nMax );
    _cMCEnergyNeutronsSim   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyNeutronsSim", "energy spectrum of neutrons in one event (generator status != 1)", nMax );
    _cMCEnergyGammasSim     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyGammasSim", "energy spectrum of gammas in one event (generator status != 1)", nMax );
    _cMCEnergyLambda0sSim   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyLambda0sSim", "energy spectrum of Lambda0s in one event (generator status != 1)", nMax );
    _cMCEnergySigma0sSim    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergySigma0sSim", "energy spectrum of Sigma0s in one event (generator status != 1)", nMax );
    _cMCEnergyXi0sSim       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyXi0sSim", "energy spectrum of Xi0s in one event (generator status != 1)", nMax );
  
    _cMCEnergyRemainingSim  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyRemainingSim", "energy spectrum of the remaining particles in one event (generator status != 1)", nMax );
    
  }

  

  if (_fillMCGen) {    
   
    _cMCEnergyGen            = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyGen", "energy spectrum of all MC particles in one event (generator status == 1)", nMax );

    _cMCEnergyHChGen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyHChGen", "energy of charged particles/hadrons (generator status == 1)", nMax );
    _cMCEnergyH0Gen          = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyH0Gen", "energy of neutral particles/hadrons (generator status == 1)", nMax );  
    _cMCEnergyGGen           = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyGGen", "energy of gammas and pi0s (generator status == 1)", nMax );  
    _cMCEnergyFractionHChGen = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyFractionHChGen", "energy fraction of charged particles/hadrons (generator status == 1)", nMax );
    _cMCEnergyFractionH0Gen  = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyFractionH0Gen", "energy fraction of neutral particles/hadrons (generator status == 1)", nMax );
    _cMCEnergyFractionGGen   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyFractionGGen", "energy fraction of gammas and pi0s (generator status == 1)", nMax );  



    _cMCEnergyElectronsGen   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyElectronsGen", "energy spectrum of the e+/- in one event (generator status == 1)", nMax );
    _cMCEnergyMuonsGen       = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyMuonsGen", "energy spectrum of the mu+/- in one event (generator status == 1)", nMax );
    _cMCEnergyTausGen        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyTausGen", "energy spectrum of the tau+/- in one event (generator status == 1)", nMax );
    _cMCEnergyNusGen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyNusGen", "energy spectrum of the neutrinos in one event (generator status == 1)", nMax );
       
    _cMCEnergyPiChGen        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPiChGen", "energy spectrum of Pi+/- in one event (generator status == 1)", nMax );
    _cMCEnergyKChGen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyKChGen", "energy spectrum of K+/- in one event (generator status == 1)", nMax );
    _cMCEnergyProtonsGen     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyProtonsGen", "energy spectrum of protons in one event (generator status == 1)", nMax );
    _cMCEnergyPi0Gen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyPi0Gen", "energy spectrum of Pi0 in one event (generator status == 1)", nMax );
    _cMCEnergyK0lGen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyK0lGen", "energy spectrum of K0l in one event (generator status == 1)", nMax );
    _cMCEnergyK0sGen         = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyK0sGen", "energy spectrum of K0s in one event (generator status == 1)", nMax );
    _cMCEnergyNeutronsGen    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyNeutronsGen", "energy spectrum of neutrons in one event (generator status == 1)", nMax );
    _cMCEnergyGammasGen      = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyGammasGen", "energy spectrum of gammas in one event (generator status == 1)", nMax );
    _cMCEnergyLambda0sGen    = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyLambda0sGen", "energy spectrum of Lambda0s in one event (generator status == 1)", nMax );
    _cMCEnergySigma0sGen     = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergySigma0sGen", "energy spectrum of Sigma0s in one event (generator status == 1)", nMax );
    _cMCEnergyXi0sGen        = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyLambda0sGen", "energy spectrum of Xi0s in one event (generator status == 1)", nMax );

    
    _cMCEnergyLostInBeamPipe = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyLostInBeamPipe", "energy of particles lost in beam pipe (generator status == 1)", nMax ); 
    _cMCEnergyRemainingGen   = AIDAProcessor::histogramFactory(this)->createCloud1D( "MCEnergyRemainingGen", "energy spectrum of the remaining particles in one event (generator status == 1)", nMax );
        
  }
   




  if (_fillSimCaloHit) {  
    
    _cNumberSimCaloHits    = AIDAProcessor::histogramFactory(this)->createCloud1D( "NumberSimCaloHits", "number of SimCaloHits per event", nMax );
    _cEnergySimCaloHitsSum = AIDAProcessor::histogramFactory(this)->createCloud1D( "EnergySimCaloHitsSum", "energy sum of the SimCaloHits per event", nMax );    
    _cEnergySimCaloHits    = AIDAProcessor::histogramFactory(this)->createCloud1D( "EnergySimCaloHits", "energy spectrum of the SimCaloHits per event", nMax );
      
  }




  if (_fillCaloHit) {  
    
    _cNumberCaloHits    = AIDAProcessor::histogramFactory(this)->createCloud1D( "NumberCaloHits", "number of CaloHits per event", nMax );
    _cEnergyCaloHitsSum = AIDAProcessor::histogramFactory(this)->createCloud1D( "EnergyCaloHitsSum", "energy sum of the CaloHits per event", nMax );    
    _cEnergyCaloHits    = AIDAProcessor::histogramFactory(this)->createCloud1D( "EnergyCaloHits", "energy spectrum of the CaloHits per event", nMax );
      
  }




  if (_fillTracks) {  

    _cNumberTracks              = AIDAProcessor::histogramFactory(this)->createCloud1D( "NumberTracks", "number of tracks per event", nMax );
    _cNumberTrackerHitsPerTrack = AIDAProcessor::histogramFactory(this)->createCloud1D( "NumberTrackerHitsPerTrack", "number of tracker hits per track", nMax );
    _cMomentumTracks            = AIDAProcessor::histogramFactory(this)->createCloud1D( "MomentumTracks", "|p| of tracks", nMax );
    _cNumberMCParticlesPerTrack = AIDAProcessor::histogramFactory(this)->createCloud1D( "NumberMCParticlesPerTrack", "number of MC particles per track", nMax );

  }



 if (_fillReconstructedParticles) {  

   _cNumberReconstructedParticles     = AIDAProcessor::histogramFactory(this)->createCloud1D( "NumberReconstructedParticles", "number of reconstructed particles per event", nMax );
   _cEnergyReconstructedParticles     = AIDAProcessor::histogramFactory(this)->createCloud1D( "EnergyReconstructedParticles", "energy spectrum of reconstructed particles", nMax );
   _cEnergySumReconstructedParticles  = AIDAProcessor::histogramFactory(this)->createCloud1D( "EnergySumReconstructedParticles", "energy sum of all reconstructed particles in an event", nMax );
   
 }



 if (_fillComparisonMCReco) {  

   _cNumberMCvsNumberReco      = AIDAProcessor::histogramFactory(this)->createCloud2D( "NumberMCvsNumberReco", "number of MC particles vs. number of reconstructed particles" , nMax );
   _cNumberMCChvsNumberRecoCh  = AIDAProcessor::histogramFactory(this)->createCloud2D( "NumberMCChvsNumberRecoCh", "number of charged MCparticles vs. number of charged reconstructed particles" , nMax );
   _cNumberMCNvsNumberRecoN    = AIDAProcessor::histogramFactory(this)->createCloud2D( "NumberMCNvsNumberRecoN", "number of neutral MCparticles vs. number of neutral reconstructed particles" , nMax );
   
   _cEnergyMCvsEnergyReco      = AIDAProcessor::histogramFactory(this)->createCloud2D( "EnergyMCvsEnergyReco", "energy of MC particles vs. energy of reconstructed particles" , nMax );
   _cEnergyMCChvsEnergyRecoCh  = AIDAProcessor::histogramFactory(this)->createCloud2D( "EnergyMCChvsEnergyRecoCh", "energy of charged MC particles vs. energy of charged reconstructed particles" , nMax );
   _cEnergyMCNvsEnergyRecoN    = AIDAProcessor::histogramFactory(this)->createCloud2D( "EnergyMCNvsEnergyRecoN", "energy of neutral MC particles vs. energy of neutral reconstructed particles" , nMax );


 }


}




void CheckPlots::fillMCGenCheckPlots(LCEvent * evt){


  if (_fillMCGen) {
    
    unsigned int NMCGen             = 0;
	  
    unsigned int NMCHChGen          = 0;
    unsigned int NMCH0Gen           = 0;
    unsigned int NMCGGen            = 0;
	  
    unsigned int NMCElectronsGen    = 0;
    unsigned int NMCMuonsGen        = 0;
    unsigned int NMCTausGen         = 0;
    unsigned int NMCNusGen          = 0;
    
    unsigned int NMCPiChGen         = 0;
    unsigned int NMCKChGen          = 0;
    unsigned int NMCProtonsGen      = 0;
    unsigned int NMCPi0Gen          = 0;
    unsigned int NMCK0lGen          = 0;
    unsigned int NMCK0sGen          = 0;
    unsigned int NMCNeutronsGen     = 0;
    unsigned int NMCGammasGen       = 0;
    unsigned int NMCLambda0sGen     = 0;
    unsigned int NMCSigma0sGen      = 0;
    unsigned int NMCXi0sGen         = 0;

    unsigned int NMCLostInBeamPipe  = 0;	  
    unsigned int NMCRemainingGen    = 0;
    
    double energyGen                = 0.0;
    double energyMCHChGen           = 0.0;
    double energyMCH0Gen            = 0.0;
    double energyMCGGen             = 0.0;
    double energyLostInBeamPipe     = 0.0;

       
    try {
      
      const std::vector< std::string >* strVec = evt->getCollectionNames() ;
      std::vector< std::string >::const_iterator name ;
    
      for( name = strVec->begin() ; name != strVec->end() ; name++) {
      
	LCCollection* col = evt->getCollection( *name ) ;

	if ( col->getTypeName() == LCIO::MCPARTICLE ) {

	  int nMCP = col->getNumberOfElements();
	  
	  for(int i = 0; i < nMCP ; ++i){
	    
	    MCParticle* mcp = dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;
	   	    
	    const double* p = mcp->getMomentum();
	    float e = mcp->getEnergy();

	    if ( mcp->getGeneratorStatus() == 1 ) {

	      double pt = hypot(p[0],p[1]);
	      double theta = atan2(pt,p[2]);

	      if ( fabs(theta) > _thetaCut ) {
		
		++NMCGen;
		energyGen += e;

                #ifdef MARLIN_USE_AIDA
		_cMCEnergyGen->fill(e);
	        #endif

		
		switch (abs(mcp->getPDG())) {
		  
		case 11  : {
		  ++NMCElectronsGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyElectronsGen->fill(e);
	          #endif

		  ++NMCHChGen;
		  energyMCHChGen += e;
		  break;
		}
		case 13  : {
		  ++NMCMuonsGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyMuonsGen->fill(e);
	          #endif

		  ++NMCHChGen;
		  energyMCHChGen += e;
		  break;
		}
		case 15  : {
		  ++NMCTausGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyTausGen->fill(e);
	          #endif

		  ++NMCHChGen;
		  energyMCHChGen += e;
		  break;
		}
		case  12 : {
		  ++NMCNusGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyNusGen->fill(e);
	          #endif
		  break;
		}
		case  14 : {
		  ++NMCNusGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyNusGen->fill(e);
	          #endif
		  break;
		}
		case  16 : {
		  ++NMCNusGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyNusGen->fill(e);
	          #endif
		  break;
		}
		  
		case 211  : {
		  ++NMCPiChGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyPiChGen->fill(e);
	          #endif

		  ++NMCHChGen;
		  energyMCHChGen += e;
		  break;
		}
		case 321  : {
		  ++NMCKChGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyKChGen->fill(e);
	          #endif

		  ++NMCHChGen;
		  energyMCHChGen += e;
		  break;
		}
		case 2212  : {
		  ++NMCProtonsGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyProtonsGen->fill(e);
	          #endif

		  ++NMCHChGen;
		  energyMCHChGen += e;
		  break;
		}
		case 111  : {
		  ++NMCPi0Gen;
		  ++NMCGGen;
		  energyMCGGen += e;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyPi0Gen->fill(e);
	          #endif
		  break;
		}
		case 130 : {
		  ++NMCK0lGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyK0lGen->fill(e);
	          #endif

		  ++NMCH0Gen;
		  energyMCH0Gen += e;
		  break;
		}
		case 310  : {
		  ++NMCK0sGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyK0sGen->fill(e);
	          #endif

		  ++NMCHChGen;
		  energyMCHChGen += e;
		  break;
		}
		case 2112 : {
		  ++NMCNeutronsGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyNeutronsGen->fill(e);
	          #endif

		  ++NMCH0Gen;
		  energyMCH0Gen += e;
		  break;
		}
		case 22 : {
		  ++NMCGammasGen;
		  ++NMCGGen;
		  energyMCGGen += e;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyGammasGen->fill(e);
	          #endif
		  break;
		}
		case 3122 : {
		  ++NMCLambda0sGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyLambda0sGen->fill(e);
	          #endif

		  ++NMCHChGen;
		  energyMCHChGen += e;
		  break;
		}
		case 3212 : {
		  ++NMCSigma0sGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergySigma0sGen->fill(e);
	          #endif

		  ++NMCHChGen;
		  energyMCHChGen += e;
		  break;
		}
		case 3322 : {
		  ++NMCXi0sGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyXi0sGen->fill(e);
	          #endif

		  ++NMCHChGen;
		  energyMCHChGen += e;
		  break;
		}

		default  : {

		  // debug
		  std::cout << "default case for for MCP (generator status == 1)" << mcp->getPDG() << std::endl;
		  
		  ++NMCRemainingGen;

                  #ifdef MARLIN_USE_AIDA
		  _cMCEnergyRemainingGen->fill(e);
	          #endif

		  ++NMCHChGen;
		  energyMCHChGen += e;
		  break;
		}
		  
		} // end of switch
		
	      }
	      else { // in beam pipe

		++NMCGen;
		energyGen += e;

		++NMCLostInBeamPipe;
		energyLostInBeamPipe += e;
		
	      }
	      
	    }
	    
	  }
	  		  
	}
	
      }
          
    }
    catch(DataNotAvailableException &e){
      std::cout << "MC particle collection not available in Check Plot processor" << std::endl ;
    };


    #ifdef MARLIN_USE_AIDA
    _cMCNumberGen            -> fill(NMCGen);
    _cMCNumberElectronsGen   -> fill(NMCElectronsGen);
    _cMCNumberMuonsGen       -> fill(NMCMuonsGen);   
    _cMCNumberTausGen        -> fill(NMCTausGen);    
    _cMCNumberNusGen         -> fill(NMCNusGen);  

    _cMCNumberHChGen         -> fill(NMCHChGen);
    _cMCNumberH0Gen          -> fill(NMCH0Gen);
    _cMCNumberGGen           -> fill(NMCGGen);
    _cMCFractionHChGen       -> fill(double(NMCHChGen)/double(NMCGen));
    _cMCFractionH0Gen        -> fill(double(NMCH0Gen)/double(NMCGen));
    _cMCFractionGGen         -> fill(double(NMCGGen)/double(NMCGen));
    _cMCEnergyHChGen         -> fill(energyMCHChGen);
    _cMCEnergyH0Gen          -> fill(energyMCH0Gen);
    _cMCEnergyGGen           -> fill(energyMCGGen);
    _cMCEnergyFractionHChGen -> fill(energyMCHChGen/energyGen);
    _cMCEnergyFractionH0Gen  -> fill(energyMCH0Gen/energyGen);
    _cMCEnergyFractionGGen   -> fill(energyMCGGen/energyGen);

    _cMCNumberPiChGen        -> fill(NMCPiChGen);
    _cMCNumberKChGen         -> fill(NMCKChGen);   
    _cMCNumberProtonsGen     -> fill(NMCProtonsGen);    
    _cMCNumberPi0Gen         -> fill(NMCPi0Gen);
    _cMCNumberK0lGen         -> fill(NMCK0lGen);
    _cMCNumberK0sGen         -> fill(NMCK0sGen);   
    _cMCNumberNeutronsGen    -> fill(NMCNeutronsGen); 
    _cMCNumberGammasGen      -> fill(NMCGammasGen);    
    _cMCNumberLambda0sGen    -> fill(NMCLambda0sGen);
    _cMCNumberSigma0sGen     -> fill(NMCSigma0sGen);
    _cMCNumberXi0sGen        -> fill(NMCXi0sGen);

    _cMCNumberLostInBeamPipe -> fill(NMCLostInBeamPipe);
    _cMCEnergyLostInBeamPipe -> fill(energyLostInBeamPipe);

    _cMCNumberRemainingGen   -> fill(NMCRemainingGen);       
    
    _cMCEnergySumGen         -> fill(energyGen);
    #endif  	


    _nMC = NMCGen - NMCLostInBeamPipe;
    _nMCCh = NMCElectronsGen + NMCMuonsGen + NMCTausGen + NMCHChGen;
    _nMCN = NMCH0Gen + NMCGGen;
    
    _energyMC = energyGen - energyLostInBeamPipe;
    _energyMCCh = energyMCHChGen;
    _energyMCN  = energyMCH0Gen + energyMCGGen;

  }

}









void CheckPlots::fillMCSimCheckPlots(LCEvent * evt){


  if (_fillMCSim) {
	  
    unsigned int NMCSim             = 0;
    
    unsigned int NMCElectronsSim    = 0;
    unsigned int NMCMuonsSim        = 0;
    unsigned int NMCTausSim         = 0;
    unsigned int NMCNusSim          = 0;
    
    unsigned int NMCPiChSim         = 0;
    unsigned int NMCKChSim          = 0;
    unsigned int NMCProtonsSim      = 0;
    unsigned int NMCPi0Sim          = 0;
    unsigned int NMCK0lSim          = 0;
    unsigned int NMCK0sSim          = 0;
    unsigned int NMCNeutronsSim     = 0;
    unsigned int NMCGammasSim       = 0;
    unsigned int NMCLambda0sSim     = 0;
    unsigned int NMCSigma0sSim      = 0;
    unsigned int NMCXi0sSim         = 0;
    
    unsigned int NMCRemainingSim    = 0;
    
    double energySim                = 0.0;
	  

    try {
      
      const std::vector< std::string >* strVec = evt->getCollectionNames() ;
      std::vector< std::string >::const_iterator name ;
      
      for( name = strVec->begin() ; name != strVec->end() ; name++) {
	
	LCCollection* col = evt->getCollection( *name ) ;
	
	if ( col->getTypeName() == LCIO::MCPARTICLE ) {
	  
	  int nMCP = col->getNumberOfElements();
	  
	  for(int i = 0; i < nMCP ; ++i){
	    
	    MCParticle* mcp = dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;
	    
	    float e = mcp->getEnergy();
	    
	    if (mcp->getGeneratorStatus() != 1 ) {
	      
	      ++NMCSim;
	      energySim += e;
		
              #ifdef MARLIN_USE_AIDA  
	      _cMCEnergySim->fill(e);
              #endif

		
	      switch (abs(mcp->getPDG())) {
		  
	      case 11  : {
		++NMCElectronsSim;

                #ifdef MARLIN_USE_AIDA  
		_cMCEnergyElectronsSim->fill(e);
                #endif
		break;
	      }
	      case 13  : {
		++NMCMuonsSim;

                #ifdef MARLIN_USE_AIDA  
		_cMCEnergyMuonsSim->fill(e);
                #endif
		break;
	      }		
	      case 15  : {
		++NMCTausSim;

                #ifdef MARLIN_USE_AIDA  
		_cMCEnergyTausSim->fill(e);
                #endif
		break;
	      }
	
	      case  12 : {
		++NMCNusSim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyNusSim->fill(e);
                #endif
		break;
	      }
	      case  14 : {
		++NMCNusSim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyNusSim->fill(e);
                #endif
		break;
	      }
	      case 16 : {
		++NMCNusSim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyNusSim->fill(e);
                #endif
		break;
	      }
		  
	      case 211  : {
		++NMCPiChSim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyPiChSim->fill(e);
                #endif
		break;
	      }
	      case 321  : {
		++NMCKChSim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyKChSim->fill(e);
                #endif
		break;
	      }
	      case 2212  : {
		++NMCProtonsSim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyProtonsSim->fill(e);
                #endif
		break;
	      }
	      case 111  : {
		++NMCPi0Sim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyPi0Sim->fill(e);
                #endif
		break;
	      }
	      case 130 : {
		++NMCK0lSim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyK0lSim->fill(e);
                #endif
		break;
	      }
	      case 310  : {
		++NMCK0sSim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyK0sSim->fill(e);
                #endif
		break;
	      }
	      case 2112 : {
		++NMCNeutronsSim;

                 #ifdef MARLIN_USE_AIDA 
		_cMCEnergyNeutronsSim->fill(e);
                #endif
		break;
	      }
	      case 22  : {
		++NMCGammasSim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyGammasSim->fill(e);
                #endif
		break;
	      }
	      case 3122 : {
		++NMCLambda0sSim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyLambda0sSim->fill(e);
                #endif
		break;
	      }
	      case 3212 : {
		++NMCSigma0sSim;
		  
                #ifdef MARLIN_USE_AIDA
		_cMCEnergySigma0sSim->fill(e);
	        #endif
		break;
	      }
	      case 3322 : {
		++NMCXi0sSim;

                #ifdef MARLIN_USE_AIDA
	        _cMCEnergyXi0sSim->fill(e);
	        #endif
		break;
	      }
	
	      default  : {

		// debug
		std::cout << "default case for for MCP (generator status != 1)" << mcp->getPDG() << std::endl;
		  
		++NMCRemainingSim;

                #ifdef MARLIN_USE_AIDA 
		_cMCEnergyRemainingSim->fill(e);
                #endif
	 	break;
	      }
		
	      } // end of switch
	      
	    }
	    
	  }
	  
	}

      }

    }
    catch(DataNotAvailableException &e){
      std::cout << "MC particle collection not available in Check Plot processor" << std::endl ;
    };
    

    #ifdef MARLIN_USE_AIDA
    _cMCNumberSim             -> fill(NMCSim);
    _cMCNumberElectronsSim    -> fill(NMCElectronsSim);
    _cMCNumberMuonsSim        -> fill(NMCMuonsSim);   
    _cMCNumberTausSim         -> fill(NMCTausSim);    
    _cMCNumberNusSim          -> fill(NMCNusSim);
    
    _cMCNumberPiChSim         -> fill(NMCPiChSim);
    _cMCNumberKChSim          -> fill(NMCKChSim);   
    _cMCNumberProtonsSim      -> fill(NMCProtonsSim);    
    _cMCNumberPi0Sim          -> fill(NMCPi0Sim);
    _cMCNumberK0lSim          -> fill(NMCK0lSim);
    _cMCNumberK0sSim          -> fill(NMCK0sSim);   
    _cMCNumberNeutronsSim     -> fill(NMCNeutronsSim); 
    _cMCNumberGammasSim       -> fill(NMCGammasSim);    
    _cMCNumberLambda0sSim     -> fill(NMCLambda0sSim);  
    _cMCNumberSigma0sSim      -> fill(NMCSigma0sSim);  
    _cMCNumberXi0sSim         -> fill(NMCXi0sSim);  

    _cMCNumberRemainingSim    -> fill(NMCRemainingSim);       
	    
    _cMCEnergySumSim          -> fill(energySim);
    #endif


	    
  }

}




void CheckPlots::fillSimCaloHitCheckPlots(LCEvent * evt) {


  if (_fillSimCaloHit) {  

    try {

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

                #ifdef MARLIN_USE_AIDA
		_cEnergySimCaloHits->fill(simEnergy);
	        #endif

		energySimSum += simEnergy;
	      }
	      
	    }
	  
	  }
	
	}
      
      }

      #ifdef MARLIN_USE_AIDA
      _cNumberSimCaloHits->fill(numberSimHits);
      _cEnergySimCaloHitsSum->fill(energySimSum);
      #endif
    
    }
    catch(DataNotAvailableException &e){std::cout << "SimCalorimeterHit collection not available in Check Plot processor" << std::endl; };
 
  }
  
}




void CheckPlots::fillCaloHitCheckPlots(LCEvent * evt) {


  if (_fillCaloHit) {  

    try {

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

                #ifdef MARLIN_USE_AIDA
		_cEnergyCaloHits->fill(Energy);
                #endif

		energySum += Energy;
	      }
	      
	    }
	    
	  }
	  
	}
	
      }

      #ifdef MARLIN_USE_AIDA
      _cNumberCaloHits->fill(numberHits);
      _cEnergyCaloHitsSum->fill(energySum);
      #endif




      //    cDifferenceEnergyCaloSumAndMCEnergy->fill(energySum-SEMC);





    }
    catch(DataNotAvailableException &e){std::cout << "CalorimeterHit collection not available in Check Plot processor" << std::endl; };
    
  }
  
}




void CheckPlots::fillTrackCheckPlots(LCEvent * evt) {


  if (_fillTracks) {

    try {
    
      std::vector< std::string >::const_iterator iter;
      const std::vector< std::string >* ColNames = evt->getCollectionNames();
    
      for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
      
	LCCollection* col = evt->getCollection( *iter ) ;
  
	if ( (col->getTypeName() == LCIO::TRACK) && (*iter == _colNameTracks) ) {

	  int nTracks = col->getNumberOfElements();
	  
          #ifdef MARLIN_USE_AIDA
	  _cNumberTracks->fill(nTracks);	
	  #endif

	  for(int j=0; j<nTracks; ++j){
	  
	    Track* track = dynamic_cast<Track*>(col->getElementAt(j));

	    int nTrackerHits = track->getTrackerHits().size();
	    
            #ifdef MARLIN_USE_AIDA
	    _cNumberTrackerHitsPerTrack->fill(nTrackerHits);	
     	    #endif

	    const double absP = MarlinUtil::getAbsMomentum(track,_bField);
	  
            #ifdef MARLIN_USE_AIDA
	    _cMomentumTracks->fill(absP);	
	    #endif

	    try {

	      LCCollection* LCRcolTracks = evt->getCollection(_colNameRelationTrackToMCP);

	      LCRelationNavigator* navTracks = new LCRelationNavigator(LCRcolTracks);
	      const LCObjectVec& relMCParticlesToTrack = navTracks->getRelatedToObjects(track); 

	      int nOfRelatedMCParticles = relMCParticlesToTrack.size();

              #ifdef MARLIN_USE_AIDA
	      _cNumberMCParticlesPerTrack->fill(nOfRelatedMCParticles);
              #endif

	    }
	    catch(DataNotAvailableException &e){std::cout << "no valid LCRelation between track and MC particle in event " << _nEvt << std::endl; };
	    
	  }
	  
	}
	
      }
      
    }
    catch(DataNotAvailableException &e){std::cout << "no valid track collection available in event " << _nEvt << std::endl; };
  }

}




void CheckPlots::fillReconstructedParticlesCheckPlots(LCEvent * evt) {


  if (_fillReconstructedParticles) {

    int nReco   = 0;
    int nRecoCh = 0;
    int nRecoN  = 0;

    double energyReco   = 0.0;
    double energyRecoCh = 0.0;
    double energyRecoN  = 0.0;
    

    try {

      int nReconstructedParticles = 0;
      double energySumReconstructedParticles = 0.0;
    
      std::vector< std::string >::const_iterator iter;
      const std::vector< std::string >* ColNames = evt->getCollectionNames();
    
      for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
      
	LCCollection* col = evt->getCollection( *iter ) ;
      
	if ( (col->getTypeName() == LCIO::RECONSTRUCTEDPARTICLE) && (*iter == _colNameReconstructedParticles) ) {

	  nReconstructedParticles += col->getNumberOfElements();
	  
	  for(int j=0; j<nReconstructedParticles; ++j){
	  
	    ReconstructedParticle* recoParticle = dynamic_cast<ReconstructedParticle*>(col->getElementAt(j));

	    double energyReconstructedParticle = recoParticle->getEnergy();
	    
	    energySumReconstructedParticles += energyReconstructedParticle;

            #ifdef MARLIN_USE_AIDA
            _cEnergyReconstructedParticles->fill(energyReconstructedParticle);
            #endif
	    

	    if ( (recoParticle->getTracks().size()) > 0 ) {
	      
	      ++nRecoCh;
	      energyRecoCh += energyReconstructedParticle;

	    }
	    else {

	      ++nRecoN;
	      energyRecoN += energyReconstructedParticle;
	      
	    }

	  }
	  
	}
	
      }
      	  
      #ifdef MARLIN_USE_AIDA
      _cNumberReconstructedParticles->fill(nReconstructedParticles);
      _cEnergySumReconstructedParticles->fill(energySumReconstructedParticles);
      #endif
	
      nReco = nReconstructedParticles;
      energyReco = energySumReconstructedParticles;

    }
    catch(DataNotAvailableException &e){std::cout << "no valid reconstructed particle collection available in event " << _nEvt << std::endl; };
  

    _nReco = nReco;
    _nRecoCh = nRecoCh;
    _nRecoN = nRecoN;
    
    _energyReco = energyReco;
    _energyRecoCh = energyRecoCh;
    _energyRecoN = energyRecoN;

  }

}




void CheckPlots::fillComparisonMCRecoPlots() {

  if (_fillComparisonMCReco) {

    if ( _fillMCGen && _fillReconstructedParticles ) {

      _cNumberMCvsNumberReco->fill(_nMC,_nReco);
      _cNumberMCChvsNumberRecoCh->fill(_nMCCh,_nRecoCh);
      _cNumberMCNvsNumberRecoN->fill(_nMCN,_nRecoN);
      
      _cEnergyMCvsEnergyReco->fill(_energyMC,_energyReco);
      _cEnergyMCChvsEnergyRecoCh->fill(_energyMCCh,_energyRecoCh);
      _cEnergyMCNvsEnergyRecoN->fill(_energyMCN,_energyRecoN);

    }

  }

}
