#include "TJjetsPFOAnalysisProcessor.h"

TJjetsPFOAnalysisProcessor aTJjetsPFOAnalysisProcessor ;


TJjetsPFOAnalysisProcessor::TJjetsPFOAnalysisProcessor() :
    Processor("TJjetsPFOAnalysisProcessor"),
    // Parameters for PFOAnalysis
    m_nRun(0),
    m_nEvt(0),
    m_nRunSum(0),
    m_nEvtSum(0),
    m_printing(0),
    m_lookForQuarksWithMotherZ(0),
    m_mcPfoSelectionRadius(500.f),
    m_mcPfoSelectionMomentum(0.01f),
    m_mcPfoSelectionLowEnergyNPCutOff(1.2f),
    m_nPfosTotal(0),
    m_nPfosNeutralHadrons(0),
    m_nPfosPhotons(0),
    m_nPfosTracks(0),
    m_pfoEnergyTotal(0.f),
    m_pfoEnergyNeutralHadrons(0.f),
    m_pfoEnergyPhotons(0.f),
    m_pfoEnergyTracks(0.f),
    m_pfoECalToEmEnergy(0.f),
    m_pfoECalToHadEnergy(0.f),
    m_pfoHCalToEmEnergy(0.f),
    m_pfoHCalToHadEnergy(0.f),
    m_pfoMuonToEnergy(0.f),
    m_pfoOtherEnergy(0.f),
    m_pfoMassTotal(0.f),
    m_nPfoTargetsTotal(0),
    m_nPfoTargetsNeutralHadrons(0),
    m_nPfoTargetsPhotons(0),
    m_nPfoTargetsTracks(0),
    m_pfoTargetsEnergyTotal(0.f),
    m_pfoTargetsEnergyNeutralHadrons(0.f),
    m_pfoTargetsEnergyPhotons(0.f),
    m_pfoTargetsEnergyTracks(0.f),
    m_mcEnergyENu(0.f),
    m_mcEnergyFwd(0.f),
    m_eQQ(-99.f),
    m_eQ1(-99.f),
    m_eQ2(-99.f),
    m_costQQ(-99.f),
    m_costQ1(-99.f),
    m_costQ2(-99.f),
    m_mQQ(-99.f),
    m_thrust(-99.f),
    m_qPdg(-99.f),
    m_pTFile(NULL),
    m_pTTree(NULL),
    m_hPfoEnergySum(NULL),
    m_hPfoEnergySumL7A(NULL)
{
    _description = "TJjetsPFOAnalysisProcessor: Basically a copy-pasted version of the PFOAnalysis processor working on individual TrueJet jets" ;




    // Standard input collections
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
            "InputAllPFOsCollection",
            "Name of the PFOs collection",
            _colAllPFOs,
            std::string("PandoraPFOs")
    );

    registerInputCollection( LCIO::MCPARTICLE,
            "MCParticleCollection",
            "Name of the MC particle collection",
            _colMC,
            std::string("MCParticle")
    );


    // Collection for TrueJet

  	registerInputCollection( LCIO::LCRELATION,
      		"RecoMCTruthLink",
      		"Name of the RecoMCTruthLink input collection"  ,
  			_recoMCTruthLink,
  			std::string("RecoMCTruthLink") ) ;

  	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
  			"TrueJets" ,
  			"Name of the TrueJetCollection output collection"  ,
  			_trueJetCollectionName ,
  			std::string("TrueJets") ) ;

  	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
  			"FinalColourNeutrals" ,
  			"Name of the FinalColourNeutralCollection output collection"  ,
  			_finalColourNeutralCollectionName ,
  			std::string("FinalColourNeutrals") ) ;

  	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
  			"InitialColourNeutrals" ,
  			"Name of the InitialColourNeutralCollection output collection"  ,
  			_initialColourNeutralCollectionName ,
  			std::string("InitialColourNeutrals") ) ;

  	 registerInputCollection( LCIO::LCRELATION,
  			"TrueJetPFOLink" ,
  			"Name of the TrueJetPFOLink output collection"  ,
  			_trueJetPFOLink,
  			std::string("TrueJetPFOLink") ) ;

  	registerInputCollection( LCIO::LCRELATION,
  			"TrueJetMCParticleLink" ,
  			"Name of the TrueJetMCParticleLink output collection"  ,
  			_trueJetMCParticleLink,
  			std::string("TrueJetMCParticleLink") ) ;

  	registerInputCollection( LCIO::LCRELATION,
  			"FinalElementonLink" ,
  			"Name of the  FinalElementonLink output collection"  ,
  			_finalElementonLink,
  			std::string("FinalElementonLink") ) ;

  	registerInputCollection( LCIO::LCRELATION,
  			"InitialElementonLink" ,
  			"Name of the  InitialElementonLink output collection"  ,
  			_initialElementonLink,
  			std::string("InitialElementonLink") ) ;

  	registerInputCollection( LCIO::LCRELATION,
  			"FinalColourNeutralLink" ,
  			"Name of the  FinalColourNeutralLink output collection"  ,
  			_finalColourNeutralLink,
  			std::string("FinalColourNeutralLink") ) ;

  	registerInputCollection( LCIO::LCRELATION,
  			"InitialColourNeutralLink" ,
  			"Name of the  InitialColourNeutralLink output collection"  ,
  			_initialColourNeutralLink,
  			std::string("InitialColourNeutralLink") ) ;



    // Collections for PFOAnalysis

    registerProcessorParameter(
        "LookForQuarksWithMotherZ",
        "Flag to look for quarks with mother Z",
        m_lookForQuarksWithMotherZ,
        int(0));

    registerProcessorParameter(
        "MCPfoSelectionRadius",
        "MC pfo selection radius",
        m_mcPfoSelectionRadius,
        float(500.f));

    registerProcessorParameter(
        "MCPfoSelectionMomentum",
        "MC pfo selection momentum",
        m_mcPfoSelectionMomentum,
        float(0.01f));

    registerProcessorParameter(
        "MCPfoSelectionLowEnergyNPCutOff",
        "MC pfo selection neutron and proton low energy cut-off",
        m_mcPfoSelectionLowEnergyNPCutOff,
        float(1.0f));

    registerProcessorParameter(
        "RootFile",
        "Name of the output root file",
        m_rootFile,
        std::string("PFOAnalysis.root"));

}



void TJjetsPFOAnalysisProcessor::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  m_nRun = 0;
  m_nEvt = 0;
  m_nRunSum = 0;
  m_nEvtSum = 0;
  this->Clear();

  m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
  makeNTuple() ;
}


void TJjetsPFOAnalysisProcessor::makeNTuple() {
  m_pTTree = new TTree("PfoAnalysisTree", "PfoAnalysisTree");
  m_pTTree->SetDirectory(m_pTFile);
  m_pTTree->Branch("run", &m_nRun, "run/I");
  m_pTTree->Branch("event", &m_nEvt, "event/I");
  m_pTTree->Branch("jet", &m_nJet, "jet/I");

  m_pTTree->Branch("jetInitElPDG", &m_jetInitElPDG, "jetInitElPDG/I");
  m_pTTree->Branch("jetFinElPDG", &m_jetInitElPDG, "jetFinElPDG/I");

  m_pTTree->Branch("nPfosTotal", &m_nPfosTotal, "nPfosTotal/I");
  m_pTTree->Branch("nPfosNeutralHadrons", &m_nPfosNeutralHadrons, "nPfosNeutralHadrons/I");
  m_pTTree->Branch("nPfosPhotons", &m_nPfosPhotons, "nPfosPhotons/I");
  m_pTTree->Branch("nPfosTracks", &m_nPfosTracks, "nPfosTracks/I");
  m_pTTree->Branch("pfoEnergyTotal", &m_pfoEnergyTotal, "pfoEnergyTotal/F");
  m_pTTree->Branch("pfoEnergyNeutralHadrons", &m_pfoEnergyNeutralHadrons, "pfoEnergyNeutralHadrons/F");
  m_pTTree->Branch("pfoEnergyPhotons", &m_pfoEnergyPhotons, "pfoEnergyPhotons/F");
  m_pTTree->Branch("pfoEnergyTracks", &m_pfoEnergyTracks, "pfoEnergyTracks/F");
  m_pTTree->Branch("pfoECalToEmEnergy", &m_pfoECalToEmEnergy, "pfoECalToEmEnergy/F");
  m_pTTree->Branch("pfoECalToHadEnergy", &m_pfoECalToHadEnergy, "pfoECalToHadEnergy/F");
  m_pTTree->Branch("pfoHCalToEmEnergy", &m_pfoHCalToEmEnergy, "pfoHCalToEmEnergy/F");
  m_pTTree->Branch("pfoHCalToHadEnergy", &m_pfoHCalToHadEnergy, "pfoHCalToHadEnergy/F");
  m_pTTree->Branch("pfoOtherEnergy", &m_pfoOtherEnergy, "pfoOtherEnergy/F");
  m_pTTree->Branch("pfoMuonToEnergy", &m_pfoMuonToEnergy, "pfoMuonToEnergy/F");
  m_pTTree->Branch("pfoMassTotal", &m_pfoMassTotal, "pfoMassTotal/F");
  m_pTTree->Branch("pfoEnergies", &m_pfoEnergies);
  m_pTTree->Branch("pfoPx", &m_pfoPx);
  m_pTTree->Branch("pfoPy", &m_pfoPy);
  m_pTTree->Branch("pfoPz", &m_pfoPz);
  m_pTTree->Branch("pfoCosTheta", &m_pfoCosTheta);
  m_pTTree->Branch("pfoTargetEnergies", &m_pfoTargetEnergies);
  m_pTTree->Branch("pfoTargetPx", &m_pfoTargetPx);
  m_pTTree->Branch("pfoTargetPy", &m_pfoTargetPy);
  m_pTTree->Branch("pfoTargetPz", &m_pfoTargetPz);
  m_pTTree->Branch("pfoTargetCosTheta", &m_pfoTargetCosTheta);
  m_pTTree->Branch("pfoPdgCodes", &m_pfoPdgCodes);
  m_pTTree->Branch("pfoTargetPdgCodes", &m_pfoTargetPdgCodes);
  m_pTTree->Branch("nPfoTargetsTotal", &m_nPfoTargetsTotal, "nPfoTargetsTotal/I");
  m_pTTree->Branch("nPfoTargetsNeutralHadrons", &m_nPfoTargetsNeutralHadrons, "nPfoTargetsNeutralHadrons/I");
  m_pTTree->Branch("nPfoTargetsPhotons", &m_nPfoTargetsPhotons, "nPfoTargetsPhotons/I");
  m_pTTree->Branch("nPfoTargetsTracks", &m_nPfoTargetsTracks, "nPfoTargetsTracks/I");
  m_pTTree->Branch("pfoTargetsEnergyTotal", &m_pfoTargetsEnergyTotal, "pfoTargetsEnergyTotal/F");
  m_pTTree->Branch("pfoTargetsEnergyNeutralHadrons", &m_pfoTargetsEnergyNeutralHadrons, "pfoTargetsEnergyNeutralHadrons/F");
  m_pTTree->Branch("pfoTargetsEnergyPhotons", &m_pfoTargetsEnergyPhotons, "pfoTargetsEnergyPhotons/F");
  m_pTTree->Branch("pfoTargetsEnergyTracks", &m_pfoTargetsEnergyTracks, "pfoTargetsEnergyTracks/F");
  m_pTTree->Branch("mcEnergyENu", &m_mcEnergyENu, "mcEnergyENu/F");
  m_pTTree->Branch("mcEnergyFwd", &m_mcEnergyFwd, "mcEnergyFwd/F");
  m_pTTree->Branch("eQQ", &m_eQQ, "eQQ/F");
  m_pTTree->Branch("eQ1", &m_eQ1, "eQ1/F");
  m_pTTree->Branch("eQ2", &m_eQ2, "eQ2/F");
  m_pTTree->Branch("costQQ", &m_costQQ, "costQQ/F");
  m_pTTree->Branch("costQ1", &m_costQ1, "costQ1/F");
  m_pTTree->Branch("costQ2", &m_costQ2, "costQ2/F");
  m_pTTree->Branch("mQQ", &m_mQQ, "mQQ/F");
  m_pTTree->Branch("thrust", &m_thrust, "thrust/F");
  m_pTTree->Branch("qPdg", &m_qPdg, "qPdg/I");

  m_hPfoEnergySum = new TH1F("fPFA", "total pfo energy", 10000, 0., 5000.);
  m_hPfoEnergySumL7A = new TH1F("fPFA_L7A", "total pfo energy < 0.7 A", 10000, 0., 5000.);
  m_hPfoEnergySum->SetDirectory(m_pTFile);
  m_hPfoEnergySumL7A->SetDirectory(m_pTFile);
}



void TJjetsPFOAnalysisProcessor::processRunHeader( LCRunHeader* /*run*/) {
  m_nRun = 0;
  m_nEvt = 0;
  ++m_nRunSum;
}





void TJjetsPFOAnalysisProcessor::processEvent( LCEvent * event ) {

  // tj is a pointer to a Trujet_Parser, with the data of this processor object:
  TrueJet_Parser* tj= this ;
  // this method gets all the collections needed + initialises a few convienent variables.
  tj->getall(event);


  m_nRun = event->getRunNumber();
  m_nEvt = event->getEventNumber();
  ++m_nEvtSum;

  if ((m_nEvtSum % 100) == 0)
    std::cout << " processed events: " << m_nEvtSum << std::endl;

  m_jets.clear();
  this->findTrueJetParticles(event);

  for (int i_jet=0; i_jet<njets(); i_jet++) {
    m_nJet = i_jet;
    this->Clear();

    if( initial_elementon(i_jet) != NULL) {
      m_jetInitElPDG = initial_elementon(i_jet)->getPDG();
    } else {
      m_jetInitElPDG = 0;
    }
    if( final_elementon(i_jet) != NULL) {
      m_jetFinElPDG = final_elementon(i_jet)->getPDG();
    } else {
      m_jetFinElPDG = 0;
    }


    this->ExtractCollections(m_jets[i_jet]);
    this->MakeQuarkVariables(m_jets[i_jet]);
    this->PerformPfoAnalysis();

    m_pTTree->Fill();
    
    delete m_jets[i_jet];
  }
}



void TJjetsPFOAnalysisProcessor::check( LCEvent * /*event*/ ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TJjetsPFOAnalysisProcessor::end(){
  if (m_printing > -1)
  {
    std::cout << "PfoAnalysis::end() " << this->name() << " processed " << m_nEvtSum << " events in " << m_nRunSum << " runs " << std::endl
    << "Rootfile: " << m_rootFile.c_str() << std::endl;
  }

  m_pTFile->cd();
  m_pTTree->Write();
  m_hPfoEnergySum->Write();
  m_hPfoEnergySumL7A->Write();

  m_pTFile->Close();
  delete m_pTFile;
}
