#include "SLDCorrection.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include "EVENT/MCParticle.h"
#include "EVENT/Vertex.h"
#include "EVENT/ReconstructedParticle.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "IMPL/ParticleIDImpl.h"
#include <IMPL/VertexImpl.h>
#include "UTIL/PIDHandler.h"
#include "marlin/VerbosityLevels.h"
#include <GeometryUtil.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TF1.h"
#include "TTree.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TRatioPlot.h"
#include "TAxis.h"
#include "TLine.h"
#include "DDMarlinCED.h"

#include "flightDirection.h"
#include "AssignParticlestoSLD.h"
#include "FindParticle.h"

using namespace lcio ;
using namespace marlin ;

SLDCorrection aSLDCorrection;

SLDCorrection::SLDCorrection() :

	Processor("SLDCorrection"),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),
	m_Bfield(0.f),
	foundFlightDirection(true),
	m_nSLDecayOfBHadron(0),
	m_nSLDecayOfCHadron(0),
	m_nSLDecayOfTauLepton(0),
	m_nSLDecayTotal(0),
	m_nSLDecayToElectron(0),
	m_nSLDecayToMuon(0),
	m_nSLDecayToTau(0),
	m_nTauNeutrino(0),
	m_nNeutrino(0),
	m_nChargedPFOwoTrack(0),
	n_SLDStatus(0),
	n_NuPxResidual(0),
	n_NuPyResidual(0),
	n_NuPzResidual(0),
	n_NuEResidual(0),
	n_NuPxNormalizedResidual(0),
	n_NuPyNormalizedResidual(0),
	n_NuPzNormalizedResidual(0),
	n_NuENormalizedResidual(0),
	n_secondaryVertex(0)
{
	_description = "SLDCorrection finds semi-leptonic decays within jets and performs a correction to 4-momentum of the jet due to the missing neutrino(s)";

	registerInputCollection(	LCIO::MCPARTICLE,
					"MCParticleCollection" ,
					"Name of the MCParticle collection"  ,
					m_mcParticleCollection,
					std::string("MCParticle")
				);

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"PfoCollection",
					"Name of input pfo collection",
					m_inputPfoCollection,
					std::string("PandoraPFOs")
				);

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"JetCollection",
					"Name of input jet collection",
					m_inputJetCollection,
					std::string("Durham_nJets")
				);

	registerInputCollection(	LCIO::VERTEX,
					"PrimaryVertex",
					"Name of Primary Vertex Collection",
					m_inputPrimaryVertex,
					std::string("PrimaryVertex")
				);

	registerInputCollection(	LCIO::VERTEX,
					"BuildUpVertex",
					"Name of BuildUp Vertex Collection",
					m_inputBuildUpVertex,
					std::string("BuildUpVertex")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"RecoMCTruthLinkCollection",
					"Name of input RecoMCTruthLink Collection",
					m_RecoMCTruthLinkCollection,
					std::string("RecoMCTruthLinkCollection")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthRecoLinkCollection",
					"Name of input MCTruthRecoLink Collection",
					m_MCTruthRecoLinkCollection,
					std::string("MCTruthRecoLinkCollection")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"ClusterMCTruthLinkCollection",
					"Name of input m_ClusterMCTruthLink Collection",
					m_ClusterMCTruthLinkCollection,
					std::string("ClusterMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthClusterLinkCollection",
					"Name of input MCTruthClusterLink Collection",
					m_MCTruthClusterLinkCollection,
					std::string("MCTruthClusterLink")
				);

	registerOutputCollection(	LCIO::VERTEX,
					"SemiLeptonicDecayVertex",
					"Name of Semi-Leptonic Decay Vertices Collection",
					m_SLDVertex,
					std::string("SemiLeptonicDecayVertex")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"SemiLeptonicDecayVertexRP",
					"Name of Semi-Leptonic Decay Vertices Reconstructed Particle Collection",
					m_SLDVertexRP,
					std::string("SemiLeptonicDecayVertexRP")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"ReconstructedNeutrino",
					"Name of Reconstructed Neutrino Collection",
					m_reconstructedNeutrino,
					std::string("ReconstructedNeutrino")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"JetSLDLinkName",
					"Name of the JetSemiLeptonicDecayLinkName output collection",
					m_JetSLDLinkName,
					std::string("JetSLDLinkName")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"SLDJetLinkName",
					"Name of the SemiLeptonicDecayJetLinkName output collection",
					m_SLDJetLinkName,
					std::string("SLDJetLinkName")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"mcNurecoNuLinkName",
					"Name of the trueNeutrino-reconstructedNeutrino output Link collection",
					m_mcNurecoNuLinkName,
					std::string("mcNurecoNuLinkName")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"recoNumcNuLinkName",
					"Name of the trueNeutrino-reconstructedNeutrino output Link collection",
					m_recoNumcNuLinkName,
					std::string("recoNumcNuLinkName")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"SLDNuLinkName",
					"Name of the NeutrinoSemiLeptonicDecayLinkName output collection",
					m_SLDNuLinkName,
					std::string("SLDNuLinkName")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"NuSLDLinkName",
					"Name of the SemiLeptonicDecayNeutrinoLinkName output collection",
					m_NuSLDLinkName,
					std::string("NuSLDLinkName")
				);

	registerProcessorParameter(	"includeBSLD",
					"do correction for semi-leptonic decays of B-Hadrons",
					m_includeBSLD,
					bool(true)
				);

	registerProcessorParameter(	"includeCSLD",
					"do correction for semi-leptonic decays of C-Hadrons",
					m_includeCSLD,
					bool(true)
				);

	registerProcessorParameter(	"includeTSLD",
					"do correction for semi-leptonic decays of Tau-Leptons",
					m_includeTSLD,
					bool(true)
				);

	registerProcessorParameter(	"cheatSLDLeptons",
					"Cheat semi-leptonic decays lepton from MCTruth",
					m_cheatSLDLeptons,
					bool(true)
				);

	registerProcessorParameter(	"cheatFlightDirection",
					"Cheat Flight direction of mother hadron",
					m_cheatFlightDirection,
					bool(true)
				);

	registerProcessorParameter(	"vertexingScenario",
					"Scenario for finding flight direction of mother hadron: 1 = default , 2 = assign jet axis , 3 = assign flight direction of leading particle in the jet",
					m_vertexingScenario,
					int(1)
				);

	registerProcessorParameter(	"cheatLepton4momentum",
					"Cheat FourMomentum of lepton in semi-leptonic decays",
					m_cheatLepton4momentum,
					bool(true)
				);

	registerProcessorParameter(	"cheatCharged4momentum",
					"Cheat FourMomentum of charged visibles in semi-leptonic decays",
					m_cheatCharged4momentum,
					bool(true)
				);

	registerProcessorParameter(	"cheatNeutral4momentum",
					"Cheat FourMomentum of neutral visibles in semi-leptonic decays",
					m_cheatNeutral4momentum,
					bool(true)
				);

	registerProcessorParameter(	"cheatPVAcharged",
					"Cheat Particle ID for Charged Decay Products",
					m_cheatPVAcharged,
					bool(true)
				);

	registerProcessorParameter(	"chargedCosAcceptanceAngleSLD4",
					"Acceptance angle for charged PFOs to be assigned to the semi-leptonic decay when lepton is in secondary vertex",
					m_chargedCosAcceptanceAngleSLD4,
					float(0.0)
				);

	registerProcessorParameter(	"chargedCosAcceptanceAngleSLD5",
					"Acceptance angle for charged PFOs to be assigned to the semi-leptonic decay when lepton is with third vertex",
					m_chargedCosAcceptanceAngleSLD5,
					float(0.0)
				);

	registerProcessorParameter(	"cheatPVAneutral",
					"Cheat Particle ID for Neutral Decay Products",
					m_cheatPVAneutral,
					bool(true)
				);

	registerProcessorParameter(	"cheatSolutionSign",
					"Cheat the signe of solution for the neutrino correction",
					m_cheatSolutionSign,
					bool(false)
				);

	registerProcessorParameter(	"neutralCosAcceptanceAngle",
					"Acceptance angle for neutral PFOs to be assigned to the semi-leptonic decay",
					m_neutralCosAcceptanceAngle,
					float(0.0)
				);

	registerProcessorParameter(	"BSLDChargedSLD4InvMassCut",
					"Cut on Invariant mass of visible charged decay products for semi-leptonic decay of B-Hadron when lepton is in secondary vertex",
					m_BSLDChargedSLD4InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"BSLDChargedSLD5InvMassCut",
					"Cut on Invariant mass of visible charged decay products for semi-leptonic decay of B-Hadron when lepton is with third vertex",
					m_BSLDChargedSLD5InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"BSLDNeutralSLD4InvMassCut",
					"Cut on Invariant mass of visible neutral decay products for semi-leptonic decay of B-Hadron when lepton is in secondary vertex",
					m_BSLDNeutralSLD4InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"BSLDNeutralSLD5InvMassCut",
					"Cut on Invariant mass of visible neutral decay products for semi-leptonic decay of B-Hadron when lepton is with third vertex",
					m_BSLDNeutralSLD5InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"CSLDChargedSLD4InvMassCut",
					"Cut on Invariant mass of visible charged decay products for semi-leptonic decay of C-Hadron when lepton is in secondary vertex",
					m_CSLDChargedSLD4InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"CSLDChargedSLD5InvMassCut",
					"Cut on Invariant mass of visible charged decay products for semi-leptonic decay of C-Hadron when lepton is with third vertex",
					m_CSLDChargedSLD5InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"CSLDNeutralSLD4InvMassCut",
					"Cut on Invariant mass of visible neutral decay products for semi-leptonic decay of C-Hadron when lepton is in secondary vertex",
					m_CSLDNeutralSLD4InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"CSLDNeutralSLD5InvMassCut",
					"Cut on Invariant mass of visible neutral decay products for semi-leptonic decay of C-Hadron when lepton is with third vertex",
					m_CSLDNeutralSLD5InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"nIterFlightDirCorrection",
					"Number of iterations for correcting flight direction of CHARGED parent hadron",
					m_nIterFlightDirCorrection,
					int(1)
				);

	registerProcessorParameter(	"BSLD4SigmaAlpha",
					"Angular uncertainty (radian) due to flight direction reconstruction for semi-leptonic decay of B-Hadron when leptin is in 2nd vertex",
					m_BSLD4SigmaAlpha,
					float(0.004)
				);

	registerProcessorParameter(	"BSLD5SigmaAlpha",
					"Angular uncertainty (radian) due to flight direction reconstruction for semi-leptonic decay of B-Hadron when leptin is with 3rd vertex",
					m_BSLD5SigmaAlpha,
					float(0.010)
				);

	registerProcessorParameter(	"BSLD4SigmaECPVA",
					"Energy uncertainty (GeV) due to Charged Particle to Vertex Association (CPVA) for semi-leptonic decay of B-Hadron when leptin is in 2nd vertex",
					m_BSLD4SigmaECPVA,
					float(0.0)
				);

	registerProcessorParameter(	"BSLD4SigmaENPVA",
					"Energy uncertainty (GeV) due to Neutral Particle to Vertex Association (NPVA) for semi-leptonic decay of B-Hadron when leptin is in 2nd vertex",
					m_BSLD4SigmaENPVA,
					float(0.0)
				);

	registerProcessorParameter(	"BSLD4SigmaAlphaCPVA",
					"Angular uncertainty (radian) due to Charged Particle to Vertex Association (NPVA) for semi-leptonic decay of B-Hadron when leptin is in 2nd vertex",
					m_BSLD4SigmaAlphaCPVA,
					float(0.0)
				);

	registerProcessorParameter(	"BSLD4SigmaAlphaNPVA",
					"Angular uncertainty (radian) due to Neutral Particle to Vertex Association (NPVA) for semi-leptonic decay of B-Hadron when leptin is in 2nd vertex",
					m_BSLD4SigmaAlphaNPVA,
					float(0.0)
				);

	registerProcessorParameter(	"BSLD5SigmaECPVA",
					"Energy uncertainty (GeV) due to Charged Particle to Vertex Association (CPVA) for semi-leptonic decay of B-Hadron when leptin is with 3rd vertex",
					m_BSLD5SigmaECPVA,
					float(0.0)
				);

	registerProcessorParameter(	"BSLD5SigmaENPVA",
					"Energy uncertainty (GeV) due to Neutral Particle to Vertex Association (NPVA) for semi-leptonic decay of B-Hadron when leptin is with 3rd vertex",
					m_BSLD5SigmaENPVA,
					float(0.0)
				);

	registerProcessorParameter(	"BSLD5SigmaAlphaCPVA",
					"Angular uncertainty (radian) due to Charged Particle to Vertex Association (NPVA) for semi-leptonic decay of B-Hadron when leptin is with 3rd vertex",
					m_BSLD5SigmaAlphaCPVA,
					float(0.0)
				);

	registerProcessorParameter(	"BSLD5SigmaAlphaNPVA",
					"Angular uncertainty (radian) due to Neutral Particle to Vertex Association (NPVA) for semi-leptonic decay of B-Hadron when leptin is with 3rd vertex",
					m_BSLD5SigmaAlphaNPVA,
					float(0.0)
				);

	registerProcessorParameter(	"sigmaENu",
					"Energy uncertainty (GeV) for reconstructed neutrino",
					m_sigmaENu,
					float(4.0)
				);

	registerProcessorParameter(	"sigmaAlphaNu",
					"Angular uncertainty (radian) for reconstructed neutrino direction",
					m_sigmaAlphaNu,
					float(0.1)
				);

	registerProcessorParameter(	"displayEvent",
					"Display recoLepton, downstraem vertex and RP and reconstructed secondary vertex in Event Display",
					m_displayEvent,
					bool(true)
				);

	registerProcessorParameter(	"fillRootTree",
					"Fill root tree to check processor performance",
					m_fillRootTree,
					bool(true)
				);

	registerProcessorParameter(	"traceEvent",
					"trace events with large discrepancy",
					m_traceEvent,
					bool(false)
				);

	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("Output.root")
				);

	registerProcessorParameter(	"BsldMode",
					"Event selection based on semi-leptonic decays of B-hadrons",
					m_BSLDMode,
					int(0)
				);

	registerProcessorParameter(	"CsldMode",
					"Event selection based on semi-leptonic decays of C-hadrons",
					m_CSLDMode,
					int(0)
				);

	registerProcessorParameter(	"TsldMode",
					"Event selection based on semi-leptonic decays of tau-leptons",
					m_TSLDMode,
					int(0)
				);

	registerProcessorParameter(	"sldMode",
					"Event selection based on semi-leptonic decays of All hadrons and tau-leptons",
					m_SLDMode,
					int(0)
				);


}


void SLDCorrection::init()
{

	streamlog_out(DEBUG) << "	init called  " << std::endl;
	m_Bfield = MarlinUtil::getBzAtOrigin();
	if ( !m_cheatPVAcharged ) m_cheatCharged4momentum = false;
	if ( !m_cheatPVAneutral ) m_cheatNeutral4momentum = false;
	printParameters();
	DDMarlinCED::init(this);

	if ( m_fillRootTree )
	{
		m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
		m_pTTree1 = new TTree("SLDCorrection", "SLDCorrection");
		m_pTTree1->SetDirectory(m_pTFile);
		m_pTTree1->Branch("event", &m_nEvt, "event/I");
		m_pTTree1->Branch("SLDFlavour", &m_SLDFlavour);
		m_pTTree1->Branch("SLDType", &m_SLDType);
		m_pTTree1->Branch("SLDLeptonID", &m_SLDLeptonID);
		m_pTTree1->Branch("leptonE_to_parentE", &m_leptonE_to_parentE);
		m_pTTree1->Branch("otherChargedE_to_parentE", &m_otherChargedE_to_parentE);
		m_pTTree1->Branch("allChargedE_to_parentE", &m_allChargedE_to_parentE);
		m_pTTree1->Branch("neutralE_to_parentE", &m_neutralE_to_parentE);
		m_pTTree1->Branch("neutrino_to_parentE", &m_neutrino_to_parentE);
		m_pTTree1->Branch("nSLDecayOfBHadron",&m_nSLDecayOfBHadron,"nSLDecayOfBHadron/I");
		m_pTTree1->Branch("nSLDecayOfCHadron",&m_nSLDecayOfCHadron,"nSLDecayOfCHadron/I");
		m_pTTree1->Branch("nSLDecayOfTauLepton",&m_nSLDecayOfTauLepton,"nSLDecayOfTauLepton/I");
		m_pTTree1->Branch("nSLDecayTotal",&m_nSLDecayTotal,"nSLDecayTotal/I");
		m_pTTree1->Branch("nSLDecayToElectron",&m_nSLDecayToElectron,"nSLDecayToElectron/I");
		m_pTTree1->Branch("nSLDecayToMuon",&m_nSLDecayToMuon,"nSLDecayToMuon/I");
		m_pTTree1->Branch("nSLDecayToTau",&m_nSLDecayToTau,"nSLDecayToTau/I");
		m_pTTree1->Branch("SLDStatus", &m_SLDStatus );
		m_pTTree1->Branch("PVAStatus", &m_PVAStatus );
		m_pTTree1->Branch("trueSolutionSign", &m_trueSolutionSign );
		m_pTTree1->Branch("recoSolutionSign", &m_recoSolutionSign );
		m_pTTree1->Branch("nTauNeutrino",&m_nTauNeutrino,"nTauNeutrino/I");
		m_pTTree1->Branch("nNeutrino",&m_nNeutrino,"nNeutrino/I");
		m_pTTree1->Branch("nSLD_chargedMCPwoTrack",&m_nSLD_chargedMCPwoTrack);
		m_pTTree1->Branch("GenStatParentHadron",&m_GenStatParentHadron);
		m_pTTree1->Branch("ChargeParentHadron",&m_ChargeParentHadron);
		m_pTTree1->Branch("foundRecoLepton",&m_foundRecoLepton);
		m_pTTree1->Branch("foundBuildUpVertex",&m_foundBuildUpVertex);
		m_pTTree1->Branch("foundRecoLeptonInBuildUpVertex",&m_foundRecoLeptonInBuildUpVertex);
		m_pTTree1->Branch("foundRecoLeptonInPrimaryVertex",&m_foundRecoLeptonInPrimaryVertex);
		m_pTTree1->Branch("lostChargedMCP_CosTheta",&m_lostChargedMCP_CosTheta);
		m_pTTree1->Branch("lostChargedMCP_Energy",&m_lostChargedMCP_Energy);
		m_pTTree1->Branch("lostChargedMCP_Pt",&m_lostChargedMCP_Pt);
		m_pTTree1->Branch("SLDecayXi", &m_SLDecayXi);
		m_pTTree1->Branch("SLDecayYi", &m_SLDecayYi);
		m_pTTree1->Branch("SLDecayZi", &m_SLDecayZi);
		m_pTTree1->Branch("SLDecayRi", &m_SLDecayRi);
		m_pTTree1->Branch("SLDecayXf", &m_SLDecayXf);
		m_pTTree1->Branch("SLDecayYf", &m_SLDecayYf);
		m_pTTree1->Branch("SLDecayZf", &m_SLDecayZf);
		m_pTTree1->Branch("SLDecayRf", &m_SLDecayRf);
		m_pTTree1->Branch("trueNuPx", &m_trueNuPx);
		m_pTTree1->Branch("trueNuPy", &m_trueNuPy);
		m_pTTree1->Branch("trueNuPz", &m_trueNuPz);
		m_pTTree1->Branch("trueNuE", &m_trueNuE);
		m_pTTree1->Branch("recoNuCloseInitialPx", &m_recoNuCloseInitialPx);
		m_pTTree1->Branch("recoNuCloseInitialPy", &m_recoNuCloseInitialPy);
		m_pTTree1->Branch("recoNuCloseInitialPz", &m_recoNuCloseInitialPz);
		m_pTTree1->Branch("recoNuCloseInitialE", &m_recoNuCloseInitialE);
		m_pTTree1->Branch("recoNuClosePx", &m_recoNuClosePx);
		m_pTTree1->Branch("recoNuClosePy", &m_recoNuClosePy);
		m_pTTree1->Branch("recoNuClosePz", &m_recoNuClosePz);
		m_pTTree1->Branch("recoNuCloseE", &m_recoNuCloseE);
		m_pTTree1->Branch("recoNuPosPx", &m_recoNuPosPx);
		m_pTTree1->Branch("recoNuPosPy", &m_recoNuPosPy);
		m_pTTree1->Branch("recoNuPosPz", &m_recoNuPosPz);
		m_pTTree1->Branch("recoNuPosE", &m_recoNuPosE);
		m_pTTree1->Branch("recoNuNegPx", &m_recoNuNegPx);
		m_pTTree1->Branch("recoNuNegPy", &m_recoNuNegPy);
		m_pTTree1->Branch("recoNuNegPz", &m_recoNuNegPz);
		m_pTTree1->Branch("recoNuNegE", &m_recoNuNegE);
		m_pTTree1->Branch("NuPxResidual", &m_NuPxResidual);
		m_pTTree1->Branch("NuPyResidual", &m_NuPyResidual);
		m_pTTree1->Branch("NuPzResidual", &m_NuPzResidual);
		m_pTTree1->Branch("NuEResidual", &m_NuEResidual);
		m_pTTree1->Branch("NuPxNormalizedResidual", &m_NuPxNormalizedResidual);
		m_pTTree1->Branch("NuPyNormalizedResidual", &m_NuPyNormalizedResidual);
		m_pTTree1->Branch("NuPzNormalizedResidual", &m_NuPzNormalizedResidual);
		m_pTTree1->Branch("NuENormalizedResidual", &m_NuENormalizedResidual);
		m_pTTree1->Branch("solutionSign", &m_solutionSign);
		m_pTTree1->Branch("true_E_vis", &m_true_E_vis);
		m_pTTree1->Branch("true_E_vis_prime", &m_true_E_vis_prime);
		m_pTTree1->Branch("true_P_vis_par", &m_true_P_vis_par);
		m_pTTree1->Branch("true_P_vis_par_prime", &m_true_P_vis_par_prime);
		m_pTTree1->Branch("true_P_vis_nor", &m_true_P_vis_nor);
		m_pTTree1->Branch("E_vis", &m_E_vis);
		m_pTTree1->Branch("E_vis_prime", &m_E_vis_prime);
		m_pTTree1->Branch("P_vis_par", &m_P_vis_par);
		m_pTTree1->Branch("P_vis_par_prime_squared", &m_P_vis_par_prime_squared);
		m_pTTree1->Branch("P_vis_par_prime", &m_P_vis_par_prime);
		m_pTTree1->Branch("P_vis_nor", &m_P_vis_nor);
		m_pTTree1->Branch("P_vis_nor_prime", &m_P_vis_nor_prime);
		m_pTTree1->Branch("flightDirectionStatus", &m_flightDirectionStatus);
		m_pTTree1->Branch("flightDirectionErrorCosAlpha", &m_FlightDirectionErrorCosAlpha);
		m_pTTree1->Branch("flightDirectionErrorSinAlpha", &m_FlightDirectionErrorSinAlpha);
		m_pTTree1->Branch("flightDirectionErrorAlpha", &m_FlightDirectionErrorAlpha);
		m_pTTree1->Branch("distRecoLeptonToDownStreamVertex", &m_distRecoLeptonToDownStreamVertex);
		m_pTTree1->Branch("dsVertexResidualX", &m_dsVertexResidualX);
		m_pTTree1->Branch("dsVertexResidualY", &m_dsVertexResidualY);
		m_pTTree1->Branch("dsVertexResidualZ", &m_dsVertexResidualZ);
		m_pTTree1->Branch("secVertexResidualX", &m_SecVertexResidualX);
		m_pTTree1->Branch("secVertexResidualY", &m_SecVertexResidualY);
		m_pTTree1->Branch("secVertexResidualZ", &m_SecVertexResidualZ);
		m_pTTree1->Branch("parentHadronMass", &m_parentHadronMass );
		m_pTTree1->Branch("parentHadronPDG", &m_parentHadronPDG );
		m_pTTree1->Branch("trueParentHadronFlightDistance", &m_trueParentHadronFlightDistance );
		m_pTTree1->Branch("recoParentHadronFlightDistance", &m_recoParentHadronFlightDistance );
		m_pTTree1->Branch("daughterHadronMass", &m_daughterHadronMass );
		m_pTTree1->Branch("daughterHadronPDG", &m_daughterHadronPDG );
		m_pTTree1->Branch("daughterHadronFlightDistance", &m_daughterHadronFlightDistance );
		m_pTTree1->Branch("jetEnergy", &m_jetEnergy );
		m_pTTree1->Branch("jetEnergyFractionCharged", &m_jetEnergyFractionCharged );
		m_pTTree1->Branch("jetEnergyFractionNeutralHadron", &m_jetEnergyFractionNeutralHadron );
		m_pTTree1->Branch("jetEnergyFractionPhoton", &m_jetEnergyFractionPhoton );
		m_pTTree1->Branch("jetEnergyFractionNeutrals", &m_jetEnergyFractionNeutrals );
		m_pTTree1->Branch("nTrueNeutralDecayProducts",&m_nTrueNeutralDecayProducts);
		m_pTTree1->Branch("nTrueAloneChargedDecayProducts",&m_nTrueAloneChargedDecayProducts);
		m_pTTree1->Branch("nTrueVertices",&m_nTrueVertices);
		m_pTTree1->Branch("nJetsVerticesDistributedIn",&m_nJetsVerticesDistributedIn);
		m_pTTree1->Branch("nRecoVerticesInJet",&m_nRecoVerticesInJet);
		m_pTTree1->Branch("nAloneChargedPFOs",&m_nAloneChargedPFOs);
		m_pTTree1->Branch("nAloneChargedPFOsFromSLD",&m_nAloneChargedPFOsFromSLD);
		m_pTTree1->Branch("distLeptonAlonePFOsNotFromSLD",&m_distLeptonAlonePFOsNotFromSLD);
		m_pTTree1->Branch("distLeptonAlonePFOsFromSLD",&m_distLeptonAlonePFOsFromSLD);
		m_pTTree1->Branch("weightPFOtoMCP_Lepton", &m_weightPFOtoMCP_Lepton );
		m_pTTree1->Branch("weightMCPtoPFO_Lepton", &m_weightMCPtoPFO_Lepton );
		m_pTTree1->Branch("weightPFOtoMCP_Neutral", &m_weightPFOtoMCP_Neutral );
		m_pTTree1->Branch("weightMCPtoPFO_Neutral", &m_weightMCPtoPFO_Neutral );
		m_pTTree1->Branch("weightPFOtoMCP_Charged", &m_weightPFOtoMCP_Charged );
		m_pTTree1->Branch("weightMCPtoPFO_Charged", &m_weightMCPtoPFO_Charged );

		m_pTTree2 = new TTree("FourMomentums", "FourMomentums");
		m_pTTree2->Branch("event", &m_nEvt, "event/I");
		m_pTTree2->Branch("SLDStatus", &m_SLDStatus );
		m_pTTree2->Branch("PVAStatus", &m_PVAStatus );
		m_pTTree2->Branch("trueSolutionSign", &m_trueSolutionSign );
		m_pTTree2->Branch("recoSolutionSign", &m_recoSolutionSign );
		m_pTTree2->Branch("SLDFlavour", &m_SLDFlavour);
		m_pTTree2->Branch("SLDType", &m_SLDType);
		m_pTTree2->Branch("SLDLeptonID", &m_SLDLeptonID);
		m_pTTree2->Branch("leptonE_to_parentE", &m_leptonE_to_parentE);
		m_pTTree2->Branch("otherChargedE_to_parentE", &m_otherChargedE_to_parentE);
		m_pTTree2->Branch("allChargedE_to_parentE", &m_allChargedE_to_parentE);
		m_pTTree2->Branch("neutralE_to_parentE", &m_neutralE_to_parentE);
		m_pTTree2->Branch("neutrino_to_parentE", &m_neutrino_to_parentE);
		m_pTTree2->Branch("visibleChargedInvMassCut", &m_visibleChargedInvMassCut );
		m_pTTree2->Branch("visibleNeutralInvMassCut", &m_visibleNeutralInvMassCut );
		m_pTTree2->Branch( "nChargedParticlesInPrimVertex" , &m_nChargedParticlesInPrimVertex );
		m_pTTree2->Branch( "nChargedParticlesNotInPrimVertex" , &m_nChargedParticlesNotInPrimVertex );
		m_pTTree2->Branch( "trueVisibleFourMomentumAtSLDVertex_Px" , &m_trueVisibleFourMomentumAtSLDVertex_Px );
		m_pTTree2->Branch( "trueVisibleFourMomentumAtSLDVertex_Py" , &m_trueVisibleFourMomentumAtSLDVertex_Py );
		m_pTTree2->Branch( "trueVisibleFourMomentumAtSLDVertex_Pz" , &m_trueVisibleFourMomentumAtSLDVertex_Pz );
		m_pTTree2->Branch( "trueVisibleFourMomentumAtSLDVertex_E" , &m_trueVisibleFourMomentumAtSLDVertex_E );
		m_pTTree2->Branch( "trueVisibleFourMomentumAtSLDVertex_M" , &m_trueVisibleFourMomentumAtSLDVertex_M );
		m_pTTree2->Branch( "truePVATrueFourMomentum_Px" , &m_truePVATrueFourMomentum_Px );
		m_pTTree2->Branch( "truePVATrueFourMomentum_Py" , &m_truePVATrueFourMomentum_Py );
		m_pTTree2->Branch( "truePVATrueFourMomentum_Pz" , &m_truePVATrueFourMomentum_Pz );
		m_pTTree2->Branch( "truePVATrueFourMomentum_E" , &m_truePVATrueFourMomentum_E );
		m_pTTree2->Branch( "truePVATrueFourMomentum_M" , &m_truePVATrueFourMomentum_M );
		m_pTTree2->Branch( "truePVARecoFourMomentum_Px" , &m_truePVARecoFourMomentum_Px );
		m_pTTree2->Branch( "truePVARecoFourMomentum_Py" , &m_truePVARecoFourMomentum_Py );
		m_pTTree2->Branch( "truePVARecoFourMomentum_Pz" , &m_truePVARecoFourMomentum_Pz );
		m_pTTree2->Branch( "truePVARecoFourMomentum_E" , &m_truePVARecoFourMomentum_E );
		m_pTTree2->Branch( "truePVARecoFourMomentum_M" , &m_truePVARecoFourMomentum_M );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_Px" , &m_recoPVARecoFourMomentum_Px );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_Py" , &m_recoPVARecoFourMomentum_Py );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_Pz" , &m_recoPVARecoFourMomentum_Pz );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_E" , &m_recoPVARecoFourMomentum_E );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_M" , &m_recoPVARecoFourMomentum_M );
		m_pTTree2->Branch( "usedVisibleFourMomentum_Px" , &m_usedVisibleFourMomentum_Px );
		m_pTTree2->Branch( "usedVisibleFourMomentum_Py" , &m_usedVisibleFourMomentum_Py );
		m_pTTree2->Branch( "usedVisibleFourMomentum_Pz" , &m_usedVisibleFourMomentum_Pz );
		m_pTTree2->Branch( "usedVisibleFourMomentum_E" , &m_usedVisibleFourMomentum_E );
		m_pTTree2->Branch( "usedVisibleFourMomentum_M" , &m_usedVisibleFourMomentum_M );
		m_pTTree2->Branch( "trueLeptonFourMomentum_Px" , &m_trueLeptonFourMomentum_Px );
		m_pTTree2->Branch( "trueLeptonFourMomentum_Py" , &m_trueLeptonFourMomentum_Py );
		m_pTTree2->Branch( "trueLeptonFourMomentum_Pz" , &m_trueLeptonFourMomentum_Pz );
		m_pTTree2->Branch( "trueLeptonFourMomentum_E" , &m_trueLeptonFourMomentum_E );
		m_pTTree2->Branch( "trueLeptonFourMomentum_M" , &m_trueLeptonFourMomentum_M );
		m_pTTree2->Branch( "recoLeptonFourMomentum_Px" , &m_recoLeptonFourMomentum_Px );
		m_pTTree2->Branch( "recoLeptonFourMomentum_Py" , &m_recoLeptonFourMomentum_Py );
		m_pTTree2->Branch( "recoLeptonFourMomentum_Pz" , &m_recoLeptonFourMomentum_Pz );
		m_pTTree2->Branch( "recoLeptonFourMomentum_E" , &m_recoLeptonFourMomentum_E );
		m_pTTree2->Branch( "recoLeptonFourMomentum_M" , &m_recoLeptonFourMomentum_M );
		m_pTTree2->Branch( "usedLeptonFourMomentum_Px" , &m_usedLeptonFourMomentum_Px );
		m_pTTree2->Branch( "usedLeptonFourMomentum_Py" , &m_usedLeptonFourMomentum_Py );
		m_pTTree2->Branch( "usedLeptonFourMomentum_Pz" , &m_usedLeptonFourMomentum_Pz );
		m_pTTree2->Branch( "usedLeptonFourMomentum_E" , &m_usedLeptonFourMomentum_E );
		m_pTTree2->Branch( "usedLeptonFourMomentum_M" , &m_usedLeptonFourMomentum_M );
		m_pTTree2->Branch( "truePVATrueChargedFourMomentum_Px" , &m_truePVATrueChargedFourMomentum_Px );
		m_pTTree2->Branch( "truePVATrueChargedFourMomentum_Py" , &m_truePVATrueChargedFourMomentum_Py );
		m_pTTree2->Branch( "truePVATrueChargedFourMomentum_Pz" , &m_truePVATrueChargedFourMomentum_Pz );
		m_pTTree2->Branch( "truePVATrueChargedFourMomentum_E" , &m_truePVATrueChargedFourMomentum_E );
		m_pTTree2->Branch( "truePVATrueChargedFourMomentum_M" , &m_truePVATrueChargedFourMomentum_M );
		m_pTTree2->Branch( "truePVARecoChargedFourMomentum_Px" , &m_truePVARecoChargedFourMomentum_Px );
		m_pTTree2->Branch( "truePVARecoChargedFourMomentum_Py" , &m_truePVARecoChargedFourMomentum_Py );
		m_pTTree2->Branch( "truePVARecoChargedFourMomentum_Pz" , &m_truePVARecoChargedFourMomentum_Pz );
		m_pTTree2->Branch( "truePVARecoChargedFourMomentum_E" , &m_truePVARecoChargedFourMomentum_E );
		m_pTTree2->Branch( "truePVARecoChargedFourMomentum_M" , &m_truePVARecoChargedFourMomentum_M );
		m_pTTree2->Branch( "recoSLDVertexChargedFourMomentum_Px" , &m_recoSLDVertexChargedFourMomentum_Px );
		m_pTTree2->Branch( "recoSLDVertexChargedFourMomentum_Py" , &m_recoSLDVertexChargedFourMomentum_Py );
		m_pTTree2->Branch( "recoSLDVertexChargedFourMomentum_Pz" , &m_recoSLDVertexChargedFourMomentum_Pz );
		m_pTTree2->Branch( "recoSLDVertexChargedFourMomentum_E" , &m_recoSLDVertexChargedFourMomentum_E );
		m_pTTree2->Branch( "recoSLDVertexChargedFourMomentum_M" , &m_recoSLDVertexChargedFourMomentum_M );
		m_pTTree2->Branch( "recoPVARecoChargedFourMomentum_Px" , &m_recoPVARecoChargedFourMomentum_Px );
		m_pTTree2->Branch( "recoPVARecoChargedFourMomentum_Py" , &m_recoPVARecoChargedFourMomentum_Py );
		m_pTTree2->Branch( "recoPVARecoChargedFourMomentum_Pz" , &m_recoPVARecoChargedFourMomentum_Pz );
		m_pTTree2->Branch( "recoPVARecoChargedFourMomentum_E" , &m_recoPVARecoChargedFourMomentum_E );
		m_pTTree2->Branch( "recoPVARecoChargedFourMomentum_M" , &m_recoPVARecoChargedFourMomentum_M );
		m_pTTree2->Branch( "usedChargedFourMomentum_Px" , &m_usedChargedFourMomentum_Px );
		m_pTTree2->Branch( "usedChargedFourMomentum_Py" , &m_usedChargedFourMomentum_Py );
		m_pTTree2->Branch( "usedChargedFourMomentum_Pz" , &m_usedChargedFourMomentum_Pz );
		m_pTTree2->Branch( "usedChargedFourMomentum_E" , &m_usedChargedFourMomentum_E );
		m_pTTree2->Branch( "usedChargedFourMomentum_M" , &m_usedChargedFourMomentum_M );
		m_pTTree2->Branch( "truePVATrueNeutralFourMomentum_Px" , &m_truePVATrueNeutralFourMomentum_Px );
		m_pTTree2->Branch( "truePVATrueNeutralFourMomentum_Py" , &m_truePVATrueNeutralFourMomentum_Py );
		m_pTTree2->Branch( "truePVATrueNeutralFourMomentum_Pz" , &m_truePVATrueNeutralFourMomentum_Pz );
		m_pTTree2->Branch( "truePVATrueNeutralFourMomentum_E" , &m_truePVATrueNeutralFourMomentum_E );
		m_pTTree2->Branch( "truePVATrueNeutralFourMomentum_M" , &m_truePVATrueNeutralFourMomentum_M );
		m_pTTree2->Branch( "truePVARecoNeutralFourMomentum_Px" , &m_truePVARecoNeutralFourMomentum_Px );
		m_pTTree2->Branch( "truePVARecoNeutralFourMomentum_Py" , &m_truePVARecoNeutralFourMomentum_Py );
		m_pTTree2->Branch( "truePVARecoNeutralFourMomentum_Pz" , &m_truePVARecoNeutralFourMomentum_Pz );
		m_pTTree2->Branch( "truePVARecoNeutralFourMomentum_E" , &m_truePVARecoNeutralFourMomentum_E );
		m_pTTree2->Branch( "truePVARecoNeutralFourMomentum_M" , &m_truePVARecoNeutralFourMomentum_M );
		m_pTTree2->Branch( "recoPVARecoNeutralFourMomentum_Px" , &m_recoPVARecoNeutralFourMomentum_Px );
		m_pTTree2->Branch( "recoPVARecoNeutralFourMomentum_Py" , &m_recoPVARecoNeutralFourMomentum_Py );
		m_pTTree2->Branch( "recoPVARecoNeutralFourMomentum_Pz" , &m_recoPVARecoNeutralFourMomentum_Pz );
		m_pTTree2->Branch( "recoPVARecoNeutralFourMomentum_E" , &m_recoPVARecoNeutralFourMomentum_E );
		m_pTTree2->Branch( "recoPVARecoNeutralFourMomentum_M" , &m_recoPVARecoNeutralFourMomentum_M );
		m_pTTree2->Branch( "usedNeutralFourMomentum_Px" , &m_usedNeutralFourMomentum_Px );
		m_pTTree2->Branch( "usedNeutralFourMomentum_Py" , &m_usedNeutralFourMomentum_Py );
		m_pTTree2->Branch( "usedNeutralFourMomentum_Pz" , &m_usedNeutralFourMomentum_Pz );
		m_pTTree2->Branch( "usedNeutralFourMomentum_E" , &m_usedNeutralFourMomentum_E );
		m_pTTree2->Branch( "usedNeutralFourMomentum_M" , &m_usedNeutralFourMomentum_M );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Px" , &m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Px );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Py" , &m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Py );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Pz" , &m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Pz );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_E" , &m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_E );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_M" , &m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_M );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Px" , &m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Px );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Py" , &m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Py );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Pz" , &m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Pz );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_E" , &m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_E );
		m_pTTree2->Branch( "recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_M" , &m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_M );
		m_pTTree2->Branch( "trueNeutrinoFourMomentum_Px" , &m_trueNeutrinoFourMomentum_Px );
		m_pTTree2->Branch( "trueNeutrinoFourMomentum_Py" , &m_trueNeutrinoFourMomentum_Py );
		m_pTTree2->Branch( "trueNeutrinoFourMomentum_Pz" , &m_trueNeutrinoFourMomentum_Pz );
		m_pTTree2->Branch( "trueNeutrinoFourMomentum_E" , &m_trueNeutrinoFourMomentum_E );
		m_pTTree2->Branch( "trueNeutrinoFourMomentum_M" , &m_trueNeutrinoFourMomentum_M );
		m_pTTree2->Branch( "recoNeutrinoFourMomentumClose_Px" , &m_recoNeutrinoFourMomentumClose_Px );
		m_pTTree2->Branch( "recoNeutrinoFourMomentumClose_Py" , &m_recoNeutrinoFourMomentumClose_Py );
		m_pTTree2->Branch( "recoNeutrinoFourMomentumClose_Pz" , &m_recoNeutrinoFourMomentumClose_Pz );
		m_pTTree2->Branch( "recoNeutrinoFourMomentumClose_E" , &m_recoNeutrinoFourMomentumClose_E );
		m_pTTree2->Branch( "recoNeutrinoFourMomentumClose_M" , &m_recoNeutrinoFourMomentumClose_M );
		m_pTTree2->Branch( "trueHadronFourMomentum_Px" , &m_trueHadronFourMomentum_Px );
		m_pTTree2->Branch( "trueHadronFourMomentum_Py" , &m_trueHadronFourMomentum_Py );
		m_pTTree2->Branch( "trueHadronFourMomentum_Pz" , &m_trueHadronFourMomentum_Pz );
		m_pTTree2->Branch( "trueHadronFourMomentum_E" , &m_trueHadronFourMomentum_E );
		m_pTTree2->Branch( "trueHadronFourMomentum_M" , &m_trueHadronFourMomentum_M );
		m_pTTree2->Branch( "recoHadronFourMomentum_Px" , &m_recoHadronFourMomentum_Px );
		m_pTTree2->Branch( "recoHadronFourMomentum_Py" , &m_recoHadronFourMomentum_Py );
		m_pTTree2->Branch( "recoHadronFourMomentum_Pz" , &m_recoHadronFourMomentum_Pz );
		m_pTTree2->Branch( "recoHadronFourMomentum_E" , &m_recoHadronFourMomentum_E );
		m_pTTree2->Branch( "recoHadronFourMomentum_M" , &m_recoHadronFourMomentum_M );
		m_pTTree2->Branch( "trueFlightDirection_X" , &m_trueFlightDirection_X );
		m_pTTree2->Branch( "trueFlightDirection_Y" , &m_trueFlightDirection_Y );
		m_pTTree2->Branch( "trueFlightDirection_Z" , &m_trueFlightDirection_Z );
		m_pTTree2->Branch( "recoFlightDirection_X" , &m_recoFlightDirection_X );
		m_pTTree2->Branch( "recoFlightDirection_Y" , &m_recoFlightDirection_Y );
		m_pTTree2->Branch( "recoFlightDirection_Z" , &m_recoFlightDirection_Z );
		m_pTTree2->Branch( "usedFlightDirection_X" , &m_usedFlightDirection_X );
		m_pTTree2->Branch( "usedFlightDirection_Y" , &m_usedFlightDirection_Y );
		m_pTTree2->Branch( "usedFlightDirection_Z" , &m_usedFlightDirection_Z );
		m_pTTree2->Branch( "neutralPVAAlpha", &m_NeutralPVAAlpha );
		m_pTTree2->Branch( "neutralPVASinAlpha", &m_NeutralPVASinAlpha );
		m_pTTree2->Branch( "neutralPVACosAlpha", &m_NeutralPVACosAlpha );
		m_pTTree2->Branch( "chargedPVAAlpha", &m_ChargedPVAAlpha );
		m_pTTree2->Branch( "chargedPVASinAlpha", &m_ChargedPVASinAlpha );
		m_pTTree2->Branch( "chargedPVACosAlpha", &m_ChargedPVACosAlpha );
		m_pTTree2->Branch( "PVAAlpha", &m_PVAAlpha );
		m_pTTree2->Branch( "PVASinAlpha", &m_PVASinAlpha );
		m_pTTree2->Branch( "PVACosAlpha", &m_PVACosAlpha );
		m_pTTree2->Branch( "Alpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral" , &m_Alpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral );
		m_pTTree2->Branch( "SinAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral" , &m_SinAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral );
		m_pTTree2->Branch( "CosAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral" , &m_CosAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral );
		m_pTTree2->Branch( "Alpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged" , &m_Alpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged );
		m_pTTree2->Branch( "SinAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged" , &m_SinAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged );
		m_pTTree2->Branch( "CosAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged" , &m_CosAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged );
		m_pTTree2->Branch( "visiblePxBeforePVA", &m_visiblePxBeforePVA );
		m_pTTree2->Branch( "visiblePyBeforePVA", &m_visiblePyBeforePVA );
		m_pTTree2->Branch( "visiblePzBeforePVA", &m_visiblePzBeforePVA );
		m_pTTree2->Branch( "visibleEnergyBeforePVA", &m_visibleEnergyBeforePVA );
		m_pTTree2->Branch( "visibleMassBeforePVA", &m_visibleMassBeforePVA );
		m_pTTree2->Branch( "visiblePxAfterVertexPVA", &m_visiblePxAfterVertexPVA );
		m_pTTree2->Branch( "visiblePyAfterVertexPVA", &m_visiblePyAfterVertexPVA );
		m_pTTree2->Branch( "visiblePzAfterVertexPVA", &m_visiblePzAfterVertexPVA );
		m_pTTree2->Branch( "visibleEnergyAfterVertexPVA", &m_visibleEnergyAfterVertexPVA );
		m_pTTree2->Branch( "visibleMassAfterVertexPVA", &m_visibleMassAfterVertexPVA );
		m_pTTree2->Branch( "visiblePxAfterChargedPVA", &m_visiblePxAfterChargedPVA );
		m_pTTree2->Branch( "visiblePyAfterChargedPVA", &m_visiblePyAfterChargedPVA );
		m_pTTree2->Branch( "visiblePzAfterChargedPVA", &m_visiblePzAfterChargedPVA );
		m_pTTree2->Branch( "visibleEnergyAfterChargedPVA", &m_visibleEnergyAfterChargedPVA );
		m_pTTree2->Branch( "visibleMassAfterChargedPVA", &m_visibleMassAfterChargedPVA );
		m_pTTree2->Branch( "visiblePxAfterNeutralPVA", &m_visiblePxAfterNeutralPVA );
		m_pTTree2->Branch( "visiblePyAfterNeutralPVA", &m_visiblePyAfterNeutralPVA );
		m_pTTree2->Branch( "visiblePzAfterNeutralPVA", &m_visiblePzAfterNeutralPVA );
		m_pTTree2->Branch( "visibleEnergyAfterNeutralPVA", &m_visibleEnergyAfterNeutralPVA );
		m_pTTree2->Branch( "visibleMassAfterNeutralPVA", &m_visibleMassAfterNeutralPVA );
		m_pTTree2->Branch( "NuPxResidual", &m_NuPxResidual);
		m_pTTree2->Branch( "NuPyResidual", &m_NuPyResidual);
		m_pTTree2->Branch( "NuPzResidual", &m_NuPzResidual);
		m_pTTree2->Branch( "NuEResidual", &m_NuEResidual);
		m_pTTree2->Branch( "NuPxNormalizedResidual", &m_NuPxNormalizedResidual);
		m_pTTree2->Branch( "NuPyNormalizedResidual", &m_NuPyNormalizedResidual);
		m_pTTree2->Branch( "NuPzNormalizedResidual", &m_NuPzNormalizedResidual);
		m_pTTree2->Branch( "NuENormalizedResidual", &m_NuENormalizedResidual);

		m_pTTree3 = new TTree("Errors", "Errors");
		m_pTTree3->Branch("event", &m_nEvt, "event/I");
		m_pTTree3->Branch("SLDStatus", &m_SLDStatus );
		m_pTTree3->Branch("PVAStatus", &m_PVAStatus );
		m_pTTree3->Branch("trueSolutionSign", &m_trueSolutionSign );
		m_pTTree3->Branch("recoSolutionSign", &m_recoSolutionSign );
		m_pTTree3->Branch("SLDFlavour", &m_SLDFlavour);
		m_pTTree3->Branch("SLDType", &m_SLDType);
		m_pTTree3->Branch("SLDLeptonID", &m_SLDLeptonID);
		m_pTTree3->Branch( "trueVisibleFourMomentumPx", &m_trueVisibleFourMomentumPx );
		m_pTTree3->Branch( "trueVisibleFourMomentumPy", &m_trueVisibleFourMomentumPy );
		m_pTTree3->Branch( "trueVisibleFourMomentumPz", &m_trueVisibleFourMomentumPz );
		m_pTTree3->Branch( "trueVisibleFourMomentumE", &m_trueVisibleFourMomentumE );
		m_pTTree3->Branch( "recoVisibleFourMomentumPx", &m_recoVisibleFourMomentumPx );
		m_pTTree3->Branch( "recoVisibleFourMomentumPy", &m_recoVisibleFourMomentumPy );
		m_pTTree3->Branch( "recoVisibleFourMomentumPz", &m_recoVisibleFourMomentumPz );
		m_pTTree3->Branch( "recoVisibleFourMomentumE", &m_recoVisibleFourMomentumE );
		m_pTTree3->Branch( "residualVisibleFourMomentumPx", &m_residualVisibleFourMomentumPx );
		m_pTTree3->Branch( "residualVisibleFourMomentumPy", &m_residualVisibleFourMomentumPy );
		m_pTTree3->Branch( "residualVisibleFourMomentumPz", &m_residualVisibleFourMomentumPz );
		m_pTTree3->Branch( "residualVisibleFourMomentumE", &m_residualVisibleFourMomentumE );
		m_pTTree3->Branch( "sigmaPxPx_Det", &m_sigmaPxPx_Det );
		m_pTTree3->Branch( "sigmaPxPy_Det", &m_sigmaPxPy_Det );
		m_pTTree3->Branch( "sigmaPyPy_Det", &m_sigmaPyPy_Det );
		m_pTTree3->Branch( "sigmaPxPz_Det", &m_sigmaPxPz_Det );
		m_pTTree3->Branch( "sigmaPyPz_Det", &m_sigmaPyPz_Det );
		m_pTTree3->Branch( "sigmaPzPz_Det", &m_sigmaPzPz_Det );
		m_pTTree3->Branch( "sigmaPxE_Det", &m_sigmaPxE_Det );
		m_pTTree3->Branch( "sigmaPyE_Det", &m_sigmaPyE_Det );
		m_pTTree3->Branch( "sigmaPzE_Det", &m_sigmaPzE_Det );
		m_pTTree3->Branch( "sigmaEE_Det", &m_sigmaEE_Det );
		m_pTTree3->Branch( "normalizedResidualVisibleFourMomentumPx", &m_normalizedResidualVisibleFourMomentumPx );
		m_pTTree3->Branch( "normalizedResidualVisibleFourMomentumPy", &m_normalizedResidualVisibleFourMomentumPy );
		m_pTTree3->Branch( "normalizedResidualVisibleFourMomentumPz", &m_normalizedResidualVisibleFourMomentumPz );
		m_pTTree3->Branch( "normalizedResidualVisibleFourMomentumE", &m_normalizedResidualVisibleFourMomentumE );
		m_pTTree3->Branch( "trueFlightDirectionUx", &m_trueFlightDirectionUx );
		m_pTTree3->Branch( "trueFlightDirectionUy", &m_trueFlightDirectionUy );
		m_pTTree3->Branch( "trueFlightDirectionUz", &m_trueFlightDirectionUz );
		m_pTTree3->Branch( "recoFlightDirectionUx", &m_recoFlightDirectionUx );
		m_pTTree3->Branch( "recoFlightDirectionUy", &m_recoFlightDirectionUy );
		m_pTTree3->Branch( "recoFlightDirectionUz", &m_recoFlightDirectionUz );
		m_pTTree3->Branch( "residualFlightDirectionUx", &m_residualFlightDirectionUx );
		m_pTTree3->Branch( "residualFlightDirectionUy", &m_residualFlightDirectionUy );
		m_pTTree3->Branch( "residualFlightDirectionUz", &m_residualFlightDirectionUz );
		m_pTTree3->Branch( "sigmaUxUx", &m_sigmaUxUx );
		m_pTTree3->Branch( "sigmaUxUy", &m_sigmaUxUy );
		m_pTTree3->Branch( "sigmaUyUy", &m_sigmaUyUy );
		m_pTTree3->Branch( "sigmaUxUz", &m_sigmaUxUz );
		m_pTTree3->Branch( "sigmaUyUz", &m_sigmaUyUz );
		m_pTTree3->Branch( "sigmaUzUz", &m_sigmaUzUz );
		m_pTTree3->Branch( "normalizedResidualFlightDirectionUx", &m_normalizedResidualFlightDirectionUx );
		m_pTTree3->Branch( "normalizedResidualFlightDirectionUy", &m_normalizedResidualFlightDirectionUy );
		m_pTTree3->Branch( "normalizedResidualFlightDirectionUz", &m_normalizedResidualFlightDirectionUz );
		m_pTTree3->Branch( "trueNeutrinoFourMomentumPx", &m_trueNeutrinoFourMomentumPx );
		m_pTTree3->Branch( "trueNeutrinoFourMomentumPy", &m_trueNeutrinoFourMomentumPy );
		m_pTTree3->Branch( "trueNeutrinoFourMomentumPz", &m_trueNeutrinoFourMomentumPz );
		m_pTTree3->Branch( "trueNeutrinoFourMomentumE", &m_trueNeutrinoFourMomentumE );
		m_pTTree3->Branch( "recoNeutrinoFourMomentumClosePx", &m_recoNeutrinoFourMomentumClosePx );
		m_pTTree3->Branch( "recoNeutrinoFourMomentumClosePy", &m_recoNeutrinoFourMomentumClosePy );
		m_pTTree3->Branch( "recoNeutrinoFourMomentumClosePz", &m_recoNeutrinoFourMomentumClosePz );
		m_pTTree3->Branch( "recoNeutrinoFourMomentumCloseE", &m_recoNeutrinoFourMomentumCloseE );
		m_pTTree3->Branch( "sigmaNeutrinoPxPx", &m_sigmaNeutrinoPxPx );
		m_pTTree3->Branch( "sigmaNeutrinoPxPy", &m_sigmaNeutrinoPxPy );
		m_pTTree3->Branch( "sigmaNeutrinoPyPy", &m_sigmaNeutrinoPyPy );
		m_pTTree3->Branch( "sigmaNeutrinoPxPz", &m_sigmaNeutrinoPxPz );
		m_pTTree3->Branch( "sigmaNeutrinoPyPz", &m_sigmaNeutrinoPyPz );
		m_pTTree3->Branch( "sigmaNeutrinoPzPz", &m_sigmaNeutrinoPzPz );
		m_pTTree3->Branch( "sigmaNeutrinoPxE", &m_sigmaNeutrinoPxE );
		m_pTTree3->Branch( "sigmaNeutrinoPyE", &m_sigmaNeutrinoPyE );
		m_pTTree3->Branch( "sigmaNeutrinoPzE", &m_sigmaNeutrinoPzE );
		m_pTTree3->Branch( "sigmaNeutrinoEE", &m_sigmaNeutrinoEE );
		m_pTTree3->Branch( "residualNeutrinoFourMomentumPx", &m_residualNeutrinoFourMomentumPx );
		m_pTTree3->Branch( "residualNeutrinoFourMomentumPy", &m_residualNeutrinoFourMomentumPy );
		m_pTTree3->Branch( "residualNeutrinoFourMomentumPz", &m_residualNeutrinoFourMomentumPz );
		m_pTTree3->Branch( "residualNeutrinoFourMomentumE", &m_residualNeutrinoFourMomentumE );
		m_pTTree3->Branch( "normalizedResidualNeutrinoFourMomentumPx", &m_normalizedResidualNeutrinoFourMomentumPx );
		m_pTTree3->Branch( "normalizedResidualNeutrinoFourMomentumPy", &m_normalizedResidualNeutrinoFourMomentumPy );
		m_pTTree3->Branch( "normalizedResidualNeutrinoFourMomentumPz", &m_normalizedResidualNeutrinoFourMomentumPz );
		m_pTTree3->Branch( "normalizedResidualNeutrinoFourMomentumE", &m_normalizedResidualNeutrinoFourMomentumE );
		m_pTTree3->Branch( "recoNeutrinoDirectionError" , &m_recoNeutrinoDirectionError );
		m_pTTree3->Branch( "PCAatLeptonX" , &m_PCAatLeptonX );
		m_pTTree3->Branch( "PCAatLeptonY" , &m_PCAatLeptonY );
		m_pTTree3->Branch( "PCAatLeptonZ" , &m_PCAatLeptonZ );
		m_pTTree3->Branch( "PCAatOtherParticleX" , &m_PCAatOtherParticleX );
		m_pTTree3->Branch( "PCAatOtherParticleY" , &m_PCAatOtherParticleY );
		m_pTTree3->Branch( "PCAatOtherParticleZ" , &m_PCAatOtherParticleZ );
		m_pTTree3->Branch( "JetAxisX" , &m_JetAxisX );
		m_pTTree3->Branch( "JetAxisY" , &m_JetAxisY );
		m_pTTree3->Branch( "JetAxisZ" , &m_JetAxisZ	);
		m_pTTree3->Branch( "PrimaryVertexX" , &m_PrimaryVertexX );
		m_pTTree3->Branch( "PrimaryVertexY" , &m_PrimaryVertexY );
		m_pTTree3->Branch( "PrimaryVertexZ" , &m_PrimaryVertexZ );
		m_pTTree3->Branch( "AngleTrueFlightDirectionJet" , &m_AngleTrueFlightDirectionJet );
		m_pTTree3->Branch( "AngleRecoFlightDirectionJet" , &m_AngleRecoFlightDirectionJet );
		m_pTTree3->Branch( "AngleLeptonJet" , &m_AngleLeptonJet );
		m_pTTree3->Branch( "AngleDSVertexJet" , &m_AngleDSVertexJet );
		m_pTTree3->Branch( "AngleRecoNeutrinoJet" , &m_AngleRecoNeutrinoJet );
		m_pTTree3->Branch( "AngleTrueNeutrinoJet" , &m_AngleTrueNeutrinoJet );
		m_pTTree3->Branch( "LeptonDistanceFromPV" , &m_LeptonDistanceFromPV );
		m_pTTree3->Branch( "DSVDistanceFromPV" , &m_DSVDistanceFromPV );
		m_pTTree3->Branch( "flightDirectionErrorCosAlpha", &m_FlightDirectionErrorCosAlpha);
		m_pTTree3->Branch( "flightDirectionErrorSinAlpha", &m_FlightDirectionErrorSinAlpha);
		m_pTTree3->Branch( "flightDirectionErrorAlpha", &m_FlightDirectionErrorAlpha);
		m_pTTree3->Branch( "distRecoLeptonToDownStreamVertex", &m_distRecoLeptonToDownStreamVertex);
		m_pTTree3->Branch( "parentHadronMass", &m_parentHadronMass );
		m_pTTree3->Branch( "parentHadronPDG", &m_parentHadronPDG );
		m_pTTree3->Branch( "trueParentHadronFlightDistance", &m_trueParentHadronFlightDistance );
		m_pTTree3->Branch( "recoParentHadronFlightDistance", &m_recoParentHadronFlightDistance );
		m_pTTree3->Branch( "daughterHadronMass", &m_daughterHadronMass );
		m_pTTree3->Branch( "daughterHadronPDG", &m_daughterHadronPDG );
		m_pTTree3->Branch( "daughterHadronFlightDistance", &m_daughterHadronFlightDistance );
		m_pTTree3->Branch( "Lepton3DImpactParameter" , &m_Lepton3DImpactParameter );
		m_pTTree3->Branch( "OtherParticle3DImpactParameter" , &m_OtherParticle3DImpactParameter );
		h_SLDStatus = new TH1I( "SLDStatus" , ";" , 7 , 0 , 7 ); n_SLDStatus = 0;
		h_SLDStatus->GetXaxis()->SetBinLabel(1,"No #font[32]{l}^{REC}");
		h_SLDStatus->GetXaxis()->SetBinLabel(2,"#font[32]{l}#notin^{}jet");
		h_SLDStatus->GetXaxis()->SetBinLabel(3,"#font[32]{l}#in^{}Vtx^{Prim.}");
		h_SLDStatus->GetXaxis()->SetBinLabel(4,"#font[32]{l}#in^{}2^{nd}Vtx");
		h_SLDStatus->GetXaxis()->SetBinLabel(5,"#font[32]{l}+3^{rd}Vtx");
		h_SLDStatus->GetXaxis()->SetBinLabel(6,"#font[32]{l}+trk^{alone}");
		h_SLDStatus->GetXaxis()->SetBinLabel(7,"other");

		h_BHadronType = new TH1F( "BHadronType" , ";" , 8 , -0.5 , 7.5 );
		h_BHadronType->GetXaxis()->SetBinLabel(1,"B^{0}");//PDG = 511
		h_BHadronType->GetXaxis()->SetBinLabel(2,"B^{#pm}");//PDG = 521
		h_BHadronType->GetXaxis()->SetBinLabel(3,"B^{0}_{s}");//PDG = 531
		h_BHadronType->GetXaxis()->SetBinLabel(4,"B^{#pm}_{c}");//PDG = 541
		h_BHadronType->GetXaxis()->SetBinLabel(5,"#Lambda^{0}_{b}");//PDG = 5122
		h_BHadronType->GetXaxis()->SetBinLabel(6,"#Xi^{#pm}_{b}");//PDG = 5132
		h_BHadronType->GetXaxis()->SetBinLabel(7,"#Xi^{0}_{b}");//PDG = 5232
		h_BHadronType->GetXaxis()->SetBinLabel(8,"#Omega^{-}_{b}");//PDG = 5332

		h_CHadronType = new TH1F( "CHadronType" , ";" , 16 , -0.5 , 15.5 );
		h_CHadronType->GetXaxis()->SetBinLabel(1,"D^{#pm}");//PDG = 411
		h_CHadronType->GetXaxis()->SetBinLabel(2,"D^{*}(2010)^{#pm}");//PDG = 413
		h_CHadronType->GetXaxis()->SetBinLabel(3,"D^{*}_{2}(2460)^{#pm}");//PDG = 415
		h_CHadronType->GetXaxis()->SetBinLabel(4,"D^{0}");//PDG = 421
		h_CHadronType->GetXaxis()->SetBinLabel(5,"D^{*}(2007)^{0}");//PDG = 423
		h_CHadronType->GetXaxis()->SetBinLabel(6,"D^{*}_{2}(2460)^{0}");//PDG = 425
		h_CHadronType->GetXaxis()->SetBinLabel(7,"D^{#pm}_{s}");//PDG = 431
		h_CHadronType->GetXaxis()->SetBinLabel(8,"D^{*#pm}_{s}");//PDG = 433
		h_CHadronType->GetXaxis()->SetBinLabel(9,"D^{*}_{s2}");//PDG = 435
		h_CHadronType->GetXaxis()->SetBinLabel(10,"#eta_{c}");//PDG = 441
		h_CHadronType->GetXaxis()->SetBinLabel(11,"J/#psi");//PDG = 443
		h_CHadronType->GetXaxis()->SetBinLabel(12,"#Lambda^{#pm}_{c}");//PDG = 4122
		h_CHadronType->GetXaxis()->SetBinLabel(13,"#Xi^{0}_{c}");//PDG = 4132
		h_CHadronType->GetXaxis()->SetBinLabel(14,"#Xi^{#pm}_{c}");//PDG = 4232
		h_CHadronType->GetXaxis()->SetBinLabel(15,"#Omega^{0}_{c}");//PDG = 4332
		h_CHadronType->GetXaxis()->SetBinLabel(16,"#Omega^{*0}_{c}");//PDG = 4334

		h_NuPxResidual = new TH1F( "PxResidual" , "; _{}p_{x,#nu}^{REC} - p_{x,#nu}^{MC} [GeV]; Normalized #SLDs / 0.1 GeV" , 200 , -10.0 , 10.0 ); n_NuPxResidual = 0;
		h_NuPyResidual = new TH1F( "PyResidual" , "; _{}p_{y,#nu}^{REC} - p_{y,#nu}^{MC} [GeV]; Normalized #SLDs / 0.1 GeV" , 200 , -10.0 , 10.0 ); n_NuPyResidual = 0;
		h_NuPzResidual = new TH1F( "PzResidual" , "; _{}p_{z,#nu}^{REC} - p_{z,#nu}^{MC}  [GeV]; Normalized #SLDs / 0.1 GeV" , 200 , -10.0 , 10.0 ); n_NuPzResidual = 0;
		h_NuEResidual = new TH1F( "EResidual" , "; _{}E_{#nu}^{REC} - E_{#nu}^{MC} [GeV]; Normalized #SLDs / 0.1 GeV" , 200 , -10.0 , 10.0 ); n_NuEResidual = 0;
		h_NuPxNormalizedResidual = new TH1F( "PxNormalizedResidual" , "; ( _{}p_{x,#nu}^{REC} - p_{x,#nu}^{MC} ) / #sigma_{p_{x,#nu}}; Normalized #SLDs / 0.1" , 200 , -10.0 , 10.0 ); n_NuPxNormalizedResidual = 0;
		h_NuPyNormalizedResidual = new TH1F( "PyNormalizedResidual" , "; ( _{}p_{y,#nu}^{REC} - p_{y,#nu}^{MC} ) / #sigma_{p_{y,#nu}}; Normalized #SLDs / 0.1" , 200 , -10.0 , 10.0 ); n_NuPyNormalizedResidual = 0;
		h_NuPzNormalizedResidual = new TH1F( "PzNormalizedResidual" , "; ( _{}p_{z,#nu}^{REC} - p_{z,#nu}^{MC} ) / #sigma_{p_{z,#nu}}; Normalized #SLDs / 0.1" , 200 , -10.0 , 10.0 ); n_NuPzNormalizedResidual = 0;
		h_NuENormalizedResidual = new TH1F( "ENormalizedResidual" , "; ( _{}E_{#nu}^{REC} - E_{#nu}^{MC} ) / #sigma_{E_{#nu}}; Normalized #SLDs / 0.1" , 200 , -10.0 , 10.0 ); n_NuENormalizedResidual = 0;
		h_recoNuPx_mcNuPx = new TH2F( "p_{x}^{#nu}" , "; _{}p_{x,#nu}^{MC} [GeV] ;  _{}p_{x,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
		h_recoNuPy_mcNuPy = new TH2F( "p_{y}^{#nu}" , "; _{}p_{y,#nu}^{MC} [GeV] ;  _{}p_{y,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
		h_recoNuPz_mcNuPz = new TH2F( "p_{z}^{#nu}" , "; _{}p_{z,#nu}^{MC} [GeV] ;  _{}p_{z,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
		h_recoNuE_mcNuE = new TH2F( "E^{#nu}" , "; _{}E_{#nu}^{MC} [GeV] ;  _{}E_{#nu}^{REC} [GeV]" , 100 , 0.0 , 100.0 , 100 , 0.0 , 100.0 );
		h_parentPx_daughtersPx = new TH2F( "p_{x} conservation" , "; _{}Px_{parentHadron}^{MC} [GeV] ;  _{}Px_{daughters}^{MC} [GeV]" , 200 , -100.0 , 100.0 , 200 , -100.0 , 100.0 );
		h_parentPy_daughtersPy = new TH2F( "p_{y} conservation" , "; _{}Py_{parentHadron}^{MC} [GeV] ;  _{}Py_{daughters}^{MC} [GeV]" , 200 , -100.0 , 100.0 , 200 , -100.0 , 100.0 );
		h_parentPz_daughtersPz = new TH2F( "p_{z} conservation" , "; _{}Pz_{parentHadron}^{MC} [GeV] ;  _{}Pz_{daughters}^{MC} [GeV]" , 200 , -100.0 , 100.0 , 200 , -100.0 , 100.0 );
		h_parentE_daughtersE = new TH2F( "E conservation" , "; _{}E_{parentHadron}^{MC} [GeV] ;  _{}E_{daughters}^{MC} [GeV]" , 100 , 0.0 , 100.0 , 100 , 0.0 , 100.0 );
		h_recoPFOLinkedToElectron_Type = new TH1I( "PFOTypeofTrueElectron" , "; PFO Type" , 8 , 0 , 8 );
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(1,"e^{#pm}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(3,"#pi^{#pm}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(4,"#gamma");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(5,"K^{0}_{S}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(6,"#Lambda");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(7,"Other");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(8,"Not Found");
		h_recoPFOLinkedToMuon_Type = new TH1I( "PFOTypeofTrueMuon" , "; PFO Type" , 8 , 0 , 8 );
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(1,"e^{#pm}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(3,"#pi^{#pm}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(4,"#gamma");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(5,"K^{0}_{S}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(6,"#Lambda");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(7,"Other");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(8,"Not Found");
		h_SLDecayFlavour = new TH1I( "SLDecayFlavour" , "; SLDecay Type" , 3 , 0 , 3 );
		h_SLDecayFlavour->GetXaxis()->SetBinLabel(1,"SLD of B-Hadron");
		h_SLDecayFlavour->GetXaxis()->SetBinLabel(2,"SLD of C-Hadron");
		h_SLDecayFlavour->GetXaxis()->SetBinLabel(3,"LD of Tau-Lepton");
		h_SLDecayModeB = new TH1I( "SLDMode B-Hadron" , "; SLDecay Mode" , 3 , 0 , 3 );
		h_SLDecayModeB->GetXaxis()->SetBinLabel(1,"X#rightarrow e#nu_{e}Y");
		h_SLDecayModeB->GetXaxis()->SetBinLabel(2,"X#rightarrow #mu#nu_{#mu}Y");
		h_SLDecayModeB->GetXaxis()->SetBinLabel(3,"X#rightarrow #tau#nu_{#tau}Y");
		h_SLDecayModeC = new TH1I( "SLDMode C-Hadron" , "; SLDecay Mode" , 3 , 0 , 3 );
		h_SLDecayModeC->GetXaxis()->SetBinLabel(1,"X#rightarrow e#nu_{e}Y");
		h_SLDecayModeC->GetXaxis()->SetBinLabel(2,"X#rightarrow #mu#nu_{#mu}Y");
		h_SLDecayModeC->GetXaxis()->SetBinLabel(3,"X#rightarrow #tau#nu_{#tau}Y");
		h_SLDecayOrder = new TH1I( "SLDecayOrder" , "; SLDecay Type" , 2 , 0 , 2 );
		h_SLDecayOrder->GetXaxis()->SetBinLabel(1,"with UpStream/DownStream SLD");
		h_SLDecayOrder->GetXaxis()->SetBinLabel(2,"without UpStream/DownStream SLD");
		h_foundVertex = new TH2I( "Viertices" , "; primary vertex ; secondary vertex" , 2 , 0 , 2 , 2 , 0 , 2 );
		h_foundVertex->GetXaxis()->SetBinLabel(1,"vertex not found");
		h_foundVertex->GetXaxis()->SetBinLabel(2,"vertex found");
		h_foundVertex->GetYaxis()->SetBinLabel(1,"vertex not found");
		h_foundVertex->GetYaxis()->SetBinLabel(2,"vertex found");
		h_secondaryVertex = new TH1F( "secondary_vertices" , ";" , 11 , 0 , 11 ); n_secondaryVertex = 0;
		h_secondaryVertex->GetXaxis()->SetBinLabel(1,"reco lep not found");
		h_secondaryVertex->GetXaxis()->SetBinLabel(2,"reco lep not in jet");
		h_secondaryVertex->GetXaxis()->SetBinLabel(3,"lep in Sec. Vtx");
		h_secondaryVertex->GetXaxis()->SetBinLabel(4,"Ch. 3^{rd} Vtx + lep");
		h_secondaryVertex->GetXaxis()->SetBinLabel(5,"N. 3^{rd} Vtx + lep");
		h_secondaryVertex->GetXaxis()->SetBinLabel(6,"jet axis");
		h_secondaryVertex->GetXaxis()->SetBinLabel(7,"lead. par. in jet (Ch.)");
		h_secondaryVertex->GetXaxis()->SetBinLabel(8,"lead. par. in jet (N.)");
		h_secondaryVertex->GetXaxis()->SetBinLabel(9,"lead. par. in jet (#gamma)");
		h_secondaryVertex->GetXaxis()->SetBinLabel(10,"SLD-lep");
		h_secondaryVertex->GetXaxis()->SetBinLabel(11,"lep in Prim. Vtx");
		h_parentHadronCharge = new TH1I( "parentHadronCharge" , "; Parent Hadron Charge" , 5 , 0 , 5 );
		h_parentHadronCharge->GetXaxis()->SetBinLabel(1,"-2");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(2,"-1");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(3,"0");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(4,"1");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(5,"2");
		h_MCPTracks = new TH1I( "chargedMCPTracks" , ";" , 2 , 0 , 2 );
		h_MCPTracks->GetXaxis()->SetBinLabel(1,"track found");
		h_MCPTracks->GetXaxis()->SetBinLabel(2,"track lost");
		h_MCPTracks_Eweighted = new TH1I( "chargedMCPTracks" , "Energy weighted;" , 2 , 0 , 2 );
		h_MCPTracks_Eweighted->GetXaxis()->SetBinLabel(1,"track found");
		h_MCPTracks_Eweighted->GetXaxis()->SetBinLabel(2,"track lost");
		h_MCPTracks_Ptweighted = new TH1I( "chargedMCPTracks" , "p_{T} weighted;" , 2 , 0 , 2 );
		h_MCPTracks_Ptweighted->GetXaxis()->SetBinLabel(1,"track found");
		h_MCPTracks_Ptweighted->GetXaxis()->SetBinLabel(2,"track lost");
		h_FlightDirectionError = new TH1F ( "Flight Direction Error" , "; cos #alpha" , 200 , -1.0 , 1.0 );
		h_distRecoLeptonToDownStreamVertex = new TH1F( "distance of RecoLepton to DownStreamVertex" , ";r^{3D} [mm]", 100 , 0.0 , 10.0 );
	}

	BHadPDGs = { 511 , 521 , 531 , 541 , 5122 , 5132 , 5232 , 5332 };
	CHadPDGs = { 411 , 413 , 415 , 421 , 423 , 425 , 431 , 433 , 435 , 441 , 443 , 4122 , 4132 , 4232 , 4332 , 4334 };
}

void SLDCorrection::Clear()
{
	m_nSLD_chargedMCPwoTrack.clear();
	m_GenStatParentHadron.clear();
	m_ChargeParentHadron.clear();
	m_nTrueNeutralDecayProducts.clear();
	m_nTrueAloneChargedDecayProducts.clear();
	m_nTrueVertices.clear();
	m_nJetsVerticesDistributedIn.clear();
	m_nRecoVerticesInJet.clear();
	m_nAloneChargedPFOs.clear();
	m_nAloneChargedPFOsFromSLD.clear();
	m_distLeptonAlonePFOsNotFromSLD.clear();
	m_distLeptonAlonePFOsFromSLD.clear();
	m_foundRecoLepton.clear();
	m_foundBuildUpVertex.clear();
	m_foundRecoLeptonInBuildUpVertex.clear();
	m_foundRecoLeptonInPrimaryVertex.clear();
	m_lostChargedMCP_CosTheta.clear();
	m_lostChargedMCP_Energy.clear();
	m_lostChargedMCP_Pt.clear();
	m_SLDFlavour.clear();
	m_SLDType.clear();
	m_SLDLeptonID.clear();
	m_leptonE_to_parentE.clear();
	m_otherChargedE_to_parentE.clear();
	m_allChargedE_to_parentE.clear();
	m_neutralE_to_parentE.clear();
	m_neutrino_to_parentE.clear();
	m_nSLDecayOfBHadron = 0;
	m_nSLDecayOfCHadron = 0;
	m_nSLDecayOfTauLepton = 0;
	m_nSLDecayTotal = 0;
	m_nSLDecayToElectron = 0;
	m_nSLDecayToMuon = 0;
	m_nSLDecayToTau = 0;
	m_nTauNeutrino = 0;
	m_nNeutrino = 0;
	m_nChargedPFOwoTrack = 0;
	m_SLDecayXi.clear();
	m_SLDecayYi.clear();
	m_SLDecayZi.clear();
	m_SLDecayRi.clear();
	m_SLDecayXf.clear();
	m_SLDecayYf.clear();
	m_SLDecayZf.clear();
	m_SLDecayRf.clear();
	m_trueNuPx.clear();
	m_trueNuPy.clear();
	m_trueNuPz.clear();
	m_trueNuE.clear();
	m_recoNuCloseInitialPx.clear();
	m_recoNuCloseInitialPy.clear();
	m_recoNuCloseInitialPz.clear();
	m_recoNuCloseInitialE.clear();
	m_recoNuClosePx.clear();
	m_recoNuClosePy.clear();
	m_recoNuClosePz.clear();
	m_recoNuCloseE.clear();
	m_recoNuPosPx.clear();
	m_recoNuPosPy.clear();
	m_recoNuPosPz.clear();
	m_recoNuPosE.clear();
	m_recoNuNegPx.clear();
	m_recoNuNegPy.clear();
	m_recoNuNegPz.clear();
	m_recoNuNegE.clear();
	m_NuPxResidual.clear();
	m_NuPyResidual.clear();
	m_NuPzResidual.clear();
	m_NuEResidual.clear();
	m_NuPxNormalizedResidual.clear();
	m_NuPyNormalizedResidual.clear();
	m_NuPzNormalizedResidual.clear();
	m_NuENormalizedResidual.clear();
	m_solutionSign.clear();
	m_true_E_vis.clear();
	m_true_E_vis_prime.clear();
	m_true_P_vis_par.clear();
	m_true_P_vis_par_prime.clear();
	m_true_P_vis_nor.clear();
	m_E_vis.clear();
	m_E_vis_prime.clear();
	m_P_vis_par.clear();
	m_P_vis_par_prime.clear();
	m_P_vis_par_prime_squared.clear();
	m_P_vis_nor.clear();
	m_P_vis_nor_prime.clear();
	m_flightDirectionStatus.clear();
	m_FlightDirectionErrorCosAlpha.clear();
	m_FlightDirectionErrorSinAlpha.clear();
	m_FlightDirectionErrorAlpha.clear();
	m_dsVertexResidualX.clear();
	m_dsVertexResidualY.clear();
	m_dsVertexResidualZ.clear();
	m_SecVertexResidualX.clear();
	m_SecVertexResidualY.clear();
	m_SecVertexResidualZ.clear();
	m_parentHadronMass.clear();
	m_parentHadronPDG.clear();
	m_trueParentHadronFlightDistance.clear();
	m_recoParentHadronFlightDistance.clear();
	m_daughterHadronMass.clear();
	m_daughterHadronPDG.clear();
	m_daughterHadronFlightDistance.clear();
	m_jetEnergy.clear();
	m_jetEnergyFractionCharged.clear();
	m_jetEnergyFractionNeutralHadron.clear();
	m_jetEnergyFractionPhoton.clear();
	m_jetEnergyFractionNeutrals.clear();
	m_SLDStatus.clear();
	m_PVAStatus.clear();
	m_trueSolutionSign.clear();
	m_recoSolutionSign.clear();
	m_nChargedParticlesInPrimVertex.clear();
	m_nChargedParticlesNotInPrimVertex.clear();
	m_weightPFOtoMCP_Lepton.clear();
	m_weightMCPtoPFO_Lepton.clear();
	m_weightPFOtoMCP_Neutral.clear();
	m_weightMCPtoPFO_Neutral.clear();
	m_weightPFOtoMCP_Charged.clear();
	m_weightMCPtoPFO_Charged.clear();
	m_distRecoLeptonToDownStreamVertex.clear();

	m_visibleChargedInvMassCut.clear();
	m_visibleNeutralInvMassCut.clear();
	m_trueVisibleFourMomentumAtSLDVertex_Px.clear();
	m_trueVisibleFourMomentumAtSLDVertex_Py.clear();
	m_trueVisibleFourMomentumAtSLDVertex_Pz.clear();
	m_trueVisibleFourMomentumAtSLDVertex_E.clear();
	m_trueVisibleFourMomentumAtSLDVertex_M.clear();
	m_truePVATrueFourMomentum_Px.clear();
	m_truePVATrueFourMomentum_Py.clear();
	m_truePVATrueFourMomentum_Pz.clear();
	m_truePVATrueFourMomentum_E.clear();
	m_truePVATrueFourMomentum_M.clear();
	m_truePVARecoFourMomentum_Px.clear();
	m_truePVARecoFourMomentum_Py.clear();
	m_truePVARecoFourMomentum_Pz.clear();
	m_truePVARecoFourMomentum_E.clear();
	m_truePVARecoFourMomentum_M.clear();
	m_recoPVARecoFourMomentum_Px.clear();
	m_recoPVARecoFourMomentum_Py.clear();
	m_recoPVARecoFourMomentum_Pz.clear();
	m_recoPVARecoFourMomentum_E.clear();
	m_recoPVARecoFourMomentum_M.clear();
	m_usedVisibleFourMomentum_Px.clear();
	m_usedVisibleFourMomentum_Py.clear();
	m_usedVisibleFourMomentum_Pz.clear();
	m_usedVisibleFourMomentum_E.clear();
	m_usedVisibleFourMomentum_M.clear();
	m_trueLeptonFourMomentum_Px.clear();
	m_trueLeptonFourMomentum_Py.clear();
	m_trueLeptonFourMomentum_Pz.clear();
	m_trueLeptonFourMomentum_E.clear();
	m_trueLeptonFourMomentum_M.clear();
	m_recoLeptonFourMomentum_Px.clear();
	m_recoLeptonFourMomentum_Py.clear();
	m_recoLeptonFourMomentum_Pz.clear();
	m_recoLeptonFourMomentum_E.clear();
	m_recoLeptonFourMomentum_M.clear();
	m_usedLeptonFourMomentum_Px.clear();
	m_usedLeptonFourMomentum_Py.clear();
	m_usedLeptonFourMomentum_Pz.clear();
	m_usedLeptonFourMomentum_E.clear();
	m_usedLeptonFourMomentum_M.clear();
	m_truePVATrueChargedFourMomentum_Px.clear();
	m_truePVATrueChargedFourMomentum_Py.clear();
	m_truePVATrueChargedFourMomentum_Pz.clear();
	m_truePVATrueChargedFourMomentum_E.clear();
	m_truePVATrueChargedFourMomentum_M.clear();
	m_truePVARecoChargedFourMomentum_Px.clear();
	m_truePVARecoChargedFourMomentum_Py.clear();
	m_truePVARecoChargedFourMomentum_Pz.clear();
	m_truePVARecoChargedFourMomentum_E.clear();
	m_truePVARecoChargedFourMomentum_M.clear();
	m_recoSLDVertexChargedFourMomentum_Px.clear();
	m_recoSLDVertexChargedFourMomentum_Py.clear();
	m_recoSLDVertexChargedFourMomentum_Pz.clear();
	m_recoSLDVertexChargedFourMomentum_E.clear();
	m_recoSLDVertexChargedFourMomentum_M.clear();
	m_recoPVARecoChargedFourMomentum_Px.clear();
	m_recoPVARecoChargedFourMomentum_Py.clear();
	m_recoPVARecoChargedFourMomentum_Pz.clear();
	m_recoPVARecoChargedFourMomentum_E.clear();
	m_recoPVARecoChargedFourMomentum_M.clear();
	m_usedChargedFourMomentum_Px.clear();
	m_usedChargedFourMomentum_Py.clear();
	m_usedChargedFourMomentum_Pz.clear();
	m_usedChargedFourMomentum_E.clear();
	m_usedChargedFourMomentum_M.clear();
	m_truePVATrueNeutralFourMomentum_Px.clear();
	m_truePVATrueNeutralFourMomentum_Py.clear();
	m_truePVATrueNeutralFourMomentum_Pz.clear();
	m_truePVATrueNeutralFourMomentum_E.clear();
	m_truePVATrueNeutralFourMomentum_M.clear();
	m_truePVARecoNeutralFourMomentum_Px.clear();
	m_truePVARecoNeutralFourMomentum_Py.clear();
	m_truePVARecoNeutralFourMomentum_Pz.clear();
	m_truePVARecoNeutralFourMomentum_E.clear();
	m_truePVARecoNeutralFourMomentum_M.clear();
	m_recoPVARecoNeutralFourMomentum_Px.clear();
	m_recoPVARecoNeutralFourMomentum_Py.clear();
	m_recoPVARecoNeutralFourMomentum_Pz.clear();
	m_recoPVARecoNeutralFourMomentum_E.clear();
	m_recoPVARecoNeutralFourMomentum_M.clear();
	m_usedNeutralFourMomentum_Px.clear();
	m_usedNeutralFourMomentum_Py.clear();
	m_usedNeutralFourMomentum_Pz.clear();
	m_usedNeutralFourMomentum_E.clear();
	m_usedNeutralFourMomentum_M.clear();
	m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Px.clear();
	m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Py.clear();
	m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Pz.clear();
	m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_E.clear();
	m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_M.clear();
	m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Px.clear();
	m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Py.clear();
	m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Pz.clear();
	m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_E.clear();
	m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_M.clear();
	m_trueNeutrinoFourMomentum_Px.clear();
	m_trueNeutrinoFourMomentum_Py.clear();
	m_trueNeutrinoFourMomentum_Pz.clear();
	m_trueNeutrinoFourMomentum_E.clear();
	m_trueNeutrinoFourMomentum_M.clear();
	m_recoNeutrinoFourMomentumClose_Px.clear();
	m_recoNeutrinoFourMomentumClose_Py.clear();
	m_recoNeutrinoFourMomentumClose_Pz.clear();
	m_recoNeutrinoFourMomentumClose_E.clear();
	m_recoNeutrinoFourMomentumClose_M.clear();
	m_trueHadronFourMomentum_Px.clear();
	m_trueHadronFourMomentum_Py.clear();
	m_trueHadronFourMomentum_Pz.clear();
	m_trueHadronFourMomentum_E.clear();
	m_trueHadronFourMomentum_M.clear();
	m_recoHadronFourMomentum_Px.clear();
	m_recoHadronFourMomentum_Py.clear();
	m_recoHadronFourMomentum_Pz.clear();
	m_recoHadronFourMomentum_E.clear();
	m_recoHadronFourMomentum_M.clear();
	m_trueFlightDirection_X.clear();
	m_trueFlightDirection_Y.clear();
	m_trueFlightDirection_Z.clear();
	m_recoFlightDirection_X.clear();
	m_recoFlightDirection_Y.clear();
	m_recoFlightDirection_Z.clear();
	m_usedFlightDirection_X.clear();
	m_usedFlightDirection_Y.clear();
	m_usedFlightDirection_Z.clear();
	m_NeutralPVAAlpha.clear();
	m_NeutralPVASinAlpha.clear();
	m_NeutralPVACosAlpha.clear();
	m_ChargedPVAAlpha.clear();
	m_ChargedPVASinAlpha.clear();
	m_ChargedPVACosAlpha.clear();
	m_PVAAlpha.clear();
	m_PVASinAlpha.clear();
	m_PVACosAlpha.clear();
	m_Alpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral.clear();
	m_SinAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral.clear();
	m_CosAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral.clear();
	m_Alpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged.clear();
	m_SinAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged.clear();
	m_CosAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged.clear();
	m_visiblePxBeforePVA.clear();
	m_visiblePyBeforePVA.clear();
	m_visiblePzBeforePVA.clear();
	m_visibleEnergyBeforePVA.clear();
	m_visibleMassBeforePVA.clear();
	m_visiblePxAfterVertexPVA.clear();
	m_visiblePyAfterVertexPVA.clear();
	m_visiblePzAfterVertexPVA.clear();
	m_visibleEnergyAfterVertexPVA.clear();
	m_visibleMassAfterVertexPVA.clear();
	m_visiblePxAfterChargedPVA.clear();
	m_visiblePyAfterChargedPVA.clear();
	m_visiblePzAfterChargedPVA.clear();
	m_visibleEnergyAfterChargedPVA.clear();
	m_visibleMassAfterChargedPVA.clear();
	m_visiblePxAfterNeutralPVA.clear();
	m_visiblePyAfterNeutralPVA.clear();
	m_visiblePzAfterNeutralPVA.clear();
	m_visibleEnergyAfterNeutralPVA.clear();
	m_visibleMassAfterNeutralPVA.clear();

	m_trueVisibleFourMomentumPx.clear();
	m_trueVisibleFourMomentumPy.clear();
	m_trueVisibleFourMomentumPz.clear();
	m_trueVisibleFourMomentumE.clear();
	m_recoVisibleFourMomentumPx.clear();
	m_recoVisibleFourMomentumPy.clear();
	m_recoVisibleFourMomentumPz.clear();
	m_recoVisibleFourMomentumE.clear();
	m_residualVisibleFourMomentumPx.clear();
	m_residualVisibleFourMomentumPy.clear();
	m_residualVisibleFourMomentumPz.clear();
	m_residualVisibleFourMomentumE.clear();
	m_sigmaPxPx_Det.clear();
	m_sigmaPxPy_Det.clear();
	m_sigmaPyPy_Det.clear();
	m_sigmaPxPz_Det.clear();
	m_sigmaPyPz_Det.clear();
	m_sigmaPzPz_Det.clear();
	m_sigmaPxE_Det.clear();
	m_sigmaPyE_Det.clear();
	m_sigmaPzE_Det.clear();
	m_sigmaEE_Det.clear();
	m_normalizedResidualVisibleFourMomentumPx.clear();
	m_normalizedResidualVisibleFourMomentumPy.clear();
	m_normalizedResidualVisibleFourMomentumPz.clear();
	m_normalizedResidualVisibleFourMomentumE.clear();
	m_trueFlightDirectionUx.clear();
	m_trueFlightDirectionUy.clear();
	m_trueFlightDirectionUz.clear();
	m_recoFlightDirectionUx.clear();
	m_recoFlightDirectionUy.clear();
	m_recoFlightDirectionUz.clear();
	m_residualFlightDirectionUx.clear();
	m_residualFlightDirectionUy.clear();
	m_residualFlightDirectionUz.clear();
	m_sigmaUxUx.clear();
	m_sigmaUxUy.clear();
	m_sigmaUyUy.clear();
	m_sigmaUxUz.clear();
	m_sigmaUyUz.clear();
	m_sigmaUzUz.clear();
	m_normalizedResidualFlightDirectionUx.clear();
	m_normalizedResidualFlightDirectionUy.clear();
	m_normalizedResidualFlightDirectionUz.clear();
	m_trueNeutrinoFourMomentumPx.clear();
	m_trueNeutrinoFourMomentumPy.clear();
	m_trueNeutrinoFourMomentumPz.clear();
	m_trueNeutrinoFourMomentumE.clear();
	m_recoNeutrinoFourMomentumClosePx.clear();
	m_recoNeutrinoFourMomentumClosePy.clear();
	m_recoNeutrinoFourMomentumClosePz.clear();
	m_recoNeutrinoFourMomentumCloseE.clear();
	m_sigmaNeutrinoPxPx.clear();
	m_sigmaNeutrinoPxPy.clear();
	m_sigmaNeutrinoPyPy.clear();
	m_sigmaNeutrinoPxPz.clear();
	m_sigmaNeutrinoPyPz.clear();
	m_sigmaNeutrinoPzPz.clear();
	m_sigmaNeutrinoPxE.clear();
	m_sigmaNeutrinoPyE.clear();
	m_sigmaNeutrinoPzE.clear();
	m_sigmaNeutrinoEE.clear();
	m_residualNeutrinoFourMomentumPx.clear();
	m_residualNeutrinoFourMomentumPy.clear();
	m_residualNeutrinoFourMomentumPz.clear();
	m_residualNeutrinoFourMomentumE.clear();
	m_normalizedResidualNeutrinoFourMomentumPx.clear();
	m_normalizedResidualNeutrinoFourMomentumPy.clear();
	m_normalizedResidualNeutrinoFourMomentumPz.clear();
	m_normalizedResidualNeutrinoFourMomentumE.clear();
	m_recoNeutrinoDirectionError.clear();
	m_PCAatLeptonX.clear();
	m_PCAatLeptonY.clear();
	m_PCAatLeptonZ.clear();
	m_PCAatOtherParticleX.clear();
	m_PCAatOtherParticleY.clear();
	m_PCAatOtherParticleZ.clear();
	m_JetAxisX.clear();
	m_JetAxisY.clear();
	m_JetAxisZ.clear();
	m_PrimaryVertexX.clear();
	m_PrimaryVertexY.clear();
	m_PrimaryVertexZ.clear();
	m_AngleTrueFlightDirectionJet.clear();
	m_AngleRecoFlightDirectionJet.clear();
	m_AngleLeptonJet.clear();
	m_AngleDSVertexJet.clear();
	m_AngleRecoNeutrinoJet.clear();
	m_AngleTrueNeutrinoJet.clear();
	m_LeptonDistanceFromPV.clear();
	m_DSVDistanceFromPV.clear();
	m_Lepton3DImpactParameter.clear();
	m_OtherParticle3DImpactParameter.clear();
}

void SLDCorrection::processRunHeader()
{
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
}

void SLDCorrection::processEvent( EVENT::LCEvent *pLCEvent )
{
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	bool b_BSLDMode = false;
	bool b_CSLDMode = false;
	bool b_TSLDMode = false;
	bool b_SLDMode = false;

	LCCollection *MCParticleCollection{};
	LCCollection *jetCollection{};
	int nTauNeutrino = 0;
	++m_nEvtSum;
	Clear();
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

	IMPL::LCCollectionVec* semiLeptonicVertex(NULL);
	semiLeptonicVertex = new IMPL::LCCollectionVec( LCIO::VERTEX );
	//semiLeptonicVertex->setSubset( true );
	IMPL::LCCollectionVec* semiLeptonicVertexRP(NULL);
	semiLeptonicVertexRP = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	//semiLeptonicVertexRP->setSubset( true );
	IMPL::LCCollectionVec* Neutrinos(NULL);
	Neutrinos = new IMPL::LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
	EVENT::LCCollection* JetSLDLink(NULL);
	//JetSLDLink = new EVENT::LCCollection( LCIO::LCRELATION );
	EVENT::LCCollection* SLDJetLink(NULL);
	//SLDJetLink = new EVENT::LCCollection( LCIO::LCRELATION );
	EVENT::LCCollection* mcNurecoNuLink(NULL);
	EVENT::LCCollection* recoNumcNuLink(NULL);
	EVENT::LCCollection* NuSLDLink(NULL);
	EVENT::LCCollection* SLDNuLink(NULL);

	LCRelationNavigator JetSLDRelNav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::VERTEX  );
	LCRelationNavigator SLDJetRelNav(LCIO::VERTEX , LCIO::RECONSTRUCTEDPARTICLE  );
	LCRelationNavigator NeutrinoSLDRelNav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::VERTEX  );
	LCRelationNavigator SLDNeutrinoRelNav(LCIO::VERTEX , LCIO::RECONSTRUCTEDPARTICLE  );
	LCRelationNavigator MCNuRecoNuRelNav(LCIO::MCPARTICLE , LCIO::RECONSTRUCTEDPARTICLE  );
	LCRelationNavigator RecoNuMCNuRelNav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE  );

	try
	{
		vtxVector BsemiLeptonicVertices{};
		vtxVector CsemiLeptonicVertices{};
		vtxVector TsemiLeptonicVertices{};
		vtxVector semiLeptonicVertices{};
		pfoVector semiLeptonicVertexRecoParticles{};
		pfoVectorVector neutrinos{};
		pfoVector jetsOfSemiLeptonicDecays{};
		mcpVector mcNeutrinos{};
		IntVector SLDStatus{};
		IntVector PVAStatus{};
		IntVector solutionSigns{};

		vtxVector tempSemiLeptonicVertices{};
		pfoVector tempSemiLeptonicVertexRecoParticles{};
		pfoVector tempNeutrinos{};
		pfoVector tempJetsOfSemiLeptonicDecays{};

		MCParticleCollection = pLCEvent->getCollection( m_mcParticleCollection );
		jetCollection = pLCEvent->getCollection( m_inputJetCollection );
		int nMCP = MCParticleCollection->getNumberOfElements();
		for ( int i_mcp = 0 ; i_mcp < nMCP ; ++i_mcp )
		{
			MCParticle *testLepton = dynamic_cast<EVENT::MCParticle*>( MCParticleCollection->getElementAt( i_mcp ) );
			if ( testLepton->isOverlay() ) continue;
			if ( abs( testLepton->getPDG() ) == 16 && ( testLepton->getGeneratorStatus() ) == 1 ) ++nTauNeutrino;
			bool primarySLDecay = false;
			bool downStreamSLDecay = false;
			bool upStreamSLDecay = false;
			bool isBHadronSLDecay = false;
			bool isCHadronSLDecay = false;
			bool isTauLeptonSLDecay = false;
			if ( ( abs( testLepton->getPDG() ) == 11 || abs( testLepton->getPDG() ) == 13 || abs( testLepton->getPDG() ) == 15 ) )
			{
				int chargedLeptonPDG;
				for ( long unsigned int i_parent = 0 ; i_parent < ( testLepton->getParents() ).size() ; ++i_parent )
				{
					MCParticle *parent = testLepton->getParents()[ i_parent ];
					primarySLDecay = hasPrimarySLDecay( parent , chargedLeptonPDG );
					if ( primarySLDecay ) downStreamSLDecay = hasDownStreamSLDecay( parent );
					if ( primarySLDecay ) upStreamSLDecay = hasUpStreamSLDecay( parent );
				}
				if ( primarySLDecay )
				{
					bool solveSLD = true;
					std::vector< TLorentzVector > recoNeutrinoFourMomentum;
					std::vector< std::vector< float > > recoNeutrinoCovMat;
					isBHadronSLDecay = checkBHadronSLDecay( testLepton );
					isCHadronSLDecay = checkCHadronSLDecay( testLepton );
					isTauLeptonSLDecay = checkTauLeptonSLDecay( testLepton );
					if ( isBHadronSLDecay )
					{
						++m_nSLDecayOfBHadron;
						m_SLDFlavour.push_back( 5 );
						if ( m_fillRootTree ) h_SLDecayFlavour->Fill( 0.5 );
					}
					else if ( isCHadronSLDecay )
					{
						++m_nSLDecayOfCHadron;
						m_SLDFlavour.push_back( 4 );
						if ( m_fillRootTree ) h_SLDecayFlavour->Fill( 1.5 );
					}
					else if ( isTauLeptonSLDecay )
					{
						++m_nSLDecayOfTauLepton;
						m_SLDFlavour.push_back( 15 );
						if ( m_fillRootTree ) h_SLDecayFlavour->Fill( 2.5 );
					}
					else
					{
						m_SLDFlavour.push_back( 0 );
					}
					if ( downStreamSLDecay || upStreamSLDecay )
					{
						m_SLDType.push_back( 0 );
						solveSLD = false;
					}
					else
					{
						m_SLDType.push_back( 1 );
					}
					if ( abs( chargedLeptonPDG ) == 11 )
					{
						++m_nSLDecayToElectron;
						if ( isBHadronSLDecay && m_fillRootTree ) h_SLDecayModeB->Fill( 0.5 );
						if ( isCHadronSLDecay && m_fillRootTree ) h_SLDecayModeC->Fill( 0.5 );
						solveSLD = solveSLD && ( chargedLeptonPDG == testLepton->getPDG() );
					}
					else if ( abs( chargedLeptonPDG ) == 13 )
					{
						++m_nSLDecayToMuon;
						if ( isBHadronSLDecay && m_fillRootTree ) h_SLDecayModeB->Fill( 1.5 );
						if ( isCHadronSLDecay && m_fillRootTree ) h_SLDecayModeC->Fill( 1.5 );
						solveSLD = solveSLD && ( chargedLeptonPDG == testLepton->getPDG() );
					}
					else if ( abs( chargedLeptonPDG ) == 15 )
					{
						++m_nSLDecayToTau;
						if ( isBHadronSLDecay && m_fillRootTree ) h_SLDecayModeB->Fill( 2.5 );
						if ( isCHadronSLDecay && m_fillRootTree ) h_SLDecayModeC->Fill( 2.5 );
						solveSLD = solveSLD && ( chargedLeptonPDG == testLepton->getPDG() );
					}
					m_SLDecayXi.push_back( testLepton->getParents()[ 0 ]->getVertex()[ 0 ] );
					m_SLDecayYi.push_back( testLepton->getParents()[ 0 ]->getVertex()[ 1 ] );
					m_SLDecayZi.push_back( testLepton->getParents()[ 0 ]->getVertex()[ 2 ] );
					m_SLDecayRi.push_back( sqrt( pow( testLepton->getParents()[ 0 ]->getVertex()[ 0 ] , 2 ) + pow( testLepton->getParents()[ 0 ]->getVertex()[ 1 ] , 2 ) + pow( testLepton->getParents()[ 0 ]->getVertex()[ 2 ] , 2 ) ) );
					m_SLDecayXf.push_back( testLepton->getParents()[ 0 ]->getEndpoint()[ 0 ] );
					m_SLDecayYf.push_back( testLepton->getParents()[ 0 ]->getEndpoint()[ 1 ] );
					m_SLDecayZf.push_back( testLepton->getParents()[ 0 ]->getEndpoint()[ 2 ] );
					m_SLDecayRf.push_back( sqrt( pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 0 ] , 2 ) + pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 1 ] , 2 ) + pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 2 ] , 2 ) ) );
					if ( isBHadronSLDecay && !m_includeBSLD ) solveSLD = false;
					if ( isCHadronSLDecay && !m_includeCSLD ) solveSLD = false;
					if ( isTauLeptonSLDecay && !m_includeTSLD ) solveSLD = false;
					if ( abs( testLepton->getPDG() ) == 15 ) solveSLD = false;
					streamlog_out(DEBUG3) << "" << std::endl;
					streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<< Found a primary semi-leptonic decay >>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					if ( downStreamSLDecay )
					{
						streamlog_out(DEBUG3) << "	There is/are downstream semi-leptonic(s) decay in primary semi-leptonic decay products" << std::endl;
					}
					if ( upStreamSLDecay )
					{
						streamlog_out(DEBUG3) << "	There is/are upstream semi-leptonic(s) decay in primary semi-leptonic decay products" << std::endl;
					}
					if ( downStreamSLDecay || upStreamSLDecay )
					{
						if ( m_fillRootTree ) h_SLDecayOrder->Fill( 0.5 );
					}
					if ( solveSLD )
					{
						m_SLDLeptonID.push_back( testLepton->getPDG() );
						if ( m_fillRootTree ) h_SLDecayOrder->Fill( 1.5 );
						streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
						streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<< There are no upstream and downstream semi-leptonic decay >>>>>>>>>>>>>>>>>" << std::endl;
						streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
						doSLDCorrection( pLCEvent , testLepton , BsemiLeptonicVertices , semiLeptonicVertexRecoParticles , jetsOfSemiLeptonicDecays , neutrinos , SLDStatus , PVAStatus , solutionSigns , mcNeutrinos );
						m_parentHadronMass.push_back( ( testLepton->getParents()[ 0 ] )->getMass() );
						m_parentHadronPDG.push_back( ( testLepton->getParents()[ 0 ] )->getPDG() );
						for ( unsigned int i_Btype = 0 ; i_Btype < BHadPDGs.size() ; ++i_Btype )
						{
							if ( abs( ( testLepton->getParents()[ 0 ] )->getPDG() ) == BHadPDGs[ i_Btype ] && m_fillRootTree ) h_BHadronType->Fill( i_Btype );
						}
						m_trueParentHadronFlightDistance.push_back( std::sqrt( std::pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 0 ] - testLepton->getParents()[ 0 ]->getVertex()[ 0 ] , 2 ) + std::pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 1 ] - testLepton->getParents()[ 0 ]->getVertex()[ 1 ] , 2 ) + std::pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 2 ] - testLepton->getParents()[ 0 ]->getVertex()[ 2 ] , 2 ) ) );
						for ( unsigned int i_d = 0 ; i_d < ( testLepton->getParents()[ 0 ] )->getDaughters().size() ; ++i_d )
						{
							MCParticle *mcDaughter = ( testLepton->getParents()[ 0 ] )->getDaughters()[ i_d ];
							int daughterPDG = std::abs( mcDaughter->getPDG() );
							if ( daughterPDG < 11 || daughterPDG > 16 )
							{
								m_daughterHadronMass.push_back( mcDaughter->getMass() );
								m_daughterHadronPDG.push_back( mcDaughter->getPDG() );
								for ( unsigned int i_Ctype = 0 ; i_Ctype < CHadPDGs.size() ; ++i_Ctype )
								{
									if ( abs( mcDaughter->getPDG() ) == CHadPDGs[ i_Ctype ] && m_fillRootTree ) h_CHadronType->Fill( i_Ctype );
								}
							}
						}
					}
				}
			}
		}
		m_nTauNeutrino = nTauNeutrino;
		m_nSLDecayTotal = m_nSLDecayOfBHadron + m_nSLDecayOfCHadron;

		if ( m_fillRootTree )
		{
			m_pTTree1->Fill();
			m_pTTree2->Fill();
			m_pTTree3->Fill();
		}
		for ( unsigned int bSLD = 0 ; bSLD < BsemiLeptonicVertices.size() ; ++bSLD ) semiLeptonicVertices.push_back( BsemiLeptonicVertices[ bSLD ] );
		for ( unsigned int cSLD = 0 ; cSLD < CsemiLeptonicVertices.size() ; ++cSLD ) semiLeptonicVertices.push_back( CsemiLeptonicVertices[ cSLD ] );
		for ( unsigned int tSLD = 0 ; tSLD < TsemiLeptonicVertices.size() ; ++tSLD ) semiLeptonicVertices.push_back( TsemiLeptonicVertices[ tSLD ] );
		if ( m_BSLDMode == 0 && m_nSLDecayOfBHadron == 0 ) b_BSLDMode = true;
		if ( m_BSLDMode == 1 ) b_BSLDMode = true;
		if ( m_BSLDMode == 2 && m_nSLDecayOfBHadron > 0 ) b_BSLDMode = true;
		if ( m_BSLDMode == 3 && BsemiLeptonicVertices.size() == m_nSLDecayOfBHadron && BsemiLeptonicVertices.size() > 0 ) b_BSLDMode = true;
		if ( m_CSLDMode == 0 && m_nSLDecayOfCHadron == 0 ) b_CSLDMode = true;
		if ( m_CSLDMode == 1 ) b_CSLDMode = true;
		if ( m_CSLDMode == 2 && m_nSLDecayOfCHadron > 0 ) b_CSLDMode = true;
		if ( m_CSLDMode == 3 && CsemiLeptonicVertices.size() == m_nSLDecayOfCHadron && CsemiLeptonicVertices.size() > 0 ) b_CSLDMode = true;
		if ( m_TSLDMode == 0 && m_nSLDecayOfTauLepton == 0 ) b_TSLDMode = true;
		if ( m_TSLDMode == 1 ) b_TSLDMode = true;
		if ( m_TSLDMode == 2 && m_nSLDecayOfTauLepton > 0 ) b_TSLDMode = true;
		if ( m_TSLDMode == 3 && TsemiLeptonicVertices.size() == m_nSLDecayOfTauLepton && TsemiLeptonicVertices.size() > 0 ) b_TSLDMode = true;
		if ( m_SLDMode == 0 && m_nSLDecayTotal == 0 ) b_SLDMode = true;
		if ( m_SLDMode == 1 ) b_SLDMode = true;
		if ( m_SLDMode == 2 && m_nSLDecayTotal > 0 ) b_SLDMode = true;
		if ( m_SLDMode == 3 && semiLeptonicVertices.size() == m_nSLDecayTotal && semiLeptonicVertices.size() > 0 ) b_SLDMode = true;
		streamlog_out(DEBUG9) << "	Number of BSLD: 	Found = " << m_nSLDecayOfBHadron << " , Corrected = " << BsemiLeptonicVertices.size() << std::endl;
		streamlog_out(DEBUG9) << "	Number of CSLD: 	Found = " << m_nSLDecayOfCHadron << " , Corrected = " << CsemiLeptonicVertices.size() << std::endl;
		streamlog_out(DEBUG9) << "	Number of TSLD: 	Found = " << m_nSLDecayOfTauLepton << " , Corrected = " << TsemiLeptonicVertices.size() << std::endl;
		streamlog_out(DEBUG9) << "	Number of SLD: 		Found = " << m_nSLDecayTotal << " , Corrected = " << semiLeptonicVertices.size() << std::endl;
		streamlog_out(DEBUG9) << "	Event selected by number of BSLD:		" << ( b_BSLDMode ? "TRUE" : "FALSE" ) << std::endl;
		streamlog_out(DEBUG9) << "	Event selected by number of CSLD:		" << ( b_CSLDMode ? "TRUE" : "FALSE" ) << std::endl;
		streamlog_out(DEBUG9) << "	Event selected by number of TSLD:		" << ( b_TSLDMode ? "TRUE" : "FALSE" ) << std::endl;
		streamlog_out(DEBUG9) << "	Event selected by number of SLD:		" << ( b_SLDMode ? "TRUE" : "FALSE" ) << std::endl;


		semiLeptonicVertex->parameters().setValue( "nBHadronSLD_found" , ( int )m_nSLDecayOfBHadron );
		semiLeptonicVertex->parameters().setValue( "nCHadronSLD_found" , ( int )m_nSLDecayOfCHadron );
		semiLeptonicVertex->parameters().setValue( "nTauLeptonSLD_found" , ( int )m_nSLDecayOfTauLepton );
		semiLeptonicVertex->parameters().setValue( "nTotalSLD_found" , ( int )m_nSLDecayTotal );
		semiLeptonicVertex->parameters().setValue( "nBHadronSLD_solved" , ( int )( BsemiLeptonicVertices.size() ) );
		semiLeptonicVertex->parameters().setValue( "nCHadronSLD_solved" , ( int )( CsemiLeptonicVertices.size() ) );
		semiLeptonicVertex->parameters().setValue( "nTauLeptonSLD_solved" , ( int )( TsemiLeptonicVertices.size() ) );
		semiLeptonicVertex->parameters().setValue( "nSolvedSLD" , ( int )( semiLeptonicVertices.size() ) );
		semiLeptonicVertex->parameters().setValues( "SLDLepStatus" , ( IntVec )SLDStatus );
		semiLeptonicVertex->parameters().setValues( "PVAStatus" , ( IntVec )PVAStatus );
		semiLeptonicVertex->parameters().setValues( "trueSolutionSign" , ( IntVec )solutionSigns );

		jetCollection->parameters().setValue( "nBHadronSLD_found" , ( int )m_nSLDecayOfBHadron );
		jetCollection->parameters().setValue( "nCHadronSLD_found" , ( int )m_nSLDecayOfCHadron );
		jetCollection->parameters().setValue( "nTauLeptonSLD_found" , ( int )m_nSLDecayOfTauLepton );
		jetCollection->parameters().setValue( "nTotalSLD_found" , ( int )m_nSLDecayTotal );
		jetCollection->parameters().setValue( "nBHadronSLD_solved" , ( int )( BsemiLeptonicVertices.size() ) );
		jetCollection->parameters().setValue( "nCHadronSLD_solved" , ( int )( CsemiLeptonicVertices.size() ) );
		jetCollection->parameters().setValue( "nTauLeptonSLD_solved" , ( int )( TsemiLeptonicVertices.size() ) );
		jetCollection->parameters().setValue( "nSolvedSLD" , ( int )( semiLeptonicVertices.size() ) );
		jetCollection->parameters().setValues( "SLDLepStatus" , ( IntVec )SLDStatus );
		jetCollection->parameters().setValues( "PVAStatus" , ( IntVec )PVAStatus );
		jetCollection->parameters().setValues( "trueSolutionSign" , ( IntVec )solutionSigns );

		for ( unsigned int i_sld = 0 ; i_sld < semiLeptonicVertices.size() ; ++i_sld )
		{
			semiLeptonicVertex->addElement( semiLeptonicVertices[ i_sld ] );
			semiLeptonicVertexRP->addElement( semiLeptonicVertexRecoParticles[ i_sld ] );
			for ( unsigned int i_nu = 0 ; i_nu < neutrinos[ i_sld].size() ; ++i_nu )
			{
				Neutrinos->addElement( neutrinos[ i_sld][ i_nu ] );
				NeutrinoSLDRelNav.addRelation( neutrinos[ i_sld][ i_nu ] , semiLeptonicVertices[ i_sld ] , 1.0 );
				SLDNeutrinoRelNav.addRelation( semiLeptonicVertices[ i_sld ] , neutrinos[ i_sld][ i_nu ] , 1.0 );
			}
			JetSLDRelNav.addRelation( jetsOfSemiLeptonicDecays[ i_sld ] , semiLeptonicVertices[ i_sld ] , 1.0 );
			SLDJetRelNav.addRelation( semiLeptonicVertices[ i_sld ] , jetsOfSemiLeptonicDecays[ i_sld ] , 1.0 );
		}
		if ( mcNeutrinos.size() == semiLeptonicVertices.size() )
		{
			for ( unsigned int i_sld = 0 ; i_sld < mcNeutrinos.size() ; ++i_sld )
			{
				MCNuRecoNuRelNav.addRelation( mcNeutrinos[ i_sld ] , neutrinos[ i_sld ][ 0 ] , 1.0 );
				RecoNuMCNuRelNav.addRelation( neutrinos[ i_sld ][ 0 ] , mcNeutrinos[ i_sld ] , 1.0 );
			}
			mcNurecoNuLink = MCNuRecoNuRelNav.createLCCollection();
			recoNumcNuLink = RecoNuMCNuRelNav.createLCCollection();
		}
		JetSLDLink = JetSLDRelNav.createLCCollection();
		SLDJetLink = SLDJetRelNav.createLCCollection();
		NuSLDLink = NeutrinoSLDRelNav.createLCCollection();
		SLDNuLink = SLDNeutrinoRelNav.createLCCollection();
		pLCEvent->addCollection( semiLeptonicVertex , m_SLDVertex );
		pLCEvent->addCollection( semiLeptonicVertexRP , m_SLDVertexRP );
		pLCEvent->addCollection( Neutrinos , m_reconstructedNeutrino );
		pLCEvent->addCollection( JetSLDLink , m_JetSLDLinkName );
		pLCEvent->addCollection( SLDJetLink , m_SLDJetLinkName );
		pLCEvent->addCollection( NuSLDLink , m_NuSLDLinkName );
		pLCEvent->addCollection( SLDNuLink , m_SLDNuLinkName );
		if ( mcNeutrinos.size() == semiLeptonicVertices.size() )
		{
			pLCEvent->addCollection( mcNurecoNuLink , m_mcNurecoNuLinkName );
			pLCEvent->addCollection( recoNumcNuLink , m_recoNumcNuLinkName );
		}

		streamlog_out(MESSAGE0) << "	Found " << semiLeptonicVertices.size() << " semi-leptonic decays with " << semiLeptonicVertexRecoParticles.size() << " associated reconstructed particles in " << jetsOfSemiLeptonicDecays.size() << " jets" << std::endl;
		for ( unsigned int i_sld = 0 ; i_sld < semiLeptonicVertices.size() ; ++i_sld )
		{
			streamlog_out(DEBUG9) << "		SLD[ " << i_sld << " ]:  Status of Charged Lepton: " << SLDStatus[ i_sld ] << std::endl;
			streamlog_out(DEBUG9) << "		SLD[ " << i_sld << " ]:  Status of Association of Particles to SLD Vertex: " << PVAStatus[ i_sld ] << std::endl;
		}
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "	Input collection not found in event " << m_nEvt << std::endl;
	}
	setReturnValue( "BSLDMode" , b_BSLDMode );
	setReturnValue( "CSLDMode" , b_CSLDMode );
	setReturnValue( "TSLDMode" , b_TSLDMode );
	setReturnValue( "SLDMode" , b_SLDMode );

}

bool SLDCorrection::hasPrimarySLDecay( MCParticle *parentHadron , int &chargedLeptonPDG )
{
	bool hasSLDecay = false;
	if ( parentHadron->getGeneratorStatus() == 2 && ( floor( abs( parentHadron->getPDG() ) / 100 ) == 5 || ( floor( abs( parentHadron->getPDG() ) / 1000 ) == 5 ) || floor( abs( parentHadron->getPDG() ) / 100 ) == 4 || ( floor( abs( parentHadron->getPDG() ) / 1000 ) == 4 ) || ( abs( parentHadron->getPDG() ) == 15 ) ) )
	{
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentHadron->getDaughters() ).size() ; ++i_daughter )
		{
			MCParticle *daughter = parentHadron->getDaughters()[ i_daughter ];
			if ( daughter->isOverlay() ) continue;
			if ( daughter->getGeneratorStatus() == 1 )
			{
				if ( abs( daughter->getPDG() ) == 11 || abs( daughter->getPDG() ) == 13 || abs( daughter->getPDG() ) == 15 )
				{
					int leptonCharge = ( int ) daughter->getCharge();
					int expectedNeutrinoPDG = -1 * ( daughter->getPDG() - leptonCharge );
					for ( long unsigned int i_NuCandidate = 0 ; i_NuCandidate < ( parentHadron->getDaughters() ).size() ; ++i_NuCandidate )
					{
						MCParticle *secondDaughter = parentHadron->getDaughters()[ i_NuCandidate ];
						if ( secondDaughter->getPDG() == expectedNeutrinoPDG && secondDaughter->getGeneratorStatus() == 1 )
						{
							hasSLDecay = true;
							chargedLeptonPDG = daughter->getPDG();
						}
					}
				}
			}
			else if( abs( daughter->getPDG() ) == 15 )
			{
				int leptonCharge = ( int ) daughter->getCharge();
				int expectedNeutrinoPDG = -1 * ( daughter->getPDG() - leptonCharge );
				for ( long unsigned int i_NuCandidate = 0 ; i_NuCandidate < ( parentHadron->getDaughters() ).size() ; ++i_NuCandidate )
				{
					MCParticle *secondDaughter = parentHadron->getDaughters()[ i_NuCandidate ];
					if ( secondDaughter->getPDG() == expectedNeutrinoPDG && secondDaughter->getGeneratorStatus() == 1 )
					{
						hasSLDecay = true;
						chargedLeptonPDG = daughter->getPDG();
					}
				}
			}
		}
	}
	return hasSLDecay;
}

bool SLDCorrection::hasDownStreamSLDecay( MCParticle *parentHadron )
{
	bool hasSLDecay = false;
	bool primarySLDecay = false;
	bool downStreamSLDecay = false;
	for ( long unsigned int i_primaryDaughter = 0 ; i_primaryDaughter < ( parentHadron->getDaughters() ).size() ; ++i_primaryDaughter )
	{
		int lepPDG;
		MCParticle *primaryDaughter = parentHadron->getDaughters()[ i_primaryDaughter ];
		primarySLDecay = primarySLDecay || hasPrimarySLDecay( primaryDaughter , lepPDG );
		downStreamSLDecay = downStreamSLDecay || hasDownStreamSLDecay( primaryDaughter );
	}
	hasSLDecay = primarySLDecay || downStreamSLDecay ;
	return hasSLDecay;
}

bool SLDCorrection::hasUpStreamSLDecay( MCParticle *parentHadron )
{
	bool hasSLDecay = false;
	bool primarySLDecay = false;
	bool upStreamSLDecay = false;
	for ( long unsigned int i_upperParent = 0 ; i_upperParent < ( parentHadron->getParents() ).size() ; ++i_upperParent )
	{
		int lepPDG;
		MCParticle *upperParent = parentHadron->getParents()[ i_upperParent ];
		primarySLDecay = primarySLDecay || hasPrimarySLDecay( upperParent , lepPDG );
		upStreamSLDecay = upStreamSLDecay || hasUpStreamSLDecay( upperParent );
	}
	hasSLDecay = primarySLDecay || upStreamSLDecay ;
	return hasSLDecay;
}

bool SLDCorrection::checkBHadronSLDecay( MCParticle *SLDLepton )
{
	bool isBHadronSLDecay = false;
	MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
	if ( floor( fabs( parentHadron->getPDG() ) / 100 ) == 5 || floor( fabs( parentHadron->getPDG() ) / 1000 ) == 5 ) isBHadronSLDecay = true;
	return isBHadronSLDecay;
}

bool SLDCorrection::checkCHadronSLDecay( MCParticle *SLDLepton )
{
	bool isCHadronSLDecay = false;
	MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
	if ( floor( fabs( parentHadron->getPDG() ) / 100 ) == 4 || floor( fabs( parentHadron->getPDG() ) / 1000 ) == 4 ) isCHadronSLDecay = true;
	return isCHadronSLDecay;
}

bool SLDCorrection::checkTauLeptonSLDecay( MCParticle *SLDLepton )
{
	bool TauLeptonSLDecay = false;
	MCParticle *parent = SLDLepton->getParents()[ 0 ];
	if ( abs( parent->getPDG() ) == 15 ) TauLeptonSLDecay = true;
	return TauLeptonSLDecay;
}

void SLDCorrection::doSLDCorrection( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton , vtxVector& semiLeptonicVertices , pfoVector& semiLeptonicVertexRecoParticles , pfoVector& jetsOfSemiLeptonicDecays , pfoVectorVector& neutrinos , IntVector &sldStatus , IntVector &pvaStatus , IntVector &solutionSigns , mcpVector &trueNeutrinos )
{
	++n_SLDStatus;
	showTrueParameters( SLDLepton );
	pfoVector neutrinosOfThisSLD{};

	VertexImpl *semiLeptonicVertex = new VertexImpl;
	ReconstructedParticleImpl *semiLeptonicVertexRecoParticle = new ReconstructedParticleImpl;
	ReconstructedParticleImpl *recoNeutrinoZero = new ReconstructedParticleImpl;
	ReconstructedParticleImpl *recoNeutrinoPos = new ReconstructedParticleImpl;
	ReconstructedParticleImpl *recoNeutrinoNeg = new ReconstructedParticleImpl;

	Vertex *SLDVertex = NULL;
	ReconstructedParticle *SLDVertexRP = NULL;

	mcpVector trueDecayProducts{};
	mcpVector trueNeutralDecayProducts{};
	mcpVector trueChargedDecayProducts{};
	pfoVector truePVAChargedDecayProducts{};
	pfoVector truePVANeutralDecayProducts{};
	pfoVector recoSLDVertexDecayProducts{};
	pfoVector recoPVAVertexDecayProducts{};
	pfoVector recoPVAChargedDecayProducts{};
	pfoVector recoPVANeutralDecayProducts{};
	pfoVector recoPVADecayProducts{};

	int nChargedParticlesInPrimVertex = 0;
	int nChargedParticlesNotInPrimVertex = 0;


	pfoVector jetVector{};
	pfoVector allPFOsInJet{};
	pfoVector chargedPFOsInJet{};
	pfoVector aloneChargedPFOsInJet{};
	pfoVector neutralPFOsInJet{};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////																		////
////			Get Primary Information from event							////
////																		////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( m_RecoMCTruthLinkCollection ) );
	LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( m_MCTruthRecoLinkCollection ) );
	LCCollection *primaryVertexCollection = pLCEvent->getCollection( m_inputPrimaryVertex );
	Vertex* primaryVertex = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
	Vertex* startVertex = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
	LCCollection *jetCollection = pLCEvent->getCollection( m_inputJetCollection );
	for ( int i_jet = 0 ; i_jet < jetCollection->getNumberOfElements(); ++i_jet)
	{
		ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( jetCollection->getElementAt( i_jet ) );
		jetVector.push_back( jet );
	}
	LCCollection *buildUpVertexCollection = pLCEvent->getCollection( m_inputBuildUpVertex );
	vtxVector buildUpVertexVector{};
	for ( int i_vtx = 0 ; i_vtx < buildUpVertexCollection->getNumberOfElements(); ++i_vtx)
	{
		Vertex* vtx = dynamic_cast<Vertex*>( buildUpVertexCollection->getElementAt( i_vtx ) );
		buildUpVertexVector.push_back( vtx );
	}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////																		////
////		Draw semi-leptonic decay with MCParticles						////
////																		////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	if ( m_displayEvent )
	{
		DDMarlinCED::newEvent( this ); // refresh
		DDMarlinCED::drawDD4hepDetector( this->_theDetector , 0 , std::vector<std::string>{} ); // draw geometry
		DDCEDPickingHandler& pHandler = DDCEDPickingHandler::getInstance();
		pHandler.update(pLCEvent);
		drawMCParticles( SLDLepton->getParents()[ 0 ] , SLDLepton->getParents()[ 0 ] );
	}
	//m_displayEvent = false;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////																		////
////	Initialize Input/Output variables for neutrino correction			////
////																		////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	MCParticle* trueNeutrino = getTrueNeutrino( SLDLepton );
	MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
	TLorentzVector trueVisibleFourMomentumAtSLDVertex( 0.0 , 0.0 , 0.0 , 0.0 );
	for ( unsigned int i_d = 0 ; i_d < parentHadron->getDaughters().size() ; ++i_d )
	{
		if ( parentHadron->getDaughters()[ i_d ] != trueNeutrino ) trueVisibleFourMomentumAtSLDVertex += TLorentzVector( parentHadron->getDaughters()[ i_d ]->getMomentum() , parentHadron->getDaughters()[ i_d ]->getEnergy() );
	}

	TLorentzVector trueNeutrinoFourMomentum( trueNeutrino->getMomentum() , trueNeutrino->getEnergy() );
	TLorentzVector trueHadronFourMomentum( parentHadron->getMomentumAtEndpoint() , parentHadron->getEnergy() );
	TLorentzVector trueLeptonFourMomentum( SLDLepton->getMomentum() , SLDLepton->getEnergy() );
	TLorentzVector truePVATrueFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector truePVARecoFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoPVARecoFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );

	TLorentzVector recoLeptonFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector truePVATrueChargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector truePVARecoChargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoSLDVertexChargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoPVARecoChargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector truePVATrueNeutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector truePVARecoNeutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoPVARecoNeutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );

	TLorentzVector recoNeutrinoFourMomentumPos( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumNeg( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumClose( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoHadronFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );

	std::vector< float > NeutrinoCovMat( 10 , 0.0 );
	std::vector< float > NeutrinoCovMatPos( 10 , 0.0 );
	std::vector< float > NeutrinoCovMatNeg( 10 , 0.0 );
	std::vector< float > NeutrinoCovMatZero( 10 , 0.0 );
	std::vector< float > CovMatrixDetector( 10 , 0.0 );
	std::vector< float > CovMatrixPVA( 10, 0.0 );
	std::vector< float > CovMatrixFlightDirection( 6, 0.0 );

	TLorentzVector leptonFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector chargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector neutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector visibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );

	std::vector<float> truePrimaryVertex{};
	truePrimaryVertex.push_back( parentHadron->getVertex()[ 0 ] );
	truePrimaryVertex.push_back( parentHadron->getVertex()[ 1 ] );
	truePrimaryVertex.push_back( parentHadron->getVertex()[ 2 ] );

	TVector3 trueFlightDirection( 0.0 , 0.0 , 0.0 );
	TVector3 recoFlightDirection( 0.0 , 0.0 , 0.0 );
	TVector3 flightDirection( 0.0 , 0.0 , 0.0 );
	double parentHadronMass = parentHadron->getMass();// Cheated for the time being
	float helicesDistance = 0.0;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////																	    ////
////	Get True Visible decay products, four-momentum & flight direction   ////
////																	    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	float weightPFOtoMCP = 0.0;
	float weightMCPtoPFO = 0.0;
	ReconstructedParticle* linkedRecoLepton = NULL;
	pfoVector tempTruePVAChargedDecayProducts{};
	pfoVector tempTruePVANeutralDecayProducts{};
	getTrueDecayProducts( parentHadron , SLDLepton , trueNeutrino , trueNeutralDecayProducts , trueChargedDecayProducts );
	getTruePVADecayProducts( parentHadron , SLDLepton , trueNeutrino , linkedRecoLepton , weightPFOtoMCP , weightMCPtoPFO , tempTruePVANeutralDecayProducts , tempTruePVAChargedDecayProducts , RecoMCParticleNav , MCParticleRecoNav );
	for ( unsigned int i_par = 0 ; i_par < tempTruePVANeutralDecayProducts.size() ; ++i_par )
	{
		if ( ( tempTruePVANeutralDecayProducts[ i_par ]->getTracks() ).size() == 0 )
		{
			truePVANeutralDecayProducts.push_back( tempTruePVANeutralDecayProducts[ i_par ] );
		}
		else
		{
			truePVAChargedDecayProducts.push_back( tempTruePVANeutralDecayProducts[ i_par ] );
		}
	}
	for ( unsigned int i_par = 0 ; i_par < tempTruePVAChargedDecayProducts.size() ; ++i_par )
	{
		if ( ( tempTruePVAChargedDecayProducts[ i_par ]->getTracks() ).size() == 0 )
		{
			truePVANeutralDecayProducts.push_back( tempTruePVAChargedDecayProducts[ i_par ] );
		}
		else
		{
			truePVAChargedDecayProducts.push_back( tempTruePVAChargedDecayProducts[ i_par ] );
		}
	}
	for ( unsigned int i_par = 0 ; i_par < truePVAChargedDecayProducts.size() ; ++i_par )
	{
		if ( isParticleInVertex( truePVAChargedDecayProducts[ i_par ] , primaryVertex ) )
		{
			++nChargedParticlesInPrimVertex;
		}
		else
		{
			++nChargedParticlesNotInPrimVertex;
		}
	}

	streamlog_out(DEBUG5) << "		True Charged Lepton of Semi-Leptonic Decay:" << std::endl;
	streamlog_out(DEBUG5) << *SLDLepton << std::endl;
	streamlog_out(DEBUG5) << "		Reco Charged Lepton of Semi-Leptonic Decay:" << std::endl;
	if ( linkedRecoLepton == NULL )
	{
		streamlog_out(DEBUG5) << "There is no reconstructed particle linked to the true charged Lepton" << std::endl;
	}
	else
	{
		streamlog_out(DEBUG5) << *linkedRecoLepton << std::endl;
	}


	streamlog_out(DEBUG5) << "" << std::endl;
	streamlog_out(DEBUG5) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG5) << "----------------- " << trueNeutralDecayProducts.size() << " Stable Neutral MCParticles Found------------------" << std::endl;
	streamlog_out(DEBUG5) << "----------------------------------------------------------------------" << std::endl;
	for ( unsigned int i_par = 0 ; i_par < trueNeutralDecayProducts.size() ; ++i_par )
	{
		truePVATrueNeutralFourMomentum += TLorentzVector( ( trueNeutralDecayProducts[ i_par ] )->getMomentum() , ( trueNeutralDecayProducts[ i_par ] )->getEnergy() );
		streamlog_out(DEBUG5) << "		Found One Stable Neutral Particle in Semi-Leptonic Decay Products:" << std::endl;
		streamlog_out(DEBUG5) << *trueNeutralDecayProducts[ i_par ] << std::endl;
	}
	streamlog_out(DEBUG5) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG5) << "--------------- " << truePVANeutralDecayProducts.size() << " Neutral Reconstructed Particle Found----------------" << std::endl;
	streamlog_out(DEBUG5) << "----------------------------------------------------------------------" << std::endl;
	for ( unsigned int i_par = 0 ; i_par < truePVANeutralDecayProducts.size() ; ++i_par )
	{
		truePVARecoNeutralFourMomentum += TLorentzVector( ( truePVANeutralDecayProducts[ i_par ] )->getMomentum() , ( truePVANeutralDecayProducts[ i_par ] )->getEnergy() );
		streamlog_out(DEBUG5) << "		Found One Neutral Reconstructed Particle in Semi-Leptonic Decay Products:" << std::endl;
		streamlog_out(DEBUG5) << *truePVANeutralDecayProducts[ i_par ] << std::endl;
	}
	streamlog_out(DEBUG5) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG5) << "----------------- " << trueChargedDecayProducts.size() << " Stable Charged MCParticles Found------------------" << std::endl;
	streamlog_out(DEBUG5) << "----------------------------------------------------------------------" << std::endl;
	for ( unsigned int i_par = 0 ; i_par < trueChargedDecayProducts.size() ; ++i_par )
	{
		truePVATrueChargedFourMomentum += TLorentzVector( ( trueChargedDecayProducts[ i_par ] )->getMomentum() , ( trueChargedDecayProducts[ i_par ] )->getEnergy() );
		streamlog_out(DEBUG5) << "		Found One Stable Charged Particle in Semi-Leptonic Decay Products:" << std::endl;
		streamlog_out(DEBUG5) << *trueChargedDecayProducts[ i_par ] << std::endl;
	}
	streamlog_out(DEBUG5) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG5) << "--------------- " << truePVAChargedDecayProducts.size() << " Charged Reconstructed Particle Found----------------" << std::endl;
	streamlog_out(DEBUG5) << "----------------------------------------------------------------------" << std::endl;
	for ( unsigned int i_par = 0 ; i_par < truePVAChargedDecayProducts.size() ; ++i_par )
	{
		truePVARecoChargedFourMomentum += TLorentzVector( ( truePVAChargedDecayProducts[ i_par ] )->getMomentum() , ( truePVAChargedDecayProducts[ i_par ] )->getEnergy() );
		streamlog_out(DEBUG5) << "		Found One Charged Reconstructed Particle in Semi-Leptonic Decay Products:" << std::endl;
		streamlog_out(DEBUG5) << *truePVAChargedDecayProducts[ i_par ] << std::endl;
	}
	truePVATrueFourMomentum = truePVATrueNeutralFourMomentum + truePVATrueChargedFourMomentum + trueLeptonFourMomentum;
	truePVARecoFourMomentum = truePVARecoNeutralFourMomentum + truePVARecoChargedFourMomentum + recoLeptonFourMomentum;
	streamlog_out(DEBUG5) << "" << std::endl;
	streamlog_out(DEBUG5) << "		True Visible Energy At the vertex of Semi-Leptonic Decay, E_vis 	= " << trueVisibleFourMomentumAtSLDVertex.E() << " GeV" << std::endl;
	streamlog_out(DEBUG5) << "		True Visible Energy from decay chain of Semi-Leptonic Decay, E_vis 	= " << truePVATrueFourMomentum.E() << " GeV" << std::endl;



	std::vector<double> trueStartVertex{};
	std::vector<double> trueSLDVertex{};
	getTrueFlightDirection( SLDLepton , trueFlightDirection , trueStartVertex , trueSLDVertex );

	m_true_E_vis.push_back( trueVisibleFourMomentumAtSLDVertex.E() );
	m_true_E_vis_prime.push_back( ( pow( parentHadron->getMass() , 2 ) + pow( trueVisibleFourMomentumAtSLDVertex.M() , 2 ) ) / ( 2 * parentHadron->getMass() ) );
	TVector3 p_vis = TVector3( trueVisibleFourMomentumAtSLDVertex.Px() , trueVisibleFourMomentumAtSLDVertex.Py() , trueVisibleFourMomentumAtSLDVertex.Pz() );
	m_true_P_vis_par.push_back( p_vis.Dot( trueFlightDirection ) );
	m_true_P_vis_nor.push_back( sqrt( pow( p_vis.Mag() , 2 ) - pow( p_vis.Dot( trueFlightDirection ) , 2 ) ) );
	m_true_P_vis_par_prime.push_back( sqrt( pow( ( pow( parentHadron->getMass() , 2 ) - pow( trueVisibleFourMomentumAtSLDVertex.M() , 2 ) ) / ( 2 * parentHadron->getMass() ) , 2 ) - pow( sqrt( pow( p_vis.Mag() , 2 ) - pow( p_vis.Dot( trueFlightDirection ) , 2 ) ) , 2 ) ) );
	streamlog_out(DEBUG8) << "		 Flight Dir of Parent Hadron:		( " << trueFlightDirection.X() << " , " << trueFlightDirection.Y() << " , " << trueFlightDirection.Z() <<" )" << std::endl;
	double _parentHadronMass = parentHadron->getMass();
	streamlog_out(DEBUG8) << "		 Parent Hadron Mass:			" << _parentHadronMass << std::endl;
	double _parentHadronEnergy = parentHadron->getEnergy();
	streamlog_out(DEBUG8) << "		 Parent Hadron Energy:			" << _parentHadronEnergy << std::endl;
	double _visibleEnergy = trueVisibleFourMomentumAtSLDVertex.E();
	streamlog_out(DEBUG8) << "		 Visible Energy:				" << _visibleEnergy << std::endl;
	double _visibleEnergyPrime = ( pow( parentHadron->getMass() , 2 ) + pow( trueVisibleFourMomentumAtSLDVertex.M() , 2 ) ) / ( 2 * parentHadron->getMass() );
	streamlog_out(DEBUG8) << "		 Visible Energy prime:			" << _visibleEnergyPrime << std::endl;
	double _visibleMomentum = p_vis.Mag();
	streamlog_out(DEBUG8) << "		 Visible Momentum:				" << _visibleMomentum << " (" << p_vis.Px() << " , " << p_vis.Py() << " , " << p_vis.Pz() << ")" << std::endl;
	double _visibleMomentumParallel = p_vis.Dot( trueFlightDirection );
	streamlog_out(DEBUG8) << "		|Visible Momentum (par)|:		" << _visibleMomentumParallel << std::endl;
	double _visibleMomentumPerpendicular = sqrt( pow( p_vis.Mag() , 2 ) - pow( p_vis.Dot( trueFlightDirection ) , 2 ) );
	streamlog_out(DEBUG8) << "		|Visible Momentum (nor)|:		" << _visibleMomentumPerpendicular << std::endl;
	double _visibleMomentumParallelPrime = sqrt( pow( ( pow( parentHadron->getMass() , 2 ) - pow( trueVisibleFourMomentumAtSLDVertex.M() , 2 ) ) / ( 2 * parentHadron->getMass() ) , 2 ) - pow( sqrt( pow( p_vis.Mag() , 2 ) - pow( p_vis.Dot( trueFlightDirection ) , 2 ) ) , 2 ) );
	streamlog_out(DEBUG8) << "		|Visible Momentum (par-prime)|:	+/-" << _visibleMomentumParallelPrime << std::endl;
	double _cheatedNeutrinoSolutionPos = ( ( _visibleEnergy * _visibleEnergyPrime ) - ( +1.0 * _visibleMomentumParallel * _visibleMomentumParallelPrime ) ) * _parentHadronMass / ( pow( _visibleEnergy , 2 ) - pow( _visibleMomentumParallel , 2 ) ) - _visibleEnergy;
	double _cheatedNeutrinoSolutionNeg = ( ( _visibleEnergy * _visibleEnergyPrime ) - ( -1.0 * _visibleMomentumParallel * _visibleMomentumParallelPrime ) ) * _parentHadronMass / ( pow( _visibleEnergy , 2 ) - pow( _visibleMomentumParallel , 2 ) ) - _visibleEnergy;
	double _trueNeutrinoEnergy = trueNeutrinoFourMomentum.E();
	streamlog_out(DEBUG8) << "		 Neutrino Solution (+):			" << _cheatedNeutrinoSolutionPos << std::endl;
	streamlog_out(DEBUG8) << "		 Neutrino Solution (-):			" << _cheatedNeutrinoSolutionNeg << std::endl;
	streamlog_out(DEBUG8) << "		 True Neutrino Energy :			" << _trueNeutrinoEnergy << std::endl;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////																	    ////
////	Identify Reconstrcuted Lepton and the jet/vertex assigned to SLD	////
////																		////
////////////////////////////////////////////////////////////////////////////////
////				If there is no Reconstrcuted Lepton						////
////								OR										////
////				Reconstrcuted Lepton is not in a jet					////
////								OR										////
////		     Reconstrcuted Lepton is in Primary Vertex					////
////																		////
////				     Neutrino Correction FAILS							////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	int SLDStatus = -999;
	if ( linkedRecoLepton == NULL )
	{
		SLDStatus = 1;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||| Reconstructed Lepton is not found ||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		if ( m_fillRootTree ) h_SLDStatus->Fill( 0.5 );
		return;
	}
	recoLeptonFourMomentum = TLorentzVector( linkedRecoLepton->getMomentum() , linkedRecoLepton->getEnergy() );

	bool recoLeptonIsInJet = false;
	ReconstructedParticle *assignedJet = NULL;
	assignedJet = getJetAssignedToParticle( linkedRecoLepton , jetVector , recoLeptonIsInJet );
	if ( !recoLeptonIsInJet )
	{
		SLDStatus = 2;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||| Reconstructed Lepton doesn't belong to any jet |||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		if ( m_fillRootTree ) h_SLDStatus->Fill( 1.5 );
		return;
	}

	if ( isParticleInVertex( linkedRecoLepton , primaryVertex ) )
	{
		SLDStatus = 3;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||| Reconstructed Lepton is in primary vertex ||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		if ( m_fillRootTree ) h_SLDStatus->Fill( 2.5 );
		return;
	}
	investigateJetEnergyContent( assignedJet );
	jetsOfSemiLeptonicDecays.push_back( assignedJet );

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////																		////
////	   Investigate Alone charged Particles in the assigned jet			////
////																		////
////			     FromSLD and NotFromSLD is Cheated						////
////																		////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	aloneChargedPFOsInJet = getParticlesWithAloneTracks( linkedRecoLepton , assignedJet , primaryVertex , buildUpVertexVector );

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////																	    ////
////	   Reconstruction of Flight Direction of Parent Hadron			    ////
////																	    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	TVector3 leptonDirection = TVector3( SLDLepton->getMomentum() ); leptonDirection.SetMag( 1.0 );
	TVector3 jetAxis = TVector3( assignedJet->getMomentum() ); jetAxis.SetMag( 1.0 );
	if ( m_displayEvent ) ced_line( primaryVertex->getPosition()[ 0 ] + 10.0 * jetAxis.X() , primaryVertex->getPosition()[ 1 ] + 10.0 * jetAxis.Y() , primaryVertex->getPosition()[ 2 ] + 10.0 * jetAxis.Z() , primaryVertex->getPosition()[ 0 ] , primaryVertex->getPosition()[ 1 ] , primaryVertex->getPosition()[ 2 ] , 2 , 1 , 0x000000 );

	vtxVector verticesInJet = getVerticesInJet( assignedJet , buildUpVertexVector );
	m_nRecoVerticesInJet.push_back( verticesInJet.size() );

	int vertexingScenario = m_vertexingScenario;
	pfoVector sortedChargedPFOs;
	if ( aloneChargedPFOsInJet.size() != 0 ) sortParticles( sortedChargedPFOs , aloneChargedPFOsInJet , jetAxis );
	if ( sortedChargedPFOs.size() == 0 && m_vertexingScenario == 4 ) vertexingScenario = 1;
	m_flightDirectionStatus.push_back( vertexingScenario );
	double hadronFlightLength;
	TVector3 daughterHadronFlightDirection;
	double daughterHadronFlightDistance = 0.0;
	std::vector<double> sldVertexPosition{};
	TVector3 PCAatLepton( 0.0 , 0.0 , 0.0 );
	TVector3 PCAatOtherParticle( 0.0 , 0.0 , 0.0 );
	double lepton3DImpactParameter = 0.0;
	double OtherParticle3DImpactParameter = 0.0;
	if ( SLDStatus < 0 ) SLDStatus = getRecoFlightDirection( linkedRecoLepton , recoFlightDirection , hadronFlightLength , primaryVertex , startVertex , SLDVertex , SLDVertexRP , assignedJet , verticesInJet , sortedChargedPFOs , helicesDistance , vertexingScenario , daughterHadronFlightDirection , daughterHadronFlightDistance , sldVertexPosition , PCAatLepton , PCAatOtherParticle , lepton3DImpactParameter , OtherParticle3DImpactParameter );
	if ( flightDirection.Dot( jetAxis ) < 0.0 ) flightDirection = -1.0 * flightDirection;
	m_daughterHadronFlightDistance.push_back( daughterHadronFlightDistance );
	m_recoParentHadronFlightDistance.push_back( hadronFlightLength );
	m_distRecoLeptonToDownStreamVertex.push_back( helicesDistance );
	m_FlightDirectionErrorCosAlpha.push_back( trueFlightDirection.Dot( recoFlightDirection ) );
	m_FlightDirectionErrorSinAlpha.push_back( std::sin( acos( trueFlightDirection.Dot( recoFlightDirection ) ) ) );
	m_FlightDirectionErrorAlpha.push_back( acos( trueFlightDirection.Dot( recoFlightDirection ) ) );
	m_Lepton3DImpactParameter.push_back( lepton3DImpactParameter );
	m_OtherParticle3DImpactParameter.push_back( OtherParticle3DImpactParameter );
	m_nChargedParticlesInPrimVertex.push_back( nChargedParticlesInPrimVertex );
	m_nChargedParticlesNotInPrimVertex.push_back( nChargedParticlesNotInPrimVertex );
	m_weightPFOtoMCP_Lepton.push_back( weightPFOtoMCP );
	m_weightMCPtoPFO_Lepton.push_back( weightMCPtoPFO );
	m_PCAatLeptonX.push_back( PCAatLepton.X() );
	m_PCAatLeptonY.push_back( PCAatLepton.Y() );
	m_PCAatLeptonZ.push_back( PCAatLepton.Z() );
	m_PCAatOtherParticleX.push_back( PCAatOtherParticle.X() );
	m_PCAatOtherParticleY.push_back( PCAatOtherParticle.Y() );
	m_PCAatOtherParticleZ.push_back( PCAatOtherParticle.Z() );
	m_JetAxisX.push_back( jetAxis.X() );
	m_JetAxisY.push_back( jetAxis.Y() );
	m_JetAxisZ.push_back( jetAxis.Z() );
	m_PrimaryVertexX.push_back( primaryVertex->getPosition()[ 0 ] );
	m_PrimaryVertexY.push_back( primaryVertex->getPosition()[ 1 ] );
	m_PrimaryVertexZ.push_back( primaryVertex->getPosition()[ 2 ] );
	TVector3 LepFromPV( PCAatLepton.X() - primaryVertex->getPosition()[ 0 ] , PCAatLepton.Y() - primaryVertex->getPosition()[ 1 ] , PCAatLepton.Z() - primaryVertex->getPosition()[ 2 ] );
	m_LeptonDistanceFromPV.push_back( LepFromPV.Mag() );
	LepFromPV.SetMag( 1.0 );
	TVector3 DSVFromPV( PCAatOtherParticle.X() - primaryVertex->getPosition()[ 0 ] , PCAatOtherParticle.Y() - primaryVertex->getPosition()[ 1 ] , PCAatOtherParticle.Z() - primaryVertex->getPosition()[ 2 ] );
	m_DSVDistanceFromPV.push_back( DSVFromPV.Mag() );
	DSVFromPV.SetMag( 1.0 );
	m_AngleLeptonJet.push_back( acos( LepFromPV.Dot( jetAxis ) ) );
	m_AngleDSVertexJet.push_back( acos( DSVFromPV.Dot( jetAxis ) ) );
	m_AngleTrueFlightDirectionJet.push_back( acos( trueFlightDirection.Dot( jetAxis ) ) );
	m_AngleRecoFlightDirectionJet.push_back( acos( recoFlightDirection.Dot( jetAxis ) ) );
	if ( m_cheatFlightDirection )
	{
		flightDirection = trueFlightDirection;
	}
	else
	{
		flightDirection = recoFlightDirection;
		//flightDirection = jetAxis;
	}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////																				////
////		Particle Association to the vertex of the semi-leptonic decay			////
////																				////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

	float chargedCosAcceptanceAngle = 0.0;
	TVector3 PVAConeAxis( 0.0 , 0.0 , 0.0 );
	if ( SLDStatus == 4 )
	{
		chargedCosAcceptanceAngle = m_chargedCosAcceptanceAngleSLD4;
		//PVAConeAxis = daughterHadronFlightDirection;
	}
	else if ( SLDStatus == 5 )
	{
		chargedCosAcceptanceAngle = m_chargedCosAcceptanceAngleSLD5;
		//PVAConeAxis = daughterHadronFlightDirection;
	}
	PVAConeAxis = ( m_cheatFlightDirection ? trueFlightDirection : recoFlightDirection );
	//PVAConeAxis = daughterHadronFlightDirection;
	PVAConeAxis.SetMag( 1.0 );

	for ( unsigned int i_par = 0 ; i_par < assignedJet->getParticles().size() ; ++i_par )
	{
		ReconstructedParticle* jetParticle = assignedJet->getParticles()[ i_par ];
		if ( jetParticle == linkedRecoLepton ) continue;
		TVector3 jetParticleMomentum = TVector3( jetParticle->getMomentum() );
		jetParticleMomentum.SetMag( 1.0 );
		if ( jetParticle->getTracks().size() == 0 && jetParticleMomentum.Dot( PVAConeAxis ) >= m_neutralCosAcceptanceAngle )
		{
			streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG2) << "-------- Added One Neutral PFO to SLDecay products candidates --------" << std::endl;
			streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG2) << *jetParticle << std::endl;
			allPFOsInJet.push_back( jetParticle );
			neutralPFOsInJet.push_back( jetParticle );
		}
		else if ( jetParticle->getTracks().size() != 0 && jetParticleMomentum.Dot( PVAConeAxis ) >= chargedCosAcceptanceAngle )
		{
			bool particleIsInAVertex = false;
			for ( unsigned int i_vtx = 0 ; i_vtx < buildUpVertexVector.size() ; ++i_vtx )
			{
				if ( !particleIsInAVertex )
				{
					Vertex *testVertex = buildUpVertexVector[ i_vtx ];
					streamlog_out(DEBUG0) << "	Looking for particle in build up vertex[ " << i_vtx << " ]: " << testVertex << std::endl;
					streamlog_out(DEBUG0) << *testVertex << std::endl;
					particleIsInAVertex = isParticleInVertex( jetParticle , testVertex );
				}
			}
			if ( !particleIsInAVertex )
			{
				streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
				streamlog_out(DEBUG2) << "-------- Added One Charged PFO to SLDecay products candidates --------" << std::endl;
				streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
				streamlog_out(DEBUG2) << *jetParticle << std::endl;
				allPFOsInJet.push_back( jetParticle );
				chargedPFOsInJet.push_back( jetParticle );
			}
		}
	}
	double InvMassCutCharged = 0.0;
	double InvMassCutNeutral = 0.0;
	TLorentzVector visibleTLV( 0.0 , 0.0 , 0.0 , 0.0 );
	if ( SLDStatus == 4 )
	{
		InvMassCutCharged = m_BSLDChargedSLD4InvMassCut;
		InvMassCutNeutral = m_BSLDNeutralSLD4InvMassCut;
		//ReconstructedParticle* sldVertexRP = SLDVerticesRP[ 0 ];
		for ( unsigned int i_par = 0 ; i_par < SLDVertexRP->getParticles().size() ; ++i_par )
		{
			ReconstructedParticle* chargedDecayProduct = SLDVertexRP->getParticles()[ i_par ];
			if ( chargedDecayProduct != linkedRecoLepton ) recoSLDVertexDecayProducts.push_back( chargedDecayProduct );
			if ( m_displayEvent ) drawReconstructedParticle( chargedDecayProduct , primaryVertex , 0x0075df , 0x2e8e04 );
		}
	}
	else if ( SLDStatus == 5 )
	{
		InvMassCutCharged = m_BSLDChargedSLD5InvMassCut;
		InvMassCutNeutral = m_BSLDNeutralSLD5InvMassCut;
		//ReconstructedParticle* sldVertexRP = SLDVerticesRP[ 0 ];
		for ( unsigned int i_par = 0 ; i_par < SLDVertexRP->getParticles().size() ; ++i_par )
		{
			ReconstructedParticle* chargedDecayProduct = SLDVertexRP->getParticles()[ i_par ];
			if ( chargedDecayProduct != linkedRecoLepton ) recoSLDVertexDecayProducts.push_back( chargedDecayProduct );
			if ( m_displayEvent ) drawReconstructedParticle( chargedDecayProduct , primaryVertex , 0x0075df , 0x2e8e04 );
		}
	}
	for ( unsigned int i_par = 0 ; i_par < recoSLDVertexDecayProducts.size() ; ++i_par )
	{
		visibleTLV += TLorentzVector( recoSLDVertexDecayProducts[ i_par ]->getMomentum() , recoSLDVertexDecayProducts[ i_par ]->getEnergy() );
	}

	vtxVector availableVerticesInJet{};
	for ( unsigned int i_vtx = 0 ; i_vtx < verticesInJet.size() ; ++i_vtx )
	{
		if ( verticesInJet[ i_vtx ] != SLDVertex ) availableVerticesInJet.push_back( verticesInJet[ i_vtx ] );
	}
	if ( m_cheatPVAneutral )
	{
		for ( unsigned int i_par = 0 ; i_par < truePVANeutralDecayProducts.size() ; ++i_par )
		{
			recoPVANeutralDecayProducts.push_back( truePVANeutralDecayProducts[ i_par ] );
			if ( m_displayEvent ) drawReconstructedParticle( truePVANeutralDecayProducts[ i_par ] , primaryVertex , 0x0075df , 0x2e8e04 );
		}
		if ( m_cheatPVAcharged )
		{
			for ( unsigned int i_par = 0 ; i_par < truePVAChargedDecayProducts.size() ; ++i_par )
			{
				recoPVAChargedDecayProducts.push_back( truePVAChargedDecayProducts[ i_par ] );
				if ( m_displayEvent ) drawReconstructedParticle( truePVAChargedDecayProducts[ i_par ] , primaryVertex , 0x0075df , 0x2e8e04 );
			}
		}
		else
		{
			m_visiblePxBeforePVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Px() );
			m_visiblePyBeforePVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Py() );
			m_visiblePzBeforePVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Pz() );
			m_visibleEnergyBeforePVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).E() );
			m_visibleMassBeforePVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).M() );
			assignVerticesToSemiLeptonicDecay( recoPVAChargedDecayProducts , availableVerticesInJet , InvMassCutCharged , PVAConeAxis , startVertex , visibleTLV );
			m_visiblePxAfterVertexPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Px() );
			m_visiblePyAfterVertexPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Py() );
			m_visiblePzAfterVertexPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Pz() );
			m_visibleEnergyAfterVertexPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).E() );
			m_visibleMassAfterVertexPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).M() );
			assignParticlesToSemiLeptonicDecay( recoPVAChargedDecayProducts , chargedPFOsInJet , InvMassCutCharged , PVAConeAxis , visibleTLV );
			m_visiblePxAfterChargedPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Px() );
			m_visiblePyAfterChargedPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Py() );
			m_visiblePzAfterChargedPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Pz() );
			m_visibleEnergyAfterChargedPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).E() );
			m_visibleMassAfterChargedPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).M() );
		}
	}
	else
	{
		if ( m_cheatPVAcharged )
		{
			for ( unsigned int i_par = 0 ; i_par < truePVAChargedDecayProducts.size() ; ++i_par )
			{
				recoPVAChargedDecayProducts.push_back( truePVAChargedDecayProducts[ i_par ] );
				visibleTLV += TLorentzVector( truePVAChargedDecayProducts[ i_par ]->getMomentum() , truePVAChargedDecayProducts[ i_par ]->getEnergy() );
				if ( m_displayEvent ) drawReconstructedParticle( truePVAChargedDecayProducts[ i_par ] , primaryVertex , 0x0075df , 0x2e8e04 );
			}
		}
		else
		{
			m_visiblePxBeforePVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Px() );
			m_visiblePyBeforePVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Py() );
			m_visiblePzBeforePVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Pz() );
			m_visibleEnergyBeforePVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).E() );
			m_visibleMassBeforePVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).M() );
			assignVerticesToSemiLeptonicDecay( recoPVAChargedDecayProducts , availableVerticesInJet , InvMassCutCharged , PVAConeAxis , startVertex , visibleTLV );
			m_visiblePxAfterVertexPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Px() );
			m_visiblePyAfterVertexPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Py() );
			m_visiblePzAfterVertexPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Pz() );
			m_visibleEnergyAfterVertexPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).E() );
			m_visibleMassAfterVertexPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).M() );
			assignParticlesToSemiLeptonicDecay( recoPVAChargedDecayProducts , chargedPFOsInJet , InvMassCutCharged , PVAConeAxis , visibleTLV );
			m_visiblePxAfterChargedPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Px() );
			m_visiblePyAfterChargedPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Py() );
			m_visiblePzAfterChargedPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Pz() );
			m_visibleEnergyAfterChargedPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).E() );
			m_visibleMassAfterChargedPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).M() );
		}
		assignParticlesToSemiLeptonicDecay( recoPVANeutralDecayProducts , neutralPFOsInJet , InvMassCutNeutral , PVAConeAxis , visibleTLV );
		m_visiblePxAfterNeutralPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Px() );
		m_visiblePyAfterNeutralPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Py() );
		m_visiblePzAfterNeutralPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).Pz() );
		m_visibleEnergyAfterNeutralPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).E() );
		m_visibleMassAfterNeutralPVA.push_back( ( visibleTLV + recoLeptonFourMomentum ).M() );
	}
	if ( m_displayEvent ) drawReconstructedParticle( linkedRecoLepton , primaryVertex , 0xfc0000 , 0x2e8e04 );

	m_visibleChargedInvMassCut.push_back( InvMassCutCharged );
	m_visibleNeutralInvMassCut.push_back( InvMassCutNeutral );

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////																		////
////		Form visible four-momentum for neutrino correction				////
////																		////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	recoPVADecayProducts.push_back( linkedRecoLepton );
	for ( unsigned int i_par = 0 ; i_par < recoSLDVertexDecayProducts.size() ; ++i_par )
	{
		recoSLDVertexChargedFourMomentum += TLorentzVector( recoSLDVertexDecayProducts[ i_par ]->getMomentum() , recoSLDVertexDecayProducts[ i_par ]->getEnergy() );
		recoPVADecayProducts.push_back( recoSLDVertexDecayProducts[ i_par ] );
	}
	for ( unsigned int i_par = 0 ; i_par < recoPVAVertexDecayProducts.size() ; ++i_par )
	{
		recoPVARecoChargedFourMomentum += TLorentzVector( recoPVAVertexDecayProducts[ i_par ]->getMomentum() , recoPVAVertexDecayProducts[ i_par ]->getEnergy() );
		recoPVADecayProducts.push_back( recoPVAVertexDecayProducts[ i_par ] );
	}
	for ( unsigned int i_par = 0 ; i_par < recoPVAChargedDecayProducts.size() ; ++i_par )
	{
		recoPVARecoChargedFourMomentum += TLorentzVector( recoPVAChargedDecayProducts[ i_par ]->getMomentum() , recoPVAChargedDecayProducts[ i_par ]->getEnergy() );
		recoPVADecayProducts.push_back( recoPVAChargedDecayProducts[ i_par ] );
	}
	for ( unsigned int i_par = 0 ; i_par < recoPVANeutralDecayProducts.size() ; ++i_par )
	{
		recoPVARecoNeutralFourMomentum += TLorentzVector( recoPVANeutralDecayProducts[ i_par ]->getMomentum() , recoPVANeutralDecayProducts[ i_par ]->getEnergy() );
		recoPVADecayProducts.push_back( recoPVANeutralDecayProducts[ i_par ] );
	}
	recoPVARecoFourMomentum = recoLeptonFourMomentum + recoSLDVertexChargedFourMomentum + recoPVARecoChargedFourMomentum + recoPVARecoNeutralFourMomentum;

	if ( m_cheatLepton4momentum )
	{
		leptonFourMomentum = trueLeptonFourMomentum;
	}
	else
	{
		leptonFourMomentum = recoLeptonFourMomentum;
	}
	if ( m_cheatPVAneutral )
	{
		if ( m_cheatNeutral4momentum )
		{
			neutralFourMomentum = truePVATrueNeutralFourMomentum;
		}
		else
		{
			neutralFourMomentum = truePVARecoNeutralFourMomentum;
		}
	}
	else
	{
		neutralFourMomentum = recoPVARecoNeutralFourMomentum;
	}
	if ( m_cheatPVAcharged )
	{
		if ( m_cheatCharged4momentum )
		{
			chargedFourMomentum = truePVATrueChargedFourMomentum;
		}
		else
		{
			chargedFourMomentum = truePVARecoChargedFourMomentum;
		}
	}
	else
	{
		chargedFourMomentum = recoSLDVertexChargedFourMomentum + recoPVARecoChargedFourMomentum;
	}
	visibleFourMomentum = leptonFourMomentum + neutralFourMomentum + chargedFourMomentum;

	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		||||||||||||||| Input to SLDecay Correction |||||||||||||||" << std::endl;
	streamlog_out(DEBUG8) << "		SLD Vertex Status: " << SLDStatus << std::endl;
	streamlog_out(DEBUG8) << "			     (  			PDG 	, Mass		, Px		, Py		, Pz		, E		, Charge	)" << std::endl;
	streamlog_out(DEBUG8) << "		Neutrino" << std::endl;
	streamlog_out(DEBUG8) << "			True:(				" << trueNeutrino->getPDG() << "	, " << trueNeutrino->getMass() << "	, " << trueNeutrino->getMomentum()[ 0 ] << "	, " << trueNeutrino->getMomentum()[ 1 ] << "	, " << trueNeutrino->getMomentum()[ 2 ] << "	, " << trueNeutrino->getEnergy() << "	, " << trueNeutrino->getCharge() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		Hadron" << std::endl;
	streamlog_out(DEBUG8) << "			True:(				" << parentHadron->getPDG() << "	, " << parentHadron->getMass() << "	, " << parentHadron->getMomentum()[ 0 ] << "	, " << parentHadron->getMomentum()[ 1 ] << "	, " << parentHadron->getMomentum()[ 2 ] << "	, " << parentHadron->getEnergy() << "	, " << parentHadron->getCharge() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		Lepton" << std::endl;
	streamlog_out(DEBUG8) << "			True:(				" << SLDLepton->getPDG() << "	, " << SLDLepton->getMass() << "	, " << SLDLepton->getMomentum()[ 0 ] << "	, " << SLDLepton->getMomentum()[ 1 ] << "	, " << SLDLepton->getMomentum()[ 2 ] << "	, " << SLDLepton->getEnergy() << "	, " << SLDLepton->getCharge() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "			Reco:(				" << linkedRecoLepton->getType() << "	, " << linkedRecoLepton->getMass() << "	, " << linkedRecoLepton->getMomentum()[ 0 ] << "	, " << linkedRecoLepton->getMomentum()[ 1 ] << "	, " << linkedRecoLepton->getMomentum()[ 2 ] << "	, " << linkedRecoLepton->getEnergy() << "	, " << linkedRecoLepton->getCharge() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "			Used:(				" << "xxx" << "	, " << leptonFourMomentum.M() << "	, " << leptonFourMomentum.Px() << "	, " << leptonFourMomentum.Py() << "	, " << leptonFourMomentum.Pz() << "	, " << leptonFourMomentum.E() << "	, " << (int)SLDLepton->getCharge() << "		)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		Charged" << std::endl;
	streamlog_out(DEBUG8) << "			TruePVATrue4Mom:(		" << "qqq" << "	, " << truePVATrueChargedFourMomentum.M() << "	, " << truePVATrueChargedFourMomentum.Px() << "	, " << truePVATrueChargedFourMomentum.Py() << "	, " << truePVATrueChargedFourMomentum.Pz() << "	, " << truePVATrueChargedFourMomentum.E() << "	, " << "qqq" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			TruePVAReco4Mom:(		" << "qqq" << "	, " << truePVARecoChargedFourMomentum.M() << "	, " << truePVARecoChargedFourMomentum.Px() << "	, " << truePVARecoChargedFourMomentum.Py() << "	, " << truePVARecoChargedFourMomentum.Pz() << "	, " << truePVARecoChargedFourMomentum.E() << "	, " << "qqq" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			RecoPVAReco4Mom:(		" << "qqq" << "	, " << ( recoPVARecoChargedFourMomentum + recoSLDVertexChargedFourMomentum ).M() << "	, " << ( recoPVARecoChargedFourMomentum + recoSLDVertexChargedFourMomentum ).Px() << "	, " << ( recoPVARecoChargedFourMomentum + recoSLDVertexChargedFourMomentum ).Py() << "	, " << ( recoPVARecoChargedFourMomentum + recoSLDVertexChargedFourMomentum ).Pz() << "	, " << ( recoPVARecoChargedFourMomentum + recoSLDVertexChargedFourMomentum ).E() << "	, " << "qqq" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			Used:(				" << "qqq" << "	, " << chargedFourMomentum.M() << "	, " << chargedFourMomentum.Px() << "	, " << chargedFourMomentum.Py() << "	, " << chargedFourMomentum.Pz() << "	, " << chargedFourMomentum.E() << "	, " << "nnn" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		Neutral" << std::endl;
	streamlog_out(DEBUG8) << "			TruePVATrue4Mom:(		" << "nnn" << "	, " << truePVATrueNeutralFourMomentum.M() << "	, " << truePVATrueNeutralFourMomentum.Px() << "	, " << truePVATrueNeutralFourMomentum.Py() << "	, " << truePVATrueNeutralFourMomentum.Pz() << "	, " << truePVATrueNeutralFourMomentum.E() << "	, " << "qqq" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			TruePVAReco4Mom:(		" << "nnn" << "	, " << truePVARecoNeutralFourMomentum.M() << "	, " << truePVARecoNeutralFourMomentum.Px() << "	, " << truePVARecoNeutralFourMomentum.Py() << "	, " << truePVARecoNeutralFourMomentum.Pz() << "	, " << truePVARecoNeutralFourMomentum.E() << "	, " << "qqq" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			RecoPVAReco4Mom:(		" << "nnn" << "	, " << recoPVARecoNeutralFourMomentum.M() << "	, " << recoPVARecoNeutralFourMomentum.Px() << "	, " << recoPVARecoNeutralFourMomentum.Py() << "	, " << recoPVARecoNeutralFourMomentum.Pz() << "	, " << recoPVARecoNeutralFourMomentum.E() << "	, " << "qqq" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			Used:(				" << "nnn" << "	, " << neutralFourMomentum.M() << "	, " << neutralFourMomentum.Px() << "	, " << neutralFourMomentum.Py() << "	, " << neutralFourMomentum.Pz() << "	, " << neutralFourMomentum.E() << "	, " << "nnn" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		Visible" << std::endl;
	streamlog_out(DEBUG8) << "			True:(				" << "nnn" << "	, " << truePVATrueFourMomentum.M() << "	, " << truePVATrueFourMomentum.Px() << "	, " << truePVATrueFourMomentum.Py() << "	, " << truePVATrueFourMomentum.Pz() << "	, " << truePVATrueFourMomentum.E() << "	, " << "nnn" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			Used:(				" << "nnn" << "	, " << visibleFourMomentum.M() << "	, " << visibleFourMomentum.Px() << "	, " << visibleFourMomentum.Py() << "	, " << visibleFourMomentum.Pz() << "	, " << visibleFourMomentum.E() << "	, " << "nnn" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		Flight Direction	(						X			, Y			, Z		)" << std::endl;
	streamlog_out(DEBUG8) << "			True:		(						" << trueFlightDirection.X() << "		, " << trueFlightDirection.Y() << "		, " << trueFlightDirection.Z() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "			Reco:		(						" << recoFlightDirection.X() << "		, " << recoFlightDirection.Y() << "		, " << recoFlightDirection.Z() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "			Used:		(						" << flightDirection.X() << "		, " << flightDirection.Y() << "		, " << flightDirection.Z() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "			Daughter:	(						" << daughterHadronFlightDirection.X() << "		, " << daughterHadronFlightDirection.Y() << "		, " << daughterHadronFlightDirection.Z() << "	)" << std::endl;


	if ( visibleFourMomentum.M() > parentHadronMass ) parentHadronMass = visibleFourMomentum.M();
	m_SLDStatus.push_back( SLDStatus );
	if ( m_fillRootTree ) h_SLDStatus->Fill( SLDStatus - 0.5 );

	recoNeutrinoFourMomentumPos = getNeutrinoFourMomentumModified( flightDirection , visibleFourMomentum , parentHadronMass , +1.0 );
	recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentumModified( flightDirection , visibleFourMomentum , parentHadronMass , -1.0 );
	//recoNeutrinoFourMomentumPos = getNeutrinoFourMomentum( flightDirection , visibleFourMomentum , parentHadronMass , +1.0 );
	//recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentum( flightDirection , visibleFourMomentum , parentHadronMass , -1.0 );
	//recoNeutrinoFourMomentumPos = getNeutrinoFourMomentumStandardMethod( flightDirection , visibleFourMomentum , parentHadronMass , +1.0 );
	//recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentumStandardMethod( flightDirection , visibleFourMomentum , parentHadronMass , -1.0 );

	TVector3 recoNeutrinoMomentumPos = recoNeutrinoFourMomentumPos.Vect();
	TVector3 recoNeutrinoMomentumNeg = recoNeutrinoFourMomentumNeg.Vect();



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////																		////
////		Estimate Errors due to Particle to Vertex Association			////
////																		////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	//std::vector< float > CovMatrixDetector( 10 , 0.0 );
	//std::vector< float > CovMatrixFlightDirection( 6, 0.0 );
	addNeutrinoCovarianceMatrix( recoNeutrinoFourMomentumPos , NeutrinoCovMatPos );
	addNeutrinoCovarianceMatrix( recoNeutrinoFourMomentumNeg , NeutrinoCovMatNeg );
	int PVAStatus = 0;
	if ( recoPVANeutralDecayProducts.size() == 0 ) //( without neutral PVA)
	{
		PVAStatus = 3;
	}
	else //( with neutral PVA)
	{
		if ( recoPVAVertexDecayProducts.size() + recoPVAChargedDecayProducts.size() == 0 ) //( without charged PVA)
		{
			PVAStatus = 2;
		}
		else //( with charged PVA)
		{
			PVAStatus = 1;
		}
	}
	if( m_cheatPVAneutral || m_cheatPVAcharged ) PVAStatus = 0;
	m_PVAStatus.push_back( PVAStatus );


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////																		////
////		Checking results and preparing solutions in LCIO format			////
////																		////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	int _trueSolutionSign;
	if ( fabs( _cheatedNeutrinoSolutionPos - _trueNeutrinoEnergy ) <= fabs( _cheatedNeutrinoSolutionNeg - _trueNeutrinoEnergy ) )
	{
		m_trueSolutionSign.push_back( 1 );
		_trueSolutionSign = 1;
	}
	else
	{
		m_trueSolutionSign.push_back( -1 );
		_trueSolutionSign = -1;
	}

	if ( _trueSolutionSign == 1 )
	{
		recoNeutrinoFourMomentumClose = recoNeutrinoFourMomentumPos;
		m_recoSolutionSign.push_back( 1 );
		for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element )
		{
			NeutrinoCovMat[ i_Element ] = NeutrinoCovMatPos[ i_Element ];
		}
	}
	else if ( _trueSolutionSign == -1 )
	{
		recoNeutrinoFourMomentumClose = recoNeutrinoFourMomentumNeg;
		m_recoSolutionSign.push_back( -1 );
		for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element )
		{
			NeutrinoCovMat[ i_Element ] = NeutrinoCovMatNeg[ i_Element ];
		}
	}

	streamlog_out(DEBUG4) << "" << std::endl;
	streamlog_out(DEBUG4) << "	Closest Neutrino 4-Momentum:			( " << recoNeutrinoFourMomentumClose.Px() << "	, " << recoNeutrinoFourMomentumClose.Py() << "	, " << recoNeutrinoFourMomentumClose.Pz() << "	, " << recoNeutrinoFourMomentumClose.E() << " )" << std::endl;
	streamlog_out(DEBUG4) << "	True Neutrino 4-Momentum:			( " << trueNeutrinoFourMomentum.Px() << "	, " << trueNeutrinoFourMomentum.Py() << "	, " << trueNeutrinoFourMomentum.Pz() << "	, " << trueNeutrinoFourMomentum.E() << " )" << std::endl;
	if ( fabs( recoNeutrinoFourMomentumClose.E() - trueNeutrinoFourMomentum.E() ) > 0.1 * trueNeutrinoFourMomentum.E() )
	{
		if ( m_traceEvent )
		{
			streamlog_out(MESSAGE) << "	Closest Neutrino 4-Momentum:			( " << recoNeutrinoFourMomentumClose.Px() << "	, " << recoNeutrinoFourMomentumClose.Py() << "	, " << recoNeutrinoFourMomentumClose.Pz() << "	, " << recoNeutrinoFourMomentumClose.E() << " )" << std::endl;
			streamlog_out(MESSAGE) << "	True Neutrino 4-Momentum:			( " << trueNeutrinoFourMomentum.Px() << "	, " << trueNeutrinoFourMomentum.Py() << "	, " << trueNeutrinoFourMomentum.Pz() << "	, " << trueNeutrinoFourMomentum.E() << " )" << std::endl;
			streamlog_out(MESSAGE) << "	!!! Big Difference between true and reco neutrino Energy : " << recoNeutrinoFourMomentumClose.E() - trueNeutrinoFourMomentum.E() << "  GeV in event: " << pLCEvent->getEventNumber() << std::endl;
			streamlog_out(MESSAGE) << "	Parent hadron that decays semi-leptonically:" << std::endl;
			streamlog_out(MESSAGE) << *parentHadron << std::endl;
			checkSLDInput( parentHadron );
		}
	}


	streamlog_out(DEBUG8) << "	CovMatNeutrino(+) :	" << NeutrinoCovMatPos[ 0 ] << std::endl;
	streamlog_out(DEBUG8) << "				" << NeutrinoCovMatPos[ 1 ] << "	,	" << NeutrinoCovMatPos[ 2 ] << std::endl;
	streamlog_out(DEBUG8) << "				" << NeutrinoCovMatPos[ 3 ] << "	,	" << NeutrinoCovMatPos[ 4 ] << "	,	" << NeutrinoCovMatPos[ 5 ] << std::endl;
	streamlog_out(DEBUG8) << "				" << NeutrinoCovMatPos[ 6 ] << "	,	" << NeutrinoCovMatPos[ 7 ] << "	,	" << NeutrinoCovMatPos[ 8 ] << "	,	" << NeutrinoCovMatPos[ 9 ] << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "	CovMatNeutrino(-) :	" << NeutrinoCovMatNeg[ 0 ] << std::endl;
	streamlog_out(DEBUG8) << "				" << NeutrinoCovMatNeg[ 1 ] << "	,	" << NeutrinoCovMatNeg[ 2 ] << std::endl;
	streamlog_out(DEBUG8) << "				" << NeutrinoCovMatNeg[ 3 ] << "	,	" << NeutrinoCovMatNeg[ 4 ] << "	,	" << NeutrinoCovMatNeg[ 5 ] << std::endl;
	streamlog_out(DEBUG8) << "				" << NeutrinoCovMatNeg[ 6 ] << "	,	" << NeutrinoCovMatNeg[ 7 ] << "	,	" << NeutrinoCovMatNeg[ 8 ] << "	,	" << NeutrinoCovMatNeg[ 9 ] << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;

	evaluateInputCovMat( truePVATrueFourMomentum , trueFlightDirection , trueNeutrinoFourMomentum , visibleFourMomentum , flightDirection , recoNeutrinoFourMomentumClose , CovMatrixDetector , CovMatrixFlightDirection , NeutrinoCovMat );
	TVector3 recoNeutrinoDirection( recoNeutrinoFourMomentumClose.Px() , recoNeutrinoFourMomentumClose.Py() , recoNeutrinoFourMomentumClose.Pz() ); recoNeutrinoDirection.SetMag( 1.0 );
	m_AngleRecoNeutrinoJet.push_back( acos( recoNeutrinoDirection.Dot( jetAxis ) ) );
	TVector3 trueNeutrinoDirection( trueNeutrinoFourMomentum.Px() , trueNeutrinoFourMomentum.Py() , trueNeutrinoFourMomentum.Pz() ); trueNeutrinoDirection.SetMag( 1.0 );
	m_AngleTrueNeutrinoJet.push_back( acos( trueNeutrinoDirection.Dot( jetAxis ) ) );

	if ( m_fillRootTree ) plotHistograms( trueNeutrinoFourMomentum , recoNeutrinoFourMomentumClose , NeutrinoCovMat );
	m_trueNuPx.push_back( trueNeutrinoFourMomentum.Px() );
	m_trueNuPy.push_back( trueNeutrinoFourMomentum.Py() );
	m_trueNuPz.push_back( trueNeutrinoFourMomentum.Pz() );
	m_trueNuE.push_back( trueNeutrinoFourMomentum.E() );
	m_recoNuClosePx.push_back( recoNeutrinoFourMomentumClose.Px() );
	m_recoNuClosePy.push_back( recoNeutrinoFourMomentumClose.Py() );
	m_recoNuClosePz.push_back( recoNeutrinoFourMomentumClose.Pz() );
	m_recoNuCloseE.push_back( recoNeutrinoFourMomentumClose.E() );
	m_recoNuPosPx.push_back( recoNeutrinoFourMomentumPos.Px() );
	m_recoNuPosPy.push_back( recoNeutrinoFourMomentumPos.Py() );
	m_recoNuPosPz.push_back( recoNeutrinoFourMomentumPos.Pz() );
	m_recoNuPosE.push_back( recoNeutrinoFourMomentumPos.E() );
	m_recoNuNegPx.push_back( recoNeutrinoFourMomentumNeg.Px() );
	m_recoNuNegPy.push_back( recoNeutrinoFourMomentumNeg.Py() );
	m_recoNuNegPz.push_back( recoNeutrinoFourMomentumNeg.Pz() );
	m_recoNuNegE.push_back( recoNeutrinoFourMomentumNeg.E() );

	//m_displayEvent = true;
	if ( m_displayEvent )
	{
		DDMarlinCED::draw( this , 1); // draw everything
	}
	recoHadronFourMomentum = visibleFourMomentum + recoNeutrinoFourMomentumClose;
	fillTrueRecoFourMomentum(	trueVisibleFourMomentumAtSLDVertex , truePVATrueFourMomentum , truePVARecoFourMomentum , recoPVARecoFourMomentum , visibleFourMomentum ,
								trueLeptonFourMomentum , recoLeptonFourMomentum , leptonFourMomentum ,
								truePVATrueChargedFourMomentum , truePVARecoChargedFourMomentum , recoSLDVertexChargedFourMomentum , recoPVARecoChargedFourMomentum , chargedFourMomentum ,
								truePVATrueNeutralFourMomentum , truePVARecoNeutralFourMomentum , recoPVARecoNeutralFourMomentum , neutralFourMomentum ,
								trueNeutrinoFourMomentum , recoNeutrinoFourMomentumClose , trueHadronFourMomentum , recoHadronFourMomentum ,
								trueFlightDirection , recoFlightDirection , flightDirection
							);

	streamlog_out(DEBUG8) << "		Semi-leptonic decay type: " << SLDStatus << std::endl;

	double MomentumPos[3]{ recoNeutrinoFourMomentumPos.Px() , recoNeutrinoFourMomentumPos.Py() , recoNeutrinoFourMomentumPos.Pz() };
	recoNeutrinoPos->setType( -1 * SLDLepton->getPDG() + SLDLepton->getCharge() );
	recoNeutrinoPos->setMomentum( MomentumPos );
	recoNeutrinoPos->setEnergy( recoNeutrinoFourMomentumPos.E() );
	recoNeutrinoPos->setCovMatrix( NeutrinoCovMatPos );
	recoNeutrinoPos->setMass( 0.0 );
	recoNeutrinoPos->setCharge( 0.0 );
	recoNeutrinoPos->setReferencePoint( linkedRecoLepton->getReferencePoint() );
	for ( unsigned int j = 0 ; j < linkedRecoLepton->getParticleIDs().size() ; ++j )
	{
		ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>( linkedRecoLepton->getParticleIDs()[ j ] );
		ParticleIDImpl* outPID = new ParticleIDImpl;
		outPID->setType( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setPDG( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setLikelihood( inPID->getLikelihood() );
		outPID->setAlgorithmType( inPID->getAlgorithmType() ) ;
		for ( unsigned int k = 0 ; k < inPID->getParameters().size() ; ++k ) outPID->addParameter( inPID->getParameters()[ k ] );
		recoNeutrinoPos->addParticleID( outPID );
	}
	recoNeutrinoPos->setParticleIDUsed( linkedRecoLepton->getParticleIDUsed() );
	recoNeutrinoPos->setGoodnessOfPID( linkedRecoLepton->getGoodnessOfPID() );
	recoNeutrinoPos->setStartVertex( semiLeptonicVertex );

	double MomentumNeg[3]{ recoNeutrinoFourMomentumNeg.Px() , recoNeutrinoFourMomentumNeg.Py() , recoNeutrinoFourMomentumNeg.Pz() };
	recoNeutrinoNeg->setType( -1 * SLDLepton->getPDG() + SLDLepton->getCharge() );
	recoNeutrinoNeg->setMomentum( MomentumNeg );
	recoNeutrinoNeg->setEnergy( recoNeutrinoFourMomentumNeg.E() );
	recoNeutrinoNeg->setCovMatrix( NeutrinoCovMatNeg );
	recoNeutrinoNeg->setMass( 0.0 );
	recoNeutrinoNeg->setCharge( 0.0 );
	recoNeutrinoNeg->setReferencePoint( linkedRecoLepton->getReferencePoint() );
	for ( unsigned int j = 0 ; j < linkedRecoLepton->getParticleIDs().size() ; ++j )
	{
		ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>( linkedRecoLepton->getParticleIDs()[ j ] );
		ParticleIDImpl* outPID = new ParticleIDImpl;
		outPID->setType( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setPDG( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setLikelihood( inPID->getLikelihood() );
		outPID->setAlgorithmType( inPID->getAlgorithmType() ) ;
		for ( unsigned int k = 0 ; k < inPID->getParameters().size() ; ++k ) outPID->addParameter( inPID->getParameters()[ k ] );
		recoNeutrinoNeg->addParticleID( outPID );
	}
	recoNeutrinoNeg->setParticleIDUsed( linkedRecoLepton->getParticleIDUsed() );
	recoNeutrinoNeg->setGoodnessOfPID( linkedRecoLepton->getGoodnessOfPID() );
	recoNeutrinoNeg->setStartVertex( semiLeptonicVertex );

	double MomentumZero[3]{ 0.0 , 0.0 , 0.0 };
	recoNeutrinoZero->setType( -1 * SLDLepton->getPDG() + SLDLepton->getCharge() );
	recoNeutrinoZero->setMomentum( MomentumZero );
	recoNeutrinoZero->setEnergy( 0.0 );
	recoNeutrinoZero->setCovMatrix( NeutrinoCovMatZero );
	recoNeutrinoZero->setMass( 0.0 );
	recoNeutrinoZero->setCharge( 0.0 );
	recoNeutrinoZero->setReferencePoint( linkedRecoLepton->getReferencePoint() );
	for ( unsigned int j = 0 ; j < linkedRecoLepton->getParticleIDs().size() ; ++j )
	{
		ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>( linkedRecoLepton->getParticleIDs()[ j ] );
		ParticleIDImpl* outPID = new ParticleIDImpl;
		outPID->setType( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setPDG( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setLikelihood( inPID->getLikelihood() );
		outPID->setAlgorithmType( inPID->getAlgorithmType() ) ;
		for ( unsigned int k = 0 ; k < inPID->getParameters().size() ; ++k ) outPID->addParameter( inPID->getParameters()[ k ] );
		recoNeutrinoZero->addParticleID( outPID );
	}
	recoNeutrinoZero->setParticleIDUsed( linkedRecoLepton->getParticleIDUsed() );
	recoNeutrinoZero->setGoodnessOfPID( linkedRecoLepton->getGoodnessOfPID() );
	recoNeutrinoZero->addParticle( recoNeutrinoPos );
	recoNeutrinoZero->addParticle( recoNeutrinoNeg );
	recoNeutrinoZero->setStartVertex( semiLeptonicVertex );

	streamlog_out(DEBUG7) << "	Creating semi-leptonic vertex LCObject " << std::endl;
	float vertexPosition[ 3 ]{};
	if ( SLDStatus == 4 || SLDStatus == 5 )
	{
		vertexPosition[ 0 ] = sldVertexPosition[ 0 ];
		vertexPosition[ 1 ] = sldVertexPosition[ 1 ];
		vertexPosition[ 2 ] = sldVertexPosition[ 2 ];
	}
	else
	{
		vertexPosition[ 0 ] = primaryVertex->getPosition()[ 0 ];
		vertexPosition[ 1 ] = primaryVertex->getPosition()[ 1 ];
		vertexPosition[ 2 ] = primaryVertex->getPosition()[ 2 ];
	}
	streamlog_out(DEBUG7) << "		Position: ( " << vertexPosition[ 0 ] << " , " << vertexPosition[ 1 ] << " , " << vertexPosition[ 2 ] << " )" << std::endl;
	semiLeptonicVertex->setPrimary( false );
	streamlog_out(DEBUG7) << "		IsPrimary? " << semiLeptonicVertex->isPrimary() << std::endl;
	if ( SLDStatus == 4 )
	{
		semiLeptonicVertex->setAlgorithmType( "LepIn2ndVertex" );
	}
	else if ( SLDStatus == 5 )
	{
		semiLeptonicVertex->setAlgorithmType( "Lep+3rdVertex" );
	}
	streamlog_out(DEBUG7) << "		Algorithm Type : " << semiLeptonicVertex->getAlgorithmType() << std::endl;
	semiLeptonicVertex->setChi2( 0.0 );
	streamlog_out(DEBUG7) << "		Chi2 : " << semiLeptonicVertex->getChi2() << std::endl;
	semiLeptonicVertex->setProbability( 0.0 );
	streamlog_out(DEBUG7) << "		Probability : " << semiLeptonicVertex->getProbability() << std::endl;
	semiLeptonicVertex->setPosition( vertexPosition );
	streamlog_out(DEBUG7) << "		Position : " << semiLeptonicVertex->getPosition()[ 0 ] << " , " << semiLeptonicVertex->getPosition()[ 1 ] << " , " << semiLeptonicVertex->getPosition()[ 2 ] << std::endl;
	semiLeptonicVertex->setAssociatedParticle( semiLeptonicVertexRecoParticle );
	streamlog_out(DEBUG7) << "		AssociatedParticle : " << semiLeptonicVertex->getAssociatedParticle() << std::endl;

	double VisMomentumSLD[3]{ visibleFourMomentum.Px() , visibleFourMomentum.Py() , visibleFourMomentum.Pz() };
	semiLeptonicVertexRecoParticle->setMomentum( VisMomentumSLD );
	semiLeptonicVertexRecoParticle->setEnergy( visibleFourMomentum.E() );
	semiLeptonicVertexRecoParticle->setCovMatrix( NeutrinoCovMatZero );
	semiLeptonicVertexRecoParticle->setMass( visibleFourMomentum.M() );
	semiLeptonicVertexRecoParticle->setCharge( 0.0 );
	semiLeptonicVertexRecoParticle->setReferencePoint( linkedRecoLepton->getReferencePoint() );
	/*
	for ( unsigned int j = 0 ; j < linkedRecoLepton->getParticleIDs().size() ; ++j )
	{
		ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>( linkedRecoLepton->getParticleIDs()[ j ] );
		ParticleIDImpl* outPID = new ParticleIDImpl;
		outPID->setType( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setPDG( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setLikelihood( inPID->getLikelihood() );
		outPID->setAlgorithmType( inPID->getAlgorithmType() ) ;
		for ( unsigned int k = 0 ; k < inPID->getParameters().size() ; ++k ) outPID->addParameter( inPID->getParameters()[ k ] );
		semiLeptonicVertexRecoParticle->addParticleID( outPID );
	}
	*/
	semiLeptonicVertexRecoParticle->addParticle( recoNeutrinoZero );
	for ( unsigned int i_par = 0 ; i_par < recoPVADecayProducts.size() ; ++i_par )
	{
		semiLeptonicVertexRecoParticle->addParticle( recoPVADecayProducts[ i_par ] );
	}

	semiLeptonicVertexRecoParticle->setStartVertex( startVertex );

	trueNeutrinos.push_back( trueNeutrino );
	semiLeptonicVertices.push_back( semiLeptonicVertex );
	semiLeptonicVertexRecoParticles.push_back( semiLeptonicVertexRecoParticle );
	neutrinosOfThisSLD.push_back( recoNeutrinoZero );
	neutrinosOfThisSLD.push_back( recoNeutrinoPos );
	neutrinosOfThisSLD.push_back( recoNeutrinoNeg );
	neutrinos.push_back( neutrinosOfThisSLD );
	sldStatus.push_back( SLDStatus );
	pvaStatus.push_back( PVAStatus );
	solutionSigns.push_back( _trueSolutionSign );
}

void SLDCorrection::checkSLDInput( MCParticle *SLDHadron )
{
	if ( SLDHadron->getDaughters().size() == 0 )
	{
		streamlog_out(MESSAGE) << "" << std::endl;
		streamlog_out(MESSAGE) << "		++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		streamlog_out(MESSAGE) << "		semi-leptonic decay has a decay product without daughter:" << std::endl;
		streamlog_out(MESSAGE) << SLDHadron << std::endl;
		streamlog_out(MESSAGE) << *SLDHadron << std::endl;
		streamlog_out(MESSAGE) << "		++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	}
	else
	{
		streamlog_out(MESSAGE) << "" << std::endl;
		streamlog_out(MESSAGE) << "		--------------------------------------------------------------------------------------------" << std::endl;
		streamlog_out(MESSAGE) << "		semi-leptonic decay has a decay product with " << SLDHadron->getDaughters().size() << " daughters:" << std::endl;
		streamlog_out(MESSAGE) << SLDHadron << "	SHOULD NOT BE COUNTED IN VISIBLE DECAY PRODUCTS!" << std::endl;
		streamlog_out(MESSAGE) << *SLDHadron << std::endl;
		streamlog_out(MESSAGE) << "		--------------------------------------------------------------------------------------------" << std::endl;
		streamlog_out(MESSAGE) << "Daughters: ";
		for ( unsigned int i_daughter = 0 ; i_daughter < SLDHadron->getDaughters().size() ; ++i_daughter )
		{
			streamlog_out(MESSAGE) << SLDHadron->getDaughters()[ i_daughter ] << "		";
		}
		streamlog_out(MESSAGE) << std::endl;
		for ( unsigned int i_daughter = 0 ; i_daughter < SLDHadron->getDaughters().size() ; ++i_daughter )
		{
			checkSLDInput( SLDHadron->getDaughters()[ i_daughter ] );
		}
	}
}

void SLDCorrection::addNeutrinoCovarianceMatrix( TLorentzVector neutrinoFourMomentum , std::vector< float > &NuCovMat )
{
	//	Obtain covariance matrix on Neutrino Four-Momentum (Px,Py,Pz,E)
	//	(Energy error: sigma_E, angular error: sigma_alpha).
	//
	//	Px = E . sin( Theta ) . cos( Phi )
	//	Py = E . sin( Theta ) . sin( Phi )
	//	Pz = E . cos( Theta )
	//	E  = E
	//
	//	dTheta / dAlpha = 1
	//	dPhi / dAlpha = 1 / sin( Theta )
	//
	//	define the jacobian as the 2x4 matrix:
	//
	//
	//
	//			DPx/DAlpha		DPy/DAlpha		DPz/DAlpha		DE/DAlpha
	//	 J =
	//			DPx/DE			DPy/DE			DPz/DE			DE/DE
	//
	//
	//
	//			E . ( cos( Theta) . cos( Phi ) - sin( Phi ) )		E . ( cos( Theta) . sin( Phi ) + cos( Phi ) )		-E . cos( Theta)		0
	//	 J =
	//			sin( Theta ) . cos( Phi )							sin( Theta ) . sin( Phi )							cos( Theta )			1
	//
	//
	//	Order in the PVA errors in the form of covariance matrix:
	//
	//			sigma_alpha ^ 2			0
	//	Cov =
	//					0			sigma_E ^ 2
	//
	//
//	CovMatrixFlightDirection.clear();
	const int rows			= 2; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_dim	= 4;

	TVector3 neutrinoMomentum = neutrinoFourMomentum.Vect();
	double energy = neutrinoFourMomentum.E();
	double theta = neutrinoMomentum.Theta();
	double phi = neutrinoMomentum.Phi();

	float sigma_alpha = m_sigmaAlphaNu;
	float sigma_E = m_sigmaENu;

	TMatrixD covMatrix( kspace_dim , kspace_dim );
	double jacobian_by_rows[rows*columns] =
	{
		energy * ( cos( theta ) * cos( phi ) - sin( phi ) )	,	energy * ( cos( theta ) * sin( phi ) - cos( phi ) )	,	-1.0 * energy * sin( theta )	, 	0.0	,
		sin( theta ) * cos( phi )							,	sin( theta ) * sin( phi )							, cos( theta )						,	1.0
	};
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG6) << "	Jacobian :	" << jacobian( 0 , 0 ) << "	,	" << jacobian( 0 , 1 ) << "	,	" << jacobian( 0 , 2 ) << "	,	" << jacobian( 0 , 3 ) << std::endl;
	streamlog_out(DEBUG6) << "			" << jacobian( 1 , 0 ) << "	,	" << jacobian( 1 , 1 ) << "	,	" << jacobian( 1 , 2 ) << "	,	" << jacobian( 1 , 3 ) << std::endl;
	double cov_matrix_by_rows[rows*rows] =
	{
		pow( sigma_alpha , 2 )	,			0.0		,
					0.0			,	pow( sigma_E , 2 )
	};
	streamlog_out(DEBUG6) << "	cov_matrix_by_rows[ 0 ] :	" << cov_matrix_by_rows[ 0 ] << std::endl;
	streamlog_out(DEBUG6) << "	cov_matrix_by_rows[ 1 ] :	" << cov_matrix_by_rows[ 1 ] << std::endl;
	streamlog_out(DEBUG6) << "	cov_matrix_by_rows[ 2 ] :	" << cov_matrix_by_rows[ 2 ] << std::endl;
	streamlog_out(DEBUG6) << "	cov_matrix_by_rows[ 3 ] :	" << cov_matrix_by_rows[ 3 ] << std::endl;
	TMatrixD CovMat(rows,rows, cov_matrix_by_rows, "C");
	streamlog_out(DEBUG6) << "	CovMat :	" << CovMat( 0 , 0 ) << std::endl;
	covMatrix.Mult( TMatrixD( jacobian , TMatrixD::kTransposeMult , CovMat ) , jacobian );
	NuCovMat[ 0 ] = covMatrix( 0 , 0 );
	NuCovMat[ 1 ] = covMatrix( 1 , 0 );
	NuCovMat[ 2 ] = covMatrix( 1 , 1 );
	NuCovMat[ 3 ] = covMatrix( 2 , 0 );
	NuCovMat[ 4 ] = covMatrix( 2 , 1 );
	NuCovMat[ 5 ] = covMatrix( 2 , 2 );
	NuCovMat[ 6 ] = covMatrix( 3 , 0 );
	NuCovMat[ 7 ] = covMatrix( 3 , 1 );
	NuCovMat[ 8 ] = covMatrix( 3 , 2 );
	NuCovMat[ 9 ] = covMatrix( 3 , 3 );
	streamlog_out(DEBUG6) << "	Neutrino Theta		= " << theta << std::endl;
	streamlog_out(DEBUG6) << "	Neutrino Phi		= " << phi << std::endl;
	streamlog_out(DEBUG6) << "	Neutrino Energy		= " << energy << std::endl;
	streamlog_out(DEBUG6) << "	sigmaAlpha			= " << sigma_alpha << std::endl;
	streamlog_out(DEBUG6) << "	sigmaEnergy			= " << sigma_E << std::endl;
	streamlog_out(DEBUG6) << "	Partial CovMat :	" << std::endl;
	streamlog_out(DEBUG6) << "				" << covMatrix( 0 , 0 ) << std::endl;
	streamlog_out(DEBUG6) << "				" << covMatrix( 1 , 0 ) << "	,	" << covMatrix( 1 , 1 ) << std::endl;
	streamlog_out(DEBUG6) << "				" << covMatrix( 2 , 0 ) << "	,	" << covMatrix( 2 , 1 ) << "	,	" << covMatrix( 2 , 2 ) << std::endl;
	streamlog_out(DEBUG6) << "				" << covMatrix( 3 , 0 ) << "	,	" << covMatrix( 3 , 1 ) << "	,	" << covMatrix( 3 , 2 ) << "	,	" << covMatrix( 3 , 3 ) << std::endl;
	streamlog_out(DEBUG6) << "	Overall CovMat :	" << std::endl;
	streamlog_out(DEBUG6) << "				" << NuCovMat[ 0 ] << std::endl;
	streamlog_out(DEBUG6) << "				" << NuCovMat[ 1 ] << "	,	" << NuCovMat[ 2 ] << std::endl;
	streamlog_out(DEBUG6) << "				" << NuCovMat[ 3 ] << "	,	" << NuCovMat[ 4 ] << "	,	" << NuCovMat[ 5 ] << std::endl;
	streamlog_out(DEBUG6) << "				" << NuCovMat[ 6 ] << "	,	" << NuCovMat[ 7 ] << "	,	" << NuCovMat[ 8 ] << "	,	" << NuCovMat[ 9 ] << std::endl;
	streamlog_out(DEBUG6) << "" << std::endl;
}

void SLDCorrection::evaluateInputCovMat(	TLorentzVector trueVisibleFourMomentum , TVector3 trueFlightDirection , TLorentzVector trueNeutrinoFourMomentum ,
											TLorentzVector visibleFourMomentum , TVector3 flightDirection , TLorentzVector recoNeutrinoFourMomentum ,
											std::vector< float > CovMatDetector , std::vector< float > CovMatFlightDirection , std::vector< float > CovMatNeutrino )
{
	m_trueVisibleFourMomentumPx.push_back( trueVisibleFourMomentum.Px() );
	m_trueVisibleFourMomentumPy.push_back( trueVisibleFourMomentum.Py() );
	m_trueVisibleFourMomentumPz.push_back( trueVisibleFourMomentum.Pz() );
	m_trueVisibleFourMomentumE.push_back( trueVisibleFourMomentum.E() );
	m_recoVisibleFourMomentumPx.push_back( visibleFourMomentum.Px() );
	m_recoVisibleFourMomentumPy.push_back( visibleFourMomentum.Py() );
	m_recoVisibleFourMomentumPz.push_back( visibleFourMomentum.Pz() );
	m_recoVisibleFourMomentumE.push_back( visibleFourMomentum.E() );
	m_residualVisibleFourMomentumPx.push_back( visibleFourMomentum.Px() - trueVisibleFourMomentum.Px() );
	m_residualVisibleFourMomentumPy.push_back( visibleFourMomentum.Py() - trueVisibleFourMomentum.Py() );
	m_residualVisibleFourMomentumPz.push_back( visibleFourMomentum.Pz() - trueVisibleFourMomentum.Pz() );
	m_residualVisibleFourMomentumE.push_back( visibleFourMomentum.E() - trueVisibleFourMomentum.E() );
	m_sigmaPxPx_Det.push_back( CovMatDetector[ 0 ] );
	m_sigmaPxPy_Det.push_back( CovMatDetector[ 1 ] );
	m_sigmaPyPy_Det.push_back( CovMatDetector[ 2 ] );
	m_sigmaPxPz_Det.push_back( CovMatDetector[ 3 ] );
	m_sigmaPyPz_Det.push_back( CovMatDetector[ 4 ] );
	m_sigmaPzPz_Det.push_back( CovMatDetector[ 5 ] );
	m_sigmaPxE_Det.push_back( CovMatDetector[ 6 ] );
	m_sigmaPyE_Det.push_back( CovMatDetector[ 7 ] );
	m_sigmaPzE_Det.push_back( CovMatDetector[ 8 ] );
	m_sigmaEE_Det.push_back( CovMatDetector[ 9 ] );
	m_normalizedResidualVisibleFourMomentumPx.push_back( ( visibleFourMomentum.Px() - trueVisibleFourMomentum.Px() ) / sqrt( CovMatDetector[ 0 ] ) );
	m_normalizedResidualVisibleFourMomentumPy.push_back( ( visibleFourMomentum.Py() - trueVisibleFourMomentum.Py() ) / sqrt( CovMatDetector[ 2 ] ) );
	m_normalizedResidualVisibleFourMomentumPz.push_back( ( visibleFourMomentum.Pz() - trueVisibleFourMomentum.Pz() ) / sqrt( CovMatDetector[ 5 ] ) );
	m_normalizedResidualVisibleFourMomentumE.push_back( ( visibleFourMomentum.E() - trueVisibleFourMomentum.E() ) / sqrt( CovMatDetector[ 9 ] ) );

	m_trueFlightDirectionUx.push_back( trueFlightDirection.X() );
	m_trueFlightDirectionUy.push_back( trueFlightDirection.Y() );
	m_trueFlightDirectionUz.push_back( trueFlightDirection.Z() );
	m_recoFlightDirectionUx.push_back( flightDirection.X() );
	m_recoFlightDirectionUy.push_back( flightDirection.Y() );
	m_recoFlightDirectionUz.push_back( flightDirection.Z() );
	m_residualFlightDirectionUx.push_back( flightDirection.X() - trueFlightDirection.X() );
	m_residualFlightDirectionUy.push_back( flightDirection.Y() - trueFlightDirection.Y() );
	m_residualFlightDirectionUz.push_back( flightDirection.Z() - trueFlightDirection.Z() );
	m_sigmaUxUx.push_back( CovMatFlightDirection[ 0 ] );
	m_sigmaUxUy.push_back( CovMatFlightDirection[ 1 ] );
	m_sigmaUyUy.push_back( CovMatFlightDirection[ 2 ] );
	m_sigmaUxUz.push_back( CovMatFlightDirection[ 3 ] );
	m_sigmaUyUz.push_back( CovMatFlightDirection[ 4 ] );
	m_sigmaUzUz.push_back( CovMatFlightDirection[ 5 ] );
	m_normalizedResidualFlightDirectionUx.push_back( ( flightDirection.X() - trueFlightDirection.X() ) / ( sqrt( CovMatFlightDirection[ 0 ] ) ) );
	m_normalizedResidualFlightDirectionUy.push_back( ( flightDirection.Y() - trueFlightDirection.Y() ) / ( sqrt( CovMatFlightDirection[ 2 ] ) ) );
	m_normalizedResidualFlightDirectionUz.push_back( ( flightDirection.Z() - trueFlightDirection.Z() ) / ( sqrt( CovMatFlightDirection[ 5 ] ) ) );

	m_trueNeutrinoFourMomentumPx.push_back( trueNeutrinoFourMomentum.Px() );
	m_trueNeutrinoFourMomentumPy.push_back( trueNeutrinoFourMomentum.Py() );
	m_trueNeutrinoFourMomentumPz.push_back( trueNeutrinoFourMomentum.Pz() );
	m_trueNeutrinoFourMomentumE.push_back( trueNeutrinoFourMomentum.E() );
	m_recoNeutrinoFourMomentumClosePx.push_back( recoNeutrinoFourMomentum.Px() );
	m_recoNeutrinoFourMomentumClosePy.push_back( recoNeutrinoFourMomentum.Py() );
	m_recoNeutrinoFourMomentumClosePz.push_back( recoNeutrinoFourMomentum.Pz() );
	m_recoNeutrinoFourMomentumCloseE.push_back( recoNeutrinoFourMomentum.E() );
	m_residualNeutrinoFourMomentumPx.push_back( recoNeutrinoFourMomentum.Px() - trueNeutrinoFourMomentum.Px() );
	m_residualNeutrinoFourMomentumPy.push_back( recoNeutrinoFourMomentum.Py() - trueNeutrinoFourMomentum.Py() );
	m_residualNeutrinoFourMomentumPz.push_back( recoNeutrinoFourMomentum.Pz() - trueNeutrinoFourMomentum.Pz() );
	m_residualNeutrinoFourMomentumE.push_back( recoNeutrinoFourMomentum.E() - trueNeutrinoFourMomentum.E() );
	m_sigmaNeutrinoPxPx.push_back( CovMatNeutrino[ 0 ] );
	m_sigmaNeutrinoPxPy.push_back( CovMatNeutrino[ 1 ] );
	m_sigmaNeutrinoPyPy.push_back( CovMatNeutrino[ 2 ] );
	m_sigmaNeutrinoPxPz.push_back( CovMatNeutrino[ 3 ] );
	m_sigmaNeutrinoPyPz.push_back( CovMatNeutrino[ 4 ] );
	m_sigmaNeutrinoPzPz.push_back( CovMatNeutrino[ 5 ] );
	m_sigmaNeutrinoPxE.push_back( CovMatNeutrino[ 6 ] );
	m_sigmaNeutrinoPyE.push_back( CovMatNeutrino[ 7 ] );
	m_sigmaNeutrinoPzE.push_back( CovMatNeutrino[ 8 ] );
	m_sigmaNeutrinoEE.push_back( CovMatNeutrino[ 9 ] );
	m_normalizedResidualNeutrinoFourMomentumPx.push_back( ( recoNeutrinoFourMomentum.Px() - trueNeutrinoFourMomentum.Px() ) / sqrt( CovMatNeutrino[ 0 ] ) );
	m_normalizedResidualNeutrinoFourMomentumPy.push_back( ( recoNeutrinoFourMomentum.Py() - trueNeutrinoFourMomentum.Py() ) / sqrt( CovMatNeutrino[ 2 ] ) );
	m_normalizedResidualNeutrinoFourMomentumPz.push_back( ( recoNeutrinoFourMomentum.Pz() - trueNeutrinoFourMomentum.Pz() ) / sqrt( CovMatNeutrino[ 5 ] ) );
	m_normalizedResidualNeutrinoFourMomentumE.push_back( ( recoNeutrinoFourMomentum.E() - trueNeutrinoFourMomentum.E() ) / sqrt( CovMatNeutrino[ 9 ] ) );
	TVector3 trueNeutrinoDirection( trueNeutrinoFourMomentum.Px() , trueNeutrinoFourMomentum.Py() , trueNeutrinoFourMomentum.Pz() ); trueNeutrinoDirection.SetMag( 1.0 );
	TVector3 recoNeutrinoDirection( recoNeutrinoFourMomentum.Px() , recoNeutrinoFourMomentum.Py() , recoNeutrinoFourMomentum.Pz() ); recoNeutrinoDirection.SetMag( 1.0 );
	m_recoNeutrinoDirectionError.push_back( acos( recoNeutrinoDirection.Dot( trueNeutrinoDirection ) ) );
}

void SLDCorrection::showTrueParameters( MCParticle *SLDLepton )
{
	TLorentzVector true4mom( 0.0 , 0.0 , 0.0 , 0.0 );
	MCParticle* parentHadron = SLDLepton->getParents()[ 0 ];
	streamlog_out(DEBUG4) << "	PARENT HADRON:" << std::endl;
	streamlog_out(DEBUG4) << *parentHadron << std::endl;
	TVector3 trueFliDir = TVector3( parentHadron->getMomentumAtEndpoint() );
	trueFliDir.SetMag( 1.0 );
	for ( unsigned int i_mcp = 0 ; i_mcp < parentHadron->getDaughters().size() ; ++i_mcp )
	{
		MCParticle* daughter = parentHadron->getDaughters()[ i_mcp ];
		streamlog_out(DEBUG4) << *daughter << std::endl;
		if ( std::fabs( daughter->getPDG() ) != 12 && std::fabs( daughter->getPDG() ) != 14 && std::fabs( daughter->getPDG() ) != 16 )
		{
			true4mom += TLorentzVector( daughter->getMomentum() , daughter->getEnergy() );
			streamlog_out(DEBUG4) << " ONE MCPARTICLE IS ADDED" << std::endl;
		}
	}
	streamlog_out(DEBUG4) << "	TRUE PARENT HADRON MASS =  			" << parentHadron->getMass() << std::endl;
	streamlog_out(DEBUG4) << "	TRUE FLIHT DIRECTION (x,y,z): 			" << trueFliDir.X() << "	,	" << trueFliDir.Y() << "	,	" << trueFliDir.Z() << std::endl;
	streamlog_out(DEBUG4) << "	TRUE VISIBLE FOUR-MOMENTUM (Px,Py,Pz,E):	" << true4mom.Px() << "	,	" << true4mom.Py() << "	,	" << true4mom.Pz() << "	,	" << true4mom.E() << std::endl;
	TVector3 truePvisPar = trueFliDir.Dot( true4mom.Vect() ) * trueFliDir;
	TVector3 truePvisNor = true4mom.Vect() - truePvisPar;
	streamlog_out(DEBUG4) << "	TRUE VISIBLE FOUR-MOMENTUM(par) (Px,Py,Pz):	" << truePvisPar.Px() << "	,	" << truePvisPar.Py() << "	,	" << truePvisPar.Pz() << std::endl;
	streamlog_out(DEBUG4) << "	TRUE VISIBLE FOUR-MOMENTUM(nor) (Px,Py,Pz):	" << truePvisNor.Px() << "	,	" << truePvisNor.Py() << "	,	" << truePvisNor.Pz() << std::endl;
}

TLorentzVector SLDCorrection::getNeutrinoFourMomentum( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double ParentHadronMass , float solutionSign )
{
	int sign = ( solutionSign != 0 ? solutionSign / abs( solutionSign ) : 1 );
	m_solutionSign.push_back( sign );
	const char *solSign = ( sign >= 0 ? "+" : "-" );
	streamlog_out(DEBUG4) << "" << std::endl;
	streamlog_out(DEBUG4) << "		--------------------------------------------" << std::endl;
	streamlog_out(DEBUG4) << "		Calculate Neutrino 4-Momentum for " << solSign << " solution" << std::endl;
	streamlog_out(DEBUG4) << "		--------------------------------------------" << std::endl;

	streamlog_out(DEBUG4) << "		Test 1, |flightDirection| = " << flightDirection.Mag() << std::endl;
	flightDirection.SetMag( 1.0 );
	streamlog_out(DEBUG4) << "		flightDirection:			( " << flightDirection.Px() << "	, " << flightDirection.Py() << "	, " << flightDirection.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG4) << "		Parent Hadron Mass =	 " << ParentHadronMass << std::endl;

	streamlog_out(DEBUG4) << "		Visible 4-Momentum:			( " << visibleFourMomentum.Px() << "	, " << visibleFourMomentum.Py() << "	, " << visibleFourMomentum.Pz() << "	, " << visibleFourMomentum.E() << " )" << std::endl;

	double visible_mass		= visibleFourMomentum.M();
	streamlog_out(DEBUG4) << "		Visible Inv Mass:	" << visible_mass << std::endl;

	double visible_E		= visibleFourMomentum.E();
	m_E_vis.push_back( visible_E );
	streamlog_out(DEBUG4) << "		Visible Energy:									" << visible_E << std::endl;

	TVector3 visible_p		= TVector3( visibleFourMomentum.Px() , visibleFourMomentum.Py() , visibleFourMomentum.Pz() );
	streamlog_out(DEBUG4) << "		Visible Momentum:			( " << visible_p.Px() << "	, " << visible_p.Py() << "	, " << visible_p.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG4) << "		(Visible Momentum).(FlightDirection):							" << visible_p.Dot( flightDirection ) << std::endl;

	TVector3 visible_p_par		= visible_p.Dot( flightDirection ) * flightDirection;
	m_P_vis_par.push_back( visible_p_par.Mag() );
	streamlog_out(DEBUG4) << "		Visible Momentum (par):			( " << visible_p_par.Px() << "	, " << visible_p_par.Py() << "	, " << visible_p_par.Pz() << "	)" << std::endl;

	TVector3 visible_p_nor		= visible_p - visible_p_par;
	m_P_vis_nor.push_back( visible_p_nor.Mag() );
	m_P_vis_nor_prime.push_back( visible_p_nor.Mag() );
	streamlog_out(DEBUG4) << "		Visible Momentum (nor):			( " << visible_p_nor.Px() << "	, " << visible_p_nor.Py() << "	, " << visible_p_nor.Pz() << "	)" << std::endl;

	double visible_E_prime		= ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass );
	m_E_vis_prime.push_back( visible_E_prime );
	streamlog_out(DEBUG4) << "		Visible Energy (prime):								" << visible_E_prime << std::endl;

	TVector3 visible_p_par_prime	= solutionSign * sqrt( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) - visible_p_nor.Mag2() ) * flightDirection;
	streamlog_out(DEBUG4) << "		Visible Momentum (par-prime):		( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;
	if ( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) < visible_p_nor.Mag2() )
	{
		visible_p_par_prime	= solutionSign * std::numeric_limits<double>::min() * flightDirection;
	}
	m_P_vis_par_prime.push_back( visible_p_par_prime.Mag() );
	streamlog_out(DEBUG4) << "		Visible Momentum (par-prime):		( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;

	double parent_hadron_E		= ( ( visibleFourMomentum.E() * ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) ) - visible_p_par.Dot( visible_p_par_prime ) ) * ParentHadronMass / ( pow( visible_mass , 2 ) + visible_p_nor.Mag2() );
	streamlog_out(DEBUG4) << "		Parent Hadron Energy =									" << parent_hadron_E << std::endl;
	TVector3 parent_hadron_p	= sqrt( pow( parent_hadron_E , 2 ) - pow( ParentHadronMass , 2 ) ) * flightDirection;
	streamlog_out(DEBUG4) << "		Parent Hadron Momentum:			( " << parent_hadron_p.Px() << "	, " << parent_hadron_p.Py() << "	, " << parent_hadron_p.Pz() << "	, " << parent_hadron_E << " )" << std::endl;

	double sigma_E_vis = 0.0;
	double sigma_E_vis_prime = 0.0;
	double sigma_p_vis_par = 0.0;
	double sigma_p_vis_par_prime = 0.0;
	double sigma_p_vis_nor = 0.0;

	double sigma_parent_hadron_E2	=	pow( ParentHadronMass , 2 ) * (
						pow( visible_E_prime / ( pow( visible_mass , 2 ) + pow( visible_p_nor.Mag() , 2 ) ) , 2 ) * pow( sigma_E_vis , 2 ) +
						pow( visible_E / ( pow( visible_mass , 2 ) + pow( visible_p_nor.Mag() , 2 ) ) , 2 ) * pow( sigma_E_vis_prime , 2 ) +
						pow( visible_p_par_prime.Mag() / ( pow( visible_mass , 2 ) + pow( visible_p_nor.Mag() , 2 ) ) , 2 ) * pow( sigma_p_vis_par , 2 ) +
						pow( visible_p_par.Mag() / ( pow( visible_mass , 2 ) + pow( visible_p_nor.Mag() , 2 ) ) , 2 ) * pow( sigma_p_vis_par_prime , 2 ) +
						pow( 2 * visible_p_nor.Mag() * ( visible_E_prime * visible_E - visible_p_par_prime.Mag() * visible_p_par.Mag() ) / pow( ( pow( visible_mass , 2 ) + pow( visible_p_nor.Mag() , 2 ) ) , 2 ) , 2 ) * pow( sigma_p_vis_nor , 2 )
					);
					streamlog_out(DEBUG4) << "		Parent Hadron Sigma_E =									" << sigma_parent_hadron_E2 << std::endl;

	double Neutrino_E		= parent_hadron_E - visible_E;
	streamlog_out(DEBUG4) << "		Neutrino Energy:									" << Neutrino_E << std::endl;

	TVector3 Neutrino_p_nor		= -1 * visible_p_nor;
	streamlog_out(DEBUG4) << "		Neutrino Momentum (nor):		( " << Neutrino_p_nor.Px() << "	, " << Neutrino_p_nor.Py() << "	, " << Neutrino_p_nor.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p_par		= sqrt( pow( Neutrino_E , 2 ) - Neutrino_p_nor.Mag2() ) * flightDirection;
	streamlog_out(DEBUG4) << "		Neutrino Momentum (par):		( " << Neutrino_p_par.Px() << "	, " << Neutrino_p_par.Py() << "	, " << Neutrino_p_par.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p		= Neutrino_p_nor + Neutrino_p_par;
	streamlog_out(DEBUG4) << "		Neutrino Momentum:			( " << Neutrino_p_par.Px() << "	, " << Neutrino_p_par.Py() << "	, " << Neutrino_p_par.Pz() << "	)" << std::endl;

	TLorentzVector Neutrino_tlv( Neutrino_p , Neutrino_E );
	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG2) << "	Reconstructed Neutrino Mass:	" << Neutrino_tlv.Mag() << std::endl;
	streamlog_out(DEBUG2) << "	Reconstructed Neutrino 4-Momentum:		( " << Neutrino_tlv.Px() << "	, " << Neutrino_tlv.Py() << "	, " << Neutrino_tlv.Pz() << "	, " << Neutrino_tlv.E() << " )" << std::endl;
	return Neutrino_tlv;
}

TLorentzVector SLDCorrection::getNeutrinoFourMomentumModified( TVector3 &flightDirection , TLorentzVector visibleFourMomentum , double ParentHadronMass , float solutionSign )
{
	int sign = ( solutionSign != 0 ? solutionSign / abs( solutionSign ) : 1 );
	m_solutionSign.push_back( sign );
	const char *solSign = ( sign >= 0 ? "+" : "-" );
	streamlog_out(DEBUG9) << "" << std::endl;
	streamlog_out(DEBUG9) << "		--------------------------------------------" << std::endl;
	streamlog_out(DEBUG9) << "		Calculate Neutrino 4-Momentum for " << solSign << " solution" << std::endl;
	streamlog_out(DEBUG9) << "		--------------------------------------------" << std::endl;

	streamlog_out(DEBUG9) << "		Test 1, |flightDirection| = " << flightDirection.Mag() << std::endl;
	flightDirection.SetMag( 1.0 );
	streamlog_out(DEBUG9) << "		flightDirection:			( " << flightDirection.Px() << "	, " << flightDirection.Py() << "	, " << flightDirection.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG9) << "		Parent Hadron Mass =	 " << ParentHadronMass << std::endl;

	streamlog_out(DEBUG9) << "		Visible 4-Momentum:			( " << visibleFourMomentum.Px() << "	, " << visibleFourMomentum.Py() << "	, " << visibleFourMomentum.Pz() << "	, " << visibleFourMomentum.E() << " )" << std::endl;

	double visible_mass		= visibleFourMomentum.M();
	streamlog_out(DEBUG9) << "		Visible Inv Mass:	" << visible_mass << std::endl;

	double visible_E		= visibleFourMomentum.E();
	m_E_vis.push_back( visible_E );
	streamlog_out(DEBUG9) << "		Visible Energy:									" << visible_E << std::endl;

	TVector3 visible_p		= TVector3( visibleFourMomentum.Px() , visibleFourMomentum.Py() , visibleFourMomentum.Pz() );
	streamlog_out(DEBUG9) << "		Visible Momentum:			( " << visible_p.Px() << "	, " << visible_p.Py() << "	, " << visible_p.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG9) << "		(Visible Momentum).(FlightDirection):							" << visible_p.Dot( flightDirection ) << std::endl;

	double P_vis_par		= visible_p.Dot( flightDirection );
	TVector3 visible_p_par		= P_vis_par * flightDirection;
	m_P_vis_par.push_back( P_vis_par );
	streamlog_out(DEBUG9) << "		|Visible Momentum (par)|:	" << P_vis_par << std::endl;
	streamlog_out(DEBUG9) << "		Visible Momentum (par):			( " << visible_p_par.Px() << "	, " << visible_p_par.Py() << "	, " << visible_p_par.Pz() << "	)" << std::endl;

	TVector3 visible_p_nor		= visible_p - visible_p_par;
	double P_vis_nor		= visible_p_nor.Mag();
	m_P_vis_nor.push_back( P_vis_nor );
	m_P_vis_nor_prime.push_back( P_vis_nor );
	streamlog_out(DEBUG9) << "		|Visible Momentum (nor)|:	" << P_vis_nor << std::endl;
	streamlog_out(DEBUG9) << "		Visible Momentum (nor):			( " << visible_p_nor.Px() << "	, " << visible_p_nor.Py() << "	, " << visible_p_nor.Pz() << "	)" << std::endl;

	double visible_E_prime		= ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass );
	m_E_vis_prime.push_back( visible_E_prime );
	streamlog_out(DEBUG9) << "		Visible Energy (prime):								" << visible_E_prime << std::endl;


	double visible_p_par_prime_squared = pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) - visible_p_nor.Mag2();
	streamlog_out(DEBUG9) << "		Visible Momentum (par-prime)2:		= " << visible_p_par_prime_squared << std::endl;
	m_P_vis_par_prime_squared.push_back( visible_p_par_prime_squared );
	if ( visible_p_par_prime_squared < 0.0 )
	{
		TVector3 visible_p_nor_direction = visible_p_nor; visible_p_nor_direction.SetMag( 1.0 );
		double visible_p_nor_corrected = ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass );
		streamlog_out(DEBUG9) << "		|corrected visible_p_nor|=		" << visible_p_nor_corrected << std::endl;
		double visible_p_par_corrected = sqrt( visible_p.Mag2() - pow( visible_p_nor_corrected , 2 ) );
		streamlog_out(DEBUG9) << "		|corrected visible_p_par|=		" << visible_p_par_corrected << std::endl;
		TVector3 visible_p_direction = visible_p; visible_p_direction.SetMag( 1.0 );
		double initial_alpha = acos( visible_p_direction.Dot( flightDirection ) );
		double new_alpha = acos( visible_p_par_corrected / visible_p.Mag() );
		double correctionAngle = initial_alpha - new_alpha;
		TVector3 newFlightDirection = flightDirection + sin( correctionAngle ) * visible_p_nor_direction;
		newFlightDirection.SetMag( 1.0 );
		visible_p_par		= visible_p.Dot( newFlightDirection ) * newFlightDirection;
		visible_p_nor		= visible_p - visible_p_par;
		flightDirection = newFlightDirection;
		streamlog_out(DEBUG9) << "		new flightDirection:		( " << flightDirection.Px() << "	, " << flightDirection.Py() << "	, " << flightDirection.Pz() << "	)" << std::endl;
		visible_p_par_prime_squared = std::numeric_limits<double>::min();
	}
	streamlog_out(DEBUG9) << "		Visible Momentum (par-prime)2:		= " << visible_p_par_prime_squared << std::endl;
	TVector3 visible_p_par_prime	= solutionSign * sqrt( visible_p_par_prime_squared ) * flightDirection;
	m_P_vis_par_prime.push_back( visible_p_par_prime.Mag() );
	streamlog_out(DEBUG9) << "		Visible Momentum (par-prime):		( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;

	streamlog_out(DEBUG9) << "		|Visible Momentum (par-prime)|:	" << visible_p_par_prime.Mag() << std::endl;

	double parent_hadron_E		= ( ( visibleFourMomentum.E() * visible_E_prime ) - visible_p_par.Dot( visible_p_par_prime ) ) * ParentHadronMass / ( pow( visible_mass , 2 ) + visible_p_nor.Mag2() );
	streamlog_out(DEBUG9) << "		Parent Hadron Energy =									" << parent_hadron_E << std::endl;
	TVector3 parent_hadron_p	= ( ParentHadronMass >= parent_hadron_E ? 0.0 : sqrt( pow( parent_hadron_E , 2 ) - pow( ParentHadronMass , 2 ) ) ) * flightDirection;
	streamlog_out(DEBUG9) << "		Parent Hadron Momentum:			( " << parent_hadron_p.Px() << "	, " << parent_hadron_p.Py() << "	, " << parent_hadron_p.Pz() << "	) = " << parent_hadron_p.Mag() << std::endl;

	TVector3 Neutrino_p_nor		= -1.0 * visible_p_nor;
	streamlog_out(DEBUG9) << "		Neutrino Momentum (nor):		( " << Neutrino_p_nor.Px() << "	, " << Neutrino_p_nor.Py() << "	, " << Neutrino_p_nor.Pz() << "	) = " << Neutrino_p_nor.Mag() << std::endl;

	TVector3 Neutrino_p_par		= parent_hadron_p - visible_p_par;
	streamlog_out(DEBUG9) << "		Neutrino Momentum (par):		( " << Neutrino_p_par.Px() << "	, " << Neutrino_p_par.Py() << "	, " << Neutrino_p_par.Pz() << "	) = " << parent_hadron_p.Mag() - visible_p_par.Mag() << std::endl;

	TVector3 Neutrino_p		= Neutrino_p_nor + Neutrino_p_par;
	streamlog_out(DEBUG9) << "		Neutrino Momentum:			( " << Neutrino_p.Px() << "	, " << Neutrino_p.Py() << "	, " << Neutrino_p.Pz() << "	)" << std::endl;

	double Neutrino_E		= Neutrino_p.Mag();
	streamlog_out(DEBUG9) << "		Neutrino Energy:									" << Neutrino_E << std::endl;

	TLorentzVector Neutrino_tlv( Neutrino_p , Neutrino_E );
	streamlog_out(DEBUG9) << "" << std::endl;
	streamlog_out(DEBUG9) << "	Reconstructed Neutrino Mass:	" << Neutrino_tlv.Mag() << std::endl;
	streamlog_out(DEBUG9) << "	Reconstructed Neutrino 4-Momentum:		( " << Neutrino_tlv.Px() << "	, " << Neutrino_tlv.Py() << "	, " << Neutrino_tlv.Pz() << "	, " << Neutrino_tlv.E() << " )" << std::endl;
	return Neutrino_tlv;
}

TLorentzVector SLDCorrection::getNeutrinoFourMomentumStandardMethod( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , float solutionSign )
{
	int sign = ( solutionSign != 0 ? solutionSign / abs( solutionSign ) : 1 );
	m_solutionSign.push_back( sign );
	const char *solSign = ( sign >= 0 ? "+" : "-" );
	streamlog_out(DEBUG4) << "" << std::endl;
	streamlog_out(DEBUG4) << "		--------------------------------------------" << std::endl;
	streamlog_out(DEBUG4) << "		Calculate Neutrino 4-Momentum for " << solSign << " solution" << std::endl;
	streamlog_out(DEBUG4) << "		--------------------------------------------" << std::endl;

	streamlog_out(DEBUG4) << "		Test 1, |flightDirection| = " << flightDirection.Mag() << std::endl;
	flightDirection.SetMag( 1.0 );
	streamlog_out(DEBUG4) << "		flightDirection:			( " << flightDirection.Px() << "	, " << flightDirection.Py() << "	, " << flightDirection.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG4) << "		Parent Hadron Mass =	 " << parentHadronMass << std::endl;

	streamlog_out(DEBUG4) << "		Visible 4-Momentum:			( " << visibleFourMomentum.Px() << "	, " << visibleFourMomentum.Py() << "	, " << visibleFourMomentum.Pz() << "	, " << visibleFourMomentum.E() << " )" << std::endl;

	double visible_mass		= visibleFourMomentum.M();
	streamlog_out(DEBUG4) << "		Visible Inv Mass:	" << visible_mass << std::endl;

	double visible_E		= visibleFourMomentum.E();
	m_E_vis.push_back( visible_E );
	streamlog_out(DEBUG8) << "		Visible Energy:									" << visible_E << std::endl;

	TVector3 visible_p		= TVector3( visibleFourMomentum.Px() , visibleFourMomentum.Py() , visibleFourMomentum.Pz() );
	streamlog_out(DEBUG4) << "		Visible Momentum:			( " << visible_p.Px() << "	, " << visible_p.Py() << "	, " << visible_p.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG4) << "		(Visible Momentum).(FlightDirection):							" << visible_p.Dot( flightDirection ) << std::endl;

	double P_vis_par		= visible_p.Dot( flightDirection );
	TVector3 visible_p_par		= P_vis_par * flightDirection;
	m_P_vis_par.push_back( P_vis_par );
	streamlog_out(DEBUG8) << "		|Visible Momentum (par)|:	" << P_vis_par << std::endl;
	streamlog_out(DEBUG4) << "		Visible Momentum (par):			( " << visible_p_par.Px() << "	, " << visible_p_par.Py() << "	, " << visible_p_par.Pz() << "	)" << std::endl;

	double P_vis_nor		= sqrt( visible_p.Mag2() - pow( P_vis_par , 2 ) );
	TVector3 visible_p_nor		= visible_p - visible_p_par;
	m_P_vis_nor.push_back( P_vis_nor );
	m_P_vis_nor_prime.push_back( P_vis_nor );
	streamlog_out(DEBUG8) << "		|Visible Momentum (nor)|:	" << P_vis_nor << std::endl;
	streamlog_out(DEBUG4) << "		Visible Momentum (nor):			( " << visible_p_nor.Px() << "	, " << visible_p_nor.Py() << "	, " << visible_p_nor.Pz() << "	)" << std::endl;
/*
	double A = P_vis_par * ( pow( parentHadronMass , 2 ) - pow( visible_mass , 2 ) - 2.0 * pow( P_vis_nor , 2 ) );
	streamlog_out(DEBUG8) << "		(A) ( M_B^2 - M_vis^2 - 2 * P_vis(nor)^2 ) * P_vis(par) =	" << A << std::endl;
	double B = 4.0 * pow( P_vis_nor , 2 ) * pow( visible_E , 2 ) - pow( pow( parentHadronMass , 2 ) - pow( visible_mass , 2 ) - 2.0 * pow( P_vis_nor , 2 ) , 2 );
	streamlog_out(DEBUG8) << "		(B) 4 * E_vis^2 * P_vis(par)^2 - ( M_B^2 - M_vis^2 - 2 * P_vis(nor)^2 )^2 =	" << B << std::endl;
*/
	double Q = pow( parentHadronMass , 2 ) - pow( visible_mass , 2 ) - 2.0 * pow( P_vis_nor , 2 );
	streamlog_out(DEBUG8) << "		(Q) M_B^2 - M_vis^2 - 2 * P_vis(nor)^2 =	" << Q << std::endl;
	double A = pow( visible_E , 2 ) - pow( P_vis_par , 2 );
	streamlog_out(DEBUG8) << "		(A) E_vis^2 - P_vis(par)^2 =	" << A << std::endl;
	double B = Q * P_vis_par;
	streamlog_out(DEBUG8) << "		(B) ( M_B^2 - M_vis^2 - 2 * P_vis(nor)^2 ) * P_vis(par) =	" << B << std::endl;
	double C = pow( P_vis_nor , 2 ) * pow( visible_E , 2 ) - pow( Q , 2 ) / 4.0;
	streamlog_out(DEBUG8) << "		(C) E_vis^2 * P_vis(nor)^2 - ( M_B^2 - M_vis^2 - 2 * P_vis(nor)^2 )^2 / 4 =	" << C << std::endl;

	double P_nu_par = ( pow( B , 2 ) >= 4.0 * A * C ? ( B + solutionSign * std::sqrt( pow( B , 2 ) - 4.0 * A * C ) ) / ( 2.0 * A ) : B / ( 2.0 * A ) );
	streamlog_out(DEBUG8) << "		|Neutrino Momentum (par)|=		" << P_nu_par << std::endl;
	TVector3 Neutrino_p_par		= P_nu_par * flightDirection;
	streamlog_out(DEBUG8) << "		Neutrino Momentum (par):		( " << Neutrino_p_par.Px() << "	, " << Neutrino_p_par.Py() << "	, " << Neutrino_p_par.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p_nor		= -1.0 * visible_p_nor;
	streamlog_out(DEBUG8) << "		Neutrino Momentum (nor):		( " << Neutrino_p_nor.Px() << "	, " << Neutrino_p_nor.Py() << "	, " << Neutrino_p_nor.Pz() << "	)" << std::endl;
	TVector3 Neutrino_p		= Neutrino_p_nor + Neutrino_p_par;
	streamlog_out(DEBUG8) << "		Neutrino Momentum:			( " << Neutrino_p.Px() << "	, " << Neutrino_p.Py() << "	, " << Neutrino_p.Pz() << "	)" << std::endl;
	double Neutrino_E		= Neutrino_p.Mag();
	streamlog_out(DEBUG8) << "		Neutrino Energy:									" << Neutrino_E << std::endl;

	TLorentzVector Neutrino_tlv( Neutrino_p , Neutrino_E );
	streamlog_out(DEBUG4) << "" << std::endl;
	streamlog_out(DEBUG4) << "	Reconstructed Neutrino Mass:	" << Neutrino_tlv.Mag() << std::endl;
	streamlog_out(DEBUG4) << "	Reconstructed Neutrino 4-Momentum:		( " << Neutrino_tlv.Px() << "	, " << Neutrino_tlv.Py() << "	, " << Neutrino_tlv.Pz() << "	, " << Neutrino_tlv.E() << " )" << std::endl;
	return Neutrino_tlv;
}

MCParticle* SLDCorrection::getTrueNeutrino( MCParticle *SLDLepton )
{
	MCParticle* trueNeutrino{};
	try
	{
		EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		int nNeutrinos = 0;
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( MotherHadron->getDaughters() ).size() ; ++i_daughter )
		{
			EVENT::MCParticle *daughter = MotherHadron->getDaughters()[ i_daughter ];
			if ( daughter->getGeneratorStatus() == 1 && ( abs( daughter->getPDG() ) == abs( SLDLepton->getPDG() ) + 1 ) )
			{
				streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
				streamlog_out(DEBUG0) << "------------------------------ Neutrino ------------------------------" << std::endl;
				streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
				streamlog_out(DEBUG0) << *daughter << std::endl;
				trueNeutrino = daughter;
			}
		}
		++nNeutrinos;
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	True Neutrino for semi-leptonic decay not found in MCParticles" << std::endl;
        }
	return trueNeutrino;
}

void SLDCorrection::fillTrueRecoFourMomentum(	TLorentzVector trueVisibleFourMomentumAtSLDVertex , TLorentzVector truePVATrueFourMomentum , TLorentzVector truePVARecoFourMomentum , TLorentzVector recoPVARecoFourMomentum , TLorentzVector visibleFourMomentum ,
												TLorentzVector trueLeptonFourMomentum , TLorentzVector recoLeptonFourMomentum , TLorentzVector leptonFourMomentum ,
												TLorentzVector truePVATrueChargedFourMomentum , TLorentzVector truePVARecoChargedFourMomentum , TLorentzVector recoSLDVertexChargedFourMomentum , TLorentzVector recoPVARecoChargedFourMomentum , TLorentzVector chargedFourMomentum ,
												TLorentzVector truePVATrueNeutralFourMomentum , TLorentzVector truePVARecoNeutralFourMomentum , TLorentzVector recoPVARecoNeutralFourMomentum , TLorentzVector neutralFourMomentum ,
												TLorentzVector trueNeutrinoFourMomentum , TLorentzVector recoNeutrinoFourMomentumClose , TLorentzVector trueHadronFourMomentum , TLorentzVector recoHadronFourMomentum ,
												TVector3 trueFlightDirection , TVector3 recoFlightDirection , TVector3 flightDirection )
{
	m_leptonE_to_parentE.push_back( trueLeptonFourMomentum.E() / trueHadronFourMomentum.E() );
	m_otherChargedE_to_parentE.push_back( truePVATrueChargedFourMomentum.E() / trueHadronFourMomentum.E() );
	m_allChargedE_to_parentE.push_back( ( trueLeptonFourMomentum + truePVATrueChargedFourMomentum ).E() / trueHadronFourMomentum.E() );
	m_neutralE_to_parentE.push_back( truePVATrueNeutralFourMomentum.E() / trueHadronFourMomentum.E() );
	m_neutrino_to_parentE.push_back( trueNeutrinoFourMomentum.E() / trueHadronFourMomentum.E() );

	m_trueVisibleFourMomentumAtSLDVertex_Px.push_back( trueVisibleFourMomentumAtSLDVertex.Px() );
	m_trueVisibleFourMomentumAtSLDVertex_Py.push_back( trueVisibleFourMomentumAtSLDVertex.Py() );
	m_trueVisibleFourMomentumAtSLDVertex_Pz.push_back( trueVisibleFourMomentumAtSLDVertex.Pz() );
	m_trueVisibleFourMomentumAtSLDVertex_E.push_back( trueVisibleFourMomentumAtSLDVertex.E() );
	m_trueVisibleFourMomentumAtSLDVertex_M.push_back( trueVisibleFourMomentumAtSLDVertex.M() );
	m_truePVATrueFourMomentum_Px.push_back( truePVATrueFourMomentum.Px() );
	m_truePVATrueFourMomentum_Py.push_back( truePVATrueFourMomentum.Py() );
	m_truePVATrueFourMomentum_Pz.push_back( truePVATrueFourMomentum.Pz() );
	m_truePVATrueFourMomentum_E.push_back( truePVATrueFourMomentum.E() );
	m_truePVATrueFourMomentum_M.push_back( truePVATrueFourMomentum.M() );
	m_truePVARecoFourMomentum_Px.push_back( truePVARecoFourMomentum.Px() );
	m_truePVARecoFourMomentum_Py.push_back( truePVARecoFourMomentum.Py() );
	m_truePVARecoFourMomentum_Pz.push_back( truePVARecoFourMomentum.Pz() );
	m_truePVARecoFourMomentum_E.push_back( truePVARecoFourMomentum.E() );
	m_truePVARecoFourMomentum_M.push_back( truePVARecoFourMomentum.M() );
	m_recoPVARecoFourMomentum_Px.push_back( recoPVARecoFourMomentum.Px() );
	m_recoPVARecoFourMomentum_Py.push_back( recoPVARecoFourMomentum.Py() );
	m_recoPVARecoFourMomentum_Pz.push_back( recoPVARecoFourMomentum.Pz() );
	m_recoPVARecoFourMomentum_E.push_back( recoPVARecoFourMomentum.E() );
	m_recoPVARecoFourMomentum_M.push_back( recoPVARecoFourMomentum.M() );
	m_usedVisibleFourMomentum_Px.push_back( visibleFourMomentum.Px() );
	m_usedVisibleFourMomentum_Py.push_back( visibleFourMomentum.Py() );
	m_usedVisibleFourMomentum_Pz.push_back( visibleFourMomentum.Pz() );
	m_usedVisibleFourMomentum_E.push_back( visibleFourMomentum.E() );
	m_usedVisibleFourMomentum_M.push_back( visibleFourMomentum.M() );
	TVector3 recoPVARecoDirection = recoPVARecoFourMomentum.Vect(); recoPVARecoDirection.SetMag( 1.0 );
	TVector3 truePVARecoDirection = truePVARecoFourMomentum.Vect(); truePVARecoDirection.SetMag( 1.0 );
	m_PVAAlpha.push_back( acos( recoPVARecoDirection.Dot( truePVARecoDirection ) ) );
	m_PVASinAlpha.push_back( sin( acos( recoPVARecoDirection.Dot( truePVARecoDirection ) ) ) );
	m_PVACosAlpha.push_back( recoPVARecoDirection.Dot( truePVARecoDirection ) );


	m_trueLeptonFourMomentum_Px.push_back( trueLeptonFourMomentum.Px() );
	m_trueLeptonFourMomentum_Py.push_back( trueLeptonFourMomentum.Py() );
	m_trueLeptonFourMomentum_Pz.push_back( trueLeptonFourMomentum.Pz() );
	m_trueLeptonFourMomentum_E.push_back( trueLeptonFourMomentum.E() );
	m_trueLeptonFourMomentum_M.push_back( trueLeptonFourMomentum.M() );
	m_recoLeptonFourMomentum_Px.push_back( recoLeptonFourMomentum.Px() );
	m_recoLeptonFourMomentum_Py.push_back( recoLeptonFourMomentum.Py() );
	m_recoLeptonFourMomentum_Pz.push_back( recoLeptonFourMomentum.Pz() );
	m_recoLeptonFourMomentum_E.push_back( recoLeptonFourMomentum.E() );
	m_recoLeptonFourMomentum_M.push_back( recoLeptonFourMomentum.M() );
	m_usedLeptonFourMomentum_Px.push_back( leptonFourMomentum.Px() );
	m_usedLeptonFourMomentum_Py.push_back( leptonFourMomentum.Py() );
	m_usedLeptonFourMomentum_Pz.push_back( leptonFourMomentum.Pz() );
	m_usedLeptonFourMomentum_E.push_back( leptonFourMomentum.E() );
	m_usedLeptonFourMomentum_M.push_back( leptonFourMomentum.M() );

	m_truePVATrueChargedFourMomentum_Px.push_back( truePVATrueChargedFourMomentum.Px() );
	m_truePVATrueChargedFourMomentum_Py.push_back( truePVATrueChargedFourMomentum.Py() );
	m_truePVATrueChargedFourMomentum_Pz.push_back( truePVATrueChargedFourMomentum.Pz() );
	m_truePVATrueChargedFourMomentum_E.push_back( truePVATrueChargedFourMomentum.E() );
	m_truePVATrueChargedFourMomentum_M.push_back( truePVATrueChargedFourMomentum.M() );
	m_truePVARecoChargedFourMomentum_Px.push_back( truePVARecoChargedFourMomentum.Px() );
	m_truePVARecoChargedFourMomentum_Py.push_back( truePVARecoChargedFourMomentum.Py() );
	m_truePVARecoChargedFourMomentum_Pz.push_back( truePVARecoChargedFourMomentum.Pz() );
	m_truePVARecoChargedFourMomentum_E.push_back( truePVARecoChargedFourMomentum.E() );
	m_truePVARecoChargedFourMomentum_M.push_back( truePVARecoChargedFourMomentum.M() );
	m_recoSLDVertexChargedFourMomentum_Px.push_back( recoSLDVertexChargedFourMomentum.Px() );
	m_recoSLDVertexChargedFourMomentum_Py.push_back( recoSLDVertexChargedFourMomentum.Py() );
	m_recoSLDVertexChargedFourMomentum_Pz.push_back( recoSLDVertexChargedFourMomentum.Pz() );
	m_recoSLDVertexChargedFourMomentum_E.push_back( recoSLDVertexChargedFourMomentum.E() );
	m_recoSLDVertexChargedFourMomentum_M.push_back( recoSLDVertexChargedFourMomentum.M() );
	m_recoPVARecoChargedFourMomentum_Px.push_back( recoPVARecoChargedFourMomentum.Px() );
	m_recoPVARecoChargedFourMomentum_Py.push_back( recoPVARecoChargedFourMomentum.Py() );
	m_recoPVARecoChargedFourMomentum_Pz.push_back( recoPVARecoChargedFourMomentum.Pz() );
	m_recoPVARecoChargedFourMomentum_E.push_back( recoPVARecoChargedFourMomentum.E() );
	m_recoPVARecoChargedFourMomentum_M.push_back( recoPVARecoChargedFourMomentum.M() );
	m_usedChargedFourMomentum_Px.push_back( chargedFourMomentum.Px() );
	m_usedChargedFourMomentum_Py.push_back( chargedFourMomentum.Py() );
	m_usedChargedFourMomentum_Pz.push_back( chargedFourMomentum.Pz() );
	m_usedChargedFourMomentum_E.push_back( chargedFourMomentum.E() );
	m_usedChargedFourMomentum_M.push_back( chargedFourMomentum.M() );
	TVector3 recoPVARecoChargedDirection = ( recoSLDVertexChargedFourMomentum + recoPVARecoChargedFourMomentum ).Vect(); recoPVARecoChargedDirection.SetMag( 1.0 );
	TVector3 truePVARecoChargedDirection = truePVARecoChargedFourMomentum.Vect(); truePVARecoChargedDirection.SetMag( 1.0 );
	if ( truePVARecoChargedFourMomentum.E() > 0.0 && ( recoSLDVertexChargedFourMomentum + recoPVARecoChargedFourMomentum ).E() > 0.0 )
	{
		m_ChargedPVAAlpha.push_back( acos( recoPVARecoChargedDirection.Dot( truePVARecoChargedDirection ) ) );
		m_ChargedPVASinAlpha.push_back( sin( acos( recoPVARecoChargedDirection.Dot( truePVARecoChargedDirection ) ) ) );
		m_ChargedPVACosAlpha.push_back( recoPVARecoChargedDirection.Dot( truePVARecoChargedDirection ) );
	}
	else
	{
		m_ChargedPVAAlpha.push_back( -100.0 );
		m_ChargedPVASinAlpha.push_back( -100.0 );
		m_ChargedPVACosAlpha.push_back( -100.0 );
	}

	m_truePVATrueNeutralFourMomentum_Px.push_back( truePVATrueNeutralFourMomentum.Px() );
	m_truePVATrueNeutralFourMomentum_Py.push_back( truePVATrueNeutralFourMomentum.Py() );
	m_truePVATrueNeutralFourMomentum_Pz.push_back( truePVATrueNeutralFourMomentum.Pz() );
	m_truePVATrueNeutralFourMomentum_E.push_back( truePVATrueNeutralFourMomentum.E() );
	m_truePVATrueNeutralFourMomentum_M.push_back( truePVATrueNeutralFourMomentum.M() );
	m_truePVARecoNeutralFourMomentum_Px.push_back( truePVARecoNeutralFourMomentum.Px() );
	m_truePVARecoNeutralFourMomentum_Py.push_back( truePVARecoNeutralFourMomentum.Py() );
	m_truePVARecoNeutralFourMomentum_Pz.push_back( truePVARecoNeutralFourMomentum.Pz() );
	m_truePVARecoNeutralFourMomentum_E.push_back( truePVARecoNeutralFourMomentum.E() );
	m_truePVARecoNeutralFourMomentum_M.push_back( truePVARecoNeutralFourMomentum.M() );
	m_recoPVARecoNeutralFourMomentum_Px.push_back( recoPVARecoNeutralFourMomentum.Px() );
	m_recoPVARecoNeutralFourMomentum_Py.push_back( recoPVARecoNeutralFourMomentum.Py() );
	m_recoPVARecoNeutralFourMomentum_Pz.push_back( recoPVARecoNeutralFourMomentum.Pz() );
	m_recoPVARecoNeutralFourMomentum_E.push_back( recoPVARecoNeutralFourMomentum.E() );
	m_recoPVARecoNeutralFourMomentum_M.push_back( recoPVARecoNeutralFourMomentum.M() );
	m_usedNeutralFourMomentum_Px.push_back( neutralFourMomentum.Px() );
	m_usedNeutralFourMomentum_Py.push_back( neutralFourMomentum.Py() );
	m_usedNeutralFourMomentum_Pz.push_back( neutralFourMomentum.Pz() );
	m_usedNeutralFourMomentum_E.push_back( neutralFourMomentum.E() );
	m_usedNeutralFourMomentum_M.push_back( neutralFourMomentum.M() );

	TLorentzVector recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum = recoPVARecoFourMomentum - recoLeptonFourMomentum - truePVARecoNeutralFourMomentum;
	TLorentzVector recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum = recoPVARecoFourMomentum - recoLeptonFourMomentum - truePVARecoChargedFourMomentum;
	m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Px.push_back( recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum.Px() );
	m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Py.push_back( recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum.Py() );
	m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Pz.push_back( recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum.Pz() );
	m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_E.push_back( recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum.E() );
	m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_M.push_back( recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum.M() );
	m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Px.push_back( recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum.Px() );
	m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Py.push_back( recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum.Py() );
	m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Pz.push_back( recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum.Pz() );
	m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_E.push_back( recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum.E() );
	m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_M.push_back( recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum.M() );


	TVector3 recoPVARecoNeutralDirection = recoPVARecoNeutralFourMomentum.Vect(); recoPVARecoNeutralDirection.SetMag( 1.0 );
	TVector3 truePVARecoNeutralDirection = truePVARecoNeutralFourMomentum.Vect(); truePVARecoNeutralDirection.SetMag( 1.0 );
	if ( recoPVARecoNeutralFourMomentum.E() > 0.0 && recoPVARecoNeutralFourMomentum.E() > 0.0 )
	{
		m_NeutralPVAAlpha.push_back( acos( recoPVARecoNeutralDirection.Dot( truePVARecoNeutralDirection ) ) );
		m_NeutralPVASinAlpha.push_back( sin( acos( recoPVARecoNeutralDirection.Dot( truePVARecoNeutralDirection ) ) ) );
		m_NeutralPVACosAlpha.push_back( recoPVARecoNeutralDirection.Dot( truePVARecoNeutralDirection ) );
	}
	else
	{
		m_NeutralPVAAlpha.push_back( -100.0 );
		m_NeutralPVASinAlpha.push_back( -100.0 );
		m_NeutralPVACosAlpha.push_back( -100.0 );
	}
	TVector3 recoPVARecoMomentum_minus_truePVARecoNeutralMomentum = recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum.Vect();
	TVector3 recoPVARecoMomentum_minus_truePVARecoNeutralMomentum_Direction = recoPVARecoMomentum_minus_truePVARecoNeutralMomentum; recoPVARecoMomentum_minus_truePVARecoNeutralMomentum_Direction.SetMag( 1.0 );
	TVector3 recoPVARecoMomentum_minus_truePVARecoChargedMomentum = recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum.Vect();
	TVector3 recoPVARecoMomentum_minus_truePVARecoChargedMomentum_Direction = recoPVARecoMomentum_minus_truePVARecoChargedMomentum; recoPVARecoMomentum_minus_truePVARecoChargedMomentum_Direction.SetMag( 1.0 );

	m_Alpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged.push_back( acos( recoPVARecoMomentum_minus_truePVARecoNeutralMomentum_Direction.Dot( recoPVARecoChargedDirection ) ) );
	m_SinAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged.push_back( sin( acos( recoPVARecoMomentum_minus_truePVARecoNeutralMomentum_Direction.Dot( recoPVARecoChargedDirection ) ) ) );
	m_CosAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged.push_back( recoPVARecoMomentum_minus_truePVARecoNeutralMomentum_Direction.Dot( recoPVARecoChargedDirection ) );
	m_Alpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral.push_back( acos( recoPVARecoMomentum_minus_truePVARecoChargedMomentum_Direction.Dot( recoPVARecoNeutralDirection ) ) );
	m_SinAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral.push_back( sin( acos( recoPVARecoMomentum_minus_truePVARecoChargedMomentum_Direction.Dot( recoPVARecoNeutralDirection ) ) ) );
	m_CosAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral.push_back( recoPVARecoMomentum_minus_truePVARecoChargedMomentum_Direction.Dot( recoPVARecoNeutralDirection ) );

	m_trueNeutrinoFourMomentum_Px.push_back( trueNeutrinoFourMomentum.Px() );
	m_trueNeutrinoFourMomentum_Py.push_back( trueNeutrinoFourMomentum.Py() );
	m_trueNeutrinoFourMomentum_Pz.push_back( trueNeutrinoFourMomentum.Pz() );
	m_trueNeutrinoFourMomentum_E.push_back( trueNeutrinoFourMomentum.E() );
	m_trueNeutrinoFourMomentum_M.push_back( trueNeutrinoFourMomentum.M() );
	m_recoNeutrinoFourMomentumClose_Px.push_back( recoNeutrinoFourMomentumClose.Px() );
	m_recoNeutrinoFourMomentumClose_Py.push_back( recoNeutrinoFourMomentumClose.Py() );
	m_recoNeutrinoFourMomentumClose_Pz.push_back( recoNeutrinoFourMomentumClose.Pz() );
	m_recoNeutrinoFourMomentumClose_E.push_back( recoNeutrinoFourMomentumClose.E() );
	m_recoNeutrinoFourMomentumClose_M.push_back( recoNeutrinoFourMomentumClose.M() );
	m_trueHadronFourMomentum_Px.push_back( trueHadronFourMomentum.Px() );
	m_trueHadronFourMomentum_Py.push_back( trueHadronFourMomentum.Py() );
	m_trueHadronFourMomentum_Pz.push_back( trueHadronFourMomentum.Pz() );
	m_trueHadronFourMomentum_E.push_back( trueHadronFourMomentum.E() );
	m_trueHadronFourMomentum_M.push_back( trueHadronFourMomentum.M() );
	m_recoHadronFourMomentum_Px.push_back( recoHadronFourMomentum.Px() );
	m_recoHadronFourMomentum_Py.push_back( recoHadronFourMomentum.Py() );
	m_recoHadronFourMomentum_Pz.push_back( recoHadronFourMomentum.Pz() );
	m_recoHadronFourMomentum_E.push_back( recoHadronFourMomentum.E() );
	m_recoHadronFourMomentum_M.push_back( recoHadronFourMomentum.M() );

	m_trueFlightDirection_X.push_back( trueFlightDirection.X() );
	m_trueFlightDirection_Y.push_back( trueFlightDirection.Y() );
	m_trueFlightDirection_Z.push_back( trueFlightDirection.Z() );
	m_recoFlightDirection_X.push_back( recoFlightDirection.X() );
	m_recoFlightDirection_Y.push_back( recoFlightDirection.Y() );
	m_recoFlightDirection_Z.push_back( recoFlightDirection.Z() );
	m_usedFlightDirection_X.push_back( flightDirection.X() );
	m_usedFlightDirection_Y.push_back( flightDirection.Y() );
	m_usedFlightDirection_Z.push_back( flightDirection.Z() );
}

void SLDCorrection::investigateJetEnergyContent( EVENT::ReconstructedParticle *assignedJet )
{
	double jetEnergy = assignedJet->getEnergy();
	double neutralHadronEnergy = 0.0;
	double chargedEnergy = 0.0;
	double photonEnergy = 0.0;
	double neutralsEnergy = 0.0;
	for ( unsigned int i_par = 0 ; i_par < assignedJet->getParticles().size() ; ++i_par )
	{
		if ( ( assignedJet->getParticles()[ i_par ] )->getTracks().size() != 0 )
		{
			chargedEnergy += assignedJet->getParticles()[ i_par ]->getEnergy();
		}
		else if ( ( assignedJet->getParticles()[ i_par ] )->getType() == 22 )
		{
			photonEnergy += assignedJet->getParticles()[ i_par ]->getEnergy();
			neutralsEnergy += assignedJet->getParticles()[ i_par ]->getEnergy();
		}
		else
		{
			neutralHadronEnergy += assignedJet->getParticles()[ i_par ]->getEnergy();
			neutralsEnergy += assignedJet->getParticles()[ i_par ]->getEnergy();
		}
	}
	m_jetEnergyFractionCharged.push_back( chargedEnergy / jetEnergy );
	m_jetEnergyFractionNeutralHadron.push_back( neutralHadronEnergy / jetEnergy );
	m_jetEnergyFractionPhoton.push_back( photonEnergy / jetEnergy );
	m_jetEnergyFractionNeutrals.push_back( neutralsEnergy / jetEnergy );
	m_jetEnergy.push_back( jetEnergy );
}

void SLDCorrection::plotHistograms( TLorentzVector trueFourMomentumNeutrino , TLorentzVector FourMomentumNuClose , std::vector<float> NeutrinoCovMat )
{
	double NuPxResidual = FourMomentumNuClose.Px() - trueFourMomentumNeutrino.Px(); m_NuPxResidual.push_back( NuPxResidual );
	double NuPyResidual = FourMomentumNuClose.Py() - trueFourMomentumNeutrino.Py(); m_NuPyResidual.push_back( NuPyResidual );
	double NuPzResidual = FourMomentumNuClose.Pz() - trueFourMomentumNeutrino.Pz(); m_NuPzResidual.push_back( NuPzResidual );
	double NuEResidual = FourMomentumNuClose.E() - trueFourMomentumNeutrino.E(); m_NuEResidual.push_back( NuEResidual );
	double NuPxNormalizedResidual = NuPxResidual / sqrt( NeutrinoCovMat[ 0 ] ); m_NuPxNormalizedResidual.push_back( NuPxNormalizedResidual );
	double NuPyNormalizedResidual = NuPyResidual / sqrt( NeutrinoCovMat[ 2 ] ); m_NuPyNormalizedResidual.push_back( NuPyNormalizedResidual );
	double NuPzNormalizedResidual = NuPzResidual / sqrt( NeutrinoCovMat[ 5 ] ); m_NuPzNormalizedResidual.push_back( NuPzNormalizedResidual );
	double NuENormalizedResidual = NuEResidual / sqrt( NeutrinoCovMat[ 9 ] ); m_NuENormalizedResidual.push_back( NuENormalizedResidual );
	h_NuPxResidual->Fill( NuPxResidual ); ++n_NuPxResidual;
	h_NuPyResidual->Fill( NuPyResidual ); ++n_NuPyResidual;
	h_NuPzResidual->Fill( NuPzResidual ); ++n_NuPzResidual;
	h_NuEResidual->Fill( NuEResidual ); ++n_NuEResidual;
	h_NuPxNormalizedResidual->Fill( NuPxNormalizedResidual ); ++n_NuPxNormalizedResidual;
	h_NuPyNormalizedResidual->Fill( NuPyNormalizedResidual ); ++n_NuPyNormalizedResidual;
	h_NuPzNormalizedResidual->Fill( NuPzNormalizedResidual ); ++n_NuPzNormalizedResidual;
	h_NuENormalizedResidual->Fill( NuENormalizedResidual ); ++n_NuENormalizedResidual;
	h_recoNuPx_mcNuPx->Fill( trueFourMomentumNeutrino.Px() , FourMomentumNuClose.Px() );
	h_recoNuPy_mcNuPy->Fill( trueFourMomentumNeutrino.Py() , FourMomentumNuClose.Py() );
	h_recoNuPz_mcNuPz->Fill( trueFourMomentumNeutrino.Pz() , FourMomentumNuClose.Pz() );
	h_recoNuE_mcNuE->Fill( trueFourMomentumNeutrino.E() , FourMomentumNuClose.E() );
}

void SLDCorrection::InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle )
{
	histogram->Scale( 1.0 / scale );
	histogram->SetLineColor( color );
	histogram->SetLineWidth( lineWidth );
	histogram->SetMarkerSize( markerSize );
	histogram->SetMarkerStyle( markerStyle );
	histogram->SetMarkerColor( color );
//	float fit_range = 4.0;
//	float fit_min = -2.0;
//	float fit_max = 2.0;
//	histogram->GetFunction("gaus")->SetLineColor( color );
	float y_max = 1.2 * histogram->GetMaximum();
	histogram->GetYaxis()->SetRangeUser(0.0, y_max);
	histogram->Write();
}

void SLDCorrection::check( EVENT::LCEvent *pLCEvent )
{
	try
	{
		streamlog_out(DEBUG4) << "	CHECK : processed events: " << m_nEvtSum << " , Number of semi-leptonic vertex: " << pLCEvent->getCollection(m_SLDVertex)->getNumberOfElements() << std::endl;
		streamlog_out(DEBUG4) << "	CHECK : processed events: " << m_nEvtSum << " , Number of semi-leptonic vertexRP: " << pLCEvent->getCollection(m_SLDVertexRP)->getNumberOfElements() << std::endl;
		streamlog_out(DEBUG4) << "	CHECK : processed events: " << m_nEvtSum << " , Number of Reconstrcuted Neutrinos: " << pLCEvent->getCollection(m_reconstructedNeutrino)->getNumberOfElements() << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "	Input/Output collection not found in event " << m_nEvt << std::endl;
        }
}

void SLDCorrection::end()
{
	if ( m_fillRootTree )
	{
		m_pTFile->cd();
		m_pTTree1->Write();
		m_pTTree2->Write();
		m_pTTree3->Write();
//		h_SLDStatus->Scale( 100.0 / n_SLDStatus );
		h_SLDStatus->GetYaxis()->SetTitle("number of SLDecays");
		h_SLDStatus->Write();
		h_BHadronType->Write();
		h_CHadronType->Write();
		InitializeHistogram( h_NuPxResidual , n_NuPxResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPyResidual , n_NuPyResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPzResidual , n_NuPzResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuEResidual , n_NuEResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPxNormalizedResidual , n_NuPxNormalizedResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPyNormalizedResidual , n_NuPyNormalizedResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPzNormalizedResidual , n_NuPzNormalizedResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuENormalizedResidual , n_NuENormalizedResidual , 4 , 1 , 1.0 , 1 );
		h_SLDecayFlavour->Write();
		h_SLDecayModeB->Write();
		h_SLDecayModeC->Write();
		h_SLDecayOrder->Write();
		h_NuPxResidual->Write();
		h_NuPyResidual->Write();
		h_NuPzResidual->Write();
		h_NuEResidual->Write();
		h_NuPxNormalizedResidual->Write();
		h_NuPyNormalizedResidual->Write();
		h_NuPzNormalizedResidual->Write();
		h_NuENormalizedResidual->Write();
		h_recoNuPx_mcNuPx->Write();
		h_recoNuPy_mcNuPy->Write();
		h_recoNuPz_mcNuPz->Write();
		h_recoNuE_mcNuE->Write();
		h_secondaryVertex->Scale( 100.0 / n_secondaryVertex );
		h_secondaryVertex->GetYaxis()->SetTitle("#SLDecay [%]");
		h_secondaryVertex->Write();
		m_pTFile->Close();
		delete m_pTFile;
	}
}
