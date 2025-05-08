#ifndef SLDCorrection_h
#define SLDCorrection_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "EVENT/LCStrVec.h"
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include "IMPL/LCCollectionVec.h"
#include "UTIL/LCRelationNavigator.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "IMPL/ParticleIDImpl.h"
#include <IMPL/VertexImpl.h>
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "DD4hep/Detector.h"
#include "SLDCorrectionTypes.h"
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TTree.h>
#include <TF1.h>

using namespace lcio ;
using namespace marlin ;

class SLDCorrection : public Processor
{
public:
	virtual Processor *newProcessor()
	{
		return new SLDCorrection;
	}
	SLDCorrection();
	virtual ~SLDCorrection() = default;
	SLDCorrection( const SLDCorrection& ) = delete;
	SLDCorrection &operator = ( const SLDCorrection& ) = delete;
	virtual void init();
	virtual void Clear();
	virtual void processRunHeader();
	virtual void processEvent( EVENT::LCEvent *pLCEvent );

	//	hasPrimarySLDecay checks if a MCParticle (potentially a B-/C-Hadron) decays semi-leptonically) true: decays semi-leptonically, false: doesn't decay semi-leptonically
	bool hasPrimarySLDecay( const MCP &parentHadron , int &chargedLeptonPDG );

	//	hasDownStreamSLDecay checks whether a B-/C-Hadron that decays semi-leptonically, has another semi-leptonic decay in its decay products or not, true: has another semi-leptonic decay, false: does not have any semi-leptonic decay in its decay products
	bool hasDownStreamSLDecay( const MCP &parentHadron );

	//	hasUpStreamSLDecay checks  whether a B-/C-Hadron that decays semi-leptonically, has another semi-leptonic decay in its up-stream, true has a semi-leptonic decay in up-stream, false, does not have a semi-leptonic decay in up-stream
	bool hasUpStreamSLDecay( const MCP &parentHadron );

	//	checkBHadronSLDecay checks the flavour of the Hadron decays semi-leptonicall, if B-Hadron: true, if not: false
	bool checkBHadronSLDecay( const MCP &parentHadron );

	//	checkCHadronSLDecay checks the flavour of the Hadron decays semi-leptonicall, if C-Hadron: true, if not: false
	bool checkCHadronSLDecay( const MCP &parentHadron );

	//	checkTauLeptonSLDecay checks the flavour of the Hadron decays semi-leptonicall, if tau-lepton: true, if not: false
	bool checkTauLeptonSLDecay( const MCP &parentHadron );

	//	doSLDCorrection prepares all inputs for nu-correction and makes output LCIO elements (SLDVertex, associated RP of SLDVertex, solutions for neutrino as Reconstructed Particle, etc)
	virtual void doSLDCorrection( EVENT::LCEvent *pLCEvent , const MCP &SLDLepton , size_t parentHadronIDx, VertexVector& semiLeptonicVertices , PFOVector& semiLeptonicVertexRecoParticles , PFOVector& jetsOfSemiLeptonicDecays , PFOVectorVector& neutrinos , IntVector &sldStatus , IntVector &pvaStatus , IntVector &solutionSigns , MCPVector &trueNeutrinos );

	//	showTrueParameters prints true input for nu-correction
	void showTrueParameters( const MCP &SLDLepton , size_t parentHadronIDx );

	//	getNeutrinoFourMomentum corrects neutrino energy using rapidity-based approach
	TLorentzVector getNeutrinoFourMomentum( const TVector3 &flightDirection , const TLorentzVector &visibleFourMomentum , const double &parentHadronMass , const float &solutionSign );

	//	getNeutrinoFourMomentumModified corrects neutrino energy using rapidity-based approach (checks P_vis_nor, P_vis_par, flightDirection, etc to avoid math divergence)
	TLorentzVector getNeutrinoFourMomentumModified( TVector3 &flightDirection , const TLorentzVector &visibleFourMomentum , const double &parentHadronMass , const float &solutionSign );

	//	getNeutrinoFourMomentumStandardMethod corrects neutrino energy satndard (p,E) conservation method (highly affected by numeric precision!)
	TLorentzVector getNeutrinoFourMomentumStandardMethod( const TVector3 &flightDirection , const TLorentzVector &visibleFourMomentum , const double &parentHadronMass , const float &solutionSign );

	//	getTrueNeutrinogets true four-momentum of neutrino
	MCP getTrueNeutrino( const MCP &SLDLepton, size_t parent_idx );

	//	fillTrueRecoFourMomentummakes performance evaluation variables stored in root tree
	void fillTrueRecoFourMomentum( const TLorentzVector &trueVisibleFourMomentumAtSLDVertex , const TLorentzVector &truePVATrueFourMomentum , const TLorentzVector &truePVARecoFourMomentum , const TLorentzVector &recoPVARecoFourMomentum , const TLorentzVector &visibleFourMomentum , const TLorentzVector &trueLeptonFourMomentum , const TLorentzVector &recoLeptonFourMomentum , const TLorentzVector &leptonFourMomentum , const TLorentzVector &truePVATrueChargedFourMomentum , const TLorentzVector &truePVARecoChargedFourMomentum , const TLorentzVector &recoSLDVertexChargedFourMomentum , const TLorentzVector &recoPVARecoChargedFourMomentum , const TLorentzVector &chargedFourMomentum , const TLorentzVector &truePVATrueNeutralFourMomentum , const TLorentzVector &truePVARecoNeutralFourMomentum , const TLorentzVector &recoPVARecoNeutralFourMomentum , const TLorentzVector &neutralFourMomentum , const TLorentzVector &trueNeutrinoFourMomentum , const TLorentzVector &recoNeutrinoFourMomentumClose , const TLorentzVector &trueHadronFourMomentum , const TLorentzVector &recoHadronFourMomentum , const TVector3 &trueFlightDirection , const TVector3 &recoFlightDirection , const TVector3 &flightDirection );

	//	addNeutrinoCovarianceMatrix Calculates Covariance matrix for each solution of neutrinos, based and energy and angular error
	virtual void addNeutrinoCovarianceMatrix( const TLorentzVector &neutrinoFourMomentum , std::vector< float > &NuCovMat );

	virtual void evaluateInputCovMat( const TLorentzVector &trueVisibleFourMomentum , const TVector3 &trueFlightDirection , const TLorentzVector &trueNeutrinoFourMomentum , const TLorentzVector &visibleFourMomentum , const TVector3 &flightDirection , const TLorentzVector &recoNeutrinoFourMomentum , const FloatVector &CovMatDetector , const FloatVector &CovMatNeutrino );

	virtual void plotHistograms( const TLorentzVector &trueFourMomentumNeutrino , const TLorentzVector &FourMomentumNuClose , const FloatVector &NeutrinoCovMat );

	virtual void InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle );

	void investigateJetEnergyContent( const RecoParticle &assignedJet );

	void checkSLDInput( const MCP &SLDHadron );

	virtual void check( EVENT::LCEvent *pLCEvent );

	virtual void end();

	dd4hep::Detector& _theDetector = dd4hep::Detector::getInstance();

private:



	std::string				m_mcParticleCollection{};
	std::string				m_inputPfoCollection{};
	std::string				m_inputJetCollection{};
	std::string				m_inputPrimaryVertex{};
	std::string				m_inputBuildUpVertex{};
	std::string				m_RecoMCTruthLinkCollection{};
	std::string				m_MCTruthRecoLinkCollection{};
	std::string				m_ClusterMCTruthLinkCollection{};
	std::string				m_MCTruthClusterLinkCollection{};
	std::string				m_SLDVertex{};
	std::string				m_SLDVertexRP{};
	std::string				m_reconstructedNeutrino{};
	std::string				m_JetSLDLinkName{};
	std::string				m_SLDJetLinkName{};
	std::string				m_mcNurecoNuLinkName{};
	std::string				m_recoNumcNuLinkName{};
	std::string				m_SLDNuLinkName{};
	std::string				m_NuSLDLinkName{};
	std::string				m_rootFile{};

	bool					m_includeBSLD = true;
	bool					m_includeCSLD = true;
	bool					m_includeTSLD = true;
	bool					m_cheatSLDLeptons = true;
	bool					m_cheatFlightDirection = true;
	int					m_vertexingScenario = 1;
	bool					m_cheatLepton4momentum = true;
	bool					m_cheatCharged4momentum = true;
	bool					m_cheatNeutral4momentum = true;
	bool					m_cheatPVAcharged = true;
	float					m_chargedCosAcceptanceAngleSLD4 = 0.0;
	float					m_chargedCosAcceptanceAngleSLD5 = 0.0;
	bool					m_cheatPVAneutral = true;
	bool					m_cheatSolutionSign = false;
	float					m_neutralCosAcceptanceAngle = 0.0;
	float					m_BSLDChargedSLD4InvMassCut = 0.0;
	float					m_BSLDChargedSLD5InvMassCut = 0.0;
	float					m_BSLDNeutralSLD4InvMassCut = 0.0;
	float					m_BSLDNeutralSLD5InvMassCut = 0.0;
	float					m_CSLDChargedSLD4InvMassCut = 0.0;
	float					m_CSLDChargedSLD5InvMassCut = 0.0;
	float					m_CSLDNeutralSLD4InvMassCut = 0.0;
	float					m_CSLDNeutralSLD5InvMassCut = 0.0;
	int					m_nIterFlightDirCorrection = 0;
	float					m_BSLD4SigmaAlpha = 0.004;
	float					m_BSLD5SigmaAlpha = 0.010;
	float					m_BSLD4SigmaECPVA = 0;
	float					m_BSLD4SigmaENPVA = 0;
	float					m_BSLD4SigmaAlphaCPVA = 0.0;
	float					m_BSLD4SigmaAlphaNPVA = 0.0;
	float					m_BSLD5SigmaECPVA = 0;
	float					m_BSLD5SigmaENPVA = 0;
	float					m_BSLD5SigmaAlphaCPVA = 0.0;
	float					m_BSLD5SigmaAlphaNPVA = 0.0;
	float					m_sigmaAlphaNu = 0.1;
	float					m_sigmaENu = 4.0;
	bool					m_displayEvent = true;
	bool					m_fillRootTree = true;
	bool					m_traceEvent = false;
	int					m_BSLDMode = 0;
	int					m_CSLDMode = 0;
	int					m_TSLDMode = 0;
	int					m_SLDMode = 0;

	int					m_nRun;
	int					m_nEvt;
	int					m_nRunSum;
	int					m_nEvtSum;
	double					m_Bfield;
	bool					foundFlightDirection;
	IntVector				m_SLDFlavour{}; //4: SLDecayOfCHadron, 5: SLDecayOfBHadron, 15: SLDecayOfTauLepton
	IntVector				m_SLDType{}; //0: SLDecay with DownStream/UpStream semi-leptonic decay(s), 1: SLDecay without DownStream/UpStream semi-leptonic decay(s)
	IntVector				m_SLDLeptonID{}; //+11/-11: SLDecay to electron/positron, +13/-13: SLDecay to muon/anti-muon, +15/-15: SLDecay to tau/anti-tau
	FloatVector				m_leptonE_to_parentE{};
	FloatVector				m_otherChargedE_to_parentE{};
	FloatVector				m_allChargedE_to_parentE{};
	FloatVector				m_neutralE_to_parentE{};
	FloatVector				m_neutrino_to_parentE{};
	unsigned int				m_nSLDecayOfBHadron;
	unsigned int				m_nSLDecayOfCHadron;
	unsigned int				m_nSLDecayOfTauLepton;
	unsigned int				m_nSLDecayTotal;
	int					m_nSLDecayToElectron;
	int					m_nSLDecayToMuon;
	int					m_nSLDecayToTau;
	int					m_nTauNeutrino;
	IntVector				BHadPDGs{};
	IntVector				CHadPDGs{};
	IntVector				m_nRecoVerticesInJet{};
	DoubleVector				m_jetEnergy{};
	DoubleVector				m_jetEnergyFractionCharged{};
	DoubleVector				m_jetEnergyFractionNeutralHadron{};
	DoubleVector				m_jetEnergyFractionPhoton{};
	DoubleVector				m_jetEnergyFractionNeutrals{};
	DoubleVector				m_SLDecayXi{};
	DoubleVector				m_SLDecayYi{};
	DoubleVector				m_SLDecayZi{};
	DoubleVector				m_SLDecayRi{};
	DoubleVector				m_SLDecayXf{};
	DoubleVector				m_SLDecayYf{};
	DoubleVector				m_SLDecayZf{};
	DoubleVector				m_SLDecayRf{};
	DoubleVector				m_trueNuPx{};
	DoubleVector				m_trueNuPy{};
	DoubleVector				m_trueNuPz{};
	DoubleVector				m_trueNuE{};
	DoubleVector				m_recoNuClosePx{};
	DoubleVector				m_recoNuClosePy{};
	DoubleVector				m_recoNuClosePz{};
	DoubleVector				m_recoNuCloseE{};
	DoubleVector				m_recoNuPosPx{};
	DoubleVector				m_recoNuPosPy{};
	DoubleVector				m_recoNuPosPz{};
	DoubleVector				m_recoNuPosE{};
	DoubleVector				m_recoNuNegPx{};
	DoubleVector				m_recoNuNegPy{};
	DoubleVector				m_recoNuNegPz{};
	DoubleVector				m_recoNuNegE{};
	DoubleVector				m_NuPxResidual{};
	DoubleVector				m_NuPyResidual{};
	DoubleVector				m_NuPzResidual{};
	DoubleVector				m_NuEResidual{};
	DoubleVector				m_NuPxNormalizedResidual{};
	DoubleVector				m_NuPyNormalizedResidual{};
	DoubleVector				m_NuPzNormalizedResidual{};
	DoubleVector				m_NuENormalizedResidual{};
	IntVector				m_solutionSign{};
	DoubleVector				m_true_E_vis{};
	DoubleVector				m_true_E_vis_prime{};
	DoubleVector				m_true_P_vis_par{};
	DoubleVector				m_true_P_vis_par_prime{};
	DoubleVector				m_true_P_vis_nor{};
	DoubleVector				m_E_vis{};
	DoubleVector				m_E_vis_prime{};
	DoubleVector				m_P_vis_par{};
	DoubleVector				m_P_vis_par_prime{};
	DoubleVector				m_P_vis_par_prime_squared{};
	DoubleVector				m_P_vis_nor{};
	DoubleVector				m_P_vis_nor_prime{};
	IntVector				m_flightDirectionStatus{};
	DoubleVector				m_FlightDirectionErrorSinAlpha{};
	DoubleVector				m_FlightDirectionErrorCosAlpha{};
	DoubleVector				m_FlightDirectionErrorAlpha{};
	DoubleVector				m_distRecoLeptonToDownStreamVertex{};
	DoubleVector				m_parentHadronMass{};
	IntVector				m_parentHadronPDG{};
	DoubleVector				m_trueParentHadronFlightDistance{};
	DoubleVector				m_recoParentHadronFlightDistance{};
	DoubleVector				m_daughterHadronMass{};
	IntVector				m_daughterHadronPDG{};
	DoubleVector				m_daughterHadronFlightDistance{};
	FloatVector				m_visibleChargedInvMassCut{};
	FloatVector				m_visibleNeutralInvMassCut{};

	DoubleVector				m_trueVisibleFourMomentumAtSLDVertex_Px{};
	DoubleVector				m_trueVisibleFourMomentumAtSLDVertex_Py{};
	DoubleVector				m_trueVisibleFourMomentumAtSLDVertex_Pz{};
	DoubleVector				m_trueVisibleFourMomentumAtSLDVertex_E{};
	DoubleVector				m_trueVisibleFourMomentumAtSLDVertex_M{};
	DoubleVector				m_truePVATrueFourMomentum_Px{};
	DoubleVector				m_truePVATrueFourMomentum_Py{};
	DoubleVector				m_truePVATrueFourMomentum_Pz{};
	DoubleVector				m_truePVATrueFourMomentum_E{};
	DoubleVector				m_truePVATrueFourMomentum_M{};
	DoubleVector				m_truePVARecoFourMomentum_Px{};
	DoubleVector				m_truePVARecoFourMomentum_Py{};
	DoubleVector				m_truePVARecoFourMomentum_Pz{};
	DoubleVector				m_truePVARecoFourMomentum_E{};
	DoubleVector				m_truePVARecoFourMomentum_M{};
	DoubleVector				m_recoPVARecoFourMomentum_Px{};
	DoubleVector				m_recoPVARecoFourMomentum_Py{};
	DoubleVector				m_recoPVARecoFourMomentum_Pz{};
	DoubleVector				m_recoPVARecoFourMomentum_E{};
	DoubleVector				m_recoPVARecoFourMomentum_M{};
	DoubleVector				m_usedVisibleFourMomentum_Px{};
	DoubleVector				m_usedVisibleFourMomentum_Py{};
	DoubleVector				m_usedVisibleFourMomentum_Pz{};
	DoubleVector				m_usedVisibleFourMomentum_E{};
	DoubleVector				m_usedVisibleFourMomentum_M{};
	DoubleVector				m_trueLeptonFourMomentum_Px{};
	DoubleVector				m_trueLeptonFourMomentum_Py{};
	DoubleVector				m_trueLeptonFourMomentum_Pz{};
	DoubleVector				m_trueLeptonFourMomentum_E{};
	DoubleVector				m_trueLeptonFourMomentum_M{};
	DoubleVector				m_recoLeptonFourMomentum_Px{};
	DoubleVector				m_recoLeptonFourMomentum_Py{};
	DoubleVector				m_recoLeptonFourMomentum_Pz{};
	DoubleVector				m_recoLeptonFourMomentum_E{};
	DoubleVector				m_recoLeptonFourMomentum_M{};
	DoubleVector				m_usedLeptonFourMomentum_Px{};
	DoubleVector				m_usedLeptonFourMomentum_Py{};
	DoubleVector				m_usedLeptonFourMomentum_Pz{};
	DoubleVector				m_usedLeptonFourMomentum_E{};
	DoubleVector				m_usedLeptonFourMomentum_M{};
	DoubleVector				m_truePVATrueChargedFourMomentum_Px{};
	DoubleVector				m_truePVATrueChargedFourMomentum_Py{};
	DoubleVector				m_truePVATrueChargedFourMomentum_Pz{};
	DoubleVector				m_truePVATrueChargedFourMomentum_E{};
	DoubleVector				m_truePVATrueChargedFourMomentum_M{};
	DoubleVector				m_truePVARecoChargedFourMomentum_Px{};
	DoubleVector				m_truePVARecoChargedFourMomentum_Py{};
	DoubleVector				m_truePVARecoChargedFourMomentum_Pz{};
	DoubleVector				m_truePVARecoChargedFourMomentum_E{};
	DoubleVector				m_truePVARecoChargedFourMomentum_M{};
	DoubleVector				m_recoSLDVertexChargedFourMomentum_Px{};
	DoubleVector				m_recoSLDVertexChargedFourMomentum_Py{};
	DoubleVector				m_recoSLDVertexChargedFourMomentum_Pz{};
	DoubleVector				m_recoSLDVertexChargedFourMomentum_E{};
	DoubleVector				m_recoSLDVertexChargedFourMomentum_M{};
	DoubleVector				m_recoPVARecoChargedFourMomentum_Px{};
	DoubleVector				m_recoPVARecoChargedFourMomentum_Py{};
	DoubleVector				m_recoPVARecoChargedFourMomentum_Pz{};
	DoubleVector				m_recoPVARecoChargedFourMomentum_E{};
	DoubleVector				m_recoPVARecoChargedFourMomentum_M{};
	DoubleVector				m_usedChargedFourMomentum_Px{};
	DoubleVector				m_usedChargedFourMomentum_Py{};
	DoubleVector				m_usedChargedFourMomentum_Pz{};
	DoubleVector				m_usedChargedFourMomentum_E{};
	DoubleVector				m_usedChargedFourMomentum_M{};
	DoubleVector				m_truePVATrueNeutralFourMomentum_Px{};
	DoubleVector				m_truePVATrueNeutralFourMomentum_Py{};
	DoubleVector				m_truePVATrueNeutralFourMomentum_Pz{};
	DoubleVector				m_truePVATrueNeutralFourMomentum_E{};
	DoubleVector				m_truePVATrueNeutralFourMomentum_M{};
	DoubleVector				m_truePVARecoNeutralFourMomentum_Px{};
	DoubleVector				m_truePVARecoNeutralFourMomentum_Py{};
	DoubleVector				m_truePVARecoNeutralFourMomentum_Pz{};
	DoubleVector				m_truePVARecoNeutralFourMomentum_E{};
	DoubleVector				m_truePVARecoNeutralFourMomentum_M{};
	DoubleVector				m_recoPVARecoNeutralFourMomentum_Px{};
	DoubleVector				m_recoPVARecoNeutralFourMomentum_Py{};
	DoubleVector				m_recoPVARecoNeutralFourMomentum_Pz{};
	DoubleVector				m_recoPVARecoNeutralFourMomentum_E{};
	DoubleVector				m_recoPVARecoNeutralFourMomentum_M{};
	DoubleVector				m_usedNeutralFourMomentum_Px{};
	DoubleVector				m_usedNeutralFourMomentum_Py{};
	DoubleVector				m_usedNeutralFourMomentum_Pz{};
	DoubleVector				m_usedNeutralFourMomentum_E{};
	DoubleVector				m_usedNeutralFourMomentum_M{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Px{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Py{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Pz{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_E{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_M{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Px{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Py{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Pz{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_E{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_M{};
	DoubleVector				m_trueNeutrinoFourMomentum_Px{};
	DoubleVector				m_trueNeutrinoFourMomentum_Py{};
	DoubleVector				m_trueNeutrinoFourMomentum_Pz{};
	DoubleVector				m_trueNeutrinoFourMomentum_E{};
	DoubleVector				m_trueNeutrinoFourMomentum_M{};
	DoubleVector				m_recoNeutrinoFourMomentumClose_Px{};
	DoubleVector				m_recoNeutrinoFourMomentumClose_Py{};
	DoubleVector				m_recoNeutrinoFourMomentumClose_Pz{};
	DoubleVector				m_recoNeutrinoFourMomentumClose_E{};
	DoubleVector				m_recoNeutrinoFourMomentumClose_M{};
	DoubleVector				m_trueHadronFourMomentum_Px{};
	DoubleVector				m_trueHadronFourMomentum_Py{};
	DoubleVector				m_trueHadronFourMomentum_Pz{};
	DoubleVector				m_trueHadronFourMomentum_E{};
	DoubleVector				m_trueHadronFourMomentum_M{};
	DoubleVector				m_recoHadronFourMomentum_Px{};
	DoubleVector				m_recoHadronFourMomentum_Py{};
	DoubleVector				m_recoHadronFourMomentum_Pz{};
	DoubleVector				m_recoHadronFourMomentum_E{};
	DoubleVector				m_recoHadronFourMomentum_M{};
	DoubleVector				m_trueFlightDirection_X{};
	DoubleVector				m_trueFlightDirection_Y{};
	DoubleVector				m_trueFlightDirection_Z{};
	DoubleVector				m_recoFlightDirection_X{};
	DoubleVector				m_recoFlightDirection_Y{};
	DoubleVector				m_recoFlightDirection_Z{};
	DoubleVector				m_usedFlightDirection_X{};
	DoubleVector				m_usedFlightDirection_Y{};
	DoubleVector				m_usedFlightDirection_Z{};

	DoubleVector				m_NeutralPVAAlpha{};
	DoubleVector				m_NeutralPVASinAlpha{};
	DoubleVector				m_NeutralPVACosAlpha{};
	DoubleVector				m_ChargedPVAAlpha{};
	DoubleVector				m_ChargedPVASinAlpha{};
	DoubleVector				m_ChargedPVACosAlpha{};
	DoubleVector				m_PVAAlpha{};
	DoubleVector				m_PVASinAlpha{};
	DoubleVector				m_PVACosAlpha{};
	DoubleVector				m_Alpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral{};
	DoubleVector				m_SinAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral{};
	DoubleVector				m_CosAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral{};
	DoubleVector				m_Alpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged{};
	DoubleVector				m_SinAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged{};
	DoubleVector				m_CosAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged{};

	DoubleVector				m_visiblePxBeforePVA{};
	DoubleVector				m_visiblePyBeforePVA{};
	DoubleVector				m_visiblePzBeforePVA{};
	DoubleVector				m_visibleEnergyBeforePVA{};
	DoubleVector				m_visibleMassBeforePVA{};
	DoubleVector				m_visiblePxAfterVertexPVA{};
	DoubleVector				m_visiblePyAfterVertexPVA{};
	DoubleVector				m_visiblePzAfterVertexPVA{};
	DoubleVector				m_visibleEnergyAfterVertexPVA{};
	DoubleVector				m_visibleMassAfterVertexPVA{};
	DoubleVector				m_visiblePxAfterChargedPVA{};
	DoubleVector				m_visiblePyAfterChargedPVA{};
	DoubleVector				m_visiblePzAfterChargedPVA{};
	DoubleVector				m_visibleEnergyAfterChargedPVA{};
	DoubleVector				m_visibleMassAfterChargedPVA{};
	DoubleVector				m_visiblePxAfterNeutralPVA{};
	DoubleVector				m_visiblePyAfterNeutralPVA{};
	DoubleVector				m_visiblePzAfterNeutralPVA{};
	DoubleVector				m_visibleEnergyAfterNeutralPVA{};
	DoubleVector				m_visibleMassAfterNeutralPVA{};

	DoubleVector				m_trueVisibleFourMomentumPx{};
	DoubleVector				m_trueVisibleFourMomentumPy{};
	DoubleVector				m_trueVisibleFourMomentumPz{};
	DoubleVector				m_trueVisibleFourMomentumE{};
	DoubleVector				m_recoVisibleFourMomentumPx{};
	DoubleVector				m_recoVisibleFourMomentumPy{};
	DoubleVector				m_recoVisibleFourMomentumPz{};
	DoubleVector				m_recoVisibleFourMomentumE{};
	DoubleVector				m_residualVisibleFourMomentumPx{};
	DoubleVector				m_residualVisibleFourMomentumPy{};
	DoubleVector				m_residualVisibleFourMomentumPz{};
	DoubleVector				m_residualVisibleFourMomentumE{};
	DoubleVector				m_sigmaPxPx_Det{};
	DoubleVector				m_sigmaPxPy_Det{};
	DoubleVector				m_sigmaPyPy_Det{};
	DoubleVector				m_sigmaPxPz_Det{};
	DoubleVector				m_sigmaPyPz_Det{};
	DoubleVector				m_sigmaPzPz_Det{};
	DoubleVector				m_sigmaPxE_Det{};
	DoubleVector				m_sigmaPyE_Det{};
	DoubleVector				m_sigmaPzE_Det{};
	DoubleVector				m_sigmaEE_Det{};
	DoubleVector				m_normalizedResidualVisibleFourMomentumPx{};
	DoubleVector				m_normalizedResidualVisibleFourMomentumPy{};
	DoubleVector				m_normalizedResidualVisibleFourMomentumPz{};
	DoubleVector				m_normalizedResidualVisibleFourMomentumE{};
	DoubleVector				m_trueFlightDirectionUx{};
	DoubleVector				m_trueFlightDirectionUy{};
	DoubleVector				m_trueFlightDirectionUz{};
	DoubleVector				m_recoFlightDirectionUx{};
	DoubleVector				m_recoFlightDirectionUy{};
	DoubleVector				m_recoFlightDirectionUz{};
	DoubleVector				m_residualFlightDirectionUx{};
	DoubleVector				m_residualFlightDirectionUy{};
	DoubleVector				m_residualFlightDirectionUz{};
	DoubleVector				m_trueNeutrinoFourMomentumPx{};
	DoubleVector				m_trueNeutrinoFourMomentumPy{};
	DoubleVector				m_trueNeutrinoFourMomentumPz{};
	DoubleVector				m_trueNeutrinoFourMomentumE{};
	DoubleVector				m_recoNeutrinoFourMomentumClosePx{};
	DoubleVector				m_recoNeutrinoFourMomentumClosePy{};
	DoubleVector				m_recoNeutrinoFourMomentumClosePz{};
	DoubleVector				m_recoNeutrinoFourMomentumCloseE{};
	DoubleVector				m_sigmaNeutrinoPxPx{};
	DoubleVector				m_sigmaNeutrinoPxPy{};
	DoubleVector				m_sigmaNeutrinoPyPy{};
	DoubleVector				m_sigmaNeutrinoPxPz{};
	DoubleVector				m_sigmaNeutrinoPyPz{};
	DoubleVector				m_sigmaNeutrinoPzPz{};
	DoubleVector				m_sigmaNeutrinoPxE{};
	DoubleVector				m_sigmaNeutrinoPyE{};
	DoubleVector				m_sigmaNeutrinoPzE{};
	DoubleVector				m_sigmaNeutrinoEE{};
	DoubleVector				m_residualNeutrinoFourMomentumPx{};
	DoubleVector				m_residualNeutrinoFourMomentumPy{};
	DoubleVector				m_residualNeutrinoFourMomentumPz{};
	DoubleVector				m_residualNeutrinoFourMomentumE{};
	DoubleVector				m_normalizedResidualNeutrinoFourMomentumPx{};
	DoubleVector				m_normalizedResidualNeutrinoFourMomentumPy{};
	DoubleVector				m_normalizedResidualNeutrinoFourMomentumPz{};
	DoubleVector				m_normalizedResidualNeutrinoFourMomentumE{};
	DoubleVector				m_recoNeutrinoDirectionError{};
	DoubleVector				m_PCAatLeptonX{};
	DoubleVector				m_PCAatLeptonY{};
	DoubleVector				m_PCAatLeptonZ{};
	DoubleVector				m_PCAatOtherParticleX{};
	DoubleVector				m_PCAatOtherParticleY{};
	DoubleVector				m_PCAatOtherParticleZ{};
	DoubleVector				m_JetAxisX{};
	DoubleVector				m_JetAxisY{};
	DoubleVector				m_JetAxisZ{};
	DoubleVector				m_PrimaryVertexX{};
	DoubleVector				m_PrimaryVertexY{};
	DoubleVector				m_PrimaryVertexZ{};
	DoubleVector				m_AngleTrueFlightDirectionJet{};
	DoubleVector				m_AngleRecoFlightDirectionJet{};
	DoubleVector				m_AngleLeptonJet{};
	DoubleVector				m_AngleDSVertexJet{};
	DoubleVector				m_AngleRecoNeutrinoJet{};
	DoubleVector				m_AngleTrueNeutrinoJet{};
	DoubleVector				m_LeptonDistanceFromPV{};
	DoubleVector				m_DSVDistanceFromPV{};
	DoubleVector				m_Lepton3DImpactParameter{};
	DoubleVector				m_OtherParticle3DImpactParameter{};
	IntVector				m_SLDStatus{};
	IntVector				m_PVAStatus{};
	IntVector				m_trueSolutionSign{};
	IntVector				m_recoSolutionSign{};
	IntVector				m_nChargedParticlesInPrimVertex{};
	IntVector				m_nChargedParticlesNotInPrimVertex{};
	FloatVector				m_weightPFOtoMCP_Lepton{};
	FloatVector				m_weightMCPtoPFO_Lepton{};
	TH1I					*h_SLDStatus{};
	TH1F					*h_BHadronType{};
	TH1F					*h_CHadronType{};
	TH1F					*h_NuPxResidual{};
	TH1F					*h_NuPyResidual{};
	TH1F					*h_NuPzResidual{};
	TH1F					*h_NuEResidual{};
	TH1F					*h_NuPxNormalizedResidual{};
	TH1F					*h_NuPyNormalizedResidual{};
	TH1F					*h_NuPzNormalizedResidual{};
	TH1F					*h_NuENormalizedResidual{};
	TH2F					*h_recoNuPx_mcNuPx{};
	TH2F					*h_recoNuPy_mcNuPy{};
	TH2F					*h_recoNuPz_mcNuPz{};
	TH2F					*h_recoNuE_mcNuE{};
	TH2F					*h_parentPx_daughtersPx{};
	TH2F					*h_parentPy_daughtersPy{};
	TH2F					*h_parentPz_daughtersPz{};
	TH2F					*h_parentE_daughtersE{};
	TH1I					*h_recoPFOLinkedToElectron_Type{};
	TH1I					*h_recoPFOLinkedToMuon_Type{};
	TH1I					*h_SLDecayFlavour{};
	TH1I					*h_SLDecayModeB{};
	TH1I					*h_SLDecayModeC{};
	TH1I					*h_SLDecayOrder{};
	TH2I					*h_foundVertex{};
	TH1F					*h_secondaryVertex{};
	TH1I					*h_parentHadronCharge{};
	TH1I					*h_MCPTracks{};
	TH1I					*h_MCPTracks_Eweighted{};
	TH1I					*h_MCPTracks_Ptweighted{};
	TH1F					*h_FlightDirectionError{};
	TH1F					*h_distRecoLeptonToDownStreamVertex{};
	TFile					*m_pTFile{};
	TTree					*m_pTTree1 = new TTree( "SLDCorrection", "SLDCorrection" );


};
#endif
