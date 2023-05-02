#ifndef CheatedMCOverlayRemoval_h
#define CheatedMCOverlayRemoval_h 1

#include "marlin/Processor.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "lcio.h"
#include "EVENT/LCStrVec.h"
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/PIDHandler.h"
#include "UTIL/LCRelationNavigator.h"
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <vector>
class TFile;
class TH1F;
class TH1I;
class TH2I;
class TTree;

using namespace lcio ;
using namespace marlin ;
class CheatedMCOverlayRemoval : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new CheatedMCOverlayRemoval;
		}
		CheatedMCOverlayRemoval();
		virtual ~CheatedMCOverlayRemoval() = default;
		CheatedMCOverlayRemoval(const CheatedMCOverlayRemoval&) = delete;
		CheatedMCOverlayRemoval& operator=(const CheatedMCOverlayRemoval&) = delete;
		virtual void init();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		virtual void check();
		EVENT::MCParticle* getLinkedMCP( EVENT::ReconstructedParticle *recoParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedMCP , bool getNeutralMCP , float &weightPFOtoMCP , float &weightMCPtoPFO );
		EVENT::ReconstructedParticle* getLinkedPFO( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedPFO , bool getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO );
		virtual void end();
		void Clear();

	private:

		typedef std::vector<int>		IntVector;
		typedef std::vector<double>		DoubleVector;
		typedef std::vector<float>		FloatVector;

		int                                     m_nRun;
		int                                     m_nEvt;
		int                                     m_nAllPFOs;
		int                                     m_nMCPs;

		std::string                             _MCParticleCollectionName{};
		std::string                             _recoParticleCollectionName{}; 
		std::string                             _recoMCTruthLink{};
		std::string                             _mcTruthRecoLink{};
		std::string				_OutputPfoCollection{};
};

#endif
