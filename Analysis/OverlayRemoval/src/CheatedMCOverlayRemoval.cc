#include "CheatedMCOverlayRemoval.h"
#include <iostream>
#include <vector>
#include <string>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

CheatedMCOverlayRemoval aCheatedMCOverlayRemoval ;

CheatedMCOverlayRemoval::CheatedMCOverlayRemoval() :

  Processor("CheatedMCOverlayRemoval"),
  m_nRun(0),
  m_nEvt(0),
  m_nAllPFOs(0),
  m_nMCPs(0)
{

	_description = "CheatedMCOverlayRemoval identifies MC particles that are overlay and removes the corresponding PFOs from the collection" ;
	// Inputs: MC-particles, Reco-particles, the link between the two

	registerInputCollection( LCIO::MCPARTICLE,
				 "MCParticleCollection" ,
				 "Name of the MCParticle collection"  ,
				 _MCParticleCollectionName ,
				 std::string("MCParticlesSkimmed")
				 );

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				 "RecoParticleCollection" ,
				 "Name of the ReconstructedParticles input collection"  ,
				 _recoParticleCollectionName ,
				 std::string("PandoraPFOs")
				 );

	registerInputCollection( LCIO::LCRELATION,
				 "RecoMCTruthLink",
				 "Name of the RecoMCTruthLink input collection"  ,
				 _recoMCTruthLink,
				 std::string("RecoMCTruthLink")
				 );

	registerInputCollection( LCIO::LCRELATION,
                                 "MCTruthRecoLink",
                                 "Name of the MCTruthRecoLink input collection"  ,
                                 _mcTruthRecoLink,
                                 std::string("MCTruthRecoLink")
                                 );

	registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
				 "OutputPfoCollection",
				 "Name of output PFO collection",
				 _OutputPfoCollection,
				 std::string("PFOsWithoutOverlay")
				 );	
}

void CheatedMCOverlayRemoval::init()
{

	streamlog_out(DEBUG) << "   init called  " << std::endl ;
	this->Clear();

}

void CheatedMCOverlayRemoval::Clear()
{
  m_nAllPFOs = 0;
  m_nMCPs = 0;
}

void CheatedMCOverlayRemoval::processRunHeader()
{
  streamlog_out(DEBUG0) << "   processRunHeader called" << std::endl ;
  m_nRun++ ;  
  streamlog_out(DEBUG0) << "   processRunHeader finished successfully" << std::endl ;
}

void CheatedMCOverlayRemoval::processEvent( LCEvent *pLCEvent )
{
  this->Clear();
  m_nRun = pLCEvent->getRunNumber();
  m_nEvt = pLCEvent->getEventNumber();
  streamlog_out(DEBUG) << "Processing event " << pLCEvent->getEventNumber() << std::endl;

  const EVENT::LCCollection *PFOs{};
  const EVENT::LCCollection *MCPs{};
  IMPL::LCCollectionVec* OutputPfoCollection(NULL);
  OutputPfoCollection = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  OutputPfoCollection->setSubset( true ); 
  IMPL::LCCollectionVec* OutputOverlayCollection(NULL);
  OutputOverlayCollection = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  OutputOverlayCollection->setSubset( true ); 
  try
    {
      PFOs = pLCEvent->getCollection( _recoParticleCollectionName );
      MCPs = pLCEvent->getCollection( _MCParticleCollectionName );

      LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( _recoMCTruthLink ) );
      LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( _mcTruthRecoLink ) );

      m_nAllPFOs = PFOs->getNumberOfElements();
      m_nMCPs = MCPs->getNumberOfElements();
      int nKeptPFOs = 0;
      int nRemovedPFOs = 0;
      
      //Create new collection without PFOs from MCPs that are overlay
      for (int i=0; i<m_nMCPs; i++) {
	MCParticle* mcp = (MCParticle*) MCPs->getElementAt(i);
	if (mcp->isOverlay()) { 
	  float weightPFOtoMCP = 0.0;
	  float weightMCPtoPFO = 0.0;
	  ReconstructedParticle* linkedPFO = getLinkedPFO( mcp , RecoMCParticleNav , MCParticleRecoNav , false , false , weightPFOtoMCP , weightMCPtoPFO );
	  if ( linkedPFO == NULL ) {
	    continue;
	  }
	  streamlog_out(DEBUG5) << "i = " << i << ": energy = " << linkedPFO->getEnergy () << std::endl;
	  OutputOverlayCollection->addElement(linkedPFO);
	  nRemovedPFOs++;
	}
      }

      for (int i=0; i<m_nAllPFOs; i++) {
	ReconstructedParticle* pfo = (ReconstructedParticle*) PFOs->getElementAt(i);
	float weightPFOtoMCP = 0.0;
	float weightMCPtoPFO = 0.0;
	MCParticle* linkedMCP = getLinkedMCP( pfo, RecoMCParticleNav , MCParticleRecoNav , false , false , weightPFOtoMCP , weightMCPtoPFO );
	if ( linkedMCP != NULL ) {
	  if (linkedMCP->isOverlay()) streamlog_out(DEBUG5) << "i = " << i << ": energy = " << pfo->getEnergy () << ", gen status = " << linkedMCP->getGeneratorStatus() << std::endl;
	  if (linkedMCP->isOverlay()) continue;
	}
	OutputPfoCollection->addElement(pfo);
	nKeptPFOs++;
      }

      streamlog_out(DEBUG5) << "In event " << m_nEvt << ": Kept PFOs = " << nKeptPFOs << " VS removed number of PFOs = " << nRemovedPFOs << " VS total number of PFOs = " << m_nAllPFOs << std::endl;
      pLCEvent->addCollection(OutputPfoCollection, _OutputPfoCollection.c_str() );
      m_nEvt++ ;
    }
  catch(...)
    {
      streamlog_out(WARNING) << "Check : Input collections not found in event " << m_nEvt << std::endl;
    }
  streamlog_out(DEBUG) << "nevt = " << m_nEvt << std::endl;
}


EVENT::MCParticle* CheatedMCOverlayRemoval::getLinkedMCP( EVENT::ReconstructedParticle *recoParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedMCP , bool getNeutralMCP , float &weightPFOtoMCP , float &weightMCPtoPFO )
{
  streamlog_out(DEBUG1) << "" << std::endl;
  streamlog_out(DEBUG1) << "Look for MCP linked to Reconstructed Particle:" << std::endl;

  MCParticle* linkedMCP{};
  bool foundlinkedMCP = false;
  double maxweightPFOtoMCP = 0.;
  double maxweightMCPtoPFO = 0.;
  int iPFOtoMCPmax = -1;
  int iMCPtoPFOmax = -1;
  const EVENT::LCObjectVec& MCPvec = RecoMCParticleNav.getRelatedToObjects( recoParticle );
  const EVENT::FloatVec&  MCPweightvec = RecoMCParticleNav.getRelatedToWeights( recoParticle );
  for ( unsigned int i_mcp = 0; i_mcp < MCPvec.size(); i_mcp++ )
    {
      double mcp_weight = 0.0;
      double trackWeight = ( int( MCPweightvec.at( i_mcp ) ) % 10000 ) / 1000.0;
      double clusterWeight = ( int( MCPweightvec.at( i_mcp ) ) / 10000 ) / 1000.0;
      if ( getChargedMCP && !getNeutralMCP )
	{
	  mcp_weight = trackWeight;
	}
      else if ( getNeutralMCP && !getChargedMCP )
	{
	  mcp_weight = clusterWeight;
	}
      else
	{
	  mcp_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
	}
      MCParticle *testMCP = (MCParticle *) MCPvec.at( i_mcp );
      if ( mcp_weight > maxweightPFOtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
	{
	  maxweightPFOtoMCP = mcp_weight;
	  iPFOtoMCPmax = i_mcp;
	  streamlog_out(DEBUG0) << "MCParticle at index: " << testMCP->id() << " has PDG: " << testMCP->getPDG() << " and PFO to MCParticle link weight is " << mcp_weight << std::endl;
	}
    }
  if ( iPFOtoMCPmax != -1 )
    {
      MCParticle *testMCP = (MCParticle *) MCPvec.at( iPFOtoMCPmax );
      const EVENT::LCObjectVec& PFOvec = MCParticleRecoNav.getRelatedToObjects( testMCP );
      const EVENT::FloatVec&  PFOweightvec = MCParticleRecoNav.getRelatedToWeights( testMCP );
      streamlog_out(DEBUG0) << "Test Visible MCParticle is linked to " << PFOvec.size() << " PFO(s)" << std::endl;
      for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
	{
	  double pfo_weight = 0.0;
	  double trackWeight = ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0;
	  double clusterWeight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
	  if ( getChargedMCP && !getNeutralMCP )
	    {
	      pfo_weight = trackWeight;
	    }
	  else if ( getNeutralMCP && !getChargedMCP )
	    {
	      pfo_weight = clusterWeight;
	    }
	  else
	    {
	      pfo_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
	    }
	  streamlog_out(DEBUG0) << "Test Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << trackWeight << " , Cluster: " << clusterWeight << ")" << std::endl;
	  ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
	  if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
	    {
	      maxweightMCPtoPFO = pfo_weight;
	      iMCPtoPFOmax = i_pfo;
	      streamlog_out(DEBUG0) << "PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
	    }
	}
      if ( iMCPtoPFOmax != -1 )
	{
	  if ( PFOvec.at( iMCPtoPFOmax ) == recoParticle )
	    {
	      linkedMCP = testMCP;
	      foundlinkedMCP = true;
	    }
	}
    }
  if( foundlinkedMCP )
    {
      streamlog_out(DEBUG1) << "Linked MCParticle to PFO found successfully " << std::endl;
      weightPFOtoMCP = maxweightPFOtoMCP;
      weightMCPtoPFO = maxweightMCPtoPFO;
      return linkedMCP;
    }
  else
    {
      streamlog_out(DEBUG1) << "Couldn't Find a MCParticle linked to PFO" << std::endl;
      return NULL;
    }

}

EVENT::ReconstructedParticle* CheatedMCOverlayRemoval::getLinkedPFO( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedPFO , bool getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO )
{
  streamlog_out(DEBUG1) << "" << std::endl;
  streamlog_out(DEBUG1) << "Look for PFO linked to visible MCParticle:" << std::endl;

  ReconstructedParticle* linkedPFO{};
  bool foundlinkedPFO = false;
  const EVENT::LCObjectVec& PFOvec = MCParticleRecoNav.getRelatedToObjects( mcParticle );
  const EVENT::FloatVec&  PFOweightvec = MCParticleRecoNav.getRelatedToWeights( mcParticle );
  streamlog_out(DEBUG0) << "Visible MCParticle is linked to " << PFOvec.size() << " PFO(s)" << std::endl;
  weightPFOtoMCP = 0.0;
  weightMCPtoPFO = 0.0;
  double maxweightPFOtoMCP = 0.;
  double maxweightMCPtoPFO = 0.;
  int iPFOtoMCPmax = -1;
  int iMCPtoPFOmax = -1;
  for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
    {
      double pfo_weight = 0.0;
      double trackWeight = ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0;
      double clusterWeight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
      if ( getChargedPFO && !getNeutralPFO )
	{
	  pfo_weight = trackWeight;
	}
      else if ( getNeutralPFO && !getChargedPFO )
	{
	  pfo_weight = clusterWeight;
	}
      else
	{
	  pfo_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
	}
      streamlog_out(DEBUG0) << "Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << trackWeight << " , Cluster: " << clusterWeight << ")" << std::endl;
      ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
      if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
	{
	  maxweightMCPtoPFO = pfo_weight;
	  iMCPtoPFOmax = i_pfo;
	  streamlog_out(DEBUG0) << "PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
	}
    }
  if ( getChargedPFO && maxweightMCPtoPFO < 0.8 )
    {
      streamlog_out(DEBUG1) << "MCParticle has link weight lower than 0.8 ( " << maxweightMCPtoPFO << " ), looking for linked PFO in clusters" << std::endl;
      for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
	{
	  double pfo_weight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
	  streamlog_out(DEBUG0) << "Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0 << " , Cluster: " << ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0 << ")" << std::endl;
	  ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
	  if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
	    {
	      maxweightMCPtoPFO = pfo_weight;
	      iMCPtoPFOmax = i_pfo;
	      streamlog_out(DEBUG0) << "PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
	    }
	}
    }
  if ( iMCPtoPFOmax != -1 )
    {
      ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( iMCPtoPFOmax );
      const EVENT::LCObjectVec& MCPvec = RecoMCParticleNav.getRelatedToObjects( testPFO );
      const EVENT::FloatVec&  MCPweightvec = RecoMCParticleNav.getRelatedToWeights( testPFO );
      for ( unsigned int i_mcp = 0; i_mcp < MCPvec.size(); i_mcp++ )
	{
	  double mcp_weight = 0.0;
	  double trackWeight = ( int( MCPweightvec.at( i_mcp ) ) % 10000 ) / 1000.0;
	  double clusterWeight = ( int( MCPweightvec.at( i_mcp ) ) / 10000 ) / 1000.0;
	  if ( getChargedPFO && !getNeutralPFO )
	    {
	      mcp_weight = trackWeight;
	    }
	  else if ( getNeutralPFO && !getChargedPFO )
	    {
	      mcp_weight = clusterWeight;
	    }
	  else
	    {
	      mcp_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
	    }
	  MCParticle *testMCP = (MCParticle *) MCPvec.at( i_mcp );
	  if ( mcp_weight > maxweightPFOtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
	    {
	      maxweightPFOtoMCP = mcp_weight;
	      iPFOtoMCPmax = i_mcp;
	      streamlog_out(DEBUG0) << "MCParticle at index: " << testMCP->id() << " has PDG: " << testMCP->getPDG() << " and PFO to MCParticle link weight is " << mcp_weight << std::endl;
	    }
	}
      if ( iPFOtoMCPmax != -1 )
	{
	  if ( MCPvec.at( iPFOtoMCPmax ) == mcParticle )
	    {
	      linkedPFO = testPFO;
	      foundlinkedPFO = true;
	    }
	}
    }

  if( foundlinkedPFO )
    {
      streamlog_out(DEBUG1) << "Linked PFO to MCParticle found successfully " << std::endl;
      weightPFOtoMCP = maxweightPFOtoMCP;
      weightMCPtoPFO = maxweightMCPtoPFO;
      return linkedPFO;
    }
  else
    {
      streamlog_out(DEBUG1) << "Couldn't Find a PFO linked to MCParticle" << std::endl;
      return NULL;
    }
}

void CheatedMCOverlayRemoval::check()
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void CheatedMCOverlayRemoval::end()
{

}



