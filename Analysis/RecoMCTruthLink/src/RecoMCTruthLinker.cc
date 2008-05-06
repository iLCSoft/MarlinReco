#include "RecoMCTruthLinker.h"
#include <iostream>

// LCIO 
#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/Cluster.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/LCRelation.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/LCRelationNavigator.h>

#include <math.h>
#include <map>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/AIDA.h>
#endif


using namespace lcio ;
using namespace marlin ;



RecoMCTruthLinker aRecoMCTruthLinker;


struct MCPKeep :  public LCIntExtension<MCPKeep> {} ;

typedef std::map< MCParticle* , int > MCPMap ;

typedef std::map< MCParticle* , double > MCPMapDouble ;


RecoMCTruthLinker::RecoMCTruthLinker() : Processor("RecoMCTruthLinker") {
  
  // modify processor description
  _description = "links RecontructedParticles to the MCParticle based on number of hits used" ;
  

  IntVec pdgVecDef ;
    
  pdgVecDef.push_back(  22 ) ;  // gamma
  pdgVecDef.push_back( 111 ) ;  // pi0
  pdgVecDef.push_back( 310 ) ;  // K0s
    

  registerProcessorParameter(  "KeepDaughtersPDG" , 
			       "PDG codes of particles of which the daughters will be "
			       "kept in the skimmmed MCParticle collection"  ,
			       _pdgVec ,
			       pdgVecDef ) ;
    
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection" , 
			   "Name of the MCParticle input collection"  ,
			   _mcParticleCollectionName ,
			   std::string("MCParticle") ) ;
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "RecoParticleCollection" , 
			   "Name of the ReconstructedParticles input collection"  ,
			   _recoParticleCollectionName ,
			   std::string("ReconstructedParticles") ) ;
  
  registerInputCollection( LCIO::LCRELATION,
			   "SimTrackerHitRelationName" , 
			   "Name of the  SimTrackerHit - TrackerHit relation"  ,
			   _trackHitRelationName ,
			   std::string("SimTrackerHitRelation") ) ;
  
  registerInputCollection( LCIO::LCRELATION,
			   "SimCalorimeterHitRelationName" , 
			   "Name of the  SimCalorimeterHit - CalorimeterHit relation"  ,
			   _caloHitRelationName ,
			   std::string("SimCalorimeterHitRelation") ) ;


  registerOutputCollection( LCIO::LCRELATION,
			    "RecoMCTruthLinkName" , 
			    "Name of the RecoMCTruthLink output collection"  ,
			    _recoMCTruthLinkName ,
			    std::string("RecoMCTruthLink") ) ;


  registerOutputCollection( LCIO::MCPARTICLE,
			    "MCParticlesSkimmedName" , 
			    "Name of the skimmed MCParticle  output collection"  ,
			    _mcParticlesSkimmedName ,
			    std::string("MCParticlesSkimmed") ) ;
    
    
}


void RecoMCTruthLinker::init() { 
  
  // usually a good idea to
  printParameters<MESSAGE>() ;
  
  for( IntVec::iterator it = _pdgVec.begin() ; it != _pdgVec.end() ; ++it ){

    if( *it < 0 ) {
	
      streamlog_out( WARNING ) << " init: negative PDG given - only abs value is used : " <<  *it  << std::endl ;
    }

    _pdgSet.insert( abs( *it )  )  ;
  }

  _nRun = 0 ;
  _nEvt = 0 ;
}


void RecoMCTruthLinker::processRunHeader( LCRunHeader* run) { 
  
  _nRun++ ;
} 


void RecoMCTruthLinker::processEvent( LCEvent * evt ) { 

  _nEvt ++ ;

  
  streamlog_out( DEBUG ) << " processEvent "  <<  evt->getEventNumber() << "  - " << evt->getRunNumber() 
			 << std::endl ; 

  LCCollection* mcpCol = evt->getCollection( _mcParticleCollectionName ) ;
  LCCollection* recoCol = evt->getCollection( _recoParticleCollectionName ) ;

  //    LCCollection* tHitRelCol = evt->getCollection( _trackHitRelationName ) ;
  LCCollection* cHitRelCol = evt->getCollection( _caloHitRelationName ) ;
  
  //    LCRelationNavigator tHitRelNav( tHitRelCol ) ;
  LCRelationNavigator cHitRelNav( cHitRelCol ) ;
  

  LCRelationNavigator truthRelNav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE  ) ;


  // loop over reconstructed particles
  int nReco = recoCol->getNumberOfElements() ;

  for(int i=0;i<nReco;++i){
      
    ReconstructedParticle* rec = dynamic_cast<ReconstructedParticle*> ( recoCol->getElementAt(i) ) ;

    float charge =  rec->getCharge()  ;

    streamlog_out( DEBUG0 ) << " rec particle:  e: " << rec->getEnergy()  
			    << " charge: " << rec->getCharge() 
			    << " px: " << rec->getMomentum()[0]
			    << " py: " << rec->getMomentum()[1]
			    << " pz: " << rec->getMomentum()[2]
			    << " mass: " << rec->getMass() 
			    <<  " *pdg " << rec->getParticleIDUsed()
			    << std::endl ;

    const ReconstructedParticleVec& recDaughters = rec->getParticles() ;
      
    if( recDaughters.size() > 0 ) {  // ignore compound particles
	
      streamlog_out( DEBUG0 ) << " compound particle -> ignored " << std::endl ;
	
      continue ;
    }
      
    if( fabs( charge ) > .01 ) {   // charged particle -> use track
	
      MCPMap mcpMap ;

      const TrackVec& trkVec = rec->getTracks() ;

      unsigned nTrk = trkVec.size() ;
	
      streamlog_out( DEBUG0 ) << " loop over " <<  nTrk << " tracks " << std::endl ;
	
      int nHit = 0 ;
	
      for( TrackVec::const_iterator trkIt = trkVec.begin() ;
	   trkIt != trkVec.end() ; ++trkIt ) {
	  
	const Track* trk = *trkIt ;
	  
	const TrackerHitVec& trkHits = trk->getTrackerHits() ;
	  
	for( TrackerHitVec::const_iterator hitIt = trkHits.begin() ;
	     hitIt != trkHits.end() ; ++hitIt ) { 
	    
	  const TrackerHit* hit = * hitIt ;
	    
	  const LCObjectVec& simHits  = hit->getRawHits() ;
	    
	  for( LCObjectVec::const_iterator objIt = simHits.begin() ;
	       objIt != simHits.end() ; ++objIt ){
	      
	    SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>( *objIt ) ;

	    MCParticle* mcp = simHit->getMCParticle() ; 

	    mcpMap[ mcp ]++ ;   // count the hit caused by this particle

	    ++nHit ; // total hit count

	  }
	}
      }

      if( nTrk == 0 ){
	  
	streamlog_out( WARNING ) << " charged particle without tracks   " <<  std::endl ;
	  
	if( nHit == 0 ){

	  streamlog_out( WARNING ) << " no tracker hits found " <<  std::endl ;
	}

	break ;  // won't find a particle 
      }


      // find the mc particle with the largest #hits 
	
      int maxHit = 0 ;
      MCParticle* theMCP = 0 ;
      for( MCPMap::iterator it = mcpMap.begin() ;
	   it != mcpMap.end() ; ++it ){
	  
	if( it->second > maxHit ){
	  maxHit = it->second ;
	  theMCP = it->first ;	    
	}
      }
	
      float  weight = float( maxHit) / nHit ;
	
      streamlog_out( DEBUG0 ) << " reco particle has " << maxHit << " hits of " << nHit 
			      << " [ " << weight << " ] " 
			      << " of MCParticle with pdg : " << theMCP->getPDG() 
			      << std::endl ;
	
	
      truthRelNav.addRelation(   rec , theMCP , weight ) ;
	



    } else { // neutral particles -> use cluster
	

      MCPMapDouble mcpEnergy ;
	
      const ClusterVec& cluVec = rec->getClusters() ;

      unsigned nClu = cluVec.size() ;
	
      streamlog_out( DEBUG0 ) << " loop over " <<  nClu << " clusters " << std::endl ;
	
      double eTot = 0 ;
	
      for( ClusterVec::const_iterator cluIt = cluVec.begin() ;
	   cluIt != cluVec.end() ; ++cluIt ) {
	  
	const Cluster* clu = *cluIt ;
	  
	const CalorimeterHitVec& cluHits = clu->getCalorimeterHits() ;
	  
	for( CalorimeterHitVec::const_iterator hitIt = cluHits.begin() ;
	     hitIt != cluHits.end() ; ++hitIt ) { 
	    
	  const CalorimeterHit* hit = * hitIt ;
	    
	  // 	     SimCalorimeterHit* simHit = dynamic_cast<SimCalorimeterHit*>( hit->getRawHit() ) ;
	  // 	     streamlog_out( DEBUG0 ) << "   raw hit ptr : " << simHit << std::endl ;


	  const LCObjectVec& simHits = cHitRelNav.getRelatedToObjects( (LCObject*) hit )  ;
	    
	  for( LCObjectVec::const_iterator objIt = simHits.begin() ;
	       objIt != simHits.end() ; ++objIt ){
	      
	    SimCalorimeterHit* simHit = dynamic_cast<SimCalorimeterHit*>( *objIt ) ;

	    for(int j=0;j<simHit->getNMCContributions() ;++j){
		
	      MCParticle* mcp = simHit->getParticleCont( j ) ;
		
	      double e  = simHit->getEnergyCont( j ) ;
		
	      mcpEnergy[ mcp ] +=  e ;
		
	      eTot += e ; 

	    }
	  }
	}
      }

      if( nClu == 0 ){
	  
	streamlog_out( WARNING ) << " neutral particle without clusters   " <<  std::endl ;
	  
	break ;  // won't find a particle 
      }  

      if( eTot == 0.0 ){ // fixme - this might happen if clusters are from Lcal/Muon only
	  
	streamlog_out( DEBUG ) << " no calorimeter energy found for " 
			       << " reco particle: e:"  << rec->getEnergy()  
			       << " charge: " << rec->getCharge() 
			       << " px: " << rec->getMomentum()[0]
			       << " py: " << rec->getMomentum()[1]
			       << " pz: " << rec->getMomentum()[2]
			       << " mass: " << rec->getMass() 
			       <<  " *pdg " << rec->getParticleIDUsed()
			       << std::endl ;
	
	  
	break ;  
      }


      // find the mc particle with the largest energy contribution 
	
      double eMax = 0 ;
      MCParticle* theMCP = 0 ;

      for( MCPMapDouble::iterator it = mcpEnergy.begin() ;
	   it != mcpEnergy.end() ; ++it ){
	  
	if( it->second > eMax  ){

	  eMax = it->second ;
	  theMCP = it->first ;	    
	}
      }
	
      float  weight = eMax / eTot ;
	

      if( theMCP == 0 ) {

	streamlog_out( ERROR ) << " reco particle has " << eMax << " GeV of " << eTot 
			       << " GeV [ " << weight << " ] true energy in " << nClu << " clusters "  
			       << std::endl ;
	  
	break ;

      }

      streamlog_out( DEBUG0 ) << " reco particle has " << eMax << " GeV of " << eTot 
			      << " GeV [ " << weight << " ] " 
			      << " of MCParticle with pdg : " << theMCP->getPDG() 
			      << std::endl ;
	
	
      truthRelNav.addRelation(   rec , theMCP , weight ) ;
	
    }
      
  }

  evt->addCollection(  truthRelNav.createLCCollection() , _recoMCTruthLinkName  ) ;

  //-------------- create skimmed MCParticle collection ------------------------

  LCCollectionVec* skimVec = new LCCollectionVec( LCIO::MCPARTICLE )  ;
    
  skimVec->setSubset( true) ;  // flag as subset 
    
    
  int nMCP  = mcpCol->getNumberOfElements() ;
    
  for(int i=0; i< nMCP ; i++){
      
    MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
      
      
    if( mcp->ext<MCPKeep>() == true ){
	
      continue ;    // particle allready in skim 
    }
      
    //       if ( mcp->isCreatedInSimulation()  && mcp->getGeneratorStatus()  != 0  ) 

    // 	streamlog_out( WARNING ) << " mcp->isCreatedInSimulation()  && mcp->getGeneratorStatus()  != 0 TRUE" 
    // 				 <<  " :" << mcp->getEnergy()  
    // 				 << " charge: " << mcp->getCharge() 
    // 				 << " genstat: " << mcp->getGeneratorStatus() 

    // 				 << " simstat: " << std::hex << mcp->getSimulatorStatus()  << std::dec 

    // 				 << " : "  << std::endl 
    // 				 << " isCreatedInSimulation :" << mcp->isCreatedInSimulation()	<< std::endl
    // 				 << " isBackscatter :" << mcp->isBackscatter()	<< std::endl
    // 				 << " vertexIsNotEndpointOfParent :" << mcp->vertexIsNotEndpointOfParent() << std::endl
    // 				 << " isDecayedInTracker :" << mcp->isDecayedInTracker() << std::endl
    // 				 << " isDecayedInCalorimeter :" << mcp->isDecayedInCalorimeter() << std::endl
    // 				 << " hasLeftDetector :" << mcp->hasLeftDetector() << std::endl
    // 				 << " isStopped :" << mcp->isStopped() << "  : " 
    // 				 << " pdg " << mcp->getPDG()
    // 				 << std::endl ;
      


    if ( ! mcp->isCreatedInSimulation() || mcp->getGeneratorStatus()  != 0 )  { 

      //FIXME: this is a workaround for a Mokka bug: the isCreatedInSimulation 
      // is set also for generated particles....

      // keep all generated particles (complete event)
	
      mcp->ext<MCPKeep>() = true  ;
	
      continue ;

    } else { // of those created in the simulation we keep those that actually are reconstructed
      // including all parents

      const LCObjectVec& recoObjects = truthRelNav.getRelatedFromObjects( (LCObject*) mcp )  ;

      if( recoObjects.size() > 0 ){

	streamlog_out( DEBUG0 ) << " keep MCParticle - e :" << mcp->getEnergy()  
				<< " charge: " << mcp->getCharge() 
				<< " px: " << mcp->getMomentum()[0]
				<< " py: " << mcp->getMomentum()[1]
				<< " pz: " << mcp->getMomentum()[2]
				<< " mass: " << mcp->getMass() 
				<< " pdg " << mcp->getPDG()
				<< std::endl ;

	keepMCParticle( mcp ) ;
      }	

    } // else

  }   // end mcp loop 
 
  // --- loop again and add daughters of particles that are in the skimmed list and have a pdg in
  //     the parameter vector 'KeepDaughtersPDG 

  for(int i=0; i< nMCP ; i++){
    
    MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    
    // keep the daughters of all decays in flight of particles in the pdg list (default: gamma, pi0, K0s) 
    if( mcp->ext<MCPKeep>() == true &&  mcp->isDecayedInTracker()  ){ //&& !mcp->isStopped()  ){  
      
      unsigned thePDG = abs( mcp->getPDG() ) ;
      
      if( _pdgSet.find( thePDG ) != _pdgSet.end()  ) {
	
	const MCParticleVec& daughters = mcp->getDaughters() ;
	
	streamlog_out( DEBUG0 ) << " keeping daughters of particle with pdg : " << mcp->getPDG() << " : " 
				<< " [" << mcp->getGeneratorStatus() << "] :";
	// 				  << " e :" << mcp->getEnergy() 
	// 				  << " isCreatedInSimulation :" << mcp->isCreatedInSimulation()	<< std::endl
	// 				  << " isBackscatter :" << mcp->isBackscatter()	<< std::endl
	// 				  << " vertexIsNotEndpointOfParent :" << mcp->vertexIsNotEndpointOfParent()	<< std::endl
	// 				  << " isDecayedInTracker :" << mcp->isDecayedInTracker()	<< std::endl
	// 				  << " isDecayedInCalorimeter :" << mcp->isDecayedInCalorimeter()	<< std::endl
	// 				  << " hasLeftDetector :" << mcp->hasLeftDetector()	<< std::endl
	// 				  << " isStopped :" << mcp->isStopped()    << "  : " 

	streamlog_message( DEBUG0 , 
			   if( mcp->getParents().size() ) ,
			   " parent pdg : " << mcp->getParents()[0]->getPDG() << " : "  ;
			   ) ;
	
	streamlog_out( DEBUG0 ) << std::endl ;
	
	//	<< std::endl ;
	
	for( MCParticleVec::const_iterator dIt = daughters.begin() ;
	     dIt != daughters.end() ; ++dIt ){
	  
	  (*dIt)->ext<MCPKeep>() = true ;
	  
	  streamlog_out( DEBUG0 ) <<  (*dIt)->getPDG() << ", " ;
	}
	
	streamlog_out( DEBUG0 ) << std::endl ;
	
      }
    }
  }
  
  
  for(int i=0; i< nMCP ; i++){
    
    MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    
    if( mcp->ext<MCPKeep>() == true ) {  
      
      skimVec->addElement( mcp ) ;
    }
  }    
  
  evt->addCollection(  skimVec , _mcParticlesSkimmedName ) ;
}

void  RecoMCTruthLinker::keepMCParticle( MCParticle* mcp ){
  
  mcp->ext<MCPKeep>() = true  ;
	

  const MCParticleVec& parents = mcp->getParents() ;
    
  streamlog_out( DEBUG0 ) << " keepMCParticle keep particle with pdg : " << mcp->getPDG() 
			  << std::endl ;
    
  for( MCParticleVec::const_iterator pIt = parents.begin() ;
       pIt != parents.end() ; ++pIt ){
      
    if(  (*pIt )->ext<MCPKeep>() != true  ) { // if parent not yet in skim 
	
      // add it
      keepMCParticle(  *pIt ) ;
    }
      
  }
}
  
void RecoMCTruthLinker::check( LCEvent * evt ) { 

  // ---- create some checkplots

#ifdef MARLIN_USE_AIDA
  
  // - define some static histo pointers 
  // FIXME: these need to become class members eventually ...

  //  static AIDA::IHistogram1D* hTrack_z0 ;
  static AIDA::ICloud1D* hmcp_etot ;
  static AIDA::ICloud1D* hmcp_e ;
  static AIDA::ICloud1D* hmcp_n ;
  static AIDA::ICloud1D* hmcp_ntot ;

  static AIDA::ICloud1D* hmcpsk_etot ;
  static AIDA::ICloud1D* hmcpsk_e ;
  static AIDA::ICloud1D* hmcpsk_n ;
  static AIDA::ICloud1D* hmcpsk_ntot ;
    
  if( isFirstEvent() ) { 
      
    hmcp_e = AIDAProcessor::histogramFactory(this)->
      createCloud1D( "hmcp_e", " energy/GeV - all " , 100 ) ; 
    //       createHistogram1D( "hmcp_e", " energy/GeV ", 100, 0. , 10. ) ; 
      
    hmcp_etot = AIDAProcessor::histogramFactory(this)->
      createCloud1D( "hmcp_etot", " total energy/GeV " , 100 ) ; 
    //       createHistogram1D( "hmcp_etot_e", " energy/GeV ", 1000, 0. , 1000. ) ; 

    hmcp_n = AIDAProcessor::histogramFactory(this)->
      createCloud1D( "hmcp_n", " # generated stable particles " , 100 ) ; 

    hmcp_ntot = AIDAProcessor::histogramFactory(this)->
      createCloud1D( "hmcp_ntot", "  total # particles " , 100 ) ; 


    hmcpsk_e = AIDAProcessor::histogramFactory(this)->
      createCloud1D( "hmcpsk_e", " energy/GeV - all " , 100 ) ; 
    //       createHistogram1D( "hmcpsk_e", " energy/GeV ", 100, 0. , 10. ) ; 
      
    hmcpsk_etot = AIDAProcessor::histogramFactory(this)->
      createCloud1D( "hmcpsk_etot", " total energy/GeV " , 100 ) ; 
    //       createHistogram1D( "hmcpsk_etot_e", " energy/GeV ", 1000, 0. , 1000. ) ; 

    hmcpsk_n = AIDAProcessor::histogramFactory(this)->
      createCloud1D( "hmcpsk_n", " # generated stable particles " , 100 ) ; 

    hmcpsk_ntot = AIDAProcessor::histogramFactory(this)->
      createCloud1D( "hmcpsk_ntot", "  total # particles " , 100 ) ; 
  }


  LCCollection* mcpCol = evt->getCollection( _mcParticleCollectionName ) ;

  int nMCP  = mcpCol->getNumberOfElements() ;
    
  double etot = 0.0 ;
  int nStable = 0 ;

  for(int i=0; i< nMCP ; i++){
      
    MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;

    if( mcp->getGeneratorStatus() == 1 ) {
	
      hmcp_e->fill(   mcp->getEnergy()  ) ;

      etot +=  mcp->getEnergy()  ;
	
      ++nStable ;
    }

  }
  hmcp_n->fill( nStable ) ;

  hmcp_ntot->fill( nMCP ) ;
    
  hmcp_etot->fill( etot ) ;
    

  // create the same histos now with the skimmed collection

  LCCollection* mcpskCol = evt->getCollection( _mcParticlesSkimmedName ) ;

  int nMCPSK  = mcpskCol->getNumberOfElements() ;
    
  etot = 0.0 ;
  nStable = 0 ;

  for(int i=0; i< nMCPSK ; i++){
      
    MCParticle* mcpsk = dynamic_cast<MCParticle*> ( mcpskCol->getElementAt( i ) ) ;

    if( mcpsk->getGeneratorStatus() == 1 ) {
	
      hmcpsk_e->fill(   mcpsk->getEnergy()  ) ;

      etot +=  mcpsk->getEnergy()  ;
	
      ++nStable ;
    }

  }
  hmcpsk_n->fill( nStable ) ;

  hmcpsk_ntot->fill( nMCPSK ) ;
    
  hmcpsk_etot->fill( etot ) ;
    

#endif

}


void RecoMCTruthLinker::end(){ 
  
  streamlog_out(DEBUG) << " processed " << _nEvt << " events in " << _nRun << " runs "
		       << std::endl ;

}







