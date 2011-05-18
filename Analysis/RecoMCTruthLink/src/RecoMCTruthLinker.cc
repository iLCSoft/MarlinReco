#include "RecoMCTruthLinker.h"
#include <iostream>

#include <cstdlib>

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
#include <UTIL/LCRelationNavigator.h>

#include "gearimpl/Vector3D.h"

#include <math.h>
#include <map>
#include <algorithm>

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

typedef std::map< MCParticle* , MCParticle* >  Remap_as_you_go;

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
  
  registerInputCollection( LCIO::TRACK,
			   "TrackCollection" , 
			   "Name of the Tracks input collection"  ,
			   _trackCollectionName ,
			   std::string("LDCTracks") ) ;

  registerInputCollection( LCIO::CLUSTER,
			   "ClusterCollection" , 
			   "Name of the Clusters input collection"  ,
			   _clusterCollectionName ,
			   std::string("PandoraClusters") ) ;
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "RecoParticleCollection" ,
                           "Name of the ReconstructedParticles input collection"  ,
                           _recoParticleCollectionName ,
                           std::string("PandoraPFOs") ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "SimTrackerHitRelationName" , 
			   "Name of the  SimTrackerHit - TrackerHit relation"  ,
			   _trackHitRelationName ,
			   std::string("TPCTracksMCP") ) ;
  
  registerInputCollection( LCIO::LCRELATION,
			   "SimCalorimeterHitRelationName" , 
			   "Name of the  SimCalorimeterHit - CalorimeterHit relation"  ,
			   _caloHitRelationName ,
			   std::string("RelationCaloHit") ) ;


  registerProcessorParameter("OutputTrackRelation",
                             "true: Create TrackMCTruthLink collection, false : dont ",
                             _OutputTrackRelation,
                             bool(true));

  registerOutputCollection( LCIO::LCRELATION,
			    "TrackMCTruthLinkName" , 
			    "Name of the trackMCTruthLink output collection"  ,
			    _trackMCTruthLinkName ,
			    std::string("TrackMCTruthLink") ) ;
 
  registerProcessorParameter("OutputClusterRelation",
                             "true: Create ClusterMCTruthLink collection, false : dont ",
                             _OutputClusterRelation,
                             bool(true));

 registerOutputCollection( LCIO::LCRELATION,
			    "ClusterMCTruthLinkName" , 
			    "Name of the clusterMCTruthLink output collection"  ,
			    _clusterMCTruthLinkName ,
			    std::string("ClusterMCTruthLink") ) ;


  registerProcessorParameter("FullRecoRelation",
                             "true: All reco <-> true relations are given, with weight = 10000*calo weight+"
                             "track weight (weights in permill). false: Only highest contributor linked,"
                             "and only to tracks, not clusters if there are any tracks",   
                             _FullRecoRelation,
                             bool(false)
			     );

  registerOutputCollection( LCIO::LCRELATION,
                            "RecoMCTruthLinkName" ,
                            "Name of the RecoMCTruthLink output collection"  ,
                            _recoMCTruthLinkName ,
                            std::string("RecoMCTruthLink") ) ;

  registerProcessorParameter("OutputCalohitRelation",
                             "true: Create CalohitMCTruthLink collection, false : dont ",
                             _OutputCalohitRelation,
                             bool(false));

  registerOutputCollection( LCIO::LCRELATION,
                            "CalohitMCTruthLinkName" ,
                            "Name of the updated calo-hit MCTruthLink output collection"  ,
                            _calohitMCTruthLinkName ,
                            std::string("CalohitMCTruthLink") ) ;

  registerOutputCollection( LCIO::MCPARTICLE,
			    "MCParticlesSkimmedName" , 
			    "Name of the skimmed MCParticle  output collection"  ,
			    _mcParticlesSkimmedName ,
			    std::string("MCParticlesSkimmed") ) ;
    
    
  registerProcessorParameter( "daughtersECutMeV" , 
			      "energy cut for daughters that are kept from KeepDaughtersPDG"  ,
			      _eCutMeV,
			      float( 10. )  
			      ) ;
  

  registerProcessorParameter( "SaveBremsstrahlungPhotons" , 
			      "save photons from Brems"  ,
			      _saveBremsstrahlungPhotons,
			      bool(false)  
			      ) ;
  
    
  registerProcessorParameter( "BremsstrahlungEnergyCut" , 
			      "energy cut for Brems that are kept"  ,
			      _bremsstrahlungEnergyCut,
			      float( 1. )  
			      ) ;
 

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

  streamlog_out( DEBUG4 ) << " processEvent "  <<  evt->getEventNumber() << "  - " << evt->getRunNumber() 
			 << std::endl ; 

  LCCollection* mcpCol = evt->getCollection( _mcParticleCollectionName ) ;

  //    LCCollection* tHitRelCol = evt->getCollection( _trackHitRelationName ) ;
  


  // find track to MCParticle relations

    LCCollection* trackCol = 0 ;
    LCCollection* ttrlcol  = 0 ;
    bool haveTracks = true  ;

    try{ trackCol = evt->getCollection( _trackCollectionName ) ;  }   catch(...){  haveTracks=false ; } 

    if( ! haveTracks ) {
      streamlog_out( MESSAGE ) << " Track collection : " << _trackCollectionName 
			       << " not found - cannot create relation " << std::endl ;


    } else {  //if ( haveTracks ) {
      
      trackLinker(  mcpCol ,  trackCol  ,  &ttrlcol );
      
      if (_OutputTrackRelation ) 
	evt->addCollection(  ttrlcol  , _trackMCTruthLinkName  ) ;
      
    }
  
  // find cluster to MCParticle relations and the updated calohit to MCParticle relations.

  LCCollection* clusterCol = 0 ;
  LCCollection* cHitRelCol = 0 ;
  
  bool haveClusters = true ;
  bool haveCaloHitRel = true ;

  LCCollection* ctrlcol = 0;
  LCCollection* chittrlcol = 0;

  try{ clusterCol = evt->getCollection( _clusterCollectionName ) ; }   catch(...){  haveClusters=  false ; } 
  
  if( ! haveClusters ) {
    streamlog_out( MESSAGE ) << " Cluster collection : " << _clusterCollectionName 
			     << " not found - cannot create relation " << std::endl ;
  }
  
  
  try{ cHitRelCol = evt->getCollection( _caloHitRelationName ) ;   }   catch(...){  haveCaloHitRel = false ; } 
  
  if( ! haveCaloHitRel ) {
    streamlog_out( MESSAGE ) << " CaloHit relation : " << _caloHitRelationName 
			     << " not found - cannot create relation " << std::endl ;
  }
  
  if( haveTracks && haveClusters && haveCaloHitRel ) {
    
    clusterLinker(  mcpCol ,  clusterCol,  cHitRelCol ,  &ctrlcol , &chittrlcol);

    if (_OutputClusterRelation ) evt->addCollection(  ctrlcol  , _clusterMCTruthLinkName  ) ;
  
    if (_OutputCalohitRelation ) evt->addCollection(  chittrlcol  , _calohitMCTruthLinkName ) ;
  }
  
  
  
  
  // combine track and cluster to MCParticle relations to the reconstructed particle
  // to MCParticle relation
  LCCollection*  particleCol = 0 ;
  bool haveRecoParticles = true ; 
  
  try { particleCol = evt->getCollection(  _recoParticleCollectionName ); }   catch(...){ haveRecoParticles = false ; } 

  if( ! haveRecoParticles ) {
    streamlog_out( MESSAGE ) << " ReconstructedParticle collection : " << _recoParticleCollectionName 
			     << " not found - cannot create relation " << std::endl ;
  }
  
  LCCollection* ptrlcol = 0;

  if( haveRecoParticles &&  haveTracks && haveClusters && haveCaloHitRel ) {

    particleLinker(   particleCol, ttrlcol,  ctrlcol, &ptrlcol);

    evt->addCollection(  ptrlcol  , _recoMCTruthLinkName  ) ;
  }

  if( haveTracks && haveClusters && haveCaloHitRel ) {

    LCCollectionVec* skimVec = new LCCollectionVec( LCIO::MCPARTICLE )  ;

    makeSkim(    mcpCol , ttrlcol,  ctrlcol , &skimVec );
    evt->addCollection(  skimVec , _mcParticlesSkimmedName ) ;
  }
  //If either collection has not been added to the event, we have to delete it now!
  //Don't delete them before, because they are used
  if(!_OutputClusterRelation) { delete ctrlcol; }
  if(!_OutputCalohitRelation) { delete chittrlcol; }

}


void RecoMCTruthLinker::trackLinker(  LCCollection* mcpCol ,  LCCollection* trackCol,  LCCollection** ttrlcol) { 

  

  LCRelationNavigator trackTruthRelNav(LCIO::TRACK , LCIO::MCPARTICLE  ) ;


  // loop over reconstructed tracks
  int nTrack = trackCol->getNumberOfElements() ;

  int ifoundch =0;

  for(int i=0;i<nTrack;++i){

    Track* trk = dynamic_cast<Track*> ( trackCol->getElementAt(i) ) ;

      // charged particle  or V0  -> analyse track. We need to find all seen hits
      // this track is made of, wich sim hits each of the seen hits came from,
      // and finally which true particle actually created each sim hit
	
    MCPMap mcpMap ;  // mcpMap is a map seen <-> true particle
      
    int nHit = 0 ;

    const TrackerHitVec& trkHits = trk->getTrackerHits() ;
	
    for( TrackerHitVec::const_iterator hitIt = trkHits.begin() ;
	     hitIt != trkHits.end() ; ++hitIt ) { 
	  
      const TrackerHit* hit = * hitIt ; // ... and a seen hit ... 
	  
      const LCObjectVec& simHits  = hit->getRawHits() ;
	  
      for( LCObjectVec::const_iterator objIt = simHits.begin() ;
	       objIt != simHits.end() ; ++objIt ){
	    
	SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>( *objIt ) ; // ...and a sim hit ...
	    
	MCParticle* mcp = simHit->getMCParticle() ; // ... and a true particle !
	    
        if ( mcp != 0 ) {
	   mcpMap[ mcp ]++ ;   // count the hit caused by this true particle
	} else {
	  streamlog_out( WARNING ) << " tracker SimHit without MCParticle ?!   " <<  std::endl ;
        }

	++nHit ; // total hit count

      }
    }

  
    if( nHit == 0 ){
	
      streamlog_out( WARNING ) << " no simulated tracker hits found " <<  std::endl ;
	
      continue ;  // won't find a particle 
    }
      
      
      // find the mc particle with the largest #hits 
      // also store all genstat=1 particles and 
      // all particles of type to be saved

    MCParticle* mother = 0;
    MCParticleVec theMCPs ;  // vector that will contain all true particles contributing
    theMCPs.reserve( 1000 ) ;
    std::vector<int> MCPhits;
    MCPhits.reserve( 1000 ) ;
    int ifound = 0;

    for( MCPMap::iterator it = mcpMap.begin() ;   // iterate trough the map, map->first is the
       it != mcpMap.end() ; ++it ){             // true particle, map->second is the number of
                                                // times this true particle got mapped, ie. the
                                                // number of hits it produced.

      mother= dynamic_cast<MCParticle*>(it->first->getParents()[0]); // mother of the true partic.

      if (it->first->getGeneratorStatus() == 1 ) {  // genstat 1 particle, ie. it is a bona fide
                                                  // creating true particle: enter it into the list,
                                                  // and note how many hits it produced.
        theMCPs.push_back(it->first);  MCPhits.push_back(it->second); ifound++;
      } else {  // not genstat 1. Wat should we do with it ?
        if ( mother != 0 ) { // if it has a parent, save it
          theMCPs.push_back(it->first);  MCPhits.push_back(it->second); ifound++;
        } else {
          streamlog_out( WARNING ) << 
            " track has hit(s) from a non-generator particle with no parents ?!! " 
           << std::endl;
        }	    
      }
    } // end of loop over map

    // finally calculate the weight of each true partic to the seen 
    // (= hits_from_this_true/ all_hits),
    // and add the weighted reltion
    for (int iii=0 ; iii<ifound ; iii++ ) {

      float  weight = float(MCPhits[iii]  )/float(nHit) ; 

    
      streamlog_out( DEBUG4 ) << " track has " << MCPhits[iii]  
                            << " hits of " << nHit 
    			<< " [ " << int(weight) << " ] " 
    			<< " of MCParticle with pdg : " << theMCPs[iii]->getPDG() 
    			<< " and genstat : " << theMCPs[iii]->getGeneratorStatus() 
                            << " id: " << theMCPs[iii]
    		      << std::endl ;
    
    
      trackTruthRelNav.addRelation(   trk , theMCPs[iii] , weight ) ;
    }	

    ifoundch=ifound;
  }  
  //  seen-true relation complete. add the collection

  streamlog_out( DEBUG4 ) << " track linking complete, create collection " << std::endl;

  *ttrlcol = trackTruthRelNav.createLCCollection() ;
}



void RecoMCTruthLinker::clusterLinker(  LCCollection* mcpCol ,  LCCollection* clusterCol, 
                                         LCCollection* cHitRelCol , 
                                         LCCollection** ctrlcol, LCCollection** chittrlcol) { 

  

  
  LCRelationNavigator clusterTruthRelNav(LCIO::CLUSTER , LCIO::MCPARTICLE  ) ;
  LCRelationNavigator chitTruthRelNav(LCIO::CALORIMETERHIT , LCIO::MCPARTICLE  ) ;
  LCRelationNavigator cHitRelNav( cHitRelCol ) ;
  



  // loop over reconstructed particles
  int nCluster = clusterCol->getNumberOfElements() ;
  streamlog_out( DEBUG4 ) <<" Number of Clusters "<< nCluster <<std::endl;

  std::vector<Cluster*> missingMC ;
  missingMC.reserve( nCluster ) ;
  int ifoundclu =0;
    // now for the clusters
    	

  Remap_as_you_go remap_as_you_go ;  // map from true tracks linked to hits to 
                                         // those that really should have been linked
  MCParticle* mother = 0;
	
	

  for(int i=0;i<nCluster;i++){

    MCPMapDouble mcpEnergy ;
    double eTot = 0 ;
    Cluster* clu = dynamic_cast<Cluster*> ( clusterCol->getElementAt(i) ) ;



	
      // We need to find all seen hits this clutser is made of, wich sim hits each 
      // of the seen hits came from, and finally which true particles actually created 
      // each sim hit. Contrary to the sim tracker hits above, a sim-calo hit can be
      // made by several true particles. They also have a signal size (energy) value.
      // In addition, the true particle creating sometimes needs to be back-tracked
      // to the particle actually entering tha calorimeter. 

	  
	  
    streamlog_out( DEBUG4 ) <<"Cluster clu = "<< clu << " (i = " << i << " )" << std::endl;

    const CalorimeterHitVec& cluHits = clu->getCalorimeterHits() ;
    double ecalohitsum=0.;	  
    for( CalorimeterHitVec::const_iterator hitIt = cluHits.begin() ;
	     hitIt != cluHits.end() ; ++hitIt ) { 
	    

      CalorimeterHit* hit = * hitIt ;  // ... a calo seen hit ...
      streamlog_out( DEBUG1 ) <<"   hit = "<< hit << " e " << hit->getEnergy()<< std::endl;
      ecalohitsum+= hit->getEnergy();	    
          // (MB: the below was commented out by Frank, and replaced by the navigator
          // below. I don't kow if this is due to a bug that might be fixed now. If so
          // getRawHit could be used, and the bit below would be simplified to work
          // the ame way as it does for the tracks above)
	  // 	     SimCalorimeterHit* simHit = 
          //             dynamic_cast<SimCalorimeterHit*>( hit->getRawHit() ) ;
	  // 	     streamlog_out( DEBUG4 ) << "   raw hit ptr : " << simHit << std::endl ;

        
      const LCObjectVec& simHits = cHitRelNav.getRelatedToObjects( (LCObject*) hit )  ;
      double ehit = 0.0; 
      int nsimhit = 0;        
      for( LCObjectVec::const_iterator objIt = simHits.begin() ;
           objIt != simHits.end() ; ++objIt ){
          
        SimCalorimeterHit* simHit = dynamic_cast<SimCalorimeterHit*>( *objIt ) ; // ... and a sim hit ....

        streamlog_out( DEBUG1 ) <<"      simhit = "<< simHit << std::endl;
        nsimhit++;
        for(int j=0;j<simHit->getNMCContributions() ;j++){
            
          MCParticle* mcp = simHit->getParticleCont( j ) ;
          double e  = simHit->getEnergyCont( j ) ;
          streamlog_out( DEBUG1 ) <<"         true contributor = "<< mcp << " e: " << e <<std::endl;
          if ( remap_as_you_go.find(mcp) != remap_as_you_go.end() ) {
            mcp=remap_as_you_go.find(mcp)->second;
          } else {
            MCParticle* this_Kid = mcp ;     // ... and a true particle !
          
	    //Particle gun particles dont have parents, but are "created in the simulation" and have genStat 0
            if ( mcp-> getGeneratorStatus() == 0 && ( mcp->getParents().empty() == false ) ) { 
                  // not from generator, find which true particle this
                  // hit should really be attributed to, by tracking back the history.

                  // Two cases to treat:
                  //   For some reason, the hit is not attributed to the
                  //   incomming particle, but to some particle in the shower.
                  //   Just track back to the incomming particle, which might (usually)
                  //   be a gen-stat 1 particle from the main vertex. It can also
                  //   be from a decay or interaction in the tracking made by
                  //   Geant (genstat=0, mother isDecayedInTracker), or a
                  //   decyed particle from the generator (genstat 2) that
                  //   hits the calo before decaying (lambdas, K^0_S)

                  //   Or: for back-scatterers the calorimiter hit is attributed to the
                  //   the last particle *even if this particle both
                  //   started and ended in the tracker* !!!! Then we back-track
                  //   untill we find a particle which at least started inside
                  //   the calorimeter, and then go on backtracking as above.
                  //   This case is triggered by the particle linked to the
                  //   hit being DecayedInTracker, hence the case where a
                  //   back-scatter actually ends in the calorimeter is
                  //   treated as the first case.


              streamlog_out( DEBUG3 ) << " hit " <<hit<<","<< j << 
                                         " not created by generator particle. backtracking ..." 
                                      << std::endl;
              streamlog_out( DEBUG2 ) <<"       "<<mcp<<" gs "<<mcp->getGeneratorStatus()<<
                                        " dint "<<mcp->isDecayedInTracker()<<
                                        " bs "<<mcp->isBackscatter()<<
                                        " npar "<<mcp->getParents().size()<<
                                        " pdg "<<mcp->getPDG()<<
                                        " "<< mcp->getVertex()[0]<<
                                        " "<<mcp->getVertex()[1]<<
                                        " "<<mcp->getVertex()[2]<<std::endl;

              mother= dynamic_cast<MCParticle*>(mcp->getParents()[0]); 

	      bool NoChange=true;
              while ( this_Kid->isDecayedInTracker()  && 
                           mother->getGeneratorStatus() == 0){ // backtrack until the particle is 
                                                               // in a calorimeter. 
		NoChange=true;//If nothing happens in the loop, we are going to break, to make sure we are not looping endlessly
                // case back-scatter

                streamlog_out( DEBUG2 ) <<"  back-scatter o mother "<<mother<<
                                          " gs "<<mother->getGeneratorStatus()<<
                                          " dint "<<mother->isDecayedInTracker()<<
                                          " bs "<<mother->isBackscatter()<<
                                          " npar "<<mother->getParents().size()<<
                                          " pdg "<<mother->getPDG()<<std::endl;
                streamlog_out( DEBUG2 ) <<"  back-scatter o this_Kid "<<this_Kid<<
                                          " gs "<<this_Kid->getGeneratorStatus()<<
                                          " dint "<<this_Kid->isDecayedInTracker()<<
                                          " bs "<<this_Kid->isBackscatter()<<
                                          " npar "<<this_Kid->getParents().size()<<
                                          " pdg "<<this_Kid->getPDG()<<std::endl;

                while ( mother!= 0 &&   mother->getParents().size()>0 && 
                        !mother->isBackscatter() && 
                        mother->getGeneratorStatus() == 0 ) { 
                                      // back-track until the parent is not a back-scatter .
                                      // generator paricles are never back-scatterers, so
                                      // once genstat != 0, we know there will be no
                                      // back-scatterers furher down the tree.
                  this_Kid=mother ;
                  mother= dynamic_cast<MCParticle*>(mother->getParents()[0]); // (assume only one...)
		  NoChange=false;

                  streamlog_out( DEBUG2 ) <<"  back-scatter i mother "<<mother<<
                                            " gs "<<mother->getGeneratorStatus()<<
                                            " dint "<<mother->isDecayedInTracker()<<
                                            " bs "<<mother->isBackscatter()<<
                                            " npar "<<mother->getParents().size()<<
                                            " pdg "<<mother->getPDG()<<std::endl;
 
                }
                if ( mother->isBackscatter()) {
                  this_Kid=mother ;
                  mother= dynamic_cast<MCParticle*>(mother->getParents()[0]); // (assume only one...)
		  NoChange=false;
                }
		//If mother genStatus==0 and it is not backscattering and there are no parents (because mother is coming from particle gun), we may never get out of this loop
		if(NoChange) break;
              }
              streamlog_out( DEBUG2 ) <<"  mother "<<mother<<
                                        " gs "<<mother->getGeneratorStatus()<<
                                        " dint "<<mother->isDecayedInTracker()<<
                                        " bs "<<mother->isBackscatter()<<
                                        " npar "<<mother->getParents().size()<<
                                        " pdg "<<mother->getPDG()<<std::endl;

              while ( mother!= 0 &&   mother->getParents().size()>0 && 
                       mother->getGeneratorStatus() ==0 &&
                       !mother->isDecayedInTracker() ) { 
                                    // back-track as long as there is a non-generator 
                                    // mother, or the mother decayed in the tracker 
                                    // (=> the kid is the particle entering the calorimeter.)

                // case shower-particle
                this_Kid=mother ;
                mother= dynamic_cast<MCParticle*>(mother->getParents()[0]); // (assume only one...)

                streamlog_out( DEBUG2 ) <<"  shower-part mother "<<mother<<
                                          " gs "<<mother->getGeneratorStatus()<<
                                          " dint "<<mother->isDecayedInTracker()<<
                                          " bs "<<mother->isBackscatter()<<
                                          " npar "<<mother->getParents().size()<<
                                          " pdg "<<mother->getPDG()<<std::endl;


              }
              if ( mother->isDecayedInTracker() ) {
                remap_as_you_go[mcp]=this_Kid;
                mcp=this_Kid; // this_Kid started before the calo, and is the one 
                              // the hit should be attributed to

                streamlog_out( DEBUG3 ) << "   attributed to " << mcp << 
                                           " because it's origin is in tracker :" <<std::endl;

              } else {
                remap_as_you_go[mcp]=mother;
                mcp=mother;   // mother started at ip, and reached the calo, and is the one 
                              // the hit should be attributed to.

                streamlog_out( DEBUG3 ) << "   attributed to " << mcp << 
                                           " because it is a generator particle : "<< std::endl;
              }

              streamlog_out( DEBUG3 ) <<"       gs "<<mcp->getGeneratorStatus()<<
                                        " dint "<<mcp->isDecayedInTracker()<<
                                        " bs "<<mcp->isBackscatter()<<
                                        " npar "<<mcp->getParents().size()<<
                                        " pdg "<<mcp->getPDG()<<
                                        " "<<mcp->getEndpoint()[0]<<
                                        " "<<mcp->getEndpoint()[1]<<
                                        " "<<mcp->getEndpoint()[2]<<std::endl;
            }
          }  // genstat if - then - else
//        if( e == 0.0 )
//          streamlog_out( WARNING ) << " zero energy in MCContribution " 
//                                   << mcp->getPDG()  << std::endl ;

          chitTruthRelNav.addRelation(  hit , mcp , float(e) ) ;
          mcpEnergy[ mcp ] +=  e ;// count the hit-energy caused by this true particle
          eTot += e ;             // total energy
          ehit+= e;
        } // mc-contributon-to-simHit loop
      } // simHit loop 
      streamlog_out( DEBUG2 )<< "   summed contributed e: " << ehit << " ratio : " << ehit/hit->getEnergy()
                             << " nsimhit " << nsimhit <<std::endl;
    } // hit loop
    streamlog_out( DEBUG3 ) << " Sum of calohit E: " <<  ecalohitsum << " cluster E " 
                            << clu->getEnergy() << " Sum of Simcalohit E: " << eTot << std::endl;
    if( eTot == 0.0 ){ // fixme - this might happen if clusters are from Lcal/Muon only
	  

		  // save reco particle in missingMC 
      missingMC.push_back( clu  ) ;
		  
      streamlog_out( DEBUG4 ) << " no calorimeter hits found for " 
			       << " cluster particle: e:"  << clu->getEnergy()  
		               //   << " charge: " << rec->getCharge() 
		               //   << " px: " << rec->getMomentum()[0]
                               //   << " py: " << rec->getMomentum()[1]
                               //   << " pz: " << rec->getMomentum()[2]
		               //   << " mass: " << rec->getMass() 
		               //   <<  " *pdg " << rec->getParticleIDUsed()
			       << std::endl ;
	
	  
      continue ;  
    }

    // At this point, we have a list of all true particles contributiong to this cluster.
    // We now need to sum up the energies each particle contributes with and do
    // further parsing.
  
  
    // find the mc particle with the largest energy contribution 
    // also store all genstat=1 particles and 
    // all particles of type to be saved
      
    double eMax = 0 ;          // energy the most contributing true particle added to 
                               // the cluster
  
    MCParticleVec theMCPs ;    // vector that will contain all true particles contributing,
                               // if they, ex officio, will be in the skimmed collection.
    theMCPs.reserve(1000);
    std::vector<double> MCPes; // energy contribution of these
    MCPes.reserve(1000);
    int ifound = 0;            // total number of contributors
  
    MCParticleVec moreMCPs ;   // vector that will contain true particles contributing, that
                               // normally wouldn't be in the skimmed collection
    moreMCPs.reserve(1000);
    std::vector<double>moreMCPes;  // energy contribution of these
    moreMCPes.reserve(1000);
    int morefound = 0;         // number of such cases.
  
    mother = 0;
  
    for( MCPMapDouble::iterator it = mcpEnergy.begin() ;  // iterate trough the map.
         it != mcpEnergy.end() ; ++it ){
      if (it->first->getGeneratorStatus() == 1 ) {  // genstat 1 particle, ie. it is a bona fide
                                                    // creating true particle: enter it into 
                                                    // the list, and note how much energy it 
                                                    // contributed.
        theMCPs.push_back(it->first);  MCPes.push_back(it->second); ifound++;
      } else { // not genstat 1. What should we do with it ?
         if (  it->first->getParents().size() != 0 ) { 
           mother= dynamic_cast<MCParticle*>(it->first->getParents()[0]);
         } else { 
           mother = 0 ; 
         }
         if ( mother != 0 ) {
           if ( mother->isDecayedInTracker() &&            // ... so the partic apeared in the 
                                                           // tracker ...
                _pdgSet.find( abs (mother->getPDG())) != _pdgSet.end() ) {  
                                                           // ... and is of a type we want to
                                                           // save (it's mother is in _pdgSet) 
                                                           // -> also a bona fide creator.
             theMCPs.push_back(it->first);  MCPes.push_back(it->second); ifound++;
           } else { // else: add to the  moreMCPs-list, further treated below
             streamlog_out( DEBUG2 ) << " case 1 for "<< it->first << 
                                        "(morefound=" << morefound << ")" << 
                                        " mother: " << mother <<
                                        " gs "  <<mother->getGeneratorStatus()<< 
                                        " dint " << mother->isDecayedInTracker() <<
                                        " bs "  << mother->isBackscatter() << 
                                        " pdg " <<mother->getPDG() <<std::endl;
             moreMCPs.push_back(it->first);  moreMCPes.push_back(it->second); morefound++;
           }
         } else { // not genstat 1, no mother ?! Also add to the  moreMCPs-list.
           streamlog_out( WARNING ) << " case 2 for "<< it->first << 
                                       "(morefound=" << morefound << ")" << std::endl;
           moreMCPs.push_back(it->first);  moreMCPes.push_back(it->second); morefound++;
         }	    
      }
      
  
      if( it->second > eMax  ){
  
        eMax = it->second ;
      }
    }
  
    // We now have two lists of contributing true particles, one containing those that will be
    // put to the skimmed list in any case, one with those that will be put there only is they
    // produced a visible signal. In addition, we have the total contribution to the shower from
    // each true particle in separate lists.
    
    if ( morefound > 0 ) {
      for (int iii=0 ; iii<morefound ; iii++ ) {
        streamlog_out( DEBUG2 ) << " iii, moreMCPes[iii], moreMCPs[iii], gs " << iii <<
                                   " "<< moreMCPes[iii] <<
                                   " "<< moreMCPs[iii]<<
                                   " "<<moreMCPs[iii]->getGeneratorStatus() << std::endl;
      }
      for (int iii=0 ; iii<ifound ; iii++ ) {
        streamlog_out( DEBUG2 ) << " iii, MCPes[iii], theMCPs[iii] " << iii <<
                                   " "<< MCPes[iii] <<" "<< theMCPs[iii] << std::endl;
      }
      streamlog_out( DEBUG3 )<< "   morefond: " << morefound <<std::endl;
    }
    
    // figure out what to do with the cases where true particles not automatically in the
    // skimmed list contributed to the cluster: Normally, the back-tracking above has
    // done the job, but sometimes one has a neutron or a photon playing pin-ball between
    // the calorimeter and the tracker, so that a back-scatter scatters back from the
    // tracker into the same cluster it came from. In that case, we should ignore the
    // back scatter (normally, of course, a back-scatterer is - if it reaches a 
    // calorimeter - in a different cluster and should be kept as an originator).
    // I couldn't figure out a more efficient way of figuring this out than the
    // almost-always-do-nothing loop below ;(
    
  
    mother = 0;
    for (int iii=0 ; iii<morefound ; iii++ ) {
      if (  moreMCPs[iii]->getParents().size() != 0 ) { 
        mother= dynamic_cast<MCParticle*>(moreMCPs[iii]->getParents()[0]); 
        streamlog_out( DEBUG2 ) << "   iii: " << iii << " mother: " << mother  <<std::endl; 
      } else { 
        mother = 0 ; 
      }
      streamlog_out( DEBUG2 ) <<"      partic vert: "<<moreMCPs[iii]->getVertex()[0]<<
                                " "<<moreMCPs[iii]->getVertex()[1]<<
                                " "<<moreMCPs[iii]->getVertex()[2]<<std::endl;
      streamlog_out( DEBUG2 ) <<"      partic endpoint: "<<moreMCPs[iii]->getEndpoint()[0]<<
                                " "<<moreMCPs[iii]->getEndpoint()[1]<<
                                " "<<moreMCPs[iii]->getEndpoint()[2]<<std::endl;
      streamlog_out( DEBUG2 ) <<"      partic status : gs "<<moreMCPs[iii]->getGeneratorStatus()<<
                                " dint "<<moreMCPs[iii]->isDecayedInTracker ()<<
                                " bs "<<moreMCPs[iii]->isBackscatter ()<< 
                                " npar "<<moreMCPs[iii]->getParents().size()<< 
                                " pdg " << moreMCPs[iii]->getPDG() << 
                                " dinc "<<moreMCPs[iii]->isDecayedInCalorimeter ()<<
                                " bye "<<moreMCPs[iii]->hasLeftDetector () <<
                                " stop "<<moreMCPs[iii]->isStopped () <<  std::endl;
  
      while ( mother!= 0 &&  mother->getGeneratorStatus() !=2 ) { // back-track to the 
                                                                  // beginning of the chain
  
        streamlog_out( DEBUG2 ) << "       mother "<< mother << 
                                   " gs " << mother->getGeneratorStatus() << 
                                   " dint " << mother->isDecayedInTracker()<< " " <<
                                   " bs "   << mother->isBackscatter()<<
                                   " pdg "  << mother->getPDG() << 
                                   " dinc " <<mother->isDecayedInCalorimeter()<< std::endl;
        streamlog_out( DEBUG2 ) <<"       mother vert: "<<mother->getVertex()[0]<<
                                  " "<<mother->getVertex()[1]<<
                                  " "<<mother->getVertex()[2]<<std::endl;
        streamlog_out( DEBUG2 ) <<"       mother endpoint: "<<mother->getEndpoint()[0]<<
                                  " "<<mother->getEndpoint()[1]<<
                                  " "<<mother->getEndpoint()[2]<<std::endl;
  
        
        // find out if this ancestor actually directly gives rise to hits in this cluster. 
        // If so, we attribute the hits of  moreMCPs[iii] to that true particle instead.
  
        // (this ought to work, but it doesn't:)
        //aa	   for (MCParticleVec::iterator itfind=theMCPs.begin();  
        //aa                               itfind != theMCPs.end() ; ++itfind) 
        //aa  MCParticle* tmcp = dynamic_cast<MCParticle*>(*itfind);
   
        for (int kkk=0 ; kkk<ifound ; kkk++){
          MCParticle* tmcp = theMCPs[kkk];
          if ( tmcp == mother ) {
            streamlog_out( DEBUG3 ) << "        found " << moreMCPs[iii] << 
                                       "(iii= "<<iii <<")" << kkk << 
                                       " to be related to "<<mother <<
                                       " add e : " <<  moreMCPes[iii] << std::endl;
            MCPes[kkk]+= moreMCPes[iii];  
            goto endwhile;
          }
        }   
        if (  mother->getParents().size() != 0 ) { 
  
           mother= dynamic_cast<MCParticle*>(mother->getParents()[0]); 
        } else { mother = 0; }
      }
      endwhile:
      if ( mother == 0 || mother->getGeneratorStatus() ==2 ) {
  
        // no other contributing true particle found among the ancestors 
        // to moreMCPs[iii], so we add it to the
        // list of true particles to be saved.
  
        streamlog_out( DEBUG3 ) << "        No relation found for "<< moreMCPs[iii] << 
                                   ". Keep it as separate originator "<< std::endl;
        theMCPs.push_back(moreMCPs[iii]);  MCPes.push_back(moreMCPes[iii]); ifound++;
      }
  
    }
    streamlog_out( DEBUG4 ) << " ifound cluster " << ifound <<  std::endl;
  
    // finally calculate the weight of each true partic to the seen 
    // (= energy_from_this_true/ total ), and add the weighted reltion.
  
    float totwgt =0.0;
    for (int iii=0 ; iii<ifound ; iii++ ) {

      float  weight = MCPes[iii]/eTot;
      
      totwgt+= weight;
      if( theMCPs[iii] == 0 ) {
  
        streamlog_out( ERROR ) << " cluster has " << MCPes[iii] << " GeV of " << eTot 
      			 << " GeV [ " << int(weight) << " ] true energy " 
      		         << " but no MCParticle " << std::endl ;
        
        continue ;
  
      }
  
      streamlog_out( DEBUG4 ) << " cluster has " <<  MCPes[iii] << " GeV of " << eTot 
      			<< " GeV [ " << int(weight) << " ] " 
      			<< " of MCParticle with pdg : " << theMCPs[iii]->getPDG() 
      			<< " and genstat : " << theMCPs[iii]->getGeneratorStatus() 
                              << " id: " << theMCPs[iii]
      		      << std::endl ;
      
      //Particle gun particles dont have parents but have genStat 0      
      if (  theMCPs[iii]->getGeneratorStatus() == 0 &&  ( theMCPs[iii]->getParents().empty() == false) ) {
        streamlog_out( DEBUG4 ) << " mother id " << theMCPs[iii]->getParents()[0] 
         << " and genstat " 
         << theMCPs[iii]->getParents()[0]->getGeneratorStatus() << std::endl ;
  
      }
      
      // add relation. 
  
      clusterTruthRelNav.addRelation(   clu , theMCPs[iii] , weight ) ;
      
    }
    ifoundclu=ifound;
  } // cluster loop
  

  // recover missing MCParticles for neutrals :
  // attach the MCParticle with smallest angle to cluster

  // (This is from the original RecoMCTrueLinker. I didn't revise it/MB)

  for( unsigned i=0 ; i < missingMC.size() ; ++i ) {
	 
    int nMCP  = mcpCol->getNumberOfElements() ;
     
    Cluster* clu = missingMC[i] ;
	 

    gear::Vector3D recP( clu->getPosition()[0] , clu->getPosition()[1] ,
                              clu->getPosition()[2] ) ;
	 
    double recTheta = recP.theta() ;
        
    double maxProd = 0.0 ;
    MCParticle* closestMCP = 0 ;
        
    for (int j=0; j< nMCP ; j++){
       	 
      MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( j ) ) ;
       	  
      if ( fabs( mcp->getCharge() ) > 0.01 ) {
        continue ;
      }

      gear::Vector3D mcpP( mcp->getMomentum()[0] , mcp->getMomentum()[1] ,
                                    mcp->getMomentum()[2] ) ;

      if ( fabs( recTheta - mcpP.theta() ) > 0.3 ) {// fixme : proc param...
        continue ;
      }

 
      double prod  = mcpP.unit().dot(  recP.unit() ) ;

      if ( prod > maxProd ) {
        maxProd = prod ;
        closestMCP = mcp ;
      }
    }
    if ( maxProd > 0. ) {
       
      streamlog_out( DEBUG4 ) << "  neutral cluster particle recovered"  
                       << clu->getEnergy()
                       //<< " px: " << rec->getMomentum()[0]
       		//<< " py: " << rec->getMomentum()[1]
       		//<< " pz: " << rec->getMomentum()[2]
       		<< " maxProd: " << maxProd 
                       << " px: " << closestMCP->getMomentum()[0]
       		<< " py: " << closestMCP->getMomentum()[1]                                           
       		<< " pz: " << closestMCP->getMomentum()[2]
       		<< std::endl ;
       	  
       clusterTruthRelNav.addRelation(   clu , closestMCP ,  maxProd ) ;
    }


  }
  //  seen-true relation complete. add the collection

  streamlog_out( DEBUG4 ) << " cluster linking complete, create collection " << std::endl;
  *ctrlcol = clusterTruthRelNav.createLCCollection() ;
  *chittrlcol = chitTruthRelNav.createLCCollection() ;
} 


void RecoMCTruthLinker::particleLinker(  LCCollection* particleCol, LCCollection* ttrlcol, 
                                          LCCollection* ctrlcol, LCCollection** ptrlcol) {

  LCRelationNavigator particleTruthRelNav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE  ) ;

  LCRelationNavigator      trackTruthRelNav = LCRelationNavigator(  ttrlcol );
  LCRelationNavigator      clusterTruthRelNav = LCRelationNavigator(  ctrlcol );

  LCObjectVec mcvec; 
  int nPart = particleCol->getNumberOfElements() ;


  static FloatVec www ;

  for(int i=0;i<nPart;++i){

    std::map< MCParticle* , int > mcmap;

    MCParticle* mcmax =0;
    float maxwgt=0. ;

    ReconstructedParticle* part = dynamic_cast<ReconstructedParticle*> ( particleCol->getElementAt(i) ) ;
    TrackVec tracks=part->getTracks();
    int ntrk =  part->getTracks().size();
    ClusterVec clusters=part->getClusters();
    int nclu =  part->getClusters().size();
    streamlog_out( DEBUG4 ) << " Treating particle " << part << " with index " << i 
                            << " it has " << ntrk << " tracks, and " << nclu << " clusters " <<std::endl;  
    int nhit[100] ;
    int nhitT = 0;
    if ( ntrk > 1 ) {

      for (int j=0 ; j < ntrk ; j++ ) {
        nhit[j] = 0;
        for ( unsigned kkk=0 ;kkk<tracks[j]->getSubdetectorHitNumbers().size(); kkk++ ) {
          nhit[j]+= tracks[j]->getSubdetectorHitNumbers()[kkk];
          nhitT+= tracks[j]->getSubdetectorHitNumbers()[kkk];
        }
        streamlog_out( DEBUG2 )  << " Track " << j <<" has " << nhit[j] << " hits " << std::endl; 
      }
      streamlog_out( DEBUG2 )  << " Total : " << nhitT << " hits " << std::endl; 
    } else {
      nhit[0]=1 ; nhitT=1;
    }
    for (int j=0 ; j < ntrk ; j++ ) {
  
      if (  tracks[j] != 0 ) {  
        mcvec = trackTruthRelNav.getRelatedToObjects( tracks[j]);
        www = trackTruthRelNav.getRelatedToWeights( tracks[j]);
        int ntp= mcvec.size();
        streamlog_out( DEBUG3 ) << "    Track " <<  tracks[j] << " with index " << j 
                                << " has " << ntp << " true particles " << std::endl;  
        if ( ntp > 0 ) {
          if ( mcvec[0] != 0 ) {

            for ( int kkk=0 ; kkk < ntp ; kkk++ ) {
              if ( !_FullRecoRelation && www[kkk]*(float(nhit[j])/float(nhitT)) > maxwgt ) {
                maxwgt= www[kkk]*(float(nhit[j])/float(nhitT)) ;
                mcmax=dynamic_cast<MCParticle*>(mcvec[kkk]);
              }
              mcmap[dynamic_cast<MCParticle*>(mcvec[kkk])] += 
                int(www[kkk]*1000.*(float(nhit[j])/float(nhitT))+0.5);
              streamlog_out( DEBUG2 ) << "    Individual track weight to " <<mcvec[kkk]<< " is " 
                                      <<  www[kkk] << ", scaled one is "
                                      <<  www[kkk]*(float(nhit[j])/float(nhitT))
                                      << " ( loop -index : " << kkk << ")"<< std::endl; 
            }
          }
        }
      }
    }
    if (  _FullRecoRelation || ntrk == 0 ) {
      double eclu[100] ;
      double ecluT = 0.;
      if ( nclu > 1 ) {
  
        for (int j=0 ; j < nclu ; j++ ) {
          eclu[j] = clusters[j]->getEnergy();
          ecluT+=clusters[j]->getEnergy();
          streamlog_out( DEBUG2 )  << " Cluster " << j <<" has energy " << eclu[j] << std::endl; 
        }
        streamlog_out( DEBUG2 )  << " Total : " << ecluT << std::endl; 
      } else {
        eclu[0]=1 ; ecluT=1;
      }
      for (int j=0 ; j < nclu  ; j++ ) {
        if ( clusters[j] != 0 ) {
          mcvec =  clusterTruthRelNav.getRelatedToObjects(clusters[j]);
          www = clusterTruthRelNav.getRelatedToWeights(clusters[j]);
          int ntp= mcvec.size();
          streamlog_out( DEBUG3 ) << "    Cluster " <<  clusters[j] << " with index " << j
                                << " has " << ntp << " true particles " << std::endl;  
          if ( ntp > 0 ) {
            if ( mcvec[0] != 0 ) {
              for ( int kkk=0 ; kkk < ntp ; kkk++ ) {
                if ( !_FullRecoRelation &&  www[kkk]*(eclu[j]/ecluT) > maxwgt ) {
                   maxwgt= www[kkk]*(eclu[j]/ecluT) ;
                   mcmax=dynamic_cast<MCParticle*>(mcvec[kkk]);
                }
                mcmap[dynamic_cast<MCParticle*>(mcvec[kkk])] += 
                     int(www[kkk]*1000.*(eclu[j]/ecluT)+0.5)*10000;
                streamlog_out( DEBUG2 ) << "    Individual cluster Weight to " <<mcvec[kkk]<< " is " 
                                      <<  www[kkk] << ", scaled one is "
                                      <<  www[kkk]*(eclu[j]/ecluT)
                                      << " ( loop -index : " << kkk << ")"<< std::endl; 
              }
            }
          }
        }
      }
    }
    if ( _FullRecoRelation ) {
      for ( std::map< MCParticle* , int >::iterator mcit = mcmap.begin() ; 
            mcit !=  mcmap.end() ; mcit++ ) { 
         // loop all MCparticles releted to the particle 
         // get the true particle
        streamlog_out( DEBUG3 ) << " particle has weight "<<mcit->second
                                << " (Track: " << int(mcit->second)%10000 
                                << " , Cluster: " << int(mcit->second)/10000 << " ) " 
     			        << " of MCParticle with pdg : " << mcit->first->getPDG() 
    			        << " and genstat : " <<  mcit->first->getGeneratorStatus() 
                                << " id: " << mcit->first 
    		                << std::endl ;

        particleTruthRelNav.addRelation(   part ,  mcit->first ,  mcit->second ) ;
      }
    } else {
      if( mcmax != NULL ) {
	//AS: FixMe: There is still something going wrong with particles from the particle gun.
	//Although there are perfect PFOs no link with the MCParticle is established.
	  particleTruthRelNav.addRelation(   part ,  mcmax ,  maxwgt ) ;
	  streamlog_out( DEBUG3 ) << " particle has weight "<<maxwgt
     	  		        << " of MCParticle with pdg : " << mcmax->getPDG() 
    	  		        << " and genstat : " <<  mcmax->getGeneratorStatus() 
                                << " id: " << mcmax 
                                << ". Particle charge and ntracks : " << part->getCharge()<<" "<<ntrk
    	  	                << std::endl ;
      } else { 
      	  streamlog_out( WARNING ) << " particle has weight "<< maxwgt
      				  << ". Particle charge and ntracks : " << part->getCharge()<<" "<<ntrk
      				  << " but no mcparticle found "
      				  << std::endl ;
      }
    }
  }
  streamlog_out( DEBUG4 ) << " particle linking complete, create collection " << std::endl;
  *ptrlcol = particleTruthRelNav.createLCCollection() ;
}
void RecoMCTruthLinker::makeSkim(   LCCollection* mcpCol ,  LCCollection* ttrlcol,  LCCollection* ctrlcol ,  LCCollectionVec** skimVec){


  LCRelationNavigator      trackTruthRelNav = LCRelationNavigator(  ctrlcol );
  LCRelationNavigator      clusterTruthRelNav = LCRelationNavigator(  ttrlcol );


  //-------------- create skimmed MCParticle collection ------------------------

//  *skimVec = new LCCollectionVec( LCIO::MCPARTICLE )  ;
  (*skimVec)->setSubset( true) ;  // flag as subset 
    
    
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

      //  truthRelNav is the one we created above, remember, so here we make sure that all
      //  true particles related to seen ones (with the logic we used there) really will
      //  be in the skimmed collection !

      const LCObjectVec& trackObjects = trackTruthRelNav.getRelatedFromObjects( (LCObject*) mcp )  ;
      const LCObjectVec& clusterObjects = clusterTruthRelNav.getRelatedFromObjects( (LCObject*) mcp )  ;

      if( trackObjects.size() > 0 || clusterObjects.size()>0){

	streamlog_out( DEBUG4 ) << " keep MCParticle - e :" << mcp->getEnergy()  
				<< " charge: " << mcp->getCharge() 
				<< " px: " << mcp->getMomentum()[0]
				<< " py: " << mcp->getMomentum()[1]
				<< " pz: " << mcp->getMomentum()[2]
				<< " mass: " << mcp->getMass() 
				<< " pdg " << mcp->getPDG()
				<< std::endl ;

        // keepMCParticles also flags all parents of a kept particle, guaranteeing that
        // the history of any contributor to a detected signal will be kept !

	keepMCParticle( mcp ) ;
      }	

    } // else

  }   // end mcp loop 
  // --- loop again and add daughters of particles that are in the skimmed list and have a pdg in
  //     the parameter vector 'KeepDaughtersPDG 

  streamlog_out( DEBUG4 ) << " First loop done. Now search for KeepDaughtersPDG:s" << std::endl;


  if(_saveBremsstrahlungPhotons){
    for(int i=0; i< nMCP ; i++){
      MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
      if( abs(mcp->getPDG()) == 22 && mcp->getEnergy() > _bremsstrahlungEnergyCut){
	if( mcp->getParents().size() ){
	  MCParticle* parent = mcp->getParents()[0];
	  if( abs(parent->getPDG()) == 11){
	    const float x = fabs(mcp->getVertex()[0]);
	    const float y = fabs(mcp->getVertex()[1]);
	    const float z = fabs(mcp->getVertex()[2]);
	    const float r= sqrt(x*x+z*z);
	    const MCParticleVec& daughters = parent->getDaughters() ;
	    const float xp = fabs(parent->getVertex()[0]);
	    const float yp = fabs(parent->getVertex()[1]);
	    const float zp = fabs(parent->getVertex()[2]);
	    const float rp= sqrt(x*x+z*z);
	    const float xpe = fabs(parent->getEndpoint()[0]);
	    const float ype = fabs(parent->getEndpoint()[1]);
	    const float zpe = fabs(parent->getEndpoint()[2]);
	    const float dx= x-xpe;
	    const float dy= y-ype;
	    const float dz= z-zpe;
	    const float dr = sqrt(dx*dx+dy*dy+dz*dz);
	    if(dr>100.){
	      //std::cout << " Brem candidate " << mcp->getPDG() << " E = " << mcp->getEnergy() << std::endl; 
	      //std::cout << mcp->getParents().size() << " parent pdg : " << parent->getPDG() << " E = " << parent->getEnergy() << " Daughters : " << daughters.size() << std::endl; 
	      
	      //std::cout << " parent " << xp << "," << yp << " " << zp << " ->  " << xpe << "," << ype << " " << zpe << std::endl; 
	      keepMCParticle( mcp ) ;
	    }
	  }
	} 
      }
    }
  }


  for(int i=0; i< nMCP ; i++){
    
    MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    
    // keep the daughters of all decays in flight of particles in the pdg list (default: gamma, pi0, K0s) 
    if( mcp->ext<MCPKeep>() == true &&  mcp->isDecayedInTracker()  ){ //&& !mcp->isStopped()  ){  
      
      unsigned thePDG = abs( mcp->getPDG() ) ;
      
      if( _pdgSet.find( thePDG ) != _pdgSet.end()  ) {
	
	const MCParticleVec& daughters = mcp->getDaughters() ;
	
	streamlog_out( DEBUG4 ) << " keeping daughters of particle with pdg : " << mcp->getPDG() << " : " 
				<< " [" << mcp->getGeneratorStatus() << "] :";
	// 				  << " e :" << mcp->getEnergy() 
	// 				  << " isCreatedInSimulation :" << mcp->isCreatedInSimulation()	<< std::endl
	// 				  << " isBackscatter :" << mcp->isBackscatter()	<< std::endl
	// 				  << " vertexIsNotEndpointOfParent :" << mcp->vertexIsNotEndpointOfParent()	<< std::endl
	// 				  << " isDecayedInTracker :" << mcp->isDecayedInTracker()	<< std::endl
	// 				  << " isDecayedInCalorimeter :" << mcp->isDecayedInCalorimeter()	<< std::endl
	// 				  << " hasLeftDetector :" << mcp->hasLeftDetector()	<< std::endl
	// 				  << " isStopped :" << mcp->isStopped()    << "  : " 

	streamlog_message( DEBUG4 , 
			   if( mcp->getParents().size() ) ,
			   " parent pdg : " << mcp->getParents()[0]->getPDG() << " : "  ;
			   ) ;
	
	streamlog_out( DEBUG4 ) << std::endl ;
	
	//	<< std::endl ;
	
	for( MCParticleVec::const_iterator dIt = daughters.begin() ;
	     dIt != daughters.end() ; ++dIt ){
	  

	  MCParticle* dau = dynamic_cast<MCParticle*>( *dIt ) ;

	  if( dau->getEnergy()*1000. >  _eCutMeV ) {
	    
	    (*dIt)->ext<MCPKeep>() = true ;
	    
	    streamlog_out( DEBUG4 ) <<  (*dIt)->getPDG() << ", " ;
	  }
	}
	
	streamlog_out( DEBUG4 ) << std::endl ;
	
      }
    }
  }
  
  streamlog_out( DEBUG4 ) << " All found, add to skimmed list " << std::endl;

  for(int i=0; i< nMCP ; i++){
    
    MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    
    if( mcp->ext<MCPKeep>() == true ) {  
      
       (*skimVec)->addElement( mcp ) ;
    }
  }    
  
  // okidoki, the skimmed collection is complete. Add it.

}

void  RecoMCTruthLinker::keepMCParticle( MCParticle* mcp ){
  
  mcp->ext<MCPKeep>() = true  ;
	

  const MCParticleVec& parents = mcp->getParents() ;
    
  streamlog_out( DEBUG4 ) << " keepMCParticle keep particle with pdg : " << mcp->getPDG() 
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
  
  streamlog_out(DEBUG) << " check " << std::endl;
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


  LCCollection* mcpCol = NULL;
  try{
    mcpCol = evt->getCollection( _mcParticleCollectionName ) ;
  } catch (DataNotAvailableException e) {
    streamlog_out(WARNING) << "RecoMCTructh::Check(): MCParticle collection \"" << _mcParticleCollectionName << "\" does not exist, skipping" << std::endl;
  }

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
  LCCollection* mcpskCol = NULL;
  try{  
    mcpskCol = evt->getCollection( _mcParticlesSkimmedName ); 
  }
  catch (DataNotAvailableException e){
    streamlog_out(WARNING) << "RecoMCTructh::Check(): MCParticleSkimmed collection \"" << _mcParticlesSkimmedName << "\" does not exist, skipping" << std::endl;
  }

  etot = 0.0 ;
  nStable = 0 ;
  int nMCPSK = 0;

  if (mcpskCol){
    nMCPSK  = mcpskCol->getNumberOfElements() ;
    for(int i=0; i< nMCPSK ; i++){
      MCParticle* mcpsk = dynamic_cast<MCParticle*> ( mcpskCol->getElementAt( i ) ) ;
      if( mcpsk->getGeneratorStatus() == 1 ) {
	hmcpsk_e->fill(   mcpsk->getEnergy()  ) ;
	etot +=  mcpsk->getEnergy()  ;
	++nStable ;
      }
    }
  }//If Skimmed Collection Exists
  //Fill this, even if MCParticles Skimmed is empty, all are 0
  hmcpsk_n->fill( nStable ) ;
  
  hmcpsk_ntot->fill( nMCPSK ) ;
    
  hmcpsk_etot->fill( etot ) ;

#endif

}


void RecoMCTruthLinker::end(){ 
  
  streamlog_out(DEBUG4) << " processed " << _nEvt << " events in " << _nRun << " runs "
			<< std::endl ;

}







