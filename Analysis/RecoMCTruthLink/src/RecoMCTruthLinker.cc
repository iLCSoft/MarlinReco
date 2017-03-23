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

#include "UTIL/LCTOOLS.h"
#include <cstdio>

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

typedef std::map< Track* , int > TrackMap ;

typedef std::map< MCParticle* , double > MCPMapDouble ;

typedef std::map< MCParticle* , MCParticle* >  Remap_as_you_go;

RecoMCTruthLinker::RecoMCTruthLinker() : Processor("RecoMCTruthLinker") {
  
  // modify processor description
  _description = "links RecontructedParticles to the MCParticle based on number of hits used" ;
  
  
  IntVec pdgVecDef ;
  
  pdgVecDef.push_back(  22 ) ;  // gamma
  pdgVecDef.push_back( 111 ) ;  // pi0
  pdgVecDef.push_back( 310 ) ;  // K0s
  
  _encoder = new UTIL::BitField64(lcio::LCTrackerCellID::encoding_string());
  
 
   //  >>>>> reconstructed or simulated thingis 
 
  
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection" , 
			   "Name of the MCParticle input collection"  ,
			   _mcParticleCollectionName ,
			   std::string("MCParticle") ) ;
  
  registerInputCollection( LCIO::TRACK,
			   "TrackCollection" , 
			   "Name of the Tracks input collection"  ,
			   _trackCollectionName ,
			   std::string("MarlinTrkTracks") ) ;
  
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
  
 // <<<<<<
 
//  >>>>> tracker hits and sim<->seen relations 

 StringVec trackerHitsInputColNamesDefault ;
     trackerHitsInputColNamesDefault.push_back("VXDCollection ");
     trackerHitsInputColNamesDefault.push_back("SITCollection"); 
     trackerHitsInputColNamesDefault.push_back("FTD_PIXELCollection"); 
     trackerHitsInputColNamesDefault.push_back("FTD_STRIPCollection"); 
     trackerHitsInputColNamesDefault.push_back("TPCCollection"); 
     trackerHitsInputColNamesDefault.push_back("SETCollection");
  
  registerInputCollections( LCIO::SIMTRACKERHIT,
			   "SimTrackerHitCollections" ,
			   "Names of the SimTrackerHits input collection"  ,
			   _simTrkHitCollectionNames ,
			    trackerHitsInputColNamesDefault ) ;

  registerProcessorParameter("UseTrackerHitRelations",
			     "true: use relations for TrackerHits, false : use getRawHits ",
			     _use_tracker_hit_relations,
			     bool(true));
  
  
  StringVec trackerHitsRelInputColNamesDefault;
  trackerHitsRelInputColNamesDefault.push_back( "VXDTrackerHitRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "SITSpacePointRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "FTDPixelTrackerHitRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "FTDSpacePointRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "TPCTrackerHitRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "SETSpacePointRelations" );
  
  
  registerInputCollections("LCRelation",
                           "TrackerHitsRelInputCollections",
                           "Name of the lcrelation collections, that link the TrackerHits to their SimTrackerHits.",
                           _colNamesTrackerHitRelations,
                           trackerHitsRelInputColNamesDefault );

  // <<<<<<
    
  //  >>>>> Calorimeter hits and sim<->seen relations 

   StringVec caloHitsInputColNamesDefault ;
     caloHitsInputColNamesDefault.push_back("BeamCalCollection");
     caloHitsInputColNamesDefault.push_back("EcalBarrelSiliconCollection");
     caloHitsInputColNamesDefault.push_back("EcalBarrelSiliconPreShowerCollection");
     caloHitsInputColNamesDefault.push_back("EcalEndcapRingCollection");
     caloHitsInputColNamesDefault.push_back("EcalEndcapRingPreShowerCollection");
     caloHitsInputColNamesDefault.push_back("EcalEndcapSiliconCollection");
     caloHitsInputColNamesDefault.push_back("EcalEndcapSiliconPreShowerCollection");
     caloHitsInputColNamesDefault.push_back("HcalBarrelRegCollection");
     caloHitsInputColNamesDefault.push_back("HcalEndCapRingsCollection");
     caloHitsInputColNamesDefault.push_back("HcalEndCapsCollection");
     caloHitsInputColNamesDefault.push_back("LHcalCollection");
     caloHitsInputColNamesDefault.push_back("LumiCalCollection");
     caloHitsInputColNamesDefault.push_back("MuonBarrelCollection");
     caloHitsInputColNamesDefault.push_back("MuonEndCapCollection");

  registerInputCollections( LCIO::SIMCALORIMETERHIT,
			   "SimCaloHitCollections" ,
			   "Names of the SimCaloHits input collections"  ,
			   _simCaloHitCollectionNames ,
			    caloHitsInputColNamesDefault ) ;

  StringVec caloHitsRelInputColNamesDefault;
  caloHitsRelInputColNamesDefault.push_back( "RelationCaloHit" );
  caloHitsRelInputColNamesDefault.push_back( "RelationLcalHit" );
  caloHitsRelInputColNamesDefault.push_back( "RelationLHcalHit" );
  caloHitsRelInputColNamesDefault.push_back( "RelationBCalHit" ) ;
  caloHitsRelInputColNamesDefault.push_back( "RelationMuonHit" );

  registerInputCollections( "LCRelation",
                            "SimCalorimeterHitRelationNames" , 
			   "Name of the  lcrelation collections, that link the SimCalorimeterHit to CalorimeterHit"  ,
			   _caloHitRelationNames ,
			    caloHitsRelInputColNamesDefault );
  
  // <<<<<<

  //  >>>> Output true<->seen realations

  registerOutputCollection( LCIO::LCRELATION,
			    "MCTruthTrackLinkName" , 
			    "Name of the trackMCTruthLink output collection"  ,
			    _mCTruthTrackLinkName ,
			    std::string("") ) ;

  
  registerOutputCollection( LCIO::LCRELATION,
			    "TrackMCTruthLinkName" , 
			    "Name of the trackMCTruthLink output collection - not created if empty()"  ,
			    _trackMCTruthLinkName ,
			    std::string("") ) ;

  
  registerOutputCollection( LCIO::LCRELATION,
			    "ClusterMCTruthLinkName" , 
			    "Name of the clusterMCTruthLink output collection - not created if empty()"  ,
			    _clusterMCTruthLinkName ,
			    std::string("") ) ;


 registerOutputCollection( LCIO::LCRELATION,
			    "MCTruthClusterLinkName" , 
			    "Name of the MCTruthClusterLink output collection"  ,
			    _mCTruthClusterLinkName ,
			    std::string("") ) ;


 registerProcessorParameter("FullRecoRelation",
                             "true: All reco <-> true relations are given, with weight = 10000*calo weight+"
                             "track weight (weights in permill). false: Only highest contributor linked,"
                             "and only to tracks, not clusters if there are any tracks",   
                             _FullRecoRelation,
                             bool(false)
                             );
  
  registerOutputCollection( LCIO::LCRELATION,
			    "RecoMCTruthLinkName" ,
			    "Name of the RecoMCTruthLink output collection - not created if empty()"  ,
			    _recoMCTruthLinkName ,
			    std::string("RecoMCTruthLink") ) ;
  
  registerOutputCollection( LCIO::LCRELATION,
                            "MCTruthRecoLinkName" ,
                            "Name of the MCTruthRecoLink output collection"  ,
                            _mCTruthRecoLinkName ,
                            std::string("") ) ;

  
   registerOutputCollection( LCIO::LCRELATION,
			    "CalohitMCTruthLinkName" ,
			    "Name of the updated calo-hit MCTruthLink output collection - not created if empty()"  ,
			    _calohitMCTruthLinkName ,
			    std::string("") ) ;
  
   // <<<<<<<
   // Output skimmed MCParticles

  registerProcessorParameter(  "KeepDaughtersPDG" , 
			       "PDG codes of particles of which the daughters will be "
			       "kept in the skimmmed MCParticle collection"  ,
			       _pdgVec ,
			       pdgVecDef ) ;
 
  registerOutputCollection( LCIO::MCPARTICLE,
			    "MCParticlesSkimmedName" , 
			    "Name of the skimmed MCParticle  output collection - not created if empty()"  ,
			    _mcParticlesSkimmedName ,
			    std::string("") ) ;
  // <<<<<<<

  // various steering parameters 
  
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
  
  registerProcessorParameter( "UsingParticleGun" , 
			      "If Using Particle Gun Ignore Gen Stat"  ,
			      _using_particle_gun,
			      bool(false)  
			      ) ;
  
  
  
  registerProcessorParameter( "BremsstrahlungEnergyCut" , 
			      "energy cut for Brems that are kept"  ,
			      _bremsstrahlungEnergyCut,
			      float( 1. )  
			      ) ;

   registerProcessorParameter( "InvertedNonDestructiveInteractionLogic" ,
                              "Work-around Mokka bug in vertex-is-not-endpoint-of-parent flag (logic inverted)"  ,
                              _invertedNonDestructiveInteractionLogic,
                              bool(false)
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

  // don't write outpur collections that have an empty name
   _OutputClusterTruthRelation =   ! _clusterMCTruthLinkName.empty() ;  
   _OutputTruthClusterRelation =  !  _mCTruthClusterLinkName.empty() ;  
   _OutputCalohitRelation =   ! _calohitMCTruthLinkName.empty() ;
   _OutputTrackTruthRelation   =  ! _trackMCTruthLinkName.empty() ;
   _OutputTruthTrackRelation   =   ! _mCTruthTrackLinkName.empty() ;
   _OutputTruthRecoRelation    = ! _mCTruthRecoLinkName.empty()   ;

  
  if( ! _use_tracker_hit_relations ){

    streamlog_out( WARNING ) << "  ====== UseTrackerHitRelations=false => not using TrackerHit-SimTrackerHit-Relations  but getRawHits() instead - \n"
			     << "         this is probably not what you want (only for backward compatibility with very old files...)" << std::endl ;
  }
   
   _nRun = 0 ;
   _nEvt = 0 ;
}


void RecoMCTruthLinker::processRunHeader( LCRunHeader* run) { 
  
  _nRun++ ;
} 


void RecoMCTruthLinker::processEvent( LCEvent * evt ) { 
  
  _nEvt ++ ;
  
  
  streamlog_out( MESSAGE0 ) << " processEvent "  <<  evt->getEventNumber() << "  - " << evt->getRunNumber() 
  << std::endl ; 

  LCCollection* mcpCol = evt->getCollection( _mcParticleCollectionName ) ;
  
  
  // find track to MCParticle relations
  
  LCCollection* trackCol = 0 ;
  LCCollection* ttrlcol  = 0 ;
  LCCollection* trtlcol = 0 ;

  bool haveTracks = true  ;
  
  try{ trackCol = evt->getCollection( _trackCollectionName ) ;  }   catch(DataNotAvailableException&){  haveTracks=false ; }
  
  if( ! haveTracks ) {
    streamlog_out( DEBUG9 ) << " Track collection : " << _trackCollectionName 
    << " not found - cannot create relation " << std::endl ;
    
    
  } else {  //if ( haveTracks ) {
    
      trackLinker( evt,  mcpCol ,  trackCol  ,  &ttrlcol , &trtlcol);

      if (_OutputTrackTruthRelation ) 
         evt->addCollection(  ttrlcol  , _trackMCTruthLinkName  ) ;
      if (_OutputTruthTrackRelation ) 
         evt->addCollection(  trtlcol  , _mCTruthTrackLinkName  ) ;
  }
  
  // find cluster to MCParticle relations and the updated calohit to MCParticle relations.
  
  LCCollection* clusterCol = 0 ;
   
  bool haveClusters = true ;
  bool haveCaloHitRel = true ;
  
  LCCollection* ctrlcol = 0;
  LCCollection* trclcol = 0;
  LCCollection* chittrlcol = 0;
  
  try{ clusterCol = evt->getCollection( _clusterCollectionName ) ; }   catch(DataNotAvailableException&){  haveClusters=  false ; } 
  
  if( ! haveClusters ) {
    streamlog_out( DEBUG9 ) << " Cluster collection : " << _clusterCollectionName 
    << " not found - cannot create relation " << std::endl ;
  }
  
  
  if( haveTracks && haveClusters ) {
    
    clusterLinker( evt, mcpCol ,  clusterCol, 
		   &ctrlcol , &trclcol , &chittrlcol );
    
    if (_OutputClusterTruthRelation )   evt->addCollection(  ctrlcol  , _clusterMCTruthLinkName  ) ;
    if (_OutputTruthClusterRelation )   evt->addCollection(  trclcol  , _mCTruthClusterLinkName );
    if (_OutputCalohitRelation ) evt->addCollection(  chittrlcol  , _calohitMCTruthLinkName ) ;
  }
  
  
  
  
  // combine track and cluster to MCParticle relations to the reconstructed particle
  // to MCParticle relation
  LCCollection*  particleCol = 0 ;
  bool haveRecoParticles = true ; 
  
  try { particleCol = evt->getCollection(  _recoParticleCollectionName ); }   catch(DataNotAvailableException&){ haveRecoParticles = false ; } 
  
  if( ! haveRecoParticles ) {
    streamlog_out( DEBUG9 ) << " ReconstructedParticle collection : " << _recoParticleCollectionName 
    << " not found - cannot create relation " << std::endl ;
  }
  
  LCCollection* ptrlcol = 0;
  LCCollection* trplcol = 0;
  if( haveRecoParticles &&  haveTracks && haveClusters && haveCaloHitRel && !_recoMCTruthLinkName.empty() ) {
    
    particleLinker(   mcpCol , particleCol, ttrlcol,  ctrlcol, trtlcol,  trclcol, &ptrlcol, &trplcol);
    
    evt->addCollection(  ptrlcol  , _recoMCTruthLinkName  ) ;
    if (_OutputTruthRecoRelation )  evt->addCollection(  trplcol  , _mCTruthRecoLinkName  ) ;
    
  }
  
  if( haveTracks && haveClusters && haveCaloHitRel && !_mcParticlesSkimmedName.empty() ) {
    
    LCCollectionVec* skimVec = new LCCollectionVec( LCIO::MCPARTICLE )  ;
    
    makeSkim(    mcpCol , ttrlcol,  ctrlcol , &skimVec );
    evt->addCollection(   skimVec , _mcParticlesSkimmedName ) ;
    if( streamlog::out.write<streamlog::DEBUG5>() ){ linkPrinter (  skimVec , particleCol, ptrlcol,  trplcol ); }
  }

  //If either collection has not been added to the event, we have to delete it now!
  //Don't delete them before, because they are used
  if(!_OutputClusterTruthRelation) { delete ctrlcol; }
  if(!_OutputTruthClusterRelation) { delete trclcol; }
  if(!_OutputTrackTruthRelation) { delete ttrlcol; }
  if(!_OutputTruthTrackRelation) { delete trtlcol; }
  if(!_OutputCalohitRelation) { delete chittrlcol; }
  if(!_OutputTruthRecoRelation ) { delete trplcol ; }
  
}
void RecoMCTruthLinker::trackLinker(   LCEvent * evt, LCCollection* mcpCol ,  LCCollection* trackCol,  
                                       LCCollection** ttrlcol,  LCCollection** trtlcol) { 


  
  // merge all the SimTrackerHit - TrackerHit relations into on collection and set up the combined navigator
  mergeTrackerHitRelations(evt);
  

  LCRelationNavigator trackTruthRelNav(LCIO::TRACK , LCIO::MCPARTICLE  ) ;
   
  // the inverse relation from MCTruth particles to tracks 
  // weight is the realtive number of hits from a given MCParticle on the track
  LCRelationNavigator truthTrackRelNav(LCIO::MCPARTICLE , LCIO::TRACK  ) ;
 
  //========== fill a map with #SimTrackerHits per MCParticle ==================
  MCPMap simHitMap ;  //  counts total simhits for every MCParticle
  for( unsigned i=0,iN=_simTrkHitCollectionNames.size() ; i<iN ; ++i){
    
    const LCCollection* col = 0 ;
    try{ col = evt->getCollection( _simTrkHitCollectionNames[i] ) ; } catch(DataNotAvailableException&) {}
    if( col )
      for( int j=0, jN= col->getNumberOfElements() ; j<jN ; ++j ) {
        
        SimTrackerHit* simHit = (SimTrackerHit*) col->getElementAt( j ) ; 
        MCParticle* mcp = simHit->getMCParticle() ;
        simHitMap[ mcp ] ++ ;
    }
  }    
  //===========================================================================

  // loop over reconstructed tracks
  int nTrack = trackCol->getNumberOfElements() ;
  
  int ifoundch =0;
  
  streamlog_out( DEBUG6 ) << " *** Sorting out Track<->MCParticle using simHit<->MCParticle." << std::endl;

  for(int i=0;i<nTrack;++i){
    
    Track* trk = dynamic_cast<Track*> ( trackCol->getElementAt(i) ) ;

    // charged particle  or V0  -> analyse track. We need to find all seen hits
    // this track is made of, wich sim hits each of the seen hits came from,
    // and finally which true particle actually created each sim hit
    
    MCPMap mcpMap ;  // mcpMap is a map seen <-> true particle
    
    int nSimHit = 0 ;
    
    const TrackerHitVec& trkHits = trk->getTrackerHits() ;
    
    for( TrackerHitVec::const_iterator hitIt = trkHits.begin() ; hitIt != trkHits.end() ; ++hitIt ) { 
      
      TrackerHit* hit = * hitIt ; // ... and a seen hit ... 
      
      const LCObjectVec& simHits  = _use_tracker_hit_relations ? *(this->getSimHits(hit)) : hit->getRawHits() ;
      MCParticle* mcp2 = 0;       
      for( LCObjectVec::const_iterator objIt = simHits.begin() ; objIt != simHits.end() ; ++objIt ) {
        

        SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>( *objIt ) ; // ...and a sim hit ...
 
        MCParticle* mcp = simHit->getMCParticle() ; // ... and a true particle !
        if ( simHits.size() > 1 ) {
          if ( mcp2 != 0 && mcp2 != mcp ) { 
  	    //  In  this case, the mcp:s will count double !!
            streamlog_out( DEBUG3 ) << " ghost/double " << mcp << " " << mcp2 << " " << hit->getCellID0() <<std::endl;
          }
          mcp2=mcp;
        }       
        if ( mcp != 0 ) {
          mcpMap[ mcp ]++ ;   // count the hit caused by this true particle
        } else {
          streamlog_out( WARNING ) << " tracker SimHit without MCParticle ?!   " <<  std::endl ;
        }
        
        ++nSimHit ; // total hit count
        
      }
    }
    
    
    if( nSimHit == 0 ){
      
      streamlog_out( WARNING ) << " No simulated tracker hits found. Set UseTrackerHitRelations to true in steering file to enable using TrackerHit relations if they are available." <<  std::endl ;
      streamlog_out( WARNING ) << trk->id() <<" " << i << " " << trk->getTrackerHits().size()  <<  std::endl ; 
 
      for( TrackerHitVec::const_iterator hitIt = trkHits.begin() ; hitIt != trkHits.end() ; ++hitIt ) { 
      
        TrackerHit* hit = * hitIt ; // ... and a seen hit ... 
	streamlog_out( WARNING ) << hit->getPosition()[0]  <<  " " <<  hit->getPosition()[1]  <<  " " <<  hit->getPosition()[2]  <<  " " << std::endl ;       
        const LCObjectVec& simHits  = *(this->getSimHits(hit)) ;
	streamlog_out( WARNING ) <<  simHits.size()  <<  std::endl ; 
      }
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
      
      
      mother = ( it->first->getParents().size()!=0  ? dynamic_cast<MCParticle*>(it->first->getParents()[0])  : 0 )  ; // mother of the true particle.
      
      if ( _using_particle_gun || it->first->getGeneratorStatus() == 1 ) {  // genstat 1 particle, ie. it is a bona fide
                                                                            // creating true particle: enter it into the list,
                                                                            // and note how many hits it produced.
        theMCPs.push_back( it->first ) ;  
        MCPhits.push_back( it->second ) ; 
        ifound++;

      } else {  // not genstat 1. Wat should we do with it ?

        if ( mother != 0 ) { // if it has a parent, save it

          theMCPs.push_back(it->first);  
          MCPhits.push_back(it->second); 
          ifound++;
          
        } else {

          streamlog_out( WARNING ) << " track has hit(s) from a non-generator particle with no parents ?!! "  << std::endl;
        }           
      }
    } // end of loop over map
    
    // hits in track
    unsigned nHit = 0;
    
    for (unsigned ihit=0; ihit<trkHits.size(); ++ihit) {
      if ( UTIL::BitSet32( trkHits[ihit]->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]) {
        nHit += 2;
      } else {
        nHit += 1;
      }
    }
    
    // finally calculate the weight of each true particle to the seen 
    // (= hits_from_this_true/ all_hits),
    // and add the weighted reltion
    for (int iii=0 ; iii<ifound ; iii++ ) {
      

      float  weight = float(MCPhits[iii]  )/float(nHit) ; 
      
      trackTruthRelNav.addRelation(   trk , theMCPs[iii] , weight ) ;
      
      int Total_SimHits_forMCP = simHitMap[ theMCPs[iii] ];
      
      float inv_weight = float(MCPhits[iii]  ) / Total_SimHits_forMCP  ;
      
      truthTrackRelNav.addRelation(   theMCPs[iii] , trk , inv_weight ) ;

      
      streamlog_out( DEBUG4 ) << "    track " << trk->id() << " has " << MCPhits[iii]  << " hits of "
      << nSimHit << " SimHits ("
      << nHit <<    " TrackerHits) "
      << " weight = " << weight << " , "
      << " inv rel weight = " << inv_weight << "  "
      << " of MCPart with pdg : " << theMCPs[iii]->getPDG()
      << " , SimHits " <<  Total_SimHits_forMCP
      << " and genstat : " << theMCPs[iii]->getGeneratorStatus()
			      << " id: " << theMCPs[iii]->id()
      << std::endl ;

      
    
    }
    
    ifoundch=ifound;
  } 
  //  seen-true relation complete. add the collection

  streamlog_out( DEBUG6 ) << " *** Sorting out Track<->MCParticle : DONE " << std::endl;
  streamlog_out( DEBUG6 ) << " *** track linking complete, create collection " << std::endl;
  
  *trtlcol = truthTrackRelNav.createLCCollection() ;
  *ttrlcol = trackTruthRelNav.createLCCollection() ;
  
  delete _navMergedTrackerHitRel ; _navMergedTrackerHitRel = 0;
  delete _navMergedCaloHitRel ; _navMergedCaloHitRel = 0;
  
}



void RecoMCTruthLinker::clusterLinker(  LCEvent * evt,  LCCollection* mcpCol ,  LCCollection* clusterCol, 
                                         LCCollection** ctrlcol, LCCollection** trclcol,LCCollection** chittrlcol) { 

  
  
  
  
  mergeCaloHitRelations(evt);
  LCRelationNavigator clusterTruthRelNav(LCIO::CLUSTER , LCIO::MCPARTICLE  ) ;
  LCRelationNavigator truthClusterRelNav( LCIO::MCPARTICLE  , LCIO::CLUSTER ) ;
  LCRelationNavigator chitTruthRelNav(LCIO::CALORIMETERHIT , LCIO::MCPARTICLE  ) ;
   
  
  Remap_as_you_go remap_as_you_go ;  // map from true tracks linked to hits to 
                                     // those that really should have been linked
  MCPMapDouble simHitMapEnergy ;  //  counts total simhits for every MCParticle


  streamlog_out( DEBUG6 ) << " *** Sorting out simHit<->MCParticle connections, and find corresponding calo hits " << std::endl;

  // This step is needed for cluster (contrary to the tracks above), because the first few branchings in the
  // shower is kept in the MCParticle collection. We therefor need to back-track to the particle actually
  // entering the calorimeter. The originator of a hit sometimes IS a generator particle, in which case
  // there is no problem. Otherwise, the criteria to keep a particle as a hit originator is:
  //  The particle is an ancestor of the one initial linked to the hit.
  //      AND
  //  [
  //   The particle is a generator particle
  //        OR
  //   The particle starts in the tracker, ends in calorimeter. That it starts in the tracker
  //    is determined by the fact that it's mother ends there. However, note the posibility of 
  //    "non-destuctive interactions", where a particles mother does *not* end at the point where
  //    the particle is created !
  //        OR
  //   The particle, or one of it's ancestors is a back-scatterer. This case is further treated in the
  //    loop over clusters, to catch the rather common case that a back-scattered particle re-enters the
  //    calorimeter close to it's start-start point (think 15 MeV charged particle!) and the hits it
  //    produces is in the same cluster as that of the initiator of the shower the backscatter come from.
  //    (This can't be detected in this first loop, since we don't know about clusters here) 
  //   ]
  //   A map is set up, so obviously the first check is wether the MCP of the hit has already been
  //   treated when looking at some previous hit, in which case one just thake that association.


  for( unsigned i=0,iN=_simCaloHitCollectionNames.size() ; i<iN ; ++i){
    const LCCollection* col = 0 ;
    try{ col = evt->getCollection( _simCaloHitCollectionNames[i] ) ; } catch(DataNotAvailableException&) {}
    if( col ) {

     streamlog_out( DEBUG6 ) << std::endl;
     streamlog_out( DEBUG6 ) << " ================= " << std::endl;
     streamlog_out( DEBUG6 ) << " Treating sim-hits in " <<  _simCaloHitCollectionNames[i] << ". It has " 
			      <<  col->getNumberOfElements() << " hits " << std::endl;


      for( int j=0, jN= col->getNumberOfElements() ; j<jN ; ++j ) {
       
        SimCalorimeterHit* simHit = (SimCalorimeterHit*) col->getElementAt( j ) ; 
        LCObjectVec caloHits  = _navMergedCaloHitRel->getRelatedFromObjects(simHit);
        if (   caloHits.size() == 0 ) { continue ;}
        if (   caloHits.size() != 1 ) { streamlog_out( WARNING ) << " Sim hit with nore than one calo hit ? " << std::endl; }
        CalorimeterHit* caloHit = dynamic_cast<CalorimeterHit*>(caloHits[0]);
        double calib_factor = caloHit->getEnergy()/simHit->getEnergy();


	streamlog_out( DEBUG4 ) << std::endl;
	streamlog_out( DEBUG4 ) << "  ================= " << std::endl;
	streamlog_out( DEBUG4 ) << "    Treating hit " << j << " sim hit id " << simHit->id() << " calo hit id " 
                                << caloHit->id() << " nb contributions " << simHit->getNMCContributions() << std::endl;
	streamlog_out( DEBUG3 ) << "      ncalo hits : " << caloHits.size() << " calib factor " << calib_factor 
                                << " position " << simHit->getPosition()[0] << " " << 
                                                               simHit->getPosition()[1] << " " << 
                                                               simHit->getPosition()[2] << " " << std::endl;

        for(int k=0;k<simHit->getNMCContributions() ;k++){
          MCParticle* mcp = simHit->getParticleCont( k ) ;

	  streamlog_out( DEBUG2 ) << " a calo hit to treat " << std::endl;

          if ( mcp == 0 ) {


	    streamlog_out( DEBUG6 ) <<" collection nb "  << i <<  " nhits " << jN << " hit " << j << " contr " << k << std::endl;
	    streamlog_out( DEBUG6 ) <<" Sim hit with no mcp. N contrib is "  <<
                  simHit->getNMCContributions() <<std::endl;
            streamlog_out( DEBUG6 ) <<"     Hit Position: " << simHit->getPosition()[0] << " " << 
                                                               simHit->getPosition()[1] << " " << 
                                                               simHit->getPosition()[2] << " " << std::endl;
            if ( simHit->getNMCContributions() > 0 ) {
              for (int l=0;l<simHit->getNMCContributions() ;l++){
                streamlog_out( DEBUG6 ) <<" Contrib " << l << " : " <<std::endl;
                streamlog_out( DEBUG6 ) <<"     Energy: " << simHit->getEnergyCont(l) <<std::endl;
                streamlog_out( DEBUG6 ) <<"     PDG:    " << simHit->getPDGCont(l) <<std::endl;
                streamlog_out( DEBUG6 ) <<"     MCPart: " << simHit->getParticleCont(l) <<std::endl;
              }
            }



            continue;
	  }

          double e  = simHit->getEnergyCont( k ) * calib_factor ;

          streamlog_out( DEBUG3 ) <<"      initial true contributor "<< k << " id " << mcp->id() 
                                  <<" (with E = " << mcp->getEnergy() << " and pdg " << mcp->getPDG() << " ) e hit: " << e <<std::endl;

          if ( remap_as_you_go.find(mcp) != remap_as_you_go.end() ) {
            // very first condition: I already know what to do from some earlier hit created by mcp.
            mcp=remap_as_you_go.find(mcp)->second;
            streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / " 
                                << caloHit->id() << " attributed to " << mcp->id() << 
                        " because it's origin mcp has already been treated : originator case 0 " <<std::endl;
          } else {
            MCParticle* mother = 0;
            MCParticle* this_Kid = mcp ;     // ... and a true particle !
            
            //Particle gun particles dont have parents, but are "created in the simulation" and have genStat 0
            if ( mcp-> getGeneratorStatus() == 0 && ( mcp->getParents().empty() == false ) ) { 

              // not from generator, find which true particle this
              // hit should really be attributed to, by tracking back the history.
              // (the other case, ie. if the if above is false is "case 1")
              
              // Two cases to treat:
              //   For some reason, the hit is not attributed to the
              //   incomming particle, but to some particle in the shower.
              //   Just track back to the incomming particle, which might (usually)
              //   be a gen-stat 1 particle from the main vertex. It can also
              //   be from a decay or interaction in the tracking made by
              //   Geant : genstat=0, mother isDecayedInTracker, (or 
              //   the production vertex is in the track, but the mother
              //   continues on - ticky case, see comment futher down), or a
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
              

              
              streamlog_out( DEBUG2 ) << "        simHit " <<simHit->id()<<","<< j << 
              " not created by generator particle. backtracking ..." 
              << std::endl;
              streamlog_out( DEBUG2 ) <<"          "<<mcp->id()<<" gs "<<mcp->getGeneratorStatus()<<
              " dint "<<mcp->isDecayedInTracker()<<
              " bs "<<mcp->isBackscatter()<<
              " ndi "<<mcp->vertexIsNotEndpointOfParent()<<
              " npar "<<mcp->getParents().size()<<
              " pdg "<<mcp->getPDG()<<
              " "<< mcp->getVertex()[0]<<
              " "<<mcp->getVertex()[1]<<
              " "<<mcp->getVertex()[2]<<
               " "<< mcp->getEndpoint()[0]<<
              " "<<mcp->getEndpoint()[1]<<
              " "<<mcp->getEndpoint()[2]<<std::endl;
             


              mother= dynamic_cast<MCParticle*>(mcp->getParents()[0]); 


              if ( !this_Kid->isBackscatter() &&  mother!= 0 &&
                     mother->getParents().size()>0 && 
                     mother->getGeneratorStatus() ==0 &&
                     !mother->isDecayedInTracker() ) { 
              
	        streamlog_out( DEBUG2 ) << "        goes into originator loop " << std::endl;
              }
	      while ( !this_Kid->isBackscatter() && mother!= 0 &&   mother->getParents().size()>0 && 
                       mother->getGeneratorStatus() ==0 &&
                       !mother->isDecayedInTracker() ) { 
                // back-track as long as there is a non-generator 
                // mother, or the mother decayed in the tracker 
                // (=> the kid is the particle entering the calorimeter.)
                
                // case shower-particle
 	        streamlog_out( DEBUG1 ) <<"          in originator loop " << std::endl;

		if ( this_Kid->vertexIsNotEndpointOfParent() != _invertedNonDestructiveInteractionLogic) {
		  MCParticle* oma=dynamic_cast<MCParticle*>(mother->getParents()[0]);
                  if ( oma->isDecayedInTracker() ) {
                    streamlog_out( DEBUG1 ) <<"          break out : gandmother "<<oma->id()<<
                                              " gs "<<oma->getGeneratorStatus()<<
                                              " dint "<<oma->isDecayedInTracker()<<
                                              " bs "<<oma->isBackscatter()<<
                                              " ndi "<<oma->vertexIsNotEndpointOfParent()<<
                                              " npar "<<oma->getParents().size()<<
                                              " pdg "<<oma->getPDG()<<std::endl;
                    break ;
                  }
		} 
                this_Kid=mother ;
                mother= dynamic_cast<MCParticle*>(mother->getParents()[0]); // (assume only one...)
                
                streamlog_out( DEBUG1 ) <<"          shower-part mother "<<mother->id()<<
                " gs "<<mother->getGeneratorStatus()<<
                " dint "<<mother->isDecayedInTracker()<<
                " bs "<<mother->isBackscatter()<<
                " ndi "<<mother->vertexIsNotEndpointOfParent()<<
                " npar "<<mother->getParents().size()<<
                " pdg "<<mother->getPDG()<<std::endl;

                
                
              }

              // Further treatment (basically determining if it this_Kid or mother that enetered the
              // calorimeter) based on why we left the while-loop 

              // here one of the while conditions is false, ie. at least one of
              // " kid is back-scatter", "no mother" , "mother has no parents", "mother is from 
              // generator", or "mother did decay in tracker" must be true. We know that this_Kid 
              // fulfills all the while-conditions except the first: obvious if at least one iteration of
              // the loop was done, since it was the mother in the previous iteration, 
              // but also true even if the loop wasn't transversed, due to the conditions
              // to at all enter this block of code (explicitly must be simulator particle with
              // mother, implicitly must have ended in the calo, since it did make calo hits.)
	      if (this_Kid->isBackscatter() ) {
                // case 2: Kid is back-scatterer. It has thus started in a calo, and entered
                //  from there into the tracking volume, and did cause hits after leaving the tracker
                //  volume again ->  this_Kid started before the calo, and is the one
                remap_as_you_go[mcp]=this_Kid;
                mcp=this_Kid; 
                streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / " 
                                << caloHit->id() << " attributed to kid " << mcp->id() << 
                        " because it's origin is a back-scatter : originator case 2 " <<std::endl;
                streamlog_out( DEBUG2 ) <<"          "<<this_Kid->id()<<
                " gs "<<this_Kid->getGeneratorStatus()<<
                " dint "<<this_Kid->isDecayedInTracker()<<
                " bs "<<this_Kid->isBackscatter()<<
                " npar "<<this_Kid->getParents().size()<<
                " pdg "<<this_Kid->getPDG()<<std::endl;
               

              } else if ( mother->isDecayedInTracker() ) { // the clear-cut case:
                remap_as_you_go[mcp]=this_Kid;
                mcp=this_Kid; // this_Kid started before the calo, and is the one 
                              // the hit should be attributed to
                

                streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / " 
                                << caloHit->id() << " attributed to kid " << mcp->id() << 
                        " because it's origin is in tracker : originator case 3 " <<std::endl;
                

              } else { // the other three cases, ie. one or sveral of "no mother", "no grand-mother", 
                       //  "generator particle" + that we know that "mother decayed in calo"
                if ( mother == 0 ) { // ... which of course implies no grand-mother, and no gen stat 
                                     // of the mother as well -> should not be possible !

                  streamlog_out( WARNING ) << "  MCparticle " << this_Kid->id() << 
                   " is a simulation particle, created in the calorimeter by nothing . "<< std::endl;


                  remap_as_you_go[mcp]=this_Kid;
                  mcp=this_Kid; // can't do better than that.

                } else { // here we know: "mother exists, but decayed in calo". In addition, two posibilities: 
                         // mother is generator particle, or there was no grand-parents. One or both 
                         // must be true here. Here it gets complicated, because what we want to know is
                         // whether this_Kid started in the tracker or not. Unluckily, we don't know that 
                         // directly, we only know where the mother ended. IF ithe mother ended in the tracker, 
                         // there is no problem, and has already been treated, but if it ended in the calo, it is 
                         // still possible that this_Kid came from a "non-destructive interaction" with the tracke-detector
                         // material. This we now try to figure out.

		  if (   this_Kid->vertexIsNotEndpointOfParent() == _invertedNonDestructiveInteractionLogic ) {
		      // This bizare condition is due to a bug in LCIO (at least for the DBD samples. this_Kid->vertexIsNotEndpointOfParent()
                      // should be true in the "non-destructive interaction", but it isn't: actually it is "false", but is "true" for the
                      // for particles that *do* originate at the end-vertex of their parent. This is a bug in Mokka.
                      // Hence the above ensures that there was NO "non-destructive interaction", and it is clear this_Kid was created at the end-point
                      // of the mother. The mother is either a generator particle (to be saved), or the "Eve" of the decay-chain (or both).
                      // So mother is the one to save and assign the hits to:
                    remap_as_you_go[mcp]=mother;
                    mcp=mother;   // mother started at ip, and reached the calo, and is the one 
                                  // the hit should be attributed to.

                    streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / " 
                                << caloHit->id() << " attributed to mother " << mcp->id() << 
                        " because it is a generator particle or started in tracker : originator case 4 " <<std::endl;

                  } else { // here we DO have a "non-destructive interaction". Unluckily, we cant directly know if this took place in the
                           // tracker (in which case we should keep this_Kid as the mcp to save and assign hits to), or not.
                           // We will play a few clean tricks to find the cases where it either certain that the 
                           // "non-destructive interation" was in the tracker, or that it was in the calo. This reduces 
                           // the number of uncertain cases to play dirty tricks with.

                           // Clean tricks to play: look at the sisters of this_Kid: with some luck one of them is a promptly decaying 
                           // particle, eg. a pi0.
                           // This sister will be flagged as decayed in calo/tracker, and from that we know for certain that the 
                           // "non-destructive interaction" was in the calo/tracker. 
                           // It can also be that one of the sisters is flagged as a back-scatter, which only happens in the calo. 
                           // If the particles grand-mother is decayed in calo, and the mother isn't from a "non-destructive interaction",
                           // the "non-destructive interaction" was in the calo.

                           // Finally, we can check the distance of the end-point of the mother (sure to be in the calo) to 
                           // the vertex of this_Kid. If this is small, this_Kid *probably* started in the calo.
                    

                    int starts_in_tracker = 0 ;
                    unsigned lll=0 ;
                    int has_pi0 = 0 ;
                    int oma_in_calo = 0;
                    double rdist =0.;
		    int has_bs = 0;

		    for ( unsigned kkk=0 ; kkk < mother->getDaughters().size() ; kkk++ ) {
                      MCParticle* sister = dynamic_cast<MCParticle*>(mother->getDaughters()[kkk]);
                      if ( sister == this_Kid ) continue;
 		      if (  abs(sister->getVertex()[0]-this_Kid->getVertex()[0]) > 0.1 ||
                            abs(sister->getVertex()[1]-this_Kid->getVertex()[1]) > 0.1 ||
                            abs(sister->getVertex()[2]-this_Kid->getVertex()[2]) > 0.1 ) continue;  // must check that it is the same vertex:
                                                                                                    // several "non-destructive interactions" can
                                                                                                    // take place (think delta-rays !)
 		      if ( sister->isBackscatter()) {
                        has_bs = 1 ; 
                        lll=kkk;
                        break ;
		      } else if ( sister->isDecayedInTracker() ) {
                        starts_in_tracker = 1 ;
                        lll=kkk;
                        break ;
                      }
                      // any pi0:s at all ? (it doesn't matter that we break at the two cases above, 
                      // because if we do, it doesn't matter if there are
                      // pi0 sisters or not !)
                      if ( sister->getPDG() == 111 ) {
                         has_pi0 = 1 ;
                      }
                    }
                    // if not already clear-cut, calculate distance vertext to mother end-point
                    if ( starts_in_tracker != 1 && has_bs != 1 && has_pi0 != 1 && oma_in_calo != 1 ) {
                      rdist=sqrt(pow(mother->getEndpoint()[0]-this_Kid->getVertex()[0],2)+
                                   pow(mother->getEndpoint()[1]-this_Kid->getVertex()[1],2)+
                                   pow(mother->getEndpoint()[2]-this_Kid->getVertex()[2],2));
                      if ( mother->getParents().size() != 0 ) {
		        MCParticle* oma=dynamic_cast<MCParticle*>(mother->getParents()[0]);
                        if ( oma->isDecayedInCalorimeter() ) {
                          oma_in_calo = 1 ;
			  streamlog_out( DEBUG1 ) << "          grandmother in calo " << std::endl;
                        }
                      }
                    }
                    streamlog_out( DEBUG1 ) << "          " <<  starts_in_tracker << " " << has_pi0 << " " <<  has_bs << " " << oma_in_calo << std::endl;
		    //                    }
                    if ( starts_in_tracker == 1 ) { // this_Kid is a clear-cut hit-originator

                      remap_as_you_go[mcp]=this_Kid;
                      mcp=this_Kid; // this_Kid started before the calo, and is the one 
                                    // the hit should be attributed to
                
                      streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / " 
                                << caloHit->id() << " attributed to kid " << mcp->id() << 
                        " because it's origin could be deduced to be in tracker : originator case 5 " <<std::endl;

		      streamlog_out( DEBUG2 ) << "        Details of case 5: "
                                             << this_Kid->getVertex()[0] << " " 
					     << this_Kid->getVertex()[1] << " " << this_Kid->getVertex()[2] << " " 
                                             << mother->getEndpoint()[0] << " " 
					     << mother->getEndpoint()[1] << " " << mother->getEndpoint()[2] << " " 
                                             << mother->getGeneratorStatus() <<  " " 
                                             << this_Kid->vertexIsNotEndpointOfParent() << std::endl;
                      streamlog_out( DEBUG1 ) << " starts in tracker " <<  std::endl;

		    } else if ( has_pi0 != 0 || has_bs != 0 || oma_in_calo != 0 ) { // clear-cut case of this_Kid starting in the calo. 
                                                                // We do know that the mother 
                                                                // is a generator particle and/or the "Eve" of the cascade, 
                                                                // so we should attribute hits to the
                                                                // mother and save it.

                      remap_as_you_go[mcp]=mother;
                      mcp=mother;   // mother started at ip, and reached the calo, and is the one 
                                    // the hit should be attributed to.

                      streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / " 
                                << caloHit->id() << " attributed to mother " << mcp->id() << 
                                 " because it's origin could be deduced to be in tracker : originator case 6 " <<std::endl;
 		      streamlog_out( DEBUG2 ) << "        Case 6 details: kid starts in calo " << " " 
					      << this_Kid->getVertex()[1] << " " << this_Kid->getVertex()[2] << " " << rdist <<  " " 
					      << has_pi0  << " " << has_bs << " " <<  oma_in_calo << std::endl;


		    } else { // un-clear case: no pi0 nor back-scatteres among the sisters to help to decide.
                             // Use distance this_Kid-startpoint to mother-endpoint. We know that the latter is
                             // in the calo, so if this is small, guess that the start point of this_Kid is
                             // also in the calo. Calos are dense, so typically in the case the "non-destructive interaction"
                             // is in the calo, one would guess  that the distance is small, ie. large distance ->
                             // unlikely that it was in the calo.

		      if ( rdist > 200. ) { // guess "non-destructive interaction" not in calo -> this_Kid is originator 
                       
                        remap_as_you_go[mcp]=this_Kid;
                        mcp=this_Kid; // this_Kid started before the calo, and is the one 
                                    // the hit should be attributed to
                        streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / " 
                                << caloHit->id() << " attributed to kid " << mcp->id() << 
                                 " because it's origin is guessed to be in tracker : originator case 7 " <<std::endl;
 		        streamlog_out( DEBUG2 ) << "        Case 7 details: guess kid starts in tracker " << " " 
						  << this_Kid->getVertex()[1] << " " << this_Kid->getVertex()[2] << " " << rdist <<  " " 
                                                  << has_pi0  << " " << std::endl;

                      } else { // guess in calo -> mother is originator
                        remap_as_you_go[mcp]=mother;
                        mcp=mother;   // mother started at ip, and reached the calo, and is the one 
                                    // the hit should be attributed to.
                        streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / " 
                                << caloHit->id() << " attributed to mother " << mcp->id() << 
                                 " because it's origin is guessed in tracker : originator case 8 " <<std::endl;
		        streamlog_out( DEBUG2 ) << "        Case 8 details:  guess kid starts in calo " << " " 
						  << this_Kid->getVertex()[1] << " " << this_Kid->getVertex()[2] << " " << rdist <<  " " 
                                                  << has_pi0  << " " << std::endl;
                        }

                    }
                  }          
                }
              }
              
            } else {
              // first case: the hit generating mcp itself already fulfills the firest
              // criterium, ie. it is a generator particle
              streamlog_out( DEBUG3 ) << "     Hit " << simHit->id() << " / " 
                                << caloHit->id() << " attributed to " << mcp->id() << 
                        " because it's origin is a generator particle : originator case 1 " <<std::endl;

              remap_as_you_go[mcp]=mcp;
	    } // genstat if - then - else
          }

          streamlog_out( DEBUG4 ) <<"    Final assignment for contribution " << k << " to " << simHit->id() << " / " 
				  << caloHit->id() << " : " << mcp->id() << std::endl;
          streamlog_out( DEBUG4 ) <<"      gs "<<mcp->getGeneratorStatus()<<
              " dint "<<mcp->isDecayedInTracker()<<
              " bs "<<mcp->isBackscatter()<<
              " npar "<<mcp->getParents().size()<<
              " pdg "<<mcp->getPDG()<<
              " "<<mcp->getEndpoint()[0]<<
              " "<<mcp->getEndpoint()[1]<<
              " "<<mcp->getEndpoint()[2]<<std::endl;
 
          // decided which mcp this sim-hit should be associated to :
	  //	already_known:          

          simHitMapEnergy[ mcp ] += e ;
          chitTruthRelNav.addRelation(  caloHit , mcp , e ) ;

        } // mc-contributon-to-simHit loop
      }
    }
  } // sim-calo-hits loop   

  streamlog_out( DEBUG6 ) << " *** Sorting out simHit<->MCParticle connections : DONE " << std::endl;
  streamlog_out( DEBUG6 ) << " *** Sorting out Cluster<->MCParticle using simHit<->MCParticle, re-assigning the latter in some rare cases." << std::endl;

  
  // loop over reconstructed particles
  int nCluster = clusterCol->getNumberOfElements() ;

  streamlog_out( DEBUG6 ) << std::endl;
  streamlog_out( DEBUG6 ) << " ================= " << std::endl;
  streamlog_out( DEBUG6 ) <<" Treating clusters. There are "<< nCluster <<" of them " << std::endl;
  
  std::vector<Cluster*> missingMC ;
  missingMC.reserve( nCluster ) ;
  int ifoundclu =0;
  // now for the clusters
  
  
  MCParticle* mother = 0;
	
  MCPMapDouble mcpEnergyTot ;
	

  for(int i=0;i<nCluster;i++){
    
    MCPMapDouble mcpEnergy ;
    double eTot = 0 ;
    Cluster* clu = dynamic_cast<Cluster*> ( clusterCol->getElementAt(i) ) ;

    
    
    
    // We need to find all seen hits this clutser is made of, which sim hits each 
    // of the seen hits came from, and finally which true particles actually created 
    // each sim hit. Contrary to the sim tracker hits above, a sim-calo hit can be
    // made by several true particles. They also have a signal size (energy) value.
    // In addition, the true particle creating sometimes needs to be back-tracked
    // to the particle actually entering tha calorimeter. 
    
    
    
    
    const CalorimeterHitVec& cluHits = clu->getCalorimeterHits() ;

    streamlog_out( DEBUG5 ) << std::endl;
    streamlog_out( DEBUG5 ) << "  ================= " << std::endl;
    streamlog_out( DEBUG5 ) <<"Cluster clu = "<< clu->id() << " (i = " << i << " ) with " << cluHits.size() << " hits " << std::endl;

    double ecalohitsum=0.;        
    double ecalohitsum_unknown=0.;	  
    int no_sim_hit = 0;
    for( CalorimeterHitVec::const_iterator hitIt = cluHits.begin() ;
       hitIt != cluHits.end() ; ++hitIt ) { 
      
      
      CalorimeterHit* hit = * hitIt ;  // ... a calo seen hit ...
      ecalohitsum+= hit->getEnergy();       
      
      
      const LCObjectVec& simHits  = *(this->getCaloHits(hit)) ;

      streamlog_out( DEBUG4 ) << std::endl;
      streamlog_out( DEBUG4 ) <<"     Treating hit = "<< hit->id() << " e " << hit->getEnergy()<< " nb sim hits : " <<
	simHits.size() << std::endl;

      double ehit = 0.0; 
      int nsimhit = 0;        
      for( LCObjectVec::const_iterator objIt = simHits.begin() ;
          objIt != simHits.end() ; ++objIt ){
        
        SimCalorimeterHit* simHit = dynamic_cast<SimCalorimeterHit*>( *objIt ) ; // ... and a sim hit ....
        
        double calib_factor = hit->getEnergy()/simHit->getEnergy();
        streamlog_out( DEBUG3 ) <<"     simhit = "<< simHit->id() << " has " << simHit->getNMCContributions() 
                                << " contributors " << std::endl;
        nsimhit++;
        for(int j=0;j<simHit->getNMCContributions() ;j++){
          
          MCParticle* mcp = simHit->getParticleCont( j ) ;
          double e  = simHit->getEnergyCont( j ) * calib_factor;
          streamlog_out( DEBUG3 ) <<"     true contributor = "<< mcp->id() << " e: " << e 
                                  <<" mapped to " <<remap_as_you_go.find(mcp)->second->id() << std::endl;
          if ( mcp == 0 ) {
            streamlog_out( DEBUG7 ) <<"      simhit = "<< simHit << " has no creator " <<std::endl;
            // streamlog_out( DEBUG7 ) <<"         true contributor = "<< mcp << " e: " << e 
            //    <<" mapped to " <<remap_as_you_go.find(mcp)->second << std::endl;
            continue ; 
          }
          mcp=remap_as_you_go.find(mcp)->second;
          mcpEnergy[ mcp ] +=  e ;// count the hit-energy caused by this true particle
          eTot += e ;             // total energy
          ehit+= e;
        } // mc-contributon-to-simHit loop
      } // simHit loop 

      streamlog_out( DEBUG4 )<< "     summed contributed e: " << ehit << " ratio : " << ehit/hit->getEnergy()
                             << " nsimhit " << nsimhit <<std::endl;
      if ( nsimhit == 0 ) {
        
        streamlog_out( DEBUG5 ) << " Warning: no simhits for calohit " << hit << 
             ". Will have to guess true particle ... " << std::endl;
        no_sim_hit = 1;
        ecalohitsum_unknown+= hit->getEnergy();	    
	streamlog_out( DEBUG5 )<< " sim-less calohit E and position " << hit->getEnergy() << " "
             << hit->getPosition()[0] << " " <<  hit->getPosition()[1] << " " <<  hit->getPosition()[2] << " " << std::endl;	    
	streamlog_out( DEBUG5 )<<  " sim-less calohit clust E and position " << clu->getEnergy() 
             << " "<< clu->getPosition()[0] << " " <<  clu->getPosition()[1] << " " <<  clu->getPosition()[2] << " " << std::endl;	    
      }
    } // hit loop
    if ( no_sim_hit == 1 ) {
      streamlog_out( DEBUG6 ) << "   Warning, there are  sim-less calohits in cluster " << clu->id() << std::endl;
    } 
    streamlog_out( DEBUG5 ) << std::endl;
    streamlog_out( DEBUG5 ) << "    Sum of calohit E: " <<  ecalohitsum << " cluster E " 
                            << clu->getEnergy() << " Sum of Simcalohit E: " << eTot 
                            << "  no_sim_hit: " <<  no_sim_hit << ". Energy from unknow source : " 
                            << ecalohitsum_unknown << std::endl;
    if( eTot == 0.0 ){ // fixme - this might happen if clusters are from Lcal/Muon only
      
      
      // save reco particle in missingMC 
      missingMC.push_back( clu  ) ;
      
      streamlog_out( DEBUG8 ) << " no calorimeter hits found for " 
			      << " cluster " << clu->id() << " e:"  << clu->getEnergy()  
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
      if ( it->first == 0 ) {   // ( if == 0, this cluster contains some (but not all) sim-hits with unknown origin.
                                //   If *all* sim-hits would have had unknown origin, we would already have "continue":ed above)
	streamlog_out( MESSAGE ) << " SimHit with unknown origin in cluster " << clu << " ( " << clu->id() << " ) " << std::endl;
	continue;
      }  
                                            
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
          } else { // else: if the mother is a BS, add to the  moreMCPs-list, further treated below,
                   //  otherwise keep as a bona fide creator.
            streamlog_out( DEBUG2 ) << " case 1 for "<< it->first->id() << 
            "(morefound=" << morefound << ")" << 
	      " mother: " << mother->id() <<
            " gs "  <<mother->getGeneratorStatus()<< 
            " dint " << mother->isDecayedInTracker() <<
            " bs "  << mother->isBackscatter() << 
            " pdg " <<mother->getPDG() <<std::endl;
            if (  mother->isBackscatter() == 1 ) {
              moreMCPs.push_back(it->first);  moreMCPes.push_back(it->second); morefound++;
            } else {
              theMCPs.push_back(it->first);  MCPes.push_back(it->second); ifound++;
            }
          }
        } else { // not genstat 1, no mother ?! Also add to the  moreMCPs-list.
          streamlog_out( DEBUG6 ) << " case 2 for "<< it->first->id() << 
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
	  " "<< moreMCPs[iii]->id()<<
        " "<<moreMCPs[iii]->getGeneratorStatus() << std::endl;
      }
      for (int iii=0 ; iii<ifound ; iii++ ) {
        streamlog_out( DEBUG2 ) << " iii, MCPes[iii], theMCPs[iii] " << iii <<
	  " "<< MCPes[iii] <<" "<< theMCPs[iii]->id() << std::endl;
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
        streamlog_out( DEBUG2 ) << "   iii: " << iii << " mother: " << mother->id()  <<std::endl; 
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
        
        streamlog_out( DEBUG2 ) << "       mother "<< mother->id() << 
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
        
        
        // find out if this ancestor is a back-scatterer actually directly giving rise to hits in this cluster. 
        // If so, we attribute the hits of  moreMCPs[iii] to that true particle instead.
        
        
        for (int kkk=0 ; kkk<ifound ; kkk++){
          MCParticle* tmcp = theMCPs[kkk];
          if ( tmcp == mother ) {
            // find hits related to moreMCPs[iii]
            streamlog_out( DEBUG3 ) << "        found " << moreMCPs[iii] << 
              "(iii= "<<iii <<")" << kkk << 
              " to be related to "<<mother->id() <<
              " add e : " <<  moreMCPes[iii] << std::endl;
            LCObjectVec hitvec = chitTruthRelNav.getRelatedFromObjects(moreMCPs[iii]);
            FloatVec evec=chitTruthRelNav.getRelatedFromWeights(moreMCPs[iii]);
	    int one_mod_done = 0;
            for ( unsigned lll=0 ; lll<hitvec.size() ; lll++ ) {
                 CalorimeterHit* hit  = dynamic_cast<CalorimeterHit*>(hitvec[lll]);
                 for( CalorimeterHitVec::const_iterator hitIt = cluHits.begin() ;
                            hitIt != cluHits.end() ; ++hitIt ) { 
                   if (  *hitIt == hit ) { // ... a calo seen hit ...
                     chitTruthRelNav.removeRelation(  hit ,moreMCPs[iii] );  
                     chitTruthRelNav.addRelation(  hit , mother , evec[lll] ) ;
                     one_mod_done = 1;
                     break ; 
                   }
                 }
             
            }
	    if ( one_mod_done == 1 ) {
              MCPes[kkk]+= moreMCPes[iii]; 
              goto endwhile;
            } else {
              streamlog_out( DEBUG3 ) << "        However, no hits related to the current cluster found ??" << std::endl;
            } 
            
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
        
        streamlog_out( DEBUG3 ) << "        No relation found for "<< moreMCPs[iii]->id() << 
        ". Keep it as separate originator "<< std::endl;
        theMCPs.push_back(moreMCPs[iii]);  MCPes.push_back(moreMCPes[iii]); ifound++;
      }
      
    }
    streamlog_out( DEBUG5 ) << " cluster " << clu->id() << " , E = " << clu->getEnergy() << " , ifound = " << ifound <<  std::endl;
    
    // finally calculate the weight of each true partic to the seen 
    // (= energy_from_this_true/ total ), and add the weighted reltion.
    
    float totwgt =0.0;
    for (int iii=0 ; iii<ifound ; iii++ ) {

      float  weight = (MCPes[iii]/eTot)*(clu->getEnergy()-ecalohitsum_unknown)/clu->getEnergy();
      mcpEnergyTot[ theMCPs[iii] ] +=   weight*clu->getEnergy();
      
      totwgt+= weight;
      if( theMCPs[iii] == 0 ) {
        
        streamlog_out( ERROR ) << " cluster " <<clu->id() << " has " << MCPes[iii] << " GeV of " << eTot 
        << " GeV [ " << weight << " ] true energy " 
        << " but no MCParticle " << std::endl ;
        
        continue ;
        
      }
      
      streamlog_out( DEBUG5 ) << " cluster " <<clu->id() << " has " <<  MCPes[iii] << " of " << eTot 
      << "  [ " << weight << " ] " 
      << " of MCParticle " << theMCPs[iii]->id() << " with pdg : " << theMCPs[iii]->getPDG() 
      << " and genstat : " << theMCPs[iii]->getGeneratorStatus() 
      << std::endl ;
      
      //Particle gun particles dont have parents but have genStat 0      
      if (  theMCPs[iii]->getGeneratorStatus() == 0 &&  ( theMCPs[iii]->getParents().empty() == false) ) {
        streamlog_out( DEBUG5 ) << " mother id " << theMCPs[iii]->getParents()[0]->id() 
        << " and genstat " 
        << theMCPs[iii]->getParents()[0]->getGeneratorStatus() << std::endl ;
        
      }
      

      
      clusterTruthRelNav.addRelation(   clu , theMCPs[iii] , weight ) ;
      // The reverse map. Differs in what the weight means: here it means:
      // "this cluster got wgt of all seen cluster energy the true produced"
      // (in the  other one it means:
      // "this true contributed wgt to the total seen energy of the cluster")
      weight=(MCPes[iii]/simHitMapEnergy[theMCPs[iii]]);
      truthClusterRelNav.addRelation(   theMCPs[iii] , clu , weight ) ;
      
    }
    ifoundclu=ifound;
  } // cluster loop
  
  streamlog_out( DEBUG6 ) << " *** Sorting out Cluster<->MCParticle : DONE" << std::endl;
  streamlog_out( DEBUG6 ) << " *** Guessing Cluster<->MCParticle in cases where MCParticle<->simHit is broken (LCAL/DBD)" << std::endl;
  
  // recover missing MCParticles for neutrals, typically LCal. Not fixed in DBD ... :
  // attach the MCParticle with smallest angle to cluster
  
  // (This is from the original RecoMCTruthLinker. I didn't revise it/MB)
  
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
      
      streamlog_out( DEBUG5 ) 
      << "  neutral cluster particle recovered"  
      << clu->getEnergy()
      << " maxProd: " << maxProd 
      << " px: " << closestMCP->getMomentum()[0]
      << " py: " << closestMCP->getMomentum()[1]                                           
      << " pz: " << closestMCP->getMomentum()[2]
      << std::endl ;
      
      clusterTruthRelNav.addRelation(   clu , closestMCP ,  1.0 ) ;
      truthClusterRelNav.addRelation(   closestMCP ,  clu , 1.0 ) ;
      // Could try to also add relation calohit <-> MCPart from this
      // fixup. However, this makes more confusion that clarification,
      // so, for now at least, leave it out
      // const CalorimeterHitVec& cluHits = clu->getCalorimeterHits() ;
      //
      // for( CalorimeterHitVec::const_iterator hitIt = cluHits.begin() ;
      //        hitIt != cluHits.end() ; ++hitIt ) { 
      //
      //  CalorimeterHit* hit = * hitIt ;  // ... a calo seen hit ...
      //  streamlog_out( DEBUG5 )
      //  <<"recuperated   hit = "<< hit << " e " << hit->getEnergy()<< std::endl;
      // chitTruthRelNav.addRelation(  hit , closestMCP , hit->getEnergy() ) ;
      // }
    }
    
    
  }

  //  seen-true relation complete. add the collection
  
  streamlog_out( DEBUG6 ) << " *** Cluster linking complete, create collection " << std::endl;
  *trclcol = truthClusterRelNav.createLCCollection() ;
  *ctrlcol = clusterTruthRelNav.createLCCollection() ;
  *chittrlcol = chitTruthRelNav.createLCCollection() ;
} 


void RecoMCTruthLinker::particleLinker(  LCCollection* mcpCol, LCCollection* particleCol, 
                                          LCCollection* ttrlcol, LCCollection* ctrlcol,
                                          LCCollection* trtlcol, LCCollection* trclcol,
                                          LCCollection** ptrlcol, LCCollection** trplcol) {

  LCRelationNavigator particleTruthRelNav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE  ) ;
  LCRelationNavigator truthParticleRelNav( LCIO::MCPARTICLE , LCIO::RECONSTRUCTEDPARTICLE ) ;

  LCRelationNavigator      trackTruthRelNav = LCRelationNavigator(  ttrlcol );
  LCRelationNavigator      clusterTruthRelNav = LCRelationNavigator(  ctrlcol );
  LCRelationNavigator      truthTrackRelNav = LCRelationNavigator(  trtlcol );
  LCRelationNavigator      truthClusterRelNav = LCRelationNavigator(  trclcol );

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
    streamlog_out( DEBUG3 ) << "      ======== " << std::endl;
    streamlog_out( DEBUG3 ) << "      Treating particle " << part->id() << " with index " << i 
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
        streamlog_out( DEBUG2 )  << "         Track " <<  tracks[j]->id() << " with index " << j <<" has " << nhit[j] << " hits " << std::endl; 
      }
      streamlog_out( DEBUG2 )  << "         Total : " << nhitT << " hits " << std::endl; 
    } else {
      nhit[0]=1 ; nhitT=1;
    }
    for (int j=0 ; j < ntrk ; j++ ) {
      
      if (  tracks[j] != 0 ) {  
        mcvec = trackTruthRelNav.getRelatedToObjects( tracks[j]);
        www = trackTruthRelNav.getRelatedToWeights( tracks[j]);
        int ntp= mcvec.size();
        streamlog_out( DEBUG3 ) << "      Track " <<  tracks[j]->id() << " with index " << j 
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
              streamlog_out( DEBUG2 ) << "        Individual track weight to " <<mcvec[kkk]<< " is " 
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
          streamlog_out( DEBUG2 )  << "        Cluster "  <<  clusters[j]->id() << " with index " << j <<" has energy " << eclu[j] << std::endl; 
        }
        streamlog_out( DEBUG2 )  << "       Total : " << ecluT << std::endl; 
      } else {
        eclu[0]=1 ; ecluT=1;
      }
      for (int j=0 ; j < nclu  ; j++ ) {
        if ( clusters[j] != 0 ) {
          mcvec =  clusterTruthRelNav.getRelatedToObjects(clusters[j]);
          www = clusterTruthRelNav.getRelatedToWeights(clusters[j]);
          int ntp= mcvec.size();
          streamlog_out( DEBUG3 ) << "    Cluster " <<  clusters[j]->id() << " with index " << j
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
                streamlog_out( DEBUG2 ) << "        Individual cluster Weight to " <<mcvec[kkk]<< " is " 
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
        streamlog_out( DEBUG4 ) << "     particle " << part->id() <<" has weight "<<mcit->second
        << " (Track: " << int(mcit->second)%10000 
        << " , Cluster: " << int(mcit->second)/10000 << " ) " 
        << " of MCParticle with pdg : " << mcit->first->getPDG() 
        << " and genstat : " <<  mcit->first->getGeneratorStatus() 
				<< " id: " << mcit->first->id() 
        << std::endl ;
        
        particleTruthRelNav.addRelation(   part ,  mcit->first ,  mcit->second ) ;
      }
    } else {
      if( mcmax != NULL ) {
        //AS: FixMe: There is still something going wrong with particles from the particle gun.
        //Although there are perfect PFOs no link with the MCParticle is established.
        particleTruthRelNav.addRelation(   part ,  mcmax ,  maxwgt ) ;
        streamlog_out( DEBUG3 ) << "      particle " << part->id() << " has weight "<<maxwgt
        << " of MCParticle with pdg : " << mcmax->getPDG() 
        << " and genstat : " <<  mcmax->getGeneratorStatus() 
				<< " id: " << mcmax->id() 
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
  // The reverse map. Differs in what the weight means: here it means:
  // "this cluster got wgt of all seen cluster energy the true produced"
  // (in the  other one it means:
  // "this true contributed wgt to the total seen energy of the cluster")
  LCObjectVec partvec; 
  static FloatVec www_clu ;
  static FloatVec www_trk ;
  LCObjectVec cluvec_t;
  LCObjectVec trkvec_t;
  ClusterVec cluvec_p;
  TrackVec trkvec_p;
  int nMCP  = mcpCol->getNumberOfElements() ;
  for(int i=0;i< nMCP;i++){
    MCParticle* mcp =dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    partvec  = particleTruthRelNav.getRelatedFromObjects(mcp);
    cluvec_t = truthClusterRelNav.getRelatedToObjects(mcp);
    www_clu  = truthClusterRelNav.getRelatedToWeights(mcp);
    trkvec_t = truthTrackRelNav.getRelatedToObjects(mcp);
    www_trk  = truthTrackRelNav.getRelatedToWeights(mcp);

    for ( unsigned j=0 ;  j<partvec.size() ; j++ ) {
      ReconstructedParticle* msp = dynamic_cast<ReconstructedParticle*>(partvec[j]) ;
      cluvec_p = msp->getClusters() ;
      trkvec_p = msp->getTracks() ;
      float  c_wgt=0. ;
      for ( unsigned k=0 ; k<cluvec_p.size() ; k++ ) {
        for ( unsigned l=0 ; l<cluvec_t.size() ; l++ ) {
          if ( cluvec_p[k] == cluvec_t[l] ) {
            c_wgt+=www_clu[l];
          }
        }
      }
      float  t_wgt=0. ;
      for ( unsigned k=0 ; k<trkvec_p.size() ; k++ ) {
        for ( unsigned l=0 ; l<trkvec_t.size() ; l++ ) {
          if ( trkvec_p[k] == trkvec_t[l] ) {
            t_wgt+=www_trk[l];
          }
        }
      }
      float wgt=int(c_wgt*1000)*10000 + int(t_wgt*1000) ;
      truthParticleRelNav.addRelation(   mcp, msp  , wgt ) ;
      streamlog_out( DEBUG4 ) << "    True Particle " << mcp->id() << " ( pdg " << mcp->getPDG()  
			      << " ) has weight " << c_wgt <<" / " << t_wgt << "( " << int(wgt) << " ) to particle "  
			      << msp->id() << "  with " << cluvec_p.size() << " clusters, and " 
                     << trkvec_p.size() << " tracks " <<  std::endl; 
    }
  }
  streamlog_out( DEBUG6 ) << " particle linking complete, create collection " << std::endl;
  *ptrlcol = particleTruthRelNav.createLCCollection() ;
  *trplcol = truthParticleRelNav.createLCCollection() ;
}
void RecoMCTruthLinker::makeSkim(   LCCollection* mcpCol ,  LCCollection* ttrlcol,  LCCollection* ctrlcol ,  LCCollectionVec** skimVec){
  
  
  LCRelationNavigator      trackTruthRelNav = LCRelationNavigator(  ttrlcol );
  LCRelationNavigator      clusterTruthRelNav = LCRelationNavigator(  ctrlcol );
  
  
  //-------------- create skimmed MCParticle collection ------------------------
  
  //  *skimVec = new LCCollectionVec( LCIO::MCPARTICLE )  ;
  (*skimVec)->setSubset( true) ;  // flag as subset 
  
  
  int nMCP  = mcpCol->getNumberOfElements() ;
  
  for(int i=0; i< nMCP ; i++){
    
    MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    
    
    if( mcp->ext<MCPKeep>() == true ){
      
      continue ;    // particle allready in skim 
    }
    
    
    
    
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
        
        streamlog_out( DEBUG5 ) << " keep MCParticle - e :" << mcp->getEnergy()  
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
  
  streamlog_out( DEBUG5 ) << " First loop done. Now search for KeepDaughtersPDG:s" << std::endl;
  
  
  if(_saveBremsstrahlungPhotons){
    for(int i=0; i< nMCP ; i++){
      MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
      if( ( abs(mcp->getPDG()) == 22 ) && ( mcp->getEnergy() > _bremsstrahlungEnergyCut ) ){
        if( mcp->getParents().size() ){
          MCParticle* parent = mcp->getParents()[0];
          if( abs(parent->getPDG()) == 11){
            const float x = mcp->getVertex()[0];
            const float y = mcp->getVertex()[1];
            const float z = mcp->getVertex()[2];
            const float xpe = parent->getEndpoint()[0];
            const float ype = parent->getEndpoint()[1];
            const float zpe = parent->getEndpoint()[2];
            const float dx  = x - xpe;
            const float dy  = y - ype;
            const float dz  = z - zpe;
            const float dr = sqrt(dx*dx+dy*dy+dz*dz);
            if(dr>100.){
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
        
        streamlog_out( DEBUG5 ) << " keeping daughters of particle with pdg : " << mcp->getPDG() << " : " 
        << " [" << mcp->getGeneratorStatus() << "] :";
        //                                << " e :" << mcp->getEnergy() 
        //                                << " isCreatedInSimulation :" << mcp->isCreatedInSimulation() << std::endl
        //                                << " isBackscatter :" << mcp->isBackscatter() << std::endl
        //                                << " vertexIsNotEndpointOfParent :" << mcp->vertexIsNotEndpointOfParent()     << std::endl
        //                                << " isDecayedInTracker :" << mcp->isDecayedInTracker()       << std::endl
        //                                << " isDecayedInCalorimeter :" << mcp->isDecayedInCalorimeter()       << std::endl
        //                                << " hasLeftDetector :" << mcp->hasLeftDetector()     << std::endl
        //                                << " isStopped :" << mcp->isStopped()    << "  : " 
        
        streamlog_message( DEBUG5 , 
                          if( mcp->getParents().size() ) ,
                          " parent pdg : " << mcp->getParents()[0]->getPDG() << " : "  ;
                          ) ;
        
        
        //      << std::endl ;
        
        for( MCParticleVec::const_iterator dIt = daughters.begin() ;
            dIt != daughters.end() ; ++dIt ){
          
          
          MCParticle* dau = dynamic_cast<MCParticle*>( *dIt ) ;
          
          if( dau->getEnergy()*1000. >  _eCutMeV ) {
            
            (*dIt)->ext<MCPKeep>() = true ;
            
            streamlog_out( DEBUG5 ) <<  (*dIt)->getPDG() << ", " ;
          }
        }
        
        streamlog_out( DEBUG5 ) << std::endl ;
        
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
  
  streamlog_out( DEBUG3 ) << " keepMCParticle keep particle with pdg : " << mcp->getPDG() 
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
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG9) << "RecoMCTructh::Check(): MCParticle collection \"" << _mcParticleCollectionName << "\" does not exist, skipping" << std::endl;
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
    streamlog_out(DEBUG9) << "RecoMCTructh::Check(): MCParticleSkimmed collection \"" << _mcParticlesSkimmedName << "\" does not exist, skipping" << std::endl;
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

const LCObjectVec* RecoMCTruthLinker::getSimHits( TrackerHit* trkhit, const FloatVec* weights ){
  
  const LCObjectVec* obj = & _navMergedTrackerHitRel->getRelatedToObjects(trkhit);
  
  if( obj->empty() == false  ) { 

    if(weights != 0 ) weights = & _navMergedTrackerHitRel->getRelatedToWeights(trkhit);
    
  }
  else {
    streamlog_out( WARNING ) << "getSimHits :  TrackerHit : " << trkhit << " has no sim hits related. CellID0 = " << 
         trkhit->getCellID0() << " pos = " << trkhit->getPosition()[0] << " " << trkhit->getPosition()[1] << " " << 
         trkhit->getPosition()[2] << std::endl ;
  }
  
  return obj;
  
}


void RecoMCTruthLinker::end(){ 
  
  streamlog_out(DEBUG6) << " processed " << _nEvt << " events in " << _nRun << " runs "
  << std::endl ;
  
}


/** helper function to get collection safely */
inline lcio::LCCollection* getCollection(lcio::LCEvent* evt, const std::string name ){
  
  if( name.size() == 0 )
    return 0 ;
  
  try{
    
    return evt->getCollection( name ) ;
    
  } catch( lcio::DataNotAvailableException& e ){
    
    streamlog_out( DEBUG2 ) << "getCollection :  DataNotAvailableException : " << name <<  std::endl ;
    
    return 0 ;
  }
}


void RecoMCTruthLinker::mergeTrackerHitRelations(LCEvent * evt){
  
  unsigned nCol = _colNamesTrackerHitRelations.size() ;
  
  //--- copy existing collections to a vector first
  std::vector<LCCollection*> colVec ;
  
  for( unsigned i=0; i < nCol ; ++i) {
    
    LCCollection* col  =  getCollection ( evt , _colNamesTrackerHitRelations[i] ) ;
    
    if( col != 0 ){ 
      
      colVec.push_back( col ) ;
    } else {
      
      streamlog_out(DEBUG2) << " mergeTrackerHitRelations: input collection missing : " << _colNamesTrackerHitRelations[i] << std::endl ;
    }
  }
  
  
   
  streamlog_out( DEBUG2 ) <<  " mergeTrackerHitRelations: copied collection parameters ... " << std::endl ;

  
  
  nCol = colVec.size() ;
  
  for( unsigned i=0; i < nCol ; ++i) {
    
    LCCollection* col  =  colVec[i] ;
    
    if( i == 0 ){
      // copy collection flags and collection parameters from first collections
      
      _mergedTrackerHitRelCol = new LCCollectionVec( col->getTypeName() )  ;

      
      _mergedTrackerHitRelCol->setFlag( col->getFlag() ) ;
 
      StringVec stringKeys ;
      col->getParameters().getStringKeys( stringKeys ) ;
      for(unsigned i=0; i< stringKeys.size() ; i++ ){
        StringVec vals ;
        col->getParameters().getStringVals(  stringKeys[i] , vals ) ;
        _mergedTrackerHitRelCol->parameters().setValues(  stringKeys[i] , vals ) ;   
      }
      StringVec intKeys ;
      col->getParameters().getIntKeys( intKeys ) ;
      for(unsigned i=0; i< intKeys.size() ; i++ ){
        IntVec vals ;
        col->getParameters().getIntVals(  intKeys[i] , vals ) ;
        _mergedTrackerHitRelCol->parameters().setValues(  intKeys[i] , vals ) ;   
      }
      StringVec floatKeys ;
      col->getParameters().getFloatKeys( floatKeys ) ;
      for(unsigned i=0; i< floatKeys.size() ; i++ ){
        FloatVec vals ;
        col->getParameters().getFloatVals(  floatKeys[i] , vals ) ;
        _mergedTrackerHitRelCol->parameters().setValues(  floatKeys[i] , vals ) ;   
      }

      
      
    }
    
    int nEle = col->getNumberOfElements() ;
    
    for(int j=0; j < nEle ; ++j){
      
      _mergedTrackerHitRelCol->addElement(  col->getElementAt(j) ) ;
      
    }
    
  }    
   
  if( nCol != 0 ) _navMergedTrackerHitRel = new LCRelationNavigator( _mergedTrackerHitRelCol );
  
}







void RecoMCTruthLinker::mergeCaloHitRelations(LCEvent * evt){
  
  unsigned nCol = _caloHitRelationNames.size() ;
  
  //--- copy existing collections to a vector first
  std::vector<LCCollection*> colVec ;
  
  for( unsigned i=0; i < nCol ; ++i) {
    
    LCCollection* col  =  getCollection ( evt , _caloHitRelationNames[i] ) ;
    
    if( col != 0 ){ 
      
      colVec.push_back( col ) ;
      
    } else {
      
      streamlog_out(DEBUG2) << " mergeCaloHitRelations: input collection missing : " << _caloHitRelationNames[i] << std::endl ;
    }
  }
  
  
   
  streamlog_out( DEBUG2 ) <<  " mergeCaloHitRelations: copied collection parameters ... " << std::endl ;

  
  
  nCol = colVec.size() ;
  
  for( unsigned i=0; i < nCol ; ++i) {
    
    LCCollection* col  =  colVec[i] ;
    
    if( i == 0 ){
      // copy collection flags and collection parameters from first collections
      
      _mergedCaloHitRelCol = new LCCollectionVec( col->getTypeName() )  ;

      
      _mergedCaloHitRelCol->setFlag( col->getFlag() ) ;
 
      StringVec stringKeys ;
      col->getParameters().getStringKeys( stringKeys ) ;
      for(unsigned i=0; i< stringKeys.size() ; i++ ){
        StringVec vals ;
        col->getParameters().getStringVals(  stringKeys[i] , vals ) ;
        _mergedCaloHitRelCol->parameters().setValues(  stringKeys[i] , vals ) ;   
      }
      StringVec intKeys ;
      col->getParameters().getIntKeys( intKeys ) ;
      for(unsigned i=0; i< intKeys.size() ; i++ ){
        IntVec vals ;
        col->getParameters().getIntVals(  intKeys[i] , vals ) ;
        _mergedCaloHitRelCol->parameters().setValues(  intKeys[i] , vals ) ;   
      }
      StringVec floatKeys ;
      col->getParameters().getFloatKeys( floatKeys ) ;
      for(unsigned i=0; i< floatKeys.size() ; i++ ){
        FloatVec vals ;
        col->getParameters().getFloatVals(  floatKeys[i] , vals ) ;
        _mergedCaloHitRelCol->parameters().setValues(  floatKeys[i] , vals ) ;   
      }

      
      
    }
    
    int nEle = col->getNumberOfElements() ;
    
    for(int j=0; j < nEle ; ++j){
      
      _mergedCaloHitRelCol->addElement(  col->getElementAt(j) ) ;
      
    }
    
  }    
   
  if( nCol != 0 ) _navMergedCaloHitRel = new LCRelationNavigator( _mergedCaloHitRelCol );
  
}

const LCObjectVec* RecoMCTruthLinker::getCaloHits( CalorimeterHit* calohit, const FloatVec* weights ){
  
  const LCObjectVec* obj = & _navMergedCaloHitRel->getRelatedToObjects(calohit);
  
  if( obj->empty() == false  ) { 

    if(weights != 0 ) weights = & _navMergedCaloHitRel->getRelatedToWeights(calohit);
    
  }
  else {
    
    streamlog_out( DEBUG5 ) << "getCaloHits :  CalorimeterHit : " << calohit << " has no sim hits related. CellID0 = " << 
       calohit->getCellID0() << " pos = " << calohit->getPosition()[0] << " " << calohit->getPosition()[1] << " " << 
       calohit->getPosition()[2] << std::endl ;
  }
  return obj;
  
}



void RecoMCTruthLinker::linkPrinter (  LCCollection* mcpCol, LCCollection* particleCol, LCCollection* ptrlcol, LCCollection* trplcol) {

  LCRelationNavigator      particleTruthRelNav = LCRelationNavigator(  ptrlcol );
  LCRelationNavigator      truthParticleRelNav  = LCRelationNavigator(  trplcol );

  if ( ! _FullRecoRelation ) {
    // Not useful in this case ...
    return;
  }
  static FloatVec www ; 
  static FloatVec www_other_way ; 

  LCObjectVec mcvec; 
  int nPart = particleCol->getNumberOfElements() ;
  LCObjectVec partvec; 

  for(int i=0;i<nPart;++i){
    
    ReconstructedParticle* part = dynamic_cast<ReconstructedParticle*> ( particleCol->getElementAt(i) ) ;
    double clu_e=0.0;
    if ( part->getClusters().size() > 0 ) {
      Cluster* clu = part->getClusters()[0];
      if ( clu != 0 ) {
        clu_e=clu->getEnergy();
      }
    }
    mcvec = particleTruthRelNav.getRelatedToObjects( part);
    www = particleTruthRelNav.getRelatedToWeights( part);
    int ntp= mcvec.size();
    streamlog_out( DEBUG5 ) << "    Particle " <<  part->id() << " (q: " <<  int(part->getCharge()) << 
      " ) has " << ntp << " true particles. E= "<< part->getEnergy() << " " << clu_e ;
    streamlog_out( DEBUG3 )<< "       Index, id, PDG and energy of contributors: " << std::endl;  
    if ( ntp > 0 ) {
      if ( mcvec[0] != 0 ) {
        double total_trk_weight = 0.0;
        double total_clu_weight = 0.0;
        double total_e_from_neutrals = 0.0;
        double true_charged_E=0.0;
        double total_ecalo_in_this_true=0.;
        double total_ecalo_neutral_in_charged=0.;
        for ( int kkk=0 ; kkk < ntp ; kkk++ ) {
          total_trk_weight+=(int(www[kkk])%10000)/1000.0;
          total_clu_weight+=(int(www[kkk])/10000)/1000.0;
          MCParticle* mcp = dynamic_cast< MCParticle*>( mcvec[kkk]);
          streamlog_out( DEBUG3 ) <<"       "<< kkk << " " << mcp->id()  << " " << mcp->getPDG()  << " " << mcp->getEnergy()  << std::endl ;
          if ( part->getCharge() != 0 ) {
            if ( mcp->getCharge() != 0 && (int(www[kkk])%10000)/1000.0 > 0.6 ) {
              true_charged_E=mcp->getEnergy();
              double ecalo_from_this_true=((int(www[kkk])/10000)/1000.0)*clu_e;
              if ( ecalo_from_this_true != 0.0 ) {
                partvec = truthParticleRelNav.getRelatedToObjects( mcp);
                www_other_way =  truthParticleRelNav.getRelatedToWeights( mcp);
                for ( unsigned jjj=0 ; jjj < partvec.size() ; jjj++ ) {
                  if ( partvec[jjj] == part ) {
                    total_ecalo_in_this_true = ecalo_from_this_true/(int(www_other_way[jjj]/10000)/1000.0) ;
                    break;
                  }
                }
              }
            } else {
              total_ecalo_neutral_in_charged+=((int(www[kkk])/10000)/1000.0)*clu_e;
            }
          } else {
            if ( mcp->getCharge() == 0 ) {
              total_e_from_neutrals+=((int(www[kkk])/10000)/1000.0)*clu_e;
	    }
	  }
        }
        if ( part->getCharge() != 0 ) {
	  streamlog_out( DEBUG5 ) << "     True E_ch " << true_charged_E << " total E_calo from true " << total_ecalo_in_this_true 
                                  << " E_calo from neutral " << total_ecalo_neutral_in_charged << " track weight " << 
          total_trk_weight << " cluster weight " << total_clu_weight  << std::endl ;
        } else {
	  streamlog_out( DEBUG5 ) << "     True E_neutral " <<  total_e_from_neutrals << " track weight " << 
          total_trk_weight << " cluster weight " << total_clu_weight  << std::endl ;
        }
      }
    }
  }  

  int nMCPart = mcpCol->getNumberOfElements() ;

  for(int i=0;i<nMCPart;++i){
    
    MCParticle* mcpart = dynamic_cast<MCParticle*> ( mcpCol->getElementAt(i) ) ;

    partvec = truthParticleRelNav.getRelatedToObjects( mcpart);
    www =  truthParticleRelNav.getRelatedToWeights( mcpart);
    int nsp= partvec.size();
    if ( nsp > 0 ) {
      streamlog_out( DEBUG5 ) << "    TrueParticle " <<  mcpart->id() << " (q: " <<  int(mcpart->getCharge()) << 
        " ) has " << nsp << " seen particles. E= "<<mcpart->getEnergy() << " pt= " <<
	sqrt(mcpart->getMomentum()[0]* mcpart->getMomentum()[0]+ mcpart->getMomentum()[1]* mcpart->getMomentum()[1]);
      streamlog_out( DEBUG3  )<< "       index, id, and enregy of reco particle: " << std::endl;  
      if ( partvec[0] != 0 ) {
        double total_trk_weight = 0.0;
        double total_clu_weight = 0.0;
        for ( int kkk=0 ; kkk < nsp ; kkk++ ) {
          total_trk_weight+=(int(www[kkk])%10000)/1000.0;
          total_clu_weight+=(int(www[kkk])/10000)/1000.0;
           ReconstructedParticle* recopart = dynamic_cast< ReconstructedParticle*>(partvec[kkk]);
	   streamlog_out( DEBUG3 )<<"       "<< kkk << " " << recopart->id()  << " " << recopart->getEnergy()    << std::endl ; 
        }
        streamlog_out( DEBUG5 ) << "    Total track weight " << total_trk_weight << " , Total cluster weight " << total_clu_weight << std::endl ;
      }
    }
  } 
}



