#include "TJjetsPFOAnalysisProcessor.h"

void TJjetsPFOAnalysisProcessor::findTrueJetParticles(LCEvent* evt) {
  for (int i_jet=0; i_jet<njets(); i_jet++) {
    ReconstructedParticleVec jet_recos = seen_partics(i_jet);

    LCRelationNavigator* relation_TrueJetFinalColourNeutral   =  new LCRelationNavigator( evt->getCollection( _finalColourNeutralLink ) );
    LCRelationNavigator* relation_TrueJetInitialColourNeutral =  new LCRelationNavigator( evt->getCollection( _initialColourNeutralLink ) );
    LCRelationNavigator* relation_TrueJetFinalElementon       =  new LCRelationNavigator( evt->getCollection( _finalElementonLink ) );
    LCRelationNavigator* relation_TrueJetInitialElementon     =  new LCRelationNavigator( evt->getCollection( _initialElementonLink ) );
    LCRelationNavigator* relation_TrueJetMCParticle           =  new LCRelationNavigator( evt->getCollection( _trueJetMCParticleLink ) );

    LCRelationNavigatorVec tj_mc_relations {
      relation_TrueJetFinalColourNeutral, relation_TrueJetInitialColourNeutral,
      relation_TrueJetFinalElementon, relation_TrueJetInitialElementon,
      relation_TrueJetMCParticle
    };

    MCParticleList unique_jet_mcs{};

    for (auto i_rel = 0; i_rel<tj_mc_relations.size(); i_rel++) {
      LCObjectVec jet_mc_objects = tj_mc_relations[i_rel]->getRelatedToObjects( getJets()->at(i_jet) );
      streamlog_out(DEBUG) << " N FinalCN " << jet_mc_objects.size() << std::endl;
      for (int i_mc=0; i_mc<jet_mc_objects.size(); i_mc++) {
        MCParticle* jet_mc = dynamic_cast<MCParticle*>(jet_mc_objects[i_mc]);
        streamlog_out(DEBUG) << " FinalCN " << jet_mc << std::endl;
        if ( jet_mc != NULL ) {
          unique_jet_mcs.insert( jet_mc );
        }
      }
    }
    MCParticleVec jet_mcs (unique_jet_mcs.begin(), unique_jet_mcs.end());

    m_jets.push_back( new JetContentPair(jet_mcs, jet_recos) );
  }
}
