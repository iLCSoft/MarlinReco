from Configurables import EDM4hep2LcioTool, EventDataSvc, Lcio2EDM4hepTool, MarlinProcessorWrapper
from k4FWCore import ApplicationMgr, IOSvc

io_svc = IOSvc("IOSvc")
io_svc.Input = "truth_link_events.edm4hep.root"
io_svc.Output = "truth_links.edm4hep.root"

edm4hep_to_lcio = EDM4hep2LcioTool("EDM4hep2Lcio")
edm4hep_to_lcio.convertAll = True
lcio_to_edm4hep = Lcio2EDM4hepTool("Lcio2EDM4hep")
lcio_to_edm4hep.convertAll = False
lcio_to_edm4hep.collNameMapping = {
    "CaloHitMCParticleLinks": "CaloHitMCParticleLinks",
    "ClusterMCParticleLinks": "ClusterMCParticleLinks",
}

reco_mc_truth_linker = MarlinProcessorWrapper("RecoMCTruthLinker")
reco_mc_truth_linker.ProcessorType = "RecoMCTruthLinker"
reco_mc_truth_linker.Parameters = {
    "MCParticleCollection": ["MCParticles"],
    "TrackCollection": ["EmptyTracks"],
    "ClusterCollection": ["Clusters"],
    "SimCaloHitCollections": ["ECalSimHits", "HCalSimHits"],
    "SimCalorimeterHitRelationNames": ["ECalLinks", "HCalLinks"],
    "CalohitMCTruthLinkName": ["CaloHitMCParticleLinks"],
    "ClusterMCTruthLinkName": ["ClusterMCParticleLinks"],
    "Verbosity": ["WARNING"],
}
reco_mc_truth_linker.EDM4hep2LcioTool = edm4hep_to_lcio
reco_mc_truth_linker.Lcio2EDM4hepTool = lcio_to_edm4hep

ApplicationMgr(
    TopAlg=[reco_mc_truth_linker],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
)
