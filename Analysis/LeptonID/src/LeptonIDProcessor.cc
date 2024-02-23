#include "LeptonIDProcessor.h"
#include "PIDutil.h"
#include <iostream>
#include <cmath>
#include <optional>

#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>

// ----- include for verbosity dependent logging ---------
#include "marlin/VerbosityLevels.h"

#include <marlin/AIDAProcessor.h>

#include <TVector3.h>

using namespace lcio;
using namespace marlin;

LeptonIDProcessor aLeptonIDProcessor;

LeptonIDProcessor::LeptonIDProcessor() : Processor("LeptonIDProcessor"), _mvareader("!Color:!Silent")
{

    // modify processor description
    _description = "The LeptonIDProcessor is used to identify electrons and muons using TMVA classifiers. It also produces the TTrees needed to train them";

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                            "PandoraPFOsName",
                            "Name of the collection holding the PID info usually PandoraPFOs",
                            _PandoraPFOsColName,
                            std::string("PandoraPFOs"));

    registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                            "InputPFOsName",
                            "Name of the PFO collection the training or evaluation should be done on",
                            _InputPFOsColName,
                            std::string("DistilledPFOs"));

    registerInputCollection(LCIO::LCRELATION,
                            "RecoMCTruthLinkName",
                            "Name of the RecoMCTruthLink collection",
                            _RecoMCTruthLinkName,
                            std::string("RecoMCTruthLink"));

    registerInputCollection(LCIO::LCRELATION,
                            "MCTruthRecoLinkName",
                            "Name of the MCTruthRecoLink collection",
                            _MCTruthRecoLinkName,
                            std::string("MCTruthRecoLink"));

    registerInputCollection(LCIO::CLUSTER,
                            "PandoraClustersName",
                            "Name of the PandoraClusters collection",
                            _PandoraClustersColName,
                            std::string("PandoraClusters"));

    registerProcessorParameter("MCTruthRecoCweightCut",
                               "Cut used on MCTruthRecoCluster weight",
                               _MCTruthRecoCweightCut,
                               0.5F);

    registerProcessorParameter("MCTruthRecoTweightCut",
                               "Cut used on MCTruthRecoTrack weight",
                               _MCTruthRecoTweightCut,
                               0.5F);

    registerProcessorParameter("RecoMCTruthCweightCut",
                               "Cut used on RecoMCTruthCluster weight",
                               _RecoMCTruthCweightCut,
                               0.8F);

    registerProcessorParameter("RecoMCTruthTweightCut",
                               "Cut used on RecoMCTruthTrack weight",
                               _RecoMCTruthTweightCut,
                               0.5F);

    registerProcessorParameter("BuildTree",
                               "Decides if the Processor builds a TTree that can be used for MVA training",
                               _buildTree,
                               true);

    registerProcessorParameter("EvalMVA",
                               "Decides if the Processor does the MVA evaluation (requires weight files)",
                               _evalMVA,
                               false);

    registerProcessorParameter("weightfile",
                               "path to the BDTG weightfile",
                               _weightfile,
                               std::string("weights/PID_multi_jet_dEdx_800t_cm50r80_BDTG.weights.xml"));

    registerProcessorParameter("mvaname",
                               "Name given to the MVA method during training, e.g. BDTG.",
                               _mvaname,
                               std::string("BDTG"));

    registerProcessorParameter("pidname",
                               "Name used to store the PID results in the ReconstructedParticle collection.",
                               _pidname,
                               std::string("LeptonID"));

    registerProcessorParameter("dEdxname",
                               "Name used to read the dEdx distance from the PandoraPFO collection parameters.",
                               _dEdxname,
                               std::string("dEdxPIDv2"));

std::vector<std::string> usedVars = {
    "e_over_p",
    "ecal_share",
    "seenEcalDep",
    "seenYokeDep",
    "shape0",
    "shape1",
    "shape8",
    "shape17",
    "shape19",
    "shape20",
    "shape21",
    "shape22",
    "shape23",
    "cluEllipsoid_r1",
    "cluEllipsoid_r2",
    "cluEllipsoid_r3",
    "cluEllipsoid_r_ave",
    "cluEllipsoid_density",
    "cluEllipsoid_eccentricity_T",
    "cluEllipsoid_eccentricity_L",
    "dEdxDist_e"
};

    registerProcessorParameter("UsedVariables",
                               "List of variables used during training",
                               _usedVariables,
                               usedVars);
}

void LeptonIDProcessor::setupTree()
{
    _tree->Branch("truePDG", &_truePDG);
    _tree->Branch("parentPDG", &_parentPDG);

    for (auto &v : _vars)
    {
        _tree->Branch(v.first.c_str(), &(v.second));
    }
}

void LeptonIDProcessor::setupMVAReader()
{

    for (const auto &varname : _usedVariables) {
        _mvareader.AddVariable(varname, &_vars[varname]);
    }
    // XXX: the reader also needs the spectator variables but does not use them...
    _mvareader.AddSpectator("abs(truePDG)", (int*)nullptr);
    _mvareader.AddSpectator("abs(parentPDG)", (int*)nullptr);
    _mvareader.AddSpectator("seenP", (float*)nullptr);
    _mvareader.AddSpectator("seenCosTheta", (float*)nullptr);

    _mvareader.BookMVA(_mvaname, _weightfile);
}

void LeptonIDProcessor::init()
{
    // usually a good idea to
    printParameters();

    streamlog_out(DEBUG) << "   init called  " << std::endl;
    if (_buildTree) {
        // if I keep this here it magically appears in the aida file :)
        _tree = new TTree("pfos", "pfos");
        setupTree();
    }

    if (_evalMVA) {
        setupMVAReader();
    }

}

void LeptonIDProcessor::clear()
{
    _truePDG = 0;
    _parentPDG = 0;

    for (auto v : _vars) {
        v.second = 0.0f;
    }
}

void LeptonIDProcessor::processRunHeader(LCRunHeader *run)
{
    (void) run;
}

void LeptonIDProcessor::processEvent(LCEvent *evt)
{
    streamlog_out(DEBUG) << "process Event called" << std::endl;
    LCCollection *PandoraPFOs = nullptr;
    LCCollection *InputPFOs = nullptr;
    LCCollection *RecoMCTruthLink = nullptr;
    LCCollection *MCTruthRecoLink = nullptr;
    LCCollection *PandoraClusters = nullptr;

    try
    {
        PandoraPFOs = evt->getCollection(_PandoraPFOsColName);
        InputPFOs = evt->getCollection(_InputPFOsColName);
        PandoraClusters = evt->getCollection(_PandoraClustersColName);
        if (_buildTree) {
            RecoMCTruthLink = evt->getCollection(_RecoMCTruthLinkName);
            MCTruthRecoLink = evt->getCollection(_MCTruthRecoLinkName);
        }
    }
    catch (lcio::DataNotAvailableException const& e)
    {
        streamlog_out(WARNING) << "atleast one collection not available" << std::endl;
        streamlog_out(WARNING) << e.what() << std::endl;
        return;
    }

    std::optional<LCRelationNavigator> RecoMCTruthNavigator;
    std::optional<LCRelationNavigator> MCTruthRecoNavigator;
    if (_buildTree) {
        RecoMCTruthNavigator = LCRelationNavigator(RecoMCTruthLink);
        MCTruthRecoNavigator = LCRelationNavigator(MCTruthRecoLink);
    }

    PIDHandler pidh(PandoraPFOs);
    int algoID = pidh.addAlgorithm(_pidname, {"mva_e", "mva_mu", "mva_had"});
    int dEdxAlgId = pidh.getAlgorithmID(_dEdxname);
    int dEdxDist_e_Id = pidh.getParameterIndex(dEdxAlgId, "electron_dEdxdistance");

    int n = InputPFOs->getNumberOfElements();

    for (int i = 0; i < n; i++)
    {
        streamlog_out(DEBUG) << "pfo loop: " << i << std::endl;
        auto *p = dynamic_cast<ReconstructedParticle *>(InputPFOs->getElementAt(i));

        // We only want to look at charged particles with exactly one track and cluster for now
        if (abs(p->getCharge()) != 1.0) {
            continue;
        }

       const auto& tracks = p->getTracks();
        int nTracks = tracks.size();
        if (nTracks != 1) {
            continue;
        }

       const auto& clusters = p->getClusters();
        int nClusters = clusters.size();
        if (nClusters != 1) {
            continue;
        }

        if (_buildTree) {
            MCParticle* mcp = getMCParticle(p, *RecoMCTruthNavigator, *MCTruthRecoNavigator, _RecoMCTruthCweightCut, _RecoMCTruthTweightCut, _MCTruthRecoCweightCut, _MCTruthRecoTweightCut);
            if (mcp == nullptr) {
                continue;
            }
            _truePDG = mcp->getPDG();

            auto parent = getMCParent(mcp);
            int parentPDG = 0;
            if (parent != nullptr) {
                parentPDG = parent->getPDG();
            }
            _parentPDG = parentPDG;
        }

        TVector3 P(p->getMomentum());
        _vars["seenP"] = P.Mag();
        _vars["seenLogP"] = std::log(P.Mag());
        _vars["seenCosTheta"] = P.CosTheta();
        _vars["seenE"] = p->getEnergy();

       const auto& dEdxParameters = pidh.getParticleID(p, dEdxAlgId).getParameters();
        if (dEdxParameters.size() > (size_t) dEdxDist_e_Id) {
            _vars["dEdxDist_e"] = dEdxParameters[dEdxDist_e_Id];
        }

        auto cluster = clusters[0];
        getShowerShapes(cluster);

        // TODO: access subdetector energies by name as for the shapes
       const auto& subEnergies = cluster->getSubdetectorEnergies();
        if (subEnergies.size() >= 2) {
            _vars["seenEcalDep"] = subEnergies[0];
            _vars["seenHcalDep"] = subEnergies[1];
            if (subEnergies[0] > 0) {
                _vars["ecal_share"] = subEnergies[0] / (subEnergies[0] + subEnergies[1]);
            }
        }
        if (subEnergies.size() >= 3) {
            _vars["seenYokeDep"] = subEnergies[2];
        }

        _vars["e_over_p"] = _vars["seenEcalDep"] / _vars["seenP"];

        auto wgtp = getWeightedPoints3D(cluster, PandoraClusters);
        getClusterShapes(wgtp);

        if (_buildTree && std::isfinite(_vars["cluEllipsoid_r3"])) {
            _tree->Fill();
        }
        if (_evalMVA) {
            auto mvares = _mvareader.EvaluateMulticlass(_mvaname);
            int pidClassID = std::distance(mvares.begin(), std::max_element(mvares.begin(), mvares.end()));

            int pdg[3] = {11, 13, 211};

            pidh.setParticleID(p, pidClassID, pdg[pidClassID], mvares[pidClassID], algoID, mvares);
        }
        clear();
    }
}

void LeptonIDProcessor::check(LCEvent *)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void LeptonIDProcessor::end()
{
}

void LeptonIDProcessor::getClusterShapes(WeightedPoints3D &wgpt)
{
    _vars["cluEllipsoid_r1"] = wgpt.getElipsoid_r1();
    _vars["cluEllipsoid_r2"] = wgpt.getElipsoid_r2();
    _vars["cluEllipsoid_r3"] = wgpt.getElipsoid_r3();
    _vars["cluEllipsoid_vol"] = wgpt.getElipsoid_vol();
    _vars["cluEllipsoid_r_ave"] = wgpt.getElipsoid_r_ave();
    if (_vars["cluEllipsoid_vol"] > 0) {
        _vars["cluEllipsoid_density"] = wgpt.getElipsoid_density();
    }
    _vars["cluEllipsoid_eccentricity_T"] = wgpt.getTransverseElipsis_eccentricity();
    _vars["cluEllipsoid_eccentricity_L"] = wgpt.getLongitudinalElipsis_eccentricity();
}

void LeptonIDProcessor::getShowerShapes(Cluster *clu)
{
   const auto& shapes = clu->getShape();
    // TODO: make sure that the size is right
    // TODO: ideally access the shapes via their name, but that first needs a fix elsewhere...
    for (size_t i = 0; i <= 23; i++)
    {
        if (std::isfinite(shapes[i])) {
            _vars["shape" + std::to_string(i)] = shapes[i];
        }
    }
}
