#include "PIDutil.h"
#include <cmath>
#include <set>
#include <TF2.h>

#include <marlin/VerbosityLevels.h>

#include <EVENT/ReconstructedParticle.h>

LCObject *getRelated(LCObject *p, LCRelationNavigator const &nav, int &_cw, int &_tw)
{
    int wcut = 500;
    auto weights = nav.getRelatedToWeights(p);
    int max_ci = -1;
    int max_ti = -1;
    int max_cw = 0;
    int max_tw = 0;
    for (size_t i = 0; i < weights.size(); i++)
    {
        float w = weights[i];
        // weights in permille
        int cw = int(w) / 10000;
        int tw = int(w) % 10000;
        if (cw > max_cw)
        {
            max_cw = cw;
            max_ci = i;
        }
        if (tw > max_tw)
        {
            max_tw = tw;
            max_ti = i;
        }
    }
    if (max_ci != max_ti)
    {
        // confusion, skip particle
        return nullptr;
    }
    if (max_cw < wcut || max_tw < wcut)
    {
        // too messy, skip
        return nullptr;
    }

    _cw = max_cw;
    _tw = max_tw;
    auto objects = nav.getRelatedToObjects(p);
    return objects[max_ci];

}

WeightedPoints3D getWeightedPoints3D(const Cluster *clu, const LCCollection *PandoraClusters)
{
    auto hits = clu->getCalorimeterHits();
    size_t n_hits = hits.size();
    if (n_hits > 0) {
        // we have a REC file and can use all details
        std::vector<double> ehits, xhits, yhits, zhits;
        ehits.reserve(n_hits);
        xhits.reserve(n_hits);
        yhits.reserve(n_hits);
        zhits.reserve(n_hits);
        for (const auto *hit: hits) {
            ehits.push_back(hit->getEnergy());
            auto pos = hit->getPosition();
            xhits.push_back(pos[0]);
            yhits.push_back(pos[1]);
            zhits.push_back(pos[2]);
        }
        return WeightedPoints3D(int(n_hits), ehits.data(), xhits.data(), yhits.data(), zhits.data());
    } else {
        // we have a DST file and have to call the more complicated constructor and some stuff will be approximate
        auto *coga = clu->getPosition();
        std::vector<double> cogv(coga, coga + 3);

        auto shapes = clu->getShape();

        StringVec shape_keys;
        PandoraClusters->getParameters().getStringVals("ClusterShapeParameters", shape_keys);
        int npoints = 0;
        double wgt_sum = 0;
        double wgt2_sum = 0;
        double wgt4_sum = 0;
        streamlog_out(DEBUG) << "shape_keys: [";
        for (unsigned k = 0; k < shape_keys.size(); k++) {
            streamlog_out(DEBUG) << shape_keys[k] << ", ";
            if (shape_keys[k] == "npoints") {
                npoints = int(shapes[k]);
            } else if (shape_keys[k] == "sum_wgt") {
                wgt_sum = shapes[k];
            } else if (shape_keys[k] == "sum_wgt^2") {
                wgt2_sum = shapes[k];
            } else if (shape_keys[k] == "sum_wgt^4") {
                wgt4_sum = shapes[k];
            }
        }
        streamlog_out(DEBUG) << "]" << std::endl;

        auto pos_error = clu->getPositionError();
        double seen_covmat[3][3];
        int n = 0;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j <= i; j++) {
                seen_covmat[j][i] = pos_error[n] * wgt_sum * wgt_sum / wgt2_sum;
                seen_covmat[i][j] = pos_error[n] * wgt_sum * wgt_sum / wgt2_sum;
                n++;
            }
        }

        std::vector<double> cov_fr_clu;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                cov_fr_clu.push_back(seen_covmat[i][j]);
            }
        }

        auto dir_error = clu->getDirectionError();
        std::vector<double> thphcov(dir_error.begin(), dir_error.end());

        return WeightedPoints3D(cogv, cov_fr_clu, thphcov, npoints, wgt_sum, wgt2_sum, wgt4_sum);
    }
}

MCParticle *getMCParent(const MCParticle *mcp)
{
    MCParticle *parent = nullptr;

    auto parents = mcp->getParents();
    if (parents.size() > 0) {
        if (parents[0]->getPDG() != mcp->getPDG())
        {
            parent = parents[0];
        } else {
            // recursively walk up the chain to the real parent
            parent = getMCParent(parents[0]);
        }
    }

    return parent;
}

MCParticle *getMCParticle(ReconstructedParticle *p, LCRelationNavigator const &RecoMCTruthNavigator, LCRelationNavigator const &MCTruthRecoNavigator)
{
    int RecoMC_cw = 0;
    int RecoMC_tw = 0;
    auto mcp_candidate = getRelated(p, RecoMCTruthNavigator, RecoMC_cw, RecoMC_tw);
    if (mcp_candidate == nullptr)
    {
        // could not find suitable mc particle, skip
        return nullptr;
    }

    // reverse check if mc particle also had most of its contributions in our reco particle
    int MCReco_cw = 0;
    int MCReco_tw = 0;
    auto reco_candidate = getRelated(mcp_candidate, MCTruthRecoNavigator, MCReco_cw, MCReco_tw);
    if (reco_candidate != p)
    {
        // p got most of its track/cluster hits from mcp but mcp contributed stronger to other reco particles, skip
        return nullptr;
    }

    // now cast mcp_candidate to a real MCParticle to use it
    return dynamic_cast<MCParticle *>(mcp_candidate);
}