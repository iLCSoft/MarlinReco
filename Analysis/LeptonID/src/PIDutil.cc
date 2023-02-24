#include "PIDutil.h"
#include <cmath>
#include <set>
#include <TF2.h>

#include <marlin/VerbosityLevels.h>

#include <EVENT/ReconstructedParticle.h>

float getCorrectdEdx(float dEdx, Track *track)
{
    // dEdx correction stolen from the dEdx correction processor
    double lambda = abs(atan(track->getTanLambda()) * 180.0 / M_PI);
    std::vector<double> par = {0.970205, 0.0007506, 4.41781e-8, 5.8222e-8};

    double f3 = 1 / (par[0] + par[1] * lambda + par[2] * pow(lambda, 2) + par[3] * pow(lambda, 3));

    float dEdx_v2 = dEdx * f3;
    return dEdx_v2;
}

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
        auto res = WeightedPoints3D(int(n_hits), ehits.data(), xhits.data(), yhits.data(), zhits.data());
        return res;
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

std::vector<float> getPosMomFromTrackState(const TrackState *ts, double B) {
    double trackPt = B * 3.0e-4 / std::fabs(ts->getOmega());
    float trackPx = trackPt * std::cos(ts->getPhi());
    float trackPy = trackPt * std::sin(ts->getPhi());
    float trackPz = trackPt * ts->getTanLambda();

    float trackXs = ts->getReferencePoint()[0] - ts->getD0() * std::sin(ts->getPhi());
    float trackYs = ts->getReferencePoint()[1] + ts->getD0() * std::cos(ts->getPhi());
    float trackZs = ts->getReferencePoint()[2] + ts->getZ0();
    return std::vector<float>{trackXs, trackYs, trackZs, trackPx, trackPy, trackPz};
}

double distFun(const double *val, const double *par)
{
    double phi = val[0];
    double v = val[1];

    double x_ref = par[0];
    double y_ref = par[1];
    double z_ref = par[2];

    double inv_Omega_minus_d_zero = par[3];

    double sin_phi_zero = par[4];
    double cos_phi_zero = par[5];

    double z_zero = par[6];
    double tan_lambda_over_Omega = par[7];
    double phi_zero = par[8];
    double inv_Omega = par[9];

    double x_third = par[10];
    double y_third = par[11];
    double z_third = par[12];

    double Px_third = par[13];
    double Py_third = par[14];
    double Pz_third = par[15];

    double x_track = x_ref + inv_Omega_minus_d_zero * sin_phi_zero - inv_Omega * sin(phi);
    double y_track = y_ref - inv_Omega_minus_d_zero * cos_phi_zero + inv_Omega * cos(phi);
    double z_track = z_ref + z_zero - tan_lambda_over_Omega * (phi - phi_zero);

    double x_line = x_third + Px_third * v;
    double y_line = y_third + Py_third * v;
    double z_line = z_third + Pz_third * v;

    double dist_x_sq = (x_track - x_line) * (x_track - x_line);
    double dist_y_sq = (y_track - y_line) * (y_track - y_line);
    double dist_z_sq = (z_track - z_line) * (z_track - z_line);

    double dist = sqrt(dist_x_sq + dist_y_sq + dist_z_sq);

    return dist;
}

double intersectTrackLine(Track *track, TVector3 momentumOfLine,
                          std::vector<double> pointOnLine,
                          Vertex *primaryVertex, TVector3 &PCAatTrack,
                          TVector3 &PCAatLine)
{
    double phi_min, phi_max, v_min, v_max;
    phi_min = track->getPhi() - 3.14159265359 / 40.0;
    phi_max = track->getPhi() + 3.14159265359 / 40.0;

    v_min = -1.0 * (momentumOfLine.Px() * (pointOnLine[0] - primaryVertex->getPosition()[0])
                    + momentumOfLine.Py() * (pointOnLine[1] - primaryVertex->getPosition()[1])
                    + momentumOfLine.Pz() * (pointOnLine[2] - primaryVertex->getPosition()[2])
                    ) / momentumOfLine.Mag2();
    v_max = 0.0;

    auto dist_f = TF2("dist", distFun, phi_min, phi_max, v_min, v_max, 16, 2);

    auto x_ref = track->getReferencePoint();
    dist_f.SetParameter(0, x_ref[0]);
    dist_f.SetParameter(1, x_ref[1]);
    dist_f.SetParameter(2, x_ref[2]);

    auto invOmega = 1.0 / track->getOmega();
    dist_f.SetParameter(3, invOmega - track->getD0());

    // TODO: check if needs correction to fit into [0, 2pi)
    float phi_zero = track->getPhi();
    dist_f.SetParameter(4, sin(phi_zero));
    dist_f.SetParameter(5, cos(phi_zero));

    dist_f.SetParameter(6, track->getZ0());
    dist_f.SetParameter(7, track->getTanLambda() * invOmega);
    dist_f.SetParameter(8, phi_zero);
    dist_f.SetParameter(9, invOmega);

    dist_f.SetParameter(10, pointOnLine[0]);
    dist_f.SetParameter(11, pointOnLine[1]);
    dist_f.SetParameter(12, pointOnLine[2]);

    dist_f.SetParameter(13, momentumOfLine[0]);
    dist_f.SetParameter(14, momentumOfLine[1]);
    dist_f.SetParameter(15, momentumOfLine[2]);

    // TODO: Yasser had an extra SetRange call, is this needed?

    double phi = phi_zero;
    double v = 0.0;

    double distance = dist_f.GetMinimumXY(phi, v);

    double xTrackPCA = x_ref[0] + (invOmega - track->getD0()) * sin(phi_zero) - sin(phi) * invOmega;
    double yTrackPCA = x_ref[1] - (invOmega - track->getD0()) * cos(phi_zero) + cos(phi) * invOmega;
    double zTrackPCA = x_ref[2] + track->getZ0() - (phi-phi_zero) * track->getTanLambda() * invOmega;

    PCAatTrack = TVector3(xTrackPCA, yTrackPCA, zTrackPCA);

    double xLine = pointOnLine[0] + momentumOfLine.Px() * v;
	double yLine = pointOnLine[1] + momentumOfLine.Py() * v;
	double zLine = pointOnLine[2] + momentumOfLine.Pz() * v;

    PCAatLine = TVector3(xLine, yLine, zLine);

    return distance;
}

const MCParticle *checkSLD(const MCParticle *mcp, const MCParticle *parent)
{
    if (mcp == nullptr || parent == nullptr)
    {
        return nullptr;
    }

    int parentPDG = abs(parent->getPDG());
    bool isCHadron = parentPDG / 100 == 4 || parentPDG / 1000 == 4;
    bool isBHadron = parentPDG / 100 == 5 || parentPDG / 1000 == 5;

    if (!(isCHadron || isBHadron))
    {
        return nullptr;
    }

    int mcPDG = mcp->getPDG();
    int mcPDGsign = (0 < mcPDG) - (mcPDG < 0);
    // nu pdg is one number higher and it should have opposite sign
    int nuPDG = -1 * mcPDGsign * (abs(mcPDG) + 1);

    auto &daughters = parent->getDaughters();
    for (const auto *d : daughters)
    {
        if (d->getGeneratorStatus() == 1 && d->getPDG() == nuPDG)
        {
            return d;
        }
    }

    return nullptr;
}

VertexVec getJetVertices(const ReconstructedParticle *jet)
{
    std::set<Vertex *> jetVertices{};
    // iterate over jet pfos & save all non prim vertices
    for (const auto *reco : jet->getParticles())
    {
        auto vtx = reco->getStartVertex();
        if (vtx != nullptr && !vtx->isPrimary())
        {
            jetVertices.insert(vtx);
        }
    }
    return VertexVec(jetVertices.begin(), jetVertices.end());
}

Vertex *getClosestVtx(const Vertex *primVtx, VertexVec jetVertices)
{
    Vertex *closestVtx = nullptr;

    const float* primVtxPos = primVtx->getPosition();

    double minDist = 1000000.0;
    for (auto *vtx : jetVertices)
    {
        const float* vtxPos = vtx->getPosition();

        double dist_x_sq = (vtxPos[0] - primVtxPos[0]) * (vtxPos[0] - primVtxPos[0]);
        double dist_y_sq = (vtxPos[1] - primVtxPos[1]) * (vtxPos[1] - primVtxPos[1]);
        double dist_z_sq = (vtxPos[2] - primVtxPos[2]) * (vtxPos[2] - primVtxPos[2]);

        double dist = sqrt(dist_x_sq + dist_y_sq + dist_z_sq);

        if (dist < minDist)
        {
            minDist = dist;
            closestVtx = vtx;
        }
    }
    return closestVtx;
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