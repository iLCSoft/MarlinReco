#include "Math/Vector3D.h"
#include <ROOT/RVec.hxx>
#include <vector>
#include <algorithm>
#include <map>
#include <random>

using namespace ROOT::Math;
using namespace ROOT::VecOps;
using namespace std;

// in mm/ns
#define SPEED_OF_LIGHT 299.792458
// const float rInner = 1804.8;

struct RandomEngines{
    RandomEngines(){
        //IMPORTANT: I need to define several engines, if I use multiple threads! I hope there is no > 30 threads
        for (int i=0; i<30; i++) engines.push_back( std::mt19937( rd() ) );
    }
    std::vector< std::mt19937 > engines; // or std::default_random_engine;
    std::random_device rd;
};

RandomEngines engines;

std::normal_distribution<float> gaus50(0., 0.05);
std::normal_distribution<float> gaus100(0., 0.1);

RVec <XYZVector> hitPos(const RVec<float>& x, const RVec<float>& y, const RVec<float>& z){
    auto constructHit = [](float x, float y, float z) { return XYZVector(x, y, z); };
    return Map(x, y, z, constructHit);
}


RVec <float> dToImpact(const RVec<XYZVector>& hits, const XYZVector& rImpact){
    auto getDistance = [&](const XYZVector& hit) { return (hit - rImpact).R(); };
    return Map(hits, getDistance);
}

RVec <float> getTimeAtSurface(const RVec<float>& tHit, const RVec<float>& dToImpact){
    auto correctForDistance = [](float tHit, float dToImpact) { return tHit - dToImpact/SPEED_OF_LIGHT; };
    return Map(tHit, dToImpact, correctForDistance);
}



RVec <float> dToLine(const RVec <XYZVector>& hits, const XYZVector& p0, const XYZVector& p){
    RVec <float> distance;
    for (const auto& hit:hits){
        float d = (hit - p0).Cross(p.Unit()).R();
        distance.push_back(d);
    }
    return distance;
}

RVec <float> smear(unsigned int slot, const RVec <float>& times, std::normal_distribution<float>& gaus){
    auto smearWithGaussian = [&](float time){ return time + gaus( engines.engines[slot] ); };
    return Map(times, smearWithGaussian);
}

RVec <bool> selectHits_OLD(const RVec<float>& dToLine, const RVec<int>& layer_hit, bool only_closest=true, int n_layers=10, float cyl_cut=9999.){
    int nHits = dToLine.size();
    RVec <bool> selected(nHits);

    if(only_closest){
        std::map<int, std::vector< std::pair<int, float> > > layer2hit;
        for (int i=0; i < nHits; ++i){
            if( dToLine[i] < cyl_cut ){
                layer2hit[ layer_hit[i] ].push_back( { i, dToLine[i] } );
            }
        }

        for (int layer=0; layer<n_layers; ++layer){
            if ( layer2hit[layer].empty() ) continue;
            int idx = (*std::min_element( layer2hit[layer].begin(), layer2hit[layer].end(), [](const auto& l, const auto& r) { return l.second < r.second; })).first;
            selected[idx] = true;
        }
    }
    else{
        for (int i=0; i < nHits; ++i){
            selected[i] = dToLine[i] < cyl_cut && layer_hit[i] < n_layers;
        }
    }
    return selected;
}

RVec <bool> selectFrankHits(const RVec<float>& dToLine, const RVec<int>& layer_hit){
    int nHits = dToLine.size();
    RVec <bool> selected(nHits);

    int maxLayer = *std::max_element( layer_hit.begin(), layer_hit.end() );

    std::map<int, std::vector< std::pair<int, float> > > layer2hit;
    for (int i=0; i < nHits; ++i){
        layer2hit[ layer_hit[i] ].push_back( { i, dToLine[i] } );
    }
    for (int layer=0; layer<=maxLayer; ++layer){
        if ( layer2hit[layer].empty() ) continue;
        int idx = (*std::min_element( layer2hit[layer].begin(), layer2hit[layer].end(), [](const auto& l, const auto& r) { return l.second < r.second; })).first;
        selected[idx] = true;
    }
    return selected;
}

RVec <bool> selectCylinderHits(const RVec<float>& dPerp, float cutR){
    int nHits = dPerp.size();
    RVec <bool> selected(nHits);
    if (nHits == 0) return selected;

    float minHitR = *std::min_element(dPerp.begin(), dPerp.end());
    float r = max(cutR, minHitR);
    bool foundHit = false;
    while (not foundHit){
        for (int i=0; i < nHits; ++i){
            if ( dPerp[i] < r ){
                selected[i] = true;
                foundHit = true;
            }
        }
        if (foundHit) return selected;
        r += 0.05;
    }
    return selected;
}

RVec <bool> selectMedianHits(const RVec<float>& tSurface, float cutTime=180.){
    RVec<float> times = tSurface;
    size_t n = times.size();
    std::sort(times.begin(), times.end());
    float median = 0.;
    if (n % 2 != 0) median = times[n / 2];
    else median = ( times[(n-1)/2] + times[n/2] )/2.0;

    auto getDtToMedian = [&](float time){ return std::abs(time - median)*1000.; }; // in ps!
    auto distancesToMedian = Map(tSurface, getDtToMedian);

    float dtMin = *std::min_element( distancesToMedian.begin(), distancesToMedian.end() );

    auto applyCut = [&](float dt){ return dt <= cutTime; };

    if ( cutTime < dtMin ) cutTime = dtMin;
    return Map(distancesToMedian, applyCut);
}

float getClosestDtoLine(const RVec<float>& dToLine){
    return *std::min_element(dToLine.begin(), dToLine.end());
}


float tofClosest(const RVec<float>& tHit, const RVec<float>& dToImpact){
    int min_idx = std::min_element(dToImpact.begin(), dToImpact.end()) - dToImpact.begin();
    return tHit[min_idx] - dToImpact[min_idx]/SPEED_OF_LIGHT;
}

float getAverage(const RVec<float>& v){
    if( v.empty() ) return 0;
    auto const count = static_cast<float>( v.size() );
    return std::reduce(v.begin(), v.end()) / count;
}

int getNHitsInLayers(const RVec<int>& layer_hit, int n_layers=10){
    return std::count_if(layer_hit.begin(), layer_hit.end(), [&](int layer){return layer < n_layers;});
}

