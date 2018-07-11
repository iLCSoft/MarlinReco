#ifndef PFOUTILITIES_H
#define PFOUTILITIES_H 1

#include <EVENT/ReconstructedParticle.h>
#include <EVENT/LCObject.h>
#include "TrackHitPair.h"

#include "lcio.h"
#include "TrackHitPair.h"
#include "HelixClass.h"
#include <string>
#include <map>
#include <set>
#include <algorithm>

#include <marlinutil/GeometryUtil.h>
#include <vector>

using namespace std;

#define FORMATTED_OUTPUT_TRACK_CLUSTER_full(out, N1, E1,E2,E3,N2,E4,N3,E5,E6,E7) \
    out <<                                                                                      \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E3        <<                                   \
    std::right << std::setw(widthInt)      <<    N2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E4        <<                                   \
    std::right << std::setw(widthInt  )    <<    N3        <<                                   \
    std::right << std::setw(widthFloat)    <<    E5        <<                                   \
    std::right << std::setw(widthFloat)    <<    E6        <<                                   \
    std::right << std::setw(widthFloat)    <<    E7  	   << std::endl

#define FORMATTED_OUTPUT_TRACK_CLUSTER(out, N1, E1,E2,E3,N2,N3) \
    out <<                                                                                      \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E3        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N2        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N3	   << std::endl

#define FORMATTED_OUTPUT_TRACK(out, N1, E1,E2,E3,N2,N3) \
    out <<                                                                                      \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E3        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N2        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N3 	   << std::endl

#define FORMATTED_OUTPUT_CLUSTER(out, N1, E1,E2,E3,N2,N3) \
    out <<                                                                                      \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E3        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N2        <<                                   \
    std::right << std::setw(widthSmallInt) <<    N3 	   << std::endl

#define FORMATTED_OUTPUT_MC(out, N1,E1) \
    out <<                                                                                      \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        << std::endl

namespace PfoUtil{

  typedef std::vector<EVENT::ReconstructedParticle*> PfoList;

  static bool PfoSortFunction(EVENT::ReconstructedParticle* lhs,EVENT::ReconstructedParticle* rhs){

  //  true if lhs goes before

   const float lhs_energy  = lhs->getEnergy();
   const float rhs_energy  = rhs->getEnergy();
   const EVENT::ClusterVec lhs_clusters = lhs->getClusters();
   const EVENT::TrackVec   lhs_tracks   = lhs->getTracks();
   const EVENT::TrackVec   rhs_tracks   = rhs->getTracks();
   const EVENT::ClusterVec rhs_clusters = rhs->getClusters();
   const int lhs_particleType  = abs(lhs->getType());
   const int rhs_particleType  = abs(rhs->getType());

   int lhs_type(0);
   int rhs_type(0);

   if(lhs_clusters.size()==0)lhs_type=1;
   if(rhs_clusters.size()==0)rhs_type=1;
   if(lhs_particleType==22)lhs_type=10;
   if(rhs_particleType==22)rhs_type=10;
   if(lhs_particleType==2112)lhs_type=20;
   if(rhs_particleType==2112)rhs_type=20;

   if(lhs_type==rhs_type)return (lhs_energy>rhs_energy);
   return (lhs_type<rhs_type);

  }

  float TimeAtEcal(const Track* pTrack, float &tof);

  // Calculate cluster times
  void GetClusterTimes(const Cluster* cluster, float &meanTime, int &nCaloHitsUsed, float &meanTimeEcal, 
                       int &nEcal, float &meanTimeHcalEndcap, int &nHcalEnd, bool correctHitTimesForTimeOfFlight);

}//namespace

#endif // PFOUTILITIES_H
