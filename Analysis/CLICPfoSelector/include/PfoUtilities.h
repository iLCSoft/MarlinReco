#ifndef PFOUTILITIES_H
#define PFOUTILITIES_H 1

#include <EVENT/ReconstructedParticle.h>
#include <EVENT/LCObject.h>
#include "TrackHitPair.h"

#include <vector>


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
/*
  float Time(const Track* pTrack, float &tof){

    std::cout << "Using PfoUtil::TimeAtEcal" << std::endl;
    float time = 0.0;
    std::cout << time << std::endl;
    return time;
  }
*/
  void GetClusterTimes(const Cluster* cluster, float &meanTime, int &nCaloHitsUsed, float &meanTimeEcal, 
                       int &nEcal, float &meanTimeHcalEndcap, int &nHcalEnd);

}//namespace

#endif // PFOUTILITIES_H
