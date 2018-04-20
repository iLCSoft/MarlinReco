#include "PfoUtilities.h"

#include "marlin/VerbosityLevels.h"
#include <CalorimeterHitType.h>

#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>

#include <IMPL/LCCollectionVec.h>

#include <marlinutil/GeometryUtil.h>
#include <CalorimeterHitType.h>
#include <set>
#include <algorithm>

float PfoUtil::TimeAtEcal(const Track* pTrack, float &tof){

  streamlog_out( MESSAGE ) << " Using PfoUtil::TimeAtEcal " << std::endl;
  float bField = MarlinUtil::getBzAtOrigin();
//  streamlog_out( MESSAGE ) << " bfield " << bField << std::endl;
  const TrackState *pTrackState = pTrack->getTrackState(TrackState::AtCalorimeter);
//  streamlog_out( MESSAGE ) << " track state " << pTrackState->getReferencePoint() << std::endl;

  const float* locationAtECal = pTrackState->getReferencePoint();
  HelixClass helix;
  helix.Initialize_Canonical(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(), pTrack->getTanLambda(), bField);

  tof = sqrt( locationAtECal[0]*locationAtECal[0]+
              locationAtECal[1]*locationAtECal[1]+
              locationAtECal[2]*locationAtECal[2])/300;
  
  float distance[3] = {0.0, 0.0, 0.0};
  float minTime = helix.getDistanceToPoint(locationAtECal, distance);
  
  const float px = helix.getMomentum()[0];
  const float py = helix.getMomentum()[1];
  const float pz = helix.getMomentum()[2];
  const float E = sqrt(px*px+py*py+pz*pz+0.139*0.139);
  minTime = minTime/300*E-tof;

  streamlog_out( MESSAGE ) << " minTime: " << minTime << std::endl;
  return minTime;

}

void PfoUtil::GetClusterTimes(const Cluster* cluster, float &meanTime, int &nCaloHitsUsed, float &meanTimeEcal, 
                              int &nEcal, float &meanTimeHcalEndcap, int &nHcalEnd, bool correctHitTimesForTimeOfFlight){

  streamlog_out( MESSAGE ) << " Using PfoUtil::GetClusterTimes " << std::endl;

  float sumTimeEnergy(0.f);
  float sumEnergy(0.f);
  float sumEnergyEcal(0.f);
  float sumTimeEnergyEcal(0.f);
  float sumEnergyHcalEndcap(0.f);
  float sumTimeEnergyHcalEndcap(0.f);
  meanTime     = std::numeric_limits<float>::max();
  meanTimeEcal = std::numeric_limits<float>::max();
  meanTimeHcalEndcap = std::numeric_limits<float>::max();
  nEcal    = 0;
  nHcalEnd = 0;
  nCaloHitsUsed = 0;

  CalorimeterHitVec hits = cluster->getCalorimeterHits();
  std::vector<float> hittimes;
  std::vector<float> tofCorrections;
  std::vector<float> deltaTimes;

  for(unsigned int ihit=0;ihit<hits.size();++ihit){

    // optionally correct hit times for straight line tof (may have already been done in another processor)
    if(correctHitTimesForTimeOfFlight){
      const float x = hits[ihit]->getPosition()[0];
      const float y = hits[ihit]->getPosition()[1];
      const float z = hits[ihit]->getPosition()[2];
      const float r = sqrt(x*x+y*y+z*z);
      const float tof = r/300.;
      tofCorrections.push_back(tof);
      hittimes.push_back(hits[ihit]->getTime()-tof);
    }else{
      hittimes.push_back(hits[ihit]->getTime());
    }
  }

  std::sort(hittimes.begin(),hittimes.end());
  
  int iMedian = static_cast<int>(hits.size()/2.);
  float medianTime = hittimes[iMedian];
  streamlog_out( MESSAGE ) << " Median time : " << medianTime << std::endl;

  for(unsigned int ihit=0;ihit<hits.size();++ihit)deltaTimes.push_back( fabs(hittimes[ihit]-medianTime)); 
  std::sort(deltaTimes.begin(),deltaTimes.end());
  
  unsigned ihit90 = 0;

  if (hits.size() > 1) {
    ihit90 = static_cast<int>((hits.size()*9)/10.);
    if(ihit90>=hits.size()-1)ihit90=hits.size()-2;
  } else {
    ihit90 = 0;
  }
  
  streamlog_out( MESSAGE ) << " hits " << hits.size() << " hit 90 = " << ihit90 << std::endl;    
  float deltaMedian = deltaTimes[ihit90]+0.1;
  streamlog_out( MESSAGE ) << " deltaCut : " << deltaMedian << std::endl;

  for(unsigned int ihit=0;ihit<hits.size();++ihit){
    CalorimeterHit *hit = hits[ihit];
    float hitTime  = hits[ihit]->getTime();
    if(correctHitTimesForTimeOfFlight)hitTime -= tofCorrections[ihit];

    if( (hitTime - medianTime) < deltaMedian){
      sumEnergy += hit->getEnergy();
      sumTimeEnergy += hit->getEnergy()*hitTime;
      nCaloHitsUsed++;
      streamlog_out( MESSAGE ) << " Using : " << hit->getEnergy() << " : " << hit->getTime() << std::endl;

      CHT ch = hit->getType();
      if(ch.is(CHT::ecal)){
        nEcal++;
        sumEnergyEcal += hit->getEnergy();
        sumTimeEnergyEcal += hit->getEnergy()*hitTime;
      }else{
//      float z = hit->getPosition()[2]; 
        if(!ch.is(CHT::barrel)){
          nHcalEnd++;
          sumEnergyHcalEndcap += hit->getEnergy();
          sumTimeEnergyHcalEndcap += hit->getEnergy()*hitTime;
        }
      }
    }else{
      streamlog_out( MESSAGE ) << " notus : " << hit->getEnergy() << " : " << hit->getTime() << std::endl;
    }
  }
  
  if (sumEnergy > 0.f)meanTime = sumTimeEnergy/sumEnergy;
  if (sumEnergyEcal > 0.f)meanTimeEcal = sumTimeEnergyEcal/sumEnergyEcal;
  if (sumEnergyHcalEndcap > 0.f)meanTimeHcalEndcap = sumTimeEnergyHcalEndcap/sumEnergyHcalEndcap;

  streamlog_out( MESSAGE ) << sumEnergy << " " << sumEnergyEcal << " " << nEcal << std::endl;

  return;


}  
