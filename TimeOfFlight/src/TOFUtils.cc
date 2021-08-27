#include "TOFUtils.h"

#include <cmath>
#include <algorithm>

#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/Operators.h>

#include "marlinutil/CalorimeterHitType.h"


namespace TOFUtils{

using namespace lcio ;


float computeFlightLength( EVENT::Track* trk){

    const TrackState* tsIp = trk->getTrackState( TrackState::AtIP ) ;
    const TrackState* tscalo = trk->getTrackState( TrackState::AtCalorimeter ) ;

    float phiIp = tsIp->getPhi() ;
    float phiCalo = tscalo->getPhi() ;

    float omega = tsIp->getOmega()  ;
    float tanL = tsIp->getTanLambda() ;

    float dPhi = std::abs(phiIp - phiCalo);

    // if dPhi > PI assume phi is flipped
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;

    return dPhi/std::abs(omega) * std::sqrt( 1. + tanL * tanL ) ;
}

float computeFlightLength(const TrackState* ts0,  const TrackState* ts1 ){
    float phi1 = ts1->getPhi() ;
    float phi0 = ts0->getPhi() ;

    float omega = ts0->getOmega()  ;
    float tanL = ts0->getTanLambda() ;

    float dPhi = std::abs(phi0 - phi1);

    // if dPhi > PI assume phi is flipped
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;

    return dPhi/std::abs(omega) * std::sqrt( 1. + tanL * tanL ) ;
}


  int layer( EVENT::CalorimeterHit* h ) {

    return CHT( h->getType() ).layer() ;
  }


  /// helper function to check if this is an Ecal hit
  bool isEcal( EVENT::CalorimeterHit* h ) {

    return CHT( h->getType() ).caloID() == CHT::ecal ;
  }


  std::string caloTypeStr(  EVENT::CalorimeterHit* h ) {

    std::stringstream s ;
    s << CHT( h->getType() ) ;
    return s.str() ;
  }


  std::string CaloHitData::toString(){

    std::stringstream s ;
    s << "  l= " << layer ;
    s << ", st= " << smearedTime ;
    s << ", tr= " << timeResolution  ;
    s << ", dIP=" << distanceFromIP  ;
    s << ", dRP= " << distanceFromReferencePoint  ;
    s << ", dSL= " << distancefromStraightline  ;
//    EVENT::CalorimeterHit* lcioHit = nullptr ;

    return s.str() ;
  }


  float computeDistanceFromLine( EVENT::CalorimeterHit* h, const dd4hep::rec::Vector3D& point,
				 const dd4hep::rec::Vector3D& unitDir) {

    dd4hep::rec::Vector3D pos( h->getPosition()[0],
			       h->getPosition()[1],
			       h->getPosition()[2] ) ;

    dd4hep::rec::Vector3D diff = pos - point ;


    return diff.cross( unitDir ).r() ;

  }


  CaloHitDataVec findHitsClosestToLine( const CaloHitLayerMap& layerMap ){

    CaloHitDataVec hitVec ;

    for( auto m : layerMap ){

      const CaloHitDataVec& chv = m.second ;

      CaloHitData* closestHit =
	*std::min_element( chv.begin() , chv.end () ,
			   [](CaloHitData* c0, CaloHitData* c1 ){ return c0->distancefromStraightline < c1->distancefromStraightline  ; }
	  ) ;

      hitVec.push_back( closestHit ) ;
    }

    return hitVec ;
  }


  std::pair<float,float> computeTOFEstimator( const CaloHitDataVec& chv ){

    const static float c_mm_per_ns = 299.792458 ;

    int   nHit = 0 ;
    float mean = 0. ;
    float meansq = 0. ;

    for( auto ch : chv ){
      float t = ch->smearedTime - ch->distanceFromReferencePoint / c_mm_per_ns ;
      mean   += t   ;
      ++nHit ;
    }
    mean /= nHit ;

    for( auto ch : chv ){
      float t = ch->smearedTime - ch->distanceFromReferencePoint / c_mm_per_ns ;
      meansq += ( t - mean ) * ( t - mean )   ;
    }

    return std::make_pair( mean, std::sqrt(meansq/nHit ) );
  }


} //namespace
