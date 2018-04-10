#include "TOFUtils.h"

#include <cmath>

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
    
    const TrackState* tsIP = trk->getTrackState( TrackState::AtIP ) ; 	
    const TrackState* tscalo = trk->getTrackState( TrackState::AtCalorimeter ) ; 	
    
    float Phicalo = tscalo->getPhi() ;
    float PhiIP = tsIP->getPhi() ;
    
    float Omega = tsIP->getOmega()  ;
    float tanL = tsIP->getTanLambda() ;
    
    float length = (PhiIP-Phicalo)*(1/Omega) * sqrt( 1 + tanL * tanL ) ;

    return length ;
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
 
}
