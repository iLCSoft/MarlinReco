#ifndef TOFUtils_h
#define TOFUtils_h 1

/******************************************************
 * Utilities for cumputing time of flight estimators.
 * 
 * @author F. gaede, DESY, 2018
 * 
 ******************************************************
 */

#include <vector>
#include <map>
#include <memory>

#include "EVENT/CalorimeterHit.h"
#include "EVENT/Track.h"
#include "DDRec/Vector3D.h"

namespace EVENT{
  class Track ;
}



namespace TOFUtils{


  /// handle for calorimeter hit meta data used for TOF estimators
  struct CaloHitData{
    CaloHitData() = delete ;
    CaloHitData(const CaloHitData&) = default ;
    CaloHitData& operator=(const CaloHitData&) = default ;
    CaloHitData(EVENT::CalorimeterHit* h) :  lcioHit(h) {} ;
    EVENT::CalorimeterHit* lcioHit = nullptr ;
    int   layer = 0 ;
    float smearedTime = 0. ;
    float timeResolution = 0. ;
    float distanceFromIP = 0. ;
    float distanceFromReferencePoint = 0. ;
    float distancefromStraightline = 0. ;
    std::string toString() ;
  } ;
  
  /// define a vector of unique pointers to extended calo hits
  typedef std::vector< std::unique_ptr<CaloHitData> > CaloHitUPtrVec ;

  typedef std::vector< CaloHitData* > CaloHitDataVec ;
  
  /// define map type for storing the hits by layer
  typedef std::map< int, CaloHitDataVec > CaloHitLayerMap ;
  

  /* compute the distance of the hit from a straight line
   * defined by l = point + lambda * unitDir
   */
  float computeDistanceFromLine( EVENT::CalorimeterHit* h, const dd4hep::rec::Vector3D& point,
				   const dd4hep::rec::Vector3D& unitDir) ;

  
  /// return vector with hits that have the shortest distancefromStraightline for every layer
  CaloHitDataVec findHitsClosestToLine( const CaloHitLayerMap& layerMap ) ;


  /// compute the flight length of the particle from the IP to the calorimeter
  float computeFlightLength( EVENT::Track* trk) ;

  /// compute the flight length of the particle between the two given track states
  float computeFlightLength(const EVENT::TrackState* ts0,  const EVENT::TrackState* ts1 ) ;

  /// helper function to get the layer of the calo hit
  int layer( EVENT::CalorimeterHit* h ) ; 


  /// helper function to check if this is an Ecal hit 
  bool isEcal( EVENT::CalorimeterHit* h ) ; 


  /// string with calo type information
  std::string caloTypeStr(  EVENT::CalorimeterHit* h ) ;


  /* compute the average TOF from all given hits - correcting for time of 
   * flight from reference point
   */
   std::pair<float,float>  computeTOFEstimator( const CaloHitDataVec& chv ) ;

  
}



#endif
