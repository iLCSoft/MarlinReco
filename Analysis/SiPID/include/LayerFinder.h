/*
 * LayerFinder.h
 *
 *  Created on: Nov 11, 2016
 *      Author: slukic
 */

#ifndef INCLUDE_LAYERFINDER_H_
#define INCLUDE_LAYERFINDER_H_

#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <UTIL/CellIDDecoder.h>
#include "LCIOTypes.h"
#include "streamlog/streamlog.h"

#include <DD4hep/Detector.h>
#include "DDRec/DetectorData.h"

typedef std::map<EVENT::LCCollection*, lcio::CellIDDecoder<lcio::TrackerHitPlane>*> CollectionMap;

/* Utility class to connect the dots:
 * - Detector element type flagword
 * - Struct extension (for sensitive thickness)
 * - Hit collection
 * - decoder (layer number decoder string)
 */
class LayerResolverBase {

public:
  LayerResolverBase() = delete ;
  LayerResolverBase(const LayerResolverBase &) = delete ;
  LayerResolverBase(const int _detTypeFlag,
                    const std::string _collectionName,
                    const std::string _detectorName,
                    double _sensThickCheatVal=-1.);

  virtual ~LayerResolverBase() ;

  const LayerResolverBase& operator=(const LayerResolverBase&) = delete;

  /* The detector type flag helps to distinguish pointers
   * to plane from petal resolver objects at runtime.
   */
  int GetDetTypeFlag() const { return detTypeFlag; }

  int SetCollection(EVENT::LCEvent *) ;
  std::string GetCollectionName() const { return collectionName; }
  std::string GetDetectorName() const { return detectorName; }
  std::string GetCollectionType() const;
  std::string GetCollectionEncoding() const;

  // Event-by-event:
  int   GetNumberOfHits() const;
  lcio::TrackerHitPlane* GetHit(int i) const;

  virtual unsigned GetNumberOfLayers() const = 0;

  lcio::CellIDDecoder<lcio::TrackerHitPlane>* GetDecoder() const { return decoder; }
  bool HasCollection() const { return static_cast<bool>(decoder); }
  int DecodeLayer(lcio::TrackerHitPlane* thit) const { return (*decoder)(thit)["layer"]; }
  int DecodeSystem(lcio::TrackerHitPlane* thit) const { return (*decoder)(thit)["system"]; }


  virtual double SensitiveThickness(int nLayer) const { return (this->*ThicknessSensitive)(nLayer); }
  virtual double SensitiveThickness(lcio::TrackerHitPlane* thit) const {
    return (this->*ThicknessSensitive)(DecodeLayer(thit));
  }

  typedef  double (LayerResolverBase::*LayerResolverFn)(int) const;
  bool CheatsSensThickness() const { return (ThicknessSensitive==&LayerResolverBase::SensitiveThicknessCheat) ; }

protected:

  LayerResolverFn ThicknessSensitive;
  virtual double SensitiveThicknessRead(int nLayer) const = 0;
  double SensitiveThicknessCheat(int) const { return sensThickCheatVal; };
  // Thickness value used for cheating
  const double sensThickCheatVal;


  // Constant in run:
  int detTypeFlag;
  std::string collectionName;
  std::string detectorName;
  // Event-to-event
  EVENT::LCCollection *collection;
  lcio::CellIDDecoder<lcio::TrackerHitPlane>* decoder;

};


template <class T> class LayerResolver : public LayerResolverBase {

public:
  LayerResolver() = delete ;
  LayerResolver(const LayerResolver<T> &lt) = delete ;
  LayerResolver(const int _detTypeFlag, T *,
                const std::string _collectionName,
                const std::string _detectorName,
                double _sensThickCheatVal=-1.);

  ~LayerResolver() {};

  unsigned GetNumberOfLayers() const { return layering->layers.size(); }

  const LayerResolver& operator=(const LayerResolver<T>&) = delete ;
  const T *GetLayering() const {return layering;};

protected:
  double SensitiveThicknessRead(int nLayer) const ;
  const T *layering;

};

typedef LayerResolver<dd4hep::rec::ZDiskPetalsData> PetalResolver;
typedef LayerResolver<dd4hep::rec::ZPlanarData> PlaneResolver;



/* LayerFinder class parses tracker hit collections in the event
 * finds where the hit belogs,
 * and uses a vector of LayerResolver objects to decode the
 * layer number and get the sensitive thickness
 */

class LayerFinder {

public:
  LayerFinder() = delete;
  // Constructor with the vector of collection names that the finder
  // will use when looking for the collections.
  LayerFinder(EVENT::StringVec _collectionNames, dd4hep::Detector&, lcio::FloatVec sensThickCheatVals);
  LayerFinder( const LayerFinder& ) = delete;

  LayerFinder& operator = ( const LayerFinder& ) = delete;

  /* Reads the decoders of whichever collections are found in the event
   * among those stored in knownDetectors. Returns zero on success and -1 if
   * some of the decoders could not be read.
   */
  int ReadCollections(EVENT::LCEvent *);

  /* Returns the sensitive thickness of the layer where the hit was recorded
   * Also returns the detector type flag in the second argument
   */
  double SensitiveThickness(lcio::TrackerHitPlane*);

  // Sensitive thickness of layer
  void ReportHandledDetectors();

  typedef std::map<int, LayerResolverBase*> ResolverMap;
  typedef ResolverMap::iterator ResolverMapIter;

//  CellIDDecoder<TrackerHitPlane>* GetDecoder() const { return decoder; }

protected:
  ResolverMap layerResolvers{};

  // Search for the TrackerHit collection that contains the hit and
  // use the collection decoder string to decode the system from CellID
  int FindSystem(lcio::TrackerHitPlane* thit) const ;

};



#endif /* INCLUDE_LAYERFINDER_H_ */
