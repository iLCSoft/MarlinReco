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

typedef std::map<EVENT::LCCollection*, CellIDDecoder<TrackerHitPlane>*> CollectionMap;

/* Utility class to connect the dots:
 * - Detector element type flagword
 * - Struct extension (for sensitive thickness)
 * - Hit collection
 * - decoder (layer number decoder string)
 */
class LayerResolverBase {

public:
  LayerResolverBase(const LayerResolverBase &);
  LayerResolverBase(const int _detTypeFlag,
                    const std::string _collectionName,
                    const std::string _detectorName,
                    double _sensThickCheatVal=-1.);

  virtual ~LayerResolverBase() ;

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
  TrackerHitPlane* GetHit(int i) const;

  virtual unsigned GetNumberOfLayers() const = 0;

  int Decode(TrackerHitPlane* thit, const char *what="layer") const { return (*decoder)(thit)[what];}
  bool HasCollection() const { return static_cast<bool>(decoder); }

  const LayerResolverBase& operator=(const LayerResolverBase&);

  virtual double SensitiveThickness(int nLayer) const { return (this->*ThicknessSensitive)(nLayer); }
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
  CellIDDecoder<TrackerHitPlane>* decoder;

  LayerResolverBase();

};


template <class T> class LayerResolver : public LayerResolverBase {

public:
  LayerResolver(const LayerResolver<T> &lt);
  LayerResolver(const int _detTypeFlag, T *,
                const std::string _collectionName,
                const std::string _detectorName,
                double _sensThickCheatVal=-1.);

  ~LayerResolver() {};

  unsigned GetNumberOfLayers() const { return layering->layers.size(); }

  const LayerResolver& operator=(const LayerResolver<T>&);
  const T *GetLayering() const {return layering;};

protected:
  double SensitiveThicknessRead(int nLayer) const ;
  const T *layering;

  LayerResolver();
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
  // Constructor with the vector of collection names that the finder
  // will use when looking for the collection.
  LayerFinder(EVENT::StringVec _collectionNames, dd4hep::Detector&, FloatVec sensThickCheatVals, int elementMask);


  /* Reads the decoders of whichever collections are found in the event
   * among those stored in knownDetectors. Returns zero on success and -1 if
   * some of the decoders could not be read.
   */
  int ReadCollections(EVENT::LCEvent *);

  /* Returns the sensitive thickness of the layer where the hit was recorded
   * Also returns the detector type flag in the second argument
   */
  double SensitiveThickness(TrackerHitPlane*, int &detTypeFlag);

  EVENT::LCCollection* GetCollection(EVENT::LCObject*);
  CellIDDecoder<TrackerHitPlane>* GetDecoder(EVENT::LCObject*);
  int GetLayer(TrackerHitPlane*);
  // Sensitive thickness of layer
  void ReportKnownDetectors();

  typedef std::map<int, LayerResolverBase*> ResolverMap;
  typedef ResolverMap::iterator ResolverMapIter;

protected:
  ResolverMap knownDetectors{};

  // Default constructor. Not very useful.
  LayerFinder();
};



#endif /* INCLUDE_LAYERFINDER_H_ */
