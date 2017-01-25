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

using namespace DD4hep;

typedef std::map<EVENT::LCCollection*, CellIDDecoder<TrackerHitPlane>*> CollectionMap;

/* Utility class to connect the dots:
 * - Detector element type flagword
 * - Struct extension (for sensitive thickness)
 * - Hit collection
 * - decoder (layer number decoder string)
 */
class LayerResolver {

public:
  LayerResolver(const LayerResolver &);
  LayerResolver(const int _detTypeFlag, const std::string _collectionName, double _sensThickCheatVal=-1.);

  virtual ~LayerResolver() ;

  /* The detector type flag helps to distinguish pointers
   * to plane from petal resolver objects at runtime.
   */
  int GetDetTypeFlag() const { return detTypeFlag; }

  int SetCollection(EVENT::LCEvent *) ;
  std::string GetCollectionName() const { return collectionName; }
  std::string GetCollectionType() const;
  std::string GetCollectionEncoding() const;

  // Event-by-event:
  int   GetNumberOfHits() const;
  TrackerHitPlane* GetHit(int i) const;

  virtual int GetNumberOfLayers() const = 0;

  int Decode(TrackerHitPlane* thit) const { return (*decoder)(thit)["layer"];}

  const LayerResolver& operator=(const LayerResolver&);

  double SensitiveThickness(int nLayer) const { return (this->*ThicknessSensitive)(nLayer); }
  typedef  double (LayerResolver::*LayerResolverFn)(int) const;
  bool CheatsSensThickness() const { return (ThicknessSensitive==&LayerResolver::SensitiveThicknessCheat) ; }

protected:

  LayerResolverFn ThicknessSensitive;
  virtual double SensitiveThicknessRead(int nLayer) const = 0;
  double SensitiveThicknessCheat(int nLayer) const { return sensThickCheatVal; };
  // Thickness value used for cheating
  const double sensThickCheatVal;


  // Constant in run:
  int detTypeFlag;
  std::string collectionName;
  // Event-to-event
  EVENT::LCCollection *collection;
  CellIDDecoder<TrackerHitPlane>* decoder;

  LayerResolver();

};


class PetalResolver : public LayerResolver{
public:
  PetalResolver(const PetalResolver &lt);
  PetalResolver(const int _detTypeFlag,
            DDRec::ZDiskPetalsData *,
            const std::string _collectionName,
            double _sensThickCheatVal=-1.);

  ~PetalResolver() {};

  int GetNumberOfLayers() const { return layering->layers.size(); }

  const PetalResolver& operator=(const PetalResolver&);
  const DDRec::ZDiskPetalsData *GetLayering() const {return layering;};

protected:
  double SensitiveThicknessRead(int nLayer) const ;
  const DDRec::ZDiskPetalsData *layering;

  PetalResolver();

};


class PlaneResolver : public LayerResolver{
public:
  PlaneResolver(const PlaneResolver &lt);
  PlaneResolver(const int _detTypeFlag,
            DDRec::ZPlanarData *,
            const std::string _collectionName,
            double _sensThickCheatVal=-1.);

  ~PlaneResolver() {};

  int GetNumberOfLayers() const { return layering->layers.size(); }

  const PlaneResolver& operator=(const PlaneResolver&);
  const DDRec::ZPlanarData *GetLayering() const {return layering;};


protected:
  double SensitiveThicknessRead(int nLayer) const ;
  const DDRec::ZPlanarData *layering;

  PlaneResolver();

};



/* LayerFinder class parses tracker hit collections in the event
 * finds where the hit belogs,
 * and uses a vector of LayerResolver objects to decode the
 * layer number and get the sensitive thickness
 */

class LayerFinder {

public:
  // Constructor with the vector of collection names that the finder
  // will use when looking for the collection.
  LayerFinder(EVENT::StringVec _collectionNames, Geometry::LCDD&, FloatVec sensThickCheatVals, int elementMask);


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

protected:
  std::vector<LayerResolver*> knownDetectors;

  // Default constructor. Not very useful.
  LayerFinder();
};



#endif /* INCLUDE_LAYERFINDER_H_ */
