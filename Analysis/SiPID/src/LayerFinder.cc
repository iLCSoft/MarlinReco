/*
 * LayerFinder.cc
 *
 *  Created on: Nov 11, 2016
 *      Author: slukic
 */

#include <DD4hep/Detector.h>
#include <DD4hep/DetType.h>
#include <LayerFinder.h>
#include "DD4hep/DD4hepUnits.h"

/******************************************
 * Implementation of base class LayerResolver
 * ************************************** */

LayerResolverBase::LayerResolverBase(const LayerResolverBase &lt) :
  ThicknessSensitive( lt.CheatsSensThickness() ? &LayerResolverBase::SensitiveThicknessCheat : NULL ),
  sensThickCheatVal(lt.SensitiveThicknessCheat(0)),
  detTypeFlag(lt.GetDetTypeFlag()),
  collectionName(lt.GetCollectionName()),
  collection(lt.collection),
  decoder(lt.decoder)
{}

LayerResolverBase::LayerResolverBase(const int _detTypeFlag,
              const std::string _collectionName, double _sensThickCheatVal) :
  ThicknessSensitive( _sensThickCheatVal>0. ? &LayerResolverBase::SensitiveThicknessCheat : NULL ),
  sensThickCheatVal(_sensThickCheatVal),
  detTypeFlag(_detTypeFlag),
  collectionName(_collectionName),
  collection(NULL),
  decoder(NULL)
{}

LayerResolverBase::LayerResolverBase() :
  ThicknessSensitive(NULL),
  sensThickCheatVal(0),
  detTypeFlag(0),
  collectionName(""),
  collection(NULL),
  decoder(NULL)
{}


LayerResolverBase::~LayerResolverBase() {
  if(decoder) delete decoder;
  decoder=NULL;
}

const LayerResolverBase& LayerResolverBase::operator =(const LayerResolverBase& lr) {
  this->detTypeFlag = lr.detTypeFlag;
  this->collectionName = lr.collectionName;
  this->decoder = lr.decoder;
  return *this;
}

std::string LayerResolverBase::GetCollectionType() const {
  if (collection) return collection->getTypeName();
  else return "NONE";
}

std::string LayerResolverBase::GetCollectionEncoding() const {
  if (collection)
    return collection->getParameters().getStringVal( lcio::LCIO::CellIDEncoding );
  else return "EMPTY";
}

// Event-by-event:
int LayerResolverBase::GetNumberOfHits() const {
  if (collection) return collection->getNumberOfElements();
  else return -1;
}

TrackerHitPlane* LayerResolverBase::GetHit(int i) const {
  if (collection) return dynamic_cast<TrackerHitPlane*>(collection->getElementAt(i));
  else return NULL;
}

int LayerResolverBase::SetCollection(EVENT::LCEvent *evt) {

  delete decoder;
  decoder = NULL;

  streamlog_out(DEBUG) << "LayerResolver::SetCollection: Looking for collection "
      << collectionName << ".\n";
  try {
    collection = evt->getCollection(collectionName);
  }
  catch(EVENT::DataNotAvailableException &dataex) {
    streamlog_out(MESSAGE) << "Collection " << collectionName
        << " not found in event #" << evt->getEventNumber() << ".\n";
    collection = NULL;
    return -1;
  }

  decoder = new CellIDDecoder<TrackerHitPlane>(collection);

  streamlog_out(DEBUG) << "Found collection of type " << collection->getTypeName()
                       << " with encoding "
                       << collection->getParameters().getStringVal( lcio::LCIO::CellIDEncoding ) << "\n";
  return 0;
}



/******************************************
 * Implementation of class SmartResolver
 * ************************************** */

template <class T>
LayerResolver<T>::LayerResolver(const LayerResolver<T> &lt) :
  LayerResolverBase(lt),
  layering(lt.GetLayering())
{
  ThicknessSensitive = (lt.CheatsSensThickness() ?
      &LayerResolver::SensitiveThicknessCheat :
      &LayerResolver::SensitiveThicknessRead );
}

template <class T>
LayerResolver<T>::LayerResolver(const int _detTypeFlag,
          T *_layering,
          const std::string _collectionName,
          double _sensThickCheatVal) :
  LayerResolverBase(_detTypeFlag, _collectionName, _sensThickCheatVal),
  layering(_layering)
{
  if(_sensThickCheatVal < 0) ThicknessSensitive = (LayerResolverFn)(&LayerResolver::SensitiveThicknessRead);
}

template <class T>
LayerResolver<T>::LayerResolver() :
  LayerResolverBase(),
  layering(NULL)
{}

template <class T>
const LayerResolver<T>& LayerResolver<T>::operator=(const LayerResolver<T>& pr)
  {
  this->sensThickCheatVal = pr.sensThickCheatVal;
  this->detTypeFlag = pr.detTypeFlag;
  this->collectionName = pr.collectionName;
  this->collection = pr.collection;
  this->decoder = pr.decoder;
  this->layering = pr.layering;

  this->ThicknessSensitive = (pr.CheatsSensThickness() ?
      &LayerResolverBase::SensitiveThicknessCheat :
      &LayerResolverBase::SensitiveThicknessRead );

  return *this;
}

template <class T>
double LayerResolver<T>::SensitiveThicknessRead(int nLayer) const {
  return layering->layers.at(nLayer).thicknessSensitive;
}



/******************************************
 * Implementation of class LayerFinder
 * ************************************** */

LayerFinder::LayerFinder(EVENT::StringVec _collectionNames, dd4hep::Detector& theDetector,
                          FloatVec sensThickCheatVals, int elementMask) :
  knownDetectors()
{

  const std::vector< dd4hep::DetElement > &detElements = theDetector.detectors("tracker", true);

  if(_collectionNames.size() != detElements.size()) {
    streamlog_out(WARNING) << "There are " <<  detElements.size() << " tracker detector elements in "
        "the geometry and " << _collectionNames.size() << " tracker hit collection names have been "
            "set in the parameters.\n";
  }

  streamlog_out(MESSAGE) << "Tracker has " << detElements.size() << " elements:\n";

  for(unsigned i=0; i<detElements.size(); i++) {

    streamlog_out(MESSAGE) << "Detector element #" << i << " of type " << detElements.at(i).type();
    streamlog_out(MESSAGE) << ", named " << detElements.at(i).name() << "\n";
    streamlog_out(MESSAGE) << " ... expects collection name " << _collectionNames[i] << "\n";
    streamlog_out(MESSAGE) << " ... has type flags " << detElements.at(i).typeFlag() << "\n";

    if(! (elementMask & 1 << i) ) {
      streamlog_out(MESSAGE) << "Turned OFF in the element mask (see ElementMask parameter).\n";
      continue;
    }

    int tf = detElements.at(i).typeFlag();

    LayerResolverBase *sr = NULL;

    try {
      streamlog_out(DEBUG) << "Trying ZPlanarData.\n";
      dd4hep::rec::ZPlanarData* layering = detElements.at(i).extension<dd4hep::rec::ZPlanarData>() ;
      sr = new PLANE_RESOLVER(tf, layering, _collectionNames[i], sensThickCheatVals.at(i));
    }
    catch ( std::exception &e) {
      streamlog_out(DEBUG) << "Caught exception " << e.what() << std::endl;
      try {
            streamlog_out(DEBUG) << "Trying ZDiskPetalsData.\n";
            dd4hep::rec::ZDiskPetalsData *layering = detElements.at(i).extension<dd4hep::rec::ZDiskPetalsData>() ;
            sr = new PETAL_RESOLVER(tf, layering, _collectionNames[i], sensThickCheatVals.at(i));
          }
      catch ( std::exception &e1) {
        streamlog_out(DEBUG) << "Caught exception " << e1.what() << std::endl;
      }
    }

    if (!sr) {
      streamlog_out(ERROR) << "Could not find a semiconductor layering extension in the detector element! Aborting.";
      exit(0);
    }

    streamlog_out(MESSAGE) << "Detector element has " << sr->GetNumberOfLayers() << " layers with thicknesses: ";

    for(unsigned iLayer = 0; iLayer < sr->GetNumberOfLayers(); iLayer++) {

      double sensThick = sr->SensitiveThickness(iLayer)/dd4hep::mm;
      streamlog_out(MESSAGE) << sensThick << " mm, ";

      if(sensThick < 1.e-6 && sensThickCheatVals.at(i) < 0.) {
        streamlog_out(ERROR) << "Detector element #" << i << ", named " << detElements.at(i).name()
            << " has sensitive thickness " << sensThick << " mm in layer " << iLayer
            << " and cheat value has not been set.\n";
        exit(0);
      }

      streamlog_out(MESSAGE) << "\n";
    }

    knownDetectors.push_back(sr);


    if(sensThickCheatVals.at(i) > 0.) {
      streamlog_out(MESSAGE) << " ... Replacing this detector sensitive thicknesses with value "
          << sensThickCheatVals.at(i)/dd4hep::mm << " mm.\n";
    }

  } // Loop over detector elements (i)

}



int LayerFinder::ReadCollections(EVENT::LCEvent *evt) {

  int nFound=0;

  streamlog_out(DEBUG5) << "CollectionFinder::ReadCollections. Event #" << evt->getEventNumber()
      << ". Number of collections to look for: " << knownDetectors.size() << "\n";

  for(std::vector<LayerResolverBase*>::iterator collit=knownDetectors.begin(); collit!=knownDetectors.end(); collit++) {

    if ((*collit)->SetCollection(evt) == 0) nFound++;

  } // End loop over known detectors
  return nFound;

}


double LayerFinder::SensitiveThickness(TrackerHitPlane* thit, int &tf) {

  streamlog_out(DEBUG5) << "CollectionFinder::SensitiveThickness() for hit ID = " << thit->getCellID0() << ".\n";

  for(std::vector<LayerResolverBase*>::iterator cit=knownDetectors.begin(); cit!=knownDetectors.end(); cit++) {

    streamlog_out(DEBUG5) << " ... Looking in collection " << (*cit)->GetCollectionName() ;
    streamlog_out(DEBUG5) << " of type "
                         << (*cit)->GetCollectionType() << " with " << (*cit)->GetNumberOfHits() << " hits.\n";

    for(int i=0; i<(*cit)->GetNumberOfHits(); i++) {

      if( thit == (*cit)->GetHit(i) ) {

        streamlog_out(DEBUG5) << " ... Found matching hit " << (*cit)->GetHit(i)->getCellID0()
                                 << " with energy " << (*cit)->GetHit(i)->getEDep()
                                 << " in collection " << (*cit)->GetCollectionName() << ".\n";

        int nLayer = (*cit)->Decode(thit);
        tf = (*cit)->GetDetTypeFlag();
        int nLayers = (*cit)->GetNumberOfLayers();

        streamlog_out(DEBUG5) << " ... Layer number is " << nLayer << ".\n";
        streamlog_out(DEBUG5) << " ... Number of layers is " << nLayers << ".\n";
        if(nLayer >= nLayers) {
          streamlog_out(ERROR) << "CollectionFinder::SensitiveThickness() -- Layer number out of range!\n";
          return 1.;
        }

        double t = (*cit)->SensitiveThickness(nLayer);
        streamlog_out(DEBUG5) << " ... Thickness is " << t/dd4hep::mm << " mm.\n";
        return t;

        // We found it so we can stop searching
        break;
      }
    }
  }
  return -1.;
}


void LayerFinder::ReportKnownDetectors() {


  streamlog_out(DEBUG5) << "LayerFinder knows the following detectors:\n";
  for(std::vector<LayerResolverBase*>::iterator dit=knownDetectors.begin(); dit!=knownDetectors.end(); dit++) {
    std::string dettype;
    if      ((*dit)->GetDetTypeFlag() & dd4hep::DetType::BARREL) dettype = "BARREL";
    else if ((*dit)->GetDetTypeFlag() & dd4hep::DetType::ENDCAP) dettype = "ENDCAP";
    streamlog_out(DEBUG5) << "Detector of type " << dettype;
    streamlog_out(DEBUG5) << " associated with collection name " << (*dit)->GetCollectionName() << "\n";
    streamlog_out(DEBUG5) << "Currently looking at collection of type " << (*dit)->GetCollectionType();
    streamlog_out(DEBUG5) << " with encoding " << (*dit)->GetCollectionEncoding() << "\n";
  }
}

