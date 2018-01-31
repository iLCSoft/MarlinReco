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

/***************************************************
 * Implementation of base class LayerResolverBase
 * *********************************************** */

LayerResolverBase::LayerResolverBase(const int _detTypeFlag,
    const std::string _collectionName, const std::string _detectorName, double _sensThickCheatVal) :
  ThicknessSensitive( _sensThickCheatVal>0. ? &LayerResolverBase::SensitiveThicknessCheat : NULL ),
  sensThickCheatVal(_sensThickCheatVal),
  detTypeFlag(_detTypeFlag),
  collectionName(_collectionName),
  detectorName(_detectorName),
  collection(NULL),
  decoder(NULL)
{}

LayerResolverBase::~LayerResolverBase() {
  if(decoder) delete decoder;
  decoder=NULL;
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

  if (decoder) delete decoder;
  decoder = NULL;

  streamlog_out(DEBUG5) << "LayerResolver::SetCollection: Looking for collection "
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

  streamlog_out(DEBUG5) << "Found collection of type \'" << collection->getTypeName()
                       << "\' with encoding \'"
                       << collection->getParameters().getStringVal( lcio::LCIO::CellIDEncoding ) << "\'\n";
  return 0;
}



/******************************************
 * Implementation of class LayerResolver
 * ************************************** */

template <class T>
LayerResolver<T>::LayerResolver(const int _detTypeFlag,
          T *_layering,
          const std::string _collectionName,
          const std::string _detectorName,
          double _sensThickCheatVal) :
  LayerResolverBase(_detTypeFlag, _collectionName, _detectorName, _sensThickCheatVal),
  layering(_layering)
{
  if(_sensThickCheatVal < 0) ThicknessSensitive = (LayerResolverFn)(&LayerResolver<T>::SensitiveThicknessRead);
}


template <class T>
double LayerResolver<T>::SensitiveThicknessRead(int nLayer) const {
  return layering->layers.at(nLayer).thicknessSensitive;
}



/******************************************
 * Implementation of class LayerFinder
 * ************************************** */

LayerFinder::LayerFinder(EVENT::StringVec _collectionNames, dd4hep::Detector& theDetector,
                          FloatVec sensThickCheatVals) :
  layerResolvers()
{

  const std::vector< dd4hep::DetElement > &detElements = theDetector.detectors("tracker", true);

  if(_collectionNames.size() != detElements.size()) {
    streamlog_out(ERROR) << "There are " <<  detElements.size() << " tracker detector elements in "
        "the geometry and " << _collectionNames.size() << " tracker hit collection names have been "
            "set in the parameters.\n";
    exit(0);
  }

  streamlog_out(MESSAGE) << "Tracker has " << detElements.size() << " elements:\n";

  for(unsigned i=0; i<detElements.size(); i++) {

    streamlog_out(MESSAGE) << "Detector element #" << i << " of type \'" << detElements.at(i).type();
    streamlog_out(MESSAGE) << "\', named \'" << detElements.at(i).name() << "\'\n";
    streamlog_out(MESSAGE) << " ... expects collection name " << _collectionNames[i] << "\n";
    streamlog_out(MESSAGE) << " ... has type flags " << detElements.at(i).typeFlag() << "\n";

    int tf = detElements.at(i).typeFlag();

    LayerResolverBase *sr = NULL;

    try {
      streamlog_out(DEBUG) << "Trying ZPlanarData.\n";
      dd4hep::rec::ZPlanarData* layering = detElements.at(i).extension<dd4hep::rec::ZPlanarData>() ;
      sr = new PlaneResolver(tf, layering, _collectionNames[i], detElements.at(i).name(), sensThickCheatVals.at(i));
    }
    catch ( std::exception &e) {
      streamlog_out(DEBUG) << "Caught exception " << e.what() << std::endl;
      try {
            streamlog_out(DEBUG) << "Trying ZDiskPetalsData.\n";
            dd4hep::rec::ZDiskPetalsData *layering = detElements.at(i).extension<dd4hep::rec::ZDiskPetalsData>() ;
            sr = new PetalResolver(tf, layering, _collectionNames[i], detElements.at(i).name(), sensThickCheatVals.at(i));
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
        streamlog_out(ERROR) << "Detector element #" << i << ", named \'" << detElements.at(i).name()
            << "\' has sensitive thickness " << sensThick << " mm in layer " << iLayer
            << " and cheat value has not been set.\n";
        exit(0);
      }

      streamlog_out(MESSAGE) << "\n";
    }

    ResolverMapIter it = layerResolvers.find(detElements.at(i).id());
    if (it == layerResolvers.end()) {
      layerResolvers[detElements.at(i).id()] = sr;
    }
    else {
      streamlog_out(ERROR) << "Detector element \'" << detElements.at(i).name()
          << "\' has ID = " << detElements.at(i).id()
          << " which is already taken by the detector element \'"
          << it->second->GetDetectorName() << "\'! Aborting.\n";
      exit(0);
    }


    if(sensThickCheatVals.at(i) > 0.) {
      streamlog_out(MESSAGE) << " ... Replacing this detector sensitive thicknesses with value "
          << sensThickCheatVals.at(i)/dd4hep::mm << " mm.\n";
    }

  } // Loop over detector elements (i)

}


int LayerFinder::ReadCollections(EVENT::LCEvent *evt) {

  for(auto collit : layerResolvers) {

    collit.second->SetCollection(evt);

  } // End loop over known detectors

  return 0;
}


double LayerFinder::SensitiveThickness(TrackerHitPlane* thit) {

  streamlog_out(DEBUG5) << "LayerFinder::SensitiveThickness() for hit ID = " << thit->getCellID0() << ".\n";

  // Decode system ID where the hit is located.
  int systemID = FindSystem(thit);

  // Check if the system is registered for handling by the processor
  ResolverMapIter resolver = layerResolvers.find(systemID);
  if (resolver == layerResolvers.end()) {
    streamlog_out(WARNING) << "Hit is located in the system with ID " << systemID
        << " which is not handled by the processor.\n";
    return -1.;
  }

  streamlog_out(DEBUG5) << "Hit is located in subdetector \'" << resolver->second->GetDetectorName()
      << "\' with system ID " << systemID << "\n";

  double t = resolver->second->SensitiveThickness(thit);

  streamlog_out(DEBUG5) << " ... Thickness is " << t/dd4hep::mm << " mm.\n";
  return t;
}


void LayerFinder::ReportHandledDetectors() {


  streamlog_out(DEBUG5) << "LayerFinder handles the following detectors:\n";
  for(auto dit : layerResolvers) {
    std::string dettype;
    if      (dit.second->GetDetTypeFlag() & dd4hep::DetType::BARREL) dettype = "BARREL";
    else if (dit.second->GetDetTypeFlag() & dd4hep::DetType::ENDCAP) dettype = "ENDCAP";
    streamlog_out(DEBUG5) << "Detector \'" << dit.second->GetDetectorName() << "\' of type \'" << dettype;
    streamlog_out(DEBUG5) << "\' associated with collection name \'" << dit.second->GetCollectionName() << "\'.\n";
    streamlog_out(DEBUG5) << "Currently looking at collection of type \'" << dit.second->GetCollectionType();
    streamlog_out(DEBUG5) << "\' with encoding \'" << dit.second->GetCollectionEncoding() << "\'\n";
  }
}


int LayerFinder::FindSystem(TrackerHitPlane* thit) const {

  for( auto cit : layerResolvers) {

    for(int i=0; i<cit.second->GetNumberOfHits(); i++) {

      if( thit == cit.second->GetHit(i) ) {

        streamlog_out(DEBUG5) << " ... Found matching hit " << cit.second->GetHit(i)->getCellID0()
                                 << " with energy " << cit.second->GetHit(i)->getEDep()
                                 << " in collection " << cit.second->GetCollectionName() << ".\n";

        return cit.second->DecodeSystem(thit);
      }
    } // Loop over hits
  } // Loop over collections

  return -1;

}
