/*
 * LayerFinder.cc
 *
 *  Created on: Nov 11, 2016
 *      Author: slukic
 */

#include <DD4hep/LCDD.h>
#include <DD4hep/DetType.h>
#include <LayerFinder.h>
#include "DD4hep/DD4hepUnits.h"

/******************************************
 * Implementation of base class LayerResolver
 * ************************************** */

LayerResolver::LayerResolver(const LayerResolver &lt) :
  ThicknessSensitive( lt.CheatsSensThickness() ? &LayerResolver::SensitiveThicknessCheat : NULL ),
  sensThickCheatVal(lt.SensitiveThicknessCheat(0)),
  detTypeFlag(lt.GetDetTypeFlag()),
  collectionName(lt.GetCollectionName()),
  collection(lt.collection),
  decoder(lt.decoder)
{}

LayerResolver::LayerResolver(const int _detTypeFlag,
              const std::string _collectionName, double _sensThickCheatVal) :
  ThicknessSensitive( _sensThickCheatVal>0. ? &LayerResolver::SensitiveThicknessCheat : NULL ),
  sensThickCheatVal(_sensThickCheatVal),
  detTypeFlag(_detTypeFlag),
  collectionName(_collectionName),
  collection(NULL),
  decoder(NULL)
{}

LayerResolver::LayerResolver() :
  ThicknessSensitive(NULL),
  sensThickCheatVal(0),
  detTypeFlag(0),
  collectionName(""),
  collection(NULL),
  decoder(NULL)
{}


LayerResolver::~LayerResolver() {
  if(decoder) delete decoder;
  decoder=NULL;
}

const LayerResolver& LayerResolver::operator =(const LayerResolver& lr) {
  this->detTypeFlag = lr.detTypeFlag;
  this->collectionName = lr.collectionName;
  this->decoder = lr.decoder;
  return *this;
}

std::string LayerResolver::GetCollectionType() const {
  if (collection) return collection->getTypeName();
  else return "NONE";
}

std::string LayerResolver::GetCollectionEncoding() const {
  if (collection)
    return collection->getParameters().getStringVal( lcio::LCIO::CellIDEncoding );
  else return "EMPTY";
}

// Event-by-event:
int LayerResolver::GetNumberOfHits() const {
  if (collection) return collection->getNumberOfElements();
  else return -1;
}

TrackerHitPlane* LayerResolver::GetHit(int i) const {
  if (collection) return dynamic_cast<TrackerHitPlane*>(collection->getElementAt(i));
  else return NULL;
}

int LayerResolver::SetCollection(EVENT::LCEvent *evt) {

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
 * Implementation of class PetalResolver
 * ************************************** */

PetalResolver::PetalResolver(const PetalResolver &lt) :
  LayerResolver(lt),
  layering(lt.GetLayering())
{
  if(ThicknessSensitive == NULL) ThicknessSensitive = (LayerResolverFn)(&PetalResolver::SensitiveThicknessRead);
}

PetalResolver::PetalResolver(const int _detTypeFlag,
          DDRec::ZDiskPetalsData *_layering,
          const std::string _collectionName,
          double _sensThickCheatVal) :
  LayerResolver(_detTypeFlag, _collectionName, _sensThickCheatVal),
  layering(_layering)
{
  if(_sensThickCheatVal < 0) ThicknessSensitive = (LayerResolverFn)(&PetalResolver::SensitiveThicknessRead);
}

PetalResolver::PetalResolver() :
  LayerResolver(),
  layering(NULL)
{}

const PetalResolver& PetalResolver::operator=(const PetalResolver& pr)
  {
  this->detTypeFlag = pr.detTypeFlag;
  this->collectionName = pr.collectionName;
  this->decoder = pr.decoder;
  this->layering = pr.layering;
  return *this;
}

double PetalResolver::SensitiveThicknessRead(int nLayer) const {
  return layering->layers.at(nLayer).thicknessSensitive;
}




/******************************************
 * Implementation of class PlaneResolver
 * ************************************** */

PlaneResolver::PlaneResolver(const PlaneResolver &lt) :
  LayerResolver(lt),
  layering(lt.GetLayering())
{
  if(!lt.CheatsSensThickness()) ThicknessSensitive = (LayerResolverFn)(&PlaneResolver::SensitiveThicknessRead);
}

PlaneResolver::PlaneResolver(const int _detTypeFlag,
          DDRec::ZPlanarData *_layering,
          const std::string _collectionName,
          double _sensThickCheatVal) :
  LayerResolver(_detTypeFlag, _collectionName, _sensThickCheatVal),
  layering(_layering)
{
  if(_sensThickCheatVal < 0) ThicknessSensitive = (LayerResolverFn)(&PlaneResolver::SensitiveThicknessRead);
}

PlaneResolver::PlaneResolver() :
  LayerResolver(),
  layering(NULL)
{}


const PlaneResolver& PlaneResolver::operator=(const PlaneResolver& pr)
  {
  this->detTypeFlag = pr.detTypeFlag;
  this->collectionName = pr.collectionName;
  this->decoder = pr.decoder;
  this->layering = pr.layering;
  return *this;
}

double PlaneResolver::SensitiveThicknessRead(int nLayer) const {
  return layering->layers.at(nLayer).thicknessSensitive;
}


LayerFinder::LayerFinder() :
knownDetectors()
{}


/******************************************
 * Implementation of class LayerFinder
 * ************************************** */

LayerFinder::LayerFinder(EVENT::StringVec _collectionNames, Geometry::LCDD& lcdd,
                          FloatVec sensThickCheatVals, int elementMask) :
  knownDetectors()
{

  const std::vector< Geometry::DetElement > &detElements = lcdd.detectors("tracker", true);

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

    if(! (elementMask & 1 << i) ) {
      streamlog_out(MESSAGE) << "Turned OFF in the element mask (see ElementMask parameter).\n";
      continue;
    }

    int tf = detElements.at(i).typeFlag();
    detElements.at(i).name();

    if (tf & DD4hep::DetType::BARREL) {

      DDRec::ZPlanarData *layering = detElements.at(i).extension<DD4hep::DDRec::ZPlanarData>() ;
      streamlog_out(MESSAGE) << " ... has " << layering->layers.size() << " layers with thicknesses: ";

      for(unsigned iLayer = 0; iLayer < layering->layers.size(); iLayer++) {

        double sensThick = layering->layers.at(i).thicknessSensitive/dd4hep::mm;
        streamlog_out(MESSAGE) << sensThick << " mm, ";

        if(sensThick < 1.e-6 && sensThickCheatVals.at(i) < 0.) {
          streamlog_out(ERROR) << "Detector element #" << i << ", named " << detElements.at(i).name()
              << " has sensitive thickness " << sensThick << " mm in layer " << iLayer
              << " and cheat value has not been set.\n";
          exit(0);
        }
      }

      streamlog_out(MESSAGE) << "\n";
      knownDetectors.push_back(new PlaneResolver(tf, layering, _collectionNames[i], sensThickCheatVals.at(i)));

    }
    else if (tf & DD4hep::DetType::ENDCAP) {

      DDRec::ZDiskPetalsData *layering = detElements.at(i).extension<DD4hep::DDRec::ZDiskPetalsData>() ;
      streamlog_out(MESSAGE) << " ... has " << layering->layers.size() << " layers with thicknesses: ";

      for(unsigned iLayer = 0; iLayer < layering->layers.size(); iLayer++) {

        double sensThick = layering->layers.at(i).thicknessSensitive/dd4hep::mm;
        streamlog_out(MESSAGE) << sensThick << " mm, ";

        if(sensThick < 1.e-6 && sensThickCheatVals.at(i) < 0.) {
          streamlog_out(ERROR) << "Detector element #" << i << ", named " << detElements.at(i).name()
              << " has sensitive thickness " << sensThick << " mm in layer " << iLayer
              << " and cheat value has not been set.\n";
          exit(0);
        }
      }


      streamlog_out(MESSAGE) << "\n";
      knownDetectors.push_back(new PetalResolver(tf, layering, _collectionNames[i], sensThickCheatVals.at(i)));
    }
    else {
      streamlog_out(WARNING) << "Detector element #" << i << " of type " << detElements.at(i).type();
      streamlog_out(WARNING) << ", named " << detElements.at(i).name() << " has an unknown type (barrel/endcap)"
          " and will not be added!\n";
    }

    if(sensThickCheatVals.at(i) > 0.) {
      streamlog_out(MESSAGE) << " ... Replacing this detector sensitive thicknesses with value "
          << sensThickCheatVals.at(i)/dd4hep::mm << " mm.\n";
    }

  } // Loop over detector elements (i)

}



int LayerFinder::ReadCollections(EVENT::LCEvent *evt) {

  int nFound=0;

  streamlog_out(DEBUG) << "CollectionFinder::ReadCollections. Event #" << evt->getEventNumber()
      << ". Number of collections to look for: " << knownDetectors.size() << "\n";

  for(std::vector<LayerResolver*>::iterator collit=knownDetectors.begin(); collit!=knownDetectors.end(); collit++) {

    if ((*collit)->SetCollection(evt) == 0) nFound++;

  } // End loop over known detectors
  return nFound;

}


double LayerFinder::SensitiveThickness(TrackerHitPlane* thit, int &tf) {

  streamlog_out(DEBUG) << "CollectionFinder::SensitiveThickness() for hit ID = " << thit->getCellID0() << ".\n";

  for(std::vector<LayerResolver*>::iterator cit=knownDetectors.begin(); cit!=knownDetectors.end(); cit++) {

    streamlog_out(DEBUG) << " ... Looking in collection " << (*cit)->GetCollectionName() ;
    streamlog_out(DEBUG) << " of type "
                         << (*cit)->GetCollectionType() << " with " << (*cit)->GetNumberOfHits() << " hits.\n";
    for(int i=0; i<(*cit)->GetNumberOfHits(); i++) {

      if( thit == (*cit)->GetHit(i) ) {

        streamlog_out(DEBUG) << " ... Found matching hit " << (*cit)->GetHit(i)->getCellID0()
                                 << " with energy " << (*cit)->GetHit(i)->getEDep()
                                 << " in collection " << (*cit)->GetCollectionName() << ".\n";

        int nLayer = (*cit)->Decode(thit);
        tf = (*cit)->GetDetTypeFlag();
        int nLayers = (*cit)->GetNumberOfLayers();

        streamlog_out(DEBUG) << " ... Layer number is " << nLayer << ".\n";
        streamlog_out(DEBUG) << " ... Number of layers is " << nLayers << ".\n";
        if(nLayer >= nLayers) {
          streamlog_out(ERROR) << "CollectionFinder::SensitiveThickness() -- Layer number out of range!\n";
          return 1.;
        }


        if (tf & DD4hep::DetType::BARREL) {
          double t = dynamic_cast<PlaneResolver*>(*cit)->SensitiveThickness(nLayer);
          streamlog_out(DEBUG) << " ... Thickness is " << t/dd4hep::mm << " mm.\n";
          return t;
        }
        if (tf & DD4hep::DetType::ENDCAP) {
          double t = dynamic_cast<PetalResolver*>(*cit)->SensitiveThickness(nLayer);
          streamlog_out(DEBUG) << " ... Thickness is " << t/dd4hep::mm << " mm.\n";
          return t;
        }
        streamlog_out(WARNING) << " ... Detector type flag not found!!!\n";

      }
    }
  }
  return -1.;
}


void LayerFinder::ReportKnownDetectors() {

  streamlog_out(DEBUG) << "CollectionFinder knows the following detectors:\n";
  for(std::vector<LayerResolver*>::iterator dit=knownDetectors.begin(); dit!=knownDetectors.end(); dit++) {
    std::string dettype;
    if      ((*dit)->GetDetTypeFlag() & DD4hep::DetType::BARREL) dettype = "BARREL";
    else if ((*dit)->GetDetTypeFlag() & DD4hep::DetType::ENDCAP) dettype = "ENDCAP";
    streamlog_out(DEBUG) << "Detector of type " << dettype;
    streamlog_out(DEBUG) << " associated with collection name " << (*dit)->GetCollectionName() << "\n";
    streamlog_out(DEBUG) << "Currently looking at collection of type " << (*dit)->GetCollectionType();
    streamlog_out(DEBUG) << " with encoding " << (*dit)->GetCollectionEncoding() << "\n";
  }
}

