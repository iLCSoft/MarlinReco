#include "DDStripSplitter.h"
#include <iostream>
#include <map>

using std::endl;

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>

#include <EVENT/SimCalorimeterHit.h>

#include "CalorimeterHitType.h"

#include "IMPL/CalorimeterHitImpl.h"
#include <IMPL/LCRelationImpl.h>

#include <UTIL/LCRelationNavigator.h>

#include <marlin/Global.h>
#include <marlin/Exceptions.h>

#include "DD4hep/DetectorSelector.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"

#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

// needed for uint
#include <sys/types.h>

#include "streamlog/streamlog.h"

#define RELATIONFROMTYPESTR "FromType"
#define RELATIONTOTYPESTR "ToType"


using namespace lcio ;
using namespace marlin ;

DDStripSplitter aDDStripSplitter ;


DDStripSplitter::DDStripSplitter() : Processor("DDStripSplitter") {
  // modify processor description
  _description = "DDStripSplitter applies SSA to a strip-based calorimeter ..." ;

  registerInputCollection( LCIO::CALORIMETERHIT,
			   "ECALcollection_evenLayers",
			   "Name of the ECAL collection in even layers",
			   _ecalCollectionEvenLayers,
			   std::string("ECalBarrelScEvenCollectionRec") );
  
  registerInputCollection( LCIO::CALORIMETERHIT,
			   "ECALcollection_oddLayers",
			   "Name of the ECAL collection in odd layers",
			   _ecalCollectionOddLayers,
			   std::string("ECalBarrelScOddCollectionRec") );

  registerInputCollection( LCIO::LCRELATION,
			   "LCRelations_evenLayers",
			   "name of the relation collection for even layer hits",
			   _inputRelationsColEven,
			   std::string("ECalBarrelScHitsEvenRecRelations") );
  
  registerInputCollection( LCIO::LCRELATION,
			   "LCRelations_oddLayers",
			   "name of the relation collection for odd layer hits",
			   _inputRelationsColOdd,
			   std::string("ECalBarrelScHitsOddRecRelations") );
  
  registerInputCollection( LCIO::MCPARTICLE,
                           "MCParticleCollection",
                           "name of MCParticle collection (used for some plots)",
                           _mcParticleCollectionName,
                           std::string("MCParticle") );

  registerProcessorParameter( "splitEcalCollection",
                              "name of output collection containing split hits",
			      _splitEcalCollection,
                              std::string("ECalBarrelSplitCollection") );

  registerProcessorParameter( "unsplitEcalCollection",
                              "name of output collection containing unsplit hits",
			      _unsplitEcalCollection,
                              std::string("ECalBarrelUnSplitCollection") );

  registerProcessorParameter( "splitEcalRelCol",
                              "name of relation collection for split hits",
			      _splitEcalRelCol,
                              std::string("ECalBarrelSplitRelations") );

  registerProcessorParameter( "stripIntersecCollName",
			      "name of (optional) output collection containing strip intersections",
			      _stripIntersecCollName,
			      std::string("stripIntersections") );

  registerProcessorParameter( "evenStripEndsCollName",
			      "name of (optional) output collection containing ends of strips in even layers",
			      _evenStripEndsCollName,
			      std::string("stripEndsEven") );

  registerProcessorParameter( "oddStripEndsCollName",
                              "name of (optional) output collection containing ends of strips in odd layers",
			      _oddStripEndsCollName,
                              std::string("stripEndsOdd") );

  registerProcessorParameter( "virtualCellsDefault",
			      "number of virtual cells per strip (used if info not found in gear file)",
			      _ecalStrip_default_nVirt,
			      int(1) );

  registerProcessorParameter( "saveIntersectionCollection",
			      "save collection with strip interactions?",
			      _saveIntersections,
			      false);

  registerProcessorParameter( "isBarrelHits",
			      "are hits in these collections in the barrel (true) or endcap (false) ?",
			      _isBarrel,
			      true);
  
  // code for layer info for cellID decoder
  registerProcessorParameter("CellIDLayerString" ,
                             "name of the part of the cellID that holds the layer" , 
                             _cellIDLayerString , 
                             std::string("layer")
                             );
  registerProcessorParameter("CellIDModuleString" ,
                             "name of the part of the cellID that holds the module" , 
                             _cellIDModuleString , 
                             std::string("module")
                             );
  registerProcessorParameter("CellIDStaveString" ,
                             "name of the part of the cellID that holds the stave" , 
                             _cellIDStaveString , 
                             std::string("stave")
                             );


}


void DDStripSplitter::init() { 

  _symmetry=-999;
  _cellSize=999;
  _nVirtual=999;
  _evenIsTransverse=-1;
  _flag_rel.setBit(LCIO::LCREL_WEIGHTED); // for the hit relations

  return;
}

void DDStripSplitter::processRunHeader( LCRunHeader*  /*run*/) { 
   return;
}



void DDStripSplitter::setupGeometry() {


  unsigned int includeFlag(0);
  unsigned int excludeFlag(0);

  if ( _isBarrel ) {
    includeFlag = ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL);
  } else {
    includeFlag = ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP);
  }
  excludeFlag = ( dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD );

  dd4hep::Detector & lcdd = dd4hep::Detector::getInstance();
  const std::vector< dd4hep::DetElement>& theDetectors = dd4hep::DetectorSelector(lcdd).detectors(  includeFlag, excludeFlag ) ;
  streamlog_out( DEBUG2 ) << " getExtension :  includeFlag: " << dd4hep::DetType( includeFlag ) << " excludeFlag: " << dd4hep::DetType( excludeFlag )
			  << "  found : " << theDetectors.size() << "  - first det: " << theDetectors.at(0).name() << std::endl ;

  if( theDetectors.size()  != 1 ){
    std::stringstream es ;
    streamlog_out (ERROR) << " getExtension: selection is not unique (or empty)  includeFlag: " << dd4hep::DetType( includeFlag ) << " excludeFlag: " << dd4hep::DetType( excludeFlag )
			  << " --- found detectors : " ;
    for( unsigned i=0, N= theDetectors.size(); i<N ; ++i ){
      streamlog_out (ERROR) << theDetectors.at(i).name() << ", " ;
    }
    throw marlin::StopProcessingException(this);
  }

  int ical=0;
  _caloGeomData = theDetectors.at(ical).extension<dd4hep::rec::LayeredCalorimeterData>();
  if ( ! _caloGeomData ) {
    streamlog_out ( WARNING ) << "could not get calorimeter geometry information!" << endl;
    throw marlin::StopProcessingException(this);
  }
  _symmetry = _caloGeomData->inner_symmetry;
  _stripLength = _caloGeomData->layers[0].cellSize0;
  _stripWidth  = _caloGeomData->layers[0].cellSize1;

  // check that layer 1 is orthogonal to this one
  float lay1_strSize0 = _caloGeomData->layers[1].cellSize0;
  float lay1_strSize1 = _caloGeomData->layers[1].cellSize1;
  if ( fabs(lay1_strSize0-_stripWidth)/_stripWidth > 0.05 || fabs(lay1_strSize1-_stripLength)/_stripLength > 0.05 ) {
    streamlog_out( ERROR ) << "doesn't look like a basic strip-based calo?? don't know how to deal with this geometry!" << endl;
    for (size_t ilay = 0 ; ilay<_caloGeomData->layers.size(); ilay++ ) {
      streamlog_out( MESSAGE ) << "strip size in layer " << ilay << " : " << _caloGeomData->layers[ilay].cellSize0 << " " << _caloGeomData->layers[ilay].cellSize1 << endl;
    }
    throw marlin::StopProcessingException(this);
  }

  _evenIsTransverse = 0;
  if ( _stripLength<_stripWidth ) {
    _evenIsTransverse = 1;
    float temp = _stripLength;
    _stripLength=_stripWidth;
    _stripWidth=temp;
  }
  _nVirtual = int(_stripLength/_stripWidth);
  _stripAspectRatio = _stripLength/_stripWidth;

  // convert from cm to mm (dd4hep units -> lcio units...)
  _stripLength*=10;
  _stripWidth*=10;

  streamlog_out( DEBUG ) << "strip length, width = " << _stripLength << " " << _stripWidth << " mm " << endl;
  if ( _stripAspectRatio < 2. ) {
    streamlog_out( ERROR ) << "this strip is very short: probably not worth using the strip splitter!" << endl;
    throw marlin::StopProcessingException(this);
  }

  return;
}


void DDStripSplitter::processEvent( LCEvent * evt ) { 

  if (_symmetry<0) setupGeometry();

  if ( _stripAspectRatio<2.0 ) {
    streamlog_out ( WARNING ) << " -- not a long enough strip, not worth splitting -- length, width, aspect ratio: " << 
      _stripLength << " " << _stripWidth << " " << _stripAspectRatio << ", doing nothing" << endl;
    return;
  }

  std::pair < TVector3, TVector3 > stripEnds;
  std::string toSplit;
  std::string toSplitRel;
  int orientation;

  //  std::map < IMPL::LCCollectionVec*, std::string > outputcolls;

  if (_saveIntersections) {
    intersectionHits = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
    intersectionHits->setFlag(intersectionHits->getFlag()|( 1 << LCIO::RCHBIT_LONG)); // store position

    stripEndsEvenCol = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
    stripEndsEvenCol->setFlag(stripEndsEvenCol->getFlag()|( 1 << LCIO::RCHBIT_LONG)); // store position

    stripEndsOddCol = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
    stripEndsOddCol->setFlag(stripEndsOddCol->getFlag()|( 1 << LCIO::RCHBIT_LONG)); // store position
  }

  IMPL::LCCollectionVec* splitStripHits   = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
  splitStripHits->setFlag(splitStripHits->getFlag()|( 1 << LCIO::RCHBIT_LONG)); // store position
  IMPL::LCCollectionVec* unSplitStripHits = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
  unSplitStripHits->setSubset();

  LCCollectionVec *splitRelCol  = new LCCollectionVec(LCIO::LCRELATION);
  splitRelCol->setFlag(_flag_rel.getFlag());
  splitRelCol->parameters().setValue( RELATIONFROMTYPESTR , LCIO::CALORIMETERHIT ) ;
  splitRelCol->parameters().setValue( RELATIONTOTYPESTR   , LCIO::SIMCALORIMETERHIT ) ;

  for (int icol=0; icol<2; icol++) { // loop over strip collections in even and odd layers (assumed to have perpendicular orientations)
    switch (icol) {
    case 0:
      if ( _evenIsTransverse==1 ) {
	orientation = TRANSVERSE;
      } else if ( _evenIsTransverse==0 ) {
	orientation = LONGITUDINAL;
      } else {
	streamlog_out( ERROR) << "strip orientation undefined!!" << endl;
	throw marlin::StopProcessingException(this);
      }
      toSplit = _ecalCollectionEvenLayers;
      toSplitRel = _inputRelationsColEven;
      break;
    case 1:
      if ( _evenIsTransverse==1 ) {
	orientation = LONGITUDINAL;
      } else if ( _evenIsTransverse==0 ) {
	orientation = TRANSVERSE;
      } else {
	streamlog_out( ERROR ) << "strip orientation undefined!!" << endl;
	throw marlin::StopProcessingException(this);
      }
      toSplit = _ecalCollectionOddLayers;
      toSplitRel = _inputRelationsColOdd;
      break;
    default:
      streamlog_out ( ERROR ) << "crazy stuff!!! abandoning event..." << endl;
      throw marlin::StopProcessingException(this);
      return;
    }

    try {
      LCCollection * col = evt->getCollection( toSplit.c_str() );
      if (!col) continue;

      LCCollection * colrel = evt->getCollection( toSplitRel.c_str() );
      if (!colrel) continue;
      LCRelationNavigator navi(colrel);

      const std::string layerCodingString(col->getParameters().getStringVal(LCIO::CellIDEncoding));
      // set up cellid encoding for new collections
      splitStripHits->parameters().setValue(LCIO::CellIDEncoding, layerCodingString);
      unSplitStripHits->parameters().setValue(LCIO::CellIDEncoding, layerCodingString);

      // get the cellid decoder for this collection
      CellIDDecoder<CalorimeterHit> id( col ) ;
      _decoder = &id;
	
      // loop over the collection's hits
      int nelem = col->getNumberOfElements();
      for (int j=0; j < nelem; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(j) );
	if (!hit) {
	  streamlog_out ( ERROR ) << "ERROR  null hit in collection " << toSplit.c_str() << " " << j << endl;
	  throw marlin::StopProcessingException(this);
	  continue;
	}
	// split the hits
	std::vector <CalorimeterHit*> splitHits = getVirtualHits(evt, hit, orientation, _isBarrel); // assume the first one (should be only one)

	// find the corresponding simhit (for the relations)
	SimCalorimeterHit* simhit =  (SimCalorimeterHit*) navi.getRelatedToObjects(hit)[0];

	// add (new) hits to collections
	if (splitHits.size()==0) { // not split, add original hit
	  unSplitStripHits->addElement(hit);
          splitRelCol->addElement( new LCRelationImpl(hit,simhit,1.0) );
	} else { // split was split, add the virtual hits
	  for (uint hh=0; hh<splitHits.size(); hh++) {
	    splitStripHits->addElement(splitHits[hh]);
	    // weight = fraction of orig hit which is in this split hit <---- to be discussed!!!!! DJ
	    float weight = splitHits[hh]->getEnergy() / hit->getEnergy(); 
	    splitRelCol->addElement( new LCRelationImpl(splitHits[hh],simhit,weight) );
	  }
	}
      } // loop over hits in collection
    } catch(DataNotAvailableException &e) {};
  } // long/trans loop

  // add the new collections to the event
  evt->addCollection( splitStripHits, _splitEcalCollection );
  evt->addCollection( unSplitStripHits, _unsplitEcalCollection );
  evt->addCollection( splitRelCol, _splitEcalRelCol );

  if (_saveIntersections) {
    evt->addCollection(intersectionHits,_stripIntersecCollName);
    evt->addCollection(stripEndsEvenCol,_evenStripEndsCollName);
    evt->addCollection(stripEndsOddCol, _oddStripEndsCollName);
  }

  return;
}

std::vector <CalorimeterHit*> DDStripSplitter::getVirtualHits(LCEvent* evt, CalorimeterHit* hit, int orientation, bool barrel ) {

  // this splits the strip into zero or more hits along its length
  // by looking at nearby hits with different orientation (trans/long or square)

  int layer  = (*_decoder)(hit)[_cellIDLayerString];
  int module = (*_decoder)(hit)[_cellIDModuleString];
  int stave  = (*_decoder)(hit)[_cellIDStaveString];

  TVector3 pp;
  pp.SetXYZ( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);

  std::vector <CalorimeterHit*> newhits;

  // get the ends of this strip
  std::pair < TVector3, TVector3 > stripEnds = getStripEnds(hit, orientation, barrel);
  TVector3 stripDir = stripEnds.first - stripEnds.second;

  if (_saveIntersections) { 
    // make new collection with a hit at each end of a strip
    //   for debugging purposes
    float ppp[3];
    ppp[0] = stripEnds.first.X();
    ppp[1] = stripEnds.first.Y();
    ppp[2] = stripEnds.first.Z();
    CalorimeterHitImpl* interhit = new CalorimeterHitImpl();
    interhit->setPosition( ppp );
    interhit->setEnergy(0.03);

    if ( layer%2 == 0 ) stripEndsEvenCol->addElement(interhit);
    else                stripEndsOddCol->addElement(interhit);

    ppp[0] = stripEnds.second.X();
    ppp[1] = stripEnds.second.Y();
    ppp[2] = stripEnds.second.Z();
    interhit = new CalorimeterHitImpl();
    interhit->setPosition( ppp );
    interhit->setEnergy(0.03);

    if ( layer%2 == 0 ) stripEndsEvenCol->addElement(interhit);
    else                stripEndsOddCol->addElement(interhit);
  }

  // decide which collections to use to split the strip
  int splitterOrientation;
  std::string splitterCol;
  if ( layer%2 == 0 ) {
    splitterCol = _ecalCollectionOddLayers;
    splitterOrientation = _evenIsTransverse ? LONGITUDINAL : TRANSVERSE ;
  } else {
    splitterCol = _ecalCollectionEvenLayers;
    splitterOrientation = _evenIsTransverse ? TRANSVERSE : LONGITUDINAL ;
  }

  std::map <int, float> virtEnergy;
  int nSplitters(0);

  // loop over splitter cols, find nearby hits
  for (int jj=0; jj<2; jj++) { // strips, cells

    try {
      LCCollection * col = evt->getCollection( splitterCol.c_str() );
      if (!col) continue;

      CellIDDecoder<CalorimeterHit> id( col ) ;
      _decoder2 = &id;

      int nelem = col->getNumberOfElements();

      for (int j=0; j < nelem; ++j) {
	CalorimeterHit * hit2 = dynamic_cast<CalorimeterHit*>(col->getElementAt(j) );
	if (!hit2) {
	  streamlog_out ( ERROR ) << "ERROR  null hit2 in collection " <<  splitterCol << " " << j << endl;
	  throw marlin::StopProcessingException(this);
	  continue;
	}
	int layer2  = (*_decoder)(hit2)[_cellIDLayerString];
	int module2 = (*_decoder)(hit2)[_cellIDModuleString];
	int stave2  = (*_decoder)(hit2)[_cellIDStaveString];
	
	int dlayer = abs(layer2-layer);
	int dstave = abs(stave2-stave);
	int dmodule = abs(module2-module);

	// are the two hits close enough to look at further?

	// if hits in same module and same stave, require that only one layer difference
	if (dmodule==0 && dstave==0 && dlayer>1) continue;

	if (barrel) {
	  dstave = min( dstave, _symmetry-dstave);
	  if ( dstave==0 && dmodule>1 ) continue; // allow same stave and +- 1 module
	  if ( dmodule==0 && dstave>1 ) continue; // or same module +- 1 stave
	  if ( dstave==0 && dlayer>1) continue;   // if in same stave, require dlayer==1
	} else { // endcap
	  dstave = min( dstave, 4-dstave);
	  if (dmodule!=0) continue; // different endcap
	  if (dstave>1) continue;   // more than 1 stave (=quarter endcap) apart
	  if (dlayer>1) continue;   // more than 1 layer apart
	}

	// simple distance check for remaining hit pairs
	float dist = sqrt( pow(hit2->getPosition()[0] - hit->getPosition()[0], 2) + 
			   pow(hit2->getPosition()[1] - hit->getPosition()[1], 2) + 
			   pow(hit2->getPosition()[2] - hit->getPosition()[2], 2) );
	
	if (dist>2*_stripLength) continue;

	// for remaining hits, check if they overlap
	TVector3 stripDir2(0,0,0);
	if (jj==0) { //strip
	  std::pair < TVector3, TVector3 > stripEnds2 = getStripEnds(hit2, splitterOrientation, barrel);
	  stripDir2 = stripEnds2.first - stripEnds2.second;
	} // leave 0 for cell
	
	  // check if strips intersect
	TVector3 intercept = stripIntersect(hit, stripDir, hit2, stripDir2);
	if (intercept.Mag()>0) { // intercept found, calculate in which virtual cell
	  nSplitters++;
	  float frac(-1);
	  for (int ii=0; ii<3; ii++) {
	    float dx = stripEnds.second[ii] - stripEnds.first[ii];
	    if (fabs(dx)>0.1) {
	      frac = (intercept[ii]-stripEnds.first[ii])/dx;
	      break;
	    }
	  }

	  if (frac>=0.0 && frac<=1.0) {
	    int segment = int(frac*_nVirtual);
	    if (segment>=0 && segment<_nVirtual) {
	      if (virtEnergy.find(segment)!=virtEnergy.end()) {
		virtEnergy[segment] += hit2->getEnergy();
	      } else {
		virtEnergy[segment] = hit2->getEnergy();
	      }
	      
	      if (_saveIntersections) {
		CalorimeterHitImpl* interhit = new CalorimeterHitImpl();
		float pos[3];
		pos[0] = intercept.X();
		pos[1] = intercept.Y();
		pos[2] = intercept.Z();
		interhit->setPosition( pos );
		interhit->setEnergy(0.1);
		intersectionHits->addElement(interhit);
	      }
	      
	    } else {
	      streamlog_out ( WARNING ) << "strange segment " << segment << " frac = " << frac << " nvirt = " << _nVirtual << endl;
	    }
	  } else {
	    streamlog_out ( WARNING ) << "strange frac " << frac << endl;
	  }
	  
	}
      }
    } catch(DataNotAvailableException &e) {};
  }


  // now create the virtual cells, and assign energy
  std::map <int, float>::iterator it;  
  float totenergy(0);
  for (it=virtEnergy.begin(); it!=virtEnergy.end(); it++) {
    totenergy+=it->second;
  }
  for (it=virtEnergy.begin(); it!=virtEnergy.end(); it++) {
    // energy of hit
    float energy = hit->getEnergy()*it->second/totenergy;
    // position of hit
    TVector3 virtualCentre = stripEnds.second-stripEnds.first;
    virtualCentre*= (it->first + 0.5)/_nVirtual;
    virtualCentre += stripEnds.first;
    // make the new hit
    CalorimeterHitImpl* newhit = new CalorimeterHitImpl();

    float pos[3];
    pos[0] = virtualCentre.X();
    pos[1] = virtualCentre.Y();
    pos[2] = virtualCentre.Z();

    newhit->setType(hit->getType());
    newhit->setPosition( pos );
    newhit->setEnergy(energy);

    newhit->setCellID0(hit->getCellID0());
    newhit->setCellID1(hit->getCellID1());

    newhits.push_back(newhit);

  }

  return newhits;
}


std::pair < TVector3, TVector3 > DDStripSplitter::getStripEnds(CalorimeterHit* hit, int orientation, bool barrel) {
  // calculate the positions of the strip ends

  TVector3 stripcentre(hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);
  TVector3 stripend1(stripcentre);
  TVector3 stripend2(stripcentre);

  int stave  = (*_decoder)(hit)[_cellIDStaveString];

  if (barrel) {
    if (orientation == TRANSVERSE) { // transverse, along z axis in barrel region
      stripend1.SetZ(stripcentre.Z() - _stripLength/2.);
      stripend2.SetZ(stripcentre.Z() + _stripLength/2.);
    } else if (orientation == LONGITUDINAL) { // longitudinal, along x-y in barrel
      float phiRotAngle = (stave-1)*2.*TMath::Pi()/_symmetry; // for new ILD models (dd4hep based)
      TVector3 stripHalfVec(_stripLength/2.,0,0);
      stripHalfVec.RotateZ(phiRotAngle - TMath::Pi()/2.);
      stripend1-=stripHalfVec;
      stripend2+=stripHalfVec;
    }
  } else { // endcap
    // is the slab horizontal?
    bool horizontalSlab = stave%2==1;
    // are the strips in the slab horizontal?
    bool horizontalStrip(true);
    if (horizontalSlab) {
      if      (orientation == LONGITUDINAL) horizontalStrip=false;
      else if (orientation == TRANSVERSE)   horizontalStrip=true;
    } else {
      if      (orientation == LONGITUDINAL) horizontalStrip=true;
      else if (orientation == TRANSVERSE)   horizontalStrip=false;
    }

    if (horizontalStrip) {
      stripend1.SetX(stripcentre.X() - _stripLength/2.);
      stripend2.SetX(stripcentre.X() + _stripLength/2.);
    } else {
      stripend1.SetY(stripcentre.Y() - _stripLength/2.);
      stripend2.SetY(stripcentre.Y() + _stripLength/2.);
    }

  }

  std::pair < TVector3, TVector3 > output (stripend1, stripend2);

  return output;
}


void DDStripSplitter::check( LCEvent *  /*evt*/ ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void DDStripSplitter::end(){
  streamlog_out ( MESSAGE ) << "DDStripSplitter::end()  " << name() << std::endl ;
}

TVector3 DDStripSplitter::stripIntersect(CalorimeterHit* hit0, TVector3 axis0, CalorimeterHit* hit1, TVector3 axis1) {
  // find intercept of hit1 with hit0
  // hit0 must be a strip
  // hit1 can be an orthogonal strip, or a cell
  // axis0,1 are direction of strip

  // centre position of cell/strip
  TVector3 stripCentre[2];
  stripCentre[0].SetXYZ( hit0->getPosition()[0], hit0->getPosition()[1], hit0->getPosition()[2] );
  stripCentre[1].SetXYZ( hit1->getPosition()[0], hit1->getPosition()[1], hit1->getPosition()[2] );

  // direction of strip long axis
  // 0,0,0 for square cell
  TVector3 stripDir[2];
  stripDir[0]=axis0;
  stripDir[1]=axis1;

  // deal with cell case
  // define it's direction as perpendicular to strip and vector cell centre to origin
  bool isStrip[2];
  for (int i=0; i<2; i++) {
    if (stripDir[i].Mag()>1e-10) {
      isStrip[i]=true;
    } else {
      isStrip[i]=false;
      stripDir[i] = stripDir[1-i].Cross(stripCentre[i]);
    }
    // ensure dir is normalised
    stripDir[i]*=1./stripDir[i].Mag();
  }

  if (!isStrip[0]) {
    streamlog_out ( ERROR ) << "ERROR from DDStripSplitter::stripIntersect, first hit should be a strip" << endl;
    throw marlin::StopProcessingException(this);
  }

  TVector3 p[2][2]; // ends of strips
  for (int j=0; j<2; j++) {
    for (int i=0; i<2; i++) {
      float ll = isStrip[j] ? _stripLength : _stripWidth*1.1; // inflate a little for cells
      p[j][i] = stripCentre[j]-TMath::Power(-1, i)*0.5*ll*stripDir[j];
    }
  }

  TVector3 pNorm[2][2];
  for (int j=0; j<2; j++) {
    for (int i=0; i<2; i++) {
      float mm = p[j][i].Mag();
      pNorm[j][i]=p[j][i];
      pNorm[j][i]*=1./mm;
    }
  }


  TVector3 inPlane[2]; // difference between points: line inside plane
  for (int j=0; j<2; j++) {
    inPlane[j]=p[j][0]-p[j][1];
  }

  // vector normal to both lines (this is normal to the plane we're interested in)
  TVector3 normal = inPlane[0].Cross(inPlane[1]);
  float mag = normal.Mag();
  normal*=1./mag;

  // point on line [0]
  TVector3 point(p[0][0]+p[0][1]); point*=0.5;

  // calculate the projected positions of ends of [1] on 
  //    the plane which contains "point" with normal "normal"
  TVector3 qPrime[2];
  for (int i=0; i<2; i++) {
    float d = ((point - p[1][i]).Dot(normal))/(pNorm[1][i].Dot(normal));
    qPrime[i] = p[1][i] + d*pNorm[1][i];
  }

  // find the intersection of lines qPrime and p (they are in same plane)
  TVector3 a(p[0][1]-p[0][0]);
  TVector3 b(qPrime[1]-qPrime[0]);
  TVector3 c(qPrime[0]-p[0][0]);
  
  float factor = (c.Cross(b)).Dot(a.Cross(b))/((a.Cross(b)).Mag2());

  TVector3 x = p[0][0] + a*factor;

  // check that two lines really intercept
  bool intersect = true;
  for (int ii=0; ii<2; ii++) {
    for (int j=0; j<3; j++) {
      float d0 = ii==0 ? x[j]-p[0][0][j] : x[j]-qPrime[0][j];
      float d1 = ii==0 ? x[j]-p[0][1][j] : x[j]-qPrime[1][j];
      if (d0*d1>1e-10) {
	intersect = false;
      }
    }
  }

  if (intersect) return x;
  else return TVector3(0,0,0);
}
