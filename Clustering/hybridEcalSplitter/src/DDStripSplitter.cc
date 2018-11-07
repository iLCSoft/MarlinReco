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

#include <UTIL/LCRelationNavigator.h>

#include <TFile.h>
#include <TH2F.h>
#include <TMath.h>
#include <TVector3.h>
#include "TLorentzVector.h"

#include <marlin/Global.h>

#include "DD4hep/DetectorSelector.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"



#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

// needed for uint
#include <sys/types.h>

using namespace lcio ;
using namespace marlin ;

DDStripSplitter aDDStripSplitter ;


DDStripSplitter::DDStripSplitter() : Processor("DDStripSplitter") {
  // modify processor description
  _description = "DDStripSplitter applies SSA to a strip-based calorimeter ..." ;

  std::vector <std::string> ecalCollectionsTranStrips;
  ecalCollectionsTranStrips.push_back(std::string("ECALScTransverseBarrel"));
  registerInputCollections( LCIO::CALORIMETERHIT,
                            "ECALcollections_tranStrips",
			    "Name of the ECAL transverse strip collections",
			    _ecalCollectionsTranStrips,
			    ecalCollectionsTranStrips);
  
  std::vector <std::string> ecalCollectionsLongStrips;
  ecalCollectionsLongStrips.push_back(std::string("ECALScLongitudinalBarrel"));
  registerInputCollections( LCIO::CALORIMETERHIT,
                            "ECALcollections_longStrips",
			    "Name of the ECAL longitudinal strip collections",
			    _ecalCollectionsLongStrips,
			    ecalCollectionsLongStrips);

  registerProcessorParameter( "virtualCellsDefault",
			      "number of virtual cells per strip (used if info not found in gear file)",
			      _ecalStrip_default_nVirt,
			      int(1) );

  registerProcessorParameter( "saveIntersectionCollection",
			      "save collection with strip interactions?",
			      _saveIntersections,
			      false);
  
  registerProcessorParameter( "saveCheckRootHistograms",
			      "save root file with debugging histograms?",
			      _makePlots,
			      false);

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

  _fout=NULL;

  _symmetry=-999;
  _cellSize=999;
  _nVirtual=999;

  if (_makePlots) { // define output root file and some histograms to put in it
    _fout = new TFile("hybridCheckHistos.root","recreate");
    h_phiModuleCheck = new TH2F("phi_HitModule","phi_HitModule", 20, -4, 16, 100, -2*TMath::Pi(), 2*TMath::Pi());

    h_phiThetaMC = new TH2F("MCphiTheta", "MCphiTheta", 100, -2*TMath::Pi(), 2*TMath::Pi(), 100, -2*TMath::Pi(), 2*TMath::Pi());

    h_stripDist_intercept   = new TH1F("dist_in","dist_in",100,0,2);
    h_stripDist_nointercept = new TH1F("dist_no","dist_no",100,0,2);

    for (int i=0; i<2; i++) {
      TString reg = i==0 ? "_barrel" : "_endcap" ;
      h_stavemodule[i] = new TH2F("stavemodule"+reg,"stavemodule"+reg,35, 0, 35, 35,0,35);;
      h_layer[i]  = new TH1F("layer"+reg, "layer"+reg, 35,0,35);;

      h_staveX[i]  = new TH2F("staveX_"+reg,  "staveX_"+reg,  100,-2000,2000,10,0,10);
      h_staveY[i]  = new TH2F("staveY_"+reg,  "staveY_"+reg,  100,-2000,2000,10,0,10);
      h_staveZ[i]  = new TH2F("staveZ_"+reg,  "staveZ_"+reg,  100,-5000,5000,10,0,10);
      h_moduleX[i] = new TH2F("moduleX_"+reg, "moduleX_"+reg, 100,-2000,2000,10,0,10);
      h_moduleY[i] = new TH2F("moduleY_"+reg, "moduleY_"+reg, 100,-2000,2000,10,0,10);
      h_moduleZ[i] = new TH2F("moduleZ_"+reg, "moduleZ_"+reg, 100,-5000,5000,10,0,10);
    }
    
    for (int ieb=0; ieb<2; ieb++) {
      for (int is=0; is<10; is++) {
	for (int im=0; im<10; im++) {
	  TString hhname = ieb==0 ? "barrel" : "endcap";
	  hhname+="_stave"; hhname+=is;
	  hhname+="_module"; hhname+=im;
	  h_cth_phi[ieb][is][im] = new TH2F(hhname+"_angle",hhname+"_angle",100,-1,1,100,-TMath::Pi(), TMath::Pi());
	  h_XY[ieb][is][im] = new TH2F(hhname+"_xy",hhname+"_xy",100,-3000,3000,100,-3000,3000);
	}
      }
    }
  }

  return;
}

void DDStripSplitter::processRunHeader( LCRunHeader*  /*run*/) { 
   return;
}



void DDStripSplitter::setupGeometry() {


  unsigned int includeFlag(0);
  unsigned int excludeFlag(0);

  includeFlag = ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL);
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
    assert(0);
  }

  int ical=0;
  _caloGeomData = theDetectors.at(ical).extension<dd4hep::rec::LayeredCalorimeterData>();
  if ( ! _caloGeomData ) {
    streamlog_out ( WARNING ) << "could not get calorimeter geometry information!" << endl;
    assert(0);
  }
  _symmetry = _caloGeomData->inner_symmetry;
  _stripLength = _caloGeomData->layers[0].cellSize0;
  _stripWidth  = _caloGeomData->layers[0].cellSize1;

  if ( _stripLength<_stripWidth ) {
    float temp = _stripLength;
    _stripLength=_stripWidth;
    _stripWidth=temp;
  }
  _nVirtual = int(_stripLength/_stripWidth);
  _stripAspectRatio = _stripLength/_stripWidth;

  return;
}


void DDStripSplitter::processEvent( LCEvent * evt ) { 

  if (_symmetry<0) setupGeometry();

  if (_makePlots) {
    // first fill some simple MC histos
    try {
      LCCollection * col = evt->getCollection("MCParticle");
      if (col->getNumberOfElements()>0) {
	MCParticle * mcp = dynamic_cast<MCParticle*>(col->getElementAt(0) );
	TVector3 mom(mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2]);
	h_phiThetaMC->Fill(mom.Phi(), mom.Theta());
      }
    }
    catch(DataNotAvailableException &e) {};
  }

  if ( _stripAspectRatio<2.0 ) {
    streamlog_out ( WARNING ) << " -- not a long enough strip, not worth splitting -- length, width, aspect ratio: " << 
      _stripLength << " " << _stripWidth << " " << _stripAspectRatio << ", doing nothing" << endl;
    return;
  }

  std::pair < TVector3, TVector3 > stripEnds;
  std::vector <std::string> * toSplit;
  //  std::vector <std::string> * stripSplitter;
  int orientation;

  std::map < IMPL::LCCollectionVec*, std::string > outputcolls;

  if (_saveIntersections) {
    intersectionHits = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
    intersectionHits->setFlag(intersectionHits->getFlag()|( 1 << LCIO::RCHBIT_LONG)); // store position

    stripEndsTransCol = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
    stripEndsTransCol->setFlag(stripEndsTransCol->getFlag()|( 1 << LCIO::RCHBIT_LONG)); // store position

    stripEndsLongCol = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
    stripEndsLongCol->setFlag(stripEndsLongCol->getFlag()|( 1 << LCIO::RCHBIT_LONG)); // store position
  }

  for (int icol=0; icol<2; icol++) { // loop over longitudinal and transverse strip collections
    switch (icol) {
    case 0:
      orientation = TRANSVERSE;
      toSplit = &_ecalCollectionsTranStrips;
      //stripSplitter = &_ecalCollectionsLongStrips;
      break;
    case 1:
      orientation = LONGITUDINAL;
      toSplit = &_ecalCollectionsLongStrips;
      //stripSplitter = &_ecalCollectionsTranStrips;
      break;
    default:
      streamlog_out ( ERROR ) << "ERROR crazy stuff!!! abandoning event..." << endl;
      return;
    }

    for (uint i=0; i<toSplit->size(); i++) { // loop over collections of this type (long/trans)
      try {
	LCCollection * col = evt->getCollection( toSplit->at(i).c_str() );
	if (!col) continue;

	const std::string layerCodingString(col->getParameters().getStringVal(LCIO::CellIDEncoding));

	// is this a barrel or endcap collection?
	TString hhname = toSplit->at(i);
	bool barrel(true);
	if (hhname.Contains("Barrel") || hhname.Contains("barrel")) {
	  barrel = true;
	} else if (hhname.Contains("Endcap") || hhname.Contains("endcap")) {
	  barrel = false;
	} else {
	  streamlog_out ( ERROR ) << "WARNING: cannot tell if collection is for barrel or endcap..." << hhname << endl;
	}

	// make new collections for split and unsplit strips
	IMPL::LCCollectionVec* splitStripHits = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
	splitStripHits->setFlag(splitStripHits->getFlag()|( 1 << LCIO::RCHBIT_LONG)); // store position
	splitStripHits->parameters().setValue(LCIO::CellIDEncoding, layerCodingString);

	IMPL::LCCollectionVec* unSplitStripHits = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
	unSplitStripHits->setSubset();
	unSplitStripHits->parameters().setValue(LCIO::CellIDEncoding, layerCodingString);
	
	// get the cellid decoder for this collection
	CellIDDecoder<CalorimeterHit> id( col ) ;
	_decoder = &id;
	
	// loop over the collection's hits
	int nelem = col->getNumberOfElements();
	for (int j=0; j < nelem; ++j) {
	  CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(j) );
	  if (!hit) {
	    streamlog_out ( ERROR ) << "ERROR  null hit in collection " << toSplit->at(i).c_str() << " " << j << endl;
	    continue;
	  }
	  // split the hits
	  std::vector <CalorimeterHit*> splitHits = getVirtualHits(evt, hit, orientation, barrel);

	  // add (new) hits to collections
	  if (splitHits.size()==0) { // not split, add original hit
	    unSplitStripHits->addElement(hit);
	  } else { // split was split, add the virtual hits
	    for (uint hh=0; hh<splitHits.size(); hh++) {
	      splitStripHits->addElement(splitHits[hh]);
	    }
	  }

	} // loop over hits in collection

	if (splitStripHits->getNumberOfElements()>0)
	  outputcolls[splitStripHits] = toSplit->at(i) + "Split";
	if (unSplitStripHits->getNumberOfElements()>0)
	  outputcolls[unSplitStripHits] = toSplit->at(i) + "UnSplit";

      } catch(DataNotAvailableException &e) {};
    } // loop over collections
  } // long/trans loop

  std::map < IMPL::LCCollectionVec*, std::string >::iterator ii;

  // add the new collections to the event
  for (ii=outputcolls.begin(); ii!=outputcolls.end(); ii++) {
    evt->addCollection(ii->first, ii->second);
  }

  if (_saveIntersections) {
    evt->addCollection(intersectionHits,"stripIntersections");
    evt->addCollection(stripEndsTransCol,"stripEndsT");
    evt->addCollection(stripEndsLongCol, "stripEndsL");
  }

  return;
}

std::vector <CalorimeterHit*> DDStripSplitter::getVirtualHits(LCEvent* evt, CalorimeterHit* hit, int orientation, bool barrel ) {

  // this splits the strip into zero or more hits along its length
  // by looking at nearby hits with different orientation (trans/long or square)

  int ieb=!barrel;

  int layer  = (*_decoder)(hit)[_cellIDLayerString];
  int module = (*_decoder)(hit)[_cellIDModuleString];
  int stave  = (*_decoder)(hit)[_cellIDStaveString];

  TVector3 pp;
  pp.SetXYZ( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);

  if (_makePlots) {
    h_stavemodule[ieb]->Fill(stave, module);
    h_layer[ieb] ->Fill(layer);

    h_staveX[ieb] ->Fill(hit->getPosition()[0], stave);
    h_staveY[ieb] ->Fill(hit->getPosition()[1], stave);
    h_staveZ[ieb] ->Fill(hit->getPosition()[2], stave);
    h_moduleX[ieb]->Fill(hit->getPosition()[0], module);
    h_moduleY[ieb]->Fill(hit->getPosition()[1], module);
    h_moduleZ[ieb]->Fill(hit->getPosition()[2], module);

    h_cth_phi[ieb][stave][module]->Fill(pp.CosTheta(), pp.Phi());
    h_XY[ieb][stave][module]->Fill(pp.X(), pp.Y());
  }

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

    if (orientation==TRANSVERSE) stripEndsTransCol->addElement(interhit);
    else if (orientation==LONGITUDINAL) stripEndsLongCol->addElement(interhit);

    ppp[0] = stripEnds.second.X();
    ppp[1] = stripEnds.second.Y();
    ppp[2] = stripEnds.second.Z();
    interhit = new CalorimeterHitImpl();
    interhit->setPosition( ppp );
    interhit->setEnergy(0.03);

    if (orientation==TRANSVERSE) stripEndsTransCol->addElement(interhit);
    else if (orientation==LONGITUDINAL) stripEndsLongCol->addElement(interhit);
  }

  // decide which collections to use to split the strip
  int splitterOrientation;
  std::vector <std::string> * splitterCols;
  if ( orientation==TRANSVERSE ) {
    splitterOrientation = LONGITUDINAL;
    splitterCols = &_ecalCollectionsLongStrips;
  } else if ( orientation==LONGITUDINAL ) {
    splitterOrientation = TRANSVERSE;
    splitterCols = &_ecalCollectionsTranStrips;
  } else {
    streamlog_out ( DEBUG ) << "no need to split this orientation";
    return newhits;
  }

  std::map <int, float> virtEnergy;
  int nSplitters(0);

  // loop over splitter cols, find nearby hits
  for (int jj=0; jj<2; jj++) { // strips, cells

    std::vector <std::string> * splitter = splitterCols;

    for (uint i=0; i<splitter->size(); i++) {
      try {
	LCCollection * col = evt->getCollection( splitter->at(i).c_str() );
	if (!col) continue;

	CellIDDecoder<CalorimeterHit> id( col ) ;
	_decoder2 = &id;

	int nelem = col->getNumberOfElements();

	for (int j=0; j < nelem; ++j) {
	  CalorimeterHit * hit2 = dynamic_cast<CalorimeterHit*>(col->getElementAt(j) );
	  if (!hit2) {
	    streamlog_out ( ERROR ) << "ERROR  null hit2 in collection " <<  splitter->at(i).c_str() << " " << j << endl;
	    continue;
	  }
	  
	  // int layer2  = (*_decoder)(hit2)["K-1"];
	  // int module2 = (*_decoder)(hit2)["M"];
	  // int stave2  = (*_decoder)(hit2)["S-1"];

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
	  if (_makePlots) {
	    if (intercept.Mag()>0) h_stripDist_intercept  ->Fill(dist/_stripLength);
	    else                   h_stripDist_nointercept->Fill(dist/_stripLength);
	  }
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
      if (_makePlots) {
	h_phiModuleCheck->Fill(1.0*stave+0.5, stripcentre.Phi());
      }
      //      float phiRotAngle = stave*2.*TMath::Pi()/_symmetry; // some definition has changed
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
  if (_makePlots && _fout) {
    _fout->Write(0);
    _fout->Close();
  }

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

  if (!isStrip[0]) streamlog_out ( ERROR ) << "ERROR from DDStripSplitter::stripIntersect, first hit should be a strip" << endl;

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
