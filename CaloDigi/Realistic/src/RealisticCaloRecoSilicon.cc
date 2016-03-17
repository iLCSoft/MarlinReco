#include "RealisticCaloRecoSilicon.h"
#include <algorithm>
#include <cassert>
#include <iostream>

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"
#include "DDRec/DetectorData.h"

using std::cout;
using std::endl;

RealisticCaloRecoSilicon aRealisticCaloRecoSilicon;

RealisticCaloRecoSilicon::RealisticCaloRecoSilicon() : RealisticCaloReco::Processor("RealisticCaloRecoSilicon") {
  _description = "Performs fist reconstruction of silicon ECAL hits";
}

DD4hep::DDRec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag=0) {


  DD4hep::DDRec::LayeredCalorimeterData * theExtension = 0;

  DD4hep::Geometry::LCDD & lcdd = DD4hep::Geometry::LCDD::getInstance();
  const std::vector< DD4hep::Geometry::DetElement>& theDetectors = DD4hep::Geometry::DetectorSelector(lcdd).detectors(  includeFlag, excludeFlag ) ;


  streamlog_out( DEBUG2 ) << " getExtension :  includeFlag: " << DD4hep::DetType( includeFlag ) << " excludeFlag: " << DD4hep::DetType( excludeFlag )
              << "  found : " << theDetectors.size() << "  - first det: " << theDetectors.at(0).name() << std::endl ;

  if( theDetectors.size()  != 1 ){

    std::stringstream es ;
    es << " getExtension: selection is not unique (or empty)  includeFlag: " << DD4hep::DetType( includeFlag ) << " excludeFlag: " << DD4hep::DetType( excludeFlag )
       << " --- found detectors : " ;
    for( unsigned i=0, N= theDetectors.size(); i<N ; ++i ){
      es << theDetectors.at(i).name() << ", " ;
    }
    throw std::runtime_error( es.str() ) ;
  }

  theExtension = theDetectors.at(0).extension<DD4hep::DDRec::LayeredCalorimeterData>();

  return theExtension;
}

void RealisticCaloRecoSilicon::init() {
  RealisticCaloReco::init();
  _caloData=0;
}

void RealisticCaloRecoSilicon::getGeometryInformation() {

  if ( _cht_caloid != CHT::ecal ) {
    streamlog_out (ERROR) << "seems not no be an ECAL collection? " << _cht_caloid << std::endl;
    assert(0);
  }

  if ( _cht_layout == CHT::barrel ) {
    try {
      _caloData = getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::ELECTROMAGNETIC | DD4hep::DetType::BARREL), 
			       ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) );
    } catch(std::exception &e) {
      streamlog_out (ERROR) << "Could not get barrel ECAL parameters from DD4hep!" << std::endl;
    }
  } else if ( _cht_layout == CHT::endcap ) {
    try {
      _caloData = getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::ELECTROMAGNETIC | DD4hep::DetType::ENDCAP), 
				( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) );
    } catch(std::exception &e) {
      streamlog_out (ERROR) << "Could not get endcap ECAL parameters from DD4hep!" << std::endl;
    }
  } else if ( _cht_layout == CHT::ring ) {
    try {
      _caloData = getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::ELECTROMAGNETIC | DD4hep::DetType::ENDCAP | DD4hep::DetType::AUXILIARY ),
			       ( DD4hep::DetType::FORWARD ) );
    } catch(std::exception &e) {
      streamlog_out (ERROR) << "Could not get ring ECAL parameters from DD4hep!" << std::endl;
    }
  }

  assert ( _caloData );

}


float RealisticCaloRecoSilicon::reconstructEnergy(const CalorimeterHit* hit) {
  // here the input energy should be in MIPs
  float energy = hit->getEnergy();
  // what layer is this hit in?
  int layer   = (*_idDecoder) (hit)[_cellIDLayerString];
  // now correct for sampling fraction
  energy *= getLayerCalib( layer );
  return energy;
}

void RealisticCaloRecoSilicon::resetGaps() {
  for (int im=0; im<MAXMOD; im++) {
    for (int is=0; is<MAXSTA; is++) {
      for (int it=0; it<MAXTOW; it++) {
	for (int il=0; il<MAXLAY; il++) {
	  for (int iw=0; iw<MAXWAF; iw++) {
	    edgeHitsModStaTowWaf[im][is][it][il][iw].clear();
	  }
	}
      }
    }
  }
  return;
}

bool RealisticCaloRecoSilicon::cellAtWaferEdge(const CalorimeterHit* hit) {
  // cell at edge of sensor, but not at corner
  // this needs to be implemented once geometry and cellIDs are fixed...
  return false;
}

bool RealisticCaloRecoSilicon::cellAtWaferCorner(const CalorimeterHit* hit) {
  // cell at corner of sensor
  // this needs to be implemented once geometry and cellIDs are fixed...
  return false;
}

void RealisticCaloRecoSilicon::prepareForGaps(const CalorimeterHit* hit) {
  // store hits at wafer edge in this event
  // store by module/stave/etx to make life easier later
  bool corner = cellAtWaferCorner(hit);
  bool edge = cellAtWaferEdge(hit);
  if ( corner || edge ) {
    int module  = (*_idDecoder) (hit)[_cellIDModuleString];
    int stave   = (*_idDecoder) (hit)[_cellIDStaveString];
    int tower   = (*_idDecoder) (hit)[_cellIDTowerString];
    int layer   = (*_idDecoder) (hit)[_cellIDLayerString];
    int wafer   = (*_idDecoder) (hit)[_cellIDWaferString];
    edgeHitsModStaTowWaf[module][stave][tower][layer][wafer].push_back(hit);
  }
  return;
}

std::vector < const CalorimeterHit* > RealisticCaloRecoSilicon::getGapNeighbours( const CalorimeterHit* mainHit, int stave, int layer ) {
  std::vector < const CalorimeterHit* > nns;
  // to be implemented: needs knowledge of geometry in simulation.
  return nns;
}


void RealisticCaloRecoSilicon::fillGaps() {
  // Gaps between modules, towers, and wafers
  // if two hits are directly opposite each other across a gap, increase their energy
  // to be implemented
  return;
}

