// -*- C++ -*-
//
// Package:    SiStripDigi
// Class:      SiStripGeomBuilder
//
//
// Original Author: Jordi Duarte Campderros
//         Created:  Wed Jul  13 16:23:11 CET 2011
//
// @autor J. Duarte Campderros, IFCA, Spain,  duarte@ifca.unican.es
//
//

// standard
#include <stdlib.h>
#include <string>

// Include Marlin
#include <streamlog/streamlog.h>

#include "SiStripGeomBuilder.h"
#include "SiStripGeomFTD.h"
// #include "SiStripGeomVXD.h"

namespace sistrip {

SiStripGeom* SiStripGeomBuilder::Build(const std::string& detector, const double& pitchFront, const double& pitchRear) {
  SiStripGeom* p = 0;

  // Building the concrete SiStripGeom
  if (detector == "FTD") {
    p = new SiStripGeomFTD("FTD", pitchFront, pitchRear);
  }
  // else if( detector == "VXD" )
  //{
  //	p = new SiStripGeomVXD();
  // }
  //  Add another type here and codify a new concrete
  //  class derived from SiStripGeom
  else {
    streamlog_out(ERROR) << "Detector type unknown: " << detector << std::endl;
    exit(-1);
  }

  return p;
}

} // namespace sistrip
