// -*- C++ -*-
//
// Package:    SiStripDigi
// Class:      SiStripGeomBuilder
// 
/**\class SiStripGeomBuilder SiStripBuilder/include/SiStripGeomBuilder.h 

//! Description: Builder to create the different subdetectors FTD, SET, VXT,...
//!              (SiStripGeom concrete instances).

//! Implementation: The Builder is called by the SiStripDigi processor acting 
                 as a client. The class calls the instances of SiStripGeom
   	         selected by the user via the steering xml file: FTD, VXD, SET....
	         Return a generic SiStripGeom.     
*/
//
// Original Author: Jordi Duarte Campderros  
//         Created:  Wed Jul  13 16:18:11 CET 2011
// 
// @autor J. Duarte Campderros, IFCA, Spain,  duarte@ifca.unican.es
//
//

#ifndef SiStripGeomBuilder_HH
#define SiStripGeomBuilder_HH

#include<string>

/**
\addtogroup SiStripDigi SiStripDigi
@{
*/

namespace sistrip
{

class SiStripGeom;

class SiStripGeomBuilder 
{

  public:
     static SiStripGeom * Build(const std::string & detector,const double & pitchFront,
		     const double & pitchRear);
 
};

}

/** @} */

#endif
