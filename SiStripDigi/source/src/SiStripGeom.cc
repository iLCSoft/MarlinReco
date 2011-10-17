#include "SiStripGeom.h"

// Include standard header files
#include "Colours.h"
#include "PhysicalConstants.h"

#include <cstdlib>
#include <iomanip>
#include <algorithm>

// Include CLHEP header files
#include <CLHEP/Vector/EulerAngles.h>
#include <CLHEP/Vector/Rotation.h>

// Include Gear header files
//#include <gear/BField.h>
//#include <gear/VXDParameters.h>
//#include <gear/VXDLayerLayout.h>
#include <gearimpl/Vector3D.h>

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Namespaces
//using namespace CLHEP;
//using namespace marlin;

namespace sistrip 
{
	SiStripGeom::SiStripGeom(const std::string & subdetector)
		:_gearType(subdetector)
	{

	}
	
	SiStripGeom::~SiStripGeom()
	{
	//	std::cout << "Deleting SiStripGeom" << std::endl;
	}
	
// GEOMETRY PROPERTIES

//
// Encode cell ID
//
/*int SiStripGeom::encodeCellID(short int layerID, short int ladderID, short int sensorID) const
{
	return layerID*LAYERCOD + ladderID*LADDERCOD + sensorID*SENSORCOD;
}*/

//
// Decode cell ID
//
/*void SiStripGeom::decodeCellID(short int & layerID, short int & ladderID, short int & sensorID, int cellID) const
{
	layerID  =  cellID / LAYERCOD;
	ladderID = (cellID - layerID*LAYERCOD) / LADDERCOD;
	sensorID = (cellID - layerID*LAYERCOD - ladderID*LADDERCOD) / SENSORCOD;
}*/

//
// Encode strip ID
//
int SiStripGeom::encodeStripID(StripType type, int stripID) const
{
	// TODO: Posible mejora del algoritmo
	int encodedStripID = 0;
	
	// Strip in R-Phi
	if (type == RPhi) 
	{
		encodedStripID = (stripID + STRIPOFF)*STRIPCODRPHI;
	}
	
	// Strip in Z
	if (type == Z) 
	{
		encodedStripID = (stripID + STRIPOFF)*STRIPCODZ;
	}

	return encodedStripID;
}

//
// Decode strip ID: FIXME TO BE FIXED---> Now the codification is done with the cellDec
//
/*std::pair<StripType,int> SiStripGeom::decodeStripID(const int & encodedStripID) const 
{
	StripType stype;
	int stripID;

	// Strip in RPhi
	if(encodedStripID/STRIPCODRPHI > 0.) 
	{
		stripID = encodedStripID/STRIPCODRPHI - STRIPOFF;
		stype    = RPhi;
		
	}
	// Strip in Z
	else if(encodedStripID/STRIPCODZ > 0.) 
	{
		stripID = encodedStripID/STRIPCODZ - STRIPOFF;
		stype    = Z;
		
	}
	else
	{
		// Error
		streamlog_out(ERROR) << "SiStripGeom::decodeStripID: "
			<< encodedStripID
			<< " - problem to identify if strips in Z or R-Phi!!!"
			<< std::endl;
		exit(0);
	}

	return std::pair<StripType,int>(stype,stripID);
}*/

//
// Get layer real ID
//
int SiStripGeom::getLayerRealID(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if(_layerRealID.size()>(unsigned short int)layerID) 
	{
		return _layerRealID[layerID];
	}
	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getLayerRealID - layerID: " << layerID << " out of range!!!"
			<< std::endl;
		exit(-1);
	}
}

//
// Transform real layer ID to C-type numbering 0 - n ...
//
short int SiStripGeom::getLayerIDCTypeNo(int realLayerID) const
{
	for (unsigned int i=0; i<_layerRealID.size(); i++) 
	{
		if(realLayerID == _layerRealID[i]) 
		{
			return i;
		}
	}
	streamlog_out(ERROR) << "SiStripGeom::getLayerIDCTypeNo - layer: " << 
		realLayerID << " not found in layer list!!!" << std::endl;
	exit(-1);
}

//
// Get layer type
//
short int SiStripGeom::getLayerType(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_layerType.size()>(unsigned short int)layerID)
	{
		return _layerType[layerID];
	}
	else 
	{
		//FIXME---> _layerType is empty!!! FIXME
		std::cout << " size " << _layerType.size() << std::endl;
		streamlog_out(ERROR) << "SiStripGeom::getLayerType - layerID: " 
			<< layerID << " out of range!!!" << std::endl;
		exit(0);
	}
}

//
// Get layer radius
//
double SiStripGeom::getLayerRadius(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_layerRadius.size()>(unsigned short int)layerID) 
	{
		return _layerRadius[layerID];
	}
	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getLayerRadius - layerID: " << layerID << " out of range!!!"
			<< std::endl;
		exit(-1);
	}
}

//
// Get semiangle layer (petals)
//
double SiStripGeom::getLayerHalfPhi(const int & layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_layerHalfPhi.size()>(unsigned short int)layerID)
	{
		return _layerHalfPhi[layerID];
	}
	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getLayerHalfPhi - layerID: " << layerID << " out of range!!!"
			<< std::endl;
		exit(0);
	}
}

//
// Get layer phi zero angle
//
double SiStripGeom::getLayerPhi0(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_layerPhi0.size()>(unsigned short int)layerID)
	{
		return _layerPhi0[layerID];
	}
	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getLayerPhi0 - layerID: " << layerID << " out of range!!!"
			<< std::endl;
		exit(0);
	}
}

//
// Get number of ladders
//
short int SiStripGeom::getNLadders(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_numberOfLadders.size()>(unsigned short int)layerID)
	{
		return _numberOfLadders[layerID];
	}
	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getNLadders - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder thickness
//
double SiStripGeom::getLadderThick(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_ladderThick.size()>(unsigned short int)layerID)
	{
		return _ladderThick[layerID];
	}
   	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getLadderThick - layerID: " << layerID << " out of range!!!"
			<< std::endl;
		exit(-1);
	}
}

// FIXME::TODO: queda pendiente desde aqui hacia abajo

//
// Get ladder width
//
double SiStripGeom::getLadderWidth(short int layerID) const
{
   if (_ladderWidth.size()>(unsigned short int)layerID) return _ladderWidth[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderWidth - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder length
//
double SiStripGeom::getLadderLength(short int layerID) const
{
   if (_ladderLength.size()>(unsigned short int)layerID) return _ladderLength[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderLength - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder offset in Y
//
double SiStripGeom::getLadderOffsetY(short int layerID) const
{
   if (_ladderOffsetY.size()>(unsigned short int)layerID) return _ladderOffsetY[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderOffsetY - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder offset in Z
//
double SiStripGeom::getLadderOffsetZ(short int layerID) const
{
   if (_ladderOffsetZ.size()>(unsigned short int)layerID) return _ladderOffsetZ[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderOffsetZ - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder rotation - phi angle (in system of units defined in PhysicalConstants.h)
//
double SiStripGeom::getLadderPhi(short int layerID, short int ladderID) const
{
   if (ladderID<getNLadders(layerID)) return (getLayerPhi0(layerID) + 2*pi/getNLadders(layerID)*ladderID);
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderPhi - ladderID: " << ladderID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder rotation - theta angle
//
/*double SiStripGeom::getLadderTheta(short int layerID) const
{
   if (_layerTheta.size()>(unsigned short int)layerID) return _layerTheta[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderTheta - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}*/

//
// Get number of sensors for given ladder
//
short int SiStripGeom::getNSensors(short int layerID) const
{
   if (_numberOfSensors.size()>(unsigned short int)layerID) return _numberOfSensors[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getNSensors - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get number of strips (in each sensor)
//
/*int SiStripGeom::getSensorNStrips(const int & layerID, const int & sensorID) const
{
   if(_sensorNStrips.size()>(unsigned short int)layerID)
   {
	   return _sensorNStrips[layerID];
   }
   else
   {
	   return 0;
   }
}

//
// Get number of RPhi strips (in each sensor)
//
int SiStripGeom::getSensorNStripsInRPhi(short int layerID) const
{
	if(_sensorNStripsInRPhi.size()>(unsigned short int)layerID)
	{
		return _sensorNStripsInRPhi[layerID];
	}
	else
	{
		return 0;
	}
}*/

//
// Get sensor pitch in Z axis for barrel-type and forward-type sensors
//
/*double SiStripGeom::getSensorPitchInZ(short int layerID) const
{
   if(_sensorPitchInZ.size()>(unsigned short int)layerID)
   {
	   return _sensorPitchInZ[layerID];
   }
   else
   {
	   return 0.;
   }
}*/


//
// Get sensor thickness
//
double SiStripGeom::getSensorThick(short int layerID) const
{
   if (_sensorThick.size()>(unsigned short int)layerID) return _sensorThick[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getSensorThick - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get sensor width
//
double SiStripGeom::getSensorWidthMax(short int layerID) const
{
	if(_sensorWidth.size()>(unsigned short int)layerID)
	{
		return _sensorWidth[layerID];
	}
	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getSensorWidth - layerID: " 
			<< layerID << " out of range!!!" << std::endl;
		exit(0);
	}
}

//
// Get sensor width2
//
double SiStripGeom::getSensorWidthMin(short int layerID) const
{
	if(_sensorWidth2.size()>(unsigned short int)layerID)
	{
		return _sensorWidth2[layerID];
	}
	else
	{
		return getSensorWidthMax(layerID);
	}
}

//
// Get sensor length
//
double SiStripGeom::getSensorLength(short int layerID) const
{
   if (_sensorLength.size()>(unsigned short int)layerID) return _sensorLength[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getSensorLength - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get gap size inbetween sensors
//
double SiStripGeom::getSensorGapInBetween(short int layerID) const
{
   if (_sensorGapInBtw.size()>(unsigned short int)layerID) return _sensorGapInBtw[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getSensorGapInBetween - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get width of sensor rim in Z (passive part of silicon)
//
double SiStripGeom::getSensorRimWidthInZ(short int layerID) const
{
  if (_sensorRimWidthInZ.size()>(unsigned short int)layerID) return _sensorRimWidthInZ[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getSensorRimWidthInZ - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   } 
}

//
// Get width of sensor rim in R-Phi (passive part of silicon)
//
double SiStripGeom::getSensorRimWidthInRPhi(short int layerID) const
{
   if (_sensorRimWidthInRPhi.size()>(unsigned short int)layerID) return _sensorRimWidthInRPhi[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getSensorRimWidthInRPhi - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}


// TRANSFORMATION METHODS - GLOBAL REF. SYSTEM

//
// Get Z-position of given strip in local ref system (in system of units defined in PhysicalConstants.h);
// strips are considered to be perpendicular to beam axis for both barrel-type and forward-type sensors.
//SUBDETECTOR DEPENDENT
/*double SiStripGeom::getStripPosInZ(short int layerID, int stripID) const
{
	// Get pitch
	double sensPitch = getSensorPitchInZ(layerID);
	
	// Calculate position
	double posZ = sensPitch*(stripID + 0.5);
	// Error
	if( (posZ<0.) || (posZ>getSensorLength(layerID)) ) 
	{
		streamlog_out(ERROR) 
			<< "SiStripGeom::getStripPosInZ - position out of sensor!!!"
			<< std::endl;
		exit(-1);
	}

	// Return Z position of given strip in local ref. system
	return posZ;
}

//
// Get R-Phi position of given strip in local ref system (in system of units defined in PhysicalConstants.h);
// strips are considered to be parallel to beam axis for barrel-type sensors and at angle alpha for forward-type
// sensors (see getSensorPitchInRPhi method).
//
double SiStripGeom::getStripPosInRPhi(short int layerID, int stripID, double posZ) const
{
	// Get pitch
	double sensPitch = getSensorPitchInRPhi(layerID, posZ);
	
	// Calculate position
	double posRPhi = sensPitch*(stripID + 0.5);
	
	// Error
	if ( (posRPhi<0.) || (posRPhi>getSensorWidth(layerID)) ) 
	{
		streamlog_out(ERROR) 
			<< "SiStripGeom::getStripPosInRPhi - position out of sensor!!!"
			<< std::endl;
		exit(-1);
	}
	
	// Return R-Phi position of given strip in local ref. system
	return posRPhi;
}

//
// Get strip ID (in Z-direction, i.e. measures the Rphi direction), 
// point is given in local ref. system
// 
//
int SiStripGeom::getStripIDInZ(short int layerID, double posZ ) const
{
	// Get pitch
	double sensPitch = getSensorPitchInZ(layerID);

	if (sensPitch == 0) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripIDInZ " 
			<< "- division by zero (sensPitch is zero)!!!"
			<< std::endl;
		exit(-1);
	}

	// Get number of strips
	int sensNStrips = getSensorNStripsInZ(layerID);

	int stripID;
	// Calculate stripID
	if (posZ <= 0.)
	{
		stripID = 0;
	}
	else 
	{
	   stripID = floor(posZ/sensPitch);
	   if (stripID >= sensNStrips)
	   {
		   stripID = sensNStrips - 1;
	   }
	}

   	// Error
	if (stripID >= sensNStrips) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripIDInZ " 
			<< "- stripID in Z greater than number of strips!!!"
			<< std::endl;
		exit(-1);
	}
	// Return stripID
	return stripID;
}
*/
//
// Get strip ID (in R-Phi), point is given in local ref. system; strips are
// considered to be parallel to beam axis for barrel-type sensors and at angle
// alpha for forward-type sensors (see getSensorPitchInRPhi method).
//
// IMPLEMENTATION DEPENDENT OF THE SUBDETECTOR

/*int SiStripGeom::getStripIDInRPhi(short int layerID, double posRPhi, double posZ ) const
{
	// Get pitch
	double sensPitch = getSensorPitchInRPhi(layerID, posZ);
	if (sensPitch == 0) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripIDInRPhi - division by zero (sensPitch is zero)!!!"
			<< std::endl;
		exit(-1);
	}
	// Get number of strips
	int sensNStrips = getSensorNStripsInRPhi(layerID);
	
   	int stripID;
	// Calculate stripID
	if (posRPhi <= 0.)
	{
		stripID = 0;
	}
	else 
	{
		stripID = floor(posRPhi/sensPitch);
		if (stripID >= sensNStrips) 
		{
			stripID = sensNStrips - 1;
		}
	}
	// Error
	if (stripID >= sensNStrips) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripIDInRPhi - stripID in RPhi greater than number of strips!!!"
			<< std::endl;
		exit(-1);
	}
	// Return stripID
	return stripID;
}*/


// PRINT METHODS

//
// Method printing general Gear parameters
//
/*
void SiStripGeom::printGearParams() const
{
   streamlog_out(MESSAGE3) << std::endl
                           << " "
	                   		<< DUNDERL
                           << DBLUE
                           << "Gear parameters:"
                           << ENDCOLOR
                           << " "
                           << std::endl  << std::endl;

   // Gear type: BField
   streamlog_out(MESSAGE3) << "  B field [T]:       "
                           << _magField/T
                           << std::endl << std::endl;

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Print general info
      for (int i=0; i<_numberOfLayers; ++i){

         streamlog_out(MESSAGE2) << std::endl;
         streamlog_out(MESSAGE2) << "  Layer: "           << _layerRealID[i]      << std::endl;
         if (_layerType[i]==pixel) {
         streamlog_out(MESSAGE2) << "  LayerType:       " << "pixel"              << std::endl;
         }
         else {
         streamlog_out(MESSAGE2) << "  LayerType:       " << "strip"              << std::endl;
         }
         streamlog_out(MESSAGE2) << "  NumberOfLadders: " << _numberOfLadders[i]  << std::endl;
         streamlog_out(MESSAGE2) << "  Radius[mm]:      " << _layerRadius[i]/mm   << std::endl;
         streamlog_out(MESSAGE2) << "  Width[mm]:       " << _ladderWidth[i]/mm   << std::endl;
         streamlog_out(MESSAGE2) << "  Length[mm]:      " << _ladderLength[i]/mm  << std::endl;
         streamlog_out(MESSAGE2) << "  Phi0:            " << _layerPhi0[i]        << std::endl;
         streamlog_out(MESSAGE2) << "  Theta:           " << _layerTheta[i]       << std::endl;
         streamlog_out(MESSAGE2) << "  OffsetY[mm]:     " << _ladderOffsetY[i]/mm << std::endl;
 //        streamlog_out(MESSAGE2) << "  OffsetZ[mm]:     " << _ladderOffsetZ[i]/mm << std::endl; OUT//
      }

   }
   // Gear type: unknown - error
   else {
   	streamlog_out(ERROR) << "Unknown gear type!"
   	                     << std::endl;

      exit(0);
   }
}
*/
//
// Method printing sensor Gear parameters (parameters: sensorID)
//
/*void SiStripGeom::printSensorParams(short int layerID) const
{

	// Gear type: VXD
   if (_gearType == "VXD") {

      // Print sensor parameters
      streamlog_out(MESSAGE2) << "    Parameters: " << _sensorThick[layerID]/um       << "um thick, "
                              << "with "            << _sensorPitchInZ[layerID]/um    << "um pitch "
                              << "and "             << _sensorNStripsInZ[layerID]     << " strips in Z"
                              << ", resp. "         << _sensorPitchInRPhi[layerID]/um << "um pitch "
                              << "and "             << _sensorNStripsInRPhi[layerID]  << " strips in R-Phi."
                              << std::endl;
   }
   // Gear type: unknown - error
   else {
   	streamlog_out(ERROR) << "Unknown gear type!"
   	                     << std::endl;

      exit(0);
   }
}
*/
} // Namespace;

