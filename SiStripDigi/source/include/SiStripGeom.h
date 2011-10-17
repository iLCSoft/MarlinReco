#ifndef SISTRIPGEOM_H
#define SISTRIPGEOM_H 1

#include<map>
#include<vector>
#include<string>


#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Vector/ThreeVector.h>
// Include Gear header files
//#include <gear/GearParameters.h>

#include "EVENT/SimTrackerHit.h"

#include "UTIL/BitField64.h"

namespace sistrip 
{

// Define constants
#define EPS 5           // in micrometers -- FIXME: Cambiado de 5 a 9 PROV
#define LAYERCOD   1000 // Const to encode and decode layers==disk
#define LADDERCOD    10 // Const to encode and decode ladders==petal
#define SENSORCOD     1 // Const to encode and decode sensors

#define STRIPOFF      1 // Const to encode and decode strip offset
#define STRIPCODRPHI  1 // Const to encode and decode strip in RPhi
#define STRIPCODZ    -1 // Const to encode and decode strip in Z

// Define strip type
enum StripType { RPhi = 0, Z = 1 };

// Define layer type
enum LayerType { pixel = 0, stripB = 1, stripF = 2};

//! Gear geometry class - holds all geometry information about silicon strip
//! sensors. The data are taken directly from Gear xml file and values are
//! saved in the system of units defined in PhysicalConstants.h. The local
//! reference system is defined as follows: X-axis is perpendicular to the
//! beam axis and to the ladder (sensor) plane; Y-axis lies in a ladder
//! (sensor) plane and is also perpendicular to the beam axis; Z-axis is
//! parallel to the beam axis(for zero theta); (0,0,0) point is positioned
//! such as X, Y, Z coordinates are always positive. Strips are considered
//! to be either perpendicular to the beam axis or parallel with the beam
//! axis (SSDs) or both (DSSDs).
//!
//! @author Z. Drasal, Charles University Prague
//!
//! Thu Jul 14 (J. Duarte)
//! Converted to abstract class used by the builder to construct the 
//! different subdetectors which going to use silicon strips (FTD,SIT,..)

class SiStripGeom
{
	public:
		
		//!Constructor
		SiStripGeom(const std::string & detector);
		//!Destructor
		virtual ~SiStripGeom();
		
		//!Method initializing class - reads Gear parameters from XML file
		//!Pure virtual method, have to be implemented by a concrete instance
		virtual void initGearParams() = 0;

		//!Pure virtual method, PROVISIONAL
		//virtual std::map<std::string,short int> cellIDDecProv(EVENT::SimTrackerHit * & simHit) = 0;

		// MAGNETIC FIELD
		//!Get magnetic field - x
		virtual double getBx() { return _magField.getX(); }
		//!Get magnetic field - y
		virtual double getBy() { return _magField.getY(); }
		//!Get magnetic field - z
		virtual double getBz() { return _magField.getZ(); }
		//!Get magnetic field - Three vector
		virtual CLHEP::Hep3Vector get3MagField() { return _magField; }

		// GEOMETRY PROPERTIES
		//!Get gear type
		virtual inline std::string getGearType() { return _gearType; }
		
		// ENCODING
		//! Returns the input cellID0 where the field sensor is put to 0
		virtual int cellID0withSensor0(const int & cellID0) const =0 ;

		//!Encode cellID
		//virtual int encodeCellID(short int layerID, short int ladderID, short int sensorID) const;
		//!Decode cellID1
		//virtual void decodeCellID(short int & layerID, short int & ladderID, short int & sensorID, int cellID) const;
		virtual std::map<std::string,int> decodeCellID(const UTIL::BitField64 & cellDec) const = 0;
		virtual std::map<std::string,int> decodeCellID(const int & cellDec) const = 0;
		
		//!Encode stripID
		virtual int encodeStripID(StripType type, int stripID) const;
		//!Decode stripID
		virtual std::pair<StripType,int> decodeStripID(const int & encStripID) const = 0;
		virtual std::pair<StripType,int> decodeStripID(const UTIL::BitField64 & cellDec) const = 0;
		//!Stores the cellID0 and cellID1 of the LCIO object to the file
		virtual void updateCanonicalCellID(const int & cellID, const int & stripType,
				const int & stripID, UTIL::BitField64 * bf) = 0;
		
		// LAYER PROPERTIES
		//!Get number of layers
		virtual short int getNLayers() const {return _numberOfLayers;}
		//!Get layer real ID
		virtual int getLayerRealID(short int layerID) const;
		
		//!Transform real layer ID to C-type numbering 0 - n ...
		virtual short int getLayerIDCTypeNo(int realLayerID) const;
		
		//!Get layer type
		virtual short int getLayerType(short int layerID) const;

		//!Get layer radius
		virtual double getLayerRadius(short int layerID) const;
	
		//!Get layer semiangle
		virtual double getLayerHalfPhi(const int & layerID) const;
		//!Get layer phi zero angle
		virtual double getLayerPhi0(short int layerID) const;

		
		// LADDER PROPERTIES
		//!Get number of ladders
		virtual short int getNLadders(short int layerID) const;
		
		//!Get ladder thickness
		virtual double getLadderThick(short int layerID) const;
		//!Get ladder width (the wider one for forward sensors)
		virtual double getLadderWidth(short int layerID) const;
		//!Get ladder length
		virtual double getLadderLength(short int layerID) const;
		//!Get ladder offset in Y
		virtual double getLadderOffsetY(short int layerID) const;
		//!Get ladder offset in Z
		virtual double getLadderOffsetZ(short int layerID) const;

		//!Get ladder rotation - phi angle
		virtual double getLadderPhi(short int layerID, short int ladderID) const;
		//!Get ladder rotation - theta angle
		virtual double getLadderTheta(short int layerID) const = 0;

		
		// SENSOR PROPERTIES
		//!Get number of sensors for given ladder
		virtual short int getNSensors(short int layerID) const;
		
		//!Get number of strips in Z axis (in each sensor)
/*		virtual int getSensorNStripsInZ(short int layerID) const;
		//!Get number of strips in R-Phi (in each sensor)
		virtual int getSensorNStripsInRPhi(short int layerID) const;
		//!Get sensor pitch in Z axis for barrel-type and forward-type sensors
		virtual double getSensorPitchInZ(short int layerID) const;
		
		//!Get sensor pitch in the R-Phi plane of the local coordinate system
		//!(posZ must be in local coordinate system)
		virtual double getSensorPitchInRPhi(short int layerID, double posZ) const = 0;*/
		//! Get number of strips
		virtual int getSensorNStrips(const int & layerID, const int & sensorID) const=0;
		//! Get sensor pitch 
		virtual double getSensorPitch(const int & layerID, const int & sensorID, 
				const double & posZ) const = 0;

		//!Get sensor thickness
		virtual double getSensorThick(short int layerID) const;
		//!Get sensor width (the wider one for forward-type sensors)
		virtual double getSensorWidthMax(short int layerID) const;
		//!Get sensor width 2 (the narrower one for forward-type sensors)
		virtual double getSensorWidthMin(short int layerID) const;
		//!Get sensor length
		virtual double getSensorLength(short int layerID) const;
		
		//!Get gap size inbetween sensors
		virtual double getSensorGapInBetween(short int layerID) const;

		//!Get width of sensor rim in Z (passive part of silicon)
		virtual double getSensorRimWidthInZ(short int layerID) const;
		
		//!Get width of sensor rim in R-Phi (passive part of silicon)
		virtual double getSensorRimWidthInRPhi(short int layerID) const;


		// TRANSFORMATION METHODS - GLOBAL REF. SYSTEM
		//!Transform given point from global ref. system to local ref. system (sensor)
		virtual CLHEP::Hep3Vector transformPointToLocal(short int layerID, 
				short int ladderID, short int sensorID, const CLHEP::Hep3Vector & point) = 0;

		//!Transform given vector from global ref. system to local ref. system (sensor)
		virtual CLHEP::Hep3Vector transformVecToLocal(short int layerID, 
				short int ladderID, short int sensorID,
				const CLHEP::Hep3Vector & vec) = 0;

		//!Transform given matrix from global ref. system to local ref. system (sensor)
		virtual CLHEP::HepMatrix transformMatxToLocal(short int layerID, 
				short int ladderID, const CLHEP::HepMatrix & matrix) = 0;
		
		
		// TRANSFORMATION METHODS - LOCAL REF. SYSTEM
		//!Transform given point from local ref. system (sensor) to global ref. system
		virtual CLHEP::Hep3Vector transformPointToGlobal(short int layerID, 
				short int ladderID, short int sensorID, const CLHEP::Hep3Vector & point) = 0;

		//!Transform given vector from local ref. system (sensor) to global ref. system
		virtual CLHEP::Hep3Vector transformVecToGlobal(short int layerID, 
				short int ladderID, short int sensorID, 
				const CLHEP::Hep3Vector & vec) = 0;

		//!Transform given diagonal matrix from local ref. system (sensor) to global ref. system
		virtual CLHEP::HepMatrix transformMatxToGlobal(short int layerID, 
				short int ladderID, const CLHEP::HepMatrix & matrix) = 0;

		
		// OTHER METHODS - GLOBAL REF. SYSTEM
		//!Get info whether the given point is inside of Si sensor
		virtual bool isPointInsideSensor (short int layerID, short int ladderID, 
				short int sensorID, const CLHEP::Hep3Vector & point) const = 0;


		// OTHER METHODS - LOCAL REF. SYSTEM
		//!Get info whether the given point is out of Si sensor
		virtual bool isPointOutOfSensor(short int layerID, 
				const CLHEP::Hep3Vector & point) const = 0;

		//!Get strip y-position
		virtual double getStripPosY(const int & layerID, const int & sensorID, 
				const int & stripID, const double & posZ) const =0;
		//!Get strip ID
		virtual int getStripID(const int & layerID, const int & sensorID, 
				const double & posRPhi, const double & posZ) const = 0;

		//! Transforming a given point to the local ref. frame of a petal 
		//! which is rotated around its center an angle stAngle
		virtual CLHEP::Hep3Vector transformPointToRotatedLocal(const int & diskID, 
				const int & sensorID, const CLHEP::Hep3Vector & point) const = 0;
		
		//! Transforming a given point from the reference system rotated around
		//! the center of a petal to its local ref. frame (inverse of 
		//! transformPointToRotatedLocal)
		virtual CLHEP::Hep3Vector transformPointFromRotatedLocal(const int & diskID, 
				const int & sensorID, const CLHEP::Hep3Vector & point) const = 0;

		//! Get director vector of a strip (inside the local Ref. system)
		//! The vector is defined to describe the strip from z_local=0
		virtual CLHEP::Hep3Vector getStripUnitVector(const int & diskID, 
				const int & sensorID, const int & stripID) const = 0;
		
		//! Get the point which crossed two strips
		virtual CLHEP::Hep3Vector getCrossLinePoint(const int & diskID, 
				const int & sensorID, 
				const int & stripIDFront, const int & stripIDRear) const = 0;

		// PRINT METHODS
		//!Method printing general Gear parameters
		virtual void printGearParams() const = 0;

		//!Method printing sensor Gear parameters
		virtual void printSensorParams(short int layerID) const = 0;
	
	protected:		
		std::string _gearType;  // Not necessary              //!< GearType
		
		// Magnetic field
		CLHEP::Hep3Vector _magField;
	
		// Geometry parameters - layers (disks)
		short int _numberOfLayers;
		
		std::vector<double> _layerHalfPhi;
		std::vector<double> _layerThickness;
		std::vector<double> _layerHalfThickness;
		std::vector<double> _layerRadius;
		std::vector<double> _layerZ; //Z-center NEW
	
		std::vector<int> _layerRealID;
		std::vector<int> _layerType;
	     
	     	std::vector<double> _layerPhi0;
	     	std::vector<double> _layerTheta;
		
		// Geometry parameters - ladders
		std::vector<int> _numberOfLadders;
		std::vector<int> _laddersInLayer;
		std::vector<double> _layerLadderLength;
		std::vector<double> _layerLadderWidth;
		std::vector<double> _layerLadderHalfWidth;
		std::vector<double> _layerActiveSiOffset;
		std::vector<double> _layerLadderGap;
		std::vector<double> _layerPhiOffset;
		std::vector<double> _ladderThick;
		std::vector<double> _ladderWidth;
		std::vector<double> _ladderLength;
		std::vector<double> _ladderOffsetY;
		std::vector<double> _ladderOffsetZ;
	 
		// Geometry - sensors
		std::vector<int> _numberOfSensors;
		std::vector<int> _sensorNStripsInFront;
		std::vector<int> _sensorNStripsInRear;
		std::vector<double> _sensorPitchInFront; 
		std::vector<double> _sensorPitchInRear;
		std::vector<double> _sensorThick;
		std::vector<double> _sensorWidth;
		std::vector<double> _sensorWidth2;
		std::vector<double> _sensorLength;
		std::vector<double> _sensorGapInBtw;
		std::vector<double> _sensorRimWidthInZ;
		std::vector<double> _sensorRimWidthInRPhi;
}; // Class

} // Namespace

#endif // SISTRIPGEOM_H
