#ifndef SISTRIPGEOMFTD_H
#define SISTRIPGEOMFTD_H 1

#include<vector>
#include<string>

// Include CLHEP header files
// Include CLHEP header files

#include "SiStripGeom.h"

// Include Gear header files
#include "gear/FTDParameters.h"
#include "gear/FTDLayerLayout.h"

// Include LCIO stuff
#include "EVENT/SimTrackerHit.h"
#include "UTIL/BitField64.h"

//class gear::FTDParameters;
//class gear::FTDLayerLayout;
/**
\addtogroup SiStripDigi SiStripDigi
@{
*/
namespace sistrip 
{

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
//! Thu Jul 14 (J. Duarte Campderros)
//! Converted to abstract class used by the builder to construct the 
//! different subdetectors which going to use silicon strips (FTD,SIT,..)

enum LayerTypeFTD { strip = sistrip::stripF };

class SiStripGeomFTD: public SiStripGeom
{
	public:
		
		//!Constructor
		SiStripGeomFTD(const std::string &, const double & pitchFront, 
				const double & pitchRear);
		//!Destructor
		virtual ~SiStripGeomFTD();
		
		//!Method initializing class - reads Gear parameters from XML file
		virtual void initGearParams();

		//! Returns the input cellID0 where the field sensor is put to 0
		virtual int cellID0withSensor0(const int & cellID0) const;

		//!Method to extract codification ID, 
		//std::map<std::string,short int> cellIDDecProv(EVENT::SimTrackerHit * & simHit);
		virtual std::map<std::string,int> decodeCellID(const UTIL::BitField64 & cellID) const;
		virtual std::map<std::string,int> decodeCellID(const int & cellID) const;

		//!Method to extract strip and strip type (cellID1)
		virtual std::pair<StripType,int> decodeStripID(const int & encodedStripID) const; 
		virtual std::pair<StripType,int> decodeStripID(const UTIL::BitField64 & cellEnc) const; 
		
		//!Stores the cellID0 and cellID1 of the LCIO object to the file
		virtual void updateCanonicalCellID(const int & cellID, const int & stripType,
				const int & stripID, UTIL::BitField64 * bf);


		// TRANSFORMATION METHODS - GLOBAL REF. SYSTEM
		//!Transform given point from global ref. system to local ref. system (sensor)
		virtual CLHEP::Hep3Vector transformPointToLocal(short int layerID, 
				short int ladderID, short int sensorID, const CLHEP::Hep3Vector & point);

		//!Transform given vector from global ref. system to local ref. system (sensor)
		virtual CLHEP::Hep3Vector transformVecToLocal(short int layerID, 
				short int ladderID, short int sensorID,
				const CLHEP::Hep3Vector & vec);

		//!Transform given matrix from global ref. system to local ref. system (sensor)
		virtual CLHEP::HepMatrix transformMatxToLocal(short int layerID, 
				short int ladderID, const CLHEP::HepMatrix & matrix);
		
		
		// TRANSFORMATION METHODS - LOCAL REF. SYSTEM
		//!Transform given point from local ref. system (sensor) to global ref. system
		virtual CLHEP::Hep3Vector transformPointToGlobal(short int layerID, 
				short int ladderID, short int sensorID, const CLHEP::Hep3Vector & point);

		//!Transform given vector from local ref. system (sensor) to global ref. system
		virtual CLHEP::Hep3Vector transformVecToGlobal(short int layerID, 
				short int ladderID, short int sensorID,
				const CLHEP::Hep3Vector & vec);

		//!Transform given diagonal matrix from local ref. system (sensor) to global ref. system
		virtual CLHEP::HepMatrix transformMatxToGlobal(short int layerID, 
				short int ladderID, const CLHEP::HepMatrix & matrix);

		
		// OTHER METHODS - GLOBAL REF. SYSTEM
		//!Get info whether the given point is inside of Si sensor
		virtual bool isPointInsideSensor (short int layerID, 
				short int ladderID, short int sensorID, const CLHEP::Hep3Vector & point) const;

		
		// OTHER METHODS - LOCAL REF. SYSTEM
		//! Get info whether the given point is out of Si sensor (parameters: layerID,
		//! space point in local ref. system)
		virtual bool isPointOutOfSensor(short int layerID, const CLHEP::Hep3Vector & point) const;
		
		//! Get sensor pitch for strips in the local ref. system.
		//! Currently it means the faced to IP sensors of the petals.
		virtual double getSensorPitch(const int & diskID, const int & sensorID,
				const double & posZ) const;
		
		//! Get Y-centered position in the local ref. system for the stripID
		virtual double getStripPosY(const int & diskID, const int & sensorID, 
				const int & stripID, const double & posZ) const;

		//! Get number of strips per sensor
		virtual int getSensorNStrips(const int & diskID, const int & sensorID) const;
		
		//! Transforming a given point to the local ref. frame of a petal 
		//! which is rotated around its center an angle stAngle
		virtual CLHEP::Hep3Vector transformPointToRotatedLocal(const int & diskID, 
				const int & sensorID, const CLHEP::Hep3Vector & point) const;
		
		//! Transforming a given point from the reference system rotated around
		//! the center of a petal to its local ref. frame (inverse of 
		//! transformPointToRotatedLocal)
		virtual CLHEP::Hep3Vector transformPointFromRotatedLocal(const int & diskID, 
				const int & sensorID, const CLHEP::Hep3Vector & point) const;

		//! Get director vector of a strip (inside the local Ref. system)
		//! The vector is defined to describe the strip from z_local=0
		virtual CLHEP::Hep3Vector getStripUnitVector(const int & diskID, 
				const int & sensorID, const int & stripID) const;
		
		//! Get the point which crosses two strips
		virtual CLHEP::Hep3Vector getCrossLinePoint(const int & diskID, 
				const int & sensorID,
				const int & stripIDFront, const int & stripIDRear) const;
		
		// OTHER METHODS - IDENTIFYING
		//! Get strip ID (in Phi), points are given in local ref. system; strips are
		//! considered to be perpendicular to beam axis. It is perpedicular to the 
		//! RPhi LRF and paralel to the z-axis of the LRF. Recall that the number
		//! of strips in the RPhi LRF plane is dependent of Z (in LRF) due to 
		//! the trapezoid shape of the sensors.
		//! That sensor are the faced to IP on each petal (in the current design)
		virtual int getStripID(const int & diskID, const int & sensorID,
				const double & posRPhi, const double & posZ ) const;

		//!Method impossed by consistent with the mother class.
		//!The FTD sensors are all perpendicular to the beam axis
		virtual double getLadderTheta(short int diskID) const { return M_PI/2.0;}
		

		// PRINT METHODS
		//!Method printing general Gear parameters
		virtual void printGearParams() const;

		//!Method printing sensor Gear parameters
		virtual void printSensorParams(short int layerID) const;

	private:
		gear::FTDParameters * _ftdParams;
		gear::FTDLayerLayout * _ftdLayer;
		
		// Pitch of the front and rear sensor in the middle of 
		// the sensor. FIXME: These data member are initialized by
		// the processors  (Digi and Cluster) and could introduced
		// an inconsistency between differents values...
		double _pitchFront;
		double _pitchRear;
		//FIXME PROVISIONAL?
		double getLadderOffsetX(const short int & layerID) const;
		//! Get the stereo angle of the strips 
		double getStereoAngle(const int & diskID,const int & sensorID) const;

		std::vector<double> _layerOuterRadius;
		std::vector<double> _layerPetalOpAngle;
		std::vector<double> _ladderZOffsetSign0;
	
}; // Class

} // Namespace

/** @} */

#endif // SISTRIPGEOMFTD_H
