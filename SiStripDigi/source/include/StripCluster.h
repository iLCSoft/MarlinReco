#ifndef STRIPCLUSTER_H
#define STRIPCLUSTER_H 1

#include "Signal.h"

// Include CLHEP header files
#include <CLHEP/Vector/ThreeVector.h>

namespace sistrip {

//! This class holds all information about strip clusters, where the strip
//! cluster is defined as a bunch of strips, where at least one strip has its
//! signal above so-called seed threshold and other strips above threshold lower
//! than seed.
//!
//! @author Z. Drasal, Charles University Prague
//!

class StripCluster 
{
	public:
		//!Constructor
		StripCluster( short int layerID  , short int ladderID , short int sensorID,
				CLHEP::Hep3Vector position, CLHEP::Hep3Vector posSigma, 
				double charge,  short int size) :
			_iLayer(layerID)   , _iLadder(ladderID) , _iSensor(sensorID),
			_position(position), _posSigma(posSigma), _charge(charge), _size(size)
	{ _stripFront = 0; _stripRear  = 0;}     
		//!Constructor with strips ID
		StripCluster( short int layerID  , short int ladderID , short int sensorID,
				CLHEP::Hep3Vector position, CLHEP::Hep3Vector posSigma, 
				double charge, short int size, 
				short int stripFront, short int stripRear) :
			_iLayer(layerID)   , _iLadder(ladderID) , _iSensor(sensorID),
			_position(position), _posSigma(posSigma), _charge(charge), 
			_size(size), _stripFront(stripFront), _stripRear(stripRear)
	{;}
		//!Destructor
		~StripCluster();
		
		// SET METHODS
		
		//!Set cluster layer ID
		inline void setLayerID(short int iLayer) {_iLayer = iLayer;}

		//!Set cluster ladder ID
		inline void setLadderID(short int iLadder) {_iLadder = iLadder;}
		
		//!Set cluster sensor ID
		inline void setSensorID(short int iSensor) {_iSensor = iSensor;}
		
		//!Set cluster position X in cm
		inline void setPosX( double posX) {_position.setX(posX);}
		
		//!Set cluster position Y in cm
		inline void setPosY( double posY) {_position.setY(posY);}
		
		//!Set cluster position Z in cm
		inline void setPosZ( double posZ) {_position.setZ(posZ);}
		
		//!Set cluster position Three vector in cm
		void set3Position( const CLHEP::Hep3Vector & position);
		
		//!Set cluster - position sigma X in cm
		inline void setPosSigmaX( double posSigmaX) {_posSigma.setX(posSigmaX);}
		
		//!Set cluster - position sigma Y in cm
		inline void setPosSigmaY( double posSigmaY) {_posSigma.setY(posSigmaY);}
		
		//!Set cluster - position sigma Z in cm
		inline void setPosSigmaZ( double posSigmaZ) {_posSigma.setZ(posSigmaZ);}
		
		//!Set cluster position Three vector in cm
		void set3PosSigma( const CLHEP::Hep3Vector & position);
		
		//!Set cluster charge
		inline void setCharge( double charge) {_charge = charge;}
		
		//!Set time when the cluster has been created by a particle in s
		inline void setTime( double time) {_time = time;}
		
		//!Set cluster size (how many strips contributed)
		inline void setSize( short int size) {_size = size;}
		
		//!Update MC truth information about SimTrackerHits, which contributed
		void updateSimHitMap(SimTrackerHitMap simHitMap);
		
		//!Set front strip ID
		inline void setStripFront(const int & stripId) {_stripFront = stripId;}
		
		//!Set rear strip ID
		inline void setStripRear(const int & stripId) {_stripRear = stripId;}
		
		// GET METHODS

//!Get cluster layer ID
	inline short int getLayerID() const {return _iLayer;}

//!Get cluster ladder ID
	inline short int getLadderID() const {return _iLadder;}

//!Get cluster sensor ID
	inline short int getSensorID() const {return _iSensor;}

//!Get cluster position X
   inline double getPosX() const {return _position.getX();}

//!Get cluster position Y
   inline double getPosY() const {return _position.getY();}

//!Get cluster position Z
   inline double getPosZ() const {return _position.getZ();}

//!Get cluster position Three vector
   inline CLHEP::Hep3Vector get3Position() const {return _position;}

//!Get cluster - position sigma X
   inline double getPosSigmaX() const {return _posSigma.getX();}

//!Get cluster - position sigma Y
   inline double getPosSigmaY() const {return _posSigma.getY();}

//!Get cluster - position sigma Z
   inline double getPosSigmaZ() const {return _posSigma.getZ();}

//!Get cluster position Three vector
   inline CLHEP::Hep3Vector get3PosSigma() const {return _posSigma;}

//!Get cluster charge
   inline double getCharge() const {return _charge;}

//!Get time when the cluster has been created by a particle
	inline double getTime() const {return _time;}

//!Get cluster size
	inline short int getSize() const {return _size;}

//!Get MC truth information about SimTrackerHits, which contributed
   inline const SimTrackerHitMap & getSimHitMap() const {return _simHitMap;}

//!Get MC truth information - total sum of individual weights
   float getSimHitWeightSum();
   
   		//!Get front strip ID
		inline int getStripFront() const { return _stripFront;}
		
		//!Get rear strip ID
		inline int getStripRear() const { return _stripRear;}


 protected:

   short int  _iLayer;          //!<ID number of a layer
   short int  _iLadder;         //!<ID number of a ladder
   short int  _iSensor;         //!<ID number of a sensor

   CLHEP::Hep3Vector _position; //!<Cluster position in cm
   CLHEP::Hep3Vector _posSigma; //!<Cluster - position sigma in cm

   double _charge;              //!<Total charge deposited
   double _time;                //!<Time when the cluster has been created by a particle
   short int _size;             //!<Cluster size

   short int _stripFront;       //!<Strip ID of the front sensor
   short int _stripRear;        //!<Strip ID of the rear sensor

   SimTrackerHitMap _simHitMap; //!< Map of SimTrkHits which contributed to the signal

}; // Class

} // Namespace

#endif // STRIPCLUSTER_H
