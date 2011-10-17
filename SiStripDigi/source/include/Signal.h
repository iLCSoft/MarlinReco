#ifndef SIGNAL_H
#define SIGNAL_H

// Include LCIO header files
#include <lcio.h>
#include <map>
#include <EVENT/SimTrackerHit.h>

namespace sistrip {

// Typedefs
typedef std::map<EVENT::SimTrackerHit *, float>  SimTrackerHitMap; // Hit, weight

//!
//! Signal class holds all information that one gets from strip, pixel, etc ...
//!
//! @author Z. Drasal, Charles University Prague
//!
class Signal {

 public:

//!Constructor
   Signal(double charge, double time) : _charge(charge), _time(time) {;}

//!Destructor
   ~Signal() {};


// SET METHODS

//!Set signal
   inline void setCharge(double charge) {_charge = charge;}

//!Update signal
   inline void updateCharge(double charge) {_charge += charge;}

//!Set time when signal was created
   inline void setTime(double time) {_time = time;}

//!Update MC truth information about SimTrackerHits, which contributed
   void updateSimHitMap(EVENT::SimTrackerHit * simHit, float weight);

//!Update MC truth information about SimTrackerHits, which contributed
   void updateSimHitMap(SimTrackerHitMap simHitMap);


// GET METHODS

//!Get signal
   inline double getCharge() const {return _charge;}

//!Get time when signal was created
   inline double getTime() const {return _time;}

//!Get MC truth information about SimTrackerHits, which contributed
   inline const SimTrackerHitMap & getSimHitMap() const {return _simHitMap;}

//!Get MC truth information - total sum of individual weights
   float getSimHitWeightSum();


 private:

   double           _charge;    //!< Signal charge
   double           _time;      //!< Time when signal has been created

   SimTrackerHitMap _simHitMap; //!< Map of SimTrkHits which contributed to the signal

}; // Class

} // Namespace

#endif // SIGNAL_H
