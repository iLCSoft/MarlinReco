#ifndef Calibration_h
#define Calibration_h 1

//C++
#include "iostream"
#include "string"

//LCIO
#include "lcio.h"
#include "UTIL/LCFixedObject.h"


#define CalibrationNINTVals 4  // N event and N hits in it
#define CalibrationNFLOATVals 0 
#define CalibrationNDOUBLEVals 4 // Energies

class Calibration : public UTIL::LCFixedObject<CalibrationNINTVals,
  CalibrationNFLOATVals,CalibrationNDOUBLEVals> {
  
public: 
  
    /** Convenient constructor.
     */
  Calibration(int nevt,  int n1,  int n2,  int n3, 
	      double en1, double  en2,double  en3, double enr );

  /** 'Copy constructor' needed to interpret LCCollection read from file/database.
   */
  Calibration(EVENT::LCObject* obj) : UTIL::LCFixedObject<CalibrationNINTVals,
							 CalibrationNFLOATVals,
							 CalibrationNDOUBLEVals>(obj) { } 
    
    /** Important for memory handling*/
  virtual ~Calibration(){};
  
  // the class interface:
  int getNEvt();
  int getNhit1();
  int getNhit2();
  int getNhit3();
  double getEnr1();
  double getEnr2();
  double getEnr3();
  double getEnr4();

    // -------- need to implement abstract methods from LCGenericObject
    const std::string getTypeName() const { 
	return std::string("Calibration");
    } 
    const std::string getDataDescription() const {
	return std::string("i:Nevent,Nhits[3],d:Energies[4]"); 
    }
    
}; // class


#endif
