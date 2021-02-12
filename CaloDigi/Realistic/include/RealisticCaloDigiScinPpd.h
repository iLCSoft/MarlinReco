#ifndef DIGITIZER_DDCCALODIGISCINT_H
#define DIGITIZER_DDCCALODIGISCINT_H 1

#include "RealisticCaloDigi.h"

using namespace lcio ;
using namespace marlin ;

/** === RealisticCaloDigiScinPpd Processor === <br>
    realistic digitisation of scint+PPD (SiPM, MPPC) calorimeter hits
    D.Jeans 02/2016.
*/

class RealisticCaloDigiScinPpd : public RealisticCaloDigi {
  
 public:
  virtual Processor*  newProcessor() { return new RealisticCaloDigiScinPpd ; }
  RealisticCaloDigiScinPpd() ;

 protected:
  int getMyUnit() const {return NPE;}
  float digitiseDetectorEnergy(float energy) const ; // apply scin+PPD specific effects
  float convertEnergy( float energy, int inputUnit ) const; // convert energy from input to output scale

  float _PPD_pe_per_mip{};         // # photoelectrons/MIP for PPD
  int   _PPD_n_pixels{};           // # pixels in PPD
  float _misCalibNpix{};           // miscalibration of # PPD pixels
  float _pixSpread{};              // relative spread of PPD pixel signal

} ;

#endif



