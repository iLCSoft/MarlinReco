#ifndef DIGITIZER_REALISTICCALODIGISILICON_H
#define DIGITIZER_REALISTICCALODIGISILICON_H 1

#include "RealisticCaloDigi.h"

using namespace lcio ;
using namespace marlin ;

/** === RealisticCaloDigiSilicon Processor === <br>
    realistic digitisation of silicon calorimeter hits
    D.Jeans 02/2016.
*/

class RealisticCaloDigiSilicon : public RealisticCaloDigi {
  
 public:
  virtual Processor*  newProcessor() { return new RealisticCaloDigiSilicon ; }
  RealisticCaloDigiSilicon() ;

 protected:
  int getMyUnit() {return MIP;}
  float convertEnergy( float energy, int inputUnit ); // convert energy from input to output (MIP) scale 
  float digitiseDetectorEnergy(float energy);         // apply silicon-specific realistic digitisation
  float _ehEnergy{};                                    // energy to create e-h pair in silicon
} ;

#endif



