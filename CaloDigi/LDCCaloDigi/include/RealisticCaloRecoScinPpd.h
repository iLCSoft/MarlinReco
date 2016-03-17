#ifndef REALISTICCALORECOSCINPPD_H
#define REALISTICCALORECOSCINPPD_H 1

#include "RealisticCaloReco.h"

/** === RealisticCaloRecoSilicon Processor === <br>
    realistic reconstruction of scint+PPD calorimeter hits
    D.Jeans 02/2016.
*/

class RealisticCaloRecoScinPpd : public RealisticCaloReco {

 public:

  virtual Processor*  newProcessor() { return new RealisticCaloRecoScinPpd ; }

  RealisticCaloRecoScinPpd();

 protected:

  // virtual void init();

  virtual float reconstructEnergy(const CalorimeterHit* hit);

  // no gap correction for now
  virtual void resetGaps() {}
  virtual void prepareForGaps(const CalorimeterHit* hit) {}
  virtual void fillGaps() {}

  float _PPD_pe_per_mip;         // # photoelectrons/MIP for MPPC
  int   _PPD_n_pixels;           // # pixels in MPPC

} ;

#endif 
