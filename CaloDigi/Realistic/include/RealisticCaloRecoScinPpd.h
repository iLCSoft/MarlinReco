#ifndef REALISTICCALORECOSCINPPD_H
#define REALISTICCALORECOSCINPPD_H 1

#include "RealisticCaloReco.h"

/**
\addtogroup CaloDigi CaloDigi
@{

\addtogroup RealisticCaloRecoScinPpd RealisticCaloRecoScinPpd
@{
Realistic reconstruction of scint+PPD calorimeter hits.
=== RealisticCaloRecoSilicon Processor === <br>
    realistic reconstruction of scint+PPD calorimeter hits
    D.Jeans 02/2016.
*/

class RealisticCaloRecoScinPpd : public RealisticCaloReco {

 public:

  virtual Processor*  newProcessor() { return new RealisticCaloRecoScinPpd ; }

  RealisticCaloRecoScinPpd();

 protected:
  virtual float reconstructEnergy(const CalorimeterHit* hit);

  float _PPD_pe_per_mip{};         // # photoelectrons/MIP for MPPC
  int   _PPD_n_pixels{};           // # pixels in MPPC
} ;

/** @} @}*/

#endif 
