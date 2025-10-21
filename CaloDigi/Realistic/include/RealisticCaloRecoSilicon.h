#ifndef REALISTICCALORECOSILICON_H
#define REALISTICCALORECOSILICON_H 1

#include "RealisticCaloReco.h"

/** === RealisticCaloRecoSilicon Processor === <br>
    realistic reconstruction of silicon calorimeter hits
    D.Jeans 02/2016.

    24 March 2016: removed gap corrections - to be put into separate processor

*/

class RealisticCaloRecoSilicon : public RealisticCaloReco {

public:
  virtual Processor* newProcessor() { return new RealisticCaloRecoSilicon; }
  RealisticCaloRecoSilicon();

protected:
  virtual void init();
  virtual float reconstructEnergy(const CalorimeterHit* hit);
};

#endif
