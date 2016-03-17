#ifndef REALISTICCALORECOSILICON_H
#define REALISTICCALORECOSILICON_H 1

#include "RealisticCaloReco.h"

#include "DDRec/DetectorData.h"

/** === RealisticCaloRecoSilicon Processor === <br>
    realistic reconstruction of silicon calorimeter hits
    D.Jeans 02/2016.
*/

class RealisticCaloRecoSilicon : public RealisticCaloReco {

 public:
  virtual Processor*  newProcessor() { return new RealisticCaloRecoSilicon ; }
  RealisticCaloRecoSilicon();

 protected:


  virtual void init();

  virtual void getGeometryInformation();

  virtual float reconstructEnergy(const CalorimeterHit* hit);
  virtual void resetGaps();
  virtual void prepareForGaps(const CalorimeterHit* hit);
  virtual void fillGaps();
  virtual bool cellAtWaferEdge(const CalorimeterHit* hit);
  virtual bool cellAtWaferCorner(const CalorimeterHit* hit);

  //virtual const CalorimeterHit* getEdgeGapNeighbour( const CalorimeterHit* mainHit, int stave, int layer );
  //virtual std::vector < const CalorimeterHit* > getCornerGapNeighbours( const CalorimeterHit* mainHit, int stave, int layer );
  virtual std::vector < const CalorimeterHit* > getGapNeighbours( const CalorimeterHit* mainHit, int stave, int layer );

  enum { MAXMOD=7, MAXSTA=8, MAXTOW=10, MAXLAY=33, MAXWAF=60 };
  std::vector < const CalorimeterHit* > edgeHitsModStaTowWaf[MAXMOD][MAXSTA][MAXTOW][MAXLAY][MAXWAF];

  DD4hep::DDRec::LayeredCalorimeterData* _caloData;

  float _cellSize;

} ;

#endif 
