#ifndef ComputeShowerShapesProcessor_hh
#define ComputeShowerShapesProcessor_hh 1

#include <marlin/Processor.h>
#include "EVENT/LCCollection.h"

#include "ClusterShapesMR.h"

using namespace lcio ;
using namespace marlin ;

class ComputeShowerShapesProcessor : public Processor{
public:
  virtual Processor*  newProcessor() { return new ComputeShowerShapesProcessor ; }
  ComputeShowerShapesProcessor();
  virtual void init( LCEvent * evt ) ;
  virtual void processRunHeader( LCRunHeader* run);
  virtual void processEvent( LCEvent * evt );
  virtual void check( LCEvent * evt );
  virtual void end();
 
private:
  marlinreco::ClusterShapes *pClusterShapes;
  std::string _PfoCollection;
  std::string _ClusterCollection;
  float _X01,_X02;
  float _Rm1,_Rm2;
  LCCollection* _PFOCol;
};

#endif
