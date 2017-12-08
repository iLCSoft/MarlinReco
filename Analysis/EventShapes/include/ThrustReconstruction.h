#ifndef ThrustReconstruction_h
#define ThrustReconstruction_h 1
#include <vector>
#include "marlin/Processor.h"
#include "lcio.h"
#include <iostream>
#include <string>
#include <IMPL/ReconstructedParticleImpl.h>

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Random/RanluxEngine.h>

namespace CLHEP{}    // declare namespace CLHEP for backward compatibility
using namespace CLHEP ;


using namespace lcio ;
using namespace marlin ;



/** Thrust processor for marlin. 
  * Calculates the thrust axis and thrust value for each event using two 
  * different algorithms: 
  * Tasso algorithm --- calculates only the principle thrust value and axis
  * Jetnet algorithm --- calculates the principle thrust value and axis
  *                                 the major thrust value and axis
  *                                 the minor thrust value and axis
 */
class ThrustReconstruction : public Processor {
  
public:
  
  virtual Processor* newProcessor() { return new ThrustReconstruction;}
  
  
  ThrustReconstruction() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void modifyRunHeader( LCRunHeader* run ) {}
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void modifyEvent( LCEvent * evt ) {} 
  
 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
protected:
  int TassoThrust();
  int JetsetThrust();
  double sign(double a,double b);
  double min(double a,double b);

  /** Input collection name.
   */

  std::string  _inputCollectionName{};
  int _typeOfThrustFinder{};

  float _principleThrustValue{};
  float _majorThrustValue{};
  float _minorThrustValue{};
  Hep3Vector _principleThrustAxis{};
  Hep3Vector _majorThrustAxis{};
  Hep3Vector _minorThrustAxis{};
  float _min{},_max{};
  LCCollection* _inParVec{};
  std::vector<Hep3Vector> _partMom{};
  std::string filename{};
  RanluxEngine myrnd{};
} ;

#endif
