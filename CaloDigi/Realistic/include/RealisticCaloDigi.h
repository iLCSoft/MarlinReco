#ifndef REALISTICCALODIGI_H
#define REALISTICCALODIGI_H 1

#include "marlin/Processor.h"
#include <IMPL/LCFlagImpl.h>
#include <EVENT/SimCalorimeterHit.h>
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;

/** === RealisticCaloDigi Processor === <br>
    Digitisation of calorimeter hits
    e.g. timing, dead cells, miscalibrations
    this is virtual class, technology-blind
    technology-specific classes can inherit from this one
    D. Jeans 02/2016, rewrite of parts of ILDCaloDigi, DDCaloDigi
 */

class RealisticCaloDigi : virtual public Processor {
  
 public:

  RealisticCaloDigi( ) ;
  RealisticCaloDigi ( const RealisticCaloDigi& ) = delete;
  RealisticCaloDigi& operator=(const RealisticCaloDigi&) = delete;

  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;

 protected:
  
  // energy scales we know about
  enum { MIP, GEVDEP, NPE };

  virtual float EnergyDigi(float energy, int id0, int id1);
  virtual std::vector < std::pair < float , float > > timingCuts( const SimCalorimeterHit * hit );

  // virtual methods to be be overloaded in tech-specific derived classes
  virtual int   getMyUnit()=0;
  virtual float digitiseDetectorEnergy(float energy) = 0 ;
  virtual float convertEnergy( float energy, int inScale ) = 0; // convert energy from input to output scale

  // general parameters

  std::vector<std::string> _inputCollections{};
  std::vector<std::string> _outputCollections{};
  std::vector<std::string> _outputRelCollections{};


  // parameters for digitization effects


  float _threshold_value{};         // hit energy threshold
  std::string _threshold_unit{};    // hit energy threshold unit

  int   _time_apply{};              // apply timing cuts?
  int   _time_correctForPropagation{}; // correct times for propagation?
  float _time_windowMin{};          // defn of timing window
  float _time_windowMax{};

  float _calib_mip{};               // MIP calibration factor (most probable energy deposit by MIP in active material of one layer)

  float _misCalib_uncorrel{};       // general miscalibration (uncorrelated between channels)
  bool  _misCalib_uncorrel_keep{};  // if true, use the same cell miscalibs over events (requires more memory)

  float _misCalib_correl{};         // general miscalibration (100% uncorrelated between channels)

  float _deadCell_fraction{};       // fraction of random dead channels
  bool  _deadCell_keep{};           // keep same cells dead between events? (requires more memory)

  float _elec_noiseMip{};           // electronics noise (as fraction of MIP)
  float _elec_rangeMip{};           // electronics dynamic range (in terms of MIPs)

  std::string _cellIDLayerString{};
  
  // internal variables

  int _threshold_iunit{};

  LCFlagImpl _flag{};
  LCFlagImpl _flag_rel{};

  float _event_correl_miscalib{};

  std::map < std::pair <int, int> , float > _cell_miscalibs{};
  std::map < std::pair <int, int> , bool  > _cell_dead{};

} ;

#endif



