#ifndef REALISTICCALODIGI_H
#define REALISTICCALODIGI_H 1

#include "marlin/Processor.h"
#include <IMPL/LCFlagImpl.h>
#include <EVENT/SimCalorimeterHit.h>
#include "lcio.h"
#include <string>
#include <vector>
#include <functional>
#include <optional>

using namespace lcio ;
using namespace marlin ;

/** === RealisticCaloDigi Processor === <br>
    Digitisation of calorimeter hits
    e.g. timing, dead cells, miscalibrations
    this is virtual class, technology-blind
    technology-specific classes can inherit from this one
    D. Jeans 02/2016, rewrite of parts of ILDCaloDigi, DDCaloDigi
    R. Ete 11/2020, rewrite of charge integration and extension of timing treatment
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
  // integration result types
  using integr_res = std::pair<float,float>;
  using integr_res_opt = std::optional<integr_res>;
  using integr_function = std::function<integr_res_opt(const EVENT::SimCalorimeterHit*)>;

  virtual float EnergyDigi(float energy, int id0, int id1);
  virtual integr_res_opt Integrate( const SimCalorimeterHit * hit ) const;
  
  integr_res_opt StandardIntegration( const SimCalorimeterHit * hit ) const ;
  integr_res_opt ROCIntegration( const SimCalorimeterHit * hit ) const ;
  float SmearTime(float time) const;

  // virtual methods to be be overloaded in tech-specific derived classes
  virtual int   getMyUnit()=0;
  virtual float digitiseDetectorEnergy(float energy) = 0 ;
  virtual float convertEnergy( float energy, int inScale ) = 0; // convert energy from input to output scale

  // general parameters

  std::vector<std::string> _inputCollections{};
  std::vector<std::string> _outputCollections{};
  std::vector<std::string> _outputRelCollections{};


  // parameters for digitization effects

  std::string _integration_method{}; // timing calculation method
  float _threshold_value{};         // hit energy threshold
  std::string _threshold_unit{};    // hit energy threshold unit

  int   _time_apply{};              // apply timing cuts?
  int   _time_correctForPropagation{}; // correct times for propagation?
  float _time_windowMin{};          // defn of timing window
  float _time_windowMax{};
  float _fast_shaper{};             // fast shaper value. unit in ns
  float _slow_shaper{};             // slow shaper value. unit in ns
  float _time_resol{};              // time resolution (unit ns)

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
  integr_function _integr_function{};

} ;

#endif



