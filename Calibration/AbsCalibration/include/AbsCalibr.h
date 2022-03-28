#ifndef AbsCalibr_h
#define AbsCalibr_h 1

#include <vector>

#include "marlin/Processor.h"
#include "lcio.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;

namespace marlin {


/**            === AbsCalibr ==== <br>
 *  Processor makes:<br>
 *
 *      Output file for Absolute Energy Calibration <br>
 *      Create collection of energies of calorimeters <br>
 *
 *    @author V. L. Morgunov, A Zhelezov (DESY/ITEP)<br>
 *
 *
 */

  class AbsCalibr : public Processor {
  
  public:
  
    virtual Processor*  newProcessor() { return new AbsCalibr ; }
    AbsCalibr() ;
    //-----------------------------------------------------------------------
    virtual void init() ;
    virtual void processRunHeader( LCRunHeader* run ) ;
    virtual void processEvent( LCEvent * evt ) ; 
    virtual void check( LCEvent * evt ) ; 
    virtual void end() ;
    //-----------------------------------------------------------------------

  protected:

    int _nRun{};
    int _nEvt{};

    enum {
      ECAL1=0,
      ECAL2,
      HCAL
    };

    vector<int> _nlayer{};
    vector<float> _coeff{};
    vector<float> _cuts{};
  } ;
} //namespace marlin
#endif
