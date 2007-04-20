#ifndef MaterialDB_h
#define MaterialDB_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <marlin/Global.h>


using namespace lcio ;
using namespace marlin ;


/**======= MaterialDB ============ <br>
 * Processor builds material database used by DELPHI fitting package <br>
 * Information about material shapes and properties are read in from GEAR steering and
 * passed to the DELPHI fitter via global C structures. Processor defines also 
 * extrapolation surfaces on which Track parameters can be evaluated <br>
 * <h4>Input</h4>
 * Processor doesn't require any LCIO collection but needs GEAR steering
 * to create material database <br> 
 * <h4>Output</h4> 
 * No LCIO collections are produced by the processor. 
 * It is meant only as material database builder. <br>
 * @param UseMaterials When this flag is set to 1 material database is built 
 * otherwise no materials are assumed to be present in detector and DELPHI fit
 * ignores effects of particle interactions with detector materials <br> 
 * (default value 1) <br>
 * @param UseExtrapolations When this flag is set to 1 extrapolation surfaces are defined.
 * Otherwise track parameters can be calculated only at the measurement point and no
 * track extrapolations are possible <br>
 * (default value 1) <br>
 * <br>
 * @author A. Raspereza (MPI Munich)
 * @version $Id: MaterialDB.h,v 1.2 2007-04-20 13:44:39 rasp Exp $ 
 */

class MaterialDB : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new MaterialDB ; }
  
  
  MaterialDB() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection name.
   */
  std::string _colName ;

  int _nRun ;
  int _nEvt ;

  std::vector<float> _ladder_dim_z;
  std::vector<float> _layerRadius;
  std::vector<float> _ladderGaps;
  std::vector<float> _zFTD;
  std::vector<float> _stripLine_final_z;
  std::vector<float> _rInFTD;
  std::vector<float> _rOutFTD;
  std::vector<float> _rSIT;
  std::vector<float> _halfZSIT;

  float _support_thickness, _si_thickness;
  float _radlen_ber, _radlen_si, _radlen_si872,_radlen_kapton;
  float _dedx_si,_dedx_ber,_dedx_si872,_dedx_kapton;
  float _dedx_copper;
  float _radlen_copper;
  float _electronicEnd_thickness;
  float _electronicEnd_length;
  float _shell_thickness;
  float _stripLine_thickness;
  float _stripLine_beampipe_r;
  float _FTDdisk_thickness;
  float _FTD_outerSupport_dR;
  float _FTD_innerSupport_dR;
  float _FTD_outerSupport_thickness;
  float _FTD_innerSupport_thickness;
  float _zFTDOuterCyllinderStart;
  float _zFTDOuterCyllinderEnd;
  float _zFTDInnerConeStart;
  float _zFTDInnerConeEnd;
  float _FTD_copper_thickness;
  float _FTD_kaptonCyl_thickness;
  
  float _VTXShell_Radius,_VTXShell_thickness,_VTXShell_HalfZ;
  float _VTXEndPlate_innerRadius;

  float _bField;

  float _beamPipeRadius,_beamPipeHalfZ,_beamPipe_thickness;
  float _SITLayer_thickness;

  int _useExtrapolations;
  int _useMaterials;

} ;

#endif



