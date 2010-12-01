#ifndef CLICCDRMaterialDB_h
#define CLICCDRMaterialDB_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <iostream>
#include <string>
#include <marlin/Global.h>


using namespace lcio ;
using namespace marlin ;


/**======= CLICCDRMaterialDB ============ <br>
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
 * @param BuildSET When this flag is set to 1 the SET detector is built. BuildSET must be
 * set to 0 for the Mokka models containing no SET device. <br> 
 * (default value 1) <br>
 * <br>
 * @author A. Raspereza (MPI Munich)
 * @version $Id: CLICCDRMaterialDB.h,v 1.7 2008-06-05 13:40:16 rasp Exp $ 
 */

class CLICCDRMaterialDB : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CLICCDRMaterialDB ; }
  
  
  CLICCDRMaterialDB() ;
  
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

  std::stringstream _errorMsg;

  std::vector<float> _zFTD;
  std::vector<float> _rInFTD;
  std::vector<float> _rOutFTD;
  std::vector<float> _dzSiFTD;
  std::vector<float> _dzSupportFTD;

  std::vector<float> _rSIT;
  std::vector<float> _halfZSIT;
  std::vector<float> _rSITSupport;
  std::vector<float> _halfZSITSupport;

  float _ladder_phi0, _ladder_distance, _ladder_thickness, _ladder_width, _ladder_length;
  float _ladder_offset, _ladder_radLength;
  float _sensitive_distance, _sensitive_thickness, _sensitive_width, _sensitive_length;
  float _sensitive_offset, _sensitive_radLength;

  float _stripLine_final_z;
  float _halfLadderGaps;

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
  
  float _VTXShell_Radius,_VTXShell_thickness,_VTXShell_HalfZ,_VTXShell_radLength;
  float _VTXEndPlate_innerRadius;

  float _bField;

  float _beamPipeRadius,_beamPipeHalfZ,_beamPipe_thickness,_beamPipe_radLength,_beamPipe_dedx;
  float _SITLayer_thickness, _SITLayerSupport_thickness;

  int _useExtrapolations;
  int _useMaterials;

  int _buildSET;

} ;

#endif



