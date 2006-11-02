#ifndef MaterialDB_h
#define MaterialDB_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <marlin/Global.h>


using namespace lcio ;
using namespace marlin ;


/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id: MaterialDB.h,v 1.1 2006-11-02 12:34:17 rasp Exp $ 
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

  int _useGearFile;

} ;

#endif



