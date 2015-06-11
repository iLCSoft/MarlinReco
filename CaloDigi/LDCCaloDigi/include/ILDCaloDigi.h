#ifndef DIGITIZER_ILDCCALODIGI_H
#define DIGITIZER_ILDCCALODIGI_H 1

#include "marlin/Processor.h"
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCFlagImpl.h>
#include "CalorimeterHitType.h"
#include "lcio.h"
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "ScintillatorPpdDigi.h"
#include "CLHEP/Random/MTwistEngine.h"

using namespace lcio ;
using namespace marlin ;

const int MAX_LAYERS = 200;
const int MAX_STAVES =  16;

/** === ILDCaloDigi Processor === <br>
 *  Simple calorimeter digitizer Processor. <br>
 *  Takes SimCalorimeterHit Collections and <br>
 *  produces CalorimeterHit Collections. <br>
 *  Simulated energy depositions in active <br>
 *  layers of calorimeters are <br>
 *  converted into physical energy. This is done <br>
 *  taking into account sampling fractions of <br>
 *  ECAL and HCAL. <br>
 *  User has to specify ECAL and HCAL SimCalorimeterHit <br>
 *  collections with processor parameters <br>
 *  HCALCollections and ECALCollections. <br>
 *  The names of the output CalorimeterHit Collections <br>
 *  are specified with processor parameters <br>
 *  ECALOutputCollection and HCALOutputCollection. <br>
 *  Conversion factors for ECAL and HCAL <br>
 *  are specified via processor parameters  <br>
 *  CalibrECAL and CalibrHCAL. <br>
 *  It should be noted that ECAL and HCAL may consist <br>
 *  of several sections with different sampling fractions. <br>
 *  To handle this situation, calibration coefficients for <br>
 *  ECAL and HCAL are passed as arrays of floats with each element <br>
 *  in this array corresponding to certain section with <br>
 *  a given sampling fraction. <br>
 *  List of layer numbers terminating each section are given through <br>
 *  processor parameters ECALLayers and HCALLayers <br>
 *  There is an option to perform digitization of <br> 
 *  both ECAL and HCAL in a digital mode. <br>
 *  Digital mode is activated by  <br>
 *  setting processor parameters <br>
 *  IfDigitalEcal / IfDigitalHcal to 1. <br>
 *  In this case CalibrECAL / CalibrHCAL will  <br>
 *  convert the number of hits into physical energy. <br>
 *  Thresholds on hit energies in ECAL and HCAL <br>
 *  are set with processor parameters <br>
 *  ECALThreshold and HCALThreshold.  <br>
 *  Relations between CalorimeterHits and SimCalorimeterHits <br>
 *  are held in the corresponding relation collection. <br>
 *  The name of this relation collection is specified <br>
 *  via processor parameter RelationOutputCollection. <br> 
 *  <h4>Input collections and prerequisites</h4>
 *  SimCalorimeterHit collections <br>
 *  <h4>Output</h4>
 *  CalorimeterHit collections for ECal and HCal. <br>
 *  Collection of relations <br>
 *  between CalorimeterHits and SimCalorimeterHits. <br> 
 *  For ECal Calorimeter hits the variable type is set to 0, <br>
 *  whereas for HCal Calorimeter hits the type is set to 1 <br>
 *  @author A. Raspereza (DESY) <br>
 *  @author M. Thomson (DESY) <br>
 *  @version $Id$ <br>
 */
class ILDCaloDigi : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ILDCaloDigi ; }
  
  
  ILDCaloDigi() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
   
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;

  virtual void fillECALGaps() ;
  
  float digitalHcalCalibCoeff(CHT::Layout,float energy );

  float analogueHcalCalibCoeff(CHT::Layout, int layer );

  float digitalEcalCalibCoeff(int layer );

  float analogueEcalCalibCoeff(int layer );

 protected:

  float ecalEnergyDigi(float energy, int id0, int id1);
  float ahcalEnergyDigi(float energy, int id0, int id1);

  float siliconDigi(float energy);
  float scintillatorDigi(float energy, bool isEcal);
  LCCollection* combineVirtualStripCells(LCCollection* col, bool isBarrel, int orientation );

  int getNumberOfVirtualCells();
  std::vector < std::pair <int, int> > & getLayerConfig();
  void checkConsistency(std::string colName, int layer);
  std::pair < int, int > getLayerProperties( std::string colName, int layer );
  int getStripOrientationFromColName( std::string colName );


  int _nRun ;
  int _nEvt ;
  
  LCFlagImpl _flag;

  std::vector<std::string> _ecalCollections;
  std::vector<std::string> _hcalCollections;
  std::vector<std::string> _outputEcalCollections;
  std::vector<std::string> _outputHcalCollections;

  std::string _outputRelCollection;

  float _thresholdEcal;
  std::string _unitThresholdEcal;
  std::vector<float> _thresholdHcal;
  std::string _unitThresholdHcal;

  int _digitalEcal;
  int _mapsEcalCorrection;
  int _digitalHcal;

  bool _ECAL_stripHits;

  std::vector<float> _calibrCoeffEcal;
  std::vector<float> _calibrCoeffHcalBarrel;
  std::vector<float> _calibrCoeffHcalEndCap;
  std::vector<float> _calibrCoeffHcalOther;

  std::vector<int> _ecalLayers;
  std::vector<int> _hcalLayers;

  int _ecalGapCorrection;
  float _ecalGapCorrectionFactor;
  float _ecalModuleGapCorrectionFactor;
  float _ecalEndcapCorrectionFactor;
  float _hcalEndcapCorrectionFactor;
  int   _hcalGapCorrection;
  float _hcalModuleGapCorrectionFactor;

  std::vector<CalorimeterHitImpl*> _calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS];
  std::vector<int> _calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS];

  float _zOfEcalEndcap;
  float _barrelPixelSizeT[MAX_LAYERS];
  float _barrelPixelSizeZ[MAX_LAYERS];
  float _endcapPixelSizeX[MAX_LAYERS];
  float _endcapPixelSizeY[MAX_LAYERS];
  float _barrelStaveDir[MAX_STAVES][2];
  
  int   _histograms;

  // timing
  int   _useEcalTiming;
  int   _ecalCorrectTimesForPropagation;
  float _ecalTimeWindowMin;
  float _ecalBarrelTimeWindowMax;
  float _ecalEndcapTimeWindowMax;
  float _ecalDeltaTimeHitResolution;
  float _ecalTimeResolution;
  bool  _ecalSimpleTimingCut;

  int   _useHcalTiming;
  int   _hcalCorrectTimesForPropagation;
  float _hcalTimeWindowMin;
  float _hcalBarrelTimeWindowMax;
  float _hcalEndcapTimeWindowMax;
  float _hcalDeltaTimeHitResolution;
  float _hcalTimeResolution;
  bool  _hcalSimpleTimingCut;
  
  ScintillatorPpdDigi* _scEcalDigi;
  ScintillatorPpdDigi* _scHcalDigi;


  // parameters for extra ECAL digitization effects
  float _calibEcalMip;                // MIP calibration factor
  int   _applyEcalDigi;               // which realistic calib to apply
  float _ecal_PPD_pe_per_mip;         // # photoelectrons/MIP for MPPC
  int   _ecal_PPD_n_pixels;           // # pixels in MPPC
  float _ehEnergy;                    // energy to create e-h pair in silicon
  float _ecal_misCalibNpix;           // miscalibration of # MPPC pixels

  float _misCalibEcal_uncorrel;       // general ECAL miscalibration (uncorrelated between channels)
  bool  _misCalibEcal_uncorrel_keep;  // if true, use the same ECAL cell miscalibs in each event (requires more memory)
  float _misCalibEcal_correl;         // general ECAL miscalibration (100% uncorrelated between channels)

  float _deadCellFractionEcal;        // fraction of random dead channels
  bool  _deadCellEcal_keep;           // keep same cells dead between events?

  float _strip_abs_length;            // absorption length along strip for non-uniformity modeling
  float _ecal_pixSpread;              // relative spread of MPPC pixel signal
  float _ecal_elec_noise;             // electronics noise (as fraction of MIP)
  float _ecalMaxDynMip;               // electronics dynamic range (in terms of MIPs)
  int _ecalStrip_default_nVirt;       // # virtual cells used in Mokka simulation of strips (if available, this is taken from gear file)
  std::string _ecal_deafult_layer_config; // ECAL layer configuration (if available, this is taken from gear file)

  // parameters for extra AHCAL digitization effects
  float _calibHcalMip;                // MIP calibration factor
  int   _applyHcalDigi;               // which realistic calib to apply
  float _hcal_PPD_pe_per_mip;         // # photoelectrons/MIP for MPPC
  int   _hcal_PPD_n_pixels;           // # pixels in MPPC
  float _hcal_misCalibNpix;           // miscalibration of # MPPC pixels

  float _misCalibHcal_uncorrel;       // general ECAL miscalibration (uncorrelated between channels)
  bool  _misCalibHcal_uncorrel_keep;  // if true, use the same AHCAL cell miscalibs in each event (requires more memory)
  float _misCalibHcal_correl;         // general ECAL miscalibration (100% uncorrelated between channels) 

  float _deadCellFractionHcal;        // fraction of random dead channels
  bool  _deadCellHcal_keep;           // keep same cells dead between events?
  float _hcal_pixSpread;              // relative spread of MPPC pixel signal
  float _hcal_elec_noise;             // electronics noise (as fraction of MIP)
  float _hcalMaxDynMip;               // electronics dynamic range (in terms of MIPs)



  // internal variables
  std::vector < std::pair <int, int> > _layerTypes;
  int   _strip_virt_cells;
  int _countWarnings;
  std::string _ecalLayout;

  float _event_correl_miscalib_ecal;
  float _event_correl_miscalib_hcal;
  
  CLHEP::MTwistEngine *_randomEngineDeadCellEcal;
  CLHEP::MTwistEngine *_randomEngineDeadCellHcal;

  std::map < std::pair <int, int> , float > _ECAL_cell_miscalibs;
  std::map < std::pair <int, int> , bool > _ECAL_cell_dead;
  std::map < std::pair <int, int> , float > _HCAL_cell_miscalibs;
  std::map < std::pair <int, int> , bool > _HCAL_cell_dead;

  enum {
    SQUARE,
    STRIP_ALIGN_ALONG_SLAB,
    STRIP_ALIGN_ACROSS_SLAB,
    SIECAL=0,
    SCECAL
  };

  std::string _cellIDLayerString ;
  std::string _cellIDModuleString ;
  std::string _cellIDStaveString ;
  std::string _cellIDIndexIString ;
  std::string _cellIDIndexJString ;
    
  TH1F* fEcal;
  TH1F* fHcal;
  TH1F* fEcalC;
  TH1F* fHcalC;
  TH1F* fEcalC1;
  TH1F* fHcalC1;
  TH1F* fEcalC2;
  TH1F* fHcalC2;
  TH2F* fHcalCvsE; 
  TH2F* fHcalLayer1;
  TH2F* fHcalLayer11;
  TH2F* fHcalLayer21;
  TH2F* fHcalLayer31;
  TH2F* fHcalLayer41;
  TH2F* fHcalLayer51;
  TH2F* fHcalLayer61;
  TH2F* fHcalLayer71;
  TH1F* fHcalRLayer1;
  TH1F* fHcalRLayer11;
  TH1F* fHcalRLayer21;
  TH1F* fHcalRLayer31;
  TH1F* fHcalRLayer41;
  TH1F* fHcalRLayer51;
  TH1F* fHcalRLayer61;
  TH1F* fHcalRLayer71;
  TH1F* fHcalRLayerNorm;

  TH1F* fEcalRLayerNorm;
  TH2F* fEcalLayer1;
  TH2F* fEcalLayer11;
  TH2F* fEcalLayer21;
  TH1F* fEcalRLayer1;
  TH1F* fEcalRLayer11;
  TH1F* fEcalRLayer21;

} ;

#endif



