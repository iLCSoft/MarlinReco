#include "MaterialDB.h"
#include <iostream>
#include <string>
#include <stdexcept>


#include <cstdlib>


#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>

#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/GearParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/BField.h>

using namespace lcio ;
using namespace marlin ;

extern "C" {
  void setmat_();
}

extern "C" {
  extern struct {
    int ncmat;
    float rcmat[100];
    float zcmax[100];
    float xrlc[100];
    int npmat;
    float zpmat[50];
    float rpmin[50];
    float rpmax[50];
    float xrlp[50];
    float xelosc[100];
    float xelosp[50];   
    float zcmin[100];
  } fkddes_; 
}

extern "C" {
  extern struct {
    int nconmat;
    float z1conmat[100];
    float z2conmat[100];
    float r1conmat[100];
    float r2conmat[100];
    float xrl1con[100];
    float xrl2con[100];
    float xel1con[100];
    float xel2con[100];
  } fkddes1_; 
}

extern "C" {
  extern struct {
    int nplmat;          //total number of ladders
    float xplmat[200];   //X coordinate of the center of the ladder
    float yplmat[200];   //Y coordinate of the center of the ladder
    float zplmat[200];   //Z coordinate of the center of the ladder
    float widplmat[200]; //width of the ladder
    float lenplmat[200]; //length of the ladder (ladder parallel to z axis)
    float phiplmat[200]; //angle between normal vector of ladder and x axis
    float xrlpl[200];    //radiation length
    float xelospl[200];  //energy loss
  } fkddes2_; 
}

extern "C" {
  extern struct {
    int nexs;
    float rzsurf[50];
    float zrmin[50];
    float zrmax[50];
    int itexts[50];
    int nexhpc;
  } fkexts_;

}

extern "C" {
  struct {
    float consb;
  } fkfild_;
}




MaterialDB aMaterialDB ;


MaterialDB::MaterialDB() : Processor("MaterialDB") {
  
  // Processor description
  _description = "Material DB builder..." ;

 
  // ***************************** //
  //                               //
  // BEAM PIPE GEOMETRY PARAMETERS //
  //                               //
  //   based on model tube00       //
  //                               //
  // ***************************** //
  //   registerProcessorParameter("BeamPipeRadius",
  // 			     "Beam Pipe Radius",
  // 			     _beamPipeRadius,
  // 			     float(10.0));

  //   registerProcessorParameter("BeamPipeHalfZ",
  // 			     "Beam Pipe Half Z",
  // 			     _beamPipeHalfZ,
  // 			     float(61.0));

  //   registerProcessorParameter("BeamPipeThickness",
  // 			     "Beam Pipe Thickness",
  // 			     _beamPipe_thickness,
  // 			     float(0.5));


  // ************************//
  //                         //
  // VXD GEOMETRY PARAMETERS //
  //                         //
  //   based on model vxd00  //
  //                         //
  // ************************//

  //   std::vector<float> layerRadius;
  //   layerRadius.push_back(15.5); 
  //   layerRadius.push_back(27.0);
  //   layerRadius.push_back(38.0);
  //   layerRadius.push_back(49.0);
  //   layerRadius.push_back(60.0);
  //   registerProcessorParameter("LayerRadii",
  // 			     "Radii of VXD layers",
  // 			     _layerRadius,
  // 			     layerRadius);
  
  //   std::vector<float> ladder_dim_z;
  //   ladder_dim_z.push_back(50.0);
  //   ladder_dim_z.push_back(125.0);
  //   ladder_dim_z.push_back(125.0);
  //   ladder_dim_z.push_back(125.0);
  //   ladder_dim_z.push_back(125.0);
  //   registerProcessorParameter("ZLadderDimension",
  //                              "Ladder z dimension",
  //                              _ladder_dim_z,
  //                              ladder_dim_z);

  //   std::vector<float> ladderGaps;
  //   ladderGaps.push_back(0.0);
  //   ladderGaps.push_back(0.04);
  //   ladderGaps.push_back(0.04);
  //   ladderGaps.push_back(0.04);
  //   ladderGaps.push_back(0.04);
  //   registerProcessorParameter("LadderGaps",
  // 			     "Ladder gaps",
  // 			     _ladderGaps,
  // 			     ladderGaps);

  //   registerProcessorParameter("SiThickness",
  //                              "Silicon Thickness",
  //                              _si_thickness,
  //                              float(0.03744));

  //   registerProcessorParameter("SupportThickness",
  //                              "Support Thickness",
  //                              _support_thickness,
  //                              float (0.28224));

  //   registerProcessorParameter("ElectronicEndThickness",
  // 			     "Electronic End Thickness",
  // 			     _electronicEnd_thickness,
  // 			     float(0.19656));

  //   registerProcessorParameter("ElectronicEndLength",
  // 			     "Electronic End Length",
  // 			     _electronicEnd_length,
  // 			     float(10.0)); 

  //   std::vector<float> stripLine_final_z;
  //   stripLine_final_z.push_back(136);
  //   stripLine_final_z.push_back(136);
  //   stripLine_final_z.push_back(140);  
  //   stripLine_final_z.push_back(145);
  //   stripLine_final_z.push_back(150);
  //   registerProcessorParameter("StripLineFinalZ",
  // 			     "Z coordinate for Stripl line End",
  // 			     _stripLine_final_z,
  // 			     stripLine_final_z);

  //   registerProcessorParameter("StripLineThickness",
  // 			     "Z coordinate for Stripl line End",
  // 			     _stripLine_thickness,
  // 			     float(0.09438));

  //   registerProcessorParameter("StripLineBeamPipeR",
  // 			     "Radius for beam pipe ",
  // 			     _stripLine_beampipe_r,
  // 			     float(23.0));
  
  //   registerProcessorParameter("ShellThickness",
  // 			     "Shell Thickness ",
  // 			     _shell_thickness,
  // 			     float(0.49392));
  
  //   _VTXShell_thickness = _shell_thickness;

  //   registerProcessorParameter("VTXShellRadius",
  // 			     "VTX Shell Radius",
  // 			     _VTXShell_Radius,
  // 			     float(65.0));

  //   registerProcessorParameter("VTXShellHalfZ",
  // 			     "VTX Shell HalfZ",
  // 			     _VTXShell_HalfZ,
  // 			     float(135.0));
  
  //   registerProcessorParameter("VTXEndPlateInnerRadius",
  // 			     "VTX Endplate Inner Radius",
  // 			     _VTXEndPlate_innerRadius,
  // 			     float(23.2));

  

  // ************************//
  //                         //
  // FTD GEOMETRY PARAMETERS //
  //                         //
  //   based on model ftd01  //
  //                         //  
  // ************************//

  //   std::vector<float> zFTD;
  //   zFTD.push_back(200);
  //   zFTD.push_back(320);
  //   zFTD.push_back(440);
  //   zFTD.push_back(550);
  //   zFTD.push_back(800);
  //   zFTD.push_back(1050);
  //   zFTD.push_back(1300);
  //   registerProcessorParameter("ZFTDDisk",
  // 			     "Z coordinates of FTD disks",
  // 			     _zFTD,
  // 			     zFTD);

  //   std::vector<float> rInFTD;
  //   rInFTD.push_back(38);
  //   rInFTD.push_back(48);
  //   rInFTD.push_back(59);
  //   rInFTD.push_back(68);
  //   rInFTD.push_back(90);
  //   rInFTD.push_back(111);
  //   rInFTD.push_back(132);
  //   registerProcessorParameter("RinFTDDisk",
  // 			     "Inner Radii of FTD disks",
  // 			     _rInFTD,
  // 			     rInFTD);

  //   std::vector<float> rOutFTD;
  //   rOutFTD.push_back(140);
  //   rOutFTD.push_back(140);
  //   rOutFTD.push_back(210);
  //   rOutFTD.push_back(270);
  //   rOutFTD.push_back(290);
  //   rOutFTD.push_back(290);
  //   rOutFTD.push_back(290);
  //   registerProcessorParameter("RoutFTDDisk",
  // 			     "Outer Radii of FTD disks",
  // 			     _rOutFTD,
  // 			     rOutFTD);

  //   registerProcessorParameter("FTDDiskThickness",
  // 			     "FTD Disk Thickness",
  // 			     _FTDdisk_thickness,
  // 			     float(0.3));
  
  //   registerProcessorParameter("FTDInnerSupportdR",
  // 			     "Inner Support Ring dR",
  // 			     _FTD_innerSupport_dR,
  // 			     float(2));

  //   registerProcessorParameter("FTDOuterSupportdR",
  // 			     "Outer Support Ring dR",
  // 			     _FTD_outerSupport_dR,
  // 			     float(10));

  //   registerProcessorParameter("FTDInnerSupportThickness",
  // 			     "Inner Support Ring Thickness",
  // 			     _FTD_innerSupport_thickness,
  // 			     float(4));

  //   registerProcessorParameter("FTDOuterSupportThickness",
  // 			     "Outer Support Ring Thickness",
  // 			     _FTD_outerSupport_thickness,
  // 			     float(4));

  //   registerProcessorParameter("zFTDOuterCylStart",
  // 			     "z FTD Outer Cyllinder Start",
  // 			     _zFTDOuterCyllinderStart,
  // 			     float(800.));

  //   registerProcessorParameter("zFTDOuterCylEnd",
  // 			     "z FTD Outer Cyllinder End",
  // 			     _zFTDOuterCyllinderEnd,
  // 			     float(1300.));
  
  //   registerProcessorParameter("zFTDInnerConeStart",
  // 			     "z FTD Inner Cone Start",
  // 			     _zFTDInnerConeStart,
  // 			     float(550.));

  //   registerProcessorParameter("zFTDInnerConeEnd",
  // 			     "z FTD Inner Cone End",
  // 			     _zFTDInnerConeEnd,
  // 			     float(1300.));
  
  //   registerProcessorParameter("FTDCopperThickness",
  // 			     "FTD Copper Thickness",
  // 			     _FTD_copper_thickness,
  // 			     float(0.08));
  
  //   registerProcessorParameter("FTDOuterCylThickness",
  // 			     "FTD Outer Cyllinder Thickness",
  // 			     _FTD_kaptonCyl_thickness,
  // 			     float(1.0));

  // SIT Geometry (model sit00)

  //   std::vector<float> rSIT;
  //   rSIT.push_back(160.);
  //   rSIT.push_back(300.);
  //   registerProcessorParameter("SITLayerRadii",
  // 			     "SIT Layer Radii",
  // 			     _rSIT,
  // 			     rSIT);

  //   std::vector<float> halfZSIT;
  //   halfZSIT.push_back(380.);
  //   halfZSIT.push_back(660.);
  //   registerProcessorParameter("SITLayerHalfZ",
  // 			     "SIT Layer HalfZ",
  // 			     _halfZSIT,
  // 			     halfZSIT);

  //   registerProcessorParameter("SITLayerThickness",
  // 			     "SIT Layer Thickness",
  // 			     _SITLayer_thickness,
  // 			     float(0.3));


  // ********************//
  // Material properties //
  //        dE/dx        //
  // ********************//

  //   registerProcessorParameter("dedxBer",
  //                              "Energy loss for Be",
  //                              _dedx_ber,
  //                              float(1.85*1.59e-3));

  //   registerProcessorParameter("dedxSi",
  //                              "Energy loss for Si",
  //                               _dedx_si,
  //                               float(2.33*1.66e-3));
  
  //   registerProcessorParameter("dedxSi872",
  //                              "Energy loss for Si",
  //                               _dedx_si872,
  //                               float(8.72*1.66e-3));
 
  //   registerProcessorParameter("dedxKapton",
  //                              "Energy loss for Kapton",
  //                               _dedx_kapton,
  //                               float(1.42*1.7e-3));

  //   registerProcessorParameter("dedxCopper",
  //                              "Energy loss for Copper",
  //                               _dedx_copper,
  //                               float(8.96*1.6e-3));

  // ********************//
  // Material properties //
  //   Rad length        //
  // ********************//
  
  //   registerProcessorParameter("BeRadiationLength",
  //                              "Be Radiation Length",
  //                              _radlen_ber,
  //                              float(35.28));

  //   registerProcessorParameter("SiRadiationLength",
  //                              "Si Radiation Length",
  //                              _radlen_si,
  //                              float(9.36));

  //   registerProcessorParameter("Si872RadiationLength",
  //                              "Si 8.72 Radiation Length",
  //                              _radlen_si872,
  //                              float(9.36*2.33/8.72));

  //   registerProcessorParameter("KaptonRadiationLength",
  //                              "Kapton Radiation Length",
  //                              _radlen_kapton,
  //                              float(28.6));

  //   registerProcessorParameter("CopperRadiationLength",
  //                              "Copper Radiation Length",
  //                              _radlen_copper,
  //                              float(1.43));
  
  registerProcessorParameter("UseExtrapolations",
			     "Use Extrapolations in Fit",
			     _useExtrapolations,
			     int(1));

  registerProcessorParameter("UseMaterials",
			     "Use material database",
			     _useMaterials,
			     int(1));

}


void MaterialDB::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  int Ncmat=0;
  int Npmat=0;
  int Nconmat=0;
  int Nexs=0;
  int Nplmat=0;

  // **************************************** //
  // ** Building Database for VTX Detector ** //
  // **************************************** //

  //--Get GEAR Parameters--
  const gear::VXDParameters& pVXDDetMain = Global::GEAR->getVXDParameters();
  const gear::VXDLayerLayout& pVXDLayerLayout = pVXDDetMain.getVXDLayerLayout();
  const gear::GearParameters& pVXDDet = Global::GEAR->getGearParameters("VXDInfra");
  const gear::GearParameters& pBeamPipe = Global::GEAR->getGearParameters("BeamPipe");


  //--Berillium beam-pipe--
  _beamPipeRadius = float(pBeamPipe.getDoubleVal("BeamPipeRadius"));
  _beamPipeHalfZ  = float(pBeamPipe.getDoubleVal("BeamPipeHalfZ"));
  _beamPipe_thickness = float(pBeamPipe.getDoubleVal("BeamPipeThickness"));
  _beamPipe_radLength = 0.1*float(pBeamPipe.getDoubleVal("BeamPipeProperties_RadLen"));
  _beamPipe_dedx = 10.0*float(pBeamPipe.getDoubleVal("BeamPipeProperties_dEdx"));

  fkddes_.rcmat[Ncmat] = 0.1*_beamPipeRadius;
  fkddes_.zcmin[Ncmat] = -0.1*_beamPipeHalfZ;
  fkddes_.zcmax[Ncmat] = 0.1*_beamPipeHalfZ;
  fkddes_.xrlc[Ncmat] = 0.1*_beamPipe_thickness / _beamPipe_radLength;
  fkddes_.xelosc[Ncmat] = 0.1*_beamPipe_thickness * _beamPipe_dedx;

  Ncmat++;


  //   fkexts_.itexts[Nexs] = 0;  
  //   fkexts_.rzsurf[Nexs] =  0;
  //   fkexts_.zrmin[Nexs] = -1000.;
  //   fkexts_.zrmax[Nexs] = 1000.;
  //   Nexs++;

  //   fkexts_.itexts[Nexs] = 1;  
  //   fkexts_.rzsurf[Nexs] = 0.0 ;
  //   fkexts_.zrmin[Nexs] = 0.;
  //   fkexts_.zrmax[Nexs] = 1000.;
  //   Nexs++;

  //   int Nb = 10;
  //   float dR = 0.1*_beamPipeRadius;

  //Extension Surface
  fkexts_.itexts[Nexs] = 0;  
  fkexts_.rzsurf[Nexs] = 0.05*_beamPipeRadius ;
  fkexts_.zrmin[Nexs] = -10000.;
  fkexts_.zrmax[Nexs] = 10000.;
  Nexs++;

  //   for (int iB=0;iB<Nb;++iB) {
  //       fkexts_.itexts[Nexs] = 0;  
  //       fkexts_.rzsurf[Nexs] = 0.01+dR*(0.5+float(iB)) ;
  //       fkexts_.zrmin[Nexs] = -10000.;
  //       fkexts_.zrmax[Nexs] = 10000.;
  //       Nexs++;
  //   }
  
  // A.R. FIXME : Introduce additional cones and tubes (Tube00.cc)
  //  Cyllinders and cones in VXD00

  //--Get additional parameters from GEAR--
  int nLayersVTX = pVXDLayerLayout.getNLayers();
  int nLadderGaps = int(pVXDDet.getDoubleVals("LadderGaps").size());
  int nStripLines = int(pVXDDet.getDoubleVals("StripLineFinalZ").size());
  _dedx_si = 10.0*float(pVXDDet.getDoubleVal("ActiveLayerProperties_dEdx"));
  _dedx_ber = 10.0*float(pVXDDet.getDoubleVal("SupportLayerProperties_dEdx"));
  _dedx_kapton = 10.0*float(pVXDDet.getDoubleVal("StripLineProperties_dEdx"));
  _radlen_kapton = 0.1*float(pVXDDet.getDoubleVal("StripLineProperties_RadLen"));
  _electronicEnd_thickness = float(pVXDDet.getDoubleVal("ElectronicEndThickness"));
  _electronicEnd_length = float(pVXDDet.getDoubleVal("ElectronicEndLength"));
  _stripLine_thickness =  float(pVXDDet.getDoubleVal("StripLineThickness"));
  _stripLine_beampipe_r = float(pVXDDet.getDoubleVal("StripLineBeamPipeRadius"));

  // Check values
  if (nLadderGaps != nLayersVTX) {
    _errorMsg << "MaterialDB Processor : vector size of LadderGaps vector ("
	      << nLadderGaps << ")  not equal to number of VXD Layers (" << nLayersVTX << ")" << std::endl;
    throw gear::Exception(_errorMsg.str());
  }
  if (nStripLines != nLayersVTX) {
    _errorMsg << "MaterialDB Processor : vector size of StripLines vector ("
	      << nStripLines << ")  not equal to number of VXD Layers (" << nLayersVTX << ")" << std::endl;
    throw gear::Exception(_errorMsg.str());
  }
  if (pVXDDetMain.getShellOuterRadius()<pVXDDetMain.getShellInnerRadius()) {
    _errorMsg << "Outer Shell Radius (" << pVXDDetMain.getShellOuterRadius() 
	      << ") is smaller than inner Shell radius (" << pVXDDetMain.getShellInnerRadius()
	      << ")" << std::endl;
    throw gear::Exception(_errorMsg.str());
  }

  _shell_thickness = float(pVXDDetMain.getShellOuterRadius()-pVXDDetMain.getShellInnerRadius());
  _VTXShell_thickness = _shell_thickness;
  _VTXShell_Radius = float(pVXDDetMain.getShellInnerRadius());
  _VTXShell_HalfZ = float(pVXDDetMain.getShellHalfLength());
  _VTXEndPlate_innerRadius = float(pVXDDet.getDoubleVal("VXDEndPlateInnerRadius"));
  _VTXShell_radLength = float(pVXDDetMain.getShellRadLength());


  //--The Ladder structure (cylinder or realistic ladder)--
  int nLadders;
  float Pi = acos(-1);
  //  float deg2rad = Pi / 180.0;

  for (int i=0; i<nLayersVTX; ++i) {
    nLadders = pVXDLayerLayout.getNLadders(i);

    _ladder_phi0 = float(pVXDLayerLayout.getPhi0(i));
    _ladder_distance = float(pVXDLayerLayout.getLadderDistance(i));
    _ladder_thickness = float(pVXDLayerLayout.getLadderThickness(i));
    _ladder_width = float(pVXDLayerLayout.getLadderWidth(i));
    _ladder_length = float (pVXDLayerLayout.getLadderLength(i));
    _ladder_offset = float (pVXDLayerLayout.getLadderOffset(i));
    _ladder_radLength = 0.1*float(pVXDLayerLayout.getLadderRadLength(i));

    _sensitive_distance = float(pVXDLayerLayout.getSensitiveDistance(i));
    _sensitive_thickness = float(pVXDLayerLayout.getSensitiveThickness(i));
    _sensitive_width = float(pVXDLayerLayout.getSensitiveWidth(i));
    _sensitive_length = float(pVXDLayerLayout.getSensitiveLength(i));
    _sensitive_offset = float (pVXDLayerLayout.getSensitiveOffset(i));
    _sensitive_radLength = 0.1*float(pVXDLayerLayout.getSensitiveRadLength(i));

    _stripLine_final_z = float(pVXDDet.getDoubleVals("StripLineFinalZ")[i]);
    _halfLadderGaps = float(pVXDDet.getDoubleVals("LadderGaps")[i]);

    if (nLadders < 3) {
      //--Cylinder--

      // Beryllium support    
      fkddes_.rcmat[Ncmat] = 0.1*(_ladder_distance + 0.5*_ladder_thickness);
      fkddes_.zcmin[Ncmat] = -0.1*(_ladder_length + _halfLadderGaps);
      fkddes_.zcmax[Ncmat] = 0.1*(_ladder_length + _halfLadderGaps);
      fkddes_.xrlc[Ncmat] = 0.1*_ladder_thickness/_ladder_radLength;
      fkddes_.xelosc[Ncmat] = 0.1*_ladder_thickness * _dedx_ber;

      //     fkddes_.rcmat[Ncmat] = 0.1*_ladder_distance;
      //     fkddes_.zcmin[Ncmat] = -0.1*_ladder_length;
      //     fkddes_.zcmax[Ncmat] = 0.1*_ladder_length;
      //     fkddes_.xrlc[Ncmat] = 0.1*(_ladder_thickness/_ladder_radLength+_si_thickness/_sensitive_radLength);
      //     fkddes_.xelosc[Ncmat] = 0.1*(_ladder_thickness*_dedx_ber+_sensitive_thickness*_dedx_si);

      // Extension Surface for Kalman filter
      fkexts_.itexts[Nexs] = 0;
      fkexts_.rzsurf[Nexs] = fkddes_.rcmat[Ncmat];
      fkexts_.zrmin[Nexs] = fkddes_.zcmin[Ncmat];
      fkexts_.zrmax[Nexs] = fkddes_.zcmax[Ncmat];

      Ncmat++;
      Nexs++;

      // Electronic Ends;
      // right-hand part
      fkddes_.rcmat[Ncmat] = 0.1*(_ladder_distance + 0.5*_electronicEnd_thickness);
      fkddes_.zcmin[Ncmat] = 0.1*(_ladder_length + _halfLadderGaps);
      fkddes_.zcmax[Ncmat] = 0.1*(_ladder_length + _halfLadderGaps+_electronicEnd_length);
      fkddes_.xrlc[Ncmat] = 0.1*_electronicEnd_thickness / _sensitive_radLength;
      fkddes_.xelosc[Ncmat] = 0.1*_electronicEnd_thickness * _dedx_si;
      
      Ncmat++;

      // left-hand part
      fkddes_.rcmat[Ncmat] = 0.1*(_ladder_distance + 0.5*_electronicEnd_thickness);
      fkddes_.zcmin[Ncmat] = -0.1*(_ladder_length + _halfLadderGaps+_electronicEnd_length);
      fkddes_.zcmax[Ncmat] = -0.1*(_ladder_length + _halfLadderGaps);
      fkddes_.xrlc[Ncmat] = 0.1*_electronicEnd_thickness / _sensitive_radLength;
      fkddes_.xelosc[Ncmat] = 0.1*_electronicEnd_thickness * _dedx_si;

      Ncmat++;

      // Active Silicon
      // right-hand part
      fkddes_.rcmat[Ncmat] = 0.1*(_sensitive_distance + 0.5*_sensitive_thickness);
      fkddes_.zcmin[Ncmat] = 0.1*_halfLadderGaps;
      fkddes_.zcmax[Ncmat] = 0.1*(_sensitive_length + _halfLadderGaps);
      fkddes_.xrlc[Ncmat] = 0.1*_sensitive_thickness / _sensitive_radLength;
      fkddes_.xelosc[Ncmat] = 0.1*_sensitive_thickness * _dedx_si;    

      Ncmat++;

      // left-hand part
      fkddes_.rcmat[Ncmat] = 0.1*(_sensitive_distance + 0.5*_sensitive_thickness);
      fkddes_.zcmin[Ncmat] = -0.1*(_sensitive_length + _halfLadderGaps);
      fkddes_.zcmax[Ncmat] = -0.1*_halfLadderGaps;
      fkddes_.xrlc[Ncmat] = 0.1*_sensitive_thickness / _sensitive_radLength;
      fkddes_.xelosc[Ncmat] = 0.1*_sensitive_thickness * _dedx_si;    

      Ncmat++;

      //     // Strip Lines
      //     // right-hand part    
      //     fkddes1_.z1conmat[Nconmat] = 0.1*(_ladder_length+_halfLadderGaps+
      // 				    _electronicEnd_length+_shell_thickness);
      //     fkddes1_.z2conmat[Nconmat] = 0.1*_stripLine_final_z;
      //     fkddes1_.r1conmat[Nconmat] = 0.1*(_ladder_distance+0.5*_stripLine_thickness);
      //     fkddes1_.r2conmat[Nconmat] = 0.1*(_stripLine_beampipe_r+0.5*_stripLine_thickness);
      //     float dR = fkddes1_.r2conmat[Nconmat]-fkddes1_.r1conmat[Nconmat];
      //     float dZ = fkddes1_.z2conmat[Nconmat]-fkddes1_.z1conmat[Nconmat];
      //     float alpha = atan(dR/dZ);
      //     float cosa  = fabs(cos(alpha));
      //     float dL = _stripLine_thickness*cosa;
      //     fkddes1_.xrl1con[Nconmat]=0.1*dL/_radlen_kapton;
      //     fkddes1_.xrl2con[Nconmat]=0.1*dL/_radlen_kapton;
      //     fkddes1_.xel1con[Nconmat]=0.1*dL*_dedx_kapton;
      //     fkddes1_.xel2con[Nconmat]=0.1*dL*_dedx_kapton;  
      //     Nconmat++;
      
      //     // left-hand part    
      //     fkddes1_.z2conmat[Nconmat] = -0.1*(_ladder_length+_halfLadderGaps+
      // 				       _electronicEnd_length+_shell_thickness);
      //     fkddes1_.z1conmat[Nconmat] = -0.1*_stripLine_final_z;
      //     fkddes1_.r2conmat[Nconmat] = 0.1*(_ladder_distance+0.5*_stripLine_thickness);
      //     fkddes1_.r1conmat[Nconmat] = 0.1*(_stripLine_beampipe_r+0.5*_stripLine_thickness);
      //     dR = fkddes1_.r2conmat[Nconmat]-fkddes1_.r1conmat[Nconmat];
      //     dZ = fkddes1_.z2conmat[Nconmat]-fkddes1_.z1conmat[Nconmat];
      //     alpha = atan(dR/dZ);
      //     cosa  = fabs(cos(alpha));
      //     dL = _stripLine_thickness*cosa;
      //     fkddes1_.xrl1con[Nconmat]=0.1*dL/_radlen_kapton;
      //     fkddes1_.xrl2con[Nconmat]=0.1*dL/_radlen_kapton;
      //     fkddes1_.xel1con[Nconmat]=0.1*dL*_dedx_kapton;
      //     fkddes1_.xel2con[Nconmat]=0.1*dL*_dedx_kapton;  
      //     Nconmat++;  

    } else {
      //--Realistic ladder structure--

      float currPhi;
      float angleLadders = 2*Pi / nLadders;
      float cosphi, sinphi;

      _ladder_distance += 0.5* _ladder_thickness;
      _sensitive_distance +=0.5* _sensitive_thickness;

      for (int j=0; j<nLadders; ++j) {

        currPhi = _ladder_phi0 + (angleLadders * j);
        cosphi = cos(currPhi);
        sinphi = sin(currPhi);

        //Beryllium support
        fkddes2_.xplmat[Nplmat] = 0.1*(_ladder_distance*cosphi - _ladder_offset*sinphi);
        fkddes2_.yplmat[Nplmat] = 0.1*(_ladder_distance*sinphi + _ladder_offset*cosphi);
        fkddes2_.zplmat[Nplmat] = 0.0; //Ladder is centered around z=0

        fkddes2_.widplmat[Nplmat] = 0.1*(_ladder_width); 
        fkddes2_.lenplmat[Nplmat] = 0.1*(2.*(_ladder_length+_halfLadderGaps));
        fkddes2_.phiplmat[Nplmat] = currPhi;

        fkddes2_.xrlpl[Nplmat] = 0.1*_ladder_thickness/_ladder_radLength;
        fkddes2_.xelospl[Nplmat] = 0.1*_ladder_thickness*_dedx_ber;

        Nplmat++;

        //Sensitive Si part
        fkddes2_.xplmat[Nplmat] = 0.1*(_sensitive_distance*cosphi - _sensitive_offset*sinphi);
        fkddes2_.yplmat[Nplmat] = 0.1*(_sensitive_distance*sinphi + _sensitive_offset*cosphi);
        fkddes2_.zplmat[Nplmat] = 0.0; //Ladder is centered around z=0

        fkddes2_.widplmat[Nplmat] = 0.1*(_sensitive_width); 
        fkddes2_.lenplmat[Nplmat] = 0.1*(2.*(_sensitive_length+_halfLadderGaps));
        fkddes2_.phiplmat[Nplmat] = currPhi;

        fkddes2_.xrlpl[Nplmat] = 0.1*_sensitive_thickness/_sensitive_radLength;
        fkddes2_.xelospl[Nplmat] = 0.1*_sensitive_thickness*_dedx_si;

        Nplmat++;
      }

      //Extension Surface (used by the Kalman filter); Radius is set to the distance from IP to the center of the support
      fkexts_.itexts[Nexs] = 0;
      fkexts_.rzsurf[Nexs] = 0.1*_ladder_distance;
      fkexts_.zrmin[Nexs] = - 0.1*_ladder_length;
      fkexts_.zrmax[Nexs] = 0.1*_ladder_length;
      Nexs++;	
    }
  }

  //--Cryostat--
  float AlRadius = float(pVXDDet.getDoubleVal("CryostatAlRadius"));
  float AlHalfLength = float(pVXDDet.getDoubleVal("CryostatAlHalfZ"));
  float AlThickness = float(pVXDDet.getDoubleVal("CryostatAlThickness"));
  float AlZEndCap = float(pVXDDet.getDoubleVal("CryostatAlZEndCap"));
  float AlRinEndCap = float(pVXDDet.getDoubleVal("CryostatAlInnerR"));
  float xrad_cryo = 0.1*float(pVXDDet.getDoubleVal("Cryostat_RadLen"));
  float dedx_cryo = 10.0*float(pVXDDet.getDoubleVal("Cryostat_dEdx"));
  

  // Al cryostat barrel
  fkddes_.rcmat[Ncmat] = 0.1*(AlRadius+0.5*AlThickness);
  fkddes_.zcmin[Ncmat] = -0.1*AlHalfLength;
  fkddes_.zcmax[Ncmat] = 0.1*AlHalfLength;
  fkddes_.xrlc[Ncmat] = 0.1*AlThickness/xrad_cryo;
  fkddes_.xelosc[Ncmat] = 0.1*AlThickness*dedx_cryo;
  Ncmat++;
  

  // Al cryostat endcaps
  fkddes_.zpmat[Npmat] = -0.1*(AlZEndCap+0.5*AlThickness);
  fkddes_.rpmin[Npmat] = 0.1*AlRinEndCap;
  fkddes_.rpmax[Npmat] = 0.1*(AlRadius+AlThickness);;
  fkddes_.xrlp[Npmat] = 0.1*AlThickness/xrad_cryo;
  fkddes_.xelosp[Npmat] = 0.1*AlThickness*dedx_cryo;
  Npmat++;

  fkddes_.zpmat[Npmat] = 0.1*(AlZEndCap+0.5*AlThickness);
  fkddes_.rpmin[Npmat] = 0.1*AlRinEndCap;
  fkddes_.rpmax[Npmat] = 0.1*(AlRadius+AlThickness);;
  fkddes_.xrlp[Npmat] = 0.1*AlThickness/xrad_cryo;
  fkddes_.xelosp[Npmat] = 0.1*AlThickness*dedx_cryo;
  Npmat++;


  //  Outer support cyllinder for VTX
  fkddes_.rcmat[Ncmat] = 0.1*(_VTXShell_Radius+0.5*_VTXShell_thickness);
  fkddes_.zcmin[Ncmat] = -0.1*_VTXShell_HalfZ;
  fkddes_.zcmax[Ncmat] = 0.1*_VTXShell_HalfZ;
  fkddes_.xrlc[Ncmat] = 0.1*_VTXShell_thickness/_VTXShell_radLength;
  fkddes_.xelosc[Ncmat] = 0.1*_VTXShell_thickness*_dedx_ber;
  Ncmat++;
  //  EndPlate support disk for VTX ; left part
  fkddes_.zpmat[Npmat] = -0.1*(_VTXShell_HalfZ+0.5*_VTXShell_thickness);
  fkddes_.rpmin[Npmat] = 0.1*_VTXEndPlate_innerRadius;
  fkddes_.rpmax[Npmat] = 0.1*(_VTXShell_Radius+_VTXShell_thickness);
  fkddes_.xrlp[Npmat] = 0.1*_VTXShell_thickness/_VTXShell_radLength;
  fkddes_.xelosp[Npmat] = 0.1*_VTXShell_thickness*_dedx_ber;  
  Npmat++;    
  //  EndPlate support disk for VTX ; right part
  fkddes_.zpmat[Npmat] = 0.1*(_VTXShell_HalfZ+0.5*_VTXShell_thickness);
  fkddes_.rpmin[Npmat] = 0.1*_VTXEndPlate_innerRadius;
  fkddes_.rpmax[Npmat] = 0.1*(_VTXShell_Radius+_VTXShell_thickness);
  fkddes_.xrlp[Npmat] = 0.1*_VTXShell_thickness/_VTXShell_radLength;
  fkddes_.xelosp[Npmat] = 0.1*_VTXShell_thickness*_dedx_ber;  
  Npmat++;    


  // **************************************** //
  // ** Building Database for FTD Detector ** //
  // **************************************** //
  const gear::GearParameters& pFTDDet = Global::GEAR->getGearParameters("FTD");


  // Check which version of the FTD this is:
  // i)  Support Rings
  // ii) Support Disks
      
  // Planar detectors in FTD

  bool has_si872 (false);
  bool has_SupportDisks (false);

  try{
    int nFTDSupport_dZ = int(pFTDDet.getDoubleVals("FTDDiskSupportThickness").size());
    nFTDSupport_dZ = 0; // avoid compiler warning 
    has_SupportDisks = true;
  }
  catch(gear::UnknownParameterException &e){}

  try{
    float _dedx_si872 = float(pFTDDet.getDoubleVal("Silicon872_dEdx"));
    _dedx_si872 = 0.; // avoid compiler warning 
    has_si872 = true;
  }
  catch(gear::UnknownParameterException &e){}


      
  // Planar detectors in FTD

  if(has_SupportDisks==true && has_si872==false) {
 
    std::cout << "build FTD according to the SupportDisk structure" << std::endl;

    int nLayersFTD = 0;
    int nFTDZ = int(pFTDDet.getDoubleVals("FTDZCoordinate").size());
    int nFTDRin = int(pFTDDet.getDoubleVals("FTDInnerRadius").size());
    int nFTDRout = int(pFTDDet.getDoubleVals("FTDOuterRadius").size());
    int nFTDSi_dZ = int(pFTDDet.getDoubleVals("FTDDiskSiThickness").size());
    int nFTDSupport_dZ = int(pFTDDet.getDoubleVals("FTDDiskSupportThickness").size());
  

    if (nFTDZ == nFTDRin && nFTDRin == nFTDRout && nFTDRout == nFTDSi_dZ && nFTDSi_dZ == nFTDSupport_dZ) {
      nLayersFTD = nFTDZ;
      _zFTD.resize(nLayersFTD);
      _rInFTD.resize(nLayersFTD);
      _rOutFTD.resize(nLayersFTD);
      _dzSiFTD.resize(nLayersFTD);
      _dzSupportFTD.resize(nLayersFTD);
    }
    else {
      _errorMsg << "Size of vectors FTDZCoordinate, FTDInnerRadius, FTDInnerRadius, FTDDiskSiThickness and FTDDiskSupportThickness are not equal --->" 
		<< " # FTDZCoordinate : " << nFTDZ
		<< " # FTDInnerRadius : " << nFTDRin
		<< " # FTDOuterRadius : " << nFTDRout 
		<< " # FTDDiskSiThickness : " << nFTDSi_dZ
		<< " # FTDDiskSupportThickness : " << nFTDSupport_dZ
		<< std::endl ;
      throw gear::Exception(_errorMsg.str());
    }
    for (int i=0;i<nLayersFTD;++i) {
      _zFTD[i] = float(pFTDDet.getDoubleVals("FTDZCoordinate")[i]);
      _rInFTD[i] = float(pFTDDet.getDoubleVals("FTDInnerRadius")[i]);
      _rOutFTD[i] = float(pFTDDet.getDoubleVals("FTDOuterRadius")[i]);
      _dzSiFTD[i] = float(pFTDDet.getDoubleVals("FTDDiskSiThickness")[i]);
      _dzSupportFTD[i] = float(pFTDDet.getDoubleVals("FTDDiskSupportThickness")[i]);
    } 

    _zFTDOuterCyllinderStart = float(pFTDDet.getDoubleVal("zFTDOuterCylinderStart"));
    _zFTDOuterCyllinderEnd = float(pFTDDet.getDoubleVal("zFTDOuterCylinderEnd"));
    _zFTDInnerConeStart = float(pFTDDet.getDoubleVal("zFTDInnerConeStart"));
    _zFTDInnerConeEnd = float(pFTDDet.getDoubleVal("zFTDInnerConeEnd"));
    _FTD_copper_thickness = float(pFTDDet.getDoubleVal("FTDCopperThickness"));
    _FTD_kaptonCyl_thickness = float(pFTDDet.getDoubleVal("FTDOuterCylinderThickness"));

    _dedx_si = 10.0*float(pFTDDet.getDoubleVal("Silicon_dEdx"));
    _dedx_kapton = 10.0*float(pFTDDet.getDoubleVal("Kapton_dEdx"));
    _dedx_copper = 10.0*float(pFTDDet.getDoubleVal("Copper_dEdx"));
    _radlen_si = 0.1*float(pFTDDet.getDoubleVal("Silicon_RadLen"));  
    _radlen_kapton = 0.1*float(pFTDDet.getDoubleVal("Kapton_RadLen"));
    _radlen_copper = 0.1*float(pFTDDet.getDoubleVal("Copper_RadLen"));

    for (int i=0;i<nLayersFTD;++i) {
      // FTD Si Disks
      float dedx = _dedx_si;
      float radlen  = _radlen_si;

      // right-hand part
      fkddes_.zpmat[Npmat] = 0.1*_zFTD[i];
      fkddes_.rpmin[Npmat] = 0.1*_rInFTD[i];
      fkddes_.rpmax[Npmat] = 0.1*_rOutFTD[i];
      fkddes_.xrlp[Npmat] = 0.1*_dzSiFTD[i]/radlen;
      fkddes_.xelosp[Npmat] = 0.1*_dzSiFTD[i]*dedx;

      fkexts_.itexts[Nexs] = 1;
      fkexts_.rzsurf[Nexs] = fkddes_.zpmat[Npmat];
      fkexts_.zrmin[Nexs] = fkddes_.rpmin[Npmat];
      fkexts_.zrmax[Nexs] = fkddes_.rpmax[Npmat];

      Nexs++;
      Npmat++;
    
      // left-hand part
      fkddes_.zpmat[Npmat] = -0.1*_zFTD[i];
      fkddes_.rpmin[Npmat] = 0.1*_rInFTD[i];
      fkddes_.rpmax[Npmat] = 0.1*_rOutFTD[i];
      fkddes_.xrlp[Npmat] = 0.1*_dzSiFTD[i]/radlen;
      fkddes_.xelosp[Npmat] = 0.1*_dzSiFTD[i]*dedx;

      fkexts_.itexts[Nexs] = 1;
      fkexts_.rzsurf[Nexs] = fkddes_.zpmat[Npmat];
      fkexts_.zrmin[Nexs] = fkddes_.rpmin[Npmat];
      fkexts_.zrmax[Nexs] = fkddes_.rpmax[Npmat];

      Nexs++;
      Npmat++;
      
      // Support Disks
      // right-hand part
      fkddes_.zpmat[Npmat] = 0.1*(_zFTD[i]+(_dzSiFTD[i]/2.0)+(_dzSupportFTD[i]/2.0));
      fkddes_.rpmin[Npmat] = 0.1*_rInFTD[i];
      fkddes_.rpmax[Npmat] = 0.1*_rOutFTD[i];
      fkddes_.xrlp[Npmat] = 0.1*(_dzSupportFTD[i]/_radlen_kapton);
      fkddes_.xelosp[Npmat] = 0.1*(_dzSupportFTD[i]*_dedx_kapton);
      Npmat++;
   
      // left-hand part
      fkddes_.zpmat[Npmat] = -0.1*(_zFTD[i]+(_dzSiFTD[i]/2.0)+(_dzSupportFTD[i]/2.0));
      fkddes_.rpmin[Npmat] = 0.1*_rInFTD[i];
      fkddes_.rpmax[Npmat] = 0.1*_rOutFTD[i];
      fkddes_.xrlp[Npmat] = 0.1*_dzSupportFTD[i]/_radlen_kapton;
      fkddes_.xelosp[Npmat] = 0.1*_dzSupportFTD[i]*_dedx_kapton;
      Npmat++;
   
    }

    // Outer support cyllinders (FTD)

    // copper cables; left part
    fkddes_.rcmat[Ncmat] = 0.1*(_rOutFTD[nLayersFTD-1]+
				_FTD_outerSupport_dR+0.5+
				0.5*_FTD_copper_thickness);
    fkddes_.zcmin[Ncmat] = -0.1*_zFTDOuterCyllinderEnd;
    fkddes_.zcmax[Ncmat] = -0.1*_zFTDOuterCyllinderStart;
    fkddes_.xrlc[Ncmat] = 0.1*_FTD_copper_thickness/_radlen_copper;
    fkddes_.xelosc[Ncmat] = 0.1*_FTD_copper_thickness*_dedx_copper;
    Ncmat++;
    // copper cables; right part
    fkddes_.rcmat[Ncmat] = 0.1*(_rOutFTD[nLayersFTD-1]+
				_FTD_outerSupport_dR+0.5+ 
				0.5*_FTD_copper_thickness);
    fkddes_.zcmin[Ncmat] = 0.1*_zFTDOuterCyllinderStart;
    fkddes_.zcmax[Ncmat] = 0.1*_zFTDOuterCyllinderEnd;
    fkddes_.xrlc[Ncmat] = 0.1*_FTD_copper_thickness/_radlen_copper;
    fkddes_.xelosc[Ncmat] = 0.1*_FTD_copper_thickness*_dedx_copper;
    Ncmat++;  
    // kapton cyllinder; left part
    fkddes_.rcmat[Ncmat] = 0.1*(_rOutFTD[nLayersFTD-1]+
				_FTD_outerSupport_dR+0.5+
				_FTD_copper_thickness+
				0.5*_FTD_kaptonCyl_thickness);
    fkddes_.zcmin[Ncmat] = -0.1*_zFTDOuterCyllinderEnd;
    fkddes_.zcmax[Ncmat] = -0.1*_zFTDOuterCyllinderStart;
    fkddes_.xrlc[Ncmat] = 0.1*_FTD_kaptonCyl_thickness/_radlen_kapton;
    fkddes_.xelosc[Ncmat] = 0.1*_FTD_kaptonCyl_thickness*_dedx_kapton;
    Ncmat++;
    // kapton cyllinder; right part
    fkddes_.rcmat[Ncmat] = 0.1*(_rOutFTD[nLayersFTD-1]+
				_FTD_outerSupport_dR+0.5+
				_FTD_copper_thickness+
				0.5*_FTD_kaptonCyl_thickness);
    fkddes_.zcmin[Ncmat] = 0.1*_zFTDOuterCyllinderStart;
    fkddes_.zcmax[Ncmat] = 0.1*_zFTDOuterCyllinderEnd;
    fkddes_.xrlc[Ncmat] = 0.1*_FTD_kaptonCyl_thickness/_radlen_kapton;
    fkddes_.xelosc[Ncmat] = 0.1*_FTD_kaptonCyl_thickness*_dedx_kapton;
    Ncmat++;

  }
  
  else if(has_si872==true && has_SupportDisks==false) {

    std::cout << "build FTD according to the SupportRing structure" << std::endl;

    int nLayersFTD = 0;
    int nFTDZ = int(pFTDDet.getDoubleVals("FTDZCoordinate").size());
    int nFTDRin = int(pFTDDet.getDoubleVals("FTDInnerRadius").size());
    int nFTDRout = int(pFTDDet.getDoubleVals("FTDOuterRadius").size());
  

    if (nFTDZ == nFTDRin && nFTDRin == nFTDRout) {
      nLayersFTD = nFTDZ;
      _zFTD.resize(nLayersFTD);
      _rInFTD.resize(nLayersFTD);
      _rOutFTD.resize(nLayersFTD);
    }
    else {
      std::cout << "Size of vectors FTDZCoordinate, FTDInnerRadius and  FTDInnerRadius are not equal --->" << std::endl;
      std::cout << "# FTDZCoordinate : " << nFTDZ << std::endl;
      std::cout << "# FTDInnerRadius : " << nFTDRin << std::endl;
      std::cout << "# FTDOuterRadius : " << nFTDRout << std::endl;
      exit(1);
    }
    for (int i=0;i<nLayersFTD;++i) {
      _zFTD[i] = float(pFTDDet.getDoubleVals("FTDZCoordinate")[i]);
      _rInFTD[i] = float(pFTDDet.getDoubleVals("FTDInnerRadius")[i]);
      _rOutFTD[i] = float(pFTDDet.getDoubleVals("FTDOuterRadius")[i]);
    } 
    _FTDdisk_thickness = float(pFTDDet.getDoubleVal("FTDDiskThickness"));
    _FTD_innerSupport_dR = float(pFTDDet.getDoubleVal("FTDInnerSupportdR"));
    _FTD_outerSupport_dR = float(pFTDDet.getDoubleVal("FTDOuterSupportdR"));
    _FTD_innerSupport_thickness = float(pFTDDet.getDoubleVal("FTDInnerSupportThickness"));
    _FTD_outerSupport_thickness = float(pFTDDet.getDoubleVal("FTDOuterSupportThickness"));
    _zFTDOuterCyllinderStart = float(pFTDDet.getDoubleVal("zFTDOuterCylinderStart"));
    _zFTDOuterCyllinderEnd = float(pFTDDet.getDoubleVal("zFTDOuterCylinderEnd"));
    _zFTDInnerConeStart = float(pFTDDet.getDoubleVal("zFTDInnerConeStart"));
    _zFTDInnerConeEnd = float(pFTDDet.getDoubleVal("zFTDInnerConeEnd"));
    _FTD_copper_thickness = float(pFTDDet.getDoubleVal("FTDCopperThickness"));
    _FTD_kaptonCyl_thickness = float(pFTDDet.getDoubleVal("FTDOuterCylinderThickness"));
    int iLast = pFTDDet.getIntVal("LastHeavyLayer");
    _dedx_si = 10.0*float(pFTDDet.getDoubleVal("Silicon_dEdx"));
    if (iLast>0) {
      _dedx_si872 = 10.0*float(pFTDDet.getDoubleVal("Silicon872_dEdx"));
      _radlen_si872 = 0.1*float(pFTDDet.getDoubleVal("Silicon872_RadLen"));
    }
    _dedx_kapton = 10.0*float(pFTDDet.getDoubleVal("Kapton_dEdx"));
    _dedx_copper = 10.0*float(pFTDDet.getDoubleVal("Copper_dEdx"));
    _radlen_si = 0.1*float(pFTDDet.getDoubleVal("Silicon_RadLen"));  
    _radlen_kapton = 0.1*float(pFTDDet.getDoubleVal("Kapton_RadLen"));
    _radlen_copper = 0.1*float(pFTDDet.getDoubleVal("Copper_RadLen"));

    for (int i=0;i<nLayersFTD;++i) {
      // FTD Si Disks
      float dedx = _dedx_si;
      float radlen  = _radlen_si;
      if (i<iLast) {
	dedx = _dedx_si872;
	radlen = _radlen_si872;
      }
      // right-hand part
      fkddes_.zpmat[Npmat] = 0.1*_zFTD[i];
      fkddes_.rpmin[Npmat] = 0.1*_rInFTD[i];
      fkddes_.rpmax[Npmat] = 0.1*_rOutFTD[i];
      fkddes_.xrlp[Npmat] = 0.1*_FTDdisk_thickness/radlen;
      fkddes_.xelosp[Npmat] = 0.1*_FTDdisk_thickness*dedx;

      fkexts_.itexts[Nexs] = 1;
      fkexts_.rzsurf[Nexs] = fkddes_.zpmat[Npmat];
      fkexts_.zrmin[Nexs] = fkddes_.rpmin[Npmat];
      fkexts_.zrmax[Nexs] = fkddes_.rpmax[Npmat];

      Nexs++;
      Npmat++;
    
      // left-hand part
      fkddes_.zpmat[Npmat] = -0.1*_zFTD[i];
      fkddes_.rpmin[Npmat] = 0.1*_rInFTD[i];
      fkddes_.rpmax[Npmat] = 0.1*_rOutFTD[i];
      fkddes_.xrlp[Npmat] = 0.1*_FTDdisk_thickness/radlen;
      fkddes_.xelosp[Npmat] = 0.1*_FTDdisk_thickness*dedx;

      fkexts_.itexts[Nexs] = 1;
      fkexts_.rzsurf[Nexs] = fkddes_.zpmat[Npmat];
      fkexts_.zrmin[Nexs] = fkddes_.rpmin[Npmat];
      fkexts_.zrmax[Nexs] = fkddes_.rpmax[Npmat];

      Nexs++;
      Npmat++;
      
      // Inner Support Rings
      // right-hand part
      fkddes_.zpmat[Npmat] = 0.1*_zFTD[i];
      fkddes_.rpmin[Npmat] = 0.1*(_rInFTD[i]-_FTD_innerSupport_dR);
      fkddes_.rpmax[Npmat] = 0.1*_rInFTD[i];
      fkddes_.xrlp[Npmat] = 0.1*_FTD_innerSupport_thickness/_radlen_kapton;
      fkddes_.xelosp[Npmat] = 0.1*_FTD_innerSupport_thickness*_dedx_kapton;
      Npmat++;
   
      // left-hand part
      fkddes_.zpmat[Npmat] = -0.1*_zFTD[i];
      fkddes_.rpmin[Npmat] = 0.1*(_rInFTD[i]-_FTD_innerSupport_dR);
      fkddes_.rpmax[Npmat] = 0.1*_rInFTD[i];
      fkddes_.xrlp[Npmat] = 0.1*_FTD_innerSupport_thickness/_radlen_kapton;
      fkddes_.xelosp[Npmat] = 0.1*_FTD_innerSupport_thickness*_dedx_kapton;
      Npmat++;
   
      // Outer Support
      // right-hand part
      fkddes_.zpmat[Npmat] = 0.1*_zFTD[i];
      fkddes_.rpmin[Npmat] = 0.1*_rOutFTD[i];
      fkddes_.rpmax[Npmat] = 0.1*(_rOutFTD[i]+_FTD_outerSupport_dR);
      fkddes_.xrlp[Npmat] = 0.1*_FTD_outerSupport_thickness/_radlen_kapton;
      fkddes_.xelosp[Npmat] = 0.1*_FTD_outerSupport_thickness*_dedx_kapton;
      Npmat++;
   
      // left-hand part
      fkddes_.zpmat[Npmat] = -0.1*_zFTD[i];
      fkddes_.rpmin[Npmat] = 0.1*_rOutFTD[i];
      fkddes_.rpmax[Npmat] = 0.1*(_rOutFTD[i]+_FTD_outerSupport_dR);
      fkddes_.xrlp[Npmat]  = 0.1*_FTD_outerSupport_thickness/_radlen_kapton;
      fkddes_.xelosp[Npmat]= 0.1*_FTD_outerSupport_thickness*_dedx_kapton;
      Npmat++;
    }
  
  

    // copper cables; left part
    fkddes_.rcmat[Ncmat] = 0.1*(_rOutFTD[nLayersFTD-1]+
				_FTD_outerSupport_dR+0.5+
				0.5*_FTD_copper_thickness);
    fkddes_.zcmin[Ncmat] = -0.1*_zFTDOuterCyllinderEnd;
    fkddes_.zcmax[Ncmat] = -0.1*_zFTDOuterCyllinderStart;
    fkddes_.xrlc[Ncmat] = 0.1*_FTD_copper_thickness/_radlen_copper;
    fkddes_.xelosc[Ncmat] = 0.1*_FTD_copper_thickness*_dedx_copper;
    Ncmat++;
  
    // copper cables; right part
    fkddes_.rcmat[Ncmat] = 0.1*(_rOutFTD[nLayersFTD-1]+
				_FTD_outerSupport_dR+0.5+ 
				0.5*_FTD_copper_thickness);
    fkddes_.zcmin[Ncmat] = 0.1*_zFTDOuterCyllinderStart;
    fkddes_.zcmax[Ncmat] = 0.1*_zFTDOuterCyllinderEnd;
    fkddes_.xrlc[Ncmat] = 0.1*_FTD_copper_thickness/_radlen_copper;
    fkddes_.xelosc[Ncmat] = 0.1*_FTD_copper_thickness*_dedx_copper;
    Ncmat++;  
  
    // kapton cyllinder; left part
    fkddes_.rcmat[Ncmat] = 0.1*(_rOutFTD[nLayersFTD-1]+
				_FTD_outerSupport_dR+0.5+
				_FTD_copper_thickness+
				0.5*_FTD_kaptonCyl_thickness);
    fkddes_.zcmin[Ncmat] = -0.1*_zFTDOuterCyllinderEnd;
    fkddes_.zcmax[Ncmat] = -0.1*_zFTDOuterCyllinderStart;
    fkddes_.xrlc[Ncmat] = 0.1*_FTD_kaptonCyl_thickness/_radlen_kapton;
    fkddes_.xelosc[Ncmat] = 0.1*_FTD_kaptonCyl_thickness*_dedx_kapton;
    Ncmat++;
  
    // kapton cyllinder; right part
    fkddes_.rcmat[Ncmat] = 0.1*(_rOutFTD[nLayersFTD-1]+
				_FTD_outerSupport_dR+0.5+
				_FTD_copper_thickness+
				0.5*_FTD_kaptonCyl_thickness);
    fkddes_.zcmin[Ncmat] = 0.1*_zFTDOuterCyllinderStart;
    fkddes_.zcmax[Ncmat] = 0.1*_zFTDOuterCyllinderEnd;
    fkddes_.xrlc[Ncmat] = 0.1*_FTD_kaptonCyl_thickness/_radlen_kapton;
    fkddes_.xelosc[Ncmat] = 0.1*_FTD_kaptonCyl_thickness*_dedx_kapton;
    Ncmat++;
  }
  else {
    _errorMsg << "MaterialDB Processor : FTD Geometery not correctly described. \n"
	      << " It is neither SupportRing or SupportDisk based."
	      <<  std::endl;
    throw gear::Exception(_errorMsg.str());
  }

  // **************************************** //
  // ** Building Database for SIT Detector ** //
  // **************************************** //
  const gear::GearParameters& pSITDet = Global::GEAR->getGearParameters("SIT");

  // SIT layers

  int nSITR = int(pSITDet.getDoubleVals("SITLayerRadius").size());
  int nSITHL = int(pSITDet.getDoubleVals("SITLayerHalfLength").size());
  int SITModel = int(pSITDet.getIntVal("SITModel"));
  int nLayersSIT = 0;

  if (nSITR == nSITHL) {
    nLayersSIT = nSITR;
    _rSIT.resize(nLayersSIT);
    _halfZSIT.resize(nLayersSIT);
    if (SITModel>0) {
      _rSITSupport.resize(nLayersSIT);
      _halfZSITSupport.resize(nLayersSIT);
    }
  }
  else {
    _errorMsg << "Size of SITLayerRadius vector (" << nSITR 
	      << ") is not equal to the size of SITHalfLength vector ("
	      << nSITHL << ")" << std::endl;
    throw gear::Exception(_errorMsg.str());
  }

  for (int iL=0;iL<nLayersSIT;++iL) {
    _rSIT[iL] = float(pSITDet.getDoubleVals("SITLayerRadius")[iL]);
    _halfZSIT[iL] = float(pSITDet.getDoubleVals("SITLayerHalfLength")[iL]);
    if (SITModel>0) {
      _rSITSupport[iL] = float(pSITDet.getDoubleVals("SITSupportLayerRadius")[iL]);
      _halfZSITSupport[iL] = float(pSITDet.getDoubleVals("SITSupportLayerHalfLength")[iL]);
    }
  }

  _radlen_si = 0.1*float(pSITDet.getDoubleVal("SITLayer_RadLen"));
  _dedx_si = 10.*float(pSITDet.getDoubleVal("SITLayer_dEdx"));
  _SITLayer_thickness =  float(pSITDet.getDoubleVal("SITLayerThickness"));

  if (SITModel>0) {
    _radlen_ber = 0.1*float(pSITDet.getDoubleVal("SITSupportLayer_RadLen"));
    _dedx_ber = 10.*float(pSITDet.getDoubleVal("SITSupportLayer_dEdx"));
    _SITLayerSupport_thickness =  float(pSITDet.getDoubleVal("SITSupportLayerThickness"));
  }


  for (int iL = 0; iL < nLayersSIT; ++iL) {
    
    float rmin = _rSIT[iL];
    if (SITModel>0) {
      if (_rSITSupport[iL]<rmin)
	rmin = _rSITSupport[iL];      
    }

    rmin = rmin - 0.01; 

    float halfZ = _halfZSIT[iL];

    fkddes_.rcmat[Ncmat] = 0.1*_rSIT[iL];
    fkddes_.zcmin[Ncmat] = -0.1*_halfZSIT[iL];
    fkddes_.zcmax[Ncmat] = 0.1*_halfZSIT[iL];
    fkddes_.xrlc[Ncmat] = 0.1*_SITLayer_thickness/_radlen_si;
    fkddes_.xelosc[Ncmat] = 0.1*_SITLayer_thickness*_dedx_si;     
    Ncmat++;

    if (SITModel>0) {
      fkddes_.rcmat[Ncmat] = 0.1*_rSITSupport[iL];
      fkddes_.zcmin[Ncmat] = -0.1*_halfZSITSupport[iL];
      fkddes_.zcmax[Ncmat] = 0.1*_halfZSITSupport[iL];
      fkddes_.xrlc[Ncmat] = 0.1*_SITLayerSupport_thickness/_radlen_ber;
      fkddes_.xelosc[Ncmat] = 0.1*_SITLayerSupport_thickness*_dedx_ber;
      Ncmat++;
    }

    fkexts_.itexts[Nexs] = 0;
    fkexts_.rzsurf[Nexs] = rmin;
    fkexts_.zrmin[Nexs]  = -halfZ;
    fkexts_.zrmax[Nexs]  = halfZ;

    Nexs++;

  } 


  // **************************************** //
  // ** Building Database for TPC Detector ** //
  // **************************************** //

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  //  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  //  const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;
  
  //   float innerrad = 0.1 * float( planeExt[0] ) ;
  //   float outerrad = 0.1 *float( planeExt[1] ) ;
  float maxdrift = 0.1 * float( gearTPC.getMaxDriftLength() );
  //  float npadrows = padLayout.getNRows() ;
  //  float tpcpixz = 0.1 * float(gearTPC.getDoubleVal("tpcPixZ")) ;
  //  float ionpoten = 0.1 * float(gearTPC.getDoubleVal("tpcIonPotential")) ;  
  //  float tpcrpres = 0.1 * float(gearTPC.getDoubleVal("tpcRPhiResMax")) ;  
  //  float tpczres = 0.1 * float(gearTPC.getDoubleVal("tpcZRes")) ;
  //  float tpcbfield = float(gearTPC.getDoubleVal("BField")) ;
  //  float xralu = 8.9;
  //  float dedxalu = 2.70*1.62e-3;
  //  float xrargon = 10971.;
  //  float dedxargon = 0.0018*1.52e-3;
  


  float RTPCINN = 0.1*float(gearTPC.getDoubleVal("tpcInnerRadius"));
  float RTPCOUT = 0.1*float(gearTPC.getDoubleVal("tpcOuterRadius"));
  float TPCTHBI = 0.1*float(gearTPC.getDoubleVal("tpcInnerWallThickness"));
  float TPCTHBO = 0.1*float(gearTPC.getDoubleVal("tpcOuterWallThickness"));
  //  float TPCACRI = innerrad;
  //  float TPCACRO = outerrad;
  float TPCHLFZ = maxdrift;
  float xralu = 0.1*float(gearTPC.getDoubleVal("TPCWallProperties_RadLen"));
  float dedxalu = 10.*float(gearTPC.getDoubleVal("TPCWallProperties_dEdx"));
  float xrargon = 0.1*float(gearTPC.getDoubleVal("TPCGasProperties_RadLen"));
  float dedxargon = 10.*float(gearTPC.getDoubleVal("TPCGasProperties_dEdx"));
  

  // inner tube    
  fkddes_.rcmat[Ncmat] = RTPCINN + TPCTHBI/2.;
  fkddes_.zcmin[Ncmat] =  -TPCHLFZ;
  fkddes_.zcmax[Ncmat] = TPCHLFZ;
  fkddes_.xrlc[Ncmat] = TPCTHBI/xralu;
  fkddes_.xelosc[Ncmat] = TPCTHBI*dedxalu;

  fkexts_.itexts[Nexs] = 0;
  fkexts_.rzsurf[Nexs] = fkddes_.rcmat[Ncmat]+0.1;
  fkexts_.zrmin[Nexs] = fkddes_.zcmin[Ncmat];
  fkexts_.zrmax[Nexs] = fkddes_.zcmax[Ncmat];

  Nexs++;
  Ncmat++;

  int ncyl = 50;
 
  // Gas volume in TPC
  float xstep = (RTPCOUT-RTPCINN-TPCTHBO-TPCTHBI)/float(ncyl);
  for (int icyl=0;icyl<ncyl;++icyl) {
    fkddes_.rcmat[Ncmat] = RTPCINN + TPCTHBI + (float(icyl)+0.5)*xstep;
    fkddes_.zcmin[Ncmat] =  -TPCHLFZ;
    fkddes_.zcmax[Ncmat] = TPCHLFZ;
    fkddes_.xrlc[Ncmat] = xstep/xrargon;
    fkddes_.xelosc[Ncmat] = xstep*dedxargon;
    Ncmat++;
  }

  // outer tube
  //  fkddes_.rcmat[Ncmat] = RTPCOUT;
  //  fkddes_.zcmin[Ncmat] =  -TPCHLFZ;
  //  fkddes_.zcmax[Ncmat] = TPCHLFZ;
  //  fkddes_.xrlc[Ncmat] = TPCTHBO/xralu;
  //  fkddes_.xelosc[Ncmat] = TPCTHBO*dedxalu;
  //  Ncmat++;
    
  // endplates
  //  fkddes_.rpmin[Npmat] = RTPCINN;
  //  fkddes_.rpmax[Npmat] = RTPCOUT;
  //  fkddes_.zpmat[Npmat] = -(TPCHLFZ-TPCTHKE/2.);
  //  fkddes_.xrlp[Npmat] = 0.35;
  //  fkddes_.xelosp[Npmat] = 0.35*xralu*dedxalu;
  //  Npmat++;

  //  fkddes_.rpmin[Npmat] = RTPCINN;
  //  fkddes_.rpmax[Npmat] = RTPCOUT;
  //  fkddes_.zpmat[Npmat] = TPCHLFZ-TPCTHKE/2.;
  //  fkddes_.xrlp[Npmat] = 0.35;
  //  fkddes_.xelosp[Npmat] = 0.35*xralu*dedxalu;
  //  Npmat++;


  // setting numbers of planar, cyllinder, conical and realistic ladder  shapes and extrapolation surfaces
  fkddes_.npmat = Npmat;
  fkddes_.ncmat = Ncmat;
  fkddes1_.nconmat = Nconmat;
  fkexts_.nexs = Nexs;

  //Set to 0 for no planar Material
  fkddes2_.nplmat = Nplmat;
  //fkddes2_.nplmat = 0;

  if (_useMaterials == 0) {
    fkddes_.npmat = 0;
    fkddes_.ncmat = 0;
    fkddes1_.nconmat = 0; 
    fkddes2_.nplmat = 0;   
  }
    
  if (_useExtrapolations == 0)
    fkexts_.nexs = 0;

  // Call to FORTRAN Routine
  setmat_();

}

void MaterialDB::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  //  _bField = Global::parameters->getFloatVal("BField");
  _bField = float(Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z());
  fkfild_.consb = 2.997924e-3*_bField;
} 

void MaterialDB::processEvent( LCEvent * evt ) { 
  _nEvt ++ ;
  _bField = float(Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z());
  fkfild_.consb = 2.997924e-3*_bField;
}



void MaterialDB::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void MaterialDB::end(){ 
  //   std::cout << "MaterialDB::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
}

