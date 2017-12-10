/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#define ROOT_DEBUG 1
/**
If you don't want to generate information of position resolution of generated TrackerHits,
set 0.
*/

#include "FPCCDClustering.h"
#include "FPCCDPixelHit.h"
#include "FPCCDData.h"

#include <iostream>

#include <EVENT/MCParticle.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <gear/ZPlanarParameters.h>
#include <gear/ZPlanarLayerLayout.h>
#include <UTIL/CellIDEncoder.h>
#include "UTIL/LCTOOLS.h"

#include <ILDCellIDEncoding.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>

#include <gsl/gsl_randist.h>

#include <cmath>
#include <algorithm>
#include <sstream>
#include <map>
#include <vector>
#include <list>
#include <utility>
#include <stack>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <TVector3.h>


#include "CLHEP/Vector/LorentzVector.h"
#include <marlin/Global.h>
#include <EVENT/Track.h>
#include <gear/GEAR.h>
#include <gear/BField.h>


using namespace lcio ;
using namespace marlin ;
using namespace std ;

FPCCDClustering aFPCCDClustering ;

// =====================================================================
FPCCDClustering::FPCCDClustering() : Processor("FPCCDClustering") {

  // modify processor description
  _description = "FPCCDClustering icteats TrackerHits from FPCCDPixelHits" ;

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "Debug",
      "Debugging option",
      _debug,
      int(0)); 

  registerProcessorParameter( "EnergyDigitizatoin",
      "switch for digitization of energy deposit",
      _energyDigitization,
      int(1)); 

  registerProcessorParameter( "RandomNoise",
      "Random noise option",
      _randomNoise,
      int(0)); 

  registerProcessorParameter( "FPCCD_PixelSize" ,
      "Pixel size of FPCCD (unit:mm) (default: 0.005)",
      _pixelSize,
      float(0.0050) ) ;

  // Start of New edit
  //<New variable is 'PiXelSizeVec'(so far this variable wasn't vector). So I've replaced PixelSizeVec with something which there were. And I
  // don't write here the place. Please search the word 'PixelSizeVec'.

  FloatVec PixelSizeVec; // This is a temporary variable. Each ladder's PixelSize are different.
  PixelSizeVec.push_back(0.005);
  PixelSizeVec.push_back(0.005);
  PixelSizeVec.push_back(0.010);
  PixelSizeVec.push_back(0.010);
  PixelSizeVec.push_back(0.010);
  PixelSizeVec.push_back(0.010);

  registerProcessorParameter( "Each_FPCCD_pixelSize(mm)",
      "Each ladder's Pixel size of FPCCD (unit:mm) (default:0.005)",
      _pixelSizeVec,
      PixelSizeVec );
  // End of New edit


  registerProcessorParameter( "PixelHeight" , 
      "Pixel Height(mm)",
      _pixelheight ,
      float(0.015));

  registerProcessorParameter( "PointResolutionRPhi" ,
      "R-Phi Resolution in VTX"  ,
      _pointResoRPhi ,
      float(0.001440)) ; // 5/sqrt(12)

  registerProcessorParameter( "PointResolutionZ" ,
      "Z Resolution in VTX" ,
      _pointResoZ ,
      float(0.001440));   // 5/sqrt(12)

  registerProcessorParameter("positionReso_ReadingFile_ON",
      "true: reading file for allocation of position reso. false: rough allocation of position reso. ",
      _positionReso_ReadingFile_ON,
      bool(true));

  registerProcessorParameter("positionReso_ReadingFile",
      "On: reading file for allocation of position reso.",
      _positionReso_ReadingFile,
    std::string("MarlinReco/v01-05/TrackDigi/FPCCDDigi/src/FPCCDClustering_ResoMap_MS_EL_ON_ver0.dat"));

  registerProcessorParameter("noiseElectronSigma",
      "gaussian sigma of number of electron of noise",
      _electronNoiseRate,
      double(50));

  registerProcessorParameter("NumberOfBitsForEnergyDeposit",
      "Number of Bits used for energy depsit",
      _nbitsForEdep,
      int(7));

  registerProcessorParameter("NumberOfElectronForStep",
      "Number of electron for 1 step",
      _electronsPerStep,
      int(25));

  registerProcessorParameter("Threshold",
      "Cell Threshold in electrons",
      _threshold,
      double(200));

  registerProcessorParameter("ElectronsPerKeV",
      "Electrons per keV",
      _electronsPerKeV,
      double(276)); // 3.62eV/1electron


  /*
     registerProcessorParameter( "NewTrackingSystem",
     "Use new tracking system",
     _new_tracking_system,
     bool(true)); 
   */

  registerProcessorParameter( "ClusterRejection_1stCut_isActive",
      "If true, 1stCut is active. false makes this inactive. minus value makes cut on the layer inactive ",
      _firstCut.isActive,
      bool(false)); 

  IntVec ZWidthVec(6);
  ZWidthVec[0] = 15; 
  ZWidthVec[1] = 15; 
  ZWidthVec[2] = 8; 
  ZWidthVec[3] = 8; 
  ZWidthVec[4] = 8; 
  ZWidthVec[5] = 8; 
  registerProcessorParameter( "ClusterRejection_1stCut_ZWidth",
      "Clusters with ZWidth >= this value are discarded. minus value makes cut on the layer inactive",
      _firstCut.ZWidth,
      ZWidthVec); 

  IntVec RPhiWidthVec(6);
  RPhiWidthVec[0] = 10; 
  RPhiWidthVec[1] = 10; 
  RPhiWidthVec[2] = 6; 
  RPhiWidthVec[3] = 6; 
  RPhiWidthVec[4] = 6; 
  RPhiWidthVec[5] = 6; 
  registerProcessorParameter( "ClusterRejection_1stCut_RPhiWidth",
      "Clusters with RPhiWidth >= this value are discarded. minus value makes cut on the layer inactive ",
      _firstCut.RPhiWidth,
      RPhiWidthVec); 

  IntVec nPixVec(6);
  nPixVec[0] = 20; 
  nPixVec[1] = 20; 
  nPixVec[2] = 15; 
  nPixVec[3] = 15; 
  nPixVec[4] = 15; 
  nPixVec[5] = 15; 
  registerProcessorParameter( "ClusterRejection_1stCut_nPix",
      "Clusters with nPix >= this value are discarded. minus value makes cut on the layer inactive",
      _firstCut.nPix,
      nPixVec); 

  registerProcessorParameter( "ClusterRejection_Mori2ndCut_isActive",
      "true: active, false: inactive",
      _m2Cut.isActive,
      bool(false)); 

  FloatVec ZParVec(6);
  ZParVec[0] = 90; 
  ZParVec[1] = 90; 
  ZParVec[2] = 280; 
  ZParVec[3] = 280; 
  ZParVec[4] = 600; 
  ZParVec[5] = 600; 
  registerProcessorParameter( "ClusterRejection_Mori2ndCut_areaZ_parameter",
      "You can choose float[0,~900]. 0 discards most of bad clusters, although some clusters are sacrified. minus value makes cut on the layer inactive",
      _m2Cut.zpar,
      ZParVec); 


  registerProcessorParameter( "ClusterRejection_Kamai2ndCut_isActive",
      "true: active, false: inactive",
      _k2Cut.isActive,
      bool(false)); 

  FloatVec BParVec(6);
  BParVec[0] = 2; 
  BParVec[1] = 2; 
  BParVec[2] = 1; 
  BParVec[3] = 1; 
  BParVec[4] = 1; 
  BParVec[5] = 1; 
  registerProcessorParameter( "ClusterRejection_Kamai2ndCut_b_parameter",
      "You can choose float[0,10]. less than 2 discards most of bad clusters, although some clusters are sacrified. minus value makes cut on the layer inactive",
      _k2Cut.bpar,
      BParVec); 

  IntVec MinZWidthVec(6);
  MinZWidthVec[0] = 4; 
  MinZWidthVec[1] = 4; 
  MinZWidthVec[2] = 2; 
  MinZWidthVec[3] = 2; 
  MinZWidthVec[4] = 2; 
  MinZWidthVec[5] = 2; 
  registerProcessorParameter( "ClusterRejection_Kamai2ndCut_Save_minimum_Cluster_Z_Width_parameter",
      "You can choose int[0,~7]. If a cluster with ZWidth <= this value, this is not discarded.",
      _k2Cut.minZWidth,
      MinZWidthVec); 

  registerProcessorParameter( "MakeRelationOftrkhitsAndsimthits",
      "If false, Clustering doesn't make LCRelation.",
      _makeRelation,
      bool(true));


  registerProcessorParameter( "RemovePixelHitsCollection",
      "Remove VTX Pixel Hits collection after processing",
      _remove_pixelhits_collection,
      bool(false)); 

  // Input collections
  registerInputCollection( LCIO::SIMTRACKERHIT,
      "VTXCollectionName" , 
      "Name of the VTX SimTrackerHit collection"  ,
      _colNameSTH ,
      std::string("VXDCollection") ) ;

  registerInputCollection( LCIO::LCGENERICOBJECT,
      "VTXPixelHitCollection" , 
      "Name of the VTX PixelHit collection"  ,
      _colNameVTX ,
      std::string("VTXPixelHits") ) ;

  // Output collections
  registerOutputCollection( LCIO::TRACKERHIT,
      "VTXHitCollection" ,
      "Name of the vxd TrackerHit output collection"  ,
      _outColNameVTX ,
      std::string("VXDTrackerHits") ) ;//old: VXDTrackerHits

  registerOutputCollection( LCIO::LCRELATION,
      "SimTrackerHitRelCollection" ,
      "Name of TrackerHit SimTrackerHit relation collection"  ,
      _outRelColNameVTX ,
      std::string("VXDTrackerHitRelations") ) ;//old: VTXTrackerHitRelations
#if ROOT_DEBUG
  registerProcessorParameter( "outputRootFileName",
      "output root file name & path (default: ./positionReso.root)",
      _rootFileName,
      std::string("./positionReso.root")); 

  registerProcessorParameter( "treeName",
      "tree name which belongs to the above output root file name & path(default: tree)",
      _treeName,
      std::string("t")); 
#endif

}


// =====================================================================
void FPCCDClustering::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  InitGeometry();

  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_default_seed = _ranSeed;

#if ROOT_DEBUG  
  _rootf = new TFile(_rootFileName.c_str(), "RECREATE");
  _tree = new TTree(_treeName.c_str()," This is tree_name ");


  _tree->Branch("nEvt",&_nEvt,"nEvt/I");
  _tree->Branch("nlinks",&_link.nlink,"nlinks/I");
  _tree->Branch("weight",_link.weight,"weight[nlinks]/F");
  _tree->Branch("simthits_x",_simthits.x,"simthits_x[nlinks]/D");
  _tree->Branch("simthits_y",_simthits.y,"simthits_y[nlinks]/D");
  _tree->Branch("simthits_z",_simthits.z,"simthits_z[nlinks]/D");
  _tree->Branch("simthits_R",_simthits.R,"simthits_R[nlinks]/D");
  _tree->Branch("simthits_pdg",_simthits.pdg,"simthits_pdg[nlinks]/I");
  _tree->Branch("simthits_mass",_simthits.mass,"simthits_mass[nlinks]/D");
  _tree->Branch("simthits_vx",_simthits.vx,"simthits_vx[nlinks]/D");
  _tree->Branch("simthits_vy",_simthits.vy,"simthits_vy[nlinks]/D");
  _tree->Branch("simthits_vz",_simthits.vz,"simthits_vz[nlinks]/D");
  _tree->Branch("simthits_px",_simthits.px,"simthits_px[nlinks]/F");
  _tree->Branch("simthits_py",_simthits.py,"simthits_py[nlinks]/F");
  _tree->Branch("simthits_pz",_simthits.pz,"simthits_pz[nlinks]/F");
  _tree->Branch("simthits_pAbs",_simthits.pAbs,"simthits_pAbs[nlinks]/D");
  _tree->Branch("simthits_theta",_simthits.theta,"simthits_theta[nlinks]/D");
  _tree->Branch("simthits_phi",_simthits.phi,"simthits_phi[nlinks]/D");
  _tree->Branch("simthits_area_theta",_simthits.area_theta,"simthits_area_theta[nlinks]/D");
  _tree->Branch("simthits_area_phi",_simthits.area_phi,"simthits_area_phi[nlinks]/D");
  _tree->Branch("simthits_layer",_simthits.layer,"simthits_layer[nlinks]/I");
  _tree->Branch("simthits_ladder",_simthits.ladder,"simthits_ladder[nlinks]/I");
  _tree->Branch("simthits_mcp_px",_simthits.mcp_px,"simthits_mcp_px[nlinks]/D");
  _tree->Branch("simthits_mcp_py",_simthits.mcp_py,"simthits_mcp_py[nlinks]/D");
  _tree->Branch("simthits_mcp_pz",_simthits.mcp_pz,"simthits_mcp_pz[nlinks]/D");
  _tree->Branch("simthits_mcp_energy",_simthits.mcp_energy,"simthits_mcp_energy[nlinks]/D");
  _tree->Branch("simthits_mcp_time",_simthits.mcp_time,"simthits_mcp_time[nlinks]/F");
  _tree->Branch("simthits_mcp_Pt",_simthits.mcp_Pt,"simthits_mcp_Pt[nlinks]/D");
  _tree->Branch("simthits_mcp_isCreatedInSimulation",_simthits.mcp_isCreatedInSimulation,"simthits_mcp_isCreatedInSimulation[nlinks]/I");
  _tree->Branch("simthits_mcp_isBackscatter",_simthits.mcp_isBackscatter,"simthits_mcp_isBackscatter[nlinks]/I");
  _tree->Branch("simthits_mcp_vertexIsNotEndpointOfParent",_simthits.mcp_vertexIsNotEndpointOfParent,"simthits_mcp_vertexIsNotEndpointOfParent[nlinks]/I");
  _tree->Branch("simthits_mcp_isDecayedInTracker",_simthits.mcp_isDecayedInTracker,"simthits_mcp_isDecayedInTracker[nlinks]/I");
  _tree->Branch("simthits_mcp_isDecayedInCalorimeter",_simthits.mcp_isDecayedInCalorimeter,"simthits_mcp_isDecayedInCalorimeter[nlinks]/I");
  _tree->Branch("simthits_mcp_hasLeftDetector",_simthits.mcp_hasLeftDetector,"simthits_mcp_hasLeftDetector[nlinks]/I");
  _tree->Branch("simthits_mcp_isStopped",_simthits.mcp_isStopped,"simthits_mcp_isStopped[nlinks]/I");
  _tree->Branch("simthits_mcp_d0",_simthits.mcp_d0,"simthits_mcp_d0[nlinks]/F");
  _tree->Branch("simthits_mcp_z0",_simthits.mcp_z0,"simthits_mcp_z0[nlinks]/F");
  _tree->Branch("simthits_mcp_phi0",_simthits.mcp_phi0,"simthits_mcp_phi0[nlinks]/F");
  _tree->Branch("simthits_mcp_omega",_simthits.mcp_omega,"simthits_mcp_omega[nlinks]/F");
  _tree->Branch("simthits_mcp_tanL",_simthits.mcp_tanL,"simthits_mcp_tanL[nlinks]/F");
  _tree->Branch("simthits_xi",_simthits.xi,"simthits_xi[nlinks]/D");
  _tree->Branch("simthits_zeta",_simthits.zeta,"simthits_zeta[nlinks]/D");
  _tree->Branch("simthits_edep",_simthits.edep,"simthits_edep[nlinks]/F");





  _tree->Branch("trkhits_x",_trkhits.x,"trkhits_x[nlinks]/D");
  _tree->Branch("trkhits_y",_trkhits.y,"trkhits_y[nlinks]/D");
  _tree->Branch("trkhits_z",_trkhits.z,"trkhits_z[nlinks]/D");
  _tree->Branch("trkhits_R",_trkhits.R,"trkhits_R[nlinks]/D");
  _tree->Branch("trkhits_area_theta",_trkhits.area_theta,"trkhits_area_theta[nlinks]/D");
  _tree->Branch("trkhits_area_phi",_trkhits.area_phi,"trkhits_area_phi[nlinks]/D");
  _tree->Branch("trkhits_CWidth_RPhi",_trkhits.CWidth_RPhi,"trkhits_CWidth_RPhi[nlinks]/i");
  _tree->Branch("trkhits_CWidth_Z",_trkhits.CWidth_Z,"trkhits_CWidth_Z[nlinks]/i");
  _tree->Branch("trkhits_tilt",_trkhits.tilt,"trkhits_tilt[nlinks]/i");
  _tree->Branch("trkhits_nPix",_trkhits.nPix,"trkhits_nPix[nlinks]/i");
  _tree->Branch("trkhits_layer",_trkhits.layer,"trkhits_layer[nlinks]/I");
  _tree->Branch("trkhits_ladder",_trkhits.ladder,"trkhits_ladder[nlinks]/I");
  _tree->Branch("trkhits_xi",_trkhits.xi,"trkhits_xi[nlinks]/D");
  _tree->Branch("trkhits_zeta",_trkhits.zeta,"trkhits_zeta[nlinks]/D");
  _tree->Branch("trkhits_edep",_trkhits.edep,"trkhits_edep[nlinks]/F");

  _tree->Branch("trkhits_tposX",_trkhits.tposX,"trkhits_tposX[nlinks]/D");
  _tree->Branch("trkhits_tposY",_trkhits.tposY,"trkhits_tposY[nlinks]/D");
  _tree->Branch("trkhits_tposZ",_trkhits.tposZ,"trkhits_tposZ[nlinks]/D");
  _tree->Branch("trkhits_cx",_trkhits.cx,"trkhits_cx[nlinks]/D");
  _tree->Branch("trkhits_cy",_trkhits.cy,"trkhits_cy[nlinks]/D");
  _tree->Branch("trkhits_cz",_trkhits.cz,"trkhits_cz[nlinks]/D");
  _tree->Branch("trkhits_dot",_trkhits.dot,"trkhits_dot[nlinks]/D");

  _tree->Branch("diff_RPhi",_diff.RPhi,"diff_RPhi[nlinks]/D");
  _tree->Branch("diff_Z",_diff.Z,"diff_Z[nlinks]/D");
#endif






  if(_positionReso_ReadingFile_ON == true){
    _fin.open( _positionReso_ReadingFile.c_str() );
    if(!_fin.is_open()){
      std::cout << "========FATAL ERROR============================" << std::endl;
      std::cout << "positionReso_ReadingFile is not opened!! exit!!" << std::endl;
      std::cout << _positionReso_ReadingFile                         << std::endl;
      std::cout << "The above path seems to be wrong.              " << std::endl;
      std::cout << "FPCCDClustering exits!                         " << std::endl;
      std::cout << "===============================================" << std::endl;
      exit(10);
    }

    float a = 0;
    for(int i = 0 ; i < 6 ; i++){
      for(int j = 1 ; j <= 49 ; j++){
        _fin >> a ;
        short int code = i * 100 + j;
        _resolutionMapRPhi.insert(std::make_pair(code,a));
      }
    }
    for(int i = 0 ; i < 6 ; i++){
      for(int j = 1 ; j <= 49 ; j++){
        _fin >> a ;
        short int code = i * 100 + j;
        _resolutionMapZ.insert(std::make_pair(code,a));
      }
    }
  }
}

// =====================================================================
void FPCCDClustering::InitGeometry() 
{ 
  // Save frequently used parameters.

  const gear::ZPlanarParameters &gearVXD = Global::GEAR->getVXDParameters();
  const gear::ZPlanarLayerLayout &layerVXD=gearVXD.getVXDLayerLayout();

  _nLayer = layerVXD.getNLayers() ;
  _geodata.resize(_nLayer);
  _maxLadder = 0;

  for(int ly=0;ly<_nLayer;ly++){
    _geodata[ly].nladder = layerVXD.getNLadders(ly);     // Number of ladders in this layer
    if( _maxLadder < _geodata[ly].nladder ) { _maxLadder = _geodata[ly].nladder; }
    _geodata[ly].rmin = layerVXD.getSensitiveDistance(ly); // Distance of sensitive area from IP
    _geodata[ly].dphi = (2*M_PI)/(double)_geodata[ly].nladder;
    _geodata[ly].phi0 = layerVXD.getPhi0(ly);  // phi offset.
    _geodata[ly].sthick = layerVXD.getSensitiveThickness(ly);
    _geodata[ly].sximin = -layerVXD.getSensitiveOffset(ly)
      -layerVXD.getSensitiveWidth(ly)/2.0;
    _geodata[ly].sximax = -layerVXD.getSensitiveOffset(ly)
      +layerVXD.getSensitiveWidth(ly)/2.0;
    _geodata[ly].hlength = layerVXD.getSensitiveLength(ly);
    _geodata[ly].cosphi.resize( _geodata[ly].nladder ) ;
    _geodata[ly].sinphi.resize( _geodata[ly].nladder ) ;
    _geodata[ly].ladder_incline.resize( _geodata[ly].nladder ) ;
    _geodata[ly].num_xi_pixel = (int)(layerVXD.getSensitiveWidth(ly)/_pixelSizeVec[ly]);
    _geodata[ly].num_zeta_pixel = (int)(2*layerVXD.getSensitiveLength(ly)/_pixelSizeVec[ly]);
    for(int ld=0;ld<_geodata[ly].nladder;ld++) {
      double phi = _geodata[ly].phi0 + _geodata[ly].dphi*ld;
      _geodata[ly].cosphi[ld] = cos(phi);
      _geodata[ly].sinphi[ld] = sin(phi);
      double incline =  phi - (M_PI/2.0);
      while( incline >  1.0*M_PI || incline < -1.0*M_PI ){ incline > 1.0*M_PI ? incline -= 2.0*M_PI : incline += 2.0*M_PI ;}
      _geodata[ly].ladder_incline[ld] = incline; 
    }
  }
} 


// =====================================================================
void FPCCDClustering::processRunHeader( LCRunHeader*  /*run*/) { 
  _nRun++ ;
} 

// =====================================================================
void FPCCDClustering::modifyEvent( LCEvent * evt )
{ 

  LCCollectionVec* pHitCol=0;
  LCCollection* STHcol=0;
  try { pHitCol = dynamic_cast<LCCollectionVec*>( evt->getCollection( _colNameVTX ) );
    if( _remove_pixelhits_collection ) {  pHitCol->setTransient(true); } }
    catch(DataNotAvailableException &e){
      if (_debug >= 1) { std::cout << "Collection " << _colNameVTX.c_str() << " is unavailable in event " 
        << _nEvt << std::endl; }
    }

    try { STHcol = evt->getCollection( _colNameSTH ) ; }
    catch(DataNotAvailableException &e){
      if (_debug >= 1) { std::cout << "Collection " << _colNameSTH.c_str() << " is unavailable in event " 
        << _nEvt << std::endl; }
    }


    if( _debug >= 1 ) {
      std::cout << " Collection =" << _colNameVTX << " nevt=" << _nEvt << std::endl;
      std::cout << " Collection =" << _colNameSTH << " nevt=" << _nEvt << std::endl;
      if( pHitCol != 0 ) { std::cout << " number of elements is " << pHitCol->getNumberOfElements() << std::endl; }
      if( STHcol != 0 ) { std::cout << " number of elements is " << STHcol->getNumberOfElements() << std::endl; }
    }

    if( pHitCol != 0 && STHcol != 0){    
      FPCCDData theData(_nLayer, _maxLadder);  // prepare object to make pixelhits
      int nhit=theData.unpackPixelHits(*pHitCol);//This unpacks LCGenericObjectVec made in FPCCDDigitizer (precisely, FPCCDData).
      //And return number of pixel hits generated by one simtrackerhit in FPCCDDigitizer.
      if( _debug >= 2 ) { LCTOOLS::dumpEvent( evt ) ;}

      if( _debug >=2 ) { std::cout << " nhit=" << nhit << std::endl; }
      if( nhit > 0 ) {  // Output Trackhit, if there are pixel hits

        LCCollectionVec* trkHitVec = new LCCollectionVec( LCIO::TRACKERHITPLANE );
        LCCollectionVec* relCol = new LCCollectionVec( LCIO::LCRELATION );

        //The flow of Flag setting. 
        LCFlagImpl lcFlag(0);
        lcFlag.setBit( LCIO::LCREL_WEIGHTED );//LCIO::LEREL_WEIGHTED = 31. See LCIO and LCFlagImpl in class list.
        relCol->setFlag( lcFlag.getFlag() );


        makeTrackerHitVec( &theData, STHcol, relCol, trkHitVec);
        evt->addCollection( trkHitVec, _outColNameVTX ) ;
        evt->addCollection( relCol, _outRelColNameVTX ) ;

        if( _debug >= 2 ) {
          std::cout << "dumpEvent FPCCDClustering" << std::endl;
          LCTOOLS::dumpEvent( evt ) ;
        }
        if( _remove_pixelhits_collection ) {
          for(int i = pHitCol->getNumberOfElements(); i>0; i--){
            pHitCol->removeElementAt(i-1);
          }
        }
      }
      if( _remove_pixelhits_collection ) {  evt->removeCollection( _colNameVTX ); }
      theData.clear();

    } // End of process when VXD has hits

#if ROOT_DEBUG

    //******************ROOT Tree Session*****************************************//

    bool Go = true;
    if( pHitCol != 0 && STHcol != 0){    

      LCCollectionVec* rootRelCol=0;
      try { rootRelCol = dynamic_cast<LCCollectionVec*>( evt->getCollection( _outRelColNameVTX ) );}
      catch(DataNotAvailableException &e){
        Go = false;
        std::cout << "miss getCollection!! at _nEvt = " << _nEvt <<  std::endl;
        if (_debug >= 1) { std::cout << "Collection " << _outRelColNameVTX.c_str() << " is unavailable in event " 
          << _nEvt << std::endl; }
      }
      /*
         _tree->Branch("nlinks",&_link.nlink,"nlinks/I");
         _tree->Branch("weight",_link.weight,"weight[nlinks]/F");
         _tree->Branch("simthits_x",_simthits.x,"simthits_x[nlinks]/D");
         _tree->Branch("simthits_y",_simthits.y,"simthits_y[nlinks]/D");
         _tree->Branch("simthits_z",_simthits.z,"simthits_z[nlinks]/D");
         _tree->Branch("simthits_pdg",_simthits.pdg,"simthits_pdg[nlinks]/I");

         _tree->Branch("trkhits_x",_trkhits.x,"trkhits_x[nlinks]/D");
         _tree->Branch("trkhits_y",_trkhits.y,"trkhits_y[nlinks]/D");
         _tree->Branch("trkhits_z",_trkhits.z,"trkhits_z[nlinks]/D");
         _tree->Branch("trkhits_CWidth_RPhi",_trkhits.CWidth_RPhi,"trkhits_CWidth_RPhi[nlinks]/D");
         _tree->Branch("trkhits_CWidth_Z",_trkhits.CWidth_Z,"trkhits_CWidth_Z[nlinks]/D");
       */

      if(Go == true){
        LCFlagImpl forRootlcFlag(0);
        forRootlcFlag.setBit( LCIO::LCREL_WEIGHTED );//LCIO::LEREL_WEIGHTED = 31. See LCIO and LCFlagImpl in class list.
        //std::cout << "check1" << std::endl;
        rootRelCol->setFlag( forRootlcFlag.getFlag() );
        //std::cout << "check2" << std::endl;

        _link.nlink = rootRelCol->getNumberOfElements(); 
        for(unsigned int ri = 0; ri < _link.nlink; ri++ ){
          LCRelation* link_p = dynamic_cast<LCRelation*>( rootRelCol->getElementAt(ri) );
          _link.weight[ri] = link_p->getWeight();
          TrackerHit* trkhit_p = dynamic_cast<TrackerHit*>( link_p->getFrom() );
          SimTrackerHit* simthit_p = dynamic_cast<SimTrackerHit*>( link_p->getTo() );

          _simthits.x[ri] = simthit_p->getPosition()[0];
          _simthits.y[ri] = simthit_p->getPosition()[1];
          _simthits.z[ri] = simthit_p->getPosition()[2];


          _simthits.pdg[ri] = simthit_p->getMCParticle()->getPDG();
          _simthits.mass[ri] = simthit_p->getMCParticle()->getMass();
          _simthits.px[ri] = simthit_p->getMomentum()[0]; 
          _simthits.py[ri] = simthit_p->getMomentum()[1]; 
          _simthits.pz[ri] = simthit_p->getMomentum()[2]; 
          _simthits.mcp_px[ri] = simthit_p->getMCParticle()->getMomentum()[0];
          _simthits.mcp_py[ri] = simthit_p->getMCParticle()->getMomentum()[1];
          _simthits.mcp_pz[ri] = simthit_p->getMCParticle()->getMomentum()[2];
          _simthits.mcp_Pt[ri] = sqrt(pow(_simthits.mcp_px[ri],2) + pow(_simthits.mcp_py[ri],2)); 
          _simthits.vx[ri] = simthit_p->getMCParticle()->getVertex()[0]; 
          _simthits.vy[ri] = simthit_p->getMCParticle()->getVertex()[1]; 
          _simthits.vz[ri] = simthit_p->getMCParticle()->getVertex()[2]; 
          _simthits.mcp_energy[ri] = simthit_p->getMCParticle()->getEnergy(); 
          _simthits.mcp_time[ri] = simthit_p->getMCParticle()->getTime(); 
          _simthits.mcp_isCreatedInSimulation[ri] = simthit_p->getMCParticle()->isCreatedInSimulation(); 
          _simthits.mcp_isBackscatter[ri] = simthit_p->getMCParticle()->isBackscatter(); 
          _simthits.mcp_vertexIsNotEndpointOfParent[ri] = simthit_p->getMCParticle()->vertexIsNotEndpointOfParent(); 
          _simthits.mcp_isDecayedInTracker[ri] = simthit_p->getMCParticle()->isDecayedInTracker(); 
          _simthits.mcp_isDecayedInCalorimeter[ri] = simthit_p->getMCParticle()->isDecayedInCalorimeter(); 
          _simthits.mcp_hasLeftDetector[ri] = simthit_p->getMCParticle()->hasLeftDetector(); 
          _simthits.mcp_isStopped[ri] = simthit_p->getMCParticle()->isStopped(); 













          _trkhits.x[ri] = trkhit_p->getPosition()[0];
          _trkhits.y[ri] = trkhit_p->getPosition()[1];
          _trkhits.z[ri] = trkhit_p->getPosition()[2];

          double signRPhi = _trkhits.x[ri] * _simthits.y[ri] - _trkhits.y[ri] * _simthits.x[ri];

          if(signRPhi > 0){
            _diff.RPhi[ri] = sqrt(pow(_trkhits.x[ri] - _simthits.x[ri],2) + pow(_trkhits.y[ri] - _simthits.y[ri],2));
          }
          else{
            _diff.RPhi[ri] = (-1)*sqrt(pow(_trkhits.x[ri] - _simthits.x[ri],2) + pow(_trkhits.y[ri] - _simthits.y[ri],2));
          }

          _diff.Z[ri] = _trkhits.z[ri] - _simthits.z[ri];


          int cellid1 = trkhit_p->getCellID1();
          unsigned int   xiwidth =  cellid1 & 0x000001ff ;
          unsigned int zetawidth = ( cellid1 >> 9 ) & 0x000001ff;
          unsigned int      nPix = ( cellid1 >> 18) & 0x000003ff;
          unsigned int      tilt = ( cellid1 >> 28) & 0x00000003;
          _trkhits.CWidth_RPhi[ri] = xiwidth;
          _trkhits.CWidth_Z[ri] = zetawidth;
          _trkhits.nPix[ri] = nPix;
          _trkhits.tilt[ri] = tilt;


          _simthits.pAbs[ri] = sqrt(pow(static_cast<double>(_simthits.px[ri]),2) + pow(static_cast<double>(_simthits.py[ri]),2) + pow(static_cast<double>(_simthits.pz[ri]),2));
          _simthits.theta[ri] = 180/M_PI* acos( static_cast<double>(_simthits.pz[ri]) / _simthits.pAbs[ri]  ); //from 0 to pi
          _simthits.phi[ri] = 180/M_PI* atan2( static_cast<double>(_simthits.py[ri]) ,static_cast<double>(_simthits.px[ri])); // from -pi to pi


          _simthits.R[ri] = sqrt(pow(_simthits.x[ri],2) + pow(_simthits.y[ri],2) + pow(_simthits.z[ri],2));
          _simthits.area_theta[ri] = 180/M_PI* acos( _simthits.z[ri] / _simthits.R[ri]  ); //from 0 to pi
          _simthits.area_phi[ri] = 180/M_PI* atan2( _simthits.y[ri] ,_simthits.x[ri]); // from -pi to pi

          _trkhits.R[ri] = sqrt(pow(_trkhits.x[ri],2) + pow(_trkhits.y[ri],2) + pow(_trkhits.z[ri],2));
          _trkhits.area_theta[ri] = 180/M_PI* acos( _trkhits.z[ri] / _trkhits.R[ri]  ); //from 0 to pi
          _trkhits.area_phi[ri] = 180/M_PI* atan2( _trkhits.y[ri] ,_trkhits.x[ri]); // from -pi to pi








          const int cellId_sim = simthit_p->getCellID0();
          UTIL::BitField64 encoder_sim( lcio::LCTrackerCellID::encoding_string() ) ; 
          encoder_sim.setValue(cellId_sim) ;
          _simthits.layer[ri]    = encoder_sim[lcio::LCTrackerCellID::layer()] ;
          _simthits.ladder[ri] = encoder_sim[lcio::LCTrackerCellID::module()] ;

          const int cellId_trk = trkhit_p->getCellID0();
          UTIL::BitField64 encoder_trk( lcio::LCTrackerCellID::encoding_string() ) ; 
          encoder_trk.setValue(cellId_trk) ;
          _trkhits.layer[ri]    = encoder_trk[lcio::LCTrackerCellID::layer()] ;
          _trkhits.ladder[ri] = encoder_trk[lcio::LCTrackerCellID::module()] ;

          double lpos[3];//xi,eta,zeta
          double gpos[3];
          gpos[0] = _trkhits.x[ri];
          gpos[1] = _trkhits.y[ri];
          gpos[2] =_trkhits.z[ri];
          int tlayer = _trkhits.layer[ri];
          int tladder = _trkhits.ladder[ri];
          double sinphi = _geodata[tlayer].sinphi[tladder];
          double cosphi = _geodata[tlayer].cosphi[tladder];
          lpos[0] = gpos[0] * sinphi - gpos[1] * cosphi + gpos[2] * 0;
          lpos[1] = gpos[0] * cosphi + gpos[1] * sinphi + gpos[2] * 0;        
          lpos[2] = gpos[0] * 0      + gpos[1] * 0      + gpos[2] * 1;
          _trkhits.xi[ri] = lpos[0];
          _trkhits.zeta[ri] = lpos[2];

          gpos[0] = _simthits.x[ri];
          gpos[1] = _simthits.y[ri];
          gpos[2] =_simthits.z[ri];
          tlayer = _simthits.layer[ri];
          tladder = _simthits.ladder[ri];
          sinphi = _geodata[tlayer].sinphi[tladder];
          cosphi = _geodata[tlayer].cosphi[tladder];
          lpos[0] = gpos[0] * sinphi - gpos[1] * cosphi + gpos[2] * 0;
          lpos[1] = gpos[0] * cosphi + gpos[1] * sinphi + gpos[2] * 0;        
          lpos[2] = gpos[0] * 0      + gpos[1] * 0      + gpos[2] * 1;
          _simthits.xi[ri] = lpos[0];
          _simthits.zeta[ri] = lpos[2];

          double par[5];
          calcTrackParameterOfMCP(simthit_p->getMCParticle(), par);
          _simthits.mcp_d0[ri] = par[0];
          _simthits.mcp_z0[ri] = par[1];
          _simthits.mcp_omega[ri] = par[2];
          _simthits.mcp_phi0[ri] = par[3];
          _simthits.mcp_tanL[ri] = par[4];


          /////////////////////Calc of Dot of cluster and global position//////////////////////

#if 0
          const double xioffset[6] = {-5.5+1.575320104,
            -5.5+1.575320104,
            -11+1.531923469,
            -11+1.531923469,
            -11+2.293817688,
            -11+2.293817688};

          const double layerdistance[6] = { 15.9575,
            18.0425,
            36.9575,
            39.0425,
            57.9575,
            60.0425};
          const double hlength[6] = {62.5,62.5,125,125,125,125};

          std::cout << "xioffset : " << std::endl;
          for(int jj = 0; jj < 6; jj++) std::cout << xioffset[jj] << "  " << std::endl;
          std::cout << "layerdistance : " << std::endl;
          for(int jj = 0; jj < 6; jj++) std::cout << layerdistance[jj] << "  " << std::endl;
          std::cout << "hlength : " << std::endl;
          for(int jj = 0; jj < 6; jj++) std::cout << hlength[jj] << "  " << std::endl;

          std::cout << " _geodata[ly].sximin  : "  <<  std::endl;
          for(int jj = 0; jj < 6; jj++) std::cout <<  _geodata[jj].sximin << "  " << std::endl;
          std::cout << " _geodata[ly].rmin  : "    << std::endl;
          for(int jj = 0; jj < 6; jj++) std::cout <<  _geodata[jj].rmin << "  " << std::endl;
          std::cout << " _geodata[ly].hlength  : " << std::endl;
          for(int jj = 0; jj < 6; jj++) std::cout <<  _geodata[jj].hlength << "  " << std::endl;
#endif

          int mly = _trkhits.layer[ri];
          int mld = _trkhits.ladder[ri];
          gear::Vector3D tpos(trkhit_p->getPosition()[0],trkhit_p->getPosition()[1],trkhit_p->getPosition()[2]);

          /*
             double TrkX = tpos.x();
             double TrkY = tpos.y();
             double TrkZ = tpos.z();
           */
          double   tposR = tpos.rho();
          double tposphi = tpos.phi();

          double TrkXi = tposR*(cos(tposphi)*_geodata[mly].sinphi[mld]-sin(tposphi)*_geodata[mly].cosphi[mld]) - _geodata[mly].sximin;
          double TrkZeta = tpos.z()+_geodata[mly].hlength;

          _trkhits.edep[ri] = trkhit_p->getEDep();
          _simthits.edep[ri] = simthit_p->getEDep();
          TVector3 tposDir(TrkXi+_geodata[mly].sximin,TrkZeta-_geodata[mly].hlength,_geodata[mly].rmin);
          _trkhits.tposX[ri] = tposDir.X();
          _trkhits.tposY[ri] = tposDir.Y();
          _trkhits.tposZ[ri] = tposDir.Z();

          double sign = 1;
          if(tilt>1) sign = -1;
          float pixelsize = _pixelSizeVec[mly];

          TVector3 ClusterDir(pixelsize*(xiwidth-1),sign*pixelsize*(zetawidth-1),_pixelheight);
          TVector3 revClusterDir(-ClusterDir.X(),-ClusterDir.Y(),ClusterDir.Z());
          //  gear::Vector3D ClusterDir(sign*((double)xiwidth-1),(double)zetawidth-1,0);
          if( revClusterDir.Unit().Dot( tposDir.Unit() ) > ClusterDir.Unit().Dot( tposDir.Unit() )){
            ClusterDir.SetX(revClusterDir.X());
            ClusterDir.SetY(revClusterDir.Y());
          }
          _trkhits.cx[ri] = ClusterDir.Unit().X();
          _trkhits.cy[ri] = ClusterDir.Unit().Y();
          _trkhits.cz[ri] = ClusterDir.Unit().Z();

          _trkhits.dot[ri] = ClusterDir.Unit().Dot( tposDir.Unit() );//global座標とclusterの内積

        }


        //std::cout << "fill" << std::endl;
        _tree->Fill();
        //std::cout << "fill clear " << std::endl;
      }
    }
#endif

    //***********************ROOT Tree Session end*********************************//

    _nEvt ++ ;
}

// ====================================================================
void FPCCDClustering::makeTrackerHitVec(FPCCDData* pHitData, LCCollection* STHcol, LCCollectionVec* relCol,  LCCollectionVec* trkHitVec)
{
  //At first, relCol and trkHitVec are void data.
  typedef std::multimap< std::pair<int,int>, SimTrackerHit*> RelMap_t;
  //std::multimap< std::pair<int,int>, SimTrackerHit*> relMap;//"multimap" allows duplication of Key and "map" doesn't.
  RelMap_t relMap;//"multimap" allows duplication of Key and "map" doesn't.

  relMap.clear();



  // Convert pixel hits to TrackerHits
  for(int layer=0;layer<_nLayer;layer++){
    int     nladder = _geodata[layer].nladder;
    double   sximin = _geodata[layer].sximin;
    double   sximax = _geodata[layer].sximax;
    double sxiwidth = sximax-sximin;
    double  hlength = _geodata[layer].hlength;

    // ----------------------------------------
    // Start clustering of hits in each laddder
    // ----------------------------------------
    for(int ladder=0;ladder<nladder;ladder++){

      FPCCDLadderHit_t ladderHit; //typedef std::map<FPCCDHitLoc_t, FPCCDPixelHit*>  FPCCDLadderHit_t;

      //typedef std::pair<unsigned int, unsigned int> FPCCDHitLoc_t;
      //typedef std::map<FPCCDHitLoc_t, FPCCDPixelHit*>  FPCCDLadderHit_t;
      // namely, std::map<std::pair<unsigned int, unsigned int>, FPCCDPixelHit*>



      PixelHitMap_t::iterator it=pHitData->itBegin(layer, ladder);      
      /*
         is defined as below in FPCCDData.h.
         typedef std::map<unsigned int, FPCCDPixelHit*> PixelHitMap_t; //key will be encorded cellID. 
         typedef std::vector< std::vector<PixelHitMap_t> > PixelDataBuf_t;
       */

      //(1) Copy all hits into local variables, ladderHit
      while( it != pHitData->itEnd(layer, ladder) ) {

        FPCCDPixelHit *aHit=(*it).second;
        // Set the noise
        aHit->setEdep( aHit->getEdep() + gsl_ran_gaussian(_rng, _electronNoiseRate/_electronsPerKeV)*1e-6);//_electronNoiseRate = 50 (by default)
        if( _energyDigitization){
          EnergyDigitizer( aHit );
        }
        if( aHit->getEdep() > ( _threshold/_electronsPerKeV ) * 1e-6 ) {  // Do the threshold cut

          int xiID=aHit->getXiID();
          int zetaID=aHit->getZetaID();
          FPCCDHitLoc_t hitloc(xiID, zetaID);
          ladderHit.insert(map<FPCCDHitLoc_t, FPCCDPixelHit*>::value_type(hitloc, aHit));
        }
        it++;
      }
      // Set the random noise hit
      if( _randomNoise && _energyDigitization){
        unsigned int noiseXi = 0;
        unsigned int noiseZeta = 0;
        double noiseEdep = 0;

        for(unsigned int i=0; i< gsl_ran_poisson( _rng, (sxiwidth/_pixelSizeVec[layer])*(2*hlength/_pixelSizeVec[layer])*3.1671*1e-5) ; i++){ // Random noise is generated by 4 sigma plobability.

          noiseXi = (unsigned int)gsl_ran_flat( _rng, 0, sxiwidth/_pixelSizeVec[layer] );
          noiseZeta = (unsigned int)gsl_ran_flat( _rng, 0, 2*hlength/_pixelSizeVec[layer] );          
          FPCCDHitLoc_t noiseHitLoc( noiseXi, noiseZeta );
          FPCCDLadderHit_t::iterator noiseIt=ladderHit.find( noiseHitLoc );

          if( noiseIt == ladderHit.end() ){
            noiseEdep = gsl_ran_gaussian_tail( _rng, _threshold, _electronNoiseRate )/(_electronsPerKeV * 1e+6 );
            FPCCDPixelHit* noiseHit = new FPCCDPixelHit(layer, ladder, noiseXi, noiseZeta, noiseEdep, FPCCDPixelHit::kBKG, 0);
            EnergyDigitizer( noiseHit);            
            ladderHit.insert( map<FPCCDHitLoc_t, FPCCDPixelHit*>::value_type( noiseHitLoc, noiseHit) );
          }          
        }
      }
      // Do clustering
      if( ladderHit.size() > 0 ) {
        FPCCDClusterVec_t clusterVec;
        //typedef std::vector<FPCCDPixelHit*> FPCCDCluster_t;
        //typedef std::vector<FPCCDCluster_t*> FPCCDClusterVec_t;
        //namely, std::vector<std::vector<FPCCDPixelHit*>>

        makeClustersInALadder(layer, ladderHit, clusterVec);
        if( _debug >= 1 ) { 
          std::cout << "Layer:" << layer << " ladder:" << ladder << " # of clusters=" << clusterVec.size() << std::endl;
          std::cout << "  # pixels in each clusters are "; 
          for(unsigned int i=0;i<clusterVec.size();i++) {
            std::cout << clusterVec[i]->size() << " " ;
            if ( i > 1000 ) { std::cout << " ... omitting the rest. " ; break; }
          }
          std::cout << endl;
        }
        // Calulate position from cluster and 
        makeTrackerHit(STHcol, layer, ladder, clusterVec, relMap, relCol, trkHitVec);
        //At this phase, relCol and trkHitVec haven't been filled and changed.

        // Now clean up clusters in clusterVec;
        for(int i=clusterVec.size()-1;i>=0;i--){ delete clusterVec[i]; }
      }
    } // End of Ladder loop
  } // End of Layer loop
}




void FPCCDClustering::makeTrackerHit(LCCollection* STHcol, int layer, int ladder, FPCCDClusterVec_t &clusterVec, std::multimap< std::pair<int,int>, SimTrackerHit*> relMap, LCCollectionVec* relCol, LCCollectionVec* trkHitVec )
{
  //trkHitVec and relCol is yet void data. 
  //This Version is different from default version at the point of last area of this function scope.
  //The way of making relation between simthits and trkhits is different.

  //typedef std::vector<FPCCDPixelHit*> FPCCDCluster_t;
  //typedef std::vector<FPCCDCluster_t*> FPCCDClusterVec_t;
  //namely, std::vector<std::vector<FPCCDPixelHit*>>



  CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( lcio::LCTrackerCellID::encoding_string() , trkHitVec ) ;
  for(unsigned int ic=0;ic<clusterVec.size();ic++) {
    FPCCDCluster_t *cluster=clusterVec[ic];

    double   xiene = 0.0; double zetaene = 0.0; double  enesum = 0.0;
    int      maxxi = 0; int      minxi = 0; int    maxzeta = 0; int    minzeta = 0;
    FPCCDPixelHit::HitQuality_t trackquality = FPCCDPixelHit::kSingle;
    /* is defined as below in FPCCDPixelHit.h.
       typedef enum { kSingle=0, kSignalOverlap=0x01, kBKGOverlap=0x02, kBKG=0x03} HitQuality_t;
     */

    FPCCDPixelHit   *maxXiHit=(*cluster)[0]; FPCCDPixelHit   *minXiHit=(*cluster)[0];
    FPCCDPixelHit *maxZetaHit=(*cluster)[0]; FPCCDPixelHit *minZetaHit=(*cluster)[0];

    unsigned int nPix=cluster->size();    
    std::list<int> orderVecOfCluster(0); //Mori added 20121203
    for(unsigned int i=0;i<nPix;i++) {
      /***Mori added 20121203***/
      int overlaidSize = (*cluster)[i]->getSizeOfOrderID();//Sometimes one pixel is generated by more than one simthits. 
      //So, usually overlaidSize is one, but sometimes 2 or more.
      for(int ios = 0; ios < overlaidSize; ios++){
        orderVecOfCluster.push_back( (*cluster)[i]->getOrderID(ios) ); //Mori added 20121204
      }
      /********end**************/

      FPCCDPixelHit *aHit=(*cluster)[i]; trackquality = aHit->getQuality(); enesum+=aHit->getEdep();
      xiene+=((double)(aHit->getXiID()+0.5))*_pixelSizeVec[layer]*aHit->getEdep();
      zetaene+=((double)(aHit->getZetaID()+0.5))*_pixelSizeVec[layer]*aHit->getEdep();

      if(i==0){
        maxxi = aHit->getXiID(); minxi = maxxi;
        maxzeta = aHit->getZetaID(); minzeta = maxzeta;
      }
      else{
        if(  aHit->getXiID() > maxxi){ maxxi = aHit->getXiID(); maxXiHit=aHit;}
        if(  aHit->getXiID() < minxi){ minxi = aHit->getXiID(); minXiHit=aHit;}
        if(aHit->getZetaID() > maxzeta){ maxzeta = aHit->getZetaID();maxZetaHit=aHit;}
        if(aHit->getZetaID() < minzeta){ minzeta = aHit->getZetaID();minZetaHit=aHit;}
      }
      FPCCDPixelHit::HitQuality_t addedQuality = aHit->getQuality();
      if( trackquality != FPCCDPixelHit::kBKGOverlap){
        if(trackquality == FPCCDPixelHit::kBKG){
          addedQuality=FPCCDPixelHit::kBKG ? trackquality=FPCCDPixelHit::kBKG : trackquality=FPCCDPixelHit::kBKGOverlap;
        }
      }
      else if(addedQuality==FPCCDPixelHit::kSignalOverlap){trackquality=FPCCDPixelHit::kSignalOverlap;       }
    }

    unsigned int xiWidth = maxxi-minxi + 1;
    unsigned int zetaWidth = maxzeta-minzeta + 1; // cluster shapes info. could be used to reject background.
    //Mori added new if statement
    if(_firstCut.isActive != true || _firstCut.ZWidth[layer] < 0 || zetaWidth < static_cast<unsigned int>(_firstCut.ZWidth[layer]) ){ 
      if(_firstCut.isActive != true || _firstCut.RPhiWidth[layer] < 0 || xiWidth < static_cast<unsigned int>(_firstCut.RPhiWidth[layer])  ){


        unsigned short int tilt=0;
        if( xiWidth==1 || zetaWidth==1 ) tilt=0;
        else{
          if(maxXiHit->getZetaID()-minXiHit->getZetaID()>0 && maxZetaHit->getXiID()-minZetaHit->getXiID()>0) tilt=1;
          else if(maxXiHit->getZetaID()-minXiHit->getZetaID()<=0 && maxZetaHit->getXiID()-minZetaHit->getXiID()<=0) tilt=2;
          else tilt=0;
        }

        double xiave = xiene/enesum;
        double zetaave = zetaene/enesum;
        double eta0 = _geodata[layer].rmin+0.5*_pixelheight +(layer%2)*(_geodata[layer].sthick-_pixelheight);

        //Mori2ndCut
        double mxi = xiave+_geodata[layer].sximin ;
        double mzeta = zetaave-_geodata[layer].hlength;
        if(_m2Cut.isActive != true || _m2Cut.zpar[layer] < 0 || ( (mxi * mzeta >= -_m2Cut.zpar[layer] || tilt != 1) && (mxi * mzeta <= _m2Cut.zpar[layer] || tilt != 2) ) ) {


          double newpos[3];
          newpos[0] = (xiave+_geodata[layer].sximin)*_geodata[layer].sinphi[ladder]
            + eta0*_geodata[layer].cosphi[ladder];
          newpos[1] =-(xiave+_geodata[layer].sximin)*_geodata[layer].cosphi[ladder]
            + eta0*_geodata[layer].sinphi[ladder];
          newpos[2] = zetaave-_geodata[layer].hlength;
          //store hit variables

          //Kamai2ndCut
          if(_k2Cut.isActive != true || int(zetaWidth) <= _k2Cut.minZWidth[layer] || float(zetaWidth) > (std::abs(newpos[2])*_pixelheight/sqrt(pow(newpos[0],2)+pow(newpos[1],2))/_pixelSizeVec[layer] - _k2Cut.bpar[layer] ) ){



            //*******************************Mori added 20121203*************************************//
            // Here finds the numbers of types of simthits that create trkhit in this scope.
            std::vector<int> typesOfOrderID;
            orderVecOfCluster.sort(); //sort() also sort order Vec 
            std::list<int>::iterator itlist;
            int tmpV = -1000;
            for(itlist = orderVecOfCluster.begin(); itlist != orderVecOfCluster.end(); itlist++){
              if(*itlist != tmpV){ 
                typesOfOrderID.push_back(*itlist);
                tmpV = *itlist;
              }
            }
            //********************************************************************************************//

            float pointResoRPhi = 0.0;
            float pointResoZ    = 0.0;
            if(_positionReso_ReadingFile_ON == true){

              short int codeRPhi = layer * 100 + xiWidth; 
              if(xiWidth <= 49){ 
                pointResoRPhi      = _resolutionMapRPhi[codeRPhi];
                if(pointResoRPhi < 0.00001){ pointResoRPhi = (sqrt((float)(xiWidth + 1.5)) - 1)*1e-3; }
              }
              else{ pointResoRPhi = (sqrt((float)(xiWidth + 1.5)) - 1)*1e-3; }

              short int codeZ    = layer * 100 + zetaWidth; 
              if(zetaWidth <= 49){ 
                pointResoZ      = _resolutionMapZ[codeZ];
                if(pointResoZ < 0.00001){ pointResoZ = (sqrt((float)(zetaWidth + 1.5)) - 1)*1e-3; }
              }
              else{ pointResoZ = (sqrt((float)(zetaWidth + 1.5)) - 1)*1e-3; }

            }
            else{
              pointResoRPhi = _pointResoRPhi;
              pointResoZ    = _pointResoZ;
            }


            TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl ;

            trkHit->setType( 100 + layer);

            cellid_encoder[ lcio::LCTrackerCellID::subdet() ] = lcio::ILDDetID::VXD ;
            cellid_encoder[ lcio::LCTrackerCellID::side()   ] = 0 ;
            cellid_encoder[ lcio::LCTrackerCellID::layer()  ] = layer ;
            cellid_encoder[ lcio::LCTrackerCellID::module() ] = ladder ;
            cellid_encoder[ lcio::LCTrackerCellID::sensor() ] = 0 ;      

            cellid_encoder.setCellID( trkHit );

            unsigned int cellid1=0;
            cellid1 = ((trackquality << 30) & 0xc0000000) |
              ((tilt << 28) & 0x30000000) |
              ((nPix << 18) & 0x0ffc0000) |
              ((zetaWidth << 9) & 0x0003fe00) |
              ((xiWidth) & 0x000001ff) ;

            trkHit->setCellID1(cellid1);
            float u_direction[2] ;
            u_direction[0] = M_PI/2.0;
            u_direction[1] = _geodata[layer].ladder_incline[ladder] ;

            float v_direction[2] ;
            v_direction[0] = 0.0 ;
            v_direction[1] = 0.0 ;

            trkHit->setU( u_direction ) ;
            trkHit->setV( v_direction ) ;


            //printf("mori check : pointResoRPhi,pointResoZ : %f,%f \n",pointResoRPhi,pointResoZ);
            trkHit->setdU( pointResoRPhi ) ;
            trkHit->setdV( pointResoZ ) ;      

            trkHit->setPosition( newpos ) ;
            trkHit->setEDep( enesum );

            //LCRelationImpl* rel = new LCRelationImpl ;
            if(_makeRelation == true){
              SimTrackerHit* p_simthit; 
              int ntypes = typesOfOrderID.size();
              for(int i = 0; i < ntypes; i++){ 
                LCRelationImpl* rel = new LCRelationImpl ;
                if(typesOfOrderID[i] >= 0){
                  p_simthit =  dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( typesOfOrderID[i] ) ) ;      
                  float nTo = static_cast<float>( ntypes ); 
                  rel->setTo( p_simthit ); 
                  rel->setFrom( trkHit ); 
                  rel->setWeight( nTo ); //Weight is defined as nFrom/nTo by Mori. ->今はnTo/nFrom
                  relCol->addElement( rel );
                }// If typesOfOrderID[i] == -1, then it is background data.
              }
            }
            trkHitVec->addElement( trkHit );
          }//end of Kamai2ndCut
        }//end of Mori2ndCut
      }//end of if statement for RPhi Cut of firstCut
    }//end of if statement for Z Cut of firstCut
  }// Repeat pixel hits in this ladder.

} 






//void FPCCDClustering::viewerOfClusterShape(FPCCDClusterVec_t &clusterVec);

// =====================================================================
void FPCCDClustering::makeClustersInALadder(int layer, FPCCDLadderHit_t &ladderHit, FPCCDClusterVec_t &clusterVec)
{
  //typedef std::pair<unsigned int, unsigned int> FPCCDHitLoc_t;
  //typedef std::map<FPCCDHitLoc_t, FPCCDPixelHit*>  FPCCDLadderHit_t;
  // namely, std::map<std::pair<unsigned int, unsigned int>, FPCCDPixelHit*>

  //typedef std::vector<FPCCDPixelHit*> FPCCDCluster_t;
  //typedef std::vector<FPCCDCluster_t*> FPCCDClusterVec_t;
  //namely, std::vector<std::vector<FPCCDPixelHit*>>


  //(2) Start pixel hit clustering 
  //    All pixels adjucent pixels are clustered, if it has a hit.  Threshold is not considered yet
  //    2013_06_07 I think thresold cut has been done before makeClustersInALadder(...).

  int   maxxi = _geodata[layer].num_xi_pixel;
  int maxzeta = _geodata[layer].num_zeta_pixel;

  FPCCDLadderHit_t::iterator ih = ladderHit.begin();
  while( ih != ladderHit.end() ) {
    FPCCDPixelHit *h=(*ih).second;
    (*ih).second=0;
    ih++;
    if( h==0 ) { continue; }
    //  std::cerr << "Found a seed pixel.  Edep=" << h->getEdep() << std::endl;
    //  Found a ssed of a cluster

    FPCCDCluster_t *cluster=new FPCCDCluster_t();
    cluster->push_back(h);
    int xi0=h->getXiID();
    int zeta0=h->getZetaID();
    typedef std::pair<int, int> Direction_t;
    std::stack<Direction_t>  direction;   
    for(int i=-1;i<=1;i++) { // Set initial search direction, neighbouring 8 pixels are searched
      for(int j=-1;j<=1;j++) {
        if( i!=0 || j!= 0 ) { direction.push(Direction_t(i,j)); }
      }
    } 
    // Now look for all directions and find pixel with hit.
    while( ! direction.empty() ) {
      Direction_t newdir=direction.top();
      int xi=xi0+newdir.first;
      int zeta=zeta0+newdir.second;
      //          std::cerr << "searching at (xi,zeta)=" << xi << "  " << zeta << std::endl;
      if( xi < 0 || xi > maxxi || zeta < 0 || zeta > maxzeta ) {
        direction.pop(); continue;
      }
      FPCCDHitLoc_t newloc(xi, zeta);
      FPCCDLadderHit_t::iterator look=ladderHit.find(newloc);//If not found, end iterator is returned.
      if( look == ladderHit.end() ) {
        direction.pop(); continue;
      }
      FPCCDPixelHit *aHit=(*look).second;
      if( aHit == 0 ) {
        direction.pop(); continue;
      }
      cluster->push_back( aHit );
      //          std::cout << " Found hit .. at xi, zeta=" << xi << " " << zeta << std::endl;
      (*look).second=0;
      for (int ix=-1;ix<=1;ix++) {
        for(int iz=-1;iz<=1;iz++) {
          if( ix != 0 || iz != 0 ) { 
            direction.push(Direction_t(newdir.first+ix, newdir.second+iz));
          }
        }
      }
    } // Complete search of all nighbours of the seed
    if( cluster->size() > 0 ) {
      //Mori added
      if(_firstCut.isActive != true || _firstCut.nPix[layer] < 0 || cluster->size() < static_cast<unsigned int>(_firstCut.nPix[layer])){
        clusterVec.push_back(cluster);  
      }
    }
  } 
}


// =====================================================================
void FPCCDClustering::check( LCEvent *  /*evt*/ ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


// =====================================================================
void FPCCDClustering::end(){ 

#if ROOT_DEBUG
  _rootf->Write();
#endif

  streamlog_out(MESSAGE) << " end()  " << name() 
    << " processed " << _nEvt << " events in " << _nRun << " runs "
    << std::endl ;

}


// =================================================================
// =================================================================
void FPCCDClustering::EnergyDigitizer( FPCCDPixelHit* aHit ){

  int nEle = (int)(aHit->getEdep()*1e+6 * _electronsPerKeV) ;
  int nStep = 1 << _nbitsForEdep ;
  int count = nEle/ _electronsPerStep ; //  int devided by int equales int. 

  if(count > nStep ) count = nStep ;

  aHit->setEdep( (count * _electronsPerStep)*1e-6 / _electronsPerKeV) ;
  return ;
}














void FPCCDClustering::calcTrackParameterOfMCP(MCParticle* pmcp, double* par){
  //KalTest Reference Manual page 17 

  double bz = Global::GEAR->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  const double* p;
  p = pmcp->getMomentum();
  const double pi = 3.14159265359;
  const double pt = sqrt(p[0]*p[0]+p[1]*p[1]) ;
  double radius = pt / (2.99792458E-4*bz) ; // for r in mm, p in GeV and Bz in Tesla
  double omega = pmcp->getCharge()/radius ;
  double tanL = p[2]/pt ;
  double oldphi0 = atan2(p[1],p[0]);
  while( oldphi0 <= -pi ){  oldphi0 += 2. * pi ; } 
  while( oldphi0 >   pi ){  oldphi0 -= 2. * pi ; } 

  double oldsinphi0 = std::sin(oldphi0);
  double oldcosphi0 = std::cos(oldphi0);
  double centx = pmcp->getVertex()[0] + 1/omega * oldcosphi0;
  double centy = pmcp->getVertex()[1] + 1/omega * oldsinphi0;
  double newphi0;
  if(omega > 0){ newphi0 = std::atan2(centy,centx); }   
  else{ newphi0 = std::atan2(-centy,-centx); }   
  double sinphi0 = std::sin(newphi0);
  double cosphi0 = std::cos(newphi0);
  double d0 = centx * cosphi0 + centy * sinphi0 - 1.0/omega ;
  double z0 = pmcp->getVertex()[2] - 1.0/omega *(newphi0 - oldphi0)*tanL;
  /*
     std::cout << "Test Pivot Transformation " << std::endl;
     std::cout << "d0 : " << d0 << std::endl;
     std::cout << "z0 : " << z0 << std::endl;
     std::cout << "omega : " << omega << std::endl;
     std::cout << "phi0 : " << newphi0 << std::endl;
     std::cout << "tanL : " << tanL << std::endl;
   */
  par[0] = d0;
  par[1] = z0;
  par[2] = omega;
  par[3] = newphi0;
  par[4] = tanL;
}
