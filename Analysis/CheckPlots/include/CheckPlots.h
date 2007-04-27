#ifndef CheckPlots_h
#define CheckPlots_h 1


#include <iostream>
#include <string>
#include <vector>

#include "marlin/Processor.h"
#include <marlin/Exceptions.h>
#include <marlin/Global.h>

#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/Track.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/LCRelationNavigator.h>

#include <gear/GEAR.h>
#include <gear/TPCParameters.h>

#include <MarlinUtil.h>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/ICloud2D.h>
#endif


using namespace lcio ;
using namespace marlin ;


/**
 *    This processor provides check plots. The plots are arranged in the following different categories:
 *    <br>
 *    1. MC particle related plots <br>
 *    2. Plots related to the simulated hits in tracking and calorimeter devices <br>
 *    3. Plots related to the hits in tracking and calorimeter devices <br>
 *    4. Plots related to tracks <br>
 *    5. Plots related to reconstructed particles <br>
 *    6. more to come <br>
 *    <br>
 *    Steering parameters: <br>
 *    FillMCGen: toggles the check plots for the MC particles with generator status == 1 ( 0 or 1 ) <br>
 *    FillMCSim : toggles the check plots for the MC particles with generator status != 1 ( 0 or 1 ) <br>
 *    FillSimCalo : produces check plots for the simulated calorimeter hits <br>
 *    SimECut : cut on the simulated energy in the cell. Only energies above this cut are taken into account. <br>
 *    FillCalo : produces check plots for the calorimeter hits <br>
 *    ECut : cut on the energy in the cell. Only energies above this cut are taken into account. <br>
 *    ThetaCut : cut in theta to assign particles which are lost in the beam pipe
 *    FillTracks : toggles the check plots for the tracks
 *    ColNameTracks : name of the LCCollection of tracks
 *    ColNameRelationTrackToMCP : name of the LCRelation collection connecting the tracks and the corresponding MC particle
 *    FillReconstructedParticles : toggles the check plots for the reconstructed paricles
 *    ColNameReconstructedParticles : name of the LCCollection of reconstructed paricles
 *
 *    @author O. Wendt (DESY)
 *    @version $Id: CheckPlots.h,v 1.3 2007-04-27 13:30:46 owendt Exp $
 *
 */
class CheckPlots : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CheckPlots ; }
    
  CheckPlots() ;
  
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  
  
 private:

  int _nRun;
  int _nEvt;

  float _bField;

  int _fillMCGen;
  float _thetaCut;
  int _fillMCSim;

  int _fillSimCaloHit;
  float _simECut;
  
  int _fillCaloHit;
  float _ECut;

  int _fillTracks;
  std::string _colNameTracks;
  std::string _colNameRelationTrackToMCP;

  int _fillReconstructedParticles;
  std::string _colNameReconstructedParticles;

  int _fillComparisonMCReco;
  int _nMC;
  int _nMCCh;
  int _nMCN;
  int _nReco;
  int _nRecoCh;
  int _nRecoN;

  double _energyMC;
  double _energyMCCh;
  double _energyMCN;
  double _energyReco;
  double _energyRecoCh;
  double _energyRecoN;



  #ifdef MARLIN_USE_AIDA
    
  // MCPs with generator status != 1
  // numbers per event
  AIDA::ICloud1D* _cMCNumberSim;
  AIDA::ICloud1D* _cMCEnergySumSim;
  
  AIDA::ICloud1D* _cMCNumberElectronsSim;
  AIDA::ICloud1D* _cMCNumberMuonsSim;
  AIDA::ICloud1D* _cMCNumberTausSim;
  AIDA::ICloud1D* _cMCNumberNusSim;
  
  AIDA::ICloud1D* _cMCNumberPiChSim;
  AIDA::ICloud1D* _cMCNumberKChSim;
  AIDA::ICloud1D* _cMCNumberProtonsSim;
  AIDA::ICloud1D* _cMCNumberPi0Sim;
  AIDA::ICloud1D* _cMCNumberK0lSim;
  AIDA::ICloud1D* _cMCNumberK0sSim;
  AIDA::ICloud1D* _cMCNumberNeutronsSim;
  AIDA::ICloud1D* _cMCNumberGammasSim;
  AIDA::ICloud1D* _cMCNumberLambda0sSim;
  AIDA::ICloud1D* _cMCNumberSigma0sSim;
  AIDA::ICloud1D* _cMCNumberXi0sSim;  

  AIDA::ICloud1D* _cMCNumberRemainingSim;
  
    

    
  // MCPs with generator status == 1
  // numbers per event
  AIDA::ICloud1D* _cMCNumberGen;
  AIDA::ICloud1D* _cMCEnergySumGen;
  
  AIDA::ICloud1D* _cMCNumberHChGen;
  AIDA::ICloud1D* _cMCNumberH0Gen;
  AIDA::ICloud1D* _cMCNumberGGen; // all gammas, i.e. gammas and pi0s
  AIDA::ICloud1D* _cMCFractionHChGen;
  AIDA::ICloud1D* _cMCFractionH0Gen;
  AIDA::ICloud1D* _cMCFractionGGen;

  AIDA::ICloud1D* _cMCNumberElectronsGen;
  AIDA::ICloud1D* _cMCNumberMuonsGen;
  AIDA::ICloud1D* _cMCNumberTausGen;
  AIDA::ICloud1D* _cMCNumberNusGen;
  
  AIDA::ICloud1D* _cMCNumberPiChGen;
  AIDA::ICloud1D* _cMCNumberKChGen;
  AIDA::ICloud1D* _cMCNumberProtonsGen;
  AIDA::ICloud1D* _cMCNumberPi0Gen;
  AIDA::ICloud1D* _cMCNumberK0lGen;
  AIDA::ICloud1D* _cMCNumberK0sGen;
  AIDA::ICloud1D* _cMCNumberNeutronsGen;
  AIDA::ICloud1D* _cMCNumberGammasGen;
  AIDA::ICloud1D* _cMCNumberLambda0sGen;
  AIDA::ICloud1D* _cMCNumberSigma0sGen;
  AIDA::ICloud1D* _cMCNumberXi0sGen;  

  AIDA::ICloud1D* _cMCNumberLostInBeamPipe;
  AIDA::ICloud1D* _cMCNumberRemainingGen;
    
    
    
  // MCPs with generator status != 1
  // numbers per single particle
  AIDA::ICloud1D* _cMCEnergySim;

  AIDA::ICloud1D* _cMCEnergyElectronsSim;
  AIDA::ICloud1D* _cMCEnergyMuonsSim;
  AIDA::ICloud1D* _cMCEnergyTausSim;
  AIDA::ICloud1D* _cMCEnergyNusSim;
  
  AIDA::ICloud1D* _cMCEnergyPiChSim;
  AIDA::ICloud1D* _cMCEnergyKChSim;
  AIDA::ICloud1D* _cMCEnergyProtonsSim;
  AIDA::ICloud1D* _cMCEnergyPi0Sim;
  AIDA::ICloud1D* _cMCEnergyK0lSim;
  AIDA::ICloud1D* _cMCEnergyK0sSim;
  AIDA::ICloud1D* _cMCEnergyNeutronsSim;
  AIDA::ICloud1D* _cMCEnergyGammasSim;
  AIDA::ICloud1D* _cMCEnergyLambda0sSim;
  AIDA::ICloud1D* _cMCEnergySigma0sSim;
  AIDA::ICloud1D* _cMCEnergyXi0sSim;  

  AIDA::ICloud1D* _cMCEnergyRemainingSim;
  


  // MCPs with generator status == 1
  // numbers per single particle
  AIDA::ICloud1D* _cMCEnergyGen;

  AIDA::ICloud1D* _cMCEnergyHChGen;
  AIDA::ICloud1D* _cMCEnergyH0Gen;
  AIDA::ICloud1D* _cMCEnergyGGen; // all gammas, i.e. gammas and pi0s
  AIDA::ICloud1D* _cMCEnergyFractionHChGen;
  AIDA::ICloud1D* _cMCEnergyFractionH0Gen;
  AIDA::ICloud1D* _cMCEnergyFractionGGen;

  AIDA::ICloud1D* _cMCEnergyElectronsGen;
  AIDA::ICloud1D* _cMCEnergyMuonsGen;
  AIDA::ICloud1D* _cMCEnergyTausGen;
  AIDA::ICloud1D* _cMCEnergyNusGen;
  
  AIDA::ICloud1D* _cMCEnergyPiChGen;
  AIDA::ICloud1D* _cMCEnergyKChGen;
  AIDA::ICloud1D* _cMCEnergyProtonsGen;
  AIDA::ICloud1D* _cMCEnergyPi0Gen;
  AIDA::ICloud1D* _cMCEnergyK0lGen;
  AIDA::ICloud1D* _cMCEnergyK0sGen;
  AIDA::ICloud1D* _cMCEnergyNeutronsGen;
  AIDA::ICloud1D* _cMCEnergyGammasGen;
  AIDA::ICloud1D* _cMCEnergyLambda0sGen;
  AIDA::ICloud1D* _cMCEnergySigma0sGen;
  AIDA::ICloud1D* _cMCEnergyXi0sGen;  

  AIDA::ICloud1D* _cMCEnergyLostInBeamPipe;
  AIDA::ICloud1D* _cMCEnergyRemainingGen;



  // SimCalorimeterHit spectra      
  AIDA::ICloud1D* _cNumberSimCaloHits;
  AIDA::ICloud1D* _cEnergySimCaloHitsSum;  
  AIDA::ICloud1D* _cEnergySimCaloHits;
      
  // CalorimeterHit spectra      
  AIDA::ICloud1D* _cNumberCaloHits;
  AIDA::ICloud1D* _cEnergyCaloHitsSum;  
  AIDA::ICloud1D* _cEnergyCaloHits;


  // clouds corresponding to the tracks
  AIDA::ICloud1D* _cNumberTracks;
  AIDA::ICloud1D* _cNumberTrackerHitsPerTrack;
  AIDA::ICloud1D* _cMomentumTracks;
  AIDA::ICloud1D* _cNumberMCParticlesPerTrack;


  // clouds corresponding to the reconstructed particles
  AIDA::ICloud1D* _cNumberReconstructedParticles;
  AIDA::ICloud1D* _cEnergyReconstructedParticles;
  AIDA::ICloud1D* _cEnergySumReconstructedParticles;


  // clouds for comparison of MC tree and reconstructed particles
  AIDA::ICloud2D* _cNumberMCvsNumberReco;
  AIDA::ICloud2D* _cNumberMCChvsNumberRecoCh;
  AIDA::ICloud2D* _cNumberMCNvsNumberRecoN;
  /*
  AIDA::ICloud2D* _cNumberMCH0vsNumberRecoH0;
  AIDA::ICloud2D* _cNumberMCGammavsNumberRecoGamma;
  */

  AIDA::ICloud2D* _cEnergyMCvsEnergyReco;
  AIDA::ICloud2D* _cEnergyMCChvsEnergyRecoCh;
  AIDA::ICloud2D* _cEnergyMCNvsEnergyRecoN;
     
  #endif



  void createClouds();

  void fillMCGenCheckPlots(LCEvent * evt);
  void fillMCSimCheckPlots(LCEvent * evt);

  //  void fillSimTrackerHitCheckPlots(LCEvent * evt);
  void fillSimCaloHitCheckPlots(LCEvent * evt);

  //  void fillTrackerHitCheckPlots(LCEvent * evt);
  void fillCaloHitCheckPlots(LCEvent * evt);

  void fillTrackCheckPlots(LCEvent * evt);

  void fillReconstructedParticlesCheckPlots(LCEvent * evt);

  void fillComparisonMCRecoPlots();
  
} ;

#endif



