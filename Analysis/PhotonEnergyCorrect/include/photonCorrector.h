#ifndef photonCorrector_h
#define photonCorrector_h 1

#include "EVENT/ReconstructedParticle.h"
#include <cassert>

class photonCorrector {
 public:
  photonCorrector() {
    _barrelendcap_costhlimit    = 0.;
    _energyLin_const            = 0.;
    _energyLin_logen            = 0.;
    _phiBarrelCorr_pos_const    = 0.;
    _phiBarrelCorr_pos_logen    = 0.;
    _phiBarrelCorr_depth        = 0.;
    _phiBarrelCorr_width1       = 0.;
    _phiBarrelCorr_width2       = 0.;
    _costhCorr_gaus1_norm_const = 0.;
    _costhCorr_gaus1_norm_logen = 0.;
    _costhCorr_gaus1_mean       = 0.;
    _costhCorr_gaus1_sigm       = 0.;
    _costhCorr_gaus2_norm_const = 0.;
    _costhCorr_gaus2_norm_logen = 0.;
    _costhCorr_gaus2_mean       = 0.;
    _costhCorr_gaus2_sigm       = 0.;
    _costhCorr_gaus3_norm       = 0.;
    _costhCorr_gaus3_mean       = 0.;
    _costhCorr_gaus3_sigm       = 0.;
    _costhCorr_endcap_scale     = 0.;
    _endcap_gaus1_norm          = 0.;
    _endcap_gaus1_mean          = 0.;
    _endcap_gaus1_sigm          = 0.;
    _endcap_gaus2_norm          = 0.;
    _endcap_gaus2_mean          = 0.;
    _endcap_gaus2_sigm          = 0.;
    _assumed_boxsize            = 0.;
    _assumed_endZ               = 0.;

    _phiBias_barrel_p0_1 = 0.;
    _phiBias_barrel_p0_2 = 0.;
    _phiBias_barrel_p0_3 = 0.;
    _phiBias_barrel_p0_4 = 0.;
    _phiBias_barrel_p1_1 = 0.;
    _phiBias_barrel_p1_2 = 0.;
    _phiBias_barrel_p2_1 = 0.;
    _phiBias_barrel_p2_2 = 0.;
    _phiBias_barrel_p3_1 = 0.;
    _phiBias_barrel_p4_1 = 0.;
    _phiBias_barrel_p5_1 = 0.;
    _phiBias_barrel_p5_2 = 0.;
    _phiBias_barrel_p5_3 = 0.;
    _phiBias_barrel_p6_1 = 0.;
    _phiBias_barrel_p6_2 = 0.;
    _phiBias_barrel_p7_1 = 0.;
    _phiBias_barrel_p7_2 = 0.;
    _phiBias_barrel_p7_3 = 0.;
    _thetaBias_barrel_p0_1 = 0;
    _thetaBias_barrel_p0_2 = 0;
    _thetaBias_barrel_p1_1 = 0;
    _thetaBias_barrel_p1_2 = 0;
    _thetaBias_endcap_p0_1 = 0.;
    _thetaBias_endcap_p0_2 = 0.;
    _thetaBias_endcap_p1_1 = 0.;
    _thetaBias_endcap_p1_2 = 0.;
    _thetaBias_endcap_p2_1 = 0.;
    _thetaBias_endcap_p2_2 = 0.;

  }
  ~photonCorrector() {}


  // the main energy corrector
  float photonEnergyCorrection( EVENT::ReconstructedParticle* rp );
  void  photonDirectionCorrection( EVENT::ReconstructedParticle* rp , float& cor_theta, float& cor_phi );

  // parameter setters
  void set_barrelendcap_limit         ( float x ) { _barrelendcap_costhlimit    = x; }
  void set_energyLin_const            ( float x ) { _energyLin_const            = x; }
  void set_energyLin_logen            ( float x ) { _energyLin_logen            = x; }
  void set_phiBarrelCorr_pos_const    ( float x ) { _phiBarrelCorr_pos_const    = x; }
  void set_phiBarrelCorr_pos_logen    ( float x ) { _phiBarrelCorr_pos_logen    = x; }
  void set_phiBarrelCorr_depth        ( float x ) { _phiBarrelCorr_depth        = x; }
  void set_phiBarrelCorr_width1       ( float x ) { _phiBarrelCorr_width1       = x; }
  void set_phiBarrelCorr_width2       ( float x ) { _phiBarrelCorr_width2       = x; }
  void set_costhCorr_gaus1_norm_const ( float x ) { _costhCorr_gaus1_norm_const = x; }
  void set_costhCorr_gaus1_norm_logen ( float x ) { _costhCorr_gaus1_norm_logen = x; }
  void set_costhCorr_gaus1_mean       ( float x ) { _costhCorr_gaus1_mean       = x; }
  void set_costhCorr_gaus1_sigm       ( float x ) { _costhCorr_gaus1_sigm       = x; }
  void set_costhCorr_gaus2_norm_const ( float x ) { _costhCorr_gaus2_norm_const = x; }
  void set_costhCorr_gaus2_norm_logen ( float x ) { _costhCorr_gaus2_norm_logen = x; }
  void set_costhCorr_gaus2_mean       ( float x ) { _costhCorr_gaus2_mean       = x; }
  void set_costhCorr_gaus2_sigm       ( float x ) { _costhCorr_gaus2_sigm       = x; }
  void set_costhCorr_gaus3_norm       ( float x ) { _costhCorr_gaus3_norm       = x; }
  void set_costhCorr_gaus3_mean       ( float x ) { _costhCorr_gaus3_mean       = x; }
  void set_costhCorr_gaus3_sigm       ( float x ) { _costhCorr_gaus3_sigm       = x; }
  void set_costhCorr_endcap_scale     ( float x ) { _costhCorr_endcap_scale     = x; }
  void set_endcap_gaus1_norm          ( float x ) { _endcap_gaus1_norm          = x; }
  void set_endcap_gaus1_mean          ( float x ) { _endcap_gaus1_mean          = x; }
  void set_endcap_gaus1_sigm          ( float x ) { _endcap_gaus1_sigm          = x; }
  void set_endcap_gaus2_norm          ( float x ) { _endcap_gaus2_norm          = x; }
  void set_endcap_gaus2_mean          ( float x ) { _endcap_gaus2_mean          = x; }
  void set_endcap_gaus2_sigm          ( float x ) { _endcap_gaus2_sigm          = x; }
  void set_assumed_boxsize            ( float x ) { _assumed_boxsize            = x; }
  void set_assumed_endZ               ( float x ) { _assumed_endZ               = x; }


  void set_phiBias_barrel_p0_1    (float x) {_phiBias_barrel_p0_1   = x;}
  void set_phiBias_barrel_p0_2    (float x) {_phiBias_barrel_p0_2   = x;}
  void set_phiBias_barrel_p0_3    (float x) {_phiBias_barrel_p0_3   = x;}
  void set_phiBias_barrel_p0_4    (float x) {_phiBias_barrel_p0_4   = x;}
  void set_phiBias_barrel_p1_1    (float x) {_phiBias_barrel_p1_1   = x;}
  void set_phiBias_barrel_p1_2    (float x) {_phiBias_barrel_p1_2   = x;}
  void set_phiBias_barrel_p2_1    (float x) {_phiBias_barrel_p2_1   = x;}
  void set_phiBias_barrel_p2_2    (float x) {_phiBias_barrel_p2_2   = x;}
  void set_phiBias_barrel_p3_1    (float x) {_phiBias_barrel_p3_1   = x;}
  void set_phiBias_barrel_p4_1    (float x) {_phiBias_barrel_p4_1   = x;}
  void set_phiBias_barrel_p5_1    (float x) {_phiBias_barrel_p5_1   = x;}
  void set_phiBias_barrel_p5_2    (float x) {_phiBias_barrel_p5_2   = x;}
  void set_phiBias_barrel_p5_3    (float x) {_phiBias_barrel_p5_3   = x;}
  void set_phiBias_barrel_p6_1    (float x) {_phiBias_barrel_p6_1   = x;}
  void set_phiBias_barrel_p6_2    (float x) {_phiBias_barrel_p6_2   = x;}
  void set_phiBias_barrel_p7_1    (float x) {_phiBias_barrel_p7_1   = x;}
  void set_phiBias_barrel_p7_2    (float x) {_phiBias_barrel_p7_2   = x;}
  void set_phiBias_barrel_p7_3    (float x) {_phiBias_barrel_p7_3   = x;}
  void set_thetaBias_barrel_p0_1  (float x) {_thetaBias_barrel_p0_1 = x;}
  void set_thetaBias_barrel_p0_2  (float x) {_thetaBias_barrel_p0_2 = x;}
  void set_thetaBias_barrel_p1_1  (float x) {_thetaBias_barrel_p1_1 = x;}
  void set_thetaBias_barrel_p1_2  (float x) {_thetaBias_barrel_p1_2 = x;}
  void set_thetaBias_endcap_p0_1  (float x) {_thetaBias_endcap_p0_1 = x;}
  void set_thetaBias_endcap_p0_2  (float x) {_thetaBias_endcap_p0_2 = x;}
  void set_thetaBias_endcap_p1_1  (float x) {_thetaBias_endcap_p1_1 = x;}
  void set_thetaBias_endcap_p1_2  (float x) {_thetaBias_endcap_p1_2 = x;}
  void set_thetaBias_endcap_p2_1  (float x) {_thetaBias_endcap_p2_1 = x;}
  void set_thetaBias_endcap_p2_2  (float x) {_thetaBias_endcap_p2_2 = x;}


  float get_barrelendcap_limit         ( ) { return _barrelendcap_costhlimit    ; }
  float get_energyLin_const            ( ) { return _energyLin_const            ; }
  float get_energyLin_logen            ( ) { return _energyLin_logen            ; }
  float get_phiBarrelCorr_pos_const    ( ) { return _phiBarrelCorr_pos_const    ; }
  float get_phiBarrelCorr_pos_logen    ( ) { return _phiBarrelCorr_pos_logen    ; }
  float get_phiBarrelCorr_depth        ( ) { return _phiBarrelCorr_depth        ; }
  float get_phiBarrelCorr_width1       ( ) { return _phiBarrelCorr_width1       ; }
  float get_phiBarrelCorr_width2       ( ) { return _phiBarrelCorr_width2       ; }
  float get_costhCorr_gaus1_norm_const ( ) { return _costhCorr_gaus1_norm_const ; }
  float get_costhCorr_gaus1_norm_logen ( ) { return _costhCorr_gaus1_norm_logen ; }
  float get_costhCorr_gaus1_mean       ( ) { return _costhCorr_gaus1_mean       ; }
  float get_costhCorr_gaus1_sigm       ( ) { return _costhCorr_gaus1_sigm       ; }
  float get_costhCorr_gaus2_norm_const ( ) { return _costhCorr_gaus2_norm_const ; }
  float get_costhCorr_gaus2_norm_logen ( ) { return _costhCorr_gaus2_norm_logen ; }
  float get_costhCorr_gaus2_mean       ( ) { return _costhCorr_gaus2_mean       ; }
  float get_costhCorr_gaus2_sigm       ( ) { return _costhCorr_gaus2_sigm       ; }
  float get_costhCorr_gaus3_norm       ( ) { return _costhCorr_gaus3_norm       ; }
  float get_costhCorr_gaus3_mean       ( ) { return _costhCorr_gaus3_mean       ; }
  float get_costhCorr_gaus3_sigm       ( ) { return _costhCorr_gaus3_sigm       ; }
  float get_costhCorr_endcap_scale     ( ) { return _costhCorr_endcap_scale     ; }
  float get_endcap_gaus1_norm          ( ) { return _endcap_gaus1_norm          ; }
  float get_endcap_gaus1_mean          ( ) { return _endcap_gaus1_mean          ; }
  float get_endcap_gaus1_sigm          ( ) { return _endcap_gaus1_sigm          ; }
  float get_endcap_gaus2_norm          ( ) { return _endcap_gaus2_norm          ; }
  float get_endcap_gaus2_mean          ( ) { return _endcap_gaus2_mean          ; }
  float get_endcap_gaus2_sigm          ( ) { return _endcap_gaus2_sigm          ; }
  float get_assumed_boxsize            ( ) { return _assumed_boxsize            ; }
  float get_assumed_endZ               ( ) { return _assumed_endZ               ; }


  float get_phiBias_barrel_p0_1    () {return _phiBias_barrel_p0_1  ;}
  float get_phiBias_barrel_p0_2    () {return _phiBias_barrel_p0_2  ;}
  float get_phiBias_barrel_p0_3    () {return _phiBias_barrel_p0_3  ;}
  float get_phiBias_barrel_p0_4    () {return _phiBias_barrel_p0_4  ;}
  float get_phiBias_barrel_p1_1    () {return _phiBias_barrel_p1_1  ;}
  float get_phiBias_barrel_p1_2    () {return _phiBias_barrel_p1_2  ;}
  float get_phiBias_barrel_p2_1    () {return _phiBias_barrel_p2_1  ;}
  float get_phiBias_barrel_p2_2    () {return _phiBias_barrel_p2_2  ;}
  float get_phiBias_barrel_p3_1    () {return _phiBias_barrel_p3_1  ;}
  float get_phiBias_barrel_p4_1    () {return _phiBias_barrel_p4_1  ;}
  float get_phiBias_barrel_p5_1    () {return _phiBias_barrel_p5_1  ;}
  float get_phiBias_barrel_p5_2    () {return _phiBias_barrel_p5_2  ;}
  float get_phiBias_barrel_p5_3    () {return _phiBias_barrel_p5_3  ;}
  float get_phiBias_barrel_p6_1    () {return _phiBias_barrel_p6_1  ;}
  float get_phiBias_barrel_p6_2    () {return _phiBias_barrel_p6_2  ;}
  float get_phiBias_barrel_p7_1    () {return _phiBias_barrel_p7_1  ;}
  float get_phiBias_barrel_p7_2    () {return _phiBias_barrel_p7_2  ;}
  float get_phiBias_barrel_p7_3    () {return _phiBias_barrel_p7_3  ;}
  float get_thetaBias_barrel_p0_1  () {return _thetaBias_barrel_p0_1;}
  float get_thetaBias_barrel_p0_2  () {return _thetaBias_barrel_p0_2;}
  float get_thetaBias_barrel_p1_1  () {return _thetaBias_barrel_p1_1;}
  float get_thetaBias_barrel_p1_2  () {return _thetaBias_barrel_p1_2;}
  float get_thetaBias_endcap_p0_1  () {return _thetaBias_endcap_p0_1;}
  float get_thetaBias_endcap_p0_2  () {return _thetaBias_endcap_p0_2;}
  float get_thetaBias_endcap_p1_1  () {return _thetaBias_endcap_p1_1;}
  float get_thetaBias_endcap_p1_2  () {return _thetaBias_endcap_p1_2;}
  float get_thetaBias_endcap_p2_1  () {return _thetaBias_endcap_p2_1;}
  float get_thetaBias_endcap_p2_2  () {return _thetaBias_endcap_p2_2;}


  void printParams();

 private:

  float energyLinearise( float en );
  float gapCompensatedEnergy( EVENT::ReconstructedParticle* rp );
  float gapCompensate_barrelPhi( float en, float phi );
  float gapCompensate_theta( float en, float costh );
  float gapCompensate_endcap( float xAcross );
  float getDistanceAcrossEndcapQuadrant( float costh, float phi );
  float getBarrelFoldedPhi( float phi);

  float getBarrelCorrectedPhi( float en, float phi );
  float getCorrectedPhi      ( float en, float costh, float phi );

  float getBarrelCorrectedTheta( float en, float theta );
  float getEndcapCorrectedTheta( float en, float theta );
  float getCorrectedTheta      ( float en, float theta );


  // -------------------

  float _barrelendcap_costhlimit{};

  float _energyLin_const{};
  float _energyLin_logen{};

  float _phiBarrelCorr_pos_const{};
  float _phiBarrelCorr_pos_logen{};
  float _phiBarrelCorr_depth{};
  float _phiBarrelCorr_width1{};
  float _phiBarrelCorr_width2{};

  float _costhCorr_gaus1_norm_const{};
  float _costhCorr_gaus1_norm_logen{};
  float _costhCorr_gaus1_mean{};
  float _costhCorr_gaus1_sigm{};

  float _costhCorr_gaus2_norm_const{};
  float _costhCorr_gaus2_norm_logen{};
  float _costhCorr_gaus2_mean{};
  float _costhCorr_gaus2_sigm{};

  float _costhCorr_gaus3_norm{};
  float _costhCorr_gaus3_mean{};
  float _costhCorr_gaus3_sigm{};
  float _costhCorr_endcap_scale{};

  float _endcap_gaus1_norm{};
  float _endcap_gaus1_mean{};
  float _endcap_gaus1_sigm{};
  float _endcap_gaus2_norm{};
  float _endcap_gaus2_mean{};
  float _endcap_gaus2_sigm{};

  float _phiBias_barrel_p0_1  {};
  float _phiBias_barrel_p0_2  {};
  float _phiBias_barrel_p0_3  {};
  float _phiBias_barrel_p0_4  {};
  float _phiBias_barrel_p1_1  {};
  float _phiBias_barrel_p1_2  {};
  float _phiBias_barrel_p2_1  {};
  float _phiBias_barrel_p2_2  {};
  float _phiBias_barrel_p3_1  {};
  float _phiBias_barrel_p4_1  {};
  float _phiBias_barrel_p5_1  {};
  float _phiBias_barrel_p5_2  {};
  float _phiBias_barrel_p5_3  {};
  float _phiBias_barrel_p6_1  {};
  float _phiBias_barrel_p6_2  {};
  float _phiBias_barrel_p7_1  {};
  float _phiBias_barrel_p7_2  {};
  float _phiBias_barrel_p7_3  {};
  float _thetaBias_barrel_p0_1{};
  float _thetaBias_barrel_p0_2{};
  float _thetaBias_barrel_p1_1{};
  float _thetaBias_barrel_p1_2{};
  float _thetaBias_endcap_p0_1{};
  float _thetaBias_endcap_p0_2{};
  float _thetaBias_endcap_p1_1{};
  float _thetaBias_endcap_p1_2{};
  float _thetaBias_endcap_p2_1{};
  float _thetaBias_endcap_p2_2{};

  float _assumed_boxsize{};
  float _assumed_endZ{};

};
#endif
