#ifndef photonCorrector_h
#define photonCorrector_h 1

#include "EVENT/ReconstructedParticle.h"
#include <cassert>

class photonCorrector {
 public:
  photonCorrector() {_initialized=false;}
  ~photonCorrector() {}

  enum { set_std, set_noInterMod };

  void setDefaultValues(int defSet) {
    if ( defSet == set_std ) {
      setDefaultValues_interModBruteForceCorr();
    } else if ( defSet == set_noInterMod ) {
      setDefaultValues_no_interModBruteForceCorr();
    } else {
      std::cout << "photonCorrector: unknown default param set! giving up!" << std::endl;
      assert(0);
    }
    return;
  }

  void setDefaultValues_interModBruteForceCorr() {
    // this one is for the case in which inter-module gaps are corrected by "brute force" at hit reco level
    // determined for ILD_l5_o1_v02 model

    _barrel_limit=0.75;
    _endcap_limit=0.8;

    _energyLin_const = 9.87661e-01;
    _energyLin_logen = 1.40676e-02;

    _phiBarrelCorr_pos_const = 0.411540;
    _phiBarrelCorr_pos_logen = 0.0149989;
    _phiBarrelCorr_depth     = -0.0923036;
    _phiBarrelCorr_width1    =  0.014018;
    _phiBarrelCorr_width2    = 0.039961;

    _costhCorr_gaus1_norm_const = -7.97130e-02;
    _costhCorr_gaus1_norm_logen = 5.85799e-02;
    _costhCorr_gaus1_mean = 0.235;
    _costhCorr_gaus1_sigm = 0.01;

    _costhCorr_gaus2_norm_const = -3.52624e-02;
    _costhCorr_gaus2_norm_logen = 3.07763e-02;
    _costhCorr_gaus2_mean = 0.588;
    _costhCorr_gaus2_sigm = 0.009;

    _costhCorr_gaus3_norm = -0.0422968;
    _costhCorr_gaus3_mean = 0.774;
    _costhCorr_gaus3_sigm = 0.009;

    _costhCorr_endcap_scale = 1.002;

    _endcap_gaus1_norm = -0.025;
    _endcap_gaus1_mean = 855.  ;
    _endcap_gaus1_sigm = 23.   ;

    _endcap_gaus2_norm = -0.07 ;
    _endcap_gaus2_mean = 1489. ;
    _endcap_gaus2_sigm = 18.   ;

    // actual values not important, but should use same defs when determining the correction and applying it
    _assumed_boxsize=400; // size of inner endcap box (ECAL plug/ring)
    _assumed_endZ = 2411.; // start of endcap in z
  }

  void setDefaultValues_no_interModBruteForceCorr() {

    // this one is for the case in which inter-module gaps are not corrected by "brute force" at hit reco level
    // determined for ILD_l5_o1_v02 model

    _barrel_limit=0.75;
    _endcap_limit=0.8;

    _energyLin_const = 9.870e-01;
    _energyLin_logen = 1.426e-02;

    _phiBarrelCorr_pos_const = 0.412249;
    _phiBarrelCorr_pos_logen = 0.0142289;
    _phiBarrelCorr_depth     = -0.0933687;
    _phiBarrelCorr_width1    =  0.01345;
    _phiBarrelCorr_width2    = 0.0408156;

    _costhCorr_gaus1_norm_const = -0.0900;
    _costhCorr_gaus1_norm_logen = 0;
    _costhCorr_gaus1_mean = 0.235;
    _costhCorr_gaus1_sigm = 0.007256;

    _costhCorr_gaus2_norm_const = -0.0369648;
    _costhCorr_gaus2_norm_logen = 0;
    _costhCorr_gaus2_mean = 0.588;
    _costhCorr_gaus2_sigm = 0.0121604;

    _costhCorr_gaus3_norm = -0.0422968;
    _costhCorr_gaus3_mean = 0.774;
    _costhCorr_gaus3_sigm = 0.009;

    _costhCorr_endcap_scale = 1.002;

    _endcap_gaus1_norm = -0.025;
    _endcap_gaus1_mean = 855.  ;
    _endcap_gaus1_sigm = 23.   ;

    _endcap_gaus2_norm = -0.07 ;
    _endcap_gaus2_mean = 1489. ;
    _endcap_gaus2_sigm = 18.   ;

    // actual values not important, but should use same defs when determining the correction and applying it
    _assumed_boxsize=400; // size of inner endcap box (ECAL plug/ring)
    _assumed_endZ = 2411.; // start of endcap in z
  }


  // the main energy corrector
  //  float correctEnergy( EVENT::ReconstructedParticle* rp );
  float photonEnergyCorrection( EVENT::ReconstructedParticle* rp );


  // parameter setters
  void set_barrel_limit               ( float x ) { _barrel_limit               = x; }
  void set_endcap_limit               ( float x ) { _endcap_limit               = x; }
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


 private:

  bool _initialized{};

  float energyLinearise( float en );
  float gapCompensatedEnergy( EVENT::ReconstructedParticle* rp );
  float gapCompensate_barrelPhi( float en, float phi );
  float gapCompensate_theta( float en, float costh );
  float gapCompensate_endcap( float xAcross );
  float getDistanceAcrossEndcapQuadrant( float costh, float phi );
  float getBarrelFoldedPhi( float phi);

  // -------------------

  float _barrel_limit{};
  float _endcap_limit{};

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

  float _assumed_boxsize{};
  float _assumed_endZ{};

};
#endif
