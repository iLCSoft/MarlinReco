#include "photonCorrectionProcessor.h"
#include <iostream>
#include <marlin/Global.h>
#include "lcio.h"

#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "DD4hep/DD4hepUnits.h" 
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"
#include "DDRec/DetectorData.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/AIDA.h>
#endif


using namespace lcio ;
using namespace marlin ;
using std::cout;
using std::endl;

photonCorrectionProcessor aphotonCorrectionProcessor ;


photonCorrectionProcessor::photonCorrectionProcessor() : Processor("photonCorrectionProcessor") {
  // processor description
  _description = "photonCorrectionProcessor applies an energy correction to photon-like PFOs" ;

  registerProcessorParameter("inputCollection", "name of input PFO collection",
                             _inputCollection, std::string("PandoraPFOs") );

  registerProcessorParameter("modifyPFOenergies" , "apply the corrected energies to the PFOs", _modifyPFOenergies, true );
  registerProcessorParameter("modifyPFOdirection", "apply the corrected direction to the PFOs", _modifyPFOdirections, true );

  registerProcessorParameter("energyLin_const"	          , "overall energy correction: constant term",                 _energyLin_const             , float(  9.870e-01) );
  registerProcessorParameter("energyLin_logen"	          , "overall energy correction: log(e) coefficient",            _energyLin_logen             , float(  1.426e-02) );
  registerProcessorParameter("phiBarrelCorr_pos_const"    , "barrel phi correction: central position (constant)",       _phiBarrelCorr_pos_const     , float( 0.412249  ) );
  registerProcessorParameter("phiBarrelCorr_pos_logen"    , "barrel phi correction: central position (log(e) coeff)",   _phiBarrelCorr_pos_logen     , float( 0.0142289 ) );
  registerProcessorParameter("phiBarrelCorr_depth"	  , "barrel phi correction: gaussian depth",                    _phiBarrelCorr_depth         , float( -0.0933687) );
  registerProcessorParameter("phiBarrelCorr_width1"	  , "barrel phi correction: gaussian width (left side)",        _phiBarrelCorr_width1        , float(  0.01345  ) );
  registerProcessorParameter("phiBarrelCorr_width2"	  , "barrel phi correction: gaussian width (right side)",       _phiBarrelCorr_width2        , float( 0.0408156 ) );
  registerProcessorParameter("costhCorr_gaus1_norm_const" , "barrel cos(theta) correction: gaus1: norm (constant)",     _costhCorr_gaus1_norm_const  , float( -0.0900   ) );
  registerProcessorParameter("costhCorr_gaus1_norm_logen" , "barrel cos(theta) correction: gaus1: norm (log(e) coeff)", _costhCorr_gaus1_norm_logen  , float( 0         ) );
  registerProcessorParameter("costhCorr_gaus1_mean"	  , "barrel cos(theta) correction: gaus1: mean",                _costhCorr_gaus1_mean        , float( 0.235     ) );
  registerProcessorParameter("costhCorr_gaus1_sigm"	  , "barrel cos(theta) correction: gaus1: sigma",               _costhCorr_gaus1_sigm        , float( 0.007256  ) );
  registerProcessorParameter("costhCorr_gaus2_norm_const" , "barrel cos(theta) correction: gaus2: norm (constant)",     _costhCorr_gaus2_norm_const  , float( -0.0369648) );
  registerProcessorParameter("costhCorr_gaus2_norm_logen" , "barrel cos(theta) correction: gaus2: norm (log(e) coeff)", _costhCorr_gaus2_norm_logen  , float( 0         ) );
  registerProcessorParameter("costhCorr_gaus2_mean"	  , "barrel cos(theta) correction: gaus2: mean",                _costhCorr_gaus2_mean        , float( 0.588     ) );
  registerProcessorParameter("costhCorr_gaus2_sigm"	  , "barrel cos(theta) correction: gaus2: sigm",                _costhCorr_gaus2_sigm        , float( 0.0121604 ) );
  registerProcessorParameter("costhCorr_gaus3_norm"	  , "barrel cos(theta) correction: gaus3: norm (constant)",     _costhCorr_gaus3_norm        , float( -0.0422968) );
  registerProcessorParameter("costhCorr_gaus3_mean"	  , "barrel cos(theta) correction: gaus3: mean",                _costhCorr_gaus3_mean        , float( 0.774     ) );
  registerProcessorParameter("costhCorr_gaus3_sigm"	  , "barrel cos(theta) correction: gaus3: sigm",                _costhCorr_gaus3_sigm        , float( 0.009     ) );
  registerProcessorParameter("costhCorr_endcap_scale"	  , "extra correction factor for endcap",                       _costhCorr_endcap_scale      , float( 1.002     ) );
  registerProcessorParameter("endcap_gaus1_norm"	  , "across endcap module correction: gaus1 norm",              _endcap_gaus1_norm           , float( -0.025    ) );
  registerProcessorParameter("endcap_gaus1_mean"	  , "across endcap module correction: gaus1 mean",              _endcap_gaus1_mean           , float( 855.      ) );
  registerProcessorParameter("endcap_gaus1_sigm"	  , "across endcap module correction: gaus1 sigma",             _endcap_gaus1_sigm           , float( 23        ) );
  registerProcessorParameter("endcap_gaus2_norm"	  , "across endcap module correction: gaus2 norm",              _endcap_gaus2_norm           , float( -0.07     ) );
  registerProcessorParameter("endcap_gaus2_mean"	  , "across endcap module correction: gaus2 mean",              _endcap_gaus2_mean           , float( 1489.     ) );
  registerProcessorParameter("endcap_gaus2_sigm"	  , "across endcap module correction: gaus2 sigma",             _endcap_gaus2_sigm           , float( 18.       ) );

  registerProcessorParameter("phiBias_barrel_p0_1"  , "phiBias barrel_par0_1",   _phiBias_barrel_p0_1   , float( 2.36517e-05  ) );
  registerProcessorParameter("phiBias_barrel_p0_2"  , "phiBias barrel_par0_2",   _phiBias_barrel_p0_2   , float(  1.32090e-04 ) );
  registerProcessorParameter("phiBias_barrel_p0_3"  , "phiBias barrel_par0_3",   _phiBias_barrel_p0_3   , float( -3.86883e+00 ) );
  registerProcessorParameter("phiBias_barrel_p0_4"  , "phiBias barrel_par0_4",   _phiBias_barrel_p0_4   , float( -1.67809e-01 ) );
  registerProcessorParameter("phiBias_barrel_p1_1"  , "phiBias barrel_par1_1",   _phiBias_barrel_p1_1   , float( 2.28614e-05 ) );
  registerProcessorParameter("phiBias_barrel_p1_2"  , "phiBias barrel_par1_2",   _phiBias_barrel_p1_2   , float( 6.03495e-05 ) );
  registerProcessorParameter("phiBias_barrel_p2_1"  , "phiBias barrel_par2_1",   _phiBias_barrel_p2_1   , float( 0.419 ) );
  registerProcessorParameter("phiBias_barrel_p2_2"  , "phiBias barrel_par2_2",   _phiBias_barrel_p2_2   , float( 0.00728 ) );
  registerProcessorParameter("phiBias_barrel_p3_1"  , "phiBias barrel_par3_1",   _phiBias_barrel_p3_1   , float( 0.025 ) );
  registerProcessorParameter("phiBias_barrel_p4_1"  , "phiBias barrel_par4_1",   _phiBias_barrel_p4_1   , float( 0.00 ) );
  registerProcessorParameter("phiBias_barrel_p5_1"  , "phiBias barrel_par5_1",   _phiBias_barrel_p5_1   , float( 2.86667e-05 ) );
  registerProcessorParameter("phiBias_barrel_p5_2"  , "phiBias barrel_par5_2",   _phiBias_barrel_p5_2   , float( 2.49371e-05 ) );
  registerProcessorParameter("phiBias_barrel_p5_3"  , "phiBias barrel_par5_3",   _phiBias_barrel_p5_3   , float( -7.71684e-06 ) );
  registerProcessorParameter("phiBias_barrel_p6_1"  , "phiBias barrel_par6_1",   _phiBias_barrel_p6_1   , float( -1.48118e-05 ) );
  registerProcessorParameter("phiBias_barrel_p6_2"  , "phiBias barrel_par6_2",   _phiBias_barrel_p6_2   , float( -5.63786e-06 ) );
  registerProcessorParameter("phiBias_barrel_p7_1"  , "phiBias barrel_par7_1",   _phiBias_barrel_p7_1   , float( -9.38376e-06 ) );
  registerProcessorParameter("phiBias_barrel_p7_2"  , "phiBias barrel_par7_2",   _phiBias_barrel_p7_2   , float( -4.96296e-06 ) );
  registerProcessorParameter("phiBias_barrel_p7_3"  , "phiBias barrel_par7_3",   _phiBias_barrel_p7_3   , float( 2.91262e-06 ) );

  registerProcessorParameter("thetaBias_barrel_p0_1", "thetaBias barrel_par0_1", _thetaBias_barrel_p0_1 , float( -0.000166568  ) );
  registerProcessorParameter("thetaBias_barrel_p0_2", "thetaBias barrel_par0_2", _thetaBias_barrel_p0_2 , float( -7.119e-05 ) );
  registerProcessorParameter("thetaBias_barrel_p1_1", "thetaBias barrel_par1_1", _thetaBias_barrel_p1_1 , float( 0.000223618 ) );
  registerProcessorParameter("thetaBias_barrel_p1_2", "thetaBias barrel_par1_2", _thetaBias_barrel_p1_2 , float( -3.95915e-05 ) );

  registerProcessorParameter("thetaBias_endcap_p0_1", "thetaBias endcap_par0_1", _thetaBias_endcap_p0_1 , float( 0.000129478 ) );
  registerProcessorParameter("thetaBias_endcap_p0_2", "thetaBias endcap_par0_2", _thetaBias_endcap_p0_2 , float( -3.73863e-05 ) );
  registerProcessorParameter("thetaBias_endcap_p1_1", "thetaBias endcap_par1_1", _thetaBias_endcap_p1_1 , float( -0.000847783 ) );
  registerProcessorParameter("thetaBias_endcap_p1_2", "thetaBias endcap_par1_2", _thetaBias_endcap_p1_2 , float( 0.000153646 ) );
  registerProcessorParameter("thetaBias_endcap_p2_1", "thetaBias endcap_par2_1", _thetaBias_endcap_p2_1 , float( 0.000806605 ) );
  registerProcessorParameter("thetaBias_endcap_p2_2", "thetaBias endcap_par2_2", _thetaBias_endcap_p2_2 , float( -0.000132608 ) );

  registerProcessorParameter("validationPlots", "produce validation plots", _validationPlots, false );
  registerProcessorParameter("nominalEnergy", "nominal photon energy (for validation plots)", _nominalEnergy, float(200) );
}


void photonCorrectionProcessor::init() {
  streamlog_out (MESSAGE) << "hello from photonCorrectionProcessor::init" << endl;

  // first get some geometry information
  dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();
  // endcap ecal
  unsigned int includeFlag = ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ENDCAP | dd4hep::DetType::ELECTROMAGNETIC );
  unsigned int excludeFlag = ( dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD );
  const std::vector< dd4hep::DetElement>& endcapEcalDetectors = dd4hep::DetectorSelector(mainDetector).detectors(  includeFlag, excludeFlag );
  if ( endcapEcalDetectors.size() == 1 ) {
    _assumed_boxsize  = endcapEcalDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>()->extent[0]/dd4hep::mm; // r-min
    _assumed_endZ     = endcapEcalDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>()->extent[2]/dd4hep::mm; // z-min 
  } else {
    streamlog_out (ERROR) << "did not find exactly one endcap ECAL! found " << endcapEcalDetectors.size() << "; refusing to continue!" << endl;
    assert(0);
  }

  // barrel ecal
  includeFlag = ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL | dd4hep::DetType::ELECTROMAGNETIC );
  excludeFlag = ( dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD );
  const std::vector< dd4hep::DetElement>& barrelEcalDetectors = dd4hep::DetectorSelector(mainDetector).detectors(  includeFlag, excludeFlag );
  if ( barrelEcalDetectors.size() == 1 ) {
    float barrelLength = barrelEcalDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>()->extent[3]/dd4hep::mm; // z-max
    float barrelInRad  = barrelEcalDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>()->extent[0]/dd4hep::mm; // r-min
    _barrelendcap_limit_costh=barrelLength/sqrt( pow(barrelLength,2) + pow(barrelInRad,2) );
  } else {
    streamlog_out (ERROR) << "did not find exactly one barrel ECAL! found " << barrelEcalDetectors.size() << "; refusing to continue!" << endl;
    assert(0);
  }

  streamlog_out (MESSAGE) << "*** ECAL endcap dimensions: minZ, minR, barrel/endcap transition " << 
    _assumed_endZ << " , " << _assumed_boxsize << " , " << _barrelendcap_limit_costh << endl;

  _photonCorrector = new photonCorrector();

  _photonCorrector->set_barrelendcap_limit         ( _barrelendcap_limit_costh   );
  _photonCorrector->set_assumed_boxsize            ( _assumed_boxsize            );
  _photonCorrector->set_assumed_endZ               ( _assumed_endZ               );

  _photonCorrector->set_energyLin_const            ( _energyLin_const            );
  _photonCorrector->set_energyLin_logen            ( _energyLin_logen            );
  _photonCorrector->set_phiBarrelCorr_pos_const    ( _phiBarrelCorr_pos_const    );
  _photonCorrector->set_phiBarrelCorr_pos_logen    ( _phiBarrelCorr_pos_logen    );
  _photonCorrector->set_phiBarrelCorr_depth        ( _phiBarrelCorr_depth        );
  _photonCorrector->set_phiBarrelCorr_width1       ( _phiBarrelCorr_width1       );
  _photonCorrector->set_phiBarrelCorr_width2       ( _phiBarrelCorr_width2       );
  _photonCorrector->set_costhCorr_gaus1_norm_const ( _costhCorr_gaus1_norm_const );
  _photonCorrector->set_costhCorr_gaus1_norm_logen ( _costhCorr_gaus1_norm_logen );
  _photonCorrector->set_costhCorr_gaus1_mean       ( _costhCorr_gaus1_mean       );
  _photonCorrector->set_costhCorr_gaus1_sigm       ( _costhCorr_gaus1_sigm       );
  _photonCorrector->set_costhCorr_gaus2_norm_const ( _costhCorr_gaus2_norm_const );
  _photonCorrector->set_costhCorr_gaus2_norm_logen ( _costhCorr_gaus2_norm_logen );
  _photonCorrector->set_costhCorr_gaus2_mean       ( _costhCorr_gaus2_mean       );
  _photonCorrector->set_costhCorr_gaus2_sigm       ( _costhCorr_gaus2_sigm       );
  _photonCorrector->set_costhCorr_gaus3_norm       ( _costhCorr_gaus3_norm       );
  _photonCorrector->set_costhCorr_gaus3_mean       ( _costhCorr_gaus3_mean       );
  _photonCorrector->set_costhCorr_gaus3_sigm       ( _costhCorr_gaus3_sigm       );
  _photonCorrector->set_costhCorr_endcap_scale     ( _costhCorr_endcap_scale     );
  _photonCorrector->set_endcap_gaus1_norm          ( _endcap_gaus1_norm          );
  _photonCorrector->set_endcap_gaus1_mean          ( _endcap_gaus1_mean          );
  _photonCorrector->set_endcap_gaus1_sigm          ( _endcap_gaus1_sigm          );
  _photonCorrector->set_endcap_gaus2_norm          ( _endcap_gaus2_norm          );
  _photonCorrector->set_endcap_gaus2_mean          ( _endcap_gaus2_mean          );
  _photonCorrector->set_endcap_gaus2_sigm          ( _endcap_gaus2_sigm          );

  _photonCorrector->set_phiBias_barrel_p0_1        (_phiBias_barrel_p0_1  );
  _photonCorrector->set_phiBias_barrel_p0_2        (_phiBias_barrel_p0_2  );
  _photonCorrector->set_phiBias_barrel_p0_3        (_phiBias_barrel_p0_3  );
  _photonCorrector->set_phiBias_barrel_p0_4        (_phiBias_barrel_p0_4  );
  _photonCorrector->set_phiBias_barrel_p1_1        (_phiBias_barrel_p1_1  );
  _photonCorrector->set_phiBias_barrel_p1_2        (_phiBias_barrel_p1_2  );
  _photonCorrector->set_phiBias_barrel_p2_1        (_phiBias_barrel_p2_1  );
  _photonCorrector->set_phiBias_barrel_p2_2        (_phiBias_barrel_p2_2  );
  _photonCorrector->set_phiBias_barrel_p3_1        (_phiBias_barrel_p3_1  );
  _photonCorrector->set_phiBias_barrel_p4_1        (_phiBias_barrel_p4_1  );
  _photonCorrector->set_phiBias_barrel_p5_1        (_phiBias_barrel_p5_1  );
  _photonCorrector->set_phiBias_barrel_p5_2        (_phiBias_barrel_p5_2  );
  _photonCorrector->set_phiBias_barrel_p5_3        (_phiBias_barrel_p5_3  );
  _photonCorrector->set_phiBias_barrel_p6_1        (_phiBias_barrel_p6_1  );
  _photonCorrector->set_phiBias_barrel_p6_2        (_phiBias_barrel_p6_2  );
  _photonCorrector->set_phiBias_barrel_p7_1        (_phiBias_barrel_p7_1  );
  _photonCorrector->set_phiBias_barrel_p7_2        (_phiBias_barrel_p7_2  );
  _photonCorrector->set_phiBias_barrel_p7_3        (_phiBias_barrel_p7_3  );
  _photonCorrector->set_thetaBias_barrel_p0_1      (_thetaBias_barrel_p0_1);
  _photonCorrector->set_thetaBias_barrel_p0_2      (_thetaBias_barrel_p0_2);
  _photonCorrector->set_thetaBias_barrel_p1_1      (_thetaBias_barrel_p1_1);
  _photonCorrector->set_thetaBias_barrel_p1_2      (_thetaBias_barrel_p1_2);
  _photonCorrector->set_thetaBias_endcap_p0_1      (_thetaBias_endcap_p0_1);
  _photonCorrector->set_thetaBias_endcap_p0_2      (_thetaBias_endcap_p0_2);
  _photonCorrector->set_thetaBias_endcap_p1_1      (_thetaBias_endcap_p1_1);
  _photonCorrector->set_thetaBias_endcap_p1_2      (_thetaBias_endcap_p1_2);
  _photonCorrector->set_thetaBias_endcap_p2_1      (_thetaBias_endcap_p2_1);
  _photonCorrector->set_thetaBias_endcap_p2_2      (_thetaBias_endcap_p2_2);


  return;
}

void photonCorrectionProcessor::processRunHeader( LCRunHeader* /* run */ ) {
  streamlog_out (MESSAGE) << "hello from photonCorrectionProcessor::processRunHeader" << endl;
}

void photonCorrectionProcessor::processEvent( LCEvent * evt ) {
  // adjust energy of all type==22 objects in the input collection

  // for validation plots: sum the pfo energies before and after correction
  float totPfoEn[2]={0};
  float totPfoGammaEn[2]={0};
  float totomom[3]={0};

  try {
    LCCollectionVec* pfocol =  dynamic_cast<LCCollectionVec*>(evt->getCollection( _inputCollection ) );
    for (auto lcobj : (*pfocol) ) {
      ReconstructedParticleImpl* pfo = dynamic_cast<ReconstructedParticleImpl*>(lcobj);

      if ( _validationPlots ) {
	for (int j=0; j<3; j++) {
	  totomom[j]+=pfo->getMomentum()[j];
	}
      }

      if ( pfo->getType() == 22 ) {

	float origEn = pfo->getEnergy();
	float origMom(0);
	for (int i=0; i<3; i++) {
	  origMom+=pow(pfo->getMomentum()[i],2);
	}
	origMom=sqrt(origMom);
	float origTheta = acos( pfo->getMomentum()[2]/origMom);
	float origPhi = atan2( pfo->getMomentum()[1], pfo->getMomentum()[0] );

	float correctedEnergy = _photonCorrector->photonEnergyCorrection( pfo );

	float correctedMom = origMom*correctedEnergy/origEn; // just scale. in principle should be same as energy for massless photon PFO

	float correctedTheta, correctedPhi;
	_photonCorrector->photonDirectionCorrection( pfo , correctedTheta, correctedPhi );

	// these are the properties we will set for the PFO
	float newEnergy = _modifyPFOenergies ? correctedEnergy : origEn;
	float newMom    = _modifyPFOenergies ? correctedMom : origMom;
	float newTheta  = _modifyPFOdirections ? correctedTheta : origTheta;
	float newPhi    = _modifyPFOdirections ? correctedPhi : origPhi;

	float newMomentum[3];
	newMomentum[0] = newMom*sin(newTheta)*cos(newPhi);
	newMomentum[1] = newMom*sin(newTheta)*sin(newPhi);
	newMomentum[2] = newMom*cos(newTheta);

	pfo->setEnergy( newEnergy );
	pfo->setMomentum( newMomentum );

	if ( _validationPlots ) {
	  totPfoEn     [0]+=pfo->getEnergy();
	  totPfoGammaEn[0]+=pfo->getEnergy();
	  totPfoEn     [1]+=correctedEnergy;
	  totPfoGammaEn[1]+=correctedEnergy;
	}

      } else {
	totPfoEn[0]+=pfo->getEnergy();
	totPfoEn[1]+=pfo->getEnergy();
      }
    }
  } catch(DataNotAvailableException &e) {};


#ifdef MARLIN_USE_AIDA

  if ( _validationPlots ) {
    static AIDA::IHistogram2D* h_sumPfoE_orig ;
    static AIDA::IHistogram2D* h_sumPfoE_corr ;
    static AIDA::IHistogram2D* h_sumGamPfoE_orig ;
    static AIDA::IHistogram2D* h_sumGamPfoE_corr ;
    if ( isFirstEvent() ) {
      h_sumPfoE_orig    =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "sumPfoE_orig"   , "energy [GeV]", 200, -1, 1, 100, 0., 1.5*_nominalEnergy );
      h_sumPfoE_corr    =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "sumPfoE_corr"   , "energy [GeV]", 200, -1, 1, 100, 0., 1.5*_nominalEnergy );
      h_sumGamPfoE_orig =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "sumGamPfoE_orig", "energy [GeV]", 200, -1, 1, 100, 0., 1.5*_nominalEnergy );
      h_sumGamPfoE_corr =  AIDAProcessor::histogramFactory(this)->createHistogram2D( "sumGamPfoE_corr", "energy [GeV]", 200, -1, 1, 100, 0., 1.5*_nominalEnergy );
    }

    float costh = totomom[2]/sqrt( pow( totomom[0], 2 ) +
				   pow( totomom[1], 2 ) +
				   pow( totomom[2], 2 ) );

    h_sumPfoE_orig   ->fill( costh, totPfoEn[0] );
    h_sumPfoE_corr   ->fill( costh, totPfoEn[1] );
    h_sumGamPfoE_orig->fill( costh, totPfoGammaEn[0] );
    h_sumGamPfoE_corr->fill( costh, totPfoGammaEn[1] );
  }

#endif

  return;
}

void photonCorrectionProcessor::check( LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void photonCorrectionProcessor::end(){
  streamlog_out (MESSAGE) << "photonCorrectionProcessor::end()  " << name() << std::endl ;
}

