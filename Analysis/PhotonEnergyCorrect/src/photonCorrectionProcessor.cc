#include "photonCorrectionProcessor.h"
#include <iostream>
#include <marlin/Global.h>
#include "lcio.h"

#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"

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

  registerProcessorParameter("useCorrectorDefaults" , "use defaults correction paramteres (as defined in photonCorrector", 
			     _useCorrectorDefaults, true );

  registerProcessorParameter("interModGaps_haveBeenCorrectedAtHitLevel", 
			     "have the inter-Module gaps already been corrected at this hit level? (affects default param set used)",
			     _interModGaps_haveBeenCorrectedAtHitLevel, true);

  registerProcessorParameter("barrel_limit"               , "_barrel_limit               parameter",  _barrel_limit                , float( -999. ) );
  registerProcessorParameter("endcap_limit"		  , "_endcap_limit               parameter",  _endcap_limit                , float( -999. ) );
  registerProcessorParameter("energyLin_const"	          , "_energyLin_const            parameter",  _energyLin_const             , float( -999. ) );
  registerProcessorParameter("energyLin_logen"	          , "_energyLin_logen            parameter",  _energyLin_logen             , float( -999. ) );
  registerProcessorParameter("phiBarrelCorr_pos_const"    , "_phiBarrelCorr_pos_const    parameter",  _phiBarrelCorr_pos_const     , float( -999. ) );
  registerProcessorParameter("phiBarrelCorr_pos_logen"    , "_phiBarrelCorr_pos_logen    parameter",  _phiBarrelCorr_pos_logen     , float( -999. ) );
  registerProcessorParameter("phiBarrelCorr_depth"	  , "_phiBarrelCorr_depth        parameter",  _phiBarrelCorr_depth         , float( -999. ) );
  registerProcessorParameter("phiBarrelCorr_width1"	  , "_phiBarrelCorr_width1       parameter",  _phiBarrelCorr_width1        , float( -999. ) );
  registerProcessorParameter("phiBarrelCorr_width2"	  , "_phiBarrelCorr_width2       parameter",  _phiBarrelCorr_width2        , float( -999. ) );
  registerProcessorParameter("costhCorr_gaus1_norm_const" , "_costhCorr_gaus1_norm_const parameter",  _costhCorr_gaus1_norm_const  , float( -999. ) );
  registerProcessorParameter("costhCorr_gaus1_norm_logen" , "_costhCorr_gaus1_norm_logen parameter",  _costhCorr_gaus1_norm_logen  , float( -999. ) );
  registerProcessorParameter("costhCorr_gaus1_mean"	  , "_costhCorr_gaus1_mean       parameter",  _costhCorr_gaus1_mean        , float( -999. ) );
  registerProcessorParameter("costhCorr_gaus1_sigm"	  , "_costhCorr_gaus1_sigm       parameter",  _costhCorr_gaus1_sigm        , float( -999. ) );
  registerProcessorParameter("costhCorr_gaus2_norm_const" , "_costhCorr_gaus2_norm_const parameter",  _costhCorr_gaus2_norm_const  , float( -999. ) );
  registerProcessorParameter("costhCorr_gaus2_norm_logen" , "_costhCorr_gaus2_norm_logen parameter",  _costhCorr_gaus2_norm_logen  , float( -999. ) );
  registerProcessorParameter("costhCorr_gaus2_mean"	  , "_costhCorr_gaus2_mean       parameter",  _costhCorr_gaus2_mean        , float( -999. ) );
  registerProcessorParameter("costhCorr_gaus2_sigm"	  , "_costhCorr_gaus2_sigm       parameter",  _costhCorr_gaus2_sigm        , float( -999. ) );
  registerProcessorParameter("costhCorr_gaus3_norm"	  , "_costhCorr_gaus3_norm       parameter",  _costhCorr_gaus3_norm        , float( -999. ) );
  registerProcessorParameter("costhCorr_gaus3_mean"	  , "_costhCorr_gaus3_mean       parameter",  _costhCorr_gaus3_mean        , float( -999. ) );
  registerProcessorParameter("costhCorr_gaus3_sigm"	  , "_costhCorr_gaus3_sigm       parameter",  _costhCorr_gaus3_sigm        , float( -999. ) );
  registerProcessorParameter("costhCorr_endcap_scale"	  , "_costhCorr_endcap_scale     parameter",  _costhCorr_endcap_scale      , float( -999. ) );
  registerProcessorParameter("endcap_gaus1_norm"	  , "_endcap_gaus1_norm          parameter",  _endcap_gaus1_norm           , float( -999. ) );
  registerProcessorParameter("endcap_gaus1_mean"	  , "_endcap_gaus1_mean          parameter",  _endcap_gaus1_mean           , float( -999. ) );
  registerProcessorParameter("endcap_gaus1_sigm"	  , "_endcap_gaus1_sigm          parameter",  _endcap_gaus1_sigm           , float( -999. ) );
  registerProcessorParameter("endcap_gaus2_norm"	  , "_endcap_gaus2_norm          parameter",  _endcap_gaus2_norm           , float( -999. ) );
  registerProcessorParameter("endcap_gaus2_mean"	  , "_endcap_gaus2_mean          parameter",  _endcap_gaus2_mean           , float( -999. ) );
  registerProcessorParameter("endcap_gaus2_sigm"	  , "_endcap_gaus2_sigm          parameter",  _endcap_gaus2_sigm           , float( -999. ) );
  registerProcessorParameter("assumed_boxsize"	          , "_assumed_boxsize            parameter",  _assumed_boxsize             , float( -999. ) );
  registerProcessorParameter("assumed_endZ"		  , "_assumed_endZ               parameter",  _assumed_endZ                , float( -999. ) );

  registerProcessorParameter("modifyPFOenergies", "apply the corrected energies to the PFOs", _modifyPFOenergies, true );
  registerProcessorParameter("validationPlots", "produce validation plots", _validationPlots, false );
  registerProcessorParameter("validationPlotFilename", "file name for validation plots", _validationPlotFilename, std::string("photonCorrectionValidation.root") );
  registerProcessorParameter("nominalEnergy", "nominal photon energy (for validation plots)", _nominalEnergy, float(200) );
  


}


void photonCorrectionProcessor::init() {
  cout << "hello from photonCorrectionProcessor::init" << endl;

  _photonCorrector = new photonCorrector();

  if ( _useCorrectorDefaults ) {
    if ( _interModGaps_haveBeenCorrectedAtHitLevel ) {
      _photonCorrector->setDefaultValues( photonCorrector::set_std );
    } else {
      _photonCorrector->setDefaultValues( photonCorrector::set_noInterMod );
    }
  } else {
    _photonCorrector->set_barrel_limit               ( _barrel_limit               );
    _photonCorrector->set_endcap_limit               ( _endcap_limit               );
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
    _photonCorrector->set_assumed_boxsize            ( _assumed_boxsize            );
    _photonCorrector->set_assumed_endZ               ( _assumed_endZ               );
  }


  if ( _validationPlots ) {
    _fout = new TFile(_validationPlotFilename.c_str(), "recreate");
    for (int i=0; i<2; i++) {
      TString dd = i==0 ? "beforeCorr" : "afterCorr" ;
      _hEnSum[0][i] = new TH2F("sumPfoEn_"+dd, "sumPfoEn_"+dd, 100, -1, 1, 500, _nominalEnergy/2., _nominalEnergy*1.5 );
      _hEnSum[1][i] = new TH2F("sumGamEn_"+dd, "sumGamEn_"+dd, 100, -1, 1, 500, _nominalEnergy/2., _nominalEnergy*1.5 );
    }
  }



  return;
}

void photonCorrectionProcessor::processRunHeader( LCRunHeader* /* run */ ) {
  cout << "hello from photonCorrectionProcessor::processRunHeader" << endl;
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

	float correctedEnergy = _photonCorrector->photonEnergyCorrection( pfo );

	if ( _validationPlots ) {
	  totPfoEn     [0]+=pfo->getEnergy();
	  totPfoGammaEn[0]+=pfo->getEnergy();
	  totPfoEn     [1]+=correctedEnergy;
	  totPfoGammaEn[1]+=correctedEnergy;
	}

	if (_modifyPFOenergies) pfo->setEnergy( correctedEnergy );


      } else {
	totPfoEn[0]+=pfo->getEnergy();
	totPfoEn[1]+=pfo->getEnergy();
      }
    }
  } catch(DataNotAvailableException &e) {};

  if ( _validationPlots ) {

    float costh = totomom[2]/sqrt( pow( totomom[0], 2 ) +
				   pow( totomom[1], 2 ) +
				   pow( totomom[2], 2 ) );

    for (int i=0; i<2; i++) {
      _hEnSum[0][i]->Fill( costh, totPfoEn[i] );
      _hEnSum[1][i]->Fill( costh, totPfoGammaEn[i] );
    }
  }

  return;
}



void photonCorrectionProcessor::check( LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void photonCorrectionProcessor::end(){
  std::cout << "photonCorrectionProcessor::end()  " << name()
            << std::endl ;

  if ( _validationPlots ) {
    _fout->Write();
    _fout->Close();
  }

}

