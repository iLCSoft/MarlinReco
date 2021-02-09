/*
 * =====================================================================================
 *
 *       Filename:  ErrorFlow.h
 *
 *    Description:  The ErrorFlow processor computes jet-specific energy uncertainty
 *		    by summing up the uncertainty of individual particles clustered
 *		    in a jet object. 
 *		    The processor updates the covariance matrix of the jet object with
 *		    the computed energy uncertainty and creates a ROOT tree with
 *		    breanches for PFO multiplicities and various terms contributing to
 *		    the total jet energy uncertainty.
 *
 *        Version:  1.0
 *        Created:  08/20/2015 03:22:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Aliakbar Ebrahimi (aliakbar.ebrahimi@desy.de), 
 *   Organization:  DESY
 *
 * =====================================================================================
 */


#ifndef ERRORFLOW_H
#define ERRORFLOW_H 1

/* #####   HEADER FILE INCLUDES   ################################################### */
#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include <TTree.h>
#include <EVENT/ReconstructedParticle.h>


using namespace lcio ;
using namespace marlin ;


/*
 * =====================================================================================
 *        Class:  ErrorFlow
 *  Description:  
 * =====================================================================================
 */
class ErrorFlow : public Processor
{
	public:
		/* ====================  LIFECYCLE     ======================================= */
		virtual Processor*  newProcessor() { return new ErrorFlow; }

		ErrorFlow ();                             /* constructor */

		/** Called at the begin of the job before anything is read.
		* Use to initialize the processor, e.g. book histograms.
		*/
		virtual void init() ;

		/** Called for every run.
		*/
		virtual void processRunHeader( LCRunHeader* run ) ;
  
		/** Called after data processing for clean up.
		*/
		virtual void end() ;

		/* ====================  ACCESSORS     ======================================= */

		/* ====================  MUTATORS      ======================================= */
		void resetVariables();

		/* ====================  OPERATORS     ======================================= */
		/* Computes confusion terms for a jet */
		double * getRelativeConfusion(	double jetEnergy,
						double chargedHadronsEnergy,
						double photonsEnergy,
						double neutralHadronsEnergy,
						double relConfTerms[]
						);

		/** Called for every event - the working horse.
		*/
		virtual void processEvent( LCEvent * evt ) ; 

		virtual void check( LCEvent * evt ) ; 
  
	protected:
		/* ====================  METHODS       ======================================= */
		// Computes energy uncertainty SQUARED for photons
		double getPhotonSigmaESqr ( double t_phEnergy );
		// Computes energy uncertainty SQUARED for neutral hadrons
		double getNeuHadSigmaESqr ( double t_neuHadEnergy );

		// Compute total momentum from 3-momentum
		double getTotalMomentum ( const double * t_threeMomentum );

		/* ====================  DATA MEMBERS  ======================================= */

		/** Input collection name.
		*/
		std::string p_inputJetCollection {};
		std::string p_outputJetCollection {};
		std::string p_inputMCTruth {};

		// Semi-leptonic correction
		bool p_semiLepCorrection {};                 /* Add semi-leptonic energy resolution to the covariance matrix */
		bool p_confusionterm {};                 /* Add uncertainty due to confusion to the covariance matrix */
		bool p_propagateConfusiontoMomentumComp {};  /* Propagate uncertainty due to confusion to the Momentum components */
		double p_semiLepSigmaCorrFactor {};          /* A correction factor to be multiplied by total lepton energy to get semi-leptonic uncertainty */
		double p_CovMatFactorPhotons {};          /* A correction factor to be multiplied to angular uncertainties of photons */
		double p_CovMatFactorNeutralHadrons {};          /* A correction factor to be multiplied to angular uncertainties of Neutral Hadrons */

		// Confusion scale factor
		double p_scaleConf {};                        /* A factor to use to scale confusion term */
	  
		// Calorimeter resolution parameters
		double p_aECAL {};                      /* Stochastic term coefficient in ECAL resolution */
		double p_cECAL {};                      /* Constant term coefficinet in ECAL resolution */
		double p_aHCAL {};                      /* Stochastic term coefficient in HCAL resolution */
		double p_cHCAL {};                      /* Constant term coefficient in HCAL resolution */

		// Confusion scale factor
		bool p_storeTree {};                        /* Enable/disable storing in a ROOT tree */
		bool p_useFullCovMatNeut {};                /* Enable/disable using full CovMat for neutral PFOs */

		// Counters for PFOs
		int numChargedPFOs {};                       /* Number of charged PFOs in a jet */
		int numPhotons {};                           /* Number of photons in a jet */
		int numNeutralPFOs {};                       /* Number of charged PFOs in a jet */

		// Energy of PFOs
		double eChargedPFOs {};                      /* Total energy of charged PFOs in a jet  */
		double ePhotons {};                          /* Total energy of photons in a jet  */
		double eNeutralPFOs {};                      /* Total energy of neutral PFOs in a jet  */
		double eJetTotal {};                         /* Sum of energy of all PFOs */
		double eLeptonsTotal {};                     /* Sum of energy of all leptons in a jet */

		// Resolution
		double absDetResSquared {};                  /* Total absolute detector resolution */
		double relConfCharged {};                    /* Relative confusion term for charged PFOs */
		double relConfPhotons {};                    /* Relative confusion term for photons */
		double relConfNeutral {};                    /* Relative confusion term for neutral PFOs */
		double relConfSquared {};                    /* Total relative confusion squared  */
		double absConfSquared {};                    /* Total absolute confusion squared  */
		double absSemiLepResSquared {};              /* Total absolute semi-leptonic resolution  */

		int p_nRun {};
		int p_nEvt {};

	private:
		/* ====================  METHODS       ======================================= */

		/* ====================  DATA MEMBERS  ======================================= */
		std::shared_ptr<TTree> tree {};
		ReconstructedParticleVec::size_type nPFOs {};

}; /* -----  end of class ErrorFlow  ----- */


#endif
