/*
 * =====================================================================================
 *
 *       Filename:  ErrorFlow.cc
 *
 *    Description:  This is the implementation of the ErrorFlow Marlin Processor.
 *					See the header file for more information.
 *
 *        Version:  1.0
 *        Created:  11/02/2016 02:00:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Aliakbar Ebrahimi (aliakbar.ebrahimi@desy.de), 
 *   Organization:  DESY
 *
 * =====================================================================================
 */



#include "ErrorFlow.h"
#include <iostream>
#include <vector>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h> 
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/MCParticle.h>
#include "TVector3.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio ;
using namespace marlin ;


ErrorFlow aErrorFlow;


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ErrorFlow::ErrorFlow()
 *  Description: This is the class constructor. The processor parameters are
 *				 defined here.
 * =====================================================================================
 */
ErrorFlow::ErrorFlow() : Processor("ErrorFlow") {

    // Processor description
    _description = "ErrorFlow processor computes total error of a jet based on the information from each sub-detector" ;

    // Register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
            "InputPFOCollection", 
            "Input collection of jets",
            p_inputJetCollection,
            std::string( "InputJets" )
		  );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
            "OutputPFOCollection", 
            "Out collection of jets with calculated covariance matrix",
            p_outputJetCollection,
            std::string( "OutputJets" )
		  );

	registerInputCollection( LCIO::LCRELATION,
		  "InputMCTruthLinkCollection" ,
		  "Input collection of MC Truth Links",
		  p_inputMCTruth,
		  std::string( "RecoMCTruthLink" )
		  );

	// Reference of default value: doi:10.1016/j.nima.2009.07.026
	registerProcessorParameter( "ECALStochastic" ,
		  "Stochastic term coefficient in ECAL resolution",
		  p_aECAL,
		  (double) 0.16 );

	// Reference of default value: doi:10.1016/j.nima.2009.07.026
	registerProcessorParameter( "ECALConstant" ,
		  "Constant term coefficient in ECAL resolution",
		  p_cECAL,
		  (double) 0.02 );

	// Reference of defautl value: 2012 JINST 7 P09017
	registerProcessorParameter( "HCALStochastic" ,
		  "Stochastic term coefficient in HCAL resolution",
		  p_aHCAL,
		  (double) 0.50 );

	// Reference of defautl value: 2012 JINST 7 P09017
	registerProcessorParameter( "HCALConstant" ,
		  "Constant term coefficient in HCAL resolution",
		  p_cHCAL,
		  (double) 0.01 );

	registerProcessorParameter( "ConfusionScaleFactor" ,
		  "Scale factor for confusion resolution term",
		  p_scaleConf,
		  (double) 1.0 );

	registerProcessorParameter( "EnableSemiLepCorrection" ,
		  "Enable/disable semi-leptonic correction resolution to be added to covariance matrix",
		  p_semiLepCorrection,
		  (bool) false );

	registerProcessorParameter( "EnableConfusionTerm" ,
		  "Enable/disable confusion term to be added to covariance matrix",
		  p_confusionterm,
		  (bool) true );

	registerProcessorParameter( "PropagateConfusion2Mom" ,
		  "Enable/disable Propagating uncertainty due to confusion to the Momentum components",
		  p_propagateConfusiontoMomentumComp,
		  (bool) true );

	registerProcessorParameter( "SemiLepSigmaCorrectionFactor" ,
		  "A correction factor to be multiplied by total lepton energy to get semi-leptonic uncertainty",
		  p_semiLepSigmaCorrFactor,
		  (double) 0.57 );

	registerProcessorParameter( "CovMatFactorPhotons" ,
		  "A correction factor to be multiplied to angular uncertainties of photons",
		  p_CovMatFactorPhotons,
		  (double) 1.0 );

	registerProcessorParameter( "CovMatFactorNeutralHadrons" ,
		  "A correction factor to be multiplied to angular uncertainties of Neutral Hadrons",
		  p_CovMatFactorNeutralHadrons,
		  (double) 1.0 );

	registerProcessorParameter( "SotreInTree" ,
		  "Enable/disable storing computed qunatities in a ROOT tree",
		  p_storeTree,
		  (bool) false );

	registerProcessorParameter( "useFullCovMatforNeutrals" ,
		  "whether use full CovMat for neutral PFOs or use only energy uncertainty",
		  p_useFullCovMatNeut,
		  (bool) true );


} /* -----  end of function ErrorFlow::ErrorFlow  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ErrorFlow::init()
 *  Description:  The initialisation function for the processor. ROOT tree and
 *				  branches are created here.
 * =====================================================================================
 */
void ErrorFlow::init() { 

    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    // usually a good idea to
    printParameters() ;

	// Create a ROOT file if enabled in the steering file
	if ( p_storeTree ) {
	   tree = std::make_shared<TTree>("ErrorFlow", "Number of photons, charged and neutral hadrons in a jet");
	   // PFO count branches
	   tree->Branch( "numPhotons" , &numPhotons, "numPhotons/I");
	   tree->Branch( "numChargedPFOs" , &numChargedPFOs, "numChargedPFOs/I");
	   tree->Branch( "numNeutralPFOs" , &numNeutralPFOs, "numNeutralPFOs/I");
	   tree->Branch( "numPFOsTotal" , &nPFOs, "numPFOsTotal/I");
	   // Energy branches
	   tree->Branch( "ePhotons" , &ePhotons, "ePhotons/D");
	   tree->Branch( "eChargedPFOs" , &eChargedPFOs, "eChargedPFOs/D");
	   tree->Branch( "eNeutralPFOs" , &eNeutralPFOs, "eNeutralPFOs/D");
	   tree->Branch( "ePFOsTotal" , &eJetTotal, "eJetTotal/D");
	   // Jet energy resolution branches
	   tree->Branch( "absDetResSquared" , &absDetResSquared, "absDetResSquared/D");
	   tree->Branch( "relConfCharged" , &relConfCharged, "relConfCharged/D");
	   tree->Branch( "relConfPhotons" , &relConfPhotons, "relConfPhotons/D");
	   tree->Branch( "relConfNeutral" , &relConfNeutral, "relConfNeutral/D");
	   tree->Branch( "relConfSquared" , &relConfSquared, "relConfSquared/D");
	   tree->Branch( "absConfSquared" , &absConfSquared, "absConfSquared/D");
	   tree->Branch( "absSemiLepResSquared" , &absSemiLepResSquared, "absSemiLepResSquared/D");
	} // end if p_storeTree

    p_nRun = 0 ;
    p_nEvt = 0 ;

} /* -----  end of function ErrorFlow::init  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ErrorFlow::processRunHeader( LCRunHeader* run )
 *  Description:  
 * =====================================================================================
 */
void ErrorFlow::processRunHeader( LCRunHeader* /*run*/ ) { 

    p_nRun++ ;

} /* -----  end of function ErrorFlow::processRunHeader  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ErrorFlow::processEvent( LCEvent* evt )
 *  Description:  This function loops over all particles of each jet in an events and
 *                1- Computes error on the energy of each particle by using covariance
 *                matrix of the particle for charged PFOs, ECAL energy resolution for
 *                photons and HCAL energy resolution for neutral particles
 *                1- Counts number of charged PFOs, photons and neutral PFOs
 *                2- Summs up energy of charged PFOs, photons and neutral PFOs and
 *                stores them in separate root branches.
 *                This function gets called for every event.
 * =====================================================================================
 */
void ErrorFlow::processEvent( LCEvent * evt ) { 

   // Debug message
   streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
	  << "   in run:  " << evt->getRunNumber() << std::endl;

   // Get input collection (exits if collection is not available)
   LCCollection *col = evt->getCollection( p_inputJetCollection );
   
   // Setup MC relation for semi leptonic correction
   LCCollection *colMCTL = evt->getCollection( p_inputMCTruth );
   LCRelationNavigator *navMCTL = new LCRelationNavigator( colMCTL );

   // Create the output collection
   LCCollectionVec* outCol= new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
   outCol->setSubset( true );

   // If the collection is available, do the rest of the work
   if( NULL != col ) {

	  // Get number of jets in the collection
	  size_t nJets = col->getNumberOfElements();
	  streamlog_out(DEBUG) << "Number of jets in the event:" << nJets << std::endl;

	  // Loop over all jets in the event
	  for( size_t iJet=0; iJet < nJets; ++iJet ) {

		 // Reset variables
		 resetVariables();

		 // Set a pointer to the i-th jet in the collection
		 ReconstructedParticleImpl *jetPtr = dynamic_cast<ReconstructedParticleImpl*>( col->getElementAt( iJet ) );
		 TVector3 jet3Momentum( jetPtr->getMomentum()[0] , jetPtr->getMomentum()[1] , jetPtr->getMomentum()[2] );
		 double jetMom = jet3Momentum.Mag();
		 double jetEnergy = jetPtr->getEnergy();
		 float jetTheta = jet3Momentum.Theta();
		 float jetPhi = jet3Momentum.Phi();
		 std::vector< float > jetCovMatrixConf( 10, 0.0 );
		 jetCovMatrixConf[ 0 ] = pow( sin( jetTheta ) , 2 ) * pow( cos( jetPhi ) , 2 );
		 jetCovMatrixConf[ 1 ] = pow( sin( jetTheta ) , 2 ) * sin( jetPhi ) * cos( jetPhi );
		 jetCovMatrixConf[ 2 ] = pow( sin( jetTheta ) , 2 ) * pow( sin( jetPhi ) , 2 );
		 jetCovMatrixConf[ 3 ] = sin( jetTheta ) * cos( jetTheta ) * cos( jetPhi );
		 jetCovMatrixConf[ 4 ] = sin( jetTheta ) * cos( jetTheta ) * sin( jetPhi );
		 jetCovMatrixConf[ 5 ] = pow( cos( jetTheta ) , 2 );
		 jetCovMatrixConf[ 6 ] = ( jetMom / jetEnergy ) * sin( jetTheta ) * cos( jetPhi );
		 jetCovMatrixConf[ 7 ] = ( jetMom / jetEnergy ) * sin( jetTheta ) * sin( jetPhi );
		 jetCovMatrixConf[ 8 ] = ( jetMom / jetEnergy ) * cos( jetTheta );
		 jetCovMatrixConf[ 9 ] = pow( jetMom , 2 ) / pow( jetEnergy , 2 );
		 
		 // Get PFOs in the jet and number of them
		 ReconstructedParticleVec jetPFOs  = jetPtr->getParticles();
		 nPFOs = jetPFOs.size();
		 streamlog_out(DEBUG) << "Number of PFOs in the jet:" << nPFOs << std::endl;

		 // Create a vector to hold covariance matrix of the jet
		 std::vector< float > jetCovMatrix( 10, 0.0 );

		 // Loop over all PFOs in the jet
		 for ( size_t iPFO = 0; iPFO < nPFOs; ++iPFO ) {
			// Set a pointer to the i-th PFO
			ReconstructedParticle *particlePtr = jetPFOs[ iPFO ];
			LCObjectVec vecMCTL = navMCTL->getRelatedToObjects( particlePtr );
			int nTrackspfo = (particlePtr->getTracks()).size();

			// Add particle energy to the sum
			eJetTotal += particlePtr->getEnergy();
			double scaleFactor = 1.0;

			if ( p_useFullCovMatNeut )
			{
				FloatVec particleCovMatrix = particlePtr->getCovMatrix();
				if ( 0 != nTrackspfo )
				{
					++numChargedPFOs;
					eChargedPFOs += particlePtr->getEnergy();
					scaleFactor = 1.0;
				}
				// particle is not charged
				// check whether its a photon
				else if ( 22 == particlePtr->getType() )
				{
					++numPhotons;
					ePhotons += particlePtr->getEnergy();
					scaleFactor = p_CovMatFactorPhotons;
				}
				else
				{   // if not a photon, then it's a neutral hadron
					++numNeutralPFOs;
					eNeutralPFOs += particlePtr->getEnergy();
					scaleFactor = p_CovMatFactorNeutralHadrons;
				} // end if-else
				for ( size_t iElement = 0; iElement < 9; ++iElement)
				{ // for all CovMat elements except energy error ( CovMat(9) := sigma_E^2 )
					if ( iElement < 6 )
					{ // for elements not involving E ( sigma_px^2 , sigma_pxpy , sigma_py^2 , sigma_pxpz , sigma_pypz , sigma_pz^2 )
						jetCovMatrix[ iElement ] += scaleFactor * scaleFactor * particleCovMatrix [ iElement ];
					}
					else
					{ // for elements involving E ( sigma_pxE , sigma_pyE , sigma_pzE )
						jetCovMatrix[ iElement ] += scaleFactor * particleCovMatrix [ iElement ];
					}
				} // end for
				// for energy error ( CovMat(9) := sigma_E^2 )
				jetCovMatrix[ 9 ] += particleCovMatrix [ 9 ];
			}
			else
			{

				// If it is a charged particle, add its covaricance matrix to jet cov matirix
				if ( 0 != nTrackspfo ) {
				   ++numChargedPFOs;
				   eChargedPFOs += particlePtr->getEnergy();
				   FloatVec particleCovMatrix = particlePtr->getCovMatrix();
				   for ( size_t iElement = 0; iElement < 10; ++iElement) {
					  jetCovMatrix[ iElement ] += particleCovMatrix [ iElement ];
				   } // end for
				} else {	// particle is not charged
				   // check whether its a photon
				   if ( 22 == particlePtr->getType() ) {
					  ++numPhotons;
					  ePhotons += particlePtr->getEnergy();
					  jetCovMatrix[ 9 ] += getPhotonSigmaESqr( particlePtr->getEnergy() );
				   } else {   // if not a photon, then it's a neutral hadron
					  ++numNeutralPFOs;
					  eNeutralPFOs += particlePtr->getEnergy();
					  jetCovMatrix[ 9 ] += getNeuHadSigmaESqr( particlePtr->getEnergy() );
				   } // end if-else
				} // end if-else charged
			}

			 /* :TODO:03/28/2016 07:09:24 PM:: Double check everything in this procedure here */
			// Compute error for semi-leptonic decays
			// If the PFO is a lepton add its energy to Total Lepton energy
			if( 0 < vecMCTL.size() ) {
			   MCParticle *mcParticle = dynamic_cast<MCParticle*>( vecMCTL[ 0 ] );
			   // Get PDG of the particle
			   int mcpPDG = mcParticle->getPDG();

			   // If it is an electron or a muon add its energy to total lepton energy
			   if ( 11 == abs( mcpPDG ) || 13 == abs( mcpPDG ) ) {
				   /* :TODO:04/03/2016 07:01:21 PM:: Why using mcParticle for the momentum? */
				  double particleMomentum = getTotalMomentum( mcParticle->getMomentum() ); 
				  if ( particleMomentum > 3.0 ) {
					 eLeptonsTotal += particlePtr->getEnergy();
				  } // end if
			   } // end if an electron or a muon
			} // end if vecMCTL.size()

		 } //end for loop over all PFOs of a jet

		 // Compute confusion terms of the current jet
		 double relConfTerms[ 3 ] = { 0.0 };                /* Confusion terms: 0: charged, 1: photons, 2: neutral */
		 getRelativeConfusion( eJetTotal, eChargedPFOs, ePhotons, eNeutralPFOs, relConfTerms );
		 // Assign each confusion term to its own variable, just to work around a problem with ROOT
		 absDetResSquared = jetCovMatrix [ 9 ];
		 relConfCharged = relConfTerms[ 0 ];
		 relConfPhotons = relConfTerms[ 1 ];
		 relConfNeutral = relConfTerms[ 2 ];

		 // Compute absolute confusion resolution squared
		 relConfSquared = pow( relConfCharged, 2.0 ) +
			pow( relConfPhotons, 2.0 ) + pow( relConfNeutral, 2.0 );
		 absConfSquared = relConfSquared * eJetTotal * eJetTotal;

		 // Compute absolute semi-leptonic resolution squared
		 absSemiLepResSquared = p_semiLepSigmaCorrFactor * p_semiLepSigmaCorrFactor
							    * eLeptonsTotal * eLeptonsTotal;

		 // Sum up all uncertainties to get total jet energy uncertainty
		 if ( p_confusionterm && p_propagateConfusiontoMomentumComp )
		 {
			for ( int iElement = 0 ; iElement < 10 ; ++iElement )
			{
				jetCovMatrix[ iElement ] = jetCovMatrix[ iElement ]
					+ ( pow( jetEnergy , 2 ) * absConfSquared * p_scaleConf * p_scaleConf / pow( jetMom , 2 ) ) * jetCovMatrixConf[ iElement ];
			}
		 }
		 else
		 {
			 jetCovMatrix[ 9 ] = jetCovMatrix[ 9 ]
				+ ( p_confusionterm ? absConfSquared * p_scaleConf * p_scaleConf : 0.0 )
				+ ( p_semiLepCorrection ? absSemiLepResSquared : 0.0 );
		}

		 // Fill the ROOT tree if enabled in the steering file
		 if ( p_storeTree ) {
			tree->Fill();
		 } // end if p_storeTree

		 // Update the convariance matrix of the jet object
		 jetPtr->setCovMatrix( jetCovMatrix );
		 outCol->addElement( jetPtr );

	  } // end for loop over all jets

	  // Update output collection with new collection of jets which has convariance matrix
	  evt->addCollection( outCol, p_outputJetCollection );

   } else {
	  streamlog_out(DEBUG) << "Error: collection does not exist!" << std::endl;
   } //end if-else collection exists

   // Collect the garbage
   delete navMCTL;

   ++p_nEvt;

} /* -----  end of function ErrorFlow::processEvent  ----- */



/*
 * ===  FUNCTION  ======================================================================
 *         Name:  ErrorFlow::resetVariables
 *  Description:  This function resets all variables used in processEvent()
 *                                function to zero.
 * =====================================================================================
 */
void ErrorFlow::resetVariables ()
{
   // Reset PFO counters and energies
   numChargedPFOs = 0;
   numPhotons = 0;
   numNeutralPFOs = 0;
   nPFOs = 0;
   ePhotons = 0.0;
   eChargedPFOs = 0.0;
   eNeutralPFOs = 0.0;
   eJetTotal = 0.0;
   relConfCharged = 0.0;
   relConfPhotons = 0.0;
   relConfNeutral = 0.0;
   relConfSquared = 0.0;
   absConfSquared = 0.0;
   absDetResSquared = 0.0;
   eLeptonsTotal = 0.0;
   absSemiLepResSquared = 0.0;

}               /* -----  end of function ErrorFlow::resetVariables  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getConfusion
 *  Description:  This fucntion computes RELATIVE confusion terms for 3 groups of particles
 *                in each jet based on its particle content. Then total confusion of the
 *                would be these 3 terms added in quadrature.
 *                The assumption for computation is that confusion terms given in 
 *                table 4 of DOI:10.1016/j.nima.2009.09.009 are for a jet which has
 *                0.62 of its energy in charged particles, 0.28 in photons and 0.10 in
 *                neutral hadrons. Then fracton of energy of each of these three categories
 *                is computed for each jet and confusion terms are scaled accordingly.
 * =====================================================================================
 */
double *ErrorFlow::getRelativeConfusion ( double t_jetEnergy, double t_chargedHadronsEnergy,
	                             double t_photonsEnergy, double t_neutralHadronsEnergy, double t_confusion[] )
{
   // Define and initialise variables
   // Average jet energy carried by: Ref. DOI:10.1016/j.nima.2009.09.009
   const double muEnCharged = 0.62;            /* charged particles */
   /* :TODO:01/24/2016 06:29:03 PM:: What to do with the neutrinos? */
   const double muEnPhotons = 0.28;            /* Photons +1% for neutrinos*/
   const double muEnNeutral = 0.10;            /* neutral particles */

   // Average confusion term for {charged, photons, neutral} for various energies
   // Ref. DOI:10.1016/j.nima.2009.09.009, values are from table 5 multipied by 1.1
   // to get gaussian values, as suggested in section 5 of the same paper.
   //  /* :TODO:01/24/2016 07:02:23 PM:: Is the multiplication actually correct? */
   const double conf45[ 3 ]  = { 1.32, 0.88, 0.99 }; /* Jet energy 45  GeV */
   const double conf100[ 3 ] = { 0.77, 1.10, 1.43 }; /* Jet energy 100 GeV */
   const double conf180[ 3 ] = { 0.55, 1.21, 1.87 }; /* Jet energy 180 GeV */
   const double conf250[ 3 ] = { 0.22, 1.43, 1.98 }; /* Jet energy 250 GeV */
	  
   double gammaChargedHadrons = ( t_chargedHadronsEnergy / t_jetEnergy ) / muEnCharged;
   double gammaPhotons = ( t_photonsEnergy / t_jetEnergy ) / muEnPhotons;
   double gammaNeutalHadrons = ( t_neutralHadronsEnergy / t_jetEnergy ) / muEnNeutral;

   // Compute confusion terms for each group based on the jet energy
   // Formula used for interpolation: Y = ( ( X - X1 )( Y2 - Y1) / ( X2 - X1) ) + Y1
   if ( 45.0 >= t_jetEnergy ) {
	  t_confusion[ 0 ] = conf45[ 0 ] * gammaChargedHadrons / 100.0;
	  t_confusion[ 1 ] = conf45[ 1 ] * gammaPhotons / 100.0;
	  t_confusion[ 2 ] = conf45[ 2 ] * gammaNeutalHadrons / 100.0;
   }
   else if ( 45 < t_jetEnergy && 100.0 >= t_jetEnergy ) {
	  t_confusion[ 0 ] = ( ( ( t_jetEnergy - 45.0 ) * ( conf100[ 0 ] - conf45[ 0 ] )
			   / ( 100.0 - 45.0 ) ) + conf45[ 0 ] ) * gammaChargedHadrons / 100.0;
	  t_confusion[ 1 ] = ( ( ( t_jetEnergy - 45.0 ) * ( conf100[ 1 ] - conf45[ 1 ] )
			   / ( 100.0 - 45.0 ) ) + conf45[ 1 ] ) * gammaPhotons / 100.0;
	  t_confusion[ 2 ] = ( ( ( t_jetEnergy - 45.0 ) * ( conf100[ 2 ] - conf45[ 2 ])
			   / ( 100.0 - 45.0 ) ) + conf45[ 2 ] ) * gammaNeutalHadrons / 100.0;
   }
   else if ( 100 < t_jetEnergy && 180.0 >= t_jetEnergy ) {
	  t_confusion[ 0 ] = ( ( ( t_jetEnergy - 100.0 ) * ( conf180[0] - conf100[0] )
			   / ( 180.0 - 100.0 ) ) + conf100[0] ) * gammaChargedHadrons / 100.0;
	  t_confusion[ 1 ] = ( ( ( t_jetEnergy - 100.0 ) * ( conf180[1] - conf100[1] )
			   / ( 180.0 - 100.0 ) ) + conf100[1] ) * gammaPhotons / 100.0;
	  t_confusion[ 2 ] = ( ( ( t_jetEnergy - 100.0 ) * ( conf180[2] - conf100[2] )
			   / ( 180.0 - 100.0 ) ) + conf100[2] ) * gammaNeutalHadrons / 100.0;
   }
   else {
	  t_confusion[ 0 ] = ( ( ( t_jetEnergy - 180.0 ) * ( conf250[ 0 ] - conf180[ 0 ] )
			   / ( 250.0 - 180.0 ) ) + conf180[ 0 ] ) * gammaChargedHadrons / 100.0;
	  t_confusion[ 1 ] = ( ( ( t_jetEnergy - 180.0 ) * ( conf250[ 1 ] - conf180[ 1 ] )
			   / ( 250.0 - 180.0 ) ) + conf180[ 1 ] ) * gammaPhotons / 100.0;
	  t_confusion[ 2 ] = ( ( ( t_jetEnergy - 180.0 ) * ( conf250[ 2 ] - conf180[ 2 ] )
			   / ( 250.0 - 180.0 ) ) + conf180[ 2 ] ) * gammaNeutalHadrons / 100.0;
   }
   for (size_t i = 0; i < 3; ++i) {
	  streamlog_out(DEBUG) << "Confusion " << i << ":" << t_confusion[i] << std::endl;
   }

   return t_confusion;

}		/* -----  end of function getConfusionSqr  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getPhotonSigmaESqr
 *  Description:  This function computes and returns energy uncertainty SQUARED for
 *				  photons using energy resolution of electromagnetic calorimeter
 *				  The following formula is used:
 *				  sigmaE^2 = aECAL^2 . E + cECAL^2 . E^2
 * =====================================================================================
 */
double ErrorFlow::getPhotonSigmaESqr ( double t_phEnergy )
{
   double sigmaEnSquared = ( p_aECAL * p_aECAL ) * t_phEnergy +
	  ( p_cECAL * p_cECAL ) * ( t_phEnergy * t_phEnergy );

   return sigmaEnSquared;

}		/* -----  end of function ErrorFlow::getPhotonSigmaESqr  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getNeuHadSigmaESqr
 *  Description:  This function computes and returns energy uncertainty SQUARED for
 *				  neutral hadrons using energy resolution of hadron calorimeter
 *				  The following formula is used:
 *				  sigmaE^2 = aHCAL^2 . E + cHCAL^2 . E^2
 * =====================================================================================
 */
double ErrorFlow::getNeuHadSigmaESqr ( double t_neuHadEnergy )
{
   double sigmaEnSquared = ( p_aHCAL * p_aHCAL ) * t_neuHadEnergy +
	  ( p_cHCAL * p_cHCAL ) * ( t_neuHadEnergy * t_neuHadEnergy );

   return sigmaEnSquared;

}		/* -----  end of function ErrorFlow::getNeuHadSigmaESqr  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ErrorFlow::getTotalMomentum
 *  Description:  Computes and returns total momentum from a 3-momentum
 * =====================================================================================
 */
double ErrorFlow::getTotalMomentum ( const double * t_threeMomentum )
{
   float moX = t_threeMomentum[ 0 ];
   float moY = t_threeMomentum[ 1 ];
   float moZ = t_threeMomentum[ 2 ];
   double totalMomentum = sqrt( moX*moX + moY*moY +moZ*moZ );

   return totalMomentum;

}		/* -----  end of function ErrorFlow::getTotalMomentum  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ErrorFlow::check( LCRunHeader* run )
 *  Description:  
 * =====================================================================================
 */
void ErrorFlow::check( LCEvent * /*evt*/ ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
} /* -----  end of function ErrorFlow::check  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ErrorFlow::end()
 *  Description:  
 * =====================================================================================
 */
void ErrorFlow::end(){ 

    //   std::cout << "MyProcessor::end()  " << name() 
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;

} /* -----  end of function ErrorFlow::end  ----- */
