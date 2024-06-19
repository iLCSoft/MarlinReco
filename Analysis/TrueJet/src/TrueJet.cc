#include "TrueJet.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

#include <EVENT/ReconstructedParticle.h>
#include <EVENT/LCRelation.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/ParticleIDImpl.h>
#include <UTIL/PIDHandler.h>
#include <UTIL/LCRelationNavigator.h>

using namespace lcio ;
using namespace marlin ;

struct MCPpyjet : LCIntExtension<MCPpyjet> {} ;
struct TJindex : LCIntExtension<TJindex> {} ;
LCRelationNavigator* reltrue =0;

TrueJet aTrueJet ;

TrueJet::TrueJet() : Processor("TrueJet") {

    // modify processor description
    _description = "TrueJet does whatever it does ..." ;

    // register steering parameters: name, description, class-variable, default value
 

   // Inputs: MC-particles, Reco-particles, the link between the two
   
 registerInputCollection( LCIO::MCPARTICLE,
                           "MCParticleCollection" , 
                           "Name of the MCParticle collection"  ,
                            _MCParticleColllectionName ,
                           std::string("MCParticlesSkimmed") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "RecoParticleCollection" ,
                           "Name of the ReconstructedParticles input collection"  ,
                           _recoParticleCollectionName ,
                           std::string("PandoraPFOs") ) ;

  registerInputCollection( LCIO::LCRELATION,
                           "RecoMCTruthLink",
                           "Name of the RecoMCTruthLink input collection"  ,
                           _recoMCTruthLink,
                           std::string("RecoMCTruthLink") ) ;

   // Outputs: True jets (as a recoparticle, will be the sum of the _reconstructed particles_
   // created by the true particles in each true jet, in the RecoMCTruthLink sense.
   // link jet-to-reco particles, link jet-to-MC-particles.

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "TrueJets" ,
                           "Name of the TrueJetCollection output collection"  ,
                           _trueJetCollectionName ,
                           std::string("TrueJets") ) ;

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "FinalColourNeutrals" ,
                           "Name of the FinalColourNeutralCollection output collection"  ,
                           _finalColourNeutralCollectionName ,
                           std::string("FinalColourNeutrals") ) ;

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "InitialColourNeutrals" ,
                           "Name of the InitialColourNeutralCollection output collection"  ,
                           _initialColourNeutralCollectionName ,
                           std::string("InitialColourNeutrals") ) ;


  registerOutputCollection( LCIO::LCRELATION,
                            "TrueJetPFOLink" , 
                            "Name of the TrueJetPFOLink output collection"  ,
                            _trueJetPFOLink,
                            std::string("TrueJetPFOLink") ) ;

  registerOutputCollection( LCIO::LCRELATION,
                            "TrueJetMCParticleLink" , 
                            "Name of the TrueJetMCParticleLink output collection"  ,
                            _trueJetMCParticleLink,
                            std::string("TrueJetMCParticleLink") ) ;

  registerOutputCollection( LCIO::LCRELATION,
                            "FinalElementonLink" , 
                            "Name of the  FinalElementonLink output collection"  ,
                            _finalElementonLink,
                            std::string("FinalElementonLink") ) ;

  registerOutputCollection( LCIO::LCRELATION,
                            "InitialElementonLink" , 
                            "Name of the  InitialElementonLink output collection"  ,
                            _initialElementonLink,
                            std::string("InitialElementonLink") ) ;

  registerOutputCollection( LCIO::LCRELATION,
                            "FinalColourNeutralLink" , 
                            "Name of the  FinalColourNeutralLink output collection"  ,
                            _finalColourNeutralLink,
                            std::string("FinalColourNeutralLink") ) ;

  registerOutputCollection( LCIO::LCRELATION,
                            "InitialColourNeutralLink" , 
                            "Name of the  InitialColourNeutralLink output collection"  ,
                            _initialColourNeutralLink,
                            std::string("InitialColourNeutralLink") ) ;
  // steering

  registerProcessorParameter("Whizard1",
                             "true: Input has MCParticles in Whizard1 convention ",
                             _whiz1,
                             bool(false)
                             );
 }



void TrueJet::init() { 

    streamlog_out(DEBUG7) << "   init called  " << std::endl ;

    // usually a good idea to
    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}


void TrueJet::processRunHeader( LCRunHeader*  /*run*/) { 

    _nRun++ ;
} 



void TrueJet::processEvent( LCEvent * event ) { 


    // For the concepts this method works with, see the comments in TrueJet.h

    evt=event;
    streamlog_out(WARNING) << " processing event: " << evt->getEventNumber() 
        << "   in run:  " << evt->getRunNumber() << std::endl ;
    streamlog_out(MESSAGE4) << " =====================================" << std::endl;


    LCCollection* mcpcol = NULL;
    try{
        mcpcol = evt->getCollection( _MCParticleColllectionName );
    }
    catch( lcio::DataNotAvailableException& e )
    {
        streamlog_out(WARNING) <<  _MCParticleColllectionName  << " collection not available" << std::endl;
        mcpcol = NULL;
	return;
    }
    LCCollection* rmclcol = NULL;
    try{
        rmclcol = evt->getCollection( _recoMCTruthLink );
    }
    catch( lcio::DataNotAvailableException& e )
    {
        streamlog_out(MESSAGE4) << _recoMCTruthLink   << " collection not available" << std::endl;
        rmclcol = NULL;
    }

    // TODO: in w2 : get polarization1 &2 ; count sum of W:s and L/R : if this is odd, then an odd number
    //   of leptons are expected ( one beam-remnant + zero or more pairs ). For now, just  skip the
    //   warning about odd number of leptons.

    if( mcpcol != NULL ){

      if ( rmclcol != NULL ) { reltrue = new LCRelationNavigator( rmclcol ); }

      // recast the MCParticles into a PYJETS look-alike, needed for the logic to work:

 
      mcp_pyjets.reserve(4000);
      getPyjets(mcpcol );
      
      if  ( _higgs_to_glue_glue ) {
        streamlog_out(DEBUG4) << " Higgs-> gluon gluon event : TrueJet will not work correctly " << std::endl;
        if ( reltrue ) delete reltrue ;
        _nEvt++ ;
        LCCollectionVec* jet_vec = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE ) ;
        evt->addCollection( jet_vec,  _trueJetCollectionName);
        return;	
      }	



      njet=0 ; current_jet=0 ; nstr=0 ;
      for ( int kk=1 ; kk <= 25 ; kk++ ) {
        boson[kk]=0;  
        elementon[kk]=0;
        fafp_last[kk]=0;
        fafp_boson[kk]=0;
        fafp[kk]=0;
        nfsr[kk]=0 ;
      }
      bool seen[4001];
      for ( int kk=1 ; kk <= 4000 ; kk++ ) {
        jet[kk]=0;
        companion[kk] = 0;
	    seen[kk] = false;
      }

        // the actual work happens in the following routines. They are each responible to 
        // find jets stemming from different sources, as indicated by their names.

      cluster() ;
    
      true_lepton() ;

      photon();
      
      isr() ;

      nstr=0;
      for ( int k_py=1 ; k_py <= nlund ; k_py++ ) {
        streamlog_out(DEBUG0) << " Checking string at " << k_py << " k[k_py][2], jet[k_py] and k[k_py][1] : " <<   
          " " << k[k_py][2] <<  " , " << jet[k_py]  <<    " , " << k[k_py][1] <<                        std::endl;
        if ( k[k_py][2] == 92 && jet[k_py] == 0 && k[k_py][1] < 30  ) {
          streamlog_out(DEBUG4) << " string found at " << k_py <<std::endl;
          nstr++;
        }
      }

      njet = 2*nstr+n_hard_lepton+2*nclu+nisr+nphot;

      if ( nstr != 0 ) {
        string();
      } 

      n_beam_jet=0;
      for (int k_py=1 ; k_py<=nlund ; k_py++) {
        if ( k[k_py][1] >= 30 ) {
          if ( n_beam_jet == 0 ) {
            n_beam_jet=1;
          }
          jet[k_py]=njet+n_beam_jet;
        }
      }

      njet=njet+n_beam_jet;

      grouping();

      // At this point, the job is done. All the rest is printouts, putting things into collections and setting up navigators.
      // Non-debuggung code is beween ******:s and =========:s

      // FIXME : determine form the event header is an odd or even number of leptons are expected
      //      if ( ! ((n_mixed==0) && njet==2*nstr+n_hard_lepton+2*nclu+nisr+n_beam_jet &&
      //	      n_jetless==0 && n_hard_lepton%2 == 0 )) {

      if ( ! ((n_mixed==0) && njet==2*nstr+n_hard_lepton+2*nclu+nisr+n_beam_jet+nphot &&
	      n_jetless==0 )) {
        if ( ! _higgs_to_glue_glue ) {
          streamlog_out(ERROR) << " inconsiency in jet finding : " << std::endl;
          streamlog_out(ERROR) << " n_mixed/ njet/ 2*nstr+n_hard_lepton+2*nclu+nisr+n_beam_jet / n_jetless / n_hard_lepton: " << std::endl;
          streamlog_out(ERROR) << n_mixed <<" / " << njet<< " / " <<2*nstr+n_hard_lepton+2*nclu+nisr+n_beam_jet<< " / " 
                               << n_jetless << " / " << n_hard_lepton << std::endl;
          streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
          streamlog_out(ERROR) << std::endl ;
        }  
      }
 
      streamlog_out(DEBUG3) << " HEPEVT relation table up with jet assignment "  <<std::endl;
      streamlog_out(DEBUG3) << "    line   status       pdg  parent  first  last      px         py   "
                              << "       pz          E             M        jet companion type" << std::endl;
      streamlog_out(DEBUG3) << "                                       daughter                       "
                              << "                                                 jet" << std::endl;

      for ( int k_py=1 ; k_py <= nlund ; k_py++){
	streamlog_out(DEBUG3) << std::setw(7) << k_py << std::setw(7) <<  k[k_py][1] << 
                 std::setw(12) <<  k[k_py][2] << std::setw(7) <<  k[k_py][3] << 
                 std::setw(7) <<  k[k_py][4] << std::setw(7) << k[k_py][5] <<
                 std::setw(12) <<  p[k_py][1] << 
                 std::setw(12) <<  p[k_py][2] << std::setw(12) <<  p[k_py][3] << 
                 std::setw(12) <<  p[k_py][4] << std::setw(12) << p[k_py][5] <<
 
                 std::setw(7) << jet[k_py]  <<
                 std::setw(7) << companion[abs(jet[k_py])] << std::setw(7) <<  type[abs(jet[k_py])]<<std::endl;
      }
      streamlog_out(DEBUG3) << "       jet     fafp-beg  qrk/lept  fafp-end fafp-boson  type   dijet-beg  dijet-end    boson" << std::endl;
      for ( int k_py=1 ; k_py <= njet ; k_py++){
	streamlog_out(DEBUG3) << std::setw(10) << k_py << std::setw(10) <<  fafp[k_py] << 
                 std::setw(10) <<  elementon[k_py] << std::setw(10) <<  fafp_last[k_py] << 
                 std::setw(10) <<  fafp_boson[k_py] << std::setw(10) << type[k_py] << std::setw(10) << dijet_begining[k_py]  <<
                 std::setw(10) << dijet_end[k_py] << std::setw(10) <<  boson[k_py]<<std::endl;
      }


      // now for py-partic i I know which true jet it belongs to.
      // mcp_pyjets tells me which MCParticle corresponds to each pyjets one.
      // I then loop all jets creating the true_jet object for each one, loop 
      // all mcparticles to find which ones
      // belong to the current jet, find if they contribute to a reco-particle, and
      // if so, add the reco energy etc. to the true_jet-object. While loopoing, the
      // links are set up, as well.
  
      //*****************************
      // create the navigators:

      LCCollection* tjrcol = 0;
      LCCollection* tjtcol = 0;
      LCRelationNavigator truejet_pfo_Nav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::RECONSTRUCTEDPARTICLE ) ;
      LCRelationNavigator truejet_truepart_Nav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE ) ;

  
  
      // create the collection-vector that will contain the true-jet objects:

      LCCollectionVec* jet_vec = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE ) ;
      auto jetPidHandler = UTIL::PIDHandler(jet_vec);
      const auto truePidID = jetPidHandler.addAlgorithm("TrueJetPID", {});

      //============================

      double str_tmom[3]={0}, str_mom[3]={0}, str_tmomS[3]={0} , str_tE=0 , str_E=0 , str_tES=0; 

      streamlog_out(DEBUG8) << "  Number of jets found : " << njet << std::endl;

      double tmomS[25][3]={0.} ;
      double tES[25]={0.};
      int pid_type[25]={0} ;
      for ( int i_jet=1; i_jet<=njet ; i_jet++ ) { // jet-loop

       //*****************************
         ReconstructedParticleImpl* true_jet = new ReconstructedParticleImpl;
      
        pid_type[i_jet] = -type[i_jet] ;
 
        jet_vec->addElement(true_jet);
       //============================
	
	streamlog_out(DEBUG1) << std::endl;
	streamlog_out(DEBUG1) << "  Following jet " << i_jet <<   " ( "<< true_jet << ")" << std::endl;
        streamlog_out(DEBUG1) << " ============================= "  <<  std::endl;
        streamlog_out(DEBUG1) << " HEPEVT relation table with jet assignment and pfo(s), if any "  <<std::endl;
        streamlog_out(DEBUG1) << "       line    mcp          status    pdg       parent    first     last       jet     pfo(s)" << std::endl;
        streamlog_out(DEBUG1) << "                                                        daughter   daugther " << std::endl;
      }
      
      ReconstructedParticleImpl* true_jet ;
      for ( int k_py=1 ; k_py<=nlund ; k_py++){ // pyjets loop
        int i_jet =  abs(jet[k_py]);
        if ( i_jet > 0 ) {

          true_jet = dynamic_cast<ReconstructedParticleImpl*>(jet_vec->getElementAt(i_jet-1));

             // OK, add jet-to-true link

	  streamlog_out(DEBUG1) << std::setw(10) << k_py << " ("  <<mcp_pyjets[k_py] << ")" << std::setw(10) <<  k[k_py][1]%30 << 
                 std::setw(10) <<  k[k_py][2] << std::setw(10) << k[k_py][3] << 
                 std::setw(10) <<  k[k_py][4] << std::setw(10) << k[k_py][5] << std::setw(10) << jet[k_py]  ;


          //*****************************
          if (  k[k_py][1] == 1 &&  k[k_py][2] == 22 && jet[k[k_py][3]] < 0   && k[k[k_py][3]][2] != 22 ) {
                // to be able to find FSRs, set weight to -ve for them
            truejet_truepart_Nav.addRelation(  true_jet, mcp_pyjets[k_py], ( jet[k_py]>0 ? -(k[k_py][1])%30 : 0.0 )) ;
          } else {
            truejet_truepart_Nav.addRelation(  true_jet, mcp_pyjets[k_py], ( jet[k_py]>0 ? k[k_py][1]%30 : 0.0 )) ;
          }
          LCObjectVec recovec; 
	  if ( rmclcol != NULL ) {
            recovec = reltrue->getRelatedFromObjects( mcp_pyjets[k_py]);
	  }  
          //============================ 

          // (maybe this will be needed: If the current assumption that also MCPartiles created in
          //  simulation can in a consitent way be put into pyjets turns out to be wrong, then
          //  getPyjets should only load gen. stat. /= 0 particles, and then code below
          //  would be needed to assign seen gen. stat. = 0 particles to jets based on their
          //  gen. stat. /= 0 ancestor ?!)
          //            if ( recovec.size() == 0 &&  k[k_py][1] == 1 ) { // not reconstructed stable particle
          //              if ( abs(k[k_py][2]) != 12 &&  abs(k[k_py][2]) != 14 &&  abs(k[k_py][2]) != 16 ) { // not neutrino
          //                //  mcp_pyjets[k_py].getDaughter ... etc . If any decendant is seen, make association
          //              }
          //            }

          if ( recovec.size() > 0 ) { // if reconstructed
            double mom2[3] ;  
            double E,q ;

             // true quanaties for this jet (all true seen)

	    if ( ! seen[k[k_py][3]] ) {
	       // add to true-of-seen only if no ancestor of this true has already been counted
              tmomS[i_jet][0]+=p[k_py][1] ;
              tmomS[i_jet][1]+=p[k_py][2] ;
              tmomS[i_jet][2]+=p[k_py][3] ;
              tES[i_jet]+=p[k_py][4] ;
            } else {
              streamlog_out(DEBUG3) << " Ancestor of true particle " << k_py << " already seen, not added again to true-of-seen "  << std::endl;
            }
	    seen[k_py] = true ;
	    streamlog_out(DEBUG1) << " recovec size is " << recovec.size() ;
            int last_winner = i_jet ;
            for ( unsigned i_reco=0 ; i_reco<recovec.size() ; i_reco++ ) { //  reco-of-this-true loop

              // add things of this reco-particle to the current jet:

              //*****************************
              ReconstructedParticle* reco_part  = dynamic_cast<ReconstructedParticle*>(recovec[i_reco]);

              if (truejet_pfo_Nav.getRelatedFromObjects(reco_part).size() == 0 ) { // only if not yet used
	        streamlog_out(DEBUG3) << " recopart " << i_reco ;
		bool split_between_jets = false;
	        int winner=i_jet , wgt_trk[26]={0} , wgt_clu[26]={0} ;

         	LCObjectVec recomctrues = reltrue->getRelatedToObjects(reco_part);
                static FloatVec recomctrueweights;
	        recomctrueweights = reltrue->getRelatedToWeights(reco_part);
	        streamlog_out(DEBUG3) << "     mctrues of this reco " << recomctrues.size() << std::endl ;
         	for ( unsigned k_mcp_of_reco=0 ; k_mcp_of_reco<recomctrues.size() ; k_mcp_of_reco++ ) {
		  MCParticle* an_mcp = dynamic_cast<MCParticle*>(recomctrues[k_mcp_of_reco]);
		  int jetoftrue=jet[an_mcp->ext<MCPpyjet>()];
	          if ( jetoftrue != i_jet ) split_between_jets = true;
		  wgt_trk[jetoftrue]+=int(recomctrueweights[k_mcp_of_reco])%10000 ;
		  wgt_clu[jetoftrue]+=int(recomctrueweights[k_mcp_of_reco])/10000 ;
		  streamlog_out(DEBUG1) << " weights for jet " << jetoftrue << " : " << int(recomctrueweights[k_mcp_of_reco]) 
                                          << " " << int(recomctrueweights[k_mcp_of_reco])%10000 
					  << " " << int(recomctrueweights[k_mcp_of_reco])/10000 << std::endl ;
                }
                if ( split_between_jets ) {
		  double wgt_trk_max=0, wgt_clu_max=0 ;
	          int imax_trk=0, imax_clu=0;

		  for ( int j_jet=1 ; j_jet <= njet ; j_jet++ ) {
		    streamlog_out(DEBUG3) << "     jetweights for " << j_jet << " " << wgt_trk[j_jet]  << " " <<  wgt_clu[j_jet] ;
		    if ( wgt_trk[j_jet] > wgt_trk_max ) {
	              imax_trk= j_jet ; wgt_trk_max = wgt_trk[j_jet] ;
         	    }  
		    if ( wgt_clu[j_jet] > wgt_clu_max ) {
		      imax_clu= j_jet ; wgt_clu_max = wgt_clu[j_jet] ;
		    }  
		  }
		  streamlog_out(DEBUG3) << "    " << imax_clu << " " << imax_trk << " " <<  wgt_clu_max <<  " " <<  wgt_trk_max ;
		  if ( imax_clu == imax_trk || wgt_trk_max >= 750 ) {
		    winner = imax_trk;
		  } else if ( wgt_trk_max == 0 ) {
		    winner = imax_clu;
		  } else {
		    streamlog_out(DEBUG3) << " was ? " << imax_clu << " " << imax_trk << " " <<  wgt_clu_max <<  " " <<  wgt_trk_max << std::endl;
		    for ( int j_jet=1 ; j_jet <= njet ; j_jet++ ) {
		      streamlog_out(DEBUG3) << " was    jetweights for " << j_jet << " " << wgt_trk[j_jet]  << " " <<  wgt_clu[j_jet]  << std::endl;
                    }
                    // probably the cluster is somewhat more reliable in this case. Anyways, this happens very rarely (for 3 PFOs in 6400 4f_had...)
                    if ( imax_clu > 0 ) {
                      winner = imax_clu;
                    } else {
                      // ... as a last resort ...
                      winner = imax_trk ;
                    }
                  }
		}

                if ( winner > 0 && winner != last_winner ) true_jet = dynamic_cast<ReconstructedParticleImpl*>(jet_vec->getElementAt(winner-1));

		streamlog_out(DEBUG3) << " and the winner is " << winner << "  ( i_jet is " << i_jet << " ) " ;	  
                if ( winner > 0 ) {
                  last_winner = winner ;

                  const double* mom =  reco_part->getMomentum(); 
                  double psq=0.;
                  for ( int xyz=0 ; xyz<3 ; xyz++) {
                    mom2[xyz]=mom[xyz]+true_jet->getMomentum()[xyz] ;
                    psq+= mom2[xyz]*mom2[xyz];
                  }
  
                  E =  reco_part->getEnergy()+true_jet->getEnergy();
                  q =  reco_part->getCharge()+true_jet->getCharge();
                  true_jet->setCharge(q); 
                  true_jet->setEnergy(E);
                  true_jet->setMass(sqrt(E*E-psq)); 
                  true_jet->setMomentum(mom2);
        	  streamlog_out(DEBUG1) << " " << reco_part ;
                  pid_type[winner] = type[winner] ;
  
                  true_jet->addParticle(reco_part);

                    // add jet-to-reco (and v.v.) link (good to have both ways easily, reco-particle member of the
                    // true jet object just gives jet->particle, not particle->jet!)
  
                  truejet_pfo_Nav.addRelation(  true_jet, reco_part , 1.0 );
                    //============================ 
	      
                } else {
                  // FIXME: a follow on error for h->glueglue: If a particle from a non-higgs jet gets reconstructed into a PFO where
                  //   the majority contribution is from a higgs jet (which are not found in  h->glueglue, and hence has jet# 0 ).
 		  if ( ! _higgs_to_glue_glue ) {
		    streamlog_out(WARNING) << " reco particle " << i_reco << " of pythia-particle " << k_py 
                                             << " has weights = 0 to ALL MCPs ??? " << std::endl;
                    streamlog_out(WARNING) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ;
		  }
                }
              } // end if not yet used
              else {
		streamlog_out(DEBUG1) << " " << reco_part << " seen more than once ! " << std::endl;
              }
            } // end reco-of-this-true loop

	  } else { // end if reconstructed
	    if ( seen[k[k_py][3]] ) { 
              seen[k_py] = true ;  // this makes sure that in any decay-chain, only the first seen particle will includes in true-of-seen
              streamlog_out(DEBUG2) << " Ancestor of true particle " << k_py << " was seen, so particle " << k_py 
                                    << " is marked as seen (even if it wasn't seen itself " << std::endl;
            }
          }
  	  streamlog_out(DEBUG1) << std::endl;

        } // end if i_jet > 0

      } // end pyjets loop

      for ( int i_jet=1; i_jet<=njet ; i_jet++ ) { // jet-loop

        true_jet = dynamic_cast<ReconstructedParticleImpl*>(jet_vec->getElementAt(i_jet-1));
        jetPidHandler.setParticleID(true_jet, pid_type[i_jet], k[elementon[i_jet]][2], 0, truePidID, {});
        true_jet->ext<TJindex>()=i_jet;

      }

      for ( int i_jet=1; i_jet<=njet ; i_jet++ ) { // jet-loop

        true_jet = dynamic_cast<ReconstructedParticleImpl*>(jet_vec->getElementAt(i_jet-1));

        const double* mom = true_jet->getMomentum(); 

	streamlog_out(DEBUG5) << std::endl;
	streamlog_out(DEBUG9) << " summary " << i_jet << " " << type[i_jet]%100 << " " << true_jet->getEnergy() << " " 
                                <<  tE[i_jet]  << " " <<  tES[i_jet] << std::endl ;
	streamlog_out(DEBUG5) << "  jet       " <<std::setw(4)<< i_jet << " (seen)         p : " 
                                << mom[0] << " " << mom[1] << " "  << mom[2] 
                                << " " << " E " << true_jet->getEnergy()<< std::endl;
	streamlog_out(DEBUG5) << "               "   <<      "  (true)         p : " 
                                << tmom[i_jet][0] << " " << tmom[i_jet][1] << " "  << tmom[i_jet][2] 
                                << " " << " E " << tE[i_jet] << std::endl;
	streamlog_out(DEBUG5) << "               "   <<      "  (true of seen) p : " 
                                << tmomS[i_jet][0] << " " << tmomS[i_jet][1] << " "  << tmomS[i_jet][2] 
                                << " " << " E " << tES[i_jet] << std::endl;
	streamlog_out(DEBUG5) << "  ancestor 1" << std::setw(4)<<elementon[i_jet]  << "                p : " 
                                << p[elementon[i_jet]][1] << " " <<  p[elementon[i_jet]][2] << " "  <<  p[elementon[i_jet]][3] 
	                        << " " << " E " << p[elementon[i_jet]][4]<< std::endl;
	streamlog_out(DEBUG5) << "  ancestor 2" << std::setw(4)<<fafp_last[i_jet] << "                p : " 
                                << p[fafp_last[i_jet]][1] << " " <<  p[fafp_last[i_jet]][2] << " "  <<  p[fafp_last[i_jet]][3] 
                                << " " << " E " << p[fafp_last[i_jet]][4]<< std::endl;
	streamlog_out(DEBUG5) << "  ancestor 3" << std::setw(4)<<fafp[i_jet] << "                p : " 
                                << p[fafp[i_jet]][1] << " " <<  p[fafp[i_jet]][2] << " "  <<  p[fafp[i_jet]][3] 
                                << " " << " E " << p[fafp[i_jet]][4]<< std::endl;
	streamlog_out(DEBUG5) << std::endl;

        streamlog_out(DEBUG5) << "   Character of jet " << i_jet << " : "<<type[i_jet]%100 <<" (1=string, 2=lepton , "
                                <<"3=cluster, 4=isr, 5=overlay, 6=M.E. photon)."  << std::endl;
        streamlog_out(DEBUG5) << "   Jet from boson ? " << (type[i_jet]>=100 ? type[i_jet]/100 : 0) 
                                << " (0 not from boson, !=0 from boson radiated of that jet)."  << std::endl;
        streamlog_out(DEBUG5) << "   No. of FSRs : " << nfsr[i_jet] << std::endl;
        streamlog_out(DEBUG5) << "   Jet has no energy ? " << (tE[i_jet]<0.001) << std::endl;

	streamlog_out(DEBUG5) << std::endl;

        str_E+= true_jet->getEnergy(); str_tE+= tE[i_jet]; str_tES+= tES[i_jet];
        for ( int i=0 ; i<3 ; i++ ) {
          str_mom[i]+=mom[i];
          str_tmom[i]+=tmom[i_jet][i];
          str_tmomS[i]+=tmomS[i_jet][i]; 
        }

        int odd_even = 0 ;
        if ( type[i_jet]%100 == 3 && n_hard_lepton%2 == 1 ) odd_even = 1;

        if ( i_jet%2 == odd_even && type[i_jet]%100 <= 3 ) {
          streamlog_out(DEBUG6) << "      Mass of jet " << i_jet-1 << " and " << i_jet << std::endl;
          double masssq = str_E*str_E - str_mom[0]*str_mom[0] - str_mom[1]*str_mom[1] - str_mom[2]*str_mom[2] ;
          if (  masssq > 0. ) {
             streamlog_out(DEBUG6) << "         Seen         " << sqrt(masssq)  << std::endl;
          }
          double masssq_t = str_tE*str_tE - str_tmom[0]*str_tmom[0] - str_tmom[1]*str_tmom[1] - str_tmom[2]*str_tmom[2] ;
          if (  masssq_t > 0. ) {
             streamlog_out(DEBUG6) << "         true         " << sqrt(masssq_t)  << std::endl;
          }
          masssq = str_tES*str_tES - str_tmomS[0]*str_tmomS[0] - str_tmomS[1]*str_tmomS[1] - str_tmomS[2]*str_tmomS[2] ;
          if (  masssq > 0. ) {
             streamlog_out(DEBUG6) << "         true of seen " << sqrt(masssq)  << std::endl;
          }
          if ( type[i_jet]%100 == 3  ) {
             // cluster - the 91 object does NOT have the mass of the following physical states,
             // so we sum up the descendants - in the vast majority there is only one
            double clu_mass=0 , clu_E= 0., clu_P[4]= {0,0,0};
            for ( int j_py=k[k[elementon[i_jet]][4]][4] ; j_py<=k[k[elementon[i_jet]][4]][5] ; j_py++ ) {
              if ( k[j_py][3] == k[elementon[i_jet]][4] )  {
                clu_E += p[j_py][4] ;    clu_P[1] += p[j_py][1] ;  clu_P[2] += p[j_py][2] ;  clu_P[3] += p[j_py][3] ;  
              }
            }
            clu_mass=sqrt(clu_E*clu_E-clu_P[1]*clu_P[1]-clu_P[2]*clu_P[2]-clu_P[3]*clu_P[3] );
            streamlog_out(DEBUG6) << "      string mass :   " <<    clu_mass << std::endl;
	  } else if ( type[i_jet]%100 == 2 ) {
            double dilept_mass=0 , dilept_E= 0., dilept_P[4]= {0,0,0};
            for ( int j_jet=i_jet-1 ; j_jet<=i_jet ; j_jet++ ) {
              dilept_E += p[elementon[j_jet]][4] ;    
              dilept_P[1] += p[elementon[j_jet]][1] ;   
              dilept_P[2] += p[elementon[j_jet]][2] ;   
              dilept_P[3] += p[elementon[j_jet]][3] ;  
            }
            dilept_mass=sqrt(dilept_E*dilept_E-dilept_P[1]*dilept_P[1]-dilept_P[2]*dilept_P[2]-dilept_P[3]*dilept_P[3] );
            streamlog_out(DEBUG6) << "      string mass :   " <<    dilept_mass << std::endl;
          } else {  
             // string, just take the number from the 92-object
            streamlog_out(DEBUG6) << "      string mass :   " <<    p[k[elementon[i_jet]][4]][5]  << std::endl;
          }
    
	  streamlog_out(DEBUG6) << std::endl;
          str_tE=0; str_tES=0; str_E=0;
          for ( int i=0 ; i<3 ; i++ ) {
            str_tmom[i]=0. ; str_tmomS[i]=0. ; str_mom[i]=0. ; 
          } 
          if (  masssq_t > 0. ) {
            if ( abs( sqrt(masssq_t)- p[k[elementon[i_jet]][4]][5] )/  p[k[elementon[i_jet]][4]][5]  > 0.007 ) {
              if ( type[i_jet]%100 == 1 && type[i_jet-1]%100 == 1 &&  nfsr[i_jet]+nfsr[i_jet-1]==0) {
                if (  abs( tE[i_jet]+tE[i_jet-1]- p[k[elementon[i_jet]][4]][4] )/  p[k[elementon[i_jet]][4]][4] > 0.001  ) {
	  	  streamlog_out(ERROR) << " bad match M (sum/initial) " << " " <<  sqrt(masssq_t) 
                                         << " / " << p[k[elementon[i_jet]][4]][5] << std::endl;
		  streamlog_out(ERROR) << "           E (sum/initial) " << " " <<  tE[i_jet]+tE[i_jet-1] 
                                         << " / " << p[k[elementon[i_jet]][4]][4] << std::endl;

                  streamlog_out(DEBUG9) << "list HEPEVT relation table up with jet assignment "  <<std::endl;
                  streamlog_out(DEBUG9) << "list    line   status       pdg  parent  first  last      px         py          pz  "
                                          <<"        E             M        jet companion type" << std::endl;
                  streamlog_out(DEBUG9) << "list                                       daughter                           "
                                          <<"                                             jet" << std::endl;
                  for ( int i_py=1 ; i_py <= nlund ; i_py++){
	            streamlog_out(DEBUG9) <<"list "<< std::setw(7) << i_py << std::setw(7) <<  k[i_py][1] << 
                      std::setw(12) <<  k[i_py][2] << std::setw(7) <<  k[i_py][3] << 
                      std::setw(7)  <<  k[i_py][4] << std::setw(7) << k[i_py][5] <<
                      std::setw(12) <<  p[i_py][1] << 
                      std::setw(12) <<  p[i_py][2] << std::setw(12) <<  p[i_py][3] << 
                      std::setw(12) <<  p[i_py][4] << std::setw(12) << p[i_py][5] <<
 
                      std::setw(7) << jet[i_py]  <<
                      std::setw(7) << companion[abs(jet[i_py])] << std::setw(7) <<  type[abs(jet[i_py])]<<std::endl;
                  }
                  streamlog_out(DEBUG9) << "list of individual particles with bad parent/kid energy "  <<std::endl;
                  for ( int i_py=1 ; i_py<=nlund ; i_py++ ) {
                    if ( jet[i_py] == i_jet || jet[i_py] == i_jet -1 ) {
                      if ( k[i_py][1] == 11 ) { 
                        double e_kid=0;
                        for ( int jj= k[i_py][4] ; jj <=  k[i_py][5] ; jj++ ) {
                          e_kid+=p[jj][4];
                        }
                        if ( (abs(e_kid - p[i_py][4])/ p[i_py][4] ) > 0.001 ) {
                           streamlog_out(DEBUG9) << i_py << " " <<  k[i_py][4]  << " " << k[i_py][5] << " "  
                                                   << e_kid << " " << p[i_py][4] << std::endl ;
                        }
                      }
                    }
                  }
                  streamlog_out(WARNING) << "       jet     fafp-beg  qrk/lept  fafp-end fafp-boson  "
                                            <<"type   dijet-beg  dijet-end    boson" << std::endl;
                  for ( int j_jet=i_jet-1 ; j_jet <= i_jet ; j_jet++){
	            streamlog_out(WARNING) << std::setw(10) << i_jet << std::setw(10) <<  fafp[j_jet] << 
                    std::setw(10) <<  elementon[j_jet] << std::setw(10) <<  fafp_last[j_jet] << 
                    std::setw(10) <<  fafp_boson[4] << std::setw(10) << type[j_jet] << std::setw(10) << dijet_begining[j_jet]  <<
	            std::setw(10) << dijet_end[j_jet] << std::setw(10) <<  boson[j_jet]<<std::endl;
                  }
                  streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                  streamlog_out(ERROR) << std::endl ; 
                } else {
		  streamlog_out(DEBUG8) << " bad match M (sum/initial) " << " " <<  sqrt(masssq_t) 
                                          << " / " << p[k[elementon[i_jet]][4]][5] << std::endl;
		  streamlog_out(DEBUG8) << "           E (sum/initial) " << " " <<  tE[i_jet]+tE[i_jet-1] 
                                          << " / " << p[k[elementon[i_jet]][4]][4] << std::endl;
                  streamlog_out(DEBUG8) << " (As the energy matches well, this is probably due to the B-field) "  <<std::endl;
                  streamlog_out(DEBUG8) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                  streamlog_out(DEBUG8) << std::endl;
                }
              }
            }
          }
        }
      } // end jet-loop



      //*****************************

        // ! fill the two color-singlet blocks
        // post-PS part

      LCCollectionVec* fafpf_vec = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE ) ;
      auto fafpfPidHandler = UTIL::PIDHandler(fafpf_vec);
      auto fafpf_mainPidId = fafpfPidHandler.addAlgorithm("TrueJet_fafpf", {});
      std::array<int, 2> fafpf_pidIds{};
      for (size_t i = 0; i < fafpf_pidIds.size(); ++i) {
          fafpf_pidIds[i] = fafpfPidHandler.addAlgorithm("TrueJet_fafpf_jet_" + std::to_string(i), {});
      }

      LCRelationNavigator FinalColourNeutral_Nav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::RECONSTRUCTEDPARTICLE ) ;
      LCRelationNavigator FinalElementon_Nav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE ) ;

      for ( int k_dj_end=1; k_dj_end<=n_dje ; k_dj_end++ ) {
        double E=0, M=0 , mom[3]={} ;  int  pdg[26]={};
        if ( type[jets_end[1][k_dj_end]]%100 == 5 ) {
          // beam-jet: No ancestor - just use sum of true-stable
          E=tE[jets_end[1][k_dj_end]];
          double psqrs=0;
          for (int jj=0 ; jj<3 ; jj++ ) {
            mom[jj] = tmom[jets_end[1][k_dj_end]][jj];
            psqrs += mom[jj]*mom[jj];
          }
          M = sqrt(E*E-psqrs);

	} else {

             // all other cases: sum of (unique) direct decendants of the
             // di-jet generator. For strings, clusters and ISR, there is only one
             // (a 92, a 91 or a 22,  respectively), so in these cases the sum is
             // actualy just copying the descendant. For leptons, however, there  are
             // two (the next generation leptons), so the sum is non-trivial.
             // (actually, sometimes also for clustes 0 see below)

          int first=nlund , last=0;
          for ( int j_jet_end=1 ; j_jet_end<=jets_end[0][k_dj_end] ; j_jet_end++ ) {

              // decide which line should be used - this (leptons) - child (strings,gammas) - grandchild (cluster)
              // NB. In the case of leptons, it might be that the first and last leptons are not next to
              // eachother.

            int this_first=0, this_last=nlund;
            if (  type[jets_end[1][k_dj_end]]%100 == 2  || type[jets_end[1][k_dj_end]]%100 == 4 || type[jets_end[1][k_dj_end]]%100 == 6 ) {

                 // lepton : need to go one step less as the di-jet generator might
                 // be stable, and hence have no descendants.

              this_first=elementon[jets_end[j_jet_end][k_dj_end]];
              this_last=elementon[jets_end[j_jet_end][k_dj_end]] ;

            } else if (  type[jets_end[1][k_dj_end]]%100 == 3 ) { 

                 // cluster: need to go one step further (the 91 object doesn't always 
                 // have the mass of the particle it goes to). Sometimes first and
                 // last will actually become different - the cluster goes to more than one
                 // hadron.

              this_first = k[k[elementon[jets_end[j_jet_end][k_dj_end]]][4]][4];
              this_last =  k[k[elementon[jets_end[j_jet_end][k_dj_end]]][5]][5];

            } else {

              if (k[elementon[jets_end[j_jet_end][k_dj_end]]][4] != 0 ) { 
                this_first = k[elementon[jets_end[j_jet_end][k_dj_end]]][4] ;
                this_last =  k[elementon[jets_end[j_jet_end][k_dj_end]]][5];
              } else {
                this_first = elementon[jets_end[j_jet_end][k_dj_end]] ;
                this_last =  elementon[jets_end[j_jet_end][k_dj_end]];
              }

            }

            if ( this_first < first ) {
              first= this_first;
            }
            if (  this_last > last ) {
              last= this_last;
            }

            pdg[j_jet_end]=k[elementon[jets_end[j_jet_end][k_dj_end]]][2];
            if ( pdg[0] == 0 ) {
              if ( k[elementon[jets_end[j_jet_end][k_dj_end]]][4] != 0 ) {
                pdg[0]=k[k[elementon[jets_end[j_jet_end][k_dj_end]]][4]][2];
              } else {
                pdg[0]=k[elementon[jets_end[j_jet_end][k_dj_end]]][2];
              }
            }
          }
          if ( last == first ) {
            E=  p[first][4];
            for ( int ll=1 ; ll<=3 ; ll++ ) {
              mom[ll-1] =    p[first][ll];
            }
          } else {
            // need double loop for non-consecutive case ( if this could only be F95 !)
            for ( int j_jet_end=1 ; j_jet_end<=jets_end[0][k_dj_end] ; j_jet_end++ ) {
              for ( int k_py=first ; k_py<=last ; k_py++ ) {
 
                if ( abs(jet[k_py]) == jets_end[j_jet_end][k_dj_end] ) {
                  E+=  p[k_py][4];
                  for ( int ll=1 ; ll<=3 ; ll++ ) {
                    mom[ll-1] +=    p[k_py][ll];
                  }
                  if ( type[jets_end[j_jet_end][k_dj_end]]%100 != 3 ) {break;}
                }

              }
            }
          }
          double psqrs=0 ;
          for (int ii=0 ; ii<3 ; ii++ ) {
            psqrs += mom[ii]*mom[ii];
          }
          M = sqrt(E*E-psqrs);
        }

        ReconstructedParticleImpl* fafpf = new ReconstructedParticleImpl;
        fafpf->setCharge(0);
        fafpf->setEnergy(E);
        fafpf->setMass(M);
        fafpf->setMomentum(mom);
        for(int j_jet_end=1 ; j_jet_end<=jets_end[0][k_dj_end] ; j_jet_end++) {
          fafpf->addParticle(dynamic_cast<ReconstructedParticle*>(jet_vec->getElementAt(jets_end[j_jet_end][k_dj_end]-1)) );
        }

        fafpfPidHandler.setParticleID(fafpf,
                                      type[jets_end[1][k_dj_end]] % 100,  // maybe flag from boson? could be 0,1, or 2 jets in the fafp that's from boson ..
                                      pdg[0],
                                      0, fafpf_mainPidId, {});

        for ( int j_jet_end=1 ; j_jet_end<=jets_end[0][k_dj_end] ; j_jet_end++ ) {
          fafpfPidHandler.setParticleID(fafpf, type[jets_end[j_jet_end][k_dj_end]], pdg[j_jet_end], 0, fafpf_pidIds[j_jet_end - 1], {});

          if ( elementon[jets_end[j_jet_end][k_dj_end]] > 0 ) {
            FinalElementon_Nav.addRelation(  jet_vec->getElementAt(jets_end[j_jet_end][k_dj_end]-1) , 
                        mcp_pyjets[elementon[jets_end[j_jet_end][k_dj_end]]], 1.0 );
          }
          FinalColourNeutral_Nav.addRelation(  jet_vec->getElementAt(jets_end[j_jet_end][k_dj_end]-1) , fafpf);
        }
        fafpf_vec->addElement(fafpf);

      } // di-jet(end) loop
   
           
        // pre-PS part

      LCCollectionVec* fafpi_vec = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE ) ;
      auto fafpiPidHandler = UTIL::PIDHandler(fafpi_vec);
      auto fafpi_mainPidId = fafpiPidHandler.addAlgorithm("TrueJet_fafpi", {});
      std::array<int, 25> fafpi_pidIds{};
      for (size_t i = 0; i < fafpi_pidIds.size(); ++i) {
          fafpi_pidIds[i] = fafpiPidHandler.addAlgorithm("TrueJet_fafpi_jet_" + std::to_string(i), {});
      }

      LCRelationNavigator InitialElementon_Nav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE ) ;
      LCRelationNavigator InitialColourNeutral_Nav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::RECONSTRUCTEDPARTICLE ) ;

      for (int k_dj_begin=1 ; k_dj_begin<=n_djb ; k_dj_begin++ ) {
        double E=0, M=0 , mom[3]={0} ; 
        int pdg[26]= {0};

        int bosonid=23 ;

             // E and P is sum of (unique) direct decendants of the
             // initial ffbar. Normally, there is only one
             // (a 94), so the sum is actualy just copying the descendant. 
             // For leptons, however, there is no 94 - they go directly to
             // the next generation leptons), so the sum is non-trivial.
             // Sometimes this also happens for quraks, ie. they go directly to
             // the genertor quarks of the string.

        int first=nlund , last=0;
        if ( type[jets_begin[1][k_dj_begin]] == 4 ) {  // i.e. ISR

          pdg[1]=k[elementon[jets_begin[1][k_dj_begin]]][2];
          first= fafp[jets_begin[1][k_dj_begin]] ;
          last= fafp[jets_begin[1][k_dj_begin]] ;
          bosonid=22;

        } else {

        for ( int j_jet_begin=1 ; j_jet_begin<=jets_begin[0][k_dj_begin] ; j_jet_begin++ ) {


          int gen_part = fafp [jets_begin[j_jet_begin][k_dj_begin]];
	  bool take_genpart = false;
	  if (  k[gen_part][4] == 0 ) {
	    take_genpart = true;
	  } else if ( k [k [gen_part] [4]][2] != k[gen_part][2] ||  k [k [gen_part] [4]][1] == 0 ) {
	    take_genpart = true;
	  }
	  if ( take_genpart ) {
            if (  gen_part < first ) {
              first=gen_part;
            }
            if (  gen_part > last ) {
              last=gen_part;
            }
          } else {
            if (  k[gen_part][4] < first ) {
	      first= k[gen_part][4] ;
	    }  
            if (  k[gen_part][5] > last ) {
	      last= k[gen_part][5] ;
	    }  
	    
	  }
	    
          pdg[j_jet_begin]=k[elementon[jets_begin[j_jet_begin][k_dj_begin]]][2];
          int last_backtrack = fafp[jets_begin[j_jet_begin][k_dj_begin]] ;
          int this_pdg = k[last_backtrack][2] ;

          while (  k[last_backtrack][2] == this_pdg ) {
            last_backtrack = k[last_backtrack][3];
            if ( last_backtrack == 0 ) break;
          }

          if ( last_backtrack > 0) {
            if (  abs(k[last_backtrack][2]) > 21 &&   abs(k[last_backtrack][2]) <= 39 ) {
                // ie. any IVB except the gluon
              bosonid =  abs(k[last_backtrack][2]) ;
              streamlog_out(DEBUG3) << " elementon " << fafp[jets_begin[j_jet_begin][k_dj_begin]] 
                                      << " has explicit mother with PDG = " << bosonid << std::endl;
            }
          }

          if ( pdg[0] == 0 ) {

            pdg[0]=abs(k[fafp[jets_begin[j_jet_begin][k_dj_begin]]][2]);

          } 
          else if ( bosonid==23 && abs(k[fafp[jets_begin[j_jet_begin][k_dj_begin]]][2]) != pdg[0] ) {

            bosonid=24;

          } 
        }
        }
        if ( k[last][2] == 92 ) { first = last ; } //  IS THIS RIGHT ? If there is a gluon radiation ?!
	                                           // should it rather be the other way last="first"
        pdg[0]=bosonid;

        if ( last == first ) {
          E=  p[first][4];
          for ( int ll=1 ; ll<=3 ; ll++ ) {
            mom[ll-1] =    p[first][ll];
          }
        } else {

          // if several, sum up but make sure that there are no intruders between first and last.
          for ( int j_jets_begin=1 ; j_jets_begin<=jets_begin[0][k_dj_begin] ; j_jets_begin++ ) {
            for ( int k_py=first ; k_py<=last ; k_py++ ) {
              if ( abs(jet[k_py]) == jets_begin[j_jets_begin][k_dj_begin] && k[k_py][3] < first ) {
                E+=  p[k_py][4];
                for ( int ll=1 ; ll<=3 ; ll++ ) {
                  mom[ll-1] +=    p[k_py][ll];
                }
              }

            }
          }
        }

        double psqrs=0;
        for ( int ii=1 ; ii<=3 ; ii++ ) {
          psqrs += mom[ii-1]*mom[ii-1];
        }
         
        M = sqrt(E*E-psqrs);


        ReconstructedParticleImpl* fafpi = new ReconstructedParticleImpl;
        fafpi->setCharge(0);
        fafpi->setEnergy(E);
        fafpi->setMass(M);
        fafpi->setMomentum(mom);

        fafpiPidHandler.setParticleID(fafpi, type[jets_begin[1][k_dj_begin]] % 100, pdg[0], 0, fafpi_mainPidId, {});

        for(int j_jet_begin=1 ; j_jet_begin<=jets_begin[0][k_dj_begin] ; j_jet_begin++) {
          fafpiPidHandler.setParticleID(fafpi, type[jets_begin[j_jet_begin][k_dj_begin]], pdg[j_jet_begin], 0, fafpi_pidIds[j_jet_begin - 1], {});
          fafpi->addParticle(dynamic_cast<ReconstructedParticle*>(jet_vec->getElementAt(jets_begin[j_jet_begin][k_dj_begin]-1)));
        }
        fafpi_vec->addElement(fafpi);
        
        for ( int j_jet_begin=1 ; j_jet_begin<=jets_begin[0][k_dj_begin] ; j_jet_begin++ ) {
          InitialElementon_Nav.addRelation(  jet_vec->getElementAt(jets_begin[j_jet_begin][k_dj_begin]-1), 
                                              mcp_pyjets[fafp[jets_begin[j_jet_begin][k_dj_begin]]], 1.0 );
          InitialColourNeutral_Nav.addRelation(  jet_vec->getElementAt(jets_begin[j_jet_begin][k_dj_begin]-1) , fafpi);
        }
   
      }  // di-jet (begining) loop

      //============================ 

      
      //*****************************
      evt->addCollection( jet_vec,  _trueJetCollectionName);
      evt->addCollection( fafpf_vec,  _finalColourNeutralCollectionName);
      evt->addCollection( fafpi_vec,  _initialColourNeutralCollectionName);

      LCCollection* tjfpcol = 0;
      LCCollection* tjlpcol = 0;
      LCCollection* tjfccol = 0;
      LCCollection* tjlccol = 0;
      tjrcol = truejet_pfo_Nav.createLCCollection() ;
      tjtcol = truejet_truepart_Nav.createLCCollection() ;
      tjfpcol = InitialElementon_Nav.createLCCollection() ;
      tjlpcol = FinalElementon_Nav.createLCCollection() ;
      tjfccol = InitialColourNeutral_Nav.createLCCollection() ;
      tjlccol = FinalColourNeutral_Nav.createLCCollection() ;
      evt->addCollection( tjrcol , _trueJetPFOLink);
      evt->addCollection( tjtcol , _trueJetMCParticleLink);
      evt->addCollection( tjfpcol , _initialElementonLink);
      evt->addCollection( tjlpcol , _finalElementonLink);
      evt->addCollection( tjfccol , _initialColourNeutralLink);
      evt->addCollection( tjlccol , _finalColourNeutralLink);
      //============================ 



      streamlog_out(DEBUG4) << " HEPEVT relation table up with jet assignment, extended status codes, and PFOs w. energy (if any) "  <<std::endl;
      streamlog_out(DEBUG4) << " line      status   pdg    first last  first last  colour  anti-     px         py          pz     "
                              <<"     E            M        jet companion type       PFO/Energy" << std::endl;
      streamlog_out(DEBUG4) << "        init ext'd           parent     daughter          colour                               "
                              <<"                                      jet" << std::endl;

      //      for ( int i_py=1 ; i_py <= nlund ; i_py++){
      for ( int i_py=1 ; i_py <= std::min(nlund,200) ; i_py++){

        int istat=0;
        if ( jet[i_py] != 0 ) {
          static FloatVec www;
          www = truejet_truepart_Nav.getRelatedFromWeights(mcp_pyjets[i_py]);
          istat=int(www[0]);
        }

        if ( k[i_py][1] > 0 &&  k[i_py][1] < 30 ) {

	  streamlog_out(DEBUG4) << std::setw(4) << i_py << std::setw(7) <<  k[i_py][1] << std::setw(5) <<  istat << 
	    std::setw(7) <<  k[i_py][2] << std::setw(7) <<  k[i_py][3] << std::setw(5) <<  k[i_py][6] <<
                 std::setw(7) <<  k[i_py][4] << std::setw(5) << k[i_py][5] <<
                 std::setw(7) <<  k[i_py][7] << std::setw(7) << k[i_py][8] <<
                 std::setw(12) <<  p[i_py][1] << 
                 std::setw(12) <<  p[i_py][2] << std::setw(12) <<  p[i_py][3] << 
                 std::setw(12) <<  p[i_py][4] << std::setw(12) << p[i_py][5] <<
 
                 std::setw(7) << jet[i_py]  <<
                 std::setw(7) << companion[abs(jet[i_py])] << std::setw(7) <<  type[abs(jet[i_py])];
        } else {
	  streamlog_out(DEBUG3) << std::setw(4) << i_py << std::setw(7) <<  k[i_py][1] << std::setw(5) <<  istat << 
                 std::setw(5) <<  k[i_py][2] << std::setw(7) <<  k[i_py][3] <<  std::setw(5) <<  k[i_py][6] <<
                 std::setw(7) <<  k[i_py][4] << std::setw(5) << k[i_py][5] <<
                 std::setw(7) <<  k[i_py][7] << std::setw(7) << k[i_py][8] <<
                 std::setw(12) <<  p[i_py][1] << 
                 std::setw(12) <<  p[i_py][2] << std::setw(12) <<  p[i_py][3] << 
                 std::setw(12) <<  p[i_py][4] << std::setw(12) << p[i_py][5] <<
 
                 std::setw(7) << jet[i_py]  <<
                 std::setw(7) << companion[abs(jet[i_py])] << std::setw(7) <<  type[abs(jet[i_py])];
        }

        if ( reltrue ) {
          LCObjectVec recovec = reltrue->getRelatedFromObjects( mcp_pyjets[i_py]);
          if ( recovec.size() > 0 ) { // if reconstructed
            for ( unsigned i_reco=0 ; i_reco<recovec.size() ; i_reco++ ) { //  reco-of-this-true loop
              if ( k[i_py][1] > 0 &&  k[i_py][1] < 30 ) {

               ReconstructedParticle* reco_part  = dynamic_cast<ReconstructedParticle*>(recovec[i_reco]);
  	       streamlog_out(DEBUG4) << " [" << std::setw(8) << reco_part->getEnergy() <<"] ";
              } else {
                ReconstructedParticle* reco_part  = dynamic_cast<ReconstructedParticle*>(recovec[i_reco]);
  	        streamlog_out(DEBUG3) << " [" << std::setw(10) << reco_part <<" /"<< std::setw(8) << reco_part->getEnergy() <<"] ";
              }
            }
          } else if ( istat == 1 ) {
            if ( k[i_py][1] > 0 &&  k[i_py][1] < 30 ) {
    	      streamlog_out(DEBUG4) << " [" << std::setw(10) << "N.A   "<<" /"<< std::setw(8) << "N.A   "<<"] ";
            } else  {
    	      streamlog_out(DEBUG3) << " [" << std::setw(10) << "N.A   "<<" /"<< std::setw(8) << "N.A   "<<"] ";
            }
          }
        }
        if ( k[i_py][1] > 0 &&  k[i_py][1] < 30 ) {
	  streamlog_out(DEBUG4) << std::endl;
        } else {
	  streamlog_out(DEBUG3) << std::endl;
        }
      }
      streamlog_out(DEBUG6) <<  std::endl;

      streamlog_out(DEBUG6) << "       jet     fafp-beg  qrk/lept  fafp-end fafp-boson  type   "
                              <<"dijet-beg  dijet-end    boson" << std::endl;
      for ( int i_jet=1 ; i_jet <= njet ; i_jet++){
	streamlog_out(DEBUG6) << std::setw(10) << i_jet << std::setw(10) <<  fafp[i_jet] << 
                 std::setw(10) <<  elementon[i_jet] << std::setw(10) <<  fafp_last[i_jet] << 
                 std::setw(10) <<  fafp_boson[i_jet] << std::setw(10) << type[i_jet] << 
                 std::setw(10) << dijet_begining[i_jet]  <<
                 std::setw(10) << dijet_end[i_jet] << std::setw(10) <<  boson[i_jet]<<std::endl;
      }
      streamlog_out(DEBUG6) <<  std::endl;
      streamlog_out(DEBUG6) << "       Colour-neutral objects at the beginning of the parton shower " << std::endl;
      streamlog_out(DEBUG6) << "    number    px          py           pz          E           M   "
                              <<"    type       PDGs        jets " << std::endl;
      int glurad=0;
      for ( unsigned i_cn_b=0 ; i_cn_b<fafpi_vec->size() ; i_cn_b++ ) {
        
        ReconstructedParticle* fafpi = dynamic_cast<ReconstructedParticle*>(fafpi_vec->at(i_cn_b));
        LCObjectVec jts = InitialColourNeutral_Nav.getRelatedFromObjects(fafpi);

        const double* mom=fafpi->getMomentum();
	streamlog_out(DEBUG6) << std::setw(7) << i_cn_b+1 << std::setw(12) << mom[0] << 
                             std::setw(12) <<  mom[1] << std::setw(12)  << mom[2] << 
                             std::setw(12)  <<fafpi->getEnergy() << std::setw(12) << fafpi->getMass() << 
                             std::setw(7)  << fafpi->getParticleIDs()[0]->getType() << "   " ;
        for ( unsigned i_pid=0 ; i_pid < fafpi->getParticleIDs().size() ; i_pid++ ) {
	  if (fafpi->getParticleIDs()[i_pid]->getType()>=100 ) glurad=1;
          if (i_pid < fafpi->getParticleIDs().size()-1 ) {
            streamlog_out(DEBUG6)  << std::setw(3)<< fafpi->getParticleIDs()[i_pid]->getPDG()  << "," ;
          } else {
            streamlog_out(DEBUG6)  << std::setw(3)<< fafpi->getParticleIDs()[i_pid]->getPDG() << "   " ;
          }
        }
        for ( unsigned j_jet_begin=0 ; j_jet_begin<jts.size() ; j_jet_begin++ ) {
           int i_jet=jts[j_jet_begin]->ext<TJindex>(); 
           if ( j_jet_begin < jts.size()-1 ) {
             streamlog_out(DEBUG6) << std::setw(3) << i_jet << "," ; 
           } else {
             streamlog_out(DEBUG6) << std::setw(3) << i_jet  ;
           }
        }
        streamlog_out(DEBUG6) << std::endl;
      }
      if ( glurad == 1 ) {
        streamlog_out(DEBUG6) << std::endl;
        streamlog_out(DEBUG6) <<    "              gluon radiation in event : "<< std::endl;
        for ( unsigned j_cn_b=0 ; j_cn_b<fafpi_vec->size() ; j_cn_b++ ) {
        
          ReconstructedParticle* fafpi = dynamic_cast<ReconstructedParticle*>(fafpi_vec->at(j_cn_b));
          LCObjectVec jts = InitialColourNeutral_Nav.getRelatedFromObjects(fafpi);

          for ( unsigned j_pid_cn_b=1 ; j_pid_cn_b < fafpi->getParticleIDs().size() ; j_pid_cn_b++ ) {
	    if (fafpi->getParticleIDs()[j_pid_cn_b]->getType()>=100 ) {
              int i_jet=jts[j_pid_cn_b-1]->ext<TJindex>(); 
              streamlog_out(DEBUG6)  << std::setw(18)<< i_jet << " is from " 
                                       << fafpi->getParticleIDs()[j_pid_cn_b]->getType()/100  << std::endl;
            }
          }
        }
      }

      streamlog_out(DEBUG6) <<  std::endl;
      streamlog_out(DEBUG6) << "       Colour-neutral objects at the end of the parton shower " << std::endl;
      streamlog_out(DEBUG6) << "    number    px          py           pz          E           M   "
                              <<"    type       PDGs        jets " << std::endl;
      for ( unsigned i_cn_e=0 ; i_cn_e<fafpf_vec->size() ; i_cn_e++ ) {
        
        ReconstructedParticle* fafpf = dynamic_cast<ReconstructedParticle*>(fafpf_vec->at(i_cn_e));
        LCObjectVec jts = FinalColourNeutral_Nav.getRelatedFromObjects(fafpf);

        const double* mom=fafpf->getMomentum();
	streamlog_out(DEBUG6) << std::setw(7) << i_cn_e+1 << std::setw(12) << mom[0] << 
                             std::setw(12) <<  mom[1] << std::setw(12)  << mom[2] << 
                             std::setw(12)  <<fafpf->getEnergy() << std::setw(12) << fafpf->getMass() << 
                             std::setw(7)  << fafpf->getParticleIDs()[0]->getType()  << "   ";
        if (  fafpf->getParticleIDs().size() < 3 ) streamlog_out(DEBUG6)  << std::setw(3) <<" ";
        for ( unsigned i_pid=0 ; i_pid < fafpf->getParticleIDs().size() ; i_pid++ ) {
          if (i_pid < fafpf->getParticleIDs().size()-1 ) {
            streamlog_out(DEBUG6)  << std::setw(3)<< fafpf->getParticleIDs()[i_pid]->getPDG()  << "," ;
          } else {
            streamlog_out(DEBUG6)  << std::setw(3)<< fafpf->getParticleIDs()[i_pid]->getPDG() << "   " ;
          }
        }
        if (  jts.size() < 2 ) streamlog_out(DEBUG6)  << std::setw(3) <<" ";
        for ( unsigned j_jet_end=0 ; j_jet_end<jts.size() ; j_jet_end++ ) {
           int i_jet=jts[j_jet_end]->ext<TJindex>(); 
           if ( j_jet_end < jts.size()-1 ) {
             streamlog_out(DEBUG6) << std::setw(3) << i_jet << "," ; 
           } else {
             streamlog_out(DEBUG6) << std::setw(3) << i_jet  ;
           }
        }
        streamlog_out(DEBUG6) << std::endl;
      }
      streamlog_out(DEBUG6) << std::endl;
      

    }  // endif( mcpcol != NULL )


    //*****************************
    // (this guy is not deleted by lcio, so ...
    if ( reltrue ) delete reltrue ;
    // ... to avoid a memory leak!)
    _nEvt ++ ;
    //============================ 
}



void TrueJet::getPyjets(LCCollection* mcpcol )
{
  const double* ptemp;
  _higgs_to_glue_glue = false;
  bool _generator_input_format = false;
  int higgs_line = 0 ;
  int k_py=0;
  int nMCP = mcpcol->getNumberOfElements()  ;
  _top_event = false ;
  first_line=1;

  if ( nMCP > 4000 ) {
    streamlog_out(WARNING)  << " More than 4000 MCParticles in event " << evt->getEventNumber() 
                              << ",   run  " << evt->getRunNumber() << std::endl ;
  }

  for(int j_mcp=0; j_mcp< std::min(nMCP,4000) ; j_mcp++){

    MCParticle* mcp = dynamic_cast<MCParticle*>( mcpcol->getElementAt( j_mcp ) ) ;
       //  (assume this is not needed, ie. that the k-vector  can be set up correctly
       //    also for created-in-simulation  particles. The issue is if decay-products
       //    of a given particle will be consecutive).
       // if (mcp->getGeneratorStatus() == 1 || mcp->getGeneratorStatus() == 2 ) {

      if (  j_mcp == 0 && mcp->getSimulatorStatus()==0 ) { _generator_input_format = true ; }
      k_py++; 
      mcp_pyjets[k_py]=mcp;
      mcp->ext<MCPpyjet>()=k_py;
      ptemp=mcp->getMomentum();
      p[k_py][1]=ptemp[0]; 
      p[k_py][2]=ptemp[1]; 
      p[k_py][3]=ptemp[2]; 
      p[k_py][4]=mcp->getEnergy();
      p[k_py][5]=mcp->getMass();
      k[k_py][1]=mcp->getGeneratorStatus();
      if ( k[k_py][1] == 2 ) {  k[k_py][1]=11 ;}
      if ( k[k_py][1] == 3 ) {  k[k_py][1]=21 ;}
      if ( k[k_py][1] == 4 ) {  k[k_py][1]=22 ;}
      if (mcp->isOverlay() ) {  k[k_py][1]=k[k_py][1]+30 ; }
      k[k_py][2]=mcp->getPDG();
      if ( abs(k[k_py][2]) == 6 ) { _top_event = true ; }
      if ( k[k_py][2] == 25 ) {
        higgs_line = k_py ;
      }

 // }
  }
  nlund=k_py;

  for (int j_py=1; j_py <= nlund ; j_py++ ) {
    k[j_py][3] = 0;
    k[j_py][4] = 0;
    k[j_py][5] = 0;
    k[j_py][7] = 0;
    k[j_py][8] = 0;
  }
  streamlog_out(DEBUG1) << " before Motherless check: first_line = " << first_line << std::endl ; 

  for(int j_mcp=0; j_mcp< nMCP ; j_mcp++){

    MCParticle* mcp = dynamic_cast<MCParticle*>( mcpcol->getElementAt( j_mcp ) ) ;
    int i_py=0;
    if ( mcp->ext<MCPpyjet>() != 0 ) {

      i_py= mcp->ext<MCPpyjet>();
      if (  mcp->getParents().size() != 0 ) {
        k[i_py][3]=mcp->getParents()[0]->ext<MCPpyjet>();
	if ( k[i_py][3] ==  higgs_line ) {
	  if ( k[i_py][2] == 21 ) {
	    _higgs_to_glue_glue = true;
	  }
	}  
      } else {
        if ( i_py > first_line+1 &&  ! mcp->isOverlay() &&  k[i_py][3] == 0 && abs(mcp->getPDG()) != 6 && abs(mcp->getPDG()) != 24 ) {
          streamlog_out(MESSAGE4) << " Motherless generator particle " << i_py << ". pdg: " << mcp->getPDG() << std::endl ; 
          streamlog_out(MESSAGE4) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
          k[i_py][3]=-1;
        }
       
      }

      k[i_py][7]=mcp->getColorFlow()[0];
      k[i_py][8]=mcp->getColorFlow()[1];

      if (  mcp->getParents().size() > 1 ) {
        k[i_py][6]=mcp->getParents()[mcp->getParents().size()-1]->ext<MCPpyjet>();
      } else {
        k[i_py][6]=0;
      }
      
      if ( mcp->getDaughters().size() != 0 ) {

           // check daughter-status. If any of the daughters are created in simulation,
           // or has simulator-status 0 (indicating that it was ignored by Geant), 
           // set k[][1] to stable of the mother, and change it to 0 for direct descendants.
           // if the current particle has  k[][1] == 0, propagate to its descendants, that way
           // all generations of a  k[][1] == 0 ancestor will have  k[][1] == 0 once
           // the loop is done.
           
           // (A further case is if all daughters have genstat 1, simstat !=0, but the sum of their energies is
           // different from that of the mother (indicating that Geant dis something un-documented with them)
           // This case can't be fixed, yet as the check should only be applied to lines after the parton-shower,
           // which we can't determine yet. Instead it is taken care of in the gouping method.)

        if ( ! _generator_input_format ) { // none of the below can happen if we are looking at the generator input MCparticle collection 
          if ( k[i_py][1]!=0 ) {
            int has_genstat0_daughter=0;
            for ( unsigned i_dau=0 ; i_dau <  mcp->getDaughters().size() ; i_dau++ ) {
              int j_dau=mcp->getDaughters()[i_dau]->ext<MCPpyjet>();
              if ( mcp->getDaughters()[i_dau]->getGeneratorStatus()==0 || (mcp->getDaughters()[i_dau]->getSimulatorStatus()==0 && k[j_dau][1]==1 )) {
                has_genstat0_daughter=1;
              }
            }
            if ( has_genstat0_daughter==1 ) {

              if ( k[i_py][1] < 30 ) {
                k[i_py][1]=1;
              }

              for ( unsigned i_dau=0 ; i_dau <  mcp->getDaughters().size() ; i_dau++ ) {
                int j_dau=mcp->getDaughters()[i_dau]->ext<MCPpyjet>();
                if ( k[j_dau][1] < 30 ) {
                  k[j_dau][1]=0;
                }
              }

            }

          } else {

            for ( unsigned i_dau=0 ; i_dau <  mcp->getDaughters().size() ; i_dau++ ) {
              int j_dau=mcp->getDaughters()[i_dau]->ext<MCPpyjet>();
              if ( k[j_dau][1] < 30 ) {
                k[j_dau][1]=0;
              }
            }
          }
        }

           // fill in k[][4] and [5]. Note that it is possible that one gets
           // k[][1]=1 but k[][4] and [5] != 0. This means that a generator
           // particle end in a simulation particle, which is OK for the
           // logic elsewhere !

        for ( unsigned i_dau=0 ; i_dau <  mcp->getDaughters().size() ; i_dau++ ) {
          int j_dau=mcp->getDaughters()[i_dau]->ext<MCPpyjet>();
          if ( j_dau < k[i_py][4] ||  k[i_py][4] == 0 ) {
            k[i_py][4]=j_dau;
          }
          if ( j_dau > k[i_py][5] ) {
            k[i_py][5]=j_dau;
          }
	  if ( abs(k[i_py][2]) == 6 ) {  // fix for top quark - possibly not needed with whiz2 ?
	    if (k[k[i_py][4]][2] != 94 ) { 
   	      k[i_py][5]= k[i_py][4]+1;
	      k[k[i_py][5]][3] = i_py;
	    }  
	  }

        }
      } 
    }  // endif ( mcp->ext<MCPpyjet>() != 0 )
  } // end MCP loop


  // At this point basically 99% of the job is done. However, there is still the CMShower lines (pdg=94).
  // These have several mothers and several daughters. The 4-momentum is conserved in total, but not for
  // the individual components one-by-one. Here we will short-circuit them, i.e. the history codes will
  // connect mother and daughters directly. This is no problem if all pdg:s are different in the group
  // of particles going in and out of the 94. The problem arrises when this is not the case (e.g. e+e-e+e-).
  // Here we need to look at the momenta of the particles and match which dauther is closest to which
  // mother. Remember, they will *not* be identical, so there is a bit of guessing.
  // Note that the mothers of a 94 are not necessarily adjecent. They *are* in the range between k[line94][3] 
  // and k[line94][6], but there might be "intruders" between these lines.


  streamlog_out(DEBUG4)  << " HEPEVT relation table before 94-fixing"  <<std::endl;
  streamlog_out(DEBUG4)  << "       line      status    pdg       parent    first     last     second             anti   " << std::endl;
  streamlog_out(DEBUG4)  << "                                             daughter daugther    parent   colour   colour " << std::endl;
  for ( int j_py=1 ; j_py<=nlund ; j_py++ ) {
    if (  k[j_py][1] != 0 &&   k[j_py][1] < 30 ) {
      streamlog_out(DEBUG4)  << std::setw(10) << j_py << std::setw(10) <<  k[j_py][1] << 
                 std::setw(10) <<  k[j_py][2] << std::setw(10) <<  k[j_py][3] << 
                 std::setw(10) <<  k[j_py][4] << std::setw(10) << k[j_py][5]  << 
                 std::setw(10) << k[j_py][6] << std::setw(10)<< k[j_py][7] << std::setw(10) << k[j_py][8] << std::endl;
    } else { 
      streamlog_out(DEBUG3)  << std::setw(10) << j_py << std::setw(10) <<  k[j_py][1] << 
                 std::setw(10) <<  k[j_py][2] << std::setw(10) <<  k[j_py][3] << 
                 std::setw(10) <<  k[j_py][4] << std::setw(10) << k[j_py][5]  << 
                 std::setw(10) << k[j_py][6] << std::setw(10)<< k[j_py][7] << std::setw(10) << k[j_py][8] << std::endl;
    }      
  }

  //  94 fix. 

  for ( int i_py=1 ; i_py<=nlund ; i_py++ ) {

    // if the daughters of the parents to a 94 is not (only) the 94,
    // check if between the 94 line and the second daughter of the
    // parent of the 94 actually has the particle *going in* to the 94 (not the
    // 94) as parent. If so, change the dautghter lines of the incomming particle
    // to these lines (and ignore the 94). If not, move those that are not daughters
    // to the 94 instead, and set also the second daughter of the parent to the
    // 94.
 
    if ( k[i_py][2] == 94 && k[i_py][1] < 30 ) {
       // A 94: treat that:
       //   case 1: 2->2 , clear-cut connect in to out dirctly
       //   case 2: 1->1+gamma : FSR, chnage 94 to 95, else unchanged.
       //   case 3: 4->4, unique flavours, incomming has only 94 as daughter, outgoind only 94 as mother: : as case 1
       //   case 4: as case 3, but flavours not unique: Check directions etc, then connect directly
       //   case 5: 4->4, some weidness in parent-daughter relations: to study
       //   case 6: n->m, n != m: to study
       //   case 7: odd->m : Does that ever happen, except for the FSR case?
       //   case 8: Anything else : Does that ever happen
       //   case X: any of the above, but with intruders between first and last mother and/or daughter: to study
       //

      int mothers[10] ={0};
      int daughters[10][10] ={0};

      int line94 = i_py ;
      if ( _whiz1 ) { stdhep_reader_bug_workaround(line94) ; }

      int first_94_mother= k[line94][3];
      int last_94_mother= k[line94][6];
      int first_94_daughter =  k[line94][4];
      int last_94_daughter =  k[line94][5];

      int i_mo=0;
      for ( int mother = first_94_mother ; mother <= last_94_mother; mother++ ) {
        if ( (abs(k[mother][2]) == 24 && abs(k[k[mother][3]][2]) !=6 ) || k[mother][2]==23 || k[mother][2]==25 ) continue;
        if ( abs(k[mother][4]) !=  line94 ) continue; // intruder, go to next line

        i_mo++;
        mothers[0]=i_mo;
        mothers[i_mo]=mother;
        int i_dau=0;	
	daughters[i_mo][0]=0;
        int pdg_mother=k[mother][2];
        float p_mother[6];
        for ( int ip=1; ip <=5 ; ip++ ) {
          p_mother[ip]=p[mother][ip];
        }

        if ( k[mother][4] == line94 && k[mother][5] == line94 ) {

            // normal case

	  if ( _whiz1 ) { // Maybe useful not only for whiz1 ? If so, also should check k[][7:8] == 0
            int n_94_mother = 0 ;
            for ( int j_py = first_94_mother ; j_py <=  last_94_mother ; j_py++ ) {
              if ( k[j_py][4] == line94 ) { n_94_mother++ ; }
            } 
	    if ( n_94_mother == 2 ) {
	      if ( k[mother][2] > 0 ) {
	         k[mother][7] = i_py ; k[mother][8] = 0 ;
	         k[k[mother][3]][7] = i_py ; k[k[mother][3]][8] = 0 ;
              } else {
	         k[mother][7] = 0 ; k[mother][8] = i_py ;
	         k[k[mother][3]][7] = 0 ; k[k[mother][3]][8] = i_py ;
              }
              int oma=mother;
              while ( oma > 3 ) {
	         k[oma][7] = k[mother][7] ; k[oma][8]=k[mother][8] ;
                 oma=k[oma][3];
              }
	    }  
	  } // endif ( _whiz1 )
 
            // find the daughter that best matches the current mother

          int best_match = 0;
          double min_dist = 1.0e20;          
          for ( int daughter=first_94_daughter ; daughter <=  last_94_daughter ; daughter++ ) { 
            if ( k[daughter][2] == pdg_mother ) {
              double dist=((p_mother[1]-p[daughter][1])*(p_mother[1]-p[daughter][1]) +
                           (p_mother[2]-p[daughter][2])*(p_mother[2]-p[daughter][2]) +
                           (p_mother[3]-p[daughter][3])*(p_mother[3]-p[daughter][3]));
              if ( dist < min_dist ) {
                best_match=daughter;
                min_dist=dist;
              }
            }
          }

          if ( best_match == 0 ) {
            streamlog_out(ERROR)  << " No best match ... " << std::endl;
          } else {

            i_dau++;
            daughters[i_mo][0]=i_dau;
            daughters[i_mo][i_dau]=best_match;

          } 

        } else {  // mother has several daughters, find the daughter that has the same pdg as the mother (actually an intruder)

          bool found_daughter = false ;
          for ( int daughter=k[mother][4]  ; daughter <=  k[mother][5] ; daughter++ ) {
            if ( k[daughter][3] == mother && daughter != line94 ) {
              found_daughter = true;
              i_dau++;
              daughters[i_mo][0]=i_dau;
              daughters[i_mo][i_dau]=daughter;
            }
          }

          if ( ! found_daughter ) {
            streamlog_out(ERROR)  << " Only daughter of the mother of a 94 is neither the 94, "
                                    <<"nor a particle with same pdg as the mother " << std::endl;
          }
        }

        if (daughters[i_mo][0] == 0 ) {
          streamlog_out(ERROR)  << " Did not find any daughters to the mother of a 94" << std::endl;
        } 
      }

        // Here all is sorted out, update the found history relations.

      for ( int j_mo=1 ; j_mo <= mothers[0] ; j_mo++ ) {

        int mother=mothers[j_mo];

        for ( int j_dau=1 ; j_dau <= daughters[j_mo][0] ; j_dau++ ) {
          int daughter=daughters[j_mo][j_dau];
          if ( j_dau == 1 ) k[mother][4]=daughter;
          k[mother][5]=daughter;
          k[daughter][3]=mother;
        }

      }
    }  // endif this is a 94 line  
  } // for all lines

  streamlog_out(DEBUG3)  << " HEPEVT relation table, after 94 fixing"  <<std::endl;
  streamlog_out(DEBUG3)  << "       line      status    pdg       parent    first     last     second             anti   " << std::endl;
  streamlog_out(DEBUG3)  << "                                             daughter daugther    parent   colour   coulour " << std::endl;

  for ( int j_py=1 ; j_py<=nlund ; j_py++ ) {
     streamlog_out(DEBUG3)  << std::setw(10) << j_py << std::setw(10) <<  k[j_py][1] << 
                 std::setw(10) <<  k[j_py][2] << std::setw(10) <<  k[j_py][3] << 
                 std::setw(10) <<  k[j_py][4] << std::setw(10) << k[j_py][5]  << std::setw(10) << k[j_py][6] << 
                 std::setw(10) <<  k[j_py][7] << std::setw(10) << k[j_py][8] << std::endl;
  }
} 

void TrueJet::stdhep_reader_bug_workaround( int line94 )
{
  int first_94_mother= k[line94][3];
  int last_94_mother= k[line94][6];
  int first_94_daughter =  k[line94][4];
  int last_94_daughter =  k[line94][5];



  bool inconsitent_pid =
    (( k[first_94_daughter][2] != k[first_94_mother][2] && k[first_94_daughter][2] != k[last_94_mother][2] ) ||  // inconsistent pdgs
     ( k[ last_94_daughter][2] != k[first_94_mother][2] && k[ last_94_daughter][2] != k[last_94_mother][2] )) ;

  bool inconsitent_E =
    ( abs((p[last_94_mother][4]+p[first_94_mother][4]) -p[line94][4])/p[line94][4] > 0.001 ) ;  // inconsitent E

  bool gluon_daughters = (  k[first_94_daughter][2] == 21 || k[last_94_daughter][2] == 21 ) ; // 94->gluons is OK as is.

  if ( ( inconsitent_pid || inconsitent_E ) && ( ! gluon_daughters )) {
    int right_mother = 0 ;
    int sought_mother = 0 ;
    int new_mother = 0; 
    bool first = true;
    while ( 1 ) {


      if ( first ) {
        right_mother = first_94_mother;
      } else {
        right_mother = last_94_mother;
      }

      if ( first ) {
        if ( inconsitent_pid ) {
          if ( k[first_94_daughter][2]*k[last_94_daughter][2] > 0 ) {

            streamlog_out(ERROR)  << " 94 DAUGHTERS wrong of 94 on line " << line94 
                                    << " : both fermions or both anti-fermions " <<std::endl;
            // Probably does not happen ...

          } else if ( k[first_94_mother][2]*k[last_94_mother][2] > 0 ) {

            streamlog_out(DEBUG4)  << " 94 MOTHERS wrong of 94 on line " << line94 
                                     << " : both fermions or both anti-fermions " <<std::endl;

            // In this case, already the particle or antiparticle nature of the daughter tells 
            // us the mother of which particle we need
            // to find in the record to correct things.

          } else {

            streamlog_out(ERROR)  << " 94 on line " << line94 << " is wrong : differnt PDGs, but both mothers and daughters are "
                              << "fermions/anti-fermion pairs ???  " <<std::endl;
            // Probably does not happen ...
       
          } 
        } 
      }
      

      if ( k[first_94_daughter][2] == k[right_mother][2] ) {  
        sought_mother = last_94_daughter; 
      } else {
        sought_mother = first_94_daughter; 
      }

      new_mother = 0 ;
      double target[5]={0};
      for ( int jjj=1 ; jjj<=4 ; jjj++ ){target[jjj]=p[line94][jjj]-p[right_mother][jjj];}
      for ( int search_line = line94-1 ; search_line >= 1 ; search_line-- ) {
        if ( k[search_line][2] == k[sought_mother][2] ) {
          if ( abs(p[search_line][1] - target[1])/ abs(target[1]) < 0.05 &&
               abs(p[search_line][2] - target[2])/ abs(target[2]) < 0.05 &&
               abs(p[search_line][3] - target[3])/ abs(target[3]) < 0.05    ) {
            new_mother=search_line; 
            break ;
          }
        }
      }
      if ( new_mother != 0 ) { break ; }
      if ( ! first ) { break ; }
      first = false;
    } // while ( 1 )


    if ( new_mother == 0 ) {
      streamlog_out(ERROR)  << " 94 on line " << line94 << " is wrong, but could not figure out un what way. " ;
    } else {
      if ( new_mother < right_mother ) {
        first_94_mother = new_mother ; last_94_mother = right_mother ;
      } else {
        first_94_mother = right_mother ; last_94_mother = new_mother ;
      }
    }

    k[line94][3] = first_94_mother ; k[line94][6] = last_94_mother ; 
    k[first_94_mother][4] = line94 ;k[first_94_mother][5] = line94 ;
    k[last_94_mother][4] = line94 ;k[last_94_mother][5] = line94 ;
         
  } // if(inconsitent_pid || inconsitent_E )
}
void TrueJet::true_lepton()
{
  int ihard_lepton_1[4011]={0} ;
  int ihard_lepton[4011]={0};
  int lept=0;

      //! count number of hard leptons. The condition is: it should be a lepton lepton (pdg 11 to 16). 
      //! In whiz2, it is hard if the status of its parent is 21, or it comes from a W or a Z.
      //! In whiz1, it is hard either if it is  comming from line 3, or if comes from line 1, goes to 
      //! a single daughter which is an electron. (The second part is needed for beam-remenant 
      //! electrons in odd-f events)

   n_hard_lepton = 0;
   int start_loop=1;

   if ( _top_event ) { // start looking for leptons after the last top - there are un-connected leptons before that to avoid!
     for ( int k_py=nlund ; k_py>= 1 ; k_py-- ) {
       if ( abs(k[k_py][2]) == 6 ) {
         start_loop = k_py ;
         break ;
       }
     } 
   }
 
   for ( int k_py= start_loop  ; k_py<= nlund ; k_py++ ) {
     if ( ( abs(k[k_py][2]) >= 11 && abs(k[k_py][2]) <= 16 ) )  {
       streamlog_out(DEBUG1) << " found lepton, k_py = " << k_py << ", k[k_py][3] = " 
                               << k[k_py][3] << ", k[3][3] = " << k[3][3] <<  std::endl ; 
       if (
	   (( ! _whiz1 ) && ((k  [k[k_py][3]][1] == 21 ) || ( (abs(k [k[k_py][3]][2]) == 23 ||   
                             abs(k [k[k_py][3]][2]) == 24 ) && k[k [k[k_py][3]][3]][1] == 21  )) ) ||
           (   _whiz1 && ((k[k_py][3] == 3 && k[3][3] > 0) || (k[k_py][3] == -1) || 
                          (k[k_py][3] == 1 && k[k_py][4] ==k[k_py][5] && abs(k[k[k_py][5]][2])==11 )))  ) {

         n_hard_lepton++;
         ihard_lepton_1[n_hard_lepton]=k_py;
         ihard_lepton[n_hard_lepton]=0;
         streamlog_out(DEBUG4) << " lepton on line " << k_py  <<" hard according to condition 1 " << std::endl;
       }
       else if (abs(k[k[k_py][3]][2]) == 25 || 
          ((abs(k[k[k_py][3]][2]) == 24 || abs(k[k[k_py][3]][2]) == 23) && 
            abs(k[k[k[k_py][3]][3]][2]) == 25)) { // lepton from higgs, or from H->ZZ*/WW*
         n_hard_lepton++;
         ihard_lepton_1[n_hard_lepton]=k_py;
         ihard_lepton[n_hard_lepton]=0;
         streamlog_out(DEBUG4) << " lepton on line " << k_py  <<" hard according to condition 2 " << std::endl;
       }
       else if ( k[k[k[k_py][3]][3]][1] == 22 && k[k_py][1] < 20 ) { // grandmother has code 22 (= incomming), 
                                                                     // lepton has code < 20
         n_hard_lepton++;
         ihard_lepton_1[n_hard_lepton]=k_py;
         ihard_lepton[n_hard_lepton]=0;
         streamlog_out(DEBUG4) << " lepton on line " << k_py  <<" hard according to condition 3 " << std::endl;
       }
       else if ( _top_event && ( abs(k[k[k_py][3]][2]) == 24) ) { // ultimate ( via only W:s ) ancestor is a top quark
	 int oma = k[k_py][3] ;
	 while ( 1 ) {
	   if ( oma < 3 ) break ;
           if  ( abs(k[k[oma][3]][2]) == 24 ) {  
             oma=k[oma][3] ;
	   } else if (  abs(k[k[oma][3]][2]) == 6 ) {
             n_hard_lepton++;
             ihard_lepton_1[n_hard_lepton]=k_py;
             ihard_lepton[n_hard_lepton]=0;
             streamlog_out(DEBUG4) << " lepton on line " << k_py  <<" hard according to condition 4 " << std::endl;
             break;
           } else {
             break ;
           }
         }
       }
     }
   }

   // here n_hard_leptons is the number of hard leptons (i.e. directly from the initial state of from a boson).
   // and ihard_lepton_1[..] contains the list of pyjets lines of these
   
   if ( n_hard_lepton == 0 ) return ;

      //! Sort the leptons in "color singlets", ie. in groups of flavour/anti-flavour. In most events this is
      //! the order they come in, but in the case of all flavour being the same, it isn't

   // Outcome of this loop is the list of hard leptons in ihard_lepton which is ordered in pairs
   // from the same boson, or if that cant't be figured out, by flavour - anti-flavour pairs.

   for (int j_lept= 1; j_lept<=n_hard_lepton; j_lept++ ) {
     if ( ihard_lepton[j_lept] == 0 ) {
       for ( int l_lept=1; l_lept<=n_hard_lepton ; l_lept++ ) {
         if ( ihard_lepton_1[l_lept] != 0 ) {
           ihard_lepton[j_lept]= ihard_lepton_1[l_lept];
           ihard_lepton_1[l_lept]= 0;
           break ;
         }
       }

       int k_py = ihard_lepton[j_lept] ;
       int n_hint = 0;
       int n_boson = 0;
       int hint_companion = 0 ;
       int boson_companion = 0 ;
       int fallback_companion = 0;
       for ( int l_lept=1; l_lept <= n_hard_lepton ; l_lept++ ) {

         int l_py = ihard_lepton_1[l_lept] ;

         if (  flavour(k[l_py][2]) == flavour(k[k_py][2]) && k[l_py][2]*k[k_py][2] <0 ) { // consider only 
                                                                                          // lepton-flavour/anti-lepton flavour groupings
           if ( ( k[l_py][2] > 0 && k[l_py][7] ==  k[k_py][8] &&  k[l_py][7] != 0 )  || 
                ( k[l_py][2] < 0 && k[l_py][8] ==  k[k_py][7] &&  k[l_py][8] != 0 ) ){// a hint on grouping (leptons from the same 94)
             streamlog_out(DEBUG4) << " hint match for lepton on lines " << l_py << " and " << k_py << std::endl;
             hint_companion =  l_lept;
             n_hint++ ;
           }
           if ( k[l_py][3] == k[k_py][3] ) {
             if (  k[k[l_py][3]][2] == 23 || k[k[l_py][3]][2] == 24 || k[k[l_py][3]][2] == 25 ) { // leptons from the same boson
               streamlog_out(DEBUG4) << " boson match for lepton on lines " << l_py << " and " << k_py << std::endl;
               boson_companion =  l_lept;
               n_boson++;
             } 
             if ( k[l_py][7] == 0 && k[l_py][8]==0 &&  k[k_py][7] == 0 && k[k_py][8]==0 ) {// if no hint nor boson parent found, 
                                                                                           // use this: at least it has the right flavour.
               fallback_companion =  l_lept ;
             }
           }
         }
       }

       if ( n_hint == 1 ) {

         ihard_lepton[j_lept+1]= ihard_lepton_1[hint_companion] ; 
         ihard_lepton_1[hint_companion]=0;

       } else if (n_boson == 1 ) {

         ihard_lepton[j_lept+1]=ihard_lepton_1[boson_companion] ; 
         ihard_lepton_1[boson_companion]=0;

       } else {

         if ( n_hint > 1 || n_boson > 1 ) {
           streamlog_out(DEBUG4)<< " Problem finding companion lepton to lepton on line " << k_py 
                                 << " : n_boson = " << n_boson << " , n_hint = " << n_hint << std::endl;
         }

         if (fallback_companion != 0 ) { // no hint nor boson parent: take any lepton with the right 
                                         // lepton flavour as companion

           ihard_lepton[j_lept+1]=ihard_lepton_1[fallback_companion];
           ihard_lepton_1[fallback_companion]=0;

         }
       }
     }
   }

   // To summerise: at this point we have in ihard_lepton[1  - n_hard_lepton ] the pyjets line numbers of
   // all found hard leptons. These can come directly from the in-state or from a boson. If they are involved in
   // a CM-shower, the ones *going in* into it are the ones in  ihard_lepton.
   
      //! assign to jets. Basically the jet is the index in the list, -ve if it later on goes into a
      //! 94, +ve otherwise.

      // outcome of this loop is the -(jet number) for each hard lepton in jet[ (each hard-lepton pyjets-line) ] and the elementon
      // of the corresponding jet elementon[jet] = hard-lepton pyjets-line. lept on exit will be the first post-CM shower
      // incarnation of the lepton. 

   if (  n_hard_lepton%2 != 0  ) {
     if ( abs(k[ihard_lepton[n_hard_lepton]][2]) != 22 ) {

       // If there is an odd number of initial leptons or photons, possibly we are looking at an
       // odd-fermion sample. In that case, assume - resonably - that the beam-remnant is the first 
       // hard lepton seen. Make it the last instead. That way it will be the one not grouped.
       //
       // TODO: Keep an eye on this : gamma in ME (not tested) ? 
       // Also: better condition for this case (best would be process-type from the run header)

       int beamrem=ihard_lepton[1];
       for ( int  j_lept=2; j_lept<=n_hard_lepton ; j_lept++ ) {
         ihard_lepton[j_lept-1]=ihard_lepton[j_lept];
       }
       ihard_lepton[n_hard_lepton]=beamrem;
     }
   }
   for ( int  j_lept=1; j_lept<=n_hard_lepton ; j_lept++ ) {

     lept = ihard_lepton[j_lept];
     streamlog_out(DEBUG4) << " assigning leptons to jets: ihard_lepton[" << j_lept << "] = " << ihard_lepton[j_lept] 
                           << " jet : " << current_jet+j_lept <<std::endl ; 
   }
   for ( int  j_lept=1; j_lept<=n_hard_lepton ; j_lept++ ) {

     lept = ihard_lepton[j_lept];
     streamlog_out(DEBUG4) << " assigning leptons to jets: ihard_lepton[" << j_lept << "] = " << ihard_lepton[j_lept] 
                           << " and has PDG " << k[lept][2] << " jet : " << current_jet+j_lept<< std::endl ; 

        // ! initially, set to -j_lept

     jet[lept]=-(current_jet+j_lept);

     fafp_last[current_jet+j_lept]=lept;  // so it goes when too specific variable names were choosen ...


     if ( _top_event ) {  // we want to back-track the lepton to the top, so that jet grouping can group all top daughters into on "c.n."
       int oma = lept ;
       int fafp_top = 0;
       while ( 1 ) {
	 if ( oma < 3 ) break ;
         if  ( abs(k[k[oma][3]][2]) == 24 ) {  
           oma=k[oma][3] ;
	 } else if (  abs(k[k[oma][3]][2]) == 6 ) {
           fafp_top = k[oma][3]; 
           oma=k[oma][3] ;
         } else {
           break ;
         }
       }
       if ( fafp_top != 0 ) { fafp_last[current_jet+j_lept]=fafp_top; }
     }

        // ! loop until either stable descendent found

     while ( k[lept][2] != 94 && k[lept][1] != 1 && k[lept][4] == k[lept][5] ) {
       lept = k[lept][4] ;
     }
   
     streamlog_out(DEBUG2)<< " after looping to stable : " << lept << " " <<  k[lept][2] << std::endl;

     while ( true ) {
       int ida=0;
       bool is_elementon = true ;
       if ( k[lept][1] != 1 && ( k[lept][4] == k[lept][5] &&  k[lept][4] != 0 )) {
         ida=k[lept][4] ; 
         if ( k[ida][2] == k[lept][2] ) {
	   is_elementon = false;
         }
       }

       if ( is_elementon ) {

         jet[lept]=current_jet+j_lept ;
         elementon[current_jet+j_lept]=lept ;
	 break;
       } else {

         lept=ida ;

       }

     }

   }  

     // ! set jet for all further descendants of hard leptons . Don't change already set assignments. Never assign a jet to a 94.

   for ( int k_py=1; k_py<=nlund ; k_py++ ) {
     if ( jet[k_py] == 0 && k[k_py][2] != 94 && abs(jet[k[k_py][3]]) > current_jet && abs(jet[k[k_py][3]]) <= current_jet+n_hard_lepton ) {

         // ! this is indeed a descendant of a hard lepton 
       if ( k[k_py][1] == 1 ) {
         jet[k_py] = abs(jet[k[k_py][3]]);
         if ( k[k_py][2] == 22 && k[k[k_py][3] ][2] == k[elementon[jet[k_py]]][2] ) {
           nfsr[jet[k_py]]++;
         } 
       } else { 
         jet[k_py] = jet[k[k_py][3]];
       } 
     }
   }
   current_jet=current_jet+n_hard_lepton ;

}

void TrueJet::cluster()
{
   int clus[4011]={0};
   int clu=0 , jet1=0, jet2=0;

   nclu=0;
   for (int k_py=1; k_py<=nlund ; k_py++ ) {
     streamlog_out(DEBUG0) << " Checking cluster at " << k_py << " k[k_py][2], and k[k_py][1] : " <<   
          " " << k[k_py][2] <<    " , " << k[k_py][1] <<                        std::endl;
     if ( k[k_py][2] == 91 && k[k_py][1] < 30  ) {
       streamlog_out(DEBUG4) << " cluster found at " << k_py <<std::endl;
       nclu++;
       clus[nclu]=k_py;
     }
   }
   if ( nclu == 0 ) return;

   for ( int i_clu=1 ; i_clu<=nclu ; i_clu++ ) {
     clu=clus[i_clu];

       //! Jet-assignment: Decide if the particle should be
       //! assigned to the first or last quark end of the cluster, by checking
       //! the sign of the projection of the hadron momentum on the difference
       //! between the two quark directions.

       //! First of all: find first and second quark going into the cluster, normally the first is
       //! the mother of the cluster, and the second is the line after, or possibly
       //! a few linese later.

     jet1=current_jet+2*(i_clu-1)+1 ; jet2= current_jet+2*(i_clu-1)+2;
     elementon[jet1]=k[clu][3];
     elementon[jet2]=elementon[jet1];
     while ( k[elementon[jet2]+1][4] == clu ) {
       elementon[jet2]=elementon[jet2]+1;
     }
     streamlog_out(DEBUG4) << " Assigning jets " << jet1 << " or " << jet2 
                             << " to CLUSTER at " << clu << std::endl;
     assign_jet(jet1,jet2,clu) ;

   }
   current_jet=current_jet+ 2*nclu;

}

void TrueJet::photon()
{
   int phots[4011]={0};

   nphot=0;
   for (int k_py=1; k_py<=nlund ; k_py++ ) {
     streamlog_out(DEBUG0) << " Checking photon at " << k_py << " k[k_py][2], and k[k_py][1] : " <<   
          " " << k[k_py][2] <<    " , " << k[k_py][1] <<                        std::endl;
     if ( k[k_py][2] == 22 && k[k_py][1] < 30  ) {
       if ( ( ! _whiz1 ) && (k  [k[k_py][3]][1] == 21 ) ) {
         nphot++;
         phots[nphot]=k_py;
         streamlog_out(DEBUG4) << " photon on line " << k_py  <<" hard according to condition 1 " << std::endl;
       }
       else if (abs(k[k[k_py][3]][2]) == 25 ) {
         nphot++;
         phots[nphot]=k_py;	 
         streamlog_out(DEBUG4) << " photon on line " << k_py  <<" hard according to condition 2 " << std::endl;
       }
     }
   }

   if ( nphot == 0 ) return;

   for ( int i_phot=1 ; i_phot<=nphot ; i_phot++ ) {

      // ! jet assignment - index in the list. -ve if "decayimg", which might be the case
      // ! for explicitly requested gammas.

     elementon[current_jet+i_phot] =  phots[i_phot];
     fafp_last[current_jet+i_phot] =  phots[i_phot];
     if ( k[ phots[i_phot]][1] == 1 ) {
       jet[ phots[i_phot]]=current_jet+i_phot;
     } else {
       jet[ phots[i_phot]]=-(current_jet+i_phot);
     } 
   }
   for ( int  k_py=1; k_py<=nlund ; k_py++ ) {
     if ( jet[k_py] == 0 && abs(jet[k[k_py][3]]) > current_jet && abs(jet[k[k_py][3]]) <= current_jet+nphot ) {
  
          //  ! this is indeed a descendant of a hard photon
  
       jet[k_py] = abs(jet[k[k_py][3]]);
     }
   }

   current_jet=current_jet+ nphot;

}

void TrueJet::isr()
{
   int iisr[4011]={0};

   nisr=0;
   for (int k_py=1; k_py<=nlund ; k_py++ ) {
     if  ( (( ! _whiz1 ) && ( (k[k_py][3] == 3 || k[k_py][3] == 4 ) &&  ( k[k_py][2]==22 ) && ( k[k_py][1] < 20 ) )) ||
           ((  _whiz1  ) && ( (k[k_py][3] == 1 || k[k_py][3] == 3 || k[k_py][3] == 0 ) &&  ( k[k_py][2]==22 ) && 
                              (k[k_py][4]==k[k_py][5]) && (k[k[k_py][4]][2] ==22) && (k[k[k_py][4]][1] != 0)  ) )){
       nisr++;
       iisr[nisr] = k_py ;
     }
   }
   if ( nisr > 0 ) {
     // ! jet assignment - index in the list. -ve if "decayimg", which might be the case
     // ! for explicitly requested gammas.

     for ( int j_isr=1 ; j_isr<=nisr ; j_isr++ ) {
       elementon[current_jet+j_isr] =  iisr[j_isr];
       fafp_last[current_jet+j_isr] =  iisr[j_isr];
       if ( k[ iisr[j_isr]][1] == 1 ) {
         jet[ iisr[j_isr]]=current_jet+j_isr;
       } else {
         jet[ iisr[j_isr]]=-(current_jet+j_isr);
       }
     }
       
     for ( int  k_py=1; k_py<=nlund ; k_py++ ) {
       if ( jet[k_py] == 0 && abs(jet[k[k_py][3]]) > current_jet && abs(jet[k[k_py][3]]) <= current_jet+nisr ) {
  
          //  ! this is indeed a descendant of an isr
  
         jet[k_py] = abs(jet[k[k_py][3]]);
       }
     }
     current_jet=current_jet+nisr;
   }
}

void TrueJet::string()
{
  int n_left=0, sstr=0,lquark=0,lstr=0, istr=0,str_nb=0,jet1=0,jet2=0,istr_previous=0,iback=0;
  int str_index[4011]={0}, rev_index[4011]={0};

     //!! Move jets away to make room for string-induced ones at the lowest indicies

   for ( int k_py=1 ; k_py<=nlund ; k_py++ ) {
     rev_index[k_py]=0;
     if ( jet[k_py] != 0 ) {
       jet[k_py]=(jet[k_py]>0 ? jet[k_py]+2*nstr : jet[k_py]-2*nstr);
     }
   }

   if ( nstr > 0 ) {
     for (int j_jet=njet-2*nstr ; j_jet>=1 ; j_jet-- ) {
       elementon[j_jet+2*nstr]=elementon[j_jet] ; 
       elementon[j_jet]=0;
       nfsr[j_jet+2*nstr]=nfsr[j_jet] ; 
       nfsr[j_jet]=0;
       fafp_last[j_jet+2*nstr]=fafp_last[j_jet];
       fafp_last[j_jet]=0;
       boson[j_jet+2*nstr]=boson[j_jet];
       boson[j_jet]=0;
       fafp_boson[j_jet+2*nstr]=fafp_boson[j_jet];
       fafp_boson[j_jet]=0;
     }
   }

   current_jet=0 ; 

        //! Decode the event record. First find the line of the first string,
        //! the last quark of the last string (lquark, at the line before the first string),
        //! the last string (lstr, the daughter of the lquark quark), and the number of strings (nstr).
  

   for ( int k_py=1 ; k_py<= nlund ; k_py++ ) {
     streamlog_out(DEBUG0) << " particle " << k_py << " with PDG " << k[k_py][2] 
                           << " is assigned to jet " << jet[k_py] << std::endl;

     if (jet[k_py]==0 && k[k_py][2] != 91  && k[k_py][2] != 94 &&  k[k_py][2] != 21 &&  k[k_py][2] != 23 &&  
           k[k_py][2] != 24 &&  k[k_py][2] != 25 && k[k_py][1] < 30 ){ // only look for fermions or gammas (not already 
                                                                       // assigned to jets because they are leptons, from 
                                                                       // clusters, isr, whatever...) in the P.S.
       n_left++;
       str_index[n_left]=k_py;
       rev_index[k_py]=n_left;
       if ( sstr == 0 &&  k[k_py][2] == 92 ) {
         sstr=n_left;
       }
     }
   }     


     // sstr is the index in str_index of the first 92 found, and lquark is the quark on the line before that.
     // This is the last quark of the last string in the event.

   lquark=str_index[sstr-1];

   if (k[k[lquark][4]][2] != 92 ) { // there might be an intruder between the lquark and the first string. This is taken care of here
     if ( ! _higgs_to_glue_glue ) {
       int lll = lquark  ;
       while ( lll >= k[str_index[sstr]][3] ) {
         if ( k[lll][4] == str_index[sstr] ) {
           break ;
         } else {
           lll-- ;
         }
       }
       if ( lll >   k[str_index[sstr]][3] ) {
         lquark = lll ;
       } else { 
         // FIXME  : this happens in H->glueglue, Need radical re-thinking for that.
         streamlog_out(ERROR) << " INFO: Non-contigous string ancestors (1)" << std::endl;
         streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ;
         streamlog_out(ERROR) << std::endl ;
       }
     }  
   }


     // lstr is the last string and is the daughter of the last quark of the last string

   lstr=k[lquark][4] ;  istr=lstr ; str_nb = nstr ;
   streamlog_out(DEBUG3) << "  sstr, lstr, lquark, current_jet, str_nb " << sstr << " " << lstr << " " 
                         << lquark  << " " << current_jet  << " " << str_nb << std::endl;

        //! now the jet-assignment. Start by looking at the last string (lstr),
        //! stop once the first has been reached. decide if the particle should be
        //! assigned to the first or last quark end of the string, by checking
        //! the sign of the projection of the hadron momentum on the difference
        //! between the two quark directions.
  
   while( istr >= str_index[sstr]) {

     streamlog_out(DEBUG1) << " calculating jet1/2 from current_jet = " << current_jet << ", str_nb = " << str_nb  << std::endl;
     streamlog_out(DEBUG1) << " istr = " << istr << ", k[istr][3] = " << k[istr][3]  << ", lquark = " << lquark << std::endl;

     jet1= current_jet+2*(str_nb-1)+1 ; jet2=current_jet+2*(str_nb-1)+2 ; 

     elementon[jet1]=k[istr][3] ; elementon[jet2]=lquark ; // elementons: the last quarks before hadronisation, the ones
                                                           // jet-numbering pivots around.

     streamlog_out(DEBUG4) << " Assigning jets " << jet1 << " or " << jet2 << " to STRING at " << istr << std::endl;

     assign_jet(jet1,jet2,istr) ;  // assigns jets *both* to ancestors and decendants of the elementons

     streamlog_out(DEBUG4) << " Assigned jets " << std::endl; 

        //! the previous string is the mother of its hadrons, eg. the last one, which is
        //! on the line before the current string:
     
        //! new string:
  
     istr_previous=istr ;
     iback=str_index[rev_index[istr]-1]  ; istr = k[iback][3] ;
     if ( istr < str_index[sstr] ) break ;   // before the py-line of the first string -> done
  
        //! back-up: end-quark of the previous string is the line before the start-quark
        //! of the present string:

     int lquark_start = str_index[rev_index[k[istr_previous][3]]-1]  ;
     lquark=str_index[rev_index[k[istr_previous][3]]-1]  ;

     if ( k[lquark][4] != istr && lquark > 0 ) {

         //!! Once again the very unusual case: need to search for the end of the new string backwards
         //!! print a message and the event for now. This one actually happens now and then: If there
         //!! is a cluster in the event, there might be a part of the parton shower after the direct parents
         //!! of one string that eventually ends up in a later _string_ (not a cluster, that would have
         //!! been taken care of),
       streamlog_out(DEBUG3) << "  Non-contigous string ancestors (2) " << std:: endl;
       while (  lquark>-1 && k[lquark][4] != istr ) {
         lquark=lquark-1;
       }

       if ( lquark <= 0 ) {
         lquark = lquark_start ;
         while (  lquark<nlund && ( k[lquark][4] != istr || lquark == k[istr][3] )) {
           lquark=lquark+1;
         }
         if ( lquark >=  nlund ) {
           streamlog_out(ERROR) << " ERROR: Quack ?! " << std::endl ;         
           streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ;
           streamlog_out(ERROR) << std::endl ; 
         }
       }
     }
 
        //! decrement the string-counter
  
     str_nb--;
  
   }
}

void TrueJet::assign_jet(int jet1,int jet2,int this_string)
{
   double dir_diff[4]={0};
   int first_p1=0, first_p2=0;
   streamlog_out(DEBUG4) << "   assign_jet: elementon[" << jet1 << "] = " << elementon[jet1] 
                                  << ", elementon[" << jet2 << "] = " << elementon[jet2] << std::endl;
  
    // elementon is the last quark before hadronisation (i.e. either of the two going into to 92 objects).
    // jet1/2 is the jet-number of that one. The output is first_p1/2 which is the first quark in the
    // chain leading up to the elementon. fafp_last is the last quark before any 94 or boson (last in the sense
    // of going back the history!)  - often the same as first_p, while nfsr is the number of FSRs in the jet, 
    // boson is the line of any boson in the P.S. yieding a qqbar pair (typically a gluon), and fafp_boson is 
    // the ultimate ancestor of the boson.
  
   first_parton( elementon[jet1], jet1 , first_p1, fafp_last[jet1], nfsr[jet1], boson[jet1], fafp_boson[jet1]);
   streamlog_out(DEBUG4) << " back in assign_jets for jet " << jet1 << " fafp_last = " <<  fafp_last[jet1] 
                         << " nfsr = " << nfsr[jet1]<< " boson = " <<  boson[jet1]
                         << " fafp_boson = " <<  fafp_boson[jet1] << std::endl;

   first_parton( elementon[jet2], jet2 , first_p2, fafp_last[jet2], nfsr[jet2], boson[jet2], fafp_boson[jet2]);
   streamlog_out(DEBUG4) << " back in assign_jets for jet " << jet2 << " fafp_last = " <<  fafp_last[jet2] 
                         << " nfsr = " << nfsr[jet2]<< " boson = " <<  boson[jet2]
                         << " fafp_boson = " <<  fafp_boson[jet2] << std::endl;

     //   !  difference between the two quark directions: one is at elementon(jet1), the
     //   !  other is elementon(jet2) = three lines of fortran, 20 lines of C++ ...
     //   !  The closest match will determine which elementon the final hadrons are assigned to
     //   !  and hence the jet-numbering.


   double absp_1=0;
   double absp_2=0;

   for ( int kk=1; kk<=3 ; kk++ ) {
     absp_1+=p[elementon[jet1]][kk]*p[elementon[jet1]][kk] ;
     absp_2+=p[elementon[jet2]][kk]*p[elementon[jet2]][kk] ;
   }

   absp_1=sqrt(absp_1); 
   absp_2=sqrt(absp_2); 
   for ( int kk=1; kk<=3 ; kk++ ) {
     dir_diff[kk]      = 
       p[elementon[jet2]][kk]/absp_2-  p[elementon[jet1]][kk]/absp_1;
   }

   for (int k_py_had=k[this_string][4]; k_py_had<=k[this_string][5] ; k_py_had++ ) { // this loops over the first generation hadrons
     if ( k[k_py_had][3] == this_string ) {
       double dot=0;
       for ( int jj=1 ; jj<=3 ; jj++ ) {
         dot+=dir_diff[jj]*p[k_py_had][jj] ;
       }
       if ( dot <= 0. ) {  // by projcting on the *difference* between the two elelementon directions,
                           // one just need to check the sign to know which is closest!
         jet[k_py_had] = jet1 ;
       } else {
         jet[k_py_had] = jet2 ;
       } 
     }
   } 

   for ( int k_py=k[elementon[jet2]][5]+1; k_py <=nlund ; k_py++ ) { // k[elementon[jet2]][5]+1 is the first
                                                                     // hadron of the string (actually any of
                                                                     // elementon[jet1/jet2][4/5] would do...)

     if ( jet[k_py] == 0 && k[k_py][1] < 30 ) { // if not already assigned, and not overlay. NB that the first
                                                // generation already was assigned above.
       jet[k_py] = abs(jet[k[k_py][3]]);        // i.e. jet is the same as it's mother - fair enough !
       if (  k[k_py][2] == 21 || abs( k[k_py][2])<=6 ) {
               //  can happen if a paricle decay is done by gluons (Ypsilon etc.)
         jet[k_py] = -jet[k_py];
       }
       if ( jet[k_py] != 0 ) streamlog_out(DEBUG1) << "     Particle " << k_py 
              << " assigned to jet " << abs(jet[k[k_py][3]]) << std::endl;
     }
   }

}

void TrueJet::first_parton(int this_partic,int this_jet,int& first_partic,int& last_94_parent,int& nfsr_here,int& info,int& info2)
{
   int this_quark=0,quark_pdg=0,boson_ancestor=0,istat=0,boson_last_94_parent=0,istat2=0;

   this_quark=this_partic;
   jet[this_quark]=-this_jet;
   quark_pdg=k[this_quark][2];
   
   streamlog_out(DEBUG4) << "   first_parton: this_quark = " << this_quark 
                         << ",  quark_pdg = " <<  quark_pdg << std::endl;
                     

   last_94_parent =  0;

      //! back-track quarks to the begining of the event record, or until the quark-chain is broken
      //! Beginning of record is hit when the parent of a quark is an electron, photon or a whizard inserted boson 
      //! and the chain is brocken if the parent of a quark is a gluon (21) or a W (24, non-resonanc-insertion).

      //! start at this_quark, which might on entry actually be something else (a gluon or a W)

   while (  k[this_quark][3] > 0 &&
            abs(k[k[this_quark][3]][2]) != 11 &&  abs(k[k[this_quark][3]][2]) != 21    &&  abs(k[k[this_quark][3]][2]) != 22 
        &&  abs(k[k[this_quark][3]][2]) != 23 &&  abs(k[k[this_quark][3]][2]) != 24    &&  abs(k[k[this_quark][3]][2]) != 25 ) {
               // (this end-condition works both for whiz1 and whiz2)

     streamlog_out(DEBUG3) << "     start of while loop:  this_quark = " << this_quark 
                           << ", its PDG = " << k[this_quark][2] 
                           << ", its parent = " << k[this_quark][3] 
                           << ", parent's PDG  = " << abs(k[k[this_quark][3]][2])  
                           << ", quark_pdg  = " << quark_pdg << std::endl;

     //! step back in the history

     this_quark=k[this_quark][3];
     if ( quark_pdg == 21 || abs(quark_pdg) == 24 ) {
          //! looking at a gluon or a W. We've alredy steped to the
          //! parent of the boson (which should be a quark), so now we might 
          //! be loking at a different flavour:
          //! 

       quark_pdg=k[this_quark][2];

     } 

	    //! jet-assignment if requested. Set to -ve to indicate that we are in the
	    //! parton shower,

     if ( this_jet != 0 ) {jet[this_quark]=-this_jet;}

	    //! Internal brems ? If the line after this quark is a photon, and this quark and
            //! the photon has the same mother, it is an FSR

     if  ( k[this_quark+1][2] == 22 &&  k[this_quark+1][3] ==  k[this_quark][3] ) {

       streamlog_out(DEBUG5) << " FSR from  " << this_quark << std::endl;
       if ( this_jet != 0 ) {

             //! assign the FST to the jet of the quark. It's a real particle, so it
             //! assigned a +ve jet number

         jet[this_quark+1]=this_jet;
         nfsr_here=nfsr_here+1;
       }
     } 

     streamlog_out(DEBUG3) << "     end of while loop: this_quark = " << this_quark 
                           << ", its PDG = " << k[this_quark][2] 
                           << ", its parent = " << k[this_quark][3] 
                           << ", parent's PDG  = " << abs(k[k[this_quark][3]][2])  
                           << ", quark_pdg  = " << quark_pdg << std::endl;
   }

     // check if a boson parent is a "real" boson or a resonance insertion, which is flagged by the
     // fact that the parent of the boson is a beam-particle (electron or photon), with status code 21,
     // Resonance insertion were never done in Whiz1, so it is no problem that statuscode 21 is not
     // assigned for Whiz1 - the condition should always evaluate to false, anyhow!
     //
     // We also stop in the parent of the boson is a higgs, both for Whiz1 and 2.
   int resonance_insertion = 0;
   streamlog_out(DEBUG3) << "     resonance_insertion check : " <<  abs(k[k[this_quark][3]][2]) << std::endl;
   if ( abs(k[k[this_quark][3]][2]) == 23  || abs(k[k[this_quark][3]][2]) == 24 ) {
     int lanc = k[k[this_quark][3]][3];
     streamlog_out(DEBUG1) << "     lanc = " << lanc << " " << k[lanc][2] << " " << k[lanc][1] << std::endl ;
     if ( ( ( abs(k[lanc][2]) == 11 || abs(k[lanc][2]) == 22 ) && k[lanc][1] == 21 ) ||
          ( abs(k[lanc][2]) == 25    ) ){ // This is not a "real" boson, but a technical resonance insert
       resonance_insertion = 1;
     }
   }

   streamlog_out(DEBUG3) << "     resonance_insertion = " << resonance_insertion << std::endl;

     //! we are here because we've reached the end of the chain. If this happened
     //! becasue the parent was a boson radiated during the PS, we recurse. (Else we arrived at the
     //! initial state and are done.)

   if ( resonance_insertion == 0 && ( abs(k[k[this_quark][3]][2]) == 21  || abs(k[k[this_quark][3]][2]) == 24 )) {
       //! recurse.

       //! jet-number set to 0 => do not assign jet numbers anymore. The untimate first
       //! partic will go into boson_ancestor, parent of the last 94 into boson_last_94_parent.

     streamlog_out(DEBUG3) << "     in while loop calling first_parton with: k[this_quark][3] = " << k[this_quark][3] << std::endl;
     first_parton(k[this_quark][3],0 ,boson_ancestor, boson_last_94_parent, nfsr_here, istat,istat2);

       //! the way this bottoms-out: second to last argument (istat in the call) is set to 0 if
       //! we didn't find a boson.. In that case ...

     if ( istat == 0 ) {

          //! ... we end up here. Here we set the argument to something different from 0
          //! (namely boson_ancestor), so whoever called me in this step will ...

       info= boson_ancestor;
       info2= boson_last_94_parent;

     } else {

         //! ... end up here, where the same information is further transmitted up
         //! the stack.

       info = istat;
       info2 = istat2;

     }
     //! ... and so on until we are back at the top of the stack.

   } else {

        //! bottom-out

     info=0;
     info2=0;

   }

   if ( last_94_parent == 0 )  {last_94_parent=this_quark;}

   first_partic = this_quark;
      streamlog_out(DEBUG4) << "   first_parton, return:  = first_partic " << first_partic 
                         << ",  last_94_parent = " <<  last_94_parent << std::endl;
                     

}

int TrueJet::flavour(int k2)
{
  // note that flavour is a bit ambigous: for quarks it is the type of the quark (6 values in the SM),
  // for leptons it is the family of the lepton (3 valuse in the SM). Here we mean "family" in both cases.
   int k2l=0;

   if ( abs(k2)>16) return 0 ;
   if ( abs(k2) > 10 ) {
     k2l=abs(k2)-11;
   } else {
     k2l=abs(k2)-1;
   }
    return k2 > 0 ? k2l/2  : k2l/2 ;
    // sign(k2l/2,k2);
}

void TrueJet::grouping(){
   int ijet=0 ;
   for ( int k_jet=1 ; k_jet<=njet ; k_jet++ ) {
     tE[k_jet]=0;
     for ( int i=0 ; i<3 ; i++ ) {
          tmom[k_jet][i]=0. ;
     } 
   }
   for ( int i_py=1 ; i_py<=nlund ; i_py++){
     if ( jet[i_py]>0 ) {
       if ( k[i_py][1] == 11 ) {
         double Ekid=0;
         for ( int jj=k[i_py][4] ; jj <= k[i_py][5] ; jj++ ) {
           Ekid+= p[jj][4];
         }
         if ( abs(Ekid-p[i_py][4])/p[i_py][4] > 0.001 ) {
           streamlog_out(DEBUG3) << " Particle " << i_py << " has energy " << p[i_py][4] << 
                 " , but the sum of its genstat1 kids is " << Ekid << std::endl;
           streamlog_out(DEBUG3) << " indicating that Geant did something fishy and un-documented in MCParticle " 
                << std::endl ;
           streamlog_out(DEBUG3) << " (Counter-meassures taken, so it should be OK.) " << std::endl ;
           streamlog_out(DEBUG3) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl;
           streamlog_out(DEBUG3) << std::endl ;
	   const double* v = mcp_pyjets[k[i_py][4]]->getVertex();
           if ( v[0]*v[0] + v[1]*v[1] > 10. ) {
             if ( k[i_py][1] < 30 ) {
	         k[i_py][1]=1;
             }
             for ( int j_dau=k[i_py][4] ; j_dau <= k[i_py][5] ; j_dau++ ) {

               if ( k[j_dau][1] < 30 ) {
		     k[j_dau][1]=0;
               }
             }
           }
         }
         
       } else if  ( k[i_py][1] == 0 ) {
         for ( unsigned i_dau=0 ; i_dau <   mcp_pyjets[i_py]->getDaughters().size()  ; i_dau++ ) {
           int j_dau=mcp_pyjets[i_py]->getDaughters()[i_dau]->ext<MCPpyjet>();
           if ( k[j_dau][1] < 30 ) {
             k[j_dau][1]=0;
           }
         }
         
       }
     }
   }
   for ( int i_py=1 ; i_py<=nlund ; i_py++){ 
     ijet=abs(jet[i_py]);

     if ( k[i_py][1]%30 == 1 ) { // true quanaties for this jet (all true stable) NB. Checks HEPEVT code,
                             // not wether there are daughters. This is right, since Mokka/DDsim does not
                             // change the HEPEVT status in MCParticle, even if an interaction with
                             // seeable daughters takes place. To check: "un-used daughters", ie.
                             // the opposite case, where the generator has decayed a particle,
                             // but the decay happens after Mokka/DDsim has destroyed the parent.
                             // Also: generator-decays of longlived charged-particles: probaly
                             // (chek that) the MCParticle info is after application of the
                             // B-field, so the direction of the momentum has changed.
                             // [Finally: the crossing-angle. To what has it been applied? 
                             //  everything it seems, so that should be OK]

         tmom[ijet][0]+=p[i_py][1] ;
         tmom[ijet][1]+=p[i_py][2] ;
         tmom[ijet][2]+=p[i_py][3] ;
         tE[ijet]+=p[i_py][4] ;

     }
   }
      //!!   Add info on jet-type to return (1=from string,2=lepton,3=from cluster,4=isr,5=overlay,6=M.E. photon.
      //!!   (+ x00=comes from boson, boson from jet x) and indication of the companion. 


   for ( int k_jet=1 ; k_jet <= njet ; k_jet++ ) {
     if ( k_jet <= 2*nstr ){
       type[k_jet]=1;
     } else if ( k_jet<= 2*nstr+2*nclu ){
       type[k_jet]=3;
     } else if ( k_jet<=  2*nstr+2*nclu+n_hard_lepton ){
       type[k_jet]=2;
     } else if ( k_jet<=  2*nstr+2*nclu+n_hard_lepton+nphot ){
       type[k_jet]=6;
     } else if ( k_jet<= 2*nstr+n_hard_lepton+2*nclu+nphot+nisr ){
       type[k_jet]=4;
     } else if ( k_jet<= njet ){
       type[k_jet]=5;
     } 
     if ( boson[k_jet] != 0 ) {
       type[k_jet]=abs(jet[boson[k_jet]])*100+type[k_jet];
     }
     if ( fafp_boson[k_jet] != 0 ) {
       fafp[k_jet]=fafp_boson[k_jet];
     } else {
       fafp[k_jet]=fafp_last[k_jet];
     }
   }

   for  (int k_jet=1 ; k_jet<=njet-nisr-n_beam_jet-nphot-n_hard_lepton%2 ; k_jet++) {
     if ( k_jet%2 == 0 ) {
       companion[k_jet]=k_jet-1;
     } else {
       companion[k_jet]=k_jet+1;
     }
   }


   // TODO: the following would be much simpler if the groups would be sets, not lists: then the
   // order of the entries would not be important, and one would not need to explictly add
   // k_jet to it's own group. In fact, except for top-events , any jets is always grouped with
   // it's companion, so all what the following code is doing is to check is the initial c.n.
   // jets connect to a jet which is not already it's companion...
   // [In top events, where the initial b:s from the t decay are each-others companions, but should
   // not be grouped - rather we want to group jets from the same top together into an "initial
   // colour-neutral" - a misnomer in this case.]

   for ( int k_jet=1 ; k_jet<= njet ; k_jet++ ) {
     group[0][k_jet]=0;
   }

   for ( int k_jet=1 ; k_jet<= njet ; k_jet++ ) {
     if ( type[k_jet] == 4 || type[k_jet] == 6  ) { // trivial grouping of photons
       group[0][k_jet]=1 ;  group[1][k_jet]=k_jet ;
     }
     
     if ( type[k_jet]%100 <= 3 ) {  // only group strings, clusters and hard leptons

       bool k_lept=(abs(k[fafp[k_jet]][2])>=11 && abs(k[fafp[k_jet]][2])<=16);

       streamlog_out(DEBUG4) << std::endl;
 
       for ( int j_jet=1 ; j_jet<= njet ; j_jet++ ) {
         if ( type[j_jet]%100 <= 3 ) {  

           bool j_lept=(abs(k[fafp[j_jet]][2])>=11 && abs(k[fafp[j_jet]][2])<=16);

	   // the reasons to group k_jet with j_jet:

           bool colour_connected =  ((  k[fafp[j_jet]][2]<0 && k[fafp[k_jet]][2]>0 && 
                                          k[fafp[j_jet]][8] == k[fafp[k_jet]][7] && k[fafp[j_jet]][8] != 0  ) ||  
                                     (  k[fafp[j_jet]][2]>0 && k[fafp[k_jet]][2]<0 && 
                                          k[fafp[j_jet]][7] == k[fafp[k_jet]][8] && k[fafp[j_jet]][7] != 0 ));
           bool no_colour = ( k[fafp[j_jet]][7]+ k[fafp[j_jet]][8]+ k[fafp[k_jet]][7]+ k[fafp[k_jet]][8] == 0 );
          
           bool boson_connected =  (k[fafp[j_jet]][3]==k[fafp[k_jet]][3] && 
                                        (  abs(k[k[fafp[k_jet]][3]][2]) == 6  || 
                                               k[k[fafp[k_jet]][3]][2]  == 23 ||   
                                           abs(k[k[fafp[k_jet]][3]][2]) == 24 || 
					   k[k[fafp[k_jet]][3]][2]  == 25 ));
           bool same_parent = ( k[fafp[j_jet]][3]==k[fafp[k_jet]][3] && k[fafp[k_jet]][3] != 0 );
           bool same_family = (  flavour( k[fafp[k_jet]][2]) == flavour( k[fafp[j_jet]][2])) ;

 	   if (k_lept == j_lept || _top_event ) { // only group quarks with quarks, leptons with leptons, except in top events,
                                                  // where we find it more useful to group all decay-products of the tops together,
                                                  // which might imply grouping quarks with leptons. Note that the concept of
                                                  // "initial colour neutral" is therefore not really apt in this case!

             if (( fafp[j_jet] == fafp[k_jet] ) ||         // the fafp:s are the same or have the same parent. NB will group k_jet with itself!
                   colour_connected || 
                   boson_connected ||
		   (no_colour && (same_parent && same_family)) || 
                   (companion[k_jet] == j_jet  && abs(  k[fafp[k_jet]][2]) != 6) ) {  

               group[0][k_jet]++;
               group[ group[0][k_jet] ][k_jet] = j_jet;

	       if ( k_jet != j_jet ) {
	         streamlog_out(DEBUG4) << " Jet " << k_jet << " was grouped with jet " << j_jet <<" .";
	         streamlog_out(DEBUG4) << " fafp of the jets are " <<  fafp[k_jet] << " (pdg " <<  k[fafp[k_jet]][2] 
                                                    <<  ") and  " <<  fafp[j_jet]  <<" (pdg " <<  k[fafp[j_jet]][2] <<" ) ." ;
	         if ( k[fafp[k_jet]][3] != 0 ) {
	           streamlog_out(DEBUG4) << " parents of the fafp of the jets are " << k[fafp[k_jet]][3]<< " and  " <<   k[fafp[j_jet]][3] <<" .";
	         }
	         streamlog_out(DEBUG4) << " colour-connection : " <<  colour_connected<< std::endl;
               }
             }
           }
	 }
       }
       streamlog_out(DEBUG4) << " group of jet " << k_jet << " : " ;
       bool compthere=false;
       for ( int i_grp=1 ; i_grp <=  group[0][k_jet] ; i_grp++ ) {
         streamlog_out(DEBUG4) <<group[i_grp][k_jet] << " " ;
         if (group[i_grp][k_jet] == companion[k_jet] || companion[k_jet] == 0 )   { compthere = true ; }
       }
       if ( ! compthere && not _top_event ) {
         streamlog_out(ERROR) << " jet " << k_jet << " : companion " << companion[k_jet] << " not in group " << std::endl;  
       }
       streamlog_out(DEBUG4) << std::endl;
       if ( group[0][k_jet] <= 1 && ( type[k_jet] != 2 || n_hard_lepton%2 == 0 ) ) {
         streamlog_out(ERROR) << " jet " << k_jet << " ungrouped " << std::endl;
       }
     }
   }

   n_djb = 0 ;       
   n_dje = 0;  

       //! figure out the different jet combinations. (If there are no bosons, it's straight forward:
       //! it's just the odd+even combinations. Also if one only cares about the final singlet
       //! the same aplies, boson or not)

       // We want to find how many colour-neutrals we have before the P.S., and which jets they
       // go into. dijet_begining[k_jet] will be the number of that c.n., and conversely, jets_begin[initial cn#]
       // is the list of the jets each initial cn goes to.
       // Likewise (but trivially) dijet_end[k_jet] is the number of the final cn the jet, belongs to -
       // jets_end[final cn#] is the list (which always, by construction, contains two adjecent jets).

       // TODO: the below would be simpler if group would be sets. Then just set-equality would be needed to check
       //  The logic below needs identical groups not only to contain the same numbers, but also that they are in the same
       //  order!

   for ( int k_jet=1; k_jet<=njet ; k_jet++ ) {

     dijet_begining[k_jet] = 0;

     if ( group[0][k_jet] > 0 )  {
       for ( int  i_djb=1; i_djb<=n_djb ; i_djb++ ) {
         int gotit=1;
         for ( int i_grpmbr=1 ; i_grpmbr<= group[0][k_jet] ; i_grpmbr++ ) {
           if  ( ! (group[i_grpmbr][k_jet]==jets_begin[i_grpmbr][i_djb] )) {
             gotit=0;
           }
         }
         if ( gotit ) {
           dijet_begining[k_jet] = i_djb;
         }
       }
     }

     if ( dijet_begining[k_jet] == 0 && group[0][k_jet] > 0 ) {
       n_djb++;
       dijet_begining[k_jet] = n_djb;
       for ( int i_grpmbr=0 ; i_grpmbr<= group[0][k_jet] ; i_grpmbr++ ){
         jets_begin[i_grpmbr][n_djb] = group[i_grpmbr][k_jet];
       }
     }

     dijet_end[k_jet] = 0;
     for ( int i_dje=1 ; i_dje<=n_dje ; i_dje++ ) {
       if ( std::min(k_jet,companion[k_jet])==jets_end[1][i_dje] &&
           std::max(k_jet,companion[k_jet])==jets_end[2][i_dje] ) {
         dijet_end[k_jet] = i_dje;
       }
     }

     if ( dijet_end[k_jet] == 0 ) {
       n_dje++;
       dijet_end[k_jet] = n_dje; 
       if ( companion[k_jet] > 0 ) {
         jets_end[1][n_dje] = std::min(k_jet,companion[k_jet]);
         jets_end[2][n_dje] = std::max(k_jet,companion[k_jet]);
         jets_end[0][n_dje] = 2 ;
       } else {
         jets_end[1][n_dje] = k_jet;
         jets_end[0][n_dje] = 1;
       }
     }

   }   


     // collect some topology and consitency information:

   nboson=0;n_jetless=0; n_0_E_jets=0; n_2jet_clu=0;
   for ( int k_jet=1 ; k_jet<=njet ; k_jet++ ) {
     if (type[k_jet]/100 != 0 ){
       nboson++;
     }
     if (tE[k_jet]<0.000001 ) {
       n_0_E_jets++;
     }
     if ( jet[k_jet]==0 && k[k_jet][1]==1 ) {
       n_jetless++;
     }
     if ( (type[k_jet]%100== 3) && ( tE[k_jet]> 0.0000001 ) &&
	 ( tE[companion[k_jet]] >  0.0000001)) {
       n_2jet_clu++;
     }
   }

   n_mixed=0;
   for ( int k_jet=1 ; k_jet<=njet-1 ; k_jet+=2 ) {
     if ((type[k_jet]%100==1) &&
	(k[fafp[k_jet]][4] !=  k[fafp[k_jet+1]][4]) &&
         ( !( ( k[fafp[k_jet]][2] > 0 && k[fafp[k_jet]][7] ==  k[fafp[k_jet+1]][8])  || 
              ( k[fafp[k_jet]][2] < 0 && k[fafp[k_jet]][8] ==  k[fafp[k_jet+1]][7] )) ) &&
          (k[fafp[k_jet]][3] !=  k[fafp[k_jet+1]][3]) ) {

       int ii =  fafp[k_jet];
       while (  ii <= nlund ) {
         if ( ii ==  elementon[k_jet] ) {
           break;
         }  else  {
           ii = k[ii][4] ;
           if ( ii >  elementon[k_jet] ) {
             n_mixed++; 
             break ;
           }
         }
       }
     }
   }  
    //   jsum=jets_summary_t(njet,n_djb,n_dje,2*nstr, n_hard_lepton, nboson,  2*nclu,nisr , n_mixed,  n_jetless, n_2jet_clu, n_0_E_jets,jt(1:njet),dijb(1:n_djb),dije(1:n_dje)) 

}

void TrueJet::check( LCEvent *  /*evt*/ ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TrueJet::end(){ 

    //   std::cout << "TrueJet::end()  " << name() 
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;

}
