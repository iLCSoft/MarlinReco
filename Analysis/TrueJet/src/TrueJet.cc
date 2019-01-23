#include "TrueJet.h"
#include <stdlib.h>
#include <math.h>
//#include <cmath>
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

    evt=event;
    streamlog_out(WARNING) << " processing event: " << evt->getEventNumber() 
        << "   in run:  " << evt->getRunNumber() << std::endl ;
    streamlog_out(MESSAGE4) << " =====================================" << std::endl;


    LCCollection* mcpcol = NULL;
    try{
        mcpcol = evt->getCollection( _MCParticleColllectionName );
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) <<  _MCParticleColllectionName  << " collection not available" << std::endl;
        mcpcol = NULL;
    }
    LCCollection* rmclcol = NULL;
    try{
        rmclcol = evt->getCollection( _recoMCTruthLink );
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) << _recoMCTruthLink   << " collection not available" << std::endl;
        rmclcol = NULL;
    }

    if( mcpcol != NULL &&  rmclcol != NULL){

      reltrue = new LCRelationNavigator( rmclcol );

      // recast the MCParticles into a PYJETS look-alike, needed for the logic to work:

 
      //std::cout << " calling getPyjets"<< std::endl;
      mcp_pyjets.reserve(4000);
      //      MCParticleVec  mcp_pyjets;  mcp_pyjets.reserve(4000);
      getPyjets(mcpcol );

      //std::cout << " return from getPyjets"<< std::endl;
	fix94();
      njet=0 ; i_jet=0 ; nstr=0 ;
      for ( int kk=1 ; kk <= 25 ; kk++ ) {
        boson[kk]=0;  
        elementon[kk]=0;
        fafp_last[kk]=0;
        fafp_boson[kk]=0;
        fafp[kk]=0;
        nfsr[kk]=0 ;
      }
      for ( int kk=1 ; kk <= 4000 ; kk++ ) {
        jet[kk]=0;
        companion[kk] = 0;
      }

      true_lepton() ;

      cluster() ;

      isr() ;


       //   nstr=COUNT( ( 92 == k(1:nlund,2) .AND. jet(1:nlund) == 0 ) )
       //   njet = 2*nstr+n_hard_lepton+2*nclu+nisr;
      nstr=0;
      for ( int kk=1 ; kk <= nlund ; kk++ ) {
        if ( k[kk][2] == 92 && jet[kk] == 0 && k[kk][1] < 30  ) {
          nstr++;
        }
      }
      njet = 2*nstr+n_hard_lepton+2*nclu+nisr;

      if ( nstr != 0 ) {
        string();
      } 

      n_beam_jet=0;
      for (int kk=1 ; kk<=nlund ; kk++) {
        if ( k[kk][1] >= 30 ) {
          if ( n_beam_jet == 0 ) {
            n_beam_jet=1;
          }
          jet[kk]=njet+n_beam_jet;
        }
      }

      njet=njet+n_beam_jet;

      grouping();

      // At this point, the job is done. All the rest is printouts, putting things into collections and setting up navigators.
      // Non-debuggung code is beween ******:s and =========:s

      if ( ! ((n_mixed==0) && njet==2*nstr+n_hard_lepton+2*nclu+nisr+n_beam_jet &&
	      n_jetless==0 && n_hard_lepton%2 == 0 )) {
        streamlog_out(ERROR) << " inconsiency in jet finding : " << std::endl;
        streamlog_out(ERROR) << " n_mixed/ njet/ 2*nstr+n_hard_lepton+2*nclu+nisr+n_beam_jet / n_jetless / n_hard_lepton: " << std::endl;
        streamlog_out(ERROR) << n_mixed <<" / " << njet<< " / " <<2*nstr+n_hard_lepton+2*nclu+nisr+n_beam_jet<< " / " << n_jetless << " / " << n_hard_lepton << std::endl;
        streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
        streamlog_out(ERROR) << std::endl ; 
      }
 
      streamlog_out(DEBUG3) << " HEPEVT relation table up with jet assignment "  <<std::endl;
      streamlog_out(DEBUG3) << "    line   status       pdg  parent  first  last      px         py          pz          E             M        jet companion type" << std::endl;
      streamlog_out(DEBUG3) << "                                       daughter                                                                        jet" << std::endl;

      for ( int i=1 ; i <= nlund ; i++){
	streamlog_out(DEBUG3) << std::setw(7) << i << std::setw(7) <<  k[i][1] << 
                 std::setw(12) <<  k[i][2] << std::setw(7) <<  k[i][3] << 
                 std::setw(7) <<  k[i][4] << std::setw(7) << k[i][5] <<
                 std::setw(12) <<  p[i][1] << 
                 std::setw(12) <<  p[i][2] << std::setw(12) <<  p[i][3] << 
                 std::setw(12) <<  p[i][4] << std::setw(12) << p[i][5] <<
 
                 std::setw(7) << jet[i]  <<
                 std::setw(7) << companion[abs(jet[i])] << std::setw(7) <<  type[abs(jet[i])]<<std::endl;
      }
       streamlog_out(DEBUG3) << "       jet     fafp-beg  qrk/lept  fafp-end fafp-boson  type   dijet-beg  dijet-end    boson" << std::endl;
      for ( int i=1 ; i <= njet ; i++){
	streamlog_out(DEBUG3) << std::setw(10) << i << std::setw(10) <<  fafp[i] << 
                 std::setw(10) <<  elementon[i] << std::setw(10) <<  fafp_last[i] << 
                 std::setw(10) <<  fafp_boson[i] << std::setw(10) << type[i] << std::setw(10) << dijet_begining[i]  <<
                 std::setw(10) << dijet_end[i] << std::setw(10) <<  boson[i]<<std::endl;
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

      //============================

      double str_tmom[3], str_mom[3], str_tmomS[3] , str_tE , str_E , str_tES; 
      str_tE=0; str_tES=0; str_E=0;
      for ( int i=0 ; i<3 ; i++ ) {
         str_tmom[i]=0. ; str_tmomS[i]=0. ; str_mom[i]=0. ; 
      } 


      streamlog_out(DEBUG8) << "  Number of jets found : " << njet << std::endl;

      double tmomS[25][3] ;
      double tES[25];
      int pid_type[25] ;
      for ( int ijet=1; ijet<=njet ; ijet++ ) { // jet-loop
        ReconstructedParticleImpl* true_jet = new ReconstructedParticleImpl;
        tES[ijet]=0;
        for ( int i=0 ; i<3 ; i++ ) {
          tmomS[ijet][i]=0 ;
        } 
      
        //*****************************
        pid_type[ijet] = -type[ijet] ;
        //============================

        jet_vec->addElement(true_jet);
	
	streamlog_out(DEBUG1) << std::endl;
	streamlog_out(DEBUG1) << "  Following jet " << ijet <<   " ( "<< true_jet << ")" << std::endl;
        streamlog_out(DEBUG1) << " ============================= "  <<  std::endl;
        streamlog_out(DEBUG1) << " HEPEVT relation table with jet assignment and pfo(s), if any "  <<std::endl;
        streamlog_out(DEBUG1) << "       line    mcp          status    pdg       parent    first     last       jet     pfo(s)" << std::endl;
        streamlog_out(DEBUG1) << "                                                        daughter   daugther " << std::endl;
      }
      
      ReconstructedParticleImpl* true_jet ;
      for ( int i=1 ; i<=nlund ; i++){ // pyjets loop
        int ijet =  abs(jet[i]);
        if ( ijet > 0 ) {

          true_jet = dynamic_cast<ReconstructedParticleImpl*>(jet_vec->getElementAt(ijet-1));

             // OK, add jet-to-true link

	  streamlog_out(DEBUG1) << std::setw(10) << i << " ("  <<mcp_pyjets[i] << ")" << std::setw(10) <<  k[i][1]%30 << 
                 std::setw(10) <<  k[i][2] << std::setw(10) << k[i][3] << 
                 std::setw(10) <<  k[i][4] << std::setw(10) << k[i][5] << std::setw(10) << jet[i]  ;


          //*****************************
          if (  k[i][1] == 1 &&  k[i][2] == 22 && jet[k[i][3]] < 0   && k[k[i][3]][2] != 22 ) {
                // to be able to find FSRs, set weight to -ve for them
            truejet_truepart_Nav.addRelation(  true_jet, mcp_pyjets[i], ( jet[i]>0 ? -(k[i][1])%30 : 0.0 )) ;
          } else {
            truejet_truepart_Nav.addRelation(  true_jet, mcp_pyjets[i], ( jet[i]>0 ? k[i][1]%30 : 0.0 )) ;
          }
          LCObjectVec recovec; 
          recovec = reltrue->getRelatedFromObjects( mcp_pyjets[i]);
          //============================ 

          // (maybe this will be needed: If the current assumption that also MCPartiles created in
          //  simulation can in a consitent way be put into pyjets turns out to be wrong, then
          //  getPyjets should only load gen. stat. /= 0 particles, and then code below
          //  would be needed to assign seen gen. stat. = 0 particles to jets based on their
          //  gen. stat. /= 0 ancestor ?!)
          //            if ( recovec.size() == 0 &&  k[i][1] == 1 ) { // not reconstructed stable particle
          //              if ( abs(k[i][2]) != 12 &&  abs(k[i][2]) != 14 &&  abs(k[i][2]) != 16 ) { // not neutrino
          //                //  mcp_pyjets[i].getDaughter ... etc . If any decendant is seen, make association
          //              }
          //            }

          if ( recovec.size() > 0 ) { // if reconstructed
            double mom2[3] ;  
            double E,q ;

             // true quanaties for this jet (all true seen)
            tmomS[ijet][0]+=p[i][1] ;
            tmomS[ijet][1]+=p[i][2] ;
            tmomS[ijet][2]+=p[i][3] ;
            tES[ijet]+=p[i][4] ;

	    streamlog_out(DEBUG1) << " recovec size is " << recovec.size() ;
            int last_winner = ijet ;
            for ( unsigned ireco=0 ; ireco<recovec.size() ; ireco++ ) { //  reco-of-this-true loop

              // add things of this reco-particle to the current jet:

              //*****************************
              ReconstructedParticle* reco_part  = dynamic_cast<ReconstructedParticle*>(recovec[ireco]);

              if (truejet_pfo_Nav.getRelatedFromObjects(reco_part).size() == 0 ) { // only if not yet used
		streamlog_out(DEBUG3) << " recopart " << ireco ;
		bool split_between_jets = false;
	        int winner , wgt_trk[26] , wgt_clu[26] ;
		winner=ijet;
	        for ( int kkk=1 ; kkk <= njet ; kkk++ ) {
		  wgt_trk[kkk] =0 ; wgt_clu[kkk] = 0 ;
	        }	
         	LCObjectVec recomctrues = reltrue->getRelatedToObjects(reco_part);
                static FloatVec recomctrueweights;
	        recomctrueweights = reltrue->getRelatedToWeights(reco_part);
	        streamlog_out(DEBUG3) << "     mctrues of this reco " << recomctrues.size() << std::endl ;
         	for ( unsigned kkk=0 ; kkk<recomctrues.size() ; kkk++ ) {
		  int jetoftrue ;
		  MCParticle* an_mcp = dynamic_cast<MCParticle*>(recomctrues[kkk]);
		  jetoftrue=jet[an_mcp->ext<MCPpyjet>()];
	          if ( jetoftrue != ijet ) split_between_jets = true;
		  wgt_trk[jetoftrue]+=int(recomctrueweights[kkk])%10000 ;
		  wgt_clu[jetoftrue]+=int(recomctrueweights[kkk])/10000 ;
		  streamlog_out(DEBUG1) << " weights : " << int(recomctrueweights[kkk]) << " " << int(recomctrueweights[kkk])%10000 
					    << " " << int(recomctrueweights[kkk])/10000 << std::endl ;
                }
                if ( split_between_jets ) {
		  double wgt_trk_max, wgt_clu_max ;
	          int imax_trk, imax_clu;

         	  wgt_trk_max=0. ; wgt_clu_max=0. ; imax_trk = 0 ; imax_clu = 0;
		  for ( int kkk=1 ; kkk <= njet ; kkk++ ) {
		    streamlog_out(DEBUG3) << "     jetweights for " << kkk << " " << wgt_trk[kkk]  << " " <<  wgt_clu[kkk] ;
		    if ( wgt_trk[kkk] > wgt_trk_max ) {
	              imax_trk= kkk ; wgt_trk_max = wgt_trk[kkk] ;
         	    }  
		    if ( wgt_clu[kkk] > wgt_clu_max ) {
		      imax_clu= kkk ; wgt_clu_max = wgt_clu[kkk] ;
		    }  
		  }
		  streamlog_out(DEBUG3) << "    " << imax_clu << " " << imax_trk << " " <<  wgt_clu_max <<  " " <<  wgt_trk_max ;
		  if ( imax_clu == imax_trk || wgt_trk_max >= 750 ) {
		    winner = imax_trk;
		  } else if ( wgt_trk_max == 0 ) {
		    winner = imax_clu;
		  } else {
		    streamlog_out(DEBUG3) << " was ? " << imax_clu << " " << imax_trk << " " <<  wgt_clu_max <<  " " <<  wgt_trk_max << std::endl;
		    for ( int kkk=1 ; kkk <= njet ; kkk++ ) {
		      streamlog_out(DEBUG3) << " was    jetweights for " << kkk << " " << wgt_trk[kkk]  << " " <<  wgt_clu[kkk]  << std::endl;
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

		streamlog_out(DEBUG3) << " and the winner is " << winner << "  ( ijet is " << ijet << " ) " ;	  
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
                  streamlog_out(WARNING) << " reco particle " << ireco << " of pythia-particle " << i << " has weights = 0 to ALL MCPs ??? " << std::endl;
                  streamlog_out(WARNING) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                }
              } // end if not yet used
              else{
		streamlog_out(DEBUG1) << " " << reco_part << " seen more than once ! " << std::endl;
              }
            } // end reco-of-this-true loop

	  } // end if reconstructed
  	  streamlog_out(DEBUG1) << std::endl;

        } // end if ijet > 0

      } // end pyjets loop
      for ( int ijet=1; ijet<=njet ; ijet++ ) { // jet-loop

        true_jet = dynamic_cast<ReconstructedParticleImpl*>(jet_vec->getElementAt(ijet-1));
        ParticleIDImpl* pid = new ParticleIDImpl; 
        pid->setPDG(k[elementon[ijet]][2]);
        pid->setType(pid_type[ijet]);
        true_jet->addParticleID(pid);
        true_jet->ext<TJindex>()=ijet; 
      }

      for ( int ijet=1; ijet<=njet ; ijet++ ) { // jet-loop

        true_jet = dynamic_cast<ReconstructedParticleImpl*>(jet_vec->getElementAt(ijet-1));



        const double* mom = true_jet->getMomentum(); 

	streamlog_out(DEBUG5) << std::endl;
	streamlog_out(DEBUG9) << " summary " << ijet << " " << type[ijet]%100 << " " << true_jet->getEnergy() << " " <<  tE[ijet]  << " " <<  tES[ijet] << std::endl ;
	streamlog_out(DEBUG5) << "  jet       " <<std::setw(4)<< ijet << " (seen)         p : " << mom[0] << " " << mom[1] << " "  << mom[2] 
                    << " " << " E " << true_jet->getEnergy()<< std::endl;
	streamlog_out(DEBUG5) << "               "   <<      "  (true)         p : " << tmom[ijet][0] << " " << tmom[ijet][1] << " "  << tmom[ijet][2] 
                    << " " << " E " << tE[ijet] << std::endl;
	streamlog_out(DEBUG5) << "               "   <<      "  (true of seen) p : " << tmomS[ijet][0] << " " << tmomS[ijet][1] << " "  << tmomS[ijet][2] 
                    << " " << " E " << tES[ijet] << std::endl;
	streamlog_out(DEBUG5) << "  ancestor 1" << std::setw(4)<<elementon[ijet]  << "                p : " << p[elementon[ijet]][1] << " " <<  p[elementon[ijet]][2] << " "  <<  p[elementon[ijet]][3] 
	            << " " << " E " << p[elementon[ijet]][4]<< std::endl;
	streamlog_out(DEBUG5) << "  ancestor 2" << std::setw(4)<<fafp_last[ijet] << "                p : " << p[fafp_last[ijet]][1] << " " <<  p[fafp_last[ijet]][2] << " "  <<  p[fafp_last[ijet]][3] 
                    << " " << " E " << p[fafp_last[ijet]][4]<< std::endl;
	streamlog_out(DEBUG5) << "  ancestor 3" << std::setw(4)<<fafp[ijet] << "                p : " << p[fafp[ijet]][1] << " " <<  p[fafp[ijet]][2] << " "  <<  p[fafp[ijet]][3] 
                    << " " << " E " << p[fafp[ijet]][4]<< std::endl;
	streamlog_out(DEBUG5) << std::endl;

        streamlog_out(DEBUG5) << "   Character of jet " << ijet << " : "<<type[ijet]%100 <<" (1=string, 2=lepton , 3=cluster, 4=isr, 5=overlay)."  << std::endl;
        streamlog_out(DEBUG5) << "   Jet from boson ? " << (type[ijet]>=100 ? type[ijet]/100 : 0) << " (0 not from boson, !=0 from boson radiated of that jet)."  << std::endl;
        streamlog_out(DEBUG5) << "   No. of FSRs : " << nfsr[ijet] << std::endl;
        streamlog_out(DEBUG5) << "   Jet has no energy ? " << (tE[ijet]<0.001) << std::endl;

	streamlog_out(DEBUG5) << std::endl;

        str_E+= true_jet->getEnergy(); str_tE+= tE[ijet]; str_tES+= tES[ijet];
        for ( int i=0 ; i<3 ; i++ ) {
          str_mom[i]+=mom[i];
          str_tmom[i]+=tmom[ijet][i];
          str_tmomS[i]+=tmomS[ijet][i]; 
        }
        if ( ijet%2 == 0 ) {
          streamlog_out(DEBUG6) << "      Mass of jet " << ijet-1 << " and " << ijet << std::endl;
          double masssq ;
          masssq = str_E*str_E - str_mom[0]*str_mom[0] - str_mom[1]*str_mom[1] - str_mom[2]*str_mom[2] ;
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
          if ( type[ijet]%100 == 3  ) {
             // cluster - the 91 object does NOT have the mass of the following physical states,
             // so we sum up the descendants - in the vast majority there is only one
            double clu_mass=0 , clu_E= 0., clu_P[4]= {0,0,0};
            for ( int jj=k[k[elementon[ijet]][4]][4] ; jj<=k[k[elementon[ijet]][4]][5] ; jj++ ) {
              clu_E += p[jj][4] ;    clu_P[1] += p[jj][1] ;  clu_P[2] += p[jj][2] ;  clu_P[3] += p[jj][3] ;  
            }
            clu_mass=sqrt(clu_E*clu_E-clu_P[1]*clu_P[1]-clu_P[2]*clu_P[2]-clu_P[3]*clu_P[3] );
            streamlog_out(DEBUG6) << "      string mass :   " <<    clu_mass << std::endl;
	  } else if ( type[ijet]%100 == 2 ) {
            double dilept_mass=0 , dilept_E= 0., dilept_P[4]= {0,0,0};
            for ( int jj=ijet-1 ; jj<=ijet ; jj++ ) {
              dilept_E += p[elementon[jj]][4] ;    dilept_P[1] += p[elementon[jj]][1] ;   dilept_P[2] += p[elementon[jj]][2] ;   dilept_P[3] += p[elementon[jj]][3] ;  
            }
            dilept_mass=sqrt(dilept_E*dilept_E-dilept_P[1]*dilept_P[1]-dilept_P[2]*dilept_P[2]-dilept_P[3]*dilept_P[3] );
            streamlog_out(DEBUG6) << "      string mass :   " <<    dilept_mass << std::endl;
          } else {  
             // string, just take the number from the 92-object
            streamlog_out(DEBUG6) << "      string mass :   " <<    p[k[elementon[ijet]][4]][5]  << std::endl;
          }
    
	  streamlog_out(DEBUG6) << std::endl;
          str_tE=0; str_tES=0; str_E=0;
          for ( int i=0 ; i<3 ; i++ ) {
            str_tmom[i]=0. ; str_tmomS[i]=0. ; str_mom[i]=0. ; 
          } 
          if (  masssq_t > 0. ) {
            if ( abs( sqrt(masssq_t)- p[k[elementon[ijet]][4]][5] )/  p[k[elementon[ijet]][4]][5]  > 0.001 ) {
              if ( type[ijet]%100 == 1 && type[ijet-1]%100 == 1 &&  nfsr[ijet]+nfsr[ijet-1]==0) {
                if (  abs( tE[ijet]+tE[ijet-1]- p[k[elementon[ijet]][4]][4] )/  p[k[elementon[ijet]][4]][4] > 0.001  ) {
	  	  streamlog_out(ERROR) << " bad match M (sum/initial) " << " " <<  sqrt(masssq_t) << " / " << p[k[elementon[ijet]][4]][5] << std::endl;
		  streamlog_out(ERROR) << "           E (sum/initial) " << " " <<  tE[ijet]+tE[ijet-1] << " / " << p[k[elementon[ijet]][4]][4] << std::endl;

                  streamlog_out(DEBUG9) << "list HEPEVT relation table up with jet assignment "  <<std::endl;
                  streamlog_out(DEBUG9) << "list    line   status       pdg  parent  first  last      px         py          pz          E             M        jet companion type" << std::endl;
                  streamlog_out(DEBUG9) << "list                                       daughter                                                                        jet" << std::endl;
                  for ( int i=1 ; i <= nlund ; i++){
	            streamlog_out(DEBUG9) <<"list "<< std::setw(7) << i << std::setw(7) <<  k[i][1] << 
                      std::setw(12) <<  k[i][2] << std::setw(7) <<  k[i][3] << 
                      std::setw(7) <<  k[i][4] << std::setw(7) << k[i][5] <<
                      std::setw(12) <<  p[i][1] << 
                      std::setw(12) <<  p[i][2] << std::setw(12) <<  p[i][3] << 
                      std::setw(12) <<  p[i][4] << std::setw(12) << p[i][5] <<
 
                      std::setw(7) << jet[i]  <<
                      std::setw(7) << companion[abs(jet[i])] << std::setw(7) <<  type[abs(jet[i])]<<std::endl;
                  }
                  streamlog_out(ERROR) << "       jet     fafp-beg  qrk/lept  fafp-end fafp-boson  type   dijet-beg  dijet-end    boson" << std::endl;
                  for ( int i=ijet-1 ; i <= ijet ; i++){
	            streamlog_out(WARNING) << std::setw(10) << i << std::setw(10) <<  fafp[i] << 
                    std::setw(10) <<  elementon[i] << std::setw(10) <<  fafp_last[i] << 
                    std::setw(10) <<  fafp_boson[4] << std::setw(10) << type[i] << std::setw(10) << dijet_begining[i]  <<
	            std::setw(10) << dijet_end[i] << std::setw(10) <<  boson[i]<<std::endl;
                  }
                  streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                  streamlog_out(ERROR) << std::endl ; 
                } else {

		  streamlog_out(WARNING) << " bad match M (sum/initial) " << " " <<  sqrt(masssq_t) << " / " << p[k[elementon[ijet]][4]][5] << std::endl;
		  streamlog_out(WARNING) << "           E (sum/initial) " << " " <<  tE[ijet]+tE[ijet-1] << " / " << p[k[elementon[ijet]][4]][4] << std::endl;
                  streamlog_out(WARNING) << " (As the energy matches well, this is probably due to the B-field) "  <<std::endl;
                  streamlog_out(WARNING) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                  streamlog_out(WARNING) << std::endl;
                }
                for ( int kk=1 ; kk<=nlund ; kk++ ) {
                  if ( jet[kk] == ijet || jet[kk] == ijet -1 ) {
                    if ( k[kk][1] == 11 ) { 
                      double e_kid=0;
                      for ( int jj= k[kk][4] ; jj <=  k[kk][5] ; jj++ ) {
                        e_kid+=p[jj][4];
                      }
                      if ( (abs(e_kid - p[kk][4])/ p[kk][4] ) > 0.001 ) {
                         streamlog_out(WARNING) << kk << " " <<  k[kk][4]  << " " << k[kk][5] << " "  << e_kid << " " << p[kk][4] << std::endl ;
                      }
                    }
                  }
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

      LCRelationNavigator FinalColourNeutral_Nav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::RECONSTRUCTEDPARTICLE ) ;
      LCRelationNavigator FinalElementon_Nav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE ) ;

      for ( int kk=1; kk<=n_dje ; kk++ ) {
        double E=0, M=0 , mom[3]={} ;  int  pdg[26]={};
        if ( type[jets_end[1][kk]]%100 == 5 ) {
          // beam-jet: No ancestor - just use sum of true-stable
          E=tE[jets_end[1][kk]];
          double psqrs=0;
          for (int jj=0 ; jj<3 ; jj++ ) {
            mom[jj] = tmom[jets_end[1][kk]][jj];
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
          for ( int jj=1 ; jj<=jets_end[0][kk] ; jj++ ) {
              // decide which line should be used - this (leptons) - child (strings,gammas) - grandchild (cluster)
              // NB. In the case of leptons, it might be that the first and last leptons are not next to
              // eachother.
            int this_first, this_last;
            if (  type[jets_end[1][kk]]%100 == 2 ) {
                // lepton: need to go one step less as the di-jet generator might
                // be stable, and hence have no descendants.
               this_first=elementon[jets_end[jj][kk]];
               this_last=elementon[jets_end[jj][kk]] ;
            } else if (  type[jets_end[1][kk]]%100 == 3 ) { 
             // cluster: need to go one step further (the 91 object doesn't always 
             // have the mass of the particle it goes to). Sometimes first and
             // last will actually become different - the cluster goes to more than one
             // hadron.
               this_first = k[k[elementon[jets_end[jj][kk]]][4]][4];
               this_last =  k[k[elementon[jets_end[jj][kk]]][5]][5];
            } else {
               this_first = k[elementon[jets_end[jj][kk]]][4] ;
               this_last =  k[elementon[jets_end[jj][kk]]][5];
            }

            if ( this_first < first ) {
              first= this_first;
            }
            if (  this_last > last ) {
              last= this_last;
            }
            pdg[jj]=k[elementon[jets_end[jj][kk]]][2];
            if ( pdg[0] == 0 ) {
              if ( k[elementon[jets_end[jj][kk]]][4] != 0 ) {
                pdg[0]=k[k[elementon[jets_end[jj][kk]]][4]][2];
              } else {
                pdg[0]=k[elementon[jets_end[jj][kk]]][2];
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
            for ( int ii=1 ; ii<=jets_end[0][kk] ; ii++ ) {
              for ( int jj=first ; jj<=last ; jj++ ) {
                if ( abs(jet[jj]) == jets_end[ii][kk] ) {
                  E+=  p[jj][4];
                  for ( int ll=1 ; ll<=3 ; ll++ ) {
                    mom[ll-1] +=    p[jj][ll];
                  }
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
        for(int jj=1 ; jj<=jets_end[0][kk] ; jj++) {
          fafpf->addParticle(dynamic_cast<ReconstructedParticle*>(jet_vec->getElementAt(jets_end[jj][kk]-1)) );
        }
        ParticleIDImpl* pid[26];
        pid[0] = new ParticleIDImpl; 
        pid[0]->setPDG(pdg[0]);
        pid[0]->setType(type[jets_end[1][kk]]%100);  // maybe flag from boson? could be 0,1, or 2 jets in the fafp that's from boson ..
        fafpf->addParticleID(pid[0]);
        for ( int jj=1 ; jj<=jets_end[0][kk] ; jj++ ) {
          pid[jj]= new ParticleIDImpl; 
          pid[jj]->setPDG(pdg[jj]);
          pid[jj]->setType(type[jets_end[jj][kk]]);
          fafpf->addParticleID(pid[jj]);
          if ( elementon[jets_end[jj][kk]] > 0 ) {
            FinalElementon_Nav.addRelation(  jet_vec->getElementAt(jets_end[jj][kk]-1) , mcp_pyjets[elementon[jets_end[jj][kk]]], 1.0 );
          }
          FinalColourNeutral_Nav.addRelation(  jet_vec->getElementAt(jets_end[jj][kk]-1) , fafpf);
        }
        fafpf_vec->addElement(fafpf);
      }
   
           
        // pre-PS part

      LCCollectionVec* fafpi_vec = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE ) ;

      LCRelationNavigator InitialElementon_Nav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE ) ;
      LCRelationNavigator InitialColourNeutral_Nav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::RECONSTRUCTEDPARTICLE ) ;

      for (int kk=1 ; kk<=n_djb ; kk++ ) {
        double E=0, M=0 , mom[3]={} ; 
        int pdg[26]= {};
        //JL assume Z as default, check for W, H later
        int bosonid=23 ;
             // E and P is sum of (unique) direct decendants of the
             // initial ffbar. Normally, there is only one
             // (a 94), so the sum is actualy just copying the descendant. 
             // For leptons, however, there is no 94 - they go directly to
             // the next generation leptons), so the sum is non-trivial.
             // Sometimes this also happens for quraks, ie. they go directly to
             // the genertor quarks of the string.
        int first=nlund , last=0;
        for ( int jj=1 ; jj<=jets_begin[0][kk] ; jj++ ) {
          if ( k [fafp [jets_begin[jj][kk]] ] [4] == 0 ) {
            if (  fafp[jets_begin[jj][kk]] < first ) {
              first= fafp[jets_begin[jj][kk]] ;
            }
          } else {
            if (  k [fafp [jets_begin[jj][kk]] ] [4] < first ) {
              first= k[fafp[jets_begin[jj][kk]]][4] ;
            }
          }
          if ( k [fafp [jets_begin[jj][kk]] ] [5] == 0 ) {
            if (  fafp[jets_begin[jj][kk]] > last ) {
              last= fafp[jets_begin[jj][kk]] ;
            }
          } else { 
            if (  k[fafp[jets_begin[jj][kk]]][5] > last ) {
              last= k[fafp[jets_begin[jj][kk]]][5] ;
            }
          }
	  //          pdg[jj]=k[fafp[jets_begin[jj][kk]]][2];
          pdg[jj]=k[elementon[jets_begin[jj][kk]]][2];
          // TODO: maybe lift this part out into a separate loop ?
          // JL: first check whether we have an explicit mother in the event listing, like for Higgs events
          // MB: only do this if the parent *is* a boson
          if ( k[fafp[jets_begin[jj][kk]]][3] > 0) {
            if (  abs(k [k[fafp[jets_begin[jj][kk]]][3]] [2]) > 21 &&  abs(k [k[fafp[jets_begin[jj][kk]]][3]] [2]) <= 39 ) {
                // ie. any IVB except the gluon
              bosonid = abs(k [k[fafp[jets_begin[jj][kk]]][3]] [2]);
              streamlog_out(DEBUG3) << " elementon " << fafp[jets_begin[jj][kk]] << " has explicit mother with PDG = " << bosonid << std::endl;
            }
          }
          if ( pdg[0] == 0 ) {
              //JL this is _not_ the mother boson pdg yet, but a temporary storage for the pdg of the (first) elementon of this boson
            pdg[0]=abs(k[fafp[jets_begin[jj][kk]]][2]);
              //JL if ID of one of the following elementons of this boson has a different pdg, assume we have a W
          } 
          else if ( bosonid==23 && abs(k[fafp[jets_begin[jj][kk]]][2]) != pdg[0] ) {
            bosonid=24;
          } 
        }
        pdg[0]=bosonid;
        if ( last == first ) {
          E=  p[first][4];
          for ( int ll=1 ; ll<=3 ; ll++ ) {
            mom[ll-1] =    p[first][ll];
          }
        } else {
          // if several, sum up but make sure that there are no intruders between first and last.
          for ( int ii=1 ; ii<=jets_begin[0][kk] ; ii++ ) {
            for ( int jj=first ; jj<=last ; jj++ ) {
             if ( abs(jet[jj]) == jets_begin[ii][kk] && k[jj][3] < first ) {
                E+=  p[jj][4];
                for ( int ll=1 ; ll<=3 ; ll++ ) {
                  mom[ll-1] +=    p[jj][ll];
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
        ParticleIDImpl* pid[26] ;
        pid[0]= new ParticleIDImpl; 
        pid[0]->setPDG(pdg[0]);
        pid[0]->setType(type[jets_begin[1][kk]]%100);  
        fafpi->addParticleID(pid[0]);
        for(int jj=1 ; jj<=jets_begin[0][kk] ; jj++) {
          pid[jj]= new ParticleIDImpl; 
          pid[jj]->setPDG(pdg[jj]);
          pid[jj]->setType(type[jets_begin[jj][kk]]);  
          fafpi->addParticleID(pid[jj]);
          fafpi->addParticle(dynamic_cast<ReconstructedParticle*>(jet_vec->getElementAt(jets_begin[jj][kk]-1)));
        }
        fafpi_vec->addElement(fafpi);
        
        for ( int jj=1 ; jj<=jets_begin[0][kk] ; jj++ ) {
          InitialElementon_Nav.addRelation(   jet_vec->getElementAt(jets_begin[jj][kk]-1), mcp_pyjets[fafp[jets_begin[jj][kk]]], 1.0 );
          InitialColourNeutral_Nav.addRelation(  jet_vec->getElementAt(jets_begin[jj][kk]-1) , fafpi);
        }
   
      }
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
      streamlog_out(DEBUG4) << "    line      status           pdg  parent  first  last      px         py          pz          E             M        jet companion type       PFO/Energy" << std::endl;
      streamlog_out(DEBUG4) << "           init ext'd                         daughter                                                                        jet" << std::endl;

      for ( int i=1 ; i <= nlund ; i++){
        int istat=0;
        if ( jet[i] != 0 ) {
          static FloatVec www;
          www = truejet_truepart_Nav.getRelatedFromWeights(mcp_pyjets[i]);
          istat=int(www[0]);
        }
	streamlog_out(DEBUG4) << std::setw(7) << i << std::setw(7) <<  k[i][1] << std::setw(7) <<  istat << 
                 std::setw(12) <<  k[i][2] << std::setw(7) <<  k[i][3] << 
                 std::setw(7) <<  k[i][4] << std::setw(7) << k[i][5] <<
                 std::setw(12) <<  p[i][1] << 
                 std::setw(12) <<  p[i][2] << std::setw(12) <<  p[i][3] << 
                 std::setw(12) <<  p[i][4] << std::setw(12) << p[i][5] <<
 
                 std::setw(7) << jet[i]  <<
                 std::setw(7) << companion[abs(jet[i])] << std::setw(7) <<  type[abs(jet[i])];
        LCObjectVec recovec; 
        recovec = reltrue->getRelatedFromObjects( mcp_pyjets[i]);
        if ( recovec.size() > 0 ) { // if reconstructed
           for ( unsigned ireco=0 ; ireco<recovec.size() ; ireco++ ) { //  reco-of-this-true loop

             ReconstructedParticle* reco_part  = dynamic_cast<ReconstructedParticle*>(recovec[ireco]);
  	     streamlog_out(DEBUG4) << " [" << std::setw(10) << reco_part <<" /"<< std::setw(11) << reco_part->getEnergy() <<"] ";
           }
        } else if ( istat == 1 ) {
  	  streamlog_out(DEBUG4) << " [" << std::setw(10) << "N.A   "<<" /"<< std::setw(11) << "N.A   "<<"] ";
        }
	streamlog_out(DEBUG4) << std::endl;
      }
      streamlog_out(DEBUG6) <<  std::endl;
      streamlog_out(DEBUG6) << "       jet     fafp-beg  qrk/lept  fafp-end fafp-boson  type   dijet-beg  dijet-end    boson" << std::endl;
      for ( int i=1 ; i <= njet ; i++){
	streamlog_out(DEBUG6) << std::setw(10) << i << std::setw(10) <<  fafp[i] << 
                 std::setw(10) <<  elementon[i] << std::setw(10) <<  fafp_last[i] << 
                 std::setw(10) <<  fafp_boson[i] << std::setw(10) << type[i] << std::setw(10) << dijet_begining[i]  <<
                 std::setw(10) << dijet_end[i] << std::setw(10) <<  boson[i]<<std::endl;
      }
      streamlog_out(DEBUG6) <<  std::endl;
      streamlog_out(DEBUG6) << "       Colour-neutral objects at the beginning of the parton shower " << std::endl;
      streamlog_out(DEBUG6) << "    number    px          py           pz          E           M       type       PDGs        jets " << std::endl;
      int glurad=0;
      for ( unsigned kk=0 ; kk<fafpi_vec->size() ; kk++ ) {
        
        ReconstructedParticle* fafpi = dynamic_cast<ReconstructedParticle*>(fafpi_vec->at(kk));
        LCObjectVec jts;
        jts = InitialColourNeutral_Nav.getRelatedFromObjects(fafpi);

        const double* mom=fafpi->getMomentum();
	streamlog_out(DEBUG6) << std::setw(7) << kk+1 << std::setw(12) << mom[0] << std::setw(12) <<  mom[1] << std::setw(12)  << mom[2] << 
                             std::setw(12)  <<fafpi->getEnergy() << std::setw(12) << fafpi->getMass() << 
                             std::setw(7)  << fafpi->getParticleIDs()[0]->getType() << "   " ;
//<< std::setw(10) << fafpi->getParticleIDs()[0]->getPDG() ;
        for ( unsigned ipid=0 ; ipid < fafpi->getParticleIDs().size() ; ipid++ ) {
	  if (fafpi->getParticleIDs()[ipid]->getType()>=100 ) glurad=1;
          if (ipid < fafpi->getParticleIDs().size()-1 ) {
            streamlog_out(DEBUG6)  << std::setw(3)<< fafpi->getParticleIDs()[ipid]->getPDG()  << "," ;
          } else {
            streamlog_out(DEBUG6)  << std::setw(3)<< fafpi->getParticleIDs()[ipid]->getPDG() << "   " ;
          }
        }
        for ( unsigned jj=0 ; jj<jts.size() ; jj++ ) {
           int ijet=jts[jj]->ext<TJindex>(); 
           if ( jj < jts.size()-1 ) {
             streamlog_out(DEBUG6) << std::setw(3) << ijet << "," ; 
           } else {
             streamlog_out(DEBUG6) << std::setw(3) << ijet  ;
           }
        }
        streamlog_out(DEBUG6) << std::endl;
      }
      if ( glurad == 1 ) {
        streamlog_out(DEBUG6) << std::endl;
        streamlog_out(DEBUG6) <<    "              gluon radiation in event : "<< std::endl;
        for ( unsigned kk=0 ; kk<fafpi_vec->size() ; kk++ ) {
        
          ReconstructedParticle* fafpi = dynamic_cast<ReconstructedParticle*>(fafpi_vec->at(kk));
          LCObjectVec jts;
          jts = InitialColourNeutral_Nav.getRelatedFromObjects(fafpi);

          for ( unsigned jj=1 ; jj < fafpi->getParticleIDs().size() ; jj++ ) {
	    if (fafpi->getParticleIDs()[jj]->getType()>=100 ) {
              int ijet=jts[jj-1]->ext<TJindex>(); 
              streamlog_out(DEBUG6)  << std::setw(18)<< ijet << " is from " << fafpi->getParticleIDs()[jj]->getType()/100  << std::endl;
            }
          }
        }
      }

      streamlog_out(DEBUG6) <<  std::endl;
      streamlog_out(DEBUG6) << "       Colour-neutral objects at the end of the parton shower " << std::endl;
      streamlog_out(DEBUG6) << "    number    px          py           pz          E           M       type       PDGs        jets " << std::endl;
      for ( unsigned kk=0 ; kk<fafpf_vec->size() ; kk++ ) {
        
        ReconstructedParticle* fafpf = dynamic_cast<ReconstructedParticle*>(fafpf_vec->at(kk));
        LCObjectVec jts;
        jts = FinalColourNeutral_Nav.getRelatedFromObjects(fafpf);

        const double* mom=fafpf->getMomentum();
	streamlog_out(DEBUG6) << std::setw(7) << kk+1 << std::setw(12) << mom[0] << std::setw(12) <<  mom[1] << std::setw(12)  << mom[2] << 
                             std::setw(12)  <<fafpf->getEnergy() << std::setw(12) << fafpf->getMass() << 
                             std::setw(7)  << fafpf->getParticleIDs()[0]->getType()  << "   ";
        if (  fafpf->getParticleIDs().size() < 3 ) streamlog_out(DEBUG6)  << std::setw(3) <<" ";
        for ( unsigned ipid=0 ; ipid < fafpf->getParticleIDs().size() ; ipid++ ) {
          if (ipid < fafpf->getParticleIDs().size()-1 ) {
            streamlog_out(DEBUG6)  << std::setw(3)<< fafpf->getParticleIDs()[ipid]->getPDG()  << "," ;
          } else {
            streamlog_out(DEBUG6)  << std::setw(3)<< fafpf->getParticleIDs()[ipid]->getPDG() << "   " ;
          }
        }
        if (  jts.size() < 2 ) streamlog_out(DEBUG6)  << std::setw(3) <<" ";
        for ( unsigned jj=0 ; jj<jts.size() ; jj++ ) {
           int ijet=jts[jj]->ext<TJindex>(); 
           if ( jj < jts.size()-1 ) {
             streamlog_out(DEBUG6) << std::setw(3) << ijet << "," ; 
           } else {
             streamlog_out(DEBUG6) << std::setw(3) << ijet  ;
           }
        }
        streamlog_out(DEBUG6) << std::endl;
      }
      streamlog_out(DEBUG6) << std::endl;

    }  // endif( mcpcol != NULL &&  rmclcol != NULL)


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

  int iii=0;
  int nMCP = mcpcol->getNumberOfElements()  ;
  //JL: added so that first_line is the same in all events!
  first_line=1;
  if ( nMCP > 4000 ) {
    streamlog_out(WARNING)  << " More than 4000 MCParticles in event " << evt->getEventNumber() << ",   run  " << evt->getRunNumber() << std::endl ;
  }
  for(int j=0; j< std::min(nMCP,4000) ; j++){

    MCParticle* mcp = dynamic_cast<MCParticle*>( mcpcol->getElementAt( j ) ) ;
    //  (assume this is not needed, ie. that the k-vector  can be set up correctly
    //    also for created-in-simulation  particles. The issue is if decay-products
    //    of a given particle will be consecutive).
    // if (mcp->getGeneratorStatus() == 1 || mcp->getGeneratorStatus() == 2 ) {

      iii++; 
      mcp_pyjets[iii]=mcp;
      mcp->ext<MCPpyjet>()=iii;
      ptemp=mcp->getMomentum();
      p[iii][1]=ptemp[0]; 
      p[iii][2]=ptemp[1]; 
      p[iii][3]=ptemp[2]; 
      p[iii][4]=mcp->getEnergy();
      p[iii][5]=mcp->getMass();
      k[iii][1]=mcp->getGeneratorStatus();
      if ( k[iii][1] == 2 ) {  k[iii][1]=11 ;}
      if ( k[iii][1] == 3 ) {  k[iii][1]=21 ;}
      if (mcp->isOverlay() ) {  k[iii][1]=k[iii][1]+30 ; }
      k[iii][2]=mcp->getPDG();

    //    }
  }
  nlund=iii;

  for (int j=1; j <= nlund ; j++ ) {
    k[j][3] = 0;
    k[j][4] = 0;
    k[j][5] = 0;
  }
  streamlog_out(DEBUG1) << " before Motherless check: first_line = " << first_line << std::endl ; 
  for(int j=0; j< nMCP ; j++){

    MCParticle* mcp = dynamic_cast<MCParticle*>( mcpcol->getElementAt( j ) ) ;
    int i;
    if ( mcp->ext<MCPpyjet>() != 0 ) {
      i= mcp->ext<MCPpyjet>();
      if (  mcp->getParents().size() != 0 ) {
        k[i][3]=mcp->getParents()[0]->ext<MCPpyjet>();
      } else {
//JL  was i > first_line+2, but due to the > sign, this skips the 3 first particles, while we just want to omit the first 2 here!        
        if ( i > first_line+1 &&  ! mcp->isOverlay() &&  k[i][3] == 0 && abs(mcp->getPDG()) != 6 && abs(mcp->getPDG()) != 24 ) {
          streamlog_out(MESSAGE4) << " Motherless generator particle " << i << ". pdg: " << mcp->getPDG() << std::endl ; 
          streamlog_out(MESSAGE4) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
//JL          k[i][3]=3;
          k[i][3]=-1;
        }
       
      }
      if (  mcp->getParents().size() == 2 ) {
        k[i][6]=mcp->getParents()[1]->ext<MCPpyjet>();
      } else {
        k[i][6]=0;
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

        if ( k[i][1]!=0 ) {
          int has_genstat0_daughter=0;
          for ( unsigned idat=0 ; idat <  mcp->getDaughters().size() ; idat++ ) {
            int jdat=mcp->getDaughters()[idat]->ext<MCPpyjet>();
            if ( mcp->getDaughters()[idat]->getGeneratorStatus()==0 || (mcp->getDaughters()[idat]->getSimulatorStatus()==0 && k[jdat][1]==1 )) {
              has_genstat0_daughter=1;
            }
          }
          if ( has_genstat0_daughter==1 ) {

            if ( k[i][1] < 30 ) {
              k[i][1]=1;
            }
            for ( unsigned idat=0 ; idat <  mcp->getDaughters().size() ; idat++ ) {
              int jdat=mcp->getDaughters()[idat]->ext<MCPpyjet>();
              if ( k[jdat][1] < 30 ) {
                k[jdat][1]=0;
              }
            }
          }
        } else {

          for ( unsigned idat=0 ; idat <  mcp->getDaughters().size() ; idat++ ) {
            int jdat=mcp->getDaughters()[idat]->ext<MCPpyjet>();
            if ( k[jdat][1] < 30 ) {
              k[jdat][1]=0;
            }
          }
        }

           // fill in k[][4] and [5]. Note that it is possible that one gets
           // k[][1]=1 but k[][4] and [5] != 0. This means that a generator
           // particle end in a simulation particle, which is OK for the
           // logic elsewhere !

        for ( unsigned idat=0 ; idat <  mcp->getDaughters().size() ; idat++ ) {
          int jdat=mcp->getDaughters()[idat]->ext<MCPpyjet>();
          if ( jdat < k[i][4] ||  k[i][4] == 0 ) {
            k[i][4]=jdat;
          }
          if ( jdat > k[i][5] ) {
            k[i][5]=jdat;
          }
        }
      } 
    }
  }

    // follows a VERY tedious part to fix wrongly assigned oarents to 94:s sometimes done in the stdhep reader....
    // (mostly (only?) a problem for 6f-leptonic)

  streamlog_out(DEBUG1) << " HEPEVT relation table before 94-fixing"  <<std::endl;
  streamlog_out(DEBUG1) << "       line      status    pdg       parent    first     last     second " << std::endl;
  streamlog_out(DEBUG1) << "                                             daughter daugther    parent " << std::endl;
  for ( int j=1 ; j<=nlund ; j++ ) {
    streamlog_out(DEBUG1) << std::setw(10) << j << std::setw(10) <<  k[j][1] << 
                 std::setw(10) <<  k[j][2] << std::setw(10) <<  k[j][3] << 
                 std::setw(10) <<  k[j][4] << std::setw(10) << k[j][5]  << std::setw(10) << k[j][6] <<std::endl;
  }
  for ( int i=1 ; i<=nlund ; i++ ) {
    if ( k[i][4] ==  k[i][5] && k[ k[i][4]][2] == 94 ) {
        // check if the parent-daughther relation for 94:s is really correct.
        // there is a fixup of missing relations in the input stdhep in stdhepreader
        // that sometimes goew astray (my bad, in fact). If this seems to be the case
        // remove the wrong one, and leave it to fix94 to get it right
      int line94=k[i][4];
      int parent1=k[line94][3] ; int parent2=k[line94][6] ;
      int kid1=k[line94][4] ; int kid2=k[line94][5] ;
      if ( parent1 != 0 && parent2 != 0 ) {
        if ( abs((p[kid1][4]+p[kid2][4])-(p[parent1][4]+p[parent2][4]))/(p[parent1][4]+p[parent2][4]) > 0.0001 ) {
          streamlog_out(WARNING) << " Inconsitent 94 object: Energy of parents = "<< p[kid1][4]+p[kid2][4] << " , of kids = " <<  p[parent1][4]+p[parent2][4] << std::endl;
          streamlog_out(WARNING) << " parents: " << parent1 << " " << parent2 << ". Kids: " << kid1 << " " << kid2  << std::endl;
          streamlog_out(WARNING) << " pdgs:    " << k[parent1][2] << " " << k[parent2][2]  << ". Kids: " << k[kid1][2]  << " " << k[kid2][2]   << std::endl;
          int evnb_print=0;
          if ( ! ((k[ parent1 ][2] == k[ kid1 ][2] &&
              k[ parent2 ][2] == k[ kid2 ][2] ) ||
             (k[ parent2 ][2] == k[ kid1 ][2] &&
              k[ parent1 ][2] == k[ kid2 ][2] ) ) ) {
               // wrong-wrong=wrong !!!
            int p1w=0;
            if ( k[ parent1 ][2] !=  k[ kid1 ][2] && k[ parent1 ][2] !=  k[ kid2 ][2] ) {
               // parent1 is wrong 
    	      streamlog_out(ERROR) << " parent 1 is wrong flavor-wise ?? " << std::endl;
              streamlog_out(ERROR) << " "<<  line94  << " " << parent1 << " " <<  parent2  << " " <<kid1 << " " <<  kid2  << std::endl;
              streamlog_out(ERROR)  << " "<<  line94  << " " << k[parent1 ][2] << " " <<  k[parent2 ][2]  << " " <<
                                            k[kid1 ][2] << " " <<  k[kid2 ][2]  << std::endl;
              streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
              streamlog_out(ERROR) << std::endl ;  if( streamlog::out.write<streamlog::ERROR>() ){ evnb_print=1;}
              k[line94][3] = 0;
              k[parent1][4]= k[parent1][5]=0;
              p1w=1;
            }
            int p2w=0 ;
            if ( k[ parent2 ][2] !=  k[ kid1 ][2] && k[ parent2 ][2]!=  k[ kid2 ][2] ) {
             //     // parent2 is wrong 
              p2w=1;
              if ( p1w == 1 ) {
    	        streamlog_out(ERROR) << " BOTH parents wrong  flavor-wise ? ! " << std::endl;
                streamlog_out(ERROR) << " " <<  line94  << " " << parent1 << " " <<  parent2  << " " <<kid1 << " " <<  kid2  << std::endl;
                streamlog_out(ERROR) << " "  <<  line94  << " " << k[parent1 ][2] << " " <<  k[parent2 ][2]  << " " <<
                                            k[kid1 ][2] << " " <<  k[kid2 ][2]  << std::endl;
                streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                streamlog_out(ERROR) << std::endl ; if( streamlog::out.write<streamlog::ERROR>() ){ evnb_print=1;}
              } else {
                streamlog_out(DEBUG3) << " parent 2 is wrong flavor-wise  " << std::endl;
                streamlog_out(DEBUG3) << " "  <<  line94  << " " << parent1 << " " <<  parent2  << " " <<kid1 << " " <<  kid2  << std::endl;
                streamlog_out(DEBUG3) << " "  <<  line94  << " " << k[parent1 ][2] << " " <<  k[parent2 ][2]  << " " <<
                                            k[kid1 ][2] << " " <<  k[kid2 ][2]  << std::endl;
                streamlog_out(DEBUG3) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                streamlog_out(DEBUG3) << std::endl ; if( streamlog::out.write<streamlog::DEBUG3>() ){ evnb_print=1;}
  
              }
              k[line94][6] = 0;
              k[parent2][4]= k[parent2][5]=0;
            } 
            if ( p1w ==0 && p2w == 0 ) {
    	      streamlog_out(DEBUG3) << "  neither parent wrong  flavor-wise ?? " << std::endl;
              streamlog_out(DEBUG3) << " "<<  line94  << " " << parent1 << " " <<  parent2  << " " <<kid1 << " " <<  kid2  << std::endl;
              streamlog_out(DEBUG3) << " " <<  line94  << " " << k[parent1 ][2] << " " <<  k[parent2 ][2]  << " " <<
                                            k[kid1 ][2] << " " <<  k[kid2 ][2]  << std::endl;
              streamlog_out(DEBUG3) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
              streamlog_out(DEBUG3) << std::endl ; if( streamlog::out.write<streamlog::DEBUG3>() ){ evnb_print=1;} 
              k[line94][6] = 0;
              k[parent2][4]= k[parent2][5]=0;
            }
          } else {
  
            // flavours OK - check direction
  
            double dirdiff1[4];
            double dirdiff2[4];
            double dirp1[4], dirk1[4], abspp1, abspk1;
            double dirp2[4], dirk2[4], abspp2, abspk2;
            abspp1= sqrt(p[ parent1][1]*p[ parent1][1]+p[ parent1][2]*p[ parent1][2]+ p[ parent1][3]*p[ parent1][3]);
            abspk1= sqrt(p[ kid1][1]*p[ kid1][1]+p[ kid1][2]*p[ kid1][2]+ p[ kid1][3]*p[ kid1][3]);
            abspp2= sqrt(p[ parent2][1]*p[ parent2][1]+p[ parent2][2]*p[ parent2][2]+ p[ parent2][3]*p[ parent2][3]);
            abspk2= sqrt(p[ kid2][1]*p[ kid2][1]+p[ kid2][2]*p[ kid2][2]+ p[ kid2][3]*p[ kid2][3]);
            dirp1[1]=p[ parent1][1]/abspp1; dirp1[2]=p[ parent1][2]/abspp1; dirp1[3]=p[ parent1][3]/abspp1;
            dirk1[1]=p[ kid1][1]/abspk1; dirk1[2]=p[ kid1][2]/abspk1; dirk1[3]=p[ kid1][3]/abspk1;
            dirp2[1]=p[ parent2][1]/abspp2; dirp2[2]=p[ parent2][2]/abspp2; dirp2[3]=p[ parent2][3]/abspp2;
            dirk2[1]=p[ kid2][1]/abspk2; dirk2[2]=p[ kid2][2]/abspk2; dirk2[3]=p[ kid2][3]/abspk2;
  
            if (  k[ parent1 ][2] ==  k[ kid1 ][2] &&  k[ parent2 ][2] ==  k[ kid2 ][2] ) {
              dirdiff1[1]=dirp1[1]-dirk1[1]; dirdiff1[2]=dirp1[2]-dirk1[2]; dirdiff1[3]=dirp1[3]-dirk1[3]; 
              dirdiff2[1]=dirp2[1]-dirk2[1]; dirdiff2[2]=dirp2[2]-dirk2[2]; dirdiff2[3]=dirp2[3]-dirk2[3]; 
              
              if ( sqrt( dirdiff1[1]* dirdiff1[1] + dirdiff1[2]* dirdiff1[2] + dirdiff1[3]* dirdiff1[3] )  > 0.001 ) {
                  // parent1 is wrong 
    	        streamlog_out(ERROR) << " parent 1 is wrong direction-wise ?? " << std::endl;
                streamlog_out(ERROR) << " "<<  line94  << " " << parent1 << " " <<  parent2  << " " <<kid1 << " " <<  kid2  << std::endl;
                streamlog_out(ERROR)  << " " << k[parent1 ][2] << " " <<  p[parent1 ][1]  << " " <<  p[parent1 ][2]  << " " <<  p[parent1 ][3]  << std::endl;
                streamlog_out(ERROR)  << " " << k[kid1 ][2] << " " <<  p[kid1 ][1]  << " " <<  p[kid1 ][2]  << " " <<  p[kid1 ][3]  << std::endl;
                streamlog_out(ERROR)  << " " << k[kid1 ][2] << " " <<  dirdiff1[1]  << " " <<  dirdiff1[2]  << " " <<  dirdiff1[3]  << std::endl;
                streamlog_out(ERROR)  << " " <<  sqrt( dirdiff1[1]* dirdiff1[1] + dirdiff1[2]* dirdiff1[2] + dirdiff1[3]* dirdiff1[3] ) <<  std::endl;
                streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                streamlog_out(ERROR) << std::endl ; if( streamlog::out.write<streamlog::ERROR>() ){ evnb_print=1;}
                k[line94][3] = 0;
                k[parent1][4]= k[parent1][5]=0;
  
              } else if ( sqrt( dirdiff2[1]* dirdiff2[1] + dirdiff2[2]* dirdiff2[2] + dirdiff2[3]* dirdiff2[3] )  > 0.001 ) {
    	        streamlog_out(DEBUG3) << " parent 2 is wrong direction-wise. " << std::endl;
                streamlog_out(DEBUG3) << " "<<  line94  << " " << parent1 << " " <<  parent2  << " " <<kid1 << " " <<  kid2  << std::endl;
                streamlog_out(DEBUG3)  << " " << k[parent2 ][2] << " " <<  p[parent2 ][1]  << " " <<  p[parent2 ][2]  << " " <<  p[parent2 ][3]  << std::endl;
                streamlog_out(DEBUG3)  << " " << k[kid2 ][2] << " " <<  p[kid2 ][1]  << " " <<  p[kid2 ][2]  << " " <<  p[kid2 ][3]  << std::endl;
                streamlog_out(DEBUG3)  << " " << k[kid2 ][2] << " " <<  dirdiff2[1]  << " " <<  dirdiff2[2]  << " " <<  dirdiff2[3]  << std::endl;
                streamlog_out(DEBUG3)  << " " <<  sqrt( dirdiff2[1]* dirdiff2[1] + dirdiff2[2]* dirdiff2[2] + dirdiff2[3]* dirdiff2[3] ) <<  std::endl;
                streamlog_out(DEBUG3) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                streamlog_out(DEBUG3) << std::endl ; if( streamlog::out.write<streamlog::DEBUG3>() ){ evnb_print=1;}
                k[line94][6] = 0;
                k[parent2][4]= k[parent2][5]=0;
              } else {
                streamlog_out(ERROR) << " However, I can't find what exactly was wrong. " << std::endl ; 
                streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                streamlog_out(ERROR) << std::endl ; evnb_print=1;
                // all OK
              } 
            } else if (  k[ parent1 ][2] ==  k[ kid2 ][2] &&  k[ parent2 ][2] ==  k[ kid1 ][2] ) {
              dirdiff1[1]=dirp1[1]-dirk2[1]; dirdiff1[2]=dirp1[2]-dirk2[2]; dirdiff1[3]=dirp1[3]-dirk2[3]; 
              dirdiff2[1]=dirp2[1]-dirk1[1]; dirdiff2[2]=dirp2[2]-dirk1[2]; dirdiff2[3]=dirp2[3]-dirk1[3]; 
              if (  sqrt( dirdiff1[1]* dirdiff1[1] + dirdiff1[2]* dirdiff1[2] + dirdiff1[3]* dirdiff1[3] )  > 0.001 ) {
                  // parent1 is wrong 
                streamlog_out(ERROR) << " parent 1 is wrong direction-wise ?? " << std::endl;
                streamlog_out(ERROR) << " "<<  line94  << " " << parent1 << " " <<  parent2  << " " <<kid1 << " " <<  kid2  << std::endl;
                streamlog_out(ERROR)  << " " << k[parent1 ][2] << " " <<  p[parent1 ][1]  << " " <<  p[parent1 ][2]  << " " <<  p[parent1 ][3]  << std::endl;
                streamlog_out(ERROR)  << " " << k[kid2 ][2] << " " <<  p[kid2 ][1]  << " " <<  p[kid2 ][2]  << " " <<  p[kid2 ][3]  << std::endl;
                streamlog_out(ERROR)  << " " << k[kid2 ][2] << " " <<  dirdiff1[1]  << " " <<  dirdiff1[2]  << " " <<  dirdiff1[3]  << std::endl;
                streamlog_out(ERROR)  << " " <<  sqrt( dirdiff1[1]* dirdiff1[1] + dirdiff1[2]* dirdiff1[2] + dirdiff1[3]* dirdiff1[3] ) <<  std::endl;
                streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                streamlog_out(ERROR) << std::endl ; if( streamlog::out.write<streamlog::ERROR>() ){ evnb_print=1;}
  
                k[line94][3] = 0;
                k[parent1][4]= k[parent1][5]=0;
  
              } else if (  sqrt( dirdiff2[1]* dirdiff2[1] + dirdiff2[2]* dirdiff2[2] + dirdiff2[3]* dirdiff2[3] )  > 0.001 ) {
                streamlog_out(DEBUG3) << " parent 2 is wrong direction-wise. " << std::endl;
                streamlog_out(DEBUG3) << " "<<  line94  << " " << parent1 << " " <<  parent2  << " " <<kid1 << " " <<  kid2  << std::endl;
                streamlog_out(DEBUG3)  << " " << k[parent2 ][2] << " " <<  p[parent2 ][1]  << " " <<  p[parent2 ][2]  << " " <<  p[parent2 ][3]  << std::endl;
                streamlog_out(DEBUG3)  << " " << k[kid1 ][2] << " " <<  p[kid1 ][1]  << " " <<  p[kid1 ][2]  << " " <<  p[kid1 ][3]  << std::endl;
                streamlog_out(DEBUG3)  << " " << k[kid1 ][2] << " " <<  dirdiff2[1]  << " " <<  dirdiff2[2]  << " " <<  dirdiff2[3]  << std::endl;
                streamlog_out(DEBUG3)  << " " <<  sqrt( dirdiff2[1]* dirdiff2[1] + dirdiff2[2]* dirdiff2[2] + dirdiff2[3]* dirdiff2[3] ) <<  std::endl;
                streamlog_out(DEBUG3) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                streamlog_out(DEBUG3) << std::endl ; if( streamlog::out.write<streamlog::DEBUG3>() ){ evnb_print=1;}
  
                k[line94][6] = 0;
                k[parent2][4]= k[parent2][5]=0;
              } else {
                streamlog_out(ERROR) << " However, I can't find what exactly was wrong. " << std::endl ; 
                streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
                streamlog_out(ERROR) << std::endl ; if( streamlog::out.write<streamlog::ERROR>() ){ evnb_print=1;}
                // all OK
              } 
            } else {
              streamlog_out(ERROR) << " How could I arrive here ??? "  << std::endl;
              streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
              streamlog_out(ERROR) << std::endl ; if( streamlog::out.write<streamlog::ERROR>() ){ evnb_print=1;}
            }
          }
          if( evnb_print==0 ) {
            streamlog_out(WARNING) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
            streamlog_out(WARNING) << std::endl ; 
          }
        }
      }
    }
  }
    // end of a VERY tedious part to fix wrongly assigned oarents to 94:s sometimes done in the stdhep reader....

  streamlog_out(DEBUG2) << " HEPEVT relation table"  <<std::endl;
  streamlog_out(DEBUG2) << "       line      status    pdg       parent    first     last     second " << std::endl;
  streamlog_out(DEBUG2) << "                                             daughter daugther    parent " << std::endl;
  for ( int j=1 ; j<=nlund ; j++ ) {
    streamlog_out(DEBUG2) << std::setw(10) << j << std::setw(10) <<  k[j][1] << 
                 std::setw(10) <<  k[j][2] << std::setw(10) <<  k[j][3] << 
                 std::setw(10) <<  k[j][4] << std::setw(10) << k[j][5]  << std::setw(10) << k[j][6] <<std::endl;
  }
} 
void TrueJet::true_lepton()
{
  int ihard_lepton_1[4011],CMshowers[4011], nCMshowers;
  int ihard_lepton[4011];
  int lept, gotit ;
      //! count number of hard leptons. The condition is: it should be a lepton lepton (pdg 11 to 16). It is
      //! hard either if it is  comming from line 3, or if comes from line 1, goes to a single daughter which
      //! is an electron. (The second part is needed for beam-remenant electrons in odd-f events)

   n_hard_lepton = 0;
   for ( int kk=1 ; kk<= nlund ; kk++ ) {
     if ( ( abs(k[kk][2]) >= 11 && abs(k[kk][2]) <= 16 ) )  {
       streamlog_out(DEBUG1) << " found lepton, kk = " << kk << ", k[kk][3] = " << k[kk][3] << ", k[3][3] = " << k[3][3] <<  std::endl ; 
       if (
//JL          ((k[kk][3] == 3) || (k[kk][3] == 1 && k[kk][4] ==k[kk][5] && abs(k[k[kk][5]][2])==11 ))) {
// make sure that line 3 has a parent, thus it is not a parentless particle in Higgs sample
//JL          ((k[kk][3] == 3 && k[3][3] > 0) || (k[kk][3] == 0 && kk > 2) || (k[kk][3] == 1 && k[kk][4] ==k[kk][5] && abs(k[k[kk][5]][2])==11 )))) {
          ((k[kk][3] == 3 && k[3][3] > 0) || (k[kk][3] == -1) || (k[kk][3] == 1 && k[kk][4] ==k[kk][5] && abs(k[k[kk][5]][2])==11 ))) {
         n_hard_lepton++;
         ihard_lepton_1[n_hard_lepton]=kk;
         ihard_lepton[n_hard_lepton]=0;
         streamlog_out(DEBUG1) << "     hard according to condition 1 " << std::endl;
       }
       //JL also count leptons as hard if their mother or grandma is a Higgs boson
       //JL necessary to protect into W bosons from semi-leptonic quark decays etc
       else if (abs(k[k[kk][3]][2]) == 25 || ((abs(k[k[kk][3]][2]) == 24 || abs(k[k[kk][3]][2]) == 23) && abs(k[k[k[kk][3]][3]][2]) == 25)) {
         n_hard_lepton++;
         ihard_lepton_1[n_hard_lepton]=kk;
         ihard_lepton[n_hard_lepton]=0;
         streamlog_out(DEBUG1) << "     hard according to condition 2 " << std::endl;
       }
     }
   }
   if ( n_hard_lepton == 0 ) return ;

      //! Sort the leptons in "color singlets", ie. in groups of flavour/anti-flavour. In most events this is
      //! the order they come in, but in the case of all flavour being the same, it isn't

   nCMshowers=0 ;
   for ( int kk=1 ; kk<=nlund ; kk++ ) {
     if ( k[kk][2] == 94 ) {
       nCMshowers++;
       CMshowers[nCMshowers]=kk;
     }
   }

   int n_94_related=0 ;
   for (int kk=1; kk <= nCMshowers ; kk++ ) {
     int line=k[CMshowers[kk]][3];
       // DO WHILE ( (.NOT. ANY (ihard_lepton_1(1:n_hard_lepton) == line )) .AND. line > 3 ) ; line=k(line,3) ; ENDDO
     gotit=0;
//JL     while ( line > 3 ) {
//JL  include line 3 if parent is -1 (means has been set above because mother is missing)
     while ( line > 3 || k[line][3] == -1) {
       for ( int jj=1 ; jj <= n_hard_lepton ; jj++ ) {
         if (  ihard_lepton_1[jj] == line ) { 
           gotit=1 ; 
         }
       }
       if ( gotit==1 ) {
         break;
       } else {
         line=k[line][3];
       }
     }   
//JL      if ( line <= 3 ) { continue ; }
//JL  exclude line 3 if parent is -1 (means has been set above because mother is missing)
     if ( line <= 3 && k[line][3] != -1) { continue ; }
       // ll = transfer(maxloc( [(1,kk=1,n_hard_lepton)] ,(ihard_lepton_1(1:n_hard_lepton) == line)),ll) 
     int ll=0;
     for ( int jj=1 ; jj <= n_hard_lepton ; jj++ ) {
       if ( ihard_lepton_1[jj] == line ) {
         ll=jj;
         break;
       }
     }
     n_94_related=n_94_related+1;
     ihard_lepton[n_94_related] =  ihard_lepton_1[ll]  ; ihard_lepton_1[ll]=0 ;
     line=k[CMshowers[kk]][6];
     gotit=0;
//JL     while ( line > 3 ) {
//JL  include line 3 if parent is -1 (means has been set above because mother is missing)
     while ( line > 3 || k[line][3] == -1) {
       for ( int jj=1 ; jj <= n_hard_lepton ; jj++ ) {
         if (  ihard_lepton_1[jj] == line ) { 
           gotit=1 ; 
         }
       }
       if ( gotit==1 ) {
         break;
       } else {
         line=k[line][3];
       }
     }   
//JL     if ( line <= 3 ) { std::cout << " oups ! " << std::endl; }
//JL  exclude line 3 if parent is -1 (means has been set above because mother is missing)
     if ( line <= 3 && k[line][3] != -1) {  std::cout << " oups ! " << std::endl; }
     for ( int jj=1 ; jj <= n_hard_lepton ; jj++ ) {
       if ( ihard_lepton_1[jj] == line ) {
         ll=jj;
         break;
       }
     }
     n_94_related=n_94_related+1;
     ihard_lepton[n_94_related] =  ihard_lepton_1[ll]  ; ihard_lepton_1[ll]=0 ;
   }


   for (int kk= n_94_related+1; kk<=n_hard_lepton; kk++ ) {
     if ( ihard_lepton[kk] == 0 ) {
       for ( int ll=1; ll<=n_hard_lepton ; ll++ ) {
         if ( ihard_lepton_1[ll] != 0 ) {
           ihard_lepton[kk]= ihard_lepton_1[ll];
           ihard_lepton_1[ll]= 0;
           break ;
         }
       }
       //flav= (k[ihard_lepton[kk]][2])>0 ? abs(k[ihard_lepton[kk]][2])-11 : -(abs(k[ihard_lepton[kk]][2])-11) ; flav=flav/2;
       for ( int ll=1; ll <= n_hard_lepton ; ll++ ) {
         if (  flavour(k[ihard_lepton_1[ll]][2]) == -flavour(k[ihard_lepton[kk]][2]) ) {
           ihard_lepton[kk+1]=ihard_lepton_1[ll] ;
           ihard_lepton_1[ll]=0 ;
           break; 
         }
       }
     }
   }
      //! assign to jets. Basically the jet is the index in the list, -ve if it later on goes into a
      //! 94, +ve otherwise.

   for ( int  kk=1; kk<=n_hard_lepton ; kk++ ) {

     lept = ihard_lepton[kk];
     streamlog_out(DEBUG1) << " assigning leptons to jets: ihard_lepton[" << kk << "] = " << ihard_lepton[kk] << " and has PDG " << k[lept][2] << std::endl ; 

        // ! initially, set to -kk

     jet[lept]=-(i_jet+kk);
     //JL this is premature - set elementon below after checking for 94 daughters of leptons (tau's!)
     //MB reverted that, the elementon should be as close as possible to the hard interaction; the check for 94:s is
     //just to check if the jet-number should be positive or negative.
     if ( k[lept][4] == k[lept][5] &&  k[lept][4] != 0 ) {
       if ( k[ k[lept][4] ][2] == 94 ) {
          elementon[i_jet+kk]=lept ;
       } else {
         elementon[i_jet+kk]=k[lept][4];  // Behold the beauty of my new word !
       }
     } else {
       elementon[i_jet+kk]=lept ;
     } 
     fafp_last[i_jet+kk]=lept;  // so it goes when too specific variable names were choosen ...

        // ! loop until either stable or 94 descendent found

     while ( k[lept][2] != 94 && k[lept][1] != 1 && k[lept][4] == k[lept][5] ) {
       lept = k[lept][4] ;
     }

     if (  k[lept][2] == 94 ) {

         // ! 94 decendant found. find which of the two daughters of the 94
         // ! corresponds to the lepton we're considering in this iteration, and move
         // ! to that line in pyjets

       if (  k[k[lept][4]][2] == k[ihard_lepton[kk]][2] ) {
         lept = k[lept][4];
       } else if (  k[k[lept][5]][2] == k[ihard_lepton[kk]][2]  ) {
         lept=k[lept][5];
       } else {
         streamlog_out(WARNING) << " WARNING: Error in hard lepton-finding " << std:: endl;
         streamlog_out(WARNING) << lept << "  goes to " << k[lept][4] << " and "  << k[lept][5] << " which are " <<  k[k[lept][4]][2] << " and " <<  k[k[lept][5]][2]<< std:: endl;
         streamlog_out(WARNING) << " ihard_lepton on " << ihard_lepton[kk] << " which is " <<  k[ihard_lepton[kk]][2] << std::endl ;
         streamlog_out(WARNING) << " This means that the wrong guess for grouping of leptons to probable boson parent were made, due to an ambigous situation, " << std::endl ;
         streamlog_out(WARNING) << " ie. jets are OK, but di-jets are not in the expected order. " << std::endl ;
         streamlog_out(WARNING) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
         streamlog_out(WARNING) << std::endl ; 
       }
     }

        // ! Here we know for sure that line lept in pyjets is a post-first-94-if-any lepton:
        // ! set jet to +ve.
     jet[lept]=i_jet+kk ;
     //JL set elementon here after 94 check!
     //MB: reverted
   }

     // ! set jet for all further descendants of hard leptons . Don't change already set assignments,
     
   for ( int kk=1; kk<=nlund ; kk++ ) {
     if ( jet[kk] == 0 && abs(jet[k[kk][3]]) > i_jet && abs(jet[k[kk][3]]) <= i_jet+n_hard_lepton ) {

         // ! this is indeed a descendant of a hard lepton 

       jet[kk] = jet[k[kk][3]];
     }
   }
   i_jet=i_jet+n_hard_lepton ;


}
void TrueJet::cluster()
{
   int clus[4011];
   int clu , jet1, jet2;

   nclu=0;
   for (int kk=1; kk<=nlund ; kk++ ) {
     if ( k[kk][2] == 91 && k[kk][1] < 30  ) {
       nclu++;
       clus[nclu]=kk;
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

     jet1=i_jet+2*(i_clu-1)+1 ; jet2= i_jet+2*(i_clu-1)+2;
     elementon[jet1]=k[clu][3];
     elementon[jet2]=elementon[jet1];
     while ( k[elementon[jet2]+1][4] == clu ) {
       elementon[jet2]=elementon[jet2]+1;
     }

     assign_jet(jet1,jet2,clu) ;

   }
   i_jet=i_jet+ 2*nclu;

}
void TrueJet::isr()
{
  int iisr[4011];
   //   nisr = COUNT([ ( (ANY(k(kk,3)==[1,3]) .AND.  (k(kk,2)==22) &
   //          .AND.(k(kk,4)==k(kk,5)) .AND. (k(k(kk,4),2) ==22)), kk=1,nlund)  ])
   // index=[(kk,kk=1,nlund+10)]
   // iisr = PACK(index(1:nlund),[ ( (ANY(k(kk,3)==[1,3]) .AND.  (k(kk,2)==22) &
   //                     .AND.(k(kk,4)==k(kk,5)) .AND. (k(k(kk,4),2) ==22) ), kk=1,nlund)  ] )
   nisr=0;
   for (int kk=1; kk<=nlund ; kk++ ) {
     if ( (k[kk][3] == 1 || k[kk][3] == 3 || k[kk][3] == 0 ) &&  ( k[kk][2]==22 ) && (k[kk][4]==k[kk][5]) && (k[k[kk][4]][2] ==22) ) {
       nisr++;
       iisr[nisr] = kk ;
     }
   }
   if ( nisr > 0 ) {
     // ! jet assignment - index in the list. -ve if "decayimg", which might be the case
     // ! for explicitly requested gammas.
     for ( int kk=1 ; kk<=nisr ; kk++ ) {
       elementon[i_jet+kk] =  iisr[kk];
       if ( k[ iisr[kk]][1] == 1 ) {
         jet[ iisr[kk]]=i_jet+kk;
       } else {
         jet[ iisr[kk]]=-(i_jet+kk);
       }
     }
       
     for ( int  kk=1; kk<=nlund ; kk++ ) {
       if ( jet[kk] == 0 && abs(jet[k[kk][3]]) > i_jet && abs(jet[k[kk][3]]) <= i_jet+nisr ) {
  
          //  ! this is indeed a descendant of an isr
  
         jet[kk] = abs(jet[k[kk][3]]);
       }
     }
     i_jet=i_jet+nisr;
   }
}
void TrueJet::string()
{
  int n_left, sstr,lquark,lstr, istr,str_nb,jet1,jet2,istr_previous,iback;
  int str_index[4011], rev_index[4011];
     //!! Move jets away to make room for string-induced ones at the lowest indicies
       // 
       //    tjet=0
       //    WHERE ( jet /=0 ) tjet=sign(abs(jet)+2*nstr,jet) ; jet=tjet
       //    elementon= EOSHIFT(elementon(1:njet),-2*nstr)
       //    nfsr = EOSHIFT(nfsr(1:njet),-2*nstr)
       //    hard = EOSHIFT(hard(1:njet),-2*nstr)
       //    fafp_last= EOSHIFT(fafp_last(1:njet),-2*nstr)
       //    boson = EOSHIFT(boson(1:njet),-2*nstr)
       //    fafp_boson = EOSHIFT(fafp_boson(1:njet),-2*nstr)

   for ( int kk=1 ; kk<=nlund ; kk++ ) {
     rev_index[kk]=0;
     if ( jet[kk] != 0 ) {
       jet[kk]=(jet[kk]>0 ? jet[kk]+2*nstr : jet[kk]-2*nstr);
     }
   }
   if ( nstr > 0 ) {
     for (int kk=njet ; kk>=1 ; kk-- ) {
       elementon[kk+2*nstr]=elementon[kk] ; 
       elementon[kk]=0;
       nfsr[kk+2*nstr]=nfsr[kk] ; 
       nfsr[kk]=0;
       fafp_last[kk+2*nstr]=fafp_last[kk];
       fafp_last[kk]=0;
       boson[kk+2*nstr]=boson[kk];
       boson[kk]=0;
       fafp_boson[kk+2*nstr]=fafp_boson[kk];
       fafp_boson[kk]=0;
     }
   }

   i_jet=0 ; 

        //! Decode the event record. First find the line of the first string,
        //! the last quark of the last string (lquark, at the line before the first string),
        //! the last string (lstr, the daughter of the lquark quark), and the number of strings (nstr).
  
      //    str_index=PACK(=[(kk,kk=1,pjetss)],(jet(1:nlund)==0 .and. k(1:nlund,2) /= 91.and. k(1:nlund,2) /= 21))
      //   
      //    n_left=COUNT ( (jet(1:nlund)==0 .and. k(1:nlund,2) /= 91.and. k(1:nlund,2) /= 21) )
      //    rev_index(1:nlund)=0 ; FORALL ( kk=1:n_left) rev_index(str_index(kk))=kk
      //   
      //    sstr = transfer(maxloc( [(1,kk=1,n_left)] ,( (92 == [ (k(str_index(kk),2),kk=1,n_left ) ] ) ) ),sstr) ; 
      //    lquark=str_index(sstr-1) ;    lstr=K(lquark,4) ;  istr=lstr ; str_nb = nstr


   n_left=0; sstr=0;
   for ( int kk=1 ; kk<= nlund ; kk++ ) {
     streamlog_out(DEBUG0) << " particle " << kk << " with PDG " << k[kk][2] << " is assigned to jet " << jet[kk] << std::endl;
     if (jet[kk]==0 && k[kk][2] != 91 &&  k[kk][2] != 21 && k[kk][1] < 30 ){
       n_left++;
       str_index[n_left]=kk;
       rev_index[kk]=n_left;
       if ( sstr == 0 &&  k[kk][2] == 92 ) {
         sstr=n_left;
       }
     }
   }     
   lquark=str_index[sstr-1];
   if (k[k[lquark][4]][2] != 92 ) {
     streamlog_out(ERROR) << " INFO: Non-contigous string ancestors (1)" << std::endl;
     streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ;
     streamlog_out(ERROR) << std::endl ; 
  }
   lstr=k[lquark][4] ;  istr=lstr ; str_nb = nstr ;
   streamlog_out(DEBUG3) << "  sstr, lstr, lquark, i_jet, str_nb " << sstr << " " << lstr << " " << lquark  << " " << i_jet  << " " << str_nb << std::endl;

        //! now the jet-assignment. Start by looking at the last string (lstr),
        //! stop once the first has been reached. decide if the particle should be
        //! assigned to the first or last quark end of the string, by checking
        //! the sign of the projection of the hadron momentum on the difference
        //! between the two quark directions.
  
   while( istr >= str_index[sstr]) {
  
     streamlog_out(DEBUG0) << "calculating jet1/2 from i_jet = " << i_jet << ", str_nb = " << str_nb  << std::endl;
     streamlog_out(DEBUG0) << "istr =  = " << istr << ", k[istr][3] = " << k[istr][3]  << ", lquark = " << lquark << std::endl;
     jet1= i_jet+2*(str_nb-1)+1 ; jet2=i_jet+2*(str_nb-1)+2 ; elementon[jet1]=k[istr][3] ; elementon[jet2]=lquark ;

     streamlog_out(DEBUG0) << "calling assign_jet with jet1, jet2 = " << jet1 << ", " << jet2 << std::endl;
     assign_jet(jet1,jet2,istr) ;


        //! the previous string is the mother of its hadrons, eg. the last one, which is
        //! on the line before the current string:
     
        //! new string:
  
     istr_previous=istr ;
     iback=str_index[rev_index[istr]-1]  ; istr = k[iback][3] ;
  
     if ( istr < str_index[sstr] ) break ;
  
        //! back-up: end-quark of the previous string is the line before the start-quark
        //! of the present string:
  
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
       if ( lquark == 0 ) { 
         streamlog_out(ERROR) << " ERROR: Quack ?! " << std::endl ;         
         streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ;
         streamlog_out(ERROR) << std::endl ; 

       }
     }
 
        //! decrement the string-counter
  
     str_nb--;
  
   }
}
void TrueJet::assign_jet(int jet1,int jet2,int this_fafp)
{
 double dir_diff[4];
 int first_p1, first_p2;
  streamlog_out(DEBUG0) << "assign_jet: elementon[" << jet1 << "] = " << elementon[jet1] 
                                  << ", elementon[" << jet2 << "] = " << elementon[jet2] << std::endl;
  
  
  first_parton( elementon[jet1], jet1 , first_p1, fafp_last[jet1], nfsr[jet1], boson[jet1], fafp_boson[jet1]);
  first_parton( elementon[jet2], jet2 , first_p2, fafp_last[jet2], nfsr[jet2], boson[jet2], fafp_boson[jet2]);

     //   !  difference between the two quark directions: one is at elementon(jet1), the
     //   !  other is qyuark(jet2) = three lines of fortran, 20 lines of C++ ...

   //  dir_diff(1:3) = p(elementon(jet2),1:3)/sqrt(sum( p(elementon(jet2),1:3)**2))- p(elementon(jet1),1:3)/ sqrt(sum(p(elementon(jet1),1:3)**2))
   //   jet(k(fafp,4):k(fafp,5)) = jet2  
   //   FORALL ( kk=k(fafp,4):k(fafp,5), (DOT_PRODUCT(dir_diff(1:3),p(kk,1:3)) <=0))  jet(kk)= jet1

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

   for (int kk=k[this_fafp][4]; kk<=k[this_fafp][5] ; kk++ ) {
     double dot=0;
     for ( int jj=1 ; jj<=3 ; jj++ ) {
       dot+=dir_diff[jj]*p[kk][jj] ;
     }
     if ( dot <= 0. ) {
       jet[kk] = jet1 ;
     } else {
       jet[kk] = jet2 ;
     } 
   } 
   for ( int kk=k[elementon[jet2]][5]+1; kk <=nlund ; kk++ ) {
     if ( jet[kk] == 0 && k[kk][1] < 30 ) {
       jet[kk] = abs(jet[k[kk][3]]);
       if (  k[kk][2] == 21 || abs( k[kk][2])<=6 ) {
               //  can happen if a paricle decay is done by gluons (Ypsilon etc.)
         jet[kk] = -jet[kk];
       }
     }
   }

}
void TrueJet::first_parton(int this_partic,int this_jet,int& first_partic,int& last_94_parent,int& nfsr_here,int& info,int& info2)

//   INTEGER, INTENT(IN)  :: this_partic
//   INTEGER, INTENT(IN)  :: this_jet
//   INTEGER, INTENT(OUT) :: first_partic
//   INTEGER, INTENT(OUT) :: last_94_parent
//   INTEGER, INTENT(INOUT) :: nfsr
//   INTEGER, INTENT(OUT) :: info
//   INTEGER, INTENT(OUT) :: info2
{
   int this_quark,this_quark_save,quark_pdg,boson_ancestor,istat,boson_last_94_parent,istat2;

   this_quark=this_partic;
   jet[this_quark]=-this_jet;
   quark_pdg=k[this_quark][2];
   
   streamlog_out(DEBUG0) << "first_parton: this_quark = " << this_quark 
                         << ",  quark_pdg = " <<  quark_pdg << std::endl;
                     

   last_94_parent =  0;

      //! back-track quarks to the begining of the event record, or until the quark-chain is broken
      //! Beginning of record is hit when the parent of a quark is an electron, and the chain is
      //! brocken if the parent of a quark is a gluon (21) or a W (24).

      //! start at this_quark, which might on entry actually be something else (a gluon or a W)

//JL   while (  abs(k[k[this_quark][3]][2]) != 11 &&  abs(k[k[this_quark][3]][2]) != 21   &&  abs(k[k[this_quark][3]][2]) != 24  ) {
// JL and what about Z?
   while (  k[this_quark][3] > 0 &&
            abs(k[k[this_quark][3]][2]) != 11 &&  abs(k[k[this_quark][3]][2]) != 21   
        &&  abs(k[k[this_quark][3]][2]) != 23 &&  abs(k[k[this_quark][3]][2]) != 24    &&  abs(k[k[this_quark][3]][2]) != 25  ) {
 
     streamlog_out(DEBUG0) << "start of while loop:  this_quark = " << this_quark 
                           << ", its PDG = " << k[this_quark][2] 
                           << ", its parent = " << k[this_quark][3] 
                           << ", parent's PDG  = " << abs(k[k[this_quark][3]][2])  
                           << ", quark_pdg  = " << quark_pdg << std::endl;
     this_quark_save =   this_quark ;

     //! step back in the history

     this_quark=k[this_quark][3];
     if ( quark_pdg == 21 || abs(quark_pdg) == 24 ) {
          //! looking at a gluon or a W. We've alredy steped to the
          //! parent of the boson (which should be a quark), so now we might 
          //! be loking at a different flavour:
          //! 

       quark_pdg=k[this_quark][2];

       if ( quark_pdg == 94 ) {


	    //! well, actually it wasn't a quark, but a 94. Need to step
	    //! back one step further and find the right parent to go on following

//JL         if ( abs(k[k[this_quark][3]][2]) == 24 ||  abs(k[k[this_quark][3]][2]) == 21 ) {
         if ( abs(k[k[this_quark][3]][2]) == 24 ||  abs(k[k[this_quark][3]][2]) == 21 ) {

	        //! first parent of the 94 was the boson => the quark is the second one:

           this_quark=k[this_quark][6];
         } else {

	       //! and v.v.

           this_quark=k[this_quark][3];
         }

	    //! get the right pdg for future stepping

         quark_pdg=k[this_quark][2];

       }
     } 
     if (   k[this_quark][2] == 94 ) {
       streamlog_out(DEBUG0) << "this_quark is a 94, its parent = " << k[this_quark][3] 
                             << ", parent's PDG  = " << k[k[this_quark][3]][2]  
                             << ", quark_pdg  = " << quark_pdg << std::endl;


	 //! the current line is (still) a 94. Select which of the two parents of it
	 //! is to be followed. It is the one with the right flavour. A special case
	 //! is if the parents are tops: then we must either already be following a top -
	 //! in which case there is no problem - or a b. In the latter case the right top
	 //! to follow is the one with the right sign,
        
       if (  k[k[this_quark][3]][2] == quark_pdg ||
               (abs(k[k[this_quark][3]][2]) == 6 &&
		(k[k[this_quark][3]][2]>0 ? 5 : -5 ) ==  quark_pdg )) {

         this_quark =  k[this_quark][3];
       } else {
         this_quark =  k[this_quark][6];
       }

       streamlog_out(DEBUG0) << "this_quark was a 94, new this_quark = " << this_quark
                             << ", its PDG  = " << k[this_quark][2] << std::endl;

            //! save this line - it might turnout to be the parent of the last 94 (looking 
            //! towards the beginning of the history)

       last_94_parent =  this_quark;
       if (  this_quark ==  this_quark_save ) {
        streamlog_out(ERROR) << " ERROR: didnt find parant of 94 ?! " << this_quark << std::endl;
        streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ;
        streamlog_out(ERROR) << std::endl ; 

       }
     }
	    //! jet-assignment if requested. Set to -ve to indicate that we are in the
	    //! oarton shower,

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
     streamlog_out(DEBUG0) << "  end of while loop: this_quark = " << this_quark 
                           << ", its PDG = " << k[this_quark][2] 
                           << ", its parent = " << k[this_quark][3] 
                           << ", parent's PDG  = " << abs(k[k[this_quark][3]][2])  
                           << ", quark_pdg  = " << quark_pdg << std::endl;
   }

     //! we are here because we've reached the end of the chain. If this happened
     //! becasue the parent was a boson, we recurse

   if ( abs(k[k[this_quark][3]][2]) == 21  || abs(k[k[this_quark][3]][2]) == 24 ) {
       //! recurse.

       //! jet-number set to 0 => do not assign jet numbers anymore. The untimate first
       //! partic will go into boson_ancestor, parent of the last 94 into boson_last_94_parent.

     streamlog_out(DEBUG0) << "in while loop calling first_parton with: k[this_quark][3] = " << k[this_quark][3] << std::endl;
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
   if ( last_94_parent == 0 )  {last_94_parent=k[this_quark][4];}

   first_partic = this_quark;


}
int TrueJet::flavour(int k2)
{
   int k2l;

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
   for ( int kk=1 ; kk<=njet ; kk++ ) {
     tE[kk]=0;
     for ( int i=0 ; i<3 ; i++ ) {
          tmom[kk][i]=0. ;
     } 
   }
   for ( int i=1 ; i<=nlund ; i++){ 
     if ( jet[i]>0 ) {
       if ( k[i][1] == 11 ) {
         double Ekid=0;
         for ( int jj=k[i][4] ; jj <= k[i][5] ; jj++ ) {
           Ekid+= p[jj][4];
         }
         if ( abs(Ekid-p[i][4])/p[i][4] > 0.001 ) {
           streamlog_out(WARNING) << " Particle " << i << " has energy " << p[i][4] << 
                 " , but the sum of its genstat1 kids is " << Ekid << std::endl;
           streamlog_out(WARNING) << " indicating that Geant did something fishy and un-documented in MCParticle " 
                << std::endl ;
           streamlog_out(WARNING) << " (Counter-meassures taken, so it should be OK.) " << std::endl ;
           streamlog_out(WARNING) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl;
           streamlog_out(WARNING) << std::endl ; 

           if ( k[i][1] < 30 ) {
             k[i][1]=1;
           }
           for ( int jj=k[i][4] ; jj <= k[i][5] ; jj++ ) {

             if ( k[jj][1] < 30 ) {
               k[jj][1]=0;
             }
           }
       
         }
         
       } else if  ( k[i][1] == 0 ) {
         for ( unsigned idat=0 ; idat <   mcp_pyjets[i]->getDaughters().size()  ; idat++ ) {
           int jdat=mcp_pyjets[i]->getDaughters()[idat]->ext<MCPpyjet>();
           if ( k[jdat][1] < 30 ) {
             k[jdat][1]=0;
           }
         }
         
       }
     }
   }
   for ( int i=1 ; i<=nlund ; i++){ 
     ijet=abs(jet[i]);

     if ( k[i][1]%30 == 1 ) { // true quanaties for this jet (all true stable) NB. Checks HEPEVT code,
                             // not wether there are daughters. This is right, since Mokka does not
                             // change the HEPEVT status in MCParticle, even if an interaction with
                             // seeable daughters takes place. To check: "un-used daughters", ie.
                             // the opposite case, where the generator has decayed a particle,
                             // but the decay happens after Mokka has destroyed the parent.
                             // Also: generator-decays of longlived charged-particles: probaly
                             // (chek that) the MCParticle info is after application of the
                             // B-field, so the direction of the momentum has changed.
                             // [Finally: the crossing-angle. To what has it been applied? 
                             //  everything it seems, so that should be OK]
       //if ( mcp_pyjets[i]->getGeneratorStatus() != 0 ) {
         tmom[ijet][0]+=p[i][1] ;
         tmom[ijet][1]+=p[i][2] ;
         tmom[ijet][2]+=p[i][3] ;
         tE[ijet]+=p[i][4] ;
	 //}
       }
    }
      //!!   Add info on jet-type to return (1=from string,2=lepton,3=from cluster,4=isr + x00=comes from boson, 
      //!!   boson from jet x) and indication of the companion. 

      //   type=[(1,kk=1,2*nstr),(2,kk=2*nstr+1, 2*nstr+n_hard_lepton),(3,kk=2*nstr+n_hard_lepton+1, 2*nstr+n_hard_lepton+2*nclu), 4,kk=2*nstr+n_hard_lepton+2*nclu+1,njet)]
      //   FORALL (kk=1:njet-nisr) companion(kk)=kk+1
      //   FORALL (kk=1:njet-nisr,( MODULO(kk,2) == 0 )) companion(kk)=kk-1
      //  ! initial color-singlett
      //    WHERE ( fafp_boson(1:njet) /= 0 ) 
      //      fafp(1:njet) = fafp_boson(1:njet)
      //    ELSEWHERE
      //      fafp(1:njet) = fafp_last(1:njet)
      //    ENDWHERE


    for ( int kk=1 ; kk <= njet ; kk++ ) {
      if ( kk <= 2*nstr ){
       type[kk]=1;
      } else if ( kk<= 2*nstr+n_hard_lepton ){
       type[kk]=2;
      } else if ( kk<=  2*nstr+n_hard_lepton+2*nclu ){
       type[kk]=3;
      } else if ( kk<= 2*nstr+n_hard_lepton+2*nclu+nisr ){
       type[kk]=4;
      } else if ( kk<= njet ){
       type[kk]=5;
      } 
      if ( boson[kk] != 0 ) {
        type[kk]=abs(jet[boson[kk]])*100+type[kk];
      }
      if ( fafp_boson[kk] != 0 ) {
        fafp[kk]=fafp_boson[kk];
        if ( k[ fafp_boson[kk] ] [4] != 0 ) {
          if ( k [k[ fafp_boson[kk] ] [4]][2] != 94 ) {
            fafp[kk]=fafp_last[kk];
          }
        }
      } else {
        fafp[kk]=fafp_last[kk];
      }

    }
    for  (int kk=1 ; kk<=njet-nisr-n_beam_jet; kk++) {
      if ( kk%2 == 0 ) {
        companion[kk]=kk-1;
      } else {
        companion[kk]=kk+1;
      }
    }


      

    n_djb = 0 ;       
    n_dje = 0;  
    //int dijet_end[26] ; 
    //int group[26][26] ; int dijet_begining[26] ; 
     // FORALL (KK=1:njet) group(0,kk) = COUNT((k(fafp(1:njet),4)==k(fafp(kk),4) .AND.k(fafp(kk),4) /= 0  ))
     // FORALL (kk=1:njet) group(1: group(0,kk) ,kk) = pack(index(1:njet),(k(fafp(1:njet),4)==k(fafp(kk),4) .AND.k(fafp(kk),4) /= 0 ))
    for ( int kk=1 ; kk<= njet ; kk++ ) {
      group[0][kk]=0;
      for ( int jj=1 ; jj<= njet ; jj++ ) {
   
        if (k[fafp[jj]][4]==k[fafp[kk]][4] && k[fafp[kk]][4] != 0  ) {
          group[0][kk]++;
          group[group[0][kk]][kk] = jj;
        }
      }
    }

    for ( int  kk=1; kk<= njet ; kk++ ){
      if ( type[kk]%100 <= 3 ) {
        //!  possible special case: the original fermions goes directly into the final string, without
        //!  passing a 94 (first condition) or there is no copy, ie. the first f-anti-f pair already is the elementon
        if ( group[0][kk] <= 1 ) {
          if ( ( k[fafp[kk]][4] == elementon[kk] && k[fafp[companion[kk]]][4] == elementon[companion[kk]] ) ||
               ( fafp[kk]  == elementon[kk] && fafp[companion[kk]] == elementon[companion[kk]]) ) {
            //! ... this is the case
           group[0][kk]=2 ;  group[0][companion[kk]]=2 ; 
           group[1][kk] = std::min(kk,companion[kk]);
           group[2][kk] = std::max(kk,companion[kk]) ;
           group[1][companion[kk]] = std::min(kk,companion[kk]);
           group[2][companion[kk]] = std::max(kk,companion[kk]) ;
          }
        }
      }
    }

      //! figure out the different jet combinations. (If there are no bosons, it's straigh forward:
      //! it's just the odd+even combinations. Also if one only cares about the final singlet
      //! the same aplies, boson or not)

    for ( int kk=1; kk<=njet ; kk++ ) {
      dijet_begining[kk] = 0;
        // DO jj=1,n_djb
        //   IF ( ALL(group(1:group(0,kk),kk)==jets_begin(1:group(0,kk),jj)) )  dijet_begining(kk) = jj
        // ENDDO
      if ( group[0][kk] > 0 )  {
        for ( int  jj=1; jj<=n_djb ; jj++ ) {
          int gotit=1;
          for ( int ii=1 ; ii<= group[0][kk] ; ii++ ) {
            if  ( ! (group[ii][kk]==jets_begin[ii][jj] )) {
              gotit=0;
            }
          }
          if ( gotit ) {
            dijet_begining[kk] = jj;
          }
        }
      }
      if ( dijet_begining[kk] == 0 && group[0][kk] > 0 ) {
        n_djb++;
        dijet_begining[kk] = n_djb;
        for ( int jj=0 ; jj<= group[0][kk] ; jj++ ){
          jets_begin[jj][n_djb] = group[jj][kk];
        }
      }
      dijet_end[kk] = 0;
      for ( int jj=1 ; jj<=n_dje ; jj++ ) {
        if ( std::min(kk,companion[kk])==jets_end[1][jj] &&
             std::max(kk,companion[kk])==jets_end[2][jj] ) {
          dijet_end[kk] = jj;
        }
      }
      if ( dijet_end[kk] == 0 ) {
        n_dje++;
        dijet_end[kk] = n_dje; 
        if ( companion[kk] > 0 ) {
          jets_end[1][n_dje] = std::min(kk,companion[kk]);
          jets_end[2][n_dje] = std::max(kk,companion[kk]);
          jets_end[0][n_dje] = 2 ;
        } else {
          jets_end[1][n_dje] = kk;
          jets_end[0][n_dje] = 1;
        }
      }
    }   

    nboson=0;n_jetless=0; n_0_E_jets=0; n_2jet_clu=0;
    for ( int kk=1 ; kk<=njet ; kk++ ) {
      if (type[kk]/100 != 0 ){
        nboson++;
      }
      if (tE[kk]<0.000001 ) {
	n_0_E_jets++;
      }
      if ( jet[kk]==0 && k[kk][1]==1 ) {
	n_jetless++;
      }
      if ( (type[kk]%100== 3) && ( tE[kk]> 0.0000001 ) &&
	   ( tE[companion[kk]] >  0.0000001)) {
	n_2jet_clu++;
      }
    }
    n_mixed=0;
    for ( int kk=1 ; kk<=njet-1 ; kk+=2 ) {  
      if ((type[kk]%100==1) &&
          ( (k[fafp[kk]][4] !=  k[fafp[kk+1]][4]) && ( k[fafp[kk]][4] != elementon[kk]  ||  k[fafp[kk+1]][4] != elementon[kk+1]  ) ) ) {
	n_mixed++;
      }
    }  
    //   jsum=jets_summary_t(njet,n_djb,n_dje,2*nstr, n_hard_lepton, nboson,  2*nclu,nisr , n_mixed,  n_jetless, n_2jet_clu, n_0_E_jets,jt(1:njet),dijb(1:n_djb),dije(1:n_dje)) 

}

void TrueJet::fix94()
{ int nodd, line , nCMshowers, nCandidate_94s;
  int odd_lines[4011];
  int  CMshowers[4011];
  int  Candidate_94s[4011];

  first_line=1;
  fix_top() ;
  //  nodd=COUNT((K(first_line:nlund,1)==11 .AND. K(first_line:nlund,4)==0 .AND. K(first_line:nlund,2)/=21))
  //  odd_lines=PACK ( index,(K(first_line:nlund,1)==11 .AND. K(first_line:nlund,4)==0 .AND. K(first_line:nlund,2)/=21 )) 
  //  nCMshowers=COUNT((K(first_line:nlund,2)==94))
  //  CMshowers=PACK ( index,(K(first_line:nlund,2)==94)) 

   nodd = 0 ;  nCMshowers=0 ;
   for ( int kk=first_line ; kk<=nlund ; kk++ ) {
     if ( k[kk][1] == 11 && k[kk][4] == 0 && k[kk][2] != 21 ) {
       nodd++;
       odd_lines[nodd]=kk;
     }
     if ( k[kk][2] == 94 ) {
       nCMshowers++;
       CMshowers[nCMshowers]=kk;
     }
   }


  for ( int kk=1 ; kk<=nodd; kk++ ) {
    line=odd_lines[kk];
    //  nCandidate_94s=COUNT([( ( ANY( k(k(CMshowers(ipi),4:5),2) == k(line,2)) ),ipi=1,nCMshowers)])
    //  Candidate_94s=PACK(CMshowers(1:nCMshowers),[( ( ANY( k(k(CMshowers(ipi),4:5),2) == k(line,2)) ),ipi=1,nCMshowers)])
     nCandidate_94s=0 ;
     for ( int ipi=1 ; ipi<=nCMshowers ; ipi++ ) {
       if (  k[k[CMshowers[ipi]][4]][2] == k[line][2] ||  k[k[CMshowers[ipi]][5]][2] == k[line][2] ) {
         nCandidate_94s++;
         Candidate_94s[nCandidate_94s]=CMshowers[ipi] ;
       }
     }
     //	std::cout << "  nCandidate_94s "<<  nCandidate_94s << std::endl;
     if ( nCandidate_94s == 1 )  {

         if ( k[Candidate_94s[1]][3] != line ) {
           k[Candidate_94s[1]][6] = line ;
         } else {
           for ( int ll=first_line ; ll<=nlund ; ll++ ) {
             if ( k[ll][4] == Candidate_94s[1] ) {
               k[Candidate_94s[1]][6] = ll;
               break;
             }
           }
         }
         k[line][4] = Candidate_94s[1] ;
         k[line][5] = Candidate_94s[1] ;
       //  ! only one candidate, should be fine. 
       //    But check, anyhow, that there is *some*
       //         ! line that has this candiadte 94 as a daugter. If not, print a warning (has
       //         ! never happend during my checks)
       // !!!! 
       //
       //       IF ( COUNT(k(first_line:nlund,4)== Candidate_94s(1) ) == 1 ) THEN
       //         k(line,4:5) = Candidate_94s(1)
       //         second_parent(Candidate_94s(1)) = line
       //       ELSE
       //         print *, ' hum 1 ', COUNT(k(first_line:nlund,4)== Candidate_94s(1) ), Candidate_94s(1)
       //       ENDIF    
     } else {
       int gotit=0; 
       int mothers=0;
       for ( int ipi=1; ipi<= nCandidate_94s ; ipi++ ){
         // IF ( COUNT(k(first_line:nlund,4)== Candidate_94s(ipi) ) == 1 ) THEN
           // mothers=PACK ( index,(K(first_line:nlund,4)==Candidate_94s(ipi)))
         double Ekids=p[Candidate_94s[ipi]][4];
         int nda=0 ;
         for ( int jj=first_line ; jj<= nlund ; jj++ ) {
           if ( k[jj][4] ==  Candidate_94s[ipi] ) {
             nda++;
             mothers=jj;
           }
         }
         if ( nda == 1 ) {

           //    ! ... this candidate 94 hasn't already two aughters, so it might be the one.

           double Eparents=p[mothers][4]+p[line][4] ;
           if ( ((k[mothers][2]  == k[k[Candidate_94s[ipi]][4]][2] &&
                 k[line][2]     == k[k[Candidate_94s[ipi]][5]][2] )  ||
                 (k[mothers][2] == k[k[Candidate_94s[ipi]][5]][2] && 
                 k[line][2]     == k[k[Candidate_94s[ipi]][4]][2] )) &&
                 ( abs(Ekids- Eparents)/Ekids < 0.0001 ) ) {

               // Right flavour, right mass : say no more

             if ( k[Candidate_94s[ipi]][3] != line ) {
               k[Candidate_94s[ipi]][6] = line ;
             } else {
               for ( int ll=first_line ; ll<=nlund ; ll++ ) {
                 if ( k[ll][4] == Candidate_94s[ipi] ) {
                   k[Candidate_94s[ipi]][6] = ll;
                   break;
                 }
               }
             }
             k[line][4] = Candidate_94s[ipi] ;
             k[line][5] = Candidate_94s[ipi] ;
             gotit=1;
             break;
           }
         }
       }

       if ( gotit == 0 ) {
         streamlog_out(ERROR) << " hum 2 " << line << " " << k[line][4]  << " " << k[line][5] << std::endl;
         streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ;
         streamlog_out(ERROR) << std::endl ; 

         // call pylist(2)

       }
     }

     if (  k[line][4] == 0 &&  ( k[line][2] != 21 ||  p[line][4] > 0.0000001) ) {
          //! ( this never happens in my tests)
       streamlog_out(ERROR) << " ERROR: line "<< line <<" still odd " << k[line][2] << " " << p[line][4] <<std::endl;
       streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
       streamlog_out(ERROR) << std::endl ; 

          //  CALL pylist(2)
          //  is_odd = .TRUE.
     }
  }
    // ! fill in second_parent for lines where all was OK already

  for ( int kk=1 ; kk <=nCMshowers ; kk++ ){
    line=CMshowers[kk] ;
    int ipi=0;
    if (  k[line][6] == 0 ) {
      ipi=k[line][3];
      for (int jj=-5 ; jj<=5 ; jj++ ) {
        if (jj !=0 ) {
            // ! look 5 lines before and after the mother of this 94
          if ( k[ipi+jj][4] ==  k[ipi][4] ) {
            k[line][6]=ipi+jj ;
            break ;
          }
        }
      }
    }
    if (   k[line][6]== 0 ) {
        //! (this never happens, so the +/- 5 lines seems to be OK)
      streamlog_out(ERROR) << " Second parent not found " << line  <<std::endl;
      streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ; 
      streamlog_out(ERROR) << std::endl ; 

    }

  }

}

void TrueJet::fix_top()
{
    int ipi;
    int this_fla[3];
    first_line=1 ;
    int first_top = 0 ; 

    //	std::cout << " enter "<< std::endl;
      //    IF ( ANY( abs(k(1:nlund,2)) == 6 ) ) THEN 


    int kkk=1 ; 
    while ( abs(k[kkk][2]) != 6 && kkk <= nlund ) {
      kkk++ ;
    }
    first_top=kkk;

    if ( first_top < nlund ) {
      if ( (k[first_top][4] == k[first_top+1][4] &&  k[first_top][5] == k[first_top+1][5]) &&
            k[first_top][4] ==  k[first_top][5]  &&  k[k[first_top][4]][2] == 94 ) {

         //! the two top's goes into a 94: Make it correct:
         //!  mother of first two tops set to 3 (from 0) and
         //!  change the mother of the W from heaven-knows-what
         //!  to the corresponding top, and change the daughter
         //! of it to th4 94 the correspondig b goes into.
        
         k[first_top][3]=3 ;  k[first_top+1][3]=3 ;

         for ( int kk=0 ; kk <= 1 ; kk++ ) {
           ipi=first_top+kk ;
           this_fla[1]=k[ipi][2]; this_fla[2]= ( k[ipi][2] > 0 ? 5 : -5 ) ;
           while (abs(k[ipi][2]) != 5 ) { // ! ie. not yet at the b quark
              // ! move forward, either to ...
             if ( k[k[ipi][4]][2] == 94 ) {
                 // ! a 94
               ipi = k[ipi][4] ;
             } else {
                // ! or a t or b quark of the right sign
                //    this_fla=[k(ipi,2), sign(5,k(ipi,2))]
               if ( k[k[ipi][4]][2] == this_fla[1] ||  k[k[ipi][4]][2] == this_fla[2]  ) {
                 ipi = k[ipi][4];
               } else if  ( k[k[ipi][5]][2] == this_fla[1] ||  k[k[ipi][5]][2] == this_fla[2]  ) {
                 ipi = k[ipi][5];
               } else {
                 streamlog_out(ERROR)<< " CA, alors ! " << std::endl ;
                 streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ;
                 streamlog_out(ERROR) << std::endl ; 
 
               }
             }
             if (  abs(k[ipi][2]) != 5 &&  abs(k[ipi][2]) != 6 &&  abs(k[ipi][2]) != 94 ) {
               streamlog_out(ERROR)<< " Where IS the quark ?  " << std::endl ;
               streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ;
               streamlog_out(ERROR) << std::endl ; 
             }
           }
             //! we should now be at the b (no pun intended), and the next line
             //! is the W
           if ( abs(k[ipi+1][2]) == 24 ) {
             k[ipi+1][3]=k[ipi][3];
             k[ipi+1][4]=k[ipi][4];
             k[ipi+1][5]=k[ipi][5];
           } else {
             streamlog_out(ERROR) << " Not a W ?! " << std::endl;
             streamlog_out(ERROR) << " Event: " << evt->getEventNumber() << ",   run:  " << evt->getRunNumber() << std::endl ;
             streamlog_out(ERROR) << std::endl ; 

           }
         }
      } else {
         //! first top didn't go into a 94: add a correct 94-structure
         //! at the end of pyjets: first tops goes into the added 94.
         //! the 94 goes into a copy - after first correcting k(...,5) to
         //! be k(...,4)+1 (ie. the W rather than the b) - of the original top, 
         //! with mother changed to the 94. The mother of the original products
         //! of the top (b and W) are changed to the new copy. The new 94 gets
         //! it's p from the sum of the two tops, and k=[ 11, 94, first_top, nlund+2, nlund+3 ]
         //! the mother of the original tops are changed from 0 to 3.
 
        for(int kk=first_top; kk<=first_top+1 ; kk++ ) {
            //! W as second daughter 
          k[kk][5]=k[kk][4]+1;
            //! W goes to same 94 as the b:
          k[k[kk][5]][4]=k[k[kk][4]][4];
          k[k[kk][5]][5]=k[k[kk][4]][5];
    
            //! copy to end
          for ( int jj=1 ; jj<=5 ; jj++ ) {
            p[nlund+(2-first_top+kk)][jj]=p[kk][jj];
            k[nlund+(2-first_top+kk)][jj]=k[kk][jj];
          }
          k[nlund+(2-first_top+kk)][6]=0;
            // ! adjust mother
          k[nlund+(2-first_top+kk)][3]=nlund+1 ;
            // ! change mother of decendants to be the copy
          k[k[kk][4]][3] = nlund+(2-first_top+kk) ;
          k[k[kk][5]][3] = nlund+(2-first_top+kk) ;
            // ! change daughters of original to be the copy
          k[kk][4] = nlund+1;
          k[kk][5] = nlund+1;
            // ! mother of original = 3
          k[kk][3] = 3;
        }
 
          // ! construct the 94-object
           //  p(nlund+1,1:4)=p(nlund+2,1:4)+p(nlund+3,1:4)
           //  p(nlund+1,5)= sqrt(p(nlund+1,4)**2-SUM(p(nlund+1,1:3)**2))
           //  k(nlund+1,1:5)= [ 11, 94, first_top, nlund+2, nlund+3 ]
        double sum_psq=0.0;
        for ( int jj=1 ; jj<=4 ; jj++ ) {
          p[nlund+1][jj]=p[nlund+2][jj]+p[nlund+3][jj];
          if ( jj <= 3 ) {
            sum_psq+= p[nlund+1][jj]* p[nlund+1][jj];
          } 
        }
        k[nlund+1][1]= 11;
        k[nlund+1][2]= 94;
        k[nlund+1][3]= first_top;
        k[nlund+1][4]= nlund+2;
        k[nlund+1][5]= nlund+3 ;
        k[nlund+1][6]= 0 ;
        p[nlund+1][5]= sqrt(p[nlund+1][4]*p[nlund+1][4]-sum_psq);
        nlund=nlund+3;
      } 
      first_line=first_top;
  
  streamlog_out(DEBUG1) << " HEPEVT relation table (after fixtop)"  <<std::endl;
  streamlog_out(DEBUG1) << "       line      status    pdg       parent    first     last     second " << std::endl;
  streamlog_out(DEBUG1) << "                                             daughter daugther    parent " << std::endl;
  for ( int j=1 ; j<=nlund ; j++ ) {
    streamlog_out(DEBUG1) << std::setw(10) << j << std::setw(10) <<  k[j][1] << 
                 std::setw(10) <<  k[j][2] << std::setw(10) <<  k[j][3] << 
                 std::setw(10) <<  k[j][4] << std::setw(10) << k[j][5]  << std::setw(10) << k[j][6] <<std::endl;
  }
    } else {
      first_line=1;
    }
}


void TrueJet::check( LCEvent *  /*evt*/ ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TrueJet::end(){ 

    //   std::cout << "TrueJet::end()  " << name() 
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;

}
