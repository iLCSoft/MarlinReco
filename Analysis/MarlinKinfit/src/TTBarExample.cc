#include "TTBarExample.h"
#include <iostream>
#include <vector>
#include <string>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>
#endif

#include "UTIL/LCRelationNavigator.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
#include "NeutrinoFitObject.h"
#include "PxConstraint.h"
#include "PyConstraint.h"
#include "PzConstraint.h"
#include "EConstraint.h"
#include "OPALFitter.h"
#include "TwoB4JPairing.h"
#include "FourJetPairing.h"
#include "MassConstraint.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace CLHEP ;


TTBarExample aTTBarExample ;


TTBarExample::TTBarExample() : Processor("TTBarExample") {
  
  // modify processor description
  _description = "TTBarExample does a 6C fit on 6jet ttbar events (Px, Py, Pz, E, M12 = M34 = M_W (for all permutations assuming the b-jets are tagged))" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "LightJetCollectionName" , 
			   "Name of the Light Jet collection"  ,
			   _lightjetcolName ,
			   std::string("DurhamLightJets") ) ;
                           
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "BJetCollectionName" , 
			   "Name of the B Jet collection"  ,
			   _bjetcolName ,
			   std::string("DurhamBJets") ) ;
                           
  registerProcessorParameter( "ECM" ,
                              "Center-of-Mass Energy",
                              _ecm,
                              (float)500.);

}


void TTBarExample::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void TTBarExample::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void TTBarExample::processEvent( LCEvent * evt ) { 

    
    message<MESSAGE>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
  // this gets called for every event 
  // usually the working horse ...

#ifdef MARLIN_USE_AIDA
  
  // define a histogram pointer
  static AIDA::IHistogram1D* hRecTopMassBest ;    
  static AIDA::IHistogram1D* hRecTopMassAll ;    
  static AIDA::IHistogram1D* hRecTopMassNoFitBest ;    
  static AIDA::IHistogram1D* hRecTopMassNoFitAll ;    
  static AIDA::IHistogram1D* hFitProbBest ;    
  static AIDA::IHistogram1D* hFitProbAll ;    
             
    message<MESSAGE>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
  
  if( isFirstEvent() ) { 
    
    hRecTopMassBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassBest", "M_W", 200, 0., 200. ) ; 
    hRecTopMassAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassAll", "M_W", 200, 0., 200. ) ; 
    hRecTopMassNoFitBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassNoFitBest", "M_W", 200, 0., 200. ) ; 
    hRecTopMassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassNoFitAll", "M_W", 200, 0., 200. ) ; 
    hFitProbBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProb", "fit probability", 100, 0., 1. ) ; 
    hFitProbAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProbAll", "fit probability", 100, 0., 1. ) ; 

  }

#endif
   
  message<MESSAGE>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
  
  
  HepLorentzVector lvec;
  
    
  // fill histogram from LCIO data :

  //////////////////   JETS ///////////////////////////
   
     LCCollection* lightjetcol = evt->getCollection( _lightjetcolName ) ;
     LCCollection* bjetcol = evt->getCollection( _bjetcolName ) ;
     if (lightjetcol != 0 && bjetcol != 0) {
  
       int nlightJETS = lightjetcol->getNumberOfElements()  ;
       message<MESSAGE>( log() 
                      << " found " << nlightJETS
                      << " light jets in event" << evt->getEventNumber() 
                      << "  in run "           << evt->getRunNumber() 
                      ) ;
                   
       int nbJETS = bjetcol->getNumberOfElements()  ;
       message<MESSAGE>( log() 
                      << " found " << nbJETS
                      << " b jets in event" << evt->getEventNumber() 
                      << "  in run "           << evt->getRunNumber() 
                      ) ;
                   
       
  // original fit objects - save for next permutation
       JetFitObject* j1 = 0;
       JetFitObject* j2 = 0;
       JetFitObject* j3 = 0;
       JetFitObject* j4 = 0;
       // these are assumed to be the b-jets
       JetFitObject* j5 = 0;
       JetFitObject* j6 = 0;
       
       double erre = 1.0;        //   100%/sqrt(E)
       double errtheta = 0.01;   //   10mrad
       double errphi = 0.01;     //   10mrad
       
       
       for(int i=0; i< nlightJETS ; i++){
         
          ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>( lightjetcol->getElementAt( i ) ) ;
               
          if (j) {
             message<MESSAGE>( log() 
                       << " found jet in event " << evt->getEventNumber() 
                       << "  in run "          << evt->getRunNumber() 
                       ) ;
             lvec = HepLorentzVector ((j->getMomentum())[0],(j->getMomentum())[1],(j->getMomentum())[2],j->getEnergy()); 
             erre *= std::sqrt(lvec.e());
             if (i == 0) {
               j1 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               message<MESSAGE>( log() 
                       << " start four-vector of first  jet: " << *j1 
                       ) ;
             }
             else if (i == 1) {
               j2 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               message<MESSAGE>( log() 
                       << " start four-vector of second  jet: " << *j2 
                       ) ;
             }
             else if (i == 2) {
               j3 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               message<MESSAGE>( log() 
                       << " start four-vector of third  jet: " << *j3 
                       ) ;
             }
             else if (i == 3) {
               j4 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               message<MESSAGE>( log() 
                       << " start four-vector of forth  jet: " << *j4 
                       ) ;
             }
           
          }
       }
       
       for(int i=0; i< nbJETS ; i++){
         
          ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>( bjetcol->getElementAt( i ) ) ;
               
          if (j) {
             message<MESSAGE>( log() 
                       << " found b-jet in event " << evt->getEventNumber() 
                       << "  in run "          << evt->getRunNumber() 
                       ) ;
             lvec = HepLorentzVector ((j->getMomentum())[0],(j->getMomentum())[1],(j->getMomentum())[2],j->getEnergy()); 
             erre *= std::sqrt(lvec.e());
             if (i == 0) {
               j5 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               message<MESSAGE>( log() 
                       << " start four-vector of first b-jet: " << *j5 
                       ) ;
             }
             else if (i == 1) {
               j6 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               message<MESSAGE>( log() 
                       << " start four-vector of second b-jet: " << *j6 
                       ) ;
             }
           
          }
       }
       
       const int NJETS = 6;
       
       // these get changed by the fit -> reset after each permutation!
       JetFitObject fitjets[NJETS] = {*j1, *j2, *j3, *j4, *j5, *j6};
 
       // these point allways to the fitjets array, which gets reset.
       JetFitObject *jets[NJETS];
       for (int i = 0; i < NJETS; ++i) jets[i] = &fitjets[i];

       TwoB4JPairing pairing (jets);
       JetFitObject *permutedjets[NJETS];
       
       double bestprob = 0.;
       double bestmass1 = 0, bestmass2 = 0;
       double beststartmass1 = 0, beststartmass2 = 0;
       double startmass1 = 0, startmass2 = 0;
 
       for (int iperm = 0; iperm < pairing.getNPerm(); iperm++) {
     
         message<MESSAGE>( log() 
                       << " ================================================= "  
                       ) ;
         message<MESSAGE>( log() 
                       << " iperm = " << iperm 
                       ) ;

         // important: (re-)set fitjets array!
         fitjets[0] = *j1;
         fitjets[1] = *j2;
         fitjets[2] = *j3;
         fitjets[3] = *j4;
         fitjets[4] = *j5;
         fitjets[5] = *j6;

         pairing.nextPermutation (permutedjets);
         for (int i = 0; i < NJETS; ++i) {
            message<MESSAGE>( log() 
                       << "start four-vector of jet " << i << ": " << *(permutedjets[i])
                       ) ;
         }              
        
         PxConstraint pxc;
         for (int i = 0; i < NJETS; ++i)
            pxc.addToFOList (*(permutedjets[i]));
        
         PyConstraint pyc;
         for (int i = 0; i < NJETS; ++i)
            pyc.addToFOList (*(permutedjets[i]));
        
         PzConstraint pzc;
         for (int i = 0; i < NJETS; ++i)
            pzc.addToFOList (*(permutedjets[i]));
            
         message<MESSAGE>( log() 
                   << "ECM = " << _ecm
                       ) ;
         EConstraint ec(_ecm);
         for (int i = 0; i < NJETS; ++i)
            ec.addToFOList (*(permutedjets[i]));
        
         message<MESSAGE>( log() 
                    << "Value of pxc before fit: " << pxc.getValue()
                    ) ;
         message<MESSAGE>( log() 
                    << "Value of pyc before fit: " << pyc.getValue()
                    ) ;
         message<MESSAGE>( log() 
                    << "Value of pzc before fit: " << pzc.getValue()
                    ) ;
         message<MESSAGE>( log() 
                    << "Value of ec before fit: " << ec.getValue()
                    ) ;
  
         MassConstraint w1(80.4);
         MassConstraint w2(80.4);
         w1.addToFOList (*(permutedjets[0]), 1);
         w1.addToFOList (*(permutedjets[1]), 1);
         w2.addToFOList (*(permutedjets[2]), 1);
         w2.addToFOList (*(permutedjets[3]), 1);
         
         message<MESSAGE>( log() 
                       << "start mass of W 1: " << w1.getMass(1)
                       ) ;
         message<MESSAGE>( log() 
                       << "start mass of W 2: " << w2.getMass(1)
                       ) ;
                       
         // this is just a cheap way to monitor the resulting top mass:
         MassConstraint t1(175.);
         t1.addToFOList (*(permutedjets[0]), 1);
         t1.addToFOList (*(permutedjets[1]), 1);
         t1.addToFOList (*(permutedjets[4]), 1);
         MassConstraint t2(175.);
         t2.addToFOList (*(permutedjets[2]), 1);
         t2.addToFOList (*(permutedjets[3]), 1);
         t2.addToFOList (*(permutedjets[5]), 1);
         startmass1 = t1.getMass(1);
         startmass2 = t2.getMass(1);
         message<MESSAGE>( log() 
                       << "start mass of top 1: " << startmass1
                       ) ;
         message<MESSAGE>( log() 
                       << "start mass of top 2: " << startmass2
                       ) ;
#ifdef MARLIN_USE_AIDA
         hRecTopMassNoFitAll->fill( startmass1 ) ;
         hRecTopMassNoFitAll->fill( startmass2 ) ;
#endif        
       
         OPALFitter fitter;
         for (int i = 0; i < NJETS; ++i)
            fitter.addFitObject (*(permutedjets[i]));
         fitter.addConstraint (pxc);
         fitter.addConstraint (pyc);
         fitter.addConstraint (pzc);
         fitter.addConstraint (ec);
         fitter.addConstraint (w1);
         fitter.addConstraint (w2);

         double prob = fitter.fit();
         message<MESSAGE>( log() 
                       << "fit probability = " << prob 
                       ) ;
         message<MESSAGE>( log() 
                       << "error code: " << fitter.getError() 
                       ) ;
         for (int i = 0; i < NJETS; ++i) {
            message<MESSAGE>( log() 
                       << "final four-vector of jet " << i << ": " << *(permutedjets[i])
                       ) ;
         }              
         
         message<MESSAGE>( log() 
                       << "final mass of W 1: " << w1.getMass(1)
                       ) ;
         message<MESSAGE>( log() 
                       << "final mass of W 2: " << w2.getMass(1)
                       ) ;
         message<MESSAGE>( log() 
                       << "final mass of top 1: " << t1.getMass(1)
                       ) ;
         message<MESSAGE>( log() 
                       << "final mass of top 2: " << t2.getMass(1)
                       ) ;
         if (fitter.getError() == 0) {
#ifdef MARLIN_USE_AIDA
           hFitProbAll->fill( prob ) ;
           hRecTopMassAll->fill( t1.getMass(1)) ;
           hRecTopMassAll->fill( t2.getMass(1)) ;
#endif       
           if (prob > bestprob) {
             bestprob = prob;
             bestmass1 = t1.getMass(1);
             bestmass2 = t2.getMass(1);
             beststartmass1 = startmass1;
             beststartmass2 = startmass2;
           }
         }
         else {
         message<MESSAGE>( log() 
                       << "FIT ERROR = " << fitter.getError() << ", not filling histograms!"
                       ) ;
         }

       }

#ifdef MARLIN_USE_AIDA
       if (bestprob > 0) {
         hFitProbBest->fill( bestprob ) ;
         hRecTopMassBest->fill( bestmass1 ) ;
         hRecTopMassBest->fill( bestmass2 ) ;
         hRecTopMassNoFitBest->fill( beststartmass1 ) ;
         hRecTopMassNoFitBest->fill( beststartmass2 ) ;

       } 
#endif       

       delete j1;
       delete j2;
       delete j3;
       delete j4;
       delete j5;
       delete j6;
     }
    


  _nEvt ++ ;
}



void TTBarExample::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TTBarExample::end(){ 
  
}

