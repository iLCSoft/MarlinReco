#include "WW5CFit.h"
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
//#include <root/TLorentzVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
#include "PxConstraint.h"
#include "PyConstraint.h"
#include "PzConstraint.h"
#include "EConstraint.h"
#include "OPALFitter.h"
#include "FourJetPairing.h"
#include "MassConstraint.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace CLHEP ;


WW5CFit aWW5CFit ;


WW5CFit::WW5CFit() : Processor("WW5CFit") {
  
  // modify processor description
  _description = "WW5CFit does a 5C fit on 4 jet events (Px, Py, Pz, E, M12 = M34 (for all three permutations))" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "JetCollectionName" , 
			   "Name of the Jet collection"  ,
			   _jetcolName ,
			   std::string("Durham2Jets") ) ;
                           
  registerProcessorParameter( "ECM" ,
                              "Center-of-Mass Energy",
                              _ecm,
                              (float)500.);

}


void WW5CFit::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void WW5CFit::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void WW5CFit::processEvent( LCEvent * evt ) { 

    
    message<MESSAGE>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
  // this gets called for every event 
  // usually the working horse ...

#ifdef MARLIN_USE_AIDA
  
  // define a histogram pointer
  static AIDA::IHistogram1D* hRecWMassBest ;    
  static AIDA::IHistogram1D* hRecWMassAll ;    
  static AIDA::IHistogram1D* hRecWMassNoFitBest ;    
  static AIDA::IHistogram1D* hRecWMassNoFitAll ;    
  static AIDA::IHistogram1D* hFitProbBest ;    
  static AIDA::IHistogram1D* hFitProbAll ;    
             
    message<MESSAGE>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
  
  if( isFirstEvent() ) { 
    
    hRecWMassBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassBest", "M_W", 200, 0., 200. ) ; 
    hRecWMassAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassAll", "M_W", 200, 0., 200. ) ; 
    hRecWMassNoFitBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassNoFitBest", "M_W", 200, 0., 200. ) ; 
    hRecWMassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassNoFitAll", "M_W", 200, 0., 200. ) ; 
    hFitProbBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProb", "fit probability", 100, 0., 1. ) ; 
    hFitProbAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProbAll", "fit probability", 100, 0., 1. ) ; 

  }
   
  message<MESSAGE>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
  
  
  HepLorentzVector lvec;
  
    
  // fill histogram from LCIO data :

  //////////////////   JETS ///////////////////////////
   
     LCCollection* jetcol = evt->getCollection( _jetcolName ) ;
     if (jetcol != 0) {
  
       int nJETS = jetcol->getNumberOfElements()  ;
       message<MESSAGE>( log() 
                      << " found " << nJETS
                      << " jets in event" << evt->getEventNumber() 
                      << "  in run "          << evt->getRunNumber() 
                      ) ;
                   
       float yminus = jetcol ->parameters().getFloatVal( "YMinus");              
       message<MESSAGE>( log() 
                      << " yminus = " << yminus
                      ) ;
       float yplus = jetcol ->parameters().getFloatVal( "YPlus");              
       message<MESSAGE>( log() 
                      << " yplus = " << yplus
                      ) ;
       
  // original fit objects - save for next permutation
       JetFitObject* j1 = 0;
       JetFitObject* j2 = 0;
       JetFitObject* j3 = 0;
       JetFitObject* j4 = 0;
       
       double erre = 1.0;        //   100%/sqrt(E)
       double errtheta = 0.01;   //   10mrad
       double errphi = 0.01;     //   10mrad
       
       
       for(int i=0; i< nJETS ; i++){
         
          ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>( jetcol->getElementAt( i ) ) ;
               
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
       
       const int NJETS = 4;
       
       // these get changed by the fit -> reset after each permutation!
       JetFitObject fitjets[NJETS] = {*j1, *j2, *j3, *j4};
 
       // these point allways to the fitjets array, which gets reset.
       JetFitObject *jets[NJETS];
       for (int i = 0; i < NJETS; ++i) jets[i] = &fitjets[i];

       FourJetPairing pairing (jets);
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

         pairing.nextPermutation (permutedjets);
         for (int i = 0; i < NJETS; ++i)
            message<MESSAGE>( log() 
                       << "start four-vector of jet " << i << ": " << *(permutedjets[i])
                       ) ;
        
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
  
         MassConstraint w(0.);
         w.addToFOList (*(permutedjets[0]), 1);
         w.addToFOList (*(permutedjets[1]), 1);
         w.addToFOList (*(permutedjets[2]), 2);
         w.addToFOList (*(permutedjets[3]), 2);
        
         startmass1 = w.getMass(1);
         startmass2 = w.getMass(2);
         message<MESSAGE>( log() 
                       << "start mass of W 1: " << startmass1
                       ) ;
         message<MESSAGE>( log() 
                       << "start mass of W 2: " << startmass2
                       ) ;
         hRecWMassNoFitAll->fill( startmass1 ) ;
         hRecWMassNoFitAll->fill( startmass2 ) ;
        
       
         OPALFitter fitter;
         for (int i = 0; i < NJETS; ++i)
            fitter.addFitObject (*(permutedjets[i]));
         fitter.addConstraint (pxc);
         fitter.addConstraint (pyc);
         fitter.addConstraint (pzc);
         fitter.addConstraint (ec);
         fitter.addConstraint (w);

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
                       << "final mass of W 1: " << w.getMass(1)
                       ) ;
         message<MESSAGE>( log() 
                       << "final mass of W 2: " << w.getMass(2)
                       ) ;
         if (fitter.getError() == 0) {
           hFitProbAll->fill( prob ) ;
           hRecWMassAll->fill( w.getMass(1)) ;
           hRecWMassAll->fill( w.getMass(2)) ;
           if (prob > bestprob && w.getMass(1) > 50 && w.getMass(1) < 110) {
             bestprob = prob;
             bestmass1 = w.getMass(1);
             bestmass2 = w.getMass(2);
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

       if (bestprob > 0) {
         hFitProbBest->fill( bestprob ) ;
         hRecWMassBest->fill( bestmass1 ) ;
         hRecWMassBest->fill( bestmass2 ) ;
         hRecWMassNoFitBest->fill( beststartmass1 ) ;
         hRecWMassNoFitBest->fill( beststartmass2 ) ;

       } 

       delete j1;
       delete j2;
       delete j3;
       delete j4;
     }
    


#endif


  _nEvt ++ ;
}



void WW5CFit::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void WW5CFit::end(){ 
  
}

