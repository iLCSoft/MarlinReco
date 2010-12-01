////////////////////////////////////////////////////////////////
// Class TopEventILC
//
// Author: Benno List, Jenny Boehme
// Last update: $Date: 2008/09/26 09:58:11 $
//          by: $Author: boehmej $
// 
// Description: class to generate and fit top pair events at ILC
//               
////////////////////////////////////////////////////////////////
#include "TopEventILC.h"

#include "JetFitObject.h"
#include "NeutrinoFitObject.h"
#include "MassConstraint.h"
#include "cernlib.h"

#include <iostream>              // - cout
#include <cmath>            
using std::cout;
using std::endl;
using std::abs;
using std::sqrt;

// constructor: 
TopEventILC::TopEventILC()
: leptonic (false),
  pxc (1, 0),
  pyc (0, 1),
  pzc (0, 0, 1),
  ec  (0, 0, 0, 1, 500)
  {
  for (int i = 0; i < NFV; ++i) fv[i] = 0;
  for (int i = 0; i < NBFO; ++i) bfo[i] = bfosmear[i] = 0;
  w1.setMass(80.4);
  w2.setMass(80.4);
  pxc.setName ("px=0");
  pyc.setName ("py=0");
  pzc.setName ("pz=0");
  ec.setName  ("E=500");
  w.setName ("top-equalmass");
  w1.setName ("w1-mass");
  w2.setName ("w2-mass");
}

//destructor: 
TopEventILC::~TopEventILC() {
  for (int i = 0; i < NFV; ++i) delete fv[i];
  for (int i = 0; i < NBFO; ++i) {
    delete bfo[i];
    delete bfosmear[i];
  }  
}

// Generate Breit-Wigner Random number
double TopEventILC::bwrandom (double r, double e0, double gamma, double emin, double emax) const {
  double a = atan (2.0*(emax - e0)/gamma);
  double b = atan (2.0*(emin - e0)/gamma);
  return e0 + 0.5*gamma*tan (r*(a - b) + b);
}


// generate four vectors
void TopEventILC::genEvent(){
  // generate 4-vectors of top-decay:
  // 0: top-top-system -> top1 top2
  // 1: top 1 -> W1 b1
  // 2: top 2 -> W2 b2
  // 3: W1 -> j11 j12
  // 4: W2 -> j21 j22
  // 5: b1
  // 6: j11
  // 7: j12
  // 8: b2
  // 9: j21
  //10: j22 or neutrino
  
  double mtop = 174;
  double gammatop = 1.4;
  double mW   = 80.4;
  double gammaW = 2.1;
  double mb   = 5.0;
  double mj   = 1.0;
  double Ecm = 500;
      
  FReal rw[4];
  ranmar (rw, 4);
  
  FourVector *toppair = fv[0] = new FourVector (Ecm, 0, 0, 0);
  double mtop1 = bwrandom (rw[0], mtop, gammatop, mtop-3*gammatop, mtop+3*gammatop);
  double mtop2 = bwrandom (rw[1], mtop, gammatop, mtop-3*gammatop, mtop+3*gammatop);
  FourVector *top1 = fv[1] = new FourVector (mtop1, 0, 0, 0);
  FourVector *top2 = fv[2] = new FourVector (mtop2, 0, 0, 0);
  
  toppair->decayto (*top1, *top2);
  // cout << "top 1: m=" << mtop1 << " = " << top1->getM() << endl;
  // cout << "top 2: m=" << mtop2 << " = " << top2->getM() << endl;
  
  double mw1 = bwrandom (rw[2], mW, gammaW, mW-3*gammaW, mW+3*gammaW);
  double mw2 = bwrandom (rw[3], mW, gammaW, mW-3*gammaW, mW+3*gammaW);
  
  FourVector *W1 = fv[3] = new FourVector (mw1, 0, 0, 0);
  FourVector *W2 = fv[4] = new FourVector (mw2, 0, 0, 0);
  FourVector *b1 = fv[5] = new FourVector (mb, 0, 0, 0);
  FourVector *b2 = fv[8] = new FourVector (mb, 0, 0, 0);
  // cout << "W 1: m=" << mw1 << " = " << W1->getM() << endl;
  // cout << "W 2: m=" << mw2 << " = " << W2->getM() << endl;
  
  top1->decayto (*W1, *b1);
  top2->decayto (*W2, *b2);
  
  FourVector *j11 = fv[6]  = new FourVector (mj, 0, 0, 0);
  FourVector *j12 = fv[7]  = new FourVector (mj, 0, 0, 0);
  FourVector *j21 = fv[9]  = new FourVector (mj, 0, 0, 0);
  FourVector *j22 = fv[10] = new FourVector (mj, 0, 0, 0);
  
  W1->decayto (*j11, *j12);
  W2->decayto (*j21, *j22);
  
  double Eresolhad = 0.35;     // 35% / sqrt (E)
  double Eresolem = 0.10;     // 10% / sqrt (E)
  double thetaResol = 0.1;  // rad
  double phiResol = 0.1;    // rad
  
  double pxtot=0;
  double pytot=0;
  double pztot=0;
  
  for (int j = 0; j < 6; ++j) {
    int i = j+5;
    double E = fv[i]->getE();
    double theta = fv[i]->getTheta();
    double phi = fv[i]->getPhi();
    double EError = (j==4 && leptonic) ? Eresolem*sqrt(E) : Eresolhad*sqrt(E);
    
    static const char *names[] = {"b1", "j11", "j12", "b2", "j21", "j22"};
    // Create fit object with true quantities for later comparisons
    bfo[j] = new JetFitObject (E, theta, phi, EError, thetaResol, phiResol, 0);
    bfo[j]->setName (names[j]);
    
    FReal randoms[3];
    rnorml (randoms, 3);
    
    // Create fit object with smeared quantities as fit input
    double ESmear = E + EError*randoms[0];
    double thetaSmear = theta + thetaResol*randoms[1];
    double phiSmear = phi + phiResol*randoms[2];
    
    
    if (j != 5 || !leptonic) {
      bfosmear[j] = new JetFitObject (ESmear, thetaSmear, phiSmear, EError, thetaResol, phiResol, 0.);
      bfosmear[j]->setName (names[j]);
      pxtot += bfosmear[j]->getPx();
      pytot += bfosmear[j]->getPy();
      pztot += bfosmear[j]->getPz();
    }
    else {
      double pxn = -pxtot;
      double pyn = -pytot;
      double pzn = -pztot;
      double en =  sqrt (pxn*pxn+pyn*pyn+pzn*pzn);
      double phi=atan2(pyn, pxn);
      double theta = acos (pzn/en);
      bfosmear[5] = new NeutrinoFitObject (en, theta, phi, 10, 0.2, 0.2);
    
      bfosmear[5]->setName ("n22");
      // cout << "Neutrino: px=" << bfosmear[5]->getPx() << ", py=" << bfosmear[5]->getPy() 
      //      << ", pz=" << bfosmear[5]->getPz() << endl;
    }
    fvsmear[i] = new FourVector (bfosmear[j]->getE(), bfosmear[j]->getPx(), bfosmear[j]->getPy(), bfosmear[j]->getPz());
    
    pxc.addToFOList (*bfosmear[j]);
    pyc.addToFOList (*bfosmear[j]);
    pzc.addToFOList (*bfosmear[j]);
    ec.addToFOList (*bfosmear[j]);
    w.addToFOList (*bfosmear[j], j<3?1:2);
      
  }
  fvsmear[3] = new FourVector (*fvsmear[6]+*fvsmear[7]);
  fvsmear[4] = new FourVector (*fvsmear[9]+*fvsmear[10]);
  fvsmear[1] = new FourVector (*fvsmear[3]+*fvsmear[5]);
  fvsmear[2] = new FourVector (*fvsmear[4]+*fvsmear[8]);
  fvsmear[0] = new FourVector (*fvsmear[1]+*fvsmear[2]);
  
  w1.addToFOList (*bfosmear[1]);
  w1.addToFOList (*bfosmear[2]);
  w2.addToFOList (*bfosmear[4]);
  w2.addToFOList (*bfosmear[5]);
  
}

// fit it!
int TopEventILC::fitEvent (BaseFitter& fitter){
  
//   for (int i = 0; i < 6; ++i) 
//     cout << "true four-vector of jet " << i << ": " << *bfo[i] << endl;
//   for (int i = 0; i < 6; ++i) 
//     cout << "initial four-vector of jet " << i << ": " << *bfosmear[i] << endl;
  
  
  // reset lists of constraints and fitobjects
  fitter.reset();
  
  int debug = 0;
  
  if (debug) {
    cout << "TopEventILC::fitEvent: ==================================================\n";
    cout << "True vectors: \n";
    for (int i = 0; i<6; ++i) {
      cout << bfo[i]->getName() << ": " << *bfo[i] << endl;
    }
    cout << "Start vectors: \n";
    for (int i = 0; i<6; ++i) {
      cout << bfosmear[i]->getName() << ": " << *bfosmear[i] << endl;
    }
    cout << "Total: \n";
    cout << "gen:   " << *fv[0] << ", m=" << fv[0]->getM() << endl;
    cout << "smear: " << *fvsmear[0] << ", m=" << fvsmear[0]->getM() << endl;
    cout << "Top1: \n";
    cout << "gen:   " << *fv[1] << ", m=" << fv[1]->getM() << endl;
    cout << "smear: " << *fvsmear[1] << ", m=" << fvsmear[1]->getM() << endl;
    cout << "Top2: \n";
    cout << "gen:   " << *fv[2] << ", m=" << fv[2]->getM() << endl;
    cout << "smear: " << *fvsmear[2] << ", m=" << fvsmear[2]->getM() << endl;
    cout << "W1: \n";
    cout << "gen:   " << *fv[3] << ", m=" << fv[3]->getM() << endl;
    cout << "smear: " << *fvsmear[3] << ", m=" << fvsmear[3]->getM() << endl;
    cout << "W2: \n";
    cout << "gen:   " << *fv[4] << ", m=" << fv[4]->getM() << endl;
    cout << "smear: " << *fvsmear[4] << ", m=" << fvsmear[4]->getM() << endl;
  }
  
   
  for (int i = 0; i < 6; i++) {
    assert (bfosmear[i]);
    fitter.addFitObject (*bfosmear[i]);
//     pxc.addToFOList (*bfosmear[i]);
//     pyc.addToFOList (*bfosmear[i]);
//     w.addToFOList (*bfosmear[i], i<3?1:2);
  }
//   
//   w1.addToFOList (*bfosmear[1]);
//   w1.addToFOList (*bfosmear[2]);
//   w2.addToFOList (*bfosmear[4]);
//   w2.addToFOList (*bfosmear[5]);
  
  fitter.addConstraint (pxc);
  fitter.addConstraint (pyc);
  fitter.addConstraint (pzc);
  fitter.addConstraint (ec);
  fitter.addConstraint (w);
  fitter.addConstraint (w1);
  fitter.addConstraint (w2);
  
  double prob = fitter.fit();
  
//   cout << "fit probability = " << prob << endl;
//   for (int i = 0; i < 6; ++i) 
//     cout << "final four-vector of jet " << i << ": " << *bfosmear[i] << endl;
//   cout << "Constraint pxc: " << pxc.getValue() << endl;
//   cout << "Constraint pxc: " << pxc.getValue() << endl;
//   cout << "Constraint w:   " << w.getValue() << ", top mass: " << w.getMass() << endl;
//   cout << "Constraint w1:  " << w1.getValue() << ", W mass: " << w1.getMass() << endl;
//   cout << "Constraint w2:  " << w2.getValue() << ", W mass: " << w2.getMass() << endl;
   cout << "fit probability = " << prob << ", top mass: " << w.getMass() << ", dof: " << fitter.getDoF() 
        << ", iterations: " << fitter.getIterations() << endl;

  for (int j = 0; j < 6; ++j) {
    int i = j+5;
    fvfinal[i] = new FourVector (bfosmear[j]->getE(), bfosmear[j]->getPx(), bfosmear[j]->getPy(), bfosmear[j]->getPz());
  }
  
  fvfinal[3] = new FourVector (*fvfinal[6]+*fvfinal[7]);
  fvfinal[4] = new FourVector (*fvfinal[9]+*fvfinal[10]);
  fvfinal[1] = new FourVector (*fvfinal[3]+*fvfinal[5]);
  fvfinal[2] = new FourVector (*fvfinal[4]+*fvfinal[8]);
  fvfinal[0] = new FourVector (*fvfinal[1]+*fvfinal[2]);
  
  if (debug) {
    cout << "===============After Fiting ===================================\n";
    cout << "Final vectors: \n";
    for (int i = 0; i<6; ++i) {
      cout << bfosmear[i]->getName() << ": " << *bfo[i] << endl;
    }
    cout << "Total: \n";
    cout << "gen:   " << *fv[0] << ", m=" << fv[0]->getM() << endl;
    cout << "final: " << *fvfinal[0] << ", m=" << fvfinal[0]->getM() << endl;
    cout << "Top1: \n";
    cout << "gen:   " << *fv[1] << ", m=" << fv[1]->getM() << endl;
    cout << "final: " << *fvfinal[1] << ", m=" << fvfinal[1]->getM() << endl;
    cout << "Top2: \n";
    cout << "gen:   " << *fv[2] << ", m=" << fv[2]->getM() << endl;
    cout << "final: " << *fvfinal[2] << ", m=" << fvfinal[2]->getM() << endl;
    cout << "W1: \n";
    cout << "gen:   " << *fv[3] << ", m=" << fv[3]->getM() << endl;
    cout << "final: " << *fvfinal[3] << ", m=" << fvfinal[3]->getM() << endl;
    cout << "W2: \n";
    cout << "gen:   " << *fv[4] << ", m=" << fv[4]->getM() << endl;
    cout << "final: " << *fvfinal[4] << ", m=" << fvfinal[4]->getM() << endl;
    cout << "================================================\n";
  }
  

   
   return fitter.getError();


}
