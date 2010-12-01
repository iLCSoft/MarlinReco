/*! \file 
 *  \brief Implements class RootTracer
 *
 * \b Changelog:
 *
 *
 * \b CVS Log messages:
 * - $Log: RootTracer.cc,v $
 * - Revision 1.1  2010/05/25 13:23:56  boehmej
 * - fixed RootTracer
 * -
 * - Revision 1.2  2009/09/01 08:34:24  blist
 * - Kinfit update
 * -
 * - Revision 1.1  2009/05/22 12:12:07  blist
 * - New tracers
 * -
 * -
 */ 
 
#include "RootTracer.h"
#include "BaseFitter.h"
#include "BaseFitObject.h"
#include "BaseHardConstraint.h"
#include "BaseSoftConstraint.h"
#include <cassert>
#include <cstring>
#include <iostream>
using namespace std;

RootTracer::RootTracer(const char* filename, const char *option)
: file (0), tree(0)
{
  file = new TFile (filename, option);
  tree = new TTree ("trace", "Fit Tracing");
  CreateBranches();
  SetBranchAddresses();
  eventnumber = 0;
}

RootTracer::~RootTracer()
{
  file->Write();
  file->Close();
}
    

void RootTracer::initialize (BaseFitter& fitter) {
  cout << "=============== Starting fit ======================\n";
  
  printFitObjects (fitter);
  printConstraints (fitter);
  
  ++eventnumber;
  istep = 1;
  isubstep = 0;
  
  
  TString name ("event");
  TString title ("Event ");
  name += eventnumber;
  title += eventnumber;
  eventtree = new TTree (name, title);
  CreateEventBranches (fitter);
  
  
  BaseTracer::initialize (fitter);
}

void RootTracer::step (BaseFitter& fitter) {
  isubstep = 1;
  cout << "--------------- Step " << istep << " --------------------\n";
  
  printFitObjects (fitter);
  printConstraints (fitter);
  
  stepnumber = istep;
  substepnumber = isubstep;
  chi2 = fitter.getChi2();
  cout << "chi2=" << chi2 << endl;
  FillParameterValues (fitter);
  
  tree->Fill();
  eventtree->Fill();
  
  ++istep;
  BaseTracer::step (fitter);
}

void RootTracer::substep (BaseFitter& fitter, int flag) {
  cout << "---- Substep " << istep << "." << isubstep << " ----\n";
  
  printFitObjects (fitter);
  printConstraints (fitter);

  ++isubstep;
  BaseTracer::substep (fitter, flag);
}

void RootTracer::finish (BaseFitter& fitter) {

  cout << "=============== Final result ======================\n";
  printFitObjects (fitter);
  printConstraints (fitter);
  
  cout << "=============== Finished fit ======================\n";


  BaseTracer::finish (fitter);
}

void RootTracer::printFitObjects (BaseFitter& fitter) {
  FitObjectContainer* fitobjects = fitter.getFitObjects();
  if (!fitobjects) return;
  cout << "Fit objects:\n";
  for (FitObjectIterator i = fitobjects->begin(); i != fitobjects->end(); ++i) {
    BaseFitObject *fo = *i;
    assert (fo);
    cout << fo->getName() << ": " << *fo << ", chi2=" << fo->getChi2() << std::endl;
  }
}
void RootTracer::printConstraints (BaseFitter& fitter) {
  ConstraintContainer* constraints = fitter.getConstraints();
  if (!constraints) return;
  cout << "(Hard) Constraints:\n";
  for (ConstraintIterator i = constraints->begin(); i != constraints->end(); ++i) {
    BaseHardConstraint *c = *i;
    assert (c);
    cout << i-constraints->begin() << " " << c->getName() << ": " << c->getValue() << "+-" << c->getError() << std::endl;
  }
  SoftConstraintContainer* softconstraints = fitter.getSoftConstraints();
  if (!softconstraints) return;
  cout << "Soft Constraints:\n";
  for (SoftConstraintIterator i = softconstraints->begin(); i != softconstraints->end(); ++i) {
    BaseSoftConstraint *c = *i;
    assert (c);
    cout << i-softconstraints->begin() << " " << c->getName() << ": " << c->getValue() << "+-" << c->getError() << std::endl;
  }
}   

void RootTracer::SetBranchAddresses()
{
    tree->SetBranchAddress("event",&eventnumber);
    tree->SetBranchAddress("step",&stepnumber);
    tree->SetBranchAddress("substep",&substepnumber);
   
}

void RootTracer::CreateBranches()
{
// create branch with addresses for own datamembers
  tree->Branch("event",      &eventnumber,     "event/I");      
  tree->Branch("step",       &stepnumber,      "step/I");      
  tree->Branch("substep",    &substepnumber,   "substep/I");      
}
void RootTracer::CreateEventBranches(BaseFitter& fitter)
{
// create branch with addresses for own datamembers
  eventtree->Branch("step",       &stepnumber,      "step/I");      
  eventtree->SetBranchAddress("step",&stepnumber);
  eventtree->Branch("substep",    &substepnumber,   "substep/I");      
  eventtree->SetBranchAddress("substep",&substepnumber);
  eventtree->Branch("chi2",    &chi2,   "chi2/D");      
  eventtree->SetBranchAddress("chi2",&chi2);
  FitObjectContainer* fitobjects = fitter.getFitObjects();
  if (!fitobjects) return;
  cout << "RootTracer: Fit objects:\n";
  for (FitObjectIterator i = fitobjects->begin(); i != fitobjects->end(); ++i) {
    BaseFitObject *fo = *i;
    assert (fo);
    TString foname = fo->getName();
    for (int ilocal = 0; ilocal < fo->getNPar(); ilocal++) {
      int iglobal = fo->getGlobalParNum(ilocal);
      if (iglobal >= 0 && iglobal < NPARMAX) {
        TString parname = "Par";
        parname += iglobal;
        parname += "_";
        parname += fo->getParamName(ilocal);
        cout << "Parameter " << iglobal << " '" << parname << "'\n";
        eventtree->Branch(parname,    parvalue+iglobal,   parname+"/D");      
        eventtree->SetBranchAddress(parname,parvalue+iglobal);
        
      }
    }
  }
  cout << "RootTracer: (Hard) Constraints:\n";
  ConstraintContainer* constraints = fitter.getConstraints();
  for (ConstraintIterator i = constraints->begin(); i != constraints->end(); ++i) {
    BaseHardConstraint *c = *i;
    assert (c);
    int iglobal = c->getGlobalNum();
    if (iglobal >= 0 && iglobal < NPARMAX) {
      TString cname = "Const";
//      cname += iglobal;
      cname += (i-constraints->begin());
      cname += "_";
      cname += c->getName();
      eventtree->Branch(cname,    parvalue+iglobal,   cname+"/D");      
      eventtree->SetBranchAddress(cname,parvalue+iglobal);
//      cout << "Constraint " << iglobal << " '" << cname << "'\n";
      cout << "Constraint " << (i-constraints->begin()) << " '" << cname << "'\n";
    }
  }
//   cout << "RootTracer: Soft Constraints:\n";
//   SoftConstraintContainer* softconstraints = fitter.getSoftConstraints();
//   for (SoftConstraintIterator i = softconstraints->begin(); i != softconstraints->end(); ++i) {
//     BaseSoftConstraint *c = *i;
//     assert (c);
//     int iglobal = c->getGlobalNum();
//     if (iglobal >= 0 && iglobal < NPARMAX) {
//       TString cname = "Const";
// //      cname += iglobal;
//       cname += (i-softconstraints->begin());
//       cname += "_";
//       cname += c->getName();
//       eventtree->Branch(cname,    parvalue+iglobal,   cname+"/D");      
//       eventtree->SetBranchAddress(cname,parvalue+iglobal);
// //      cout << "Constraint " << iglobal << " '" << cname << "'\n";
//       cout << "Constraint " << (i-softconstraints->begin()) << " '" << cname << "'\n";
//     }
//   }
}

void RootTracer::FillParameterValues(BaseFitter& fitter) {
  FitObjectContainer* fitobjects = fitter.getFitObjects();
  if (!fitobjects) return;
  for (FitObjectIterator i = fitobjects->begin(); i != fitobjects->end(); ++i) {
    BaseFitObject *fo = *i;
    assert (fo);
    cout << "fo " << fo->getName() << ":\n";
    for (int ilocal = 0; ilocal < fo->getNPar(); ilocal++) {
      int iglobal = fo->getGlobalParNum(ilocal);
      if (iglobal >= 0 && iglobal < NPARMAX) {
        parvalue[iglobal] =  fo->getParam(ilocal);   
        cout << "Par " << ilocal << ": global " << iglobal << " = " << parvalue[iglobal] << endl;
      }
    }
  }
  ConstraintContainer* constraints = fitter.getConstraints();
  for (ConstraintIterator i = constraints->begin(); i != constraints->end(); ++i) {
    BaseHardConstraint *c = *i;
    assert (c);
    int iglobal = c->getGlobalNum();
    if (iglobal >= 0 && iglobal < NPARMAX) {
      parvalue[iglobal] =  c->getValue();   
      cout << "Const " << c->getName() << ": global " << iglobal << " = " << parvalue[iglobal] << endl;
    }
  }
}
