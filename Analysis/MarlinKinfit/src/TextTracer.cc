/*! \file 
 *  \brief Implements class TextTracer
 *
 * \b Changelog:
 *
 *
 * \b CVS Log messages:
 * - $Log: TextTracer.cc,v $
 * - Revision 1.1  2009/09/01 09:48:13  blist
 * - Added tracer mechanism, added access to fit covariance matrix
 * -
 * -
 */ 
 
#include "TextTracer.h"
#include "BaseFitter.h"
#include "BaseFitObject.h"
#include "BaseConstraint.h"
#include "BaseHardConstraint.h"
#include "BaseSoftConstraint.h"
#include <cassert>
#include <cstring>

TextTracer::TextTracer(std::ostream& os_)
: os (os_)
{}

TextTracer::~TextTracer()
{}
    

void TextTracer::initialize (BaseFitter& fitter) {
  os << "=============== Starting fit ======================\n";
  
  printFitObjects (fitter);
  printConstraints (fitter);
  
  istep = 1;
  isubstep = 0;
  
  BaseTracer::initialize (fitter);
}

void TextTracer::step (BaseFitter& fitter) {
  isubstep = 1;
  os << "--------------- Step " << istep << " --------------------\n";
  
  printFitObjects (fitter);
  printConstraints (fitter);
  
  ++istep;
  BaseTracer::step (fitter);
}

void TextTracer::substep (BaseFitter& fitter, int flag) {
  os << "---- Substep " << istep << "." << isubstep << " ----\n";
  
  printFitObjects (fitter);
  printConstraints (fitter);

  ++isubstep;
  BaseTracer::substep (fitter, flag);
}

void TextTracer::finish (BaseFitter& fitter) {

  os << "=============== Final result ======================\n";
  printFitObjects (fitter);
  printConstraints (fitter);
  
  os << "=============== Finished fit ======================\n";


  BaseTracer::finish (fitter);
}

void TextTracer::printFitObjects (BaseFitter& fitter) {
  FitObjectContainer* fitobjects = fitter.getFitObjects();
  if (!fitobjects) return;
  os << "Fit objects:\n";
  for (FitObjectIterator i = fitobjects->begin(); i != fitobjects->end(); ++i) {
    BaseFitObject *fo = *i;
    assert (fo);
    os << fo->getName() << ": " << *fo << ", chi2=" << fo->getChi2() << std::endl;
  }
}
void TextTracer::printConstraints (BaseFitter& fitter) {
  ConstraintContainer* constraints = fitter.getConstraints();
  if (constraints && constraints->size() > 0) {
    os << "Hard Constraints:\n";
    for (ConstraintIterator i = constraints->begin(); i != constraints->end(); ++i) {
      BaseConstraint *c = *i;
      assert (c);
      os << i-constraints->begin() << " " << c->getName() << ": " << c->getValue() << "+-" << c->getError() << std::endl;
    }
  }
  SoftConstraintContainer* softConstraints = fitter.getSoftConstraints();
  if (softConstraints && softConstraints->size() > 0) {
    os << "Soft Constraints:\n";
    for (SoftConstraintIterator i = softConstraints->begin(); i != softConstraints->end(); ++i) {
      BaseConstraint *c = *i;
      assert (c);
      os << i-softConstraints->begin() << " " << c->getName() << ": " << c->getValue() << "+-" << c->getError() << std::endl;
    }
  }
}
