/*! \file 
 *  \brief Implements class BaseTracer
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: BaseTracer.cc,v $
 * - Revision 1.1  2009/09/01 09:48:13  blist
 * - Added tracer mechanism, added access to fit covariance matrix
 * -
 * -
 *
 */ 

#include "BaseTracer.h"

BaseTracer::BaseTracer(): next (0) {}

BaseTracer::~BaseTracer() {}

void BaseTracer::initialize (BaseFitter& fitter) {
  if (next) next->initialize (fitter);
}

void BaseTracer::step (BaseFitter& fitter) {
  if (next) next->step (fitter);
}

void BaseTracer::substep (BaseFitter& fitter, int flag) {
  if (next) next->substep (fitter, flag);
}

void BaseTracer::finish (BaseFitter& fitter) {
  if (next) next->finish (fitter);
}

void BaseTracer::setNextTracer (BaseTracer *next_) {
  next = next_;
}

void BaseTracer::setNextTracer (BaseTracer& next_) {
  next = &next_;
}

BaseTracer *BaseTracer::getNextTracer () {
  return next;
}
