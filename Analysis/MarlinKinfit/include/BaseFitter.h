/*! \file 
 *  \brief Declares class BaseFitter
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.3  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 *
 */ 

#ifndef __BASEFITTER_H
#define __BASEFITTER_H

#include<vector>

class BaseFitObject;
class BaseConstraint;

//  Class BaseConstraint:
/// Abstract base class for fitting engines of kinematic fits
/**
 *
 * Author: Jenny List, Benno List
 * Last update: $Date: 2007-10-30 15:51:14 $
 *          by: $Author: gaede $
 *
 */
class BaseFitter {
  public:
    BaseFitter() {};
    virtual ~BaseFitter() {};
    virtual void addFitObject (BaseFitObject* fitobject_) {
      fitobjects.push_back(fitobject_);
    };
    virtual void addFitObject (BaseFitObject& fitobject_) {
      fitobjects.push_back(&fitobject_);
    };
    virtual void addConstraint (BaseConstraint* constraint_) {
      constraints.push_back(constraint_);
    };
    virtual void addConstraint (BaseConstraint& constraint_) {
      constraints.push_back(&constraint_);
    };
    virtual std::vector<BaseFitObject*>* getFitObjects() {return &fitobjects;};
    virtual std::vector<BaseConstraint*>* getConstraints() {return &constraints;};
    virtual double fit() = 0;
    virtual int getError() const = 0;
    virtual double getProbability() const = 0;
    virtual double getChi2() const = 0;
    virtual int   getIterations() const = 0;
  
    virtual void reset() {
      fitobjects.resize(0);
      constraints.resize(0);
    }  
  
  protected:
    virtual bool initialize() = 0;
    
    typedef std::vector <BaseFitObject *> FitObjectContainer;
    typedef std::vector <BaseConstraint *> ConstraintContainer;
    
    typedef FitObjectContainer::iterator FitObjectIterator;
    
    FitObjectContainer fitobjects;
    ConstraintContainer constraints;

};

#endif // __BASEFITTER_H
