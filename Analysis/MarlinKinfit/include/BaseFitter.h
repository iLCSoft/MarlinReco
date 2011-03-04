/*! \file 
 *  \brief Declares class BaseFitter
 *
 * \b Changelog:
 * - 12.2.08 BL: Add soft constraints
 *
 * \b CVS Log messages:
 * - $Log: BaseFitter.h,v $
 * - Revision 1.5  2011/03/03 15:03:02  blist
 * - Latest version, with NewFitterGSL
 * -
 * - Revision 1.4  2009/09/01 09:48:12  blist
 * - Added tracer mechanism, added access to fit covariance matrix
 * -
 * - Revision 1.3  2008/02/13 12:37:37  blist
 * - new file BaseFitter.cc
 * -
 * - Revision 1.2  2008/02/12 16:43:25  blist
 * - First Version of Soft Constraints
 * -
 * - Revision 1.1  2008/02/12 10:19:05  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.6  2008/01/30 21:48:02  blist
 * - Newton Fitter still doesnt work :-(
 * -
 * - Revision 1.5  2008/01/30 09:14:53  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.4  2008/01/29 17:18:35  blist
 * - added method getDoF
 * -
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
#include<string>
#include<map>

class BaseFitObject;
class BaseConstraint;
class BaseHardConstraint;
class BaseSoftConstraint;
class BaseTracer;

//  Class BaseConstraint:
/// Abstract base class for fitting engines of kinematic fits
/**
 *
 * Author: Jenny List, Benno List
 * Last update: $Date: 2011/03/03 15:03:02 $
 *          by: $Author: blist $
 *
 */
class BaseFitter {
  public:
    BaseFitter();
    virtual ~BaseFitter();
    virtual void addFitObject (BaseFitObject* fitobject_);
    virtual void addFitObject (BaseFitObject& fitobject_);
    virtual void addConstraint (BaseConstraint* constraint_);
    virtual void addConstraint (BaseConstraint& constraint_);
    virtual void addHardConstraint (BaseHardConstraint* constraint_);
    virtual void addHardConstraint (BaseHardConstraint& constraint_);
    virtual void addSoftConstraint (BaseSoftConstraint* constraint_);
    virtual void addSoftConstraint (BaseSoftConstraint& constraint_);
    virtual std::vector<BaseFitObject*>* getFitObjects();
    virtual std::vector<BaseHardConstraint*>* getConstraints();
    virtual std::vector<BaseSoftConstraint*>* getSoftConstraints();
    virtual double fit() = 0;
    virtual int getError() const = 0;
    virtual double getProbability() const = 0;
    virtual double getChi2() const = 0;
    virtual int    getDoF() const = 0;
    virtual int   getIterations() const = 0;
  
    virtual void reset();
    virtual bool initialize() = 0;
    
    virtual BaseTracer *getTracer();
    virtual const BaseTracer *getTracer() const;
    virtual void setTracer(BaseTracer *newTracer
                          );
    virtual void setTracer(BaseTracer& newTracer
                          );
    
    virtual const double *getGlobalCovarianceMatrix (int& idim ///< 1st dimension of global covariance matrix
                                                          ) const;                 
    virtual double *getGlobalCovarianceMatrix (int& idim ///< 1st dimension of global covariance matrix
                                              );                 
  
  protected:
    /// Copy constructor disabled
    BaseFitter (const BaseFitter& rhs);
    /// Assignment disabled
    BaseFitter& operator= (const BaseFitter& rhs);
    
    
    typedef std::vector <BaseFitObject *> FitObjectContainer;
    typedef std::vector <BaseHardConstraint *> ConstraintContainer;
    typedef std::vector <BaseSoftConstraint *> SoftConstraintContainer;
    
    typedef FitObjectContainer::iterator FitObjectIterator;
    typedef ConstraintContainer::iterator ConstraintIterator;
    typedef SoftConstraintContainer::iterator SoftConstraintIterator;
    
    FitObjectContainer      fitobjects;
    ConstraintContainer     constraints;
    SoftConstraintContainer softconstraints;    
    
    int     covDim;   ///< dimension of global covariance matrix
    double *cov;      ///< global covariance matrix of last fit problem
    bool    covValid; ///< Flag whether global covariance is valid

#ifndef FIT_TRACEOFF    
    BaseTracer *tracer;
#endif 
  public:  
#ifndef FIT_TRACEOFF    
    std::map<std::string, double> traceValues;
#endif 

};

#endif // __BASEFITTER_H
