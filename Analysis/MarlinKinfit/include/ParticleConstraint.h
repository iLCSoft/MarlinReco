/*! \file 
 *  \brief Declares class ParticleConstraint
 *
 * \b Changelog:
 * - 17.11.04 BL: First version (refactured from BaseConstraint)
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.4  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 * - Revision 1.3  2007/09/13 13:09:04  blist
 * - Better docs for ParticleConstraint
 * -
 * - Revision 1.2  2007/09/13 08:09:51  blist
 * - Updated 2nd derivatives for px,py,pz,E constraints, improved header documentation
 * -
 *
 */ 

#ifndef __PARTICLECONSTRAINT_H
#define __PARTICLECONSTRAINT_H

#include<vector>
#include<cassert>
#include "BaseConstraint.h"

class ParticleFitObject;

//  Class ParticleConstraint:
/// Abstract base class for constraints of kinematic fits
/**
 * This class defines the minimal functionality any constraint class must provide. 
 * First of all a constraint should know on with particles (or FitObject) it is applied. 
 * Where as for example a constraint on the total transvese momentum takes into
 * account all particles in the event, an invariant mass constraint usually applies only
 * to a subset of particles. 
 *
 * The particle list is implemented as a vector containing pointers to objects derived 
 * from ParticleFitObject and can be either set a whole (setFOList) or enlarged by adding
 * a single ParticleFitObject (addToFOList).
 *
 * From the four--momenta of all concerned fit objects the constraint has to be able 
 * to calculate its current value (getValue). Constraints should be formulated such that 
 * a value of zero corresponds to a perfectly fulfilled constraint.
 *
 * In order to find a solution to the constrained minimisation problem, fit algorithms
 * usually need the first order derivatives of the constraint with respect to the fit
 * parameters. Since many constraints can be most easily expressed in terms of E, px, py, pz,
 * the constraints supply their derivatives w.r.t. these parameters. If a FitObject uses
 * a different parametrisation, it is its own task to provide the additional derivatives
 * of  E, px, py, pz w.r.t. the parameters of the FitObject. Thus it is easily possible
 * to use FitObjects with different kinds of parametrisations under the same constraint.
 * Some fit algorithms also need the second derivatives of the constraint, 
 * i.e. the NewtonFitter. 
 *
 * First and second order derivatives of each constraint can be added directly to the 
 * global covariance matrix containing the derivatives of all constraints w.r.t. to all
 * parameters (add1stDerivativesToMatrix, add2ndDerivativesToMatrix). This requires the
 * constraint to know its position in the overall list of constraints (globalNum). 
 * 
 *
 * Author: Jenny Böhme, Benno List
 * $Date: 2007-10-30 15:51:14 $
 * $Author: gaede $
 *
 */

class ParticleConstraint: public BaseConstraint {
  public:
    /// Creates an empty ParticleConstraint object
    inline ParticleConstraint();
    /// Virtual destructor
    virtual ~ParticleConstraint() {};
    
    /// Adds several ParticleFitObject objects to the list
    virtual void setFOList(std::vector <ParticleFitObject*> *fitobjects_ ///< A list of BaseFitObject objects
                          ){
      for (int i = 0; i < (int) fitobjects_->size(); i++) {
        fitobjects.push_back ((*fitobjects_)[i]);
        flags.push_back (1);
      }  
    }; 
    /// Adds one ParticleFitObject objects to the list
    virtual void addToFOList(ParticleFitObject& fitobject, int flag = 1
                             ){
      fitobjects.push_back (&fitobject);
      flags.push_back (flag);
    }; 
    /// Returns the value of the constraint
    virtual double getValue() const = 0;
    
    /// Get first order derivatives. 
    /// Call this with a predefined array "der" with the necessary number of entries!
    virtual void getDerivatives(int idim,      ///< First dimension of the array
                                double der[]   ///< Array of derivatives, at least idim x idim 
                               ) const = 0;
    /// Adds first order derivatives to global covariance matrix M
    virtual void add1stDerivativesToMatrix(int idim,       ///< First dimension of the array
                                           double *M       ///< Covariance matrix, at least idim x idim 
                                          ) const 
    {assert (false);}
    /// Adds second order derivatives, multiplied by lambda, to global covariance matrix M
    virtual void add2ndDerivativesToMatrix(int idim,      ///< First dimension of the array
                                           double *M,     ///< Covariance matrix, at least idim x idim 
                                           double lambda  ///< Factor for derivatives
                                          ) const
    {assert (false);}
    
    /// Accesses position of constraint in global constraint list
    virtual int  getGlobalNum() const 
    {return globalNum;}
    /// Sets position of constraint in global constraint list
    virtual void setGlobalNum (int iglobal                ///< Global constraint number
                              ) 
    {globalNum = iglobal;}
    
    /// Invalidates any cached values for the next event
    virtual void invalidateCache() const 
    {}
  
  protected:
    /// Vector of pointers to ParticleFitObjects 
    typedef std::vector <ParticleFitObject *> FitObjectContainer;    
    /// Iterator through vector of pointers to ParticleFitObjects 
    typedef FitObjectContainer::iterator FitObjectIterator;
    /// Constant iterator through vector of pointers to ParticleFitObjects 
    typedef FitObjectContainer::const_iterator ConstFitObjectIterator;
    ///  The FitObjectContainer
    FitObjectContainer fitobjects;
    ///  The derivatives
    std::vector <double> derivatives;
    ///  The flags can be used to divide the FitObjectContainer into several subsets 
    ///  used for example to implement an equal mass constraint (see MassConstraint). 
    std::vector <int> flags;
    
    /// Position of constraint in global constraint list
    int globalNum;

};

ParticleConstraint::ParticleConstraint() {
  invalidateCache();
}

#endif // __PARTICLECONSTRAINT_H
