/*! \file 
 *  \brief Declares class BaseHardConstraint
 *
 * \b Changelog:
 * - 12.2.08 First version
 *
 * \b CVS Log messages:
 * - $Log: BaseHardConstraint.h,v $
 * - Revision 1.2  2011/03/03 15:03:02  blist
 * - Latest version, with NewFitterGSL
 * -
 * - Revision 1.1  2008/02/12 16:43:25  blist
 * - First Version of Soft Constraints
 * -
 * -
 *
 */ 

#ifndef __BASEHARDCONSTRAINT_H
#define __BASEHARDCONSTRAINT_H

#include "BaseConstraint.h"

class BaseFitObject;

//  Class BasehardConstraint:
/// Abstract base class for constraints of kinematic fits
/**
 * This class defines the minimal functionality any constraint class must provide. 
 * First of all a constraint should know on with particles (or FitObject) it is applied. 
 * Where as for example a constraint on the total transvese momentum takes into
 * account all particles in the event, an invariant mass constraint usually applies only
 * to a subset of particles. 
 *
 * The particle list is implemented as a vector containing pointers to objects derived 
 * from BaseFitObject and can be either set a whole (setFOList) or enlarged by adding
 * a single BaseFitObject (addToFOList).
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
 * Author: Jenny List, Benno List
 * Last update: $Date: 2011/03/03 15:03:02 $
 *          by: $Author: blist $
 *
 */

class BaseHardConstraint: public BaseConstraint {
  public:
    
    /// Virtual destructor
    virtual ~BaseHardConstraint();
  
    /// Adds first order derivatives to global covariance matrix M
    virtual void add1stDerivativesToMatrix (double *M,      ///< Global covariance matrix, dimension at least idim x idim
                                            int idim        ///< First dimension of array der
                                            ) const = 0;
       
    /// Adds second order derivatives to global covariance matrix M
    virtual void add2ndDerivativesToMatrix (double *M,      ///< Global covariance matrix, dimension at least idim x idim
                                            int idim,       ///< First dimension of array der
                                            double lambda   ///< Lagrange multiplier for this constraint
                                            ) const = 0;
    /// Add lambda times derivatives of chi squared to global derivative vector
    virtual void addToGlobalChi2DerVector (double *y,   ///< Vector of chi2 derivatives
                                           int idim,    ///< Vector size 
                                           double lambda //< The lambda value
                                           ) const = 0;
    /// Calculate directional derivative 
    virtual double dirDer                 (double *p,   ///< Vector of direction
                                           double *w,   ///< Work vector
                                           int idim,    ///< Vector size 
                                           double mu=1  ///< optional multiplier
                                          );
   
    /// Calculate directional derivative for abs(c)
    virtual double dirDerAbs              (double *p,   ///< Vector of direction
                                           double *w,   ///< Work vector
                                           int idim,    ///< Vector size 
                                           double mu=1  ///< optional multiplier
                                          );
   
    /// Accesses position of constraint in global constraint list
    virtual int  getGlobalNum() const = 0;
    /// Sets position of constraint in global constraint list
    virtual void setGlobalNum (int iglobal                  ///< The global constraint number
                              ) = 0;
                                 
};


#endif // __BASEHARDCONSTRAINT_H
