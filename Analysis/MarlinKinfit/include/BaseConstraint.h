/*! \file 
 *  \brief Declares class BaseConstraint
 *
 * \b Changelog:
 * - 6.6.04 BL: First doxygen docu
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.5  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 *
 */ 

#ifndef __BASECONTRAINT_H
#define __BASECONTRAINT_H

class BaseFitObject;

//  Class BaseConstraint:
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
 * Last update: $Date: 2007-10-30 15:51:14 $
 *          by: $Author: gaede $
 *
 */

class BaseConstraint {
  public:
    /// Creates an empty BaseConstraint object
    BaseConstraint() {};
    /// Virtual destructor
    virtual ~BaseConstraint() {};
    
    /// Returns the value of the constraint
    virtual double getValue() const = 0;
    
    /// Get first order derivatives. 
    /// Call this with a predefined array "der" with the necessary number of entries!
    virtual void getDerivatives (int idim,                  ///< First dimension of array der
                                 double der[]               ///< Array of derivatives, dimension at least idim x idim
                                ) const = 0;
    /// Adds first order derivatives to global covariance matrix M
    virtual void add1stDerivativesToMatrix (int idim,       ///< First dimension of array der
                                            double *M       ///< Global covariance matrix, dimension at least idim x idim
                                           ) const = 0;
    /// Adds second order derivatives to global covariance matrix M
    virtual void add2ndDerivativesToMatrix (int idim,       ///< First dimension of array der
                                            double *M,      ///< Global covariance matrix, dimension at least idim x idim
                                            double lambda   ///< Lagrange multiplier for this constraint
                                            ) const= 0;
    
    /// Accesses position of constraint in global constraint list
    virtual int  getGlobalNum() const = 0;
    /// Sets position of constraint in global constraint list
    virtual void setGlobalNum (int iglobal                  ///< The global constraint number
                              ) = 0;
};


#endif // __BASECONTRAINT_H
