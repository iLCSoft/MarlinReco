/*! \file 
 *  \brief Declares class BaseConstraint
 *
 * \b Changelog:
 * - 6.6.04 BL: First doxygen docu
 * - 12.2.08 BL: Part of functionality transferred to BaseHardConstraint
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.3  2008/02/18 09:59:34  blist
 * - MomentumConstraint and SoftGaussMomentumCOnstraint added; PConstraint is obsolete
 * -
 * - Revision 1.2  2008/02/12 16:43:25  blist
 * - First Version of Soft Constraints
 * -
 * - Revision 1.1  2008/02/12 10:19:04  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.8  2008/02/04 17:30:52  blist
 * - NewtonFitter works now!
 * -
 * - Revision 1.7  2008/01/30 09:14:53  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.6  2008/01/29 17:19:21  blist
 * - no change
 * -
 * - Revision 1.5  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 *
 */ 

#ifndef __BASECONSTRAINT_H
#define __BASECONSTRAINT_H

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
 * Last update: $Date: 2008-11-24 11:01:01 $
 *          by: $Author: beckmann $
 *
 */

class BaseConstraint {
  public:
    /// Creates an empty BaseConstraint object
    BaseConstraint();
    
    /// Copy constructor
    BaseConstraint (const BaseConstraint& rhs              ///< right hand side
                   );
    /// Assignment               
    BaseConstraint& operator= (const BaseConstraint& rhs   ///< right hand side
                              );
    
    /// Virtual destructor
    virtual ~BaseConstraint();
    
    /// Returns the value of the constraint function
    virtual double getValue() const = 0;
    
    /// Returns the error on the value of the constraint
    virtual double getError() const;
    
    /// Returns the name of the constraint
    virtual const char*getName() const;
    /// Set object's name
    virtual void setName (const char * name_);
    
    /// Get first order derivatives of the constraint function
    /// Call this with a predefined array "der" with the necessary number of entries!
    virtual void getDerivatives (int idim,                  ///< First dimension of array der
                                 double der[]               ///< Array of derivatives, dimension at least idim x idim
                                ) const = 0;
                                 
    protected:
      char *name;  
};


#endif // __BASECONSTRAINT_H
