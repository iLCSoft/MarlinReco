/*! \file 
 *  \brief Declares class EConstraint
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
 * - Revision 1.2  2007/09/13 08:09:50  blist
 * - Updated 2nd derivatives for px,py,pz,E constraints, improved header documentation
 * -
 *
 */ 
 
#ifndef __ECONSTRAINT_H
#define __ECONSTRAINT_H

#include "ParticleConstraint.h"

class ParticleFitObject;

//  Class EConstraint:
/// Sum (E) = ECM constraint
/**
 *
 * Author: Jenny List
 * Last update: $Date: 2007-10-30 15:51:14 $
 *          by: $Author: gaede $
 */
class EConstraint : public ParticleConstraint {
  public:
    /// Constructor
    EConstraint (double ecm_     ///< The center of mass energy
                );
                
    /// Virtual destructor
    virtual ~EConstraint();
    
    /// Returns the value of the constraint
    virtual double getValue() const;
    
    /// Get first order derivatives. 
    /// Call this with a predefined array "der" with the necessary number of entries!
    virtual void getDerivatives (int idim,                  ///< First dimension of array der
                                 double der[]               ///< Array of derivatives, dimension at least idim x idim
                                ) const;
    /// Adds first order derivatives to global covariance matrix M
    virtual void add1stDerivativesToMatrix (int idim,       ///< First dimension of array der
                                            double *M       ///< Global covariance matrix, dimension at least idim x idim
                                           ) const;
    /// Adds second order derivatives to global covariance matrix M
    virtual void add2ndDerivativesToMatrix (int idim,       ///< First dimension of array der
                                            double *M,      ///< Global covariance matrix, dimension at least idim x idim
                                            double lambda   ///< Lagrange multiplier for this constraint
                                            ) const;
  
    /// Accessor for the center of mass energy
    inline double getEcm() const {return ecm;};
    
  protected:
     
    double ecm;  ///< The center of mass energy

};

#endif // __ECONSTRAINT_H
