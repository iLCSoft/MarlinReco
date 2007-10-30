/*! \file 
 *  \brief Declares class MassConstraint
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.4  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 *
 */ 

#ifndef __MASSCONSTRAINT_H
#define __MASSCONSTRAINT_H

#include "ParticleConstraint.h"

class ParticleFitObject;

//  Class MassConstraint:
/// Implements constraint 0 = mass1 - mass2 - m
/**
 * This class implements different mass constraints:
 * - the invariant mass of several objects should be m
 * - the difference of the invariant masses between two 
 *   sets of objects should be m (normally m=0 in this case).
 *
 * Author: Jenny List, Benno List
 * Last update: $Date: 2007-10-30 15:51:14 $
 *          by: $Author: gaede $
 *
 */

class MassConstraint : public ParticleConstraint {
  public:
  
    /// Constructor
    MassConstraint (double mass_ = 0.   ///< The mass difference between object sets 1 and 2
                   );
    /// Virtual destructor             
    virtual ~MassConstraint();
    
    /// Returns the value of the constraint
    virtual double getValue() const;
    
    /// Get first order derivatives. 
    /// Call this with a predefined array "der" with the necessary number of entries!
    virtual void getDerivatives(int idim,      ///< First dimension of the array
                                double der[]   ///< Array of derivatives, at least idim x idim 
                               ) const;
    
                               
    virtual double getMass (int flag = 1);
    virtual void setMass (double mass_) {mass = mass_;};
    
    /// Adds first order derivatives to global covariance matrix M
    virtual void add1stDerivativesToMatrix(int idim,       ///< First dimension of the array
                                           double *M       ///< Covariance matrix, at least idim x idim 
                                          ) const;
    /// Adds second order derivatives, multiplied by lambda, to global covariance matrix M
    virtual void add2ndDerivativesToMatrix(int idim,      ///< First dimension of the array
                                           double *M,     ///< Covariance matrix, at least idim x idim 
                                           double lambda  ///< Factor for derivatives
                                          ) const;
  
  protected:
    double mass;   ///< The mass difference between object sets 1 and 2

};

#endif // __MASSCONSTRAINT_H
