/*! \file 
 *  \brief Declares class SoftBWMassConstraint
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.1  2008/10/21 08:36:50  blist
 * - Added classes SoftBWParticleConstraint, SoftBWMassConstraint
 * -
 *
 */ 

#ifndef __SOFTBWMASSCONSTRAINT_H
#define __SOFTBWMASSCONSTRAINT_H

#include "SoftBWParticleConstraint.h"

class ParticleFitObject;

//  Class SoftBWMassConstraint:
/// Implements constraint 0 = mass1 - mass2 - m
/**
 * This class implements different mass constraints:
 * - the invariant mass of several objects should be m
 * - the difference of the invariant masses between two 
 *   sets of objects should be m (normally m=0 in this case).
 *
 * Author: Jenny List, Benno List
 * Last update: $Date: 2008-11-24 11:01:01 $
 *          by: $Author: beckmann $
 *
 */

class SoftBWMassConstraint : public SoftBWParticleConstraint {
  public:
  
    /// Constructor
    SoftBWMassConstraint (double gamma_,    ///< The Gamma value
                          double mass_ = 0.   ///< The mass difference between object sets 1 and 2
                   );
    /// Virtual destructor             
    virtual ~SoftBWMassConstraint();
    
    /// Returns the value of the constraint function
    virtual double getValue() const;
    
    /// Get first order derivatives. 
    /// Call this with a predefined array "der" with the necessary number of entries!
    virtual void getDerivatives(int idim,      ///< First dimension of the array
                                double der[]   ///< Array of derivatives, at least idim x idim 
                               ) const;
    
    /// Get the actual invariant mass of the fit objects with a given flag                         
    virtual double getMass (int flag = 1       ///< The flag
                           );
    
    /// Sets the target mass of the constraint
    virtual void setMass (double mass_           ///< The new mass
                         );
    
  
  protected:
    double mass;   ///< The mass difference between object sets 1 and 2

  
    /// Second derivatives with respect to the 4-vectors of Fit objects i and j; result false if all derivatives are zero 
    virtual bool secondDerivatives (int i,                        ///< number of 1st FitObject
                                    int j,                        ///< number of 2nd FitObject
                                    double *derivatives           ///< The result 4x4 matrix 
                                   ) const;
    /// First derivatives with respect to the 4-vector of Fit objects i; result false if all derivatives are zero 
    virtual bool firstDerivatives (int i,                        ///< number of 1st FitObject
                                   double *derivatives           ///< The result 4-vector
                                  ) const;
};

#endif // __SOFTBWMASSCONSTRAINT_H
