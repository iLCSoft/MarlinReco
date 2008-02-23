/*! \file 
 *  \brief Declares class MomentumConstraint
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.1  2008/02/18 09:59:34  blist
 * - MomentumConstraint and SoftGaussMomentumCOnstraint added; PConstraint is obsolete
 * -
 * -
 *
 */ 
 
 
////////////////////////////////////////////////////////////////
// Class MomentumConstraint
//
// Author: Benno List, Jenny List
// Last update: $Date: 2008-02-23 11:18:39 $
//          by: $Author: listj $
// 
// Description: constraint 
// efact*sum(E_i) + pxfact*sum(p_x,i)+pyfact*sum(p_y,i)+pzfact*sum(p_z,i)-value = 0
//               
////////////////////////////////////////////////////////////////

#ifndef __MOMENTUMCONSTRAINT_H
#define __MOMENTUMCONSTRAINT_H

#include "ParticleConstraint.h"

class ParticleFitObject;

//  Class PConstraint:
/// Implements a constraint of the form efact*sum(E)+pxfact*sum(px)+pyfact*sum(py)+pzfact*sum(pz)=value
/**
 *
 * Author: Jenny List, Benno List
 * Last update: $Date: 2008-02-23 11:18:39 $
 *          by: $Author: listj $
 *
 */
class MomentumConstraint : public ParticleConstraint {
  public:
    MomentumConstraint (double efact_=0,      ///< Factor for energy sum
                        double pxfact_=0,     ///< Factor for px sum
                        double pyfact_=0,     ///< Factor for py sum
                        double pzfact_=0,     ///< Factor for pz sum
                        double value_ = 0     ///< Target value of sum
                       );
    virtual ~MomentumConstraint();
    virtual double getValue() const;
    /// Get first order derivatives. 
    /// Call this with a predefined array "der" with the necessary number of entries!
    virtual void getDerivatives(int idim,      ///< First dimension of the array
                                double der[]   ///< Array of derivatives, at least idim x idim 
                               ) const;
                                   
    virtual void addToGlobalDerMatrix (double lambda, int idim, double *M) const;
    
    virtual void invalidateCache() const;
  
  protected:
    void updateCache() const;
  
  
    double efact;
    double pxfact;
    double pyfact;
    double pzfact;
    double value;
    
    mutable bool cachevalid;
    mutable int  nparams;
  
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

#endif // __MOMENTUMCONSTRAINT_H
