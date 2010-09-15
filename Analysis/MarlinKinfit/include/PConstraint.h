/*! \file 
 *  \brief Declares class PConstraint
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: PConstraint.h,v $
 * - Revision 1.2  2008/02/12 11:03:32  blist
 * - Added doxygen configuration
 * -
 *
 */ 
 
 
////////////////////////////////////////////////////////////////
// Class PConstraint
//
// Author: Jenny Boehme, Bennmo List
// Last update: $Date: 2008/02/12 11:03:32 $
//          by: $Author: blist $
// 
// Description: constraint 
// a*sum(px)+b*sum(py)+c*sum(pz)+d*sum(E)=e
//               
////////////////////////////////////////////////////////////////

#ifndef __PCONSTRAINT_H
#define __PCONSTRAINT_H

#include "ParticleConstraint.h"

class ParticleFitObject;

//  Class PConstraint:
/// Implements a constraint of the form pxfact*sum(px)+pyfact*sum(py)+pzfact*sum(pz)+efact*sum(E)=value
/**
 *
 * Author: Jenny List, Benno List
 * Last update: $Date: 2008/02/12 11:03:32 $
 *          by: $Author: blist $
 *
 */
class PConstraint : public ParticleConstraint {
  public:
    PConstraint (double pxfact_=1, double pyfact_=0, double pzfact_=0,
                 double efact_=0, double value_ = 0);
    virtual ~PConstraint();
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
  
  
    double pxfact;
    double pyfact;
    double pzfact;
    double efact;
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

#endif // __PCONSTRAINT_H
