////////////////////////////////////////////////////////////////
// Class PxConstraint
//
// Author: Jenny Boehme
// Last update: $Date: 2007-10-30 15:51:14 $
//          by: $Author: gaede $
// 
// Description: sum (px) = 0 constraint
//               
////////////////////////////////////////////////////////////////

#ifndef __PXCONSTRAINT_H
#define __PXCONSTRAINT_H

#include "ParticleConstraint.h"

class ParticleFitObject;

class PxConstraint : public ParticleConstraint {
  public:
    PxConstraint ();
    virtual ~PxConstraint();
    virtual double getValue() const;
    virtual void getDerivatives (int idim, double der[]) const;
    virtual void add1stDerivativesToMatrix(int idim, double *M) const;
    virtual void add2ndDerivativesToMatrix(int idim, double *M, double lambda) const;
  
  protected:

};

#endif // __PXCONSTRAINT_H
